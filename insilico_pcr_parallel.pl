#!/usr/bin/env perl
## insilico_pcr_parallel.pl
## Fast (Perl) parallel in-silico PCR: one worker per primer pair (fork)
## Outputs:
##   <outdir>/<safe_pair_id>_amplicons.fasta
##   <outdir>/amplicons_summary.csv
##
## Usage:
##   perl insilico_pcr_parallel.pl --fasta_dir FASTA --primers primers.csv --out_dir OUT --max_mm 2 --min_len 50 --max_len 2000 --workers 0
##
## Notes:
##  - primer CSV must contain headers with 'forward' and 'reverse' (case-insensitive).
##  - each worker processes a single primer pair across all FASTA files.
##  - default workers = number of CPUs (detects via nproc / getconf / environment).
##  - all code uses core Perl only.

use strict;
use warnings;
use Getopt::Long;
use POSIX ":sys_wait_h";
use File::Basename;
use File::Spec;
use Fcntl qw(:flock);

# -------------- iupac map -> set of allowed bases --------------
my %IUPAC = (
    A => { map { $_ => 1 } qw(A) },
    C => { map { $_ => 1 } qw(C) },
    G => { map { $_ => 1 } qw(G) },
    T => { map { $_ => 1 } qw(T) },
    U => { map { $_ => 1 } qw(T) },
    R => { map { $_ => 1 } qw(A G) },
    Y => { map { $_ => 1 } qw(C T) },
    S => { map { $_ => 1 } qw(G C) },
    W => { map { $_ => 1 } qw(A T) },
    K => { map { $_ => 1 } qw(G T) },
    M => { map { $_ => 1 } qw(A C) },
    B => { map { $_ => 1 } qw(C G T) },
    D => { map { $_ => 1 } qw(A G T) },
    H => { map { $_ => 1 } qw(A C T) },
    V => { map { $_ => 1 } qw(A C G) },
    N => { map { $_ => 1 } qw(A C G T) },
);

# -------------- options and defaults --------------
my ($fasta_dir, $primers_csv, $out_dir, $max_mm, $min_len, $max_len, $workers, $help);
$max_mm  = 2;
$min_len = 20;
$max_len = 20000;
$workers = 0;   # 0 => auto-detect -> use number of CPUs
GetOptions(
    'fasta_dir=s' => \$fasta_dir,
    'primers=s'   => \$primers_csv,
    'out_dir=s'   => \$out_dir,
    'max_mm=i'    => \$max_mm,
    'min_len=i'   => \$min_len,
    'max_len=i'   => \$max_len,
    'workers=i'   => \$workers,
    'help|h'      => \$help,
) or usage();

usage() if $help;
for my $req ($fasta_dir, $primers_csv, $out_dir) {
    usage() unless defined $req;
}

# -------------- helpers --------------
sub usage {
    print STDERR <<"USAGE";
Usage:
  perl $0 --fasta_dir FASTA --primers primers.csv --out_dir OUT [--max_mm 2] [--min_len 20] [--max_len 20000] [--workers 0]
Notes:
  --workers 0 => auto-detect (all CPUs). Set >0 to limit.
Primer CSV must have headers containing 'forward' and 'reverse' (case-insensitive).
USAGE
    exit 1;
}

sub safe_name {
    my $s = shift;
    $s =~ s/[^A-Za-z0-9_.-]+/_/g;
    return $s;
}

sub revcomp {
    my $s = shift;
    $s = uc($s);
    $s =~ tr/ACGTU/TGCAA/;
    $s = reverse $s;
    # above tr: U->A then reverse; we want U -> A converted to T earlier; handle missing mapping
    $s =~ tr/U/T/;  # in case U remained
    return $s;
}

# Convert primer string into arrayref of allowed bases hashrefs (fast lookup)
sub primer_allowed_array {
    my $p = shift;
    $p = uc($p);
    my @arr;
    foreach my $ch (split //, $p) {
        if (exists $IUPAC{$ch}) {
            push @arr, $IUPAC{$ch};
        } else {
            # unknown -> allow all
            push @arr, $IUPAC{'N'};
        }
    }
    return \@arr;
}

# FASTA reader that yields (id, seq) for each record
sub read_fasta_file {
    my ($path) = @_;
    open my $fh, '<', $path or die "Can't open FASTA $path: $!";
    my $id;
    local $/ = "\n";
    my $seq = '';
    while (my $line = <$fh>) {
        chomp $line;
        next if $line eq '';
        if ($line =~ /^>(.*)/) {
            if (defined $id) {
                yield_fasta($id, $seq);
            }
            $id = (split /\s+/, $1)[0];
            $seq = '';
        } else {
            $seq .= uc($line);
        }
    }
    if (defined $id) {
        yield_fasta($id, $seq);
    }
    close $fh;
}

# Because we cannot yield easily, implement callback based reader below for speed
sub process_fasta_with_callback {
    my ($path, $cb) = @_;  # cb->(seq_id, seq)
    open my $fh, '<', $path or die "Can't open FASTA $path: $!";
    my $id;
    local $/ = "\n";
    my $seq = '';
    while (my $line = <$fh>) {
        chomp $line;
        next if $line eq '';
        if ($line =~ /^>(.*)/) {
            if (defined $id) {
                $cb->($id, $seq);
            }
            $id = (split /\s+/, $1)[0];
            $seq = '';
        } else {
            $seq .= uc($line);
        }
    }
    if (defined $id) {
        $cb->($id, $seq);
    }
    close $fh;
}

# approximate match by sliding window with early-exit mismatches
# returns arrayref of [start0, end0, mismatches]
sub find_approx_matches {
    my ($seq, $primer_arr, $max_mm) = @_;
    my $L = scalar @$primer_arr;
    my $N = length $seq;
    return [] if $L == 0 || $N < $L;
    my @hits;
    # slide
    for (my $i = 0; $i <= $N - $L; $i++) {
        my $mism = 0;
        # inner loop: compare and early exit
        for (my $k = 0; $k < $L; $k++) {
            my $base = substr($seq, $i+$k, 1);
            unless (exists $primer_arr->[$k]{$base}) {
                $mism++;
                last if $mism > $max_mm;
            }
        }
        push @hits, [$i, $i+$L, $mism] if $mism <= $max_mm;
    }
    return \@hits;
}

# -------------- read primers CSV --------------
open my $pcfh, '<', $primers_csv or die "Can't open primers CSV $primers_csv: $!";
my $hdr = <$pcfh>;
chomp $hdr;
$hdr =~ s/\r//g;
my @cols = map { lc $_ } split /,/, $hdr;
my ($fwd_idx, $rev_idx, $id_idx);
for my $i (0..$#cols) {
    $fwd_idx = $i if $cols[$i] =~ /forward|fwd|left|primer1/;
    $rev_idx = $i if $cols[$i] =~ /reverse|rev|right|primer2/;
    $id_idx  = $i if $cols[$i] =~ /pair|id|pair_id|name/;
}
die "Could not detect forward/reverse columns in primer CSV header." unless defined $fwd_idx && defined $rev_idx;

my @primer_pairs;
while (my $line = <$pcfh>) {
    chomp $line;
    $line =~ s/\r//g;
    next if $line =~ /^\s*$/;
    my @f = split /,/, $line, scalar(@cols);
    my $fwd = $f[$fwd_idx] // '';
    my $rev = $f[$rev_idx] // '';
    next unless $fwd && $rev;
    my $pid = $id_idx ? ($f[$id_idx] // '') : '';
    $pid = 'pair_' . (scalar(@primer_pairs)+1) if $pid eq '';
    push @primer_pairs, {
        pair_id => $pid,
        forward => uc($fwd),
        reverse => uc($rev),
    };
}
close $pcfh;

die "No primer pairs found in $primers_csv" unless @primer_pairs;

# -------------- gather fasta files --------------
opendir my $dh, $fasta_dir or die "Can't open dir $fasta_dir: $!";
my @fasta_files = grep { -f File::Spec->catfile($fasta_dir, $_) && /\.(fa|fasta|fna|ffn)$/i } readdir $dh;
closedir $dh;
@fasta_files = map { File::Spec->catfile($fasta_dir, $_) } sort @fasta_files;
die "No FASTA files found in $fasta_dir" unless @fasta_files;

# -------------- workers --------------
sub detect_cpus {
    # try environment
    if ($ENV{NSLOTS}) { return int($ENV{NSLOTS}); }
    # try nproc
    my $n = `nproc 2>/dev/null`; chomp $n;
    return int($n) if $n =~ /^\d+$/ && $n > 0;
    # try getconf
    $n = `getconf _NPROCESSORS_ONLN 2>/dev/null`; chomp $n;
    return int($n) if $n =~ /^\d+$/ && $n > 0;
    # fallback
    return 1;
}

my $cpu_count = detect_cpus();
$workers = $cpu_count if !$workers || $workers < 1;

print STDERR "Detected $cpu_count CPU(s). Using up to $workers worker(s).\n";

# -------------- main: fork one process per primer pair (but limit concurrent forks by $workers) --------------
mkdir $out_dir unless -d $out_dir;

my @children;
my $active = 0;
my @pending = @primer_pairs;

# store summary rows from children in temp files; main will combine
my $tmpdir = $out_dir . "/.tmp_insilico";
mkdir $tmpdir unless -d $tmpdir;

while (@pending || $active > 0) {
    if (@pending && $active < $workers) {
        my $pair = shift @pending;
        my $pid = fork();
        if (!defined $pid) {
            die "fork failed: $!";
        }
        if ($pid == 0) {
            # child process: process this primer pair
            my $pair_id = $pair->{pair_id};
            my $fwd = $pair->{forward};
            my $rev = $pair->{reverse};
            my $revrc = revcomp($rev);
            my $fwd_arr = primer_allowed_array($fwd);
            my $rev_arr = primer_allowed_array($revrc);

            my @hits;  # each = hashref for a found amplicon

            foreach my $fasta_path (@fasta_files) {
                my $sample_name = fileparse($fasta_path, qr/\.[^.]*/);
                # read FASTA and check
                process_fasta_with_callback($fasta_path, sub {
                    my ($seq_id, $seq) = @_;
                    # find forward hits and reverse hits
                    my $f_hits = find_approx_matches($seq, $fwd_arr, $max_mm);
                    my $r_hits = find_approx_matches($seq, $rev_arr, $max_mm);
                    return unless @$f_hits && @$r_hits;
                    foreach my $fhit (@$f_hits) {
                        foreach my $rhit (@$r_hits) {
                            my ($fs, $fe, $fmm) = @$fhit;
                            my ($rs, $re, $rmm) = @$rhit;
                            if ($fs <= $rs) {
                                my $amp_len = $re - $fs;
                                next if $amp_len < $min_len || $amp_len > $max_len;
                                my $amp_seq = substr($seq, $fs, $amp_len);
                                push @hits, {
                                    pair_id => $pair_id,
                                    sample_name => $sample_name,
                                    fasta_file  => $fasta_path,
                                    seq_id => $seq_id,
                                    fwd_start => $fs+1,
                                    fwd_end   => $fe,
                                    rev_start => $rs+1,
                                    rev_end   => $re,
                                    fwd_mm    => $fmm,
                                    rev_mm    => $rmm,
                                    amp_len   => $amp_len,
                                    amp_seq   => $amp_seq,
                                    fwd_pr    => $fwd,
                                    rev_pr    => $rev,
                                };
                            }
                        }
                    }
                });
            }

            # write multi-fasta output for this pair
            my $safe = safe_name($pair->{pair_id});
            my $outf = File::Spec->catfile($out_dir, "${safe}_amplicons.fasta");
            open my $ofh, '>', $outf or die "Can't write $outf: $!";
            my %counter;
            for my $hit (@hits) {
                my $sample = $hit->{sample_name};
                $counter{$sample}++;
                my $idx = $counter{$sample};
                # header: <pair>_<sample>_hit<N>|len=<len>|fwd=<fwd>|rev=<rev>
                my $header = $pair->{pair_id} . "_" . $sample . "_hit" . $idx . "|len=" . $hit->{amp_len}
                             . "|fwd=" . $hit->{fwd_pr} . "|rev=" . $hit->{rev_pr};
                print $ofh ">$header\n";
                # wrap sequence to 80 chars
                my $seq = $hit->{amp_seq};
                $seq =~ s/(.{1,80})/$1\n/g;
                print $ofh $seq;
            }
            close $ofh;

            # write summary temp file
            my $tmp_summary = File::Spec->catfile($tmpdir, "${safe}_summary.tmp");
            open my $sfh, '>', $tmp_summary or die "Can't write $tmp_summary: $!";
            for my $h (@hits) {
                print $sfh join(",", map { defined $h->{$_} ? $h->{$_} : "" } 
                    qw(pair_id fasta_file sample_name seq_id fwd_start fwd_end rev_start rev_end fwd_mm rev_mm amp_len fwd_pr rev_pr)
                ), "\n";
            }
            close $sfh;

            exit 0;  # child done
        } else {
            # parent
            push @children, $pid;
            $active++;
        }
    } else {
        # wait for any child
        my $wpid = wait();
        if ($wpid > 0) {
            $active-- if $active > 0;
        } else {
            # no children left
            last;
        }
    }
}

# wait for any remaining children
1 while (wait() != -1);

# -------------- combine summaries --------------
my $summary_csv = File::Spec->catfile($out_dir, "amplicons_summary.csv");
open my $sumfh, '>', $summary_csv or die "Can't write $summary_csv: $!";
print $sumfh join(",", qw(pair_id fasta_file sample_name seq_id fwd_start fwd_end rev_start rev_end fwd_mm rev_mm amp_len fwd_pr rev_pr)), "\n";
opendir my $tdh, $tmpdir or die "Can't open tmpdir $tmpdir: $!";
for my $tmpfile (sort readdir $tdh) {
    next unless $tmpfile =~ /\.tmp$/;
    my $path = File::Spec->catfile($tmpdir, $tmpfile);
    open my $tfh, '<', $path or next;
    while (my $line = <$tfh>) {
        print $sumfh $line;
    }
    close $tfh;
    unlink $path;
}
closedir $tdh;
close $sumfh;

# remove tmpdir
rmdir $tmpdir;



## commit

print STDERR "Done. Outputs in $out_dir\n";
exit 0;

