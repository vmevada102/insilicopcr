#!/usr/bin/env perl
# insilico_pcr_parallel_modified.pl
# Parallel in-silico PCR across multiple genome FASTA files in a directory.
# Reads primer pairs from a CSV and writes one multifasta per primer pair plus a combined summary CSV.
#
# Improvements:
# - Processes multiple genome FASTA files from a directory (--fasta_dir)
# - Reads primer pairs from a CSV (--primers)
# - Optional k-mer prefilter (--kmer N) to quickly skip sequences that lack a short exact anchor
# - Fork-based worker pool (one worker per primer pair, limited by --workers)
# - Core-Perl only (no CPAN dependencies)
#
# Usage:
#  perl insilico_pcr_parallel_modified.pl --fasta_dir FASTA_DIR --primers primers.csv --out_dir OUT_DIR \
#       --max_mm 2 --min_len 50 --max_len 2000 --kmer 6 --workers 0
#
# Notes:
#  - --workers 0 => auto-detect number of CPUs (default)
#  - --kmer 0 => disable prefilter
#  - Primer CSV must have headers that contain 'forward' and 'reverse' (case-insensitive)
#  - Output:
#      OUT_DIR/<safe_pair_id>_amplicons.fasta
#      OUT_DIR/amplicons_summary.csv
#
use strict;
use warnings;
use Getopt::Long;
use POSIX ":sys_wait_h";
use File::Basename;
use File::Spec;
use Fcntl qw(:flock);

# -------------------------
# IUPAC map (hashrefs for fast lookup)
# -------------------------
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

# -------------------------
# Options and defaults
# -------------------------
my ($fasta_dir, $primers_csv, $out_dir, $max_mm, $min_len, $max_len, $workers, $kmer, $help);
$max_mm  = 2;
$min_len = 20;
$max_len = 20000;
$workers = 0;   # 0 => auto-detect => use all CPUs
$kmer    = 0;   # 0 => disable k-mer prefilter
GetOptions(
    'fasta_dir=s' => \$fasta_dir,
    'primers=s'   => \$primers_csv,
    'out_dir=s'   => \$out_dir,
    'max_mm=i'    => \$max_mm,
    'min_len=i'   => \$min_len,
    'max_len=i'   => \$max_len,
    'workers=i'   => \$workers,
    'kmer=i'      => \$kmer,
    'help|h'      => \$help,
) or usage();

usage() if $help;
for my $req ($fasta_dir, $primers_csv, $out_dir) { usage() unless defined $req; }

sub usage {
    print STDERR <<"USAGE";
Usage:
  perl $0 --fasta_dir FASTA_DIR --primers primers.csv --out_dir OUT_DIR [--max_mm 2] [--min_len 20] [--max_len 20000] [--workers 0] [--kmer 0]
Notes:
  --workers 0 => auto-detect (all CPUs). Set >0 to limit.
  --kmer N => use first N bases of forward primer as an exact k-mer prefilter (0 = disable).
Primer CSV must have headers containing 'forward' and 'reverse' (case-insensitive).
USAGE
    exit 1;
}

# -------------------------
# Utility functions
# -------------------------
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
    $s =~ tr/U/T/;
    return $s;
}

sub primer_allowed_array {
    my $p = shift;
    $p = uc($p);
    my @arr;
    foreach my $ch (split //, $p) {
        if (exists $IUPAC{$ch}) {
            push @arr, $IUPAC{$ch};
        } else {
            push @arr, $IUPAC{'N'};
        }
    }
    return \@arr;
}

# FASTA reader with callback for memory-efficiency
sub process_fasta_with_callback {
    my ($path, $cb) = @_;
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

# Sliding-window approximate matching with early-exit mismatch counting
# Returns arrayref of [start0, end0, mismatches]
sub find_approx_matches {
    my ($seq, $primer_arr, $max_mm) = @_;
    my $L = scalar @$primer_arr;
    my $N = length $seq;
    return [] if $L == 0 || $N < $L;
    my @hits;
    for (my $i = 0; $i <= $N - $L; $i++) {
        my $mism = 0;
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

# -------------------------
# Read primer CSV
# -------------------------
open my $pcfh, '<', $primers_csv or die "Can't open primers CSV $primers_csv: $!";
my $hdr = <$pcfh>;
die "Primer CSV appears empty" unless defined $hdr;
chomp $hdr; $hdr =~ s/\r//g;
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
    chomp $line; $line =~ s/\r//g;
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

# -------------------------
# Collect FASTA files (multiple genomes) from fasta_dir
# -------------------------
opendir my $dh, $fasta_dir or die "Can't open dir $fasta_dir: $!";
my @fasta_files = grep { -f File::Spec->catfile($fasta_dir, $_) && /\.(fa|fasta|fna|ffn)$/i } readdir $dh;
closedir $dh;
@fasta_files = map { File::Spec->catfile($fasta_dir, $_) } sort @fasta_files;
die "No FASTA files found in $fasta_dir" unless @fasta_files;

# -------------------------
# CPU detection & workers
# -------------------------
sub detect_cpus {
    if ($ENV{NSLOTS}) { return int($ENV{NSLOTS}); }
    my $n = `nproc 2>/dev/null`; chomp $n;
    return int($n) if $n =~ /^\d+$/ && $n > 0;
    $n = `getconf _NPROCESSORS_ONLN 2>/dev/null`; chomp $n;
    return int($n) if $n =~ /^\d+$/ && $n > 0;
    return 1;
}

my $cpu_count = detect_cpus();
$workers = $cpu_count if !$workers || $workers < 1;
print STDERR "Detected $cpu_count CPU(s). Using up to $workers worker(s).\n";

# -------------------------
# Prepare output directories
# -------------------------
mkdir $out_dir unless -d $out_dir;
my $tmpdir = $out_dir . "/.tmp_insilico";
mkdir $tmpdir unless -d $tmpdir;

# -------------------------
# Main: fork-limited worker pool, each child handles one primer pair across all genomes
# -------------------------
my @children;
my $active = 0;
my @pending = @primer_pairs;

while (@pending || $active > 0) {
    if (@pending && $active < $workers) {
        my $pair = shift @pending;
        my $pid = fork();
        if (!defined $pid) { die "fork failed: $!"; }
        if ($pid == 0) {
            # Child process: process this primer pair
            my $pair_id = $pair->{pair_id};
            my $fwd = $pair->{forward};
            my $rev = $pair->{reverse};
            my $revrc = revcomp($rev);
            my $fwd_arr = primer_allowed_array($fwd);
            my $rev_arr = primer_allowed_array($revrc);

            # optional k-mer prefilter (exact substring of the forward primer)
            my $kmer_seq = ($kmer && length($fwd) >= $kmer) ? substr($fwd,0,$kmer) : '';

            my @hits;
            foreach my $fasta_path (@fasta_files) {
                my $sample_name = fileparse($fasta_path, qr/\.[^.]*/);
                process_fasta_with_callback($fasta_path, sub {
                    my ($seq_id, $seq) = @_;
                    # k-mer prefilter: skip sequence if kmer specified and not present
                    if ($kmer_seq ne '' && index($seq, $kmer_seq) == -1) { return; }
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

            # write multifasta for this pair
            my $safe = safe_name($pair->{pair_id});
            my $outf = File::Spec->catfile($out_dir, "${safe}_amplicons.fasta");
            open my $ofh, '>', $outf or die "Can't write $outf: $!";
            my %counter;
            for my $hit (@hits) {
                my $sample = $hit->{sample_name};
                $counter{$sample}++;
                my $idx = $counter{$sample};
                my $header = $pair->{pair_id} . "_" . $sample . "_hit" . $idx . "|len=" . $hit->{amp_len}
                             . "|fwd=" . $hit->{fwd_pr} . "|rev=" . $hit->{rev_pr};
                print $ofh ">$header\n";
                my $seq = $hit->{amp_seq};
                $seq =~ s/(.{1,80})/$1\n/g;
                print $ofh $seq;
            }
            close $ofh;

            # write summary temp file for this pair
            my $tmp_summary = File::Spec->catfile($tmpdir, "${safe}_summary.tmp");
            open my $sfh, '>', $tmp_summary or die "Can't write $tmp_summary: $!";
            for my $h (@hits) {
                print $sfh join(",", map { defined $h->{$_} ? $h->{$_} : "" }
                    qw(pair_id fasta_file sample_name seq_id fwd_start fwd_end rev_start rev_end fwd_mm rev_mm amp_len fwd_pr rev_pr)
                ), "\n";
            }
            close $sfh;
            exit 0;
        } else {
            # parent
            push @children, $pid;
            $active++;
        }
    } else {
        # wait for child
        my $wpid = wait();
        if ($wpid > 0) { $active-- if $active > 0; } else { last; }
    }
}

# wait for any remaining children
1 while (wait() != -1);

# -------------------------
# Combine temporary summaries into final CSV
# -------------------------
my $summary_csv = File::Spec->catfile($out_dir, "amplicons_summary.csv");
open my $sumfh, '>', $summary_csv or die "Can't write $summary_csv: $!";
print $sumfh join(",", qw(pair_id fasta_file sample_name seq_id fwd_start fwd_end rev_start rev_end fwd_mm rev_mm amp_len fwd_pr rev_pr)), "\n";
opendir my $tdh, $tmpdir or die "Can't open tmpdir $tmpdir: $!";
for my $tmpfile (sort readdir $tdh) {
    next unless $tmpfile =~ /\.tmp$/;
    my $path = File::Spec->catfile($tmpdir, $tmpfile);
    open my $tfh, '<', $path or next;
    while (my $line = <$tfh>) { print $sumfh $line; }
    close $tfh;
    unlink $path;
}
closedir $tdh;
close $sumfh;
rmdir $tmpdir;

print STDERR "Done. Outputs in $out_dir\n";
exit 0;

