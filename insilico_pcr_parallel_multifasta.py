#!/usr/bin/env python3
"""
insilico_pcr_parallel_multifasta.py

In-silico PCR over all FASTA files in a directory using primer pairs from a CSV.

Outputs:
 - <out_dir>/<safe_pair_id>_amplicons.fasta  (multi-FASTA, one file per primer pair)
 - <out_dir>/amplicons_summary.csv            (summary table of all amplicons)

Header format for each amplicon record:
    <pair_id>_<sample_name>_hit<N>|len=<amplicon_length>|fwd=<forward_primer>|rev=<reverse_primer>

This script uses multiprocessing (ProcessPoolExecutor) and by default will spawn as many
workers as CPUs available (os.cpu_count()).
"""

from pathlib import Path
import argparse
import csv
import re
from typing import List, Tuple, Dict
from concurrent.futures import ProcessPoolExecutor, as_completed
import os

# -------------------------
# FASTA reader
# -------------------------
def read_fasta(path: Path):
    """Yield (seq_id, sequence) for each record in a FASTA file (uppercased)."""
    with open(path, 'r') as fh:
        seq_id = None
        seq_parts = []
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if seq_id is not None:
                    yield seq_id, ''.join(seq_parts).upper()
                seq_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if seq_id is not None:
            yield seq_id, ''.join(seq_parts).upper()

# -------------------------
# IUPAC & matching helpers
# -------------------------
IUPAC = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'}, 'U': {'T'},
    'R': {'A','G'}, 'Y': {'C','T'}, 'S': {'G','C'}, 'W': {'A','T'},
    'K': {'G','T'}, 'M': {'A','C'}, 'B': {'C','G','T'}, 'D': {'A','G','T'},
    'H': {'A','C','T'}, 'V': {'A','C','G'}, 'N': {'A','C','G','T'}
}

def revcomp(seq: str) -> str:
    comp = {'A':'T','T':'A','G':'C','C':'G','U':'A','N':'N'}
    return ''.join(comp.get(b, 'N') for b in seq[::-1].upper())

def primer_allowed_set(primer: str):
    return [IUPAC.get(ch, {'A','C','G','T'}) for ch in primer.upper()]

def count_mismatches(seg: str, allowed_sets: List[set]) -> int:
    mism = 0
    for b, allow in zip(seg, allowed_sets):
        if b not in allow:
            mism += 1
    return mism

def find_approx_matches(seq: str, primer: str, max_mismatch: int) -> List[Tuple[int,int,int]]:
    L = len(primer)
    if L == 0 or len(seq) < L:
        return []
    allowed = primer_allowed_set(primer)
    hits = []
    # naive sliding window (fast enough for bacterial-size sequences)
    for i in range(0, len(seq) - L + 1):
        seg = seq[i:i+L]
        mm = count_mismatches(seg, allowed)
        if mm <= max_mismatch:
            hits.append((i, i+L, mm))
    return hits

# -------------------------
# CSV parser for primer pairs
# -------------------------
def parse_primer_csv(path: Path):
    """Return list of dicts: {'pair_id','forward','reverse'}"""
    pairs = []
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError("Primer CSV is empty or malformed.")
        headers = [h.lower() for h in reader.fieldnames]
        fwd_col = None
        rev_col = None
        id_col = None
        for h in headers:
            if any(k in h for k in ['forward','fwd','left','primer1','primer_fwd']):
                fwd_col = h
            if any(k in h for k in ['reverse','rev','right','primer2','primer_rev']):
                rev_col = h
            if any(k in h for k in ['pair','id','pair_id','name']):
                id_col = h
        if fwd_col is None or rev_col is None:
            raise ValueError("Could not detect forward/reverse columns in primer CSV. "
                             "Ensure columns like 'forward' and 'reverse' (or 'fwd','rev') exist.")
        for idx, row in enumerate(reader, start=1):
            fwd = row.get(fwd_col, '').strip()
            rev = row.get(rev_col, '').strip()
            if not fwd or not rev:
                continue
            pid = row.get(id_col, '').strip() if id_col else ''
            if not pid:
                pid = f"pair_{idx}"
            pairs.append({'pair_id': pid, 'forward': fwd.upper(), 'reverse': rev.upper()})
    return pairs

# -------------------------
# Utilities
# -------------------------
def safe_name(s: str) -> str:
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', s)

def wrap_seq(s: str, width: int = 80) -> str:
    return '\n'.join(s[i:i+width] for i in range(0, len(s), width))

# -------------------------
# Worker: process a single primer pair across all FASTA files
# -------------------------
def process_pair(pair: Dict, fasta_paths: List[str], max_mismatch: int, min_len: int, max_len: int):
    """
    Runs in a worker process.
    Returns: (pair_id, list_of_amplicons)
    Each amplicon is dict with keys:
      sample_name, seq_id, fwd_start, fwd_end, rev_start, rev_end, fwd_mm, rev_mm, amplicon_seq, fwd_primer, rev_primer
    """
    pair_id = pair['pair_id']
    fwd_pr = pair['forward']
    rev_pr = pair['reverse']
    revrc = revcomp(rev_pr)

    result_amplicons = []

    for fasta_path in fasta_paths:
        fasta_p = Path(fasta_path)
        sample_name = fasta_p.stem
        try:
            for seq_id, seq in read_fasta(fasta_p):
                seq = seq.upper()
                f_hits = find_approx_matches(seq, fwd_pr, max_mismatch)
                r_hits = find_approx_matches(seq, revrc, max_mismatch)
                for fstart0, fend0, fmism in f_hits:
                    for rstart0, rend0, rmism in r_hits:
                        # require forward bind left-of reverse bind on + strand
                        if fstart0 <= rstart0:
                            amp_len = rend0 - fstart0
                            if amp_len < min_len or amp_len > max_len:
                                continue
                            amplicon_seq = seq[fstart0:rend0]
                            result_amplicons.append({
                                'pair_id': pair_id,
                                'sample_name': sample_name,
                                'fasta_file': fasta_p.name,
                                'seq_id': seq_id,
                                'fwd_start_1based': fstart0 + 1,
                                'fwd_end_1based': fend0,
                                'rev_start_1based': rstart0 + 1,
                                'rev_end_1based': rend0,
                                'fwd_mismatches': fmism,
                                'rev_mismatches': rmism,
                                'amplicon_length': amp_len,
                                'amplicon_sequence': amplicon_seq,
                                'fwd_primer': fwd_pr,
                                'rev_primer': rev_pr
                            })
        except Exception as e:
            # don't crash worker on one bad fasta; return what we have plus note
            result_amplicons.append({
                'pair_id': pair_id,
                'error': f"Error reading {fasta_path}: {e}"
            })
    return (pair_id, result_amplicons)

# -------------------------
# Main runner (master process)
# -------------------------
def run_parallel(fasta_dir: Path, primer_csv: Path, out_dir: Path,
                 max_mismatch: int = 2, min_len: int = 20, max_len: int = 20000,
                 workers: int = None):
    primers = parse_primer_csv(primer_csv)
    fasta_paths = sorted([str(p) for p in fasta_dir.iterdir() if p.is_file() and p.suffix.lower() in ('.fa', '.fasta', '.fna', '.ffn')])
    if not fasta_paths:
        raise SystemExit(f"No FASTA files found in {fasta_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)

    # default workers -> all CPUs
    if workers is None:
        workers = os.cpu_count() or 1

    # Submit one task per primer pair
    tasks = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(process_pair, pair, fasta_paths, max_mismatch, min_len, max_len): pair for pair in primers}
        for fut in as_completed(futures):
            pair = futures[fut]
            try:
                pair_id, amplicons = fut.result()
            except Exception as e:
                print(f"[ERROR] pair {pair['pair_id']} failed: {e}")
                continue

            # Filter and collect only valid amplicon dicts
            valid_amplicons = [a for a in amplicons if 'amplicon_sequence' in a]

            # Write per-pair multifasta
            safe = safe_name(pair_id)
            out_fa = out_dir / f"{safe}_amplicons.fasta"
            with open(out_fa, 'w') as fh:
                # group hits by sample_name to number them per sample (optional)
                counters = {}
                for hit in valid_amplicons:
                    sample = hit['sample_name']
                    counters.setdefault(sample, 0)
                    counters[sample] += 1
                    hit_idx = counters[sample]

                    # header: <pair>_<sample>_hit<N>|len=<amp_len>|fwd=<forward>|rev=<reverse>
                    header = f"{pair_id}_{sample}_hit{hit_idx}|len={hit['amplicon_length']}|fwd={hit['fwd_primer']}|rev={hit['rev_primer']}"
                    fh.write(f">{header}\n")
                    fh.write(wrap_seq(hit['amplicon_sequence']) + "\n")
            print(f"Wrote {len(valid_amplicons)} amplicons -> {out_fa}")

            # Append to a combined summary CSV data structure
            # We'll write summary after loop to ensure ordering / single file.
            tasks.append((pair_id, valid_amplicons))

    # After all pairs processed, produce a single summary CSV
    summary_rows = []
    for pair_id, amps in tasks:
        for a in amps:
            summary_rows.append({
                'pair_id': pair_id,
                'fasta_file': a.get('fasta_file', ''),
                'sample_name': a.get('sample_name', ''),
                'seq_id': a.get('seq_id', ''),
                'fwd_start_1based': a.get('fwd_start_1based', ''),
                'fwd_end_1based': a.get('fwd_end_1based', ''),
                'rev_start_1based': a.get('rev_start_1based', ''),
                'rev_end_1based': a.get('rev_end_1based', ''),
                'fwd_mismatches': a.get('fwd_mismatches', ''),
                'rev_mismatches': a.get('rev_mismatches', ''),
                'amplicon_length': a.get('amplicon_length', ''),
                'fwd_primer': a.get('fwd_primer', ''),
                'rev_primer': a.get('rev_primer', '')
            })

    summary_csv = out_dir / "amplicons_summary.csv"
    fieldnames = ['pair_id','fasta_file','sample_name','seq_id',
                  'fwd_start_1based','fwd_end_1based','rev_start_1based','rev_end_1based',
                  'fwd_mismatches','rev_mismatches','amplicon_length','fwd_primer','rev_primer']
    with open(summary_csv, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in summary_rows:
            writer.writerow(r)
    print(f"Wrote summary -> {summary_csv}")

# -------------------------
# CLI
# -------------------------
def main():
    ap = argparse.ArgumentParser(description="Parallel in-silico PCR -> multifasta per primer pair")
    ap.add_argument('--fasta_dir', required=True, help='Directory with FASTA files (.fa/.fasta/.fna)')
    ap.add_argument('--primers', required=True, help='CSV file with primer pairs (columns: forward, reverse, optional id)')
    ap.add_argument('--out_dir', required=True, help='Directory to write per-pair multifasta files and summary CSV')
    ap.add_argument('--max_mismatch', type=int, default=2, help='Maximum mismatches allowed per primer (default 2)')
    ap.add_argument('--min_len', type=int, default=20, help='Minimum amplicon length (default 20)')
    ap.add_argument('--max_len', type=int, default=20000, help='Maximum amplicon length (default 20000)')
    ap.add_argument('--workers', type=int, default=None, help='Number of worker processes to use (default = number of CPUs)')
    args = ap.parse_args()

    run_parallel(Path(args.fasta_dir), Path(args.primers), Path(args.out_dir),
                 max_mismatch=args.max_mismatch, min_len=args.min_len, max_len=args.max_len,
                 workers=args.workers)

if __name__ == '__main__':
    main()

