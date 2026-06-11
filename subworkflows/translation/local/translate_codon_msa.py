#!/usr/bin/env python3
"""
Translate a codon CDS multiple-sequence alignment to a protein MSA.

Steps:
  1. Filter sequences with excessive ambiguity (N + gap fraction > threshold)
  2. Translate each retained sequence codon by codon
  3. Remove protein columns where gap/X fraction > gap_threshold
  4. Write PROT_RAW/<gene>.fa

Position mapping is trivial: protein column i = original codon column i.
No tracking file is needed.
"""

import argparse
import os
import sys

# ---------------------------------------------------------------------------
# Standard genetic code (NCBI table 1)
# ---------------------------------------------------------------------------
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

_GAP_CHARS = set('-X')


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def read_fasta(path):
    seqs = []
    header, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    seqs.append((header, ''.join(buf)))
                header = line[1:]
                buf = []
            else:
                buf.append(line.upper())
    if header is not None:
        seqs.append((header, ''.join(buf)))
    return seqs


def write_fasta(path, seqs, wrap=60):
    with open(path, 'w') as fh:
        for header, seq in seqs:
            fh.write(f'>{header}\n')
            for i in range(0, len(seq), wrap):
                fh.write(seq[i:i + wrap] + '\n')


# ---------------------------------------------------------------------------
# Sequence quality filter (replaces chi2 for this branch)
# ---------------------------------------------------------------------------

def filter_sequences(seqs, seq_ambiguity_threshold=0.7, min_seqs=4, log_fh=None):
    """
    Skip sequences where (N_count + gap_count) / total_nt > seq_ambiguity_threshold.
    Operates on the raw codon sequence string before translation.
    Returns (kept_seqs, n_skipped).
    Raises ValueError with a message if kept < min_seqs.
    """
    kept, n_skipped = [], 0
    for header, seq in seqs:
        total = len(seq)
        if total == 0:
            n_skipped += 1
            if log_fh:
                log_fh.write(f'{header}\t0.0\tempty_sequence\n')
            continue
        amb = seq.count('N') + seq.count('-')
        frac = amb / total
        if frac > seq_ambiguity_threshold:
            n_skipped += 1
            if log_fh:
                log_fh.write(f'{header}\t{frac:.4f}\tambiguity_above_threshold\n')
        else:
            kept.append((header, seq))

    if len(kept) < min_seqs:
        raise ValueError(
            f'Only {len(kept)} sequences remain after ambiguity filter '
            f'(threshold={seq_ambiguity_threshold}, min_seqs={min_seqs})'
        )
    return kept, n_skipped


# ---------------------------------------------------------------------------
# Translation
# ---------------------------------------------------------------------------

def translate_sequence(codon_seq):
    """
    Translate a codon-aligned sequence to protein.
      '---'  (all-gap codon)    → '-'
      partial gap or any 'N'    → 'X'
      stop codon                → 'X' (warn once via return value)
      valid codon               → amino acid from CODON_TABLE
    Returns (protein_str, found_stop: bool).
    Length = len(codon_seq) // 3.
    """
    aa_chars = []
    found_stop = False
    n = len(codon_seq)
    for i in range(0, n - 2, 3):
        codon = codon_seq[i:i + 3]
        if codon == '---':
            aa_chars.append('-')
        elif '-' in codon or 'N' in codon:
            aa_chars.append('X')
        else:
            aa = CODON_TABLE.get(codon, 'X')
            if aa == '*':
                aa_chars.append('X')
                found_stop = True
            else:
                aa_chars.append(aa)
    return ''.join(aa_chars), found_stop


# ---------------------------------------------------------------------------
# Column-level gap filter
# ---------------------------------------------------------------------------

def remove_high_gap_columns(prot_seqs, gap_threshold=0.9):
    """
    Remove protein alignment columns where (gap + X) fraction > gap_threshold.
    Returns new list of (header, trimmed_seq) with only retained columns.
    """
    if not prot_seqs:
        return prot_seqs
    n_seqs = len(prot_seqs)
    n_cols = len(prot_seqs[0][1])

    keep = []
    for col in range(n_cols):
        n_gap = sum(1 for _, seq in prot_seqs if seq[col] in _GAP_CHARS)
        if n_gap / n_seqs <= gap_threshold:
            keep.append(col)

    trimmed = [(h, ''.join(seq[c] for c in keep)) for h, seq in prot_seqs]
    return trimmed, len(keep)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Translate codon CDS MSA to protein MSA (raw variability branch).')
    parser.add_argument('--input', required=True,
                        help='Codon MSA FASTA (e.g. DEDUP/<gene>.fa)')
    parser.add_argument('--out_dir', required=True,
                        help='Base output dir; protein written to <out_dir>/PROT_RAW/<gene>.fa')
    parser.add_argument('--gap_threshold', type=float, default=0.9,
                        help='Column-level: drop columns with (gap+X)/n_seqs > threshold (default 0.9)')
    parser.add_argument('--seq_ambiguity_threshold', type=float, default=0.7,
                        help='Sequence-level: drop sequences with (N+gap)/len > threshold (default 0.7)')
    parser.add_argument('--min_seqs', type=int, default=4,
                        help='Skip gene if fewer than this many sequences remain (default 4)')
    parser.add_argument('--log_file', default=None,
                        help='Optional TSV log of filtered sequences')
    args = parser.parse_args()

    import re as _re
    gene = _re.sub(r'(\.(dedup|filtered))?\.fa$', '', os.path.basename(args.input))
    prot_raw_dir = os.path.join(args.out_dir, 'PROT_RAW')
    os.makedirs(prot_raw_dir, exist_ok=True)
    out_path = os.path.join(prot_raw_dir, f'{gene}.fa')

    seqs = read_fasta(args.input)
    n_input = len(seqs)
    if n_input == 0:
        print(f'{gene}\t0\t0\t0\t0', flush=True)
        return

    log_fh = None
    if args.log_file:
        os.makedirs(os.path.dirname(args.log_file) or '.', exist_ok=True)
        log_fh = open(args.log_file, 'w')
        log_fh.write('header\tambiguity_fraction\treason\n')

    try:
        kept_seqs, n_skipped = filter_sequences(
            seqs,
            seq_ambiguity_threshold=args.seq_ambiguity_threshold,
            min_seqs=args.min_seqs,
            log_fh=log_fh,
        )
    except ValueError as e:
        print(f'SKIP {gene}: {e}', file=sys.stderr)
        print(f'{gene}\t{n_input}\t{n_input - len([s for s in seqs])}\t0\t0', flush=True)
        if log_fh:
            log_fh.close()
        return
    finally:
        if log_fh:
            log_fh.close()

    n_filtered = n_input - n_skipped
    n_cols_raw = len(kept_seqs[0][1]) // 3 if kept_seqs else 0

    prot_seqs = []
    warned_stop = False
    for header, codon_seq in kept_seqs:
        prot_str, found_stop = translate_sequence(codon_seq)
        if found_stop and not warned_stop:
            print(f'WARN {gene}: stop codon(s) in mid-CDS; set to X', file=sys.stderr)
            warned_stop = True
        prot_seqs.append((header, prot_str))

    prot_seqs, n_cols_after = remove_high_gap_columns(prot_seqs, args.gap_threshold)

    if not prot_seqs or n_cols_after == 0:
        print(f'WARN {gene}: no columns remain after gap filter', file=sys.stderr)
        print(f'{gene}\t{n_input}\t{n_filtered}\t{n_cols_raw}\t0', flush=True)
        return

    write_fasta(out_path, prot_seqs)
    print(f'{gene}\t{n_input}\t{n_filtered}\t{n_cols_raw}\t{n_cols_after}', flush=True)


if __name__ == '__main__':
    main()
