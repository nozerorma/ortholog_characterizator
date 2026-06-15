#!/usr/bin/env python3
"""
Reconstruct the BMGE-trimmed CDS alignment by projecting the original codon
alignment onto the columns BMGE selected.

BMGE -t CODON produces a correct protein alignment (-oaaf) but its nucleotide
output (-of) is truncated to protein-length characters — a known BMGE bug
present through at least the bioconda 1.12 package.

The previous approach walked the trimmed protein alignment and consumed one
codon per non-gap residue from each sequence's *gap-stripped* original CDS.
That is only correct when BMGE removes gap-only columns: as soon as BMGE trims
an informative column, the per-sequence residue→codon counting desynchronises
and every downstream codon is frame-shifted to the wrong codon. (Measured: a
gene trimmed from 399→310 codon columns produced only 123/310 correct codons.)

The correct reconstruction uses BMGE's own column-selection record: the TRACK
HTML (-oh output) lists exactly which 1-based codon columns were `selected:`.
Slicing each sequence of the *aligned* input CDS at those columns reproduces a
CDS alignment that translates back to the BMGE protein alignment exactly
(verified 310/310, 654/654, 1861/1861 on real genes).

Usage:
    backtranslate_from_prot.py --cds input_cds.fa --track track.html --out trim.fa
                               [--prot prot.fa]   # optional sanity check
"""

import argparse
import re
import sys


def read_fasta(path):
    header, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts)
                header = line[1:]
                parts = []
            elif line:
                parts.append(line)
    if header is not None:
        yield header, ''.join(parts)


def parse_selected_columns(track_html):
    """Return sorted list of 1-based codon columns BMGE kept (`selected:` rows)."""
    selected = []
    with open(track_html) as fh:
        for line in fh:
            plain = re.sub(r'<[^>]*>', ' ', line)
            tokens = plain.split()
            for i, tok in enumerate(tokens):
                if tok == 'selected:':
                    for t in tokens[i + 1:]:
                        if re.match(r'^\d+-\d+$', t):
                            a, b = t.split('-')
                            selected.extend(range(int(a), int(b) + 1))
                        elif t.isdigit():
                            selected.append(int(t))
                    break
    return sorted(set(selected))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cds',   required=True, help='Original ALIGNED CDS (BMGE input)')
    ap.add_argument('--track', required=True, help='BMGE TRACK HTML (-oh output)')
    ap.add_argument('--out',   required=True, help='Output trimmed CDS alignment')
    ap.add_argument('--prot',  default=None,
                    help='Optional BMGE protein alignment (-oaaf) for a sanity check')
    args = ap.parse_args()

    selected = parse_selected_columns(args.track)
    if not selected:
        print(f'ERROR: no selected columns parsed from {args.track}', file=sys.stderr)
        sys.exit(1)

    # 0-based nucleotide slice positions for each kept codon column.
    starts = [(c - 1) * 3 for c in selected]

    records = list(read_fasta(args.cds))
    if not records:
        print(f'ERROR: no sequences in {args.cds}', file=sys.stderr)
        sys.exit(1)

    aln_len = len(records[0][1])
    max_needed = starts[-1] + 3 if starts else 0
    if max_needed > aln_len:
        print(f'ERROR: track references codon column beyond alignment '
              f'({max_needed} nt needed, alignment is {aln_len} nt) — '
              f'track/CDS mismatch', file=sys.stderr)
        sys.exit(1)

    with open(args.out, 'w') as fh:
        for header, seq in records:
            trimmed = ''.join(seq[s:s + 3] for s in starts)
            fh.write(f'>{header}\n{trimmed}\n')

    # Optional sanity check: translated trimmed CDS == protein alignment.
    if args.prot:
        try:
            from Bio.Seq import Seq
            prot = {h.split()[0]: s for h, s in read_fasta(args.prot)}
            cds_out = {h.split()[0]: s for h, s in read_fasta(args.out)}
            mism = 0
            for sid, p in prot.items():
                c = cds_out.get(sid, '')
                tr = ''.join('-' if c[i:i + 3] == '---'
                             else str(Seq(c[i:i + 3]).translate())
                             for i in range(0, len(c), 3))
                mism += sum(1 for a, b in zip(p, tr) if a != b)
            if mism:
                print(f'WARN: {mism} residue mismatches between trimmed CDS and '
                      f'protein alignment in {args.out}', file=sys.stderr)
        except ImportError:
            pass


if __name__ == '__main__':
    main()
