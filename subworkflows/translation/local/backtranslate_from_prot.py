#!/usr/bin/env python3
"""
Back-translate a BMGE protein alignment to a CDS alignment.

BMGE -t CODON produces a correct protein alignment (-oaaf) but the
nucleotide output (-of) is truncated to protein-length characters — a
known BMGE bug present through at least the bioconda 1.12 package.

This script reconstructs the proper CDS alignment by walking the
protein alignment and consuming one codon per non-gap amino acid from
each sequence's raw (gap-stripped) original CDS.

Usage:
    backtranslate_from_prot.py --prot prot.fa --cds input_cds.fa --out trim.fa
"""

import argparse
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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--prot', required=True, help='BMGE protein alignment (-oaaf output)')
    ap.add_argument('--cds',  required=True, help='Original CDS alignment (BMGE input)')
    ap.add_argument('--out',  required=True, help='Output back-translated CDS alignment')
    args = ap.parse_args()

    # Strip gaps from the original CDS to get each sequence's raw codons
    raw_cds = {}
    for header, seq in read_fasta(args.cds):
        sid = header.split()[0]
        raw_cds[sid] = seq.replace('-', '').replace('.', '').upper()

    with open(args.out, 'w') as fh:
        for header, prot in read_fasta(args.prot):
            sid = header.split()[0]
            if sid not in raw_cds:
                print(f'ERROR: sequence "{sid}" not found in CDS input', file=sys.stderr)
                sys.exit(1)

            cds = raw_cds[sid]
            codon_idx = 0
            out = []

            for aa in prot:
                if aa == '-':
                    out.append('---')
                else:
                    start = codon_idx * 3
                    codon = cds[start:start + 3]
                    if len(codon) != 3:
                        print(
                            f'ERROR: ran out of codons for "{sid}" at protein position {codon_idx} '
                            f'(raw CDS length {len(cds)}, non-gap protein length '
                            f'{sum(1 for c in prot if c != "-")})',
                            file=sys.stderr
                        )
                        sys.exit(1)
                    out.append(codon)
                    codon_idx += 1

            fh.write(f'>{header}\n{"".join(out)}\n')


if __name__ == '__main__':
    main()
