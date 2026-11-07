#!/usr/bin/env python3
# full_length_query.py — Full-length human protein of a codon alignment, as the hmmscan query.
# ortholog_characterizator | subworkflows/variability/local/

"""
FullLengthQuery: translates the human CDS of the codon alignment that enters BMGE, with its gaps removed.

Pfam is scanned on this sequence and not on the BMGE-trimmed one: trimming removes columns that are
interior to the protein, so the trimmed human sequence joins stretches that are not contiguous in the
protein and can split, shorten or lose domains. Each residue keeps the codon column (ori_codon_col of
the MAP table) it came from, which is what projects a hit onto the trimmed alignment.

Called by:  HMMSCAN_BATCH (CLI) and map_domain_variability.py (import full_query)
Inputs:     --manifest  TSV, one gene per line: query name, CDS alignment file name (files in --cds_dir)
            --cds_dir   directory holding those files (default: current directory)
Outputs:    --out_fa    one record per gene, named as the manifest query name; genes whose alignment has
                        no human record are skipped with a warning

Translation: standard code. A codon with an ambiguity character or a partial gap becomes X (X consumes
a query position, as before); an internal stop becomes X and a final stop is dropped.
"""

import argparse
import os
import re
import sys

BASES = 'TCAG'
AMINO = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON = {a + b + c: AMINO[16 * i + 4 * j + k]
         for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}

HUMAN = re.compile(r'hg38|Homo_sapiens', re.IGNORECASE)


def read_fasta(path):
    """Read a FASTA file into [(header, upper-cased sequence)]."""
    seqs, header, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    seqs.append((header, ''.join(buf)))
                header, buf = line[1:], []
            else:
                buf.append(line.upper())
    if header is not None:
        seqs.append((header, ''.join(buf)))
    return seqs


def human_row(path):
    """Aligned CDS of the human record (the first when several match), or None."""
    rows = [seq for header, seq in read_fasta(path) if HUMAN.search(header)]
    if len(rows) > 1:
        print(f'WARN {path}: {len(rows)} human records, the first is used', file=sys.stderr)
    return rows[0] if rows else None


def full_query(path):
    """(protein, columns, n_cols): the human protein, the 1-based codon column of each residue and the number of
    codon columns of the alignment; None without a human record."""
    row = human_row(path)
    if row is None:
        return None
    if len(row) % 3:
        raise ValueError(f'{path}: aligned CDS length is not a multiple of 3')
    protein, columns = [], []
    for c in range(len(row) // 3):
        codon = row[3 * c:3 * c + 3]
        if codon == '---':
            continue
        protein.append(CODON.get(codon, 'X'))
        columns.append(c + 1)
    if protein and protein[-1] == '*':
        protein.pop()
        columns.pop()
    return ''.join('X' if aa == '*' else aa for aa in protein), columns, len(row) // 3


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--manifest', required=True, help='TSV: query name, CDS alignment file name')
    ap.add_argument('--cds_dir', default='.')
    ap.add_argument('--out_fa', required=True)
    args = ap.parse_args()

    written = 0
    with open(args.manifest) as man, open(args.out_fa, 'w') as out:
        for line in man:
            if not line.strip():
                continue
            name, cds = line.rstrip('\n').split('\t')[:2]
            query = full_query(os.path.join(args.cds_dir, cds))
            if query is None:
                print(f'WARN {name}: no human record in {cds}, gene skipped', file=sys.stderr)
                continue
            out.write(f'>{name}\n{query[0]}\n')
            written += 1
    print(f'{written} queries written to {args.out_fa}', flush=True)


if __name__ == '__main__':
    main()
