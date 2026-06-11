#!/usr/bin/env python3
"""
Filter a FASTA alignment to only sequences whose IDs appear as leaves
in a Newick tree.

Fixes HyPhy "number of tree tips != number of sequences" errors that
occur when the FASTA contains species absent from the species tree.
Those species pass through BMGE and dedup untouched because
prune_tree_for_gene.py only prunes the tree to match the FASTA — it
does not prune the FASTA to match the tree.

Usage:
    filter_fasta_to_tree.py --fasta gene.fa --tree gene.pruned.nwk --out filtered.fa
"""

import argparse
import re
import sys


def get_leaf_names(newick_str):
    s = newick_str.strip()
    # Remove branch lengths
    s = re.sub(r':[0-9.eE+\-]+', '', s)
    # Remove internal node labels (tokens immediately after a closing paren)
    s = re.sub(r'\)[^(),;\n]*', ')', s)
    # Everything left that is not punctuation is a leaf name
    return set(re.findall(r'[^\s(),;]+', s))


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
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--tree',  required=True)
    ap.add_argument('--out',   required=True)
    args = ap.parse_args()

    with open(args.tree) as fh:
        leaves = get_leaf_names(fh.read())

    kept, dropped = 0, 0
    with open(args.out, 'w') as fh:
        for header, seq in read_fasta(args.fasta):
            sid = header.split()[0]
            if sid in leaves:
                fh.write(f'>{header}\n{seq}\n')
                kept += 1
            else:
                dropped += 1

    if dropped:
        print(f'[FILTER] {args.fasta}: kept {kept}, dropped {dropped} '
              f'(not in tree)', file=sys.stderr)

    if kept < 4:
        print(f'ERROR: only {kept} sequences remain after filtering '
              f'(min 4 required)', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
