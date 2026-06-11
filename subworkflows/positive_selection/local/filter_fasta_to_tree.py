#!/usr/bin/env python3
"""
Filter a FASTA alignment to sequences present in a Newick tree, renaming
them to their canonical tree-leaf name when the FASTA uses a synonym.

Name resolution order:
  1. Direct string match (FASTA name == tree leaf name).
  2. Taxid-based fuzzy match: look up the FASTA name in the taxid TSV,
     find the taxid, then check whether any of that taxid's names appears
     as a tree leaf. If so, rename the sequence to that tree-leaf name.

This fixes HyPhy "number of tree tips != number of sequences" errors
caused by synonym/spelling differences between the CDS alignments and
the species tree (e.g. Galagoides_demidovii in FASTA vs
Galagoides_demidoff in tree).

Usage:
    filter_fasta_to_tree.py --fasta gene.fa --tree gene.pruned.nwk --out filtered.fa
                            [--taxid-tsv taxid_species_family.tsv]
"""

import argparse
import csv
import re
import sys


# ── Taxid lookup ──────────────────────────────────────────────────────────────

def build_taxid_lookup(taxid_tsv):
    """
    Returns (name_to_taxid, taxid_to_names) using ALL rows in the TSV
    regardless of name_class or rank.
    """
    name_to_taxid = {}
    taxid_to_names = {}
    with open(taxid_tsv) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            sp  = row['species'].strip()
            tid = row['tax_id'].strip()
            if sp and tid:
                name_to_taxid[sp] = tid
                taxid_to_names.setdefault(tid, set()).add(sp)
    return name_to_taxid, taxid_to_names


def resolve_to_tree_name(query, name_to_taxid, taxid_to_names, tree_leaves):
    """Return the tree-leaf name for query, or None if unresolvable."""
    if query in tree_leaves:
        return query
    tid = name_to_taxid.get(query)
    if tid:
        for candidate in taxid_to_names.get(tid, set()):
            if candidate in tree_leaves:
                return candidate
    return None


# ── Newick leaf extraction ────────────────────────────────────────────────────

def get_leaf_names(newick_str):
    s = re.sub(r':[0-9.eE+\-]+', '', newick_str.strip())
    s = re.sub(r'\)[^(),;\n]*', ')', s)
    return set(re.findall(r'[^\s(),;]+', s))


# ── FASTA reader ──────────────────────────────────────────────────────────────

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


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--fasta',     required=True)
    ap.add_argument('--tree',      required=True)
    ap.add_argument('--out',       required=True)
    ap.add_argument('--taxid-tsv', default=None,
                    help='taxid_species_family.tsv for synonym → tree-name resolution')
    args = ap.parse_args()

    with open(args.tree) as fh:
        leaves = get_leaf_names(fh.read())

    name_to_taxid = taxid_to_names = None
    if args.taxid_tsv:
        name_to_taxid, taxid_to_names = build_taxid_lookup(args.taxid_tsv)

    kept, dropped, renamed = 0, 0, 0
    with open(args.out, 'w') as fh:
        for header, seq in read_fasta(args.fasta):
            sid = header.split()[0]

            if name_to_taxid:
                tree_name = resolve_to_tree_name(sid, name_to_taxid, taxid_to_names, leaves)
            else:
                tree_name = sid if sid in leaves else None

            if tree_name is None:
                dropped += 1
                continue

            if tree_name != sid:
                # Rename: replace first token in header with canonical tree name
                header = tree_name + (header[len(sid):] if header.startswith(sid) else '')
                renamed += 1

            fh.write(f'>{header}\n{seq}\n')
            kept += 1

    if dropped:
        print(f'[FILTER] {args.fasta}: kept {kept}, dropped {dropped} (no tree match)',
              file=sys.stderr)
    if renamed:
        print(f'[FILTER] {args.fasta}: renamed {renamed} sequences to canonical tree names',
              file=sys.stderr)

    if kept < 4:
        print(f'ERROR: only {kept} sequences remain after filtering (min 4 required)',
              file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
