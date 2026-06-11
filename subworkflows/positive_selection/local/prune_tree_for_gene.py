#!/usr/bin/env python3
"""
Prune a Newick species tree to only the taxa present in a FASTA alignment.

Reads leaf labels from the tree and sequence headers from the FASTA; writes
a pruned Newick retaining only leaves that appear in both. Uses a minimal
recursive Newick parser — no external libraries required.

Usage:
  python prune_tree_for_gene.py --tree species.nwk --fasta gene.fa --out pruned.nwk
"""

import argparse
import csv
import os
import re
import sys


# ---------------------------------------------------------------------------
# Minimal Newick parser / pruner
# ---------------------------------------------------------------------------

class Node:
    __slots__ = ('name', 'length', 'children')

    def __init__(self, name='', length=None, children=None):
        self.name = name
        self.length = length
        self.children = children or []

    def is_leaf(self):
        return not self.children


def parse_newick(s):
    s = s.strip().rstrip(';')
    pos = [0]

    def parse_node():
        node = Node()
        if pos[0] < len(s) and s[pos[0]] == '(':
            pos[0] += 1  # consume '('
            node.children.append(parse_node())
            while pos[0] < len(s) and s[pos[0]] == ',':
                pos[0] += 1
                node.children.append(parse_node())
            if pos[0] < len(s) and s[pos[0]] == ')':
                pos[0] += 1
        # label
        label_start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in (',', ')', ':', ';'):
            pos[0] += 1
        node.name = s[label_start:pos[0]]
        # branch length
        if pos[0] < len(s) and s[pos[0]] == ':':
            pos[0] += 1
            bl_start = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in (',', ')', ';'):
                pos[0] += 1
            try:
                node.length = float(s[bl_start:pos[0]])
            except ValueError:
                node.length = None
        return node

    return parse_node()


def to_newick(node):
    if node.is_leaf():
        s = node.name
    else:
        s = '(' + ','.join(to_newick(c) for c in node.children) + ')' + node.name
    if node.length is not None:
        s += f':{node.length}'
    return s


def prune(node, keep):
    """
    Recursively prune the tree, retaining only leaves in `keep`.
    Returns None if this subtree has no retained leaves.
    Collapses single-child internal nodes.
    """
    if node.is_leaf():
        return node if node.name in keep else None

    surviving = []
    for child in node.children:
        pruned = prune(child, keep)
        if pruned is not None:
            surviving.append(pruned)

    if not surviving:
        return None
    if len(surviving) == 1:
        # collapse: pass the child up but keep its own branch length if internal
        surviving[0].length = node.length
        return surviving[0]
    node.children = surviving
    return node


# ---------------------------------------------------------------------------
# Taxid-based name resolution (synonym → tree-leaf canonical name)
# ---------------------------------------------------------------------------

def build_taxid_lookup(taxid_tsv):
    """
    Returns (name_to_taxid, taxid_to_names) using ALL rows in the TSV
    regardless of name_class or rank — no assumptions about which name is
    'scientific' or 'synonym'.
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


def get_all_leaf_names(node):
    if node.is_leaf():
        return {node.name}
    result = set()
    for child in node.children:
        result |= get_all_leaf_names(child)
    return result


def resolve_fasta_to_tree_names(fasta_headers, tree, name_to_taxid, taxid_to_names):
    """
    Return the set of tree-leaf names that correspond to the given FASTA
    headers, via direct match or taxid-based fuzzy match.
    """
    tree_leaves = get_all_leaf_names(tree)
    keep = set()
    for fasta_name in fasta_headers:
        if fasta_name in tree_leaves:
            keep.add(fasta_name)
            continue
        tid = name_to_taxid.get(fasta_name)
        if tid:
            for candidate in taxid_to_names.get(tid, set()):
                if candidate in tree_leaves:
                    keep.add(candidate)
                    break
    return keep


# ---------------------------------------------------------------------------
# FASTA header reader
# ---------------------------------------------------------------------------

def read_fasta_headers(path):
    headers = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                headers.add(line[1:].rstrip().split()[0])
    return headers


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Prune a species tree to taxa present in a gene alignment.')
    parser.add_argument('--tree', required=True, help='Input Newick tree')
    parser.add_argument('--fasta', required=True, help='Gene alignment FASTA')
    parser.add_argument('--out', required=True, help='Output pruned Newick')
    parser.add_argument('--min_taxa', type=int, default=4,
                        help='Minimum taxa after pruning; exits with code 1 if below (default: 4)')
    parser.add_argument('--taxid-tsv', default=None,
                        help='taxid_species_family.tsv for synonym → tree-name resolution')
    args = parser.parse_args()

    with open(args.tree) as fh:
        nwk = fh.read().strip()

    tree = parse_newick(nwk)
    fasta_headers = read_fasta_headers(args.fasta)

    if args.taxid_tsv:
        name_to_taxid, taxid_to_names = build_taxid_lookup(args.taxid_tsv)
        keep = resolve_fasta_to_tree_names(fasta_headers, tree, name_to_taxid, taxid_to_names)
    else:
        keep = fasta_headers

    pruned = prune(tree, keep)

    if pruned is None:
        print(f'ERROR: no taxa from {args.fasta} found in tree', file=sys.stderr)
        sys.exit(1)

    # Count leaves
    def count_leaves(n):
        return 1 if n.is_leaf() else sum(count_leaves(c) for c in n.children)

    n_leaves = count_leaves(pruned)
    if n_leaves < args.min_taxa:
        print(f'ERROR: only {n_leaves} taxa after pruning (min {args.min_taxa})', file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    with open(args.out, 'w') as fh:
        fh.write(to_newick(pruned) + ';\n')

    print(f'{os.path.basename(args.fasta).replace(".fa","")}: {n_leaves} taxa → {args.out}',
          flush=True)


if __name__ == '__main__':
    main()
