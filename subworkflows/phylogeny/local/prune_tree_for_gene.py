#!/usr/bin/env python3
"""
Prune a Newick species tree to only the taxa present in a FASTA alignment, or
(--root-only) reroot a tree on an outgroup clade. Uses a minimal recursive
Newick parser — no external libraries required.

Usage:
  python prune_tree_for_gene.py --tree species.nwk --fasta gene.fa --out pruned.nwk
  python prune_tree_for_gene.py --root-only --tree species.nwk --output rooted.nwk \\
      --outgroup taxonA,taxonB
"""

import argparse
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


def collapse_singletons(node):
    """
    Remove redundant single-child (degree-2) internal nodes, summing the two
    edge lengths into the surviving edge. These show up after rerooting
    because the original Newick's own arbitrary root is usually itself a
    degree-2 artifact (its two "children" are really one continuous edge,
    just split in two by where the string happened to start).
    """
    node.children = [collapse_singletons(c) for c in node.children]
    if len(node.children) == 1:
        child = node.children[0]
        combined = (node.length or 0.0) + (child.length or 0.0) if (
            node.length is not None or child.length is not None) else None
        child.length = combined
        return child
    return node


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
# Outgroup rerooting
# ---------------------------------------------------------------------------

def build_adjacency(root):
    """Undirected adjacency: adj[id(node)] = [(neighbor_node, edge_length), ...].

    The parsed tree is treated purely as a source of nodes/edges here -- which
    of two adjacent nodes the parser happened to call "parent" is an artifact
    of the Newick string's own arbitrary root and is not meaningful for an
    unrooted tree (ASTRAL/IQ-TREE output). Outgroup rerooting below works
    entirely off this undirected view, not off parse-tree ancestry.
    """
    adj = {}

    def rec(node):
        adj.setdefault(id(node), [])
        for c in node.children:
            length = c.length if c.length is not None else 0.0
            adj[id(node)].append((c, length))
            adj.setdefault(id(c), [])
            adj[id(c)].append((node, length))
            rec(c)

    rec(root)
    return adj


def leafset_away(node, came_from, adj, memo):
    """Leaf names reachable from `node` without crossing back to `came_from`."""
    key = (id(node), id(came_from) if came_from is not None else None)
    if key in memo:
        return memo[key]
    neighbors = [n for n, _ in adj[id(node)] if n is not came_from]
    if not neighbors:
        s = frozenset([node.name])
    else:
        s = frozenset().union(*(leafset_away(n, node, adj, memo) for n in neighbors))
    memo[key] = s
    return s


def find_best_split(root, targets):
    """
    Over every directed edge (u -> v) of the undirected tree, find the one
    whose "away from u" side is the smallest superset of `targets`. Returns
    (u, v, length) for that edge, or None if no such side exists (shouldn't
    happen for a non-empty `targets` that's a subset of the tree's leaves).
    """
    adj = build_adjacency(root)
    memo = {}
    best = None  # (size, u, v, length)

    id_to_node = {}

    def collect(node):
        id_to_node[id(node)] = node
        for c in node.children:
            collect(c)

    collect(root)

    for u_id, u in id_to_node.items():
        for v, length in adj[u_id]:
            side = leafset_away(v, u, adj, memo)
            if targets <= side:
                if best is None or len(side) < best[0]:
                    best = (len(side), u, v, length)
    if best is None:
        return None
    _, u, v, length = best
    return u, v, length


def reroot_between(root, u, v, length):
    """
    Reroot the undirected tree so a new root sits on the edge (u, v), with
    v's "away from u" subtree on one side and u's "away from v" subtree
    (i.e. the rest of the tree) on the other. `length` (the original edge's)
    is split evenly between the two new pendant edges -- the common
    convention (e.g. ete3's set_outgroup default).
    """
    adj = build_adjacency(root)

    def build_subtree(current, came_from):
        children = []
        for nbr, nbr_length in adj[id(current)]:
            if nbr is came_from:
                continue
            child = build_subtree(nbr, current)
            child.length = nbr_length
            children.append(child)
        return Node(name=current.name, length=None, children=children)

    v_side = build_subtree(v, u)
    u_side = build_subtree(u, v)

    half = length / 2.0
    v_side.length = half
    u_side.length = half

    return Node(name='', length=None, children=[v_side, u_side])


def all_leaf_names(node):
    if node.is_leaf():
        return {node.name}
    out = set()
    for c in node.children:
        out |= all_leaf_names(c)
    return out


def root_on_outgroup(root, outgroup_names):
    """
    Reroot on the smallest side of the best-matching edge for `outgroup_names`
    (see find_best_split) -- i.e. the outgroup's own bipartition, found on the
    tree's undirected topology rather than the Newick string's arbitrary root.
    Raises ValueError if none of the outgroup taxa are leaves of the tree.
    """
    leaves = all_leaf_names(root)
    present = outgroup_names & leaves
    if not present:
        raise ValueError(
            f'none of the outgroup taxa {sorted(outgroup_names)} are leaves of this tree'
        )
    missing = outgroup_names - leaves
    if missing:
        print(f'WARNING: outgroup taxa not in tree, ignored: {sorted(missing)}',
              file=sys.stderr)
    split = find_best_split(root, present)
    if split is None:
        raise ValueError(f'could not find a tree edge separating {sorted(present)}')
    u, v, length = split
    return reroot_between(root, u, v, length)


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
        description='Prune a species tree to taxa present in a gene alignment, '
                    'or (--root-only) reroot a tree on an outgroup clade.')
    parser.add_argument('--root-only', action='store_true',
                        help='Reroot --tree on --outgroup (or leave unrooted if no '
                             '--outgroup given) and write it to --output. Skips the '
                             'FASTA-based pruning entirely; --fasta/--out/--min_taxa '
                             'are not used in this mode.')
    parser.add_argument('--tree', required=True, help='Input Newick tree')
    parser.add_argument('--fasta', help='Gene alignment FASTA (prune mode)')
    parser.add_argument('--out', help='Output pruned Newick (prune mode)')
    parser.add_argument('--output', help='Output rerooted Newick (--root-only mode)')
    parser.add_argument('--outgroup', help='Comma-separated outgroup taxon names '
                        '(--root-only mode; rerooted on their MRCA). If omitted, '
                        'the tree is written out unchanged (no midpoint rooting).')
    parser.add_argument('--min_taxa', type=int, default=4,
                        help='Minimum taxa after pruning; exits with code 1 if below (default: 4)')
    args = parser.parse_args()

    with open(args.tree) as fh:
        nwk = fh.read().strip()

    tree = parse_newick(nwk)

    if args.root_only:
        if not args.output:
            parser.error('--root-only requires --output')
        if args.outgroup:
            names = {t.strip() for t in args.outgroup.split(',') if t.strip()}
            tree = root_on_outgroup(tree, names)
            tree = collapse_singletons(tree)
        else:
            print('WARNING: --root-only with no --outgroup — writing the tree '
                  'unrooted (midpoint rooting is not implemented)', file=sys.stderr)
        os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
        with open(args.output, 'w') as fh:
            fh.write(to_newick(tree) + ';\n')
        return

    if not args.fasta or not args.out:
        parser.error('--fasta and --out are required unless --root-only is set')

    keep = read_fasta_headers(args.fasta)

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
