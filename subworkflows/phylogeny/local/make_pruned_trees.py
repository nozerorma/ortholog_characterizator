#!/usr/bin/env python3
"""
Produce two pruned species trees from the assembly-name reference tree.

1. speciesTree_assemblyname_pruned.nh
   Assembly-name tips, one per canonical group.
   The kept assembly is chosen by the same priority used in dedup
   (max n_genes → scaffold N50 → contig N50), falling back to any assembly
   in the tree if the best-quality one is absent.

2. speciesTree_speciesname_pruned.nh
   Same topology, tips relabeled to canonical_group name.

Groups with no assembly in the reference tree are excluded from both outputs
and listed in a report.

Usage:
    python scripts/local/make_pruned_trees.py \
        --tree          external_phylogenies/speciesTree_assemblyname.nh \
        --dedup_map     taxa_info/species_dedup_map.tsv \
        --metadata      taxa_info/assemblies_and_species_mammals.tsv \
        --gene_counts   <gene_counts.tsv> \
        --out_asm       external_phylogenies/speciesTree_assemblyname_pruned.nh \
        --out_sp        external_phylogenies/speciesTree_speciesname_pruned.nh \
        --report        external_phylogenies/pruned_tree_report.tsv
"""

import argparse
import csv
import sys
from collections import defaultdict
from Bio import Phylo
from copy import deepcopy


def load_gene_counts(path: str) -> dict[str, int]:
    counts = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            counts[row["assembly_name"].strip()] = int(row["n_genes"])
    return counts


def load_metadata(path: str) -> dict[str, dict]:
    meta = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            asm = row["Assembly name"].strip()
            if asm:
                meta[asm] = row
    return meta


def load_dedup_map(path: str) -> dict[str, list[str]]:
    """canonical_group -> [assembly_names]"""
    groups: dict[str, list[str]] = defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            groups[row["canonical_group"].strip()].append(row["assembly_name"].strip())
    return groups


def int_field(row: dict, field: str) -> int:
    try:
        return int(str(row[field]).replace(",", "").strip())
    except (ValueError, KeyError, TypeError):
        return 0


def select_assembly(candidates: list[str], tree_tips: set[str],
                    gene_counts: dict, meta: dict) -> str | None:
    """
    From candidates, pick the one in the tree with highest quality.
    Quality: n_genes desc → scaffold N50 desc → contig N50 desc.
    Returns None if no candidate is in the tree.
    """
    in_tree = [a for a in candidates if a in tree_tips]
    if not in_tree:
        return None
    return min(in_tree, key=lambda a: (
        -gene_counts.get(a, 0),
        -int_field(meta.get(a, {}), "scaffold N50 (bp)"),
        -int_field(meta.get(a, {}), "contig N50 (bp)"),
    ))


def prune_tree(tree, keep_tips: set[str]):
    """Remove all terminals not in keep_tips, then collapse internal singletons."""
    for tip in tree.get_terminals():
        if tip.name not in keep_tips:
            tree.prune(tip.name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree",        required=True)
    parser.add_argument("--dedup_map",   required=True)
    parser.add_argument("--metadata",    required=True)
    parser.add_argument("--gene_counts", required=True)
    parser.add_argument("--out_asm",     required=True,
                        help="Output: assembly-name pruned tree")
    parser.add_argument("--out_sp",      required=True,
                        help="Output: species-name pruned tree")
    parser.add_argument("--report",      default=None,
                        help="TSV report of per-group decisions")
    args = parser.parse_args()

    gene_counts = load_gene_counts(args.gene_counts)
    meta        = load_metadata(args.metadata)
    groups      = load_dedup_map(args.dedup_map)

    tree = Phylo.read(args.tree, "newick")
    tree_tips = {c.name for c in tree.get_terminals()}

    print(f"Tree tips      : {len(tree_tips)}", file=sys.stderr)
    print(f"Canonical groups: {len(groups)}", file=sys.stderr)

    # Select one assembly per group
    selection: dict[str, str] = {}   # canonical_group -> chosen assembly
    absent: list[dict] = []

    for group, candidates in sorted(groups.items()):
        chosen = select_assembly(candidates, tree_tips, gene_counts, meta)
        if chosen is None:
            absent.append({"canonical_group": group,
                           "assemblies": ",".join(sorted(candidates))})
        else:
            selection[group] = chosen

    print(f"Groups with tree coverage : {len(selection)}", file=sys.stderr)
    print(f"Groups absent from tree   : {len(absent)}", file=sys.stderr)

    keep_asms = set(selection.values())

    # ── Assembly-name pruned tree ─────────────────────────────────────────────
    asm_tree = deepcopy(tree)
    prune_tree(asm_tree, keep_asms)
    Phylo.write(asm_tree, args.out_asm, "newick")
    remaining = len(asm_tree.get_terminals())
    print(f"Assembly tree  : {remaining} tips → {args.out_asm}", file=sys.stderr)

    # ── Species-name pruned tree ──────────────────────────────────────────────
    asm_to_group = {asm: grp for grp, asm in selection.items()}
    sp_tree = deepcopy(asm_tree)
    for tip in sp_tree.get_terminals():
        tip.name = asm_to_group.get(tip.name, tip.name)
    Phylo.write(sp_tree, args.out_sp, "newick")
    print(f"Species tree   : {remaining} tips → {args.out_sp}", file=sys.stderr)

    # ── Report ────────────────────────────────────────────────────────────────
    if args.report:
        with open(args.report, "w", newline="") as fh:
            writer = csv.DictWriter(fh, delimiter="\t",
                                    fieldnames=["canonical_group", "chosen_assembly",
                                                "n_genes", "scaffold_N50", "status"])
            writer.writeheader()
            for grp, asm in sorted(selection.items()):
                m = meta.get(asm, {})
                writer.writerow({
                    "canonical_group": grp,
                    "chosen_assembly": asm,
                    "n_genes":         gene_counts.get(asm, 0),
                    "scaffold_N50":    int_field(m, "scaffold N50 (bp)"),
                    "status":          "in_tree",
                })
            for row in absent:
                writer.writerow({
                    "canonical_group": row["canonical_group"],
                    "chosen_assembly": row["assemblies"],
                    "n_genes":         "",
                    "scaffold_N50":    "",
                    "status":          "absent_from_tree",
                })
        print(f"Report         : {args.report}", file=sys.stderr)


if __name__ == "__main__":
    main()
