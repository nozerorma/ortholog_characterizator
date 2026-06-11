#!/usr/bin/env python3

#
# scripts/4b.prune_reference_tree.py
#
# General-purpose tree pruner and re-rooter using dendropy.
# Handles Newick (e.g. Kuderna primate tree with species-name tips) and
# Nexus (e.g. PHYLACINE with numeric TRANSLATE tips) input formats.
#
# Usage:
#   python scripts/4b.prune_reference_tree.py \
#       --tree <input_tree> \
#       [--schema newick|nexus] \
#       --keep <taxa_list.txt> \
#       [--name_map <map.tsv>] \
#       [--outgroup <tip1,tip2,...>] \
#       --out_dir <output_dir> \
#       --prefix <output_prefix>
#
# Arguments:
#   --tree      Input tree file (Newick or Nexus)
#   --schema    Tree format: newick or nexus (default: auto-detect from extension)
#   --keep      Text file with tip labels to keep (one per line; uses output labels,
#               i.e. after applying name_map)
#   --name_map  Optional TSV (col1=tip_label_in_tree, col2=desired_output_label)
#   --outgroup  Comma-separated outgroup tip labels for rooting (output labels)
#   --out_dir   Output directory
#   --prefix    Output file prefix
#
# Outputs:
#   <out_dir>/<prefix>.rooted.nwk   — rooted on outgroup (or midpoint if absent)
#   <out_dir>/<prefix>.unrooted.nwk — unrooted Newick
#

import argparse
import os
import subprocess
import sys
import tempfile
import textwrap

try:
    import dendropy
except ImportError:
    dendropy = None


def detect_schema(path):
    """Auto-detect tree format from file extension."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".nex", ".nexus"):
        return "nexus"
    elif ext in (".nwk", ".newick", ".tree", ".tre"):
        return "newick"
    else:
        # Try to peek at the file content
        with open(path) as fh:
            first_line = fh.readline().strip()
        if first_line.startswith("#NEXUS") or first_line.startswith("#nexus"):
            return "nexus"
        return "newick"


def load_name_map(path):
    """Load TSV mapping tree tip labels -> desired output labels."""
    name_map = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                name_map[parts[0].strip()] = parts[1].strip()
    return name_map


def load_keep_set(path):
    """Load set of tip labels to retain (output/final labels)."""
    keep = set()
    with open(path) as fh:
        for line in fh:
            label = line.strip()
            if label:
                keep.add(label)
    return keep


def prune_with_r(
    tree_path,
    schema,
    keep_path,
    name_map_path,
    outgroup,
    out_dir,
    prefix,
):
    """Fallback pruner using ape/phangorn when dendropy is unavailable."""
    r_code = textwrap.dedent(
        r"""
        args <- commandArgs(trailingOnly = TRUE)
        tree_path <- args[[1]]
        schema <- args[[2]]
        keep_path <- args[[3]]
        name_map_path <- args[[4]]
        outgroup_str <- args[[5]]
        out_dir <- args[[6]]
        prefix <- args[[7]]

        suppressPackageStartupMessages({
          library(ape)
          library(phangorn)
        })

        if (schema == "nexus") {
          tree <- read.nexus(tree_path)
        } else {
          tree <- read.tree(tree_path)
        }

        keep_set <- scan(keep_path, what = character(), quiet = TRUE)
        cat(sprintf("Keeping %d taxa from: %s\n", length(keep_set), keep_path))

        name_map <- c()
        if (nzchar(name_map_path)) {
          map_df <- read.delim(
            name_map_path,
            header = FALSE,
            sep = "\t",
            quote = "",
            comment.char = "#",
            stringsAsFactors = FALSE
          )
          if (nrow(map_df) > 0) {
            name_map <- stats::setNames(map_df[[2]], map_df[[1]])
          }
          cat(sprintf("Loaded %d entries from name_map: %s\n", length(name_map), name_map_path))
        }

        final_labels <- tree$tip.label
        if (length(name_map) > 0) {
          mapped_idx <- tree$tip.label %in% names(name_map)
          final_labels[mapped_idx] <- unname(name_map[tree$tip.label[mapped_idx]])
        }

        keep_idx <- final_labels %in% keep_set
        cat(sprintf("Pruning %d taxa not in keep set ...\n", sum(!keep_idx)))
        tree <- drop.tip(tree, tree$tip.label[!keep_idx])

        if (length(tree$tip.label) == 0) {
          stop("No taxa remaining after pruning. Check --keep labels.")
        }

        if (length(name_map) > 0) {
          mapped_idx <- tree$tip.label %in% names(name_map)
          tree$tip.label[mapped_idx] <- unname(name_map[tree$tip.label[mapped_idx]])
        }

        cat(sprintf("%d taxa remain after pruning.\n", length(tree$tip.label)))

        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        unrooted_path <- file.path(out_dir, paste0(prefix, ".unrooted.nwk"))
        rooted_path <- file.path(out_dir, paste0(prefix, ".rooted.nwk"))

        write.tree(tree, file = unrooted_path)
        cat(sprintf("Unrooted tree written -> %s\n", unrooted_path))

        if (nzchar(outgroup_str)) {
          outgroups <- trimws(strsplit(outgroup_str, ",", fixed = TRUE)[[1]])
          present_outgroup <- outgroups[outgroups %in% tree$tip.label]
          if (length(present_outgroup) > 0) {
            cat(sprintf("Rooting on outgroup: %s\n", paste(present_outgroup, collapse = ", ")))
            rooted <- root(tree, outgroup = present_outgroup, resolve.root = TRUE)
          } else {
            cat(sprintf(
              "WARNING: None of the outgroup taxa (%s) found in pruned tree. Falling back to midpoint rooting.\n",
              paste(outgroups, collapse = ", ")
            ))
            rooted <- midpoint(tree)
          }
        } else {
          cat("No outgroup specified. Using midpoint rooting.\n")
          rooted <- midpoint(tree)
        }

        write.tree(rooted, file = rooted_path)
        cat(sprintf("Rooted tree written   -> %s\n", rooted_path))
        """
    )

    with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False) as handle:
        handle.write(r_code)
        script_path = handle.name

    try:
        cmd = [
            "Rscript",
            script_path,
            tree_path,
            schema,
            keep_path,
            name_map_path or "",
            outgroup or "",
            out_dir,
            prefix,
        ]
        completed = subprocess.run(cmd, check=False)
        if completed.returncode != 0:
            sys.exit(completed.returncode)
    finally:
        try:
            os.unlink(script_path)
        except OSError:
            pass


def main():
    parser = argparse.ArgumentParser(
        description="Prune and re-root a reference species tree."
    )
    parser.add_argument("--tree", required=True, help="Input tree file")
    parser.add_argument(
        "--schema",
        choices=["newick", "nexus"],
        default=None,
        help="Tree format (default: auto-detect from extension)",
    )
    parser.add_argument(
        "--keep", required=True, help="File with tip labels to keep (one per line)"
    )
    parser.add_argument(
        "--name_map",
        default=None,
        help="Optional TSV: col1=tree_label, col2=output_label",
    )
    parser.add_argument(
        "--outgroup",
        default=None,
        help="Comma-separated outgroup tip labels (output labels) for rooting",
    )
    parser.add_argument("--out_dir", required=True, help="Output directory")
    parser.add_argument("--prefix", required=True, help="Output file prefix")
    args = parser.parse_args()

    # ── Determine schema ──────────────────────────────────────────────────────
    schema = args.schema or detect_schema(args.tree)
    print(f"Loading tree: {args.tree} (schema={schema})")

    if dendropy is None:
        print("dendropy not available; falling back to R ape/phangorn pruner.")
        prune_with_r(
            tree_path=args.tree,
            schema=schema,
            keep_path=args.keep,
            name_map_path=args.name_map,
            outgroup=args.outgroup,
            out_dir=args.out_dir,
            prefix=args.prefix,
        )
        return

    # ── Load tree ─────────────────────────────────────────────────────────────
    tree = dendropy.Tree.get(path=args.tree, schema=schema, preserve_underscores=True)

    # ── Load name_map if provided ─────────────────────────────────────────────
    name_map = {}
    if args.name_map:
        name_map = load_name_map(args.name_map)
        print(f"Loaded {len(name_map)} entries from name_map: {args.name_map}")

    # ── Build a mapping: current_taxon_label -> final_output_label ───────────
    # For taxa not in name_map, the output label = current label
    def get_output_label(tip_label):
        return name_map.get(tip_label, tip_label)

    # ── Load keep set (uses output/final labels) ──────────────────────────────
    keep_set = load_keep_set(args.keep)
    print(f"Keeping {len(keep_set)} taxa from: {args.keep}")

    # ── Prune: retain only taxa whose final label is in keep_set ─────────────
    if not keep_set:
        print("ERROR: keep set is empty — no taxa to retain. Check --keep labels and reconciliation output.", file=sys.stderr)
        sys.exit(1)

    taxa_to_prune = []
    for taxon in tree.taxon_namespace:
        final_label = get_output_label(taxon.label)
        if final_label not in keep_set:
            taxa_to_prune.append(taxon)

    print(f"Pruning {len(taxa_to_prune)} taxa not in keep set ...")
    tree.prune_taxa(taxa_to_prune)
    tree.purge_taxon_namespace()

    # ── Rename tips using name_map ────────────────────────────────────────────
    if name_map:
        for taxon in tree.taxon_namespace:
            if taxon.label in name_map:
                taxon.label = name_map[taxon.label]

    # ── Validate: how many taxa remain ───────────────────────────────────────
    remaining = [t.label for t in tree.taxon_namespace]
    print(f"{len(remaining)} taxa remain after pruning.")
    if len(remaining) == 0:
        print("ERROR: No taxa remaining after pruning. Check --keep labels.", file=sys.stderr)
        sys.exit(1)

    # ── Output paths ──────────────────────────────────────────────────────────
    os.makedirs(args.out_dir, exist_ok=True)
    unrooted_path = os.path.join(args.out_dir, f"{args.prefix}.unrooted.nwk")
    rooted_path = os.path.join(args.out_dir, f"{args.prefix}.rooted.nwk")

    # ── Write unrooted ────────────────────────────────────────────────────────
    tree.write(
        path=unrooted_path,
        schema="newick",
        suppress_rooting=True,
        suppress_internal_node_labels=True,
    )
    print(f"Unrooted tree written -> {unrooted_path}")

    # ── Root tree ─────────────────────────────────────────────────────────────
    if args.outgroup:
        outgroup_labels = [t.strip() for t in args.outgroup.split(",")]
        # Only use outgroup taxa that are actually present in the pruned tree
        present_labels = set(t.label for t in tree.taxon_namespace)
        present_outgroup = [t for t in outgroup_labels if t in present_labels]

        if present_outgroup:
            print(f"Rooting on outgroup: {present_outgroup}")
            if len(present_outgroup) == 1:
                mrca_node = tree.find_node_with_taxon_label(present_outgroup[0])
            else:
                mrca_node = tree.mrca(taxon_labels=present_outgroup)
            tree.reroot_at_edge(mrca_node.edge, update_bipartitions=False)
        else:
            print(
                f"WARNING: None of the outgroup taxa ({outgroup_labels}) found in pruned tree. "
                "Falling back to midpoint rooting."
            )
            tree.reroot_at_midpoint(update_bipartitions=False)
    else:
        print("No outgroup specified. Using midpoint rooting.")
        tree.reroot_at_midpoint(update_bipartitions=False)

    # ── Write rooted ──────────────────────────────────────────────────────────
    tree.write(
        path=rooted_path,
        schema="newick",
        suppress_rooting=False,
        suppress_internal_node_labels=True,
    )
    print(f"Rooted tree written   -> {rooted_path}")


if __name__ == "__main__":
    main()
