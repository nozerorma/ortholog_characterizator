#!/usr/bin/env python3

#
# scripts/local/run_codeml_per_gene.py
#
# Per-gene codeML wrapper: prune a master species tree to the taxa present
# in a single alignment, convert to PHYLIP, write a codeML control file,
# run codeML, parse the output tree, and write the gene tree to out_dir.
#
# Usage:
#   python scripts/local/run_codeml_per_gene.py \
#       <ali.fa> <species_tree.nwk> <out_dir> <seqtype> <wag_dat>
#
# Arguments:
#   ali.fa           FASTA alignment (protein or codon)
#   species_tree.nwk Rooted master Newick tree (tip labels = FASTA headers)
#   out_dir          Output directory for gene trees
#   seqtype          1=codon, 2=protein
#   wag_dat          Path to WAG rate matrix file (needed for seqtype=2)
#
# Output:
#   <out_dir>/<gene_name>.nwk   — branch-length tree from codeML
#

import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap

import dendropy


# ── Helpers ───────────────────────────────────────────────────────────────────

def parse_fasta_headers(fasta_path):
    """Return list of sequence names from a FASTA file."""
    names = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].strip())
    return names


def fasta_to_dict(fasta_path):
    """Return dict {name: sequence} from a FASTA file."""
    seqs = {}
    current = None
    buf = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(buf)
                current = line[1:].strip()
                buf = []
            else:
                buf.append(line)
    if current is not None:
        seqs[current] = "".join(buf)
    return seqs


def make_phylip_name(name, used, max_len=10):
    """
    Return a PHYLIP-safe name (exactly max_len chars, padded with spaces).
    If name fits, just right-pad it. If it's too long, truncate and
    disambiguate with a counter suffix to ensure uniqueness.
    """
    if len(name) <= max_len:
        return name.ljust(max_len)
    # Truncate: keep first (max_len-2) chars + 2-digit counter
    base = name[:max_len - 2]
    counter = 0
    candidate = base + str(counter).zfill(2)
    while candidate in used:
        counter += 1
        candidate = base + str(counter).zfill(2)
    return candidate


def write_phylip(seqs, out_path):
    """
    Write sequences in strict sequential PHYLIP format.
    codeML requires taxon name padded to exactly 10 chars then >=1 space.
    Returns name_map {original_name: phylip_name}.
    """
    used_names = set()
    name_map = {}
    records = list(seqs.items())

    n_taxa = len(records)
    if n_taxa == 0:
        raise ValueError("No sequences to write.")
    seq_len = len(records[0][1])

    lines = [f" {n_taxa} {seq_len}"]
    for orig_name, seq in records:
        phy_name = make_phylip_name(orig_name, used_names, max_len=10)
        used_names.add(phy_name.strip())
        name_map[orig_name] = phy_name.strip()
        lines.append(f"{phy_name} {seq}")

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    return name_map


def prune_tree(tree, keep_labels):
    """Prune dendropy Tree to only the taxa in keep_labels set."""
    taxa_to_prune = [t for t in tree.taxon_namespace if t.label not in keep_labels]
    tree.prune_taxa(taxa_to_prune)
    tree.purge_taxon_namespace()


def rename_tree_tips(tree, name_map):
    """Rename tree tips using name_map {old: new}."""
    for taxon in tree.taxon_namespace:
        if taxon.label in name_map:
            taxon.label = name_map[taxon.label]


def write_codeml_ctl(ctl_path, seq_file, tree_file, out_file,
                     seqtype, wag_dat, work_dir):
    """Write a codeML control file."""
    model_line = "       model = 2" if seqtype == 2 else "       model = 0"
    aarat_line = (
        f"  aaRatefile = {wag_dat}" if seqtype == 2 else "* aaRatefile = (not used for codons)"
    )

    ctl = textwrap.dedent(f"""\
      seqfile = {seq_file}
     treefile = {tree_file}
      outfile = {out_file}
        noisy = 0
      verbose = 0
      runmode = -2
      seqtype = {seqtype}
{model_line}
{aarat_line}
    fix_alpha = 0
        alpha = 0.5
        ncatG = 4
        clock = 0
  fix_blength = 0
       method = 1
""")

    with open(ctl_path, "w") as fh:
        fh.write(ctl)


def extract_codeml_tree(codeml_out_path, name_unmap):
    """
    Parse the final optimized tree string from codeML output.
    codeML writes the branch-length tree near the end of the output file.
    Returns a cleaned Newick string with original tip names.
    name_unmap: {phylip_name: original_name}
    """
    newick_pattern = re.compile(r"\(.*\);")
    last_newick = None

    with open(codeml_out_path) as fh:
        for line in fh:
            stripped = line.strip()
            if newick_pattern.search(stripped):
                last_newick = stripped

    if last_newick is None:
        raise RuntimeError(f"No Newick tree found in codeML output: {codeml_out_path}")

    # Restore original tip names in the Newick string
    # Replace phylip names (they appear before ':' or ')' or ',')
    for phy_name, orig_name in sorted(name_unmap.items(), key=lambda x: -len(x[0])):
        last_newick = re.sub(
            r'(?<![A-Za-z0-9_])' + re.escape(phy_name) + r'(?=[:),])',
            orig_name,
            last_newick,
        )

    # Strip any node labels (digits / bootstrap values on internal nodes)
    last_newick = re.sub(r'\)[\d.]+(?=[:),;])', ')', last_newick)

    return last_newick


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 6:
        print(
            "Usage: run_codeml_per_gene.py <ali.fa> <species_tree.nwk> "
            "<out_dir> <seqtype> <wag_dat>",
            file=sys.stderr,
        )
        sys.exit(1)

    ali_fa       = sys.argv[1]
    species_tree = sys.argv[2]
    out_dir      = sys.argv[3]
    seqtype      = int(sys.argv[4])
    wag_dat      = sys.argv[5]

    gene_name = os.path.splitext(os.path.basename(ali_fa))[0]
    out_nwk   = os.path.join(out_dir, f"{gene_name}.nwk")

    # Skip if already done
    if os.path.exists(out_nwk):
        print(f"[{gene_name}] Output exists. Skipping.")
        return

    os.makedirs(out_dir, exist_ok=True)

    # ── Parse FASTA ───────────────────────────────────────────────────────────
    seqs = fasta_to_dict(ali_fa)
    taxa_in_gene = set(seqs.keys())
    if len(taxa_in_gene) < 3:
        print(f"[{gene_name}] Fewer than 3 taxa ({len(taxa_in_gene)}). Skipping.")
        return

    # ── Prune species tree ────────────────────────────────────────────────────
    master_tree = dendropy.Tree.get(
        path=species_tree, schema="newick", preserve_underscores=True
    )
    prune_tree(master_tree, taxa_in_gene)

    remaining = set(t.label for t in master_tree.taxon_namespace)
    if len(remaining) < 3:
        print(
            f"[{gene_name}] Fewer than 3 taxa in tree after pruning "
            f"({len(remaining)}). Skipping."
        )
        return

    # Keep only seqs that survived pruning
    seqs = {k: v for k, v in seqs.items() if k in remaining}

    # ── Work in a temp directory ──────────────────────────────────────────────
    work_dir = tempfile.mkdtemp(prefix=f"codeml_{gene_name}_")
    try:
        phy_path  = os.path.join(work_dir, f"{gene_name}.phy")
        tree_path = os.path.join(work_dir, f"{gene_name}.nwk")
        ctl_path  = os.path.join(work_dir, f"{gene_name}.ctl")
        out_path  = os.path.join(work_dir, f"{gene_name}.codeml.out")

        # Write PHYLIP
        name_map = write_phylip(seqs, phy_path)
        name_unmap = {v: k for k, v in name_map.items()}

        # Rename tree tips to PHYLIP names and write
        rename_tree_tips(master_tree, name_map)
        master_tree.write(
            path=tree_path,
            schema="newick",
            suppress_rooting=False,
            suppress_internal_node_labels=True,
        )

        # Write codeML control file
        write_codeml_ctl(
            ctl_path=ctl_path,
            seq_file=phy_path,
            tree_file=tree_path,
            out_file=out_path,
            seqtype=seqtype,
            wag_dat=wag_dat,
            work_dir=work_dir,
        )

        # Run codeML
        result = subprocess.run(
            ["codeml", os.path.basename(ctl_path)],
            cwd=work_dir,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(
                f"[{gene_name}] codeML failed (exit {result.returncode}):\n"
                f"{result.stderr[:500]}",
                file=sys.stderr,
            )
            return

        if not os.path.exists(out_path):
            print(
                f"[{gene_name}] codeML did not produce output file: {out_path}",
                file=sys.stderr,
            )
            return

        # Parse and clean tree from codeML output
        cleaned_nwk = extract_codeml_tree(out_path, name_unmap)

        # Write final gene tree
        with open(out_nwk, "w") as fh:
            fh.write(cleaned_nwk + "\n")

        print(f"[{gene_name}] Gene tree written -> {out_nwk}")

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
