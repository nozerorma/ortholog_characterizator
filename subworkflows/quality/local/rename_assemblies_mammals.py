#!/usr/bin/env python3
"""
rename_assemblies.py — Rename alignment sequence headers from assembly IDs
to species names, and strip the '.codonified' suffix from filenames.

Runs after assembly-quality deduplication and before BMGE trimming.
At this point each gene alignment contains one representative assembly per
species, so a 1-to-1 assembly → species mapping is safe.

Species names are taken from the metadata TSV (the same file used by the
dedup step).  Spaces in species names are converted to underscores so the
names are valid Newick tokens and compatible with downstream tree tools.

Input filenames:  TRANSCRIPT#GENE.codonified.fa   (or any .fa)
Output filenames: TRANSCRIPT#GENE.fa              (.codonified stripped)

Usage
-----
python rename_assemblies.py \\
    --in_dir   <dedup_output_dir>  \\
    --out_dir  <renamed_output_dir> \\
    --metadata <assemblies_and_species.tsv>
    [--assembly_col  "Assembly name"]   # column header for short assembly ID
    [--species_col   "Species"]         # column header for species name
    [--dry_run]                         # print what would happen, write nothing
"""

import argparse
import sys
from pathlib import Path


def load_canonical_map(taxids_path: str) -> dict[str, str]:
    """
    Build species_name → canonical_name from the taxid TSV.
    Scientific name rows define the canonical: subspecies map to their parent
    binomial; species map to themselves.
    Synonym rows are resolved via their tax_id so assemblies whose metadata
    species name is a synonym (e.g. Microcebus_sp._3_GT-2019) still get
    remapped to the accepted canonical name (e.g. Microcebus_jonahi).
    """
    import csv
    canon: dict[str, str] = {}      # species_name → canonical
    tid_to_canon: dict[str, str] = {}  # tax_id → canonical (for synonym resolution)

    # Pass 1: scientific names define canonical groups
    with open(taxids_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sp         = row["species"].strip()
            tid        = row["tax_id"].strip()
            rank       = row.get("rank", "").strip().lower()
            name_class = row.get("name_class", "scientific name").strip().lower()
            if not sp or name_class != "scientific name":
                continue
            parts = sp.split("_")
            canonical = "_".join(parts[:2]) if (rank == "subspecies" or len(parts) >= 3) else sp
            canon[sp] = canonical
            tid_to_canon[tid] = canonical

    # Pass 2: synonyms resolve to the same canonical as their tax_id's scientific name
    with open(taxids_path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sp         = row["species"].strip()
            tid        = row["tax_id"].strip()
            name_class = row.get("name_class", "scientific name").strip().lower()
            if not sp or name_class == "scientific name":
                continue
            if tid in tid_to_canon:
                canon[sp] = tid_to_canon[tid]

    return canon


def load_assembly_map(tsv_path: str, assembly_col: str, species_col: str,
                      canonical_map: dict[str, str] | None = None) -> dict[str, str]:
    """
    Build assembly_id → species_name map from the metadata TSV.
    Spaces in species names are replaced with underscores.
    If canonical_map is provided, subspecies names are normalised to their
    parent binomial (e.g. Bison_bison_bison → Bison_bison).
    """
    mapping: dict[str, str] = {}
    with open(tsv_path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            asm_idx = header.index(assembly_col)
            sp_idx  = header.index(species_col)
        except ValueError as e:
            print(f"ERROR: column not found in {tsv_path}: {e}", file=sys.stderr)
            print(f"  Available columns: {header}", file=sys.stderr)
            sys.exit(1)

        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(asm_idx, sp_idx):
                continue
            asm_id  = parts[asm_idx].strip()
            species = parts[sp_idx].strip().replace(" ", "_")
            if asm_id:
                if canonical_map:
                    species = canonical_map.get(species, species)
                mapping[asm_id] = species

    return mapping


def rename_fasta(in_path: Path, out_path: Path,
                 asm_map: dict[str, str],
                 unknown_ids: set[str],
                 dry_run: bool) -> int:
    """
    Write a copy of in_path to out_path with FASTA headers renamed.
    Returns number of sequences renamed.
    """
    lines_out: list[str] = []
    n_renamed = 0
    n_unknown = 0

    with open(in_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                asm_id = line[1:].split()[0].strip()
                if asm_id in asm_map:
                    lines_out.append(f">{asm_map[asm_id]}\n")
                    n_renamed += 1
                else:
                    lines_out.append(line)
                    unknown_ids.add(asm_id)
                    n_unknown += 1
            else:
                lines_out.append(line)

    if not dry_run:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("".join(lines_out), encoding="utf-8")

    return n_renamed


def strip_codonified(stem: str) -> str:
    """Remove '.codonified' from filename stem if present."""
    return stem.replace(".codonified", "")


def main():
    parser = argparse.ArgumentParser(
        description="Rename alignment headers from assembly IDs to species names")
    parser.add_argument("--in_dir",       default=None,
                        help="Input directory of .fa files (batch mode)")
    parser.add_argument("--out_dir",      default=None,
                        help="Output directory (batch mode)")
    parser.add_argument("--in_file",      default=None,
                        help="Single input .fa file (single-file mode)")
    parser.add_argument("--out_file",     default=None,
                        help="Single output .fa file (single-file mode)")
    parser.add_argument("--metadata",     required=True,
                        help="Metadata TSV (same file used by dedup)")
    parser.add_argument("--taxids",       default=None,
                        help="taxid_species_family_mammals.tsv — when provided, subspecies "
                             "headers are normalised to their canonical binomial species name")
    parser.add_argument("--assembly_col", default="Assembly name",
                        help="Column containing short assembly ID (default: 'Assembly name')")
    parser.add_argument("--species_col",  default="Species",
                        help="Column containing species name (default: 'Species')")
    parser.add_argument("--dry_run",      action="store_true",
                        help="Print what would happen without writing any files")
    args = parser.parse_args()

    # Validate mode
    single_mode = args.in_file is not None
    batch_mode  = args.in_dir  is not None
    if single_mode == batch_mode:
        print("ERROR: specify either --in_file/--out_file OR --in_dir/--out_dir", file=sys.stderr)
        sys.exit(1)
    if single_mode and not args.out_file:
        print("ERROR: --in_file requires --out_file", file=sys.stderr)
        sys.exit(1)
    if batch_mode and not args.out_dir:
        print("ERROR: --in_dir requires --out_dir", file=sys.stderr)
        sys.exit(1)

    canonical_map = None
    if args.taxids:
        print(f"Loading canonical species map from {args.taxids} …", file=sys.stderr)
        canonical_map = load_canonical_map(args.taxids)
        print(f"  {len(canonical_map):,} species → canonical entries loaded.", file=sys.stderr)

    print(f"Loading assembly → species map from {args.metadata} …", file=sys.stderr)
    asm_map = load_assembly_map(args.metadata, args.assembly_col, args.species_col,
                                canonical_map=canonical_map)
    print(f"  {len(asm_map):,} assembly entries loaded.", file=sys.stderr)

    unknown_ids: set[str] = set()
    n_files = 0
    n_total_renamed = 0

    if single_mode:
        in_path  = Path(args.in_file)
        out_path = Path(args.out_file)
        # Strip .codonified from output filename if caller didn't already
        if ".codonified" in out_path.stem:
            out_path = out_path.parent / f"{strip_codonified(out_path.stem)}.fa"
        fa_files = [(in_path, out_path)]
    else:
        in_dir  = Path(args.in_dir)
        out_dir = Path(args.out_dir)
        raw = sorted(in_dir.glob("*.fa"))
        if not raw:
            print(f"ERROR: no .fa files found in {in_dir}", file=sys.stderr)
            sys.exit(1)
        fa_files = [(p, out_dir / f"{strip_codonified(p.stem)}.fa") for p in raw]

    for in_path, out_path in fa_files:
        n_renamed = rename_fasta(in_path, out_path, asm_map, unknown_ids, args.dry_run)
        n_files      += 1
        n_total_renamed += n_renamed

        if args.dry_run:
            print(f"  [dry] {in_path.name} → {out_path.name}  ({n_renamed} sequences renamed)")

    print(f"Processed {n_files:,} files — {n_total_renamed:,} sequence headers renamed.",
          file=sys.stderr)
    if unknown_ids:
        print(f"WARNING: {len(unknown_ids)} assembly IDs not found in metadata "
              f"(headers left unchanged):", file=sys.stderr)
        for uid in sorted(unknown_ids)[:20]:
            print(f"  {uid}", file=sys.stderr)
        if len(unknown_ids) > 20:
            print(f"  … and {len(unknown_ids)-20} more", file=sys.stderr)

    if not args.dry_run:
        dest = args.out_file if single_mode else str(out_dir)
        print(f"Output written to: {dest}", file=sys.stderr)


if __name__ == "__main__":
    main()
