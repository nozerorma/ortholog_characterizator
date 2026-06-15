#!/usr/bin/env python3
"""
Deduplicate a codon MSA by selecting the best assembly per species/subspecies
group using a precomputed per-assembly gene-count table.

Header format expected: >vs_<Assembly name>   (e.g. >vs_HLcarPer3)

Grouping:
  Assemblies are grouped by 'canonical_group' from the dedup map (--dedup_map).
  This collapses both replicate assemblies and subspecies-level assemblies into
  a single representative per canonical species group.
  Assemblies absent from the dedup map are kept without deduplication.

Selection priority (higher is better):
  1. n_genes  — number of gene alignments the assembly appears in
  2. scaffold N50 (bp)
  3. contig N50 (bp)

Assemblies with n_genes == 0 are removed unconditionally before selection and
listed in a dedicated zero-annotation report.
"""

import argparse
import csv
import os
import sys
from pathlib import Path
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def load_gene_counts(tsv_path: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    with open(tsv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            counts[row["assembly_name"].strip()] = int(row["n_genes"])
    return counts


def load_dedup_map(tsv_path: str) -> dict[str, str]:
    """Return assembly_name -> canonical_group."""
    mapping: dict[str, str] = {}
    with open(tsv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mapping[row["assembly_name"].strip()] = row["canonical_group"].strip()
    return mapping


def load_metadata(tsv_path: str) -> dict[str, dict]:
    path = Path(tsv_path)
    if not path.is_file():
        if str(tsv_path).strip().lower() in {"", "true", "false"}:
            fallback = Path(__file__).resolve().parents[3] / "taxa_info" / "assemblies_and_species_mammals.tsv"
            if fallback.is_file():
                print(
                    f"WARNING: --metadata was {tsv_path!r}; using {fallback}",
                    file=sys.stderr,
                )
                path = fallback
        if not path.is_file():
            raise FileNotFoundError(f"Metadata TSV not found: {tsv_path}")

    meta: dict[str, dict] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            asm = row["Assembly name"].strip()
            if asm:
                meta[asm] = row
    return meta


def _int_field(meta_row: dict, field: str) -> int:
    try:
        return int(str(meta_row[field]).replace(",", "").strip())
    except (ValueError, KeyError, TypeError):
        return 0


def quality_key(asm: str, gene_counts: dict[str, int], meta: dict[str, dict]):
    """Lower tuple = better quality (we use min-sort convention)."""
    n_genes = gene_counts.get(asm, 0)
    row = meta.get(asm)
    scaffold_n50 = _int_field(row, "scaffold N50 (bp)") if row else 0
    contig_n50 = _int_field(row, "contig N50 (bp)") if row else 0
    return (-n_genes, -scaffold_n50, -contig_n50)


def _parse_asm(record_id: str) -> str:
    asm = record_id
    if asm.startswith("vs_"):
        asm = asm[3:]
    return asm.split()[0]


def dedup_alignment(alignment, gene_counts: dict, dedup_map: dict, meta: dict):
    """
    Returns:
        deduplicated MultipleSeqAlignment
        list of dedup report rows  (groups with >1 assembly resolved)
        list of zero-annotation report rows  (assemblies removed for n_genes==0)
        list of no-map report rows  (assemblies not in dedup_map, kept as-is)
    """
    zero_rows: list[dict] = []
    no_map_rows: list[dict] = []

    groups: dict[str, list[tuple]] = {}  # canonical_group -> [(record, asm)]
    no_map_records: list = []

    for record in alignment:
        asm = _parse_asm(record.id)
        n_genes = gene_counts.get(asm, 0)

        if n_genes == 0:
            zero_rows.append({"sequence_id": record.id, "assembly_name": asm, "n_genes": 0})
            continue

        group = dedup_map.get(asm)
        if group is None:
            no_map_records.append(record)
            no_map_rows.append({"sequence_id": record.id, "assembly_name": asm})
            continue

        groups.setdefault(group, []).append((record, asm))

    kept: list = []
    dedup_rows: list[dict] = []

    for group, candidates in groups.items():
        best_rec, best_asm = min(candidates, key=lambda x: quality_key(x[1], gene_counts, meta))
        kept.append(best_rec)
        if len(candidates) > 1:
            best_meta = meta.get(best_asm, {})
            dedup_rows.append({
                "canonical_group":    group,
                "n_assemblies":       len(candidates),
                "kept_assembly":      best_asm,
                "kept_n_genes":       gene_counts.get(best_asm, 0),
                "kept_scaffold_N50":  _int_field(best_meta, "scaffold N50 (bp)"),
                "kept_contig_N50":    _int_field(best_meta, "contig N50 (bp)"),
                "discarded":          ",".join(
                    asm for _, asm in candidates if asm != best_asm
                ),
            })

    kept.extend(no_map_records)
    return MultipleSeqAlignment(kept), dedup_rows, zero_rows, no_map_rows


def main():
    parser = argparse.ArgumentParser(
        description="Deduplicate a codon MSA by keeping the best assembly per species group."
    )
    parser.add_argument("-i", "--input",       required=True,
                        help="Input codon alignment (.fa)")
    parser.add_argument("-o", "--output",      required=True,
                        help="Output deduplicated alignment (.fa)")
    parser.add_argument("-m", "--metadata",    required=True,
                        help="assemblies_and_species_mammals.tsv")
    parser.add_argument("--gene_counts",       required=True,
                        help="TSV from count_genes.py (assembly_name, n_genes)")
    parser.add_argument("--dedup_map",         required=True,
                        help="TSV mapping assembly_name -> canonical_group (species_dedup_map.tsv)")
    parser.add_argument("--report",            default=None,
                        help="TSV report of resolved multi-assembly groups")
    parser.add_argument("--zero_annot_report", default=None,
                        help="TSV report of assemblies removed for zero annotations")
    parser.add_argument("--no_map_report",     default=None,
                        help="TSV report of assemblies absent from dedup map (kept as-is)")
    args = parser.parse_args()

    gene_counts = load_gene_counts(args.gene_counts)
    dedup_map = load_dedup_map(args.dedup_map)
    meta = load_metadata(args.metadata)
    alignment = AlignIO.read(args.input, "fasta")

    deduped, dedup_rows, zero_rows, no_map_rows = dedup_alignment(
        alignment, gene_counts, dedup_map, meta
    )

    AlignIO.write(deduped, args.output, "fasta")
    print(f"Wrote {len(deduped)} sequences to {args.output}")

    def _write_tsv(path, fields, rows):
        if not path:
            return
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    if dedup_rows:
        print(f"Resolved {len(dedup_rows)} groups with multiple assemblies")
        _write_tsv(args.report, [
            "canonical_group", "n_assemblies", "kept_assembly",
            "kept_n_genes", "kept_scaffold_N50", "kept_contig_N50", "discarded",
        ], dedup_rows)

    if zero_rows:
        print(f"Removed {len(zero_rows)} assemblies with zero gene annotations")
        _write_tsv(args.zero_annot_report,
                   ["sequence_id", "assembly_name", "n_genes"], zero_rows)

    if no_map_rows:
        print(f"WARNING: {len(no_map_rows)} assemblies not found in dedup map — kept as-is")
        _write_tsv(args.no_map_report,
                   ["sequence_id", "assembly_name"], no_map_rows)

    if not dedup_rows and not zero_rows:
        print("No duplicate or zero-annotation assemblies found.")


if __name__ == "__main__":
    main()
