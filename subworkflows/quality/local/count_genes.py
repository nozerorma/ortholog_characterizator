#!/usr/bin/env python3
"""
Count how many gene alignment files each assembly appears in.

Scans a directory of FASTA alignment files and outputs a TSV with two columns:
  assembly_name   n_genes

Usage:
    python count_genes.py -d <cds_dir> -o <out.tsv>
"""

import argparse
import csv
import os
from collections import defaultdict


def count_genes(fasta_dir: str) -> dict[str, int]:
    counts: dict[str, int] = defaultdict(int)
    for fname in sorted(os.listdir(fasta_dir)):
        if not (fname.endswith(".fa") or fname.endswith(".fasta")):
            continue
        fpath = os.path.join(fasta_dir, fname)
        with open(fpath) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                asm = line[1:].strip().split()[0]
                if asm.startswith("vs_"):
                    asm = asm[3:]
                counts[asm] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(
        description="Count gene alignment files per assembly across a CDS directory."
    )
    parser.add_argument("-d", "--dir", required=True,
                        help="Directory of FASTA alignment files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output TSV (assembly_name, n_genes)")
    args = parser.parse_args()

    counts = count_genes(args.dir)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["assembly_name", "n_genes"])
        for asm, n in sorted(counts.items()):
            writer.writerow([asm, n])

    print(f"Wrote gene counts for {len(counts)} assemblies to {args.output}")


if __name__ == "__main__":
    main()
