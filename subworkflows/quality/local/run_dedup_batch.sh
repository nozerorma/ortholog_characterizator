#!/usr/bin/env bash
# Batch deduplication wrapper.
# Reads a manifest (gene_id <TAB> fasta_filename) and runs
# dedup_by_assembly_quality.py for each gene as a background subshell.
#
# Script interface (dedup_by_assembly_quality.py):
#   -i/--input        input codon alignment
#   -o/--output       output deduplicated alignment
#   -m/--metadata     assemblies_and_species_mammals.tsv  (required)
#   --gene_counts     TSV from count_genes.py    (required)
#   --dedup_map       species_dedup_map.tsv               (required)
#   --report          optional report TSV

set -euo pipefail

manifest=""
workers=8
gene_counts=""
metadata_tsv=""      # assemblies_and_species_mammals.tsv  →  -m/--metadata
dedup_map=""         # species_dedup_map.tsv               →  --dedup_map
scripts_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)     manifest="$2";     shift 2 ;;
        --workers)      workers="$2";      shift 2 ;;
        --gene-counts)  gene_counts="$2";  shift 2 ;;
        --metadata)     metadata_tsv="$2"; shift 2 ;;
        --dedup-map)    dedup_map="$2";    shift 2 ;;
        --scripts-dir)  scripts_dir="$2";  shift 2 ;;
        *) echo "[DEDUP] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$manifest"     ]] && { echo "[DEDUP] --manifest is required"     >&2; exit 1; }
[[ -z "$gene_counts"  ]] && { echo "[DEDUP] --gene-counts is required"  >&2; exit 1; }
[[ -z "$metadata_tsv" ]] && { echo "[DEDUP] --metadata is required"     >&2; exit 1; }
[[ -z "$dedup_map"    ]] && { echo "[DEDUP] --dedup-map is required"    >&2; exit 1; }
[[ -z "$scripts_dir"  ]] && { echo "[DEDUP] --scripts-dir is required"  >&2; exit 1; }

PYTHON="${PYTHON:-python3}"
DEDUP_PY="${scripts_dir}/dedup_by_assembly_quality.py"

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        wait -n 2>/dev/null || true
    done
}

while IFS=$'\t' read -r gene_id fasta_name; do
    [[ -z "$gene_id" || "$gene_id" == \#* ]] && continue
    [[ ! -f "$fasta_name" ]] && { echo "[DEDUP] WARN: $fasta_name not found, skipping $gene_id" >&2; continue; }

    wait_for_slot

    (
        out_fa="${gene_id}.dedup.fa"
        log_f="${gene_id}.dedup.log"
        report_f="${gene_id}.dedup_report.tsv"

        if "$PYTHON" "$DEDUP_PY" \
                --input       "$fasta_name" \
                --output      "$out_fa" \
                --metadata    "$metadata_tsv" \
                --gene_counts "$gene_counts" \
                --dedup_map   "$dedup_map" \
                --report      "$report_f" \
                > "$log_f" 2>&1; then
            [[ ! -s "$out_fa" ]] && rm -f "$out_fa"
            echo "[DEDUP] OK: $gene_id"
        else
            echo "[DEDUP] WARN: dedup failed for $gene_id" >&2
            rm -f "$out_fa"
        fi
    ) &

done < "$manifest"

while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
    wait -n 2>/dev/null || true
done

echo "[DEDUP] Batch complete"
