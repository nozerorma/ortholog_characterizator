#!/usr/bin/env bash
# Batch variability computation wrapper.
# For each gene in the manifest:
#   compute_variability.py  → ./PROT_VAR_RAW/<gene>.entropy.tsv etc.
#
# Script interface:
#   compute_variability.py: --prot_ali <fa> --taxid_tsv <tsv> (required)
#     --family_order_tsv --out_dir <dir> --var_subdir PROT_VAR_RAW
#     --alpha --beta --gamma
#     → writes <out_dir>/PROT_VAR_RAW/<gene>.entropy.tsv
#                          <gene>.clade_entropy.tsv
#                          <gene>.var.fa (variable-only protein)

set -uo pipefail

manifest=""
workers=8
alpha=1.0
beta=1.0
gamma=1.0
taxid_tsv=""
family_order_tsv=""
scripts_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)         manifest="$2";         shift 2 ;;
        --workers)          workers="$2";           shift 2 ;;
        --alpha)            alpha="$2";             shift 2 ;;
        --beta)             beta="$2";              shift 2 ;;
        --gamma)            gamma="$2";             shift 2 ;;
        --taxid-tsv)        taxid_tsv="$2";         shift 2 ;;
        --family-order-tsv) family_order_tsv="$2";  shift 2 ;;
        --scripts-dir)      scripts_dir="$2";       shift 2 ;;
        *) echo "[VAR] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$manifest"    ]] && { echo "[VAR] --manifest is required"    >&2; exit 1; }
[[ -z "$taxid_tsv"   ]] && { echo "[VAR] --taxid-tsv is required (compute_variability.py requires --taxid_tsv)" >&2; exit 1; }
[[ -z "$scripts_dir" ]] && { echo "[VAR] --scripts-dir is required" >&2; exit 1; }

PYTHON="${PYTHON:-python3}"
VAR_PY="${scripts_dir}/compute_variability.py"

mkdir -p PROT_VAR_RAW

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        wait -n 2>/dev/null || true
    done
}

while IFS=$'\t' read -r gene_id fasta_name; do
    [[ -z "$gene_id" || "$gene_id" == \#* ]] && continue
    [[ ! -f "$fasta_name" ]] && { echo "[VAR] WARN: $fasta_name not found, skipping $gene_id" >&2; continue; }

    wait_for_slot

    (
        # compute variability (script writes to ./PROT_VAR_RAW/<gene>.*)
        var_args=(
            --prot_ali    "$fasta_name"
            --taxid_tsv   "$taxid_tsv"
            --out_dir     "."
            --var_subdir  "PROT_VAR_RAW"
            --alpha       "$alpha"
            --beta        "$beta"
            --gamma       "$gamma"
        )
        [[ -n "$family_order_tsv" ]] && var_args+=(--family_order_tsv "$family_order_tsv")

        if "$PYTHON" "$VAR_PY" "${var_args[@]}" 2>&1; then
            echo "[VAR] OK: $gene_id"
        else
            echo "[VAR] WARN: variability failed for $gene_id" >&2
        fi
    ) &

done < "$manifest"

while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
    wait -n 2>/dev/null || true
done

echo "[VAR] Batch complete"
