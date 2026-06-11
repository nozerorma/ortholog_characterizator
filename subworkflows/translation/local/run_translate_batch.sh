#!/usr/bin/env bash
# Batch manual translation wrapper.
# Runs translate_codon_msa.py (gap-column removal + codon→AA) per gene.
#
# Script interface (translate_codon_msa.py):
#   --input        required: codon MSA FASTA
#   --out_dir      required: base output dir; protein written to <out_dir>/PROT_RAW/<gene>.fa
#   --gap_threshold, --seq_ambiguity_threshold, --min_seqs (optional)
#
# This wrapper uses the work dir ('.') as --out_dir, so each gene writes to
# ./PROT_RAW/<gene>.fa.  The Nextflow process collects 'PROT_RAW/*.fa'.

set -uo pipefail

manifest=""
workers=6
scripts_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)    manifest="$2";    shift 2 ;;
        --workers)     workers="$2";     shift 2 ;;
        --scripts-dir) scripts_dir="$2"; shift 2 ;;
        *) echo "[TRANSLATE] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$manifest"    ]] && { echo "[TRANSLATE] --manifest is required"    >&2; exit 1; }
[[ -z "$scripts_dir" ]] && { echo "[TRANSLATE] --scripts-dir is required" >&2; exit 1; }

PYTHON="${PYTHON:-python3}"
TRANS_PY="${scripts_dir}/translate_codon_msa.py"

mkdir -p PROT_RAW

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        wait -n 2>/dev/null || true
    done
}

while IFS=$'\t' read -r gene_id fasta_name; do
    [[ -z "$gene_id" || "$gene_id" == \#* ]] && continue
    [[ ! -f "$fasta_name" ]] && { echo "[TRANSLATE] WARN: $fasta_name not found, skipping $gene_id" >&2; continue; }

    wait_for_slot

    (
        # Script writes to <out_dir>/PROT_RAW/<gene>.fa
        if "$PYTHON" "$TRANS_PY" \
                --input   "$fasta_name" \
                --out_dir "." \
                2>&1; then
            out_fa="PROT_RAW/${gene_id}.fa"
            [[ ! -s "$out_fa" ]] && rm -f "$out_fa"
            echo "[TRANSLATE] OK: $gene_id"
        else
            echo "[TRANSLATE] WARN: translate failed for $gene_id" >&2
        fi
    ) &

done < "$manifest"

while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
    wait -n 2>/dev/null || true
done

echo "[TRANSLATE] Batch complete"
