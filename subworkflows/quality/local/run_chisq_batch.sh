#!/usr/bin/env bash
# Batch chi-square quality filter wrapper.
# Reads a manifest (gene_id <TAB> fasta_filename) and runs chisq-test.py
# for each gene as a background subshell.
#
# Script interface (chisq-test.py):
#   -i/--input                required
#   -o/--output               required
#   -ot/--output_test         required (chi2 stats TSV)
#   -or/--output_removed      required (removed sequences TSV)
#   --max_iter                default 3
#   --ambiguity_threshold     default 0.8
#   --min_seq_fraction        default 0.1

set -uo pipefail

manifest=""
workers=8
max_iter=3
ambiguity_threshold=0.8
min_seq_fraction=0.1
min_seqs=4
scripts_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)            manifest="$2";            shift 2 ;;
        --workers)             workers="$2";             shift 2 ;;
        --max-iter)            max_iter="$2";            shift 2 ;;
        --ambiguity-thresh)    ambiguity_threshold="$2"; shift 2 ;;
        --min-seq-frac)        min_seq_fraction="$2";    shift 2 ;;
        --min-seqs)            min_seqs="$2";            shift 2 ;;
        --scripts-dir)         scripts_dir="$2";         shift 2 ;;
        *) echo "[CHI2] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$manifest"    ]] && { echo "[CHI2] --manifest is required"    >&2; exit 1; }
[[ -z "$scripts_dir" ]] && { echo "[CHI2] --scripts-dir is required" >&2; exit 1; }

PYTHON="${PYTHON:-python3}"
CHI2_PY="${scripts_dir}/chisq-test.py"

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        wait -n 2>/dev/null || true
    done
}

while IFS=$'\t' read -r gene_id fasta_name; do
    [[ -z "$gene_id" || "$gene_id" == \#* ]] && continue
    [[ ! -f "$fasta_name" ]] && { echo "[CHI2] WARN: $fasta_name not found, skipping $gene_id" >&2; continue; }

    wait_for_slot

    (
        out_fa="${gene_id}.filtered.fa"
        stats_f="${gene_id}.chi2.txt"
        removed_f="${gene_id}.removed.txt"

        if "$PYTHON" "$CHI2_PY" \
                --input               "$fasta_name" \
                --output              "$out_fa" \
                --output_test         "$stats_f" \
                --output_removed      "$removed_f" \
                --max_iter            "$max_iter" \
                --ambiguity_threshold "$ambiguity_threshold" \
                --min_seq_fraction    "$min_seq_fraction" \
                --ambiguous_dir       "ambiguous" \
                --low_seq_dir         "low_seq" \
                2>&1; then

            if [[ -f "$out_fa" ]]; then
                n_seqs=$(grep -c '^>' "$out_fa" 2>/dev/null || echo 0)
                if [[ "$n_seqs" -lt "$min_seqs" ]]; then
                    echo "[CHI2] WARN: $gene_id only $n_seqs seqs (< $min_seqs), moving to low_seq/" >&2
                    mkdir -p low_seq
                    mv "$out_fa" "low_seq/${out_fa}"
                elif [[ ! -s "$out_fa" ]]; then
                    rm -f "$out_fa"
                else
                    echo "[CHI2] OK: $gene_id ($n_seqs seqs)"
                fi
            else
                echo "[CHI2] ROUTED: $gene_id (ambiguous or low_seq)"
            fi
        else
            echo "[CHI2] WARN: chi2 filter failed for $gene_id" >&2
            rm -f "$out_fa"
        fi
    ) &

done < "$manifest"

while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
    wait -n 2>/dev/null || true
done

echo "[CHI2] Batch complete"
