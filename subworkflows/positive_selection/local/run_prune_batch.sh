#!/usr/bin/env bash
# Batch species-tree pruning for positive selection.
# Prunes the species tree to taxa present in each gene's alignment.

set -uo pipefail

manifest=""
workers=4
species_tree=""
scripts_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)     manifest="$2";      shift 2 ;;
        --workers)      workers="$2";       shift 2 ;;
        --species-tree) species_tree="$2";  shift 2 ;;
        --scripts-dir)  scripts_dir="$2";   shift 2 ;;
        *) echo "[PRUNE] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$manifest"     ]] && { echo "[PRUNE] --manifest is required"      >&2; exit 1; }
[[ -z "$species_tree" ]] && { echo "[PRUNE] --species-tree is required"  >&2; exit 1; }
[[ -z "$scripts_dir"  ]] && { echo "[PRUNE] --scripts-dir is required"   >&2; exit 1; }

PYTHON="${PYTHON:-python3}"
PRUNE_PY="${scripts_dir}/prune_tree_for_gene.py"

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        wait -n 2>/dev/null || true
    done
}

while IFS=$'\t' read -r gene_id fasta_name; do
    [[ -z "$gene_id" || "$gene_id" == \#* ]] && continue
    [[ ! -f "$fasta_name" ]] && { echo "[PRUNE] WARN: $fasta_name not found, skipping $gene_id" >&2; continue; }

    wait_for_slot

    (
        out_tree="${gene_id}.pruned.nwk"
        # prune_tree_for_gene.py: --tree, --fasta, --out, [--min_taxa]
        if "$PYTHON" "$PRUNE_PY" \
                --fasta "$fasta_name" \
                --tree  "$species_tree" \
                --out   "$out_tree" 2>&1; then
            [[ ! -s "$out_tree" ]] && rm -f "$out_tree"
            echo "[PRUNE] OK: $gene_id"
        else
            echo "[PRUNE] WARN: prune failed for $gene_id" >&2
            rm -f "$out_tree"
        fi
    ) &

done < "$manifest"

while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
    wait -n 2>/dev/null || true
done

echo "[PRUNE] Batch complete"
