#!/usr/bin/env bash
# Batch HyPhy analyses wrapper.
# Reads a manifest (gene_id <TAB> fasta_filename <TAB> tree_filename) and
# runs each requested HyPhy method per gene as a background subshell.
#
# Methods: comma-separated, e.g. "FUBAR,FEL" or "MEME"

set -uo pipefail

manifest=""
workers=4
cpu_each=4
methods="FUBAR,FEL"
scripts_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)    manifest="$2";    shift 2 ;;
        --workers)     workers="$2";     shift 2 ;;
        --cpu-each)    cpu_each="$2";    shift 2 ;;
        --methods)     methods="$2";     shift 2 ;;
        --scripts-dir) scripts_dir="$2"; shift 2 ;;
        *) echo "[HYPHY] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$manifest"    ]] && { echo "[HYPHY] --manifest is required"    >&2; exit 1; }

# Parse methods into array
IFS=',' read -ra METHOD_LIST <<< "$methods"

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        wait -n 2>/dev/null || true
    done
}

run_hyphy_gene() {
    local gene_id="$1"
    local fasta="$2"
    local tree="$3"

    # Cap BLAS/OpenMP threads to the per-worker CPU allocation.
    # HyPhy's CPU argument only controls its own scheduler; library-level
    # threads (OpenBLAS, OpenMP, MKL) detect hardware CPU count via sysconf()
    # and ignore CPU entirely unless these vars are set.
    export OMP_NUM_THREADS="${cpu_each}"
    export MKL_NUM_THREADS="${cpu_each}"
    export OPENBLAS_NUM_THREADS="${cpu_each}"
    export BLAS_NUM_THREADS="${cpu_each}"

    for method in "${METHOD_LIST[@]}"; do
        method_lower="${method,,}"
        out_json="${gene_id}.${method}.json"

        # Resume: skip if JSON already exists and is non-empty
        [[ -s "$out_json" ]] && { echo "[HYPHY] SKIP (exists): ${gene_id}.${method}"; continue; }

        if hyphy "${method_lower}" \
                --alignment "$fasta" \
                --tree      "$tree" \
                --output    "$out_json" \
                CPU="${cpu_each}" \
                > "${gene_id}.${method}.log" 2>&1; then
            [[ ! -s "$out_json" ]] && rm -f "$out_json"
            echo "[HYPHY] OK: ${gene_id} / ${method}"
        else
            echo "[HYPHY] WARN: ${method} failed for ${gene_id}" >&2
            rm -f "$out_json"
        fi
    done
}

while IFS=$'\t' read -r gene_id fasta_name tree_name; do
    [[ -z "$gene_id" || "$gene_id" == \#* ]] && continue
    [[ ! -f "$fasta_name" ]] && { echo "[HYPHY] WARN: $fasta_name not found, skipping $gene_id" >&2; continue; }
    [[ ! -f "$tree_name"  ]] && { echo "[HYPHY] WARN: $tree_name not found, skipping $gene_id"  >&2; continue; }

    wait_for_slot
    ( run_hyphy_gene "$gene_id" "$fasta_name" "$tree_name" ) &

done < "$manifest"

while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
    wait -n 2>/dev/null || true
done

echo "[HYPHY] Batch complete"
