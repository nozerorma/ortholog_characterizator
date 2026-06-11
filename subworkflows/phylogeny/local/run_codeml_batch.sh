#!/usr/bin/env bash
# Batch CodeML fixed-topology gene tree estimation.
# Calls run_codeml_per_gene.py for each gene in the manifest.
#
# run_codeml_per_gene.py usage (positional):
#   <ali.fa> <species_tree.nwk> <out_dir> <seqtype> <wag_dat>
#
# Required params:
#   --ref-tree   master species tree (required)
#   --seqtype    1=codon, 2=protein (default: 1)
#   --wag-dat    path to WAG rate matrix (required for seqtype=2, optional for 1)

set -uo pipefail

manifest=""
workers=4
ref_tree=""
seqtype=1
wag_dat=""
scripts_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)    manifest="$2";    shift 2 ;;
        --workers)     workers="$2";     shift 2 ;;
        --ref-tree)    ref_tree="$2";    shift 2 ;;
        --seqtype)     seqtype="$2";     shift 2 ;;
        --wag-dat)     wag_dat="$2";     shift 2 ;;
        --scripts-dir) scripts_dir="$2"; shift 2 ;;
        # cpu-each accepted but ignored (codeml is single-threaded)
        --cpu-each)    shift 2 ;;
        *) echo "[CODEML] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$manifest"  ]] && { echo "[CODEML] --manifest is required"  >&2; exit 1; }
[[ -z "$ref_tree"  ]] && { echo "[CODEML] --ref-tree is required"  >&2; exit 1; }
[[ -z "$scripts_dir"]] && { echo "[CODEML] --scripts-dir is required" >&2; exit 1; }

# wag_dat defaults to empty string (codeml wrapper handles missing wag for seqtype=1)
[[ -z "$wag_dat" ]] && wag_dat="/dev/null"

PYTHON="${PYTHON:-python3}"
CODEML_PY="${scripts_dir}/run_codeml_per_gene.py"

mkdir -p codeml_out

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        wait -n 2>/dev/null || true
    done
}

while IFS=$'\t' read -r gene_id fasta_name; do
    [[ -z "$gene_id" || "$gene_id" == \#* ]] && continue
    [[ ! -f "$fasta_name" ]] && { echo "[CODEML] WARN: $fasta_name not found, skipping $gene_id" >&2; continue; }

    wait_for_slot

    (
        # run_codeml_per_gene.py writes <out_dir>/<gene_name>.nwk
        if "$PYTHON" "$CODEML_PY" \
                "$fasta_name" \
                "$ref_tree" \
                "codeml_out" \
                "$seqtype" \
                "$wag_dat" \
                > "${gene_id}.codeml.log" 2>&1; then

            # Rename output to our naming convention: <gene>.codeml.nwk
            out_nwk="codeml_out/${gene_id}.nwk"
            [[ -f "$out_nwk" ]] && cp "$out_nwk" "${gene_id}.codeml.nwk"
            echo "[CODEML] OK: $gene_id"
        else
            echo "[CODEML] WARN: codeml failed for $gene_id" >&2
            rm -f "${gene_id}.codeml.nwk"
        fi
    ) &

done < "$manifest"

while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
    wait -n 2>/dev/null || true
done

echo "[CODEML] Batch complete"
