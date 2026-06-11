#!/usr/bin/env bash
# Batch BMGE position-map generation.
# Reads a manifest (gene_id <TAB> track_html <TAB> cds_fasta) and runs
# parse_bmge_track.py for each gene.
#
# Script interface (parse_bmge_track.py):
#   --track_html      required: BMGE TRACK HTML
#   --cds_fasta       required: original codon FASTA used as BMGE input
#   --out_tsv         required: output TSV path
#   --gff             optional: hg38 GFF3 for genomic coordinate mapping
#   --gene_equiv_tsv  optional: gene→protein_id TSV for GFF lookup
#   --ref_species     optional: default 'Homo_sapiens'

set -uo pipefail

manifest=""
workers=6
gff=""
gene_equiv=""
scripts_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)    manifest="$2";    shift 2 ;;
        --workers)     workers="$2";     shift 2 ;;
        --gff)         gff="$2";         shift 2 ;;
        --gene-equiv)  gene_equiv="$2";  shift 2 ;;
        --scripts-dir) scripts_dir="$2"; shift 2 ;;
        *) echo "[MAP] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$manifest"    ]] && { echo "[MAP] --manifest is required"    >&2; exit 1; }
[[ -z "$scripts_dir" ]] && { echo "[MAP] --scripts-dir is required" >&2; exit 1; }

PYTHON="${PYTHON:-python3}"
PARSE_PY="${scripts_dir}/parse_bmge_track.py"

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        wait -n 2>/dev/null || true
    done
}

# Manifest columns: gene_id <TAB> track_html <TAB> cds_fasta
while IFS=$'\t' read -r gene_id track_html cds_fasta; do
    [[ -z "$gene_id" || "$gene_id" == \#* ]] && continue
    [[ ! -f "$track_html" ]] && { echo "[MAP] WARN: $track_html not found, skipping $gene_id" >&2; continue; }
    [[ ! -f "$cds_fasta"  ]] && { echo "[MAP] WARN: $cds_fasta not found, skipping $gene_id" >&2; continue; }

    wait_for_slot

    (
        out_tsv="${gene_id}.map.tsv"
        args=(
            --track_html  "$track_html"
            --cds_fasta   "$cds_fasta"
            --out_tsv     "$out_tsv"
        )
        [[ -n "$gff"        ]] && args+=(--gff             "$gff")
        [[ -n "$gene_equiv" ]] && args+=(--gene_equiv_tsv  "$gene_equiv")

        if "$PYTHON" "$PARSE_PY" "${args[@]}" 2>&1; then
            [[ ! -s "$out_tsv" ]] && rm -f "$out_tsv"
            echo "[MAP] OK: $gene_id"
        else
            echo "[MAP] WARN: map failed for $gene_id" >&2
            rm -f "$out_tsv"
        fi
    ) &

done < "$manifest"

while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
    wait -n 2>/dev/null || true
done

echo "[MAP] Batch complete"
