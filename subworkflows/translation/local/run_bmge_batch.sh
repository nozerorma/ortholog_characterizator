#!/usr/bin/env bash
# Batch BMGE trimming + protein projection wrapper.
# For each gene in the manifest: runs BMGE v1.12 in CODON mode.
# Handles the MST timeout bug (exit 124): appends a sentinel sequence,
# re-runs, then strips the sentinel from all outputs.

set -uo pipefail

manifest=""
workers=6
gap_thresh=0.0
entropy_thresh=1.0
matrix="BLOSUM62"
gff=""
gene_equiv=""
scripts_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)       manifest="$2";       shift 2 ;;
        --workers)        workers="$2";        shift 2 ;;
        --gap-thresh)     gap_thresh="$2";     shift 2 ;;
        --entropy-thresh) entropy_thresh="$2"; shift 2 ;;
        --matrix)         matrix="$2";         shift 2 ;;
        --gff)            gff="$2";            shift 2 ;;
        --gene-equiv)     gene_equiv="$2";     shift 2 ;;
        --scripts-dir)    scripts_dir="$2";    shift 2 ;;
        *) echo "[BMGE] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$manifest"    ]] && { echo "[BMGE] --manifest is required"    >&2; exit 1; }
[[ -z "$scripts_dir" ]] && { echo "[BMGE] --scripts-dir is required" >&2; exit 1; }

PYTHON="${PYTHON:-python3}"
PARSE_TRACK="${scripts_dir}/parse_bmge_track.py"
BACKTRANSLATE="${scripts_dir}/backtranslate_from_prot.py"
BMGE_TIMEOUT=300   # seconds before assuming MST bug

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        wait -n 2>/dev/null || true
    done
}

run_bmge_for_gene() {
    local gene_id="$1"
    local fasta="$2"
    local trim_out="${gene_id}.trim.fa"
    local track_out="${gene_id}.track.html"
    local prot_out="${gene_id}.prot.fa"
    local tmp_fa="${gene_id}.bmge_tmp.fa"

    # Helper: strip sentinel from a FASTA file
    strip_sentinel() {
        local f="$1"
        [[ ! -f "$f" ]] && return
        "$PYTHON" - "$f" <<'PYEOF'
import sys
path = sys.argv[1]
lines = open(path).readlines()
out = []
skip = False
for line in lines:
    if line.startswith('>') and 'BMGE_SENTINEL_DUPE' in line:
        skip = True
    elif line.startswith('>'):
        skip = False
    if not skip:
        out.append(line)
with open(path, 'w') as fh:
    fh.writelines(out)
PYEOF
    }

    # Run BMGE (with timeout guard).
    # NOTE: -of is intentionally omitted — BMGE 1.12 (bioconda) has a known
    # bug where -t CODON truncates the FASTA nucleotide output to protein
    # length. The correct CDS is reconstructed via back-translation below.
    run_once() {
        local in_fa="$1"
        timeout "$BMGE_TIMEOUT" bmge \
            -i    "$in_fa" \
            -t    CODON \
            -m    "$matrix" \
            -g    "$gap_thresh" \
            -h    "$entropy_thresh" \
            -oh   "$track_out" \
            -oaaf "$prot_out"
    }

    cp "$fasta" "$tmp_fa"
    run_once "$tmp_fa"
    exit_code=$?

    if [[ $exit_code -eq 124 ]]; then
        # MST bug workaround: append duplicate of first sequence with sentinel tag
        echo "[BMGE] WARN: timeout for $gene_id, applying MST workaround" >&2
        first_seq=$(awk '/^>/{if(p)exit; p=1} p' "$tmp_fa")
        first_hdr=$(grep '^>' "$tmp_fa" | head -1 | sed 's/^>//')
        printf ">%s_BMGE_SENTINEL_DUPE\n%s\n" "$first_hdr" "$first_seq" >> "$tmp_fa"

        run_once "$tmp_fa"
        exit_code=$?

        # Strip sentinel from all outputs regardless of exit code
        strip_sentinel "$prot_out"
        # track HTML: strip sentinel entry (best-effort)
        [[ -f "$track_out" ]] && sed -i '/BMGE_SENTINEL_DUPE/d' "$track_out" || true
    fi

    rm -f "$tmp_fa"

    if [[ $exit_code -ne 0 ]]; then
        echo "[BMGE] WARN: BMGE failed for $gene_id (exit $exit_code), skipping" >&2
        rm -f "$trim_out" "$track_out" "$prot_out"
        return
    fi

    # Drop empty outputs
    for f in "$track_out" "$prot_out"; do
        [[ ! -s "$f" ]] && rm -f "$f"
    done

    # Back-translate protein alignment → correct CDS alignment.
    # BMGE's -of is broken in CODON mode (truncates to protein length);
    # we reconstruct the CDS by mapping the protein alignment back to
    # the original gap-stripped codons.
    if [[ -s "$prot_out" ]]; then
        if ! "$PYTHON" "$BACKTRANSLATE" \
                --prot "$prot_out" \
                --cds  "$fasta" \
                --out  "$trim_out" 2>&1; then
            echo "[BMGE] WARN: back-translation failed for $gene_id, skipping" >&2
            rm -f "$trim_out" "$track_out" "$prot_out"
            return
        fi
        [[ ! -s "$trim_out" ]] && rm -f "$trim_out"
    fi

    echo "[BMGE] OK: $gene_id"
}

while IFS=$'\t' read -r gene_id fasta_name; do
    [[ -z "$gene_id" || "$gene_id" == \#* ]] && continue
    [[ ! -f "$fasta_name" ]] && { echo "[BMGE] WARN: $fasta_name not found, skipping $gene_id" >&2; continue; }

    wait_for_slot
    ( run_bmge_for_gene "$gene_id" "$fasta_name" ) &

done < "$manifest"

while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
    wait -n 2>/dev/null || true
done

echo "[BMGE] Batch complete"
