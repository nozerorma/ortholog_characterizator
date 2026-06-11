#!/bin/bash
# Per-gene HyPhy positive selection: FUBAR (site-level Bayesian) + FEL (site-level ML).
# MEME is run as a post-hoc step on FUBAR hits only (see 9.run_positive_selection.sh).
# Resume-safe: skips analyses whose JSON output already exists.
#
# Usage:
#   bash scripts/local/run_hyphy_gene.sh <cds_ali> <out_dir> <species_tree> [analyses]
#
# Arguments:
#   cds_ali      Codon CDS alignment FASTA (DEDUP/<gene>.fa)
#   out_dir      Base output directory (PSEL/ subdirs created here)
#   species_tree Species tree Newick (speciesTree_speciesname_pruned.nh)
#   analyses     Comma-separated subset: FUBAR,FEL,MEME,SLAC (default: FUBAR,FEL)

set -euo pipefail

CDS_ALI=${1}
OUT_DIR=${2}
SPECIES_TREE=${3}
ANALYSES=${4:-"FUBAR"}

gene=$(basename "$CDS_ALI" .fa)
PSEL_DIR="${OUT_DIR}/PSEL/${gene}"
TREE_DIR="${OUT_DIR}/PSEL_TREES"

mkdir -p "$PSEL_DIR" "$TREE_DIR"

# ── Prune species tree to alignment taxa ──────────────────────────────────────
PRUNED_TREE="${TREE_DIR}/${gene}.nwk"
if [ ! -f "$PRUNED_TREE" ]; then
    python scripts/local/prune_tree_for_gene.py \
        --tree "$SPECIES_TREE" \
        --fasta "$CDS_ALI" \
        --out "$PRUNED_TREE" || {
            echo "WARN ${gene}: tree pruning failed — skipping positive selection" >&2
            exit 0
        }
fi

run_analysis() {
    local method=$1
    local out_json="${PSEL_DIR}/${gene}.${method}.json"

    [ -f "$out_json" ] && { echo "SKIP ${gene} ${method}: output exists"; return 0; }

    echo "[$(date +'%H:%M:%S')] ${gene}: running ${method}..."

    # Cap BLAS/OpenMP threads to the per-worker CPU allocation.
    # HyPhy's CPU argument only controls its own scheduler; library-level
    # threads (OpenBLAS, OpenMP, MKL) detect hardware CPU count via sysconf()
    # and ignore CPU entirely unless these vars are set.
    export OMP_NUM_THREADS=2
    export MKL_NUM_THREADS=2
    export OPENBLAS_NUM_THREADS=2
    export BLAS_NUM_THREADS=2

    case "$method" in
        FUBAR)
            conda run -n bmge-tools hyphy FUBAR \
                --alignment "$CDS_ALI" \
                --tree "$PRUNED_TREE" \
                --output "$out_json" \
                --grid 20 \
                CPU=2 \
                > "${PSEL_DIR}/${gene}.FUBAR.log" 2>&1
            ;;
        FEL)
            conda run -n bmge-tools hyphy FEL \
                --alignment "$CDS_ALI" \
                --tree "$PRUNED_TREE" \
                --output "$out_json" \
                --branches All \
                CPU=2 \
                > "${PSEL_DIR}/${gene}.FEL.log" 2>&1
            ;;
        MEME)
            conda run -n bmge-tools hyphy MEME \
                --alignment "$CDS_ALI" \
                --tree "$PRUNED_TREE" \
                --output "$out_json" \
                --branches All \
                CPU=2 \
                > "${PSEL_DIR}/${gene}.MEME.log" 2>&1
            ;;
        SLAC)
            conda run -n bmge-tools hyphy SLAC \
                --alignment "$CDS_ALI" \
                --tree "$PRUNED_TREE" \
                --output "$out_json" \
                --branches All \
                --samples 0 \
                CPU=2 \
                > "${PSEL_DIR}/${gene}.SLAC.log" 2>&1
            ;;
        *)
            echo "WARN ${gene}: unknown analysis '${method}'" >&2
            ;;
    esac

    if [ $? -eq 0 ] && [ -s "$out_json" ]; then
        echo "[$(date +'%H:%M:%S')] ${gene}: ${method} done → ${out_json}"
    else
        echo "ERROR ${gene}: ${method} failed" >&2
    fi
}

IFS=',' read -ra METHODS <<< "$ANALYSES"
for method in "${METHODS[@]}"; do
    run_analysis "$method"
done
