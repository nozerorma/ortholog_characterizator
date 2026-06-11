#!/usr/bin/env bash
# run.sh — launcher for the generate-alis Nextflow pipeline.
#
# Edit the CONFIGURATION block below, then:
#   bash nextflow_rewrite/run.sh
#
# The script resolves all input paths relative to the repository root so it
# can be run from anywhere.  Every parameter maps 1-to-1 to a Nextflow param;
# see main.nf --help for the full list.
#
# ── Workflow stages ────────────────────────────────────────────────────────────
# The pipeline is split into independent stages that can be toggled on/off:
#
#   quality            deduplication + sequence filtering (chi2 or threshold)
#   translation        BMGE codon-aware translation + optional trimming
#   variability        Valdar entropy per position + optional Pfam domain mapping
#   phylogeny          per-gene trees (IQ-TREE or CodeML fixed-topology)
#   positive_selection HyPhy FUBAR/FEL/MEME per gene
#
# Skipped stages require a fallback input dir (--dedup_dir / --prot_dir).
#
# ── Filter modes ──────────────────────────────────────────────────────────────
# filter_mode=chi2       iterative chi-square outlier removal (default)
#                        removes sequences with significantly more ambiguous
#                        codons than the rest; up to chi2_max_iter rounds
# filter_mode=threshold  simple coverage cutoff, no statistics
#                        drops sequences where ambiguity > (1 - threshold_min_seq_frac)
#                        e.g. threshold_min_seq_frac=0.2 → drop if >80% gaps/N
# filter_mode=none       skip filtering entirely (dedup output goes straight to translation)
#
# ── Translation (BMGE) ────────────────────────────────────────────────────────
# BMGE v1.12 is the only translation path. It handles codon-aware trimming and
# projects retained codons to amino acids. Trimming is controlled by:
#   bmge_gap_threshold     (default 0.0 — only all-gap columns removed)
#   bmge_entropy_threshold (default 1.0 — no entropy trimming; use 0.5 for tree-only)

set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION  (edit before running)
# ══════════════════════════════════════════════════════════════════════════════

# ── Input CDS alignments ──────────────────────────────────────────────────────
CDS_DIR="/path/to/cds_alignments"          # directory of codon MSA .fa files

# ── Reference / metadata files ───────────────────────────────────────────────
ASSEMBLIES_TSV="${ROOT}/taxa_info/assemblies_and_species_mammals.tsv"
SPECIES_DEDUP_MAP="${ROOT}/taxa_info/species_dedup_map.tsv"
TAXID_TSV="${ROOT}/taxa_info/taxid_species_family_mammals.tsv"
FAMILY_ORDER_TSV="${ROOT}/taxa_info/family_to_order.tsv"

# ── Output ────────────────────────────────────────────────────────────────────
OUT_DIR="${ROOT}/results/$(date +%Y%m%d_%H%M%S)"

# ── Workflow stages to run ────────────────────────────────────────────────────
RUN_QUALITY=true
RUN_TRANSLATION=true
RUN_VARIABILITY=false
RUN_PHYLOGENY=false
RUN_POSITIVE_SELECTION=false

# ── Filter mode: chi2 | threshold | none ─────────────────────────────────────
FILTER_MODE="chi2"

# ── Optional: restrict to a subset of genes ───────────────────────────────────
GENE_LIST=""    # path to plain-text file (one gene name per line); leave empty to run all

# ══════════════════════════════════════════════════════════════════════════════

NF_DIR="${ROOT}/nextflow_rewrite"

# Build optional args
OPTIONAL_ARGS=()
[[ -n "$GENE_LIST"       ]] && OPTIONAL_ARGS+=(--gene_list        "$GENE_LIST")
[[ -n "$FAMILY_ORDER_TSV" && -f "$FAMILY_ORDER_TSV" ]] \
                              && OPTIONAL_ARGS+=(--family_order_tsv "$FAMILY_ORDER_TSV")

echo "════════════════════════════════════════════════════════"
echo " generate-alis Nextflow pipeline"
echo " CDS dir     : ${CDS_DIR}"
echo " Output      : ${OUT_DIR}"
echo " Stages      : quality=${RUN_QUALITY}  translation=${RUN_TRANSLATION}"
echo "               variability=${RUN_VARIABILITY}  phylogeny=${RUN_PHYLOGENY}"
echo "               positive_selection=${RUN_POSITIVE_SELECTION}"
echo " Filter mode : ${FILTER_MODE}"
echo "════════════════════════════════════════════════════════"
echo ""

RUN_LABEL="$(date +%Y%m%d_%H%M%S)"
[[ "$RUN_QUALITY"            == "true" ]] && RUN_LABEL="${RUN_LABEL}_quality"
[[ "$RUN_TRANSLATION"        == "true" ]] && RUN_LABEL="${RUN_LABEL}_translation"
[[ "$RUN_VARIABILITY"        == "true" ]] && RUN_LABEL="${RUN_LABEL}_variability"
[[ "$RUN_PHYLOGENY"          == "true" ]] && RUN_LABEL="${RUN_LABEL}_phylogeny"
[[ "$RUN_POSITIVE_SELECTION" == "true" ]] && RUN_LABEL="${RUN_LABEL}_psel"

nextflow run "${NF_DIR}/main.nf" \
    -name "$RUN_LABEL" \
    -with-tower \
    --cds_dir               "$CDS_DIR" \
    --out_dir               "$OUT_DIR" \
    --assemblies_tsv        "$ASSEMBLIES_TSV" \
    --species_dedup_map     "$SPECIES_DEDUP_MAP" \
    --taxid_tsv             "$TAXID_TSV" \
    --run_quality           "$RUN_QUALITY" \
    --run_translation       "$RUN_TRANSLATION" \
    --run_variability       "$RUN_VARIABILITY" \
    --run_phylogeny         "$RUN_PHYLOGENY" \
    --run_positive_selection "$RUN_POSITIVE_SELECTION" \
    --filter_mode           "$FILTER_MODE" \
    "${OPTIONAL_ARGS[@]}" \
    -profile local \
    -resume
