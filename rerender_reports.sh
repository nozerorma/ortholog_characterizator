#!/usr/bin/env bash
# rerender_reports.sh — re-render all three HTML reports from existing pipeline outputs.
#
# Usage:
#   bash rerender_reports.sh --results_dir /path/to/results/20260517_235350 [options]
#
# Options:
#   --results_dir   PATH   Required. Pipeline output directory.
#   --psel_tree     PATH   Species tree (Newick). Required for phylogenetic trees.
#   --taxid_tsv     PATH   taxid_species_family.tsv (default: auto-detect in taxa_info/)
#   --assemblies    PATH   assemblies_and_species.tsv (default: auto-detect)
#   --pfam_desc     PATH   pfam_descriptions.tsv (default: auto-detect)
#   --family_order_tsv PATH  family_to_order.tsv (default: auto-detect)
#   --sel_var_tsv   PATH   sel_variability.tsv from positive-selection run (optional)
#   --out_dir       PATH   Where to write the HTML files (default: results_dir/reports)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Defaults ──────────────────────────────────────────────────────────────────
RESULTS_DIR=""
PSEL_TREE=""
OUT_REPORT_DIR=""
SEL_VAR_TSV=""

# Auto-detect species tree
REPO_ROOT_EARLY="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
for _tcand in \
    "${REPO_ROOT_EARLY}/external_phylogenies/speciesTree_speciesname_pruned.nh" \
    "${REPO_ROOT_EARLY}/external_phylogenies/speciesTree_speciesname_pruned.nwk"; do
    [[ -f "$_tcand" ]] && { PSEL_TREE="$_tcand"; break; }
done

# Auto-detect shared resource paths relative to the pipeline root
REPO_ROOT="$(dirname "$SCRIPT_DIR")"   # one level up from nextflow_rewrite/
TAXA_DIR="${REPO_ROOT}/taxa_info"
TAXID_TSV="${TAXA_DIR}/taxid_species_family_mammals.tsv"
ASSEMBLIES_TSV=""
for _cand in \
    "${TAXA_DIR}/assemblies_and_species_mammals.tsv" \
    "${TAXA_DIR}/assemblies_and_species_mammals_merged.tsv" \
    "${TAXA_DIR}/assemblies_and_species.tsv"; do
    [[ -f "$_cand" ]] && { ASSEMBLIES_TSV="$_cand"; break; }
done
PFAM_DESC_TSV="${TAXA_DIR}/pfam_descriptions.tsv"
FAMILY_ORDER_TSV="${TAXA_DIR}/family_to_order.tsv"

# ── Argument parsing ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --results_dir)      RESULTS_DIR="$2";       shift 2 ;;
        --psel_tree)        PSEL_TREE="$2";         shift 2 ;;
        --taxid_tsv)        TAXID_TSV="$2";         shift 2 ;;
        --assemblies)       ASSEMBLIES_TSV="$2";    shift 2 ;;
        --pfam_desc)        PFAM_DESC_TSV="$2";     shift 2 ;;
        --family_order_tsv) FAMILY_ORDER_TSV="$2";  shift 2 ;;
        --sel_var_tsv)      SEL_VAR_TSV="$2";       shift 2 ;;
        --out_dir)          OUT_REPORT_DIR="$2";    shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

[[ -z "$RESULTS_DIR" ]] && { echo "ERROR: --results_dir is required"; exit 1; }
[[ ! -d "$RESULTS_DIR" ]] && { echo "ERROR: $RESULTS_DIR does not exist"; exit 1; }

OUT_REPORT_DIR="${OUT_REPORT_DIR:-${RESULTS_DIR}/reports}"
mkdir -p "$OUT_REPORT_DIR"

WORK_DIR="${OUT_REPORT_DIR}/.rerender_tmp"
rm -rf "$WORK_DIR"
mkdir -p "$WORK_DIR"
trap 'rm -rf "$WORK_DIR"' EXIT

# ── Auto-detect filter mode ───────────────────────────────────────────────────
_FILTER_MODE="none"
_PREFILTER_DIR=""
for _fm in threshold chi2; do
    if [[ -d "${RESULTS_DIR}/PREFILTER/${_fm}" ]]; then
        _FILTER_MODE="$_fm"
        _PREFILTER_DIR="${RESULTS_DIR}/PREFILTER/${_fm}"
        break
    fi
done

# ── Auto-detect translation mode and TRIM dir ─────────────────────────────────
_TRANSLATION_MODE="manual"
_TRIM_DIR=""
if [[ -d "${RESULTS_DIR}/TRIM/bmge" ]]; then
    _TRANSLATION_MODE="bmge"
    _TRIM_DIR="${RESULTS_DIR}/TRIM/bmge"
elif [[ -d "${RESULTS_DIR}/TRIM/manual" ]]; then
    _TRANSLATION_MODE="manual"
    _TRIM_DIR="${RESULTS_DIR}/TRIM/manual"
fi

# ── Auto-detect PROT dir ──────────────────────────────────────────────────────
_PROT_DIR="${RESULTS_DIR}/PROT/${_TRANSLATION_MODE}"
[[ ! -d "$_PROT_DIR" ]] && _PROT_DIR="${RESULTS_DIR}/PROT"

# ── Auto-detect sel_variability.tsv if not supplied ───────────────────────────
if [[ -z "$SEL_VAR_TSV" ]]; then
    _sv_cand="${RESULTS_DIR}/reports/positive_selection/sel_variability.tsv"
    [[ -f "$_sv_cand" ]] && SEL_VAR_TSV="$_sv_cand"
fi

echo "Results dir      : $RESULTS_DIR"
echo "Output dir       : $OUT_REPORT_DIR"
echo "Work dir         : $WORK_DIR"
echo "Assemblies       : ${ASSEMBLIES_TSV:-NOT FOUND (taxonomy/assembly plots will be skipped)}"
echo "Taxid TSV        : ${TAXID_TSV}"
echo "Family→Order     : ${FAMILY_ORDER_TSV:-NOT FOUND (order breakdown will be skipped)}"
echo "Species tree     : ${PSEL_TREE:-NOT FOUND (phylogenetic trees will be skipped)}"
echo "Translation mode : $_TRANSLATION_MODE"
echo "Filter mode      : $_FILTER_MODE"
echo "Trim dir         : ${_TRIM_DIR:-not found}"
echo "Sel-var TSV      : ${SEL_VAR_TSV:-not found (psel section in variability report will be skipped)}"
echo ""

# ── 1. Post-translation report ─────────────────────────────────────────────────
echo ">>> Rendering post_translation_report ..."
mkdir -p "$WORK_DIR/post_translation_report"

QUALITY_TSV="${RESULTS_DIR}/reports/translation/quality/quality_data_quality_per_gene.csv"
METADATA_TSV="${RESULTS_DIR}/reports/translation/metadata/metadata.combined_summary.tsv"

(
    cd "$WORK_DIR/post_translation_report"
    cp "${SCRIPT_DIR}/subworkflows/translation/local/post_translation_report.qmd" .
    [[ -f "$QUALITY_TSV"  ]] && cp "$QUALITY_TSV"  . || true
    [[ -f "$METADATA_TSV" ]] && cp "$METADATA_TSV" . || true

    PROT_DIR="$_PROT_DIR" \
    QUALITY_TSV="$(basename "$QUALITY_TSV")" \
    METADATA_TSV="$(basename "$METADATA_TSV")" \
    ASSEMBLIES_TSV="$ASSEMBLIES_TSV" \
    TAXID_TSV="$TAXID_TSV" \
    FAMILY_ORDER_TSV="$FAMILY_ORDER_TSV" \
    DEDUP_DIR="${RESULTS_DIR}/DEDUP" \
    PREFILTER_DIR="$_PREFILTER_DIR" \
    FILTER_MODE="$_FILTER_MODE" \
    TRIM_DIR="$_TRIM_DIR" \
    TRANSLATION_MODE="$_TRANSLATION_MODE" \
    SPECIES_TREE="${PSEL_TREE}" \
    quarto render post_translation_report.qmd \
        --output post_translation_report.html \
        --output-dir .
)
cp "$WORK_DIR/post_translation_report/post_translation_report.html" "$OUT_REPORT_DIR/"
echo "    → $OUT_REPORT_DIR/post_translation_report.html"

# ── 2. Variability report ──────────────────────────────────────────────────────
echo ">>> Rendering variability_report ..."
mkdir -p "$WORK_DIR/variability_report"

VAR_DIR="${RESULTS_DIR}/reports/variability"
GENE_VAR_TSV="${VAR_DIR}/gene_variability.tsv"
CLADE_VAR_TSV="${VAR_DIR}/clade_variability.tsv"
CLADE_OUT_TSV="${VAR_DIR}/clade_outliers.tsv"
DOMAIN_TSV="${VAR_DIR}/domain_variability.tsv"
PROT_VAR_RAW="${RESULTS_DIR}/PROT_VAR_RAW"

(
    cd "$WORK_DIR/variability_report"
    cp "${SCRIPT_DIR}/subworkflows/variability/local/variability_report.qmd" .
    ln -sf "${PROT_VAR_RAW}"/*.clade_entropy.tsv . 2>/dev/null || true

    GENE_VARIABILITY_TSV="$GENE_VAR_TSV" \
    CLADE_VARIABILITY_TSV="$CLADE_VAR_TSV" \
    CLADE_OUTLIERS_TSV="$CLADE_OUT_TSV" \
    DOMAIN_TSV="${DOMAIN_TSV}" \
    PFAM_DESC_TSV="$PFAM_DESC_TSV" \
    ENTROPY_DIR="$(pwd)" \
    TAXID_TSV="$TAXID_TSV" \
    FAMILY_ORDER_TSV="$FAMILY_ORDER_TSV" \
    SPECIES_TREE="${PSEL_TREE}" \
    SEL_VARIABILITY_TSV="${SEL_VAR_TSV}" \
    quarto render variability_report.qmd \
        --output variability_report.html \
        --output-dir .
)
cp "$WORK_DIR/variability_report/variability_report.html" "$OUT_REPORT_DIR/"
echo "    → $OUT_REPORT_DIR/variability_report.html"

# ── 3. Positive selection report ───────────────────────────────────────────────
echo ">>> Rendering positive_selection_report ..."
mkdir -p "$WORK_DIR/positive_selection_report"

PSEL_DIR="${RESULTS_DIR}/reports/positive_selection"
FUBAR_TSV="${PSEL_DIR}/fubar_results.tsv"
FUBAR_SITES_TSV="${PSEL_DIR}/fubar_sites.tsv"
FEL_TSV="${PSEL_DIR}/fel_results.tsv"
MEME_TSV="${PSEL_DIR}/meme_results.tsv"

(
    cd "$WORK_DIR/positive_selection_report"
    cp "${SCRIPT_DIR}/subworkflows/positive_selection/local/positive_selection_report.qmd" .

    FUBAR_TSV="$FUBAR_TSV" \
    FUBAR_SITES_TSV="$FUBAR_SITES_TSV" \
    FEL_TSV="${FEL_TSV}" \
    MEME_TSV="${MEME_TSV}" \
    GENE_VARIABILITY_TSV="${GENE_VAR_TSV}" \
    TAXID_TSV="$TAXID_TSV" \
    FAMILY_ORDER_TSV="$FAMILY_ORDER_TSV" \
    SPECIES_TREE="${PSEL_TREE}" \
    HYPHY_METHODS="FUBAR" \
    FUBAR_THRESHOLD="0.9" \
    quarto render positive_selection_report.qmd \
        --output positive_selection_report.html \
        --output-dir .
)
cp "$WORK_DIR/positive_selection_report/positive_selection_report.html" "$OUT_REPORT_DIR/"
echo "    → $OUT_REPORT_DIR/positive_selection_report.html"

echo ""
echo "Done. Reports written to: $OUT_REPORT_DIR"
