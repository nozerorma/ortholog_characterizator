#!/bin/bash
#
# generate-alis Nextflow Pipeline SLURM Launcher — PHYLOGENY
#
# ─── PURPOSE ──────────────────────────────────────────────────────────────────
# Runs quality + translation (with BMGE trimming) + gene trees + species tree.
# BMGE trimming is enabled here (-g 0.2 -h 0.5) because high-entropy and
# highly-gapped columns add noise to ML tree inference. Do NOT use these
# settings for variability, positive selection, or CAAS (use submit_bmge.sh).
#
# ─── USAGE ────────────────────────────────────────────────────────────────────
#   bash nextflow_rewrite/submit_phylo.sh
#

mkdir -p Slurm

# ============================================================
# CONFIGURATION  (edit before submitting)
# ============================================================
ROOT="/data/samanthafs/scratch/lab_anavarro/mramon/4.Generate_alignments_from_codons"
NF_DIR="${ROOT}/nextflow_rewrite"
SLURM_LOG_DIR="${NF_DIR}/Slurm"

# ── Input CDS alignments ──────────────────────────────────────────────────────
CDS_DIR="${ROOT}/alignments/Mammals_cds"

# ── Reference / metadata files ───────────────────────────────────────────────
ASSEMBLIES_TSV="${ROOT}/taxa_info/assemblies_and_species_mammals.tsv"
SPECIES_DEDUP_MAP="${ROOT}/taxa_info/species_dedup_map.tsv"
TAXID_TSV="${ROOT}/taxa_info/taxid_species_family_mammals.tsv"
FAMILY_ORDER_TSV="${ROOT}/taxa_info/family_to_order.tsv"

# ── Output ────────────────────────────────────────────────────────────────────
OUT_DIR="${ROOT}/results_phylo/$(date +%Y%m%d_%H%M%S)"
WORK_DIR="/data/samanthafs/scratch/lab_anavarro/mramon/3.Work_dirs/generate_alis_nf/work"

# ── Pre-computed inputs (optional) ────────────────────────────────────────────
# Set PROT_DIR to skip quality + translation and feed existing protein
# alignments directly into the phylogeny workflow.
#
# Example (reuse yesterday's translation run):
#   PROT_DIR="${ROOT}/results/20250101_120000/PROT/bmge"
PROT_DIR=""

# ── Workflow stages ───────────────────────────────────────────────────────────
RUN_QUALITY=true
RUN_TRANSLATION=true
RUN_VARIABILITY=false
RUN_PHYLOGENY=true
RUN_POSITIVE_SELECTION=false
RUN_UCR=false             # UCR analysis requires RUN_VARIABILITY=true; keep false here

# Auto-disable quality + translation when pre-computed protein alignments are provided
if [[ -n "$PROT_DIR" ]]; then
    RUN_QUALITY=false
    RUN_TRANSLATION=false
fi

# ── Quality filter ───────────────────────────────────────────────────────────
FILTER_MODE="chi2"            # chi2 | threshold | none
# Chi2 filter params (leave empty to use defaults from quality.config)
CHI2_MAX_ITER="10"              # max iterative removal rounds         (default: 3)
CHI2_AMBIGUITY_THRESH="0.8"  # max ambiguous-codon fraction across whole alignment  (default: 0.8)
CHI2_MIN_SEQ_FRAC="0.1"      # min fraction of sequences that must survive filtering (default: 0.1)
# Threshold filter param (used only when FILTER_MODE=threshold)
THRESHOLD_MIN_SEQ_FRAC="0.2"     # min non-gap/N fraction/seq           (default: 0.2)
MIN_SEQS=""                   # drop genes with fewer seqs after filter (default: 4)

# ── BMGE trimming (active for tree inference) ─────────────────────────────────
# Removes gapped and high-entropy columns that add noise to ML inference.
# These are BMGE paper defaults — do not use for CAAS/variability/psel.
BMGE_GAP_THRESH="0.2"         # -g: max gap rate allowed per column; 0.2 = remove columns with >20% gaps
BMGE_ENTROPY_THRESH="0.5"     # -h: remove columns with entropy >0.5 (BMGE default)

# ── Test mode ─────────────────────────────────────────────────────────────────
TEST_MODE=true
TEST_N_GENES=100

# ── Optional: restrict to a gene list ────────────────────────────────────────
GENE_LIST=""

# ── IQ-TREE ───────────────────────────────────────────────────────────────────
IQTREE_MSET="LG"              # restrict ModelFinder to this model family; "" for all
SPECIES_TREE_METHOD="astral"  # astral | backbone

# ── Pre-build optional Nextflow args ─────────────────────────────────────────
NF_OPTIONAL=""
[[ -n "$GENE_LIST"              ]] && NF_OPTIONAL+=" --gene_list $GENE_LIST"
[[ -n "$CHI2_MAX_ITER"          ]] && NF_OPTIONAL+=" --chi2_max_iter $CHI2_MAX_ITER"
[[ -n "$CHI2_AMBIGUITY_THRESH"  ]] && NF_OPTIONAL+=" --chi2_ambiguity_thresh $CHI2_AMBIGUITY_THRESH"
[[ -n "$CHI2_MIN_SEQ_FRAC"      ]] && NF_OPTIONAL+=" --chi2_min_seq_frac $CHI2_MIN_SEQ_FRAC"
[[ -n "$THRESHOLD_MIN_SEQ_FRAC" ]] && NF_OPTIONAL+=" --threshold_min_seq_frac $THRESHOLD_MIN_SEQ_FRAC"
[[ -n "$MIN_SEQS"               ]] && NF_OPTIONAL+=" --min_seqs $MIN_SEQS"
[[ -f "$FAMILY_ORDER_TSV"    ]] && NF_OPTIONAL+=" --family_order_tsv $FAMILY_ORDER_TSV"
[[ "$TEST_MODE" == "true"    ]] && NF_OPTIONAL+=" --test_mode true --test_n_genes $TEST_N_GENES"
[[ -n "$IQTREE_MSET"         ]] && NF_OPTIONAL+=" --iqtree_mset $IQTREE_MSET"
[[ -n "$SPECIES_TREE_METHOD" ]] && NF_OPTIONAL+=" --species_tree_method $SPECIES_TREE_METHOD"
# Pre-computed protein alignments
[[ -n "$PROT_DIR"            ]] && NF_OPTIONAL+=" --prot_dir $PROT_DIR"

# ============================================================

submit_job() {
    sbatch --parsable <<EOF
#!/bin/bash
#SBATCH --job-name=generate-alis-phylo
#SBATCH --partition=haswell
#SBATCH -t 96:00:00
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH -e ${SLURM_LOG_DIR}/slurm-%j.err
#SBATCH -o ${SLURM_LOG_DIR}/slurm-%j.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=miguel.ramon@upf.edu

module purge; module load modulepath/haswell
module load Miniconda3

source ~/.bashrc
conda deactivate
conda activate bmge

module load Nextflow

echo "======================================================"
echo " JOB ID  : \${SLURM_JOB_ID}"
echo " NODE    : \$(hostname)"
echo " STARTED : \$(date)"
echo "======================================================"

nextflow run "${NF_DIR}/main.nf" \\
    -name "generate-alis-phylo-$(date +%Y%m%d-%H%M%S)" \\
    -with-tower \\
    --cds_dir                "${CDS_DIR}" \\
    --out_dir                "${OUT_DIR}" \\
    --assemblies_tsv         "${ASSEMBLIES_TSV}" \\
    --species_dedup_map      "${SPECIES_DEDUP_MAP}" \\
    --taxid_tsv              "${TAXID_TSV}" \\
    --run_quality            "${RUN_QUALITY}" \\
    --run_translation        "${RUN_TRANSLATION}" \\
    --run_variability        "${RUN_VARIABILITY}" \\
    --run_phylogeny          "${RUN_PHYLOGENY}" \\
    --run_positive_selection "${RUN_POSITIVE_SELECTION}" \\
    --run_ucr                "${RUN_UCR}" \\
    --filter_mode            "${FILTER_MODE}" \\
    --bmge_gap_threshold     "${BMGE_GAP_THRESH}" \\
    --bmge_entropy_threshold "${BMGE_ENTROPY_THRESH}" \\
    ${NF_OPTIONAL} \\
    -profile slurm \\
    -work-dir "${WORK_DIR}" \\
    -resume

echo "[\$(date)] Pipeline complete"
EOF
}

# Main
job_id=$(submit_job)
echo "Submitted job : ${job_id}"
echo ""
echo "Monitor with  : squeue -u \$USER"
echo "Logs in       : ${SLURM_LOG_DIR}/"

mkdir -p SBATCH_runner
mv slurm-*.out SBATCH_runner 2>/dev/null || true
mv slurm-*.err SBATCH_runner 2>/dev/null || true
