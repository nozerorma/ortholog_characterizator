#!/bin/bash
#
# generate-alis Nextflow Pipeline SLURM Launcher
#
# ─── PURPOSE ──────────────────────────────────────────────────────────────────
# Submits a single SLURM job that runs the Nextflow master process.
# The master process submits all per-gene child jobs via Nextflow's SLURM
# executor (-profile slurm). Only this one job needs to be submitted manually.
# Resources here are for the orchestrator only; heavy compute is in
# conf/resources.config.
#
# ─── USAGE ────────────────────────────────────────────────────────────────────
#   bash nextflow_rewrite/submit_bmge.sh
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)
# File: nextflow_rewrite/submit_bmge.sh

mkdir -p Slurm

# ============================================================
# CONFIGURATION  (edit before submitting)
# ============================================================
ROOT="/data/samanthafs/scratch/lab_anavarro/mramon/4.Generate_alignments_from_codons"
NF_DIR="${ROOT}/nextflow_rewrite"
SLURM_LOG_DIR="${NF_DIR}/Slurm"

# ── Input CDS alignments ──────────────────────────────────────────────────────
CDS_DIR="${ROOT}/alignments/Primates_cds"

# ── Reference / metadata files ───────────────────────────────────────────────
# PPGA does not deduplicate assemblies, so the assembly metadata TSV is not required.
ASSEMBLIES_TSV=""
SPECIES_DEDUP_MAP="${ROOT}/taxa_info/species_dedup_map.tsv"
TAXID_TSV="${ROOT}/taxa_info/taxid_species_family_primates.tsv"
FAMILY_ORDER_TSV="${ROOT}/taxa_info/family_to_order.tsv"

# ── Output ────────────────────────────────────────────────────────────────────
OUT_DIR="${ROOT}/results_ppga/$(date +%Y%m%d_%H%M%S)"
WORK_DIR="/data/samanthafs/scratch/lab_anavarro/mramon/3.Work_dirs/generate_alis_nf/work_ppga"

# ── Pre-computed inputs (optional) ────────────────────────────────────────────
# Set any of these to skip quality + translation and feed existing outputs
# directly into variability / positive-selection.
#   PROT_DIR  → protein alignments (PROT/bmge/*.fa)   → variability
#   TRIM_DIR  → BMGE-trimmed CDS (TRIM/bmge/*.fa)     → positive selection
#   MAP_DIR   → BMGE position maps (MAP/*.map.tsv)    → coordinate back-mapping
#   DEDUP_DIR → deduplicated CDS (DEDUP/*.fa)         → re-run translation only
#
# Example (reuse yesterday's translation run):
#   PROT_DIR="${ROOT}/results/20250101_120000/PROT/bmge"
#   TRIM_DIR="${ROOT}/results/20250101_120000/TRIM/bmge"
#   MAP_DIR="${ROOT}/results/20250101_120000/MAP"
PROT_DIR=""
TRIM_DIR=""
MAP_DIR=""
DEDUP_DIR=""

# ── Workflow stages to run ────────────────────────────────────────────────────
RUN_QUALITY=true
RUN_TRANSLATION=true
RUN_VARIABILITY=true
RUN_PHYLOGENY=false
RUN_POSITIVE_SELECTION=true
RUN_UCR=true             # Ultra-Conserved Region analysis (requires RUN_VARIABILITY=true)

# PPGA keeps the raw assembly headers through quality filtering; no assembly-level deduplication.
RUN_DEDUP=false

# Auto-disable quality + translation when pre-computed inputs are provided
if [[ -n "$PROT_DIR" || -n "$TRIM_DIR" || -n "$DEDUP_DIR" ]]; then
    RUN_QUALITY=false
    RUN_TRANSLATION=false
fi

# ── Quality filter ───────────────────────────────────────────────────────────
FILTER_MODE="chi2"            # chi2 | threshold | none
# Chi2 filter params (leave empty to use defaults from quality.config)
CHI2_MAX_ITER="0"              # max iterative removal rounds         (default: 3)
CHI2_AMBIGUITY_THRESH="0.8"  # max ambiguous-codon fraction across whole alignment  (default: 0.8)
CHI2_MIN_SEQ_FRAC="0.1"      # min fraction of sequences that must survive filtering (default: 0.1)
# Threshold filter param (used only when FILTER_MODE=threshold)
THRESHOLD_MIN_SEQ_FRAC="0.2"     # min non-gap/N fraction/seq           (default: 0.2)
MIN_SEQS="10"                   # drop genes with fewer seqs after filter (default: 4)

# ── BMGE trimming ─────────────────────────────────────────────────────────────
# Default: passthrough (no trimming). Trimming removes alignment info that
# positive selection, CAAS and variability analyses depend on.
# For tree-only runs, consider: GAP_THRESH=0.1  ENTROPY_THRESH=0.5
BMGE_GAP_THRESH="0.1"         # -g: max gap rate allowed per column; 0.0 = remove any column with even one gap (maximum trimming)
BMGE_ENTROPY_THRESH="0.5"     # -h: max entropy score; 1.0 = no trimming

# ── Test mode: run on a random subset of TEST_N_GENES genes ──────────────────
TEST_MODE=true
TEST_N_GENES=100

# ── Concurrency: per-gene batch sizing (scales with dataset size) ────────────
# Heavy per-gene stages run batched: each batch is ONE SLURM job that processes
# BATCH_SIZE genes, BATCH_WORKERS at a time in parallel. Small batches in test mode
# → quick turnaround and fine-grained parallelism on a few hundred genes; large
# batches in full mode → ~16k genes do not spawn thousands of scheduler jobs.
# WORKERS is capped per stage by its CPU allocation in conf/resources.config
# (batch_light=8, batch_heavy=16 cpus).
#
# HyPhy is sized SEPARATELY: FUBAR+BUSTED (first pass) and MEME (second pass) are
# far heavier per gene than the alignment/variability stages, and MEME is the
# pipeline bottleneck — so it gets its own batch/worker knobs.
if [[ "$TEST_MODE" == "true" ]]; then
    BATCH_SIZE=16;       BATCH_WORKERS=4        # dedup / filter / BMGE / variability
    PSEL_BATCH_SIZE=8;   PSEL_BATCH_WORKERS=4   # HyPhy (FUBAR/BUSTED/MEME)
else
    BATCH_SIZE=96;       BATCH_WORKERS=6
    PSEL_BATCH_SIZE=48;  PSEL_BATCH_WORKERS=8
fi

# ── Optional: restrict to a gene list ────────────────────────────────────────
GENE_LIST=""    # path to plain-text file (one gene name per line); leave empty for all

# ── Positive selection ────────────────────────────────────────────────────────
PSEL_SPECIES_TREE="${ROOT}/external_phylogenies/primates_233_subst.tree"          # subst/site branch lengths (NOT the fossil-cal timetree: its MY branch lengths overflow HyPhy MEME's rate*length domain) — required when RUN_POSITIVE_SELECTION=true
HYPHY_METHODS="FUBAR,BUSTED"  # first-pass methods. BUSTED is the gene-level gate that decides
                              # which genes get site-level MEME (2nd pass). FEL/SLAC also available.
FUBAR_THRESHOLD="0.9"         # posterior probability threshold for positive sites

# ── Variability / Pfam ────────────────────────────────────────────────────────
RUN_HMMSCAN=true              # true | false — run hmmscan against Pfam-A.hmm
PFAM_HMM="${ROOT}/pfam/Pfam-A.hmm"           # path to Pfam-A.hmm (used when RUN_HMMSCAN=true)
PFAM_DESC_TSV="${ROOT}/pfam/pfam_descriptions.tsv"  # pfam ID → description table

# ── Genomic coordinate mapping ────────────────────────────────────────────────
ALI_REF="/homes/users/mramon/scratch/4.Generate_alignments_from_codons/ali_ref/ppga"
GFF_FILE="${ALI_REF}/Homo_sapiens.sorted.gff"
GENE_EQUIV_TSV="${ALI_REF}/proteiID_gene_equivalences.txt"

# ── UCR analysis (active only when RUN_UCR=true and RUN_VARIABILITY=true) ─────
# Detection thresholds — leave empty to keep pipeline defaults (see variability.config).
UCR_ABS_THRESHOLD=""          # absolute method C_trident floor           (default: 0.9)
UCR_ABS_MIN_LEN=""            # absolute method min consecutive positions  (default: 5)
UCR_REL_ZSCORE=""             # relative method within-gene z-score cutoff (default: 1.5)
UCR_WINDOW_SIZE=""            # sliding-window size in aa                  (default: 10)
UCR_WINDOW_THRESHOLD=""       # sliding-window mean C_trident floor        (default: 0.85)
UCR_FLANK_SIZE=""             # flanking aa on each side of every window   (default: 5)
# Optional: site-level positive-selection data from a prior pipeline run.
# Provide these when RUN_POSITIVE_SELECTION=false but psel was run previously.
FUBAR_SITES_TSV=""            # fubar_sites.tsv from a prior positive-selection run
FEL_RESULTS_TSV=""            # fel_results.tsv  from a prior positive-selection run

# ── Phylogeny ─────────────────────────────────────────────────────────────────
GENE_TREE_METHOD="iqtree"     # iqtree | codeml
SPECIES_TREE_METHOD="astral"  # astral | backbone
IQTREE_MSET="LG"              # restrict ModelFinder to this model family; set "" for all

# ── Pre-build optional Nextflow args as a single string ──────────────────────
NF_OPTIONAL=""
[[ -n "$GENE_LIST"              ]] && NF_OPTIONAL+=" --gene_list $GENE_LIST"
[[ -n "$CHI2_MAX_ITER"          ]] && NF_OPTIONAL+=" --chi2_max_iter $CHI2_MAX_ITER"
[[ -n "$CHI2_AMBIGUITY_THRESH"  ]] && NF_OPTIONAL+=" --chi2_ambiguity_thresh $CHI2_AMBIGUITY_THRESH"
[[ -n "$CHI2_MIN_SEQ_FRAC"      ]] && NF_OPTIONAL+=" --chi2_min_seq_frac $CHI2_MIN_SEQ_FRAC"
[[ -n "$THRESHOLD_MIN_SEQ_FRAC" ]] && NF_OPTIONAL+=" --threshold_min_seq_frac $THRESHOLD_MIN_SEQ_FRAC"
[[ -n "$MIN_SEQS"               ]] && NF_OPTIONAL+=" --min_seqs $MIN_SEQS"
[[ -f "$FAMILY_ORDER_TSV"       ]] && NF_OPTIONAL+=" --family_order_tsv $FAMILY_ORDER_TSV"
[[ "$TEST_MODE" == "true"       ]] && NF_OPTIONAL+=" --test_mode true --test_n_genes $TEST_N_GENES"
# Batch concurrency — one global pair scales every batched stage (overrides conf/*.config)
NF_OPTIONAL+=" --run_dedup $RUN_DEDUP"
NF_OPTIONAL+=" --dedup_batch_size $BATCH_SIZE --dedup_batch_workers $BATCH_WORKERS"
NF_OPTIONAL+=" --filter_batch_size $BATCH_SIZE --filter_batch_workers $BATCH_WORKERS"
NF_OPTIONAL+=" --translation_batch_size $BATCH_SIZE --translation_batch_workers $BATCH_WORKERS"
NF_OPTIONAL+=" --variability_batch_size $BATCH_SIZE --variability_batch_workers $BATCH_WORKERS"
NF_OPTIONAL+=" --psel_batch_size $PSEL_BATCH_SIZE --psel_batch_workers $PSEL_BATCH_WORKERS"
[[ -n "$PSEL_SPECIES_TREE"      ]] && NF_OPTIONAL+=" --psel_species_tree $PSEL_SPECIES_TREE"
[[ -n "$IQTREE_MSET"            ]] && NF_OPTIONAL+=" --iqtree_mset $IQTREE_MSET"
[[ -n "$HYPHY_METHODS"          ]] && NF_OPTIONAL+=" --hyphy_methods $HYPHY_METHODS"
[[ -n "$FUBAR_THRESHOLD"        ]] && NF_OPTIONAL+=" --fubar_threshold $FUBAR_THRESHOLD"
[[ -n "$GENE_TREE_METHOD"       ]] && NF_OPTIONAL+=" --gene_tree_method $GENE_TREE_METHOD"
[[ -n "$SPECIES_TREE_METHOD"    ]] && NF_OPTIONAL+=" --species_tree_method $SPECIES_TREE_METHOD"
[[ "$RUN_HMMSCAN" == "true"     ]] && NF_OPTIONAL+=" --run_hmmscan true"
[[ "$RUN_HMMSCAN" == "true" && -n "$PFAM_HMM"      ]] && NF_OPTIONAL+=" --pfam_hmm $PFAM_HMM"
[[ -f "$PFAM_DESC_TSV"          ]] && NF_OPTIONAL+=" --pfam_descriptions_tsv $PFAM_DESC_TSV"
[[ -n "$GFF_FILE"               ]] && NF_OPTIONAL+=" --gff_file $GFF_FILE"
[[ -n "$GENE_EQUIV_TSV"         ]] && NF_OPTIONAL+=" --gene_equiv_tsv $GENE_EQUIV_TSV"
# DEDUP-specific metadata is only needed when deduplication is enabled.
[[ "$RUN_DEDUP" == "true"      ]] && NF_OPTIONAL+=" --assemblies_tsv $ASSEMBLIES_TSV"
[[ "$RUN_DEDUP" == "true"      ]] && NF_OPTIONAL+=" --species_dedup_map $SPECIES_DEDUP_MAP"
# UCR thresholds (only passed when non-empty, so defaults in variability.config apply otherwise)
[[ -n "$UCR_ABS_THRESHOLD"      ]] && NF_OPTIONAL+=" --ucr_abs_threshold $UCR_ABS_THRESHOLD"
[[ -n "$UCR_ABS_MIN_LEN"        ]] && NF_OPTIONAL+=" --ucr_abs_min_len $UCR_ABS_MIN_LEN"
[[ -n "$UCR_REL_ZSCORE"         ]] && NF_OPTIONAL+=" --ucr_rel_zscore $UCR_REL_ZSCORE"
[[ -n "$UCR_WINDOW_SIZE"        ]] && NF_OPTIONAL+=" --ucr_window_size $UCR_WINDOW_SIZE"
[[ -n "$UCR_WINDOW_THRESHOLD"   ]] && NF_OPTIONAL+=" --ucr_window_threshold $UCR_WINDOW_THRESHOLD"
[[ -n "$UCR_FLANK_SIZE"         ]] && NF_OPTIONAL+=" --ucr_flank_size $UCR_FLANK_SIZE"
[[ -n "$FUBAR_SITES_TSV"        ]] && NF_OPTIONAL+=" --fubar_sites_tsv $FUBAR_SITES_TSV"
[[ -n "$FEL_RESULTS_TSV"        ]] && NF_OPTIONAL+=" --fel_results_tsv $FEL_RESULTS_TSV"
# Pre-computed input dirs
[[ -n "$PROT_DIR"               ]] && NF_OPTIONAL+=" --prot_dir $PROT_DIR"
[[ -n "$TRIM_DIR"               ]] && NF_OPTIONAL+=" --trim_dir $TRIM_DIR"
[[ -n "$MAP_DIR"                ]] && NF_OPTIONAL+=" --map_dir $MAP_DIR"
[[ -n "$DEDUP_DIR"              ]] && NF_OPTIONAL+=" --dedup_dir $DEDUP_DIR"

# ============================================================

submit_job() {
    sbatch --parsable <<EOF
#!/bin/bash
#SBATCH --job-name=generate-alis
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
    -name "generate-alis-ppga-$(date +%Y%m%d-%H%M%S)" \\
    -with-tower \\
    --cds_dir                "${CDS_DIR}" \\
    --out_dir                "${OUT_DIR}" \\
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
    -work-dir "${WORK_DIR}"

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
mv slurm-*.out SBATCH_runner
mv slurm-*.err SBATCH_runner
