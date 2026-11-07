# Ortholog Characterizator

Nextflow pipeline for ortholog CDS characterization, variability analysis, and positive selection detection. Tested in primate ortholog codon CDS set and in [TOGA 2 Multi Codon Alignments](https://genome.senckenberg.de/download/TOGA2/MultiCodonAlignments/)

Takes a directory of codon-level multiple sequence alignments (one `.fa` per gene) and runs them through up to six independent stages that can be toggled on or off.

---

## Pipeline overview

```
CDS alignments (.fa)
        │
        ▼
┌───────────────┐
│    QUALITY    │  deduplicate assemblies + filter low-quality sequences
└───────┬───────┘
        │
        ▼
┌───────────────┐
│  TRANSLATION  │  BMGE codon-aware trimming → amino acid alignments
└───────┬───────┘
        │
        ├──────────────────────────┬──────────────────────────┐
        ▼                          ▼                          ▼
┌───────────────┐        ┌──────────────────┐      ┌──────────────────────┐
│  VARIABILITY  │        │    PHYLOGENY      │      │  POSITIVE SELECTION  │
│ Valdar entropy│        │ IQ-TREE / CodeML  │      │   HyPhy FUBAR/FEL    │
│ + Pfam domains│        │ gene & species    │      │   per gene           │
│      │        │        │ trees             │      └──────────────────────┘
│      ▼        │        └──────────────────┘
│     UCR       │
│  ultra-cons.  │
│  regions      │
└───────────────┘
```

Each stage emits an HTML report.

---

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 21.04
- One of:
  - **Local**: all tools installed and available on `$PATH`
  - **Container** (recommended): Singularity or Apptainer image `apptainer/generate-alis-latest.img`
  - **SLURM cluster**: use `-profile slurm` — jobs are submitted automatically

Tools used internally: BMGE 1.12, IQ-TREE, PAML/CodeML, HyPhy, ASTRAL, hmmscan, Python 3, R/Quarto.

---

## Quick start (local)

1. Edit the `CONFIGURATION` block at the top of [run.sh](run.sh):

```bash
CDS_DIR="/path/to/cds_alignments"   # directory of codon MSA .fa files
```

2. Run:

```bash
bash nextflow_rewrite/run.sh
```

Results land in `results/YYYYMMDD_HHMMSS/`.

---

## SLURM cluster

Three launcher scripts are provided — pick the one that matches your run type:

| Script | Purpose |
|---|---|
| [submit_bmge.sh](submit_bmge.sh) | Variability + positive selection (default BMGE passthrough) |
| [submit_phylo.sh](submit_phylo.sh) | Quality + translation + gene/species trees (BMGE trimming enabled) |

Edit the `CONFIGURATION` block in the chosen script, then:

```bash
bash nextflow_rewrite/submit_bmge.sh   # or submit_phylo.sh
```

This submits a single lightweight orchestrator job. All heavy per-gene jobs are dispatched automatically by Nextflow. Logs go to `nextflow_rewrite/Slurm/`.

> **Note:** The launcher scripts pass `-with-tower`, which requires a `TOWER_ACCESS_TOKEN` environment variable set in your session. Remove the flag if you are not using Nextflow Tower.

---

## Input files

| Parameter | Description |
|---|---|
| `--cds_dir` | Directory of codon CDS alignments (one `.fa` per gene) — **required** |
| `--out_dir` | Output directory (default: `results`) |
| `--assemblies_tsv` | Assembly → species mapping table |
| `--taxid_tsv` | NCBI taxid → species/family table |
| `--species_dedup_map` | Species deduplication map |
| `--family_order_tsv` | Family → order table (optional) |
| `--gene_list` | Plain-text file (one gene name per line) to restrict the run to a subset of genes (optional) |

Reference TSVs are expected under `taxa_info/` relative to the repo root.

---

## Workflow stages

All stages default to **off** except `quality` and `translation`.

| Flag | Default | Description |
|---|---|---|
| `--run_quality` | `true` | Deduplication + sequence filtering |
| `--run_translation` | `true` | BMGE trimming + protein projection |
| `--run_variability` | `false` | Per-position Valdar entropy + optional Pfam mapping |
| `--run_ucr` | `false` | Ultra-Conserved Region detection (requires `run_variability=true`) |
| `--run_phylogeny` | `false` | Per-gene trees + species tree |
| `--run_positive_selection` | `false` | HyPhy FUBAR/FEL/MEME per gene |

Stages can be chained or run independently. If an upstream stage is skipped, point to its pre-computed output with `--dedup_dir`, `--prot_dir`, `--trim_dir`, or `--map_dir`.

---

## Key parameters

### Quality

| Parameter | Default | Description |
|---|---|---|
| `--filter_mode` | `chi2` | `chi2` iterative outlier removal, `threshold` fixed cutoff, `none` skip |
| `--chi2_max_iter` | `3` | Max chi-square iterations |
| `--chi2_ambiguity_thresh` | `0.8` | Stop iterating below this fraction of ambiguous codons |
| `--chi2_min_seq_frac` | `0.1` | Min fraction of sequences that must survive chi2 filtering |
| `--threshold_min_seq_frac` | `0.2` | Min non-gap/N fraction per sequence (threshold mode only) |
| `--min_seqs` | `4` | Minimum sequences per alignment after filtering |

### Translation (BMGE)

| Parameter | Default | Description |
|---|---|---|
| `--bmge_gap_threshold` | `0.0` | Remove columns gapped in > this fraction of sequences (0.0 = only all-gap columns) |
| `--bmge_entropy_threshold` | `1.0` | Max per-column entropy (1.0 = no trimming; use 0.5 for tree-only runs) |
| `--gff_file` | — | GFF annotation file (optional) |
| `--gene_equiv_tsv` | — | Gene-name equivalence table (optional) |

### Variability

| Parameter | Default | Description |
|---|---|---|
| `--valdar_alpha` | `1.0` | Valdar entropy α weight |
| `--valdar_beta` | `1.0` | Valdar entropy β weight |
| `--valdar_gamma` | `1.0` | Valdar entropy γ weight |
| `--run_hmmscan` | `false` | Run hmmscan to annotate Pfam domains |
| `--pfam_hmm` | — | Path to `Pfam-A.hmm` (required when `run_hmmscan=true`) |
| `--hmmscan_domtblout` | — | Precomputed `domtblout`; it must come from `hmmscan --cut_ga` on the full-length human queries |
| `--pfam_clans_tsv` | — | `Pfam-A.clans.tsv`, adds clan columns to the domain table |

Pfam is scanned on the full-length human protein, translated from the unfiltered CDS alignment (`--cds_dir`), with the Pfam gathering thresholds (`hmmscan --cut_ga`, the curated criterion of Pfam). BMGE removes columns that are interior to the protein, and a query made of the trimmed human sequence joins stretches that are not contiguous in it. Each hit is projected onto the trimmed alignment through the MAP tables, so the Pfam step needs `--cds_dir` and the MAP tables (from the translation stage or `--map_dir`) besides the protein alignments. `reports/variability/domain_variability.tsv` keeps the hmmscan endpoints (`ali_start`, `ali_end`, on the full-length protein), the exact trimmed columns (`alignment_positions`) and the residues of the domain lost to trimming (`n_removed`).

### UCR analysis

Requires `--run_variability true`. Detects ultra-conserved regions using two complementary methods (absolute C\_trident floor and within-gene relative z-score), plus a sliding-window approach.

| Parameter | Default | Description |
|---|---|---|
| `--ucr_abs_threshold` | `0.9` | Absolute C\_trident floor for the absolute method |
| `--ucr_abs_min_len` | `5` | Minimum number of consecutive conserved positions |
| `--ucr_rel_zscore` | `1.5` | Within-gene z-score cutoff for the relative method |
| `--ucr_window_size` | `10` | Sliding window size (amino acids) |
| `--ucr_window_threshold` | `0.85` | Mean C\_trident floor within a sliding window |
| `--ucr_flank_size` | `5` | Flanking amino acids reported on each side of a UCR |
| `--fubar_sites_tsv` | — | `fubar_sites.tsv` from a prior positive-selection run (enriches UCR output without re-running psel) |
| `--fel_results_tsv` | — | `fel_results.tsv` from a prior positive-selection run |

### Positive selection

| Parameter | Default | Description |
|---|---|---|
| `--psel_species_tree` | — | Species tree in Newick format — **required** |
| `--hyphy_methods` | `FUBAR,FEL` | Comma-separated list: `FUBAR`, `FEL`, `MEME`, `SLAC` |
| `--fubar_threshold` | `0.9` | Posterior probability threshold for positive sites |
| `--sel_variability_tsv` | — | `sel_variability.tsv` from a prior psel run (enriches the variability report without re-running psel) |
| `--variability_tsv` | — | `gene_variability.tsv` from a prior variability run (enriches the psel report without re-running variability) |

### Phylogeny

| Parameter | Default | Description |
|---|---|---|
| `--gene_tree_method` | `iqtree` | `iqtree` or `codeml` |
| `--species_tree_method` | `astral` | `astral` or `backbone` |
| `--backbone_tree` | — | Newick backbone tree (required for `backbone` mode) |
| `--outgroup_taxa` | — | Comma-separated assembly names to use as outgroup |

---

## Pre-computed inputs (skipping stages)

When an upstream stage has already been run, point downstream stages to the existing outputs to avoid recomputing:

| Parameter | Skips | Expected contents |
|---|---|---|
| `--dedup_dir` | quality | Deduplicated `.fa` files |
| `--prot_dir` | quality + translation | Protein `.fa` files (BMGE output) |
| `--trim_dir` | quality + translation | BMGE-trimmed CDS `.fa` files |
| `--map_dir` | quality + translation | BMGE position-map `.map.tsv` files |

---

## Test mode

To validate the pipeline on a small random subset before a full run:

```bash
nextflow run main.nf --test_mode true --test_n_genes 100 --test_seed 42 [other params]
```

---

## Output structure

```
results/
├── pipeline_info/          # Nextflow execution report, timeline, trace, DAG
├── quality/
│   ├── dedup/              # Deduplicated alignments
│   ├── filtered/           # Quality-filtered alignments
│   └── reports/            # HTML deduplication + filtering report
├── translation/
│   ├── prot/               # Amino acid alignments
│   ├── trim/               # BMGE-trimmed CDS alignments
│   └── reports/            # HTML translation quality report
├── variability/
│   ├── entropy/            # Per-gene and per-clade entropy TSVs
│   ├── ucr/                # UCR windows and positions TSVs (when run_ucr=true)
│   └── reports/            # HTML variability report
├── phylogeny/
│   ├── gene_trees/         # Per-gene trees
│   └── species_tree/       # Species-level tree
└── positive_selection/
    ├── hyphy/              # Per-gene HyPhy JSON results
    └── reports/            # HTML positive selection report
```

---

## Profiles

| Profile | Use |
|---|---|
| `local` | Workstation / laptop |
| `slurm` | SLURM cluster (Marvin @ UPF-CSIC by default) |
| `singularity` | Singularity container |
| `apptainer` | Apptainer container |

Combine profiles with commas: `-profile slurm,singularity`

---

## Re-running / resuming

Nextflow caches completed tasks. Add `-resume` (already included in the launcher scripts) to skip steps that finished successfully in a previous run.

---

## Re-running Pfam only

```bash
nextflow run main.nf -entry PFAM_ONLY --cds_dir <CDS alignments> --entropy_dir <PROT_VAR_RAW> --map_dir <MAP> \
    --pfam_hmm <Pfam-A.hmm> --out_dir <new results dir> -profile slurm
```

This rewrites `HMMSCAN/hmmscan.domtblout` and `reports/variability/domain_variability.tsv` from the entropy tables and MAP tables of a finished run. Entropy, UCR detection, positive selection and the reports are not run again. The `pfam_domain` column of `ucr_positions.tsv` is taken from the Pfam table when UCR is aggregated, so it is updated only when that step is run again.

---

## Re-rendering reports

To regenerate the HTML reports from an existing results directory without re-running the pipeline:

```bash
bash nextflow_rewrite/rerender_reports.sh \
    --results_dir /path/to/results/20260517_235350 \
    [--psel_tree   /path/to/species.nwk] \
    [--sel_var_tsv /path/to/sel_variability.tsv] \
    [--out_dir     /path/to/output_reports]
```

The script auto-detects translation mode, filter mode, and reference TSVs from the repository layout. Reports are written to `results_dir/reports/` by default.
