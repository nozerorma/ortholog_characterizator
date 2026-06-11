#!/usr/bin/env nextflow

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * generate-alis: Mammalian CDS alignment generation, variability, and positive selection
 * Author: Miguel Ramon (miguel.ramon@upf.edu)
 *
 * File: main.nf
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

nextflow.enable.dsl = 2

include { QUALITY            } from './workflows/quality'
include { TRANSLATION        } from './workflows/translation'
include { PHYLOGENY          } from './workflows/phylogeny'
include { VARIABILITY        } from './workflows/variability'
include { POSITIVE_SELECTION } from './workflows/positive_selection'
include { UCR_ANALYSIS       } from './subworkflows/variability/ucr'
include { TRANSLATION_REPORT } from './subworkflows/translation/report'
include { VARIABILITY_REPORT } from './subworkflows/variability/report'
include { PSEL_REPORT        } from './subworkflows/positive_selection/report'

// ── Help message ─────────────────────────────────────────────────────────────

def helpMessage() {
    log.info """
    ╔═══════════════════════════════════════════════════════════╗
    ║              generate-alis Nextflow pipeline              ║
    ╚═══════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf [options] -profile local|slurm

    Test mode:
        --test_mode              true/false  (default: false)
        --test_n_genes           int         (default: 100)
        --test_seed              int         (default: 42; set null for random)

    Required:
        --cds_dir             Directory of codon CDS alignments (.fa)
        --out_dir             Output directory (default: results)
        --assemblies_tsv      assemblies_and_species_mammals.tsv
        --taxid_tsv           taxid_species_family_mammals.tsv
        --species_dedup_map   species_dedup_map.tsv

    Workflow toggles (all default false except quality + translation):
        --run_quality            true/false  (default: true)
        --run_translation        true/false  (default: true)
        --run_phylogeny          true/false  (default: false)
        --run_variability        true/false  (default: false)
        --run_positive_selection true/false  (default: false)

    Quality options:
        --filter_mode            chi2|threshold|none  (default: chi2)
        --chi2_max_iter          int   (default: 3)
        --chi2_ambiguity_thresh  float (default: 0.8)
        --chi2_min_seq_frac      float (default: 0.1)
        --threshold_min_seq_frac float (default: 0.2; drop seqs with ambiguity > 1-value)
        --min_seqs               int   (default: 4)

    Translation options:
        --bmge_gap_threshold     float (default: 0.0; only all-gap columns removed)
        --bmge_entropy_threshold float (default: 1.0 = no trimming; set 0.5 for tree-only runs)
        --gff_file               path  (optional)
        --gene_equiv_tsv         path  (optional)

    Phylogeny options:
        --gene_tree_method       iqtree|codeml  (default: iqtree)
        --species_tree_method    astral|backbone (default: astral)
        --backbone_tree          path  (backbone mode only)
        --outgroup_taxa          str   comma-separated assembly names

    Variability options:
        --valdar_alpha/beta/gamma  float (default: 1.0)
        --run_hmmscan              true/false (default: false)
        --pfam_hmm                 path

    UCR analysis options (requires --run_variability true):
        --run_ucr                  true/false (default: false)
        --ucr_abs_threshold        float (default: 0.9)   absolute C_trident floor
        --ucr_abs_min_len          int   (default: 5)     min consecutive positions
        --ucr_rel_zscore           float (default: 1.5)   within-gene z-score cutoff
        --ucr_window_size          int   (default: 10)    sliding window size (aa)
        --ucr_window_threshold     float (default: 0.85)  sliding window C_trident floor
        --ucr_flank_size           int   (default: 5)     flanking aa on each side
        --fubar_sites_tsv          path  (optional: fubar_sites.tsv from prior psel run)
        --fel_results_tsv          path  (optional: fel_results.tsv from prior psel run)

    Positive selection options:
        --hyphy_methods          FUBAR,FEL[,SLAC]  (default: FUBAR,FEL)
        --psel_species_tree      path  (required)
        --fubar_threshold        float (default: 0.9)
        --sel_variability_tsv    path  (optional: sel_variability.tsv from a prior positive-selection run,
                                        used to enrich the variability report when psel is not run concurrently)
        --variability_tsv        path  (optional: gene_variability.tsv from a prior variability run,
                                        used to enrich the psel report when variability is not run concurrently)

    Profiles:
        -profile local     Local execution
        -profile slurm     SLURM cluster
        -profile singularity  Use Singularity container
        -profile apptainer    Use Apptainer container
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ── Input channel ─────────────────────────────────────────────────────────────

def make_cds_channel() {
    assert params.cds_dir : "ERROR: --cds_dir is required"
    def ch = Channel.fromPath("${params.cds_dir}/*.fa", checkIfExists: true)
    if (params.gene_list) {
        def genes = file(params.gene_list).readLines().collect { it.trim() }.toSet()
        ch = ch.filter { f -> genes.contains(f.baseName) }
    }
    if (params.test_mode) {
        def n    = params.test_n_genes as int
        def seed = params.test_seed != null ? params.test_seed as long : new Random().nextLong()
        log.warn "TEST MODE: sampling ${n} random alignments (seed=${seed})"
        ch = ch.randomSample(n, seed)
    }
    return ch.map { f -> tuple(f.baseName, f) }
}

def make_fallback_channel(String dir_param, String param_name, String glob = '*.fa') {
    assert params[dir_param] : "ERROR: upstream workflow was skipped — set --${dir_param} to point at existing ${param_name}"
    def ch = Channel.fromPath("${params[dir_param]}/${glob}", checkIfExists: true)
    if (params.gene_list) {
        def genes = file(params.gene_list).readLines().collect { it.trim() }.toSet()
        ch = ch.filter { f -> genes.contains(f.baseName.replace('.map', '')) }
    }
    if (params.test_mode) {
        def n    = params.test_n_genes as int
        def seed = params.test_seed != null ? params.test_seed as long : new Random().nextLong()
        log.warn "TEST MODE (fallback ${param_name}): sampling ${n} random alignments (seed=${seed})"
        ch = ch.randomSample(n, seed)
    }
    return ch.map { f -> tuple(f.baseName.replace('.map', ''), f) }
}

// ── Main workflow ─────────────────────────────────────────────────────────────

workflow {

    ch_dedup = Channel.empty()
    ch_trim  = Channel.empty()
    ch_prot  = Channel.empty()
    ch_map   = Channel.empty()

    // ── Quality ───────────────────────────────────────────────────────────
    ch_filtered = Channel.empty()
    if (params.run_quality) {
        QUALITY( make_cds_channel() )
        ch_dedup    = QUALITY.out.dedup
        ch_filtered = QUALITY.out.filtered
    }

    // ── Translation ───────────────────────────────────────────────────────
    if (params.run_translation) {
        def trans_in = params.run_quality
            ? ch_filtered
            : make_fallback_channel('dedup_dir', 'DEDUP')
        TRANSLATION( trans_in )
        ch_prot = TRANSLATION.out.prot
        ch_trim = TRANSLATION.out.trim
        ch_map  = TRANSLATION.out.map_tsv
    }

    // ── Phylogeny ─────────────────────────────────────────────────────────
    if (params.run_phylogeny) {
        def phylo_in = params.run_translation
            ? ch_prot
            : make_fallback_channel('prot_dir', 'PROT/bmge')
        PHYLOGENY( phylo_in )
    }

    // ── Variability ───────────────────────────────────────────────────────
    if (params.run_variability) {
        def var_in = params.run_translation
            ? ch_prot
            : make_fallback_channel('prot_dir', 'PROT/bmge')

        VARIABILITY( var_in )
    }

    // ── Positive selection ────────────────────────────────────────────────
    if (params.run_positive_selection) {
        def psel_in = params.run_translation
            ? ch_trim
            : make_fallback_channel('trim_dir', 'TRIM/bmge')

        def map_in = params.run_translation
            ? ch_map
            : (params.map_dir ? make_fallback_channel('map_dir', 'MAP', '*.map.tsv') : Channel.empty())

        def gene_var_ch = params.run_variability
            ? VARIABILITY.out.gene_variability_tsv
            : (params.variability_tsv
               ? Channel.value(file(params.variability_tsv))
               : Channel.value(file('NO_VARIABILITY_TSV')))

        POSITIVE_SELECTION( psel_in, map_in, gene_var_ch )
    }

    // ── UCR analysis ──────────────────────────────────────────────────────
    if (params.run_variability && params.run_ucr) {
        assert params.taxid_tsv : "ERROR: --taxid_tsv is required for UCR analysis"

        def ucr_prot_ch = params.run_translation
            ? TRANSLATION.out.prot
            : make_fallback_channel('prot_dir', 'PROT/bmge')

        def ucr_map_ch = params.run_translation
            ? TRANSLATION.out.map_tsv
            : (params.map_dir
               ? make_fallback_channel('map_dir', 'MAP', '*.map.tsv')
               : Channel.empty())

        // Positive-selection data: from concurrent run or prior-run params
        def ucr_fubar_ch = params.run_positive_selection
            ? POSITIVE_SELECTION.out.fubar_sites
            : (params.fubar_sites_tsv
               ? Channel.value(file(params.fubar_sites_tsv))
               : Channel.value(file('NO_FUBAR')))

        def ucr_fel_ch = params.run_positive_selection
            ? POSITIVE_SELECTION.out.fel_tsv
            : (params.fel_results_tsv
               ? Channel.value(file(params.fel_results_tsv))
               : Channel.value(file('NO_FEL')))

        UCR_ANALYSIS(
            VARIABILITY.out.entropy_tsvs,
            ucr_prot_ch,
            ucr_map_ch,
            VARIABILITY.out.domain_tsv,
            ucr_fubar_ch,
            ucr_fel_ch
        )
    }

    // ── Reports (all run last — require outputs from multiple workflows) ──────
    if (params.run_translation) {
        TRANSLATION_REPORT(
            TRANSLATION.out.prot.collect { it[1] },
            TRANSLATION.out.trim.collect { it[1] }.ifEmpty([]),
            TRANSLATION.out.map_tsv.collect { it[1] }.ifEmpty([])
        )
    }

    if (params.run_variability) {
        def sel_var_report_ch = params.run_positive_selection
            ? POSITIVE_SELECTION.out.sel_variability_tsv
            : (params.sel_variability_tsv
               ? Channel.value(file(params.sel_variability_tsv))
               : Channel.empty())

        def ucr_windows_report_ch = (params.run_variability && params.run_ucr)
            ? UCR_ANALYSIS.out.ucr_windows_tsv
            : Channel.empty()

        def ucr_positions_report_ch = (params.run_variability && params.run_ucr)
            ? UCR_ANALYSIS.out.ucr_positions_tsv
            : Channel.empty()

        VARIABILITY_REPORT(
            VARIABILITY.out.gene_variability_tsv,
            VARIABILITY.out.clade_variability_tsv,
            VARIABILITY.out.clade_outliers_tsv,
            VARIABILITY.out.entropy_tsvs.collect { it[1] }
                .mix( VARIABILITY.out.clade_tsvs.collect { it[1] } )
                .collect(),
            VARIABILITY.out.prot_files,
            VARIABILITY.out.domain_tsv,
            sel_var_report_ch,
            ucr_windows_report_ch,
            ucr_positions_report_ch
        )
    }

    if (params.run_positive_selection) {
        def gene_var_report_ch = params.run_variability
            ? VARIABILITY.out.gene_variability_tsv
            : (params.variability_tsv
               ? Channel.value(file(params.variability_tsv))
               : Channel.empty())

        PSEL_REPORT(
            POSITIVE_SELECTION.out.fubar_tsv,
            POSITIVE_SELECTION.out.fubar_sites,
            POSITIVE_SELECTION.out.fel_tsv,
            POSITIVE_SELECTION.out.meme_tsv,
            gene_var_report_ch.ifEmpty(file('NO_GENE_VAR'))
        )
    }
}
