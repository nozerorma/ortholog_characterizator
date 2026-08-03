#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ── Processes ─────────────────────────────────────────────────────────────────

process DETECT_UCR {
    tag   "${gene}"
    label 'process_single'

    input:
    tuple val(gene), path(entropy_tsv)

    output:
    tuple val(gene), path("${gene}.ucr.tsv"), optional: true, emit: ucr_tsv

    script:
    """
    python3 ${projectDir}/subworkflows/variability/local/detect_ucr.py \\
        --entropy_tsv       "${entropy_tsv}" \\
        --gene              "${gene}"        \\
        --out_tsv           "${gene}.ucr.tsv" \\
        --abs_threshold     ${params.ucr_abs_threshold}     \\
        --abs_min_len       ${params.ucr_abs_min_len}       \\
        --abs_max_gap       ${params.ucr_abs_max_gap}       \\
        --rel_zscore        ${params.ucr_rel_zscore}        \\
        --rel_min_len       ${params.ucr_rel_min_len}       \\
        --window_size       ${params.ucr_window_size}       \\
        --window_threshold  ${params.ucr_window_threshold}  \\
        --window_max_pos_gap ${params.ucr_window_max_pos_gap} \\
        --flank_size        ${params.ucr_flank_size}
    """
}

process AGGREGATE_UCR {
    label 'process_single'
    publishDir "${params.out_dir}/reports/variability", mode: 'copy'

    input:
    path ucr_files      // *.ucr.tsv — one per gene (may include NO_UCR sentinel)
    path entropy_files  // *.entropy.tsv — all genes
    path prot_files     // *.fa raw protein MSA — all genes
    path map_files      // *.map.tsv for position mapping (may include NO_MAPS sentinel)
    path taxid_tsv
    path domain_tsv     // domain_variability.tsv or NO_DOMAIN_TSV sentinel
    path fubar_sites    // fubar_sites.tsv or NO_FUBAR sentinel

    output:
    path "ucr_windows.tsv",             emit: ucr_windows
    path "ucr_positions.tsv",           emit: ucr_positions
    path "ucr_clade_variability.tsv",   optional: true, emit: ucr_clade_variability
    path "ucr_selection_summary.tsv",   optional: true, emit: ucr_selection_summary

    script:
    def fam_arg   = params.family_order_tsv
                    ? "--family_order_tsv '${params.family_order_tsv}'"
                    : ''
    def alpha_arg = "--alpha ${params.valdar_alpha} --beta ${params.valdar_beta} --gamma ${params.valdar_gamma}"
    def dom_arg   = (domain_tsv.name  != 'NO_DOMAIN_TSV') ? "--domain_tsv '${domain_tsv}'"   : ''
    def fubar_arg = (fubar_sites.name != 'NO_FUBAR')       ? "--fubar_sites '${fubar_sites}'" : ''
    """
    mkdir -p UCR_RAW PROT_VAR_RAW PROT_RAW MAP

    # Stage files into expected subdirectories
    for f in *.ucr.tsv;     do [[ -f "\$f" ]] && mv "\$f" UCR_RAW/      || true; done
    for f in *.entropy.tsv; do [[ -f "\$f" ]] && mv "\$f" PROT_VAR_RAW/ || true; done
    for f in *.fa;          do [[ -f "\$f" ]] && mv "\$f" PROT_RAW/      || true; done
    for f in *.map.tsv;     do [[ -f "\$f" ]] && mv "\$f" MAP/            || true; done

    python3 ${projectDir}/subworkflows/variability/local/aggregate_ucr.py \\
        --ucr_dir     UCR_RAW      \\
        --entropy_dir PROT_VAR_RAW \\
        --prot_dir    PROT_RAW     \\
        --map_dir     MAP          \\
        --taxid_tsv   "${taxid_tsv}" \\
        --out_dir     .            \\
        --flank_size  ${params.ucr_flank_size} \\
        ${alpha_arg}  \\
        ${fam_arg}    \\
        ${dom_arg}    \\
        ${fubar_arg}
    """
}

// ── Subworkflow ───────────────────────────────────────────────────────────────

workflow UCR_ANALYSIS {

    take:
    entropy_ch    // tuple(gene_id, entropy_tsv)  — from ENTROPY.out.entropy_tsv
    prot_ch       // tuple(gene_id, prot_fasta)   — raw translated protein MSA
    map_ch        // tuple(gene_id, map_tsv)       — BMGE position maps (may be empty)
    domain_tsv    // single file or sentinel file('NO_DOMAIN_TSV')
    fubar_sites   // single file or sentinel file('NO_FUBAR')

    main:

    DETECT_UCR( entropy_ch )

    // Collect per-gene UCR files; fall back to a sentinel when none detected
    ucr_files_ch = DETECT_UCR.out.ucr_tsv
        .collect { it[1] }
        .ifEmpty( [file('NO_UCR')] )

    AGGREGATE_UCR(
        ucr_files_ch,
        entropy_ch.collect { it[1] },
        prot_ch.collect    { it[1] },
        map_ch.collect     { it[1] }.ifEmpty( [file('NO_MAPS')] ),
        file(params.taxid_tsv),
        domain_tsv,
        fubar_sites
    )

    emit:
    ucr_windows_tsv         = AGGREGATE_UCR.out.ucr_windows
    ucr_positions_tsv       = AGGREGATE_UCR.out.ucr_positions
    ucr_clade_variability   = AGGREGATE_UCR.out.ucr_clade_variability
    ucr_selection_summary   = AGGREGATE_UCR.out.ucr_selection_summary
}
