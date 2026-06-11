#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ── Processes ─────────────────────────────────────────────────────────────────

process GENERATE_METADATA {
    label 'process_single'
    publishDir "${params.out_dir}/reports/translation/metadata", mode: 'copy'

    input:
    path prot_files   // all prot FA files staged in work dir

    output:
    path "metadata.*", emit: metadata_files

    script:
    // generate_metadata.py takes positional args: portion_dirs... output_prefix
    // We pass '.' as the portion directory (all prot FAs are staged here)
    """
    python3 ${projectDir}/subworkflows/translation/local/generate_metadata.py \\
        . metadata
    """
}

process QUALITY_REPORT_DATA {
    label 'process_single'
    publishDir "${params.out_dir}/reports/translation/quality", mode: 'copy'

    input:
    path prot_files   // all prot FA files staged in work dir

    output:
    path "quality_data_quality_per_gene.csv", emit: quality_tsv
    path "quality_data_quality_summary.csv",  emit: quality_summary, optional: true

    script:
    // quality_report.py: --align_dir, --out_dir, --output, [--gene_list]
    """
    python3 ${projectDir}/subworkflows/translation/local/quality_report.py \\
        --align_dir . \\
        --out_dir   . \\
        --output    quality_data
    """
}

process RENDER_TRANSLATION_REPORT {
    label 'process_single'
    publishDir "${params.out_dir}/reports", mode: 'copy'

    input:
    path prot_files
    path metadata_files
    path quality_tsv
    path trim_files
    path map_files

    output:
    path "post_translation_report.html", emit: report_html

    script:
    def filter_mode  = params.filter_mode        ?: 'chi2'
    def assemblies   = params.assemblies_tsv     ?: ''
    def taxid        = params.taxid_tsv          ?: ''
    def fam_order    = params.family_order_tsv   ?: ''
    def sp_tree      = params.psel_species_tree  ?: ''
    def dedup_dir     = "${params.out_dir}/DEDUP"
    def prefilter_dir = (filter_mode != 'none')
                        ? "${params.out_dir}/PREFILTER/${filter_mode}"
                        : ''
    def trim_dir      = "${params.out_dir}/TRIM/bmge"
    // Locate the metadata TSV (generate_metadata.py writes metadata.<something>.tsv)
    """
    cp ${projectDir}/subworkflows/translation/local/post_translation_report.qmd .
    META_TSV=\$(ls metadata.*.tsv 2>/dev/null | head -1 || echo "")
    PROT_DIR=\$(pwd) \\
    METADATA_TSV="\${META_TSV}" \\
    QUALITY_TSV="${quality_tsv}" \\
    ASSEMBLIES_TSV="${assemblies}" \\
    TAXID_TSV="${taxid}" \\
    FAMILY_ORDER_TSV="${fam_order}" \\
    SPECIES_TREE="${sp_tree}" \\
    TRANSLATION_MODE="bmge" \\
    DEDUP_DIR="${dedup_dir}" \\
    PREFILTER_DIR="${prefilter_dir}" \\
    FILTER_MODE="${filter_mode}" \\
    TRIM_DIR="${trim_dir}" \\
    quarto render post_translation_report.qmd \\
        --output post_translation_report.html \\
        --output-dir .
    """
}

// ── Subworkflow ───────────────────────────────────────────────────────────────

workflow TRANSLATION_REPORT {

    take:
    prot_files_ch   // collected list of prot FA files
    trim_files_ch   // collected list of trim FA files (may be empty)
    map_files_ch    // collected list of map TSV files (may be empty)

    main:

    GENERATE_METADATA( prot_files_ch )
    QUALITY_REPORT_DATA( prot_files_ch )

    RENDER_TRANSLATION_REPORT(
        prot_files_ch,
        GENERATE_METADATA.out.metadata_files.collect(),
        QUALITY_REPORT_DATA.out.quality_tsv,
        trim_files_ch,
        map_files_ch
    )

    emit:
    report_html     = RENDER_TRANSLATION_REPORT.out.report_html
    quality_tsv     = QUALITY_REPORT_DATA.out.quality_tsv
    quality_summary = QUALITY_REPORT_DATA.out.quality_summary
    metadata_files  = GENERATE_METADATA.out.metadata_files
}
