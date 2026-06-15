#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process RENDER_PSEL_REPORT {
    label 'process_single'
    publishDir "${params.out_dir}/reports", mode: 'copy'

    input:
    path fubar_tsv
    path fubar_sites_tsv
    path busted_tsv
    path meme_tsv
    path gene_variability_tsv

    output:
    path "positive_selection_report.html", emit: report_html

    script:
    def taxid     = params.taxid_tsv         ? params.taxid_tsv         : ''
    def fam_order = params.family_order_tsv  ? params.family_order_tsv  : ''
    def sp_tree   = params.psel_species_tree ? params.psel_species_tree : ''
    def methods   = params.hyphy_methods
    def threshold = params.fubar_threshold
    """
    cp ${projectDir}/subworkflows/positive_selection/local/positive_selection_report.qmd .
    FUBAR_TSV="${fubar_tsv}" \\
    FUBAR_SITES_TSV="${fubar_sites_tsv}" \\
    BUSTED_TSV="${busted_tsv}" \\
    MEME_TSV="${meme_tsv}" \\
    GENE_VARIABILITY_TSV="${gene_variability_tsv}" \\
    TAXID_TSV="${taxid}" \\
    FAMILY_ORDER_TSV="${fam_order}" \\
    SPECIES_TREE="${sp_tree}" \\
    HYPHY_METHODS="${methods}" \\
    FUBAR_THRESHOLD="${threshold}" \\
    FUBAR_FDR="${params.fubar_fdr}" \\
    BUSTED_FDR="${params.busted_fdr}" \\
    MEME_FDR="${params.meme_fdr}" \\
    quarto render positive_selection_report.qmd \\
        --output positive_selection_report.html \\
        --output-dir .
    """
}

workflow PSEL_REPORT {

    take:
    fubar_tsv
    fubar_sites_tsv
    busted_tsv
    meme_tsv
    gene_variability_tsv

    main:

    RENDER_PSEL_REPORT(
        fubar_tsv,
        fubar_sites_tsv,
        busted_tsv,
        meme_tsv,
        gene_variability_tsv
    )

    emit:
    report_html = RENDER_PSEL_REPORT.out.report_html
}
