#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process RENDER_VARIABILITY_REPORT {
    label 'process_single'
    publishDir "${params.out_dir}/reports", mode: 'copy'

    input:
    path gene_variability_tsv
    path clade_variability_tsv
    path clade_outliers_tsv
    path entropy_files
    path prot_files
    path domain_tsv
    path sel_variability_tsv
    path meme_tsv
    path ucr_windows_tsv
    path ucr_positions_tsv
    path ucr_clade_variability_tsv

    output:
    path "variability_report.html", emit: report_html

    script:
    def taxid     = params.taxid_tsv              ? params.taxid_tsv              : ''
    def fam_order = params.family_order_tsv       ? params.family_order_tsv       : ''
    def sp_tree   = params.psel_species_tree      ? params.psel_species_tree      : ''
    def pfam_desc = params.pfam_descriptions_tsv  ? params.pfam_descriptions_tsv  : ''
    """
    cp ${projectDir}/subworkflows/variability/local/variability_report.qmd ./

    GENE_VARIABILITY_TSV="${gene_variability_tsv}" \\
    CLADE_VARIABILITY_TSV="${clade_variability_tsv}" \\
    CLADE_OUTLIERS_TSV="${clade_outliers_tsv}" \\
    ENTROPY_DIR=\$(pwd) \\
    PROT_DIR=\$(pwd) \\
    DOMAIN_TSV="${domain_tsv}" \\
    PFAM_DESC_TSV="${pfam_desc}" \\
    TAXID_TSV="${taxid}" \\
    FAMILY_ORDER_TSV="${fam_order}" \\
    SPECIES_TREE="${sp_tree}" \\
    SEL_VARIABILITY_TSV="${sel_variability_tsv}" \\
    MEME_TSV="${meme_tsv}" \\
    UCR_WINDOWS_TSV="${ucr_windows_tsv}" \\
    UCR_POSITIONS_TSV="${ucr_positions_tsv}" \\
    UCR_CLADE_VARIABILITY_TSV="${ucr_clade_variability_tsv}" \\
    quarto render variability_report.qmd \\
        --output variability_report.html \\
        --output-dir .
    """
}

workflow VARIABILITY_REPORT {

    take:
    gene_variability_tsv
    clade_variability_tsv
    clade_outliers_tsv
    entropy_files_ch
    prot_files_ch
    domain_tsv_ch
    sel_variability_tsv_ch
    meme_tsv_ch
    ucr_windows_ch
    ucr_positions_ch
    ucr_clade_variability_ch

    main:

    RENDER_VARIABILITY_REPORT(
        gene_variability_tsv,
        clade_variability_tsv.ifEmpty(file('NO_CLADE_VAR')),
        clade_outliers_tsv.ifEmpty(file('NO_CLADE_OUTLIERS')),
        entropy_files_ch.ifEmpty([]).collect(),
        prot_files_ch.ifEmpty([]).collect(),
        domain_tsv_ch,
        sel_variability_tsv_ch.ifEmpty(file('NO_SEL_VAR')),
        meme_tsv_ch.ifEmpty(file('NO_MEME')),
        ucr_windows_ch.ifEmpty(file('NO_UCR_WINDOWS')),
        ucr_positions_ch.ifEmpty(file('NO_UCR_POSITIONS')),
        ucr_clade_variability_ch.ifEmpty(file('NO_UCR_CLADE_VAR'))
    )

    emit:
    report_html = RENDER_VARIABILITY_REPORT.out.report_html
}
