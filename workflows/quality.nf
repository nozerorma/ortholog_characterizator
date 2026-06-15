#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENE_COUNTS } from '../subworkflows/quality/count_genes'
include { DEDUP  } from '../subworkflows/quality/dedup'
include { FILTER } from '../subworkflows/quality/filter'

workflow QUALITY {

    take:
    cds_ch   // tuple(gene_id, fasta)

    main:

    GENE_COUNTS( cds_ch )
    gene_counts_out = GENE_COUNTS.out.gene_counts

    if (params.run_dedup) {
        DEDUP( cds_ch, gene_counts_out )
        dedup_out = DEDUP.out.dedup
    } else {
        dedup_out = cds_ch
    }

    if (params.filter_mode != 'none') {
        FILTER( dedup_out )
        filtered_out = FILTER.out.filtered
    } else {
        filtered_out = dedup_out
    }

    emit:
    dedup    = dedup_out      // deduplicated alignments (before quality filter)
    filtered = filtered_out   // quality-filtered alignments (fed to translation)
    gene_counts = gene_counts_out
}
