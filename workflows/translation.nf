#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { TRIM_TRANSLATE } from '../subworkflows/translation/trim_translate'

workflow TRANSLATION {

    take:
    filtered_ch   // tuple(gene_id, fasta) — quality-filtered alignments

    main:

    TRIM_TRANSLATE( filtered_ch )

    emit:
    prot    = TRIM_TRANSLATE.out.prot
    trim    = TRIM_TRANSLATE.out.trim
    map_tsv = TRIM_TRANSLATE.out.map_tsv
}
