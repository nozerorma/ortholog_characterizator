#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENTROPY } from '../subworkflows/variability/entropy'
include { PFAM    } from '../subworkflows/variability/pfam'

workflow VARIABILITY {

    take:
    prot_ch   // tuple(gene_id, prot_fasta) — translated protein alignments

    main:

    ENTROPY( prot_ch )

    if (params.run_hmmscan || params.hmmscan_domtblout) {
        PFAM(
            prot_ch,
            ENTROPY.out.entropy_tsv.collect { it[1] }
        )
        domain_tsv = PFAM.out.domain_tsv
    } else {
        domain_tsv = Channel.value(file('NO_DOMAIN_TSV'))
    }

    emit:
    gene_variability_tsv  = ENTROPY.out.gene_variability_tsv
    clade_variability_tsv = ENTROPY.out.clade_variability_tsv
    clade_outliers_tsv    = ENTROPY.out.clade_outliers_tsv
    entropy_tsvs          = ENTROPY.out.entropy_tsv
    clade_tsvs            = ENTROPY.out.clade_tsv
    prot_files            = prot_ch.collect { it[1] }
    domain_tsv            = domain_tsv
}

