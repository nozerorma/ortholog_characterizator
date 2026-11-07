#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ENTROPY } from '../subworkflows/variability/entropy'
include { PFAM    } from '../subworkflows/variability/pfam'

workflow VARIABILITY {

    take:
    prot_ch   // tuple(gene_id, prot_fasta) — translated protein alignments
    cds_ch    // tuple(gene, codon alignment) — unfiltered CDS set; Pfam is scanned on its full-length human protein
    map_ch    // tuple(gene, map_tsv) — BMGE position maps; project the Pfam hits onto the trimmed alignment

    main:

    ENTROPY( prot_ch )

    if (params.run_hmmscan || params.hmmscan_domtblout) {
        PFAM(
            ENTROPY.out.entropy_tsv,
            cds_ch,
            map_ch
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

