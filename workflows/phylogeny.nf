#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENE_TREES   } from '../subworkflows/phylogeny/gene_trees'
include { SPECIES_TREE } from '../subworkflows/phylogeny/species_tree'

workflow PHYLOGENY {

    take:
    prot_ch   // tuple(gene_id, prot_fasta)

    main:

    GENE_TREES( prot_ch )

    SPECIES_TREE( GENE_TREES.out.gene_trees.collect { it[1] } )

    emit:
    gene_trees   = GENE_TREES.out.gene_trees
    species_tree = SPECIES_TREE.out.species_tree
}
