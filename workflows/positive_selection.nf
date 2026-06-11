#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HYPHY } from '../subworkflows/positive_selection/hyphy'

workflow POSITIVE_SELECTION {

    take:
    trim_ch         // tuple(gene_id, fasta) — BMGE-processed CDS alignments
    map_tsv_ch      // tuple(gene_id, map_tsv) — BMGE position maps (or empty channel)
    gene_var_tsv_ch // path — gene_variability.tsv for sel_variability join (or NO_ sentinel)

    main:

    assert params.psel_species_tree : "ERROR: --psel_species_tree is required for the positive_selection workflow"

    HYPHY( trim_ch, map_tsv_ch, gene_var_tsv_ch )

    emit:
    fubar_tsv           = HYPHY.out.fubar_tsv
    fubar_sites         = HYPHY.out.fubar_sites   // per-site posteriors for UCR cross-reference
    fel_tsv             = HYPHY.out.fel_tsv        // per-site FEL results for UCR cross-reference
    meme_tsv            = HYPHY.out.meme_tsv
    sel_variability_tsv = HYPHY.out.sel_variability_tsv
}
