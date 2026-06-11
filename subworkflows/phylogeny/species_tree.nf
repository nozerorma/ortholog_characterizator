#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ── Processes ─────────────────────────────────────────────────────────────────

process PRUNE_BACKBONE {
    label 'process_medium'
    publishDir "${params.out_dir}/SPECIES_TREE/backbone", mode: 'copy'

    input:
    path gene_tree_files
    path backbone_tree

    output:
    path "backbone_pruned.nwk", emit: pruned_backbone

    script:
    """
    python3 ${projectDir}/subworkflows/phylogeny/local/4b.prune_reference_tree.py \\
        --backbone  "${backbone_tree}" \\
        --gene-trees . \\
        --output    backbone_pruned.nwk
    """
}

process ASTRAL_DENOVO {
    label 'process_long'
    publishDir "${params.out_dir}/SPECIES_TREE/astral", mode: 'copy'

    input:
    path gene_tree_files

    output:
    path "species_tree.nwk",        emit: species_tree
    path "species_tree.rooted.nwk", emit: species_tree_rooted
    path "gene_trees.nwk",          emit: gene_trees_cat

    script:
    def outgroup = params.outgroup_taxa ? "--outgroup '${params.outgroup_taxa}'" : ''
    """
    cat *.treefile *.codeml.nwk 2>/dev/null > gene_trees.nwk || true
    [[ -s gene_trees.nwk ]] || { echo "ERROR: no gene trees found" >&2; exit 1; }

    astral4 -i gene_trees.nwk -o species_tree.nwk 2> astral.log \\
        || astral -i gene_trees.nwk -o species_tree.nwk 2> astral.log

    # Root on outgroup taxa or midpoint
    python3 ${projectDir}/subworkflows/phylogeny/local/prune_tree_for_gene.py \\
        --root-only \\
        --tree    species_tree.nwk \\
        --output  species_tree.rooted.nwk \\
        ${outgroup}
    """
}

process ASTRAL_BACKBONE {
    label 'process_long'
    publishDir "${params.out_dir}/SPECIES_TREE/astral_backbone", mode: 'copy'

    input:
    path gene_tree_files
    path pruned_backbone

    output:
    path "species_tree.nwk",        emit: species_tree
    path "species_tree.rooted.nwk", emit: species_tree_rooted
    path "gene_trees.nwk",          emit: gene_trees_cat

    script:
    def outgroup = params.outgroup_taxa ? "--outgroup '${params.outgroup_taxa}'" : ''
    """
    cat *.treefile *.codeml.nwk 2>/dev/null > gene_trees.nwk || true
    [[ -s gene_trees.nwk ]] || { echo "ERROR: no gene trees found" >&2; exit 1; }

    astral4 -i gene_trees.nwk -c "${pruned_backbone}" -o species_tree.nwk 2> astral.log \\
        || astral -i gene_trees.nwk -c "${pruned_backbone}" -o species_tree.nwk 2> astral.log

    python3 ${projectDir}/subworkflows/phylogeny/local/prune_tree_for_gene.py \\
        --root-only \\
        --tree    species_tree.nwk \\
        --output  species_tree.rooted.nwk \\
        ${outgroup}
    """
}

// ── Subworkflow ───────────────────────────────────────────────────────────────

workflow SPECIES_TREE {

    take:
    gene_tree_files_ch   // collected list of gene tree files

    main:

    if (params.species_tree_method == 'backbone') {
        assert params.backbone_tree : "ERROR: --backbone_tree is required for backbone species tree mode"
        backbone_ch = Channel.value(file(params.backbone_tree))
        PRUNE_BACKBONE( gene_tree_files_ch, backbone_ch )
        ASTRAL_BACKBONE( gene_tree_files_ch, PRUNE_BACKBONE.out.pruned_backbone )
        species_tree_out = ASTRAL_BACKBONE.out.species_tree
        species_tree_rooted = ASTRAL_BACKBONE.out.species_tree_rooted
    } else {
        ASTRAL_DENOVO( gene_tree_files_ch )
        species_tree_out = ASTRAL_DENOVO.out.species_tree
        species_tree_rooted = ASTRAL_DENOVO.out.species_tree_rooted
    }

    emit:
    species_tree        = species_tree_out
    species_tree_rooted = species_tree_rooted
}
