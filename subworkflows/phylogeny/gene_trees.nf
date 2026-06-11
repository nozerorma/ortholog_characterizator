#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ── Processes ─────────────────────────────────────────────────────────────────

process IQTREE_BATCH {
    label 'process_batch_heavy'
    publishDir "${params.out_dir}/TREES", mode: 'copy', pattern: '*.treefile'
    publishDir "${params.out_dir}/TREES", mode: 'copy', pattern: '*.iqtree'

    input:
    tuple val(batch_id), val(manifest_text), path(fastas)

    output:
    tuple val(batch_id), path("*.treefile"), optional: true, emit: treefiles
    tuple val(batch_id), path("*.iqtree"),   optional: true, emit: iqtree_logs

    script:
    def workers   = params.genetree_batch_workers as int
    def cpu_each  = Math.max(1, (task.cpus as int).intdiv(workers))
    def mset_arg  = params.iqtree_mset ? "--mset ${params.iqtree_mset}" : ""
    """
    cat > manifest.tsv << 'MANIFEST_EOF'
${manifest_text}
MANIFEST_EOF

    bash ${projectDir}/subworkflows/phylogeny/local/run_iqtree_batch.sh \\
        --manifest    manifest.tsv \\
        --workers     ${workers} \\
        --cpu-each    ${cpu_each} \\
        --model       ${params.iqtree_model} \\
        ${mset_arg} \\
        --bootstrap   ${params.iqtree_bootstrap} \\
        --scripts-dir "${projectDir}/subworkflows/phylogeny/local"
    """
}

// ── Subworkflow ───────────────────────────────────────────────────────────────

workflow GENE_TREES {

    take:
    prot_ch   // tuple(gene_id, fasta)

    main:

    def batch_size = (params.genetree_batch_size ?: 50) as int

    batches_ch = prot_ch
        .collate(batch_size)
        .map { batch ->
            def bid           = "gtree_${batch.collect { it[0] }.take(2).join('_')}_${batch.hashCode().abs()}"
            def manifest_text = batch.collect { g, f -> "${g}\t${f.name}" }.join('\n')
            def files         = batch.collect { it[1] }
            tuple(bid, manifest_text, files)
        }

    IQTREE_BATCH( batches_ch )

    gene_trees_ch = IQTREE_BATCH.out.treefiles
        .transpose()
        .map { bid, f -> tuple(f.name.replace('.treefile', ''), f) }

    emit:
    gene_trees = gene_trees_ch   // tuple(gene_id, treefile)
}
