#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ── Processes ─────────────────────────────────────────────────────────────────

process COUNT_GENES {
    label 'process_medium'
    publishDir "${params.out_dir}/quality/gene_counts", mode: 'copy'

    input:
    path fastas  // collected fasta files staged in work dir

    output:
    path "gene_counts.tsv", emit: gene_counts

    script:
    """
    python3 ${projectDir}/subworkflows/quality/local/count_genes.py \\
        --dir    "." \\
        --output gene_counts.tsv
    """
}

// ── Subworkflow ───────────────────────────────────────────────────────────────

workflow GENE_COUNTS {

    take:
    cds_ch   // tuple(gene_id, fasta)

    main:

    // Collect only the (potentially sampled) fasta files so gene counting sees
    // the same gene pool used downstream, not the full CDS directory.
    cds_files_ch = cds_ch.map { gene_id, fa -> fa }.collect()
    COUNT_GENES( cds_files_ch )

    emit:
    gene_counts = COUNT_GENES.out.gene_counts
}
