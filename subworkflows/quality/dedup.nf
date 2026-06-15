#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ── Processes ─────────────────────────────────────────────────────────────────

process DEDUP_BATCH {
    label 'process_batch_light'
    publishDir "${params.out_dir}/quality/dedup_logs", mode: 'copy', pattern: '*.dedup.log',        saveAs: { it.replace('.dedup.log', '.log') }
    publishDir "${params.out_dir}/quality/dedup_reports", mode: 'copy', pattern: '*.dedup_report.tsv', saveAs: { it.replace('.dedup_report.tsv', '_dedup_report.tsv') }

    input:
    tuple val(batch_id), val(manifest_text), path(fastas), path(gene_counts)

    output:
    tuple val(batch_id), path("*.dedup.fa"),           optional: true, emit: dedup_fa
    tuple val(batch_id), path("*.dedup.log"),          optional: true, emit: dedup_log
    tuple val(batch_id), path("*.dedup_report.tsv"),   optional: true, emit: dedup_report

    script:
    assert params.assemblies_tsv    : "ERROR: --assemblies_tsv is required for deduplication"
    assert params.species_dedup_map : "ERROR: --species_dedup_map is required for deduplication"
    def workers = params.dedup_batch_workers as int
    """
    cat > manifest.tsv << 'MANIFEST_EOF'
${manifest_text}
MANIFEST_EOF

    bash ${projectDir}/subworkflows/quality/local/run_dedup_batch.sh \\
        --manifest    manifest.tsv \\
        --workers     ${workers} \\
        --gene-counts "${gene_counts}" \\
        --metadata    "${params.assemblies_tsv}" \\
        --dedup-map   "${params.species_dedup_map}" \\
        --scripts-dir "${projectDir}/subworkflows/quality/local"
    """
}

process RENAME_ASSEMBLIES {
    label 'process_batch_light'
    publishDir "${params.out_dir}/DEDUP", mode: 'copy', pattern: 'renamed/*.fa', saveAs: { it.replace('renamed/', '').replace('.dedup.fa', '.fa') }

    input:
    tuple val(batch_id), path(dedup_fastas)

    output:
    tuple val(batch_id), path("renamed/*.fa"), optional: true, emit: renamed_fa

    script:
    def taxids_arg = params.taxid_tsv ? "--taxids '${params.taxid_tsv}'" : ''
    """
    mkdir -p renamed
    python3 ${projectDir}/subworkflows/quality/local/rename_assemblies_mammals.py \\
        --in_dir   . \\
        --out_dir  renamed \\
        --metadata "${params.assemblies_tsv}" \\
        ${taxids_arg}
    """
}

// ── Subworkflow ───────────────────────────────────────────────────────────────

workflow DEDUP {

    take:
    cds_ch   // tuple(gene_id, fasta)
    gene_counts_ch

    main:

    // Build batches: collate N genes, create manifest + file list
    def batch_size = (params.dedup_batch_size ?: 500) as int

    batches_ch = cds_ch
        .collate(batch_size)
        .map { batch ->
            def batch_id      = "dedup_batch_${batch.collect { it[0] }.join('_').take(32)}_${batch.hashCode().abs()}"
            def manifest_text = batch.collect { gene_id, fa -> "${gene_id}\t${fa.name}" }.join('\n')
            def fastas        = batch.collect { it[1] }
            tuple(batch_id, manifest_text, fastas)
        }

    // Combine each batch with the precomputed gene_counts file
    batches_with_counts = batches_ch.combine(gene_counts_ch)

    DEDUP_BATCH( batches_with_counts )

    RENAME_ASSEMBLIES( DEDUP_BATCH.out.dedup_fa )

    // Reconstruct per-gene channel from batch outputs
    dedup_per_gene = RENAME_ASSEMBLIES.out.renamed_fa
        .transpose()
        .map { batch_id, fa ->
            def gene_id = fa.name.replace('.dedup.fa', '')
            tuple(gene_id, fa)
        }

    emit:
    dedup = dedup_per_gene   // tuple(gene_id, fasta) — headers are species names
}
