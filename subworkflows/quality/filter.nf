#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ── Processes ─────────────────────────────────────────────────────────────────

process CHI2_FILTER_BATCH {
    label 'process_batch_light'
    publishDir "${params.out_dir}/PREFILTER/chi2",           mode: 'copy', pattern: '*.filtered.fa',           saveAs: { it.replace('.filtered.fa', '.fa') }
    publishDir "${params.out_dir}/PREFILTER/chi2_ambiguous", mode: 'copy', pattern: 'ambiguous/*.filtered.fa', saveAs: { it.replace('.filtered.fa', '.fa') }
    publishDir "${params.out_dir}/PREFILTER/chi2_low_seq",   mode: 'copy', pattern: 'low_seq/*.filtered.fa',   saveAs: { it.replace('.filtered.fa', '.fa') }
    publishDir "${params.out_dir}/quality/chi2_tests",       mode: 'copy', pattern: '*.chi2.txt'
    publishDir "${params.out_dir}/quality/chi2_removed",     mode: 'copy', pattern: '*.removed.txt'

    input:
    tuple val(batch_id), val(manifest_text), path(fastas)

    output:
    tuple val(batch_id), path("*.filtered.fa"),            optional: true, emit: filtered_fa
    tuple val(batch_id), path("ambiguous/*.filtered.fa"),  optional: true, emit: filtered_fa_ambiguous
    tuple val(batch_id), path("low_seq/*.filtered.fa"),    optional: true, emit: filtered_fa_low_seq
    tuple val(batch_id), path("*.chi2.txt"),               optional: true, emit: chi2_stats
    tuple val(batch_id), path("*.removed.txt"),            optional: true, emit: removed

    script:
    def workers = params.filter_batch_workers as int
    """
    cat > manifest.tsv << 'MANIFEST_EOF'
${manifest_text}
MANIFEST_EOF

    bash ${projectDir}/subworkflows/quality/local/run_chisq_batch.sh \\
        --manifest          manifest.tsv \\
        --workers           ${workers} \\
        --max-iter          ${params.chi2_max_iter} \\
        --ambiguity-thresh  ${params.chi2_ambiguity_thresh} \\
        --min-seq-frac      ${params.chi2_min_seq_frac} \\
        --min-seqs          ${params.min_seqs} \\
        --scripts-dir       "${projectDir}/subworkflows/quality/local"
    """
}

process THRESHOLD_FILTER_BATCH {
    label 'process_batch_light'
    publishDir "${params.out_dir}/PREFILTER/threshold", mode: 'copy', pattern: '*.filtered.fa', saveAs: { it.replace('.filtered.fa', '.fa') }

    input:
    tuple val(batch_id), val(manifest_text), path(fastas)

    output:
    tuple val(batch_id), path("*.filtered.fa"), optional: true, emit: filtered_fa

    script:
    def workers = params.filter_batch_workers as int
    """
    cat > manifest.tsv << 'MANIFEST_EOF'
${manifest_text}
MANIFEST_EOF

    bash ${projectDir}/subworkflows/quality/local/run_threshold_filter_batch.sh \\
        --manifest      manifest.tsv \\
        --workers       ${workers} \\
        --min-seq-frac  ${params.threshold_min_seq_frac} \\
        --min-seqs      ${params.min_seqs} \\
        --scripts-dir   "${projectDir}/subworkflows/quality/local"
    """
}

// ── Subworkflow ───────────────────────────────────────────────────────────────

workflow FILTER {

    take:
    dedup_ch   // tuple(gene_id, fasta)

    main:

    def batch_size = (params.filter_batch_size ?: 500) as int

    batches_ch = dedup_ch
        .collate(batch_size)
        .map { batch ->
            def batch_id      = "filter_batch_${batch.collect { it[0] }.join('_').take(32)}_${batch.hashCode().abs()}"
            def manifest_text = batch.collect { gene_id, fa -> "${gene_id}\t${fa.name}" }.join('\n')
            def fastas        = batch.collect { it[1] }
            tuple(batch_id, manifest_text, fastas)
        }

    if (params.filter_mode == 'chi2') {
        CHI2_FILTER_BATCH( batches_ch )
        raw_out = CHI2_FILTER_BATCH.out.filtered_fa
    } else {
        THRESHOLD_FILTER_BATCH( batches_ch )
        raw_out = THRESHOLD_FILTER_BATCH.out.filtered_fa
    }

    filtered_per_gene = raw_out
        .transpose()
        .map { batch_id, fa ->
            def gene_id = fa.name.replace('.filtered.fa', '')
            tuple(gene_id, fa)
        }

    emit:
    filtered = filtered_per_gene   // tuple(gene_id, fasta)
}
