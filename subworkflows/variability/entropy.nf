#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ── Processes ─────────────────────────────────────────────────────────────────

process VARIABILITY_BATCH {
    label 'process_batch_light'
    // compute_variability.py writes to PROT_VAR_RAW/<gene>.*
    publishDir "${params.out_dir}/PROT_VAR_RAW", mode: 'copy', pattern: 'PROT_VAR_RAW/*.tsv', saveAs: { it.replace('PROT_VAR_RAW/', '') }
    publishDir "${params.out_dir}/PROT_VAR_RAW", mode: 'copy', pattern: 'PROT_VAR_RAW/*.fa',  saveAs: { it.replace('PROT_VAR_RAW/', '') }

    input:
    tuple val(batch_id), val(manifest_text), path(fastas)

    output:
    tuple val(batch_id), path("PROT_VAR_RAW/*.entropy.tsv"),       optional: true, emit: entropy_tsv
    tuple val(batch_id), path("PROT_VAR_RAW/*.clade_entropy.tsv"), optional: true, emit: clade_tsv

    script:
    assert params.taxid_tsv : "ERROR: --taxid_tsv is required for the variability workflow (compute_variability.py)"
    def workers   = params.variability_batch_workers as int
    def fam_arg   = params.family_order_tsv ? "--family-order-tsv '${params.family_order_tsv}'" : ''
    def producerHash = ['compute_variability.py', 'run_variability_batch.sh'].collect { name ->
        java.security.MessageDigest.getInstance('SHA-256')
            .digest(file("${projectDir}/subworkflows/variability/local/${name}").bytes)
            .encodeHex().toString() }.join(' ')
    """
    # Entropy producer SHA256: ${producerHash}
    cat > manifest.tsv << 'MANIFEST_EOF'
${manifest_text}
MANIFEST_EOF

    bash ${projectDir}/subworkflows/variability/local/run_variability_batch.sh \\
        --manifest    manifest.tsv \\
        --workers     ${workers} \\
        --alpha       ${params.valdar_alpha} \\
        --beta        ${params.valdar_beta} \\
        --gamma       ${params.valdar_gamma} \\
        --taxid-tsv   "${params.taxid_tsv}" \\
        ${fam_arg} \\
        --scripts-dir "${projectDir}/subworkflows/variability/local"
    """
}

process AGGREGATE_VARIABILITY {
    label 'process_single'
    publishDir "${params.out_dir}/reports/variability", mode: 'copy'

    input:
    path tsv_files   // all *.entropy.tsv and *.clade_entropy.tsv staged in work dir

    output:
    path "gene_variability.tsv",   emit: gene_variability
    path "clade_variability.tsv",  optional: true, emit: clade_variability
    path "clade_outliers.tsv",     optional: true, emit: clade_outliers

    script:
    // aggregate_variability.py: --var_dir (PROT_VAR_RAW dir), --out_dir, optional --pfam_desc_tsv
    def pfam_desc = params.pfam_descriptions_tsv ? "--pfam_desc_tsv '${params.pfam_descriptions_tsv}'" : ''
    def z_thresh  = params.zscore_threshold
    """
    python3 ${projectDir}/subworkflows/variability/local/aggregate_variability.py \\
        --var_dir  . \\
        --out_dir  . \\
        --outlier_z_threshold ${z_thresh} \\
        ${pfam_desc}
    """
}

// ── Subworkflow ───────────────────────────────────────────────────────────────

workflow ENTROPY {

    take:
    prot_ch   // tuple(gene_id, prot_fasta)

    main:

    def batch_size = (params.variability_batch_size ?: 500) as int

    batches_ch = prot_ch
        .collate(batch_size)
        .map { batch ->
            def bid           = "var_${batch.collect { it[0] }.take(2).join('_')}_${batch.hashCode().abs()}"
            def manifest_text = batch.collect { g, f -> "${g}\t${f.name}" }.join('\n')
            def files         = batch.collect { it[1] }
            tuple(bid, manifest_text, files)
        }

    VARIABILITY_BATCH( batches_ch )

    entropy_ch = VARIABILITY_BATCH.out.entropy_tsv
        .transpose()
        .map { bid, f -> tuple(f.name.replace('.entropy.tsv',''), f) }

    clade_ch = VARIABILITY_BATCH.out.clade_tsv
        .transpose()
        .map { bid, f -> tuple(f.name.replace('.clade_entropy.tsv',''), f) }

    AGGREGATE_VARIABILITY(
        VARIABILITY_BATCH.out.entropy_tsv.collect { it[1] }.flatten()
            .mix( VARIABILITY_BATCH.out.clade_tsv.collect { it[1] }.flatten() )
            .collect()
    )

    emit:
    entropy_tsv           = entropy_ch
    clade_tsv             = clade_ch
    gene_variability_tsv  = AGGREGATE_VARIABILITY.out.gene_variability
    clade_variability_tsv = AGGREGATE_VARIABILITY.out.clade_variability
    clade_outliers_tsv    = AGGREGATE_VARIABILITY.out.clade_outliers
}
