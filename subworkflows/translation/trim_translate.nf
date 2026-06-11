#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ── Processes ─────────────────────────────────────────────────────────────────

process BMGE_BATCH {
    label 'process_batch_light'
    publishDir "${params.out_dir}/TRIM/bmge",  mode: 'copy', pattern: '*.trim.fa',    saveAs: { it.replace('.trim.fa', '.fa') }
    publishDir "${params.out_dir}/TRACK/bmge", mode: 'copy', pattern: '*.track.html', saveAs: { it.replace('.track.html', '.html') }
    publishDir "${params.out_dir}/PROT/bmge",  mode: 'copy', pattern: '*.prot.fa',    saveAs: { it.replace('.prot.fa', '.fa') }

    input:
    tuple val(batch_id), val(manifest_text), path(fastas)

    output:
    tuple val(batch_id), path("*.trim.fa"),    optional: true, emit: trim_fa
    tuple val(batch_id), path("*.track.html"), optional: true, emit: track_html
    tuple val(batch_id), path("*.prot.fa"),    optional: true, emit: prot_fa

    script:
    def workers    = params.translation_batch_workers as int
    def gff_arg    = params.gff_file       ? "--gff '${params.gff_file}'"              : ''
    def equiv_arg  = params.gene_equiv_tsv ? "--gene-equiv '${params.gene_equiv_tsv}'" : ''
    """
    cat > manifest.tsv << 'MANIFEST_EOF'
${manifest_text}
MANIFEST_EOF

    bash ${projectDir}/subworkflows/translation/local/run_bmge_batch.sh \\
        --manifest        manifest.tsv \\
        --workers         ${workers} \\
        --gap-thresh      ${params.bmge_gap_threshold} \\
        --entropy-thresh  ${params.bmge_entropy_threshold} \\
        --matrix          ${params.bmge_matrix} \\
        ${gff_arg} ${equiv_arg} \\
        --scripts-dir "${projectDir}/subworkflows/translation/local"
    """
}

process BMGE_MAP_BATCH {
    label 'process_batch_light'
    publishDir "${params.out_dir}/MAP", mode: 'copy', pattern: '*.map.tsv'

    input:
    // manifest has 3 cols: gene_id <TAB> track_html <TAB> cds_fasta
    tuple val(batch_id), val(manifest_text), path(track_htmls), path(cds_fastas)

    output:
    tuple val(batch_id), path("*.map.tsv"), optional: true, emit: map_tsv

    script:
    def workers   = params.translation_batch_workers as int
    def gff_arg   = params.gff_file       ? "--gff '${params.gff_file}'"              : ''
    def equiv_arg = params.gene_equiv_tsv ? "--gene-equiv '${params.gene_equiv_tsv}'" : ''
    """
    cat > manifest.tsv << 'MANIFEST_EOF'
${manifest_text}
MANIFEST_EOF

    bash ${projectDir}/subworkflows/translation/local/run_bmge_map_batch.sh \\
        --manifest    manifest.tsv \\
        --workers     ${workers} \\
        ${gff_arg} ${equiv_arg} \\
        --scripts-dir "${projectDir}/subworkflows/translation/local"
    """
}

// ── Helpers ───────────────────────────────────────────────────────────────────

def make_batches(ch, prefix, batch_size) {
    return ch
        .collate(batch_size)
        .map { batch ->
            def bid           = "${prefix}_${batch.collect { it[0] }.take(2).join('_')}_${batch.hashCode().abs()}"
            def manifest_text = batch.collect { gene_id, fa -> "${gene_id}\t${fa.name}" }.join('\n')
            def files         = batch.collect { it[1] }
            tuple(bid, manifest_text, files)
        }
}

// ── Subworkflow ───────────────────────────────────────────────────────────────

workflow TRIM_TRANSLATE {

    take:
    filtered_ch   // tuple(gene_id, fasta)

    main:

    def batch_size = (params.translation_batch_size ?: 200) as int

    batches_ch = make_batches(filtered_ch, 'bmge', batch_size)
    BMGE_BATCH( batches_ch )

    prot_out = BMGE_BATCH.out.prot_fa
        .transpose()
        .map { bid, fa -> tuple(fa.name.replace('.prot.fa',''), fa) }

    trim_out = BMGE_BATCH.out.trim_fa
        .transpose()
        .map { bid, fa -> tuple(fa.name.replace('.trim.fa',''), fa) }

    // Position mapping: join track HTML with original CDS fasta per gene
    track_per_gene = BMGE_BATCH.out.track_html
        .transpose()
        .map { bid, f -> tuple(f.name.replace('.track.html',''), f) }

    track_with_cds = track_per_gene.join(filtered_ch, by: 0)

    map_batches = track_with_cds
        .collate(batch_size)
        .map { batch ->
            def bid           = "map_${batch.hashCode().abs()}"
            def manifest_text = batch.collect { g, tr, fa -> "${g}\t${tr.name}\t${fa.name}" }.join('\n')
            def tracks        = batch.collect { it[1] }
            def fastas        = batch.collect { it[2] }
            tuple(bid, manifest_text, tracks, fastas)
        }

    BMGE_MAP_BATCH( map_batches )

    map_out = BMGE_MAP_BATCH.out.map_tsv
        .transpose()
        .map { bid, f -> tuple(f.name.replace('.map.tsv',''), f) }

    emit:
    prot    = prot_out
    trim    = trim_out
    map_tsv = map_out
}
