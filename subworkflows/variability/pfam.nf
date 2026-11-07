#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ── Processes ─────────────────────────────────────────────────────────────────

process HMMSCAN_BATCH {
    label 'process_hmmscan'
    tag   "${batch_id}"

    input:
    tuple val(batch_id), val(manifest_text), path(cds_files)
    path pfam_hmm_files  // HMM + pressed binary index files (.h3f .h3i .h3m .h3p)

    output:
    path "hmmscan_${batch_id}.domtblout", emit: domtblout

    script:
    def cpus     = params.hmmscan_cpus as int
    def hmm_name = file(params.pfam_hmm).name
    def queryHash = java.security.MessageDigest.getInstance('SHA-256')
        .digest(file("${projectDir}/subworkflows/variability/local/full_length_query.py").bytes)
        .encodeHex().toString()
    """
    # Full-length query producer SHA256: ${queryHash}
    cat > manifest.tsv << 'MANIFEST_EOF'
${manifest_text}
MANIFEST_EOF

    [[ -f "${hmm_name}.h3i" ]] || hmmpress "${hmm_name}"

    # One query per gene: the human protein translated from the unfiltered CDS alignment, gaps removed
    python3 ${projectDir}/subworkflows/variability/local/full_length_query.py \\
        --manifest manifest.tsv --out_fa batch_proteins.fa

    # Pfam gathering thresholds, the curated per-family criterion of Pfam (exclusive with -E / -T)
    hmmscan --cut_ga --domtblout hmmscan_${batch_id}.domtblout --cpu ${cpus} \\
        "${hmm_name}" batch_proteins.fa > hmmscan_${batch_id}.log
    """
}

process MERGE_DOMTBLOUT {
    label 'process_single'
    publishDir "${params.out_dir}/HMMSCAN", mode: 'copy'

    input:
    path domtblout_files

    output:
    path "hmmscan.domtblout", emit: domtblout

    script:
    """
    first=\$(ls *.domtblout | sort | head -1)
    grep '^#' "\$first" | head -3 > hmmscan.domtblout
    for f in \$(ls *.domtblout | sort); do
        grep -v '^#' "\$f" >> hmmscan.domtblout || true
    done
    # the option settings record the reporting criterion (--cut_ga) for the exporter
    grep -m1 '^# Option settings' "\$first" >> hmmscan.domtblout || true
    echo "# [ok]" >> hmmscan.domtblout
    """
}

process MAP_DOMAIN_VARIABILITY {
    label 'process_single'
    publishDir "${params.out_dir}/reports/variability", mode: 'copy'

    input:
    tuple val(manifest_text), path(entropy_files), path(map_files), path(cds_files)
    path domtblout

    output:
    path "domain_variability.tsv", emit: domain_variability

    script:
    def sources = ['map_domain_variability.py', 'full_length_query.py'].collect { name ->
        java.security.MessageDigest.getInstance('SHA-256')
            .digest(file("${projectDir}/subworkflows/variability/local/${name}").bytes)
            .encodeHex().toString() }.join(' ')
    def clansArg = params.pfam_clans_tsv ? "--pfam_clans_tsv '${params.pfam_clans_tsv}'" : ''
    """
    # Pfam coordinate producer SHA256: ${sources}
    cat > manifest.tsv << 'MANIFEST_EOF'
${manifest_text}
MANIFEST_EOF

    # Entropy tables, MAP tables and codon alignments are staged side by side in the work dir
    python3 ${projectDir}/subworkflows/variability/local/map_domain_variability.py \\
        --manifest          manifest.tsv \\
        --var_dir           . \\
        --map_dir           . \\
        --cds_dir           . \\
        --alignment_coordinate_space trimmed_alignment_1based \\
        --hmmscan_domtblout "${domtblout}" \\
        --out_tsv           domain_variability.tsv \\
        ${clansArg}
    """
}

// ── Subworkflow ───────────────────────────────────────────────────────────────

workflow PFAM {

    take:
    entropy_ch       // tuple(gene_id, entropy_tsv); gene_id is <gene>, or <gene>.prot for tables of older runs
    cds_ch           // tuple(gene, codon alignment of the unfiltered CDS set)
    map_ch           // tuple(gene, BMGE position map)

    main:

    // gene_id keeps the id of the entropy tables, which every downstream join uses; CDS and MAP files
    // are named after the bare gene, so the .prot label that older entropy tables carry is dropped to join.
    genes_ch = entropy_ch
        .map { id, f -> tuple(id.replaceAll(/\.prot$/, ''), id, f) }
        .join(cds_ch)
        .join(map_ch)
        .map { gene, id, entropy, cds, mapf -> tuple(id, entropy, cds, mapf) }

    if (params.hmmscan_domtblout) {
        domtblout_ch = Channel.value(file(params.hmmscan_domtblout))
    } else {
        assert params.pfam_hmm : "ERROR: --pfam_hmm is required when --run_hmmscan is true"

        def batch_size = (params.variability_batch_size ?: 500) as int

        pfam_hmm_ch = Channel.fromPath("${params.pfam_hmm}*", checkIfExists: true).collect()

        batches_ch = genes_ch
            .collate(batch_size)
            .map { batch ->
                def bid           = "hmm_${batch.collect { it[0] }.take(2).join('_')}_${batch.hashCode().abs()}"
                def manifest_text = batch.collect { id, entropy, cds, mapf -> "${id}\t${cds.name}" }.join('\n')
                tuple(bid, manifest_text, batch.collect { it[2] })
            }

        HMMSCAN_BATCH( batches_ch, pfam_hmm_ch )
        MERGE_DOMTBLOUT( HMMSCAN_BATCH.out.domtblout.collect() )
        domtblout_ch = MERGE_DOMTBLOUT.out.domtblout
    }

    all_genes_ch = genes_ch
        .toList()
        .map { genes ->
            tuple(genes.collect { id, entropy, cds, mapf -> "${id}\t${entropy.name}\t${mapf.name}\t${cds.name}" }.join('\n'),
                  genes.collect { it[1] },
                  genes.collect { it[3] },
                  genes.collect { it[2] })
        }

    MAP_DOMAIN_VARIABILITY( all_genes_ch, domtblout_ch )

    emit:
    domain_tsv = MAP_DOMAIN_VARIABILITY.out.domain_variability
    // Note: map_domain_variability.py produces one domain_variability.tsv.
    // A separate gene_pfam.tsv is not produced by this script;
    // use aggregate_variability.py with --pfam_tsv if needed.
}
