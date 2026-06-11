#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ── Processes ─────────────────────────────────────────────────────────────────

process HMMSCAN_BATCH {
    label 'process_hmmscan'
    tag   "${batch_id}"

    input:
    tuple val(batch_id), path(prot_files)
    path pfam_hmm_files  // HMM + pressed binary index files (.h3f .h3i .h3m .h3p)

    output:
    path "hmmscan_${batch_id}.domtblout", emit: domtblout

    script:
    def cpus     = params.hmmscan_cpus as int
    def hmm_name = file(params.pfam_hmm).name
    """
    [[ -f "${hmm_name}.h3i" ]] || hmmpress "${hmm_name}"

    # Extract one representative per gene: human ortholog if present, else most complete sequence
    python3 - <<'PYEOF'
import os, glob, re

out = open("batch_proteins.fa", "w")
for fa in sorted(glob.glob("*.fa")):
    seqs = []
    hdr, buf = None, []
    with open(fa) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if hdr is not None:
                    seqs.append((hdr, ''.join(buf)))
                hdr, buf = line[1:], []
            else:
                buf.append(line.upper())
        if hdr is not None:
            seqs.append((hdr, ''.join(buf)))
    if not seqs:
        continue
    gene = os.path.basename(fa)[:-3]  # strip .fa
    # Prefer human ortholog
    human = next((h, s) for h, s in seqs
                 if re.search(r'hg38|Homo_sapiens', h, re.IGNORECASE)) \
            if any(re.search(r'hg38|Homo_sapiens', h, re.IGNORECASE) for h, s in seqs) \
            else None
    if human:
        _, rep_seq = human
    else:
        # Most complete: fewest gap characters
        _, rep_seq = min(seqs, key=lambda x: x[1].count('-'))
    ungapped = rep_seq.replace('-', '')
    out.write(f'>{gene}\\n{ungapped}\\n')
out.close()
PYEOF

    hmmscan --domtblout hmmscan_${batch_id}.domtblout --cpu ${cpus} -E 1e-3 \\
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
    echo "# [ok]" >> hmmscan.domtblout
    """
}

process MAP_DOMAIN_VARIABILITY {
    label 'process_single'
    publishDir "${params.out_dir}/reports/variability", mode: 'copy'

    input:
    path entropy_files   // *.entropy.tsv files staged in PROT_VAR_RAW subdir structure
    path prot_files      // *.fa protein files staged in work dir
    path domtblout

    output:
    path "domain_variability.tsv", emit: domain_variability

    script:
    // map_domain_variability.py:
    //   --var_dir      PROT_VAR_RAW/ (entropy TSVs)
    //   --prot_raw_dir PROT_RAW/     (protein FASTAs)
    //   --hmmscan_domtblout
    //   --out_tsv
    """
    mkdir -p PROT_VAR_RAW PROT_RAW
    # Move staged entropy TSVs and protein FASTAs into expected subdirs
    for f in *.entropy.tsv; do [[ -f "\$f" ]] && mv "\$f" PROT_VAR_RAW/ || true; done
    for f in *.fa;           do [[ -f "\$f" ]] && mv "\$f" PROT_RAW/     || true; done

    python3 ${projectDir}/subworkflows/variability/local/map_domain_variability.py \\
        --var_dir           PROT_VAR_RAW \\
        --prot_raw_dir      PROT_RAW \\
        --hmmscan_domtblout "${domtblout}" \\
        --out_tsv           domain_variability.tsv
    """
}

// ── Subworkflow ───────────────────────────────────────────────────────────────

workflow PFAM {

    take:
    prot_ch          // tuple(gene_id, prot_fasta)
    entropy_files_ch // collected list of entropy TSV files

    main:

    if (params.hmmscan_domtblout) {
        domtblout_ch = Channel.value(file(params.hmmscan_domtblout))
    } else {
        assert params.pfam_hmm : "ERROR: --pfam_hmm is required when --run_hmmscan is true"

        def batch_size = (params.variability_batch_size ?: 500) as int

        pfam_hmm_ch = Channel.fromPath("${params.pfam_hmm}*", checkIfExists: true).collect()

        batches_ch = prot_ch
            .collate(batch_size)
            .map { batch ->
                def bid   = "hmm_${batch.collect { it[0] }.take(2).join('_')}_${batch.hashCode().abs()}"
                def files = batch.collect { it[1] }
                tuple(bid, files)
            }

        HMMSCAN_BATCH( batches_ch, pfam_hmm_ch )
        MERGE_DOMTBLOUT( HMMSCAN_BATCH.out.domtblout.collect() )
        domtblout_ch = MERGE_DOMTBLOUT.out.domtblout
    }

    MAP_DOMAIN_VARIABILITY(
        entropy_files_ch,
        prot_ch.collect { it[1] },
        domtblout_ch
    )

    emit:
    domain_tsv = MAP_DOMAIN_VARIABILITY.out.domain_variability
    // Note: map_domain_variability.py produces one domain_variability.tsv.
    // A separate gene_pfam.tsv is not produced by this script;
    // use aggregate_variability.py with --pfam_tsv if needed.
}
