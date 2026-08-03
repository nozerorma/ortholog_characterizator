#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ── Processes ─────────────────────────────────────────────────────────────────

process DEDUP_IDENTICAL_BATCH {
    label 'process_batch_light'

    input:
    tuple val(batch_id), val(manifest_text), path(fastas)

    output:
    tuple val(batch_id), path("*.dedup.fa"), optional: true, emit: dedup_fa

    script:
    def workers = params.psel_batch_workers as int
    """
    cat > manifest.tsv << 'MANIFEST_EOF'
${manifest_text}
MANIFEST_EOF

    python3 - ${workers} <<'PYEOF'
import sys
import multiprocessing
from Bio import SeqIO

def dedup_one(args):
    gene_id, fasta = args
    seen, records = set(), []
    for r in SeqIO.parse(fasta, "fasta"):
        s = str(r.seq).upper()
        if s not in seen:
            seen.add(s)
            records.append(r)
    if len(records) >= 2:
        SeqIO.write(records, f"{gene_id}.dedup.fa", "fasta")

workers = int(sys.argv[1])
jobs = []
with open("manifest.tsv") as fh:
    for line in fh:
        line = line.strip()
        if line:
            gene_id, fasta = line.split("\\t")
            jobs.append((gene_id, fasta))

with multiprocessing.get_context('fork').Pool(workers) as pool:
    pool.map(dedup_one, jobs)
PYEOF
    """
}

process PRUNE_TREE_BATCH {
    label 'process_batch_light'
    publishDir "${params.out_dir}/PSEL_TREES", mode: 'copy', pattern: '*.pruned.nwk'

    input:
    tuple val(batch_id), val(manifest_text), path(fastas)
    path species_tree

    output:
    tuple val(batch_id), path("*.pruned.nwk"), optional: true, emit: pruned_trees

    script:
    // prune_tree_for_gene.py: --tree, --fasta, --out, [--min_taxa]
    def workers    = params.psel_batch_workers as int
    def taxid_arg  = params.taxid_tsv ? "--taxid-tsv '${params.taxid_tsv}'" : ''
    """
    cat > manifest.tsv << 'MANIFEST_EOF'
${manifest_text}
MANIFEST_EOF

    bash ${projectDir}/subworkflows/positive_selection/local/run_prune_batch.sh \\
        --manifest     manifest.tsv \\
        --workers      ${workers} \\
        --species-tree "${species_tree}" \\
        --scripts-dir  "${projectDir}/subworkflows/positive_selection/local" \\
        ${taxid_arg}
    """
}

process HYPHY_BATCH {
    label 'process_batch_heavy'
    // run_hyphy_batch.sh writes <gene>.<METHOD>.json in work dir
    publishDir "${params.out_dir}/PSEL", mode: 'copy', pattern: '*.*.json'

    input:
    tuple val(batch_id), val(manifest_text), path(fastas), path(trees)

    output:
    tuple val(batch_id), path("*.FUBAR.json"), optional: true, emit: fubar_json
    tuple val(batch_id), path("*.FEL.json"),   optional: true, emit: fel_json
    tuple val(batch_id), path("*.SLAC.json"),  optional: true, emit: slac_json

    script:
    def workers   = params.psel_batch_workers as int
    def cpu_each  = Math.max(1, (task.cpus as int).intdiv(workers))
    def taxid_arg = params.taxid_tsv ? "--taxid-tsv '${params.taxid_tsv}'" : ''
    """
    cat > manifest.tsv << 'MANIFEST_EOF'
${manifest_text}
MANIFEST_EOF

    bash ${projectDir}/subworkflows/positive_selection/local/run_hyphy_batch.sh \\
        --manifest    manifest.tsv \\
        --workers     ${workers} \\
        --cpu-each    ${cpu_each} \\
        --methods     "${params.hyphy_methods}" \\
        --scripts-dir "${projectDir}/subworkflows/positive_selection/local" \\
        ${taxid_arg}
    """
}

process AGGREGATE_FIRST_PASS {
    label 'process_single'
    publishDir "${params.out_dir}/reports/positive_selection", mode: 'copy'

    input:
    path json_files   // all FUBAR/FEL JSON files collected into work dir
    path map_files    // all *.map.tsv files (may be empty list)

    output:
    path "fubar_results.tsv",        emit: fubar_tsv
    path "fubar_sites.tsv",          optional: true, emit: fubar_sites
    path "fel_results.tsv",          optional: true, emit: fel_tsv
    path "meme_candidate_genes.txt", emit: meme_candidates

    script:
    def methods   = params.hyphy_methods
    def threshold = params.fubar_threshold
    def var_arg   = params.variability_tsv ? "--variability_tsv '${params.variability_tsv}'" : ''
    """
    # Arrange flat JSONs into PSEL/<gene>/ structure expected by aggregate script
    python3 - <<'PYEOF'
import os, glob, re
for f in glob.glob("*.FUBAR.json") + glob.glob("*.FEL.json") + glob.glob("*.SLAC.json"):
    m = re.match(r'^(.+?)\\.(?:FUBAR|FEL|SLAC)\\.json\$', f)
    if m:
        gene = m.group(1)
        os.makedirs(f"PSEL/{gene}", exist_ok=True)
        os.rename(f, f"PSEL/{gene}/{f}")
PYEOF

    # Arrange map TSVs into MAP/ so aggregate script can look them up by gene name
    mkdir -p MAP
    for f in *.map.tsv; do
        [[ -f "\$f" ]] && mv "\$f" MAP/ || true
    done
    map_arg=""
    [[ -n "\$(ls MAP/*.map.tsv 2>/dev/null)" ]] && map_arg="--map_dir MAP"

    python3 ${projectDir}/subworkflows/positive_selection/local/aggregate_positive_selection.py \\
        --psel_dir         PSEL \\
        --analyses         "${methods}" \\
        --out_dir          . \\
        --fubar_threshold  ${threshold} \\
        \${map_arg} \\
        ${var_arg}

    # Write meme candidate genes (genes with positive FUBAR sites)
    python3 - <<'PYEOF'
import pandas as pd, sys, os
fubar_f = "fubar_results.tsv"
meme_f  = "meme_candidate_genes.txt"
min_sites = int(${params.meme_min_sites})
if os.path.exists(fubar_f):
    df = pd.read_csv(fubar_f, sep='\\t')
    col = next((c for c in df.columns if 'n_pos' in c.lower() or 'pos_sites' in c.lower()), None)
    gene_col = next((c for c in df.columns if 'gene' in c.lower()), df.columns[0])
    if col:
        candidates = df[df[col] >= min_sites][gene_col].tolist()
        with open(meme_f, 'w') as fh:
            fh.write('\\n'.join(str(g) for g in candidates) + '\\n')
    else:
        open(meme_f, 'w').close()
else:
    open(meme_f, 'w').close()
PYEOF
    """
}

process HYPHY_MEME_BATCH {
    label 'process_batch_heavy'
    publishDir "${params.out_dir}/PSEL", mode: 'copy', pattern: '*.MEME.json'

    input:
    tuple val(batch_id), val(manifest_text), path(fastas), path(trees)

    output:
    tuple val(batch_id), path("*.MEME.json"), optional: true, emit: meme_json

    script:
    def workers  = params.psel_batch_workers as int
    def cpu_each = Math.max(1, (task.cpus as int).intdiv(workers))
    """
    cat > manifest.tsv << 'MANIFEST_EOF'
${manifest_text}
MANIFEST_EOF

    bash ${projectDir}/subworkflows/positive_selection/local/run_hyphy_batch.sh \\
        --manifest    manifest.tsv \\
        --workers     ${workers} \\
        --cpu-each    ${cpu_each} \\
        --methods     "MEME" \\
        --scripts-dir "${projectDir}/subworkflows/positive_selection/local"
    """
}

process AGGREGATE_FINAL {
    label 'process_single'
    publishDir "${params.out_dir}/reports/positive_selection", mode: 'copy'

    input:
    path fubar_tsv
    path meme_json_files   // all MEME JSONs collected
    path map_files         // all *.map.tsv files (may be empty list)
    path gene_var_tsv      // gene_variability.tsv (or NO_ sentinel if unavailable)

    output:
    path "meme_results.tsv",    optional: true, emit: meme_tsv
    path "sel_variability.tsv", optional: true, emit: sel_variability

    script:
    def var_arg = (gene_var_tsv && !gene_var_tsv.name.startsWith('NO_'))
                      ? "--variability_tsv '${gene_var_tsv}'"
                      : ''
    """
    # Arrange MEME JSONs into PSEL/<gene>/ subdirs
    python3 - <<'PYEOF'
import os, glob, re
for f in glob.glob("*.MEME.json"):
    m = re.match(r'^(.+?)\\.MEME\\.json\$', f)
    if m:
        gene = m.group(1)
        os.makedirs(f"PSEL/{gene}", exist_ok=True)
        os.rename(f, f"PSEL/{gene}/{f}")
PYEOF

    mkdir -p MAP
    for f in *.map.tsv; do
        [[ -f "\$f" ]] && mv "\$f" MAP/ || true
    done
    map_arg=""
    [[ -n "\$(ls MAP/*.map.tsv 2>/dev/null)" ]] && map_arg="--map_dir MAP"

    python3 ${projectDir}/subworkflows/positive_selection/local/aggregate_positive_selection.py \\
        --psel_dir        PSEL \\
        --analyses        "MEME" \\
        --out_dir         . \\
        --fubar_tsv       ${fubar_tsv} \\
        \${map_arg} \\
        ${var_arg}

    [[ -f sel_variability.tsv ]] || touch sel_variability.tsv
    """
}

// ── Helper ─────────────────────────────────────────────────────────────────────

def make_psel_batches(ch, prefix, batch_size) {
    return ch
        .collate(batch_size)
        .map { batch ->
            def bid           = "${prefix}_${batch.collect { it[0] }.take(2).join('_')}_${batch.hashCode().abs()}"
            def manifest_text = batch.collect { g, fa, tr -> "${g}\t${fa.name}\t${tr.name}" }.join('\n')
            def fastas        = batch.collect { it[1] }
            def trees         = batch.collect { it[2] }
            tuple(bid, manifest_text, fastas, trees)
        }
}

// ── Subworkflow ───────────────────────────────────────────────────────────────

workflow HYPHY {

    take:
    trim_ch         // tuple(gene_id, fasta) — BMGE-processed CDS alignments
    map_tsv_ch      // tuple(gene_id, map_tsv) — BMGE position maps (or empty channel)
    gene_var_tsv_ch // path — gene_variability.tsv for sel_variability join (or NO_ sentinel)

    main:

    def batch_size  = (params.psel_batch_size ?: 100) as int
    species_tree_ch = Channel.value(file(params.psel_species_tree))

    // Step 0: remove sequences that are identical at nucleotide level
    // (identical sequences cause zero-length branches and HyPhy warnings)
    dedup_batches_ch = trim_ch
        .collate(batch_size)
        .map { batch ->
            def bid           = "dedup_${batch.collect { it[0] }.take(2).join('_')}_${batch.hashCode().abs()}"
            def manifest_text = batch.collect { g, f -> "${g}\t${f.name}" }.join('\n')
            def files         = batch.collect { it[1] }
            tuple(bid, manifest_text, files)
        }

    DEDUP_IDENTICAL_BATCH( dedup_batches_ch )

    cds_ch = DEDUP_IDENTICAL_BATCH.out.dedup_fa
        .transpose()
        .map { bid, f -> tuple(f.name.replace('.dedup.fa',''), f) }

    // Step 1: prune species tree for each gene (batched, 2-col manifest)
    prune_batches_ch = cds_ch
        .collate(batch_size)
        .map { batch ->
            def bid           = "prune_${batch.collect { it[0] }.take(2).join('_')}_${batch.hashCode().abs()}"
            def manifest_text = batch.collect { g, f -> "${g}\t${f.name}" }.join('\n')
            def files         = batch.collect { it[1] }
            tuple(bid, manifest_text, files)
        }

    PRUNE_TREE_BATCH( prune_batches_ch, species_tree_ch )

    pruned_per_gene = PRUNE_TREE_BATCH.out.pruned_trees
        .transpose()
        .map { bid, f -> tuple(f.name.replace('.pruned.nwk',''), f) }

    // Join per-gene fasta with its pruned tree → tuple(gene_id, fasta, tree)
    fasta_with_tree = cds_ch.join(pruned_per_gene, by: 0)

    // Split fasta_with_tree so both HYPHY and MEME can read from it independently
    fasta_with_tree
        .multiMap { gene_id, fa, tree ->
            hyphy: tuple(gene_id, fa, tree)
            meme:  tuple(gene_id, fa, tree)
        }
        .set { fwt }

    // Step 2: first-pass HyPhy (FUBAR + FEL [+ SLAC])
    hyphy_batches_ch = make_psel_batches(fwt.hyphy, 'hyphy', batch_size)
    HYPHY_BATCH( hyphy_batches_ch )

    // Collect all JSON outputs for aggregation
    all_first_pass_json = HYPHY_BATCH.out.fubar_json.collect { it[1] }.flatten()
        .mix(HYPHY_BATCH.out.fel_json.collect { it[1] }.flatten())
        .mix(HYPHY_BATCH.out.slac_json.collect { it[1] }.flatten().ifEmpty([]))
        .collect()

    // Collect all map TSVs as a flat list of files (empty list if no maps provided)
    all_map_files = map_tsv_ch
        .map { gene_id, f -> f }
        .collect()
        .ifEmpty([])

    // Step 3: aggregate first pass (with optional position maps)
    AGGREGATE_FIRST_PASS( all_first_pass_json, all_map_files )

    // Step 4: conditional MEME on FUBAR-positive genes
    meme_candidates_ch = AGGREGATE_FIRST_PASS.out.meme_candidates
        .splitText()
        .map { it.trim() }
        .filter { it != '' }

    meme_input = fwt.meme
        .join(meme_candidates_ch.map { g -> tuple(g, true) }, by: 0)
        .map { gene_id, fa, tree, _flag -> tuple(gene_id, fa, tree) }

    meme_batches_ch = make_psel_batches(meme_input, 'meme', batch_size)
    HYPHY_MEME_BATCH( meme_batches_ch )

    all_meme_json = HYPHY_MEME_BATCH.out.meme_json
        .collect { it[1] }.flatten().ifEmpty([]).collect()

    AGGREGATE_FINAL(
        AGGREGATE_FIRST_PASS.out.fubar_tsv,
        all_meme_json,
        all_map_files,
        gene_var_tsv_ch
    )

    emit:
    fubar_tsv           = AGGREGATE_FIRST_PASS.out.fubar_tsv
    fubar_sites         = AGGREGATE_FIRST_PASS.out.fubar_sites.ifEmpty(file('NO_FUBAR_SITES'))
    fel_tsv             = AGGREGATE_FIRST_PASS.out.fel_tsv.ifEmpty(file('NO_FEL'))
    meme_tsv            = AGGREGATE_FINAL.out.meme_tsv.ifEmpty(file('NO_MEME'))
    sel_variability_tsv = AGGREGATE_FINAL.out.sel_variability.ifEmpty(file('NO_SEL_VAR'))
}
