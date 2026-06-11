#!/usr/bin/env bash
# Batch IQ-TREE gene tree inference.
# Strips all-gap/all-N sequences before running IQ-TREE on each gene.

set -uo pipefail

manifest=""
workers=4
cpu_each=4
model="MFP"
mset=""
bootstrap=1000
scripts_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)    manifest="$2";    shift 2 ;;
        --workers)     workers="$2";     shift 2 ;;
        --cpu-each)    cpu_each="$2";    shift 2 ;;
        --model)       model="$2";       shift 2 ;;
        --mset)        mset="$2";        shift 2 ;;
        --bootstrap)   bootstrap="$2";   shift 2 ;;
        --scripts-dir) scripts_dir="$2"; shift 2 ;;
        *) echo "[IQTREE] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$manifest" ]] && { echo "[IQTREE] --manifest is required" >&2; exit 1; }

PYTHON="${PYTHON:-python3}"

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        wait -n 2>/dev/null || true
    done
}

strip_bad_seqs() {
    # Remove sequences that are all gaps or all N/X
    local in_fa="$1"
    local out_fa="$2"
    "$PYTHON" - "$in_fa" "$out_fa" <<'PYEOF'
import sys
in_f, out_f = sys.argv[1], sys.argv[2]
seqs = []
hdr, buf = None, []
with open(in_f) as fh:
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

bad = set('-NX')
kept = [(h, s) for h, s in seqs if any(c not in bad for c in s)]
with open(out_f, 'w') as fh:
    for h, s in kept:
        fh.write(f'>{h}\n{s}\n')
PYEOF
}

while IFS=$'\t' read -r gene_id fasta_name; do
    [[ -z "$gene_id" || "$gene_id" == \#* ]] && continue
    [[ ! -f "$fasta_name" ]] && { echo "[IQTREE] WARN: $fasta_name not found, skipping $gene_id" >&2; continue; }

    wait_for_slot

    (
        clean_fa="${gene_id}.clean.fa"
        strip_bad_seqs "$fasta_name" "$clean_fa"

        n_seqs=$(grep -c '^>' "$clean_fa" 2>/dev/null || echo 0)
        if [[ "$n_seqs" -lt 4 ]]; then
            echo "[IQTREE] WARN: $gene_id has only $n_seqs usable seqs, skipping" >&2
            rm -f "$clean_fa"
            exit 0
        fi

        mset_arg=""
        [[ -n "$mset" ]] && mset_arg="-mset $mset"

        if iqtree -s "$clean_fa" \
                -m "$model" \
                $mset_arg \
                -B "$bootstrap" \
                -bnni \
                -T AUTO -ntmax "$cpu_each" \
                --prefix "$gene_id" \
                --quiet 2>"${gene_id}.iqtree.log"; then
            # Rename .treefile to our naming convention
            [[ -f "${gene_id}.treefile" ]] || { echo "[IQTREE] WARN: no treefile for $gene_id" >&2; }
            echo "[IQTREE] OK: $gene_id"
        else
            echo "[IQTREE] WARN: IQ-TREE failed for $gene_id" >&2
            rm -f "${gene_id}.treefile"
        fi
        rm -f "$clean_fa"
    ) &

done < "$manifest"

while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
    wait -n 2>/dev/null || true
done

echo "[IQTREE] Batch complete"
