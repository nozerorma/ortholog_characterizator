#!/usr/bin/env bash
# Batch threshold-based sequence filter.
#
# Drops any sequence whose non-gap, non-N bases are < min_seq_frac of the
# alignment length.  In other words:
#
#   ambiguity_fraction = (gaps + N) / alignment_length
#   keep if  (1 - ambiguity_fraction) >= min_seq_frac
#         i.e.  ambiguity_fraction     <= (1 - min_seq_frac)
#
# Examples:
#   min_seq_frac=0.1  → drop sequences with >90% gaps/N  (very permissive)
#   min_seq_frac=0.2  → drop sequences with >80% gaps/N
#   min_seq_frac=0.5  → drop sequences with >50% gaps/N
#   min_seq_frac=0.9  → drop sequences with >10% gaps/N  (very strict)

set -uo pipefail

manifest=""
workers=8
min_seq_frac=0.2
min_seqs=4
scripts_dir=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --manifest)     manifest="$2";     shift 2 ;;
        --workers)      workers="$2";       shift 2 ;;
        --min-seq-frac) min_seq_frac="$2";  shift 2 ;;
        --min-seqs)     min_seqs="$2";      shift 2 ;;
        --scripts-dir)  scripts_dir="$2";   shift 2 ;;
        *) echo "[THRESH] Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -z "$manifest" ]] && { echo "[THRESH] --manifest is required" >&2; exit 1; }

PYTHON="${PYTHON:-python3}"

wait_for_slot() {
    while [[ "$(jobs -pr | wc -l | tr -d ' ')" -ge "$workers" ]]; do
        wait -n 2>/dev/null || true
    done
}

while IFS=$'\t' read -r gene_id fasta_name; do
    [[ -z "$gene_id" || "$gene_id" == \#* ]] && continue
    [[ ! -f "$fasta_name" ]] && { echo "[THRESH] WARN: $fasta_name not found" >&2; continue; }

    wait_for_slot

    (
        out_fa="${gene_id}.filtered.fa"

        # Inline Python: keep sequences where non_gap_bases/aln_len >= min_seq_frac
        # i.e. drop sequences whose ambiguity fraction > (1 - min_seq_frac)
        "$PYTHON" - "$fasta_name" "$out_fa" "$min_seq_frac" "$min_seqs" <<'PYEOF'
import sys, re

in_fa, out_fa, min_frac, min_seqs = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4])

seqs = []
header = None
buf = []
with open(in_fa) as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith('>'):
            if header is not None:
                seqs.append((header, ''.join(buf)))
            header = line[1:]
            buf = []
        else:
            buf.append(line.upper())
    if header is not None:
        seqs.append((header, ''.join(buf)))

if not seqs:
    sys.exit(0)

aln_len = len(seqs[0][1])
kept = []
for hdr, seq in seqs:
    non_gap = sum(1 for c in seq if c not in '-N')
    if aln_len > 0 and non_gap / aln_len >= min_frac:
        kept.append((hdr, seq))

if len(kept) < min_seqs:
    sys.exit(0)

with open(out_fa, 'w') as fh:
    for hdr, seq in kept:
        fh.write(f'>{hdr}\n{seq}\n')
PYEOF

        [[ ! -s "$out_fa" ]] && rm -f "$out_fa"
        echo "[THRESH] OK: $gene_id"
    ) &

done < "$manifest"

while [[ "$(jobs -pr | wc -l | tr -d ' ')" -gt 0 ]]; do
    wait -n 2>/dev/null || true
done

echo "[THRESH] Batch complete"
