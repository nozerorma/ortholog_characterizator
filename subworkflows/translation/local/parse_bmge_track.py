#!/usr/bin/env python3
"""
Parse a BMGE TRACK HTML file and produce a plain-text TSV position map.

Output columns (one row per original alignment column):
  ori_codon_col   1-based column in original (pre-BMGE) CDS alignment
  status          selected | removed
  trim_codon_col  1-based rank among selected columns (= TRIM alignment col); NA if removed
  prot_ali_col    same as trim_codon_col (protein MSA column); NA if removed
  hg38_aa_pos     1-based AA position in the Homo_sapiens reference protein; NA if gapped or absent
  hg38_nt_coord   chr:genomic_position of first nt of this codon; NA if GFF not provided

Usage:
  python parse_bmge_track.py \\
      --track_html TRACK/A1BG.html \\
      --cds_fasta  DEDUP/A1BG.fa \\
      --out_tsv    MAP/A1BG.tsv \\
      [--gff       /path/to/hg38.gff3] \\
      [--protein_id ENSP00000263100] \\
      [--ref_species Homo_sapiens]
"""

import argparse
import os
import sys


# ---------------------------------------------------------------------------
# Gene → protein-ID equivalence table
# ---------------------------------------------------------------------------

def load_gene_equiv(path):
    """
    Load a two-column TSV (gene_name<TAB>protein_id) into a dict.
    Lines starting with '#' and blank lines are skipped.
    """
    equiv = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                equiv[parts[0].strip()] = parts[1].strip()
    return equiv


# ---------------------------------------------------------------------------
# Parse BMGE TRACK HTML
# ---------------------------------------------------------------------------

def parse_track_html(html_path):
    """
    Extract selected and removed column indices from a BMGE HTML tracking file.
    Lines of interest contain 'selected:' or 'removed:' as tokens.
    Indices are 1-based; supports both individual integers and 'a-b' ranges.
    Returns (selected: list[int], removed: list[int]), both sorted ascending.
    """
    selected, removed = [], []

    def _expand(tokens, start_idx):
        result = []
        for tok in tokens[start_idx:]:
            tok = tok.strip()
            if not tok:
                continue
            if '-' in tok:
                parts = tok.split('-')
                try:
                    a, b = int(parts[0]), int(parts[1])
                    result.extend(range(a, b + 1))
                except (ValueError, IndexError):
                    pass
            else:
                try:
                    result.append(int(tok))
                except ValueError:
                    pass
        return result

    with open(html_path) as fh:
        for line in fh:
            # Strip HTML tags for tokenisation
            plain = ''
            in_tag = False
            for ch in line:
                if ch == '<':
                    in_tag = True
                elif ch == '>':
                    in_tag = False
                elif not in_tag:
                    plain += ch
            tokens = plain.split()
            for i, tok in enumerate(tokens):
                if tok == 'selected:':
                    selected.extend(_expand(tokens, i + 1))
                    break
                elif tok == 'removed:':
                    removed.extend(_expand(tokens, i + 1))
                    break

    return sorted(selected), sorted(removed)


# ---------------------------------------------------------------------------
# Reference sequence helpers
# ---------------------------------------------------------------------------

def get_ref_aligned_seq(fasta_path, ref_patterns=None):
    """
    Return the uppercase codon string (with gaps) for the reference species.
    Tries each pattern in ref_patterns as a case-insensitive substring of the header.
    Returns None if not found.
    """
    if ref_patterns is None:
        ref_patterns = ['Homo_sapiens', 'hg38']

    seqs = {}
    header, buf = None, []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    seqs[header] = ''.join(buf).upper()
                header = line[1:]
                buf = []
            else:
                buf.append(line.upper())
    if header is not None:
        seqs[header] = ''.join(buf).upper()

    for pattern in ref_patterns:
        pattern_lower = pattern.lower()
        for h, seq in seqs.items():
            if pattern_lower in h.lower():
                return seq
    return None


def codon_cols_to_aa_pos(ref_aligned_seq):
    """
    Map each 1-based codon column → 1-based AA position in the reference protein.
    Gapped codon (first nt == '-') → None.
    Returns dict[codon_col_1based, aa_pos_1based | None].
    """
    col_to_aa = {}
    n_cols = len(ref_aligned_seq) // 3
    aa_counter = 0
    for col in range(n_cols):
        codon_start = col * 3
        first_nt = ref_aligned_seq[codon_start]
        if first_nt == '-':
            col_to_aa[col + 1] = None  # gapped in reference
        else:
            aa_counter += 1
            col_to_aa[col + 1] = aa_counter
    return col_to_aa


# ---------------------------------------------------------------------------
# GFF-based genomic coordinate mapping
# ---------------------------------------------------------------------------

def load_cds_positions(gff_path, protein_id, with_strand=False):
    """
    Parse GFF3/GFF2 for CDS features whose attribute field contains protein_id.
    Build the strand-aware ordered list of genomic nt positions for the CDS.
    Returns (cds_positions: list[int], chrom: str).
    Logic mirrors Extract_genomic_coordinates.py function F.
    """
    starts, ends, chrom, strand = [], [], None, None

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            if fields[2] != 'CDS':
                continue
            if protein_id not in fields[8]:
                continue
            chrom = fields[0]
            starts.append(int(fields[3]))
            ends.append(int(fields[4]))
            strand = fields[6]

    if not starts:
        return ([], chrom or '', strand) if with_strand else ([], chrom or '')

    cds_positions = []
    if strand == '+':
        for i in range(len(starts)):
            for j in range(starts[i], ends[i] + 1):
                cds_positions.append(j)
    else:  # '-' or unknown → reverse
        for i in range(len(starts)):
            curr_end = ends[-(i + 1)]
            curr_start = starts[-(i + 1)]
            for j in range(curr_end, curr_start - 1, -1):
                cds_positions.append(j)

    return (cds_positions, chrom, strand) if with_strand else (cds_positions, chrom)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Convert BMGE TRACK HTML to a plain TSV position map.')
    parser.add_argument('--track_html', required=True,
                        help='BMGE TRACK HTML file (TRACK/<gene>.html)')
    parser.add_argument('--cds_fasta', required=True,
                        help='Codon MSA FASTA used as BMGE input (DEDUP/<gene>.fa)')
    parser.add_argument('--out_tsv', required=True,
                        help='Output TSV path (MAP/<gene>.tsv)')
    parser.add_argument('--gff', default=None,
                        help='GFF3 file for hg38 genomic coordinate mapping (optional)')
    parser.add_argument('--protein_id', default=None,
                        help='Protein/transcript ID for GFF lookup; overrides --gene_equiv_tsv and gene name')
    parser.add_argument('--gene_equiv_tsv', default=None,
                        help='Two-column TSV (gene_name<TAB>protein_id) for GFF lookup; '
                             'e.g. proteiID_gene_equivalences.txt')
    parser.add_argument('--ref_species', default='Homo_sapiens',
                        help='Reference species header pattern (default: Homo_sapiens)')
    args = parser.parse_args()

    gene = os.path.basename(args.out_tsv).removesuffix('.map.tsv').removesuffix('.tsv')

    # 1. Parse TRACK HTML
    selected, removed = parse_track_html(args.track_html)
    if not selected and not removed:
        print(f'ERROR {gene}: no selected/removed positions found in {args.track_html}',
              file=sys.stderr)
        sys.exit(1)

    selected_set = set(selected)
    # trim_rank: ori_col → 1-based rank among selected columns
    trim_rank = {col: rank for rank, col in enumerate(sorted(selected), 1)}
    all_cols = sorted(selected_set | set(removed))

    # Sanity check
    expected = list(range(all_cols[0], all_cols[-1] + 1))
    if all_cols != expected:
        missing = set(expected) - set(all_cols)
        print(f'WARN {gene}: non-contiguous columns in TRACK HTML; {len(missing)} gaps',
              file=sys.stderr)

    # 2. Reference sequence → AA position mapping
    ref_seq = get_ref_aligned_seq(args.cds_fasta, [args.ref_species, 'hg38'])
    col_to_aa = {}
    if ref_seq is None:
        print(f'WARN {gene}: reference species "{args.ref_species}" not found in {args.cds_fasta}; '
              f'hg38_aa_pos will be NA', file=sys.stderr)
    else:
        col_to_aa = codon_cols_to_aa_pos(ref_seq)

    # 3. GFF genomic coordinates
    cds_positions, chrom, strand = [], '', None
    protein_id = args.protein_id
    if not protein_id and args.gene_equiv_tsv:
        equiv=load_gene_equiv(args.gene_equiv_tsv)
        protein_id=equiv.get(gene) or equiv.get(gene.removesuffix('.Homo_sapiens'))
    if args.gff:
        if protein_id:
            cds_positions, chrom, strand = load_cds_positions(args.gff, protein_id, with_strand=True)
        if not cds_positions:
            print(f'WARN {gene}: protein_id "{protein_id}" not found in {args.gff}; '
                  f'hg38_nt_coord will be NA', file=sys.stderr)

    # 4. Write TSV
    os.makedirs(os.path.dirname(args.out_tsv) or '.', exist_ok=True)
    with open(args.out_tsv, 'w') as fh:
        fh.write('ori_codon_col\tstatus\ttrim_codon_col\tprot_ali_col\t'
                 'hg38_aa_pos\thg38_nt_coord\thg38_nt_coord_2\thg38_nt_coord_3\tstrand\thuman_protein_id\thg38_codon\n')
        for col in all_cols:
            status = 'selected' if col in selected_set else 'removed'
            if status == 'selected':
                rank = trim_rank[col]
                trim_col = str(rank)
                prot_col = str(rank)
            else:
                trim_col = 'NA'
                prot_col = 'NA'

            aa_pos = col_to_aa.get(col)
            hg38_aa = str(aa_pos) if aa_pos is not None else 'NA'

            hg38_nt = 'NA'
            nt2, nt3 = 'NA', 'NA'
            coding = ref_seq[(col-1)*3:col*3] if ref_seq is not None else ''
            if len(coding)!=3 or set(coding)-set('ACGTN'):coding='NA'
            if cds_positions and aa_pos is not None:
                nt_idx = (aa_pos - 1) * 3  # 0-based index of first nt of this codon
                if 0 <= nt_idx and nt_idx+2 < len(cds_positions):
                    hg38_nt, nt2, nt3 = [f'{chrom}:{p}' for p in cds_positions[nt_idx:nt_idx+3]]

            fh.write(f"{col}\t{status}\t{trim_col}\t{prot_col}\t{hg38_aa}\t{hg38_nt}\t{nt2}\t{nt3}\t{strand or 'NA'}\t{protein_id or 'NA'}\t{coding}\n")

    n_sel = len(selected)
    n_rem = len(removed)
    print(f'{gene}: {n_sel} selected, {n_rem} removed → {args.out_tsv}', flush=True)


if __name__ == '__main__':
    main()
