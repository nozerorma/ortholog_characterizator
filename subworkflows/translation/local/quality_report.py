#!/usr/bin/env python3
"""
scripts/quality_report.py

Quality assessment report across pipeline stages for codon MSA processing.

Computes per-gene statistics at up to four stages:
  raw        — original codon alignment  (ALIGN_DIR)
  prefilter  — after outlier+chi2 filter (OUT_DIR/PREFILTER/), when -q was used
  trim       — BMGE-filtered codon alignment (OUT_DIR/TRIM/)
  prot       — BMGE amino-acid projection   (OUT_DIR/PROT/)

Metrics reported per stage:
  - n_seq, aln_len
  - mean/median/Q1/Q3 gap % (fraction of all positions that are -/N/X/?)
  - mean/median/Q3/IQR ambiguous-codon % (codon-level, matches outlier-test.py)
  - n_allgap  : sequences that are 100% gap/N
  - n_highgap : sequences with >50% gap/N
  - completeness_pct : % sequences with ≤50% gap/N

Additional derived metrics:
  - trim_cols_removed     : alignment columns removed by BMGE
  - trim_pct_cols_kept    : % columns retained after BMGE
  - prefilter_n_removed   : sequences removed by outlier+chi2 filters

Outputs:
  <prefix>_quality_per_gene.csv   — one row per gene
  <prefix>_quality_summary.csv    — overall distribution across all genes

Usage:
    python3 scripts/quality_report.py \\
        --align_dir <raw_codon_dir> \\
        --out_dir   <base_output_dir> \\
        --output    <path/prefix>    \\
        [--gene_list <file>]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

# ── Sequence-level statistics ─────────────────────────────────────────────────

_AMBIG = frozenset('-NnXx?')


def _gap_frac(seq: str) -> float:
    """Fraction of characters that are ambiguous/gap."""
    n = len(seq)
    return sum(1 for c in seq if c in _AMBIG) / n if n else 1.0


def _ambig_codon_prop(seq: str) -> float:
    """
    Fraction of triplet codons containing at least one ambiguous character.
    Matches the logic used in outlier-test.py and chisq-test.py.
    """
    seq = seq.upper()
    n_codons = len(seq) // 3
    if not n_codons:
        return 1.0
    bad = sum(
        1 for i in range(0, n_codons * 3, 3)
        if any(seq[i + j] in _AMBIG for j in range(3))
    )
    return bad / n_codons


def codon_stats(fasta_path) -> dict | None:
    """Statistics for a codon (nucleotide) FASTA alignment."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return None

    n_seq = len(records)
    aln_len = len(records[0].seq)

    gap_fracs = np.array([_gap_frac(str(r.seq)) for r in records])
    ambig_props = np.array([_ambig_codon_prop(str(r.seq)) for r in records])

    def pct(arr): return np.round(arr * 100, 2)

    return {
        'n_seq':                    n_seq,
        'aln_len':                  aln_len,
        # Gap distribution (all-character level)
        'mean_gap_pct':             float(pct(gap_fracs.mean())),
        'median_gap_pct':           float(pct(np.median(gap_fracs))),
        'q1_gap_pct':               float(pct(np.percentile(gap_fracs, 25))),
        'q3_gap_pct':               float(pct(np.percentile(gap_fracs, 75))),
        'n_allgap':                 int((gap_fracs == 1.0).sum()),
        'n_highgap':                int((gap_fracs > 0.5).sum()),
        'completeness_pct':         float(round((gap_fracs <= 0.5).mean() * 100, 2)),
        # Ambiguous-codon distribution (codon level)
        'mean_ambig_codon_pct':     float(pct(ambig_props.mean())),
        'median_ambig_codon_pct':   float(pct(np.median(ambig_props))),
        'q3_ambig_codon_pct':       float(pct(np.percentile(ambig_props, 75))),
        'iqr_ambig_codon_pct':      float(pct(
            np.percentile(ambig_props, 75) - np.percentile(ambig_props, 25)
        )),
    }


def prot_stats(fasta_path) -> dict | None:
    """Statistics for a protein (amino-acid) FASTA alignment."""
    _PROT_AMBIG = frozenset('-Xx?')
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return None

    n_seq = len(records)
    aln_len = len(records[0].seq)
    gap_fracs = np.array([
        sum(1 for c in str(r.seq) if c in _PROT_AMBIG) / len(str(r.seq))
        if len(str(r.seq)) else 1.0
        for r in records
    ])

    return {
        'n_seq':            n_seq,
        'aln_len':          aln_len,
        'mean_gap_pct':     float(round(gap_fracs.mean() * 100, 2)),
        'median_gap_pct':   float(round(np.median(gap_fracs) * 100, 2)),
        'q1_gap_pct':       float(round(np.percentile(gap_fracs, 25) * 100, 2)),
        'q3_gap_pct':       float(round(np.percentile(gap_fracs, 75) * 100, 2)),
        'n_allgap':         int((gap_fracs == 1.0).sum()),
        'n_highgap':        int((gap_fracs > 0.5).sum()),
        'completeness_pct': float(round((gap_fracs <= 0.5).mean() * 100, 2)),
    }


# ── Prefilter removal counts ──────────────────────────────────────────────────

def count_removed(removed_dir: Path, gene_stem: str) -> int:
    """
    Count rows across all removed-species report CSVs that match this gene.
    Reports are written by 1.outlier_test.sh and 1.chisq_test.sh under
    PREFILTER/TESTS/Removed-species/<gene_stem>-summary.txt.
    """
    total = 0
    for fname in removed_dir.glob(f"{gene_stem}*"):
        try:
            df = pd.read_csv(fname)
            total += len(df)
        except Exception:
            pass
    return total


# ── Per-gene aggregation ──────────────────────────────────────────────────────

def _prefix(d: dict | None, stage: str) -> dict:
    if not d:
        return {}
    return {f"{stage}_{k}": v for k, v in d.items()}


_STAGE_SUFFIXES = ('.prot', '.trim', '.dedup', '.filtered')


def _base_gene(name: str) -> str:
    """Canonical gene id from a stage FASTA filename, ignoring stage extensions
    (e.g. 'ENST….codonified.prot.fa' and 'ENST….codonified.fa' → same id)."""
    stem = name[:-3] if name.endswith('.fa') else name
    for suf in _STAGE_SUFFIXES:
        if stem.endswith(suf):
            return stem[:-len(suf)]
    return stem


def codon_cols_by_gene(dirpath: str | None) -> dict:
    """Map base gene id → codon-column count for a stage directory. Reads only
    the first record per file (aligned sequences share one length)."""
    p = Path(dirpath) if dirpath else None
    if not p or not p.is_dir():
        return {}
    out = {}
    for fa in sorted(p.glob('*.fa')):
        rec = next(SeqIO.parse(fa, 'fasta'), None)
        if rec is not None:
            out[_base_gene(fa.name)] = len(rec.seq) // 3
    return out


def process_gene(fa_path: Path, out_dir: Path, has_prefilter: bool) -> dict:
    """Collect statistics for one gene across all pipeline stages."""
    gene_stem = fa_path.stem      # e.g. ENST00000001008.6#FKBP4.codonified
    gene_name = fa_path.name      # full filename with extension
    row = {'gene': gene_stem}

    # ── Stage 1: raw ─────────────────────────────────────────────────────────
    raw = codon_stats(fa_path)
    row.update(_prefix(raw, 'raw'))

    # ── Stage 2: prefilter (optional) ────────────────────────────────────────
    if has_prefilter:
        pre_path = out_dir / 'PREFILTER' / gene_name
        pre = codon_stats(pre_path) if pre_path.exists() else None
        row.update(_prefix(pre, 'prefilter'))
        removed_dir = out_dir / 'PREFILTER' / 'TESTS' / 'Removed-species'
        row['prefilter_n_removed'] = (
            count_removed(removed_dir, gene_stem) if removed_dir.exists() else None
        )
        if raw and pre:
            row['prefilter_pct_seqs_kept'] = round(pre['n_seq'] / raw['n_seq'] * 100, 2)

    # ── Stage 3: TRIM (post-BMGE codon) ──────────────────────────────────────
    trim_path = out_dir / 'TRIM' / gene_name
    trim = codon_stats(trim_path) if trim_path.exists() else None
    row.update(_prefix(trim, 'trim'))
    if raw and trim:
        row['trim_cols_removed']   = raw['aln_len'] - trim['aln_len']
        row['trim_pct_cols_kept']  = round(trim['aln_len'] / raw['aln_len'] * 100, 2)

    # ── Stage 4: PROT (post-BMGE protein) ────────────────────────────────────
    prot_path = out_dir / 'PROT' / gene_name
    prot = prot_stats(prot_path) if prot_path.exists() else None
    row.update(_prefix(prot, 'prot'))

    return row


# ── Summary across all genes ──────────────────────────────────────────────────

def _dist(df: pd.DataFrame, col: str, label: str) -> dict:
    """Distribution stats for a numeric column."""
    s = df[col].dropna() if col in df.columns else pd.Series(dtype=float)
    if s.empty:
        return {}
    return {
        f'{label}_mean':   round(float(s.mean()),             2),
        f'{label}_median': round(float(s.median()),           2),
        f'{label}_std':    round(float(s.std()),              2),
        f'{label}_min':    round(float(s.min()),              2),
        f'{label}_max':    round(float(s.max()),              2),
        f'{label}_q25':    round(float(s.quantile(0.25)),     2),
        f'{label}_q75':    round(float(s.quantile(0.75)),     2),
    }


def build_summary(df: pd.DataFrame, has_prefilter: bool) -> dict:
    """Build the overall summary statistics across all genes."""
    def safe_int_sum(col):
        return int(df[col].dropna().astype(int).gt(0).sum()) if col in df.columns else 0

    summ = {
        'total_genes':              len(df),
        'genes_with_trim':          int(df['trim_n_seq'].notna().sum()) if 'trim_n_seq' in df.columns else 0,
        'genes_with_prot':          int(df['prot_n_seq'].notna().sum()) if 'prot_n_seq' in df.columns else 0,
        'genes_with_allgap_raw':    safe_int_sum('raw_n_allgap'),
        'genes_with_allgap_trim':   safe_int_sum('trim_n_allgap'),
    }

    # Distribution of key metrics across genes
    metrics = [
        ('raw_n_seq',                  'raw_n_seq'),
        ('raw_aln_len',                'raw_aln_len'),
        ('raw_mean_gap_pct',           'raw_gap_pct'),
        ('raw_completeness_pct',       'raw_completeness_pct'),
        ('raw_mean_ambig_codon_pct',   'raw_ambig_codon_pct'),
        ('trim_n_seq',                 'trim_n_seq'),
        ('trim_aln_len',               'trim_aln_len'),
        ('trim_mean_gap_pct',          'trim_gap_pct'),
        ('trim_completeness_pct',      'trim_completeness_pct'),
        ('trim_cols_removed',          'trim_cols_removed'),
        ('trim_pct_cols_kept',         'trim_pct_cols_kept'),
        ('prot_mean_gap_pct',          'prot_gap_pct'),
        ('prot_completeness_pct',      'prot_completeness_pct'),
    ]
    if has_prefilter:
        metrics += [
            ('prefilter_n_removed',       'prefilter_n_removed'),
            ('prefilter_pct_seqs_kept',   'prefilter_pct_seqs_kept'),
        ]

    for col, label in metrics:
        summ.update(_dist(df, col, label))

    return summ


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Quality report across pipeline stages for codon MSA processing."
    )
    parser.add_argument('--align_dir', required=True,
                        help='Directory of raw codon alignment .fa files')
    parser.add_argument('--out_dir', required=True,
                        help='Base output directory (contains TRIM/, PROT/, PREFILTER/)')
    parser.add_argument('--output', required=True,
                        help='Output prefix (e.g. reports/quality)')
    parser.add_argument('--gene_list', default=None,
                        help='Optional file: one gene/filename fragment per line')
    parser.add_argument('--trim_dir', default=None,
                        help='Explicit BMGE-trimmed codon dir (overrides out_dir/TRIM lookup '
                             'for the trim-column derivation)')
    parser.add_argument('--pretrim_dir', default=None,
                        help='BMGE input codon dir (PREFILTER, else DEDUP) used as the baseline '
                             'for trim_cols_removed / trim_pct_cols_kept')
    args = parser.parse_args()

    align_dir = Path(args.align_dir)
    out_dir   = Path(args.out_dir)
    has_prefilter = (out_dir / 'PREFILTER').exists()

    fa_files = sorted(align_dir.glob('*.fa'))
    if args.gene_list:
        genes = set(Path(args.gene_list).read_text().split())
        fa_files = [f for f in fa_files if any(g in f.name for g in genes)]

    if not fa_files:
        print(f"No .fa files found in {align_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Quality report: processing {len(fa_files)} alignments "
          f"(prefilter={'yes' if has_prefilter else 'no'}) ...", flush=True)

    rows = [process_gene(fa, out_dir, has_prefilter) for fa in fa_files]

    per_gene_df = pd.DataFrame(rows)

    # ── BMGE trim columns from explicit stage dirs ───────────────────────────
    # The per-stage subdirs are usually not staged alongside this process, so
    # the out_dir/TRIM lookup in process_gene yields nothing. When the published
    # TRIM and BMGE-input (PREFILTER/DEDUP) dirs are passed explicitly, derive
    # the codon-level trim columns here by matching base gene id.
    if args.trim_dir and args.pretrim_dir and not per_gene_df.empty:
        trim_len = codon_cols_by_gene(args.trim_dir)
        pre_len  = codon_cols_by_gene(args.pretrim_dir)
        if trim_len and pre_len:
            base = per_gene_df['gene'].map(_base_gene)
            trim_s = base.map(trim_len)
            pre_s  = base.map(pre_len)
            per_gene_df['trim_aln_len']      = trim_s
            per_gene_df['trim_cols_removed'] = (pre_s - trim_s)
            per_gene_df['trim_pct_cols_kept'] = (trim_s / pre_s * 100).round(2)
            n_ok = int(per_gene_df['trim_pct_cols_kept'].notna().sum())
            print(f"  BMGE trim columns derived for {n_ok}/{len(per_gene_df)} genes "
                  f"(trim={args.trim_dir}, baseline={args.pretrim_dir})", flush=True)
        else:
            print(f"  WARN: could not read trim_dir ({len(trim_len)} genes) / "
                  f"pretrim_dir ({len(pre_len)} genes) — trim columns skipped", file=sys.stderr)
    summ        = build_summary(per_gene_df, has_prefilter)
    summary_df  = pd.DataFrame([summ])

    out_prefix = Path(args.output)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    per_gene_path = Path(f"{args.output}_quality_per_gene.csv")
    summary_path  = Path(f"{args.output}_quality_summary.csv")

    per_gene_df.to_csv(per_gene_path, index=False)
    summary_df.to_csv(summary_path,  index=False)

    # ── Console summary ───────────────────────────────────────────────────────
    print(f"\n── Quality Report ──────────────────────────────────────────")
    print(f"  Per-gene report  : {per_gene_path}")
    print(f"  Summary report   : {summary_path}")
    print(f"\n  Total genes                   : {summ['total_genes']}")
    print(f"  Genes with TRIM output        : {summ['genes_with_trim']}")
    print(f"  Genes with PROT output        : {summ['genes_with_prot']}")
    print(f"  Genes with all-gap seqs (raw) : {summ['genes_with_allgap_raw']}")
    print(f"  Genes with all-gap seqs (trim): {summ['genes_with_allgap_trim']}")
    if 'raw_completeness_pct_mean' in summ:
        print(f"\n  Raw  — mean completeness      : {summ['raw_completeness_pct_mean']}%")
        print(f"         mean ambig-codon pct    : {summ.get('raw_ambig_codon_pct_mean', 'N/A')}%")
    if 'trim_completeness_pct_mean' in summ:
        print(f"  Trim — mean completeness      : {summ['trim_completeness_pct_mean']}%")
        print(f"         mean cols kept          : {summ.get('trim_pct_cols_kept_mean', 'N/A')}%")
    if has_prefilter and 'prefilter_n_removed_mean' in summ:
        print(f"  Pre-filter — mean seqs removed: {summ['prefilter_n_removed_mean']}")
    print(f"────────────────────────────────────────────────────────────\n")


if __name__ == '__main__':
    main()
