#!/usr/bin/env python3
"""
Aggregate HyPhy positive selection results (FUBAR + FEL + MEME) across genes.

Reads PSEL/<gene>/<gene>.{FUBAR,FEL,MEME,SLAC}.json files, extracts key
statistics, and optionally joins with gene_variability.tsv.

Outputs (written to --out_dir):
  fubar_results.tsv   One row per gene: n_pos_sites, max_posterior, mean_posterior
  fubar_sites.tsv     One row per gene×site with posteriors (all sites)
  fel_results.tsv     One row per gene×site with dN/dS and p-values
  meme_results.tsv    One row per gene×site with p-value and beta+
  sel_variability.tsv Join of fubar_results × gene_variability (if --fubar_tsv or FUBAR in --analyses,
                      and --variability_tsv provided)

Usage:
  python aggregate_positive_selection.py \
      --psel_dir OUT_DIR/PSEL \
      --analyses FUBAR,FEL \
      --out_dir  OUT_DIR/reports/positive_selection \
      [--fubar_threshold 0.9] \
      [--fubar_tsv OUT_DIR/reports/positive_selection/fubar_results.tsv] \
      [--variability_tsv OUT_DIR/reports/variability/gene_variability.tsv]
"""

import argparse
import json
import os
import sys


# ---------------------------------------------------------------------------
# BMGE position-map loader
# ---------------------------------------------------------------------------

def load_gene_map(map_dir, gene):
    """
    Load a BMGE map TSV for *gene* from *map_dir*.
    Returns dict: trim_codon_col(int) -> (ori_codon_col, hg38_aa_pos, hg38_nt_coord).
    Keys are 1-based trim alignment positions matching HyPhy site numbering.
    Returns empty dict if the file is absent.
    """
    path = os.path.join(map_dir, f'{gene}.map.tsv')
    if not os.path.isfile(path):
        return {}
    result = {}
    with open(path) as fh:
        fh.readline()  # skip header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            ori_col, status, trim_col, _prot_col, hg38_aa, hg38_nt = parts[:6]
            if status == 'selected' and trim_col != 'NA':
                try:
                    result[int(trim_col)] = (ori_col, hg38_aa, hg38_nt)
                except ValueError:
                    pass
    return result


# ---------------------------------------------------------------------------
# FUBAR parser
# ---------------------------------------------------------------------------

def _headers(d):
    """Normalise MLE headers: handles both ['name', ...] and [['name', 'desc'], ...]."""
    h = d['MLE']['headers']
    return [x[0] if isinstance(x, list) else x for x in h]


def parse_fubar(json_path, threshold=0.9):
    """
    Extract per-site posteriors from a FUBAR JSON.
    Returns (gene_summary_dict, list_of_site_dicts).
    gene_summary_dict keys: n_pos_sites, max_posterior, mean_posterior
    """
    try:
        with open(json_path) as fh:
            d = json.load(fh)
    except Exception as e:
        print(f'WARN {json_path}: {e}', file=sys.stderr)
        return None, []

    try:
        headers = _headers(d)
        content = d['MLE']['content']['0']

        # Locate columns (names vary slightly across HyPhy versions)
        def _idx(candidates):
            for c in candidates:
                if c in headers:
                    return headers.index(c)
            raise ValueError(f'None of {candidates} found in {headers}')

        idx_alpha   = _idx(['alpha'])
        idx_beta    = _idx(['beta'])
        idx_pos     = _idx(['Prob[alpha<beta]', 'P[alpha<beta]'])
        idx_neg     = _idx(['Prob[alpha>beta]', 'P[alpha>beta]'])

        site_rows = []
        pos_posteriors = []
        neg_posteriors = []
        for i, row in enumerate(content, 1):
            alpha    = float(row[idx_alpha])
            beta     = float(row[idx_beta])
            prob_pos = float(row[idx_pos])
            prob_neg = float(row[idx_neg])
            pos_posteriors.append(prob_pos)
            neg_posteriors.append(prob_neg)
            site_rows.append({
                'site':       i,
                'dS':         round(alpha,    4),
                'dN':         round(beta,     4),
                'prob_pos':   round(prob_pos, 6),
                'prob_neg':   round(prob_neg, 6),
                'is_pos_hit': prob_pos >= threshold,
                'is_neg_hit': prob_neg >= threshold,
            })

        n_pos = sum(1 for p in pos_posteriors if p >= threshold)
        n_neg = sum(1 for p in neg_posteriors if p >= threshold)
        summary = {
            'n_pos_sites':    n_pos,
            'n_neg_sites':    n_neg,
            'max_posterior':  round(max(pos_posteriors), 6) if pos_posteriors else 0.0,
            'mean_posterior': round(sum(pos_posteriors) / len(pos_posteriors), 6) if pos_posteriors else 0.0,
        }
        return summary, site_rows

    except (KeyError, ValueError, TypeError, IndexError) as e:
        print(f'WARN {json_path}: FUBAR parse error — {e}', file=sys.stderr)
        return None, []


# ---------------------------------------------------------------------------
# FEL parser
# ---------------------------------------------------------------------------

def parse_fel(json_path):
    """Extract per-site dS, dN, p-value from a FEL JSON."""
    try:
        with open(json_path) as fh:
            d = json.load(fh)
    except Exception as e:
        print(f'WARN {json_path}: {e}', file=sys.stderr)
        return []

    rows = []
    try:
        headers = _headers(d)
        content = d['MLE']['content']['0']

        try:
            idx_alpha = headers.index('alpha')
            idx_beta  = headers.index('beta')
            idx_pval  = headers.index('p-value')
        except ValueError:
            idx_alpha, idx_beta, idx_pval = 0, 1, 4

        for i, site in enumerate(content, 1):
            rows.append({
                'site':  i,
                'dS':    round(float(site[idx_alpha]), 4),
                'dN':    round(float(site[idx_beta]),  4),
                'p_fel': round(float(site[idx_pval]),  6),
            })
    except (KeyError, TypeError, IndexError) as e:
        print(f'WARN {json_path}: FEL parse error — {e}', file=sys.stderr)

    return rows


# ---------------------------------------------------------------------------
# MEME parser
# ---------------------------------------------------------------------------

def parse_meme(json_path):
    """
    Extract per-site episodic positive selection results from a MEME JSON.
    Key columns: alpha (dS), beta+ (dN on selected branches), p+ (proportion),
    p-value (site-level LRT p-value).
    """
    try:
        with open(json_path) as fh:
            d = json.load(fh)
    except Exception as e:
        print(f'WARN {json_path}: {e}', file=sys.stderr)
        return []

    rows = []
    try:
        headers = _headers(d)
        content = d['MLE']['content']['0']

        def _idx(candidates):
            for c in candidates:
                if c in headers:
                    return headers.index(c)
            raise ValueError(f'None of {candidates} in {headers}')

        # HyPhy ≥2.5 uses HTML entities in MLE headers (e.g. &alpha;, &beta;<sup>+</sup>)
        idx_alpha = _idx(['alpha', '&alpha;'])
        idx_betap = _idx(['beta+', '&beta;<sup>+</sup>'])
        idx_pp    = _idx(['p+', 'p<sup>+</sup>', 'fraction of branches under selection'])
        idx_pval  = _idx(['p-value'])

        for i, site in enumerate(content, 1):
            rows.append({
                'site':   i,
                'dS':     round(float(site[idx_alpha]), 4),
                'beta+':  round(float(site[idx_betap]), 4),
                'p+':     round(float(site[idx_pp]),    4),
                'p_meme': round(float(site[idx_pval]),  6),
            })
    except (KeyError, ValueError, TypeError, IndexError) as e:
        print(f'WARN {json_path}: MEME parse error — {e}', file=sys.stderr)

    return rows


# ---------------------------------------------------------------------------
# SLAC parser
# ---------------------------------------------------------------------------

def parse_slac(json_path):
    """Extract per-site dS, dN, p+ from a SLAC JSON."""
    try:
        with open(json_path) as fh:
            d = json.load(fh)
    except Exception as e:
        print(f'WARN {json_path}: {e}', file=sys.stderr)
        return []

    rows = []
    try:
        by_site = d['MLE']['content']['0']
        for i, site in enumerate(by_site, 1):
            # SLAC columns: ES, EN, S, N, dS, dN, dS-dN, p+, p-
            rows.append({
                'site':    i,
                'dS':      round(float(site[4]), 4),
                'dN':      round(float(site[5]), 4),
                'p_slac+': round(float(site[7]), 6),
            })
    except (KeyError, TypeError, IndexError) as e:
        print(f'WARN {json_path}: SLAC parse error — {e}', file=sys.stderr)

    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Aggregate HyPhy positive selection results.')
    parser.add_argument('--psel_dir', required=True,
                        help='Directory containing PSEL/<gene>/ subdirectories')
    parser.add_argument('--analyses', default='FUBAR,FEL',
                        help='Comma-separated analyses to aggregate (default: FUBAR,FEL)')
    parser.add_argument('--out_dir', required=True, help='Output directory for TSVs')
    parser.add_argument('--fubar_threshold', type=float, default=0.9,
                        help='Posterior probability threshold for FUBAR positive site (default: 0.9)')
    parser.add_argument('--fubar_tsv', default=None,
                        help='Pre-computed fubar_results.tsv to use for sel_variability join '
                             '(skips re-parsing FUBAR JSONs; used in AGGREGATE_FINAL)')
    parser.add_argument('--variability_tsv', default=None,
                        help='gene_variability.tsv to join with FUBAR gene-level results')
    parser.add_argument('--map_dir', default=None,
                        help='Directory of *.map.tsv BMGE position maps; when provided, '
                             'site-level outputs gain ori_codon_col, hg38_aa_pos, hg38_nt_coord columns')
    args = parser.parse_args()

    analyses = [a.strip().upper() for a in args.analyses.split(',')]
    os.makedirs(args.out_dir, exist_ok=True)

    genes = sorted(
        d for d in os.listdir(args.psel_dir)
        if os.path.isdir(os.path.join(args.psel_dir, d))
    )

    # ── FUBAR ─────────────────────────────────────────────────────────────────
    if 'FUBAR' in analyses:
        fubar_rows = []
        fubar_site_path = os.path.join(args.out_dir, 'fubar_sites.tsv')
        n_site_rows = 0

        map_cols_header = '\tori_codon_col\thg38_aa_pos\thg38_nt_coord' if args.map_dir else ''
        with open(fubar_site_path, 'w') as sfh:
            sfh.write(f'gene\tsite\tdS\tdN\tprob_pos\tprob_neg\tis_pos_hit\tis_neg_hit{map_cols_header}\n')
            for gene in genes:
                jp = os.path.join(args.psel_dir, gene, f'{gene}.FUBAR.json')
                if not os.path.isfile(jp):
                    continue
                summary, sites = parse_fubar(jp, threshold=args.fubar_threshold)
                if summary is not None:
                    fubar_rows.append({'gene': gene, **summary})
                gene_map = load_gene_map(args.map_dir, gene) if args.map_dir else {}
                for s in sites:
                    map_cols = ''
                    if args.map_dir:
                        ori_col, hg38_aa, hg38_nt = gene_map.get(s['site'], ('NA', 'NA', 'NA'))
                        map_cols = f'\t{ori_col}\t{hg38_aa}\t{hg38_nt}'
                    sfh.write(f"{gene}\t{s['site']}\t{s['dS']}\t{s['dN']}\t"
                              f"{s['prob_pos']}\t{s['prob_neg']}\t"
                              f"{int(s['is_pos_hit'])}\t{int(s['is_neg_hit'])}{map_cols}\n")
                    n_site_rows += 1

        fubar_path = os.path.join(args.out_dir, 'fubar_results.tsv')
        with open(fubar_path, 'w') as fh:
            fh.write('gene\tn_pos_sites\tn_neg_sites\tmax_posterior\tmean_posterior\n')
            for r in sorted(fubar_rows, key=lambda x: -x['n_pos_sites']):
                fh.write(f"{r['gene']}\t{r['n_pos_sites']}\t{r['n_neg_sites']}\t"
                         f"{r['max_posterior']}\t{r['mean_posterior']}\n")
        print(f'Written: {fubar_path} ({len(fubar_rows)} genes)', flush=True)
        print(f'Written: {fubar_site_path} ({n_site_rows} gene×site rows)', flush=True)

        # ── Join with variability ────────────────────────────────────────────
        if args.variability_tsv and os.path.isfile(args.variability_tsv):
            var_data = {}
            with open(args.variability_tsv) as fh:
                header = fh.readline().rstrip().split('\t')
                for line in fh:
                    parts = line.rstrip().split('\t')
                    if parts:
                        var_data[parts[0]] = dict(zip(header[1:], parts[1:]))

            sv_path = os.path.join(args.out_dir, 'sel_variability.tsv')
            var_cols = ['n_seqs', 'n_positions', 'mean_variability',
                        'mean_gap_fraction', 'pct_high_gap']
            with open(sv_path, 'w') as fh:
                fh.write('gene\tn_pos_sites\tn_neg_sites\tmax_posterior\tmean_posterior\t' +
                         '\t'.join(var_cols) + '\n')
                for r in sorted(fubar_rows, key=lambda x: -x['n_pos_sites']):
                    vr = var_data.get(r['gene'], {})
                    vcols = '\t'.join(vr.get(c, 'NA') for c in var_cols)
                    fh.write(f"{r['gene']}\t{r['n_pos_sites']}\t{r['n_neg_sites']}\t"
                             f"{r['max_posterior']}\t{r['mean_posterior']}\t{vcols}\n")
            print(f'Written: {sv_path}', flush=True)

    # ── FEL ───────────────────────────────────────────────────────────────────
    if 'FEL' in analyses:
        fel_path = os.path.join(args.out_dir, 'fel_results.tsv')
        n_rows = 0
        map_cols_header = '\tori_codon_col\thg38_aa_pos\thg38_nt_coord' if args.map_dir else ''
        with open(fel_path, 'w') as fh:
            fh.write(f'gene\tsite\tdS\tdN\tp_fel{map_cols_header}\n')
            for gene in genes:
                jp = os.path.join(args.psel_dir, gene, f'{gene}.FEL.json')
                if not os.path.isfile(jp):
                    continue
                gene_map = load_gene_map(args.map_dir, gene) if args.map_dir else {}
                for row in parse_fel(jp):
                    map_cols = ''
                    if args.map_dir:
                        ori_col, hg38_aa, hg38_nt = gene_map.get(row['site'], ('NA', 'NA', 'NA'))
                        map_cols = f'\t{ori_col}\t{hg38_aa}\t{hg38_nt}'
                    fh.write(f"{gene}\t{row['site']}\t{row['dS']}\t{row['dN']}\t{row['p_fel']}{map_cols}\n")
                    n_rows += 1
        print(f'Written: {fel_path} ({n_rows} gene×site rows)', flush=True)

    # ── MEME ──────────────────────────────────────────────────────────────────
    if 'MEME' in analyses:
        meme_path = os.path.join(args.out_dir, 'meme_results.tsv')
        n_rows = 0
        map_cols_header = '\tori_codon_col\thg38_aa_pos\thg38_nt_coord' if args.map_dir else ''
        with open(meme_path, 'w') as fh:
            fh.write(f'gene\tsite\tdS\tbeta+\tp+\tp_meme{map_cols_header}\n')
            for gene in genes:
                jp = os.path.join(args.psel_dir, gene, f'{gene}.MEME.json')
                if not os.path.isfile(jp):
                    continue
                gene_map = load_gene_map(args.map_dir, gene) if args.map_dir else {}
                for row in parse_meme(jp):
                    map_cols = ''
                    if args.map_dir:
                        ori_col, hg38_aa, hg38_nt = gene_map.get(row['site'], ('NA', 'NA', 'NA'))
                        map_cols = f'\t{ori_col}\t{hg38_aa}\t{hg38_nt}'
                    fh.write(f"{gene}\t{row['site']}\t{row['dS']}\t{row['beta+']}\t"
                             f"{row['p+']}\t{row['p_meme']}{map_cols}\n")
                    n_rows += 1
        print(f'Written: {meme_path} ({n_rows} gene×site rows)', flush=True)

    # ── SLAC ──────────────────────────────────────────────────────────────────
    if 'SLAC' in analyses:
        slac_path = os.path.join(args.out_dir, 'slac_results.tsv')
        n_rows = 0
        map_cols_header = '\tori_codon_col\thg38_aa_pos\thg38_nt_coord' if args.map_dir else ''
        with open(slac_path, 'w') as fh:
            fh.write(f'gene\tsite\tdS\tdN\tp_slac+{map_cols_header}\n')
            for gene in genes:
                jp = os.path.join(args.psel_dir, gene, f'{gene}.SLAC.json')
                if not os.path.isfile(jp):
                    continue
                gene_map = load_gene_map(args.map_dir, gene) if args.map_dir else {}
                for row in parse_slac(jp):
                    map_cols = ''
                    if args.map_dir:
                        ori_col, hg38_aa, hg38_nt = gene_map.get(row['site'], ('NA', 'NA', 'NA'))
                        map_cols = f'\t{ori_col}\t{hg38_aa}\t{hg38_nt}'
                    fh.write(f"{gene}\t{row['site']}\t{row['dS']}\t{row['dN']}\t{row['p_slac+']}{map_cols}\n")
                    n_rows += 1
        print(f'Written: {slac_path} ({n_rows} gene×site rows)', flush=True)

    # ── sel_variability from pre-computed FUBAR tsv ───────────────────────────
    # Used in AGGREGATE_FINAL where FUBAR JSONs are not available but
    # fubar_results.tsv from AGGREGATE_FIRST_PASS is passed in via --fubar_tsv.
    if ('FUBAR' not in analyses
            and args.fubar_tsv and os.path.isfile(args.fubar_tsv)
            and os.path.getsize(args.fubar_tsv) > 0
            and args.variability_tsv and os.path.isfile(args.variability_tsv)):
        var_data = {}
        with open(args.variability_tsv) as fh:
            header = fh.readline().rstrip().split('\t')
            for line in fh:
                parts = line.rstrip().split('\t')
                if parts:
                    var_data[parts[0]] = dict(zip(header[1:], parts[1:]))

        fubar_rows = []
        with open(args.fubar_tsv) as fh:
            hdr = fh.readline().rstrip().split('\t')
            for line in fh:
                parts = line.rstrip().split('\t')
                if parts:
                    fubar_rows.append(dict(zip(hdr, parts)))

        sv_path = os.path.join(args.out_dir, 'sel_variability.tsv')
        var_cols = ['n_seqs', 'n_positions', 'mean_variability',
                    'mean_gap_fraction', 'pct_high_gap']
        with open(sv_path, 'w') as fh:
            fh.write('gene\tn_pos_sites\tn_neg_sites\tmax_posterior\tmean_posterior\t' +
                     '\t'.join(var_cols) + '\n')
            for r in sorted(fubar_rows, key=lambda x: -int(x.get('n_pos_sites', 0))):
                vr = var_data.get(r['gene'], {})
                vcols = '\t'.join(vr.get(c, 'NA') for c in var_cols)
                fh.write(f"{r['gene']}\t{r['n_pos_sites']}\t{r['n_neg_sites']}\t"
                         f"{r['max_posterior']}\t{r['mean_posterior']}\t{vcols}\n")
        print(f'Written: {sv_path}', flush=True)


if __name__ == '__main__':
    main()
