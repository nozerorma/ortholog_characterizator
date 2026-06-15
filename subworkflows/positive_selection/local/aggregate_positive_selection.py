#!/usr/bin/env python3
"""
Aggregate HyPhy positive selection results (FUBAR + FEL + MEME) across genes.

Reads PSEL/<gene>/<gene>.{FUBAR,FEL,MEME,SLAC}.json files, extracts key
statistics, and optionally joins with gene_variability.tsv.

Outputs (written to --out_dir):
  fubar_results.tsv   One row per gene: n_pos_sites (posterior>=threshold),
                      n_pos_sites_fdr (empirical Bayesian FDR), max/mean_posterior
  fubar_sites.tsv     One row per gene×site; adds is_pos_hit_fdr (Bayesian FDR call)
  fel_results.tsv     One row per gene×site with dN/dS and p-values
  busted_results.tsv  One row per gene: busted_pval, busted_qval (BH-FDR), is_busted_sig
  meme_results.tsv    One row per gene×site; adds testable (dS>0 & beta+ in-bound),
                      meme_qval (BH-FDR over testable sites), is_meme_sig
  meme_candidate_genes.txt  Gate for the MEME second pass (with --write_candidates):
                      BUSTED-significant genes, or the FUBAR-FDR fallback
  sel_variability.tsv Join of fubar_results × gene_variability (if --fubar_tsv or FUBAR in --analyses,
                      and --variability_tsv provided)

Multiple-testing model: FUBAR is Bayesian (empirical Bayesian FDR, not BH);
BUSTED and MEME are frequentist (Benjamini-Hochberg). MEME sites with dS==0 or
beta+ at the optimiser boundary are excluded before FDR as degenerate.

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
# BUSTED results loader (for sel_variability join)
# ---------------------------------------------------------------------------

def load_busted_lookup(busted_tsv):
    """Return {gene: (busted_qval, is_busted_sig)} from a busted_results.tsv.
    Empty dict if the file is absent/empty (so the join degrades gracefully)."""
    if not busted_tsv or not os.path.isfile(busted_tsv) or os.path.getsize(busted_tsv) == 0:
        return {}
    out = {}
    with open(busted_tsv) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        gi = hdr.index('gene') if 'gene' in hdr else 0
        qi = hdr.index('busted_qval') if 'busted_qval' in hdr else None
        si = hdr.index('is_busted_sig') if 'is_busted_sig' in hdr else None
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) <= gi:
                continue
            q = p[qi] if (qi is not None and qi < len(p)) else 'NA'
            s = p[si] if (si is not None and si < len(p)) else 'NA'
            out[p[gi]] = (q, s)
    return out


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
# BUSTED parser (gene-level episodic selection test)
# ---------------------------------------------------------------------------

def parse_busted(json_path):
    """Return the gene-level LRT p-value for episodic diversifying selection
    from a BUSTED JSON, or None on failure."""
    try:
        with open(json_path) as fh:
            d = json.load(fh)
    except Exception as e:
        print(f'WARN {json_path}: {e}', file=sys.stderr)
        return None
    try:
        # HyPhy BUSTED stores the LRT p-value under d['test results']['p-value']
        return float(d['test results']['p-value'])
    except (KeyError, TypeError, ValueError) as e:
        print(f'WARN {json_path}: BUSTED parse error — {e}', file=sys.stderr)
        return None


# ---------------------------------------------------------------------------
# Multiple-testing helpers
# ---------------------------------------------------------------------------

def bh_fdr(pvals):
    """Benjamini-Hochberg q-values for *pvals*, returned aligned to input order.

    Frequentist FDR — use for p-value methods (BUSTED gene tests, MEME sites).
    """
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [1.0] * n
    prev = 1.0
    # walk ranks high→low so the step-up monotonicity holds
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * n / rank)
        q[i] = prev
    return q


def bayesian_fdr_cutoff(posteriors, target_fdr):
    """Empirical Bayesian FDR cutoff for FUBAR-style posteriors.

    FUBAR is Bayesian and has no p-values, so BH is inappropriate. For a set of
    called sites the expected FDR is mean(1 - posterior). Calling the highest-
    posterior sites first, that running mean is non-decreasing, so the largest
    set with realised FDR <= target is a prefix. Returns the posterior cutoff t;
    sites with posterior >= t are called. Returns 1.0 (no calls) if even the top
    site fails the target.
    """
    if not posteriors:
        return 1.0
    s = sorted(posteriors, reverse=True)
    cum = 0.0
    cutoff = 1.0
    for k, p in enumerate(s, 1):
        cum += (1.0 - p)
        if cum / k <= target_fdr:
            cutoff = p
        else:
            break
    return cutoff


# MEME per-site degeneracy guards: dS collapsing to 0 makes beta+ > dS trivially
# true (saturation artifact), and beta+ pinned at the optimiser bound is
# unidentifiable. Both are dropped before FDR so they cannot manufacture hits.
MEME_BETA_BOUNDARY = 1e5


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
    parser.add_argument('--fubar_fdr', type=float, default=0.05,
                        help='Target empirical Bayesian FDR for FUBAR site calls (default: 0.05)')
    parser.add_argument('--meme_fdr', type=float, default=0.05,
                        help='Target BH-FDR for MEME site calls, after dS>0 / boundary filtering (default: 0.05)')
    parser.add_argument('--busted_fdr', type=float, default=0.05,
                        help='Target BH-FDR for the gene-level BUSTED gate (default: 0.05)')
    parser.add_argument('--meme_min_sites', type=int, default=1,
                        help='Fallback MEME gate when BUSTED is absent: min FUBAR-FDR sites per gene (default: 1)')
    parser.add_argument('--write_candidates', action='store_true',
                        help='Write meme_candidate_genes.txt (BUSTED-FDR genes, or FUBAR fallback)')
    parser.add_argument('--fubar_tsv', default=None,
                        help='Pre-computed fubar_results.tsv to use for sel_variability join '
                             '(skips re-parsing FUBAR JSONs; used in AGGREGATE_FINAL)')
    parser.add_argument('--variability_tsv', default=None,
                        help='gene_variability.tsv to join with FUBAR gene-level results')
    parser.add_argument('--busted_tsv', default=None,
                        help='busted_results.tsv to fold into the sel_variability join '
                             '(adds busted_qval, is_busted_sig per gene)')
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
        # Buffer site rows so the genome-wide Bayesian-FDR cutoff can be computed
        # before writing the per-site is_pos_hit_fdr column. At ~16k genes this is
        # a few GB on the aggregation node; acceptable for a once-per-run step.
        buffered = []          # (gene, site, dS, dN, prob_pos, prob_neg, is_pos_hit, is_neg_hit, map_cols)
        all_pos_post = []      # genome-wide prob_pos for the FDR cutoff
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
                buffered.append((gene, s['site'], s['dS'], s['dN'], s['prob_pos'],
                                 s['prob_neg'], int(s['is_pos_hit']), int(s['is_neg_hit']),
                                 map_cols))
                all_pos_post.append(s['prob_pos'])

        # Genome-wide empirical Bayesian FDR cutoff (Bayesian, not BH).
        fdr_cutoff = bayesian_fdr_cutoff(all_pos_post, args.fubar_fdr)
        n_pos_fdr_per_gene = {}

        fubar_site_path = os.path.join(args.out_dir, 'fubar_sites.tsv')
        map_cols_header = '\tori_codon_col\thg38_aa_pos\thg38_nt_coord' if args.map_dir else ''
        with open(fubar_site_path, 'w') as sfh:
            sfh.write(f'gene\tsite\tdS\tdN\tprob_pos\tprob_neg\tis_pos_hit\tis_pos_hit_fdr\tis_neg_hit{map_cols_header}\n')
            for (gene, site, dS, dN, pp, pn, is_pos, is_neg, map_cols) in buffered:
                is_pos_fdr = int(pp >= fdr_cutoff)
                if is_pos_fdr:
                    n_pos_fdr_per_gene[gene] = n_pos_fdr_per_gene.get(gene, 0) + 1
                sfh.write(f"{gene}\t{site}\t{dS}\t{dN}\t{pp}\t{pn}\t"
                          f"{is_pos}\t{is_pos_fdr}\t{is_neg}{map_cols}\n")

        fubar_path = os.path.join(args.out_dir, 'fubar_results.tsv')
        with open(fubar_path, 'w') as fh:
            fh.write('gene\tn_pos_sites\tn_pos_sites_fdr\tn_neg_sites\tmax_posterior\tmean_posterior\n')
            for r in sorted(fubar_rows, key=lambda x: -x['n_pos_sites']):
                n_fdr = n_pos_fdr_per_gene.get(r['gene'], 0)
                r['n_pos_sites_fdr'] = n_fdr
                fh.write(f"{r['gene']}\t{r['n_pos_sites']}\t{n_fdr}\t{r['n_neg_sites']}\t"
                         f"{r['max_posterior']}\t{r['mean_posterior']}\n")
        n_genes_fdr = sum(1 for r in fubar_rows if r['n_pos_sites_fdr'] > 0)
        print(f'Written: {fubar_path} ({len(fubar_rows)} genes)', flush=True)
        print(f'Written: {fubar_site_path} ({len(buffered)} gene×site rows)', flush=True)
        print(f'FUBAR Bayesian FDR {args.fubar_fdr}: posterior cutoff {fdr_cutoff:.4f}, '
              f'{sum(n_pos_fdr_per_gene.values())} sites in {n_genes_fdr} genes', flush=True)

        # ── Join with variability ────────────────────────────────────────────
        if args.variability_tsv and os.path.isfile(args.variability_tsv):
            var_data = {}
            with open(args.variability_tsv) as fh:
                header = fh.readline().rstrip().split('\t')
                for line in fh:
                    parts = line.rstrip().split('\t')
                    if parts:
                        key = parts[0].removesuffix('.prot')
                        var_data[key] = dict(zip(header[1:], parts[1:]))

            busted_lookup = load_busted_lookup(args.busted_tsv)
            sv_path = os.path.join(args.out_dir, 'sel_variability.tsv')
            var_cols = ['n_seqs', 'n_positions', 'mean_variability',
                        'mean_gap_fraction', 'pct_high_gap']
            with open(sv_path, 'w') as fh:
                fh.write('gene\tn_pos_sites\tn_pos_sites_fdr\tn_neg_sites\tmax_posterior\tmean_posterior\t'
                         'busted_qval\tis_busted_sig\t' + '\t'.join(var_cols) + '\n')
                for r in sorted(fubar_rows, key=lambda x: -x['n_pos_sites']):
                    vr = var_data.get(r['gene'], {})
                    vcols = '\t'.join(vr.get(c, 'NA') for c in var_cols)
                    bq, bs = busted_lookup.get(r['gene'], ('NA', 'NA'))
                    fh.write(f"{r['gene']}\t{r['n_pos_sites']}\t{r.get('n_pos_sites_fdr', 'NA')}\t"
                             f"{r['n_neg_sites']}\t{r['max_posterior']}\t{r['mean_posterior']}\t"
                             f"{bq}\t{bs}\t{vcols}\n")
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
        # Buffer all sites: BH-FDR needs the full p-value set, and we only test
        # sites that survive the dS>0 / beta+ boundary guards (saturation and
        # unidentifiable estimates would otherwise manufacture significant calls).
        meme_rows = []   # dicts with gene, site, dS, beta+, p+, p_meme, map_cols, testable
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
                row['gene'] = gene
                row['map_cols'] = map_cols
                row['testable'] = (row['dS'] > 0.0 and row['beta+'] < MEME_BETA_BOUNDARY)
                meme_rows.append(row)

        # BH-FDR across the testable sites only.
        testable = [r for r in meme_rows if r['testable']]
        qvals = bh_fdr([r['p_meme'] for r in testable])
        for r, q in zip(testable, qvals):
            r['meme_qval'] = round(q, 6)
        for r in meme_rows:
            if not r['testable']:
                r['meme_qval'] = 'NA'
            r['is_meme_sig'] = int(r['testable'] and r['meme_qval'] != 'NA'
                                   and r['meme_qval'] <= args.meme_fdr)

        meme_path = os.path.join(args.out_dir, 'meme_results.tsv')
        map_cols_header = '\tori_codon_col\thg38_aa_pos\thg38_nt_coord' if args.map_dir else ''
        with open(meme_path, 'w') as fh:
            fh.write(f'gene\tsite\tdS\tbeta+\tp+\tp_meme\ttestable\tmeme_qval\tis_meme_sig{map_cols_header}\n')
            for r in meme_rows:
                fh.write(f"{r['gene']}\t{r['site']}\t{r['dS']}\t{r['beta+']}\t"
                         f"{r['p+']}\t{r['p_meme']}\t{int(r['testable'])}\t"
                         f"{r['meme_qval']}\t{r['is_meme_sig']}{r['map_cols']}\n")
        n_sig = sum(r['is_meme_sig'] for r in meme_rows)
        n_drop = sum(1 for r in meme_rows if not r['testable'])
        print(f'Written: {meme_path} ({len(meme_rows)} gene×site rows; '
              f'{len(testable)} testable, {n_drop} dropped by guards; '
              f'{n_sig} significant at BH-FDR {args.meme_fdr})', flush=True)

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

    # ── BUSTED (gene-level gate) ───────────────────────────────────────────────
    busted_sig_genes = None
    if 'BUSTED' in analyses:
        busted = []   # (gene, pval)
        for gene in genes:
            jp = os.path.join(args.psel_dir, gene, f'{gene}.BUSTED.json')
            if not os.path.isfile(jp):
                continue
            pv = parse_busted(jp)
            if pv is not None:
                busted.append((gene, pv))

        qvals = bh_fdr([p for _, p in busted])
        busted_path = os.path.join(args.out_dir, 'busted_results.tsv')
        busted_sig_genes = []
        with open(busted_path, 'w') as fh:
            fh.write('gene\tbusted_pval\tbusted_qval\tis_busted_sig\n')
            for (gene, pv), q in sorted(zip(busted, qvals), key=lambda x: x[1]):
                sig = int(q <= args.busted_fdr)
                if sig:
                    busted_sig_genes.append(gene)
                fh.write(f'{gene}\t{pv:.6g}\t{q:.6f}\t{sig}\n')
        print(f'Written: {busted_path} ({len(busted)} genes; '
              f'{len(busted_sig_genes)} significant at BH-FDR {args.busted_fdr})', flush=True)

    # ── MEME candidate genes (gate for the second pass) ────────────────────────
    if args.write_candidates:
        cand_path = os.path.join(args.out_dir, 'meme_candidate_genes.txt')
        if busted_sig_genes is not None:
            candidates = busted_sig_genes
            gate = f'BUSTED BH-FDR <= {args.busted_fdr}'
        else:
            # Fallback when BUSTED was not run: genes with >= meme_min_sites
            # FUBAR-FDR positive sites (n_pos_sites_fdr from the FUBAR block).
            candidates = [r['gene'] for r in fubar_rows
                          if r.get('n_pos_sites_fdr', 0) >= args.meme_min_sites] \
                         if 'FUBAR' in analyses else []
            gate = f'FUBAR n_pos_sites_fdr >= {args.meme_min_sites}'
        with open(cand_path, 'w') as fh:
            if candidates:
                fh.write('\n'.join(str(g) for g in candidates) + '\n')
        print(f'Written: {cand_path} ({len(candidates)} candidates; gate: {gate})', flush=True)

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
                    key = parts[0].removesuffix('.prot')
                    var_data[key] = dict(zip(header[1:], parts[1:]))

        fubar_rows = []
        with open(args.fubar_tsv) as fh:
            hdr = fh.readline().rstrip().split('\t')
            for line in fh:
                parts = line.rstrip().split('\t')
                if parts:
                    fubar_rows.append(dict(zip(hdr, parts)))

        busted_lookup = load_busted_lookup(args.busted_tsv)
        sv_path = os.path.join(args.out_dir, 'sel_variability.tsv')
        var_cols = ['n_seqs', 'n_positions', 'mean_variability',
                    'mean_gap_fraction', 'pct_high_gap']
        with open(sv_path, 'w') as fh:
            fh.write('gene\tn_pos_sites\tn_pos_sites_fdr\tn_neg_sites\tmax_posterior\tmean_posterior\t'
                     'busted_qval\tis_busted_sig\t' + '\t'.join(var_cols) + '\n')
            for r in sorted(fubar_rows, key=lambda x: -int(x.get('n_pos_sites', 0))):
                vr = var_data.get(r['gene'], {})
                vcols = '\t'.join(vr.get(c, 'NA') for c in var_cols)
                bq, bs = busted_lookup.get(r['gene'], ('NA', 'NA'))
                fh.write(f"{r['gene']}\t{r['n_pos_sites']}\t{r.get('n_pos_sites_fdr', 'NA')}\t"
                         f"{r['n_neg_sites']}\t{r['max_posterior']}\t{r['mean_posterior']}\t"
                         f"{bq}\t{bs}\t{vcols}\n")
        print(f'Written: {sv_path}', flush=True)


if __name__ == '__main__':
    main()
