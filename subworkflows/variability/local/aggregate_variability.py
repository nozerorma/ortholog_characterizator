#!/usr/bin/env python3
"""
Aggregate per-gene variability outputs from PROT_VAR_RAW/ into summary tables.

Reads:
  <var_dir>/*.entropy.tsv        per-position table (columns: gene, position,
                                  t, r, g, C_trident, variability, n_seqs)
  <var_dir>/*.clade_entropy.tsv  per-clade table (columns: gene, clade_level,
                                  clade_name, n_seqs, mean_variability, max_variability)

Writes to <out_dir>/:
  gene_variability.tsv      one row per gene (includes gap signal)
  clade_variability.tsv     one row per gene×clade (full table)
  clade_outliers.tsv        gene×clade pairs with |z| > outlier_z_threshold
  order_summary.tsv         per-order aggregation across genes
  family_summary.tsv        per-family (taxonomic) aggregation across genes
  pfam_summary.tsv          per-Pfam-domain aggregation (requires --pfam_tsv)
  pfam_order_summary.tsv    Pfam × order breakdown
  pfam_family_summary.tsv   Pfam × family breakdown
"""

import argparse
import glob
import os
import statistics
from collections import defaultdict

HIGH_GAP_THRESHOLD = 0.3  # column g > this is considered "high gap"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--var_dir', required=True, help='Path to PROT_VAR_RAW/')
    parser.add_argument('--out_dir', required=True, help='Output directory for reports')
    parser.add_argument('--pfam_tsv', default=None,
                        help='gene→Pfam mapping TSV (columns: gene, pfam_id); '
                             'if provided, writes pfam_summary.tsv')
    parser.add_argument('--pfam_desc_tsv', default=None,
                        help='Pfam ID→description TSV (columns: pfam_id, pfam_description); '
                             'if provided, adds pfam_description column to pfam_summary.tsv')
    parser.add_argument('--outlier_z_threshold', type=float, default=1.5,
                        help='|z-score| cutoff for clade_outliers.tsv (default 1.5)')
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # ── Per-gene summary from *.entropy.tsv ──────────────────────────────────
    # Column layout (0-based): gene, position, t, r, g, C_trident, variability, n_seqs
    COL_G   = 4
    COL_VAR = 6
    COL_N   = 7

    gene_rows = []
    for f in sorted(glob.glob(os.path.join(args.var_dir, '*.entropy.tsv'))):
        gene = os.path.basename(f).replace('.entropy.tsv', '')
        n_seqs = 0
        n_positions = 0
        n_variable = 0
        n_high_gap = 0
        variabilities = []
        gap_fractions = []
        with open(f) as fh:
            next(fh)  # header
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 8:
                    continue
                var = float(parts[COL_VAR])
                g   = float(parts[COL_G])
                n_seqs = int(parts[COL_N])
                n_positions += 1
                variabilities.append(var)
                gap_fractions.append(g)
                if var > 0.0:
                    n_variable += 1
                if g > HIGH_GAP_THRESHOLD:
                    n_high_gap += 1

        if n_positions == 0:
            continue
        mean_var     = sum(variabilities) / n_positions
        max_var      = max(variabilities)
        pct_variable = 100.0 * n_variable / n_positions
        mean_gap     = sum(gap_fractions) / n_positions
        pct_high_gap = 100.0 * n_high_gap / n_positions
        gene_rows.append((gene, n_seqs, n_positions, n_variable,
                          round(pct_variable, 2), round(mean_var, 6), round(max_var, 6),
                          round(mean_gap, 6), round(pct_high_gap, 2)))

    gene_out = os.path.join(args.out_dir, 'gene_variability.tsv')
    with open(gene_out, 'w') as fh:
        fh.write('gene\tn_seqs\tn_positions\tn_variable\tpct_variable\t'
                 'mean_variability\tmax_variability\tmean_gap_fraction\tpct_high_gap\n')
        for row in gene_rows:
            fh.write('\t'.join(str(x) for x in row) + '\n')
    print(f'Written: {gene_out} ({len(gene_rows)} genes)')

    # Build lookup dicts from gene_rows
    gene_mean     = {row[0]: row[5] for row in gene_rows}   # gene → mean_variability
    gene_mean_gap = {row[0]: row[7] for row in gene_rows}   # gene → mean_gap_fraction

    # ── Full clade variability table from *.clade_entropy.tsv ───────────────
    clade_rows = []
    for f in sorted(glob.glob(os.path.join(args.var_dir, '*.clade_entropy.tsv'))):
        with open(f) as fh:
            next(fh)
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 6:
                    clade_rows.append(parts[:6])

    clade_out = os.path.join(args.out_dir, 'clade_variability.tsv')
    with open(clade_out, 'w') as fh:
        fh.write('gene\tclade_level\tclade_name\tn_seqs\tmean_variability\tmax_variability\n')
        for row in clade_rows:
            fh.write('\t'.join(row) + '\n')
    print(f'Written: {clade_out} ({len(clade_rows)} gene×clade rows)')

    # ── Clade-specific outliers ──────────────────────────────────────────────
    # For each gene, z-score its per-clade mean_variability against the gene's
    # distribution across all clades at that level.
    outlier_rows = []
    for level in ('order', 'family'):
        gene_clade_data = defaultdict(list)  # gene → [(clade_name, n_seqs, mean_v)]
        for row in clade_rows:
            if row[1] == level:
                gene_clade_data[row[0]].append((row[2], int(row[3]), float(row[4])))

        for gene, entries in sorted(gene_clade_data.items()):
            vals = [e[2] for e in entries]
            if len(vals) < 2:
                continue
            gene_mean_v = statistics.mean(vals)
            gene_sd_v   = statistics.stdev(vals)
            if gene_sd_v == 0:
                continue
            for clade_name, n_seqs, mean_v in entries:
                z = (mean_v - gene_mean_v) / gene_sd_v
                if abs(z) > args.outlier_z_threshold:
                    outlier_rows.append((gene, level, clade_name, n_seqs,
                                         round(mean_v, 6), round(gene_mean_v, 6),
                                         round(gene_sd_v, 6), round(z, 4)))

    outlier_out = os.path.join(args.out_dir, 'clade_outliers.tsv')
    with open(outlier_out, 'w') as fh:
        fh.write('gene\tclade_level\tclade_name\tn_seqs\tmean_variability\t'
                 'gene_mean\tgene_sd\tz_score\n')
        for row in outlier_rows:
            fh.write('\t'.join(str(x) for x in row) + '\n')
    print(f'Written: {outlier_out} ({len(outlier_rows)} outlier gene×clade pairs)')

    # ── Per-level summaries (order / family) across genes ───────────────────
    # col layout in clade_rows: gene, clade_level, clade_name, n_seqs, mean_var, max_var
    for level in ('order', 'family'):
        level_vars      = defaultdict(list)
        level_gap       = defaultdict(list)
        level_n_species = defaultdict(int)  # max n_seqs seen across genes for this clade
        for row in clade_rows:
            if row[1] == level:
                clade_name = row[2]
                level_vars[clade_name].append(float(row[4]))
                level_n_species[clade_name] = max(level_n_species[clade_name], int(row[3]))
                if row[0] in gene_mean_gap:
                    level_gap[clade_name].append(gene_mean_gap[row[0]])

        summary_rows = []
        for clade_name, vals in sorted(level_vars.items()):
            n_genes  = len(vals)
            n_species = level_n_species[clade_name]
            mean_v   = sum(vals) / n_genes
            median_v = statistics.median(vals)
            sd_v     = statistics.stdev(vals) if n_genes > 1 else 0.0
            gaps     = level_gap[clade_name]
            mean_gap_v = sum(gaps) / len(gaps) if gaps else 0.0
            summary_rows.append((clade_name, n_genes, n_species,
                                  round(mean_v, 6), round(median_v, 6), round(sd_v, 6),
                                  round(mean_gap_v, 6)))

        out_path = os.path.join(args.out_dir, f'{level}_summary.tsv')
        with open(out_path, 'w') as fh:
            fh.write(f'{level}\tn_genes\tn_species\tmean_variability\tmedian_variability\t'
                     f'sd_variability\tmean_gap_fraction\n')
            for row in summary_rows:
                fh.write('\t'.join(str(x) for x in row) + '\n')
        print(f'Written: {out_path} ({len(summary_rows)} {level}s)')

    # ── Pfam-domain summary across genes ────────────────────────────────────
    if args.pfam_tsv:
        # pfam_id → human-readable description
        pfam_desc = {}
        if args.pfam_desc_tsv:
            with open(args.pfam_desc_tsv) as fh:
                next(fh)
                for line in fh:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) >= 2:
                        pfam_desc[parts[0]] = parts[1]

        # gene → set of pfam_ids; also build pfam_total_genes (full TSV, not just run genes)
        gene_to_pfams   = defaultdict(set)
        pfam_total_genes = defaultdict(int)  # pfam_id → total distinct genes in TSV
        with open(args.pfam_tsv) as fh:
            next(fh)  # header
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2 and parts[1]:
                    gene_to_pfams[parts[0]].add(parts[1])
                    pfam_total_genes[parts[1]] += 1
        # pfam_total_genes currently counts rows; we want distinct genes per pfam
        # Rebuild correctly: pfam → set of genes in the full TSV
        pfam_all_genes = defaultdict(set)
        with open(args.pfam_tsv) as fh:
            next(fh)
            for line in fh:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2 and parts[1]:
                    pfam_all_genes[parts[1]].add(parts[0])
        pfam_total_genes = {pfam_id: len(genes) for pfam_id, genes in pfam_all_genes.items()}

        # pfam_id → list of per-gene mean variabilities (current run only)
        pfam_vars  = defaultdict(list)
        pfam_genes = defaultdict(list)
        for gene, mean_v in gene_mean.items():
            for pfam_id in gene_to_pfams.get(gene, []):
                pfam_vars[pfam_id].append(mean_v)
                pfam_genes[pfam_id].append(gene)

        pfam_rows = []
        for pfam_id, vals in sorted(pfam_vars.items()):
            n_genes       = len(vals)
            n_genes_total = pfam_total_genes.get(pfam_id, n_genes)
            pct_covered   = round(100.0 * n_genes / n_genes_total, 2)
            mean_v        = sum(vals) / n_genes
            median_v      = statistics.median(vals)
            sd_v          = statistics.stdev(vals) if n_genes > 1 else 0.0
            genes_str     = ','.join(sorted(pfam_genes[pfam_id]))
            desc          = pfam_desc.get(pfam_id, '')
            pfam_rows.append((pfam_id, desc, n_genes, n_genes_total, pct_covered,
                              round(mean_v, 6), round(median_v, 6), round(sd_v, 6),
                              genes_str))

        pfam_out = os.path.join(args.out_dir, 'pfam_summary.tsv')
        with open(pfam_out, 'w') as fh:
            fh.write('pfam_id\tpfam_description\tn_genes\tn_genes_total\tpct_covered\t'
                     'mean_variability\tmedian_variability\tsd_variability\tgenes\n')
            for row in pfam_rows:
                fh.write('\t'.join(str(x) for x in row) + '\n')
        print(f'Written: {pfam_out} ({len(pfam_rows)} Pfam domains)')

        # ── Pfam × clade breakdowns (order and family) ───────────────────────
        for level in ('order', 'family'):
            pfam_clade_vars = defaultdict(lambda: defaultdict(list))
            pfam_clade_n    = defaultdict(lambda: defaultdict(int))
            for row in clade_rows:
                if row[1] != level:
                    continue
                row_gene   = row[0]
                clade_name = row[2]
                n_seqs     = int(row[3])
                mean_v     = float(row[4])
                for pfam_id in gene_to_pfams.get(row_gene, []):
                    pfam_clade_vars[pfam_id][clade_name].append(mean_v)
                    pfam_clade_n[pfam_id][clade_name] = max(
                        pfam_clade_n[pfam_id][clade_name], n_seqs)

            pfam_clade_rows = []
            for pfam_id in sorted(pfam_clade_vars):
                desc          = pfam_desc.get(pfam_id, '')
                n_genes_total = pfam_total_genes.get(pfam_id, 0)
                for clade_name in sorted(pfam_clade_vars[pfam_id]):
                    vals        = pfam_clade_vars[pfam_id][clade_name]
                    n_genes     = len(vals)
                    n_species   = pfam_clade_n[pfam_id][clade_name]
                    pct_covered = round(100.0 * n_genes / n_genes_total, 2) if n_genes_total else 0.0
                    mean_v      = sum(vals) / n_genes
                    median_v    = statistics.median(vals)
                    sd_v        = statistics.stdev(vals) if n_genes > 1 else 0.0
                    pfam_clade_rows.append((pfam_id, desc, clade_name,
                                            n_genes, n_genes_total, pct_covered, n_species,
                                            round(mean_v, 6), round(median_v, 6), round(sd_v, 6)))

            pfam_clade_out = os.path.join(args.out_dir, f'pfam_{level}_summary.tsv')
            with open(pfam_clade_out, 'w') as fh:
                fh.write(f'pfam_id\tpfam_description\t{level}\t'
                         f'n_genes\tn_genes_total\tpct_covered\tn_species\t'
                         f'mean_variability\tmedian_variability\tsd_variability\n')
                for row in pfam_clade_rows:
                    fh.write('\t'.join(str(x) for x in row) + '\n')
            print(f'Written: {pfam_clade_out} ({len(pfam_clade_rows)} pfam×{level} rows)')


if __name__ == '__main__':
    main()
