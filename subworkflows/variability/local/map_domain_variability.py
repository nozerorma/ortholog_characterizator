#!/usr/bin/env python3
"""
Map Pfam domain hits (from hmmscan) onto the trimmed alignment and its per-position Valdar variability.

The hmmscan query is the full-length human protein of the codon alignment that enters BMGE
(full_length_query.py), not the trimmed protein. A hit residue is projected to the trimmed alignment
through the MAP table: residue -> codon column of the unfiltered alignment (ori_codon_col) ->
prot_ali_col. Residues in columns that BMGE removed have no trimmed column; they are counted in
n_removed and do not enter alignment_positions or the variability statistics.

Run hmmscan with the Pfam gathering thresholds (--cut_ga), the curated criterion of Pfam:
  hmmscan --cut_ga --domtblout hmmscan.domtblout Pfam-A.hmm queries.fa
A domtblout without that option needs --i_evalue_threshold, applied to the independent E-value.

Inputs: --manifest  TSV, one gene per line: query name, entropy file, MAP file, CDS alignment file
                    (names relative to --var_dir, --map_dir and --cds_dir; the query name is the gene id
                    of the entropy tables: <gene>, or <gene>.prot in tables of older runs)
Output: domain_variability.tsv
  gene, pfam_id, target_name, description, clan_acc, clan_name, domain_instance, ali_start, ali_end,
  n_positions, n_analyzed, mean_variability, max_variability, mean_gap_fraction,
  alignment_positions, coordinate_space, n_residues, n_removed

Coordinate system:
  - ali_start/ali_end = hmmscan 1-based positions in the full-length human protein
  - n_residues = domain length in residues; n_removed = residues in columns removed by trimming
  - alignment_positions = exact 0-based trimmed MSA columns; coordinate_space = CAAS_0based
  - n_positions = number of those columns
  - entropy.tsv position column = 1-based trimmed MSA column (0-based + 1)
"""

import argparse
import csv
import math
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from full_length_query import full_query  # noqa: E402


# ---------------------------------------------------------------------------
# hmmscan domtblout parser
# ---------------------------------------------------------------------------

def scan_options(domtblout_path):
    """Option string hmmscan recorded in the domtblout comments ('' when absent)."""
    with open(domtblout_path) as fh:
        for line in fh:
            if line.startswith('# Option settings:'):
                return line.split(':', 1)[1].strip()
    return ''


def parse_domtblout(domtblout_path, i_evalue_threshold=None):
    """
    Parse hmmscan --domtblout output.
    Relevant columns (0-indexed, space-separated):
      0:  target_name  (HMM NAME field, e.g. 'Ig_2')
      1:  accession    (HMM ACC field, e.g. 'PF13895.13' → strip to 'PF13895')
      3:  query_name   (gene name)
      5:  qlen         (query length)
      11: c-Evalue, 12: i-Evalue
      17: ali_from     (1-based start in query sequence)
      18: ali_to       (1-based end in query sequence)
      22+: description of target (joined remainder after fixed fields)
    i_evalue_threshold, when given, keeps the hits whose independent E-value is at most that value.
    Returns list of dicts; domain_instance is per-(gene, pfam_id) counter.
    """
    hits = []
    domain_counter = defaultdict(int)

    with open(domtblout_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 23:
                continue
            try:
                i_evalue = float(fields[12])
                ali_start = int(fields[17])
                ali_end   = int(fields[18])
                query_length = int(fields[5])
            except ValueError:
                continue
            if i_evalue_threshold is not None and i_evalue > i_evalue_threshold:
                continue

            pfam_id = fields[1].split('.')[0]  # strip version: PF13895.13 → PF13895
            gene    = fields[3]
            domain_counter[(gene, pfam_id)] += 1
            hits.append({
                'pfam_id':         pfam_id,
                'target_name':     fields[0],
                'description':     ' '.join(fields[22:]),
                'gene':            gene,
                'ali_start':       ali_start,
                'ali_end':         ali_end,
                'i_evalue':        i_evalue,
                'query_length':    query_length,
                'domain_instance': domain_counter[(gene, pfam_id)],
            })

    return hits


def load_clans(clans_tsv_path):
    """
    Load Pfam-A.clans.tsv → dict: pfam_acc → (clan_acc, clan_name).
    File columns: pfam_acc, clan_acc, clan_name, pfam_name, pfam_description
    Domains not in any clan (orphans) will be absent from the dict.
    """
    mapping = {}
    with open(clans_tsv_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            pfam_acc, clan_acc, clan_name = parts[0], parts[1], parts[2]
            mapping[pfam_acc] = (clan_acc, clan_name)
    return mapping


# ---------------------------------------------------------------------------
# MAP table
# ---------------------------------------------------------------------------

def load_map(map_path):
    """
    Read <gene>.map.tsv (BMGE position map).
    Returns {ori_codon_col: trimmed 0-based column, or None when BMGE removed the column}.
    The retained columns must be exactly 0..n-1, the productive alignment.
    """
    out = {}
    with open(map_path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row['status'] == 'selected':
                out[int(row['ori_codon_col'])] = int(row['prot_ali_col']) - 1
            else:
                out[int(row['ori_codon_col'])] = None
    kept = sorted(c for c in out.values() if c is not None)
    if kept != list(range(len(kept))):
        raise ValueError(f'{map_path}: retained columns are not 1..n of the productive alignment')
    return out


# ---------------------------------------------------------------------------
# Entropy TSV loader
# ---------------------------------------------------------------------------

def load_entropy_tsv(tsv_path):
    """
    Read <gene>.entropy.tsv.
    Returns dict: {position_1based: {'t': float, 'r': float, 'g': float,
                                      'C_trident': float, 'variability': float}}
    position is the 1-based MSA column as written by compute_variability.py.
    """
    data = {}
    with open(tsv_path) as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 7:
                continue
            try:
                pos = int(parts[1])
                data[pos] = {
                    't': float(parts[2]),
                    'r': float(parts[3]),
                    'g': float(parts[4]),
                    'C_trident': float(parts[5]),
                    'variability': float(parts[6]),
                }
            except (ValueError, IndexError):
                continue
    return data


# ---------------------------------------------------------------------------
# Domain statistics
# ---------------------------------------------------------------------------

def compute_domain_stats(ali_cols_0based, entropy_data):
    """
    Collect entropy rows for the given 0-based alignment columns.
    entropy_data uses 1-based keys → convert: key = ali_col_0based + 1.
    Returns dict with n_positions, n_analyzed, mean_variability,
    max_variability, mean_gap_fraction.
    """
    variabilities = []
    gap_fracs = []
    for col0 in ali_cols_0based:
        pos1 = col0 + 1
        if pos1 in entropy_data:
            row = entropy_data[pos1]
            variabilities.append(row['variability'])
            gap_fracs.append(row['g'])

    n_positions = len(ali_cols_0based)
    n_analyzed = len(variabilities)
    if n_analyzed == 0:
        return {
            'n_positions': n_positions, 'n_analyzed': 0,
            'mean_variability': float('nan'), 'max_variability': float('nan'),
            'mean_gap_fraction': float('nan'),
        }
    return {
        'n_positions': n_positions,
        'n_analyzed': n_analyzed,
        'mean_variability': sum(variabilities) / n_analyzed,
        'max_variability': max(variabilities),
        'mean_gap_fraction': sum(gap_fracs) / n_analyzed,
    }


# ---------------------------------------------------------------------------
# Projection of one gene
# ---------------------------------------------------------------------------

def project_hit(hit, ori_cols, trimmed_col):
    """(retained trimmed columns, residues in removed columns) of a hit on the full-length protein."""
    if hit['query_length'] != len(ori_cols):
        raise ValueError(f'{hit["gene"]}: hmmscan query does not match the full-length human protein')
    if not 1 <= hit['ali_start'] <= hit['ali_end'] <= len(ori_cols):
        raise ValueError(f'{hit["gene"]}: hmmscan endpoints outside the full-length human protein')
    cols, removed = [], 0
    for residue in range(hit['ali_start'], hit['ali_end'] + 1):
        col = trimmed_col.get(ori_cols[residue - 1], 'missing')
        if col == 'missing':
            raise ValueError(f'{hit["gene"]}: codon column {ori_cols[residue - 1]} is not in the MAP table')
        if col is None:
            removed += 1
        else:
            cols.append(col)
    return cols, removed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Map Pfam domain hits of the full-length human protein onto the trimmed alignment.')
    parser.add_argument('--manifest', required=True,
                        help='TSV: query name, entropy file, MAP file, CDS alignment file')
    parser.add_argument('--var_dir', required=True,
                        help='Directory containing the entropy TSV files of the manifest')
    parser.add_argument('--map_dir', required=True,
                        help='Directory containing the MAP TSV files of the manifest')
    parser.add_argument('--cds_dir', required=True,
                        help='Directory containing the codon alignments of the manifest')
    parser.add_argument('--hmmscan_domtblout', required=True,
                        help='hmmscan --domtblout output file')
    parser.add_argument('--out_tsv', required=True,
                        help='Output TSV path (domain_variability.tsv)')
    parser.add_argument('--alignment_coordinate_space', required=True,
                        choices=['trimmed_alignment_1based'],
                        help='Coordinate space of the entropy tables: the productive trimmed alignment')
    parser.add_argument('--i_evalue_threshold', type=float, default=None,
                        help='Independent E-value cutoff; only for a domtblout produced without --cut_ga')
    parser.add_argument('--pfam_clans_tsv', default=None,
                        help='Pfam-A.clans.tsv: maps pfam_acc → clan_acc, clan_name')
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out_tsv) or '.', exist_ok=True)

    cut_ga = '--cut_ga' in scan_options(args.hmmscan_domtblout).split()
    if cut_ga and args.i_evalue_threshold is not None:
        raise ValueError('--i_evalue_threshold applies only to a domtblout produced without --cut_ga')
    if not cut_ga and args.i_evalue_threshold is None:
        raise ValueError('the domtblout does not record --cut_ga: pass --i_evalue_threshold to state the criterion')

    clan_map = {}
    if args.pfam_clans_tsv and os.path.isfile(args.pfam_clans_tsv):
        clan_map = load_clans(args.pfam_clans_tsv)
        print(f'Loaded {len(clan_map)} clan mappings from {args.pfam_clans_tsv}', flush=True)

    hits = parse_domtblout(args.hmmscan_domtblout, args.i_evalue_threshold)
    if not hits:
        print(f'WARNING: no domain hits found in {args.hmmscan_domtblout}', file=sys.stderr)

    files = {}
    with open(args.manifest) as fh:
        for line in fh:
            if line.strip():
                name, entropy, map_file, cds = line.rstrip('\n').split('\t')[:4]
                files[name] = (entropy, map_file, cds)

    hits_by_gene = defaultdict(list)
    for h in hits:
        hits_by_gene[h['gene']].append(h)

    rows = []
    for gene, gene_hits in sorted(hits_by_gene.items()):
        if gene not in files:
            print(f'WARN {gene}: not in the manifest, hits skipped', file=sys.stderr)
            continue
        entropy_file, map_file, cds_file = files[gene]

        query = full_query(os.path.join(args.cds_dir, cds_file))
        if query is None:
            print(f'WARN {gene}: no human record in {cds_file}', file=sys.stderr)
            continue
        _, ori_cols, n_cds_cols = query

        trimmed_col = load_map(os.path.join(args.map_dir, map_file))
        if n_cds_cols != len(trimmed_col):
            raise ValueError(f'{gene}: MAP table and codon alignment have different columns')
        n_kept = sum(c is not None for c in trimmed_col.values())
        entropy_data = load_entropy_tsv(os.path.join(args.var_dir, entropy_file))
        if sorted(entropy_data) != list(range(1, n_kept + 1)):
            raise ValueError(f'{gene}: entropy table and MAP table have different columns')

        for hit in gene_hits:
            ali_cols, n_removed = project_hit(hit, ori_cols, trimmed_col)
            if not ali_cols:
                print(f'WARN {gene} {hit["pfam_id"]}: domain range {hit["ali_start"]}-{hit["ali_end"]} '
                      'lies entirely in columns removed by trimming', file=sys.stderr)
                continue

            stats = compute_domain_stats(ali_cols, entropy_data)
            mean_v = stats['mean_variability']
            max_v  = stats['max_variability']
            mean_g = stats['mean_gap_fraction']

            clan_acc, clan_name = clan_map.get(hit['pfam_id'], ('', ''))

            rows.append((
                gene,
                hit['pfam_id'],
                hit['target_name'],
                hit['description'],
                clan_acc,
                clan_name,
                hit['domain_instance'],
                hit['ali_start'], hit['ali_end'],
                stats['n_positions'],
                stats['n_analyzed'],
                f'{mean_v:.6f}' if not math.isnan(mean_v) else 'NA',
                f'{max_v:.6f}'  if not math.isnan(max_v)  else 'NA',
                f'{mean_g:.6f}' if not math.isnan(mean_g) else 'NA',
                ','.join(map(str, ali_cols)),
                'CAAS_0based',
                hit['ali_end'] - hit['ali_start'] + 1,
                n_removed,
            ))

    with open(args.out_tsv, 'w') as fh:
        fh.write('gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\t'
                 'domain_instance\tali_start\tali_end\t'
                 'n_positions\tn_analyzed\tmean_variability\tmax_variability\t'
                 'mean_gap_fraction\talignment_positions\tcoordinate_space\t'
                 'n_residues\tn_removed\n')
        for row in rows:
            fh.write('\t'.join(str(x) for x in row) + '\n')

    print(f'Written: {args.out_tsv} ({len(rows)} domain hits)', flush=True)


if __name__ == '__main__':
    main()
