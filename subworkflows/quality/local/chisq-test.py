#!/usr/bin/env python3
"""
Iterative chi-square test to remove sequences with excess ambiguous codons.

Runs up to --max_iter rounds. Each round removes sequences whose ambiguous-codon
proportion is significantly higher than the rest of the alignment (one-sided chi2,
p < 0.05). Iteration stops early when no sequences are removed in a round.

After filtering, the alignment is routed:
  low_seq   — remaining sequences < --min_seq_fraction of the pre-filter count
  ambiguous — overall alignment ambiguity > --ambiguity_threshold
  normal    — proceeds to the next pipeline stage (BMGE)

Output location:
  normal    → -o  (as given)
  low_seq   → --low_seq_dir/<basename of -o>
  ambiguous → --ambiguous_dir/<basename of -o>
"""

import argparse
import os
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from scipy.stats import chi2_contingency


# ---------------------------------------------------------------------------
# Codon counting
# ---------------------------------------------------------------------------
AMBIGUOUS_CHARS = {"N", "-", "X"}


def count_codon_types(aln: MultipleSeqAlignment) -> pd.DataFrame:
    well, amb = [], []
    for record in aln:
        seq = str(record.seq)
        w = a = 0
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            if any(c in AMBIGUOUS_CHARS for c in codon):
                a += 1
            else:
                w += 1
        well.append(w)
        amb.append(a)
    return pd.DataFrame(
        {"well_aligned": well, "ambiguous": amb},
        index=[r.id for r in aln],
    )


def overall_ambiguity(aln: MultipleSeqAlignment) -> float:
    total = ambiguous = 0
    for record in aln:
        seq = str(record.seq)
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            total += 1
            if any(c in AMBIGUOUS_CHARS for c in codon):
                ambiguous += 1
    return ambiguous / total if total > 0 else 0.0


# ---------------------------------------------------------------------------
# Chi-square test
# ---------------------------------------------------------------------------
def chi2test(df: pd.DataFrame) -> pd.DataFrame:
    """One-sided leave-one-out chi2. Sequences below background get p=1.0."""
    df = df.copy()
    df["total"] = df["well_aligned"] + df["ambiguous"]
    df["ambiguous_proportion"] = df["ambiguous"] / df["total"]

    total_well = df["well_aligned"].sum()
    total_amb = df["ambiguous"].sum()
    p_values = []

    for _, row in df.iterrows():
        w_i, a_i = row["well_aligned"], row["ambiguous"]
        n_i = w_i + a_i
        w_rest = total_well - w_i
        a_rest = total_amb - a_i
        n_rest = w_rest + a_rest

        p_i = a_i / n_i
        p_rest = a_rest / n_rest if n_rest > 0 else np.nan

        if np.isnan(p_rest) or p_i <= p_rest:
            p_values.append(1.0)
            continue

        table = np.array([[w_i, a_i], [w_rest, a_rest]])
        _, p, _, _ = chi2_contingency(table, correction=False)
        p_values.append(p)

    return pd.DataFrame(
        {"Ambiguous_Proportion": df["ambiguous_proportion"], "p-value": p_values},
        index=df.index,
    )


def filter_one_round(aln: MultipleSeqAlignment, threshold: float = 0.05):
    """
    Returns filtered alignment and a list of removed-sequence dicts.
    Only removes sequences with significantly MORE ambiguous codons than the rest.
    """
    counts_df = count_codon_types(aln)
    results = chi2test(counts_df)

    kept, removed = [], []
    for record in aln:
        p = results.at[record.id, "p-value"]
        if isinstance(p, pd.Series):
            p = p.iloc[0]
        if p < threshold:
            removed.append({
                "species": record.id,
                "Ambiguous_Proportion": results.at[record.id, "Ambiguous_Proportion"],
                "p-value": p,
            })
        else:
            kept.append(record)

    return MultipleSeqAlignment(kept), removed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Iterative chi-square filtering of codon alignments."
    )
    parser.add_argument("-i",  "--input",    required=True,  help="Input alignment (.fa)")
    parser.add_argument("-o",  "--output",   required=True,  help="Output alignment (.fa) for normal fate")
    parser.add_argument("-ot", "--output_test",    required=True,
                        help="TSV with chi2 results from each iteration")
    parser.add_argument("-or", "--output_removed", required=True,
                        help="TSV with all removed sequences across iterations")
    parser.add_argument("--max_iter",            type=int,   default=3,
                        help="Maximum number of chi2 iterations (default: 3)")
    parser.add_argument("--ambiguity_threshold", type=float, default=0.8,
                        help="Overall alignment ambiguity above which output goes to ambiguous/ (default: 0.8)")
    parser.add_argument("--min_seq_fraction",    type=float, default=0.1,
                        help="If remaining sequences < this fraction of input, goes to low_seq/ (default: 0.1)")
    parser.add_argument("--ambiguous_dir", default=None,
                        help="Directory for alignments flagged as ambiguous")
    parser.add_argument("--low_seq_dir",   default=None,
                        help="Directory for alignments flagged as low-sequence-count")
    args = parser.parse_args()

    alignment = AlignIO.read(args.input, "fasta")
    original_n = len(alignment)
    current_aln = alignment
    all_removed: list[dict] = []
    all_test_rows: list[dict] = []

    for iteration in range(1, args.max_iter + 1):
        if len(current_aln) < 2:
            break

        counts_df = count_codon_types(current_aln)
        results_df = chi2test(counts_df)

        for seq_id in results_df.index:
            all_test_rows.append({
                "iteration": iteration,
                "species": seq_id,
                "Ambiguous_Proportion": results_df.at[seq_id, "Ambiguous_Proportion"],
                "p-value": results_df.at[seq_id, "p-value"],
            })

        filtered_aln, removed = filter_one_round(current_aln)
        all_removed.extend(removed)

        current_aln = filtered_aln
        if not removed:
            break

    # Determine fate
    remaining_n = len(current_aln)
    if remaining_n < args.min_seq_fraction * original_n:
        fate = "low_seq"
        out_dir = args.low_seq_dir
    elif overall_ambiguity(current_aln) > args.ambiguity_threshold:
        fate = "ambiguous"
        out_dir = args.ambiguous_dir
    else:
        fate = "normal"
        out_dir = None

    # Write filtered alignment to the appropriate location
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, os.path.basename(args.output))
    else:
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        out_path = args.output

    AlignIO.write(current_aln, out_path, "fasta")

    # Write chi2 test results (all iterations)
    os.makedirs(os.path.dirname(args.output_test) or ".", exist_ok=True)
    test_df = pd.DataFrame(all_test_rows,
                           columns=["iteration", "species", "Ambiguous_Proportion", "p-value"])
    test_df.to_csv(args.output_test, index=False)

    # Write removed sequences report
    if all_removed:
        os.makedirs(os.path.dirname(args.output_removed) or ".", exist_ok=True)
        removed_df = pd.DataFrame(all_removed)
        removed_df.to_csv(args.output_removed, index=False)
        print(f"Removed {len(all_removed)} sequences across iterations")

    print(f"Chi2 test results saved to {args.output_test}")
    print(f"Filtered alignment ({remaining_n}/{original_n} sequences) saved to {out_path}")
    print(f"FATE: {fate}")


if __name__ == "__main__":
    main()
