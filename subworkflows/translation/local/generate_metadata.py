import os
import pandas as pd
import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from collections import Counter
import argparse

def process_directory(dir_path):
    # Initialize variables
    species_counter = Counter()
    file_species_counts = []

    # Iterate over all .fa files in the directory
    for msa in os.listdir(dir_path):
        # Exclude files in subdirectories
        if os.path.isdir(f"{dir_path}/{msa}"):
            continue
        # If the file is a .fa file and is not empty
        if msa.endswith(".fa") and os.path.getsize(f"{dir_path}/{msa}") > 0:
            alignment = AlignIO.read(f"{dir_path}/{msa}", "fasta")
            species_set = set()

            for record in alignment:
                species = record.id
                species_set.add(species)
                species_counter[species] += 1

            file_species_counts.append(len(species_set))

    return species_counter, file_species_counts

def combine_portions(portion_dirs):
    combined_species_counter = Counter()
    combined_species_counts = []
    total_proteins = 0  # Initialize outside the loop

    for portion_dir in portion_dirs:
        # Accumulate total proteins across all directories
        total_proteins += len([msa for msa in os.listdir(portion_dir) if msa.endswith(".fa")])
        species_counter, file_species_counts = process_directory(portion_dir)
        combined_species_counter.update(species_counter)
        combined_species_counts.extend(file_species_counts)

    species_counts = np.array(combined_species_counts)

    # Add debugging output
    print(f"Species counts: {species_counts}")
    print(f"Total proteins: {total_proteins}")

    # Handle empty species_counts
    if total_proteins == 0:
        print("No data to process. Returning empty results.")
        metadata = {
            'Total Proteins': 0,
            'Sp/Protein Global Mean': np.nan,
            'Sp/Protein Global Median': np.nan,
            'Sp/Protein Global SD': np.nan,
            'Sp/Protein Quantile 25': np.nan,
            'Sp/Protein Quantile 50': np.nan,
            'Sp/Protein Quantile 75': np.nan,
            'MSA Above Mean n.Sp': 0,
            '% MSA Above Mean n.Sp': 0,
            'MSA Above Median n.Sp': 0,
            '% MSA Above Median n.Sp': 0,
            'MSA Above Q25 n.Sp': 0,
            '% MSA Above Q25 n.Sp': 0,
            'MSA Above Q50 n.Sp': 0,
            '% MSA Above Q50 n.Sp': 0,
            'MSA Above Q75 n.Sp': 0,
            '% MSA Above Q75 n.Sp': 0
        }
        species_df = pd.DataFrame(columns=['Occurrences'])
        return metadata, species_df

    # Calculate statistics on combined data
    global_mean = np.mean(species_counts)
    global_median = np.median(species_counts)
    global_sd = np.std(species_counts)
    quantile_25 = np.percentile(species_counts, 25)
    quantile_50 = np.percentile(species_counts, 50)
    quantile_75 = np.percentile(species_counts, 75)

    # Count genes surpassing thresholds
    genes_above_mean = np.sum(species_counts > global_mean)
    genes_above_median = np.sum(species_counts > global_median)
    genes_above_q25 = np.sum(species_counts > quantile_25)
    genes_above_q50 = np.sum(species_counts > quantile_50)
    genes_above_q75 = np.sum(species_counts > quantile_75)

    # Create a DataFrame for species counts
    species_df = pd.DataFrame.from_dict(combined_species_counter, orient='index', columns=['Occurrences'])

    # Generate the metadata report
    metadata = {
        'Total Proteins': total_proteins,
        'Sp/Protein Global Mean': global_mean,
        'Sp/Protein Global Median': global_median,
        'Sp/Protein Global SD': global_sd,
        'Sp/Protein Quantile 25': quantile_25,
        'Sp/Protein Quantile 50': quantile_50,
        'Sp/Protein Quantile 75': quantile_75,
        'MSA Above Mean n.Sp': genes_above_mean,
        '% MSA Above Mean n.Sp': genes_above_mean / total_proteins * 100,
        'MSA Above Median n.Sp': genes_above_median,
        '% MSA Above Median n.Sp': genes_above_median / total_proteins * 100,
        'MSA Above Q25 n.Sp': genes_above_q25,
        '% MSA Above Q25 n.Sp': genes_above_q25 / total_proteins * 100,
        'MSA Above Q50 n.Sp': genes_above_q50,
        '% MSA Above Q50 n.Sp': genes_above_q50 / total_proteins * 100,
        'MSA Above Q75 n.Sp': genes_above_q75,
        '% MSA Above Q75 n.Sp': genes_above_q75 / total_proteins * 100
    }

    return metadata, species_df

# Lets copy the MSA alignments which have a number of species greater than the global mean

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine portions of metadata and species data.")
    parser.add_argument("portion_dirs", nargs='+', help="List of directories containing portions to combine.")
    parser.add_argument("output_prefix", help="Prefix for the combined output files.")

    args = parser.parse_args()

    metadata, species_df = combine_portions(args.portion_dirs)

    # Save the results
    metadata_df = pd.DataFrame(metadata, index=[0])
    metadata_df.to_csv(f"{args.output_prefix}.combined_summary.tsv", index=False, sep='\t')
    species_df.to_csv(f"{args.output_prefix}.species_occurrences.tsv", sep='\t')

    print(f"Combined metadata and species occurrence files have been generated.")
