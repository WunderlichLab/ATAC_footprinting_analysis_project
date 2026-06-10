"""
Script Name: TFBS_per_enhancer.py
Description: Quantifies and visualizes the distribution of Transcription Factor Binding
             Sites (TFBS) across collected Cis-Regulatory Modules (CRMs)/enhancers.
             Handles duplicate entry filtration, map generation, and summary statistics.
Author: Christian Mei
"""

import matplotlib.pyplot as plt
import numpy as np


def tfbs_per_CRM(file_path: str) -> dict:
    """
    Parses a genomic mapping file, removes duplicate records, and counts the
    total number of TFBS associated with each unique CRM coordinate set.

    Parameters:
        file_path (str): Path to the tab-delimited intersection file containing
                         TFBS and CRM features.

    Returns:
        dict: A mapping dictionary where keys are CRM coordinate tuples
              (start, end, chromosome/identifier) and values are integer TFBS counts.
    """
    tfbs_crm_map = {}
    no_dups_list = []

    # --- Data Ingestion & Type Casting ---
    with open(file_path, 'r') as file:
        raw_list = []
        for line in file:
            split_line = line.strip().split()
            processed_line = []

            # Cast coordinates to integers (indices 1, 2, 5, 6) and leave labels as strings
            for index, item in enumerate(split_line):
                if index in [1, 2, 5, 6]:
                    processed_line.append(int(item))
                else:
                    processed_line.append(item)
            raw_list.append(processed_line)

    # --- Deduplication ---
    # Filter out redundant lines to ensure exact intersection occurrences are unique
    for item in raw_list:
        if item not in no_dups_list:
            no_dups_list.append(item)
    print("Total unique combinations: ", len(no_dups_list))

    # --- Frequency Mapping ---
    # Aggregate TFBS counts grouped by unique genomic CRM locations
    for item in no_dups_list:
        # Define CRM unique identifier tuple using target coordinate indices
        crm_id: tuple = (item[1], item[2], item[3])
        if crm_id not in tfbs_crm_map:
            tfbs_crm_map[crm_id] = 1
        else:
            tfbs_crm_map[crm_id] += 1

    return tfbs_crm_map


def tfbs_count_dist(tfbs_CRM_map: dict) -> None:
    """
    Calculates summary statistics and displays a publication-quality histogram
    of the TFBS density distribution across all processed CRMs.

    Parameters:
        tfbs_CRM_map (dict): CRM coordinate keys mapped to their respective TFBS counts.

    Returns:
        None
    """
    # Extract absolute frequency values from the dictionary mapping
    tfbs_count = list(tfbs_CRM_map.values())

    # --- Statistical Calculations ---
    print("Mean: ", np.mean(tfbs_count))
    print("Median: ", np.median(tfbs_count))

    # --- Data Visualization Construction ---
    plt.figure(figsize=(8, 6))  # Established explicit publication layout dimensions

    # Plot frequency distribution with sharp edge contrast boundaries
    plt.hist(tfbs_count, bins=20, color='#db5856', edgecolor='white', linewidth=0.8)

    # Axis typography setup
    plt.xlabel("TFBS per Known CRE", fontsize=12, labelpad=10)
    plt.ylabel("Frequency", fontsize=12, labelpad=10)
    plt.title("TFBS Distribution per Known CRE (Stage 4-6)", fontsize=14, pad=15, fontweight='semibold')

    # Plot population median threshold vertical indicator line
    plt.axvline(np.median(tfbs_count), color='black', linestyle='dashed', linewidth=1.2,
                label=f'Median: {np.median(tfbs_count)}')

    # Layout and frame optimization
    plt.legend(loc='upper right', frameon=False)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()

    plt.show()
    return None


#----CODE EXECUTION AND EXPORTING----
# Define the file paths
tfbs_crm_file = "/Users/christianmei/Documents/Boston_University/Fall_2022_Research/TF_prediction_results/Hypotheses/uncharacterized_crm_hypothesis/CRM_TFBS_map.txt"
output_counts_path = "/Users/christianmei/Documents/Boston_University/Fall_2022_Research/TF_prediction_results/Hypotheses/uncharacterized_crm_hypothesis/cluster_medium_H3K4me1/2026_06_08_tfbs_counts_per_crm.csv"

# 1. Run the mapping function and capture the dictionary variable
crm_mapping_data = tfbs_per_CRM(tfbs_crm_file)

# 2. Pass the data to your existing histogram function
tfbs_count_dist(crm_mapping_data)

# 3. Write out only the counts (the values) to a single-column CSV
with open(output_counts_path, "w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)

    # Write a header for the single data column
    writer.writerow(["TFBS_Count"])

    # Loop strictly through the map values to isolate the raw distribution counts
    for tf_count in crm_mapping_data.values():
        writer.writerow([tf_count])

    # Flush memory buffer to guarantee save execution complete
    csvfile.flush()

print(f"Distribution count table successfully saved to: {output_counts_path}")
