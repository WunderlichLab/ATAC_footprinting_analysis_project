"""
Script Name: percent_conservation.py
Description: Evaluates and correlates nucleotide conservation percentages between individual
             TFBS tracks and their parent macro-cluster regions. Generates publication-ready
             scatter plots, distribution histograms, filters candidate regions at a 70% threshold,
             and exports the processed dataset to a structured CSV.
Author: Christian Mei
"""

import csv
import matplotlib.pyplot as plt

# --- GLOBAL FILE PATH DEFINITIONS ---
tfbs_cluster_file = "/Users/christianmei/Documents/Boston_University/Fall_2022_Research/TF_prediction_results/Hypotheses/uncharacterized_crm_hypothesis/cluster_medium_H3K4me1/complete_cluster_TFBS_conservation_map.bed"
region_cluster_file = "/Users/christianmei/Documents/Boston_University/Fall_2022_Research/TF_prediction_results/Hypotheses/uncharacterized_crm_hypothesis/cluster_medium_H3K4me1/complete_cluster_region_conservation_map.bed"


# --- DATA PARSING & PROCESSING FUNCTIONS ---

def extract_TFBS(file_path: str) -> dict:
    """
    Parses TFBS intersection files to calculate individual feature footprints and
    associated PhastCons conservation lengths grouped by unique cluster IDs.

    Parameters:
        file_path (str): Path to the TFBS conservation BED file.

    Returns:
        dict: Format -> clusterID: [total_tfbs_bp, conserved_tfbs_bp]
    """
    cluster_dict = {}

    with open(file_path, 'r') as file:
        for line in file:
            split_line = line.strip().split()

            # Compute absolute span of the current individual TFBS footprint
            total_bp = int(split_line[2]) - int(split_line[1])
            conserved_bp = 0

            # Quantify overlapping conserved base pairs if intersection coordinates exist
            if len(split_line) >= 8:
                conserved_bp = int(split_line[6]) - int(split_line[5])

            cluster_id = split_line[3]

            # Accumulate raw base pair calculations under the corresponding cluster entry
            if cluster_id in cluster_dict:
                cluster_dict[cluster_id][0] += total_bp
                cluster_dict[cluster_id][1] += conserved_bp
            else:
                cluster_dict[cluster_id] = [total_bp, conserved_bp]

    return cluster_dict


def merge_intervals(intervals: list) -> list:
    """
    Collapses overlapping or adjacent coordinate ranges to ensure accurate
    base pair calculations across complex, parsed genomic intervals.

    Parameters:
        intervals (list): Unsorted list of (start, end) integer tuples.

    Returns:
        list: Sorted list of consolidated, non-overlapping coordinate tuples.
    """
    if not intervals:
        return []

    # Sort intervals hierarchically based on structural start position
    intervals.sort(key=lambda x: x[0])
    merged = [intervals[0]]

    for current_start, current_end in intervals[1:]:
        last_end = merged[-1][1]

        # Resolve overlaps by updating the boundary limit of the current block
        if current_start <= last_end:
            merged[-1] = (merged[-1][0], max(last_end, current_end))
        else:
            merged.append((current_start, current_end))
    return merged


def calculate_total_bp(intervals: list) -> int:
    """Calculates cumulative nucleotide distance across processed coordinate groups."""
    return sum(end - start for start, end in intervals)


def extract_region(file_path: str) -> dict:
    """
    Parses macro-cluster regional maps, applying interval deduplication
    to resolve actual global conservation ratios across complete spans.

    Parameters:
        file_path (str): Path to the macro-cluster region conservation BED file.

    Returns:
        dict: Format -> clusterID: [total_region_bp, conserved_region_bp]
    """
    cluster_dict = {}

    with open(file_path, 'r') as file:
        for line in file:
            split_line = line.strip().split()
            cluster_id = split_line[3]
            total_interval = (int(split_line[1]), int(split_line[2]))
            conserved_interval = (int(split_line[5]), int(split_line[6])) if len(split_line) >= 7 else None

            # Catalog individual feature intervals prior to interval reduction
            if cluster_id not in cluster_dict:
                cluster_dict[cluster_id] = {
                    'total': [total_interval],
                    'conserved': [] if conserved_interval is None else [conserved_interval]
                }
            else:
                cluster_dict[cluster_id]['total'].append(total_interval)
                if conserved_interval:
                    cluster_dict[cluster_id]['conserved'].append(conserved_interval)

    # Perform interval reduction and summarize absolute base pair lengths
    for cluster_id, intervals in cluster_dict.items():
        merged_total = merge_intervals(intervals['total'])
        merged_conserved = merge_intervals(intervals['conserved'])

        total_bp = calculate_total_bp(merged_total)
        conserved_bp = calculate_total_bp(merged_conserved)

        cluster_dict[cluster_id] = [total_bp, conserved_bp]

    return cluster_dict


# --- DATA PROCESSING PIPELINE ---

tfbs_dict = extract_TFBS(tfbs_cluster_file)
region_dict = extract_region(region_cluster_file)

tfbs_percentage_array = []
region_percentage_array = []
key_array = []
conservation_map = {}

# Compute percentage metrics for regional spans versus functional footprints
for key, values in tfbs_dict.items():
    tfbs_percent_conserved = (values[1] / values[0]) * 100
    tfbs_percentage_array.append(tfbs_percent_conserved)
    key_array.append(key)

    if key in region_dict:
        region_percent_conserved = (region_dict[key][1] / region_dict[key][0]) * 100
        region_percentage_array.append(region_percent_conserved)

    conservation_map[key] = (tfbs_percent_conserved, region_percent_conserved)

print("Total number of successfully mapped clusters:", len(key_array))
print("Conservation mapping outputs:", conservation_map)

# --- INTEGRITY CONTROL CHECK ---

tfbs_keys = set(tfbs_dict.keys())
region_keys = set(region_dict.keys())

unmatched_tfbs_keys = tfbs_keys - region_keys
unmatched_region_keys = region_keys - tfbs_keys

print("Unmatched TFBS keys:", unmatched_tfbs_keys)
print("Unmatched Region keys:", unmatched_region_keys)

# --- DATA VISUALIZATION CONSTRUCTION ---

# Figure 1: Covariation Scatter Plot (Footprint vs Region Conservation Profile)
plt.figure(figsize=(6.5, 5.5))
plt.scatter(region_percentage_array, tfbs_percentage_array, c="#e69138", alpha=0.75, edgecolors='none', s=45)
plt.xlabel("% Cluster Region bp Conserved", fontsize=11, labelpad=8)
plt.ylabel("% TFBS bp Conserved", fontsize=11, labelpad=8)
plt.title("cCRE Conservation Profile Covariation", fontsize=13, pad=12, fontweight='semibold')
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()
plt.show()

# Figure 2: Regional Conservation Frequency Distribution Profile
plt.figure(figsize=(6.5, 5.0))
plt.hist(region_percentage_array, color="#e69138", edgecolor="white", bins=15, linewidth=0.8)
plt.xlabel("% Conservation", fontsize=11, labelpad=8)
plt.ylabel("Frequency", fontsize=11, labelpad=8)
plt.title("cCRE Macro-Cluster Conservation Frequency Distribution", fontsize=13, pad=12, fontweight='semibold')

# Demarcate candidate regulatory hub selection triage threshold (>= 70%)
plt.axvline(70, color='black', linestyle='dashed', linewidth=1.2, label='Selection Threshold (70%)')
plt.legend(loc='upper left', frameon=False)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()
plt.show()

# --- CANDIDATE SELECTION FILTRATION ---

# Filter out highly stringent candidate windows matching project guidelines
high_conservation_cluster: list = []
for key, values in conservation_map.items():
    if values[1] >= 70:
        high_conservation_cluster.append(key)

print("Validated Hyper-Conserved Cluster Candidates (>= 70%):", high_conservation_cluster)
print("All processed region conservation arrays:", region_percentage_array)

# --- DATA EXPORT MANAGEMENT ---

output_table_path = "/Users/christianmei/Documents/Boston_University/Fall_2022_Research/TF_prediction_results/Hypotheses/uncharacterized_crm_hypothesis/cluster_medium_H3K4me1/2026_06_08_cluster_conservation_summary.csv"

# Write statistical output summaries directly to disk
with open(output_table_path, "w", newline="", encoding="utf-8") as csvfile:
    writer = csv.writer(csvfile)

    # Establish descriptive header schemas for reproducibility
    writer.writerow(["Cluster_ID", "TFBS_Conservation_Percentage", "Region_Conservation_Percentage"])

    # Loop through the compiled dictionary, unpacking structural tuple records
    for cluster_id, (tfbs_percent, region_percent) in conservation_map.items():
        writer.writerow([cluster_id, tfbs_percent, region_percent])

    # Flush the memory stream buffer to safely verify transaction integrity on local media
    csvfile.flush()

print(f"Table successfully saved to: {output_table_path}")