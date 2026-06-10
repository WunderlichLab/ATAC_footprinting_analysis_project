import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

# --- DATA LOADING & PREPARATION ---

# Define the file path to the BED file containing coordinates of TFBS overlapping with CREs using bedtools intersect
spacing_file_path = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/D1-D4_Rep1+2_TFBS/D1+D4_Rep1_2_TFBS_CRM_spacing_input.bed"

# Load the tab-delimited BED file into a pandas DataFrame without a header row
df = pd.read_csv(spacing_file_path, sep='\t', header=None)

# Assign descriptive column names to mapping coordinates of both TFBS and CREs
df.columns = ['TFBS_chr', 'TFBS_start', 'TFBS_end', 'CRE_chr', 'CRE_start', 'CRE_end', 'CRE']

# Sort the dataset first by CRE ID, then chronologically by the starting coordinate of each TFBS
df.sort_values(by=['CRE', 'TFBS_start'], inplace=True)


# --- SPACING CALCULATION LOGIC ---

def calculate_spacing(group):
    """
    Calculates the genomic distance between consecutive TFBS features within a single CRE group.
    Overlapping features are assigned a value of -1.
    """
    # Ensure the group is strictly ordered by genomic start position and reset index for positional indexing
    group = group.sort_values('TFBS_start').reset_index(drop=True)

    # Initialize list; the first TFBS in a CRE has no upstream neighbor to calculate spacing against
    spacing_list = [None]

    # Iterate through the group starting from the second TFBS feature
    for i in range(1, len(group)):
        prev_end = group.loc[i - 1, 'TFBS_end']
        curr_start = group.loc[i, 'TFBS_start']

        # Determine if features overlap or calculate the gap size between them
        if curr_start <= prev_end:
            spacing = -1  # Flag indicating overlapping TFBS features
        else:
            spacing = curr_start - prev_end  # Calculate base pair distance between features

        spacing_list.append(spacing)

    # Append the calculated metrics as a new column to the group DataFrame
    group['TFBS_spacing'] = spacing_list
    return group


# Apply the spacing calculation independently to each subset grouped by unique CRE IDs
df_with_spacing = df.groupby('CRE', group_keys=False).apply(calculate_spacing)

# --- STATISTICAL ANALYSIS ---

# Drop the 'None' entries corresponding to the first TFBS elements of each CRE
valid_spacing = df_with_spacing['TFBS_spacing'].dropna()

# Filter and isolate purely numerical values into a NumPy array for mathematical operations
only_spacing = np.array([x for x in valid_spacing if isinstance(x, (int, float))])

# Calculate summary statistics across all valid spacing distances
spacing_mean = np.mean(only_spacing)
spacing_median = np.median(only_spacing)
spacing_range = max(only_spacing)

# Output summary metrics to the console
print(f"Mean TFBS Spacing: {spacing_mean}")
print(f"Median TFBS Spacing: {spacing_median}")
print(f"Max TFBS Spacing: {spacing_range}")

# --- DATA TRANSFORMATION & VISUALIZATION ---

# Transform the distribution using a log scale to handle data skewness
only_spacing_modded: list = []
modded_value: int = 0
for i in only_spacing:
    # Offset by +2 to shift the -1 overlap flags to +1, preventing log domain errors at or below zero
    modded_value = math.log(2 + i)
    only_spacing_modded.append(modded_value)
    modded_value = 0

# Generate and display a box plot of the log-transformed spacing values
plt.boxplot(only_spacing_modded)
plt.xlabel('TFBS in known CREs')
plt.ylabel('Log Spacing (bp)')
plt.show()

# --- COMMENTED OUT HISTOGRAM ANALYSIS ---
'''
# Alternative visualization utilizing a log-scaled histogram to view raw spacing frequencies
plt.hist(only_spacing, bins = 10, color = '#cc0000', edgecolor = 'white', log = True)
plt.title('TFBS Spacing')
plt.xlabel('TFBS spacing (bp)')
# Mark the median distribution point with a vertical dashed line
plt.axvline(np.median(only_spacing), color='k', linestyle='dashed', linewidth=1)
plt.ylabel('Frequency')
plt.show()
'''
# DATA EXPORT
# Define the target directory and specific filename requested
output_spacing_path = "/Users/christianmei/Documents/Boston_University/Fall_2022_Research/TF_prediction_results/Hypotheses/uncharacterized_crm_hypothesis/cluster_medium_H3K4me1/2026_06_08_TFBS_CRM_spacing.csv"

# Convert the 1D numpy array into a structured single-column DataFrame
output_df = pd.DataFrame({"TFBS_Spacing": only_spacing})

# Save the data directly to a clean comma-separated CSV file without row indices
output_df.to_csv(output_spacing_path, index=False, sep=",")

print(f"Spacing values successfully saved to: {output_spacing_path}")
