"""
Script Name: hCRE_conservation.py
Description: Computes the global nucleotide conservation percentage for high confidence candidate 
             Cis-Regulatory Modules (cCRMs) by aggregating overlapping PhastCons track segments 
             and normalizing by total cCRM locus size. This is to make sure that my final hCREs still maintain
             high conservation after expanding their genomic region to occupy the entire ATAC-seq peak
Author: Christian Mei
"""

# Absolute path to the bedtools intersect mapping file (cCRM coordinates cross-referenced with PhastCons blocks)
ccrm_conservation_map_file = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Hypotheses/uncharacterized_crm_hypothesis/candidate_enhancer_conservation_map.txt"

# Initialize tracking dictionary: {(cCRM_ID, cCRM_size): accumulated_conserved_bp}
ccrm_conservation: dict = {}

# --- PARSING AND AGGREGATION ---
with open(ccrm_conservation_map_file, 'r') as file:
    for line in file:
        split_line = line.strip().split()
        conserved_bp = 0

        # If an intersection exists, columns 5 and 6 contain the coordinates of the overlapping PhastCons segment
        if len(split_line) >= 9:
            conserved_bp = int(split_line[6]) - int(split_line[5])

        # Extract structural attributes of the target cCRM locus
        cCRM_ID = split_line[3]
        cCRM_size: int = int(split_line[2]) - int(split_line[1])
        cCRM_key = (cCRM_ID, cCRM_size)
        
        # Accumulate total conserved base pairs mapped to each unique cCRM locus
        if cCRM_key in ccrm_conservation:
            ccrm_conservation[cCRM_key] += conserved_bp
        else:
            ccrm_conservation[cCRM_key] = conserved_bp

print("Raw aggregated conservation metrics:", ccrm_conservation)


# --- METRIC NORMALIZATION & OUTPUT ---
# Calculate the final percentage of conservation relative to the absolute length of each cCRM
for key, value in ccrm_conservation.items():
    percentage = (value / key[1]) * 100
    print(f"Cluster: {key[0]} has a conservation of {percentage:.2f}%")# This script is to measure the % conservation for the region as a whole for my 5 candidate enhancers
# File for this study is a map based on bedtools intersects between cCRM regions and overlapping Phastcons chunks

ccrm_conservation_map_file = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Hypotheses/uncharacterized_crm_hypothesis/candidate_enhancer_conservation_map.txt"
ccrm_conservation: dict = {}
with open(ccrm_conservation_map_file, 'r') as file:
    for line in file:
        split_line = line.strip().split()
        conserved_bp = 0

        if len(split_line) >= 9:
            conserved_bp = int(split_line[6]) - int(split_line[5])

        cCRM_ID = split_line[3]
        cCRM_size : int = int(split_line[2]) - int(split_line[1])
        cCRM_key = (cCRM_ID, cCRM_size)
        if cCRM_key in ccrm_conservation:
            ccrm_conservation[cCRM_key] += conserved_bp
        else:
            ccrm_conservation[cCRM_key] = conserved_bp
print(ccrm_conservation)

for key, value in ccrm_conservation.items():
    percentage = (value / key[1])*100
    print("Cluster:", key[0], "has a conservation of", percentage, "%")
