#!/bin/bash

# ==============================================================================
# Script Name: cluster_mapping.sh
# Description: Identifies high-density TFBS clusters, filters them via epigenetic
#              histone mark intersections (H3K4me1/H3K27ac), maps evolutionary
#              conservation tracks (PhastCons), and compiles complete genomic maps.
# Author: Christian Mei
# ==============================================================================

# --- INPUT VALIDATION ---
# Ensure exactly two inputs are provided (TFBS and PhastCons BED coordinates)
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <TFBS BED file> <Conserved regions BED file>"
    exit 1
fi

# Assign command line arguments to variables
TFBS_bed="$1"
phastcons_bed="$2"

# --- EPIGENETIC HISTONE MARK REFERENCE PATHS ---
H3K4me1_bed="/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/hypothesis_1/Histone_marks/H3K4me1_c14a.bed"
H3K27ac_bed="/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/hypothesis_1/Histone_marks/H3K27ac_c14a.bed"
H3K4me3_bed="/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/hypothesis_1/Histone_marks/H3K4me3_c14a.bed"

# --- PIPELINE DIRECTORIES & FILE PATHS ---
base_dir="/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/hypothesis_1/cluster_stringent_histone"

TFBS_cluster_raw="$base_dir/TFBS_cluster_raw.bed"
TFBS_cluster="$base_dir/TFBS_cluster.bed"
merged_TFBS="$base_dir/merged_TFBS.bed"
cluster_regions="$base_dir/cluster_regions.bed"
cluster_regions_filtered="$base_dir/cluster_regions_filtered.bed"
merged_TFBS_filtered="$base_dir/merged_TFBS_filtered.bed"
merged_phastcons="$base_dir/merged_phastcons.bed"
conserved_cluster_TFBS="$base_dir/conserved_cluster_TFBS.bed"
conserved_cluster_regions="$base_dir/conserved_cluster_regions.bed"
not_conserved_cluster_TFBS="$base_dir/not_conserved_cluster_TFBS.bed"
not_conserved_cluster_regions="$base_dir/not_conserved_cluster_regions.bed"
TFBS_phastcons="$base_dir/TFBS_phastcons.bed"
cluster_phastcons="$base_dir/cluster_phastcons.bed"
complete_cluster_TFBS_conservation_map="$base_dir/complete_cluster_TFBS_conservation_map.bed"
complete_cluster_region_conservation_map="$base_dir/complete_cluster_region_conservation_map.bed"


# ==============================================================================
# PIPELINE EXECUTION STEPS
# ==============================================================================

# --- STEP 1: INITIAL TFBS CLUSTERING ---
# Sort by chromosome then coordinates, cluster features within 18 bp max gap distance,
# and isolate tracking attributes (column 18 contains assigned Cluster ID)
sort -k1,1 -k2,2n "$TFBS_bed" | bedtools cluster -d 18 | \
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$18}' > "$TFBS_cluster_raw"


# --- STEP 2: HIGH-DENSITY WINDOW FILTRATION ---
# Isolate highly stringent regulatory hubs containing >= 7 unique binding features
awk '{cluster_counts[$5]++} END {for (cluster in cluster_counts) if (cluster_counts[cluster] >=7) print cluster}' "$TFBS_cluster_raw" > "valid_clusters.txt"
while read -r cluster; do
    awk -v cluster_id="$cluster" '$5 == cluster_id' "$TFBS_cluster_raw" >> "$TFBS_cluster"
done < "valid_clusters.txt"
rm "valid_clusters.txt"


# --- STEP 3: DEDUPLICATE HOMOLOGOUS OVERLAPPING FEATURES ---
# Merge directly overlapping TFBS tracks within defined clusters while maintaining structural ID metadata
sort -k1,1 -k2,2n "$TFBS_cluster" | bedtools merge -c 5 -o min > "$merged_TFBS"


# --- STEP 4: DEFINE MACRO-CLUSTER INTEGRATION REGIONS ---
# Aggregate clustered elements into a single contiguous genomic locus block using identical distance constraints
sort -k1,1 -k2,2n "$TFBS_cluster" | bedtools merge -d 18 -c 5 -o min > "$cluster_regions"


# --- STEP 5: EPIGENETIC ACTIVE ENHANCER SIGNATURE SELECTION ---
# Sequential intersection filters to retain loci containing overlapping active enhancer marks (H3K4me1 AND H3K27ac)
bedtools intersect -wa -a "$cluster_regions" -b "$H3K4me1_bed" | \
bedtools intersect -wa -a - -b "$H3K27ac_bed" > "$cluster_regions_filtered"


# --- STEP 6: BACK-MAP INDIVIDUAL HIGH-STRINGENCY FEATURES ---
# Cross-reference individual non-overlapping TFBS entries belonging to the validated histone-filtered clusters
bedtools intersect -wa -a "$merged_TFBS" -b "$cluster_regions_filtered" > "$merged_TFBS_filtered"


# --- STEP 7: PHASTCONS TRACK RE-SEGMENTATION ---
# Merge overlapping conservation intervals across the 15-species comparative genome tracks
bedtools merge -i "$phastcons_bed" -c 5 -o min > "$merged_phastcons"

# Restrict the footprint of conservation regions directly to our active, clustered coordinates
bedtools intersect -a "$merged_phastcons" -b "$merged_TFBS_filtered" > "$TFBS_phastcons"
bedtools intersect -a "$merged_phastcons" -b "$cluster_regions_filtered" > "$cluster_phastcons"


# --- STEP 8: RESOLVE EVOLUTIONARILY CONSERVED COMPONENT MAPS ---
# Match elements containing intersection overlaps between conservation tracks and active sites
bedtools intersect -wa -wb -a "$merged_TFBS_filtered" -b "$TFBS_phastcons" > "$conserved_cluster_TFBS"
bedtools intersect -wa -wb -a "$cluster_regions_filtered" -b "$cluster_phastcons" > "$conserved_cluster_regions"


# --- STEP 9: RESOLVE NON-CONSERVED SUBSETS & INTEGRATE MASTER MAPS ---
# Isolate fast-evolving/non-conserved elements using inverse filters (-v) and concatenate outputs
bedtools intersect -v -a "$merged_TFBS_filtered" -b "$TFBS_phastcons" > "$not_conserved_cluster_TFBS"
bedtools intersect -v -a "$cluster_regions_filtered" -b "$cluster_phastcons" > "$not_conserved_cluster_regions"

# Consolidate conservation layers to establish final multi-species reference maps
cat "$conserved_cluster_TFBS" "$not_conserved_cluster_TFBS" > "$complete_cluster_TFBS_conservation_map"
cat "$conserved_cluster_regions" "$not_conserved_cluster_regions" > "$complete_cluster_region_conservation_map"


# --- PIPELINE TERMINATION ---
echo "Pipeline completed. Outputs are located at:"
echo "Output 1: $complete_cluster_TFBS_conservation_map"
echo "Output 2: $complete_cluster_region_conservation_map"