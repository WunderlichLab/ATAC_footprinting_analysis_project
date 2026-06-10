#!/bin/bash

# ==============================================================================
# Script Name: bed_intersect.sh
# Description: Intersects Transcription Factor Binding Sites (TFBS) with targeted
#              genomic regions of interest using BEDtools, followed by unique coordinate
#              sorting and feature quantification.
# Author: Christian Mei
# ==============================================================================

# --- INPUT VALIDATION ---
# Check if the correct number of arguments is provided (Requires File A and File B)
if [ "$#" -ne 2 ]; then
    echo "Usage: ./bed_intersect.sh <fileA.bed> <fileB.bed>"
    exit 1
fi

# --- VARIABLE ASSIGNMENT ---
# Capture genomic file paths supplied as positional command-line arguments
fileA="$1"  # Input Query: Predicted TFBS regions (e.g., individual TOBIAS output)
fileB="$2"  # Input Target: Genomic intervals of interest (e.g., CREs/CRMs, enhancers)

# Define the absolute local file path destination for the processed intersection results
outputFile="/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Script_outputs/intersect_output.bed"


# --- GENOMIC INTERSECTION & SORTING ---
# bedtools intersect flags:
#   -wa : Write the original entry from A for each overlapping interval found in B.
#   -v  : Write the entries from A that don't overlap with B. Manually replace -wa if this is the desired optionc
#   -a  : Query genomic features file (TFBS).
#   -b  : Target genomic features database file (CRE intervals).
# sort flags:
#   -u  : Unique mode; discards duplicate lines after evaluating keys.
#   -k  : Defines specific keys for sorting (k1,1 = Chromosome; k2,2 = Start; k3,3 = End).
bedtools intersect -wa -a "$fileA" -b "$fileB" | sort -u -k1,1 -k2,2 -k3,3 > "$outputFile"


# --- QUANTIFICATION & SUMMARY ---
# Calculate the total number of lines contained in the generated output BED file
lineCount=$(wc -l < "$outputFile")

# Print final feature yield metrics to standard output
echo "Number of lines in the output file: $lineCount"
