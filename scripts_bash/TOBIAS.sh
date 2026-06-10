#!/bin/bash

# ==============================================================================
# Script Name: TOBIAS.sh
# Description: Automates the TOBIAS (Transcription Factor Occupancy Prediction By
#              ATAC-seq Signal) framework across all subdirectories. Performs
#              ATAC-seq bias correction, footprint scoring, and motif detection.
# Author: Christian Mei
# ==============================================================================

# Iterate through every immediate subdirectory within the current working directory
for dir in */; do

    # Enter the target subdirectory; only proceed if the directory shift succeeds
    cd "$dir" && {

        # --- PRE-PROCESSING ---
        # Purge existing BAM index files (.bai) to avoid index-mismatch calculation bugs
        find . -type f -name "*.bai" -delete

        # --- STEP 1: ATAC-seq BIAS CORRECTION ---
        # Corrects for Tn5 transposase sequence insertion bias to isolate true footprints
        TOBIAS ATACorrect \
            --bam "$(find . -name "*_chr.bam*")" \
            --genome /Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/TOBIAS_Bozek_Rep1/final_dm3.fa \
            --peaks "$(find . -name "*chr.bed")" \
            --blacklist /Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/TOBIAS_Bozek_Rep2/dm3-blacklist.v2.bed \
            --outdir ./ATACorrect_test \
            --cores 8

        # --- STEP 2: FOOTPRINT SCORE CALCULATION ---
        # Computes continuous footprint scores across candidate regulatory regions using bias-corrected signals
        TOBIAS FootprintScores \
            --signal "$(find . -name "*_corrected.bw*")" \
            --regions "$(find . -name "*chr.bed")" \
            --output ./footprints.bw \
            --cores 8

        # --- STEP 3: MOTIF DETECTION & BINDING SITE PREDICTION ---
        # Integrates MEME motifs with footprint depth to identify significantly bound genomic sites
        TOBIAS BINDetect \
            --motifs /Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/TOBIAS_Bozek_Rep2/motif2.meme \
            --motif_pvalue 0.001 \
            --signals footprints.bw \
            --genome /Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/TOBIAS_Bozek_Rep1/final_dm3.fa \
            --peaks "$(find . -name "*chr.bed")" \
            --outdir ./BINDetect_output \
            --cond_names "$(basename "$dir")" \
            --cores 8

        # --- PIPELINE CLEANUP ---
        # Revert back to the parent directory to safely execute the next iteration loop
        cd - > /dev/null
    }
done