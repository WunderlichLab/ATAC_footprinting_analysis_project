#!/bin/bash

# ==============================================================================
# Script Name: subdir_loop.sh
# Description: Mounts external media, iterates through subdirectories to convert
#              GFF3 features to BED format, reheaders BAM alignments, and
#              standardizes chromosome naming conventions to UCSC format (e.g., 'chr2L').
#              This is necessary for correct TOBIAS processing and downstream analysis
# Author: Christian Mei
# ==============================================================================

# --- ENVIRONMENT PREPARATION ---
# Mount target external volume containing raw genomic datasets
diskutil mount /Volumes/MEI_DROS

# Navigate into the target project data directory
cd /Volumes/MEI_DROS/Drosophila_Materials_Drive/Drosophila_Dataset || exit 1


# --- PROCESSING LOOP ---
# Iterate through every immediate subdirectory
for dir in */; do
    echo "Processing directory: $dir"

    # Enter target subdirectory; use logical '&&' to prevent execution in parent directory on failure
    cd "$dir" && {

        # --- STEP 1: GFF3 TO UCSC BED CONVERSION ---
        for file_gff in *.gff3; do
            # Guard against empty glob matches if no GFF3 files exist
            [ -e "$file_gff" ] || continue

            # Extract standard file prefix name
            base_name=$(basename "$file_gff" .gff3)

            # Convert raw GFF3 features to 6-column BED intervals
            gff2bed < "$file_gff" > "${base_name}.bed"

            # Map canonical FlyBase chromosome prefixes into standard UCSC nomenclature
            sed -e 's/^Y/chrY/'   -e 's/^X/chrX/'   \
                -e 's/^2L/chr2L/' -e 's/^2R/chr2R/' \
                -e 's/^3L/chr3L/' -e 's/^3R/chr3R/' \
                -e 's/^4/chr4/'   -e 's/^M/chrM/'   \
                "${base_name}.bed" > "${base_name}.chr.bed"
        done

        # --- STEP 2: BAM REHEADER & SANITIZATION ---
        for file_bam in *.bam; do
            # Guard against empty glob matches if no BAM files exist
            [ -e "$file_bam" ] || continue

            base_bam=$(basename "$file_bam" .bam)

            # Stream, modify, and rewrite the BAM header structural sequence tags (@SQ SN:)
            samtools view -H "$file_bam" | \
                sed -e 's/SN:2L/SN:chr2L/' -e 's/SN:2R/SN:chr2R/' \
                    -e 's/SN:3L/SN:chr3L/' -e 's/SN:3R/SN:chr3R/' \
                    -e 's/SN:4/SN:chr4/'   -e 's/SN:X/SN:chrX/'   \
                    -e 's/SN:Y/SN:chrY/'   -e 's/SN:M/SN:chrM/'   | \
                samtools reheader - "$file_bam" > "./${base_bam}.chr.bam"

            # Purge deprecated BAM index (.bai) tracks to avoid tracking mismatches
            rm -f *.bai
        done

        # --- REVERT CONTEXT ---
        # Return to the parent root directory to advance to the next loop iteration safely
        cd ..
    }
done