#!/bin/bash

# ==============================================================================
# Script Name: bound_merge.sh
# Description: Iterates through biological replicate directories (_Rep1), identifies
#              motif-specific significantly bound TFBS tracks from TOBIAS BINDetect,
#              and consolidates them into a master BED file per replicate condition.
# Author: Christian Mei
# ==============================================================================

# Iterate through target experimental directories appending the specific replicate suffix
# Will need to manually change the suffix if Rep2 is the focus
for folder in *_Rep1; do

    # Verify the target structural match is a valid directory
    if [[ -d "$folder" ]]; then

        # Define the absolute expected subpath to TOBIAS classification outputs
        bindetect_folder="$folder/BINDetect_output"

        if [[ -d "$bindetect_folder" ]]; then
            # Initialize array to accumulate structural file paths matching target criteria
            bound_files=()

            # --- FILE HARVESTING ---
            # Recursively capture file strings ending in 'bound.bed', utilizing null-delimiters
            # (-print0) to safely prevent whitespace parsing syntax breaks in file names.
            while IFS= read -r -d '' file; do
                # Isolate verified positive footprint intervals by explicitly filtering out unbound genomic subsets
                if [[ "$file" != *"unbound.bed" ]]; then
                    bound_files+=("$file")
                fi
            done < <(find "$bindetect_folder" -name "*bound.bed" -type f -print0)

            # --- DATA CONSOLIDATION ---
            # Execute concatenation only if the array contains captured target files
            if [ ${#bound_files[@]} -gt 0 ]; then
                # Establish descriptive output file path named after the biological replicate root
                output_file="${folder}_bound.bed"

                # Stream and aggregate the payload of all discovered BED arrays into the unified target destination
                cat "${bound_files[@]}" > "$output_file"

                echo "Compiled file created: $output_file"
            else
                echo "No matching bound.bed features found within: $bindetect_folder"
            fi
        fi
    fi
done