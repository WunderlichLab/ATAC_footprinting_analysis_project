"""
Script Name: RAMPAGE_to_promoter.py
Description: Converts raw TSS coordinates (from RAMPAGE or GFF data) into flanking promoter 
             windows by calculating peak centers and applying strand-specific offsets.
Author: Christian Mei
"""

# --- GLOBAL FILE PATH DEFINITIONS ---
# RAMPAGE Transcription Start Site Cluster dataset
rampage_file: str = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Genome/RAMPAGE_promoter/GSE89299_Dmel2_RAMPAGE_TSC.bed"

# Reference Dmel r5.57 GFF-derived annotated TSS dataset
r5_57_gff_tss_file: str = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Genome/dmel_r5.57_gff_based_genome_features/polished_genomic_features_no_chrU/dmel-all-r5.57-TSS_chr.bed"


def read_file(input_file: str) -> list:
    """
    Parses a standard BED file, automatically casting genomic coordinates 
    (start and end positions) into integers for downstream arithmetic.

    Parameters:
        input_file (str): Path to the tab-delimited genomic features file.

    Returns:
        list: A nested list containing parsed rows with casted coordinate fields.
    """
    txt_list = []
    with open(input_file, 'r') as file:
        for line in file:
            split_line = line.strip().split()
            
            # Explicitly cast structural coordinate indices to integers
            split_line[1] = int(split_line[1])
            split_line[2] = int(split_line[2])
            txt_list.append(split_line)
    return txt_list


# --- DATA INGESTION ---
rampage_file_list = read_file(rampage_file)
r5_57_gff_tss_list = read_file(r5_57_gff_tss_file)


def tss_to_promoter_300_60(tss_list: list) -> list:
    """
    Generates promoter windows defined as 300 bp upstream (-300) to 
    60 bp downstream (+60) relative to the calculated TSS peak center.

    Parameters:
        tss_list (list): Nested list containing parsed BED file elements.

    Returns:
        list: A standard 5-column BED-formatted list [chr, start, end, id, strand].
    """
    promoter_list: list = []

    for item in tss_list:
        # Resolve spatial peak range into a unified zero-point anchor center
        promoter_center = int((item[1] + item[2]) / 2)
        
        # Calculate flanking window dimensions based on sense/antisense strand orientation
        if item[5] == "+":
            promoter_start = promoter_center - 300
            promoter_end = promoter_center + 60
            promoter_list.append([item[0], promoter_start, promoter_end, item[3], item[5]])
        elif item[5] == "-":
            promoter_start = promoter_center - 60
            promoter_end = promoter_center + 300
            promoter_list.append([item[0], promoter_start, promoter_end, item[3], item[5]])
            
    return promoter_list


def tss_to_promoter_250_50(tss_list: list) -> list:
    """
    Generates promoter windows defined as 250 bp upstream (-250) to 
    50 bp downstream (+50) relative to the calculated TSS peak center.

    Parameters:
        tss_list (list): Nested list containing parsed BED file elements.

    Returns:
        list: A standard 5-column BED-formatted list [chr, start, end, id, strand].
    """
    promoter_list: list = []

    for item in tss_list:
        # Resolve spatial peak range into a unified zero-point anchor center
        promoter_center = int((item[1] + item[2]) / 2)
        
        # Calculate flanking window dimensions based on sense/antisense strand orientation
        if item[5] == "+":
            promoter_start = promoter_center - 250
            promoter_end = promoter_center + 50
            promoter_list.append([item[0], promoter_start, promoter_end, item[3], item[5]])
        elif item[5] == "-":
            promoter_start = promoter_center - 50
            promoter_end = promoter_center + 250
            promoter_list.append([item[0], promoter_start, promoter_end, item[3], item[5]])
            
    return promoter_list


def export_to_bed(data_list: list, output_file: str) -> None:
    """
    Writes processed promoter coordinate structures to disk as a standard tab-delimited BED file.

    Parameters:
        data_list (list): Nested list containing structured chromosome interval features.
        output_file (str): Output destination path.

    Returns:
        None
    """
    with open(output_file, 'w') as file:
        for entry in data_list:
            # Stringify all primitives and join with standard tab delimiters
            line = "\t".join(map(str, entry))
            file.write(line + "\n")
    return None


# --- PIPELINE EXECUTION ---
if __name__ == "__main__":
    # 1. Generate and export -300/+60 bp promoters from RAMPAGE peaks (Commented by default)
    # rampage_out_300_60 = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Genome/RAMPAGE_promoter/GSE89299_Dmel2_RAMPAGE_promoter_300_60.bed"
    # export_to_bed(tss_to_promoter_300_60(rampage_file_list), rampage_out_300_60)

    # 2. Generate and export -250/+50 bp promoters from RAMPAGE peaks
    rampage_out_250_50 = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Genome/RAMPAGE_promoter/GSE89299_Dmel2_RAMPAGE_250_50_promoter.bed"
    export_to_bed(tss_to_promoter_250_50(rampage_file_list), rampage_out_250_50)

    # 3. Generate and export -250/+50 bp promoters from reference GFF-derived TSS annotations
    gff_out_250_50 = "/Users/christianmei/Desktop/Fall_2022_Research/TF_prediction_results/Genome/dmel_r5.57_gff_based_genome_features/polished_genomic_features_no_chrU/dmel-all-r5.57-TSS_promoter.bed"
    export_to_bed(tss_to_promoter_250_50(r5_57_gff_tss_list), gff_out_250_50)
