# Pipeline for "Novel Drosophila cis-regulatory elements can be uncovered by footprinting transcription factor binding sites in ATAC-seq data"

This repository contains a complete pipeline for processing genomic datasets in *Drosophila melanogaster* to detect transcription factor (TF) footprints in previously published ATAC-seq data. We perform a genome-wide mapping of these TF footprints and even search for novel cis-regulatory elements from TF footprint clusters filtered by cluster size, active enahncer epigenetic signatures, evolutionary conservation across 15 species (PhastCons 15).

---

## Workflow Architecture

```mermaid
graph TD
    %% --- STYLE DEFINITIONS ---
    classDef input fill:#f9f,stroke:#333,stroke-width:2px,color:#000;
    classDef script fill:#bbf,stroke:#333,stroke-width:2px,color:#000;
    classDef process fill:#fff,stroke:#333,stroke-width:1px,color:#000;
    classDef output fill:#bfb,stroke:#333,stroke-width:2px,color:#000;
    classDef final fill:#fbc,stroke:#333,stroke-width:3px,color:#000;

    %% --- TITLE STYLE via HTML header tags ---
    subgraph Phase_1 ["<b><big>Phase 1: Active TF Footprint Isolation</big></b>"]
        direction TB
        A1["Raw ATAC-seq BAMs (*_chr.bam)"]:::input --> B1["Footprint Detection<br>TOBIAS v. 0.13.3 <br>(TOBIAS.sh)"]:::script
        A2["Genome Reference (final_dm3.fa)"]:::input --> B1
        A3["ATAC-seq Peaks (*chr.bed)"]:::input --> B1
        A4["dm3-blacklist.v2.bed"]:::input --> B1
        
        B1 --> C1["ATACorrect Bias Correction"]:::process
        C1 --> C2["FootprintScores Calculation"]:::process
        C2 --> C3["BINDetect Motif Binding (p < 0.001)"]:::process
        
        A5["MEME Motifs (motif2.meme)"]:::input --> C3
        C3 --> D1["Individual *bound.bed Tracks"]:::output
        
        D1 --> E1["compile_bound_tfbs.sh<br>(bound_merge)"]:::script
        E1 -.->|from right| F1["Consolidated Replicate Footprints<br>(*_bound.bed)"]:::output
    end

    subgraph Phase_2 ["<b><big>Phase 2: TF Footprint Genomic Mapping</big></b>"]
        direction TB
        G1["Raw FlyBase GFF3 Tracks"]:::input --> H1["Reformatting Drosophila Features<br>(subdir_loop.sh)"]:::script
        H1 --> I1["Standardized UCSC BED Features<br>(chr* Chromosome Layout)"]:::output
        
        I1 --> J1["Downstream Extraction<br>(Bash/Awk Parsing)"]:::process
        J1 --> K1["Polished Reference Genomic Features"]:::output
        
        G2["Online RAMPAGE Data"]:::input --> H2["tss_to_promoter_generator.py"]:::script
        H2 --> I2["Promoter Windows Generated<br>(-250/+50 bp)"]:::output
        
        K1 --> L1["BEDtools Merge Framework"]:::process
        I2 --> L1
        
        L1 --> M1["Collapsed Macro-Genomic Interest Regions"]:::output
        F1 --> N1["bed_intersect.sh<br>(bedtools intersect -wa)"]:::script
        M1 --> N1
        
        N1 --> O1["Unique Intersected TFBS Loci<br>(intersect_output.bed)"]:::output
    end

    subgraph Phase_3 ["<b><big>Phase 3: TF Footprint Clustering</big></b>"]
        direction TB
        O1 --> P1["TF Footprint information from CRE's <br> with at least 1 footprint <br> CRM_TFBS_map.txt"]:::input
        
        P1 --> Q1["tfbs_crm_distribution.py<br>(TFBS_per_CRM)"]:::script
        Q1 --> R1["TFBS Counts per Known Enhancer<br>(Median = 14)"]:::output
        
        P1 --> Q2["calculate_spacing Python Function"]:::script
        Q2 --> R2["TFBS Spacing Distribution Analysis<br>(-1 Overlap Flag / Log Transformation)"]:::output
        
        R1 --> S1["TF Footprint Cluster Search<br>(cluster_mapping.sh)"]:::script
        R2 --> S1
        
        A6["Histone Modification Tracks<br>(H3K4me1, H3K27ac, H3K4me3)"]:::input --> S1
        A7["PhastCons 15-Species Conservation Data"]:::input --> S1
        
        S1 --> T1["High-Density Clustering (d = 18 bp)"]:::process
        T1 --> T2["Stringency Filtration (>= 7 TFBS)"]:::process
        T2 --> T3["Active Enhancer Filtering<br>(H3K4me1 + H3K27ac)"]:::process
    end

    subgraph Phase_4 ["<b><big>Phase 4: Conservation Filtering & hCRE Curation</big></b>"]
        direction TB
        T3 --> U1["complete_cluster_region_conservation_map.bed"]:::output
        T3 --> U2["complete_cluster_TFBS_conservation_map.bed"]:::output
        
        U1 --> V1["Cluster Conservation Filtering<br> >=70% Conservation<br>(conservation_percentage_restored.py)"]:::script
        U2 --> V1
        
        V1 --> W1["Cluster Overlap Statistics & CSV Export"]:::output
        W1 --> X1["Manual Curation (Target Gene Expression Patterns)"]:::process
        
        X1 --> Y1["High Confidence cCREs (hCREs)"]:::output
        
        Y1 --> Z1["calculate_ccrm_conservation.py<br>(hCRE_conservation)"]:::script
        Y1 --> Z2["cCRE_hCRE_size.qmd"]:::script
        
        Z1 --> AA1["Final Quantified Evolutionary<br>Conservation Metrics"]:::final
        Z2 --> AA2["Locus Dimension & Size<br>Frequency Distributions"]:::final
    end

    linkStyle default stroke-width:4px;
```

---
## Adapting this Pipeline

To adapt this pipeline for your own use, you will need to change parts of each script manually:
- Path files
- What kind of BEDtools intersect option you want to use in bed_intersect.sh
- TF footprint clustering thresholds in cluster_mapping.sh


