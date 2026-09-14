# RNA_Profile
RNAseq and RPFseq based DE/TE analysis

# RNA-seq Differential Expression analysis (DE) and Translation Efficiency Analysis (TE) Pipeline

## Overview

This repository contains R scripts used for differential expression (DE), translation efficiency (TE), gene set enrichment analysis (GSEA), and data visualization in bulk RNA-seq and ribosome profiling studies.

The pipeline utilizes widely adopted bioinformatics tools and statistical methods, including:

* DESeq2
* RibosomeProfiling
* clusterProfiler
* fgsea
* EnhancedVolcano
* ggplot2

## Repository Structure

```text
RNAseq-DE-TE-Pipeline/
│
├── README.md
├── LICENSE
│
├── scripts/
│   ├── 01_pipeLine_DESeq2_SingleVariableInput.R
│   ├── 02_pipeline_TE_basedOnDE_byCondtions.R
│   ├── 03_TE_customGeneListLabel_Plot.R
│   ├── 04_pipeline_TE_basedOnDE_byCondtions.R
│   └── 05_gsea.R
│
├── Input/
│   └── count_tables.txt
│   └── sample_metadata_example.txt
│
└── docs/
    └── Main algorithms; Reference 
```

## Requirements

### R

Tested with:

```r
R >= 4.3
```

### Packages

```r
DESeq2
ggplot2
EnhancedVolcano
clusterProfiler
fgsea
dplyr
tidyr
readr
```

## Input Files

The pipeline expects:

1. Gene count matrix
2. Sample metadata table
3. Optional TE input files derived from RNA-seq and Ribo-seq quantification

Example (please check the specific RNA and RPF related tables in input folder):

```text
input/
├── counts_matrix.txt
└── sample_metadata.txt
```

## Workflow

### Step 1: Differential Expression

```r
source("scripts/01_DESeq2_DE.R")
```

### Step 2: Translation Efficiency Analysis

```r
source("scripts/02_TE_analysis.R")
```

### Step 3: Functional Enrichment

```r
source("scripts/03_GSEA.R")
```

### Step 4: Visualization

```r
source("scripts/04_visualization.R")
```

### Step 5: Export Results

```r
source("scripts/05_export_results.R")
```

## Notes

This repository contains analysis and visualization scripts only.

Raw sequencing data, reference genomes, alignments, and large intermediate files are intentionally excluded.

Pre-processing steps performed using standard community software (e.g., FASTX, STAR, featureCounts, RSEM) should be described in the associated manuscript Methods section.

## Citation

If you use this repository, please cite the associated publication.

