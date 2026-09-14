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
│   ├── 04_DEs_Interaction.R
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
source("scripts/pipeLine_DESeq2_SingleVariableInput.R")
```

### Step 2: Translation Efficiency Analysis

```r
source("scripts/pipeline_TE_basedOnDE_byCondtions.R")
```

### Step 3: TE Interactive Coefficient Analysis

```r
source("scripts/TE_customGeneListLabel_Plot.R")
```

### Step 4: DEs Interaction Analysis

```r
source("scripts/DEs_Interaction.R")
```

### Step 5: Functional Enrichment

```r
source("scripts/gsea.R")
```

## Notes

This repository contains analysis and visualization scripts only.

Raw sequencing data, reference genomes, alignments, and large intermediate files are intentionally excluded.

Pre-processing steps performed using standard community software (e.g., FASTX, STAR, featureCounts, RSEM) should be described in the associated manuscript Methods section.

The inputs and parameters settings might need adjustment for the specific projects and applications

Example visualization: certain code demonstrates the general structure of an RNA-versus-RPF correlation plot, which means that some example code is simplified and does not represent the complete analysis implementation used in the study

## Citation

If you use this repository, please cite the associated publication.

