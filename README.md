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

## Installation

### 1. Install R

Install R (version 4.3 or later) from [CRAN](https://cran.r-project.org/).

### 2. Install Required R Packages

Install the required CRAN and Bioconductor packages using the following commands:

```r
# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install CRAN packages
install.packages(c(
    "ggplot2",
    "dplyr",
    "tidyr",
    "readr"
))

# Install Bioconductor packages
BiocManager::install(c(
    "DESeq2",
    "RibosomeProfiling",
    "EnhancedVolcano",
    "clusterProfiler",
    "fgsea"
))
```

### 3. Load Required Packages

```r
library(DESeq2)
library(ggplot2)
library(RibosomeProfiling)
library(EnhancedVolcano)
library(clusterProfiler)
library(fgsea)
library(dplyr)
library(tidyr)
library(readr)
```

### Installation Time

The installation of R and all required packages typically takes approximately 15–30 minutes on a standard desktop computer with a stable internet connection. Installation time may vary depending on computer specifications, operating system, and internet speed.


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

## Demo and Reproducibility

### 1. Instructions to Run on Data

This repository provides R scripts for downstream differential expression (DE), translation efficiency (TE), and functional enrichment analyses using preprocessed RNA-seq and ribosome profiling (RPF) data.

**Input data:**
- Gene-level RNA-seq count matrix.
- Gene-level RPF count matrix.
- Sample metadata table containing sample identifiers and experimental conditions.
- Additional input files or gene lists, when required by individual analyses.

Example input files are provided in the `Input/` directory.

**Running the analysis:**

1. Clone or download the repository.
2. Install R and the required R packages (see Installation).
3. Open R or RStudio and set the working directory to the repository root.
4. Check and adjust the input file paths, sample information, and analysis parameters in the corresponding R scripts.
5. Run the relevant analysis scripts:

```r
# 1. Differential expression analysis
source("scripts/01_pipeLine_DESeq2_SingleVariableInput.R")

# 2. Translation efficiency analysis
source("scripts/02_pipeline_TE_basedOnDE_byCondtions.R")

# 3. TE visualization with custom gene labels
source("scripts/03_TE_customGeneListLabel_Plot.R")

# 4. Differential expression interaction analysis
source("scripts/04_DEs_Interaction.R")

# 5. Gene set enrichment analysis
source("scripts/05_gsea.R")
```

The scripts can be used for the corresponding analyses independently, provided that their required input files and parameters are supplied.

### 2. Expected Output

The analysis generates results and visualizations, including:

- Differential expression results from RNA-seq and RPF data.
- Translational efficiency (TE) analysis results.
- TE comparison and visualization plots.
- Differential expression interaction analysis results.
- Gene set enrichment analysis (GSEA) results and enrichment plots.
- PCA, hierarchical clustering, heatmaps, and other statistical visualizations, where applicable.

Output files and figures are generated according to the settings specified in the individual R scripts.

### 3. Expected Run Time

The downstream analysis demo typically takes approximately 10 minutes to complete on a standard desktop computer, depending on computer specifications and dataset size.

## Notes

This repository contains analysis and visualization scripts only.

Raw sequencing data, reference genomes, alignments, and large intermediate files are intentionally excluded.

Pre-processing steps performed using standard community software (e.g., FASTX, STAR, featureCounts, RSEM) should be described in the associated manuscript Methods section.

The inputs and parameters settings might need corresponding adjustments for the specific projects and applications

Example visualization: certain code demonstrates the general structure of an RNA-versus-RPF correlation plot, which means that some example code is simplified and does not represent the complete analysis implementation used in the study

## Citation

If you use this repository, please cite the associated publication.

