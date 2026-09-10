# RNA_Profile
RNAseq and RPFseq based DE/TE analysis

# RNA-seq Differential expression analysis (DE) and Translation Efficiency Analysis (TE)

This repository contains the R scripts used for differential expression
(DE), translation efficiency (TE), pathway analysis, and visualization
associated with this study.

## Analysis

The repository includes:

1. RNA-seq differential expression analysis
2. Translation efficiency analysis using RNA-seq and Ribo-seq data
3. PCA and other quality-control visualizations
4. Heatmap and volcano plot generation
5. Functional and pathway analysis
6. Export of analysis results

## Repository structure

- `scripts/` — analysis and visualization scripts
- `config/` — sample metadata required for the analysis
- `input/` — description of required input data
- `results/` — selected analysis results and figures
- `docs/` — additional methodological notes

## Data availability

Raw sequencing data and large intermediate files are not included
in this repository. Data availability is described in the associated
publication.

## Reproducibility

The R package environment used for the analysis is recorded in
`renv.lock`.

## Software

The analysis uses standard R/Bioconductor packages. Package versions
are recorded in `renv.lock`.
