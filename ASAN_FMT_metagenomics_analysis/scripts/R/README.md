# Downstream R analysis scripts

This folder contains the downstream R script used for microbiome and virome statistical analyses and figure generation.

## Main script

- `downstream_microbiome_virome_analysis.R`

## Analyses included

- Alpha-diversity analysis
- Beta-diversity analysis
- Bray-Curtis PCoA
- PERMANOVA
- ANOSIM
- Taxon filtering by prevalence and mean relative abundance
- MaAsLin2 differential abundance analysis
- Selected-taxon abundance plots
- S-plot using log2 fold change and CLR-Spearman correlation
- CLR-Spearman taxon-taxon network analysis
- CLR-Spearman taxon-clinical correlation heatmap
- BH-FDR q-value calculation for correlation analyses

## Required input files

The script expects feature tables and metadata files with the following general format:

- Feature table: taxa/features as rows and samples as columns, or samples as rows and taxa/features as columns
- Metadata table: one row per sample, including `SampleID` and `group`

Example file names used in the script:

- `otu_table.csv`
- `otu_table_longitudinal.csv`
- `otu_table_endpoint.csv`
- `metadata.csv`
- `metadata_longitudinal.csv`
- `metadata_endpoint.csv`

The file names can be modified in the `CONFIG` section of the R script.

## Required R packages

- dplyr
- tidyr
- readr
- tibble
- ggplot2
- vegan
- ape
- Hmisc
- igraph
- pheatmap
- grid
- gridExtra
- MaAsLin2
- officer
- rvg
- ggrepel
- ragg

Optional packages:

- compositions
- devEMF

## Notes

Raw sequencing data are not included in this repository. Raw shotgun metagenomic sequencing data are available through the NCBI Sequence Read Archive under the BioProject accession described in the manuscript.

This public script is a generalized version of the downstream analysis workflow and does not include local computer paths.
