# DLS FMT GK rat metagenomics analysis

This folder contains study-specific analysis scripts and bioinformatic workflow information used for the manuscript:

"Fecal microbiota transplantation from duodenal light stimulation–conditioned donors is associated with intestinal incretin-related remodeling and microbiome changes in diabetic Goto–Kakizaki rats."

## Raw sequencing data

The raw shotgun metagenomic sequencing data generated in this study have been deposited in the NCBI Sequence Read Archive under BioProject accession number PRJNA1478298.

## Workflow

Initial preprocessing and host-read prefiltering of shotgun metagenomic sequencing data were performed using a Nextflow workflow composed of standard nf-core modules, including fastp, Cutadapt, Bowtie2, and MultiQC.

Bacterial taxonomic profiling was performed using a MetaPhlAn-based Python wrapper script.

Viral read-level classification was performed using a custom MEGABLAST-based Python script against a custom viral reference database, followed by accession-to-taxid mapping and taxonomic aggregation.

Downstream microbiome and virome statistical analyses and figure generation were performed using R.

## Repository contents

- `workflow/`: Nextflow preprocessing run summary, preprocessing parameters, software versions, execution trace, pipeline DAG, and execution report.
- `scripts/python/`: Python scripts used for bacteriome taxonomic profiling and virome read-level assignment.
- `scripts/R/`: R script used for downstream microbiome and virome statistical analyses and figure generation.

## Analysis overview

Initial preprocessing and host-read prefiltering were performed using a Nextflow workflow composed of standard nf-core modules, including fastp, Cutadapt, Bowtie2, and MultiQC.

Bacterial taxonomic profiling was performed using a MetaPhlAn-based Python wrapper script.

Viral read-level classification was performed using a custom MEGABLAST-based Python script, followed by accession-to-taxid mapping and taxonomic aggregation.

Downstream diversity analyses, differential abundance analyses, host–microbe correlation analyses, bacteriome–virome network analyses, S-plots, heatmaps, and figure generation were performed using R.

## Notes

Raw sequencing files, large reference databases, intermediate work directories, and local server-specific paths are not included in this repository.

This repository provides study-specific preprocessing information, execution files, Python scripts, and downstream R scripts to support reproducibility of the analyses described in the manuscript.
