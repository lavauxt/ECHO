<p align="center">
  <img src="assets/logo.png" alt="ECHO logo" width="350"/>
</p>

<h1 align="center">ECHO</h1>

<p align="center">
  Exome Coverage and HMM Optimizer
</p>

<p align="center">
  <img src="https://img.shields.io/badge/License-GPLv3-blue.svg" alt="GPL-3 License"/>
</p>

# ECHO: Exome Coverage and HMM Optimizer

ECHO is a modular R package designed to streamline the ExomeDepth pipeline for Copy Number Variant (CNV) detection in exome sequencing data. It automates the end-to-end process from raw BAM files to finalized CNV reports and visualizations.

## Features
* Automated Pipeline: End-to-end execution of BAM coverage, QC, calling, and plotting.
* Parallelized: Multi-core support for high-throughput CNV calling.
* YAML-Driven: Centralized configuration management for reproducible analysis.
* Quality Control: Built-in metrics to detect low-coverage samples or exons.
* Visualization: Automated PDF generation for every detected CNV.

## Prerequisites
Before running the pipeline, ensure you have:
* Reference Genome: A FASTA file.
* Target sorted BED: A 4-column tab-delimited file (chrom, start, end, gene).
* BAMs: Coordinate-sorted BAM files with corresponding .bai index files.

## Installation
You can install the development version from GitHub:

# install.packages("remotes")
remotes::install_github("lavauxt/ECHO")

## Quick Start

ECHO uses a centralized config.yaml file to manage all pipeline parameters.

### 1. Create a config.yaml
Create a configuration file to define your input data and settings:

# config.yaml
```r
input:
  bamdir: "./data/bams"       # Directory containing .bam files
  bed: "./targets.bed"        # Path to BED file
  fasta: "./reference.fasta"  # Path to reference FASTA
  rbams: null                 # Optional: TSV with reference BAMS

output:
  dir: "./ECHO_Results"       # Output directory
  prefix: "cohort_01"         # Prefix for output files

settings:
  modechrom: "A"              # 'A' (Autosomal), 'XX', or 'XY'
  cores: 4                    # Number of cores for parallel processing
  min_corr: 0.98              # Min correlation for QC
  min_cov: 100                # Min median read depth for QC
  transition_probability: 0.0001
```

### 2. Run the Pipeline
Execute the full analysis with a single command:
```r
library(ECHO)
run_echo_pipeline("config.yaml")
```
## Pipeline Workflow

* Coverage (run_bam_coverage): Scans BAMs and extracts read counts for target regions.
* Metrics (run_metrics): Validates data quality; identifies samples/exons failing correlation or coverage thresholds.
* CNV Calling (run_cnv_calling): Performs HMM-based CNV detection using parallel processing.
* Visualization (generate_plots): Generates PDF visualizations of the CNV calls within the context of the observed versus expected read ratios.

## Configuration Reference

| Section | Parameter | Description |
| :--- | :--- | :--- |
| input | bamdir | Path to directory containing BAM files. |
| input | bed | Path to target BED file. |
| input | fasta | Path to reference FASTA. |
| settings | cores | Number of CPU cores to utilize. |
| settings | min_corr | Threshold for sample-wise correlation check. |
| settings | min_cov | Median read depth threshold for QC. |
| settings | transition_probability | HMM sensitivity parameter for ExomeDepth. |


## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).