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

ECHO is a modular R package that streamlines ExomeDepth-based copy-number variant (CNV) detection from exome sequencing data. It automates BAM coverage extraction, QC, CNV calling, plotting, HTML reporting, and optional VCF export.

## Features

* End-to-end BAM coverage, QC, CNV calling, plotting, and reporting.
* YAML-driven configuration for reproducible analyses.
* QC metrics for low coverage, low total reads, sample correlation, and exon variability.
* Interactive HTML report with per-sample CNV tables, dynamic plots, warnings, log, session information, and a collapsible confidence-scoring explanation.
* Per-sample output directories for individual reports and CNV plot PDFs.
* Optional PDF report generation.
* Multi-sample VCF export for cohort-level CNV calls.

## Prerequisites

Before running the pipeline, prepare:

* A reference FASTA file, with index files as required by Bioconductor/Rsamtools.
* A target BED file with at least 4 tab-delimited columns: chromosome, start, end, gene.
* Coordinate-sorted BAM files with corresponding `.bai` index files.

## Installation

```r
# install.packages("remotes")
remotes::install_github("lavauxt/ECHO")
```

## Quick start

Create a `config.yaml` file:

```yaml
input:
  bamdir: "./data/bams"
  bed: "./targets.bed"
  fasta: "./reference.fasta"
  rbams: null

output:
  dir: "./ECHO_Results"
  prefix: "cohort_01"

settings:
  modechrom: "A"
  cores: 4
  min_corr: 0.98
  min_cov: 100
  transition_probability: 0.0001
  pdf_output: false
```

Run the pipeline:

```r
library(ECHO)
echo(config_path = "config.yaml", vcf_output = "./ECHO_Results/cohort_01.vcf")
```

## Pipeline workflow

The main pipeline performs these steps:

1. Coverage extraction with `run_bam_coverage()`.
2. QC metrics with `run_qc_metrics()`.
3. CNV calling with `run_cnv_calling()`.
4. CNV plot generation with `generate_plots()`.
5. Interactive HTML reporting with `generate_report()`.
6. Optional multi-sample VCF export with `export_cnvs_to_vcf()`.

## Outputs

Typical outputs are written under `output.dir`:

| Output | Description |
| :--- | :--- |
| `ECHO_coverage.Rdata` | Coverage-count data. |
| `QC_metrics.tsv` | QC flags and metrics. |
| `CNV_calls.tsv` | CNV calls with confidence levels. |
| `ECHO_summary.RData` | Summary object used by plots, reports, and VCF export. |
| `ECHO_global_report.html` | Interactive cohort-level HTML report. |
| `Plots/<sample>/` | Per-sample CNV plot PDFs. |
| `<sample>/ECHO_report_<sample>.html` | Per-sample HTML report when non-global reporting is used. |
| `*.vcf` | Optional multi-sample VCF with all samples as genotype columns. |

## Report contents

The global HTML report includes:

* A cohort summary and QC overview.
* Per-sample CNV tables with `Sample`, `Chr`, `Gene`, `Start`, `End`, `Type`, `Number of exons`, `Fold change`, and `Confidence` columns.
* Dynamic CNV plots; affected exon hover text includes confidence score.
* A collapsible explanation of confidence-score thresholds.
* Collapsible pipeline warnings, full log, and session information sections.

Set `settings.pdf_output: true` to generate a PDF copy of the HTML report when `pagedown` is available.

## VCF export

When `vcf_output` is provided to `echo()`, ECHO writes one multi-sample VCF. Each record represents a distinct CNV event and each sample column contains genotype-style fields:

* `GT`: genotype (`0/0` for no event, `0/1` for event detected)
* `CN`: estimated copy-number state (`1` for DEL, `3` for DUP, `2` for reference)
* `FR`: observed/expected fold change rounded to two decimals

## Configuration reference

| Section | Parameter | Description |
| :--- | :--- | :--- |
| `input` | `bamdir` | Directory containing BAM files. |
| `input` | `bamfiles` | Optional explicit BAM file list. |
| `input` | `bed` | Target BED file. |
| `input` | `fasta` | Reference FASTA file. |
| `input` | `rbams` | Optional TSV of external reference BAMs. |
| `output` | `dir` | Output directory. |
| `output` | `prefix` | Prefix for generated files. |
| `settings` | `modechrom` | Chromosome mode: `A`, `X`, or `Y`. |
| `settings` | `cores` | Number of CPU cores. |
| `settings` | `min_corr` | Minimum sample correlation threshold. |
| `settings` | `min_cov` | Median depth threshold for QC. |
| `settings` | `transition_probability` | ExomeDepth HMM transition probability. |
| `settings` | `pdf_output` | Generate PDF copy of HTML report; default `false`. |

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).
