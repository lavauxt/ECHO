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
* Multi-sample VCF export for cohort-level CNV calls, with optional per-sample VCF files.
* BED file preprocessing with gene-name extraction, exon numbering, and panel-restriction support.
* Sex chromosome modes for chrX and chrY calling with gender-aware reference selection.

## Prerequisites

Before running the pipeline, prepare:

* A reference FASTA file with a `.fai` index (created automatically via Rsamtools if missing).
* A target BED file with at least 4 tab-delimited columns: chromosome, start, end, gene. A 5th column for exon numbers is optional.
* Coordinate-sorted BAM files with corresponding `.bai` index files (e.g. `sample.bam` + `sample.bam.bai` or `sample.bai`).

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
  min_corr: 0.98
  min_cov: 100
  min_total_reads: 300000
  max_exon_cv: 0.5
  transition_probability: 0.0001
  pdf_output: false
```

Run the pipeline:

```r
library(ECHO)
echo(config_path = "config.yaml")
```

To write a combined VCF, pass an explicit output path:

```r
echo(config_path = "config.yaml", vcf_output = "./ECHO_Results/cohort_01.vcf")
```

## Pipeline workflow

The main pipeline performs these steps:

1. Optional BED preprocessing with `process_bed_file()`.
2. Coverage extraction with `run_bam_coverage()`.
3. QC metrics with `run_qc_metrics()`.
4. CNV calling with `run_cnv_calling()`.
5. CNV plot generation with `generate_plots()` (disabled by default; enable with `plots = TRUE`).
6. Interactive HTML reporting with `generate_report()`.
7. Optional multi-sample VCF export with `export_cnvs_to_vcf()`.

## Outputs

All output files are written under `output.dir` and prefixed with `ECHO_{prefix}_`:

| Output | Description |
| :--- | :--- |
| `ECHO_{prefix}_coverage.Rdata` | Coverage-count data. |
| `ECHO_{prefix}_QC_metrics.tsv` | QC flags and metrics. |
| `ECHO_{prefix}_CNV_calls.tsv` | CNV calls with confidence levels. |
| `ECHO_{prefix}_summary.RData` | Summary object used by plots, reports, and VCF export. |
| `ECHO_{prefix}_report.html` | Interactive cohort-level HTML report. |
| `ECHO_{prefix}_targets.bed` | Preprocessed BED file (only when `bed_process != "NO"`). |
| `Plots/{sample}/ECHO_{prefix}_{sample}_{gene}_{n}.pdf` | Per-sample CNV plot PDFs (when `plots = TRUE`). |
| `{sample}/ECHO_report_{sample}.html` | Per-sample HTML report (legacy per-sample mode). |
| `ECHO_{prefix}_CNVs.vcf` | Default multi-sample VCF (when `vcf_output` is not set to `NULL`). |

## Report contents

The global HTML report includes:

* A cohort summary and QC overview.
* Per-sample CNV tables with `Sample`, `Chr`, `Gene`, `Start`, `End`, `Type`, `Number of exons`, `Fold change`, and `Confidence` columns.
* Dynamic CNV plots; affected exon hover text includes confidence score.
* A collapsible explanation of confidence-score thresholds.
* Collapsible pipeline warnings, full log, and session information sections.

Set `settings.pdf_output: true` to generate a PDF copy of the HTML report when `pagedown` is available.

## `echo()` function reference

```r
echo(
  config_path         = NULL,    # path to config.yaml
  vcf_output          = NULL,    # VCF output path; NULL disables VCF export
  save_ed_objects     = FALSE,   # save full ExomeDepth objects per sample
  report              = TRUE,    # generate HTML report
  plots               = FALSE,   # generate per-CNV PDF plots
  vcf_per_sample      = FALSE,   # write one VCF per sample in addition to combined
  sample_name_delim   = "\\.",   # regex delimiter to split BAM filenames
  sample_name_keep    = "1",     # which parts to keep, e.g. "1", "1-2", "2-3"
  sample_name_collapse = NULL,   # separator for rejoining parts (default: first char of delim)
  custom_sample_names = NULL,    # override all sample names (length must match BAMs)
  log_file            = NULL,    # log file path; default inside output dir
  gene_field_index    = 1,       # field index for gene name after splitting (legacy)
  sample_table        = NULL,    # path to TSV with sample_name and gender columns (required for modes X/Y)
  ref_bams            = NULL,    # path to TSV with a `bam` column of external reference BAMs
  panel_files         = NULL,    # vector of BED paths or a file listing BED paths for target restriction
  ...                            # named overrides when not using config_path
)
```

`vcf_output = NULL` disables VCF export entirely. When `vcf_output` is not provided but a `config_path` is used, the VCF is written to `ECHO_{prefix}_CNVs.vcf` inside the output directory by default.

## VCF export

When `vcf_output` is provided, ECHO writes one multi-sample VCF. Each record represents a distinct CNV event and each sample column contains genotype-style fields:

* `GT`: genotype (`0/0` for no event, `0/1` for event detected)
* `CN`: estimated copy-number state (`1` for DEL, `3` for DUP, `2` for reference) — this is an approximation
* `FR`: observed/expected fold change rounded to two decimals

Set `vcf_per_sample = TRUE` to additionally write one VCF per sample under `Plots/{sample}/`.

## Sex chromosome analysis

Set `settings.modechrom` to `"X"` or `"Y"` to restrict calling to those chromosomes. A sample table is required:

```yaml
settings:
  modechrom: "X"
```

```r
echo(config_path = "config.yaml", sample_table = "./samples.tsv")
```

The sample table must be a tab-delimited file with columns `sample_name` and `gender` (`M`/`F` or `male`/`female`). In mode `X`, references are restricted to samples of the same sex. In mode `Y`, only male samples are processed.

## BED preprocessing

The `bed_preprocess` config section controls how the input BED file is cleaned and annotated before coverage extraction. Set `bed_process` to `"NO"` to skip preprocessing entirely.

```yaml
bed_preprocess:
  bed_process: "STANDARD"      # "STANDARD", "REGEN", or "NO"
  exon_sep: "[_()]"            # regex splitting column-4 name into tokens
  gene_name_keep: "3"          # 1-based token index to extract as gene name (supports "1", "1-2", "1,3")
  gene_name_collapse: "_"      # separator when rejoining multiple kept tokens
  auto_exon_number: false      # if true, assign sequential exon numbers by genomic order
  region_numbering_mode: "file_order"  # "file_order" or "bed_text"
  customexon: true             # add a Custom.Exon column to the output BED
  bed_zero_based: true         # set true when input BED uses 0-based start coordinates
  skip_invalid_intervals: true # remove intervals that cannot be extracted from the FASTA
```

**STANDARD** mode extracts gene names from column 4 using `exon_sep`. When `panel_files`
is provided, targets are restricted to regions overlapping those panel BED files.

**REGEN** mode re-annotates targets against the UCSC knownGene transcript database. It
requires `genome_version: "hg19"` or `"hg38"` and the corresponding Bioconductor
annotation packages (`TxDb.Hsapiens.UCSC.hg19.knownGene` or
`TxDb.Hsapiens.UCSC.hg38.knownGene`, plus `org.Hs.eg.db`). Only `NM_` transcripts are
retained.

## Configuration reference

### `input` section

| Parameter | Description |
| :--- | :--- |
| `bamdir` | Directory containing BAM and BAI files. |
| `bamfiles` | Optional explicit file listing BAM paths (one per line). |
| `bed` | Target BED file (≥ 4 columns: chr, start, end, gene). |
| `fasta` | Reference FASTA file (`.fai` index created automatically if absent). |
| `rbams` | Optional TSV of external reference BAMs with a `bam` column. |

### `output` section

| Parameter | Description |
| :--- | :--- |
| `dir` | Output directory (created if missing). |
| `prefix` | Prefix for all generated file names. |

### `settings` section

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `modechrom` | `"A"` | Chromosome scope: `"A"` (autosomes), `"X"`, or `"Y"`. |
| `min_corr` | `0.98` | Minimum pairwise correlation to pass QC. |
| `min_cov` | `100` | Minimum median depth (reads) per sample to pass QC. |
| `min_total_reads` | `300000` | Minimum total mapped reads per sample to pass QC. |
| `max_exon_cv` | `0.5` | Maximum per-exon coefficient of variation across samples. |
| `transition_probability` | `0.0001` | ExomeDepth HMM state-transition probability. |
| `expected_CNV_length` | `50000` | Expected CNV length in bp (Viterbi scaling). |
| `n_bins_reduced` | `10000` | Bins sub-sampled for reference-set selection. |
| `phi_bins` | `1` | Over-dispersion bins for ExomeDepth model. |
| `formula` | `"cbind(test, reference) ~ 1"` | Regression formula for `select.reference.set`. |
| `pdf_output` | `false` | Generate a PDF copy of the HTML report via `pagedown`. |

### Confidence scoring thresholds

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `score_high_corr` | `0.985` | Correlation floor for HIGH confidence. |
| `score_med_corr` | `0.95` | Correlation floor for MEDIUM confidence. |
| `score_high_refs` | `3` | Minimum reference samples for HIGH. |
| `score_med_refs` | `2` | Minimum reference samples for MEDIUM. |
| `score_low_ratio_low` | `0.75` | Lower bound of the borderline ratio zone (→ LOW). |
| `score_low_ratio_high` | `1.25` | Upper bound of the borderline ratio zone (→ LOW). |
| `score_high_ratio_low` | `0.70` | Ratio below this (with good metrics) → HIGH. |
| `score_high_ratio_high` | `1.30` | Ratio above this (with good metrics) → HIGH. |
| `score_med_ratio_low` | `0.60` | Ratio below `score_low_ratio_low` and above this → MEDIUM. |
| `score_med_ratio_high` | `1.40` | Ratio above `score_low_ratio_high` and below this → MEDIUM. |
| `score_low_confidence_genes` | see below | Gene symbols always assigned LOW confidence regardless of metrics. |

Default `score_low_confidence_genes`: `PMS2`, `SMN1`, `CYP2D6`, `HBA1`, `HBA2`, `STRC`, `CYP21A2`, `GBA1`, `CFTR`.

### `bed_preprocess` section

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `bed_process` | `"NO"` | Preprocessing mode: `"STANDARD"`, `"REGEN"`, or `"NO"`. |
| `exon_sep` | `"_"` | Regex pattern defining token boundaries in column-4 names. |
| `gene_name_keep` | `null` | Token index(es) to keep as gene name (`"1"`, `"1-2"`, `"1,3"`). |
| `gene_name_collapse` | `"_"` | Separator for rejoining multiple tokens. |
| `auto_exon_number` | `true` | Assign sequential exon numbers by genomic coordinate order. |
| `region_numbering_mode` | `"bed_text"` | Numbering strategy when `auto_exon_number = false`: `"bed_text"` or `"file_order"`. |
| `customexon` | `false` | Write a `Custom.Exon` column in the output BED. |
| `bed_zero_based` | `true` | Whether input BED start coordinates are 0-based. |
| `skip_invalid_intervals` | `true` | Remove intervals that fail FASTA extraction. |

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).
