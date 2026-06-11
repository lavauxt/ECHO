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

ECHO is an R package for ExomeDepth-based copy-number variant (CNV) detection from exome sequencing data. It covers BAM coverage extraction, BED preprocessing, QC, CNV calling, PDF plotting, HTML reporting, and VCF export.

## Features

ECHO provides an end-to-end CNV workflow driven by YAML configuration. It includes QC metrics for low coverage, low total reads, sample correlation, and exon variability; BED preprocessing with gene-name extraction and exon numbering; optional CNV PDF plots; an interactive HTML report with per-sample tables, dynamic plots, warnings, log, session information, and confidence-scoring rules; and VCF export for cohort-level CNV calls, with optional per-sample VCFs.

## Prerequisites

Before running the pipeline, prepare:

* A reference FASTA file with a `.fai` index. If the index is missing, ECHO tries to create it with `Rsamtools`.
* A target BED file with at least 4 tab-delimited columns: chromosome, start, end, gene/name. A 5th exon-number column is optional.
* Coordinate-sorted BAM files with corresponding index files.

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

To export a combined VCF to an explicit path:

```r
echo(config_path = "config.yaml", vcf_output = "./ECHO_Results/cohort_01.vcf")
```

## Pipeline workflow

The main pipeline performs these steps:

1. Optional BED preprocessing with `process_bed_file()`.
2. Coverage extraction with `run_bam_coverage()`.
3. QC metrics with `run_qc_metrics()`.
4. CNV calling with `run_cnv_calling()`.
5. Optional CNV PDF plotting with `generate_plots()`.
6. Optional HTML report generation with `generate_report()`.
7. Optional VCF export with `export_cnvs_to_vcf()`.

## Outputs

All outputs are written under `output.dir`.

| Output | Description |
| :--- | :--- |
| `ECHO_{prefix}_pipeline.log` | Pipeline log with warnings and session info header. |
| `ECHO_{prefix}_coverage.Rdata` | Coverage-count data plus BED/sample metadata. |
| `ECHO_{prefix}_QC_metrics.tsv` | QC flags and metrics. |
| `ECHO_{prefix}_CNV_calls.tsv` | CNV calls with confidence labels. |
| `ECHO_{prefix}_summary.RData` | Summary object used by plots, report, and VCF export. |
| `ECHO_{prefix}_targets.bed` | Preprocessed BED file when `bed_process != "NO"`. |
| `ECHO_{prefix}_report.html` | Global HTML report. |
| `Plots/{sample}/ECHO_{prefix}_{sample}_{gene}_{n}.pdf` | Per-CNV PDF plots when `plots = TRUE`. |
| `ECHO_{prefix}_CNVs.vcf` | Combined multi-sample VCF when VCF export is enabled and no custom `vcf_output` path is supplied. |
| `Plots/{sample}/ECHO_{prefix}_{sample}.vcf` | Per-sample VCF when `vcf_per_sample = TRUE`. |

Legacy non-global HTML reporting is still supported internally: sample-specific HTML reports are written to `{output.dir}/{sample}/ECHO_report_{sample}.html` when `generate_report(..., global = FALSE)` is used.

## Report contents

The global HTML report includes a cohort summary, QC overview, per-sample CNV tables, interactive CNV plots, a collapsible explanation of confidence-score rules, and collapsible warnings, full log, and session-info sections.

Per-sample CNV tables currently display `Sample`, `Chr`, `Gene`, `Start`, `End`, `Type`, `Number of exons`, `Fold change`, and `Confidence`. Affected-exon hover text in the interactive plots includes the confidence label.

Set `settings.pdf_output: true` to also generate a PDF copy of the HTML report when `pagedown` is available.

## `echo()` function reference

```r
echo(
  config_path          = NULL,
  vcf_output           = NULL,
  save_ed_objects      = FALSE,
  report               = TRUE,
  plots                = FALSE,
  vcf_per_sample       = FALSE,
  sample_name_delim    = "\\.",
  sample_name_keep     = "1",
  sample_name_collapse = NULL,
  custom_sample_names  = NULL,
  log_file             = NULL,
  gene_field_index     = 1,
  sample_table         = NULL,
  ref_bams             = NULL,
  panel_files          = NULL,
  ...
)
```

`report = FALSE` skips HTML report generation. `plots = TRUE` enables per-CNV PDF output. `vcf_output = FALSE` disables combined VCF export; `vcf_per_sample = TRUE` additionally writes one VCF per sample.

## VCF export

When VCF export is enabled, ECHO writes one multi-sample VCF unless `vcf_output = FALSE`. Each VCF record represents one distinct CNV event, and each sample field contains:

* `GT`: genotype (`0/0` for no event, `0/1` for event detected)
* `CN`: estimated copy-number state (`1` for DEL, `2` for normal, `3` for DUP)
* `FR`: observed/expected fold change rounded to two decimals

If `vcf_per_sample = TRUE`, additional per-sample VCFs are written under `Plots/{sample}/`.

## Sex chromosome analysis

Set `settings.modechrom` to `"X"` or `"Y"` to restrict calling to those chromosomes. A sample table is required:

```yaml
settings:
  modechrom: "X"
```

```r
echo(config_path = "config.yaml", sample_table = "./samples.tsv")
```

The sample table must be tab-delimited and contain `sample_name` and `gender` columns. In `X` mode, references are restricted to samples of the same sex. In `Y` mode, only male samples are processed.

## BED preprocessing

The `bed_preprocess` section controls how the input BED is cleaned and annotated before coverage extraction. Set `bed_process: "NO"` to skip preprocessing.

```yaml
bed_preprocess:
  bed_process: "STANDARD"      # "STANDARD", "REGEN", or "NO"
  exon_sep: "[_()]"            # regex splitting column-4 names into tokens
  gene_name_keep: "3"          # 1-based token index or range: "1", "1-2", "1,3"
  gene_name_collapse: "_"      # separator when rejoining multiple kept tokens
  auto_exon_number: true        # assign sequential exon numbers by genomic order
  region_numbering_mode: "bed_text"  # used when auto_exon_number = false
  customexon: false             # include processed exon number as BED column 5
  bed_zero_based: true          # input BED starts are 0-based
  skip_invalid_intervals: true  # drop intervals that fail FASTA extraction
```

In `STANDARD` mode, gene names are parsed from column 4. If `panel_files` is provided, targets are restricted to regions overlapping those panel BED files. In `REGEN` mode, targets are re-annotated against UCSC knownGene transcript annotations and mapped to gene symbols using `org.Hs.eg.db` when available.

When a processed BED includes a 5th column, `run_bam_coverage()` imports it as `exon_number`, and the report/plotting code uses that field preferentially for exon labels.

## Configuration reference

### `input`

| Parameter | Description |
| :--- | :--- |
| `bamdir` | Directory containing BAM files. |
| `bamfiles` | Optional file listing BAM paths, one per line. |
| `bed` | Target BED file. |
| `fasta` | Reference FASTA file. |
| `rbams` | Optional TSV of external reference BAMs with a `bam` column. |

### `output`

| Parameter | Description |
| :--- | :--- |
| `dir` | Output directory. |
| `prefix` | Prefix used in generated file names. |

### `settings`

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `modechrom` | `"A"` | Chromosome scope: autosomes (`"A"`), `"X"`, or `"Y"`. |
| `min_corr` | `0.98` | Minimum pairwise correlation to pass QC. |
| `min_cov` | `100` | Minimum median depth per sample. |
| `min_total_reads` | `300000` | Minimum total reads per sample. |
| `max_exon_cv` | `0.5` | Maximum exon-level coefficient of variation across samples. |
| `transition_probability` | `0.0001` | ExomeDepth HMM transition probability. |
| `expected_CNV_length` | `50000` | Expected CNV length in bp. |
| `n_bins_reduced` | `10000` | Number of bins sub-sampled for reference-set selection. |
| `phi_bins` | `1` | Over-dispersion bins for ExomeDepth. |
| `formula` | `"cbind(test, reference) ~ 1"` | Formula passed to `select.reference.set()`. |
| `pdf_output` | `false` | Also generate PDF from the HTML report via `pagedown`. |

### Confidence scoring

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `score_high_corr` | `0.985` | Correlation threshold for HIGH confidence. |
| `score_med_corr` | `0.95` | Correlation threshold for MEDIUM confidence. |
| `score_high_refs` | `3` | Minimum number of reference samples for HIGH. |
| `score_med_refs` | `2` | Minimum number of reference samples for MEDIUM. |
| `score_low_ratio_low` | `0.75` | Lower bound of the borderline fold-change zone. |
| `score_low_ratio_high` | `1.25` | Upper bound of the borderline fold-change zone. |
| `score_high_ratio_low` | `0.70` | Lower HIGH-confidence fold-change threshold. |
| `score_high_ratio_high` | `1.30` | Upper HIGH-confidence fold-change threshold. |
| `score_med_ratio_low` | `0.60` | Lower MEDIUM-confidence fold-change threshold. |
| `score_med_ratio_high` | `1.40` | Upper MEDIUM-confidence fold-change threshold. |
| `score_low_confidence_genes` | package default | Genes always labeled LOW confidence regardless of metrics. |

Default low-confidence genes: `PMS2`, `SMN1`, `CYP2D6`, `HBA1`, `HBA2`, `STRC`, `CYP21A2`, `GBA1`, `CFTR`.

### `bed_preprocess`

| Parameter | Default | Description |
| :--- | :--- | :--- |
| `bed_process` | `"NO"` | Preprocessing mode: `"STANDARD"`, `"REGEN"`, or `"NO"`. |
| `exon_sep` | `"_"` | Regex/pattern used to split column-4 names. |
| `gene_name_keep` | `null` | Token index(es) to retain as gene name. |
| `gene_name_collapse` | `"_"` | Separator used when rejoining kept tokens. |
| `auto_exon_number` | `true` | Assign sequential exon numbers by genomic order. |
| `region_numbering_mode` | `"bed_text"` | Numbering strategy when `auto_exon_number = false`: `"bed_text"` or `"file_order"`. |
| `customexon` | `false` | Write processed exon numbering as BED column 5. |
| `bed_zero_based` | `true` | Whether input BED start coordinates are 0-based. |
| `skip_invalid_intervals` | `true` | Remove intervals that fail FASTA extraction. |

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).
