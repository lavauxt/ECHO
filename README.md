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

## Recent changes (report/plot parity with CANOPE)

`plots.R`, `report.R`, `ECHO_global_report.Rmd`, and `utils.R` were brought
up to parity with a sibling pipeline (CANOPE) that shares this
reporting/plotting architecture and had already been through several
rounds of bug-fixing on exactly this code. Three bugs and two new features:

**Bugs fixed:**

1. **Report couldn't be found.** `generate_report()` looked for the
   templates at `system.file("rmarkdown/ECHO_*.Rmd", package = "ECHO")`,
   but they ship at `inst/rmd/`. `system.file()` silently returns `""` for
   a subdirectory that doesn't exist rather than erroring, so this always
   failed post-install. Fixed to look under `rmd/`. Verified by building
   and installing a real minimal ECHO package and confirming
   `system.file("rmd/ECHO_global_report.Rmd", package = "ECHO")` resolves
   correctly.
2. **Report generation failed even once found.** `rmarkdown::render()`
   does not resolve a relative `output_file` against the caller's working
   directory — it resolves it against the directory containing the
   *input* Rmd (i.e. somewhere inside the installed package). A perfectly
   valid, already-created `output_dir` would fail with "The directory
   '...' does not exist." Fixed by building `output_file` as an absolute
   path once `output_dir` is confirmed to exist. Also added
   `knit_root_dir` to the `render()` calls, since the same default
   working-directory behaviour also affects code *inside* the Rmd (the
   Pipeline Log section reads `params$log_file` behind a `file.exists()`
   guard, so this was silently swallowed rather than erroring). Verified
   end-to-end: called `generate_report()` with a relative `output_dir`
   from a working directory unrelated to both the package install
   location and the data files, and confirmed the report renders with the
   log content intact.
3. **QC warning boxes showed as literal, escaped code instead of styled
   alerts.** In a `results='asis'` chunk, raw HTML lines indented 4+
   spaces get parsed by Pandoc as an indented code block, not raw HTML —
   and the missing-chromosome warning's `cat(sprintf('\n      <div
   ...'))` had each line indented to match the surrounding R code. Fixed
   by emitting each line via its own `cat()` call starting at column 0
   (matching the already-working style used elsewhere in the same file).
   Verified against a real render: the alert boxes now appear as actual
   `<div>` elements, not escaped text in a `<pre><code>` block.

**New: Z-score panel.** Both the static PDF (`plots.R`) and the
interactive report now include a third panel per CNV call, plotting the
test sample's z-score against each individual reference sample's own
z-score, all in units of the beta-binomial model's own implied standard
deviation. This isn't a mechanical copy of CANOPE's version (CANOPE's
model is a per-target HMM with its own mean/variance, ECHO's is a
beta-binomial ratio test) — it's the same *principle* applied to ECHO's
actual model: reuse the already-fitted, genome-wide dispersion parameter
(`rho`, ExomeDepth's `phi`) that the CI ribbon is already built from,
rather than computing a second, separate, noisier variance from just the
handful of local reference samples. The variance of a Beta-Binomial(n, p,
rho) is `n·p·(1-p)·(1+(n-1)·rho)` — verified directly against 5,000
simulated draws from that exact distribution (empirical/analytic variance
ratio: 1.035; resulting z-scores: mean 0.008, SD 1.017, 5.2% exceeding the
95% threshold vs. 5% expected).

**New: background-calibration flag.** Also ported from CANOPE: for each
call, checks whether the fraction of *non-called* exons in the plotted
window that fall outside the 95% interval is statistically higher than
the ~5% a well-calibrated interval implies (`check_background_calibration()`
in `utils.R`, a one-sided binomial test against a 5% null). When flagged,
the PDF's ratio-panel subtitle and the report get a note explaining what
it can mean (signal beyond the called boundary, a weak reference match, or
a technical/batch effect) — it's a prompt for manual review, not a
correction; it doesn't change the call, interval, or confidence score.

**Testing note:** the z-score/variance formula was verified analytically
(above) and the alert-box and calibration-flag *display* mechanisms were
each verified directly against real renders (a forced-`TRUE` calibration
flag renders as a proper blockquote, not escaped text — same underlying
Pandoc mechanism as bug #3). A full synthetic `generate_plots()` /
`generate_report()` round-trip was run end-to-end (10 samples, an injected
deletion, a real installed ECHO package) and confirmed the report and PDF
both render with the new z-score panel; the calibration flag was not
observed firing naturally in that specific synthetic dataset, since a
quick synthetic reference panel doesn't easily reproduce the kind of
inter-sample technical variability that trips it in real data — this
does not affect the confidence of the two direct checks described above.
Actual BAM/ExomeDepth calling wasn't exercised (no Bioconductor network
access in the environment this was written in) — the fixes and new
features sit entirely in the plotting/reporting layer, downstream of
`run_cnv_calling()`, and were tested against a hand-built summary RData
with the same structure `run_cnv_calling()` produces.

*(Minor housekeeping note: `ECHO_global_report.Rmd` had Windows-style
CRLF line endings; the editing process used here normalises edited
regions to LF, so a diff against the previous version will show line-ending
changes beyond the functional edits above. Content is otherwise
unchanged outside the sections described.)*

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
| `Plots/{sample}/ECHO_{prefix}_{sample}_{gene}_{n}.pdf` | Per-CNV PDF plots when `plots = TRUE` (coverage, gene location, ratio/CI, and z-score panels). |
| `ECHO_{prefix}_CNVs.vcf` | Combined multi-sample VCF when VCF export is enabled and no custom `vcf_output` path is supplied. |
| `Plots/{sample}/ECHO_{prefix}_{sample}.vcf` | Per-sample VCF when `vcf_per_sample = TRUE`. |

Legacy non-global HTML reporting is still supported internally: sample-specific HTML reports are written to `{output.dir}/{sample}/ECHO_report_{sample}.html` when `generate_report(..., global = FALSE)` is used.

## Report contents

The global HTML report includes a cohort summary, QC overview, per-sample CNV tables, interactive CNV plots, a collapsible explanation of confidence-score rules, and collapsible warnings, full log, and session-info sections.

Per-sample CNV tables currently display `Sample`, `Chr`, `Gene`, `Start`, `End`, `Type`, `Number of exons`, `Fold change`, and `Confidence`. Affected-exon hover text in the interactive plots includes the confidence label.

Each CNV call gets three interactive plots: **read ratio** (Observed/Expected with a 95% beta-binomial predictive interval), **coverage** (log-scale, test sample vs. each individual reference sample), and **z-score vs. references** (test sample and each individual reference sample's deviation from the model's expected value, in units of the model's own implied standard deviation — a real CNV should stand out clearly from a tight cluster of near-zero reference lines). The static PDF plots (`generate_plots()`) show the same three panels plus a gene-location tile.

When a call's non-called ("background") exons deviate from the 95% interval more often than chance would suggest, both the PDF and the report show a **background-calibration flag** — a note (not a correction) that the real alteration may extend beyond the called boundary, the reference set may be a weak match for that sample/region, or there may be a technical/batch effect at that locus. It doesn't change the call, interval, or confidence score; it's a prompt to take a closer look.

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
