<p align="center">
  <img src="assets/logo.png" alt="ECHO logo" width="350"/>
</p>

<div align="center">

  # ECHO
  **Exome/Panel CNV Calling on an ExomeDepth Beta-Binomial Engine**

  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
  [![R Version](https://img.shields.io/badge/R-%3E%3D%204.0.0-brightgreen.svg)](https://cran.r-project.org/)
</div>

---

## Overview

**ECHO** is an R package for detecting Copy Number Variations (CNVs) in
targeted panel and whole-exome sequencing (WES) data, built on top of the
`ExomeDepth` beta-binomial model.

ECHO handles the full pipeline around that model — coverage extraction, GC
handling, per-sample/per-exon QC, reference selection (via ExomeDepth's own
`select.reference.set()`), confidence scoring, VCF export, static PDF plots,
and an interactive HTML report — so you go from BAMs + a BED file to a
scored, plotted, reported set of CNV calls with a single `echo("config.yaml")`
call.

ECHO has a sibling pipeline, **CANOPE**, which calls CNVs with a custom
Hidden Markov Model instead. See [Related project: CANOPE](#related-project-canope)
below.

## Key Features

* **YAML-Driven Pipeline:** A single `echo("config.yaml")` command handles
  BED preprocessing, BAM coverage extraction, QC, ExomeDepth CNV calling,
  confidence scoring, VCF export, and HTML reporting.
* **Flexible Reference Sequence Source:** Coverage/GC extraction works from
  a local FASTA (`fasta_source = "file"`) or directly from an installed
  Bioconductor **BSgenome** package (`fasta_source = "bsgenome"`), with the
  BSgenome-derived FASTA cached across runs.
* **Discrete, ExomeDepth-native Coverage Extraction:** counts come from
  `ExomeDepth::getBamCounts()`, which reports genuine per-target read/
  fragment counts (not a summed depth track) — the specific quantity
  ExomeDepth's beta-binomial model is calibrated against. (CANOPE offers a
  fast megadepth-based backend as an alternative to its own extraction; that
  option isn't offered here, since megadepth reports summed/averaged
  per-base depth, a different statistical quantity that would silently
  invalidate ExomeDepth's Bayes-factor and confidence-interval math — see
  [Related project: CANOPE](#related-project-canope).)
* **Pre-Calling QC Exclusion:** `sample_qc` drops outlier samples (robust
  z-score on cross-target noise) and `exon_qc` drops problematic exons (high
  cross-sample MAD, low mean coverage, or GC content outside
  `gc_extreme_filter`) *before* any sample is handed to ExomeDepth.
* **Terminal-Exon Padding:** `pad_terminal_exons` (bp, default `0` = off)
  extends the outward-facing edge of each gene's first and last exon, where
  there's no neighbouring exon to carry the signal if coverage thins out
  right at the boundary. Padding is applied "if possible" — it never
  overlaps a neighbouring target or crosses a contig boundary, clamping
  short instead. Internal exons are untouched.
* **Gene-Boundary Spacing in Plots:** `plot_gene_gap` (default `1`) inserts
  blank x-axis space between a gene's last exon and the next gene's first
  exon in every PDF and HTML-report panel, so a window that spans more
  than one gene doesn't read as one continuous, unbroken feature.
* **Off-Target/Filler Region Handling:** BED panels sometimes include
  normalization "backbone" probes that are never a real gene exon (often
  named things like `HorsROI`). Rows whose parsed gene name matches
  `off_target_pattern` (default `"^HorsROI"`; set to `null`/`""` to
  disable) are handled per `off_target_handling`: `"na"` (default — kept,
  but excluded from exon numbering/gene grouping instead of being numbered
  as if the filler label were a gene), `"remove"` (dropped entirely), or
  `"merge"` (attached to the nearest neighbouring real gene, numbered as
  one of that gene's own exons).
* **Off-Target/Filler Region Handling in Plots:** `off_target_gene_pattern`
  (regex, default `"^HorsROI"`) flags BED rows like off-target padding
  regions that carry their own unrelated numbering. Flagged rows render as
  a fixed neutral grey tile with no text label (PDF gene track) and a
  blanked axis tick (every panel, PDF and report), instead of being
  interleaved with a gene's real exon numbers with no visual distinction.
* **Confidence Scoring:** HIGH/MEDIUM/LOW tiers from correlation, reference
  count, and read-ratio, with a built-in list of pseudogene/homology-prone
  genes always scored LOW.
* **Native Visualization:** Per-call static PDF plots (coverage, log-ratio
  with a beta-binomial predictive interval, and a z-score-vs-references
  panel) plus a self-contained, interactive HTML report — see
  [Interactive HTML report](#interactive-html-report).
* **Multi-sample Interactive Report:** one HTML file for the whole run, in
  the style of the [SPARKS](https://github.com/lavauxt/SPARKS) report: tabs,
  `plotly` figures, sortable/filterable tables with CSV export, a sample
  selector, light/dark theme. No Shiny, no server, no `rmarkdown`/`pandoc`.
* **Warnings That Cannot Be Missed:** QC thresholds are drawn on every
  metric, failing samples are flagged, and the report states explicitly
  which chromosome classes (**autosomes** / **gonosomes chrX, chrY**) were
  *not* analysed by the run — see [Chromosome scope](#chromosome-scope-autosomes-vs-gonosomes).
* **Flexible BED Annotation:** `STANDARD`/`REGEN`/`NO` preprocessing modes.
  `STANDARD` parses BED names or panel files; `REGEN` annotates intervals by
  largest overlap with UCSC knownGene exons (hg19 or hg38), preferring NM_
  transcripts when present. Unmatched intervals remain as `Unknown`.

## Installation

```R
# Install devtools if you haven't already
if (!require("devtools")) install.packages("devtools")

# Install ECHO
devtools::install_github("lavauxt/ECHO")
```

### Dependencies
Core: `ExomeDepth`, `data.table`, `matrixStats`, `yaml`, `ggplot2`,
`plotly`, `jsonlite`, `htmltools` (HTML report), and the Bioconductor
packages `Biostrings`, `Rsamtools`, `GenomeInfoDb`, `GenomicRanges`,
`IRanges`, `BiocGenerics`.

Optional (only needed for specific features): `BSgenome.Hsapiens.UCSC.*`
(`fasta_source = "bsgenome"`), `GenomicFeatures` + the matching
`TxDb.Hsapiens.UCSC.hg19.knownGene` or `TxDb.Hsapiens.UCSC.hg38.knownGene`
package (BED `REGEN` mode), and `org.Hs.eg.db` for gene symbols (Entrez IDs
are used when it is unavailable), and `ggrepel` (nicer PCA labels).

The HTML report no longer needs `rmarkdown`, `knitr`, `DT`, `kableExtra`,
`tidyr`, `dplyr` or `pagedown`.

## Quick Start

### 1. Data Preparation
BAM files (indexed) and a BED file of target regions. A reference FASTA
(indexed) or an installed BSgenome package.

### 2. Create a Configuration File

```yaml
input:
  bamdir: "data/bams/"
  bed: "data/targets.bed"
  fasta: "data/reference.fa"          # or omit and use fasta_source: bsgenome
  fasta_source: "file"                  # or "bsgenome"

output:
  dir: "result"
  prefix: "ECHO"

settings:
  modechrom: "A"                         # "A", "X", or "Y"
  min_corr: 0.98
  min_cov: 100
  min_total_reads: 5000000
  max_exon_cv: 0.5
  transition_probability: 1e-4
  expected_CNV_length: 50000
  sample_qc: true
  exon_qc: true
  qc_zscore: 3
  exon_mad_quantile: 0.90
  gc_extreme_filter: [0.15, 0.85]
  min_exon_mean: 20
  pad_terminal_exons: 0                  # e.g. 10 to pad each gene's first/last exon by 10bp
  plot_gene_gap: 1                       # blank x-axis space between adjacent genes in plots
  report_plotly: "inline"                # "inline" (self-contained, +3.6 MB) or "cdn" (needs internet to draw plots)
  report_max_profiles: 500               # CNV calls whose coverage profile is embedded in the report
```

```yaml
bed_preprocess:
  off_target_pattern: "^HorsROI"         # regex matching filler/backbone intervals, not real exons
  off_target_handling: "na"              # "na" | "remove" | "merge"
  off_target_gene_pattern: "^HorsROI"    # grey out/blank-label off-target filler regions in plots

genome_version: "hg19"
bed_process: "STANDARD"                  # or "REGEN" / "NO"
```

### 3. Run the Pipeline

```R
library(ECHO)

echo(config_path = "config.yaml")
```

### 4. Visualizing Results
The interactive HTML report (`ECHO_<prefix>_report.html`) is generated by
default (`report = TRUE`). Static per-call PDF plots are generated when
`plots = TRUE` (or `settings$plots: true` in the config).

#### Interactive HTML report

One self-contained file (about 4 MB, plotly.js is embedded): open it in any
browser, offline. It is a small client-side app — tabs are built lazily,
the URL hash (`#tab=cnvs&sample=S01`) makes a view shareable, and *Print* to
PDF expands every tab.

| Tab | Content |
|---|---|
| **Overview** | Warnings, key figures, CNV calls per sample (click a bar to filter), CNV landscape (log2 ratio along the targets) and the samples needing attention. |
| **QC & warnings** | Chromosome scope of the run, targets per chromosome, median depth / total reads / best correlation per sample **against their thresholds** (failing samples in red), and the filterable table of every QC flag (severity and type filters). |
| **Samples** | One row per sample (status OK / Warning / Fail / Excluded, depth, reads, correlation, references, overdispersion φ, calls by confidence and type, `sample_table` columns), with a detail card. |
| **CNV calls** | Filterable table (confidence, type, chromosome, free text, CSV export). Click a call for its evidence: gene track, log coverage of the test sample vs its scaled references, observed/expected with the 95 % beta-binomial interval, and z-scores vs references. |
| **Cohort** | PCA of the coverage profiles (coloured by gender when a `sample_table` is given), sample–sample Spearman correlation heat-map, depth distribution of the targets. |
| **Methods** | Settings, the confidence rules as actually applied, sample table, external reference BAMs. |
| **Run** | Warnings/errors of the pipeline log, configuration used, log tail, session info. |

**Several samples.** A *Sample* selector in the header filters every tab to
one sample (or *All samples*); with `generate_report(global = FALSE)` each
sample with a call gets its own `…/<sample>/ECHO_report_<sample>.html`. Those
single-sample files only contain that sample's data (no PCA, no correlation
matrix, no log, no sample table, no config): the only other sample names that
appear are its own reference samples, which the coverage plots show.

**Warnings.** Orange/red banners at the top of the Overview, QC and CNV tabs
(and a counter in the header) report: chromosome classes that are not
analysed or absent from the BED, samples excluded from calling, samples with
low correlation / median depth / total reads, small cohorts, low-depth or
highly variable targets, calls whose background exons are badly calibrated,
and warnings found in the pipeline log. The CNV table's *Calibration* column
shows `check` for calls with a poorly calibrated reference set.

The coverage profile of a call is embedded for the first
`report_max_profiles` calls (default 500, HIGH and MEDIUM first) to keep the
file small.

`generate_report(pdf_output = TRUE)` is ignored (the knit-time PDF export is
gone); use the browser's *Print* dialog.

#### Chromosome scope (autosomes vs gonosomes)

ECHO calls CNVs on **one chromosome class per run**, selected by
`modechrom`: `"A"` autosomes (chr1–22), `"X"` chrX or `"Y"` chrY. Whatever the
choice, part of the BED is therefore *not* analysed, and a report with no call
on chrX does not mean chrX is normal. ECHO now says so explicitly, in the QC
table (`Type = "Analysis scope"`) and in the report:

* *Autosomes are NOT analysed in this run* (`modechrom` `X`/`Y`),
* *Gonosomes (chrX / chrY) are NOT analysed in this run* (`modechrom` `A`),
* *Chromosome class missing from the BED file* (`Type = "Missing Chromosome"`),
* a note about non-canonical contigs (chrM, alt/decoy), which are never analysed.

To cover the whole genome, run ECHO once per class (`modechrom: "A"`, then
`"X"`, then `"Y"` with a `sample_table`).

## Input File Formats

**BED File**
Standard 0-based BED file (`bed_zero_based: true`, the default; set it to
`false` for 1-based closed intervals). If `bed_process` is `"STANDARD"` or
`"REGEN"`, it's automatically annotated with gene/exon information first;
the preprocessed BED is always written 0-based, and ECHO takes care of the
conversion (an input declared 1-based is no longer shifted twice).

**Sample table** (optional, required for `modechrom` `X`/`Y`)
Tab-separated, with `sample_name` and `gender` columns (any other column is
shown in the report's *Samples* tab).

## Output Format (`ECHO_<prefix>_CNV_calls.tsv`)

* `CNV.ID`, `Sample`: Row identifier and the sample the call was made in.
* `Chromosome`, `Start`, `End`: Genomic coordinates of the CNV.
* `Type`: `deletion` or `duplication` (ExomeDepth's own labels).
* `Nexons`, `Start.p`, `End.p`: Number of exons spanned and their within-gene
  start/end index.
* `BF`: ExomeDepth Bayes factor for the call.
* `Reads.expected`, `Reads.observed`, `Reads.ratio`: Observed vs.
  reference-implied read counts and their ratio.
* `Gene`: Gene name(s) overlapping the call.
* `Confidence`: `HIGH`/`MEDIUM`/`LOW` (see `score_cnv_confidence()`).
* `Correlation`, `N.comp`, `Comparator.name`: Correlation with the chosen
  reference set, how many reference samples were used, and their names.
* `Genomic.ID`: `chr:start-end` convenience string.

## Additional Features

- **BED preprocessing** (`echo_process_bed.R`) — `STANDARD`/`REGEN`/`NO`
  modes, gene/exon annotation from panel files or RefSeq, auto-set to
  `REGEN` whenever a BED file is supplied and no mode is configured.
  Exon-number extraction from the BED name (e.g. `..._ex12_...`) requires
  a genuine `ex`/`exon` token boundary (not a substring of a larger word)
  and takes the *last* such token in the name rather than the first, so
  compound names like `5utr1-ex1` or `utr2-ex10-extra` resolve to the
  right exon. Off-target/filler intervals (e.g. `HorsROI` backbone probes)
  are caught by name via `off_target_pattern` before any exon numbering
  happens — see `handle_off_target_regions()` below.
- **Off-target/filler region handling** (`echo_utils.R::handle_off_target_regions()`) —
  rows whose parsed gene matches `off_target_pattern` (default
  `"^HorsROI"`) are, per `off_target_handling`: kept with `gene = NA` and
  excluded from numbering (`"na"`, default), dropped (`"remove"`), or
  attached to the nearest neighbouring real gene on the same chromosome
  and numbered as one of its exons (`"merge"`). Without this, such rows
  get numbered as if the filler label were a gene in its own right, and
  show up interleaved with the real gene's exons in plots.
- **Sample/exon QC exclusion** (`echo_utils.R`) — `sample_qc` drops outlier
  samples (robust z-score on cross-target noise) and `exon_qc` drops
  problematic exons (high cross-sample MAD, low mean coverage, or GC
  content outside `gc_extreme_filter`) before calling. Both default `TRUE`.
- **Terminal-exon padding** (`echo_utils.R::pad_gene_terminal_exons()`) —
  extends the outward edge of each gene's first/last exon by
  `pad_terminal_exons` bp (default `0`, off) before coverage extraction, so
  a thin sliver of low coverage right at a gene boundary is less likely to
  drag the whole exon's count down. Applied "if possible": clamped short of
  the requested amount rather than overlapping a neighbouring target or
  crossing a contig boundary. Internal exon boundaries are left alone.
- **Gene-boundary spacing in plots** (`echo_utils.R::compute_gene_gap_positions()`) —
  shared by both the PDF panels (`echo_plots.R`) and the HTML report
  (`echo_report.R`), so a CNV window spanning more than one gene
  shows real blank space, in every panel, exactly between a gene's last
  exon and the next gene's first exon — lines and ribbons break there
  rather than visually bridging the gap. Controlled by `plot_gene_gap`
  (default `1`; `0` disables it).
- **QC metrics table** (`echo_metrics.R`) — per-sample and per-exon flags,
  written to TSV and shown in the report. Columns: `Sample`, `Exon`, `Type`,
  `Details`, `Gene`, `Severity` (`error` / `warning` / `info`). `Type` is one
  of `Whole sample` (low correlation, median depth or total reads),
  `Whole exon` (low median depth), `Exon variability` (high CV),
  `Cohort size` (fewer than 3 samples), `Missing Chromosome` (class absent
  from the BED) or `Analysis scope` (class present but not analysed by this
  `modechrom`). `run_qc_metrics()` takes the run's `modechrom` for this.
- **Confidence scoring** (`echo_cnv_calls.R::score_cnv_confidence()`) —
  HIGH/MEDIUM/LOW based on correlation, reference count, and read ratio,
  with a built-in list of pseudogene/homology-prone genes always scored LOW.
- **PCA of coverage profiles** (`echo_utils.R::plot_coverage_pca()`).
- **VCF export** (`echo_vcf_export.R`) — multi-sample or per-sample.
- **Interactive HTML report** (`echo_report.R::generate_report()`,
  `inst/report/echo_report.{css,js}`) — see
  [Interactive HTML report](#interactive-html-report). The old R Markdown
  templates (`inst/rmd/`) are no longer used and can be deleted.
- **Background-calibration diagnostic** (`echo_utils.R::check_background_calibration()`) —
  flags a call when the *non-called* exons in its plotted window fall
  outside the modelled 95% predictive interval more often than a
  well-calibrated interval implies (a sign the real alteration may extend
  further, or that the reference match is poor for this region).

## Changes in 0.2.0 (review and fixes)

**Behaviour changes worth knowing about**

* **Coordinates.** ECHO passed 1-based starts to `ExomeDepth::getBamCounts()`,
  which expects 0-based starts and adds 1 itself: every target was counted
  one base to the right. Fixed. With `bed_zero_based: false` the preprocessed
  BED (always 0-based) was shifted a second time. Fixed. Counts, and
  therefore calls, can change slightly compared with 0.1.0.
* **Row alignment.** `getBamCounts()` re-sorts overlapping targets; counts are
  now matched back to the BED by interval, and the run stops if they cannot be.
* **`remove_terminal_only`.** A call covering *every* exon of a one- or two-exon
  gene was always "terminal-only" and silently dropped; whole-gene events are
  now kept.
* **QC.** `run_qc_metrics()` no longer drops perfectly correlated pairs
  (a duplicated sample gave `NA`); unannotated targets are labelled
  `chr:start-end` instead of `NA:NA`; new `Severity` column and scope rows.
* **Mode X.** Reference samples are intersected with the samples that
  survived QC.

**Other fixes:** coloured PCA (gender) always failed; `output$dir` given as an
absolute path; missing `settings` entries gave `NULL` arguments;
`genome_version` under `input:` was ignored; the BAM-index check only accepted
one of the two `.bai` conventions; cached FASTA appended without a readable
`.fai`; `NA`/unknown genes in calls, plots and VCF; very long gene lists
overflowing file names; every sub-function message is now written to the
pipeline log (the report's warnings are built from it).

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).
