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
* **Same QC, report and logging as CANOPE:** one shared code base for sample / target QC, metrics, recurrence flags, VCF, palette, logger and HTML report (see [Code shared with CANOPE](#code-shared-with-canope)).
* **Faster re-runs:** `echo(reuse_coverage = TRUE)` skips the BAM counting when BAMs, BED, FASTA and counting settings are unchanged; `echo_targets()` wraps this in an opt-in `{targets}` pipeline (see [Re-running faster with targets](#re-running-faster-with-targets-opt-in)).
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
  than one gene doesn't read as one continuous, unbroken feature. A run of
  consecutive off-target (HorsROI) intervals forms one block, spaced the
  same way.
* **Off-Target/Filler Region Handling:** BED panels sometimes include
  normalization "backbone" probes that are never a real gene exon (often
  named things like `HorsROI`). Rows whose parsed gene name matches
  `off_target_pattern` (default `"^HorsROI"`; set to `null`/`""` to
  disable) are kept by default (`keep_hors_roi: true`), labelled `HorsROI`
  in the processed BED and given **no exon number**, so they never take or
  shift a real gene's exon numbers. `keep_hors_roi: false` (or
  `off_target_handling: "remove"` in `handle_off_target_regions()`) drops
  them instead.
* **Off-Target Regions in Plots:** the same `off_target_pattern` drives the
  plots (PDF and HTML report). Consecutive HorsROI intervals form one grey
  block labelled `HorsROI` in the gene track (an isolated one gets a
  narrow block with a vertical label), separated from the neighbouring
  genes by `plot_gene_gap`, and their x-axis ticks and hover text read
  `HorsROI` instead of being blank.
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

Working from a source checkout instead: `NAMESPACE` is generated by roxygen, so run `devtools::document()` once after pulling (it picks up the new exports `echo_targets()` and the shared helpers) and then `devtools::load_all()` or `devtools::install()`.

### Dependencies
Core: `ExomeDepth`, `data.table`, `matrixStats`, `yaml`, `ggplot2`,
`plotly`, `jsonlite`, `htmltools` (HTML report), and the Bioconductor
packages `Biostrings`, `Rsamtools`, `GenomeInfoDb`, `GenomicRanges`,
`IRanges`, `BiocGenerics`.

Optional (only needed for specific features): `BSgenome.Hsapiens.UCSC.*`
(`fasta_source = "bsgenome"`), `GenomicFeatures` + the matching
`TxDb.Hsapiens.UCSC.hg19.knownGene` or `TxDb.Hsapiens.UCSC.hg38.knownGene`
package (BED `REGEN` mode), and `org.Hs.eg.db` for gene symbols (Entrez IDs
are used when it is unavailable), `ggrepel` (nicer PCA labels) and `targets` (>= 1.3.0, only for `echo_targets()`).

The HTML report no longer needs `rmarkdown`, `knitr`, `DT`, `kableExtra`,
`tidyr`, `dplyr`, `pagedown` or `vcfR`.

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
  min_total_reads: 300000
  max_exon_cv: 0.5
  transition_probability: 1e-4
  expected_CNV_length: 50000
  sample_qc: true
  exon_qc: true
  qc_zscore: 3
  outlier_sides: "upper"                 # "upper" (noisy only) | "both"
  recurrence_action: "flag"              # "flag" | "downgrade" | "remove" | "off"
  log_level: "INFO"                      # "DEBUG" adds per-sample reference details
  exon_mad_quantile: 0.90
  gc_extreme_filter: [0.15, 0.85]
  min_exon_mean: 20
  pad_terminal_exons: 0                  # e.g. 10 to pad each gene's first/last exon by 10bp
  plot_gene_gap: 1                       # blank x-axis space between adjacent genes in plots
  report_plotly: "inline"                # "inline" (self-contained, +3.6 MB) or "cdn" (needs internet to draw plots)
  report_max_profiles: 500               # CNV calls whose coverage profile is embedded in the report
```

Each target's x-axis label is the exon label from column 5 of the BED when
present (otherwise ECHO's own numbering), identically in the PDF plots and
in the HTML report.

```yaml
bed_preprocess:
  off_target_pattern: "^HorsROI"         # regex matching filler/backbone intervals (BED + plots)
  keep_hors_roi: true                    # false drops them from the processed BED

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
browser, offline. Its CSS and JavaScript are embedded in the package code
(`R/cnv_report_assets.R`, shared with CANOPE), so no asset file has to be installed. It is a small client-side app — tabs are built lazily,
the URL hash (`#tab=cnvs&sample=S01`) makes a view shareable, and *Print* to
PDF expands every tab.

| Tab | What it shows |
| :--- | :--- |
| **Overview** | Run status tiles, all red (fail) alerts, chromosome scope, calls per sample, sample health |
| **QC & warnings** | Sub-tabs *Samples*, *Chromosomes* (analysed or not, depth relative to the autosomes, sex check), *Call landscape*, *Coverage*, *Reference quality*, *Sample variability*, *Cohort PCA*. Metrics are drawn against their thresholds; every flag is listed in a filterable table (level, category, sample, free text) |
| **CNV calls** | One table of all calls with a column chooser, page size, Flags filter (homology gene, background calibration, terminal exon only, single target, sparse span, recurrent), CSV and BED export, genome-wide map, and the call viewer: coverage vs. its references, observed/expected ratio with its interval, z-scores; the panels zoom together |
| **Methods** | The parameters and thresholds actually used, the confidence rules, glossary |
| **Run** | Run information, session info, searchable pipeline log (level filter, copy, download) |

**Interactivity** (identical in the ECHO and CANOPE reports): `1`-`9` switch tab, `/` focuses the search box, `j` / `k` walk through the calls, `?` lists the shortcuts. The view (tab, sample, filters) is kept in the URL hash, so a link reopens the same view; `Copy locus` and `Open in IGV` (desktop IGV on port 60151) are on every call; plots with 800+ points switch to WebGL (`scattergl`) automatically. Colours are colour-blind safe (Okabe-Ito) and come from one palette (`cnv_palette()`) used by the PDF plots, the CSS and the JavaScript, in light and dark mode.

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

The scope is computed on the full BED (read from the coverage RData, passed
as `coverage_rdata`), not on the targets kept for calling, so a class that
is in the BED but not analysed is never reported as missing.

## Input File Formats

**BED File**
Standard 0-based BED file (`bed_zero_based: true`, the default; set it to
`false` for 1-based closed intervals). If `bed_process` is `"STANDARD"` or
`"REGEN"`, it's automatically annotated with gene/exon information first;
the preprocessed BED is always written 0-based, and ECHO takes care of the
conversion (an input declared 1-based is no longer shifted twice).

`bed_process: "NO"` keeps the names already in the BED (parsed, then
numbered); `"STANDARD"` re-annotates from panel files when given.

**Sample table** (optional, required for `modechrom` `X`/`Y`)
Tab-separated, with `sample_name` and `gender` columns (any other column is
shown in the report's *Samples* tab).

## Output Format (`ECHO_<prefix>_CNV_calls.tsv`)

* `CNV.ID`, `Sample`: Row identifier and the sample the call was made in.
* `Chromosome`, `Start`, `End`: Genomic coordinates of the CNV.
* `Type`: `deletion` or `duplication` (ExomeDepth's own labels).
* `Nexons`, `Start.p`, `End.p`: Number of exons spanned and their within-gene
  start/end index (`NA` when the call starts or ends on a HorsROI interval,
  which has no exon number).
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
  happens, labelled `HorsROI` and left unnumbered.
- **Off-target/filler region handling** (`cnv_shared_core.R::handle_off_target_regions()`) —
  standalone helper: rows whose name matches `off_target_pattern` are set
  to `gene = NA` (`"na"`), dropped (`"remove"`) or attached to the nearest
  real gene on the same chromosome (`"merge"`). Inside the pipeline,
  retained intervals always end up labelled `HorsROI` without an exon
  number; `assign_exon_numbers_per_gene()` never numbers them, and
  `pad_gene_terminal_exons()` never pads them.
- **Sample/exon QC exclusion** (`cnv_shared_core.R::cnv_apply_qc()`, shared with CANOPE) — `sample_qc` drops outlier
  samples (robust z-score on cross-target noise) and `exon_qc` drops
  problematic exons (high cross-sample MAD, low mean coverage, or GC
  content outside `gc_extreme_filter`) before calling. Both default `TRUE`.
- **Terminal-exon padding** (`cnv_shared_core.R::pad_gene_terminal_exons()`) —
  extends the outward edge of each gene's first/last exon by
  `pad_terminal_exons` bp (default `0`, off) before coverage extraction, so
  a thin sliver of low coverage right at a gene boundary is less likely to
  drag the whole exon's count down. Applied "if possible": clamped short of
  the requested amount rather than overlapping a neighbouring target or
  crossing a contig boundary. Internal exon boundaries are left alone.
- **Gene-boundary spacing in plots** (`cnv_shared_core.R::compute_gene_gap_positions()`) —
  shared by both the PDF panels (`echo_plots.R`) and the HTML report
  (`echo_report.R`), so a CNV window spanning more than one gene
  shows real blank space, in every panel, exactly between a gene's last
  exon and the next gene's first exon — lines and ribbons break there
  rather than visually bridging the gap. Controlled by `plot_gene_gap`
  (default `1`; `0` disables it). A one-target block gets a short flat
  interval band so its 95% interval stays visible.
- **Shared coverage window** (`echo_plots.R::prepare_plot_data()`) — one
  computation of a call's window (coverage, observed/expected, interval,
  z-scores, labels) feeds both the PDF panels and the HTML report, so the
  two always agree. PDF panels share the same plot area, so gene tiles sit
  exactly above their points.
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
- **PCA of coverage profiles** (`cnv_shared_core.R::plot_coverage_pca()`).
- **VCF export** (`echo_vcf_export.R`) — multi-sample or per-sample.
- **Interactive HTML report** (`echo_report.R::generate_report()`,
  front end in `cnv_report_assets.R`, payload in `cnv_shared_report.R`) — see
  [Interactive HTML report](#interactive-html-report). `inst/report/` and
  the old R Markdown templates (`inst/rmd/`) are no longer used and can be
  deleted.
- **Background-calibration diagnostic** (`cnv_shared_core.R::check_background_calibration()`) —
  flags a call when the *non-called* exons in its plotted window fall
  outside the modelled 95% predictive interval more often than a
  well-calibrated interval implies (a sign the real alteration may extend
  further, or that the reference match is poor for this region).

## Code shared with CANOPE

`R/cnv_shared_core.R`, `R/cnv_shared_report.R`, `R/cnv_shared_plots.R` and `R/cnv_report_assets.R` are **byte-identical** in ECHO and CANOPE. They hold the sample / target QC (`cnv_apply_qc()`), the QC metrics table (`cnv_qc_metrics()`), recurrence flagging, the logger (`cnv_logger()`), the VCF writer, the PCA, the colour palette, the PDF plot building blocks, and the whole HTML report (payload builder, CSS, JavaScript). Only the adapters (`echo_report.R`, `canope_report.R`) and the calling code differ, so the two pipelines report the same metrics with the same names, thresholds and colours. Edit a shared file in one package and copy it to the other; `Rscript tools/check_shared.R` (in the repository root, next to both packages) lists any file that drifted and `--fix` copies it.

## Outlier samples, noisy targets and recurrent calls

One routine, `cnv_apply_qc()`, takes these decisions for **both** pipelines (and for the QC-metrics table, so the table describes exactly what was excluded):

| Option | Default | Effect |
| :--- | :--- | :--- |
| `sample_qc` | `TRUE` | Exclude outlier samples before calling: robust z-score of each sample's noise (MAD of `log2(count + 0.5)` across targets) |
| `qc_zscore` | `3` | z above which a sample is an outlier |
| `outlier_sides` | `"upper"` | `"upper"`: only unusually *noisy* samples are excluded (a cleaner-than-average sample is the best kind of reference). `"both"` also excludes unusually clean ones |
| `exon_qc` | `TRUE` | Exclude problematic targets: cross-sample MAD above `exon_mad_quantile`, mean count below `min_exon_mean`, GC outside `gc_extreme_filter` (**`exon_qc: false` also switches the GC filter off**) |
| `exon_mad_quantile`, `min_exon_mean`, `gc_extreme_filter` | `0.90`, `20`, `[0.15, 0.85]` | Thresholds of the target filter |
| `recurrence_action` | `"flag"` | A call carried by `recurrence_min_samples` samples **and** `recurrence_min_fraction` of the cohort with `recurrence_min_overlap` reciprocal overlap is almost always a noisy target. `"flag"` marks it (`Carriers`, `Recurrent`), `"downgrade"` also forces LOW confidence, `"remove"` drops it, `"off"` does nothing |
| `recurrence_min_samples`, `recurrence_min_fraction`, `recurrence_min_overlap` | `3`, `0.25`, `0.5` | Recurrence thresholds |
| `remove_terminal_only`, `penalize_terminal_only` | `FALSE` | Drop / force LOW the calls that touch only a gene's first or last exon |

All of them are settable from the YAML config (under `settings:`) and as arguments of ``echo()` / `run_cnv_calling()``. The report's QC tab lists every excluded sample and target with its reason.

## Why ECHO and CANOPE results differ

Both pipelines are fed the same BAMs and BED and still disagree; in comparisons the causes were, in order of impact:

1. **Different target sets.** A target the reference FASTA cannot supply (contig missing from the FASTA, interval out of range) is dropped, and the two pipelines drop different ones (ECHO at counting time, CANOPE at GC / target QC). Both now record what they dropped, judge the chromosome scope on the *input* BED and raise an alert when whole chromosomes are lost.
2. **Different QC settings.** With `exon_qc: false` a noisy target is called as a single-exon DEL/DUP in many samples; the recurrence flag now exposes that pattern in both reports.
3. **Different outlier rules.** Two-sided on raw counts versus one-sided after GC correction excluded different samples. Both use `cnv_apply_qc()` now (`outlier_sides: upper`).
4. **Different counts.** `ExomeDepth::getBamCounts()` counts each fragment once; CANOPE's `summarizeOverlaps()` counts both mates (about 2x the reads), and BED processing / padding differ if the settings differ. Compare ratios, not raw counts.
5. **Different models.** ECHO: beta-binomial against a comparator set chosen by ExomeDepth. CANOPE: negative-binomial HMM with distance-decayed transitions against NNLS-weighted, correlation-ranked references; their confidence scores are computed from different evidence.

Calls both engines make are the ones to trust first; a call made by one engine only deserves a look at the call viewer (coverage of the sample against its references).

## Logging

Both pipelines write the same log through `cnv_logger()`: one line per event, `[LEVEL] YYYY-MM-DD HH:MM:SS message`, on the console and in ``ECHO_<prefix>_pipeline.log` in the output directory (`log_file` overrides the path)`; the file is truncated at the start of each run and messages / warnings of every called function are routed into it. Set `log_level: DEBUG` (`INFO` default, `WARNING`, `ERROR`) in the config for the per-sample reference selection details. The log is shown in the report's *Run* tab, where it can be searched, filtered by level, copied and downloaded.

## Re-running faster with {targets} (opt-in)

A `_targets.R` for ECHO tracks the config and the BAM files (size and modification time) and runs `echo(reuse_coverage = TRUE)`. Changing the calling / QC / report settings re-runs the calling and the report on the cached coverage table; changing a BAM or the BED counts again. Targets: `echo_config`, `echo_bams`, `echo_run` (all files written) and `echo_calls` (the call table).

```r
# _targets.R
library(targets)
list(echo_targets("config.yaml"))
```

then `targets::tar_make()`. `echo_run` is `echo(reuse_coverage = TRUE)`, so even when the config changes the BAM counting is repeated only if a counting setting changed (padding, sample naming, genome, BED processing) - a fingerprint stored next to the coverage RData decides.  The {targets} package is only a *suggested* dependency (>= 1.3.0); the layer's plain functions are tested without it, the DAG itself has not been exercised end-to-end.

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).
