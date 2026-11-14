# ============================================================================
# ECHO report generator (unified).
#
# Builds the unified payload from an ECHO summary RData and delegates HTML
# rendering to the shared report_html() function in report_assets.R.  The
# resulting HTML is structurally identical to the CANOPE report -- same 8 tabs,
# same CSS, same JS, same combined labels -- only the embedded JSON data differs.
#
# Drop this file + report_assets.R into the ECHO package (inst/report/ or R/).
# ============================================================================

#' Generate the interactive HTML report (unified).
#'
#' Writes a single self-contained HTML file.  Delegates rendering to the shared
#' \code{\link{report_html}} so the output is structurally identical to the
#' CANOPE report.
#'
#' @param summary_rdata Path to the ECHO_*_summary.RData.
#' @param qc_metrics_file Path to the QC metrics TSV (optional).
#' @param output_dir Output directory (created if missing).
#' @param settings List of pipeline settings.
#' @param config Original config list (optional; shown on the Run tab).
#' @param global Logical. TRUE: one report for all samples. FALSE: per-sample.
#' @param sample_name Optional sample name(s) for per-sample reports.
#' @param sample_table Optional path to the sample table TSV.
#' @param ref_bams Optional path to the external reference BAM list.
#' @param log_file Optional path to the pipeline log.
#' @param pdf_output Deprecated and ignored.
#' @param prefix Optional output prefix.
#' @param coverage_rdata Optional path to the coverage RData.
#' @return Invisibly, the path(s) of the HTML file(s) written.
#' @export
generate_report <- function(summary_rdata, qc_metrics_file, output_dir,
                            settings, config = NULL, global = TRUE, sample_name = NULL,
                            sample_table = NULL, ref_bams = NULL, log_file = NULL,
                            pdf_output = FALSE, prefix = NULL, coverage_rdata = NULL) {

  if (!file.exists(summary_rdata)) stop("[ERROR] Summary RData not found: ", summary_rdata)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) stop("[ERROR] Failed to create output directory: ", output_dir)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  settings <- settings %||% list()

  off_pattern <- (config$off_target_pattern %||% settings$off_target_pattern) %||% "^HorsROI"
  rd <- load_report_inputs(summary_rdata, qc_metrics_file, coverage_rdata, sample_table, ref_bams, off_pattern)

  report_emit(
    build = function(focus) build_report_payload(rd, settings, config, log_file, prefix, focus = focus),
    output_dir = output_dir,
    global_name = if (!is.null(prefix) && nzchar(prefix)) paste0("ECHO_", prefix, "_report.html") else "ECHO_global_report.html",
    sample_name = function(safe) paste0("ECHO_report_", safe, ".html"),
    samples = if (global) NULL else (sample_name %||% unique(as.character(rd$cnv_calls$Sample))),
    settings = settings, global = global)
}

# ----------------------------------------------------------------------------
# Inputs (unchanged from original ECHO report)
# ----------------------------------------------------------------------------

load_report_inputs <- function(summary_rdata, qc_metrics_file, coverage_rdata,
                               sample_table, ref_bams, off_target_pattern = "^HorsROI") {
  env <- new.env()
  load(summary_rdata, envir = env)
  need <- c("counts", "bed_file", "sample_names")
  miss <- setdiff(need, ls(env))
  if (length(miss)) stop("[ERROR] Summary RData lacks: ", paste(miss, collapse = ", "))

  cnv_calls <- env$cnv_calls
  if (is.null(cnv_calls)) cnv_calls <- data.frame()
  bed_file <- normalize_target_annotation(env$bed_file, off_target_pattern)

  analysed <- as.character(env$sample_names)
  all_samples <- analysed
  all_chromosomes <- bed_file$chromosome
  bed_dropped <- env$bed_dropped
  if (!is.null(coverage_rdata) && file.exists(coverage_rdata)) {
    cov <- tryCatch(load_rdata(coverage_rdata, required = "sample_names"), error = function(e) NULL)
    if (length(cov$sample_names)) all_samples <- unique(c(as.character(cov$sample_names), analysed))
    # The INPUT chromosomes: the cleaned BED plus the intervals that were dropped
    # while extracting counts (an unreadable FASTA interval, an unknown contig...).
    # Judging the scope on the cleaned BED alone hid the loss of whole chromosomes.
    full_chr <- c(cov$bed_file$chromosome %||% cov$bed_file$Chr, cov$bed_dropped$chromosome)
    if (length(full_chr)) all_chromosomes <- as.character(full_chr)
    if (is.null(bed_dropped)) bed_dropped <- cov$bed_dropped
  }
  all_samples <- all_samples[all_samples %in% colnames(env$counts)]
  analysed <- analysed[analysed %in% all_samples]

  qc <- report_read_qc(qc_metrics_file)
  sdf <- report_read_table(sample_table, "sample_name")
  gender_map <- env$gender_map
  if (is.null(sdf) && !is.null(gender_map) && length(gender_map))
    sdf <- data.frame(sample_name = names(gender_map), gender = as.character(gender_map), stringsAsFactors = FALSE)
  rb <- report_read_table(ref_bams, "bam")

  # Per-chromosome depth and sex check, on ALL input chromosomes (the coverage
  # RData holds them even when the run only analysed the autosomes) -- the same
  # helpers CANOPE uses, so the two reports show the same metric.
  cnt_all <- env$counts; chr_all <- bed_file$chromosome
  if (exists("cov", inherits = FALSE) && is.list(cov) && !is.null(cov$counts) && !is.null(cov$bed_file)) {
    ca <- cov$bed_file$chromosome %||% cov$bed_file$Chr
    if (!is.null(ca) && nrow(cov$counts) == length(ca)) { cnt_all <- cov$counts; chr_all <- as.character(ca) }
  }
  chrom_depth <- tryCatch(cnv_chrom_depth(cnt_all[, intersect(all_samples, colnames(cnt_all)), drop = FALSE], chr_all),
                          error = function(e) NULL)
  sex_check <- if (!is.null(chrom_depth)) cnv_sex_check(chrom_depth, sdf, analysed = analysed) else NULL

  list(cnv_calls = cnv_calls, counts = env$counts, bed_file = bed_file,
       models = env$models, refs = env$refs, analysed = analysed, all_samples = all_samples,
       all_chromosomes = all_chromosomes, qc = qc, sample_df = sdf, ref_bams = rb,
       qc_info = env$qc_info, bed_dropped = bed_dropped, chrom_depth = chrom_depth, sex_check = sex_check,
       modechrom = if (exists("modechrom", envir = env, inherits = FALSE)) env$modechrom else NULL)
}

# ----------------------------------------------------------------------------
# Scoring table shown on the Methods tab
# ----------------------------------------------------------------------------

report_scoring <- function(settings) {
  f <- formals(score_cnv_confidence)
  g <- function(nm) settings[[nm]] %||% eval(f[[nm]])
  forced <- settings$score_low_confidence_genes
  if (is.null(forced)) forced <- eval(f$score_low_confidence_genes)

  high_row <- I(c("HIGH",
    sprintf(">= %.3f", g("score_high_corr")),
    sprintf(">= %d", as.integer(g("score_high_refs"))),
    sprintf("%.2f-%.2f", g("score_high_ratio_low"), g("score_high_ratio_high"))))

  med_row <- I(c("MEDIUM",
    sprintf(">= %.3f", g("score_med_corr")),
    sprintf(">= %d", as.integer(g("score_med_refs"))),
    # MEDIUM ranges: [med_low, low_low] on the loss side, [low_high, med_high] on the gain side
    sprintf("%.2f-%.2f / %.2f-%.2f",
            g("score_med_ratio_low"), g("score_low_ratio_low"),
            g("score_low_ratio_high"), g("score_med_ratio_high"))))

  low_row <- I(c("LOW", "any", "any", "otherwise"))

  list(rows = list(high_row, med_row, low_row),
    forced_low = I(as.character(forced)),
    high_refs = as.integer(g("score_high_refs")), med_refs = as.integer(g("score_med_refs")),
    terminal = if (isTRUE(settings$remove_terminal_only)) "Terminal-only calls removed." else if (isTRUE(settings$penalize_terminal_only)) "Terminal-only calls downgraded to LOW." else NULL)
}

# ----------------------------------------------------------------------------
# Adapter: everything ECHO-specific the shared payload builder needs
# ----------------------------------------------------------------------------

echo_report_adapter <- function(rd, settings, config, prefix = NULL) {
  qi <- rd$qc_info
  modechrom <- toupper(as.character(rd$modechrom %||% settings$modechrom %||% "A")[1])
  if (!modechrom %in% c("A", "X", "Y")) modechrom <- "A"
  counts <- rd$counts; bed_file <- rd$bed_file
  scoring <- report_scoring(settings)
  mode_label <- switch(modechrom, X = "chrX", Y = "chrY", "Autosomes (A)")
  genome <- (config$genome_version %||% settings$genome_version) %||% ""
  num <- function(nm, d) as.numeric(settings[[nm]] %||% d)

  thresholds <- list(
    qc_min_cov = num("min_cov", 100), min_cov = num("min_cov", 100),
    qc_min_total_reads = num("min_total_reads", 3e5), min_total_reads = num("min_total_reads", 3e5),
    qc_min_corr = num("min_corr", 0.98), min_corr = num("min_corr", 0.98),
    qc_max_exon_cv = num("max_exon_cv", 0.5), max_exon_cv = num("max_exon_cv", 0.5),
    # the exact rule that excludes noisy samples (cnv_apply_qc): robust z of the
    # per-sample MAD of log2 coverage, on the side(s) given by outlier_sides
    qc_zscore = as.numeric(qi$params$qc_zscore %||% settings$qc_zscore %||% 3),
    outlier_sides = as.character(qi$params$outlier_sides %||% settings$outlier_sides %||% "upper"))

  labels <- NULL; calls_h <- NULL        # built lazily, only when there are calls
  list(
    source = "ECHO", pkg = "ECHO", title = "ECHO report",
    meta = list(engine = "ExomeDepth", engine_label = "ECHO (ExomeDepth beta-binomial)",
                genome_build = if (nzchar(genome)) genome else "unknown", genome_version = genome,
                modechrom = modechrom, mode_label = mode_label, modechrom_label = mode_label),
    modechrom = modechrom,
    analysed_class = switch(modechrom, X = "X", Y = "Y", "autosome"),
    scope = function(input_chr, analysed_chr) echo_scope(input_chr, modechrom, analysed_chromosomes = analysed_chr),
    scope_hints = ECHO_SCOPE_HINTS,
    chrom_depth = rd$chrom_depth, sex_check = rd$sex_check,
    thresholds = thresholds,
    scoring = scoring,
    conf_rules = list(high_q = NA_integer_, med_q = NA_integer_,
                      high_refs = scoring$high_refs %||% 10, med_refs = scoring$med_refs %||% 5,
                      high_targ = 1, med_targ = 1, genes = I(as.character(scoring$forced_low %||% character()))),
    cols = list(sample = "Sample", type = "Type", confidence = "Confidence", id = "CNV.ID"),
    window = function(calls, i) {
      if (is.null(labels)) { labels <<- target_axis_labels(bed_file, counts); calls_h <<- harmonise_chr_prefix(bed_file, calls) }
      pd <- prepare_plot_data(calls_h[i, , drop = FALSE], counts, bed_file, rd$models, rd$refs,
                              gene_gap = settings$plot_gene_gap %||% 1, labels = labels)
      if (is.null(pd)) NULL else window_from_plot_data(pd)
    },
    call_table = function(calls, call_type, calib, has_plot) {
      pcol <- function(nm, default = NA) if (nm %in% colnames(calls)) calls[[nm]] else rep(default, nrow(calls))
      sp <- pcol("Start.p"); ep <- pcol("End.p")
      data.frame(
        id = calls$CNV.ID, sample = as.character(calls$Sample), chr = as.character(calls$Chromosome),
        start = calls$Start, end = calls$End, gene = as.character(pcol("Gene", "Unknown")),
        type = call_type, interval = as.character(pcol("Genomic.ID", "")),
        exons = ifelse(is.na(sp) | is.na(ep), "", ifelse(sp == ep, as.character(sp), paste0(sp, "-", ep))),
        nexons = pcol("Nexons"), num_targ = pcol("Nexons"),
        kb = round((as.numeric(calls$End) - as.numeric(calls$Start) + 1) / 1000, 2),
        ratio = pcol("Reads.ratio"), bf = pcol("BF"),
        expected = pcol("Reads.expected"), observed = pcol("Reads.observed"),
        confidence = as.character(pcol("Confidence", "LOW")), correlation = pcol("Correlation"),
        n_comp = pcol("N.comp"), comparators = as.character(pcol("Comparator.name", "")),
        calibration = calib, genomic = as.character(pcol("Genomic.ID", "")),
        gidx = as.numeric(pcol("global_start", NA)),
        carriers = as.integer(pcol("Carriers", NA_integer_)), recurrent = pcol("Recurrent", FALSE) %in% TRUE,
        single_target = as.integer(pcol("Nexons", NA)) == 1L,
        has_plot = has_plot, stringsAsFactors = FALSE, check.names = FALSE)
    })
}

#' Build the unified report payload of an ECHO run (see \code{report_build_payload})
#' @noRd
build_report_payload <- function(rd, settings, config, log_file, prefix, focus = NULL) {
  report_build_payload(rd, settings, config, log_file, prefix, focus,
                       ad = echo_report_adapter(rd, settings, config, prefix))
}
