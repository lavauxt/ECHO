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

  if (global) {
    out_file <- file.path(output_dir, if (!is.null(prefix) && nzchar(prefix))
      paste0("ECHO_", prefix, "_report.html") else "ECHO_global_report.html")
    payload <- build_report_payload(rd, settings, config, log_file, prefix, focus = NULL)
    report_html(payload, settings, out_file)
    return(invisible(out_file))
  }

  samples <- sample_name %||% unique(as.character(rd$cnv_calls$Sample))
  if (!length(samples)) { message("[INFO] No CNV calls -- skipping per-sample reports."); return(invisible(NULL)) }
  out_files <- character()
  for (s in samples) {
    safe_s <- paste0(substr(sanitize_filename(s), 1L, 80L), "_", substr(sprintf("%08x", as.integer(sum(utf8ToInt(enc2utf8(s)) * seq_along(utf8ToInt(enc2utf8(s)))) %% 2147483647)), 1L, 8L))
    sample_dir <- file.path(output_dir, safe_s)
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(sample_dir, paste0("ECHO_report_", safe_s, ".html"))
    payload <- build_report_payload(rd, settings, config, log_file, prefix, focus = s)
    report_html(payload, settings, out_file)
    out_files <- c(out_files, out_file)
  }
  invisible(out_files)
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
  if (!is.null(coverage_rdata) && file.exists(coverage_rdata)) {
    cov <- tryCatch(load_rdata(coverage_rdata, required = "sample_names"), error = function(e) NULL)
    if (length(cov$sample_names)) all_samples <- unique(c(as.character(cov$sample_names), analysed))
    full_chr <- cov$bed_file$chromosome %||% cov$bed_file$Chr
    if (length(full_chr)) all_chromosomes <- as.character(full_chr)
  }
  all_samples <- all_samples[all_samples %in% colnames(env$counts)]
  analysed <- analysed[analysed %in% all_samples]

  qc <- NULL
  if (!is.null(qc_metrics_file) && file.exists(qc_metrics_file) && file.size(qc_metrics_file) > 0) {
    qc <- tryCatch(data.table::fread(qc_metrics_file, header = TRUE, fill = TRUE, data.table = FALSE, colClasses = "character"),
                   error = function(e) NULL)
  }
  if (is.null(qc) || !all(c("Sample", "Type", "Details") %in% colnames(qc))) {
    qc <- NULL
  } else {
    for (cn in c("Exon", "Gene")) if (!cn %in% colnames(qc)) qc[[cn]] <- ""
    if (!"Severity" %in% colnames(qc)) qc$Severity <- ifelse(qc$Type %in% c("Whole exon", "Exon variability"), "info", "warning")
    qc <- qc[, c("Sample", "Exon", "Type", "Details", "Gene", "Severity"), drop = FALSE]
    qc[is.na(qc)] <- ""
  }

  sdf <- NULL
  if (!is.null(sample_table) && file.exists(sample_table)) {
    sdf <- tryCatch(utils::read.table(sample_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = ""),
                    error = function(e) NULL)
    if (is.null(sdf) || !"sample_name" %in% colnames(sdf)) sdf <- NULL
  }
  gender_map <- env$gender_map
  if (is.null(sdf) && !is.null(gender_map) && length(gender_map))
    sdf <- data.frame(sample_name = names(gender_map), gender = as.character(gender_map), stringsAsFactors = FALSE)

  rb <- NULL
  if (!is.null(ref_bams) && is.character(ref_bams) && file.exists(ref_bams)) {
    rb <- tryCatch(utils::read.table(ref_bams, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = ""),
                   error = function(e) NULL)
    if (is.null(rb) || !"bam" %in% colnames(rb)) rb <- NULL
  }

  list(cnv_calls = cnv_calls, counts = env$counts, bed_file = bed_file,
       models = env$models, refs = env$refs, analysed = analysed, all_samples = all_samples,
       all_chromosomes = all_chromosomes, qc = qc, sample_df = sdf, ref_bams = rb)
}

# ----------------------------------------------------------------------------
# Payload builder -- produces the UNIFIED schema expected by report_assets.R
# ----------------------------------------------------------------------------

rep_arr <- function(x) I(unname(x))

build_report_payload <- function(rd, settings, config, log_file, prefix, focus = NULL) {
  modechrom <- toupper(as.character(settings$modechrom %||% "A")[1])
  if (!modechrom %in% c("A", "X", "Y")) modechrom <- "A"
  counts <- rd$counts; bed_file <- rd$bed_file; calls <- rd$cnv_calls
  has_calls <- is.data.frame(calls) && nrow(calls) > 0
  all_samples <- rd$all_samples; analysed <- rd$analysed

  thresholds <- list(
    qc_min_cov = as.numeric(settings$min_cov %||% 100),
    min_cov = as.numeric(settings$min_cov %||% 100),
    qc_min_total_reads = as.numeric(settings$min_total_reads %||% 3e5),
    min_total_reads = as.numeric(settings$min_total_reads %||% 3e5),
    qc_min_corr = as.numeric(settings$min_corr %||% 0.98),
    min_corr = as.numeric(settings$min_corr %||% 0.98),
    qc_max_exon_cv = as.numeric(settings$max_exon_cv %||% 0.5),
    max_exon_cv = as.numeric(settings$max_exon_cv %||% 0.5),
    # qc_zscore drives the sample_qc outlier filter in run_cnv_calling()
    # (detect_outlier_samples(): robust z of per-sample MAD of log2 coverage,
    #  default threshold 3).  Surface it in the thresholds block so the
    #  "Sample variability" metric bar can draw the dashed cut-off line and
    #  the exclusion reason can cite the exact value a sample exceeded.
    qc_zscore = as.numeric(settings$qc_zscore %||% 3)
  )

  # ---- chromosome scope (ECHO) ----
  scope <- chromosome_scope(rd$all_chromosomes %||% bed_file$chromosome, modechrom)
  scope_rows <- chromosome_scope_warnings(scope, modechrom)
  qc <- rd$qc
  if (is.null(qc)) qc <- data.frame(Sample = character(), Exon = character(), Type = character(), Details = character(), Gene = character(), Severity = character(), stringsAsFactors = FALSE)
  if (nrow(scope_rows)) {
    have <- paste(qc$Type, qc$Details)
    add <- scope_rows[!paste(scope_rows$Type, scope_rows$Details) %in% have, , drop = FALSE]
    if (nrow(add)) qc <- rbind(qc, data.frame(Sample = "All", Exon = "All", Type = add$Type, Details = add$Details, Gene = "All", Severity = add$Severity, stringsAsFactors = FALSE))
  }
  qc$Severity[!qc$Severity %in% c("error", "warning", "info", "fail")] <- "warning"

  # ---- pipeline log (read early so samples_df can carry exclusion reasons) ----
  log_info <- read_report_log(log_file)

  # chr tables for the genome-wide map
  chr_df <- chrom_table(bed_file$chromosome, modechrom)
  chr_all_df <- chrom_table(rd$all_chromosomes %||% bed_file$chromosome, modechrom)

  # ---- per-sample metrics ----
  mat <- as.matrix(counts[, all_samples, drop = FALSE]); storage.mode(mat) <- "double"
  med_depth <- unname(col_medians(mat)); total_reads <- colSums(mat, na.rm = TRUE)
  cor_mat <- NULL; max_corr <- rep(NA_real_, length(all_samples))
  if (length(all_samples) >= 2) {
    cor_mat <- tryCatch(stats::cor(mat, method = "spearman", use = "pairwise.complete.obs"), error = function(e) NULL)
    if (!is.null(cor_mat)) { cm <- cor_mat; diag(cm) <- NA_real_; max_corr <- apply(cm, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)) }
  }

  # ---- sample_qc outlier metric (must match detect_outlier_samples() in
  #      echo_utils.R exactly, so the value shown in the report is the same
  #      value that triggered the exclusion in run_cnv_calling()) ----
  # noise_score = MAD of log2(counts + 0.5) across all targets, per sample.
  # noise_z     = robust z of noise_score across the cohort:
  #               (sample_noise - cohort_median) / cohort_MAD.
  # |noise_z| > qc_zscore (default 3) => outlier => excluded from calling.
  noise_score <- rep(NA_real_, length(all_samples))
  noise_z     <- rep(NA_real_, length(all_samples))
  if (length(all_samples) >= 2L && nrow(mat) > 0L) {
    log_mat     <- log2(mat + 0.5)
    noise_score <- apply(log_mat, 2, stats::mad, na.rm = TRUE)
    cohort_med  <- stats::median(noise_score, na.rm = TRUE)
    cohort_mad  <- stats::mad(noise_score,     na.rm = TRUE)
    if (is.finite(cohort_mad) && cohort_mad > 0)
      noise_z <- (noise_score - cohort_med) / cohort_mad
    else noise_z <- rep(0, length(noise_score))
    names(noise_score) <- all_samples
    names(noise_z)     <- all_samples
  }

  models <- rd$models %||% list(); refs <- rd$refs %||% list()
  n_refs <- vapply(all_samples, function(s) length(refs[[s]] %||% character()), integer(1))

  call_type <- if (has_calls) { tt <- tolower(as.character(calls$Type)); ifelse(grepl("^del", tt), "DEL", ifelse(grepl("^dup", tt), "DUP", toupper(tt))) } else character()
  conf_count <- function(s, lv) if (has_calls) sum(calls$Sample == s & calls$Confidence == lv, na.rm = TRUE) else 0L
  type_count <- function(s, ty) if (has_calls) sum(calls$Sample == s & call_type == ty, na.rm = TRUE) else 0L

  qc_s <- qc[qc$Sample %in% all_samples, , drop = FALSE]
  status <- vapply(all_samples, function(s) {
    sv <- qc_s$Severity[qc_s$Sample == s]
    if (!s %in% analysed) "Excluded" else if ("error" %in% sv || "fail" %in% sv) "Fail" else if ("warning" %in% sv) "Warning" else "OK"
  }, character(1))
  # A sample that was not analysed is a FAILURE, not a warning: the cohort
  # severity filter and the Sample QC status tiles must treat excluded
  # samples as failures.
  level <- vapply(all_samples, function(s) {
    if (!s %in% analysed) "fail"
    else { sv <- qc_s$Severity[qc_s$Sample == s]; if ("error" %in% sv || "fail" %in% sv) "fail" else if ("warning" %in% sv) "warn" else "ok" }
  }, character(1))
  # Why each not-analysed sample was excluded: recover the pipeline log line
  # verbatim, then enrich it with observable QC context (cohort ranges,
  # references=0, threshold checks) so the report explains WHY the sample
  # was flagged as an outlier instead of just paraphrasing the log message.
  excl_reasons <- recover_exclusion_reasons(all_samples, log_info)
  # Cohort context used by enrich_exclusion_reason() below. Computed once
  # here because the per-sample metrics (med_depth / total_reads / max_corr /
  # n_refs / noise_z) are already materialised up-front.  noise_z is the
  # exact metric the sample_qc step uses to flag outliers, so it gets a
  # real threshold + direction.
  cohort_ctx <- list(
    med_depth   = list(values = med_depth,   thr = thresholds$min_cov,           dir = "below", label = "median depth"),
    total_reads = list(values = total_reads, thr = thresholds$min_total_reads,   dir = "below", label = "total reads"),
    max_corr    = list(values = max_corr,    thr = thresholds$min_corr,          dir = "below", label = "max correlation"),
    n_refs      = list(values = n_refs,      thr = NA_real_,                     dir = NA,      label = "references"),
    noise_z     = list(values = noise_z,     thr = thresholds$qc_zscore,         dir = "above", label = "noise z (MAD of log2 coverage)")
  )
  reason_col <- vapply(seq_along(all_samples), function(i) {
    s <- all_samples[i]
    if (s %in% analysed) "" else enrich_exclusion_reason(s, excl_reasons[[s]], i, all_samples, cohort_ctx, thresholds)
  }, character(1))
  flags_txt <- vapply(all_samples, function(s) {
    k <- which(qc_s$Sample == s)
    txt <- if (length(k)) paste0(qc_s$Type[k], ": ", qc_s$Details[k], collapse = " | ") else ""
    if (!s %in% analysed) txt <- paste0("Excluded from CNV calling", if (nzchar(txt)) " | " else "", txt)
    txt
  }, character(1))

  samples_df <- data.frame(
    sample = all_samples, status = status, level = level, Analysed = ifelse(all_samples %in% analysed, "Yes", "No"),
    `Median depth` = med_depth, `Total reads` = total_reads, `Max correlation` = max_corr,
    References = n_refs,
    n_cnvs = vapply(all_samples, function(s) if (has_calls) sum(calls$Sample == s) else 0, numeric(1)),
    n_high = vapply(all_samples, conf_count, numeric(1), lv = "HIGH"),
    n_medium = vapply(all_samples, conf_count, numeric(1), lv = "MEDIUM"),
    n_low = vapply(all_samples, conf_count, numeric(1), lv = "LOW"),
    n_del = vapply(all_samples, type_count, numeric(1), ty = "DEL"),
    n_dup = vapply(all_samples, type_count, numeric(1), ty = "DUP"),
    median_depth = med_depth, total_reads = total_reads, max_corr = max_corr, n_refs = n_refs,
    sex_inferred = NA_character_, noise_z = noise_z, ref_cor = NA_real_, n_targets = NA_real_, n_dropped = NA_real_,
    reason = reason_col,
    check.names = FALSE, stringsAsFactors = FALSE)
  sdf <- rd$sample_df
  if (!is.null(sdf)) {
    m <- match(all_samples, as.character(sdf$sample_name))
    extra <- setdiff(colnames(sdf), "sample_name"); extra <- utils::head(extra, 8)
    for (cn in extra) { nm <- sub("^gender$", "Gender", cn); v <- sdf[[cn]][m]; samples_df[[nm]] <- if (is.numeric(v)) v else ifelse(is.na(v), "", as.character(v)) }
  }

  # ---- qc_flags: a Whole-sample / fail row for each not-analysed sample ----
  # Documents WHY the sample was excluded from calling (with the enriched
  # reason that includes QC context + cohort ranges), so the Sample QC flags
  # card and the global severity filter treat it as a failure.
  excl_samples <- all_samples[!all_samples %in% analysed]
  if (length(excl_samples)) {
    # Look up the enriched reason by sample name (reason_col is aligned with
    # all_samples, so match() gives us the right index).
    excl_idx <- match(excl_samples, all_samples)
    excl_reason_enriched <- reason_col[excl_idx]
    new_rows <- data.frame(
      Sample = excl_samples, Exon = "All", Type = "Whole sample",
      Details = sprintf("Not analysed: %s", excl_reason_enriched),
      Gene = "All", Severity = "fail", stringsAsFactors = FALSE)
    # idempotent: drop new entries whose (Sample, Type, Details) already exist in qc.
    existing <- paste(qc$Sample, qc$Type, qc$Details)
    dup <- paste(new_rows$Sample, new_rows$Type, new_rows$Details) %in% existing
    if (any(dup)) new_rows <- new_rows[!dup, , drop = FALSE]
    if (nrow(new_rows)) qc <- rbind(qc, new_rows)
  }

  # ---- CNV table (lowercase unified column names) ----
  max_profiles <- as.integer(settings$report_max_profiles %||% 500L)
  windows <- list(); calib <- rep("", if (has_calls) nrow(calls) else 0L)
  cnv_df <- NULL; windows_limited <- FALSE
  if (has_calls) {
    labels <- target_axis_labels(bed_file, counts)
    calls_h <- harmonise_chr_prefix(bed_file, calls)
    rank <- match(calls$Confidence, c("HIGH", "MEDIUM", "LOW"), nomatch = 4L)
    todo <- order(rank, calls$CNV.ID)
    if (length(todo) > max_profiles) { windows_limited <- TRUE; todo <- todo[seq_len(max_profiles)] }
    for (i in todo) {
      w <- tryCatch({
        pd <- prepare_plot_data(calls_h[i, , drop = FALSE], counts, bed_file,
                                models, refs,
                                gene_gap = settings$plot_gene_gap %||% 1,
                                labels = labels)
        if (is.null(pd)) NULL else window_from_plot_data(pd)
      }, error = function(e) NULL)
      if (!is.null(w)) {
        windows[[as.character(calls$CNV.ID[i])]] <- w$window
        calib[i] <- if (isTRUE(w$flag)) "check" else "ok"
      }
    }
    pcol <- function(nm, default = NA) if (nm %in% colnames(calls)) calls[[nm]] else rep(default, nrow(calls))
    sp <- pcol("Start.p"); ep <- pcol("End.p")
    cnv_df <- data.frame(
      id = calls$CNV.ID, sample = as.character(calls$Sample), chr = as.character(calls$Chromosome),
      start = calls$Start, end = calls$End, gene = as.character(pcol("Gene", "Unknown")),
      type = call_type, interval = as.character(pcol("Genomic.ID", "")),
      exons = ifelse(is.na(sp) | is.na(ep), "", ifelse(sp == ep, as.character(sp), paste0(sp, "-", ep))),
      nexons = pcol("Nexons"), num_targ = pcol("Nexons"),
      ratio = pcol("Reads.ratio"), bf = pcol("BF"),
      expected = pcol("Reads.expected"), observed = pcol("Reads.observed"),
      confidence = as.character(pcol("Confidence", "LOW")), correlation = pcol("Correlation"),
      n_comp = pcol("N.comp"), comparators = as.character(pcol("Comparator.name", "")),
      calibration = calib, genomic = as.character(pcol("Genomic.ID", "")),
      gidx = as.numeric(pcol("global_start", NA)),
      has_plot = sapply(calls$CNV.ID, function(id) id %in% names(windows)),
      stringsAsFactors = FALSE, check.names = FALSE)
    cnv_df$gene[is.na(cnv_df$gene) | !nzchar(cnv_df$gene)] <- "Unknown"
  }

  # ---- cohort ----
  pca <- report_pca(mat, all_samples, rd$sample_df)
  corr <- NULL
  if (!is.null(cor_mat) && length(all_samples) >= 2 && length(all_samples) <= 60) {
    ord <- tryCatch(stats::hclust(stats::as.dist(1 - cor_mat))$order, error = function(e) seq_along(all_samples))
    cm <- cor_mat[ord, ord, drop = FALSE]
    corr <- list(names = I(all_samples[ord]), z = lapply(seq_len(nrow(cm)), function(i) I(round(unname(cm[i, ]), 3))))
  }
  exon_med <- row_medians(mat); lm <- log10(pmax(exon_med[is.finite(exon_med)], 1))
  depth_hist <- NULL
  if (length(lm) > 1) { br <- seq(0, max(1, ceiling(max(lm) * 10) / 10), length.out = 41); hh <- graphics::hist(lm, breaks = br, plot = FALSE, include.lowest = TRUE); depth_hist <- list(mid = I(round(hh$mids, 3)), counts = I(hh$counts)) }

  # ---- alerts, settings, scoring ----
  # (log_info was read earlier, just after qc, so samples_df can carry the
  #  exclusion reason recovered from the log.)
  alerts <- build_report_alerts(qc, scope_rows, analysed, all_samples, cnv_df, calib, modechrom, thresholds, log_info, qc_available = !is.null(rd$qc))
  settings_kv <- report_settings_kv(settings, modechrom)
  scoring <- report_scoring(settings)

  # ---- derive CANOPE-style warnings from ECHO alerts (for unified display) ----
  warns <- list()
  for (a in alerts) {
    lvl <- switch(a$severity, error = "fail", fail = "fail", warning = "warn", warn = "warn", info = "info", "info")
    warns[[length(warns) + 1]] <- list(level = lvl, category = a$kind, text = paste0(a$title, ": ", a$text), sample = if (length(a$samples)) a$samples[1] else NULL)
  }

  # ---- derive chromosomes block from ECHO scope ----
  # status: one row per chromosome CLASS (what chromStrip() in the report JS
  # reads: class / label / in_input / analysed / reason).
  chr_status_df <- data.frame(
    class = scope$class, label = scope$label, in_input = scope$present, analysed = scope$analysed,
    reason = ifelse(scope$analysed, "", ifelse(scope$present, "not analysed in this mode", "not in BED")),
    stringsAsFactors = FALSE)
  # inventory: ONE ROW PER CHROMOSOME (columns chromosome / targets_in_bed /
  # targets_analysed / analysed), as the report's "Targets per chromosome"
  # card expects. chromosome_scope() returns one row per CLASS (autosome, X,
  # Y, other) and has NO 'chromosome' column, so the previous scope$chromosome
  # was NULL and data.frame() aborted the whole report with
  # "arguments imply differing number of rows: 0, 4". Build the inventory
  # from the per-chromosome table of the input BED (chr_all_df) instead.
  chr_inv_df <- data.frame(
    chromosome = chr_all_df$chromosome, class = chr_all_df$class,
    targets_in_bed = chr_all_df$targets,
    targets_analysed = ifelse(chr_all_df$analysed, chr_all_df$targets, 0L),
    analysed = chr_all_df$analysed, stringsAsFactors = FALSE)
  # Natural order (1-22, X, Y, M) for the inventory table and its bar chart.
  chr_inv_df <- chr_inv_df[order(factor(chr_inv_df$chromosome,
                                        levels = chrom_sort_levels(chr_inv_df$chromosome))),
                           , drop = FALSE]
  # extent: x-axis layout of the genome-wide map. Each call is plotted at its
  # gidx (global target index within the ANALYSED bed), so the bands must be
  # derived from chr_df (the analysed targets, in BED order), not from the
  # per-class scope. The report JS lays each band at [off, off + len] and
  # advances off by len + 1 per chromosome; with len = targets - 1 every band
  # therefore sits exactly over its chromosome's gidx range. unname(): the
  # columns of chrom_table() may carry names, and jsonlite would then emit an
  # object instead of the array the JS expects.
  chr_extent <- list(chr = I(unname(as.character(chr_df$chromosome))),
                     len = I(unname(pmax(as.numeric(chr_df$targets) - 1, 0))))

  scope_incomplete <- any(scope$class %in% c("autosome", "X", "Y") & !scope$analysed)
  mode_label <- switch(modechrom, X = "chrX", Y = "chrY", "Autosomes (A)")

  meta <- list(
    source = "ECHO",
    title = paste0("ECHO report", if (!is.null(prefix) && nzchar(prefix)) paste0(" - ", prefix) else "", if (!is.null(focus)) paste0(" - ", focus) else ""),
    prefix = prefix %||% "", generated = format(Sys.time(), "%Y-%m-%d %H:%M"),
    version = tryCatch(as.character(utils::packageVersion("ECHO")), error = function(e) ""),
    engine = "ExomeDepth", engine_label = "ECHO (ExomeDepth beta-binomial)",
    genome_build = (config$genome_version %||% settings$genome_version) %||% "unknown",
    modechrom = modechrom, mode_label = mode_label, modechrom_label = mode_label,
    genome_version = (config$genome_version %||% settings$genome_version) %||% "",
    legacy_workspace = FALSE,
    focus = if (is.null(focus)) NULL else focus,
    n_targets = nrow(bed_file), n_samples = length(all_samples), n_analysed = length(analysed),
    n_calls = if (has_calls) nrow(calls) else 0L,
    scope_incomplete = scope_incomplete, qc_available = !is.null(rd$qc),
    windows_limited = windows_limited, windows_limit = max_profiles)

  # ---- conf_rules from scoring ----
  conf_rules <- list(
    high_q = NA_integer_, med_q = NA_integer_,
    high_refs = scoring$high_refs %||% 10, med_refs = scoring$med_refs %||% 5,
    high_targ = 1, med_targ = 1,
    genes = I(as.character(scoring$forced_low %||% character())))

  payload <- list(
    meta = meta, thresholds = thresholds, conf_rules = conf_rules, params = as.list(settings),
    summary = list(
      samples_requested = length(all_samples), samples_called = length(analysed),
      samples_with_calls = length(unique(if (has_calls) calls$Sample else character())),
      calls = if (has_calls) nrow(calls) else 0L, targets = nrow(bed_file),
      both_chromosome_classes = if (scope_incomplete) FALSE else TRUE),
    samples = rep_table(samples_df), flags = list(), sample_attrs = list(),
    warnings = warns, alerts = alerts,
    scope = rep_table(scope), chr = rep_table(chr_df), chr_all = rep_table(chr_all_df),
    chromosomes = list(status = rep_table(chr_status_df), messages = I(character()), inventory = rep_table(chr_inv_df), depth = NULL, extent = chr_extent, sex = NULL),
    qc = NULL, qc_flags = rep_table(qc),
    corr = corr, pca = pca, depth_hist = depth_hist,
    calls = if (is.null(cnv_df)) rep_empty_table(c("id","sample","chr","start","end","gene","type","confidence","has_plot")) else rep_table(cnv_df),
    plots = windows, windows = windows,
    refs = lapply(refs[all_samples[all_samples %in% names(refs)]], function(x) I(as.character(x))),
    settings = lapply(seq_along(settings_kv), function(i) I(c(names(settings_kv)[i], settings_kv[[i]]))),
    scoring = scoring,
    sample_table = if (is.null(rd$sample_df)) NULL else rep_table(utils::head(rd$sample_df, 2000)),
    ref_bams = if (is.null(rd$ref_bams)) NULL else rep_table(utils::head(rd$ref_bams, 2000)),
    log = log_info, config_yaml = report_config_yaml(config))

  if (!is.null(focus)) payload <- restrict_payload_to_sample(payload, focus)
  payload
}

# ----------------------------------------------------------------------------
# Helpers (carried over from the original ECHO report)
# ----------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a
`%nz%` <- function(a, b) if (is.null(a) || !length(a) || (length(a) == 1 && is.na(a))) b else a

chr_scope_label <- function(modechrom) switch(modechrom, X = "chrX", Y = "chrY", "autosomes")

chrom_table <- function(chr_vec, modechrom) {
  chr_vec <- as.character(chr_vec)
  if (!length(chr_vec)) return(data.frame(chromosome = character(), targets = integer(), first = integer(), last = integer(), class = character(), analysed = logical(), stringsAsFactors = FALSE))
  ch <- unique(chr_vec); ii <- split(seq_along(chr_vec), factor(chr_vec, levels = ch))
  out <- data.frame(chromosome = ch, targets = lengths(ii), first = vapply(ii, min, integer(1)), last = vapply(ii, max, integer(1)), class = chromosome_class(ch), stringsAsFactors = FALSE, row.names = NULL)
  out$analysed <- out$class == switch(modechrom, X = "X", Y = "Y", "autosome")
  out
}

window_from_plot_data <- function(pd) {
  r <- function(x, d = 3) I(round(unname(x), d))
  list(flag = isTRUE(pd$bg_calibration$flag),
       window = list(px = r(pd$px), xlab = I(pd$xlab), gene = I(pd$gene), off = I(pd$off),
                     gene_group = I(as.integer(pd$gene_group)), affected = I(pd$affected),
                     test_log = r(pd$test_log), ref_names = I(names(pd$ref_log)), ref_cov = unname(lapply(pd$ref_log, r)),
                     ratio = r(pd$ratio), lo = r(pd$lo), hi = r(pd$hi),
                     observed = r(pd$observed, 1), expected = r(pd$expected, 1),
                     z = r(pd$z), ref_z = unname(lapply(pd$ref_z, r)),
                     subtitle = pd$ci_subtitle))
}

report_pca <- function(mat, sample_names, sample_df) {
  if (length(sample_names) < 3) return(NULL)
  tryCatch({
    lm <- log2(mat + 1)
    rv <- apply(lm, 1, stats::var, na.rm = TRUE)
    keep <- is.finite(rv) & rv > 0
    if (sum(keep) < 2) return(NULL)
    p <- stats::prcomp(t(lm[keep, , drop = FALSE]), scale. = TRUE, center = TRUE)
    ve <- summary(p)$importance[2, ] * 100
    grp <- rep("All samples", length(sample_names))
    if (!is.null(sample_df) && "gender" %in% colnames(sample_df)) {
      g <- as.character(sample_df$gender[
        match(sample_names, as.character(sample_df$sample_name))])
      grp <- ifelse(is.na(g) | !nzchar(g), "unknown", g)
    }
    list(names = I(sample_names),
         pc1 = I(round(unname(p$x[, 1]), 3)),
         pc2 = I(round(unname(p$x[, 2]), 3)),
         var1 = round(ve[1], 1), var2 = round(ve[2], 1),
         group = I(grp), n_rows = sum(keep),
         n_samples = length(sample_names))
  }, error = function(e) {
    message("[WARNING] PCA skipped: ", conditionMessage(e))
    NULL
  })
}

build_report_alerts <- function(qc, scope_rows, analysed, all_samples,
                                 cnv_df, calib, modechrom, thresholds, log_info,
                                 qc_available = TRUE) {
  out <- list()

  add <- function(kind, severity, title, text, tab = "qc", samples = NULL) {
    sm <- if (length(samples)) I(as.character(samples)) else NULL
    out[[length(out) + 1L]] <<- list(
      kind = kind, severity = severity, title = title,
      text = text, tab = tab, samples = sm)
  }

  for (i in seq_len(nrow(scope_rows))) {
    d <- scope_rows$Details[i]
    if (scope_rows$Type[i] == "Missing Chromosome") {
      title <- "Chromosome class missing from the BED file"
    } else if (grepl("^Autosomes", d)) {
      title <- "Autosomes are NOT analysed in this run"
    } else if (grepl("^Gonosomes", d)) {
      title <- "Gonosomes (chrX / chrY) are NOT analysed in this run"
    } else {
      title <- "Non-canonical contigs"
    }
    add("scope", scope_rows$Severity[i], title, d, "qc")
  }

  k <- which(qc$Type == "Cohort size")
  for (i in k) add("cohort", qc$Severity[i], "Small cohort", qc$Details[i], "qc")

  if (!qc_available) {
    add("qc", "info", "No QC metrics file", "The QC metrics TSV was not found.", "qc")
    if (length(analysed) < 3 && !length(k))
      add("cohort", "warning", "Small cohort",
          sprintf("Only %d sample(s).", length(analysed)), "qc")
  }

  excl <- setdiff(all_samples, analysed)
  if (length(excl)) {
    excl_txt <- paste(utils::head(excl, 10), collapse = ", ")
    if (length(excl) > 10) excl_txt <- paste0(excl_txt, ", ...")
    # Not-analysed samples are a failure, not a warning, so the global
    # severity filter and the header alert banner treat them as failures.
    add("sample", "fail",
        paste0(length(excl), " sample(s) excluded"),
        excl_txt, "samples", excl)
  }

  ws <- qc[qc$Type == "Whole sample" & qc$Sample %in% all_samples, , drop = FALSE]
  if (nrow(ws)) {
    kind <- sub(":.*$", "", ws$Details)
    for (kd in unique(kind)) {
      rows <- ws[kind == kd, , drop = FALSE]
      sev <- if ("error" %in% rows$Severity) "error" else "warning"
      add("sample", sev,
          paste0(kd, " in ", length(unique(rows$Sample)), " sample(s)"),
          paste(rows$Sample, collapse = ", "), "qc", unique(rows$Sample))
    }
  }

  ne <- sum(qc$Type == "Whole exon")
  if (ne) add("qc", "info", paste0(ne, " target(s) with low depth"), "", "qc")

  nv <- sum(qc$Type == "Exon variability")
  if (nv) add("qc", "info", paste0(nv, " target(s) high variability"), "", "qc")

  if (is.null(cnv_df) || nrow(cnv_df) == 0) {
    add("calls", "info", "No CNV called", "", "cnvs")
  } else {
    nchk <- sum(calib == "check")
    if (nchk)
      add("calls", "warning",
          paste0(nchk, " call(s) calibration check"), "", "cnvs",
          unique(cnv_df$sample[calib == "check"]))
  }

  nw <- length(log_info$warnings)
  if (nw) add("run", "info", paste0(nw, " log warning(s)"), "", "run")

  out
}

read_report_log <- function(log_file, tail_n = 400L) {
  si <- I(utils::capture.output(utils::sessionInfo()))
  if (is.null(log_file) || !is.character(log_file) || !file.exists(log_file))
    return(list(warnings = I(character()), tail = I(character()),
                session = si, truncated = FALSE))
  ln <- tryCatch(readLines(log_file, warn = FALSE, encoding = "UTF-8"),
                error = function(e) character())
  ln <- iconv(ln, "UTF-8", "UTF-8", sub = "?")
  w <- unique(ln[grepl("^\\[(WARNING|ERROR)\\]", ln)])
  list(warnings = I(utils::head(w, 300)),
       tail = I(utils::tail(ln, tail_n)),
       session = si,
       truncated = length(ln) > tail_n)
}

# Recover the "why" for each excluded sample from the pipeline log lines.
# Matches both shapes produced by the ECHO/CANOPE pipelines:
#   "...excluding N outlier sample(s) from calling (used neither as test nor
#    reference): SGT2501756"
#   "N sample(s) excluded: SGT2501756"
# An explicit samples$reason / samples$exclusion_reason column always wins;
# this helper is called BEFORE the samples$reason column is materialised, so
# callers can fold its output straight into samples_df.
recover_exclusion_reasons <- function(sample_names, log_info) {
  out <- stats::setNames(rep(NA_character_, length(sample_names)), sample_names)
  log_lines <- c(as.character(log_info$warnings %||% character()),
                as.character(log_info$tail     %||% character()))
  pat1 <- "excluding[[:space:]]+[0-9]+[[:space:]]+outlier sample\\(s\\)[[:space:]]+from calling[[:space:]]*\\(([^)]*)\\)[[:space:]]*:[[:space:]]*(.+)$"
  pat2 <- "([0-9]+)[[:space:]]+sample\\(s\\)[[:space:]]+excluded[[:space:]]*:[[:space:]]*(.+)$"
  for (ln in log_lines) {
    m <- regmatches(ln, regexec(pat1, ln, perl = TRUE))[[1]]
    if (length(m) >= 3) {
      why <- trimws(m[2])
      for (nm in strsplit(m[3], ",[[:space:]]*")[[1]]) {
        nm <- trimws(nm)
        if (nzchar(nm) && is.na(out[nm])) out[nm] <- sprintf("Outlier sample excluded from calling (%s)", why)
      }
      next
    }
    m <- regmatches(ln, regexec(pat2, ln, perl = TRUE))[[1]]
    if (length(m) >= 3) {
      for (nm in strsplit(m[3], ",[[:space:]]*")[[1]]) {
        nm <- trimws(nm)
        if (nzchar(nm) && is.na(out[nm])) out[nm] <- "Sample excluded from calling (sample QC failure)"
      }
    }
  }
  as.list(out)
}

# Build a SHORT tag explaining why a not-analysed sample was excluded.
# The full plain-English explanation (formula, analogy, threshold) lives in
# the Methods tab of the report (see "How outlier samples are detected" in
# report_assets_refactored_v13.R), so the sample table, cohort snapshot,
# qc_flags and header alert only need a one-line tag.
#
# Tag priority:
#   1. noise_z outlier  -> "Outlier: noise z = 3.42 (> 3)"
#   2. failed threshold -> "Low median depth: 78 (< 100)" / "High noise z: ..."
#   3. fallback         -> the bare log-line reason, or a generic message.
#
# `idx` is the 1-based position of `sample` in `all_samples`; `cohort_ctx`
# is the list built in build_report_payload() with one entry per metric:
#   list(values = <numeric vec>, thr = <num|NA>, dir = "below"|"above"|NA,
#        label = <string>)
enrich_exclusion_reason <- function(sample, base_reason, idx, all_samples,
                                    cohort_ctx, thresholds) {
  fmt_short <- function(v, thr) {
    if (is.na(v) || !is.finite(v)) return("n/a")
    # integers stay integer; ratios show 2 decimals; everything else 3 sig figs.
    if (is.finite(thr) && thr == round(thr) && v == round(v)) return(format(v, scientific = FALSE))
    sprintf("%.2f", v)
  }
  # 1. noise_z outlier (the exact metric detect_outlier_samples() uses).
  nz_ctx <- cohort_ctx$noise_z
  if (!is.null(nz_ctx)) {
    nz_val <- nz_ctx$values[idx]; nz_thr <- nz_ctx$thr
    if (is.finite(nz_val) && is.finite(nz_thr) && abs(nz_val) > nz_thr)
      return(sprintf("Outlier: noise z = %s (> %s)",
                     fmt_short(nz_val, nz_thr), fmt_short(nz_thr, nz_thr)))
  }
  # 2. Any standard per-sample threshold that failed.
  for (nm in names(cohort_ctx)) {
    ctx <- cohort_ctx[[nm]]
    val <- ctx$values[idx]; thr <- ctx$thr; dir <- ctx$dir
    if (is.na(thr) || !is.finite(thr) || is.na(dir)) next
    if (is.na(val) || !is.finite(val)) next
    if (dir == "below" && val < thr)
      return(sprintf("Low %s: %s (< %s)", ctx$label, fmt_short(val, thr), fmt_short(thr, thr)))
    if (dir == "above" && abs(val) > thr)
      return(sprintf("High %s: %s (> %s)", ctx$label, fmt_short(val, thr), fmt_short(thr, thr)))
  }
  # 3. Fallback.
  if (!is.null(base_reason) && nzchar(base_reason)) base_reason
  else "Excluded from calling by sample QC"
}

report_settings_kv <- function(settings, modechrom) {
  keep <- settings[vapply(settings, function(v) is.atomic(v) && !is.null(v), logical(1))]
  vals <- lapply(keep, function(v) paste(format(v, scientific = FALSE, trim = TRUE), collapse = ", "))
  if (is.null(keep$modechrom)) vals <- c(list(modechrom = modechrom), vals)
  vals
}

report_config_yaml <- function(config) {
  if (is.null(config)) return(NULL)
  tryCatch({
    y <- yaml::as.yaml(config)
    ln <- strsplit(y, "\n", fixed = TRUE)[[1]]
    paste(utils::head(ln, 400), collapse = "\n")
  }, error = function(e) NULL)
}

report_scoring <- function(settings) {
  f <- formals(score_cnv_confidence)
  g <- function(nm) settings[[nm]] %||% f[[nm]]
  forced <- settings$score_low_confidence_genes
  if (is.null(forced)) forced <- eval(f$score_low_confidence_genes)

  high_row <- I(c("HIGH",
    sprintf(">= %.3f", g("score_high_corr")),
    sprintf(">= %d", as.integer(g("score_high_refs"))),
    sprintf("%.2f-%.2f", g("score_high_ratio_low"), g("score_high_ratio_high"))))

  med_row <- I(c("MEDIUM",
    sprintf(">= %.3f", g("score_med_corr")),
    sprintf(">= %d", as.integer(g("score_med_refs"))),
    sprintf("%.2f-%.2f / %.2f-%.2f",
            g("score_low_ratio_low"), g("score_med_ratio_low"),
            g("score_low_ratio_high"), g("score_med_ratio_high"))))

  low_row <- I(c("LOW", "any", "any", "otherwise"))

  list(rows = list(high_row, med_row, low_row),
    forced_low = I(as.character(forced)),
    high_refs = as.integer(g("score_high_refs")), med_refs = as.integer(g("score_med_refs")),
    terminal = if (isTRUE(settings$remove_terminal_only)) "Terminal-only calls removed." else if (isTRUE(settings$penalize_terminal_only)) "Terminal-only calls downgraded to LOW." else NULL)
}

restrict_payload_to_sample <- function(p, s) {
  rows <- function(tb, keep) {
    tb$data <- lapply(tb$data, function(v) I(unname(v[keep])))
    tb$n <- sum(keep)
    tb
  }
  p$samples <- rows(p$samples, p$samples$data$sample == s)
  p$calls <- rows(p$calls, p$calls$data$sample == s)
  p$qc_flags <- rows(p$qc_flags, p$qc_flags$data$Sample %in% c(s, "All"))
  ids <- as.character(p$calls$data$id)
  p$plots <- p$plots[names(p$plots) %in% ids]
  p$windows <- p$plots
  p$refs <- p$refs[names(p$refs) == s]

  keep_alert <- function(a) {
    a$kind != "run" && (is.null(a$samples) || s %in% a$samples)
  }
  trim_alert <- function(a) {
    if (length(a$samples) > 1) {
      a$samples <- I(s)
      a$text <- "See QC and warnings tab."
    }
    a
  }
  p$alerts <- lapply(Filter(keep_alert, p$alerts), trim_alert)

  p$pca <- NULL
  p$corr <- NULL
  p$sample_table <- NULL
  p$ref_bams <- NULL
  p$config_yaml <- NULL
  p$meta$n_samples <- 1L
  p$meta$n_analysed <- 1L
  p$meta$n_calls <- p$calls$n
  p$meta$samples_requested <- 1L
  p$meta$samples_called <- if (p$calls$n > 0L) 1L else 0L
  p$meta$samples_with_calls <- if (p$calls$n > 0L) 1L else 0L
  p$meta$focus <- s
  # Summary counters are derived from the filtered payload, not copied from
  # the original cohort report.
  if (is.list(p$summary)) {
    for (nm in intersect(c("samples_requested", "samples_analysed", "samples_called",
                           "samples_with_calls", "n_samples", "n_analysed"), names(p$summary)))
      p$summary[[nm]] <- if (nm %in% c("samples_called", "samples_with_calls")) as.integer(p$calls$n > 0L) else 1L
    if ("calls" %in% names(p$summary)) p$summary$calls <- p$calls$n
    if ("n_calls" %in% names(p$summary)) p$summary$n_calls <- p$calls$n
  }
  p
}

col_medians <- function(m) apply(m, 2, function(x) stats::median(x, na.rm = TRUE))
row_medians <- function(m) apply(m, 1, function(x) stats::median(x, na.rm = TRUE))
load_rdata <- function(path, required = NULL) { e <- new.env(); load(path, envir = e); as.list(e) }
sanitize_filename <- function(s) gsub("[^A-Za-z0-9._-]", "_", as.character(s))
