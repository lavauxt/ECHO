#' Generate the interactive HTML report
#'
#' Writes a single self-contained HTML file (no server, no Shiny, no network
#' access needed once generated). The page is a small client-side application in
#' the style of the SPARKS unified report: tabs (Overview, QC & warnings,
#' Samples, CNV calls, Cohort, Methods, Run), \pkg{plotly} figures, sortable and
#' filterable tables with CSV export, a light/dark/auto theme and \code{#tab=}
#' deep links. A \emph{sample selector} in the header filters every tab, so one
#' report serves a whole multi-sample run; clicking a CNV row shows its coverage
#' evidence (log coverage vs references, observed/expected with the
#' beta-binomial interval, z-scores) computed from the same model as the PDF
#' plots.
#'
#' \strong{Warnings.} The report surfaces
#' \itemize{
#'   \item the chromosome scope of the run: chromosome classes (autosomes, chrX,
#'     chrY) that are \emph{not analysed} by the selected \code{modechrom}, and
#'     classes missing from the BED;
#'   \item sample-level QC failures (low correlation, depth, total reads),
#'     exon-level QC, cohort size, and samples excluded from calling;
#'   \item calls whose background exons are poorly calibrated.
#' }
#' Every metric is drawn against its threshold, and failing samples are
#' coloured.
#'
#' @param summary_rdata Path to the \code{ECHO_*_summary.RData} written by
#'   \code{run_cnv_calling()}.
#' @param qc_metrics_file Path to the QC metrics TSV (optional).
#' @param output_dir Output directory (created if missing).
#' @param settings List of pipeline settings (\code{cfg$settings}).
#' @param config Original config list (optional; shown on the Run tab).
#' @param global Logical. \code{TRUE} (default): one report for all samples.
#'   \code{FALSE}: one single-sample report per sample (or for
#'   \code{sample_name}); each file only contains that sample's calls, QC rows
#'   and references -- no other sample is named in it apart from its own
#'   reference samples.
#' @param sample_name Optional sample name(s) (only used when \code{global = FALSE}).
#' @param sample_table Optional path to the sample table (TSV with
#'   \code{sample_name} and \code{gender}).
#' @param ref_bams Optional path to the external reference BAM list (TSV with a
#'   \code{bam} column).
#' @param log_file Optional path to the pipeline log (warnings and a tail are shown).
#' @param pdf_output Deprecated and ignored (the old knit-time PDF export is
#'   gone). Use the browser's Print dialog: every tab is expanded for printing.
#' @param prefix Optional output prefix: the global report is named
#'   \code{ECHO_<prefix>_report.html}.
#' @param coverage_rdata Optional path to the coverage RData (output of
#'   \code{run_bam_coverage()}); used to list samples that were excluded
#'   before CNV calling.
#'
#' @details Optional \code{settings}: \code{report_plotly} (\code{"inline"},
#'   default, or \code{"cdn"} for a ~3.6 MB smaller file that needs internet to
#'   draw plots) and \code{report_max_profiles} (default \code{500}): maximum
#'   number of calls for which the per-call coverage profile is embedded
#'   (HIGH and MEDIUM confidence calls first).
#'
#' @return Invisibly, the path(s) of the HTML file(s) written.
#' @export
generate_report <- function(summary_rdata, qc_metrics_file, output_dir,
                            settings, config = NULL, global = TRUE, sample_name = NULL,
                            sample_table = NULL, ref_bams = NULL, log_file = NULL,
                            pdf_output = FALSE, prefix = NULL, coverage_rdata = NULL) {
  stop_if_not_file(summary_rdata, paste0("[ERROR] Summary RData not found: ", summary_rdata))
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("[ERROR] Package 'jsonlite' is required for report generation.")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) stop("[ERROR] Failed to create output directory: ", output_dir)
  }
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  settings   <- settings %||% list()
  if (isTRUE(pdf_output))
    message("[WARNING] pdf_output is no longer supported by the interactive report; ",
            "use the browser's Print dialog (all tabs are expanded when printing).")

  rd <- load_report_inputs(summary_rdata, qc_metrics_file, coverage_rdata, sample_table, ref_bams)

  if (global) {
    out_file <- file.path(output_dir, if (!is.null(prefix) && nzchar(prefix))
      paste0("ECHO_", prefix, "_report.html") else "ECHO_global_report.html")
    payload <- build_report_payload(rd, settings, config, log_file, prefix, focus = NULL)
    write_report_html(payload, out_file, settings)
    message("[INFO] Global report written: ", out_file)
    return(invisible(out_file))
  }

  samples <- sample_name %||% unique(as.character(rd$cnv_calls$Sample))
  if (!length(samples)) {
    message("[INFO] No CNV calls -- skipping per-sample reports.")
    return(invisible(NULL))
  }
  out_files <- character()
  for (s in samples) {
    sample_dir <- file.path(output_dir, sanitize_filename(s))
    dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(sample_dir, paste0("ECHO_report_", sanitize_filename(s), ".html"))
    payload  <- build_report_payload(rd, settings, config, log_file, prefix, focus = s)
    write_report_html(payload, out_file, settings)
    out_files <- c(out_files, out_file)
    message("[INFO] Report written: ", out_file)
  }
  invisible(out_files)
}

# ----------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------

#' Load and normalise everything the report needs
#' @noRd
load_report_inputs <- function(summary_rdata, qc_metrics_file, coverage_rdata,
                               sample_table, ref_bams) {
  env <- new.env()
  load(summary_rdata, envir = env)
  need <- c("counts", "bed_file", "sample_names")
  miss <- setdiff(need, ls(env))
  if (length(miss)) stop("[ERROR] Summary RData lacks: ", paste(miss, collapse = ", "))

  cnv_calls <- env$cnv_calls
  if (is.null(cnv_calls)) cnv_calls <- data.frame()
  bed_file <- as.data.frame(env$bed_file, stringsAsFactors = FALSE)
  if (!"chromosome" %in% colnames(bed_file)) {
    chr_col <- intersect(colnames(bed_file), c("chr", "chrom", "seqnames"))[1]
    if (is.na(chr_col)) stop("[ERROR] bed_file has no chromosome column.")
    bed_file$chromosome <- bed_file[[chr_col]]
  }
  if (!"gene" %in% colnames(bed_file)) bed_file$gene <- NA_character_

  analysed <- as.character(env$sample_names)
  all_samples <- analysed
  if (!is.null(coverage_rdata) && file.exists(coverage_rdata)) {
    cov_names <- tryCatch(load_rdata(coverage_rdata, required = "sample_names")$sample_names,
                          error = function(e) NULL)
    if (length(cov_names)) all_samples <- unique(c(as.character(cov_names), analysed))
  }
  all_samples <- all_samples[all_samples %in% colnames(env$counts)]
  analysed    <- analysed[analysed %in% all_samples]

  qc <- NULL
  if (!is.null(qc_metrics_file) && file.exists(qc_metrics_file) && file.size(qc_metrics_file) > 0) {
    qc <- tryCatch(data.table::fread(qc_metrics_file, header = TRUE, fill = TRUE, data.table = FALSE,
                                     colClasses = "character"),
                   error = function(e) NULL)
  }
  if (is.null(qc) || !all(c("Sample", "Type", "Details") %in% colnames(qc))) {
    if (!is.null(qc_metrics_file) && file.exists(qc_metrics_file))
      message("[WARNING] QC metrics file is empty or unreadable; the QC section only shows what can be recomputed.")
    qc <- NULL
  } else {
    for (cn in c("Exon", "Gene")) if (!cn %in% colnames(qc)) qc[[cn]] <- ""
    if (!"Severity" %in% colnames(qc)) {
      # QC files written before the Severity column existed
      qc$Severity <- ifelse(qc$Type %in% c("Whole exon", "Exon variability"), "info", "warning")
    }
    qc <- qc[, c("Sample", "Exon", "Type", "Details", "Gene", "Severity"), drop = FALSE]
    qc[is.na(qc)] <- ""
  }

  sdf <- NULL
  if (!is.null(sample_table) && file.exists(sample_table)) {
    sdf <- tryCatch(utils::read.table(sample_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                                      check.names = FALSE, quote = "", comment.char = ""),
                    error = function(e) NULL)
    if (is.null(sdf) || !"sample_name" %in% colnames(sdf)) {
      message("[WARNING] sample_table is unreadable or lacks a 'sample_name' column; ignored by the report.")
      sdf <- NULL
    }
  }
  gender_map <- env$gender_map
  if (is.null(sdf) && !is.null(gender_map) && length(gender_map)) {
    sdf <- data.frame(sample_name = names(gender_map), gender = as.character(gender_map),
                      stringsAsFactors = FALSE)
  }

  rb <- NULL
  if (!is.null(ref_bams) && is.character(ref_bams) && file.exists(ref_bams)) {
    rb <- tryCatch(utils::read.table(ref_bams, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                                     check.names = FALSE, quote = "", comment.char = ""),
                   error = function(e) NULL)
    if (is.null(rb) || !"bam" %in% colnames(rb)) rb <- NULL
  }

  list(cnv_calls = cnv_calls, counts = env$counts, bed_file = bed_file,
       models = env$models, refs = env$refs, analysed = analysed, all_samples = all_samples,
       qc = qc, sample_df = sdf, ref_bams = rb)
}

# ----------------------------------------------------------------------------
# Payload
# ----------------------------------------------------------------------------

rep_arr   <- function(x) I(unname(x))

#' data.frame -> column-oriented table the JS app consumes
#' @noRd
rep_table <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  data <- lapply(df, function(v) {
    if (is.factor(v)) v <- as.character(v)
    if (is.double(v)) v <- round(v, 4)
    I(unname(v))
  })
  list(cols = I(colnames(df)), n = nrow(df), data = data)
}

rep_empty_table <- function(cols) {
  rep_table(stats::setNames(rep(list(character()), length(cols)), cols))
}

`%nz%` <- function(a, b) if (is.null(a) || !length(a) || (length(a) == 1 && is.na(a))) b else a

#' @noRd
chr_scope_label <- function(modechrom) {
  switch(modechrom, X = "chrX", Y = "chrY", "autosomes")
}

#' Build the full JSON payload of the report
#' @noRd
build_report_payload <- function(rd, settings, config, log_file, prefix, focus = NULL) {
  modechrom <- toupper(as.character(settings$modechrom %||% "A")[1])
  if (!modechrom %in% c("A", "X", "Y")) modechrom <- "A"
  counts   <- rd$counts
  bed_file <- rd$bed_file
  calls    <- rd$cnv_calls
  has_calls <- is.data.frame(calls) && nrow(calls) > 0
  all_samples <- rd$all_samples
  analysed    <- rd$analysed

  thresholds <- list(
    min_corr        = settings$min_corr %||% 0.98,
    min_cov         = settings$min_cov %||% 100,
    min_total_reads = settings$min_total_reads %||% 300000,
    max_exon_cv     = settings$max_exon_cv %||% 0.5
  )

  # ---- chromosome scope --------------------------------------------------
  scope      <- chromosome_scope(bed_file$chromosome, modechrom)
  scope_rows <- chromosome_scope_warnings(scope, modechrom)
  qc <- rd$qc
  if (is.null(qc)) {
    qc <- data.frame(Sample = character(), Exon = character(), Type = character(),
                     Details = character(), Gene = character(), Severity = character(),
                     stringsAsFactors = FALSE)
  }
  # Make sure the scope warnings are present even when the QC file predates them
  # (or was not supplied): the report is the place users look first.
  if (nrow(scope_rows)) {
    have <- paste(qc$Type, qc$Details)
    add  <- scope_rows[!paste(scope_rows$Type, scope_rows$Details) %in% have, , drop = FALSE]
    if (nrow(add)) {
      qc <- rbind(qc, data.frame(Sample = "All", Exon = "All", Type = add$Type, Details = add$Details,
                                 Gene = "All", Severity = add$Severity, stringsAsFactors = FALSE))
    }
  }
  qc$Severity[!qc$Severity %in% c("error", "warning", "info")] <- "warning"

  chr_vec <- as.character(bed_file$chromosome)
  chr_ord <- unique(chr_vec)
  chr_df <- do.call(rbind, lapply(chr_ord, function(ch) {
    ii <- which(chr_vec == ch)
    data.frame(chromosome = ch, targets = length(ii), first = min(ii), last = max(ii),
               class = chromosome_class(ch), stringsAsFactors = FALSE)
  }))
  chr_df$analysed <- chr_df$class == switch(modechrom, X = "X", Y = "Y", "autosome")

  # ---- per-sample metrics -------------------------------------------------
  mat <- as.matrix(counts[, all_samples, drop = FALSE])
  storage.mode(mat) <- "double"
  med_depth   <- if (requireNamespace("matrixStats", quietly = TRUE))
    matrixStats::colMedians(mat, na.rm = TRUE) else apply(mat, 2, stats::median, na.rm = TRUE)
  total_reads <- colSums(mat, na.rm = TRUE)
  cor_mat <- NULL
  max_corr <- rep(NA_real_, length(all_samples))
  if (length(all_samples) >= 2) {
    cor_mat <- tryCatch(stats::cor(mat, method = "spearman", use = "pairwise.complete.obs"),
                        error = function(e) NULL)
    if (!is.null(cor_mat)) {
      cm <- cor_mat; diag(cm) <- NA_real_
      max_corr <- apply(cm, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
    }
  }

  models <- rd$models %||% list()
  refs   <- rd$refs %||% list()
  phi <- vapply(all_samples, function(s) {
    m <- models[[s]]
    if (is.null(m) || !length(m) || is.na(m[1])) NA_real_ else as.numeric(m[1])
  }, numeric(1))
  n_refs <- vapply(all_samples, function(s) length(refs[[s]] %||% character()), integer(1))

  # ExomeDepth labels calls "deletion" / "duplication"; the report uses DEL / DUP.
  call_type <- if (has_calls) {
    tt <- tolower(as.character(calls$Type))
    ifelse(grepl("^del", tt), "DEL", ifelse(grepl("^dup", tt), "DUP", toupper(tt)))
  } else character()
  conf_count <- function(s, lv) if (has_calls) sum(calls$Sample == s & calls$Confidence == lv, na.rm = TRUE) else 0L
  type_count <- function(s, ty) if (has_calls) sum(calls$Sample == s & call_type == ty, na.rm = TRUE) else 0L

  qc_s <- qc[qc$Sample %in% all_samples, , drop = FALSE]
  status <- vapply(all_samples, function(s) {
    sv <- qc_s$Severity[qc_s$Sample == s]
    if (!s %in% analysed) "Excluded" else if ("error" %in% sv) "Fail" else if ("warning" %in% sv) "Warning" else "OK"
  }, character(1))
  flags_txt <- vapply(all_samples, function(s) {
    k <- which(qc_s$Sample == s)
    txt <- if (length(k)) paste0(qc_s$Type[k], ": ", qc_s$Details[k], collapse = " | ") else ""
    if (!s %in% analysed) txt <- paste0("Excluded from CNV calling (sample outlier QC or not eligible for modechrom ",
                                        modechrom, ")", if (nzchar(txt)) " | " else "", txt)
    txt
  }, character(1))

  samples_df <- data.frame(
    Sample = all_samples, Status = status, Analysed = ifelse(all_samples %in% analysed, "Yes", "No"),
    `Median depth` = med_depth, `Total reads` = total_reads, `Max correlation` = max_corr,
    `References` = n_refs, Phi = phi,
    `CNV calls` = vapply(all_samples, function(s) if (has_calls) sum(calls$Sample == s) else 0, numeric(1)),
    HIGH = vapply(all_samples, conf_count, numeric(1), lv = "HIGH"),
    MEDIUM = vapply(all_samples, conf_count, numeric(1), lv = "MEDIUM"),
    LOW = vapply(all_samples, conf_count, numeric(1), lv = "LOW"),
    DEL = vapply(all_samples, type_count, numeric(1), ty = "DEL"),
    DUP = vapply(all_samples, type_count, numeric(1), ty = "DUP"),
    check.names = FALSE, stringsAsFactors = FALSE)
  sdf <- rd$sample_df
  if (!is.null(sdf)) {
    m <- match(all_samples, as.character(sdf$sample_name))
    extra <- setdiff(colnames(sdf), "sample_name")
    extra <- utils::head(extra, 8)
    for (cn in extra) {
      nm <- if (cn %in% colnames(samples_df)) paste0(cn, " (table)") else sub("^gender$", "Gender", cn)
      v  <- sdf[[cn]][m]
      samples_df[[nm]] <- if (is.numeric(v)) v else ifelse(is.na(v), "", as.character(v))
    }
  }
  samples_df[["QC flags"]] <- flags_txt

  # ---- CNV table and per-call windows ---------------------------------------
  max_profiles <- as.integer(settings$report_max_profiles %||% 500L)
  windows <- list(); calib <- rep("", if (has_calls) nrow(calls) else 0L)
  cnv_df <- NULL
  windows_limited <- FALSE
  if (has_calls) {
    exon_index <- resolve_exon_index(bed_file)
    calls_h <- harmonise_chr_prefix(bed_file, calls)
    rank <- match(calls$Confidence, c("HIGH", "MEDIUM", "LOW"), nomatch = 4L)
    todo <- order(rank, calls$CNV.ID)
    if (length(todo) > max_profiles) { windows_limited <- TRUE; todo <- todo[seq_len(max_profiles)] }
    for (i in todo) {
      w <- tryCatch({
        pd <- prepare_plot_data(calls_h[i, , drop = FALSE], counts, bed_file, exon_index, models, refs,
                                gene_gap = settings$plot_gene_gap %||% 1)
        if (is.null(pd)) NULL else window_from_plot_data(pd, bed_file, exon_index)
      }, error = function(e) {
        message("[WARNING] Report: coverage profile unavailable for CNV.ID ", calls$CNV.ID[i], ": ", conditionMessage(e))
        NULL
      })
      if (!is.null(w)) {
        windows[[as.character(calls$CNV.ID[i])]] <- w$window
        calib[i] <- if (isTRUE(w$flag)) "check" else "ok"
      }
    }
    pcol <- function(nm, default = NA) if (nm %in% colnames(calls)) calls[[nm]] else rep(default, nrow(calls))
    sp <- pcol("Start.p"); ep <- pcol("End.p")
    cnv_df <- data.frame(
      ID = calls$CNV.ID, Sample = as.character(calls$Sample), Chr = as.character(calls$Chromosome),
      Start = calls$Start, End = calls$End, Gene = as.character(pcol("Gene", "Unknown")),
      Type = call_type,
      Exons = ifelse(is.na(sp) | is.na(ep), "", ifelse(sp == ep, as.character(sp), paste0(sp, "-", ep))),
      Nexons = pcol("Nexons"), Ratio = pcol("Reads.ratio"), BF = pcol("BF"),
      Expected = pcol("Reads.expected"), Observed = pcol("Reads.observed"),
      Confidence = as.character(pcol("Confidence", "LOW")), Correlation = pcol("Correlation"),
      N.comp = pcol("N.comp"), Comparators = as.character(pcol("Comparator.name", "")),
      Calibration = calib, Genomic = as.character(pcol("Genomic.ID", "")),
      gidx = as.numeric(pcol("global_start", NA)), stringsAsFactors = FALSE, check.names = FALSE)
    cnv_df$Gene[is.na(cnv_df$Gene) | !nzchar(cnv_df$Gene)] <- "Unknown"
  }

  # ---- cohort plots -----------------------------------------------------------
  pca <- report_pca(mat, all_samples, rd$sample_df)
  corr <- NULL
  if (!is.null(cor_mat) && length(all_samples) >= 2 && length(all_samples) <= 60) {
    ord <- tryCatch(stats::hclust(stats::as.dist(1 - cor_mat))$order, error = function(e) seq_along(all_samples))
    cm <- cor_mat[ord, ord, drop = FALSE]
    corr <- list(names = I(all_samples[ord]),
                 z = lapply(seq_len(nrow(cm)), function(i) I(round(unname(cm[i, ]), 3))))
  }
  exon_med <- if (requireNamespace("matrixStats", quietly = TRUE))
    matrixStats::rowMedians(mat, na.rm = TRUE) else apply(mat, 1, stats::median, na.rm = TRUE)
  lm <- log10(pmax(exon_med[is.finite(exon_med)], 1))
  depth_hist <- NULL
  if (length(lm) > 1) {
    br <- seq(0, max(1, ceiling(max(lm) * 10) / 10), length.out = 41)
    hh <- graphics::hist(lm, breaks = br, plot = FALSE, include.lowest = TRUE)
    depth_hist <- list(mid = I(round(hh$mids, 3)), counts = I(hh$counts))
  }

  # ---- alerts ---------------------------------------------------------------------
  alerts <- build_report_alerts(qc, scope_rows, analysed, all_samples, samples_df, cnv_df, calib,
                                modechrom, thresholds, log_file, qc_available = !is.null(rd$qc))

  # ---- log, settings, scoring ----------------------------------------------------------
  log_info <- read_report_log(log_file)
  settings_kv <- report_settings_kv(settings, modechrom)
  scoring <- report_scoring(settings)

  scope_incomplete <- any(scope$class %in% c("autosome", "X", "Y") & !scope$analysed)
  meta <- list(
    title = paste0("ECHO report", if (!is.null(prefix) && nzchar(prefix)) paste0(" - ", prefix) else "",
                   if (!is.null(focus)) paste0(" - ", focus) else ""),
    prefix = prefix %nz% "", generated = format(Sys.time(), "%Y-%m-%d %H:%M"),
    version = tryCatch(as.character(utils::packageVersion("ECHO")), error = function(e) ""),
    modechrom = modechrom, modechrom_label = chr_scope_label(modechrom),
    genome_version = (config$genome_version %||% settings$genome_version) %nz% "",
    n_targets = nrow(bed_file), n_samples = length(all_samples), n_analysed = length(analysed),
    n_calls = if (has_calls) nrow(calls) else 0L,
    focus = if (is.null(focus)) NULL else focus,
    scope_incomplete = scope_incomplete, qc_available = !is.null(rd$qc),
    off_target_pattern = settings$off_target_gene_pattern %nz% "",
    windows_limited = windows_limited, windows_limit = max_profiles
  )

  payload <- list(
    meta = meta, thresholds = thresholds,
    scope = rep_table(scope), chr = rep_table(chr_df),
    alerts = alerts, qc = rep_table(qc), samples = rep_table(samples_df),
    cnvs = if (is.null(cnv_df)) rep_empty_table(c("ID", "Sample", "Chr", "Start", "End", "Gene", "Type",
                                                  "Confidence")) else rep_table(cnv_df),
    windows = windows, refs = lapply(refs[all_samples[all_samples %in% names(refs)]], function(x) I(as.character(x))),
    pca = pca, corr = corr, depth_hist = depth_hist,
    settings = lapply(seq_along(settings_kv), function(i) I(c(names(settings_kv)[i], settings_kv[[i]]))),
    scoring = scoring,
    sample_table = if (is.null(rd$sample_df)) NULL else rep_table(utils::head(rd$sample_df, 2000)),
    ref_bams = if (is.null(rd$ref_bams)) NULL else rep_table(utils::head(rd$ref_bams, 2000)),
    log = log_info,
    config_yaml = report_config_yaml(config),
    session = I(utils::capture.output(utils::sessionInfo()))
  )

  if (!is.null(focus)) payload <- restrict_payload_to_sample(payload, focus)
  payload
}

#' One call's plotting window -> compact JSON-ready list
#' @noRd
window_from_plot_data <- function(pd, bed_file, exon_index) {
  er <- pd$exon_range
  cd <- pd$cov_data
  ref_names <- names(pd$ref_z_list)
  test <- cd$coverage[cd$group == "Test sample"]
  ref_cov <- lapply(ref_names, function(r) I(round(cd$coverage[cd$group == r], 3)))
  genes <- as.character(bed_file$gene[er]); genes[is.na(genes)] <- ""
  list(
    flag = isTRUE(pd$bg_calibration$flag),
    window = list(
      px = I(round(pd$px, 3)), xlab = I(exon_index[er]), gene = I(genes),
      gene_group = I(as.integer(pd$gene_group)),
      affected = I(er %in% (pd$idx_start:pd$idx_end)),
      test = I(round(test, 3)), ref_names = I(ref_names), ref_cov = ref_cov,
      ratio = I(round(pd$ci_data$ratio, 3)), lo = I(round(pd$ci_data$lo, 3)), hi = I(round(pd$ci_data$hi, 3)),
      observed = I(round(pd$observed, 1)), expected = I(round(pd$expected, 1)),
      z = I(round(pd$z_data$z, 3)),
      ref_z = unname(lapply(pd$ref_z_list, function(v) I(round(unname(v), 3)))),
      subtitle = pd$ci_subtitle
    )
  )
}

#' PCA scores of the coverage profiles
#' @noRd
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
      g <- as.character(sample_df$gender[match(sample_names, as.character(sample_df$sample_name))])
      grp <- ifelse(is.na(g) | !nzchar(g), "unknown", g)
    }
    list(names = I(sample_names), pc1 = I(round(unname(p$x[, 1]), 3)), pc2 = I(round(unname(p$x[, 2]), 3)),
         var1 = round(ve[1], 1), var2 = round(ve[2], 1), group = I(grp),
         n_rows = sum(keep), n_samples = length(sample_names))
  }, error = function(e) {
    message("[WARNING] Report: PCA skipped: ", conditionMessage(e))
    NULL
  })
}

#' Alerts shown at the top of the tabs
#' @noRd
build_report_alerts <- function(qc, scope_rows, analysed, all_samples, samples_df, cnv_df, calib,
                                modechrom, thresholds, log_file, qc_available = TRUE) {
  out <- list()
  add <- function(kind, severity, title, text, tab = "qc", samples = NULL) {
    out[[length(out) + 1L]] <<- list(kind = kind, severity = severity, title = title, text = text,
                                     tab = tab, samples = if (length(samples)) I(as.character(samples)) else NULL)
  }
  # chromosome scope / missing classes
  for (i in seq_len(nrow(scope_rows))) {
    d <- scope_rows$Details[i]
    title <- if (scope_rows$Type[i] == "Missing Chromosome") {
      "Chromosome class missing from the BED file"
    } else if (grepl("^Autosomes", d)) {
      "Autosomes are NOT analysed in this run"
    } else if (grepl("^Gonosomes", d)) {
      "Gonosomes (chrX / chrY) are NOT analysed in this run"
    } else "Non-canonical contigs"
    add("scope", scope_rows$Severity[i], title, d, "qc")
  }
  # cohort size
  k <- which(qc$Type == "Cohort size")
  for (i in k) add("cohort", qc$Severity[i], "Small cohort", qc$Details[i], "qc")
  if (!qc_available) {
    add("qc", "info", "No QC metrics file", paste0("The QC metrics TSV of this run was not found: sample and exon ",
        "QC flags are missing from this report (chromosome-scope warnings are still computed)."), "qc")
    if (length(analysed) < 3 && !length(k)) {
      add("cohort", "warning", "Small cohort",
          sprintf(paste0("Only %d sample(s): ExomeDepth needs several comparable samples to build a reference ",
                         "set; expect few or unreliable calls."), length(analysed)), "qc")
    }
  }
  # excluded samples
  excl <- setdiff(all_samples, analysed)
  if (length(excl)) {
    add("sample", "warning", paste0(length(excl), " sample(s) excluded from CNV calling"),
        paste0(paste(utils::head(excl, 10), collapse = ", "), if (length(excl) > 10) ", ..." else "",
               ": flagged as outliers by the sample QC, or not eligible for modechrom ", modechrom, "."),
        "samples", excl)
  }
  # whole-sample QC, grouped by kind of failure
  ws <- qc[qc$Type == "Whole sample" & qc$Sample %in% all_samples, , drop = FALSE]
  if (nrow(ws)) {
    kind <- sub(":.*$", "", ws$Details)
    for (kd in unique(kind)) {
      rows <- ws[kind == kd, , drop = FALSE]
      sev <- if ("error" %in% rows$Severity) "error" else "warning"
      detail <- paste0(rows$Sample, " (", sub("^[^:]*:\\s*", "", sub("\\s*\\(.*$", "", rows$Details)), ")")
      add("sample", sev, paste0(kd, " in ", length(unique(rows$Sample)), " sample(s)"),
          paste0(paste(utils::head(detail, 12), collapse = ", "), if (length(detail) > 12) ", ..." else ""),
          "qc", unique(rows$Sample))
    }
  }
  # exon-level
  ne <- sum(qc$Type == "Whole exon")
  if (ne) add("qc", "info", paste0(ne, " target(s) with low median depth"),
              paste0("Median depth below min_cov = ", thresholds$min_cov, ". Calls on them are unreliable."), "qc")
  nv <- sum(qc$Type == "Exon variability")
  if (nv) add("qc", "info", paste0(nv, " target(s) with high variability"),
              paste0("Coefficient of variation above max_exon_cv = ", thresholds$max_exon_cv, "."), "qc")
  # calls
  if (is.null(cnv_df) || nrow(cnv_df) == 0) {
    add("calls", "info", "No CNV called", "ExomeDepth reported no CNV for the analysed chromosome class.", "cnvs")
  } else {
    nchk <- sum(calib == "check")
    if (nchk) add("calls", "warning", paste0(nchk, " call(s) with a poorly calibrated background"),
                  "More background exons than expected lie outside the 95% interval: check the reference set or the region (see the CNV calls tab, 'Calibration').",
                  "cnvs", unique(cnv_df$Sample[calib == "check"]))
  }
  # log
  lg <- read_report_log(log_file)
  nw <- length(lg$warnings)
  if (nw) add("run", "info", paste0(nw, " warning/error line(s) in the pipeline log"), "See the Run tab.", "run")
  out
}

#' @noRd
read_report_log <- function(log_file, tail_n = 400L) {
  if (is.null(log_file) || !is.character(log_file) || !file.exists(log_file)) {
    return(list(warnings = I(character()), tail = I(character()), truncated = FALSE))
  }
  ln <- tryCatch(readLines(log_file, warn = FALSE, encoding = "UTF-8"), error = function(e) character())
  ln <- iconv(ln, "UTF-8", "UTF-8", sub = "?")
  w  <- unique(ln[grepl("^\\[(WARNING|ERROR)\\]", ln)])
  list(warnings = I(utils::head(w, 300)), tail = I(utils::tail(ln, tail_n)), truncated = length(ln) > tail_n)
}

#' @noRd
report_settings_kv <- function(settings, modechrom) {
  keep <- settings[vapply(settings, function(v) is.atomic(v) && !is.null(v), logical(1))]
  vals <- lapply(keep, function(v) paste(format(v, scientific = FALSE, trim = TRUE), collapse = ", "))
  if (is.null(keep$modechrom)) vals <- c(list(modechrom = modechrom), vals)
  vals
}

#' @noRd
report_config_yaml <- function(config) {
  if (is.null(config)) return(NULL)
  tryCatch({
    y <- yaml::as.yaml(config)
    ln <- strsplit(y, "\n", fixed = TRUE)[[1]]
    paste(utils::head(ln, 400), collapse = "\n")
  }, error = function(e) NULL)
}

#' Confidence rules as actually applied by score_cnv_confidence()
#' @noRd
report_scoring <- function(settings) {
  f <- formals(score_cnv_confidence)
  g <- function(nm) settings[[nm]] %||% f[[nm]]
  rng <- function(lo, hi) sprintf("<= %.2f or >= %.2f", lo, hi)
  forced <- settings$score_low_confidence_genes
  if (is.null(forced)) forced <- eval(f$score_low_confidence_genes)
  list(
    rows = list(
      I(c("HIGH", sprintf(">= %.3f", g("score_high_corr")), sprintf(">= %d", as.integer(g("score_high_refs"))),
          rng(g("score_high_ratio_low"), g("score_high_ratio_high")))),
      I(c("MEDIUM", sprintf(">= %.3f", g("score_med_corr")), sprintf(">= %d", as.integer(g("score_med_refs"))),
          rng(g("score_low_ratio_low"), g("score_low_ratio_high")))),
      I(c("LOW", "any", "any", "otherwise"))),
    forced_low = I(as.character(forced)),
    terminal = if (isTRUE(settings$remove_terminal_only)) {
      "Terminal-only calls (a single end exon of a multi-exon gene) are removed (remove_terminal_only)."
    } else if (isTRUE(settings$penalize_terminal_only)) {
      "Terminal-only calls are downgraded to LOW (penalize_terminal_only)."
    } else NULL
  )
}

#' Keep a single sample in a payload (per-sample reports)
#'
#' Nothing about the other samples is kept apart from the reference samples
#' this sample was compared with (they appear in the coverage profiles).
#' @noRd
restrict_payload_to_sample <- function(p, s) {
  rows <- function(tb, keep) {
    tb$data <- lapply(tb$data, function(v) I(unname(v[keep])))
    tb$n <- sum(keep)
    tb
  }
  p$samples <- rows(p$samples, p$samples$data$Sample == s)
  p$cnvs    <- rows(p$cnvs, p$cnvs$data$Sample == s)
  p$qc      <- rows(p$qc, p$qc$data$Sample %in% c(s, "All"))
  ids <- as.character(p$cnvs$data$ID)
  p$windows <- p$windows[names(p$windows) %in% ids]
  p$refs <- p$refs[names(p$refs) == s]
  # Alerts that name several samples are reduced to a generic sentence, and
  # everything that could name other samples (log, config, sample table) is dropped.
  p$alerts <- lapply(Filter(function(a) a$kind != "run" && (is.null(a$samples) || s %in% a$samples), p$alerts),
                     function(a) {
                       if (length(a$samples) > 1) {
                         a$samples <- I(s)
                         a$title <- sub(" in [0-9]+ sample\\(s\\)$", "", a$title)
                         a$title <- sub("^[0-9]+ call\\(s\\)", "Call(s)", a$title)
                         a$text <- if (a$kind == "calls") a$text else "This sample is among the affected ones; see the QC & warnings tab for its own values."
                       }
                       a
                     })
  p$pca <- NULL; p$corr <- NULL
  p$sample_table <- NULL
  p$ref_bams <- NULL
  p$config_yaml <- NULL
  p$log <- list(warnings = I(character()), tail = I(character()), truncated = FALSE)
  p$session <- I(character())
  p$meta$n_samples <- 1L
  p
}

# ----------------------------------------------------------------------------
# HTML assembly
# ----------------------------------------------------------------------------

#' @noRd
report_asset_path <- function(name) {
  p <- system.file("report", name, package = "ECHO")
  if (!nzchar(p)) {
    for (cand in c(file.path("inst", "report", name), file.path("..", "..", "inst", "report", name))) {
      if (file.exists(cand)) { p <- cand; break }
    }
  }
  if (!nzchar(p) || !file.exists(p)) stop("[ERROR] Report asset '", name, "' not found; reinstall ECHO.")
  p
}

#' @noRd
read_utf8 <- function(path) {
  raw <- readBin(path, "raw", file.size(path))
  txt <- rawToChar(raw)
  Encoding(txt) <- "UTF-8"
  txt
}

#' Inert JSON for a <script type="application/json"> element
#' @noRd
report_json <- function(x) {
  json <- as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "null",
                                        digits = 6, force = TRUE))
  json <- gsub("</", "<\\/", json, fixed = TRUE)
  gsub("<!--", "\\u003c!--", json, fixed = TRUE)
}

#' Location of the plotly.js bundle shipped with the plotly R package
#' @noRd
report_plotly_js <- function() {
  deps <- tryCatch(suppressWarnings(suppressMessages(htmltools::findDependencies(plotly::plot_ly()))),
                   error = function(e) list())
  dep <- Filter(function(d) identical(d$name, "plotly-main"), deps)
  if (!length(dep)) return(NULL)
  dep <- dep[[1]]
  path <- file.path(dep$src$file, dep$script)
  if (!file.exists(path)) path <- system.file(dep$src$file, dep$script, package = dep$package %||% "plotly")
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  gsub("</script", "<\\/script", read_utf8(path), fixed = TRUE)
}

#' @noRd
write_report_html <- function(payload, out_file, settings = list()) {
  css <- read_utf8(report_asset_path("echo_report.css"))
  js  <- read_utf8(report_asset_path("echo_report.js"))
  mode <- tolower(as.character(settings$report_plotly %||% "inline")[1])
  plotly_tag <- NULL
  if (mode != "cdn") {
    src <- report_plotly_js()
    if (is.null(src)) {
      message("[WARNING] plotly.js bundle not found; the report will load it from a CDN (needs internet).")
    } else {
      plotly_tag <- paste0("<script>\n", src, "\n</script>\n")
    }
  }
  if (is.null(plotly_tag)) {
    plotly_tag <- "<script src=\"https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js\"></script>\n"
  }
  title <- htmltools::htmlEscape(payload$meta$title)
  html <- paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>", title, "</title>\n<style>\n", css, "\n</style>\n</head>\n<body>\n",
    "<header class=\"top\"><div class=\"top-inner\">\n",
    "<div class=\"title-row\"><h1 id=\"title\">", title, "</h1><span id=\"badges\" class=\"title-row\"></span>",
    "<span class=\"spacer\"></span><span id=\"sample-pick\" class=\"sample-pick\"></span>",
    "<span class=\"muted small\" id=\"generated\"></span>",
    "<button class=\"btn\" id=\"theme-btn\" type=\"button\" aria-label=\"Colour theme\">Theme</button></div>\n",
    "<nav class=\"tabs\" id=\"tabs\" role=\"tablist\" aria-label=\"Report sections\"></nav>\n",
    "</div></header>\n<main id=\"panels\"></main>\n",
    "<noscript><p style=\"padding:16px\">This report needs JavaScript enabled.</p></noscript>\n",
    plotly_tag,
    "<script id=\"echo-data\" type=\"application/json\">", report_json(payload), "</script>\n",
    "<script>\n", js, "\n</script>\n</body>\n</html>\n")
  con <- file(out_file, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(enc2utf8(html)), con)
  invisible(out_file)
}
