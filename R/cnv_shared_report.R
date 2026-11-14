# ============================================================================
# Report building blocks shared by ECHO and CANOPE.
#
# THIS FILE IS BYTE-IDENTICAL IN BOTH PACKAGES (a test enforces it). Edit it in
# one place and copy it to the other. It builds the pieces of the unified JSON
# payload that must be computed the same way whatever the caller: the pipeline
# log summary, PCA / correlation / depth histogram, the per-sample table with
# its flags and exclusion reasons, the alert list and the per-sample report
# restriction. The pipeline-specific parts (call table, per-call plot windows,
# scoring rules) stay in echo_report.R / canope_report.R.
# ============================================================================

rep_arr <- function(x) I(unname(x))

# ---------------------------------------------------------------------------
# Pipeline log
# ---------------------------------------------------------------------------

#' Summarise a pipeline log for the report
#'
#' Reads the \code{[LEVEL] YYYY-mm-dd HH:MM:SS text} lines written by
#' \code{\link{cnv_logger}} (older logs without timestamps are read too).
#'
#' @param log_file Path of the log (\code{NULL} / missing file allowed).
#' @param tail_n Number of trailing lines embedded for the "Run & log" tab.
#' @return List: \code{warnings} (unique WARNING/ERROR lines), \code{tail},
#'   \code{session}, \code{truncated}, \code{counts} (lines per level),
#'   \code{steps} (the "Step k/N" lines), \code{seconds} (wall time between the
#'   first and the last timestamp, or \code{NULL}).
#' @noRd
report_log <- function(log_file, tail_n = 1000L) {
  si <- I(utils::capture.output(print(utils::sessionInfo())))
  empty <- list(warnings = I(character()), tail = I(character()), session = si, truncated = FALSE,
                counts = list(DEBUG = 0L, INFO = 0L, WARNING = 0L, ERROR = 0L),
                steps = I(character()), seconds = NULL)
  if (is.null(log_file) || !is.character(log_file) || !length(log_file) || !file.exists(log_file))
    return(empty)
  ln <- tryCatch(readLines(log_file, warn = FALSE, encoding = "UTF-8"), error = function(e) character())
  ln <- iconv(ln, "UTF-8", "UTF-8", sub = "?")
  ln <- ln[nzchar(trimws(ln))]
  lvl <- sub("^\\[(DEBUG|INFO|WARNING|ERROR)\\].*$", "\\1", ln)
  lvl[!lvl %in% c("DEBUG", "INFO", "WARNING", "ERROR")] <- NA_character_
  counts <- lapply(c(DEBUG = "DEBUG", INFO = "INFO", WARNING = "WARNING", ERROR = "ERROR"),
                   function(k) sum(lvl == k, na.rm = TRUE))
  ts <- regmatches(ln, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}", ln))
  secs <- if (length(ts) >= 2L) {
    t <- as.POSIXct(ts, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    if (all(is.na(t))) NULL else as.numeric(difftime(max(t, na.rm = TRUE), min(t, na.rm = TRUE), units = "secs"))
  } else NULL
  list(warnings = I(utils::head(unique(ln[lvl %in% c("WARNING", "ERROR")]), 300L)),
       tail = I(utils::tail(ln, tail_n)),
       session = si, truncated = length(ln) > tail_n, counts = counts,
       steps = I(utils::head(ln[grepl("Step [0-9]+/[0-9?]+ completed", ln)], 60L)),
       seconds = secs)
}

# ---------------------------------------------------------------------------
# Cohort views
# ---------------------------------------------------------------------------

#' PCA block of the payload
#' @param mat Numeric matrix, targets x samples (raw counts).
#' @param group Optional group label per sample (e.g. sex) used for colouring.
#' @param max_pcs Number of principal components embedded.
#' @noRd
report_pca <- function(mat, group = NULL, max_pcs = 5L) {
  nm <- colnames(mat)
  if (length(nm) < 3L) return(NULL)
  tryCatch({
    lm <- log2(pmax(mat, 0) + 1)
    rv <- apply(lm, 1, stats::var, na.rm = TRUE)
    keep <- is.finite(rv) & rv > 0
    if (sum(keep) < 2L) return(NULL)
    p  <- stats::prcomp(t(lm[keep, , drop = FALSE]), scale. = TRUE, center = TRUE)
    k  <- min(max_pcs, ncol(p$x))
    ve <- summary(p)$importance[2, seq_len(k)] * 100
    grp <- if (is.null(group)) rep("All samples", length(nm)) else ifelse(is.na(group) | !nzchar(as.character(group)), "unknown", as.character(group))
    list(names = rep_arr(nm), samples = rep_arr(nm),
         pc1 = rep_arr(round(p$x[, 1], 3)), pc2 = rep_arr(round(p$x[, 2], 3)),
         var1 = round(ve[1], 1), var2 = round(ve[2], 1),
         pcs = stats::setNames(lapply(seq_len(k), function(i) rep_arr(round(p$x[, i], 3))), paste0("PC", seq_len(k))),
         var = rep_arr(round(ve, 1)), group = rep_arr(grp),
         n_rows = sum(keep), n_targets = sum(keep), n_samples = length(nm))
  }, error = function(e) { message("[WARNING] PCA skipped: ", conditionMessage(e)); NULL })
}

#' Clustered sample-to-sample correlation block
#' @param cor_mat Spearman correlation matrix (named).
#' @param max_n Largest cohort drawn (the matrix grows quadratically).
#' @noRd
report_corr <- function(cor_mat, max_n = 150L) {
  if (is.null(cor_mat) || nrow(cor_mat) < 2L || nrow(cor_mat) > max_n) return(NULL)
  ord <- if (nrow(cor_mat) >= 3L && !anyNA(cor_mat))
    tryCatch(stats::hclust(stats::as.dist(1 - cor_mat), method = "average")$order, error = function(e) seq_len(nrow(cor_mat)))
  else seq_len(nrow(cor_mat))
  cm <- cor_mat[ord, ord, drop = FALSE]
  list(names = rep_arr(rownames(cm)), z = lapply(seq_len(nrow(cm)), function(i) rep_arr(round(cm[i, ], 3))))
}

#' Histogram of the median depth of every target
#' @noRd
report_depth_hist <- function(mat) {
  if (!ncol(mat)) return(NULL)
  em <- row_medians(mat)
  lm <- log10(pmax(em[is.finite(em)], 1))
  if (length(lm) < 2L) return(NULL)
  br <- seq(0, max(1, ceiling(max(lm) * 10) / 10), length.out = 41L)
  hh <- graphics::hist(lm, breaks = br, plot = FALSE, include.lowest = TRUE)
  list(mid = rep_arr(round(hh$mids, 3)), counts = rep_arr(hh$counts))
}

#' Spearman correlation of the samples and each sample's best partner
#' @noRd
report_correlation <- function(mat) {
  cm <- NULL; best <- rep(NA_real_, ncol(mat))
  if (ncol(mat) >= 2L) {
    cm <- tryCatch(suppressWarnings(stats::cor(mat, method = "spearman", use = "pairwise.complete.obs")),
                   error = function(e) NULL)
    if (!is.null(cm)) {
      c2 <- cm; diag(c2) <- NA_real_   # by position: a duplicated sample (r = 1) stays visible
      best <- apply(c2, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
    }
  }
  list(matrix = cm, best = unname(best))
}

# ---------------------------------------------------------------------------
# Sample table
# ---------------------------------------------------------------------------

#' Why was a sample excluded?
#'
#' Priority: the reason recorded by \code{cnv_apply_qc()} (exact), then the
#' reason recovered from the pipeline log (older workspaces), then the first
#' per-sample threshold the sample fails, then a generic text.
#' @param excluded Data frame (\code{sample}, \code{reason}) or \code{NULL}.
#' @noRd
report_exclusion_reasons <- function(samples, excluded, log_info = NULL, metrics = NULL,
                                     thresholds = list(), status_msg = NULL) {
  logged <- if (!is.null(log_info)) recover_exclusion_reasons(samples, log_info) else list()
  vapply(samples, function(s) {
    if (!is.null(excluded) && nrow(excluded) && s %in% excluded$sample)
      return(excluded$reason[match(s, excluded$sample)])
    if (!is.null(status_msg) && !is.na(status_msg[s]) && nzchar(status_msg[s])) return(status_msg[[s]])
    if (!is.null(logged[[s]]) && !is.na(logged[[s]])) return(logged[[s]])
    if (!is.null(metrics)) {
      i <- match(s, metrics$sample)
      chk <- list(list("median_depth", thresholds$min_cov, "below", "median depth"),
                  list("total_reads", thresholds$min_total_reads, "below", "total reads"),
                  list("max_corr", thresholds$min_corr, "below", "max correlation"))
      for (k in chk) {
        v <- metrics[[k[[1]]]][i]
        if (is.finite(v) && is.finite(k[[2]] %||% NA_real_) && v < k[[2]])
          return(sprintf("Low %s: %s (< %s)", k[[4]], fmt_num(v, 2), fmt_num(k[[2]], 2)))
      }
    }
    "Excluded from calling by sample QC"
  }, character(1))
}

# Recover the "why" for each excluded sample from the pipeline log lines
# (workspaces written before qc_info existed).
recover_exclusion_reasons <- function(sample_names, log_info) {
  out <- stats::setNames(rep(NA_character_, length(sample_names)), sample_names)
  log_lines <- c(as.character(log_info$warnings %||% character()), as.character(log_info$tail %||% character()))
  pat1 <- "excluding[[:space:]]+[0-9]+[[:space:]]+outlier sample\\(s\\)[[:space:]]+from calling[[:space:]]*\\(([^)]*)\\)[[:space:]]*:[[:space:]]*(.+)$"
  pat2 <- "([0-9]+)[[:space:]]+sample\\(s\\)[[:space:]]+excluded[[:space:]]*:[[:space:]]*(.+)$"
  for (ln in log_lines) {
    m <- regmatches(ln, regexec(pat1, ln, perl = TRUE))[[1]]
    if (length(m) >= 3) {
      for (nm in strsplit(m[3], ",[[:space:]]*(?![^(]*\\))", perl = TRUE)[[1]]) {
        nm <- trimws(sub("[[:space:]]*\\(.*$", "", nm))
        if (nzchar(nm) && nm %in% names(out) && is.na(out[nm]))
          out[nm] <- sprintf("Outlier sample excluded from calling (%s)", trimws(m[2]))
      }
      next
    }
    m <- regmatches(ln, regexec(pat2, ln, perl = TRUE))[[1]]
    if (length(m) >= 3)
      for (nm in strsplit(m[3], ",[[:space:]]*")[[1]]) {
        nm <- trimws(nm)
        if (nzchar(nm) && nm %in% names(out) && is.na(out[nm])) out[nm] <- "Sample excluded from calling (sample QC failure)"
      }
  }
  as.list(out)
}

#' The per-sample table of the report (one schema for both pipelines)
#'
#' @param mat Raw counts (targets x samples) of every candidate sample.
#' @param analysed Samples that were actually called.
#' @param noise \code{sample_qc} table from \code{cnv_apply_qc()} (\code{sample},
#'   \code{noise_score}, \code{robust_z}); computed here if \code{NULL}.
#' @param calls Standardised calls: data frame with \code{sample}, \code{type}
#'   (DEL/DUP), \code{confidence}, optionally \code{recurrent} (logical).
#' @param n_refs,ref_cor,n_targets,n_dropped Named vectors (by sample) from
#'   the caller; \code{NA} when the pipeline does not provide them.
#' @param status Named character vector overriding the status of a sample
#'   (\code{"failed"}); the default is \code{"called"} / \code{"excluded"}.
#' @param metrics \code{sample_metrics} from \code{cnv_apply_qc()} (preferred
#'   source of depth / reads / correlation / noise).
#' @param qc The QC flag table (\code{Sample}, \code{Type}, \code{Details},
#'   \code{Severity}), used for the per-sample level.
#' @param extra_flags Optional named list (by sample) of character vectors
#'   appended to the \code{flags} text.
#' @return Data frame, lowercase column names only.
#' @noRd
report_sample_table <- function(mat, analysed, noise = NULL, calls = NULL, n_refs = NULL,
                                ref_cor = NULL, n_targets = NULL, n_dropped = NULL,
                                status = NULL, qc = NULL, qc_zscore = 3, outlier_sides = "upper",
                                sample_df = NULL, reasons = NULL, sex_inferred = NULL,
                                extra_level = NULL, metrics = NULL) {
  all_samples <- colnames(mat)
  nS <- length(all_samples)
  get <- function(v, default = NA_real_) {
    if (is.null(v)) return(rep(default, nS))
    out <- unname(v[all_samples]); out[is.na(names(v[all_samples]))] <- default; out
  }
  # The numbers that decided the QC (computed once on the raw counts of the
  # mode-filtered targets, see cnv_apply_qc()); recomputed from `mat` only for
  # workspaces written before they were stored.
  if (is.null(metrics) && nS >= 1L && nrow(mat) > 0L) metrics <- cnv_sample_metrics(mat, qc_zscore, outlier_sides)
  mi <- match(all_samples, metrics$sample)
  med_depth   <- metrics$median_depth[mi]
  total_reads <- metrics$total_reads[mi]
  max_corr    <- metrics$max_corr[mi]
  if (is.null(noise) && !is.null(metrics)) noise <- data.frame(sample = metrics$sample, noise_score = metrics$noise_score, robust_z = metrics$robust_z)
  noise_z <- if (is.null(noise)) rep(NA_real_, nS) else noise$robust_z[match(all_samples, noise$sample)]
  noise_s <- if (is.null(noise)) rep(NA_real_, nS) else noise$noise_score[match(all_samples, noise$sample)]

  has_calls <- is.data.frame(calls) && nrow(calls) > 0L
  cnt <- function(f) vapply(all_samples, function(s) if (has_calls) sum(calls$sample == s & f(calls), na.rm = TRUE) else 0, numeric(1))
  conf_n <- function(lv) cnt(function(d) toupper(as.character(d$confidence)) == lv)
  type_n <- function(ty) cnt(function(d) d$type == ty)
  rec_n  <- if (has_calls && "recurrent" %in% names(calls)) cnt(function(d) d$recurrent %in% TRUE) else rep(0, nS)

  st <- ifelse(all_samples %in% analysed, "called", "excluded")
  if (!is.null(status)) { o <- status[all_samples]; st[!is.na(o)] <- o[!is.na(o)] }
  qc_s <- if (is.null(qc)) data.frame(Sample = character(), Severity = character()) else qc[qc$Sample %in% all_samples, , drop = FALSE]
  level <- vapply(seq_len(nS), function(i) {
    s <- all_samples[i]; sv <- qc_s$Severity[qc_s$Sample == s]
    if (st[i] %in% c("excluded", "failed") || "error" %in% sv || "fail" %in% sv) "fail"
    else if ("warning" %in% sv || isTRUE(extra_level[[s]] == "warn")) "warn" else "ok"
  }, character(1))

  df <- data.frame(
    sample = all_samples, status = st, level = level,
    Analysed = ifelse(st == "called", "Yes", "No"),
    n_cnvs = cnt(function(d) rep(TRUE, nrow(d))),
    n_del = type_n("DEL"), n_dup = type_n("DUP"),
    n_high = conf_n("HIGH"), n_medium = conf_n("MEDIUM"), n_low = conf_n("LOW"), n_recurrent = rec_n,
    median_depth = round(med_depth, 1), total_reads = round(total_reads), max_corr = round(max_corr, 4),
    noise_score = round(noise_s, 4), noise_z = round(noise_z, 2),
    n_refs = get(n_refs), ref_cor = round(get(ref_cor), 3), n_targets = get(n_targets), n_dropped = get(n_dropped),
    sex_inferred = if (is.null(sex_inferred)) rep(NA_character_, nS) else unname(sex_inferred[all_samples]),
    reason = if (is.null(reasons)) rep("", nS) else unname(reasons[all_samples]),
    stringsAsFactors = FALSE, row.names = NULL)
  if (!is.null(sample_df) && "sample_name" %in% names(sample_df)) {
    m <- match(all_samples, as.character(sample_df$sample_name))
    for (cn in utils::head(setdiff(colnames(sample_df), "sample_name"), 8L)) {
      nm <- sub("^gender$", "Gender", cn); v <- sample_df[[cn]][m]
      df[[nm]] <- if (is.numeric(v)) v else ifelse(is.na(v), "", as.character(v))
    }
  }
  df
}

# ---------------------------------------------------------------------------
# Chromosome tables
# ---------------------------------------------------------------------------

#' Per-chromosome table of a target vector (rows in BED order)
#' @param analysed_class Chromosome class analysed by the run (\code{"autosome"}, \code{"X"}, \code{"Y"}).
#' @noRd
report_chrom_table <- function(chr_vec, analysed_class = "autosome") {
  chr_vec <- as.character(chr_vec)
  if (!length(chr_vec))
    return(data.frame(chromosome = character(), targets = integer(), first = integer(), last = integer(),
                      class = character(), analysed = logical(), stringsAsFactors = FALSE))
  ch <- unique(chr_vec); ii <- split(seq_along(chr_vec), factor(chr_vec, levels = ch))
  out <- data.frame(chromosome = ch, targets = as.integer(lengths(ii)),
                    first = vapply(ii, min, integer(1)), last = vapply(ii, max, integer(1)),
                    class = chromosome_class(ch), stringsAsFactors = FALSE, row.names = NULL)
  out$analysed <- out$class %in% analysed_class
  out
}

#' Targets per chromosome: in the input BED vs reaching the caller
#'
#' A chromosome that is in the BED but has no analysed target (\code{lost}) is
#' the signature of a silent upstream failure and is reported loudly.
#' @param input_chr Chromosome of every INPUT target (before any drop / filter).
#' @param analysed_chr Chromosome of the targets that were called.
#' @param analysed_class Classes the run analyses.
#' @noRd
report_chrom_inventory <- function(input_chr, analysed_chr, analysed_class = "autosome") {
  input_chr <- as.character(input_chr); analysed_chr <- as.character(analysed_chr)
  ch <- unique(input_chr)
  ch <- ch[order(factor(ch, levels = chrom_sort_levels(ch)))]
  cls <- chromosome_class(ch)
  n_in  <- vapply(ch, function(k) sum(input_chr == k), numeric(1))
  n_ana <- vapply(ch, function(k) sum(analysed_chr == k), numeric(1))
  in_scope <- cls %in% analysed_class
  data.frame(chromosome = ch, class = cls, targets_in_bed = as.integer(n_in),
             targets_analysed = as.integer(n_ana), analysed = in_scope & n_ana > 0,
             lost = in_scope & n_ana == 0, stringsAsFactors = FALSE, row.names = NULL)
}

# ---------------------------------------------------------------------------
# Alerts
# ---------------------------------------------------------------------------

#' The alert list shown in the header banner and on the Overview tab
#'
#' @param qc QC flag table (\code{Sample}, \code{Exon}, \code{Type}, \code{Details},
#'   \code{Gene}, \code{Severity}).
#' @param scope_rows Output of \code{chromosome_scope_warnings()}.
#' @param calls Standardised calls (see \code{report_sample_table}) with
#'   \code{calibration} (\code{"check"}/\code{"ok"}) and \code{recurrent}.
#' @param extra Additional pipeline-specific alerts (list).
#' @noRd
report_alerts <- function(qc, scope_rows, analysed, all_samples, calls = NULL, log_info = NULL,
                          qc_available = TRUE, extra = NULL, dropped_intervals = 0L,
                          recurrence = NULL, lost_chromosomes = character()) {
  out <- list()
  add <- function(kind, severity, title, text, tab = "qc", samples = NULL)
    out[[length(out) + 1L]] <<- list(kind = kind, severity = severity, title = title, text = text, tab = tab,
                                     samples = if (length(samples)) rep_arr(as.character(samples)) else NULL)

  for (i in seq_len(nrow(scope_rows))) {
    d <- scope_rows$Details[i]
    title <- if (scope_rows$Type[i] == "Missing Chromosome") "Chromosome class missing from the BED file"
      else if (grepl("^Autosomes", d)) "Autosomes are NOT analysed in this run"
      else if (grepl("^Gonosomes", d)) "Gonosomes (chrX / chrY) are NOT analysed in this run"
      else "Non-canonical contigs"
    add("scope", scope_rows$Severity[i], title, d, "qc")
  }
  if (length(lost_chromosomes))
    add("scope", "error", sprintf("%d chromosome(s) in the BED have NO analysed target", length(lost_chromosomes)),
        paste0(paste(lost_chromosomes, collapse = ", "), ": every interval was dropped or removed by QC, so no CNV can be called there. ",
               "See 'Target removed' in the QC tab for the reason."), "qc")
  if (dropped_intervals > 0L)
    add("scope", "warning", sprintf("%d BED interval(s) dropped before counting", dropped_intervals),
        "Intervals that could not be extracted (unknown contig, coordinates outside the reference, unreadable sequence) are listed under 'Target removed' in the QC tab.", "qc")
  k <- which(qc$Type == "Cohort size")
  for (i in k) add("cohort", qc$Severity[i], "Small cohort", qc$Details[i], "qc")
  if (!qc_available) {
    add("qc", "info", "No QC metrics file", "The QC metrics TSV was not found.", "qc")
    if (length(analysed) < 3 && !length(k))
      add("cohort", "warning", "Small cohort", sprintf("Only %d sample(s).", length(analysed)), "qc")
  }
  excl <- setdiff(all_samples, analysed)
  if (length(excl)) {
    txt <- paste(utils::head(excl, 10), collapse = ", "); if (length(excl) > 10) txt <- paste0(txt, ", ...")
    add("sample", "fail", paste0(length(excl), " sample(s) excluded"), txt, "samples", excl)
  }
  ws <- qc[qc$Type == "Whole sample" & qc$Sample %in% all_samples & !qc$Sample %in% excl, , drop = FALSE]
  if (nrow(ws)) {
    kind <- sub(":.*$", "", ws$Details)
    for (kd in unique(kind)) {
      rows <- ws[kind == kd, , drop = FALSE]
      add("sample", if (any(rows$Severity %in% c("error", "fail"))) "error" else "warning",
          paste0(kd, " in ", length(unique(rows$Sample)), " sample(s)"),
          paste(unique(rows$Sample), collapse = ", "), "qc", unique(rows$Sample))
    }
  }
  nt <- sum(qc$Type == "Target removed")
  if (nt) add("qc", "info", paste0(nt, " target(s) removed before calling"),
              "Low mean count, high variability or extreme GC (see 'Target removed' in the QC tab).", "qc")
  ne <- sum(qc$Type == "Whole exon"); if (ne) add("qc", "info", paste0(ne, " target(s) with low depth"), "", "qc")
  nv <- sum(qc$Type == "Exon variability"); if (nv) add("qc", "info", paste0(nv, " target(s) with high variability"), "", "qc")

  if (is.null(calls) || !nrow(calls)) {
    add("calls", "info", "No CNV called", "", "cnvs")
  } else {
    chk <- calls$calibration %in% "check"
    if (any(chk)) add("calls", "warning", paste0(sum(chk), " call(s) calibration check"), "", "cnvs", unique(calls$sample[chk]))
    rc <- calls$recurrent %in% TRUE
    if (any(rc)) add("calls", "warning", paste0(sum(rc), " recurrent call(s)"),
                     sprintf("Same call in >= %s sample(s): a hallmark of a noisy target rather than a germline CNV%s.",
                             recurrence$min_samples %||% 3,
                             if (identical(recurrence$action, "downgrade")) " (confidence downgraded to LOW)" else ""),
                     "cnvs", unique(calls$sample[rc]))
  }
  for (a in extra) out[[length(out) + 1L]] <- a
  nw <- length(log_info$warnings)
  if (nw) add("run", "info", paste0(nw, " log warning(s)"), "", "run")
  out
}

#' Alerts in the legacy "warnings" shape (level / category / text / sample)
#' @noRd
report_alerts_to_warnings <- function(alerts) {
  lapply(alerts, function(a) list(
    level = switch(a$severity, error = "fail", fail = "fail", warning = "warn", warn = "warn", "info"),
    category = a$kind, text = paste0(a$title, if (nzchar(a$text %||% "")) paste0(": ", a$text) else ""),
    sample = if (length(a$samples)) a$samples[1] else NULL))
}

# ---------------------------------------------------------------------------
# Settings / config / per-sample restriction
# ---------------------------------------------------------------------------

report_settings_kv <- function(settings, modechrom = NULL) {
  keep <- settings[vapply(settings, function(v) is.atomic(v) && !is.null(v), logical(1))]
  vals <- lapply(keep, function(v) paste(format(v, scientific = FALSE, trim = TRUE), collapse = ", "))
  if (!is.null(modechrom) && is.null(keep$modechrom)) vals <- c(list(modechrom = modechrom), vals)
  vals
}

report_config_yaml <- function(config) {
  if (is.null(config)) return(NULL)
  tryCatch(paste(utils::head(strsplit(yaml::as.yaml(config), "\n", fixed = TRUE)[[1]], 400), collapse = "\n"),
           error = function(e) NULL)
}

# Per-target QC rows (low depth, high variability, removed targets) can number in
# the tens of thousands on a whole-exome panel; the report keeps the first
# `max_rows` of them (in file order) and always keeps every sample-level row.
CNV_QC_TARGET_TYPES <- c("Whole exon", "Exon variability", "Target removed")
report_cap_qc <- function(qc, max_rows = 5000L) {
  is_t <- qc$Type %in% CNV_QC_TARGET_TYPES
  if (sum(is_t) <= max_rows) return(list(qc = qc, total = nrow(qc), capped = FALSE))
  keep <- !is_t; keep[utils::head(which(is_t), max_rows)] <- TRUE
  list(qc = qc[keep, , drop = FALSE], total = nrow(qc), capped = TRUE)
}

# A file-system-safe, collision-free name for a sample: sanitised prefix + a short
# hash of the full name (two samples that only differ by characters removed in the
# sanitised prefix still get different folders).
report_safe_name <- function(s, width = 80L) {
  u <- utf8ToInt(enc2utf8(as.character(s)))
  h <- sprintf("%08x", as.integer(sum(u * seq_along(u)) %% 2147483647))
  paste0(substr(sanitize_filename(s), 1L, width), "_", h)
}

# Write one global report, or one report per sample, for either pipeline.
#   build(focus)  -> payload (focus = NULL for the global one)
#   global_name   -> file name of the global report
#   sample_name(safe) -> file name of the report of one sample (written in its own folder)
report_emit <- function(build, output_dir, global_name, sample_name, samples = NULL,
                        settings = list(), global = TRUE) {
  if (global) {
    out <- file.path(output_dir, global_name)
    report_html(build(NULL), settings, out)
    return(invisible(out))
  }
  if (!length(samples)) { message("[INFO] No samples to report on -- skipping per-sample reports."); return(invisible(NULL)) }
  out <- character()
  for (s in samples) {
    safe <- report_safe_name(s)
    dir <- file.path(output_dir, safe); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    f <- file.path(dir, sample_name(safe))
    report_html(build(s), settings, f)
    out <- c(out, f)
  }
  invisible(out)
}

restrict_payload_to_sample <- function(p, s) {
  rows <- function(tb, keep) { tb$data <- lapply(tb$data, function(v) I(unname(v[keep]))); tb$n <- sum(keep); tb }
  p$samples  <- rows(p$samples, p$samples$data$sample == s)
  p$calls    <- rows(p$calls, p$calls$data$sample == s)
  p$qc_flags <- rows(p$qc_flags, p$qc_flags$data$Sample %in% c(s, "All"))
  ids <- as.character(p$calls$data$id)
  p$plots <- p$plots[names(p$plots) %in% ids]; p$windows <- p$plots
  p$refs <- p$refs[names(p$refs) == s]
  keep_alert <- function(a) a$kind != "run" && (is.null(a$samples) || s %in% a$samples)
  trim_alert <- function(a) { if (length(a$samples) > 1) { a$samples <- I(s); a$text <- "See QC and warnings tab." }; a }
  p$alerts <- lapply(Filter(keep_alert, p$alerts), trim_alert)
  p$warnings <- report_alerts_to_warnings(p$alerts)
  p$pca <- NULL; p$corr <- NULL; p$sample_table <- NULL; p$ref_bams <- NULL; p$config_yaml <- NULL
  nC <- p$calls$n
  p$meta[c("n_samples", "n_analysed", "samples_requested")] <- 1L
  p$meta$n_calls <- nC
  p$meta$samples_called <- p$meta$samples_with_calls <- as.integer(nC > 0L)
  p$meta$focus <- s
  if (is.list(p$summary)) {
    for (nm in intersect(c("samples_requested", "samples_analysed", "samples_called", "samples_with_calls",
                           "n_samples", "n_analysed"), names(p$summary)))
      p$summary[[nm]] <- if (nm %in% c("samples_called", "samples_with_calls")) as.integer(nC > 0L) else 1L
    if ("calls" %in% names(p$summary))   p$summary$calls <- nC
    if ("n_calls" %in% names(p$summary)) p$summary$n_calls <- nC
  }
  p
}

# ---------------------------------------------------------------------------
# Report inputs read the same way by both pipelines
# ---------------------------------------------------------------------------

#' Read the QC-metrics TSV of a run into the six-column QC table (or NULL)
#' @noRd
report_read_qc <- function(qc_metrics_file) {
  if (is.null(qc_metrics_file) || !is.character(qc_metrics_file) || !file.exists(qc_metrics_file) ||
      file.size(qc_metrics_file) == 0) return(NULL)
  qc <- tryCatch(data.table::fread(qc_metrics_file, header = TRUE, fill = TRUE, data.table = FALSE, colClasses = "character"),
                 error = function(e) NULL)
  if (is.null(qc) || !all(c("Sample", "Type", "Details") %in% colnames(qc))) return(NULL)
  for (cn in c("Exon", "Gene")) if (!cn %in% colnames(qc)) qc[[cn]] <- ""
  if (!"Severity" %in% colnames(qc))
    qc$Severity <- ifelse(qc$Type %in% c("Whole exon", "Exon variability", "Target removed"), "info", "warning")
  qc <- qc[, c("Sample", "Exon", "Type", "Details", "Gene", "Severity"), drop = FALSE]
  qc[is.na(qc)] <- ""
  qc
}

#' Read a tab-separated side table (sample sheet, reference BAM list); NULL when unusable
#' @param required Column that must be present.
#' @noRd
report_read_table <- function(path, required) {
  if (is.null(path) || !is.character(path) || !length(path) || !file.exists(path)) return(NULL)
  tb <- tryCatch(utils::read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE,
                                   quote = "", comment.char = ""), error = function(e) NULL)
  if (is.null(tb) || !required %in% colnames(tb)) NULL else tb
}

#' QC rows ("Whole sample") for the sex-check mismatches, in the shape of the QC table
#' @noRd
report_sex_flags <- function(sx, analysed, modechrom) {
  if (!is.data.frame(sx) || !nrow(sx)) return(NULL)
  row <- function(s, d) data.frame(Sample = s, Exon = "All", Type = "Whole sample", Details = d, Gene = "All",
                                   Severity = "warning", stringsAsFactors = FALSE)
  mm <- sx[sx$mismatch %in% TRUE, , drop = FALSE]
  wm <- sx[sx$mode_mismatch %in% TRUE & sx$sample %in% analysed, , drop = FALSE]
  out <- c(lapply(seq_len(nrow(mm)), function(i)
             row(mm$sample[i], sprintf("Sex mismatch: annotated %s, chrX/chrY depth looks %s", mm$annotated_sex[i], mm$inferred_sex[i]))),
           lapply(seq_len(nrow(wm)), function(i)
             row(wm$sample[i], sprintf("Sex vs modechrom: looks %s but this run uses modechrom = '%s'", wm$inferred_sex[i], modechrom))))
  if (length(out)) do.call(rbind, out) else NULL
}

#' Chromosome-depth matrix (rows chromosomes, columns samples) in the report shape
#' @noRd
report_chrom_depth <- function(chrom_depth) {
  if (is.null(chrom_depth) || is.null(dim(chrom_depth)) || !ncol(chrom_depth)) return(NULL)
  chrom_depth <- chrom_depth[order(factor(rownames(chrom_depth), levels = chrom_sort_levels(rownames(chrom_depth)))), , drop = FALSE]
  list(chromosomes = rep_arr(rownames(chrom_depth)), samples = rep_arr(colnames(chrom_depth)),
       z = lapply(seq_len(nrow(chrom_depth)), function(i) rep_arr(round(as.numeric(chrom_depth[i, ]), 3))))
}

# ---------------------------------------------------------------------------
# One payload builder for both pipelines
# ---------------------------------------------------------------------------

#' Top-level keys of the unified report payload (the shared JS reads exactly these)
#' @noRd
REPORT_PAYLOAD_KEYS <- c("meta", "thresholds", "conf_rules", "params", "summary", "samples", "flags",
                         "sample_attrs", "warnings", "alerts", "scope", "chr", "chr_all", "chromosomes",
                         "qc", "qc_flags", "corr", "pca", "depth_hist", "calls", "plots", "windows",
                         "refs", "settings", "scoring", "sample_table", "ref_bams", "log", "config_yaml")

#' Serialise a plot window (\code{pd} of either engine) for the HTML report
#' @noRd
window_from_plot_data <- function(pd) {
  r <- function(x, d = 3) I(round(unname(x), d))
  list(flag = isTRUE(pd$bg_calibration$flag),
       window = list(px = r(pd$px), xlab = I(pd$xlab), gene = I(pd$gene), off = I(pd$off),
                     gene_group = I(as.integer(pd$gene_group)), affected = I(pd$affected),
                     test_log = r(pd$test_log), ref_names = I(names(pd$ref_log)), ref_cov = unname(lapply(pd$ref_log, r)),
                     ratio = r(pd$ratio), lo = r(pd$lo), hi = r(pd$hi),
                     observed = r(pd$observed, 1), expected = r(pd$expected, 1),
                     z = r(pd$z), ref_z = unname(lapply(pd$ref_z, r)),
                     baseline = pd$baseline %||% 1, labels = pd$labels,
                     subtitle = pd$ci_subtitle))
}

#' Build the unified report payload (everything both pipelines have in common)
#'
#' The pipeline-specific parts are supplied by an \emph{adapter} \code{ad}, a
#' list made by \code{echo_report_adapter()} / \code{canope_report_adapter()}:
#' \describe{
#'   \item{\code{source, pkg, title}}{Names used in the meta block.}
#'   \item{\code{meta}}{Extra meta fields (engine, genome build, mode label...).}
#'   \item{\code{modechrom, analysed_class}}{The run's chromosome mode and class.}
#'   \item{\code{scope(input_chr, analysed_chr), scope_hints}}{Chromosome scope.}
#'   \item{\code{thresholds, conf_rules, scoring}}{Rules shown in the report.}
#'   \item{\code{cols}}{Names of the \code{sample, type, confidence, id} columns of the call table.}
#'   \item{\code{window(calls, i)}}{Plot window of call \code{i} (see \code{window_from_plot_data}).}
#'   \item{\code{call_table(calls, call_type, calib, has_plot)}}{The lowercase unified call table.}
#'   \item{\code{order(calls)}}{Priority of the calls whose window is embedded.}
#'   \item{\code{ref_cor, n_targets, failure_msg, extra_qc, extra_alerts, chrom_depth, sex_check, settings_kv, params}}{Optional; \\code{chrom_depth} / \\code{sex_check} are the raw matrix and table of \\code{cnv_chrom_depth()} / \\code{cnv_sex_check()}.}
#' }
#' @param rd Inputs: \code{cnv_calls, counts, bed_file, models, refs, analysed,
#'   all_samples, all_chromosomes, qc, sample_df, ref_bams, qc_info, bed_dropped,
#'   cohort} (optional raw-count matrix used for the cohort views).
#' @noRd
report_build_payload <- function(rd, settings, config, log_file, prefix, focus = NULL, ad) {
  qi <- rd$qc_info
  counts <- rd$counts; bed_file <- rd$bed_file; calls <- rd$cnv_calls
  has_calls <- is.data.frame(calls) && nrow(calls) > 0
  all_samples <- rd$all_samples; analysed <- rd$analysed
  rec_cfg <- qi$recurrence
  thresholds <- ad$thresholds
  cc <- ad$cols

  # ---- chromosome scope: judged on the INPUT chromosomes (incl. dropped intervals) ----
  input_chr  <- as.character(rd$all_chromosomes %||% bed_file$chromosome)
  scope      <- ad$scope(input_chr, bed_file$chromosome)
  scope_rows <- chromosome_scope_warnings(scope, "modechrom", ad$scope_hints)
  inv        <- report_chrom_inventory(input_chr, bed_file$chromosome, ad$analysed_class)
  chr_df     <- report_chrom_table(bed_file$chromosome, ad$analysed_class)
  chr_all_df <- report_chrom_table(input_chr, ad$analysed_class)

  qc <- rd$qc
  if (is.null(qc)) qc <- data.frame(Sample = character(), Exon = character(), Type = character(),
                                    Details = character(), Gene = character(), Severity = character(), stringsAsFactors = FALSE)
  if (nrow(scope_rows)) {
    add <- scope_rows[!paste(scope_rows$Type, scope_rows$Details) %in% paste(qc$Type, qc$Details), , drop = FALSE]
    if (nrow(add)) qc <- rbind(qc, data.frame(Sample = "All", Exon = "All", Type = add$Type, Details = add$Details,
                                              Gene = "All", Severity = add$Severity, stringsAsFactors = FALSE))
  }
  if (!is.null(ad$extra_qc) && nrow(ad$extra_qc)) {
    ex <- ad$extra_qc[!paste(ad$extra_qc$Sample, ad$extra_qc$Type, ad$extra_qc$Details) %in% paste(qc$Sample, qc$Type, qc$Details), , drop = FALSE]
    if (nrow(ex)) qc <- rbind(qc, ex[, names(qc), drop = FALSE])
  }
  sx <- ad$sex_check
  sx_rows <- report_sex_flags(sx, rd$analysed, ad$modechrom)
  if (!is.null(sx_rows)) {
    sx_rows <- sx_rows[!paste(sx_rows$Sample, sx_rows$Type, sx_rows$Details) %in% paste(qc$Sample, qc$Type, qc$Details), , drop = FALSE]
    if (nrow(sx_rows)) qc <- rbind(qc, sx_rows[, names(qc), drop = FALSE])
  }
  if (is.null(ad$sex_inferred) && is.data.frame(sx) && nrow(sx)) ad$sex_inferred <- stats::setNames(as.character(sx$inferred_sex), sx$sample)
  qc$Severity[!qc$Severity %in% c("error", "warning", "info", "fail")] <- "warning"
  log_info <- report_log(log_file)

  # ---- samples ----
  mat <- rd$cohort
  if (is.null(mat)) mat <- as.matrix(counts[, all_samples, drop = FALSE])
  mat <- mat[, all_samples, drop = FALSE]; storage.mode(mat) <- "double"
  models <- rd$models %||% list(); refs <- rd$refs %||% list()
  n_refs <- vapply(all_samples, function(s) length(refs[[s]] %||% character()), numeric(1))
  # a sample that was meant to be called but has no model / no reference set FAILED
  # (it is not the same as "no CNV found")
  failed <- stats::setNames(rep(NA_character_, length(all_samples)), all_samples)
  failed[analysed[vapply(analysed, function(s) is.null(models[[s]]) || n_refs[[s]] == 0, logical(1))]] <- "failed"
  is_failed <- failed[all_samples] %in% "failed"
  fail_msg <- stats::setNames(ifelse(is_failed, "No reference set could be built", NA_character_), all_samples)
  # the pipeline's own reason (e.g. the error message of the caller) wins
  given <- ad$failure_msg[intersect(names(ad$failure_msg), all_samples[is_failed])]
  fail_msg[names(given)] <- unlist(given)
  excl_df <- qi$excluded_samples
  sample_metrics <- if (is.null(qi$sample_metrics)) NULL else qi$sample_metrics[qi$sample_metrics$sample %in% all_samples, , drop = FALSE]

  # call-level view shared by the sample table and the alerts
  call_type <- if (has_calls) { tt <- tolower(as.character(calls[[cc$type]])); ifelse(grepl("^del", tt), "DEL", ifelse(grepl("^dup", tt), "DUP", toupper(tt))) } else character()
  rec <- if (has_calls && "Recurrent" %in% names(calls)) calls$Recurrent %in% TRUE else rep(FALSE, if (has_calls) nrow(calls) else 0L)
  call_view <- if (has_calls) data.frame(sample = as.character(calls[[cc$sample]]), type = call_type,
                                         confidence = as.character(calls[[cc$confidence]]), recurrent = rec, stringsAsFactors = FALSE) else NULL
  sdf <- rd$sample_df
  base_tab <- report_sample_table(mat, analysed, calls = call_view, n_refs = n_refs, status = failed, qc = qc,
                                  qc_zscore = thresholds$qc_zscore, outlier_sides = thresholds$outlier_sides,
                                  sample_df = NULL, metrics = sample_metrics)
  reasons <- report_exclusion_reasons(all_samples, excl_df, log_info, base_tab, thresholds, status_msg = fail_msg)
  # the column answers "why was this sample NOT analysed": empty for the samples that were
  reasons[all_samples %in% analysed & !is_failed] <- ""
  n_tg <- ad$n_targets %||% stats::setNames(rep(nrow(counts), length(all_samples)), all_samples)
  samples_df <- report_sample_table(mat, analysed, calls = call_view, n_refs = n_refs, status = failed, qc = qc,
                                    qc_zscore = thresholds$qc_zscore, outlier_sides = thresholds$outlier_sides,
                                    sample_df = sdf, reasons = reasons, metrics = sample_metrics,
                                    ref_cor = ad$ref_cor, n_targets = n_tg,
                                    sex_inferred = ad$sex_inferred, extra_level = ad$extra_level)
  n_mode <- qi$n_targets$mode %||% nrow(counts)
  samples_df$n_dropped <- ifelse(is.na(samples_df$n_targets), NA_real_, n_mode - samples_df$n_targets)

  # why each excluded sample is a failure, as a QC row too
  excl_samples <- all_samples[!all_samples %in% analysed]
  if (length(excl_samples)) {
    new_rows <- data.frame(Sample = excl_samples, Exon = "All", Type = "Whole sample",
                           Details = sprintf("Not analysed: %s", reasons[excl_samples]), Gene = "All",
                           Severity = "fail", stringsAsFactors = FALSE)
    new_rows <- new_rows[!paste(new_rows$Sample, new_rows$Type, new_rows$Details) %in% paste(qc$Sample, qc$Type, qc$Details), , drop = FALSE]
    if (nrow(new_rows)) qc <- rbind(qc, new_rows)
  }

  # ---- CNV table and per-call plot windows ----
  max_profiles <- as.integer(settings$report_max_profiles %||% 500L)
  windows <- list(); calib <- rep("", if (has_calls) nrow(calls) else 0L)
  cnv_df <- NULL; windows_limited <- FALSE
  if (has_calls) {
    todo <- if (is.function(ad$order)) ad$order(calls) else
      order(match(calls[[cc$confidence]], c("HIGH", "MEDIUM", "LOW"), nomatch = 4L), calls[[cc$id]])
    if (length(todo) > max_profiles) { windows_limited <- TRUE; todo <- todo[seq_len(max_profiles)] }
    for (i in todo) {
      w <- tryCatch(ad$window(calls, i), error = function(e) NULL)
      if (!is.null(w)) { windows[[as.character(calls[[cc$id]][i])]] <- w$window; calib[i] <- if (isTRUE(w$flag)) "check" else "ok" }
    }
    has_plot <- as.character(calls[[cc$id]]) %in% names(windows)
    cnv_df <- ad$call_table(calls, call_type, calib, has_plot)
    cnv_df$gene[is.na(cnv_df$gene) | !nzchar(cnv_df$gene)] <- "Unknown"
    # flags every pipeline shows the same way (the adapter may supply its own)
    low_genes <- toupper(as.character(unlist(ad$conf_rules$genes)))
    if (!"homology_gene" %in% names(cnv_df))
      cnv_df$homology_gene <- vapply(strsplit(toupper(as.character(cnv_df$gene)), "[;, ]+"),
                                     function(g) any(g %in% low_genes), logical(1))
    if (!"bg_flag" %in% names(cnv_df)) cnv_df$bg_flag <- cnv_df$calibration %in% "check"
    if (!"terminal_only" %in% names(cnv_df))
      cnv_df$terminal_only <- if ("terminal_only" %in% names(calls)) calls$terminal_only %in% TRUE else rep(FALSE, nrow(cnv_df))
    # reference quality of the test sample = median correlation over its calls (when the
    # pipeline gave no per-sample value)
    if (is.null(ad$ref_cor) && "correlation" %in% names(cnv_df)) {
      rc <- tapply(as.numeric(cnv_df$correlation), as.character(cnv_df$sample), function(x) stats::median(x, na.rm = TRUE))
      samples_df$ref_cor <- round(as.numeric(rc[samples_df$sample]), 3)
    }
  }

  # ---- cohort views ----
  cm <- report_correlation(mat)
  grp <- if (!is.null(sdf) && "gender" %in% colnames(sdf)) as.character(sdf$gender[match(all_samples, as.character(sdf$sample_name))]) else NULL
  pca <- report_pca(mat, grp)
  corr <- report_corr(cm$matrix, 60L)
  depth_hist <- report_depth_hist(mat)

  # ---- alerts, settings ----
  lost <- inv$chromosome[inv$lost]
  n_dropped_iv <- if (is.null(rd$bed_dropped)) 0L else nrow(rd$bed_dropped)
  alert_calls <- if (is.null(cnv_df)) NULL else data.frame(sample = cnv_df$sample, calibration = cnv_df$calibration,
                                                           recurrent = cnv_df$recurrent, stringsAsFactors = FALSE)
  extra <- c(if (any(is_failed))
               list(list(kind = "sample", severity = "fail",
                         title = sprintf("%d sample(s) FAILED", sum(is_failed)),
                         text = paste(all_samples[is_failed], collapse = ", "),
                         tab = "samples", samples = rep_arr(all_samples[is_failed]))),
             ad$extra_alerts)
  alerts <- report_alerts(qc, scope_rows, analysed, all_samples, alert_calls, log_info,
                          qc_available = !is.null(rd$qc), dropped_intervals = n_dropped_iv,
                          recurrence = rec_cfg, lost_chromosomes = lost, extra = extra)
  warns <- report_alerts_to_warnings(alerts)
  settings_kv <- ad$settings_kv %||% report_settings_kv(settings, ad$modechrom)

  chr_status_df <- data.frame(class = scope$class, label = scope$label, in_input = scope$present, analysed = scope$analysed,
                              reason = scope$reason, stringsAsFactors = FALSE)
  chr_inv_df <- inv[, c("chromosome", "class", "targets_in_bed", "targets_analysed", "analysed", "lost")]
  chr_extent <- list(chr = I(unname(as.character(chr_df$chromosome))),
                     len = I(unname(pmax(as.numeric(chr_df$targets) - 1, 0))))
  scope_incomplete <- any(scope$class %in% c("autosome", "X", "Y") & !scope$analysed) || length(lost) > 0L

  meta <- utils::modifyList(list(
    source = ad$source,
    title = paste0(ad$title, if (!is.null(prefix) && nzchar(prefix)) paste0(" - ", prefix) else "", if (!is.null(focus)) paste0(" - ", focus) else ""),
    prefix = prefix %||% "", generated = format(Sys.time(), "%Y-%m-%d %H:%M"),
    version = tryCatch(as.character(utils::packageVersion(ad$pkg)), error = function(e) ""),
    legacy_workspace = is.null(qi), palette = cnv_palette(),
    focus = if (is.null(focus)) NULL else focus,
    n_targets = nrow(bed_file), n_targets_input = length(input_chr), n_samples = length(all_samples), n_analysed = length(analysed),
    n_calls = if (has_calls) nrow(calls) else 0L,
    scope_incomplete = scope_incomplete, qc_available = !is.null(rd$qc),
    windows_limited = windows_limited, windows_limit = max_profiles), ad$meta %||% list())

  qc_cap <- report_cap_qc(qc, as.integer(settings$report_max_qc_rows %||% 5000L))
  meta$qc_capped <- qc_cap$capped; meta$qc_rows_total <- qc_cap$total

  payload <- list(
    meta = meta, thresholds = thresholds, conf_rules = ad$conf_rules, params = ad$params %||% as.list(settings),
    summary = list(
      samples_requested = length(all_samples), samples_called = sum(samples_df$status == "called"),
      samples_with_calls = length(unique(if (has_calls) calls[[cc$sample]] else character())),
      calls = if (has_calls) nrow(calls) else 0L, targets = nrow(bed_file),
      both_chromosome_classes = !scope_incomplete),
    samples = rep_table(samples_df), flags = list(), sample_attrs = list(),
    warnings = warns, alerts = alerts,
    scope = rep_table(scope), chr = rep_table(chr_df), chr_all = rep_table(chr_all_df),
    chromosomes = list(status = rep_table(chr_status_df), messages = I(character()), inventory = rep_table(chr_inv_df),
                       depth = report_chrom_depth(ad$chrom_depth), extent = chr_extent,
                       sex = if (is.data.frame(sx) && nrow(sx)) rep_table(sx) else NULL),
    qc = NULL, qc_flags = rep_table(qc_cap$qc),
    corr = corr, pca = pca, depth_hist = depth_hist,
    calls = if (is.null(cnv_df)) rep_empty_table(c("id", "sample", "chr", "start", "end", "gene", "type", "confidence", "has_plot")) else rep_table(cnv_df),
    plots = windows, windows = windows,
    refs = lapply(refs[all_samples[all_samples %in% names(refs)]], function(x) I(as.character(x))),
    settings = lapply(seq_along(settings_kv), function(i) I(c(names(settings_kv)[i], settings_kv[[i]]))),
    scoring = ad$scoring,
    sample_table = if (is.null(sdf)) NULL else rep_table(utils::head(sdf, 2000)),
    ref_bams = if (is.null(rd$ref_bams)) NULL else rep_table(utils::head(rd$ref_bams, 2000)),
    log = log_info, config_yaml = report_config_yaml(config))

  if (!is.null(focus)) payload <- restrict_payload_to_sample(payload, focus)
  payload
}
