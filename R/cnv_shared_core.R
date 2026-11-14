# ============================================================================
# cnv_shared_core.R  --  code shared BYTE-FOR-BYTE by ECHO and CANOPE
#
# This file is identical in both packages (each package keeps its own copy and
# its own namespace). tests/testthat/test-shared-sync.R fails if the two copies
# ever drift apart, so a fix made here must be copied to the sibling package.
#
# Contents
#   1. small utilities           (%||%, sanitize_filename, load_rdata, ...)
#   2. BED / target helpers      (exon numbering, off-target, padding, labels)
#   3. chromosome scope          (classes, what a run did / did not analyse)
#   4. logging                   (timestamped file + console logger)
#   5. sample / target QC        (outlier detection, QC metrics table)
#   6. call post-processing      (terminal-only flag, cohort recurrence flag)
#   7. palette                   (one colour scheme for PDF plots and HTML)
#   8. VCF writer                (plain-text VCF 4.2 building blocks)
#
# Nothing in here may refer to a package name, an engine or a call-table
# column name that only one pipeline has.
# ============================================================================

#' @importFrom data.table := .N .I
NULL

# ---------------------------------------------------------------------------
# 1. Small utilities
# ---------------------------------------------------------------------------

#' Null coalescing operator
#' @param a Primary value.
#' @param b Fallback used when \code{a} is \code{NULL}.
#' @return \code{a} unless it is \code{NULL}, otherwise \code{b}.
#' @export
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Sanitize a string for use as a file or directory name
#'
#' Keeps letters, digits, \code{.}, \code{_} and \code{-}; everything else
#' becomes \code{_}. Leading dots are replaced so the result can never be
#' \code{"."}, \code{".."} or a hidden file.
#' @param name Character vector.
#' @return Character vector.
#' @export
sanitize_filename <- function(name) {
  x <- gsub("[^[:alnum:]._-]", "_", as.character(name))
  sub("^\\.+", "_", x)
}


#' Stop if a file does not exist
#' @param path File path.
#' @param msg Error message.
#' @export
stop_if_not_file <- function(path, msg) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !file.exists(path))
    stop(msg, call. = FALSE)
  invisible(path)
}

#' Load an RData file into an isolated environment
#' @param path Path to the .RData file.
#' @param required Optional character vector of object names that must exist.
#' @return A named list of the objects in the file.
#' @export
load_rdata <- function(path, required = NULL) {
  stop_if_not_file(path, paste0("[ERROR] RData file not found: ", path))
  env <- new.env()
  load(path, envir = env)
  if (!is.null(required)) {
    miss <- setdiff(required, ls(env))
    if (length(miss))
      stop("[ERROR] Missing variables in RData: ", paste(miss, collapse = ", "),
           "\nAvailable objects: ", paste(ls(env), collapse = ", "), call. = FALSE)
  }
  as.list(env)
}

#' Column / row medians (matrixStats when installed)
#' @param m Numeric matrix.
#' @noRd
col_medians <- function(m) {
  m <- as.matrix(m)
  if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::colMedians(m, na.rm = TRUE)
  else apply(m, 2, stats::median, na.rm = TRUE)
}
#' @noRd
row_medians <- function(m) {
  m <- as.matrix(m)
  if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::rowMedians(m, na.rm = TRUE)
  else apply(m, 1, stats::median, na.rm = TRUE)
}

#' Check that a BAM file has an index on disk
#' @param bam Path to a BAM file.
#' @return Logical.
#' @noRd
bam_has_index <- function(bam) {
  # Rsamtools accepts both "x.bam.bai" and "x.bai".
  file.exists(paste0(bam, ".bai")) ||
    file.exists(sub("\\.bam$", ".bai", bam, ignore.case = TRUE))
}

#' Normalise a chromosome name to a reference naming style (vectorised)
#' @param chr_vec Character vector of chromosome names.
#' @param ref_chromosomes Reference chromosome names (decides the "chr" prefix).
#' @return Character vector.
#' @noRd
normalize_chromosome_vec <- function(chr_vec, ref_chromosomes) {
  chr_vec <- as.character(chr_vec)
  if (any(grepl("^chr", ref_chromosomes))) {
    ifelse(grepl("^chr", chr_vec), chr_vec, paste0("chr", chr_vec))
  } else {
    sub("^chr", "", chr_vec)
  }
}

#' Round-trip safe JSON number formatter used by report builders
#' @noRd
fmt_num <- function(x, digits = 2L) {
  ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "n/a")
}

# ---------------------------------------------------------------------------
# 2. BED / target helpers
# ---------------------------------------------------------------------------

#' Column names of a BED-like data frame
#'
#' Understands the lower-case (\code{chromosome/start/end/gene}), capitalised
#' (\code{Chr/Start/End/Gene}) and upper-case (\code{CHROM/START/END/GENE})
#' spellings used by the two pipelines.
#' @param df Data frame.
#' @return Named list: \code{chrom}, \code{start}, \code{end}, \code{gene}.
#' @noRd
bed_cols <- function(df) {
  nm <- names(df)
  pick <- function(cands, default) {
    hit <- cands[cands %in% nm]
    if (length(hit)) hit[1L] else default
  }
  list(chrom = pick(c("chromosome", "Chromosome", "Chr", "CHROM", "CHR", "chrom"), "chromosome"),
       start = pick(c("start", "Start", "START"), "start"),
       end   = pick(c("end", "End", "END"), "end"),
       gene  = pick(c("gene", "Gene", "GENE"), "gene"))
}

#' Natural chromosome order: 1-22, X, Y, M (with and without "chr"), then the rest
#' @noRd
chrom_sort_levels <- function(chrs = character()) {
  base <- c(as.character(1:22), "X", "Y", "M")
  unique(c(paste0("chr", base), base, unique(as.character(chrs))))
}

#' Label given to retained off-target (HorsROI) targets
#' @noRd
OFF_TARGET_LABEL <- "HorsROI"

#' Which gene labels denote off-target (HorsROI) targets
#' @param genes Character vector of gene labels.
#' @param pattern Optional extra regex (the run's \code{off_target_pattern}).
#' @return Logical vector.
#' @noRd
is_off_target_gene <- function(genes, pattern = "^HorsROI") {
  g   <- as.character(genes)
  hit <- !is.na(g) & toupper(g) == toupper(OFF_TARGET_LABEL)
  if (length(pattern) == 1L && !is.na(pattern) && nzchar(pattern))
    hit <- hit | (!is.na(g) & grepl(pattern, g, perl = TRUE))
  hit
}

#' Genes covered by target index ranges
#' @param genes Gene label of every target.
#' @param lo,hi Target index ranges (one pair per call).
#' @param drop_unknown Also drop \code{"Unknown"} labels.
#' @return One space-separated gene string per range (\code{"Unknown"} if none).
#' @noRd
genes_in_targets <- function(genes, lo, hi, drop_unknown = FALSE) {
  n <- length(genes)
  vapply(seq_along(lo), function(j) {
    a <- suppressWarnings(as.integer(lo[j])); b <- suppressWarnings(as.integer(hi[j]))
    if (is.na(a) || is.na(b) || b < a || a > n || b < 1L) return("Unknown")
    g <- unique(as.character(genes[seq.int(max(1L, a), min(b, n))]))
    g <- g[!is.na(g) & nzchar(g) & (!drop_unknown | g != "Unknown")]
    if (length(g)) paste(g, collapse = " ") else "Unknown"
  }, character(1))
}

#' Reclassify off-target / filler BED intervals before exon numbering
#'
#' @param df BED-like data frame.
#' @param pattern Regex identifying off-target intervals (matched against
#'   \code{OriginalName} when present, else the gene column).
#' @param handling \code{"na"} (gene set to NA, interval kept), \code{"remove"}
#'   or \code{"merge"} (attach to the nearest real gene on the chromosome).
#' @param verbose Emit a message.
#' @export
handle_off_target_regions <- function(df, pattern = "^HorsROI",
                                      handling = c("na", "remove", "merge"),
                                      verbose = TRUE) {
  handling <- match.arg(handling)
  if (is.null(pattern) || !nzchar(pattern)) return(df)
  if (is.null(df) || nrow(df) == 0) return(df)
  bc <- bed_cols(df)
  chrom_col <- bc$chrom; start_col <- bc$start; end_col <- bc$end; gene_col <- bc$gene
  if (!gene_col %in% names(df)) return(df)

  if ("OriginalName" %in% names(df)) {
    name_vals <- as.character(df[["OriginalName"]])
    is_off <- !is.na(name_vals) & grepl(pattern, name_vals, perl = TRUE)
  } else {
    gene_vals <- as.character(df[[gene_col]])
    is_off <- !is.na(gene_vals) & grepl(pattern, gene_vals, perl = TRUE)
  }

  n_off <- sum(is_off)
  if (n_off == 0) return(df)

  if (handling == "remove") {
    out <- df[!is_off, , drop = FALSE]
    if (verbose) message(sprintf(
      "[INFO] handle_off_target_regions: removed %d off-target interval(s) matching /%s/.", n_off, pattern))
    return(out)
  }

  if (handling == "na") {
    df[[gene_col]][is_off] <- NA_character_
    if (verbose) message(sprintf(
      "[INFO] handle_off_target_regions: set gene = NA for %d off-target interval(s) matching /%s/ (kept, excluded from exon numbering).",
      n_off, pattern))
    return(df)
  }

  # handling == "merge"
  ord     <- order(df[[chrom_col]], df[[start_col]], df[[end_col]])
  n       <- length(ord)
  chrom_s <- as.character(df[[chrom_col]])[ord]
  start_s <- df[[start_col]][ord]
  end_s   <- df[[end_col]][ord]
  gene_s  <- as.character(df[[gene_col]])[ord]
  off_s   <- is_off[ord]
  prev_gene <- character(n); prev_end <- numeric(n)
  g <- NA_character_; e <- NA_real_; last_chrom <- NA_character_
  for (i in seq_len(n)) {
    if (!identical(chrom_s[i], last_chrom)) { g <- NA_character_; e <- NA_real_; last_chrom <- chrom_s[i] }
    prev_gene[i] <- g; prev_end[i] <- e
    if (!off_s[i]) { g <- gene_s[i]; e <- end_s[i] }
  }
  next_gene <- character(n); next_start <- numeric(n)
  g <- NA_character_; s <- NA_real_; last_chrom <- NA_character_
  for (i in rev(seq_len(n))) {
    if (!identical(chrom_s[i], last_chrom)) { g <- NA_character_; s <- NA_real_; last_chrom <- chrom_s[i] }
    next_gene[i] <- g; next_start[i] <- s
    if (!off_s[i]) { g <- gene_s[i]; s <- start_s[i] }
  }
  new_gene_s  <- gene_s
  n_merged    <- 0L
  merged_flag <- logical(n)
  for (i in which(off_s)) {
    has_prev <- !is.na(prev_gene[i])
    has_next <- !is.na(next_gene[i])
    if (has_prev && has_next) {
      d_prev <- start_s[i] - prev_end[i]
      d_next <- next_start[i] - end_s[i]
      new_gene_s[i] <- if (d_prev <= d_next) prev_gene[i] else next_gene[i]
      n_merged <- n_merged + 1L; merged_flag[i] <- TRUE
    } else if (has_prev) {
      new_gene_s[i] <- prev_gene[i]; n_merged <- n_merged + 1L; merged_flag[i] <- TRUE
    } else if (has_next) {
      new_gene_s[i] <- next_gene[i]; n_merged <- n_merged + 1L; merged_flag[i] <- TRUE
    }
  }
  df[[gene_col]][ord] <- new_gene_s
  # Drop the off-target marker so hors_roi_rows() no longer re-flags merged rows.
  if ("OriginalName" %in% names(df) && any(merged_flag))
    df[["OriginalName"]][ord[merged_flag]] <- new_gene_s[merged_flag]
  if (verbose) message(sprintf(
    "[INFO] handle_off_target_regions: merged %d/%d off-target interval(s) matching /%s/ into their nearest neighbouring gene (%d had no real gene on their chromosome to attach to).",
    n_merged, n_off, pattern, n_off - n_merged))
  df
}

#' Assign sequential exon numbers within each chromosome/gene locus
#'
#' Numbers exons 1..n per gene in genomic order, whatever the row order of the
#' input. Duplicate intervals share the number of their first occurrence. The
#' returned data frame has exactly \code{nrow(bed_file)} rows in the input
#' order, so \code{exon_number} can be bound back by position. Targets with no
#' real gene (NA, \code{""}, \code{"."}, off-target) get \code{NA}.
#' @param bed_file BED-like data frame.
#' @return \code{bed_file} with an integer \code{exon_number} column.
#' @export
assign_exon_numbers_per_gene <- function(bed_file) {
  bc <- bed_cols(bed_file)
  chrom_col <- bc$chrom; start_col <- bc$start; end_col <- bc$end; gene_col <- bc$gene
  stopifnot(all(c(chrom_col, start_col, end_col, gene_col) %in% names(bed_file)))
  n_in <- nrow(bed_file)
  dt <- data.table::as.data.table(bed_file)
  dt[, .orig_row := .I]
  data.table::setnames(dt, c(chrom_col, start_col, end_col, gene_col),
                       c("chrom", "start", "end", "gene"))
  dt[, start := as.numeric(start)]
  dt[, end   := as.numeric(end)]
  dt[, chrom := as.character(chrom)]
  dt[, gene  := as.character(gene)]
  dt[, .key := paste(chrom, start, end, gene, sep = "␟")]
  dup_rows <- duplicated(dt, by = ".key")
  if (any(dup_rows)) {
    warning(sprintf(
      "Found %d duplicate row(s); duplicates share the exon_number of their first occurrence.",
      sum(dup_rows)), immediate. = TRUE)
  }
  dt_unique <- dt[!dup_rows]
  dt_unique[, .chrom_fac := factor(chrom, levels = chrom_sort_levels(chrom))]
  data.table::setorder(dt_unique, .chrom_fac, start, end)
  dt_unique[, .chrom_fac := NULL]

  valid_gene_rows <- !is.na(dt_unique$gene) & !dt_unique$gene %in% c("", ".", "NA") &
    !is_off_target_gene(dt_unique$gene, pattern = NULL)
  dt_unique[, exon_number := NA_integer_]
  if (any(valid_gene_rows))
    dt_unique[valid_gene_rows, exon_number := seq_len(.N), by = .(chrom, gene)]
  dt[, exon_number := dt_unique$exon_number[match(.key, dt_unique$.key)]]
  data.table::setorder(dt, .orig_row)
  dt[, c(".orig_row", ".key") := NULL]
  data.table::setnames(dt, c("chrom", "start", "end", "gene"),
                       c(chrom_col, start_col, end_col, gene_col))
  out <- as.data.frame(dt)
  stopifnot(nrow(out) == n_in)
  out
}

#' Exon index used for plotting (existing \code{exon_number} or recomputed)
#' @param bed_file BED-like data frame.
#' @export
compute_exon_index <- function(bed_file) {
  if ("exon_number" %in% names(bed_file) && !all(is.na(bed_file$exon_number)))
    return(bed_file$exon_number)
  assign_exon_numbers_per_gene(bed_file)$exon_number
}

#' Human-readable label of every target
#'
#' \code{GENE:exon} for annotated targets, \code{chr:start-end} for the others
#' (off-target / unannotated) instead of a meaningless \code{NA:NA}. Both
#' pipelines use this in the QC table, the plots and the report.
#' @param bed_file BED-like data frame (an \code{exon_number} column is
#'   computed when absent).
#' @return Character vector, one label per row.
#' @export
target_labels <- function(bed_file) {
  bc <- bed_cols(bed_file)
  exon <- if ("exon_number" %in% names(bed_file)) bed_file$exon_number
          else assign_exon_numbers_per_gene(bed_file)$exon_number
  gene <- as.character(bed_file[[bc$gene]])
  ok   <- !is.na(gene) & !gene %in% c("", ".", "NA") & !is_off_target_gene(gene, NULL) & !is.na(exon)
  ifelse(ok, paste0(gene, ":", exon),
         paste0(bed_file[[bc$chrom]], ":", format(bed_file[[bc$start]], scientific = FALSE, trim = TRUE),
                "-", format(bed_file[[bc$end]], scientific = FALSE, trim = TRUE)))
}

#' Best overlap per query interval
#'
#' For every query range that overlaps at least one subject range, picks the
#' subject with the widest overlap (ties go to the first subject, so sort the
#' subjects for a reproducible result). Shared by the BED annotation of both
#' pipelines.
#' @param gr1,gr2 \code{GRanges} (query, subject).
#' @param hits Output of \code{GenomicRanges::findOverlaps(gr1, gr2)} (may be
#'   \code{NULL}).
#' @return List with integer vectors \code{q} (query index) and \code{s}
#'   (subject index), one entry per annotated query.
#' @export
select_best_hits <- function(gr1, gr2, hits) {
  if (is.null(hits) || length(hits) == 0L) return(list(q = integer(0), s = integer(0)))
  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  w  <- IRanges::width(GenomicRanges::pintersect(gr1[qh], gr2[sh]))
  d  <- data.frame(q = qh, s = sh, w = w)
  d  <- d[order(d$q, -d$w, d$s), , drop = FALSE]
  d  <- d[!duplicated(d$q), , drop = FALSE]
  list(q = d$q, s = d$s)
}

#' Pad the outer edge of each gene's first and last exon
#'
#' @param bed_file BED-like data frame.
#' @param padding Bases to add (0 disables).
#' @param chr_lengths Optional named vector of contig lengths (clamp at the end).
#' @param verbose Emit a message.
#' @export
pad_gene_terminal_exons <- function(bed_file, padding = 0, chr_lengths = NULL, verbose = TRUE) {
  if (is.null(padding) || length(padding) != 1 || is.na(padding) || padding <= 0)
    return(bed_file)
  padding <- as.integer(round(padding))
  bc <- bed_cols(bed_file)
  chrom_col <- bc$chrom; start_col <- bc$start; end_col <- bc$end; gene_col <- bc$gene
  stopifnot(all(c(chrom_col, start_col, end_col, gene_col) %in% names(bed_file)))
  n_in <- nrow(bed_file)
  if (n_in == 0) return(bed_file)
  numbered <- assign_exon_numbers_per_gene(bed_file)
  original_exon_number <- if ("exon_number" %in% names(bed_file)) bed_file$exon_number else NULL
  dt <- data.table::as.data.table(numbered)
  dt[, .orig_row := .I]
  data.table::setnames(dt, c(chrom_col, start_col, end_col, gene_col),
                       c("chrom", "start", "end", "gene"))
  dt[, chrom := as.character(chrom)]
  dt[, is_first := exon_number == 1L]
  dt[, is_last  := exon_number == max(exon_number), by = .(chrom, gene)]
  no_gene <- is.na(dt$exon_number) | is.na(dt$gene) | dt$gene %in% c("", ".", "Unknown")
  dt[no_gene, c("is_first", "is_last") := FALSE]
  dt[, .chrom_fac := factor(chrom, levels = chrom_sort_levels(chrom))]
  data.table::setorder(dt, .chrom_fac, start, end)
  n          <- nrow(dt)
  chrom_id   <- as.integer(dt$.chrom_fac)
  start_v    <- dt$start
  end_v      <- dt$end
  is_first_v <- dt$is_first
  is_last_v  <- dt$is_last
  right_extend <- integer(n)
  left_extend  <- integer(n)
  if (n > 1) {
    same_chr_pair <- chrom_id[-n] == chrom_id[-1]
    gap        <- pmax(start_v[-1] - end_v[-n] - 1L, 0L)
    want_left  <- ifelse(is_last_v[-n],  padding, 0L)
    want_right <- ifelse(is_first_v[-1], padding, 0L)
    demand     <- want_left + want_right
    grant_left  <- integer(n - 1L)
    grant_right <- integer(n - 1L)
    has_demand  <- same_chr_pair & demand > 0
    fits        <- has_demand & demand <= gap
    tight       <- has_demand & demand > gap
    grant_left[fits]  <- want_left[fits]
    grant_right[fits] <- want_right[fits]
    grant_left[tight]  <- as.integer(floor(gap[tight] * want_left[tight] / demand[tight]))
    grant_right[tight] <- gap[tight] - grant_left[tight]
    right_extend[-n] <- grant_left
    left_extend[-1]  <- grant_right
  }
  has_prev <- c(FALSE, if (n > 1) chrom_id[-1] == chrom_id[-n] else logical(0))
  has_next <- c(if (n > 1) chrom_id[-n] == chrom_id[-1] else logical(0), FALSE)
  no_prev_first <- is_first_v & !has_prev
  if (any(no_prev_first))
    left_extend[no_prev_first] <- pmin(padding, pmax(start_v[no_prev_first] - 1L, 0L))
  no_next_last <- is_last_v & !has_next
  if (any(no_next_last)) {
    chr_len_here <- if (!is.null(chr_lengths)) {
      unname(chr_lengths[as.character(dt$chrom[no_next_last])])
    } else rep(NA_real_, sum(no_next_last))
    avail <- ifelse(!is.na(chr_len_here), chr_len_here - end_v[no_next_last], 1e12)
    right_extend[no_next_last] <- as.integer(pmin(padding, pmax(avail, 0)))
  }
  n_padded_start <- sum(left_extend > 0L)
  n_padded_end   <- sum(right_extend > 0L)
  n_clamped      <- sum(is_first_v & left_extend < padding) + sum(is_last_v & right_extend < padding)
  dt[, start := start_v - left_extend]
  dt[, end   := end_v   + right_extend]
  data.table::setorder(dt, .orig_row)
  dt[, c(".orig_row", ".chrom_fac", "is_first", "is_last", "exon_number") := NULL]
  data.table::setnames(dt, c("chrom", "start", "end", "gene"),
                       c(chrom_col, start_col, end_col, gene_col))
  out <- as.data.frame(dt)
  if (!is.null(original_exon_number)) {
    out$exon_number <- original_exon_number
    out <- out[, names(bed_file), drop = FALSE]
  }
  stopifnot(nrow(out) == n_in)
  if (verbose) message(sprintf(
    "[INFO] pad_gene_terminal_exons: requested %d bp padding -- extended %d gene start(s) and %d gene end(s); %d side(s) received less than the full request.",
    padding, n_padded_start, n_padded_end, n_clamped))
  out
}

#' Gap-inserted x-axis positions for CNV window plots
#'
#' A gap of \code{gap} units separates consecutive genes. When \code{bed_file}
#' has a logical \code{off_target} column a run of consecutive off-target
#' intervals forms one block, separated from its neighbours.
#' @param bed_file BED-like data frame.
#' @param exon_range Integer row indices to place.
#' @param gap Gap size in x units.
#' @export
compute_gene_gap_positions <- function(bed_file, exon_range, gap = 1) {
  n <- length(exon_range)
  if (n == 0) {
    return(data.frame(idx = integer(0), px = numeric(0),
                      gene_break = logical(0), gene_group = integer(0)))
  }
  gap <- if (is.null(gap) || is.na(gap) || gap < 0) 1 else gap
  genes <- as.character(bed_file[[bed_cols(bed_file)$gene]][exon_range])
  genes[is.na(genes)] <- ""
  off <- if ("off_target" %in% names(bed_file)) as.logical(bed_file$off_target[exon_range]) else rep(FALSE, n)
  off[is.na(off)] <- FALSE
  both_off   <- off[-1] & off[-n]
  gene_break <- c(FALSE, (genes[-1] != genes[-n] & !both_off) | off[-1] != off[-n])
  step       <- ifelse(gene_break, 1 + gap, 1)
  step[1]    <- 0
  px         <- cumsum(step) + 1
  data.frame(idx = exon_range, px = px, gene_break = gene_break,
             gene_group = cumsum(gene_break) + 1L, stringsAsFactors = FALSE)
}

#' Flag background-exon calibration issues
#'
#' Under a well-calibrated model ~5\% of unaffected targets fall outside the
#' 95\% band; a binomial test flags windows where clearly more do.
#' @param ratio,lo,hi Observed ratio and its band, per target.
#' @param is_affected Logical, targets inside the call.
#' @param min_n Minimum background targets for a verdict.
#' @export
check_background_calibration <- function(ratio, lo, hi, is_affected, min_n = 5) {
  bg   <- !is_affected
  n_bg <- sum(bg, na.rm = TRUE)
  if (n_bg == 0)
    return(list(n_background = 0L, n_outside = 0L, pct_outside = NA_real_, flag = FALSE))
  outside <- (ratio[bg] < lo[bg]) | (ratio[bg] > hi[bg])
  outside[is.na(outside)] <- FALSE
  n_outside   <- sum(outside)
  pct_outside <- 100 * n_outside / n_bg
  flag <- n_bg >= min_n &&
    stats::pbinom(n_outside - 1L, size = n_bg, prob = 0.05, lower.tail = FALSE) < 0.05
  list(n_background = n_bg, n_outside = n_outside, pct_outside = pct_outside, flag = flag)
}

#' PCA of sample coverage profiles (static PDF)
#'
#' @param counts Data frame / matrix with one column per sample.
#' @param sample_names Sample columns to use.
#' @param output_pdf Output path (\code{NULL} prints the plot).
#' @param color_by Optional grouping vector, one value per sample.
#' @param scale Scale rows before PCA.
#' @export
plot_coverage_pca <- function(counts, sample_names, output_pdf = NULL,
                              color_by = NULL, scale = TRUE) {
  if (length(sample_names) < 3) {
    warning("PCA needs at least 3 samples (got ", length(sample_names), "); skipped.")
    return(invisible(NULL))
  }
  count_mat <- as.matrix(counts[, sample_names, drop = FALSE])
  log_mat   <- log2(count_mat + 1)
  row_var   <- apply(log_mat, 1, stats::var, na.rm = TRUE)
  keep      <- row_var > 0 & !is.na(row_var)
  if (sum(keep) < 2) {
    warning("Insufficient variation for PCA (fewer than 2 informative rows).")
    return(invisible(NULL))
  }
  log_mat <- log_mat[keep, , drop = FALSE]
  pca     <- stats::prcomp(t(log_mat), scale. = scale, center = TRUE)
  var_exp <- summary(pca)$importance[2, ] * 100
  pca_df  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Sample = sample_names)
  # The grouping column must exist BEFORE ggplot() is called (the plot keeps its
  # own copy of the data).
  has_group <- !is.null(color_by) && length(color_by) == nrow(pca_df)
  if (has_group) pca_df$Group <- as.factor(ifelse(is.na(color_by) | !nzchar(as.character(color_by)),
                                                  "unknown", as.character(color_by)))
  pal <- cnv_palette()
  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, label = Sample))
  if (has_group) {
    p <- p + ggplot2::aes(colour = Group)
    group_sizes <- table(pca_df$Group)
    if (sum(group_sizes >= 3) >= 1) {
      ell_df <- pca_df[pca_df$Group %in% names(group_sizes)[group_sizes >= 3], , drop = FALSE]
      p <- p + ggplot2::stat_ellipse(data = ell_df, ggplot2::aes(x = PC1, y = PC2, colour = Group),
                                     inherit.aes = FALSE, type = "norm", linetype = "dashed")
    }
  }
  p <- p +
    (if (has_group) ggplot2::geom_point(size = 3) else ggplot2::geom_point(size = 3, colour = pal$accent)) +
    ggplot2::labs(x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
                  y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
                  title = "PCA of sample coverage profiles") +
    cnv_theme()
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(size = 3, max.overlaps = 15)
  } else {
    p <- p + ggplot2::geom_text(check_overlap = TRUE, size = 3, hjust = -0.2, vjust = 0.5)
  }
  if (has_group) {
    p <- p + ggplot2::scale_colour_manual(values = cnv_group_colours(levels(pca_df$Group))) +
      ggplot2::guides(colour = ggplot2::guide_legend(title = "Group"))
  } else {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  if (!is.null(output_pdf)) {
    dir.create(dirname(output_pdf), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_pdf, p, width = 8, height = 6)
    message("[INFO] PCA plot saved: ", output_pdf)
  } else {
    print(p)
  }
  invisible(list(pca = pca, var_exp = var_exp))
}

# ---------------------------------------------------------------------------
# 3. Chromosome scope
# ---------------------------------------------------------------------------

#' Classify chromosome names as autosome / X / Y / MT / other
#'
#' Autosomes are \code{1}-\code{22} (with or without a \code{chr} prefix);
#' \code{23}/\code{24} are read as X/Y. Everything else (alt/decoy contigs) is
#' \code{"other"}; chrM is \code{"MT"} when \code{separate_mt} is TRUE and
#' \code{"other"} otherwise.
#' @param chr Character vector of chromosome names.
#' @param separate_mt Report chrM as its own class.
#' @return Character vector.
#' @export
chromosome_class <- function(chr, separate_mt = FALSE) {
  s   <- toupper(sub("^chr", "", trimws(as.character(chr)), ignore.case = TRUE))
  out <- rep("other", length(s))
  out[grepl("^([1-9]|1[0-9]|2[0-2])$", s)] <- "autosome"
  out[s %in% c("X", "23")] <- "X"
  out[s %in% c("Y", "24")] <- "Y"
  if (separate_mt) out[s %in% c("M", "MT")] <- "MT"
  out[is.na(chr)] <- NA_character_
  out
}

#' Which chromosome classes are in the input, and which a run analysed
#'
#' @param chromosomes Chromosome of every INPUT target (before any filtering).
#' @param analysed_classes Classes the run's settings analyse (subset of
#'   \code{"autosome"}, \code{"X"}, \code{"Y"}).
#' @param analysed_chromosomes Optional chromosomes of the targets that really
#'   reached the caller; a class that the settings analyse but with no
#'   surviving target is reported as removed by QC.
#' @param reasons Optional named character vector (names = class) with the
#'   reason a present class is not analysed.
#' @return Data frame: \code{class}, \code{label}, \code{targets},
#'   \code{present}, \code{analysed}, \code{reason}.
#' @export
chromosome_scope <- function(chromosomes, analysed_classes = "autosome",
                             analysed_chromosomes = NULL, reasons = NULL) {
  cls  <- chromosome_class(chromosomes)
  keys <- c("autosome", "X", "Y", "other")
  scope <- data.frame(
    class   = keys,
    label   = c("Autosomes (chr1-22)", "chrX", "chrY", "Other contigs"),
    targets = vapply(keys, function(k) sum(cls == k, na.rm = TRUE), integer(1), USE.NAMES = FALSE),
    stringsAsFactors = FALSE)
  scope$present <- scope$targets > 0L
  surviving <- if (is.null(analysed_chromosomes)) NULL else chromosome_class(analysed_chromosomes)
  scope$analysed <- scope$class %in% analysed_classes & scope$present &
    (is.null(surviving) | scope$class %in% surviving)
  scope$reason <- ifelse(scope$analysed, "",
    ifelse(!scope$present, "not in the BED / count table",
      ifelse(scope$class %in% analysed_classes, "all targets were removed by QC / coverage filters",
             "not analysed in this mode")))
  if (!is.null(reasons)) {
    hit <- !scope$analysed & scope$present & scope$class %in% names(reasons)
    scope$reason[hit] <- unname(reasons[scope$class[hit]])
  }
  scope
}

#' Warnings about the chromosome scope of a run
#'
#' Reported in the QC table, the log and the report alerts so a run that left
#' the autosomes or the gonosomes out is never silent.
#' @param scope Output of \code{chromosome_scope()}.
#' @param setting Name of the mode setting, for the hint text.
#' @param hints Named character vector (names = class) of how to analyse a
#'   class that this run skipped.
#' @return Data frame \code{Type}, \code{Details}, \code{Severity}.
#' @export
chromosome_scope_warnings <- function(scope, setting = "modechrom", hints = NULL) {
  rows <- list()
  add  <- function(type, details, severity)
    rows[[length(rows) + 1L]] <<- data.frame(Type = type, Details = details,
                                             Severity = severity, stringsAsFactors = FALSE)
  get  <- function(k) scope[scope$class == k, , drop = FALSE]
  auto <- get("autosome"); x <- get("X"); y <- get("Y"); other <- get("other")
  hint <- function(k) if (!is.null(hints) && k %in% names(hints)) hints[[k]] else paste("another", setting)

  if (!auto$present) add("Missing Chromosome", "Autosomes (chr1-22) not present in BED file", "error")
  if (!x$present)    add("Missing Chromosome", "chrX not present in BED file", "warning")
  if (!y$present)    add("Missing Chromosome", "chrY not present in BED file", "warning")

  skipped <- scope[scope$class %in% c("autosome", "X", "Y") & scope$present & !scope$analysed, , drop = FALSE]
  if ("autosome" %in% skipped$class) {
    add("Analysis scope",
        sprintf("Autosomes (chr1-22) NOT analysed in this run (%d targets; %s): no autosomal CNV is reported. Re-run with %s to analyse them.",
                auto$targets, auto$reason, hint("autosome")), "warning")
  }
  gon <- skipped[skipped$class %in% c("X", "Y"), , drop = FALSE]
  if (nrow(gon)) {
    add("Analysis scope",
        sprintf("Gonosomes NOT analysed in this run: %s. No CNV can be detected on them; re-run with %s.",
                paste(sprintf("%s (%d targets; %s)", gon$label, gon$targets, gon$reason), collapse = " and "),
                paste(unique(vapply(gon$class, hint, character(1))), collapse = " and ")),
        "warning")
  }
  if (other$present) {
    add("Analysis scope",
        sprintf("%d target(s) lie on non-canonical contigs (e.g. chrM, alt/decoy sequences); they are never analysed.",
                other$targets), "info")
  }
  if (!length(rows))
    return(data.frame(Type = character(), Details = character(), Severity = character(),
                      stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# 4. Logging
# ---------------------------------------------------------------------------

#' Create a pipeline logger
#'
#' One logger per run. Every line -- the pipeline's own and everything the
#' called functions emit with \code{message()} or \code{warning()} -- is written
#' to the log file as \code{[LEVEL] YYYY-mm-dd HH:MM:SS text} and echoed to the
#' console with the same format. The log file is what the HTML report's "Run"
#' tab reads, so it must contain the whole run (and only that run: the file is
#' truncated at start unless \code{append = TRUE}).
#'
#' @param log_file Log path (\code{NULL} logs to the console only).
#' @param tag Short name used in the banner (\code{"ECHO"} / \code{"CANOPE"}).
#' @param level Minimum level: \code{"DEBUG"}, \code{"INFO"}, \code{"WARNING"}
#'   or \code{"ERROR"}.
#' @param console Echo lines to the console.
#' @param append Keep the existing log content instead of truncating.
#' @param total_steps Number of \code{$step()} calls expected (for "Step k/N").
#' @return An object of class \code{cnv_logger}: a list of closures
#'   \code{log(msg, level)}, \code{info}, \code{warn}, \code{error},
#'   \code{debug}, \code{capture(expr)}, \code{step(name, expr, on_error)},
#'   \code{header(...)}, \code{footer()}, \code{counts()}, \code{path}.
#' @export
cnv_logger <- function(log_file = NULL, tag = "CNV", level = "INFO",
                       console = TRUE, append = FALSE, total_steps = NA_integer_) {
  lv_rank <- c(DEBUG = 10L, INFO = 20L, WARNING = 30L, ERROR = 40L)
  level   <- toupper(level)
  if (!level %in% names(lv_rank)) level <- "INFO"
  min_rank <- lv_rank[[level]]
  st <- new.env(parent = emptyenv())
  st$busy   <- FALSE
  st$step   <- 0L
  st$total  <- total_steps
  st$t0     <- Sys.time()
  st$n      <- c(DEBUG = 0L, INFO = 0L, WARNING = 0L, ERROR = 0L)
  if (!is.null(log_file)) {
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    if (!append || !file.exists(log_file)) cat("", file = log_file)
  }
  stamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  emit <- function(msg, type = "INFO") {
    type <- toupper(type)
    if (type == "WARN") type <- "WARNING"
    if (!type %in% names(lv_rank)) type <- "INFO"
    st$n[[type]] <- st$n[[type]] + 1L
    if (lv_rank[[type]] < min_rank) return(invisible(NULL))
    ln <- strsplit(paste(as.character(msg), collapse = " "), "\n", fixed = TRUE)[[1]]
    if (!length(ln)) ln <- ""
    lines <- sprintf("[%s] %s %s", type, stamp(), ln)
    if (!is.null(log_file)) cat(lines, file = log_file, append = TRUE, sep = "\n")
    if (console) {
      st$busy <- TRUE
      on.exit(st$busy <- FALSE, add = TRUE)
      for (l in lines) message(l)
    }
    invisible(NULL)
  }

  # Normalise a captured message: pick the level out of a leading "[LEVEL]" tag
  # and drop a timestamp that the emitter already embedded.
  parse_tagged <- function(txt, default = "INFO") {
    txt <- sub("\n+$", "", txt)
    m <- regmatches(txt, regexec("^\\[(DEBUG|INFO|WARNING|WARN|ERROR)\\]\\s*(.*)$", txt))[[1]]
    if (length(m) == 3L) { type <- m[2]; body <- m[3] } else { type <- default; body <- txt }
    body <- sub("^[0-9]{4}-[0-9]{2}-[0-9]{2}[ T][0-9]{2}:[0-9]{2}:[0-9]{2}(\\.[0-9]+)?\\s+", "", body)
    list(type = type, body = body)
  }

  capture <- function(expr) {
    old <- options(warn = 1L)
    on.exit(options(old), add = TRUE)
    withCallingHandlers(
      expr,
      message = function(m) {
        if (st$busy) return(invisible(NULL))
        txt <- conditionMessage(m)
        if (!nzchar(trimws(txt))) { invokeRestart("muffleMessage") }
        p <- parse_tagged(txt)
        if (nzchar(trimws(p$body))) emit(p$body, p$type)
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        txt <- conditionMessage(w)
        if (!grepl("sequence levels not in the other", txt, fixed = TRUE)) {
          p <- parse_tagged(txt, "WARNING")
          emit(p$body, if (p$type == "INFO") "WARNING" else p$type)
        }
        invokeRestart("muffleWarning")
      })
  }

  step <- function(name, expr, on_error = c("stop", "warn")) {
    on_error <- match.arg(on_error)
    st$step <- st$step + 1L
    k   <- st$step
    tot <- if (is.na(st$total)) "?" else as.character(st$total)
    emit(sprintf("Step %d/%s: %s ...", k, tot, name))
    t0 <- Sys.time()
    ok <- TRUE
    res <- tryCatch(expr, error = function(e) {
      ok <<- FALSE
      emit(sprintf("Step %d/%s failed (%s): %s", k, tot, name, conditionMessage(e)),
           if (on_error == "stop") "ERROR" else "WARNING")
      if (on_error == "stop") stop(e)
      NULL
    })
    if (ok) emit(sprintf("Step %d/%s completed in %.1fs", k, tot,
                         as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    invisible(res)
  }

  header <- function(config = NULL, packages = NULL) {
    emit(sprintf("%s pipeline started", tag))
    if (!is.null(config)) emit(paste("Configuration:", config))
    pk <- unique(c(tag, packages))
    ver <- vapply(pk, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_), "")
    ver <- ver[!is.na(ver)]
    if (length(ver)) emit(paste0("Packages: ", paste(sprintf("%s %s", names(ver), ver), collapse = ", ")))
    emit(paste0(R.version.string, " on ", R.version$platform))
    if (!is.null(log_file)) emit(paste("Log file:", basename(log_file)))
    invisible(NULL)
  }

  footer <- function(ok = TRUE) {
    secs <- as.numeric(difftime(Sys.time(), st$t0, units = "secs"))
    emit(sprintf("%s finished %s in %s (%d warning(s), %d error(s))", tag,
                 if (ok) "successfully" else "WITH ERRORS",
                 if (secs >= 120) sprintf("%.1f min", secs / 60) else sprintf("%.1f s", secs),
                 st$n[["WARNING"]], st$n[["ERROR"]]),
         if (ok) "INFO" else "ERROR")
    invisible(NULL)
  }

  structure(list(
    path = log_file, tag = tag, log = emit,
    debug = function(...) emit(paste0(..., collapse = ""), "DEBUG"),
    info  = function(...) emit(paste0(..., collapse = ""), "INFO"),
    warn  = function(...) emit(paste0(..., collapse = ""), "WARNING"),
    error = function(...) emit(paste0(..., collapse = ""), "ERROR"),
    capture = capture, step = step, header = header, footer = footer,
    set_total = function(n) st$total <- as.integer(n),
    counts = function() st$n),
    class = "cnv_logger")
}

# ---------------------------------------------------------------------------
# 5. Sample and target QC
# ---------------------------------------------------------------------------

#' Detect outlier samples (robust z-score of coverage noise)
#'
#' The noise of a sample is the MAD of \code{log2(count + pseudocount)} across
#' targets; its robust z-score is \code{(noise - median) / MAD} across the
#' cohort. With \code{sides = "upper"} (default) only samples that are
#' unusually \emph{noisy} are outliers -- a sample that is much cleaner than the
#' cohort is the best kind of reference and is never excluded. \code{"both"}
#' also flags unusually clean samples.
#'
#' @param counts Numeric matrix / data frame, targets x samples.
#' @param pseudocount Added before the log2.
#' @param z_threshold Robust z above which a sample is an outlier.
#' @param sides \code{"upper"} (default) or \code{"both"}.
#' @return Data frame: \code{sample}, \code{noise_score}, \code{robust_z},
#'   \code{is_outlier}.
#' @export
detect_outlier_samples <- function(counts, pseudocount = 0.5, z_threshold = 3,
                                   sides = c("upper", "both")) {
  sides  <- match.arg(sides)
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  log_counts   <- log2(pmax(counts, 0) + pseudocount)
  sample_noise <- apply(log_counts, 2, stats::mad, na.rm = TRUE)
  med  <- stats::median(sample_noise, na.rm = TRUE)
  mad0 <- stats::mad(sample_noise, na.rm = TRUE)
  z <- if (length(sample_noise) > 1L && is.finite(mad0) && mad0 > 0)
    (sample_noise - med) / mad0 else rep(0, length(sample_noise))
  is_outlier <- if (sides == "upper") z > z_threshold else abs(z) > z_threshold
  is_outlier[is.na(is_outlier)] <- FALSE
  data.frame(sample = colnames(counts), noise_score = unname(sample_noise),
             robust_z = unname(z), is_outlier = is_outlier, stringsAsFactors = FALSE)
}

#' Detect problematic targets
#'
#' A target is problematic when its cross-sample MAD (computed on
#' \code{log2(count + 1)}, so it does not grow with depth) is in the top
#' \code{1 - mad_quantile}, when its mean count is below \code{min_mean}, or
#' when its GC fraction is outside \code{[gc_min, gc_max]}. At least one target
#' per chromosome is always kept so no chromosome vanishes.
#'
#' @param count_matrix Numeric matrix, targets x samples.
#' @param chromosomes Optional chromosome label per target.
#' @param mad_quantile Quantile of the target MAD above which a target is flagged.
#' @param min_mean Minimum mean count.
#' @param gc Optional GC content per target (0-1 or 0-100).
#' @param gc_min,gc_max GC bounds (fractions).
#' @return Data frame: \code{exon} (row index), \code{mean}, \code{mad},
#'   \code{low_mean}, \code{high_mad}, \code{gc_extreme}, \code{problematic},
#'   \code{reason}.
#' @export
detect_problematic_exons <- function(count_matrix, chromosomes = NULL,
                                     mad_quantile = 0.90, min_mean = 20,
                                     gc = NULL, gc_min = 0.10, gc_max = 0.90) {
  count_matrix <- as.matrix(count_matrix)
  storage.mode(count_matrix) <- "double"
  exon_mean <- rowMeans(count_matrix, na.rm = TRUE)
  exon_mad  <- apply(log2(pmax(count_matrix, 0) + 1), 1, stats::mad, na.rm = TRUE)
  mad_thresh <- stats::quantile(exon_mad, probs = mad_quantile, na.rm = TRUE)
  low_mean   <- exon_mean < min_mean | !is.finite(exon_mean)
  high_mad   <- exon_mad > mad_thresh
  high_mad[is.na(high_mad)] <- FALSE
  gc_extreme <- gc_missing <- rep(FALSE, nrow(count_matrix))
  if (!is.null(gc)) {
    if (length(gc) != nrow(count_matrix))
      stop("'gc' must have one value per row of count_matrix (", nrow(count_matrix),
           " rows, got ", length(gc), ").")
    gc_val     <- if (suppressWarnings(max(gc, na.rm = TRUE)) > 1) gc / 100 else gc
    gc_missing <- !is.finite(gc_val)           # interval absent from the FASTA / GC failed
    gc_extreme <- !gc_missing & (gc_val < gc_min | gc_val > gc_max)
  }
  problematic <- low_mean | high_mad | gc_extreme | gc_missing
  if (!is.null(chromosomes)) {
    for (chr in unique(chromosomes)) {
      idx <- which(chromosomes == chr)
      if (length(idx) > 0 && all(problematic[idx])) {
        keep <- idx[which.min(exon_mad[idx])]
        problematic[keep] <- FALSE
        low_mean[keep] <- high_mad[keep] <- gc_extreme[keep] <- gc_missing[keep] <- FALSE
      }
    }
  }
  reason <- vapply(seq_along(problematic), function(i) {
    r <- c(if (low_mean[i]) sprintf("low mean count (%s < %s)", fmt_num(exon_mean[i], 1), format(min_mean)),
           if (high_mad[i]) sprintf("high variability (MAD %s > %s)", fmt_num(exon_mad[i], 2), fmt_num(mad_thresh, 2)),
           if (gc_extreme[i]) sprintf("extreme GC (outside %s-%s)", format(gc_min), format(gc_max)),
           if (gc_missing[i]) "GC content unavailable (interval not found in the reference sequence)")
    paste(r, collapse = "; ")
  }, character(1))
  data.frame(exon = seq_len(nrow(count_matrix)), mean = exon_mean, mad = exon_mad,
             low_mean = low_mean, high_mad = high_mad, gc_extreme = gc_extreme,
             problematic = problematic, reason = reason, stringsAsFactors = FALSE)
}

#' Per-sample coverage metrics (one definition for QC table and report)
#'
#' @param counts Raw count matrix, targets x samples.
#' @param z_threshold,sides Passed to \code{detect_outlier_samples()}.
#' @return Data frame: \code{sample}, \code{median_depth}, \code{total_reads},
#'   \code{max_corr} (best Spearman correlation with another sample, diagonal
#'   excluded by position), \code{noise_score}, \code{robust_z}, \code{is_outlier}.
#' @export
cnv_sample_metrics <- function(counts, z_threshold = 3, sides = "upper") {
  mat <- as.matrix(counts); storage.mode(mat) <- "double"
  nm  <- colnames(mat)
  out <- data.frame(sample = nm, median_depth = unname(col_medians(mat)),
                    total_reads = unname(colSums(mat, na.rm = TRUE)),
                    max_corr = NA_real_, noise_score = NA_real_, robust_z = NA_real_,
                    is_outlier = FALSE, stringsAsFactors = FALSE)
  if (length(nm) >= 2L) {
    cm <- suppressWarnings(stats::cor(mat, method = "spearman", use = "pairwise.complete.obs"))
    diag(cm) <- NA_real_
    out$max_corr <- unname(apply(cm, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)))
    nz <- detect_outlier_samples(mat, z_threshold = z_threshold, sides = sides)
    out$noise_score <- nz$noise_score; out$robust_z <- nz$robust_z; out$is_outlier <- nz$is_outlier
  }
  out
}

#' Decide which samples and targets are removed before calling
#'
#' The ONE place where the sample-outlier and problematic-target rules are
#' applied, so ECHO and CANOPE (and the QC-metrics step of each) always agree.
#' Order: (1) outlier samples among all candidate samples (test + reference),
#' (2) problematic targets, judged on the samples that remain.
#'
#' @param counts Numeric matrix, targets x candidate samples (column names).
#' @param chromosomes Chromosome of every target.
#' @param gc GC fraction per target (\code{NULL} skips the GC rule).
#' @param sample_qc,exon_qc Switch the two filters on/off.
#' @param qc_zscore,outlier_sides Sample-outlier threshold and sidedness.
#' @param exon_mad_quantile,min_exon_mean,gc_range Target rules.
#' @return A list: \code{sample_qc} (data frame from
#'   \code{detect_outlier_samples}, or \code{NULL}), \code{excluded_samples}
#'   (data frame \code{sample}, \code{reason}), \code{kept_samples},
#'   \code{target_qc} (data frame from \code{detect_problematic_exons}, or
#'   \code{NULL}), \code{keep_targets} (row indices to keep) and
#'   \code{sample_metrics} (see \code{cnv_sample_metrics()}).
#' @export
cnv_apply_qc <- function(counts, chromosomes = NULL, gc = NULL,
                         sample_qc = TRUE, exon_qc = TRUE,
                         qc_zscore = 3, outlier_sides = "upper",
                         exon_mad_quantile = 0.90, min_exon_mean = 20,
                         gc_range = c(0.15, 0.85)) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  kept  <- colnames(counts)
  excl  <- data.frame(sample = character(), reason = character(), stringsAsFactors = FALSE)
  sqc   <- if (ncol(counts) >= 2L) detect_outlier_samples(counts, z_threshold = qc_zscore, sides = outlier_sides) else NULL
  if (isTRUE(sample_qc) && !is.null(sqc) && any(sqc$is_outlier)) {
    o <- sqc[sqc$is_outlier, , drop = FALSE]
    excl <- data.frame(sample = o$sample,
                       reason = sprintf("noisy outlier (robust z = %.2f > %g)", o$robust_z, qc_zscore),
                       stringsAsFactors = FALSE)
    kept <- setdiff(kept, o$sample)
  }
  tqc  <- NULL
  keep <- seq_len(nrow(counts))
  if (isTRUE(exon_qc) && length(kept) >= 1L) {
    tqc <- detect_problematic_exons(
      counts[, kept, drop = FALSE], chromosomes = chromosomes,
      mad_quantile = exon_mad_quantile, min_mean = min_exon_mean, gc = gc,
      gc_min = if (!is.null(gc_range)) gc_range[1] else 0.10,
      gc_max = if (!is.null(gc_range)) gc_range[2] else 0.90)
    keep <- which(!tqc$problematic)
  }
  list(sample_qc = sqc, excluded_samples = excl, kept_samples = kept,
       target_qc = tqc, keep_targets = keep,
       sample_metrics = cnv_sample_metrics(counts, qc_zscore, outlier_sides))
}

#' Build the QC metrics table (one schema for both pipelines)
#'
#' Columns: \code{Sample}, \code{Exon}, \code{Type}, \code{Details},
#' \code{Gene}, \code{Severity} (\code{"error"}, \code{"warning"} or
#' \code{"info"}). \code{Type} is one of \code{"Cohort size"},
#' \code{"Whole sample"}, \code{"Whole exon"}, \code{"Exon variability"},
#' \code{"Target removed"}, \code{"Missing Chromosome"} or
#' \code{"Analysis scope"}.
#'
#' All metrics are computed on the RAW count matrix of every candidate sample
#' (before GC correction / sample / target filters), so a sample that was
#' excluded still shows the numbers that got it excluded.
#'
#' @param counts Numeric matrix (targets x samples, column names = samples).
#' @param targets BED-like data frame aligned with the rows of \code{counts}
#'   (chromosome/start/end/gene[/exon_number]).
#' @param min_corr,min_cov,min_total_reads,max_exon_cv QC thresholds.
#' @param min_samples Cohort size below which a "Cohort size" warning is added.
#' @param qc Result of \code{cnv_apply_qc()} (adds the excluded-sample /
#'   noisy-sample rows and the removed-target rows).
#' @param qc_zscore Threshold used by \code{qc} (for the Details text).
#' @param sample_qc \code{TRUE} when outliers were actually excluded; if
#'   \code{FALSE} they are reported as warnings only.
#' @param removed_targets Optional data frame (\code{label}, \code{gene},
#'   \code{reason}) of targets dropped before \code{counts} was built (e.g. by
#'   the GC filter or an unreadable reference interval).
#' @param scope_rows Output of \code{chromosome_scope_warnings()}.
#' @return Data frame with the columns above (zero rows when nothing is flagged).
#' @export
cnv_qc_metrics <- function(counts, targets,
                           min_corr = 0.98, min_cov = 100,
                           min_total_reads = 300000, max_exon_cv = 0.5,
                           min_samples = 3L, qc = NULL, qc_zscore = 3,
                           sample_qc = TRUE, removed_targets = NULL,
                           scope_rows = NULL) {
  mat <- as.matrix(counts)
  storage.mode(mat) <- "double"
  sample_names <- colnames(mat)
  labels <- target_labels(targets)
  bc     <- bed_cols(targets)
  gene   <- as.character(targets[[bc$gene]])
  gene   <- ifelse(is.na(gene) | gene %in% c("", ".", "NA"), "Unknown", gene)

  out <- list()
  add <- function(sample, exon, type, details, gene, severity) {
    n <- max(length(sample), length(exon), length(details))
    if (n == 0L) return(invisible(NULL))
    out[[length(out) + 1L]] <<- data.frame(
      Sample = rep_len(as.character(sample), n), Exon = rep_len(as.character(exon), n),
      Type = rep_len(type, n), Details = rep_len(as.character(details), n),
      Gene = rep_len(as.character(gene), n), Severity = rep_len(severity, n),
      stringsAsFactors = FALSE)
  }

  if (length(sample_names) < min_samples) {
    msg <- sprintf("Only %d sample(s): at least %d comparable samples are needed to build reference sets; expect few or unreliable calls.",
                   length(sample_names), min_samples)
    message("[WARNING] ", msg)
    add("All", "All", "Cohort size", msg, "All", "warning")
  }

  if (!is.null(qc) && nrow(qc$excluded_samples)) {
    ex <- qc$excluded_samples
    add(ex$sample, "All", "Whole sample", paste("Excluded from analysis:", ex$reason), "All", "error")
  } else if (!sample_qc && !is.null(qc$sample_qc) && any(qc$sample_qc$is_outlier)) {
    o <- qc$sample_qc[qc$sample_qc$is_outlier, , drop = FALSE]
    add(o$sample, "All", "Whole sample",
        sprintf("Noisy outlier (robust z = %.2f > %g); kept because sample_qc = FALSE", o$robust_z, qc_zscore),
        "All", "warning")
  }

  total_reads   <- colSums(mat, na.rm = TRUE)
  sample_median <- col_medians(mat)
  max_corr <- rep(NA_real_, length(sample_names))
  if (length(sample_names) >= 2L) {
    cm <- stats::cor(mat, method = "spearman", use = "pairwise.complete.obs")
    diag(cm) <- NA_real_   # by position: a duplicated sample (r = 1) must stay visible
    max_corr <- apply(cm, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
  }
  low_corr <- if (length(sample_names) >= 2L) which(is.na(max_corr) | max_corr < min_corr) else integer(0)
  if (length(low_corr))
    add(sample_names[low_corr], "All", "Whole sample",
        ifelse(is.na(max_corr[low_corr]), "Undefined correlation: no usable pairwise comparison",
               sprintf("Low correlation: %.3f (< %s)", max_corr[low_corr], format(min_corr))),
        "All", "warning")
  low_med <- which(sample_median < min_cov)
  if (length(low_med))
    add(sample_names[low_med], "All", "Whole sample",
        sprintf("Low median depth: %.1f (< %s)", sample_median[low_med], format(min_cov)), "All", "warning")
  low_reads <- which(total_reads < min_total_reads)
  if (length(low_reads))
    add(sample_names[low_reads], "All", "Whole sample",
        sprintf("Low total reads: %s (< %s)", format(round(total_reads[low_reads]), big.mark = ",", trim = TRUE),
                format(min_total_reads, big.mark = ",", scientific = FALSE, trim = TRUE)),
        "All", "warning")

  exon_median <- row_medians(mat)
  fail_exon   <- which(exon_median < min_cov)
  if (length(fail_exon))
    add("All", labels[fail_exon], "Whole exon",
        sprintf("Low median depth: %.1f (< %s)", exon_median[fail_exon], format(min_cov)),
        gene[fail_exon], "info")

  exon_mean <- rowMeans(mat, na.rm = TRUE)
  exon_sd   <- if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::rowSds(mat, na.rm = TRUE)
               else apply(mat, 1, stats::sd, na.rm = TRUE)
  exon_cv   <- exon_sd / pmax(exon_mean, 1)
  high_cv   <- which(exon_cv > max_exon_cv & !is.na(exon_cv))
  if (length(high_cv))
    add("All", labels[high_cv], "Exon variability",
        sprintf("High CV (> %s): %.3f", format(max_exon_cv), exon_cv[high_cv]), gene[high_cv], "info")

  if (!is.null(qc$target_qc) && any(qc$target_qc$problematic)) {
    rm <- which(qc$target_qc$problematic)
    add("All", labels[rm], "Target removed",
        paste("Excluded before calling:", qc$target_qc$reason[rm]), gene[rm], "info")
  }
  if (!is.null(removed_targets) && nrow(removed_targets))
    add("All", removed_targets$label, "Target removed",
        paste("Excluded before calling:", removed_targets$reason),
        removed_targets$gene %||% "Unknown", "info")

  if (!is.null(scope_rows) && nrow(scope_rows))
    add("All", "All", scope_rows$Type, scope_rows$Details, "All", scope_rows$Severity)

  if (!length(out))
    return(data.frame(Sample = character(), Exon = character(), Type = character(),
                      Details = character(), Gene = character(), Severity = character(),
                      stringsAsFactors = FALSE))
  do.call(rbind, out)
}

#' Write a QC metrics table to TSV
#' @param metrics Data frame from \code{cnv_qc_metrics()}.
#' @param output_file Output path.
#' @return \code{metrics}, invisibly.
#' @export
cnv_write_qc <- function(metrics, output_file) {
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(metrics, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(metrics)
}

# ---------------------------------------------------------------------------
# 6. Call post-processing
# ---------------------------------------------------------------------------

#' Flag CNV calls that only touch terminal exons of a gene
#'
#' A call covering every exon of the gene is a whole-gene event, not a boundary
#' artefact, and is never flagged.
#' @param cnv_calls Data frame with integer columns \code{global_start} and
#'   \code{global_end} (target indices into \code{bed_file}).
#' @param bed_file BED-like data frame (\code{exon_number} computed if absent).
#' @return \code{cnv_calls} with a logical \code{terminal_only} column.
#' @export
flag_terminal_only_calls <- function(cnv_calls, bed_file) {
  if (nrow(cnv_calls) == 0) return(cnv_calls)
  if (!"exon_number" %in% colnames(bed_file)) bed_file <- assign_exon_numbers_per_gene(bed_file)
  bc <- bed_cols(bed_file)
  cnv_calls$terminal_only <- FALSE
  gene_all <- as.character(bed_file[[bc$gene]])
  chr_all  <- as.character(bed_file[[bc$chrom]])
  for (i in seq_len(nrow(cnv_calls))) {
    s <- as.integer(cnv_calls$global_start[i]); e <- as.integer(cnv_calls$global_end[i])
    if (is.na(s) || is.na(e) || s < 1L || e > nrow(bed_file) || e < s) next
    affected <- bed_file$exon_number[s:e]
    genes <- unique(gene_all[s:e])
    if (length(genes) != 1L || is.na(genes) || !nzchar(genes)) next
    rows <- which(gene_all == genes & chr_all == chr_all[s])
    gene_exons <- bed_file$exon_number[rows]
    gene_exons <- gene_exons[!is.na(gene_exons)]
    if (!length(gene_exons) || all(is.na(affected))) next
    covers_whole_gene <- all(gene_exons %in% affected)
    if (!covers_whole_gene && all(affected %in% c(min(gene_exons), max(gene_exons))))
      cnv_calls$terminal_only[i] <- TRUE
  }
  cnv_calls
}

#' Cohort-level recurrence of CNV calls
#'
#' Calls that recur in many unrelated samples of the same run are far more
#' likely to be a systematic artefact (a noisy target, a capture problem) than
#' real germline CNVs. For each call this counts the samples that carry a
#' same-type call with at least \code{min_overlap} reciprocal overlap.
#'
#' @param chr,start,end,type,sample Vectors describing the calls.
#' @param n_samples Number of analysed samples (denominator of the fraction).
#' @param min_samples Minimum carriers for a call to be "recurrent".
#' @param min_fraction Minimum fraction of analysed samples carrying it.
#' @param min_overlap Minimum reciprocal overlap between two calls.
#' @return Data frame: \code{n_carriers}, \code{frac_carriers}, \code{recurrent}.
#' @export
cnv_recurrence <- function(chr, start, end, type, sample, n_samples,
                           min_samples = 3L, min_fraction = 0.25, min_overlap = 0.5) {
  n <- length(chr)
  res <- data.frame(n_carriers = rep(1L, n), frac_carriers = rep(NA_real_, n),
                    recurrent = rep(FALSE, n), stringsAsFactors = FALSE)
  if (n == 0L) return(res)
  chr <- as.character(chr); type <- toupper(as.character(type)); sample <- as.character(sample)
  start <- as.numeric(start); end <- as.numeric(end)
  grp <- paste(chr, type, sep = "\r")
  for (g in unique(grp)) {
    idx <- which(grp == g)
    if (length(idx) < 2L) next
    s <- start[idx]; e <- end[idx]; sm <- sample[idx]
    w <- pmax(e - s + 1, 1)
    for (a in seq_along(idx)) {
      ov <- pmin(e[a], e) - pmax(s[a], s) + 1
      hit <- ov > 0 & ov / w[a] >= min_overlap & ov / w >= min_overlap
      res$n_carriers[idx[a]] <- length(unique(sm[hit]))
    }
  }
  res$frac_carriers <- res$n_carriers / max(1, n_samples)
  res$recurrent <- res$n_carriers >= min_samples & res$frac_carriers >= min_fraction
  res
}

#' Apply the recurrence policy to a call table
#'
#' Adds the columns \code{Carriers} (number of samples carrying the same call)
#' and \code{Recurrent} (logical) -- always, whatever \code{action} is, so the
#' output schema is stable -- and, depending on \code{action}, downgrades
#' (\code{"downgrade"}: confidence -> LOW) or removes (\code{"remove"}) the
#' recurrent calls.
#' @param calls Call table.
#' @param cols Named list with the column names of the pipeline's call table:
#'   \code{sample}, \code{chr}, \code{start}, \code{end}, \code{type},
#'   \code{confidence}.
#' @param n_samples Number of analysed samples.
#' @param action \code{"flag"} (default), \code{"downgrade"}, \code{"remove"} or \code{"off"}.
#' @inheritParams cnv_recurrence
#' @return \code{calls} with columns \code{Carriers} and \code{Recurrent}.
#' @export
cnv_apply_recurrence <- function(calls, cols, n_samples, action = c("flag", "downgrade", "remove", "off"),
                                 min_samples = 3L, min_fraction = 0.25, min_overlap = 0.5) {
  action <- match.arg(action)
  if (is.null(calls)) return(calls)
  # Always present, so the output schema does not depend on the setting.
  calls$Carriers  <- rep(NA_integer_, nrow(calls))
  calls$Recurrent <- rep(FALSE, nrow(calls))
  if (!nrow(calls) || action == "off") return(calls)
  rec <- cnv_recurrence(calls[[cols$chr]], calls[[cols$start]], calls[[cols$end]],
                        calls[[cols$type]], calls[[cols$sample]], n_samples,
                        min_samples, min_fraction, min_overlap)
  calls$Carriers  <- rec$n_carriers
  calls$Recurrent <- rec$recurrent
  nrec <- sum(rec$recurrent)
  if (nrec) {
    message(sprintf("[WARNING] %d call(s) recur in >= %d sample(s) (>= %.0f%% of the cohort): likely systematic artefacts%s.",
                    nrec, min_samples, 100 * min_fraction,
                    switch(action, downgrade = "; confidence downgraded to LOW",
                           remove = "; removed (recurrence_action = 'remove')", "; flagged only")))
    if (action == "downgrade" && !is.null(cols$confidence) && cols$confidence %in% names(calls))
      calls[[cols$confidence]][rec$recurrent] <- "LOW"
    if (action == "remove") calls <- calls[!rec$recurrent, , drop = FALSE]
  }
  calls
}

# ---------------------------------------------------------------------------
# 7. Palette (PDF plots and HTML report share it)
# ---------------------------------------------------------------------------

#' The colour scheme shared by the static plots and the HTML report
#'
#' Colour-blind-safe (Okabe-Ito derived) and consistent: deletions are always
#' vermilion, duplications always blue, HIGH/MEDIUM/LOW confidence green /
#' amber / grey, severities red / amber / blue.
#' @return Named list of hex colours (light and dark variants for the report).
#' @export
cnv_palette <- function() {
  list(
    del = "#D55E00", dup = "#0072B2", normal = "#7A7A75",
    high = "#1B9E77", medium = "#E69F00", low = "#8C8C86", na = "#C3C2B7",
    fail = "#B3261E", warn = "#B36B00", ok = "#0B6B3A", info = "#2A78D6",
    accent = "#2A78D6", test = "#3B4A5A", reference = "#B4B4AE", band = "#E4E4DF",
    affected = "#D55E00", grid = "#E1E0D9", text = "#1B1B1A",
    dark = list(del = "#F08A4B", dup = "#56B4E9", high = "#3CC59A", medium = "#F2B134",
                low = "#8C8C86", fail = "#FF8A80", warn = "#FAB219", ok = "#5FD38D",
                info = "#6AA6F0", accent = "#6AA6F0", test = "#A9BCD0", reference = "#6F6E68", band = "#2C2C2A")
  )
}

#' ggplot2 theme shared by every static plot
#' @param base_size Base font size.
#' @export
cnv_theme <- function(base_size = 11) {
  pal <- cnv_palette()
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", colour = pal$text, hjust = 0),
      plot.subtitle    = ggplot2::element_text(colour = "#52514E", size = base_size - 1),
      axis.title       = ggplot2::element_text(colour = "#52514E"),
      axis.text        = ggplot2::element_text(colour = "#52514E"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = pal$grid, linewidth = 0.3),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(8, 12, 8, 8))
}

#' Colours for an arbitrary set of group labels (PCA groups, genes, ...)
#' @param groups Character vector of group labels.
#' @return Named character vector of colours.
#' @export
cnv_group_colours <- function(groups) {
  groups <- unique(as.character(groups))
  base <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#7A7A75")
  col <- rep_len(base, length(groups))
  col[tolower(groups) %in% c("unknown", "na", "")] <- "#A6A6A0"
  stats::setNames(col, groups)
}

# ---------------------------------------------------------------------------
# 8. VCF writer
# ---------------------------------------------------------------------------

#' Percent-encode the characters that would corrupt a VCF INFO value
#' @param x Character vector.
#' @noRd
vcf_info_escape <- function(x) {
  x <- enc2utf8(as.character(x))
  for (p in list(c("%", "%25"), c(";", "%3B"), c("=", "%3D"), c(",", "%2C"),
                 c("\t", "%09"), c("\n", "%0A"), c("\r", "%0D"), c(" ", "%20")))
    x <- gsub(p[1], p[2], x, fixed = TRUE)
  x
}

#' Sort key for VCF records: 1-22, X, Y, M, then everything else
#' @noRd
vcf_chrom_order <- function(chr) {
  stripped <- sub("^chr", "", as.character(chr), ignore.case = TRUE)
  num <- suppressWarnings(as.integer(stripped))
  special <- match(toupper(stripped), c("X", "Y", "M", "MT"))
  ifelse(!is.na(num), num, ifelse(!is.na(special), 22 + special, Inf))
}

#' Write a plain-text VCF 4.2 file
#'
#' Shared by both pipelines: each builds its meta lines, records and per-sample
#' fields; this writes them. Always LF line endings and UTF-8, a valid
#' header-only file when there are no records, and every field converted with
#' \code{as.character()} -- never through \code{apply()} on a data frame, which
#' would pad numeric columns such as POS with spaces.
#'
#' @param path Output file.
#' @param meta Character vector of \code{##} lines.
#' @param samples Sample names (VCF columns after FORMAT).
#' @param records Data frame with \code{CHROM}, \code{POS}, \code{ID},
#'   \code{ALT}, \code{INFO} (zero rows allowed).
#' @param format FORMAT string (e.g. \code{"GT:CN:FR"}).
#' @param sample_fields Character matrix, one row per record, columns named by
#'   \code{samples}.
#' @return \code{path}, invisibly.
#' @export
vcf_write <- function(path, meta, samples, records, format = "GT", sample_fields = NULL) {
  n <- nrow(records)
  cols <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  if (length(samples)) cols <- c(cols, "FORMAT", samples)
  lines <- c(meta, paste(cols, collapse = "\t"))
  if (n > 0L) {
    body <- list(as.character(records$CHROM), format(records$POS, scientific = FALSE, trim = TRUE),
                 as.character(records$ID), rep("N", n), as.character(records$ALT),
                 rep(".", n), rep("PASS", n), as.character(records$INFO))
    if (length(samples)) {
      body <- c(body, list(rep(format, n)),
                lapply(samples, function(s) as.character(sample_fields[, s])))
    }
    lines <- c(lines, do.call(paste, c(body, sep = "\t")))
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(enc2utf8(paste0(paste(lines, collapse = "\n"), "\n"))), con)
  message("[INFO] Wrote ", n, " CNV record(s) to ", path)
  invisible(path)
}

# ---------------------------------------------------------------------------
# Annotation helpers (both pipelines)
# ---------------------------------------------------------------------------

#' Extract gene symbols from structured BED target annotations
#'
#' Recognises RefSeq transcript prefixes, exon/UTR labels, score annotations,
#' and trailing genomic coordinates. Underscores and hyphens inside gene
#' symbols are preserved (for example, SCG5_GREM1 and BCR-ABL1).
#' Multiple comma-separated annotations are reduced to their unique symbols.
#' @param x Character vector of BED annotation names.
#' @param collapse Separator used when an interval maps to distinct symbols.
#' @return Character vector; NA where no structured annotation was recognised.
#' @noRd
extract_bed_gene_symbols <- function(x, collapse = ";") {
    if (is.null(x)) return(character(0))

    extract_one <- function(annotation) {
        annotation <- trimws(as.character(annotation))
        if (is.na(annotation) || !nzchar(annotation) || annotation == ".") {
            return(NA_character_)
        }

        has_coordinates <- grepl(
            "_(?:chr[^_]+|[0-9]+|X|Y|M|MT)_[0-9]+_[0-9]+$",
            annotation, perl = TRUE, ignore.case = TRUE
        )
        has_refseq <- grepl("^(?:NM|NR|XM|XR)_[0-9]+_", annotation,
                            perl = TRUE, ignore.case = TRUE)
        has_ensembl <- grepl("^ENST[0-9]+(?:\\.[0-9]+)?_", annotation,
                             perl = TRUE, ignore.case = TRUE)
        has_exon_or_utr <- grepl(
            "_(?:(?:5|3)?utr[0-9]*(?:[-_]ex[0-9]+)?|ex(?:on)?[_-]?[0-9]+)(?:$|[_(-])",
            annotation, perl = TRUE, ignore.case = TRUE
        )

        # Do not guess at custom formats. The caller can use its legacy,
        # explicitly configured field selection for those names.
        if (!(has_coordinates || has_refseq || has_ensembl || has_exon_or_utr)) {
            return(NA_character_)
        }

        # Remove the terminal coordinate triplet before symbol extraction.
        annotation <- sub(
            "_(?:chr[^_]+|[0-9]+|X|Y|M|MT)_[0-9]+_[0-9]+$",
            "", annotation, perl = TRUE, ignore.case = TRUE
        )
        # Scores such as (251/4972) are metadata, not part of the symbol.
        annotation <- gsub("\\s*\\([^()]*\\)", "", annotation, perl = TRUE)

        # Strip transcript accessions, not arbitrary underscore-delimited
        # prefixes. This is what protects fused symbols such as SCG5_GREM1.
        annotation <- sub("^(?:NM|NR|XM|XR)_[0-9]+_", "", annotation,
                          perl = TRUE, ignore.case = TRUE)
        annotation <- sub("^ENST[0-9]+(?:\\.[0-9]+)?_", "", annotation,
                          perl = TRUE, ignore.case = TRUE)

        # Remove structural exon/UTR suffixes while leaving the complete
        # preceding symbol (including internal underscores and hyphens).
        annotation <- sub(
            "_(?:(?:5|3)?utr[0-9]*(?:[-_]ex[0-9]+)?|ex(?:on)?[_-]?[0-9]+)(?:[_-].*)?$",
            "", annotation, perl = TRUE, ignore.case = TRUE
        )
        annotation <- trimws(annotation)
        if (!nzchar(annotation) || annotation == ".") return(NA_character_)
        annotation
    }

    vapply(as.character(x), function(value) {
        if (is.na(value) || !nzchar(trimws(value)) || trimws(value) == ".") {
            return(NA_character_)
        }
        annotations <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
        genes <- vapply(annotations, extract_one, character(1))
        genes <- genes[!is.na(genes) & nzchar(genes)]
        if (!length(genes)) return(NA_character_)
        paste(unique(genes), collapse = collapse)
    }, character(1), USE.NAMES = FALSE)
}

#' Filter a data frame by its chromosome column
#'
#' Accepts \code{CHR}, \code{chromosome} or \code{Chromosome} and ignores a
#' \code{chr} prefix mismatch (\code{chrX} matches \code{X}).
#' @param df Input data frame.
#' @param include Optional vector of chromosomes to include.
#' @param exclude Optional vector of chromosomes to exclude.
#' @return Filtered data frame (unchanged if empty).
#' @export
filter_chromosomes <- function(df, include = NULL, exclude = NULL) {
  if (is.null(df) || nrow(df) == 0) return(df)
  chrom_col <- intersect(c("chromosome", "Chromosome", "CHR"), colnames(df))[1]
  if (is.na(chrom_col)) stop("Data frame must contain a 'chromosome', 'Chromosome' or 'CHR' column.")
  norm <- function(x) unique(c(x, sub("^chr", "", x), paste0("chr", sub("^chr", "", x))))
  if (!is.null(include)) df <- df[df[[chrom_col]] %in% norm(include), , drop = FALSE]
  if (!is.null(exclude)) df <- df[!df[[chrom_col]] %in% norm(exclude), , drop = FALSE]
  df
}

# ---------------------------------------------------------------------------
# Per-chromosome depth and sex check (both pipelines)
# ---------------------------------------------------------------------------

#' Median target depth of every chromosome relative to the sample's autosomal median
#'
#' Computed on the counts of ALL input chromosomes (before the mode filter), so a
#' run restricted to the autosomes still reports chrX / chrY depth and the sex
#' check has something to work from.
#'
#' @param counts Numeric matrix / data frame, rows = targets, columns = samples.
#' @param chromosomes Chromosome of each row of \code{counts}.
#' @return Matrix (rows chromosomes in input order, columns samples), or \code{NULL}
#'   when no usable autosomal baseline exists.
#' @noRd
cnv_chrom_depth <- function(counts, chromosomes) {
  chromosomes <- as.character(chromosomes)
  if (is.null(dim(counts)) || nrow(counts) != length(chromosomes) || !ncol(counts)) return(NULL)
  auto <- chromosome_class(chromosomes) == "autosome"
  if (!any(auto, na.rm = TRUE)) return(NULL)
  lv <- unique(chromosomes)
  fac <- factor(chromosomes, levels = lv)
  cd <- vapply(colnames(counts), function(sn) {
    v <- as.numeric(if (is.data.frame(counts)) counts[[sn]] else counts[, sn]); base <- stats::median(v[auto], na.rm = TRUE)
    vapply(split(v, fac), function(x) if (is.finite(base) && base > 0) stats::median(x, na.rm = TRUE) / base else NA_real_, numeric(1))
  }, numeric(length(lv)))
  if (is.null(dim(cd))) cd <- matrix(cd, nrow = length(lv), dimnames = list(lv, colnames(counts)))
  cd
}

#' Normalise a sex / gender annotation to "F" / "M"
#' @param x Character vector (F, M, female, male, XX, XY ...).
#' @return \code{"F"}, \code{"M"} or \code{NA}.
#' @noRd
cnv_norm_sex <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ifelse(x %in% c("F", "FEMALE", "XX", "W"), "F",
         ifelse(x %in% c("M", "MALE", "XY"), "M", NA_character_))
}

#' Infer the sex of every sample from chrX / chrY depth
#'
#' A female has ~1x chrX and ~0x chrY, a male ~0.5x chrX and a clearly non-zero
#' chrY. Only meant to flag an annotation / mode mismatch, not to replace a
#' proper sex check.
#'
#' @param chrom_depth Matrix from \code{cnv_chrom_depth()}.
#' @param x_threshold chrX ratio below which a sample looks male.
#' @param y_threshold chrY ratio at/above which a sample looks male.
#' @return Data frame \code{sample, x_ratio, y_ratio, inferred_sex}
#'   (\code{"F"}, \code{"M"} or \code{"unclear"}), or \code{NULL} without chrX.
#' @noRd
cnv_infer_sex <- function(chrom_depth, x_threshold = 0.75, y_threshold = 0.1) {
  if (is.null(chrom_depth) || is.null(dim(chrom_depth))) return(NULL)
  cls <- chromosome_class(rownames(chrom_depth))
  if (!any(cls == "X", na.rm = TRUE)) return(NULL)
  xr <- as.numeric(chrom_depth[which(cls == "X")[1L], ])
  yr <- if (any(cls == "Y", na.rm = TRUE)) as.numeric(chrom_depth[which(cls == "Y")[1L], ]) else rep(NA_real_, length(xr))
  male   <- !is.na(xr) & xr < x_threshold  & (is.na(yr) | yr >= y_threshold)
  female <- !is.na(xr) & xr >= x_threshold & (is.na(yr) | yr <  y_threshold)
  data.frame(sample = colnames(chrom_depth), x_ratio = xr, y_ratio = yr,
             inferred_sex = ifelse(male, "M", ifelse(female, "F", "unclear")),
             stringsAsFactors = FALSE)
}

#' Sex check table: inferred sex vs the annotation and vs the sex the mode expects
#'
#' @param chrom_depth Matrix from \code{cnv_chrom_depth()}.
#' @param sample_table Optional data frame with \code{sample_name} and \code{gender}.
#' @param expected_sex \code{"F"}, \code{"M"} or \code{NA}: the sex the run mode is meant for.
#' @param analysed Samples that are actually called (mode mismatches only matter for them).
#' @param lg Logger from \code{cnv_logger()} (may be \code{NULL}).
#' @return \code{cnv_infer_sex()} table plus \code{annotated_sex}, \code{mismatch}
#'   and \code{mode_mismatch}; \code{NULL} without chrX.
#' @noRd
cnv_sex_check <- function(chrom_depth, sample_table = NULL, expected_sex = NA_character_,
                          analysed = character(), mode = NULL, lg = NULL) {
  sx <- cnv_infer_sex(chrom_depth)
  if (is.null(sx)) return(NULL)
  annotated <- if (is.data.frame(sample_table) && all(c("sample_name", "gender") %in% colnames(sample_table)))
    cnv_norm_sex(sample_table$gender[match(sx$sample, sample_table$sample_name)]) else rep(NA_character_, nrow(sx))
  sx$annotated_sex <- annotated
  sx$mismatch <- !is.na(annotated) & sx$inferred_sex %in% c("F", "M") & annotated != sx$inferred_sex
  sx$mode_mismatch <- if (is.na(expected_sex)) FALSE else {
    known <- ifelse(!is.na(annotated), annotated, ifelse(sx$inferred_sex %in% c("F", "M"), sx$inferred_sex, NA_character_))
    !is.na(known) & known != expected_sex
  }
  if (!is.null(lg)) {
    if (any(sx$mismatch))
      lg$warn(sprintf("Annotated sex disagrees with chrX/chrY coverage for: %s", paste(sx$sample[sx$mismatch], collapse = ", ")))
    bad <- sx$mode_mismatch & sx$sample %in% analysed
    if (any(bad))
      lg$warn(sprintf("modechrom = '%s' but these analysed samples look like the other sex: %s",
                      mode %||% "?", paste(sx$sample[bad], collapse = ", ")))
  }
  sx
}

# ---------------------------------------------------------------------------
# Helpers of the opt-in {targets} layers (echo_targets() / canope_targets())
# ---------------------------------------------------------------------------

#' Stop with a clear message when the optional {targets} package is missing
#' @noRd
cnv_require_targets <- function() {
  if (!requireNamespace("targets", quietly = TRUE))
    stop("The 'targets' package (>= 1.3.0) is required for the targets layer. ",
         "Install it with: install.packages('targets')", call. = FALSE)
  invisible(TRUE)
}

#' Files below \code{dir} written at or after \code{since}
#'
#' A \code{format = "file"} target must return the paths it produced; the
#' pipelines write many files (tables, VCF, reports, logs), so the layer reports
#' what the run touched instead of guessing names.
#' @noRd
cnv_files_since <- function(dir, since, slack = 2) {
  f <- list.files(dir, recursive = TRUE, full.names = TRUE)
  f <- f[!dir.exists(f)]
  f[as.numeric(file.mtime(f)) >= as.numeric(since) - slack]
}
