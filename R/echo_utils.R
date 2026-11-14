
#' Quantile function of the Beta-Binomial distribution
#'
#' NOTE: uses the large-`size` approximation `qbeta(p, a, b) * size`.
#' For thin panels with small per-exon read counts the returned intervals
#' may deviate from the true beta-binomial quantiles. Use
#' \code{VGAM::qbetabinom.ab()} for exact quantiles if accuracy in that
#' regime matters.
#' @noRd
qbetabinom <- function(p, size, rho, prob) {
a <- prob * (1 - rho) / rho
b <- (1 - prob) * (1 - rho) / rho
qbeta(p, a, b) * size
}

#' @import data.table
NULL


#' Convert global exon indices to within-gene indices
#' @export
add_within_gene_indices <- function(cnv_calls, bed_file) {
if (nrow(cnv_calls) == 0) return(cnv_calls)
# Guard against double-conversion: once start.p has been rewritten to
# within-gene indices, it no longer equals global_start.
if ("global_start" %in% names(cnv_calls) &&
    any(cnv_calls$global_start != as.integer(cnv_calls$start.p), na.rm = TRUE)) {
  stop("[ERROR] add_within_gene_indices(): 'start.p' has already been ",
       "converted to within-gene indices; refusing to convert again.")
}
bed_sorted <- assign_exon_numbers_per_gene(bed_file)
exon_numbers <- bed_sorted$exon_number
cnv_calls$global_start <- as.integer(cnv_calls$start.p)
cnv_calls$global_end   <- as.integer(cnv_calls$end.p)
cnv_calls$start.p <- exon_numbers[cnv_calls$global_start]
cnv_calls$end.p   <- exon_numbers[cnv_calls$global_end]
cnv_calls
}


#' Load and validate YAML configuration
#' @export
load_config  <- function(yaml_path) {
if (!file.exists(yaml_path)) {
found_path  <- system.file(yaml_path, package = "ECHO")
if (found_path == "") {
stop("[ERROR] Config file not found at: ", yaml_path)
}
yaml_path  <- found_path
}
cfg  <- yaml::read_yaml(yaml_path)
required  <- c("input", "output")
if (!all(required %in% names(cfg))) {
stop("[ERROR] YAML must contain 'input' and 'output' sections")
}
return(cfg)
}


#' Validate a BED file
#'
#' @param bed_path Path to the BED file.
#' @param verbose Print messages.
#' @param bed_zero_based Logical. With the default \code{TRUE} (UCSC BED), an
#'   interval is invalid if \code{end <= start} or \code{start < 0}. With
#'   \code{FALSE} (1-based, closed intervals) a one-base interval
#'   \code{start == end} is legitimate and only \code{end < start} or
#'   \code{start < 1} is invalid.
#' @export
validate_bed_regions  <- function(bed_path, verbose = TRUE, bed_zero_based = TRUE) {
bed  <- tryCatch(data.table::fread(bed_path, header = FALSE, sep = "\t"),
error = function(e) NULL)
if (is.null(bed) || nrow(bed) == 0) {
issues <- list(unreadable_file = bed_path)
if (verbose) message("[WARNING] Could not read a non-empty BED file for validation: ", bed_path)
return(invisible(issues))
}
if (ncol(bed)  < 3) {
issues <- list(insufficient_columns = ncol(bed))
if (verbose) message("[WARNING] BED file has fewer than 3 columns.")
return(invisible(issues))
}
data.table::setnames(bed, 1:3, c("chrom", "start", "end"))
if (ncol(bed)  >= 4) data.table::setnames(bed, 4, "name") else bed$name  <- NA_character_
issues <- list()
if (verbose) message("[INFO] Validating BED file: ", bed_path, " (", nrow(bed), " rows)")
chrom_chr <- trimws(as.character(bed$chrom))
bad_chrom <- which(is.na(chrom_chr) | !nzchar(chrom_chr) | chrom_chr == ".")
if (length(bad_chrom)) {
    issues$invalid_chromosome <- as.data.frame(bed[bad_chrom, ])
    if (verbose) message("[WARNING] ", length(bad_chrom), " BED row(s) have a missing or blank chromosome.")
}
start_num <- suppressWarnings(as.numeric(as.character(bed$start)))
end_num <- suppressWarnings(as.numeric(as.character(bed$end)))
bad_order <- if (isFALSE(bed_zero_based)) {
    end_num < start_num | start_num < 1
} else {
    end_num <= start_num | start_num < 0
}
bad_coords <- which(!is.finite(start_num) | !is.finite(end_num) |
                    start_num != floor(start_num) | end_num != floor(end_num) |
                    bad_order)
if (length(bad_coords)) {
    issues$invalid_coordinates <- as.data.frame(bed[bad_coords, ])
    if (verbose) message("[WARNING] ", length(bad_coords),
                         if (isFALSE(bed_zero_based))
                           " row(s) have non-integer coordinates, end < start, or a start below 1."
                         else
                           " row(s) have non-integer coordinates, end <= start, or a negative start.")
}
dup_rows <- which(duplicated(bed[, c("chrom", "start", "end")]))
if (length(dup_rows)) {
    issues$duplicate_regions <- as.data.frame(bed[dup_rows, ])
    if (verbose) message("[WARNING] ", length(dup_rows), " duplicate chrom/start/end row(s) found.")
}
if (verbose && length(issues) == 0) message("[INFO] BED validation found no issues.")
invisible(issues)
}




#' Select the targets and samples a run analyses (ECHO)
#'
#' The single place where \code{modechrom} and the sample table are applied, so
#' \code{\link{run_cnv_calling}} and \code{\link{run_qc_metrics}} always work on
#' exactly the same targets and samples.
#'
#' @param objs List with \code{counts}, \code{sample_names}, \code{bed_file}
#'   and optionally \code{reference_sample_names} (see
#'   \code{\link{load_rdata}}).
#' @param modechrom \code{"A"} (autosomes), \code{"X"} or \code{"Y"}.
#' @param sample_table Path to a TSV with \code{sample_name} and \code{gender};
#'   required for \code{"X"} / \code{"Y"}.
#' @return List: \code{modechrom}, \code{counts}, \code{bed_file} (restricted to
#'   the mode), \code{sample_names} (samples to call), \code{reference_sample_names},
#'   \code{candidate} (all samples that may serve as test or reference) and
#'   \code{gender_map} (\code{NULL} in mode A).
#' @noRd
echo_mode_data <- function(objs, modechrom = "A", sample_table = NULL) {
  modechrom <- toupper(as.character(modechrom %||% "A")[1])
  if (!modechrom %in% c("A", "X", "Y"))
    stop("[ERROR] modechrom must be 'A' (autosomes), 'X' or 'Y' (got '", modechrom, "').")
  counts       <- objs$counts
  bed_file     <- objs$bed_file
  sample_names <- as.character(objs$sample_names)
  reference_sample_names <- as.character(objs$reference_sample_names %||% character(0))
  candidate <- unique(c(sample_names, reference_sample_names))
  missing_count_cols <- setdiff(candidate, colnames(counts))
  if (length(missing_count_cols))
    stop("[ERROR] Sample columns missing from coverage counts: ",
         paste(missing_count_cols, collapse = ", "))

  mode_class <- switch(modechrom, X = "X", Y = "Y", "autosome")
  message("[INFO] ", switch(modechrom, X = "Mode X: restricting to chromosome X",
                            Y = "Mode Y: restricting to chromosome Y",
                            "Mode A: autosomes only (excluding X and Y)"))
  keep_idx <- which(chromosome_class(bed_file$chromosome) == mode_class)
  if (length(keep_idx) == 0)
    stop("[ERROR] No ", switch(modechrom, X = "chrX", Y = "chrY", "autosomal"), " intervals found")
  bed_file <- bed_file[keep_idx, , drop = FALSE]
  counts   <- counts[keep_idx, , drop = FALSE]
  rownames(bed_file) <- rownames(counts) <- NULL

  gender_map <- NULL
  if (modechrom %in% c("X", "Y")) {
    if (is.null(sample_table) || !file.exists(sample_table))
      stop("[ERROR] sample_table required for mode ", modechrom)
    gender_df <- utils::read.table(sample_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (!all(c("sample_name", "gender") %in% colnames(gender_df)))
      stop("[ERROR] sample_table must have 'sample_name' and 'gender'")
    gender_df$gender <- toupper(substr(gender_df$gender, 1, 1))
    gender_map <- stats::setNames(gender_df$gender, gender_df$sample_name)
    missing <- setdiff(candidate, names(gender_map))
    if (length(missing) > 0)
      stop("[ERROR] Missing gender for analysis/reference samples: ", paste(missing, collapse = ", "))
  }
  if (modechrom == "Y") {
    male_samples <- names(gender_map[gender_map == "M"])
    sample_names <- intersect(sample_names, male_samples)
    reference_sample_names <- intersect(reference_sample_names, male_samples)
    candidate <- unique(c(sample_names, reference_sample_names))
    if (length(sample_names) == 0) stop("[ERROR] No male analysis samples for mode Y")
    message("[INFO] Mode Y: keeping ", length(sample_names), " male analysis sample(s) and ",
            length(reference_sample_names), " male external reference sample(s)")
  }
  list(modechrom = modechrom, counts = counts, bed_file = bed_file, sample_names = sample_names,
       reference_sample_names = reference_sample_names, candidate = candidate, gender_map = gender_map)
}

#' Chromosome scope of an ECHO run
#'
#' @param chromosomes Chromosome of every target of the BED \emph{before} the
#'   mode filter (and before any interval was dropped).
#' @param modechrom \code{"A"}, \code{"X"} or \code{"Y"}.
#' @param analysed_chromosomes Chromosomes of the targets that reached the caller.
#' @return See \code{\link{chromosome_scope}}.
#' @noRd
echo_scope <- function(chromosomes, modechrom = "A", analysed_chromosomes = NULL) {
  modechrom <- toupper(as.character(modechrom %||% "A")[1])
  cls <- switch(modechrom, X = "X", Y = "Y", "autosome")
  chromosome_scope(chromosomes, analysed_classes = cls, analysed_chromosomes = analysed_chromosomes,
                   reasons = c(autosome = sprintf("not analysed with modechrom = '%s'", modechrom),
                               X = sprintf("not analysed with modechrom = '%s'", modechrom),
                               Y = sprintf("not analysed with modechrom = '%s'", modechrom)))
}

# Hints shown when a chromosome class was skipped by the selected modechrom.
ECHO_SCOPE_HINTS <- c(autosome = "modechrom = 'A'", X = "modechrom = 'X'", Y = "modechrom = 'Y'")


# ---------------------------------------------------------------------------
# Coverage reuse (echo(reuse_coverage = TRUE)) and the targets layer
# ---------------------------------------------------------------------------

#' BAM files an ECHO run reads: the test BAMs and the external reference BAMs
#'
#' Resolves \code{bamfiles} (vector, single BAM or text list), \code{bamdir} and
#' \code{rbams} (directory or TSV with a \code{bam} column) the way
#' \code{run_bam_coverage()} does.
#' @noRd
echo_input_bams <- function(bamfiles = NULL, bamdir = NULL, rbams = NULL) {
  out <- character()
  if (length(bamfiles)) {
    bf <- as.character(bamfiles)
    out <- c(out, if (length(bf) > 1L || grepl("\\.bam$", bf, ignore.case = TRUE)) bf
                  else if (file.exists(bf)) trimws(readLines(bf, warn = FALSE)))
  } else if (!is.null(bamdir) && dir.exists(bamdir)) {
    out <- c(out, list.files(bamdir, pattern = "\\.bam$", full.names = TRUE))
  }
  if (length(rbams) == 1L && !is.na(rbams)) {
    if (dir.exists(rbams)) {
      out <- c(out, list.files(rbams, pattern = "\\.bam$", full.names = TRUE))
    } else if (file.exists(rbams)) {
      tb <- tryCatch(utils::read.delim(rbams, header = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(tb) && "bam" %in% colnames(tb)) {
        rb <- trimws(as.character(tb$bam)); rb <- rb[!is.na(rb) & nzchar(rb)]
        rel <- !grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", rb)
        rb[rel] <- file.path(dirname(normalizePath(rbams)), rb[rel])
        out <- c(out, rb)
      }
    }
  }
  unique(out[nzchar(out)])
}

#' Fingerprint of everything that decides the coverage RData
#'
#' MD5 of the size and modification time of every input file plus the counting
#' settings: equal fingerprints mean the BAM counting would give the same table.
#' @noRd
echo_coverage_fingerprint <- function(files, settings = list()) {
  files <- sort(unique(as.character(unlist(files))))
  files <- files[!is.na(files) & file.exists(files)]
  fi <- file.info(files)
  txt <- c(sprintf("%s|%s|%s", normalizePath(files, mustWork = FALSE), fi$size, as.numeric(fi$mtime)),
           vapply(names(settings), function(n) paste0(n, "=", paste(deparse(settings[[n]]), collapse = "")), ""),
           paste0("ECHO=", tryCatch(as.character(utils::packageVersion("ECHO")), error = function(e) "dev")))
  tmp <- tempfile(); on.exit(unlink(tmp), add = TRUE)
  writeLines(txt, tmp)
  unname(tools::md5sum(tmp))
}
