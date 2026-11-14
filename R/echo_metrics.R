#' Execute QC Metrics Detection
#'
#' Computes the per-sample, per-target and chromosome-scope QC metrics of a run
#' and writes the flagged rows to a TSV file. It applies exactly the same
#' mode / sample-outlier / target filters as \code{\link{run_cnv_calling}}
#' (through \code{\link{cnv_apply_qc}}), so the QC table explains every sample
#' and target that calling will exclude. All metrics are computed on the RAW
#' counts of every candidate sample.
#'
#' Output columns: \code{Sample}, \code{Exon}, \code{Type}, \code{Details},
#' \code{Gene} and \code{Severity} (\code{"error"}, \code{"warning"} or
#' \code{"info"}). \code{Type} is one of \code{"Cohort size"},
#' \code{"Whole sample"}, \code{"Whole exon"}, \code{"Exon variability"},
#' \code{"Target removed"}, \code{"Missing Chromosome"} (class absent from the
#' BED) or \code{"Analysis scope"} (class present but not analysed by this run).
#'
#' @param rdata_file  Character string. Path to the coverage RData.
#' @param min_corr    Numeric. Minimum correlation. Default \code{0.98}.
#' @param min_cov     Numeric. Minimum median read depth. Default \code{100}.
#' @param min_total_reads Numeric. Minimum total reads per sample. Default \code{300000}.
#' @param max_exon_cv Numeric. Maximum allowed coefficient of variation per exon.
#' @param output_file Character string. Output path. Default \code{"./QC_metrics.tsv"}.
#' @param modechrom Character. \code{"A"}, \code{"X"} or \code{"Y"}; the
#'   chromosome mode the CNV calling will use. Default \code{"A"}.
#' @param sample_table Optional sample sheet (required for modes X / Y).
#' @param sample_qc,exon_qc,qc_zscore,outlier_sides,exon_mad_quantile,gc_extreme_filter,min_exon_mean
#'   The sample / target filters, as in \code{\link{run_cnv_calling}}.
#'
#' @return Invisibly returns the metrics data frame.
#' @export
run_qc_metrics <- function(rdata_file, min_corr = 0.98, min_cov = 100,
                           min_total_reads = 300000, max_exon_cv = 0.5,
                           output_file = "./QC_metrics.tsv", modechrom = "A",
                           sample_table = NULL,
                           sample_qc = TRUE, exon_qc = TRUE, qc_zscore = 3,
                           outlier_sides = "upper", exon_mad_quantile = 0.90,
                           gc_extreme_filter = c(0.15, 0.85), min_exon_mean = 20) {
  message("[INFO] BEGIN QC Metrics")
  min_corr        <- min_corr %||% 0.98
  min_cov         <- min_cov %||% 100
  min_total_reads <- min_total_reads %||% 300000
  max_exon_cv     <- max_exon_cv %||% 0.5

  objs <- load_rdata(rdata_file, required = c("counts", "sample_names", "bed_file"))
  md   <- echo_mode_data(objs, modechrom, sample_table)
  counts   <- md$counts
  bed_file <- assign_exon_numbers_per_gene(md$bed_file)
  gc <- if ("GC" %in% colnames(counts)) as.numeric(counts$GC) else NULL

  qc <- cnv_apply_qc(counts[, md$candidate, drop = FALSE], chromosomes = counts$chromosome, gc = gc,
                     sample_qc = sample_qc, exon_qc = exon_qc, qc_zscore = qc_zscore,
                     outlier_sides = outlier_sides, exon_mad_quantile = exon_mad_quantile,
                     min_exon_mean = min_exon_mean, gc_range = gc_extreme_filter)

  # Targets that never reached the counts (unreadable FASTA interval, ...)
  dropped <- objs$bed_dropped
  removed_targets <- NULL
  if (!is.null(dropped) && nrow(dropped)) {
    removed_targets <- data.frame(
      label = sprintf("%s:%s-%s", dropped$chromosome, dropped$start, dropped$end),
      gene = ifelse(is.na(dropped$gene) | dropped$gene == "", "Unknown", dropped$gene),
      reason = dropped$reason, stringsAsFactors = FALSE)
  }

  # Chromosome scope is judged on the INPUT chromosomes (including the intervals
  # dropped while extracting counts) so a whole chromosome lost there is visible.
  input_chr <- c(as.character(objs$bed_file$chromosome), as.character(dropped$chromosome))
  scope <- echo_scope(input_chr, md$modechrom, analysed_chromosomes = bed_file$chromosome[qc$keep_targets])
  scope_rows <- chromosome_scope_warnings(scope, "modechrom", ECHO_SCOPE_HINTS)
  for (i in seq_len(nrow(scope_rows)))
    if (scope_rows$Severity[i] != "info") message("[WARNING] ", scope_rows$Details[i])

  final_metrics <- cnv_qc_metrics(
    counts[, md$candidate, drop = FALSE], bed_file, min_corr = min_corr, min_cov = min_cov,
    min_total_reads = min_total_reads, max_exon_cv = max_exon_cv,
    qc = qc, qc_zscore = qc_zscore, sample_qc = sample_qc,
    removed_targets = removed_targets, scope_rows = scope_rows)
  cnv_write_qc(final_metrics, output_file)
  message("[INFO] QC metrics written to: ", output_file, " (", nrow(final_metrics), " flag(s))")
  message("[INFO] END QC Metrics")
  invisible(final_metrics)
}
