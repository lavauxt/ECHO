#' Execute QC Metrics Detection
#'
#' Checks coverage for low sample correlation, low median depth per sample,
#' low total read count, and high per‑exon coefficient of variation.
#' Writes flagged rows to a TSV file.
#'
#' @param rdata_file  Character string. Path to the coverage RData.
#' @param min_corr    Numeric. Minimum correlation. Default \code{0.98}.
#' @param min_cov     Numeric. Minimum median read depth. Default \code{100}.
#' @param min_total_reads Numeric. Minimum total reads per sample. Default \code{5e6}.
#' @param max_exon_cv Numeric. Maximum allowed coefficient of variation per exon.
#'   Exons with CV > this value are flagged. Default \code{0.5}.
#' @param output_file Character string. Output path. Default \code{"./QC_metrics.tsv"}.
#'
#' @return Invisibly returns the metrics data frame.
#' @export
run_qc_metrics <- function(rdata_file,
                           min_corr = 0.98,
                           min_cov = 100,
                           min_total_reads = 5e6,
                           max_exon_cv = 0.5,
                           output_file = "./QC_metrics.tsv") {
    message("[INFO] BEGIN QC Metrics")
    objs <- load_rdata(rdata_file, required = c("counts", "sample_names", "bed_file"))
    counts       <- objs$counts
    sample_names <- objs$sample_names
    bed_file     <- objs$bed_file

    message("[INFO] Computing QC Metrics")
    bed_file$exon_in_gene <- ave(seq_len(nrow(bed_file)), bed_file$gene, FUN = seq_along)

    dt_counts   <- data.table::as.data.table(counts[, sample_names, drop = FALSE])
    sample_median <- sapply(dt_counts, median, na.rm = TRUE)
    total_reads   <- colSums(dt_counts, na.rm = TRUE)

    if (ncol(dt_counts) < 2) {
        max_corr <- setNames(rep(NA_real_, ncol(dt_counts)), names(dt_counts))
        message("[WARNING] Correlation QC skipped: fewer than 2 samples")
    } else {
        corr_matrix <- stats::cor(dt_counts, use = "pairwise.complete.obs")
        max_corr <- apply(corr_matrix, 1, function(x) {
            others <- x[x != 1]
            if (!length(others)) return(NA_real_)
            max(others, na.rm = TRUE)
        })
    }

    add_metric <- function(m_list, samples, exons, type, details, genes) {
        m_list$Sample  <- c(m_list$Sample,  samples)
        m_list$Exon    <- c(m_list$Exon,    exons)
        m_list$Type    <- c(m_list$Type,    rep(type, length(samples)))
        m_list$Details <- c(m_list$Details, details)
        m_list$Gene    <- c(m_list$Gene,    genes)
        m_list
    }

    m <- list(Sample = character(), Exon = character(), Type = character(),
              Details = character(), Gene = character())

    low_corr <- which(!is.na(max_corr) & max_corr < min_corr)
    if (length(low_corr)) {
        details <- paste("Low correlation:", round(max_corr[low_corr], 2))
        m <- add_metric(m, sample_names[low_corr], rep("All", length(low_corr)),
                        "Whole sample", details, rep("All", length(low_corr)))
    }

    low_med <- which(sample_median < min_cov)
    if (length(low_med)) {
        details <- paste("Low median depth:", round(sample_median[low_med], 2))
        m <- add_metric(m, sample_names[low_med], rep("All", length(low_med)),
                        "Whole sample", details, rep("All", length(low_med)))
    }

    low_reads <- which(total_reads < min_total_reads)
    if (length(low_reads)) {
        details <- paste("Low total reads:", round(total_reads[low_reads], 0))
        m <- add_metric(m, sample_names[low_reads], rep("All", length(low_reads)),
                        "Whole sample", details, rep("All", length(low_reads)))
    }

    exon_median <- apply(dt_counts, 1, median, na.rm = TRUE)
    fail_exon <- which(exon_median < min_cov)
    if (length(fail_exon)) {
        exon_labels <- paste0(bed_file$gene[fail_exon], ":", bed_file$exon_in_gene[fail_exon])
        details <- paste("Low median depth:", round(exon_median[fail_exon], 2))
        m <- add_metric(m, rep("All", length(fail_exon)), exon_labels,
                        "Whole exon", details, bed_file$gene[fail_exon])
    }

    # --- FIX: Efficient CV calculation using matrixStats if available ---
    if (requireNamespace("matrixStats", quietly = TRUE)) {
        exon_sd <- matrixStats::rowSds(as.matrix(dt_counts), na.rm = TRUE)
        exon_mean <- rowMeans(dt_counts, na.rm = TRUE)
        exon_cv <- exon_sd / exon_mean
    } else {
        # Fallback to apply (slower but works)
        exon_cv <- apply(dt_counts, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
    }
    high_cv <- which(exon_cv > max_exon_cv & !is.na(exon_cv))
    if (length(high_cv)) {
        exon_labels <- paste0(bed_file$gene[high_cv], ":", bed_file$exon_in_gene[high_cv])
        details <- paste("High CV (>", max_exon_cv, "):", round(exon_cv[high_cv], 3))
        m <- add_metric(m, rep("All", length(high_cv)), exon_labels,
                        "Exon variability", details, bed_file$gene[high_cv])
    }

    final_metrics <- data.frame(m, stringsAsFactors = FALSE)

    if (nrow(final_metrics) == 0) {
        final_metrics <- data.frame(
            Sample = character(), Exon = character(), Type = character(),
            Details = character(), Gene = character(), stringsAsFactors = FALSE
        )
    }

    utils::write.table(final_metrics, file = output_file, quote = FALSE, row.names = FALSE, sep = "\t")
    message("[INFO] END QC Metrics")
    invisible(final_metrics)
}