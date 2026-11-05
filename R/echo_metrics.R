#' Execute QC Metrics Detection
#'
#' Checks coverage for low sample correlation, low median depth per sample,
#' low total read count, and high per-exon coefficient of variation, and
#' reports the chromosome scope of the run: a missing autosome/chrX/chrY class
#' in the BED, and the chromosome classes (autosomes / gonosomes) that the
#' selected \code{modechrom} does NOT analyse. Writes the flagged rows to a TSV
#' file.
#'
#' Output columns: \code{Sample}, \code{Exon}, \code{Type}, \code{Details},
#' \code{Gene} and \code{Severity} (\code{"error"}, \code{"warning"} or
#' \code{"info"}). \code{Type} is one of \code{"Whole sample"},
#' \code{"Whole exon"}, \code{"Exon variability"}, \code{"Cohort size"},
#' \code{"Missing Chromosome"} (class absent from the BED) or
#' \code{"Analysis scope"} (class present but not analysed by this run).
#'
#' @param rdata_file  Character string. Path to the coverage RData.
#' @param min_corr    Numeric. Minimum correlation. Default \code{0.98}.
#' @param min_cov     Numeric. Minimum median read depth. Default \code{100}.
#' @param min_total_reads Numeric. Minimum total reads per sample. Default \code{300000}.
#' @param max_exon_cv Numeric. Maximum allowed coefficient of variation per exon.
#' @param output_file Character string. Output path. Default \code{"./QC_metrics.tsv"}.
#' @param modechrom Character. The chromosome mode the CNV calling will use
#'   (\code{"A"}, \code{"X"} or \code{"Y"}); drives the autosome/gonosome scope
#'   warnings. Default \code{"A"}.
#'
#' @return Invisibly returns the metrics data frame.
#' @export
run_qc_metrics <- function(rdata_file, min_corr = 0.98, min_cov = 100,
                           min_total_reads = 300000, max_exon_cv = 0.5,
                           output_file = "./QC_metrics.tsv", modechrom = "A") {
    message("[INFO] BEGIN QC Metrics")
    min_corr        <- min_corr %||% 0.98
    min_cov         <- min_cov %||% 100
    min_total_reads <- min_total_reads %||% 300000
    max_exon_cv     <- max_exon_cv %||% 0.5
    modechrom       <- toupper(as.character(modechrom %||% "A")[1])

    objs         <- load_rdata(rdata_file, required = c("counts", "sample_names", "bed_file"))
    counts       <- objs$counts
    sample_names <- objs$sample_names
    bed_file     <- objs$bed_file

    message("[INFO] Computing QC Metrics")

    bed_file <- assign_exon_numbers_per_gene(bed_file)
    exon_in_gene <- bed_file$exon_number
    # "GENE:3" for annotated targets, "chr:start-end" for unannotated ones
    # (instead of the meaningless "NA:NA").
    chr_col <- intersect(colnames(bed_file), c("chromosome", "chr", "chrom", "seqnames"))[1]
    if (is.na(chr_col)) chr_col <- colnames(bed_file)[1]
    has_gene <- !is.na(bed_file$gene) & !bed_file$gene %in% c("", ".") & !is.na(exon_in_gene)
    exon_label_all <- ifelse(has_gene, paste0(bed_file$gene, ":", exon_in_gene),
                             paste0(bed_file[[chr_col]], ":", bed_file$start, "-", bed_file$end))
    gene_all <- ifelse(is.na(bed_file$gene) | bed_file$gene == "", "Unknown", bed_file$gene)

    mat           <- as.matrix(counts[, sample_names, drop = FALSE])
    sample_median <- if (requireNamespace("matrixStats", quietly = TRUE))
        matrixStats::colMedians(mat, na.rm = TRUE) else apply(mat, 2, median, na.rm = TRUE)
    names(sample_median) <- sample_names
    total_reads <- colSums(mat, na.rm = TRUE)

    if (ncol(mat) < 2) {
        max_corr <- setNames(rep(NA_real_, ncol(mat)), colnames(mat))
        message("[WARNING] Correlation QC skipped: fewer than 2 samples")
    } else {
        corr_matrix <- stats::cor(mat, method = "spearman", use = "pairwise.complete.obs")
        # Drop the diagonal by position. The previous `x[x != 1]` also dropped
        # every genuinely perfect (r = 1) pair, e.g. a duplicated sample, and
        # returned NA for it.
        diag(corr_matrix) <- NA_real_
        max_corr <- apply(corr_matrix, 1, function(x) {
            if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
        })
        names(max_corr) <- rownames(corr_matrix)
    }

    m <- list(Sample = character(), Exon = character(), Type = character(),
              Details = character(), Gene = character(), Severity = character())
    add_metric <- function(samples, exons, type, details, genes, severity = "warning") {
        n <- length(samples)
        m$Sample   <<- c(m$Sample, samples)
        m$Exon     <<- c(m$Exon, exons)
        m$Type     <<- c(m$Type, rep_len(type, n))
        m$Details  <<- c(m$Details, details)
        m$Gene     <<- c(m$Gene, genes)
        m$Severity <<- c(m$Severity, rep_len(severity, n))
    }

    if (length(sample_names) < 3) {
        msg <- sprintf(paste0("Only %d sample(s): ExomeDepth needs several comparable samples to build ",
                              "a reference set; expect few or unreliable calls."), length(sample_names))
        message("[WARNING] ", msg)
        add_metric("All", "All", "Cohort size", msg, "All")
    }

    # With >=2 samples, an undefined maximum correlation usually means this
    # sample has no usable pairwise comparison (e.g. constant/all-missing data).
    # A one-sample cohort is different: correlation QC is globally skipped.
    low_corr <- if (ncol(mat) >= 2L) {
        which(is.na(max_corr) | max_corr < min_corr)
    } else integer(0)
    if (length(low_corr))
        add_metric(sample_names[low_corr], rep("All", length(low_corr)),
                   "Whole sample", ifelse(is.na(max_corr[low_corr]),
                       "Undefined correlation: no usable pairwise comparison",
                       paste0("Low correlation: ", round(max_corr[low_corr], 2),
                              " (< ", min_corr, ")")),
                   rep("All", length(low_corr)))

    low_med <- which(sample_median < min_cov)
    if (length(low_med))
        add_metric(sample_names[low_med], rep("All", length(low_med)),
                   "Whole sample", paste0("Low median depth: ", round(sample_median[low_med], 2),
                                          " (< ", min_cov, ")"),
                   rep("All", length(low_med)))

    low_reads <- which(total_reads < min_total_reads)
    if (length(low_reads))
        add_metric(sample_names[low_reads], rep("All", length(low_reads)),
                   "Whole sample", paste0("Low total reads: ", round(total_reads[low_reads], 0),
                                          " (< ", format(min_total_reads, scientific = FALSE), ")"),
                   rep("All", length(low_reads)))

    exon_median <- if (requireNamespace("matrixStats", quietly = TRUE))
        matrixStats::rowMedians(mat, na.rm = TRUE) else apply(mat, 1, median, na.rm = TRUE)
    fail_exon   <- which(exon_median < min_cov)
    if (length(fail_exon)) {
        add_metric(rep("All", length(fail_exon)), exon_label_all[fail_exon],
                   "Whole exon", paste0("Low median depth: ", round(exon_median[fail_exon], 2)),
                   gene_all[fail_exon], severity = "info")
    }

    if (requireNamespace("matrixStats", quietly = TRUE)) {
        exon_sd   <- matrixStats::rowSds(mat, na.rm = TRUE)
        exon_mean <- rowMeans(mat, na.rm = TRUE)
        exon_cv   <- exon_sd / exon_mean
    } else {
        exon_cv <- apply(mat, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
    }
    high_cv <- which(exon_cv > max_exon_cv & !is.na(exon_cv))
    if (length(high_cv)) {
        add_metric(rep("All", length(high_cv)), exon_label_all[high_cv],
                   "Exon variability",
                   paste0("High CV (>", max_exon_cv, "): ", round(exon_cv[high_cv], 3)),
                   gene_all[high_cv], severity = "info")
    }

    # Chromosome scope: classes missing from the BED, and classes this run
    # (modechrom) does not analyse -- autosomes vs gonosomes.
    scope <- chromosome_scope(bed_file[[chr_col]], modechrom)
    scope_rows <- chromosome_scope_warnings(scope, modechrom)
    if (nrow(scope_rows)) {
        for (i in seq_len(nrow(scope_rows))) {
            if (scope_rows$Severity[i] != "info")
                message("[WARNING] ", scope_rows$Details[i])
        }
        add_metric(rep("All", nrow(scope_rows)), rep("All", nrow(scope_rows)),
                   scope_rows$Type, scope_rows$Details, rep("All", nrow(scope_rows)),
                   severity = scope_rows$Severity)
    }

    final_metrics <- data.frame(m, stringsAsFactors = FALSE)
    data.table::fwrite(final_metrics, file = output_file, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message("[INFO] END QC Metrics")
    invisible(final_metrics)
}
