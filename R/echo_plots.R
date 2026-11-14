#' Coverage window of one CNV call (shared by the PDF plots and the HTML report)
#'
#' @param call_row One CNV call (needs \code{Sample}, \code{Chromosome},
#'   \code{global_start}, \code{global_end}).
#' @param bed_file Output of \code{normalize_target_annotation()}.
#' @param labels Output of \code{target_axis_labels()}.
#' @return A list of per-target vectors over the window, or \code{NULL} when
#'   the sample has no references / model or the window is empty.
#' @noRd
prepare_plot_data <- function(call_row, counts, bed_file, models, refs, gene_gap = 1,
                              labels = target_axis_labels(bed_file, counts)) {
    sample      <- as.character(call_row$Sample)
    ref_samples <- refs[[sample]]
    model       <- models[[sample]]
    if (!length(ref_samples) || !length(model) || is.na(model[1])) return(NULL)
    idx_start  <- as.numeric(call_row$global_start)
    idx_end    <- as.numeric(call_row$global_end)
    target_chr <- normalize_chromosome_vec(call_row$Chromosome, bed_file$chromosome)

    full    <- seq(max(1, idx_start - 5), min(nrow(bed_file), idx_end + 5))
    bed_chr <- normalize_chromosome_vec(bed_file$chromosome[full], bed_file$chromosome)
    single_chr <- length(unique(bed_chr)) == 1
    prev <- FALSE; new_chr <- ""; er <- full
    if (!single_chr) {
        prev    <- bed_chr[1] != target_chr
        new_chr <- bed_file$chromosome[if (prev) full[1] else utils::tail(full, 1)]
        er      <- full[bed_chr == target_chr]
    }
    if (length(er) == 0) return(NULL)

    gp       <- compute_gene_gap_positions(bed_file, er, gap = gene_gap)
    affected <- er %in% (idx_start:idx_end)
    cnt      <- function(s) as.numeric(counts[er, s])

    test_counts <- cnt(sample)
    test_median <- stats::median(test_counts)
    # max(., 1): a reference with zero median coverage in this window would
    # otherwise give an infinite scaling factor and NaN plots.
    ref_scaled <- lapply(ref_samples, function(r) { v <- cnt(r); v * test_median / max(stats::median(v), 1) })
    names(ref_scaled) <- ref_samples

    ref_counts <- rowSums(counts[er, ref_samples, drop = FALSE])
    totals     <- test_counts + ref_counts
    expected   <- pmax(ref_counts * (test_median / max(stats::median(ref_counts), 1)), 1)
    p_expected <- pmin(pmax(expected / totals, 1e-6), 1 - 1e-6)
    rho        <- max(model[1], 0.005)
    ratio      <- test_counts / expected
    lo <- qbetabinom(0.025, totals, rho, p_expected) / expected
    hi <- qbetabinom(0.975, totals, rho, p_expected) / expected

    # z-scores in units of the beta-binomial model's own SD
    sd_z  <- sqrt(pmax(totals * p_expected * (1 - p_expected) * (1 + (totals - 1) * rho), 1))
    ref_z <- lapply(ref_scaled, function(v) (v - expected) / sd_z)

    bg_calib <- check_background_calibration(ratio, lo, hi, affected)
    ci_subtitle <- if (isTRUE(bg_calib$flag)) {
        sprintf("%d%% of background exons outside CI (%d/%d) \u2014 check region/reference match",
                round(bg_calib$pct_outside), bg_calib$n_outside, bg_calib$n_background)
    } else NULL

    list(sample = sample, exon_range = er, px = gp$px, gene_group = gp$gene_group,
         gene = as.character(bed_file$gene[er]),
         off = as.logical(bed_file$off_target[er] %||% rep(FALSE, length(er))),
         xlab = labels[er], affected = affected,
         test_log = log(pmax(test_counts, 1)),
         ref_log = lapply(ref_scaled, function(v) log(pmax(v, 1))),
         ratio = ratio, lo = lo, hi = hi, observed = test_counts, expected = expected,
         z = (test_counts - expected) / sd_z, ref_z = ref_z,
         bg_calibration = bg_calib, ci_subtitle = ci_subtitle,
         call_type = toupper(substr(as.character(call_row$Type), 1, 3)),
         labels = list(cov = "log coverage", ratio = "Observed / expected",
                       z = "Z-score (beta-binomial SD)"),
         single_chr = single_chr, prev = prev, new_chr = new_chr,
         idx_start = idx_start, idx_end = idx_end)
}

#' Generate CNV Detection Plots
#'
#' Creates PDF plots summarising CNV calls.
#'
#' @param rdata_file Character string. Path to summary RData.
#' @param output_dir Character string. Output directory for PDF files.
#' @param modechrom Chromosome filter.
#' @param prefix Filename prefix.
#' @param log_file Optional path to log file.
#' @param gene_gap Numeric >= 0. Extra x-axis units inserted between a
#'   gene's last exon and the next gene's first exon, in every panel, so
#'   the two are visually separated instead of sitting flush like ordinary
#'   neighbouring exons. A run of consecutive off-target (HorsROI) intervals is
#'   one block, separated the same way. \code{0} disables the extra spacing. See
#'   \code{\link{compute_gene_gap_positions}}. Default \code{1}.
#' @param off_target_pattern Regex identifying off-target targets (the run's
#'   \code{off_target_pattern}); targets labelled \code{"HorsROI"} always are.
#' @return Invisibly returns the number of PDFs written.
#' @export
generate_plots <- function(rdata_file, output_dir = "./plots", modechrom = "A",
                           prefix = NULL, log_file = NULL, gene_gap = 1,
                           off_target_pattern = "^HorsROI") {
    lg <- cnv_logger(log_file, tag = "ECHO", append = TRUE)
    log_msg <- function(msg, type = "INFO") lg$log(msg, type)
    log_msg("BEGIN plot generation")
    objs      <- load_rdata(rdata_file, required = c("cnv_calls", "counts", "bed_file", "models", "refs"))
    cnv_calls <- objs$cnv_calls
    counts    <- objs$counts
    models    <- objs$models
    refs      <- objs$refs

    if (is.null(cnv_calls) || nrow(cnv_calls) == 0) {
        log_msg("No CNV calls to plot.")
        return(invisible(0))
    }
    bed_file <- normalize_target_annotation(objs$bed_file, off_target_pattern)
    labels   <- target_axis_labels(bed_file, counts)

    prefix_str <- if (is.null(prefix) || prefix == "") format(Sys.time(), "%Y%m%d-%H%M%S") else prefix
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    inc      <- switch(modechrom, X = "chrX", Y = "chrY", NULL)
    exc      <- if (modechrom == "A") c("chrX", "chrY") else NULL
    cnv_plot <- filter_chromosomes(cnv_calls, include = inc, exclude = exc)
    if (nrow(cnv_plot) == 0) {
        log_msg("No CNV calls remain after chromosome filtering.")
        return(invisible(0))
    }
    cnv_plot <- harmonise_chr_prefix(bed_file, cnv_plot)

    # Older call tables may carry a stale/raw Gene label. Rebuild it from the
    # normalized BED annotations whenever global target indices are saved.
    if (all(c("global_start", "global_end") %in% names(cnv_plot))) {
        cnv_plot$Gene <- genes_in_targets(bed_file$gene, cnv_plot$global_start, cnv_plot$global_end,
                                          drop_unknown = TRUE)
    } else if ("Gene" %in% names(cnv_plot)) {
        g <- as.character(cnv_plot$Gene)
        parsed <- extract_bed_gene_symbols(g, collapse = ";")
        ok <- !is.na(parsed) & nzchar(parsed)
        g[ok] <- parsed[ok]
        cnv_plot$Gene <- g
    }

    n_written <- 0L
    for (i in seq_len(nrow(cnv_plot))) {
        tryCatch({
            call_row <- cnv_plot[i, ]
            sample   <- call_row$Sample
            skip <- if (!length(refs[[sample]])) "no reference samples"
                    else if (!length(models[[sample]]) || is.na(models[[sample]][1])) "missing model parameters"
            if (is.null(skip)) {
                pd <- prepare_plot_data(call_row, counts, bed_file, models, refs, gene_gap = gene_gap, labels = labels)
                if (is.null(pd)) skip <- "empty exon window"
            }
            if (!is.null(skip)) {
                log_msg(paste0("Skipping plot for ", sample, " (call ", i, "): ", skip), "WARNING")
                next
            }
            # Calls spanning many genes would otherwise exceed the 255-char filename limit.
            gene_str   <- substr(sanitize_filename(as.character(call_row$Gene)), 1, 60)
            sample_dir <- file.path(output_dir, sanitize_filename(pd$sample))
            dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
            file_path <- file.path(sample_dir,
                                   paste0("ECHO_", prefix_str, "_", pd$sample, "_", gene_str, "_", i, ".pdf"))
            save_cnv_pdf(list(create_coverage_plot(pd), create_gene_tile_plot(pd),
                              create_ci_plot(pd), create_zscore_plot(pd)), file_path)
            n_written <- n_written + 1L
            log_msg(paste("Saved plot:", file_path))
        }, error = function(e) {
            log_msg(paste("Failed to plot call", i, ":", e$message), "ERROR")
        })
    }
    log_msg(paste0("END plot generation (", n_written, " PDF(s) written)"))
    invisible(n_written)
}
