harmonise_chr_prefix <- function(ref_df, target_df) {
    has_chr_ref <- any(grepl("^chr", ref_df$chromosome))
    has_chr_tgt <- any(grepl("^chr", target_df$Chromosome))
    if (has_chr_ref && !has_chr_tgt) target_df$Chromosome <- paste0("chr", target_df$Chromosome)
    else if (!has_chr_ref && has_chr_tgt) target_df$Chromosome <- sub("^chr", "", target_df$Chromosome)
    target_df
}

prepare_plot_data <- function(call_row, counts, bed_file, exon_index, models, refs, gene_gap = 1) {
    sample     <- call_row$Sample
    idx_start  <- as.numeric(call_row$global_start)
    idx_end    <- as.numeric(call_row$global_end)
    target_chr <- normalize_chromosome_vec(call_row$Chromosome, bed_file$chromosome)
    ref_samples <- refs[[sample]]
    if (is.null(ref_samples) || !length(ref_samples)) return(NULL)
    if (is.null(models[[sample]]) || length(models[[sample]]) < 1 || is.na(models[[sample]][1])) return(NULL)

    exon_range_full <- seq(max(1, idx_start - 5), min(nrow(bed_file), idx_end + 5))
    bed_chr   <- normalize_chromosome_vec(bed_file$chromosome[exon_range_full], bed_file$chromosome)
    single_chr <- length(unique(bed_chr)) == 1
    prev <- FALSE; new_chr <- ""
    if (!single_chr) {
        prev    <- bed_chr[1] != target_chr
        new_chr <- if (prev) bed_file$chromosome[exon_range_full[1]] else bed_file$chromosome[utils::tail(exon_range_full, 1)]
        exon_range <- exon_range_full[bed_chr == target_chr]
    } else {
        exon_range <- exon_range_full
    }
    if (length(exon_range) == 0) return(NULL)

    # Gap-inserted x-axis positions, so every panel shows a visual break
    # between a gene's last exon and the next gene's first exon rather
    # than plotting them as if they were plain neighbouring exons. See
    # compute_gene_gap_positions() for the full rationale.
    gap_pos    <- compute_gene_gap_positions(bed_file, exon_range, gap = gene_gap)
    px         <- gap_pos$px
    gene_group <- gap_pos$gene_group

    test_median <- median(counts[exon_range, sample])
    ref_median  <- median(rowSums(counts[exon_range, ref_samples, drop = FALSE]))
    test_log    <- log(pmax(counts[exon_range, sample], 1))

    cov_list <- lapply(ref_samples, function(r) {
        scaling <- test_median / median(counts[exon_range, r])
        r_log   <- log(pmax(counts[exon_range, r] * scaling, 1))
        data.frame(exon_idx = exon_range, px = px, gene_group = gene_group,
                   coverage = r_log, group = r,
                   color_group = "Reference samples", stringsAsFactors = FALSE)
    })
    cov_list[[length(cov_list) + 1]] <- data.frame(
        exon_idx = exon_range, px = px, gene_group = gene_group, coverage = test_log,
        group = "Test sample", color_group = "Test sample", stringsAsFactors = FALSE)
    cov_data <- do.call(rbind, cov_list)

    pt_data <- data.frame(
        exon_idx    = exon_range,
        px          = px,
        gene_group  = gene_group,
        coverage    = test_log,
        color_group = ifelse(exon_range %in% (idx_start:idx_end), "Affected exon(s)", "Test sample"))
    pt_data <- pt_data[pt_data$color_group == "Affected exon(s)", ]

    ref_counts  <- rowSums(counts[exon_range, ref_samples, drop = FALSE])
    test_counts <- counts[exon_range, sample]
    totals      <- test_counts + ref_counts
    expected    <- ref_counts * (test_median / ref_median)
    expected_safe <- pmax(expected, 1)
    p_expected  <- expected_safe / totals
    p_expected  <- pmin(pmax(p_expected, 1e-6), 1 - 1e-6)
    ratio <- test_counts / expected_safe
    rho   <- models[[sample]][1]

    mins <- vapply(seq_along(exon_range), function(i) {
        qbetabinom(0.025, totals[i], max(rho, 0.005), p_expected[i]) / expected_safe[i]
    }, numeric(1))
    maxs <- vapply(seq_along(exon_range), function(i) {
        qbetabinom(0.975, totals[i], max(rho, 0.005), p_expected[i]) / expected_safe[i]
    }, numeric(1))

    ci_data <- data.frame(
        exon        = exon_range,
        px          = px,
        gene_group  = gene_group,
        ratio       = ratio,
        lo          = mins,
        hi          = maxs,
        is_affected = factor(exon_range %in% (idx_start:idx_end),
                             levels = c(FALSE, TRUE), labels = c("Observed", "Affected")))

    # New: z-score panel, matching CANOPE's design.
    #
    # ECHO's calling model is fundamentally different from CANOPE's (a
    # beta-binomial test/expected ratio, not a per-target NB(mean, var)
    # HMM), so this isn't a copy-paste of CANOPE's z-score formula — it's
    # the same *principle* (reuse the model's own already-validated
    # variance; don't invent a new, noisy small-sample one) applied to
    # ECHO's actual model. `rho` here is ExomeDepth's fitted overdispersion
    # parameter (genome-wide, not just this local window), the same
    # quantity the CI ribbon above is built from. The variance of a
    # Beta-Binomial(n, p, rho) is n*p*(1-p)*(1+(n-1)*rho) — standard result,
    # and consistent with the qbetabinom() parameterisation used above
    # (a = p(1-rho)/rho, b = (1-p)(1-rho)/rho  =>  a+b+1 = 1/rho).
    var_betabinom <- totals * p_expected * (1 - p_expected) *
        (1 + (totals - 1) * max(rho, 0.005))
    var_betabinom <- pmax(var_betabinom, 1)
    sd_z <- sqrt(var_betabinom)

    test_z <- (test_counts - expected_safe) / sd_z
    ref_z_list <- lapply(ref_samples, function(r) {
        scaling_r    <- test_median / median(counts[exon_range, r])
        ref_r_scaled <- counts[exon_range, r] * scaling_r
        (ref_r_scaled - expected_safe) / sd_z
    })
    names(ref_z_list) <- ref_samples

    z_data <- data.frame(
        exon        = exon_range,
        px          = px,
        gene_group  = gene_group,
        z           = test_z,
        is_affected = ci_data$is_affected)

    bg_calib <- check_background_calibration(ratio, mins, maxs, exon_range %in% (idx_start:idx_end))
    
    # Create the subtitle for CI plot
    ci_subtitle <- if (isTRUE(bg_calib$flag)) {
        sprintf("%d%% of background exons outside CI (%d/%d) \u2014 check region/reference match",
                round(bg_calib$pct_outside),
                bg_calib$n_outside, bg_calib$n_background)
    } else NULL

    list(cov_data = cov_data, pt_data = pt_data, ci_data = ci_data,
         z_data = z_data, ref_z_list = ref_z_list, bg_calibration = bg_calib,
         exon_range = exon_range, px = px, gene_group = gene_group,
         single_chr = single_chr, prev = prev,
         new_chr = new_chr, sample = sample, idx_start = idx_start, idx_end = idx_end,
         ci_subtitle = ci_subtitle)
}

save_cnv_pdf <- function(p_cov, p_genes, p_ci, p_zscore, file_path) {
    grDevices::pdf(file_path, useDingbats = FALSE, width = 8, height = 13)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(8, 1)))
    print(p_cov,    vp = grid::viewport(layout.pos.row = 1:3, layout.pos.col = 1))
    print(p_genes,  vp = grid::viewport(layout.pos.row = 4,   layout.pos.col = 1))
    print(p_ci,     vp = grid::viewport(layout.pos.row = 5:6, layout.pos.col = 1))
    print(p_zscore, vp = grid::viewport(layout.pos.row = 7:8, layout.pos.col = 1))
    grDevices::dev.off()
    invisible(NULL)
}

apply_xaxis_formatting <- function(p, single_chr, prev, exon_range, exon_index, px) {
    if (length(exon_range) == 0) return(p)
    min_p <- min(px); max_p <- max(px)
    if (single_chr) {
        b <- px; l <- exon_index[exon_range]; lim <- NULL
    } else if (prev) {
        b <- c((min_p - 6):(min_p - 1), px)
        l <- c(rep("", 6), exon_index[exon_range])
        lim <- c(min_p - 6.75, max_p)
    } else {
        b <- c(px, (max_p + 1):(max_p + 6))
        l <- c(exon_index[exon_range], rep("", 6))
        lim <- c(min_p, max_p + 6.75)
    }
    p + ggplot2::scale_x_continuous(breaks = b, labels = l, limits = lim) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, size = 6, hjust = 1))
}

create_coverage_plot <- function(cov_data, pt_data, single_chr, prev, exon_range, exon_index, px, sample_name) {
    cov_data$color_group <- ifelse(cov_data$color_group == "Test sample",
                                   paste0("Test sample (", sample_name, ")"),
                                   cov_data$color_group)
    cols <- c("Reference samples" = "gray",
              setNames("blue", paste0("Test sample (", sample_name, ")")),
              "Affected exon(s)" = "red")
    p_cov <- ggplot2::ggplot() +
        ggplot2::geom_point(data = subset(cov_data, color_group == "Reference samples"),
                            ggplot2::aes(x = px, y = coverage, color = color_group),
                            size = 2.5, alpha = 0.7) +
        ggplot2::geom_point(data = subset(cov_data, grepl("Test sample", color_group)),
                            ggplot2::aes(x = px, y = coverage, color = color_group),
                            size = 2.5) +
        ggplot2::geom_point(data = pt_data,
                            ggplot2::aes(x = px, y = coverage, color = color_group),
                            size = 3.5) +
        ggplot2::scale_colour_manual(
            values = cols,
            guide  = ggplot2::guide_legend(override.aes = list(size = 4), nrow = 1, title = NULL)) +
        ggplot2::labs(y = "Log (Coverage)", x = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank(),
                       legend.key = ggplot2::element_rect(fill = "white", colour = NA))
    apply_xaxis_formatting(p_cov, single_chr, prev, exon_range, exon_index, px)
}

create_gene_tile_plot <- function(bed_file, exon_range, single_chr, prev, new_chr, px) {
    if (length(exon_range) == 0) return(ggplot2::ggplot())
    temp       <- cbind(row = seq_len(nrow(bed_file)), bed_file)[exon_range, ]
    gene_names <- unique(bed_file$gene[exon_range])
    gene_names <- gene_names[!is.na(gene_names) & gene_names != ""]
    n_genes    <- length(gene_names)
    if (n_genes == 0) return(ggplot2::ggplot())
    if (n_genes == 1)      pal <- c("darkblue")
    else if (n_genes == 2) pal <- c("steelblue", "purple4")
    else {
        if (requireNamespace("RColorBrewer", quietly = TRUE)) {
            pal <- RColorBrewer::brewer.pal(min(n_genes, 9), "Purples")
            if (n_genes > 9) pal <- colorRampPalette(pal)(n_genes)
            pal <- rev(pal)
        } else if (requireNamespace("scales", quietly = TRUE)) {
            pal <- scales::hue_pal()(n_genes)
        } else pal <- rainbow(n_genes)
    }
    names(pal) <- gene_names
    # mid/width are computed on the gap-inserted px scale (not raw
    # exon_range indices) so each gene's tile lines up with that gene's
    # points/lines in the other panels, and adjacent tiles get real blank
    # space between them at a gene boundary instead of sitting flush.
    gene_tiles <- data.frame(
        gene  = gene_names,
        mid   = as.numeric(sapply(gene_names, function(g) mean(px[temp$gene == g]))),
        width = as.numeric(sapply(gene_names, function(g) sum(temp$gene == g))) - 0.5,
        y     = 1, stringsAsFactors = FALSE)
    if (!single_chr) {
        gene_tiles <- rbind(gene_tiles, data.frame(
            gene  = new_chr,
            mid   = ifelse(prev, min(px) - 5, max(px) + 5),
            width = 3.5, y = 1, stringsAsFactors = FALSE))
        pal <- c(pal, setNames("gray50", new_chr))
    }
    ggplot2::ggplot(gene_tiles, ggplot2::aes(x = mid, y = y, fill = gene, width = width, label = gene)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label = gene), size = 4.5, fontface = "bold", color = "white") +
        ggplot2::scale_fill_manual(values = pal) +
        ggplot2::theme_bw(base_family = "sans") +
        ggplot2::theme(legend.position = "none", panel.grid = ggplot2::element_blank(),
                       axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
        ggplot2::labs(x = "", y = "")
}

create_ci_plot <- function(ci_data, single_chr, prev, exon_range, exon_index, px, subtitle = NULL) {
    p <- ggplot2::ggplot(ci_data, ggplot2::aes(x = px)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, group = gene_group),
                             fill = "grey80", colour = NA) +
        ggplot2::geom_point(ggplot2::aes(y = ratio, color = is_affected), size = 3.5) +
        ggplot2::scale_color_manual(
            values = c("Observed" = "blue", "Affected" = "red"),
            guide  = ggplot2::guide_legend(override.aes = list(shape = 19, size = 3))) +
        ggplot2::labs(x = "", y = "Observed / Expected",
                      subtitle = subtitle) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none", legend.title = ggplot2::element_blank())
    apply_xaxis_formatting(p, single_chr, prev, exon_range, exon_index, px)
}

#' Z-Score Panel vs Reference Samples
#'
#' Ported from CANOPE. Shows the test sample's z-score (blue/red, by
#' affected status) against each individual reference sample's own z-score
#' (small gray points) at every exon in the window, all measured in units of
#' the beta-binomial model's own implied SD (see \code{prepare_plot_data()}
#' for the variance derivation) — so a well-behaved region should show the
#' reference points clustered near zero with the test sample clearly
#' separated only where a real CNV is present.
#' @noRd
create_zscore_plot <- function(z_data, ref_z_list, single_chr, prev, exon_range, exon_index) {
    px         <- z_data$px
    gene_group <- z_data$gene_group
    ref_list <- lapply(names(ref_z_list), function(r) {
        data.frame(px = px, gene_group = gene_group, z = ref_z_list[[r]],
                  sample = r, stringsAsFactors = FALSE)
    })
    ref_df <- do.call(rbind, ref_list)

    z_lim <- suppressWarnings(max(abs(c(ref_df$z, z_data$z)), na.rm = TRUE) * 1.1)
    if (!is.finite(z_lim) || z_lim <= 0) z_lim <- 1

    p <- ggplot2::ggplot() +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
        # Points only -- no line joins reference or test z-scores across exons.
        ggplot2::geom_point(data = ref_df,
                            ggplot2::aes(x = px, y = z), colour = "grey60", size = 1.5, alpha = 0.7) +
        ggplot2::geom_point(data = z_data, ggplot2::aes(x = px, y = z, color = is_affected),
                            size = 3) +
        ggplot2::scale_color_manual(values = c("Observed" = "blue", "Affected" = "red")) +
        ggplot2::coord_cartesian(ylim = c(-z_lim, z_lim)) +
        ggplot2::labs(y = "Z-score vs references", x = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")
    apply_xaxis_formatting(p, single_chr, prev, exon_range, exon_index, px)
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
#'   neighbouring exons. \code{0} disables the extra spacing. See
#'   \code{\link{compute_gene_gap_positions}}. Default \code{1}.
#' @return Invisibly returns the number of PDFs written.
#' @export
generate_plots <- function(rdata_file, output_dir = "./plots", modechrom = "A",
                           prefix = NULL, log_file = NULL, gene_gap = 1) {
    log_msg <- function(msg, type = "INFO") {
        if (!is.null(log_file)) cat(paste0("[", type, "] ", Sys.time(), " ", msg, "\n"),
                                    file = log_file, append = TRUE)
        message(msg)
    }
    log_msg("[INFO] BEGIN plot generation")
    objs      <- load_rdata(rdata_file, required = c("cnv_calls", "counts", "bed_file", "models", "refs"))
    cnv_calls <- objs$cnv_calls
    counts    <- objs$counts
    bed_file  <- objs$bed_file
    models    <- objs$models
    refs      <- objs$refs

    if (is.null(cnv_calls) || nrow(cnv_calls) == 0) {
        log_msg("[INFO] No CNV calls to plot.")
        return(invisible(0))
    }

    prefix_str <- if (is.null(prefix) || prefix == "") format(Sys.time(), "%Y%m%d-%H%M%S") else prefix
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    inc      <- switch(modechrom, X = "chrX", Y = "chrY", NULL)
    exc      <- if (modechrom == "A") c("chrX", "chrY") else NULL
    cnv_plot <- filter_chromosomes(cnv_calls, include = inc, exclude = exc)
    if (nrow(cnv_plot) == 0) {
        log_msg("[INFO] No CNV calls remain after chromosome filtering.")
        return(invisible(0))
    }

    cnv_plot <- harmonise_chr_prefix(bed_file, cnv_plot)

    exon_col <- grep("^(exon_number|custom\\.exon|exonnum)$",
                     colnames(bed_file), ignore.case = TRUE, value = TRUE)
    if (length(exon_col) == 0) {
        # fall back to the legacy "exon" column only if exon_number is absent
        exon_col <- grep("^exon$", colnames(bed_file), ignore.case = TRUE, value = TRUE)
    }
    if (length(exon_col) > 0) {
        exon_index <- bed_file[[exon_col[1]]]
    } else if (ncol(bed_file) >= 5) {
        exon_index <- bed_file[[5]]
    } else {
        exon_index <- compute_exon_index(bed_file)
    }
    exon_index <- ifelse(is.na(exon_index), "", as.character(exon_index))

    n_written <- 0L

    for (i in seq_len(nrow(cnv_plot))) {
        tryCatch({
            call_row <- cnv_plot[i, ]
            sample   <- call_row$Sample
            if (is.null(refs[[sample]]) || !length(refs[[sample]])) {
                log_msg(paste("[WARNING] Skipping plot for", sample, "(call", i, "): no reference samples"), "WARNING")
                next
            }
            if (is.null(models[[sample]]) || length(models[[sample]]) < 1 || is.na(models[[sample]][1])) {
                log_msg(paste("[WARNING] Skipping plot for", sample, "(call", i, "): missing model parameters"), "WARNING")
                next
            }
            plot_data <- prepare_plot_data(call_row, counts, bed_file, exon_index, models, refs, gene_gap = gene_gap)
            if (is.null(plot_data)) {
                log_msg(paste("[WARNING] Skipping plot for", sample, "(call", i, "): empty exon window"), "WARNING")
                next
            }
            gene_str   <- sanitize_filename(as.character(call_row$Gene))
            sample_dir <- file.path(output_dir, sanitize_filename(plot_data$sample))
            dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
            file_path <- file.path(sample_dir,
                                   paste0("ECHO_", prefix_str, "_", plot_data$sample, "_", gene_str, "_", i, ".pdf"))

            p_cov   <- create_coverage_plot(plot_data$cov_data, plot_data$pt_data,
                                            plot_data$single_chr, plot_data$prev,
                                            plot_data$exon_range, exon_index, plot_data$px,
                                            sample_name = plot_data$sample)
            p_genes <- create_gene_tile_plot(bed_file, plot_data$exon_range,
                                             plot_data$single_chr, plot_data$prev, plot_data$new_chr,
                                             plot_data$px)
            p_ci    <- create_ci_plot(plot_data$ci_data, plot_data$single_chr, plot_data$prev,
                                      plot_data$exon_range, exon_index, plot_data$px, 
                                      subtitle = plot_data$ci_subtitle)
            p_zscore <- create_zscore_plot(plot_data$z_data, plot_data$ref_z_list,
                                           plot_data$single_chr, plot_data$prev,
                                           plot_data$exon_range, exon_index)
            save_cnv_pdf(p_cov, p_genes, p_ci, p_zscore, file_path)
            n_written <- n_written + 1L
            log_msg(paste("[INFO] Saved plot:", file_path))
        }, error = function(e) {
            if (!is.null(grDevices::dev.list())) grDevices::dev.off()
            log_msg(paste("[ERROR] Failed to plot call", i, ":", e$message), "ERROR")
        })
    }
    log_msg(paste("[INFO] END plot generation (", n_written, " PDF(s) written)"))
    invisible(n_written)
}