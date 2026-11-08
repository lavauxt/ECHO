#' Normalise the target annotation shared by the PDF plots and the HTML report
#'
#' Guarantees the columns \code{chromosome}, \code{gene} (parsed symbol;
#' missing labels become \code{"Unknown"}), \code{exon_number} when the BED has
#' one, and a logical \code{off_target} flag for retained HorsROI targets. Also
#' covers older summary files where the fourth BED column was saved as a raw
#' annotation because preprocessing was skipped.
#'
#' @param bed_file Target annotation (one row per count-matrix row).
#' @param off_target_pattern The run's \code{off_target_pattern}.
#' @return The normalised data frame.
#' @noRd
normalize_target_annotation <- function(bed_file, off_target_pattern = "^HorsROI") {
    bed_file <- as.data.frame(bed_file, stringsAsFactors = FALSE)
    find_col <- function(candidates) {
        hit <- match(tolower(candidates), tolower(colnames(bed_file)))
        hit <- hit[!is.na(hit)]
        if (length(hit)) colnames(bed_file)[hit[1]] else NA_character_
    }
    chr_col <- find_col(c("chromosome", "chr", "chrom", "seqnames"))
    if (is.na(chr_col)) stop("[ERROR] bed_file has no chromosome column.")
    bed_file$chromosome <- as.character(bed_file[[chr_col]])

    gene_col <- find_col(c("gene", "gene_symbol", "symbol"))
    genes <- if (is.na(gene_col)) rep(NA_character_, nrow(bed_file)) else trimws(as.character(bed_file[[gene_col]]))
    parsed <- extract_bed_gene_symbols(genes, collapse = ";")
    ok <- !is.na(parsed) & nzchar(parsed)
    genes[ok] <- parsed[ok]
    genes[is.na(genes) | !nzchar(genes) | genes %in% c(".", "NA")] <- "Unknown"
    bed_file$gene <- genes

    exon_col <- find_col(c("exon_number", "custom.exon", "exonnum", "exon"))
    if (!is.na(exon_col)) bed_file$exon_number <- bed_file[[exon_col]]
    bed_file$off_target <- is_off_target_gene(genes, off_target_pattern)
    bed_file
}

harmonise_chr_prefix <- function(ref_df, target_df) {
    has_chr_ref <- any(grepl("^chr", ref_df$chromosome))
    has_chr_tgt <- any(grepl("^chr", target_df$Chromosome))
    if (has_chr_ref && !has_chr_tgt) target_df$Chromosome <- paste0("chr", target_df$Chromosome)
    else if (!has_chr_ref && has_chr_tgt) target_df$Chromosome <- sub("^chr", "", target_df$Chromosome)
    target_df
}

#' Per-target exon label used on the x axes of the PDF plots and the HTML report
#'
#' @param bed_file Target annotation (one row per count-matrix row).
#' @return Character vector, \code{""} where a target has no exon number.
#' @noRd
resolve_exon_index <- function(bed_file) {
    exon_col <- grep("^(exon_label|exon_number|custom\\.exon|exonnum)$",
                     colnames(bed_file), ignore.case = TRUE, value = TRUE)
    if (length(exon_col) == 0) {
        exon_col <- grep("^exon$", colnames(bed_file), ignore.case = TRUE, value = TRUE)
    }
    exon_index <- if (length(exon_col) > 0) {
        bed_file[[exon_col[1]]]
    } else if (ncol(bed_file) >= 5) {
        bed_file[[5]]
    } else {
        compute_exon_index(bed_file)
    }
    ifelse(is.na(exon_index), "", as.character(exon_index))
}

#' X-axis label of every target, shared by the PDF plots and the HTML report
#'
#' User-curated exon labels (BED column 5, kept as \code{counts$exon_label})
#' win over the derived numbering; off-target targets are labelled with their
#' gene label (\code{"HorsROI"}) instead of an empty string.
#'
#' @param bed_file Output of \code{normalize_target_annotation()}.
#' @param counts Optional count data frame aligned with \code{bed_file}.
#' @noRd
target_axis_labels <- function(bed_file, counts = NULL) {
    lab <- if (!is.null(counts) && "exon_label" %in% colnames(counts) &&
               length(counts$exon_label) == nrow(bed_file)) {
        as.character(counts$exon_label)
    } else {
        resolve_exon_index(bed_file)
    }
    lab[is.na(lab)] <- ""
    off <- bed_file$off_target %||% rep(FALSE, nrow(bed_file))
    lab[off] <- as.character(bed_file$gene[off])
    lab
}

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
         single_chr = single_chr, prev = prev, new_chr = new_chr,
         idx_start = idx_start, idx_end = idx_end)
}

# ---------------------------------------------------------------------------
# PDF panels
# ---------------------------------------------------------------------------

#' X range shared by every panel (the 6 empty slots show the next chromosome)
#' @noRd
pdf_x_limits <- function(pd) {
    mn <- min(pd$px); mx <- max(pd$px)
    if (pd$single_chr) c(mn, mx) else if (pd$prev) c(mn - 6.75, mx) else c(mn, mx + 6.75)
}

apply_xaxis_formatting <- function(p, pd) {
    px <- pd$px; mn <- min(px); mx <- max(px)
    if (pd$single_chr) {
        b <- px; l <- pd$xlab; lim <- NULL
    } else if (pd$prev) {
        b <- c((mn - 6):(mn - 1), px); l <- c(rep("", 6), pd$xlab); lim <- pdf_x_limits(pd)
    } else {
        b <- c(px, (mx + 1):(mx + 6)); l <- c(pd$xlab, rep("", 6)); lim <- pdf_x_limits(pd)
    }
    p + ggplot2::scale_x_continuous(breaks = b, labels = l, limits = lim) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, size = 6, hjust = 1))
}

#' Test-sample points (all, then affected on top) and reference points
#' @noRd
window_points <- function(pd, y, ref_y) {
    test <- data.frame(px = pd$px, y = y, affected = pd$affected)
    list(ref = data.frame(px = rep(pd$px, length(ref_y)), y = unlist(ref_y, use.names = FALSE)),
         test = test, hit = test[test$affected, , drop = FALSE])
}

# Subtle call-boundary shading shared by the PDF evidence panels.
affected_bg <- function(pd) {
    hit <- pd$px[!is.na(pd$affected) & pd$affected & is.finite(pd$px)]
    if (!length(hit)) return(NULL)
    data.frame(xmin = min(hit) - 0.4, xmax = max(hit) + 0.4,
               ymin = -Inf, ymax = Inf)
}

create_coverage_plot <- function(pd) {
    d <- window_points(pd, pd$test_log, pd$ref_log)
    test_lab <- paste0("Test sample (", pd$sample, ")")
    d$ref$color_group <- rep("Reference samples", nrow(d$ref))
    d$test$color_group <- rep(test_lab, nrow(d$test))
    d$hit$color_group <- rep("Affected exon(s)", nrow(d$hit))
    d$test$shape_group <- factor(d$test$affected, levels = c(FALSE, TRUE),
                                 labels = c("Observed", "Affected"))
    cols <- c("Reference samples" = "gray", stats::setNames("blue", test_lab), "Affected exon(s)" = "red")
    bg <- affected_bg(pd)
    p_cov <- ggplot2::ggplot(mapping = ggplot2::aes(x = px, y = y, color = color_group))
    if (!is.null(bg)) {
        p_cov <- p_cov + ggplot2::geom_rect(data = bg,
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#E34948", alpha = 0.07)
    }
    p_cov <- p_cov +
        ggplot2::geom_point(data = d$ref, size = 2.5, alpha = 0.7) +
        ggplot2::geom_point(data = d$test, ggplot2::aes(shape = shape_group), size = 2.5) +
        ggplot2::geom_point(data = d$hit, shape = 17, size = 3.5) +
        ggplot2::scale_shape_manual(values = c("Observed" = 16, "Affected" = 17), guide = "none") +
        ggplot2::scale_colour_manual(
            values = cols,
            guide  = ggplot2::guide_legend(override.aes = list(size = 4), nrow = 1, title = NULL)) +
        ggplot2::labs(y = "Log (Coverage)", x = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank(),
                       legend.key = ggplot2::element_rect(fill = "white", colour = NA))
    apply_xaxis_formatting(p_cov, pd)
}

#' Fill colours of the gene tiles
#' @noRd
gene_palette <- function(n) {
    if (n == 0) return(character())
    if (n == 1) return("darkblue")
    if (n == 2) return(c("steelblue", "purple4"))
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
        pal <- RColorBrewer::brewer.pal(min(n, 9), "Purples")
        if (n > 9) pal <- grDevices::colorRampPalette(pal)(n)
        return(rev(pal))
    }
    if (requireNamespace("scales", quietly = TRUE)) return(scales::hue_pal()(n))
    grDevices::rainbow(n)
}

#' Readable label colour on a tile fill
#' @noRd
label_colour <- function(fill) {
    rgb <- grDevices::col2rgb(fill) / 255
    ifelse(0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ] > 0.6, "grey15", "white")
}

#' Narrow blocks get a vertical label (same rule in the PDF and the HTML report)
#' @noRd
vertical_block_label <- function(n_targets, label) {
    n_targets == 1 | (n_targets == 2 & nchar(label) > 5)
}

#' Gene track: one tile per contiguous block (gene run or off-target interval)
#'
#' Tiles are built per block, not per gene name: a gene split by a HorsROI
#' interval, or several HorsROI intervals, used to be merged into one tile
#' centred on the average position, overlapping its neighbours.
#' @noRd
create_gene_tile_plot <- function(pd) {
    runs  <- split(seq_along(pd$px), pd$gene_group)
    tiles <- do.call(rbind, lapply(runs, function(ii) {
        a <- min(pd$px[ii]); b <- max(pd$px[ii])
        data.frame(gene = pd$gene[ii[1]], off = pd$off[ii[1]], mid = (a + b) / 2,
                   width = b - a + 0.5, n = length(ii), stringsAsFactors = FALSE)
    }))
    tiles <- tiles[!is.na(tiles$gene) & nzchar(tiles$gene), , drop = FALSE]
    if (!nrow(tiles)) return(ggplot2::ggplot())
    genes <- unique(tiles$gene[!tiles$off])
    pal   <- stats::setNames(gene_palette(length(genes)), genes)
    tiles$fill <- ifelse(tiles$off, "grey75", unname(pal[tiles$gene]))
    if (!pd$single_chr) {
        tiles <- rbind(tiles, data.frame(
            gene = pd$new_chr, off = FALSE, mid = if (pd$prev) min(pd$px) - 5 else max(pd$px) + 5,
            width = 3.5, n = 4L, fill = "gray50", stringsAsFactors = FALSE))
    }
    tiles$colour <- label_colour(tiles$fill)
    vert <- vertical_block_label(tiles$n, tiles$gene)
    tiles$angle <- ifelse(vert, 90, 0)
    tiles$size  <- ifelse(vert, 3, 4.5)
    ggplot2::ggplot(tiles, ggplot2::aes(x = mid, y = 1, width = width)) +
        ggplot2::geom_tile(ggplot2::aes(fill = fill)) +
        ggplot2::geom_text(ggplot2::aes(label = gene, angle = angle, size = size, colour = colour),
                           fontface = "bold") +
        ggplot2::scale_fill_identity() + ggplot2::scale_colour_identity() + ggplot2::scale_size_identity() +
        ggplot2::coord_cartesian(xlim = pdf_x_limits(pd)) +
        ggplot2::theme_bw(base_family = "sans") +
        ggplot2::theme(legend.position = "none", panel.grid = ggplot2::element_blank(),
                       axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank()) +
        ggplot2::labs(x = "", y = "")
}

#' 95% interval band per block; a one-target block gets a short flat band
#' (a ribbon through a single point is invisible).
#' @noRd
interval_bands <- function(pd, half = 0.3) {
    d <- data.frame(px = pd$px, lo = pd$lo, hi = pd$hi, gene_group = pd$gene_group)
    single <- stats::ave(d$px, d$gene_group, FUN = length) == 1
    if (!any(single)) return(d)
    s <- d[single, , drop = FALSE]
    d <- rbind(d[!single, , drop = FALSE], within(s, px <- px - half), within(s, px <- px + half))
    d[order(d$gene_group, d$px), , drop = FALSE]
}

create_ci_plot <- function(pd) {
    ci <- data.frame(px = pd$px, gene_group = pd$gene_group, ratio = pd$ratio, lo = pd$lo, hi = pd$hi,
                     is_affected = factor(pd$affected, levels = c(FALSE, TRUE), labels = c("Observed", "Affected")))
    bg <- affected_bg(pd)
    p <- ggplot2::ggplot(ci, ggplot2::aes(x = px)) +
        # Draw the confidence interval first. The affected-region shading must
        # be added afterwards; otherwise the ribbon is painted over it and the
        # red CNV overlay disappears in the CI panel.
        ggplot2::geom_ribbon(data = interval_bands(pd), ggplot2::aes(ymin = lo, ymax = hi, group = gene_group),
                             fill = "grey80", colour = NA)
    if (!is.null(bg)) {
        p <- p + ggplot2::geom_rect(data = bg,
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#E34948", alpha = 0.18)
    }
    p <- p +
        ggplot2::geom_point(ggplot2::aes(y = ratio, color = is_affected, shape = is_affected), size = 3.5) +
        ggplot2::scale_shape_manual(values = c("Observed" = 16, "Affected" = 17), guide = "none") +
        ggplot2::scale_color_manual(
            values = c("Observed" = "blue", "Affected" = "red"),
            guide  = ggplot2::guide_legend(override.aes = list(shape = c(16, 17), size = 3))) +
        ggplot2::labs(x = "", y = "Observed / Expected", subtitle = pd$ci_subtitle) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none", legend.title = ggplot2::element_blank())
    apply_xaxis_formatting(p, pd)
}

#' Z-Score Panel vs Reference Samples
#'
#' Ported from CANOPE. Shows the test sample's z-score (blue/red, by
#' affected status) against each individual reference sample's own z-score
#' (small gray points) at every exon in the window, all measured in units of
#' the beta-binomial model's own implied SD (see \code{prepare_plot_data()}),
#' so a well-behaved region shows the reference points clustered near zero
#' with the test sample clearly separated only where a real CNV is present.
#' @noRd
create_zscore_plot <- function(pd) {
    d <- window_points(pd, pd$z, pd$ref_z)
    d$test$is_affected <- factor(d$test$affected, levels = c(FALSE, TRUE), labels = c("Observed", "Affected"))
    z_lim <- suppressWarnings(max(abs(c(d$ref$y, d$test$y)), na.rm = TRUE) * 1.1)
    if (!is.finite(z_lim)) z_lim <- 0
    # Keep the +/-2 and +/-3 visual reference guides in view even for quiet calls.
    z_lim <- max(z_lim, 3.2)
    bg <- affected_bg(pd)
    p <- ggplot2::ggplot()
    if (!is.null(bg)) {
        p <- p + ggplot2::geom_rect(data = bg,
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#E34948", alpha = 0.07)
    }
    p <- p +
        ggplot2::geom_hline(yintercept = c(-3, 3), linetype = "dashed", colour = "#B3261E", alpha = 0.7, size = 0.4) +
        ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted", colour = "#B36B00", alpha = 0.8, size = 0.4) +
        ggplot2::geom_hline(yintercept = 0, linetype = "solid", colour = "black", size = 0.3) +
        # Points only -- no line joins reference or test z-scores across exons.
        ggplot2::geom_point(data = d$ref, ggplot2::aes(x = px, y = y), colour = "grey60", size = 1.5, alpha = 0.7) +
        ggplot2::geom_point(data = d$test, ggplot2::aes(x = px, y = y, color = is_affected, shape = is_affected), size = 3) +
        ggplot2::scale_shape_manual(values = c("Observed" = 16, "Affected" = 17), guide = "none") +
        ggplot2::scale_color_manual(values = c("Observed" = "blue", "Affected" = "red")) +
        ggplot2::coord_cartesian(ylim = c(-z_lim, z_lim)) +
        ggplot2::labs(y = "Z-score vs references", x = NULL) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none")
    apply_xaxis_formatting(p, pd)
}

#' Write the four panels of one call to a PDF, with aligned plot areas
#' @noRd
save_cnv_pdf <- function(plots, file_path) {
    grobs <- lapply(plots, ggplot2::ggplotGrob)
    # Same left/right margins in every panel so the gene tiles sit above their points.
    if (length(unique(vapply(grobs, function(g) length(g$widths), integer(1)))) == 1L) {
        w <- do.call(grid::unit.pmax, lapply(grobs, function(g) g$widths))
        grobs <- lapply(grobs, function(g) { g$widths <- w; g })
    }
    rows <- list(1:3, 4, 5:6, 7:8)
    grDevices::pdf(file_path, useDingbats = FALSE, width = 8, height = 13)
    on.exit(grDevices::dev.off(), add = TRUE)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(8, 1)))
    for (k in seq_along(grobs)) {
        grid::pushViewport(grid::viewport(layout.pos.row = rows[[k]], layout.pos.col = 1))
        grid::grid.draw(grobs[[k]])
        grid::upViewport()
    }
    invisible(NULL)
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
    log_msg <- function(msg, type = "INFO") {
        if (!is.null(log_file)) cat(paste0("[", type, "] ", Sys.time(), " ", msg, "\n"),
                                    file = log_file, append = TRUE)
        message("[", type, "] ", msg)
    }
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
