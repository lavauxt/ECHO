# ---------------------------------------------------------------------------
# cnv_shared_plots.R -- static (PDF) CNV window panels, IDENTICAL in ECHO and
# CANOPE (kept byte-identical by tests/testthat/test-shared-sync.R).
#
# Every panel is drawn from one "window" object \code{pd} that each pipeline
# builds from its own model (ECHO: beta-binomial, CANOPE: negative binomial):
#
#   sample, call_type ("DEL"/"DUP"), px, gene_group, gene, off, xlab, affected,
#   test_log, ref_log (list), ratio, lo, hi, observed, expected, z, ref_z (list),
#   ci_subtitle, single_chr, prev, new_chr, labels (list: cov, ratio, z)
#
# Colours come from cnv_palette(): the SAME values the HTML report uses.
# ---------------------------------------------------------------------------

#' Colour of the affected exons of a window: vermilion for a deletion, blue for a duplication
#' @noRd
affected_colour <- function(pd) {
  pal <- cnv_palette()
  tp <- tolower(as.character(pd$call_type %||% ""))
  if (identical(tp, "dup")) pal$dup else if (identical(tp, "del")) pal$del else pal$affected
}

#' Axis labels of the three evidence panels (engines word them differently)
#' @noRd
panel_labels <- function(pd) {
  lab <- pd$labels %||% list()
  list(cov = lab$cov %||% "log2 coverage", ratio = lab$ratio %||% "Observed / expected",
       z = lab$z %||% "Z-score vs references")
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
    pal <- cnv_palette(); aff <- affected_colour(pd)
    cols <- c("Reference samples" = pal$reference, stats::setNames(pal$test, test_lab), "Affected exon(s)" = aff)
    bg <- affected_bg(pd)
    p_cov <- ggplot2::ggplot(mapping = ggplot2::aes(x = px, y = y, color = color_group))
    if (!is.null(bg)) {
        p_cov <- p_cov + ggplot2::geom_rect(data = bg,
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = affected_colour(pd), alpha = 0.08)
    }
    p_cov <- p_cov +
        ggplot2::geom_point(data = d$ref, size = 2.5, alpha = 0.7) +
        ggplot2::geom_point(data = d$test, ggplot2::aes(shape = shape_group), size = 2.5) +
        ggplot2::geom_point(data = d$hit, shape = 17, size = 3.5) +
        ggplot2::scale_shape_manual(values = c("Observed" = 16, "Affected" = 17), guide = "none") +
        ggplot2::scale_colour_manual(
            values = cols,
            guide  = ggplot2::guide_legend(override.aes = list(size = 4), nrow = 1, title = NULL)) +
        ggplot2::labs(y = panel_labels(pd)$cov, x = NULL) +
        cnv_theme(10) +
        ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank(),
                       legend.key = ggplot2::element_rect(fill = "white", colour = NA))
    apply_xaxis_formatting(p_cov, pd)
}

#' Fill colours of the gene tiles (the shared group palette, softened)
#' @noRd
gene_palette <- function(n) {
    if (n == 0) return(character())
    unname(cnv_group_colours(seq_len(n)))
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
        cnv_theme(10) +
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
                             fill = cnv_palette()$band, colour = NA)
    if (!is.null(bg)) {
        p <- p + ggplot2::geom_rect(data = bg,
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = affected_colour(pd), alpha = 0.16)
    }
    p <- p +
        # the "no change" line: ratio 1 for observed/expected, 0 for a log2 ratio
        ggplot2::geom_hline(yintercept = pd$baseline %||% 1, linetype = "dashed",
                            colour = cnv_palette()$normal, linewidth = 0.4) +
        ggplot2::geom_point(ggplot2::aes(y = ratio, color = is_affected, shape = is_affected), size = 3.5) +
        ggplot2::scale_shape_manual(values = c("Observed" = 16, "Affected" = 17), guide = "none") +
        ggplot2::scale_color_manual(
            values = c("Observed" = cnv_palette()$test, "Affected" = affected_colour(pd)),
            guide  = ggplot2::guide_legend(override.aes = list(shape = c(16, 17), size = 3))) +
        ggplot2::labs(x = "", y = panel_labels(pd)$ratio, subtitle = pd$ci_subtitle) +
        cnv_theme(10) +
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
            inherit.aes = FALSE, fill = affected_colour(pd), alpha = 0.08)
    }
    p <- p +
        ggplot2::geom_hline(yintercept = c(-3, 3), linetype = "dashed", colour = cnv_palette()$fail, alpha = 0.7, linewidth = 0.4) +
        ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dotted", colour = cnv_palette()$warn, alpha = 0.8, linewidth = 0.4) +
        ggplot2::geom_hline(yintercept = 0, linetype = "solid", colour = cnv_palette()$text, linewidth = 0.3) +
        # Points only -- no line joins reference or test z-scores across exons.
        ggplot2::geom_point(data = d$ref, ggplot2::aes(x = px, y = y), colour = cnv_palette()$reference, size = 1.5, alpha = 0.7) +
        ggplot2::geom_point(data = d$test, ggplot2::aes(x = px, y = y, color = is_affected, shape = is_affected), size = 3) +
        ggplot2::scale_shape_manual(values = c("Observed" = 16, "Affected" = 17), guide = "none") +
        ggplot2::scale_color_manual(values = c("Observed" = cnv_palette()$test, "Affected" = affected_colour(pd))) +
        ggplot2::coord_cartesian(ylim = c(-z_lim, z_lim)) +
        ggplot2::labs(y = panel_labels(pd)$z, x = NULL) +
        cnv_theme(10) +
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

# ---------------------------------------------------------------------------
# Target annotation shared by the PDF plots and the HTML report
# ---------------------------------------------------------------------------

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

#' Make the chromosome naming of a call table match the BED (\code{chr1} vs \code{1})
#' @param ref_df Target annotation with a \code{chromosome} column.
#' @param target_df Call table.
#' @param col Chromosome column of \code{target_df} (auto-detected).
#' @noRd
harmonise_chr_prefix <- function(ref_df, target_df, col = intersect(c("Chromosome", "CHR", "chromosome"), names(target_df))[1]) {
    if (is.na(col) || !nrow(target_df)) return(target_df)
    has_chr_ref <- any(grepl("^chr", ref_df$chromosome))
    has_chr_tgt <- any(grepl("^chr", target_df[[col]]))
    if (has_chr_ref && !has_chr_tgt) target_df[[col]] <- paste0("chr", target_df[[col]])
    else if (!has_chr_ref && has_chr_tgt) target_df[[col]] <- sub("^chr", "", target_df[[col]])
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
