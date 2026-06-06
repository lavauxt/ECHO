# Local copy of qbetabinom (from ExomeDepth, with permission)
qbetabinom <- function(p, size, rho, prob) {
    a <- prob * (1 - rho) / rho
    b <- (1 - prob) * (1 - rho) / rho
    qbeta(p, a, b) * size
}

#' Harmonise chromosome prefix between BED and CNV data
#'
#' @param ref_df Reference data frame with a \code{chromosome} column (BED file).
#' @param target_df CNV data frame with a \code{Chromosome} column.
#'
#' @return \code{target_df} with \code{Chromosome} prefixes aligned to \code{ref_df}.
#' @noRd
harmonise_chr_prefix <- function(ref_df, target_df) {
  has_chr_ref <- any(grepl("^chr", ref_df$chromosome))
  has_chr_tgt <- any(grepl("^chr", target_df$Chromosome))
  
  if (has_chr_ref && !has_chr_tgt) {
    target_df$Chromosome <- paste0("chr", target_df$Chromosome)
  } else if (!has_chr_ref && has_chr_tgt) {
    target_df$Chromosome <- sub("^chr", "", target_df$Chromosome)
  }
  target_df
}

#' Build exon index labels for plotting
#'
#' @param bed_file BED annotation data frame.
#'
#' @return Character or numeric vector of exon labels, one per row of \code{bed_file}.
#' @noRd
compute_exon_index <- function(bed_file) {
  exon_col <- which(colnames(bed_file) == "exon")
  if (length(exon_col) == 1) {
    bed_file[[exon_col]]
  } else {
    stats::ave(seq_len(nrow(bed_file)), bed_file$gene, FUN = seq_along)
  }
}

#' Prepare coverage and ratio data for a single CNV call
#'
#' @param call_row One-row CNV call data frame.
#' @param counts Coverage count matrix.
#' @param bed_file BED annotation data frame.
#' @param exon_index Exon labels from \code{compute_exon_index()}.
#' @param models Per-sample model parameters from CNV calling.
#' @param refs Per-sample reference sample names from CNV calling.
#'
#' @return Named list of plot inputs, or \code{NULL} if the region is empty.
#' @noRd
prepare_plot_data <- function(call_row, counts, bed_file, exon_index, models, refs) {
  sample    <- call_row$Sample
  idx_start <- as.numeric(call_row$global_start)
  idx_end   <- as.numeric(call_row$global_end)
  target_chr <- normalize_chromosome_vec(call_row$Chromosome, bed_file$chromosome)

  ref_samples <- refs[[sample]]
  if (is.null(ref_samples) || !length(ref_samples)) {
    return(NULL)
  }
  if (is.null(models[[sample]]) || length(models[[sample]]) < 1) {
    return(NULL)
  }

  exon_range_full <- seq(max(1, idx_start - 5), min(nrow(bed_file), idx_end + 5))
  bed_chr <- normalize_chromosome_vec(bed_file$chromosome[exon_range_full], bed_file$chromosome)
  single_chr <- length(unique(bed_chr)) == 1
  prev <- FALSE
  new_chr <- ""
  if (!single_chr) {
    prev <- bed_chr[1] != target_chr
    new_chr <- if (prev) bed_file$chromosome[exon_range_full[1]] else
               bed_file$chromosome[utils::tail(exon_range_full, 1)]
    exon_range <- exon_range_full[bed_chr == target_chr]
  } else {
    exon_range <- exon_range_full
  }
  if (length(exon_range) == 0) return(NULL)

  test_median <- median(counts[exon_range, sample])
  ref_median  <- median(rowSums(counts[exon_range, ref_samples, drop = FALSE]))
  test_log    <- log(pmax(counts[exon_range, sample], 1))
  
  cov_list <- lapply(ref_samples, function(r) {
    scaling <- test_median / median(counts[exon_range, r])
    r_log <- log(pmax(counts[exon_range, r] * scaling, 1))
    data.frame(exon_idx = exon_range, coverage = r_log, group = r,
               color_group = "Reference samples", stringsAsFactors = FALSE)
  })
  cov_list[[length(cov_list) + 1]] <- data.frame(
    exon_idx = exon_range, coverage = test_log, group = "Test sample",
    color_group = "Test sample", stringsAsFactors = FALSE
  )
  cov_data <- do.call(rbind, cov_list)
  
  pt_data <- data.frame(
    exon_idx = exon_range, coverage = test_log,
    color_group = ifelse(exon_range %in% (idx_start:idx_end), "Affected exon(s)", "Test sample")
  )
  pt_data <- pt_data[pt_data$color_group == "Affected exon(s)", ]
  

  ref_counts <- rowSums(counts[exon_range, ref_samples, drop = FALSE])
  test_counts <- counts[exon_range, sample]
  totals <- test_counts + ref_counts
  expected <- ref_counts * (test_median / ref_median)
  expected_safe <- pmax(expected, 1)
  p_expected <- expected_safe / totals
  p_expected <- pmin(pmax(p_expected, 1e-6), 1 - 1e-6)
  ratio <- test_counts / expected_safe
  
  rho <- models[[sample]][1]

  mins <- vapply(seq_along(exon_range), function(i) {
    qbetabinom(0.025, totals[i], max(rho, 0.005), p_expected[i]) /
      expected_safe[i]
  }, numeric(1))
  maxs <- vapply(seq_along(exon_range), function(i) {
    qbetabinom(0.975, totals[i], max(rho, 0.005), p_expected[i]) /
      expected_safe[i]
  }, numeric(1))
  
  ci_data <- data.frame(
    exon = exon_range,
    ratio = ratio,
    lo = mins,
    hi = maxs,
    is_affected = factor(exon_range %in% (idx_start:idx_end),
                         levels = c(FALSE, TRUE),
                         labels = c("Observed", "Affected"))
  )
  
  list(
    cov_data = cov_data,
    pt_data  = pt_data,
    ci_data  = ci_data,
    exon_range = exon_range,
    single_chr = single_chr,
    prev       = prev,
    new_chr    = new_chr,
    sample     = sample,
    idx_start  = idx_start,
    idx_end    = idx_end
  )
}

#' Save a three-panel CNV plot to PDF
#'
#' @param p_cov Coverage ggplot object.
#' @param p_genes Gene-tile ggplot object.
#' @param p_ci Observed/expected ratio ggplot object.
#' @param file_path Output PDF path.
#'
#' @return Invisibly returns \code{NULL}.
#' @noRd
save_cnv_pdf <- function(p_cov, p_genes, p_ci, file_path) {
  grDevices::pdf(file_path, useDingbats = FALSE, width = 8, height = 10)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(6, 1)))
  print(p_cov,   vp = grid::viewport(layout.pos.row = 1:3, layout.pos.col = 1))
  print(p_genes, vp = grid::viewport(layout.pos.row = 4,   layout.pos.col = 1))
  print(p_ci,    vp = grid::viewport(layout.pos.row = 5:6, layout.pos.col = 1))
  grDevices::dev.off()
  invisible(NULL)
}

#' Apply consistent x-axis formatting for exon plots
#'
#' @param p ggplot object.
#' @param single_chr Logical; whether the region spans a single chromosome.
#' @param prev Logical; whether the preceding chromosome should be shown.
#' @param exon_range Integer vector of exon row indices.
#' @param exon_index Exon labels for \code{exon_range}.
#'
#' @return A ggplot object with formatted x-axis.
#' @noRd
apply_xaxis_formatting <- function(p, single_chr, prev, exon_range, exon_index) {
  if (length(exon_range) == 0) return(p)
  min_e <- min(exon_range)
  max_e <- max(exon_range)
  if (single_chr) {
    b <- exon_range
    l <- exon_index[exon_range]
    lim <- NULL
  } else if (prev) {
    b <- (min_e - 6):max_e
    l <- c(rep("", 6), exon_index[exon_range])
    lim <- c(min_e - 6.75, max_e)
  } else {
    b <- min_e:(max_e + 6)
    l <- c(exon_index[exon_range], rep("", 6))
    lim <- c(min_e, max_e + 6.75)
  }
  p +
    ggplot2::scale_x_continuous(breaks = b, labels = l, limits = lim) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, size = 6, hjust = 1))
}

#' Create coverage plot for CNV visualization
#'
#' @param cov_data Long-format coverage data frame.
#' @param pt_data Highlighted affected-exon points data frame.
#' @param single_chr Logical; single-chromosome region flag.
#' @param prev Logical; show preceding chromosome context.
#' @param exon_range Integer vector of exon row indices.
#' @param exon_index Exon labels for the x-axis.
#' @param sample_name Character string to label the test sample.
#'
#' @return A ggplot object.
#' @noRd
create_coverage_plot <- function(cov_data, pt_data, single_chr, prev, exon_range, exon_index, sample_name) {
  cov_data$color_group <- ifelse(cov_data$color_group == "Test sample",
                                 paste0("Test sample (", sample_name, ")"),
                                 cov_data$color_group)

  cols <- c("Reference samples" = "gray",
            setNames("blue", paste0("Test sample (", sample_name, ")")),
            "Affected exon(s)" = "red")

  p_cov <- ggplot2::ggplot() +
    ggplot2::geom_point(data = subset(cov_data, color_group == "Reference samples"),
                        ggplot2::aes(x = exon_idx, y = coverage, color = color_group),
                        size = 2.5, alpha = 0.7) +
    ggplot2::geom_point(data = subset(cov_data, grepl("Test sample", color_group)),
                        ggplot2::aes(x = exon_idx, y = coverage, color = color_group),
                        size = 2.5) +
    ggplot2::geom_point(data = pt_data,
                        ggplot2::aes(x = exon_idx, y = coverage, color = color_group),
                        size = 3.5) +
    ggplot2::scale_colour_manual(values = cols,
                                 guide = ggplot2::guide_legend(override.aes = list(size = 4),
                                                               nrow = 1, title = NULL)) +
    ggplot2::labs(y = "Log (Coverage)", x = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top",
                   legend.title = ggplot2::element_blank(),
                   legend.key = ggplot2::element_rect(fill = "white", colour = NA),
                   panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
                   panel.grid.minor = ggplot2::element_line(color = "grey95", linewidth = 0.15))

  apply_xaxis_formatting(p_cov, single_chr, prev, exon_range, exon_index)
}

#' Create gene tile annotation plot
#'
#' @param bed_file BED annotation data frame.
#' @param exon_range Integer vector of exon row indices.
#' @param single_chr Logical; single-chromosome region flag.
#' @param prev Logical; show preceding chromosome context.
#' @param new_chr Character label for an adjacent chromosome tile.
#'
#' @return A ggplot object.
#' @noRd
create_gene_tile_plot <- function(bed_file, exon_range, single_chr, prev, new_chr) {
  if (length(exon_range) == 0) return(ggplot2::ggplot())
  temp <- cbind(row = seq_len(nrow(bed_file)), bed_file)[exon_range, ]
  gene_names <- unique(bed_file$gene[exon_range])

  n_genes <- length(gene_names)
  if (n_genes == 1) {
    pal <- c("darkblue")
  } else if (n_genes == 2) {
    pal <- c("steelblue", "purple4")
  } else {
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      pal <- RColorBrewer::brewer.pal(min(n_genes, 9), "Purples")
      if (n_genes > 9) pal <- colorRampPalette(pal)(n_genes)
      pal <- rev(pal)
    } else if (requireNamespace("scales", quietly = TRUE)) {
      pal <- scales::hue_pal()(n_genes)
    } else {
      # Fallback to base R colours
      pal <- rainbow(n_genes)
    }
  }
  names(pal) <- gene_names

  gene_tiles <- data.frame(
    gene  = gene_names,
    mid   = as.numeric(sapply(gene_names, function(g) mean(exon_range[temp$gene == g]))),
    width = as.numeric(sapply(gene_names, function(g) sum(temp$gene == g))) - 0.5,
    y     = 1, stringsAsFactors = FALSE
  )
  if (!single_chr) {
    gene_tiles <- rbind(gene_tiles, data.frame(
      gene = new_chr, mid = ifelse(prev, min(exon_range) - 5, max(exon_range) + 5),
      width = 3.5, y = 1, stringsAsFactors = FALSE
    ))
    pal <- c(pal, setNames("gray50", new_chr))
  }

  ggplot2::ggplot(gene_tiles, ggplot2::aes(x = mid, y = y, fill = gene, width = width, label = gene)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = gene), size = 4.5, fontface = "bold", color = "white") +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::theme_bw(base_family = "sans") +
    ggplot2::theme(legend.position = "none",
                   panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   plot.margin = ggplot2::unit(c(.5, .5, .5, .5), "cm"),
                   text = ggplot2::element_text(family = "sans", size = 10)) +
    ggplot2::labs(x = "", y = "")
}

#' Create observed/expected ratio plot with confidence intervals
#'
#' @param ci_data Data frame with \code{exon}, \code{ratio}, \code{lo}, \code{hi},
#'   and \code{is_affected}.
#' @param single_chr Logical; single-chromosome region flag.
#' @param prev Logical; show preceding chromosome context.
#' @param exon_range Integer vector of exon row indices.
#' @param exon_index Exon labels for the x-axis.
#'
#' @return A ggplot object.
#' @noRd
create_ci_plot <- function(ci_data, single_chr, prev, exon_range, exon_index) {
  p <- ggplot2::ggplot(ci_data, ggplot2::aes(x = exon)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey80", colour = NA) +
    ggplot2::geom_point(ggplot2::aes(y = ratio, color = is_affected), size = 3.5) +
    ggplot2::scale_color_manual(values = c("Observed" = "blue", "Affected" = "red"),
                                guide = ggplot2::guide_legend(override.aes = list(shape = 19, size = 3))) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none", legend.title = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
                   panel.grid.minor = ggplot2::element_line(color = "grey95", linewidth = 0.15)) +
    ggplot2::labs(x = "", y = "Observed/Expected")
  apply_xaxis_formatting(p, single_chr, prev, exon_range, exon_index)
}

#' Generate CNV Detection Plots
#'
#' Creates PDF plots summarising CNV calls from ExomeDepth‑style workflows.
#' Each PDF contains three panels: coverage comparison, gene annotation tiles,
#' and observed/expected ratio with beta‑binomial confidence intervals.
#'
#' @param rdata_file Character string. Path to summary RData.
#' @param output_dir Character string. Output directory for PDF files.
#' @param modechrom Chromosome filter.
#' @param prefix Filename prefix.
#' @param log_file Optional path to log file.
#'
#' @return Invisibly returns the number of PDFs written.
#' @export
generate_plots <- function(rdata_file,
                           output_dir = "./plots",
                           modechrom = "A",
                           prefix = NULL,
                           log_file = NULL) {
    
    log_msg <- function(msg, type = "INFO") {
        if (!is.null(log_file)) {
            timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
            cat(paste0("[", type, "] ", timestamp, " ", msg, "\n"), file = log_file, append = TRUE)
        }
        message(msg)
    }
    
    log_msg("[INFO] BEGIN plot generation")
    
    objs <- load_rdata(
        rdata_file,
        required = c("cnv_calls", "counts", "bed_file", "models", "refs")
    )
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

    inc <- if (modechrom %in% c("XX", "XY")) "chrX" else NULL
    exc <- if (modechrom == "A") c("chrX", "chrY") else NULL
    cnv_plot <- filter_chromosomes(cnv_calls, include = inc, exclude = exc)

    if (nrow(cnv_plot) == 0) {
        log_msg("[INFO] No CNV calls remain after chromosome filtering.")
        return(invisible(0))
    }

    cnv_plot <- harmonise_chr_prefix(bed_file, cnv_plot)
    exon_index <- compute_exon_index(bed_file)
    n_written <- 0L

    for (i in seq_len(nrow(cnv_plot))) {
        tryCatch({
            call_row <- cnv_plot[i, ]
            sample <- call_row$Sample
            if (is.null(refs[[sample]]) || !length(refs[[sample]])) {
                log_msg(paste("[WARNING] Skipping plot for", sample, "(call", i, "): no reference samples"), "WARNING")
                next
            }
            if (is.null(models[[sample]])) {
                log_msg(paste("[WARNING] Skipping plot for", sample, "(call", i, "): missing model parameters"), "WARNING")
                next
            }
            plot_data <- prepare_plot_data(call_row, counts, bed_file, exon_index, models, refs)
            if (is.null(plot_data)) {
                log_msg(paste("[WARNING] Skipping plot for", sample, "(call", i, "): empty exon window"), "WARNING")
                next
            }

            gene_str <- sanitize_filename(as.character(call_row$Gene))
            file_path <- file.path(output_dir, paste0(prefix_str, ".", plot_data$sample, ".", gene_str, "_", i, ".pdf"))

            # Only create PDF if we have data
            p_cov <- create_coverage_plot(plot_data$cov_data, plot_data$pt_data,
                                          plot_data$single_chr, plot_data$prev,
                                          plot_data$exon_range, exon_index,
                                          sample_name = plot_data$sample)
            p_genes <- create_gene_tile_plot(bed_file, plot_data$exon_range,
                                             plot_data$single_chr, plot_data$prev,
                                             plot_data$new_chr)
            p_ci    <- create_ci_plot(plot_data$ci_data, plot_data$single_chr,
                                      plot_data$prev, plot_data$exon_range, exon_index)

            save_cnv_pdf(p_cov, p_genes, p_ci, file_path)
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