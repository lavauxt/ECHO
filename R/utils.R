#' Stop if value is NULL or empty
#'
#' @param val Object to validate.
#' @param msg Error message to display if validation fails.
#' @export
stop_if_missing <- function(val, msg) {
    if (is.null(val) || length(val) == 0) stop(msg)
}

#' Stop if file does not exist
#'
#' @param path Path to the file.
#' @param msg Error message to display if file is missing.
#' @export
stop_if_not_file <- function(path, msg) {
    if (is.null(path) || !file.exists(path)) stop(msg)
}

#' Load an RData file into an isolated environment
#'
#' @param path Path to the \code{.RData} file.
#' @param required Optional character vector of required object names.
#' @return A named list of objects from the file.
#' @noRd
load_rdata <- function(path, required = NULL) {
    stop_if_not_file(path, paste0("[ERROR] RData file not found: ", path))
    env <- new.env()
    load(path, envir = env)
    if (!is.null(required)) {
        missing <- setdiff(required, ls(env))
        if (length(missing) > 0) {
            stop("[ERROR] Missing variables in RData: ", paste(missing, collapse = ", "),
                 "\nAvailable objects: ", paste(ls(env), collapse = ", "))
        }
    }
    as.list(env)
}

#' Check that a BAM file has an index on disk
#'
#' @param bam Path to a BAM file.
#' @return Logical.
#' @noRd
bam_has_index <- function(bam) {
    file.exists(paste0(bam, ".bai")) || file.exists(paste0(bam, ".bam.bai"))
}

#' Normalise a chromosome name to match a reference naming style (vectorised)
#'
#' @param chr_vec Character vector of chromosome names.
#' @param ref_chromosomes Character vector of reference chromosome names (e.g. BED).
#' @return Character vector with prefix aligned to reference style.
#' @noRd
normalize_chromosome_vec <- function(chr_vec, ref_chromosomes) {
    has_chr_ref <- any(grepl("^chr", ref_chromosomes))
    if (has_chr_ref) {
        ifelse(grepl("^chr", chr_vec), chr_vec, paste0("chr", chr_vec))
    } else {
        ifelse(grepl("^chr", chr_vec), sub("^chr", "", chr_vec), chr_vec)
    }
}

#' Filter a data frame by chromosome column
#'
#' @param df Input data frame.
#' @param include Optional vector of chromosomes to include.
#' @param exclude Optional vector of chromosomes to exclude.
#'
#' @return Filtered data frame (empty if input is empty).
#' @export
filter_chromosomes <- function(df, include = NULL, exclude = NULL) {
    if (is.null(df) || nrow(df) == 0) return(df)
    chrom_col <- intersect(c("chromosome", "Chromosome"), colnames(df))[1]
    if (is.na(chrom_col)) {
        stop("Data frame must contain a 'chromosome' or 'Chromosome' column.")
    }
    norm <- function(x) {
        unique(c(x, sub("^chr", "", x), paste0("chr", sub("^chr", "", x))))
    }
    if (!is.null(include)) {
        df <- df[df[[chrom_col]] %in% norm(include), ]
    }
    if (!is.null(exclude)) {
        df <- df[!df[[chrom_col]] %in% norm(exclude), ]
    }
    df
}

#' Null coalescing operator
#'
#' Returns the first argument if not NULL, otherwise the second.
#'
#' @param a Primary value.
#' @param b Fallback value.
#' @export
`%||%` <- function(a, b) {
    if (!is.null(a)) a else b
}

#' @import data.table
NULL

#' Sanitize a string for use as a filename
#'
#' @param name Character string to sanitize.
#' @export
sanitize_filename <- function(name) {
    gsub("[^[:alnum:]]", "_", name)
}

#' Convert global exon indices to within‑gene indices (using start‑order)
#'
#' For a set of CNV calls (each with global start.p and end.p), this function
#' replaces those columns with per‑gene exon numbers. The original global indices
#' are preserved as `global_start` and `global_end`.
#'
#' @param cnv_calls Data frame of CNV calls (from ExomeDepth).
#' @param bed_file BED annotation data frame (must have columns `gene`, `start`).
#'
#' @return The same data frame with updated `start.p` / `end.p` (per‑gene) and
#'   added `global_start` / `global_end` columns.
#' @export
add_within_gene_indices <- function(cnv_calls, bed_file) {
    if (nrow(cnv_calls) == 0) return(cnv_calls)
    cnv_calls$global_start <- as.numeric(cnv_calls$start.p)
    cnv_calls$global_end   <- as.numeric(cnv_calls$end.p)
    exon_in_gene <- ave(seq_len(nrow(bed_file)), bed_file$gene,
                        FUN = function(idx) rank(bed_file$start[idx], ties.method = "first"))
    cnv_calls$start.p <- exon_in_gene[cnv_calls$global_start]
    cnv_calls$end.p   <- exon_in_gene[cnv_calls$global_end]
    cnv_calls
}

#' Compute exon index for plotting (consistent with add_within_gene_indices)
#'
#' @param bed_file BED annotation data frame.
#' @return Integer vector of per‑gene exon numbers, one per row.
#' @export
compute_exon_index <- function(bed_file) {
    if ("exon_number" %in% colnames(bed_file) && any(!is.na(bed_file$exon_number))) {
        return(as.integer(bed_file$exon_number))
    }

    if ("exon" %in% colnames(bed_file) && any(!is.na(bed_file$exon))) {
        parsed_exon <- suppressWarnings(as.integer(gsub("[^0-9]", "", bed_file$exon)))
        if (any(!is.na(parsed_exon))) return(parsed_exon)
    }

    ave(seq_len(nrow(bed_file)), bed_file$gene,
        FUN = function(idx) rank(bed_file$start[idx], ties.method = "first"))
}

#' Load and validate YAML configuration
#'
#' @param yaml_path Character string. Path to a YAML file, or a filename resolved
#'   via \code{system.file()} when not found locally.
#'
#' @return Named list with \code{input}, \code{output}, and optional \code{settings}
#'   sections from the YAML file.
#' @export
load_config <- function(yaml_path) {
    if (!file.exists(yaml_path)) {
        found_path <- system.file(yaml_path, package = "ECHO")
        if (found_path == "") {
            stop("[ERROR] Config file not found at: ", yaml_path)
        }
        yaml_path <- found_path
    }
    cfg <- yaml::read_yaml(yaml_path)
    required <- c("input", "output")
    if (!all(required %in% names(cfg))) {
        stop("[ERROR] YAML must contain 'input' and 'output' sections")
    }
    return(cfg)
}

#' Stop if value is NULL or empty
#'
#' @param val Object to validate.
#' @param msg Error message to display if validation fails.
#' @export
stop_if_missing <- function(val, msg) {
    if (is.null(val) || length(val) == 0) stop(msg)
}

#' Stop if file does not exist
#'
#' @param path Path to the file.
#' @param msg Error message to display if file is missing.
#' @export
stop_if_not_file <- function(path, msg) {
    if (is.null(path) || !file.exists(path)) stop(msg)
}

#' Load an RData file into an isolated environment
#'
#' @param path Path to the \code{.RData} file.
#' @param required Optional character vector of required object names.
#' @return A named list of objects from the file.
#' @noRd
load_rdata <- function(path, required = NULL) {
    stop_if_not_file(path, paste0("[ERROR] RData file not found: ", path))
    env <- new.env()
    load(path, envir = env)
    if (!is.null(required)) {
        missing <- setdiff(required, ls(env))
        if (length(missing) > 0) {
            stop("[ERROR] Missing variables in RData: ", paste(missing, collapse = ", "),
                 "\nAvailable objects: ", paste(ls(env), collapse = ", "))
        }
    }
    as.list(env)
}

#' Check that a BAM file has an index on disk
#'
#' @param bam Path to a BAM file.
#' @return Logical.
#' @noRd
bam_has_index <- function(bam) {
    file.exists(paste0(bam, ".bai")) || file.exists(paste0(bam, ".bam.bai"))
}

#' Normalise a chromosome name to match a reference naming style (vectorised)
#'
#' @param chr_vec Character vector of chromosome names.
#' @param ref_chromosomes Character vector of reference chromosome names (e.g. BED).
#' @return Character vector with prefix aligned to reference style.
#' @noRd
normalize_chromosome_vec <- function(chr_vec, ref_chromosomes) {
    has_chr_ref <- any(grepl("^chr", ref_chromosomes))
    if (has_chr_ref) {
        ifelse(grepl("^chr", chr_vec), chr_vec, paste0("chr", chr_vec))
    } else {
        ifelse(grepl("^chr", chr_vec), sub("^chr", "", chr_vec), chr_vec)
    }
}

#' Filter a data frame by chromosome column
#'
#' @param df Input data frame.
#' @param include Optional vector of chromosomes to include.
#' @param exclude Optional vector of chromosomes to exclude.
#'
#' @return Filtered data frame (empty if input is empty).
#' @export
filter_chromosomes <- function(df, include = NULL, exclude = NULL) {
    if (is.null(df) || nrow(df) == 0) return(df)
    chrom_col <- intersect(c("chromosome", "Chromosome"), colnames(df))[1]
    if (is.na(chrom_col)) {
        stop("Data frame must contain a 'chromosome' or 'Chromosome' column.")
    }
    norm <- function(x) {
        unique(c(x, sub("^chr", "", x), paste0("chr", sub("^chr", "", x))))
    }
    if (!is.null(include)) {
        df <- df[df[[chrom_col]] %in% norm(include), ]
    }
    if (!is.null(exclude)) {
        df <- df[!df[[chrom_col]] %in% norm(exclude), ]
    }
    df
}

#' Null coalescing operator
#'
#' Returns the first argument if not NULL, otherwise the second.
#'
#' @param a Primary value.
#' @param b Fallback value.
#' @export
`%||%` <- function(a, b) {
    if (!is.null(a)) a else b
}

#' @import data.table
NULL

#' Sanitize a string for use as a filename
#'
#' @param name Character string to sanitize.
#' @export
sanitize_filename <- function(name) {
    gsub("[^[:alnum:]]", "_", name)
}

#' Assign sequential exon numbers within each gene
#'
#' Sorts rows by chromosome, start, end and then within each gene assigns
#' 1..n based on genomic order. Duplicate intervals (same chrom, start, end)
#' are removed and a warning is issued.
#'
#' @param bed_file data.frame with columns 'chromosome', 'start', 'end', 'gene'
#' @return A modified data.frame with added column 'exon_number'
#' @noRd
assign_exon_numbers_per_gene <- function(bed_file) {
  # Ensure we have a standard chromosome column name
  chrom_col <- if ("chromosome" %in% names(bed_file)) "chromosome" else "Chr"
  start_col <- if ("start" %in% names(bed_file)) "start" else "Start"
  end_col   <- if ("end"   %in% names(bed_file)) "end"   else "End"
  gene_col  <- if ("gene"  %in% names(bed_file)) "gene"  else "Gene"
  
  # Check for required columns
  stopifnot(all(c(chrom_col, start_col, end_col, gene_col) %in% names(bed_file)))
  
  # Convert to data.table for efficiency and easier grouping
  dt <- data.table::as.data.table(bed_file)
  setnames(dt, c(chrom_col, start_col, end_col, gene_col),
           c("chrom", "start", "end", "gene"))
  
  # Remove exact duplicates (same chrom, start, end, gene)
  dup_rows <- duplicated(dt, by = c("chrom", "start", "end", "gene"))
  if (any(dup_rows)) {
    warning(sprintf("Removed %d duplicate rows (identical chromosome, start, end, gene).",
                    sum(dup_rows)), immediate. = TRUE)
    dt <- dt[!dup_rows]
  }
  
  # Sort by chromosome (ensuring natural order), then start, then end
  chrom_levels <- c(paste0("chr", c(1:22, "X", "Y", "M")),
                    c(as.character(1:22), "X", "Y", "M"))
  dt[, chrom_fac := factor(chrom, levels = unique(c(chrom_levels, unique(chrom))))]
  setorder(dt, chrom_fac, start, end)
  dt[, chrom_fac := NULL]
  
  # Assign exon numbers per gene
  dt[, exon_number := seq_len(.N), by = "gene"]
  
  # Restore original column names
  setnames(dt, c("chrom", "start", "end", "gene"),
           c(chrom_col, start_col, end_col, gene_col))
  as.data.frame(dt)
}

#' Convert global exon indices to within‑gene indices (using start‑order)
#'
#' For a set of CNV calls (each with global start.p and end.p), this function
#' replaces those columns with per‑gene exon numbers. The original global indices
#' are preserved as `global_start` and `global_end`.
#'
#' @param cnv_calls Data frame of CNV calls (from ExomeDepth).
#' @param bed_file BED annotation data frame (must have columns `gene`, `start`).
#'
#' @return The same data frame with updated `start.p` / `end.p` (per‑gene) and
#'   added `global_start` / `global_end` columns.
#' @export
add_within_gene_indices <- function(cnv_calls, bed_file) {
  if (nrow(cnv_calls) == 0) return(cnv_calls)
  
  # Make a stable copy and ensure bed_file has exon numbers
  bed_sorted <- assign_exon_numbers_per_gene(bed_file)
  exon_numbers <- bed_sorted$exon_number
  
  cnv_calls$global_start <- as.integer(cnv_calls$start.p)
  cnv_calls$global_end   <- as.integer(cnv_calls$end.p)
  
  cnv_calls$start.p <- exon_numbers[cnv_calls$global_start]
  cnv_calls$end.p   <- exon_numbers[cnv_calls$global_end]
  
  cnv_calls
}

#' Compute exon index for plotting (consistent with add_within_gene_indices)
#'
#' @param bed_file BED annotation data frame.
#' @return Integer vector of per‑gene exon numbers, one per row.
#' @export
compute_exon_index <- function(bed_file) {
  # Reuse the same numbering logic
  if ("exon_number" %in% names(bed_file) &&
      all(!is.na(bed_file$exon_number)) &&
      length(unique(bed_file$exon_number)) > 0) {
    return(bed_file$exon_number)
  }
  bed_numbered <- assign_exon_numbers_per_gene(bed_file)
  bed_numbered$exon_number
}

#' Plot PCA of Sample Coverage Profiles
#'
#' @param counts Data frame of read counts (rows = exons, columns = samples).
#' @param sample_names Character vector of sample names.
#' @param output_pdf Path to save the PCA plot (or NULL for direct display).
#' @param color_by Optional vector of sample groups (e.g., gender, batch).
#' @param scale Logical, whether to scale the data (recommended).
#' @return Invisibly returns a list with PCA object and variance explained.
#' @export
plot_coverage_pca <- function(counts, sample_names, output_pdf = NULL,
                              color_by = NULL, scale = TRUE) {
  # Extract only count columns
  count_mat <- as.matrix(counts[, sample_names, drop = FALSE])
  # Log transform and remove zero-variance rows
  log_mat <- log2(count_mat + 1)
  row_var <- apply(log_mat, 1, var, na.rm = TRUE)
  keep <- row_var > 0 & !is.na(row_var)
  if (sum(keep) < 2) {
    warning("Insufficient variation for PCA (fewer than 2 informative rows).")
    return(invisible(NULL))
  }
  log_mat <- log_mat[keep, ]
  
  pca <- prcomp(t(log_mat), scale. = scale, center = TRUE)
  var_exp <- summary(pca)$importance[2, ] * 100
  pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Sample = sample_names)
  
  # Base plot
  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, label = Sample))
  
  if (!is.null(color_by) && length(color_by) == nrow(pca_df)) {
    pca_df$Group <- as.factor(color_by)
    p <- p + ggplot2::aes(colour = Group)
    # Add dashed ellipses (always drawn when grouping exists)
    if (length(unique(color_by)) >= 2) {
      p <- p + ggplot2::stat_ellipse(ggplot2::aes(colour = Group), 
                                     type = "norm", linetype = "dashed")
    }
  } else {
    p <- p + ggplot2::aes(colour = "Sample")
  }
  
  p <- p + 
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
                  y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
                  title = "PCA of Sample Coverage Profiles") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      # Axis titles in blue
      axis.title.x = ggplot2::element_text(color = "blue", face = "bold"),
      axis.title.y = ggplot2::element_text(color = "blue", face = "bold"),
      # Axis tick labels (numbers) in blue
      axis.text.x = ggplot2::element_text(color = "blue"),
      axis.text.y = ggplot2::element_text(color = "blue"),
      # Plot title in blue
      plot.title = ggplot2::element_text(color = "blue", face = "bold", hjust = 0.5),
      # Legend text (if present) in dark blue
      legend.text = ggplot2::element_text(color = "darkblue"),
      legend.title = ggplot2::element_text(color = "darkblue")
    )
  
  # Handle labels without overlapping
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(size = 3, max.overlaps = 15)
  } else {
    p <- p + ggplot2::geom_text(check_overlap = TRUE, size = 3, hjust = -0.2, vjust = 0.5)
  }
  
  # Remove colour legend if no grouping was provided
  if (is.null(color_by)) {
    p <- p + ggplot2::theme(legend.position = "none")
  } else {
    p <- p + ggplot2::guides(colour = ggplot2::guide_legend(title = "Group"))
  }
  
  if (!is.null(output_pdf)) {
    ggplot2::ggsave(output_pdf, p, width = 8, height = 6)
    message("[INFO] PCA plot saved: ", output_pdf)
  } else {
    print(p)
  }
  invisible(list(pca = pca, var_exp = var_exp))
}