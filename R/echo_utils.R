# =============================================================================
# utils.R — ECHO pipeline utilities
# =============================================================================

#' Stop if value is NULL or empty
#' @param val Object to validate.
#' @param msg Error message if validation fails.
#' @export
stop_if_missing <- function(val, msg) {
    if (is.null(val) || length(val) == 0) stop(msg)
}

#' Stop if file does not exist
#' @param path Path to the file.
#' @param msg Error message if file is missing.
#' @export
stop_if_not_file <- function(path, msg) {
    if (is.null(path) || !file.exists(path)) stop(msg)
}

#' Load an RData file into an isolated environment
#' @param path Path to the .RData file.
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
#' @param bam Path to a BAM file.
#' @return Logical.
#' @noRd
bam_has_index <- function(bam) {
    file.exists(paste0(bam, ".bai")) || file.exists(paste0(bam, ".bam.bai"))
}

#' Normalise a chromosome name to match a reference naming style (vectorised)
#' @param chr_vec Character vector of chromosome names.
#' @param ref_chromosomes Character vector of reference chromosome names.
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
#' @param df Input data frame.
#' @param include Optional vector of chromosomes to include.
#' @param exclude Optional vector of chromosomes to exclude.
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

#' Quantile function of the Beta-Binomial distribution
#'
#' Local copy of \code{ExomeDepth:::qbetabinom} (from ExomeDepth v1.1.15,
#' with permission), kept here as the single shared copy so that both
#' \code{plots.R} and \code{ECHO_global_report.Rmd} depend on ECHO's own
#' maintained code rather than an unexported ExomeDepth internal accessed
#' via \code{:::}, which offers no stability guarantee across ExomeDepth
#' versions.
#'
#' @param p Numeric vector of probabilities.
#' @param size Numeric vector of trial totals.
#' @param rho Numeric overdispersion parameter.
#' @param prob Numeric vector of expected success probabilities.
#' @return Numeric vector of quantiles.
#' @noRd
qbetabinom <- function(p, size, rho, prob) {
    a <- prob * (1 - rho) / rho
    b <- (1 - prob) * (1 - rho) / rho
    qbeta(p, a, b) * size
}

#' Null coalescing operator
#' @param a Primary value.
#' @param b Fallback value.
#' @export
`%||%` <- function(a, b) {
    if (!is.null(a)) a else b
}

#' @import data.table
NULL

#' Sanitize a string for use as a filename
#' @param name Character string to sanitize.
#' @export
sanitize_filename <- function(name) {
    gsub("[^[:alnum:]]", "_", name)
}

#' Assign sequential exon numbers within each gene
#'
#' Sorts rows by chromosome, start, end and then within each gene assigns
#' 1..n based on genomic order. Duplicate intervals (same chrom, start, end,
#' gene) share the \code{exon_number} of their first occurrence but are
#' \strong{not} dropped from the output.
#'
#' BUGFIX: this previously did `dt <- dt[!dup_rows]` whenever duplicates were
#' found, silently shrinking the row count -- directly contradicting this
#' function's own contract (every caller binds `exon_number` back onto
#' `bed_file`/`counts` purely by position: `counts$exon_number <-
#' bed_file$exon_number` in bam_coverage.R, `exon_in_gene[fail_exon]` in
#' metrics.R, `exon_numbers[cnv_calls$global_start]` in
#' add_within_gene_indices()). Any duplicate row shifted every downstream
#' target's exon/gene label by one and, in `add_within_gene_indices()`,
#' misaligned the CNV-to-exon index lookup outright. Fixed (matching
#' CANOPE's `assign_exon_numbers_per_gene()`) to number only the unique
#' (chrom, start, end, gene) combinations internally, then map the result
#' back onto *every* original row -- so `nrow(output) == nrow(bed_file)`
#' always holds, regardless of duplicates.
#'
#' @param bed_file data.frame with columns chromosome/Chr, start/Start, end/End, gene/Gene
#' @return A modified data.frame (original row order and row count preserved)
#'   with added column exon_number
#' @noRd
assign_exon_numbers_per_gene <- function(bed_file) {
    chrom_col <- if ("chromosome" %in% names(bed_file)) "chromosome" else "Chr"
    start_col <- if ("start"      %in% names(bed_file)) "start"      else "Start"
    end_col   <- if ("end"        %in% names(bed_file)) "end"        else "End"
    gene_col  <- if ("gene"       %in% names(bed_file)) "gene"       else "Gene"

    stopifnot(all(c(chrom_col, start_col, end_col, gene_col) %in% names(bed_file)))

    n_in <- nrow(bed_file)

    dt <- data.table::as.data.table(bed_file)
    dt[, .orig_row := .I]  # remember incoming row order so we never permute it
    data.table::setnames(dt, c(chrom_col, start_col, end_col, gene_col),
                         c("chrom", "start", "end", "gene"))
    dt[, .key := paste(chrom, start, end, gene, sep = "\r")]

    dup_rows <- duplicated(dt, by = ".key")
    if (any(dup_rows)) {
        warning(sprintf(
            "Found %d duplicate row(s) (identical chromosome, start, end, gene); duplicates share the exon_number of their first occurrence (rows are NOT dropped).",
            sum(dup_rows)), immediate. = TRUE)
    }

    # Number only the unique (chrom, start, end, gene) combinations -- a
    # duplicate row must not get its own exon_number bumped from the count,
    # or a gene with a repeated interval would appear to have more exons
    # than it does.
    dt_unique <- dt[!dup_rows]

    # IMPORTANT: bed_file/counts elsewhere in the pipeline are ordered by
    # FASTA contig order (see bam_coverage.R), which need not match this
    # hardcoded karyotype order. Any caller downstream relies on this
    # function returning rows in the SAME order they came in so that
    # positional assignments like `counts$exon_number <-
    # bed_file$exon_number` or `exon_numbers[cnv_calls$global_start]` stay
    # aligned. So we compute the numbering on a sorted *copy* of the row
    # order and map the result back onto the original row order, rather
    # than sorting dt itself.
    chrom_levels <- c(paste0("chr", c(1:22, "X", "Y", "M")),
                      c(as.character(1:22), "X", "Y", "M"))
    dt_unique[, .chrom_fac := factor(chrom, levels = unique(c(chrom_levels, unique(chrom))))]
    data.table::setorder(dt_unique, .chrom_fac, start, end)
    dt_unique[, .chrom_fac := NULL]

    dt_unique[, exon_number := seq_len(.N), by = "gene"]

    # Map exon_number back onto every original row (including duplicates)
    # by key, then restore the original row order. Guarantees
    # nrow(output) == nrow(bed_file) always.
    dt[, exon_number := dt_unique$exon_number[match(.key, dt_unique$.key)]]

    data.table::setorder(dt, .orig_row)  # restore original (input) row order
    dt[, c(".orig_row", ".key") := NULL]

    data.table::setnames(dt, c("chrom", "start", "end", "gene"),
                         c(chrom_col, start_col, end_col, gene_col))
    out <- as.data.frame(dt)
    stopifnot(nrow(out) == n_in)
    out
}

#' Convert global exon indices to within-gene indices (using start-order)
#'
#' @param cnv_calls Data frame of CNV calls (from ExomeDepth).
#' @param bed_file BED annotation data frame (must have columns gene, start).
#' @return The same data frame with updated start.p / end.p (per-gene) and
#'   added global_start / global_end columns.
#' @export
add_within_gene_indices <- function(cnv_calls, bed_file) {
    if (nrow(cnv_calls) == 0) return(cnv_calls)

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
#' @return Integer vector of per-gene exon numbers, one per row.
#' @export
compute_exon_index <- function(bed_file) {
    if ("exon_number" %in% names(bed_file) &&
        all(!is.na(bed_file$exon_number)) &&
        length(unique(bed_file$exon_number)) > 0) {
        return(bed_file$exon_number)
    }
    bed_numbered <- assign_exon_numbers_per_gene(bed_file)
    bed_numbered$exon_number
}

#' Load and validate YAML configuration
#'
#' @param yaml_path Character string. Path to a YAML file, or a filename
#'   resolved via system.file() when not found locally.
#' @return Named list with input, output, and optional settings sections.
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
    count_mat <- as.matrix(counts[, sample_names, drop = FALSE])
    log_mat   <- log2(count_mat + 1)
    row_var   <- apply(log_mat, 1, var, na.rm = TRUE)
    keep      <- row_var > 0 & !is.na(row_var)
    if (sum(keep) < 2) {
        warning("Insufficient variation for PCA (fewer than 2 informative rows).")
        return(invisible(NULL))
    }
    log_mat <- log_mat[keep, ]

    pca     <- prcomp(t(log_mat), scale. = scale, center = TRUE)
    var_exp <- summary(pca)$importance[2, ] * 100
    pca_df  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Sample = sample_names)

    p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, label = Sample))

    if (!is.null(color_by) && length(color_by) == nrow(pca_df)) {
        pca_df$Group <- as.factor(color_by)
        p <- p + ggplot2::aes(colour = Group)
        if (length(unique(color_by)) >= 2) {
            p <- p + ggplot2::stat_ellipse(ggplot2::aes(colour = Group),
                                           type = "norm", linetype = "dashed")
        }
    }
    p <- p +
        ggplot2::geom_point(size = 3,
                            colour = if (is.null(color_by)) "steelblue" else NULL) +
        ggplot2::labs(x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
                      y = paste0("PC2 (", round(var_exp[2], 1), "%)"),
                      title = "PCA of Sample Coverage Profiles") +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.title.x  = ggplot2::element_text(color = "blue",  face = "bold"),
            axis.title.y  = ggplot2::element_text(color = "blue",  face = "bold"),
            axis.text.x   = ggplot2::element_text(color = "blue"),
            axis.text.y   = ggplot2::element_text(color = "blue"),
            plot.title    = ggplot2::element_text(color = "blue",  face = "bold", hjust = 0.5),
            legend.text   = ggplot2::element_text(color = "darkblue"),
            legend.title  = ggplot2::element_text(color = "darkblue")
        )

    if (requireNamespace("ggrepel", quietly = TRUE)) {
        p <- p + ggrepel::geom_text_repel(size = 3, max.overlaps = 15)
    } else {
        p <- p + ggplot2::geom_text(check_overlap = TRUE, size = 3, hjust = -0.2, vjust = 0.5)
    }

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

#' Validate a BED file for common structural and coordinate errors
#'
#' Lightweight, non-fatal sanity checks run on a BED file before it enters
#' the ECHO pipeline: malformed/inverted coordinates, duplicate regions,
#' and (when present) a consistency check between the BED's own
#' start/end columns and any \code{"_chrN_start_end"}-style coordinate
#' suffix embedded in the name column (a convention used by some
#' capture-panel design exports, e.g. \code{"..._chr13_32889589_32889829"}).
#' That last check compares the embedded 1-based coordinate against
#' \code{start + 1} (the standard 0-based -> 1-based conversion). If most
#' rows share the exact same non-zero offset, that is reported as a likely
#' systematic authoring bug in the BED file (e.g. the export script
#' subtracting one base too many) rather than N independent bad rows.
#'
#' This function never stops the pipeline — only warns — since some
#' findings (e.g. multi-segment merged exons) are expected quirks rather
#' than fatal errors.
#'
#' @param bed_path Character string. Path to a BED file.
#' @param verbose Logical. Print a summary of findings. Default TRUE.
#' @return Invisibly, a named list of diagnostic data frames (empty list
#'   if no issues were found).
#' @export
validate_bed_regions <- function(bed_path, verbose = TRUE) {
    bed <- tryCatch(data.table::fread(bed_path, header = FALSE, sep = "\t"),
                     error = function(e) NULL)
    if (is.null(bed) || nrow(bed) == 0) {
        if (verbose) message("[WARNING] Could not read BED file for validation: ", bed_path)
        return(invisible(list()))
    }
    if (ncol(bed) < 3) {
        if (verbose) message("[WARNING] BED file has fewer than 3 columns; skipping validation.")
        return(invisible(list()))
    }
    data.table::setnames(bed, 1:3, c("chrom", "start", "end"))
    if (ncol(bed) >= 4) data.table::setnames(bed, 4, "name") else bed$name <- NA_character_

    issues <- list()
    if (verbose) message("[INFO] Validating BED file: ", bed_path, " (", nrow(bed), " rows)")

    bad_coords <- which(!is.na(bed$start) & !is.na(bed$end) & (bed$end <= bed$start | bed$start < 0))
    if (length(bad_coords)) {
        issues$invalid_coordinates <- as.data.frame(bed[bad_coords, ])
        if (verbose) message("[WARNING] ", length(bad_coords),
                              " row(s) have end <= start or a negative start (e.g. row ",
                              bad_coords[1], ").")
    }

    dup_rows <- which(duplicated(bed[, c("chrom", "start", "end")]))
    if (length(dup_rows)) {
        issues$duplicate_regions <- as.data.frame(bed[dup_rows, ])
        if (verbose) message("[WARNING] ", length(dup_rows), " duplicate chrom/start/end row(s) found.")
    }

    # Optional: cross-check any "_chrN_start_end"-style coordinate embedded
    # in the name column against the BED's own columns.
    pat <- "chr\\w+_([0-9]+)_([0-9]+)"
    has_pattern <- grepl(pat, bed$name)
    if (any(has_pattern)) {
        sub_bed <- bed[has_pattern, ]
        all_matches <- regmatches(sub_bed$name, gregexpr(pat, sub_bed$name))
        embedded_start <- vapply(all_matches, function(m) {
            if (!length(m)) return(NA_real_)
            min(as.numeric(sub(paste0("^", pat, "$"), "\\1", m)))
        }, numeric(1))
        embedded_end <- vapply(all_matches, function(m) {
            if (!length(m)) return(NA_real_)
            max(as.numeric(sub(paste0("^", pat, "$"), "\\2", m)))
        }, numeric(1))
        expected_start <- sub_bed$start + 1
        start_diff <- embedded_start - expected_start
        end_diff   <- embedded_end - sub_bed$end
        off <- which(start_diff != 0 | end_diff != 0)
        if (length(off)) {
            issues$name_coordinate_mismatch <- as.data.frame(sub_bed[off, ])
            pct <- round(100 * length(off) / nrow(sub_bed), 1)
            if (verbose) {
                message("[WARNING] ", length(off), " of ", nrow(sub_bed), " (", pct, "%) row(s) with a ",
                        "name-embedded coordinate disagree with the BED start/end columns.")
                tab_start <- table(start_diff[off])
                dominant_shift <- as.numeric(names(tab_start)[which.max(tab_start)])
                dominant_frac  <- max(tab_start) / length(off)
                if (dominant_frac > 0.9 && dominant_shift != 0) {
                    message("[WARNING]   ", round(dominant_frac * 100, 1), "% of the mismatched rows share the ",
                            "exact same offset (", dominant_shift, " bp on the start coordinate). A single, ",
                            "uniform shift across nearly every row usually means the script that generated ",
                            "this BED file has a systematic off-by-one/two error (its 'Start' column is ",
                            abs(dominant_shift), " bp too ", ifelse(dominant_shift > 0, "small", "large"),
                            " relative to the position encoded in the name column) rather than isolated bad rows.")
                }
            }
        }
    }

    if (verbose && length(issues) == 0) message("[INFO] BED validation found no issues.")
    invisible(issues)
}

#' Flag Background-Exon Calibration Issues for a CNV Call Window
#'
#' Ported from CANOPE (a sibling pipeline sharing this reporting/plotting
#' architecture). For a call's plotted window, tests whether the fraction
#' of *non-called* ("background") exons whose observed ratio falls outside
#' the modelled 95% predictive interval is statistically higher than the
#' ~5% a well-calibrated interval implies (one-sided binomial test against
#' a 5% null). A high fraction here doesn't necessarily mean the interval
#' itself is wrong for this sample/region in general — in practice it's
#' more often a sign that the real alteration extends beyond the exons
#' that got called, that the reference set is a poor match for this
#' specific sample/region, or a technical/batch difference between the
#' test sample and its references at this locus. It's a per-call
#' diagnostic flag, not an automatic correction — it doesn't change the
#' call, the interval, or the confidence score.
#'
#' @param ratio Numeric vector of Observed/Expected for every exon in the
#'   plotted window (background and affected together).
#' @param lo,hi Numeric vectors (same length as \code{ratio}) giving the
#'   95% predictive interval bounds at each exon.
#' @param is_affected Logical vector (same length); \code{TRUE} for exons
#'   already called as part of this CNV — excluded from the check, since
#'   those are expected to sit outside the interval.
#' @param min_n Minimum number of background exons required before
#'   flagging (default 5) — below this the percentage is too noisy on its
#'   own to test meaningfully.
#' @return A list with \code{n_background}, \code{n_outside},
#'   \code{pct_outside}, and \code{flag} — \code{flag} is \code{TRUE} when a
#'   one-sided binomial test of \code{n_outside} against a 5% null rate is
#'   significant at p < 0.05.
#' @export
check_background_calibration <- function(ratio, lo, hi, is_affected, min_n = 5) {
    bg <- !is_affected
    n_bg <- sum(bg, na.rm = TRUE)
    if (n_bg == 0) {
        return(list(n_background = 0L, n_outside = 0L, pct_outside = NA_real_, flag = FALSE))
    }
    outside <- (ratio[bg] < lo[bg]) | (ratio[bg] > hi[bg])
    outside[is.na(outside)] <- FALSE
    n_outside <- sum(outside)
    pct_outside <- 100 * n_outside / n_bg
    flag <- n_bg >= min_n &&
        stats::pbinom(n_outside - 1L, size = n_bg, prob = 0.05, lower.tail = FALSE) < 0.05
    list(n_background = n_bg, n_outside = n_outside, pct_outside = pct_outside, flag = flag)
}

# =============================================================================
# Pre-calling sample/exon QC exclusion -- ported from CANOPE for feature
# parity. CANOPE's run_canope() has always excluded noisy samples/exons
# before HMM calling via sample_qc/exon_qc; ECHO had no equivalent, so every
# sample and exon -- however noisy -- went straight into ExomeDepth. These
# two functions are generic (they operate on a plain numeric matrix, not on
# any HMM-specific state), so they're reused here verbatim; only the wiring
# in run_cnv_calling() below is ECHO-specific.
# =============================================================================

#' Detect Outlier Samples (robust z-score on MAD noise)
#'
#' @param counts      Numeric matrix or data frame (targets x samples).
#' @param pseudocount Numeric.
#' @param z_threshold Numeric. Robust z-score above which a sample is an outlier.
#'
#' @return Data frame: sample, noise_score, robust_z, is_outlier.
#' @export
detect_outlier_samples <- function(counts, pseudocount = 0.5, z_threshold = 3) {
    counts <- as.matrix(counts)
    log_counts   <- log2(counts + pseudocount)
    sample_noise <- apply(log_counts, 2, stats::mad, na.rm = TRUE)

    med  <- median(sample_noise, na.rm = TRUE)
    mad0 <- stats::mad(sample_noise, na.rm = TRUE)

    # Guard against mad0 == 0 (all samples identical noise) or a single sample
    z_scores  <- if (length(sample_noise) > 1L && is.finite(mad0) && mad0 > 0)
        (sample_noise - med) / mad0 else rep(0, length(sample_noise))

    is_outlier <- abs(z_scores) > z_threshold

    data.frame(
        sample       = names(sample_noise),
        noise_score  = sample_noise,
        robust_z     = z_scores,
        is_outlier   = is_outlier,
        stringsAsFactors = FALSE
    )
}

#' Detect Problematic Exons
#'
#' Flags exons with high MAD (cross-sample noise), low mean coverage,
#' extreme GC content, or non-finite values. Ensures at least one exon per
#' chromosome is retained to prevent empty chromosomes.
#'
#' @param count_matrix Numeric matrix (targets x samples).
#' @param chromosomes  Optional factor/character vector of chromosome labels.
#' @param mad_quantile Numeric. Top quantile of exon MAD to flag.
#' @param min_mean     Numeric. Minimum mean coverage to retain.
#' @param gc           Optional numeric vector of GC content (0-1 or 0-100).
#' @param gc_min       Numeric. Minimum GC fraction.
#' @param gc_max       Numeric. Maximum GC fraction.
#'
#' @return Data frame: exon (row index), mean, mad, problematic.
#' @export
detect_problematic_exons <- function(
    count_matrix,
    chromosomes  = NULL,
    mad_quantile = 0.90,
    min_mean     = 20,
    gc           = NULL,
    gc_min       = 0.10,
    gc_max       = 0.90
) {
    count_matrix <- as.matrix(count_matrix)
    exon_mean    <- rowMeans(count_matrix, na.rm = TRUE)
    exon_mad     <- apply(count_matrix, 1, stats::mad, na.rm = TRUE)
    mad_thresh   <- stats::quantile(exon_mad, probs = mad_quantile, na.rm = TRUE)

    problematic  <- (exon_mad > mad_thresh | exon_mean < min_mean |
                       !is.finite(exon_mean))

    if (!is.null(gc)) {
        if (length(gc) != nrow(count_matrix))
            stop("'gc' must have one value per row of count_matrix (", nrow(count_matrix),
                 " rows, got ", length(gc), ").")
        gc_val       <- if (max(gc, na.rm = TRUE) > 1) gc / 100 else gc
        problematic  <- problematic | gc_val < gc_min | gc_val > gc_max |
            !is.finite(gc_val)
    }

    # Guarantee at least one exon per chromosome
    if (!is.null(chromosomes)) {
        for (chr in unique(chromosomes)) {
            idx <- which(chromosomes == chr)
            if (length(idx) > 0 && all(problematic[idx]))
                problematic[idx[which.min(exon_mad[idx])]] <- FALSE
        }
    }

    data.frame(
        exon        = seq_len(nrow(count_matrix)),
        mean        = exon_mean,
        mad         = exon_mad,
        problematic = problematic
    )
}
