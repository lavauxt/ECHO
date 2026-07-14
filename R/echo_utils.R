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

#' Reclassify off-target / filler BED intervals before exon numbering
#'
#' Some target panels include intervals that were never a real gene exon at
#' all: normalization/backbone probes placed off-target for coverage
#' calibration, commonly named things like \code{"HorsROI"} ("hors ROI" is
#' French for "outside the region of interest"), \code{"OffTarget"},
#' \code{"Backbone"}, and so on. Left alone, the BED-name parser extracts
#' whatever the first token of that name happens to be (e.g.
#' \code{"HorsROI"}) as if it were a gene symbol, and
#' \code{\link{assign_exon_numbers_per_gene}} then numbers it 1..n exactly
#' like a real gene with its own exons -- so a plot window that happens to
#' straddle one of these intervals shows it interleaved with the real
#' gene's exons under its own (fake) "gene" tile. This function catches
#' those rows, by name, \strong{before} any numbering happens, and lets the
#' caller choose what should happen to them.
#'
#' @param df data.frame with chromosome/Chr, start/Start, end/End,
#'   gene/Gene columns (1-based coordinates; genomic order not required --
#'   \code{handling = "merge"} sorts internally).
#' @param pattern Character. A regular expression (matched against the
#'   \code{gene} column, case-sensitively, consistent with the exon-name
#'   parser) identifying off-target/filler rows -- e.g. \code{"^HorsROI"},
#'   or \code{"^(HorsROI|OffTarget|Backbone)$"} for a panel using several
#'   such labels. \code{NULL} or \code{""} disables this feature entirely
#'   (returns \code{df} unchanged). Default \code{"^HorsROI"}.
#' @param handling One of:
#'   \itemize{
#'     \item \code{"na"} (default) -- keep the interval (it still gets
#'       coverage extracted and still contributes background signal) but
#'       set \code{gene} to \code{NA} on those rows, so
#'       \code{assign_exon_numbers_per_gene()} and everything downstream
#'       (plots, VCF gene column, confidence scoring) leaves them
#'       un-numbered and out of any gene-based grouping, instead of
#'       numbering the filler label as if it were a gene.
#'     \item \code{"remove"} -- drop those rows entirely.
#'     \item \code{"merge"} -- reassign \code{gene} to the nearest
#'       neighbouring \emph{real} gene on the same chromosome (ties go to
#'       the preceding gene), so the interval is treated as one of that
#'       gene's own exons and gets numbered like any other exon. A row
#'       with no real gene anywhere on its chromosome is left unchanged.
#'   }
#' @param verbose Logical. Print a one-line summary. Default \code{TRUE}.
#' @return The same data.frame: column set and row order preserved for
#'   \code{"na"}/\code{"merge"}; row count reduced for \code{"remove"}.
#' @export
handle_off_target_regions <- function(df, pattern = "^HorsROI",
                                      handling = c("na", "remove", "merge"),
                                      verbose = TRUE) {
    handling <- match.arg(handling)
    if (is.null(pattern) || !nzchar(pattern)) return(df)
    if (is.null(df) || nrow(df) == 0) return(df)

    chrom_col <- if ("chromosome" %in% names(df)) "chromosome" else "Chr"
    start_col <- if ("start"      %in% names(df)) "start"      else "Start"
    end_col   <- if ("end"        %in% names(df)) "end"        else "End"
    gene_col  <- if ("gene"       %in% names(df)) "gene"       else "Gene"
    if (!gene_col %in% names(df)) return(df)

    gene_vals <- as.character(df[[gene_col]])
    is_off <- !is.na(gene_vals) & grepl(pattern, gene_vals, perl = TRUE)
    n_off <- sum(is_off)
    if (n_off == 0) return(df)

    if (handling == "remove") {
        out <- df[!is_off, , drop = FALSE]
        if (verbose) {
            message(sprintf("[INFO] handle_off_target_regions: removed %d off-target interval(s) matching /%s/.",
                            n_off, pattern))
        }
        return(out)
    }

    if (handling == "na") {
        df[[gene_col]][is_off] <- NA_character_
        if (verbose) {
            message(sprintf(
                "[INFO] handle_off_target_regions: set gene = NA for %d off-target interval(s) matching /%s/ (kept, excluded from exon numbering).",
                n_off, pattern))
        }
        return(df)
    }

    # handling == "merge": walk the off-target rows in genomic order and
    # attach each one to whichever real gene -- the previous one or the
    # next one, on the same chromosome -- sits closer. A simple forward/
    # backward carry-forward pass (O(n), two linear scans) rather than a
    # per-row search, since a panel BED can run into the tens of
    # thousands of rows.
    ord     <- order(df[[chrom_col]], df[[start_col]], df[[end_col]])
    n       <- length(ord)
    chrom_s <- as.character(df[[chrom_col]])[ord]
    start_s <- df[[start_col]][ord]
    end_s   <- df[[end_col]][ord]
    gene_s  <- gene_vals[ord]
    off_s   <- is_off[ord]

    prev_gene <- character(n); prev_end   <- numeric(n)
    g <- NA_character_; e <- NA_real_; last_chrom <- NA_character_
    for (i in seq_len(n)) {
        if (!identical(chrom_s[i], last_chrom)) { g <- NA_character_; e <- NA_real_; last_chrom <- chrom_s[i] }
        prev_gene[i] <- g; prev_end[i] <- e
        if (!off_s[i]) { g <- gene_s[i]; e <- end_s[i] }
    }

    next_gene <- character(n); next_start <- numeric(n)
    g <- NA_character_; s <- NA_real_; last_chrom <- NA_character_
    for (i in rev(seq_len(n))) {
        if (!identical(chrom_s[i], last_chrom)) { g <- NA_character_; s <- NA_real_; last_chrom <- chrom_s[i] }
        next_gene[i] <- g; next_start[i] <- s
        if (!off_s[i]) { g <- gene_s[i]; s <- start_s[i] }
    }

    new_gene_s <- gene_s
    n_merged   <- 0L
    for (i in which(off_s)) {
        has_prev <- !is.na(prev_gene[i])
        has_next <- !is.na(next_gene[i])
        if (has_prev && has_next) {
            d_prev <- start_s[i] - prev_end[i]
            d_next <- next_start[i] - end_s[i]
            new_gene_s[i] <- if (d_prev <= d_next) prev_gene[i] else next_gene[i]
            n_merged <- n_merged + 1L
        } else if (has_prev) {
            new_gene_s[i] <- prev_gene[i]; n_merged <- n_merged + 1L
        } else if (has_next) {
            new_gene_s[i] <- next_gene[i]; n_merged <- n_merged + 1L
        } # else: no real gene anywhere on this chromosome -- leave as-is
    }

    df[[gene_col]][ord] <- new_gene_s
    if (verbose) {
        message(sprintf(
            "[INFO] handle_off_target_regions: merged %d/%d off-target interval(s) matching /%s/ into their nearest neighbouring gene (%d had no real gene on their chromosome to attach to).",
            n_merged, n_off, pattern, n_off - n_merged))
    }
    df
}

#' Assign sequential exon numbers within each gene
#'
#' Sorts rows by chromosome, start, end and then within each gene assigns
#' 1..n based on genomic order. Duplicate intervals (same chrom, start, end,
#' gene) share the \code{exon_number} of their first occurrence but are
#' \strong{not} dropped from the output.
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

#' Pad the outer edge of each gene's first and last exon
#'
#' Capture-based coverage often drops off right at the true edge of a
#' target interval (probe/bait tiling is rarely perfect exactly at the
#' boundary, and reads whose alignment barely spans the edge get soft-
#' clipped or excluded). For an internal exon this is usually harmless --
#' its neighbours carry the signal -- but for a gene's *first* (lowest-
#' coordinate) or *last* (highest-coordinate) exon there is no such
#' neighbour on the outward side, so a thin sliver of low/zero coverage
#' right at that edge can pull the whole exon's count down and trigger a
#' spurious QC flag or CNV call. This function extends only the outward-
#' facing edge of those two terminal exons per gene (both edges, for a
#' single-exon gene) by \code{padding} bases, giving that boundary the same
#' kind of margin an internal exon already has "for free" from its
#' neighbour.
#'
#' "First"/"last" follows the same purely-genomic (chrom, start, end)
#' ordering that \code{\link{assign_exon_numbers_per_gene}} uses everywhere
#' else in ECHO -- i.e. lowest-coordinate exon is #1 regardless of strand.
#' Internal exons (and the inward-facing edge of a terminal exon) are left
#' untouched: padding an internal boundary would just eat into the intron
#' between two already-covered exons for no benefit.
#'
#' Padding is applied on a best-effort basis ("if possible"): the function
#' never creates an overlap with whatever interval sits next to it on the
#' same chromosome (a neighbouring exon of the same gene or of a different
#' one), and never pushes a coordinate below 1 or past the contig length
#' (when \code{chr_lengths} is supplied). Where the available gap is
#' narrower than \code{padding}, that side is extended only as far as the
#' gap allows.
#'
#' Rows whose \code{gene} value is missing/empty or the literal placeholder
#' \code{"Unknown"} are left unpadded, since "first/last exon of a gene"
#' isn't meaningful for unannotated regions.
#'
#' @param bed_file data.frame with columns chromosome/Chr, start/Start,
#'   end/End, gene/Gene (1-based, inclusive coordinates).
#' @param padding Integer >= 0. Bases to add to the outward edge of each
#'   gene's first and last exon. \code{0} (the default) disables padding
#'   and returns \code{bed_file} unchanged.
#' @param chr_lengths Optional named numeric vector (names = chromosome,
#'   values = contig length) used to cap the last exon's End at the contig
#'   boundary. If \code{NULL}, no contig-length clamp is applied (only the
#'   neighbouring-interval clamp).
#' @param verbose Logical. Print a one-line summary of how many
#'   starts/ends were extended and how many sides were clamped short of
#'   the requested padding. Default \code{TRUE}.
#' @return The same data.frame (original row order and row count
#'   preserved), with Start/End adjusted for terminal-exon rows only.
#' @export
pad_gene_terminal_exons <- function(bed_file, padding = 0, chr_lengths = NULL, verbose = TRUE) {
    if (is.null(padding) || length(padding) != 1 || is.na(padding) || padding <= 0) {
        return(bed_file)
    }
    padding <- as.integer(round(padding))

    chrom_col <- if ("chromosome" %in% names(bed_file)) "chromosome" else "Chr"
    start_col <- if ("start"      %in% names(bed_file)) "start"      else "Start"
    end_col   <- if ("end"        %in% names(bed_file)) "end"        else "End"
    gene_col  <- if ("gene"       %in% names(bed_file)) "gene"       else "Gene"
    stopifnot(all(c(chrom_col, start_col, end_col, gene_col) %in% names(bed_file)))

    n_in <- nrow(bed_file)
    if (n_in == 0) return(bed_file)

    # Reuse the pipeline's own per-gene ordering so "first"/"last" here
    # always agrees with exon_number everywhere else in ECHO.
    numbered <- assign_exon_numbers_per_gene(bed_file)

    dt <- data.table::as.data.table(numbered)
    dt[, .orig_row := .I]  # remember incoming row order
    data.table::setnames(dt, c(chrom_col, start_col, end_col, gene_col),
                         c("chrom", "start", "end", "gene"))

    dt[, is_first := exon_number == 1L]
    dt[, is_last  := exon_number == max(exon_number), by = "gene"]
    no_gene <- is.na(dt$gene) | dt$gene %in% c("", ".", "Unknown")
    dt[no_gene, c("is_first", "is_last") := FALSE]

    # Sort a copy by genomic position (per chromosome) so each terminal
    # exon can see its nearest neighbour on either side -- regardless of
    # which gene that neighbour belongs to -- and never be padded into it.
    chrom_levels <- c(paste0("chr", c(1:22, "X", "Y", "M")),
                      c(as.character(1:22), "X", "Y", "M"))
    dt[, .chrom_fac := factor(chrom, levels = unique(c(chrom_levels, unique(chrom))))]
    data.table::setorder(dt, .chrom_fac, start, end)

    n        <- nrow(dt)
    chrom_id <- as.integer(dt$.chrom_fac)
    start_v  <- dt$start
    end_v    <- dt$end
    is_first_v <- dt$is_first
    is_last_v  <- dt$is_last

    right_extend <- integer(n)  # applies to is_last rows: bp added to end
    left_extend  <- integer(n)  # applies to is_first rows: bp subtracted from start

    # Gap k (k = 1..n-1) sits between sorted row k and row k+1. Both may
    # want a share of it at once -- row k if it's a last exon growing
    # rightward, row k+1 if it's a first exon growing leftward (this is
    # the one place two *different* genes' terminal exons can compete for
    # the same free space). Give each what it asks for if the gap is big
    # enough for both; otherwise split the gap between them so neither
    # padded interval ever crosses into the other's.
    if (n > 1) {
        same_chr_pair <- chrom_id[-n] == chrom_id[-1]
        gap        <- pmax(start_v[-1] - end_v[-n] - 1L, 0L)
        want_left  <- ifelse(is_last_v[-n],  padding, 0L)  # row k wants to grow right
        want_right <- ifelse(is_first_v[-1], padding, 0L)  # row k+1 wants to grow left
        demand     <- want_left + want_right

        grant_left  <- integer(n - 1L)
        grant_right <- integer(n - 1L)
        has_demand  <- same_chr_pair & demand > 0
        fits        <- has_demand & demand <= gap
        tight       <- has_demand & demand > gap

        grant_left[fits]  <- want_left[fits]
        grant_right[fits] <- want_right[fits]
        # Proportional split of a too-small gap, rounded down so the two
        # grants can never sum to more than the gap itself.
        grant_left[tight]  <- as.integer(floor(gap[tight] * want_left[tight] / demand[tight]))
        grant_right[tight] <- gap[tight] - grant_left[tight]

        right_extend[-n] <- grant_left
        left_extend[-1]  <- grant_right
    }

    # Rows at a chromosome boundary (no same-chromosome neighbour on the
    # relevant side) have no interval to compete with there, so they fall
    # back to the contig start (position 1) / contig length instead.
    has_prev <- c(FALSE, if (n > 1) chrom_id[-1] == chrom_id[-n] else logical(0))
    has_next <- c(if (n > 1) chrom_id[-n] == chrom_id[-1] else logical(0), FALSE)

    no_prev_first <- is_first_v & !has_prev
    if (any(no_prev_first)) {
        left_extend[no_prev_first] <- pmin(padding, pmax(start_v[no_prev_first] - 1L, 0L))
    }

    no_next_last <- is_last_v & !has_next
    if (any(no_next_last)) {
        chr_len_here <- if (!is.null(chr_lengths)) {
            unname(chr_lengths[as.character(dt$chrom[no_next_last])])
        } else {
            rep(NA_real_, sum(no_next_last))
        }
        avail <- ifelse(!is.na(chr_len_here), chr_len_here - end_v[no_next_last], Inf)
        right_extend[no_next_last] <- pmin(padding, pmax(avail, 0))
    }

    n_padded_start <- sum(left_extend > 0L)
    n_padded_end   <- sum(right_extend > 0L)
    n_clamped      <- sum(is_first_v & left_extend  < padding) +
                       sum(is_last_v  & right_extend < padding)

    dt[, start := start_v - left_extend]
    dt[, end   := end_v   + right_extend]

    data.table::setorder(dt, .orig_row)  # restore original (input) row order
    dt[, c(".orig_row", ".chrom_fac", "is_first", "is_last", "exon_number") := NULL]
    data.table::setnames(dt, c("chrom", "start", "end", "gene"),
                         c(chrom_col, start_col, end_col, gene_col))
    out <- as.data.frame(dt)
    stopifnot(nrow(out) == n_in)

    if (verbose) {
        message(sprintf(
            "[INFO] pad_gene_terminal_exons: requested %d bp padding -- extended %d gene start(s) and %d gene end(s); %d side(s) received less than the full request (shared gap with a neighbouring interval, or a contig boundary).",
            padding, n_padded_start, n_padded_end, n_clamped))
    }
    out
}

#' Compute gap-inserted x-axis positions for CNV window plots
#'
#' All of ECHO's per-call plots (the four PDF panels in \code{plots.R} and
#' the three interactive panels in \code{ECHO_global_report.Rmd}) lay a
#' window of exons out along a single x-axis. Plotted at plain 1..n integer
#' positions, a gene boundary inside that window looks identical to an
#' ordinary intron between two exons of the *same* gene -- there's nothing
#' to tell a reader "these two points belong to different genes" other than
#' the tile-track colour (PDF only; the HTML report has no tile track at
#' all). This function computes an alternative x-position (\code{px}) for
#' each exon in the window that inserts \code{gap} extra, unlabelled axis
#' units wherever the \code{gene} column changes between consecutive exons
#' -- i.e. between a gene's last exon and the next gene's first exon --
#' while keeping ordinary within-gene spacing at a plain 1 unit. The result
#' is blank visual space at every gene boundary, in every panel, with no
#' change to which exons are shown or how they're labelled.
#'
#' It also returns \code{gene_group}, a per-position integer that increments
#' at every such boundary. Passing this as the \code{group} aesthetic on a
#' \code{geom_ribbon()} keeps that visual gap genuinely blank -- otherwise
#' ggplot draws a single connected ribbon straight across it, right through
#' the empty space \code{px} just created. (Points are never joined by a
#' line in these plots, so this only matters for the ribbon.)
#'
#' @param bed_file data.frame with a gene/Gene column. \code{exon_range}
#'   values are row indices into this data.frame.
#' @param exon_range Integer vector of \code{bed_file} row indices, in the
#'   order they'll be plotted along the x-axis (ascending genomic order,
#'   as produced by \code{prepare_plot_data()}/\code{get_cnv_plot_data()}).
#' @param gap Numeric >= 0. Extra x-axis units inserted at each gene
#'   boundary. \code{0} falls back to plain 1..n spacing (no visual gap,
#'   but \code{gene_group} is still computed correctly). Default \code{1}.
#' @return data.frame with one row per element of \code{exon_range}:
#'   \code{idx} (the original \code{bed_file} row index, i.e. the input
#'   \code{exon_range} value), \code{px} (the x-axis position to plot at),
#'   \code{gene_break} (TRUE at the first exon of a new gene, i.e.
#'   immediately after a gap) and \code{gene_group} (integer, constant
#'   within a gene's run of exons, incrementing at each \code{gene_break}).
#' @export
compute_gene_gap_positions <- function(bed_file, exon_range, gap = 1) {
    if (length(exon_range) == 0) {
        return(data.frame(idx = integer(0), px = numeric(0),
                          gene_break = logical(0), gene_group = integer(0)))
    }
    gap <- if (is.null(gap) || is.na(gap) || gap < 0) 1 else gap

    gene_col <- if ("gene" %in% names(bed_file)) "gene" else "Gene"
    genes    <- as.character(bed_file[[gene_col]][exon_range])
    genes[is.na(genes)] <- ""

    gene_break    <- c(FALSE, genes[-1] != genes[-length(genes)])
    step          <- ifelse(gene_break, 1 + gap, 1)
    step[1]       <- 0
    px            <- cumsum(step) + 1

    data.frame(idx = exon_range, px = px, gene_break = gene_break,
              gene_group = cumsum(gene_break) + 1L, stringsAsFactors = FALSE)
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
