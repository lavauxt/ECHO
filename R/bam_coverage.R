#' Extract BAM Coverage
#'
#' Reads BAM files and computes exon-level read counts using
#' \code{ExomeDepth::getBamCounts}.
#'
#' @param bamfiles Character string. File containing BAM file paths (one per line).
#' @param bamdir Character string. Directory containing BAM files.
#' @param bed Character string. Path to target BED file.
#' @param fasta Character string. Path to reference genome FASTA file.
#' @param rbams Character string. Optional TSV file with a \code{bam} column.
#' @param data_out Character string. Output \code{.RData} file path.
#' @param include.chr Logical. Passed to \code{ExomeDepth::getBamCounts}.
#' @param verbose Logical. Print progress messages.
#' @param sample_name_delim Character string. Delimiter to split filenames.
#' @param sample_name_keep Character string. Which parts to keep.
#' @param custom_sample_names Optional character vector.
#' @param bed_zero_based Logical. If TRUE, BED start coordinates are 0‑based.
#' @param skip_invalid_intervals Logical. If TRUE, intervals that cannot be
#'   extracted from FASTA are silently removed (with a warning). Default TRUE.
#'
#' @return Invisibly \code{NULL}. Coverage object saved to disk.
#' @export
run_bam_coverage <- function(
    bamfiles = NULL,
    bamdir = NULL,
    bed,
    fasta,
    rbams = NULL,
    data_out = "./ECHO_coverage.Rdata",
    include.chr = FALSE,
    verbose = TRUE,
    sample_name_delim = "\\.",
    sample_name_keep = "1",
    custom_sample_names = NULL,
    bed_zero_based = TRUE,
    skip_invalid_intervals = TRUE
) {

    if (verbose) message("[INFO] ", Sys.time(), " BEGIN bam coverage calculation")

    if (is.null(bamfiles) && is.null(bamdir))
        stop("[ERROR] Either bamfiles or bamdir must be provided")

    stop_if_not_file(bed, "[ERROR] BED file not found")
    stop_if_not_file(fasta, "[ERROR] FASTA file not found")

    # -------------------------------------------------------------------------
    # Ensure FASTA index exists
    # -------------------------------------------------------------------------
    if (!file.exists(paste0(fasta, ".fai"))) {
        if (requireNamespace("Rsamtools", quietly = TRUE)) {
            message("[INFO] Creating FASTA index: ", fasta, ".fai")
            Rsamtools::indexFa(fasta)
        } else {
            stop("[ERROR] FASTA index missing and Rsamtools not available. ",
                 "Please run: samtools faidx ", fasta)
        }
    }

    # -------------------------------------------------------------------------
    # Get FASTA chromosome info
    # -------------------------------------------------------------------------
    fa <- Rsamtools::FaFile(fasta)
    seq_info <- Rsamtools::seqinfo(fa)
    fasta_chroms <- trimws(as.character(GenomeInfoDb::seqnames(seq_info)))
    fasta_lengths <- GenomeInfoDb::seqlengths(seq_info)
    names(fasta_lengths) <- fasta_chroms

    if (verbose) {
        message("[INFO] FASTA contains ", length(fasta_chroms), " sequences")
        message("[INFO] First few FASTA seqnames: ",
                paste(head(fasta_chroms, 5), collapse = ", "))
    }

    # -------------------------------------------------------------------------
    # Collect BAM files
    # -------------------------------------------------------------------------
    bams <- character()
    if (!is.null(bamfiles) && file.exists(bamfiles))
        bams <- c(bams, readLines(bamfiles))
    if (!is.null(bamdir) && dir.exists(bamdir))
        bams <- c(bams, list.files(bamdir, pattern = "\\.bam$", full.names = TRUE))
    if (!is.null(rbams) && file.exists(rbams))
        bams <- c(utils::read.csv(rbams, header = TRUE, sep = "\t")$bam, bams)

    bams <- unique(bams[nzchar(bams)])
    if (!length(bams)) stop("[ERROR] No BAM files found")

    missing_idx <- !vapply(bams, bam_has_index, logical(1))
    if (any(missing_idx))
        stop("[ERROR] Missing BAM index for: ", paste(bams[missing_idx], collapse = ", "))

    # -------------------------------------------------------------------------
    # Sample name extraction
    # -------------------------------------------------------------------------
    if (!is.null(custom_sample_names)) {
        if (length(custom_sample_names) != length(bams))
            stop("[ERROR] custom_sample_names length (", length(custom_sample_names),
                 ") does not match number of BAMs (", length(bams), ")")
        sample_names <- custom_sample_names
        if (verbose) message("[INFO] Using custom sample names")
    } else {
        base_names <- tools::file_path_sans_ext(basename(bams))
        keep_parts <- function(x, delim, keep_spec) {
            parts <- strsplit(x, delim)[[1]]
            if (grepl("-", keep_spec)) {
                idx <- as.numeric(strsplit(keep_spec, "-")[[1]])
                parts <- parts[idx[1]:idx[2]]
            } else {
                idx <- as.numeric(keep_spec)
                parts <- parts[idx]
            }
            paste(parts, collapse = delim)
        }
        sample_names <- vapply(base_names, function(fname) {
            keep_parts(fname, sample_name_delim, sample_name_keep)
        }, character(1))
        names(sample_names) <- NULL
        if (verbose) {
            message("[INFO] Extracted sample names using delimiter: ", sample_name_delim,
                    ", keeping part(s): ", sample_name_keep)
            n_show <- min(3, length(sample_names))
            for (i in seq_len(n_show))
                message("  ", bams[i], " -> ", sample_names[i])
            if (length(sample_names) > n_show) message("  ...")
        }
    }

    if (any(duplicated(sample_names))) {
        dup_names <- unique(sample_names[duplicated(sample_names)])
        stop("[ERROR] Duplicate sample names after extraction: ", paste(dup_names, collapse = ", "),
             "\nTry a different delimiter or part specification, or provide custom_sample_names")
    }

    # -------------------------------------------------------------------------
    # Read BED file
    # -------------------------------------------------------------------------
    bed_file <- data.table::fread(bed, header = FALSE)
    colnames(bed_file) <- if (ncol(bed_file) == 5) {
        c("chromosome", "start", "end", "gene", "exon_number")
    } else {
        c("chromosome", "start", "end", "gene")
    }

    # Trim whitespace
    bed_file$chromosome <- trimws(as.character(bed_file$chromosome))

    # -------------------------------------------------------------------------
    # Convert coordinates if BED is 0‑based
    # -------------------------------------------------------------------------
    if (bed_zero_based) {
        bed_file$start <- bed_file$start + 1
        if (verbose) message("[INFO] Converted BED start coordinates from 0‑based to 1‑based")
    }

    # -------------------------------------------------------------------------
    # Harmonise chromosome names with FASTA
    # -------------------------------------------------------------------------
    map_chrom <- function(bed_chr) {
        if (bed_chr %in% fasta_chroms) return(bed_chr)
        with_chr <- ifelse(grepl("^chr", bed_chr), bed_chr, paste0("chr", bed_chr))
        if (with_chr %in% fasta_chroms) return(with_chr)
        without_chr <- sub("^chr", "", bed_chr)
        if (without_chr %in% fasta_chroms) return(without_chr)
        return(NA_character_)
    }

    bed_file$original_chrom <- bed_file$chromosome
    bed_file$chromosome <- vapply(bed_file$chromosome, map_chrom, character(1))

    # Remove unmappable rows
    removed_rows <- is.na(bed_file$chromosome)
    if (any(removed_rows)) {
        removed <- unique(bed_file$original_chrom[removed_rows])
        if (verbose) message("[WARNING] Removed ", sum(removed_rows),
                             " rows with chromosomes not in FASTA: ", paste(removed, collapse = ", "))
        bed_file <- bed_file[!removed_rows, ]
    }

    # Validate numeric coordinates
    bed_file$start <- as.numeric(bed_file$start)
    bed_file$end <- as.numeric(bed_file$end)
    valid_coords <- which(bed_file$start > 0 & bed_file$end >= bed_file$start)
    if (length(valid_coords) < nrow(bed_file)) {
        message("[WARNING] Removing ", nrow(bed_file) - length(valid_coords), " rows with invalid coordinates")
        bed_file <- bed_file[valid_coords, ]
    }

    # Remove out‑of‑bound intervals
    chr_len <- fasta_lengths[bed_file$chromosome]
    out_of_bounds <- which(bed_file$end > chr_len | bed_file$start < 1)
    if (length(out_of_bounds) > 0) {
        message("[WARNING] Removing ", length(out_of_bounds), " rows outside chromosome boundaries")
        bed_file <- bed_file[-out_of_bounds, ]
    }

    # -------------------------------------------------------------------------
    # Pre‑validate each interval with Rsamtools::scanFa (skip problematic ones)
    # -------------------------------------------------------------------------
    if (skip_invalid_intervals && nrow(bed_file) > 0) {
        if (verbose) message("[INFO] Testing FASTA extraction for ", nrow(bed_file), " intervals...")
        keep <- logical(nrow(bed_file))
        for (i in seq_len(nrow(bed_file))) {
            gr <- GenomicRanges::GRanges(
                seqnames = bed_file$chromosome[i],
                ranges = IRanges::IRanges(start = bed_file$start[i], end = bed_file$end[i])
            )
            ok <- tryCatch({
                seq <- Rsamtools::scanFa(fa, param = gr)
                length(seq) == 1 && BiocGenerics::width(seq)[1] > 0
            }, error = function(e) FALSE)
            keep[i] <- ok
        }
        if (!all(keep)) {
            removed <- which(!keep)
            message("[WARNING] Removing ", length(removed), " intervals that cannot be extracted from FASTA.")
            if (verbose) {
                for (r in utils::head(removed, 10)) {
                    message("   Row ", r, ": ", bed_file$chromosome[r], ":",
                            bed_file$start[r], "-", bed_file$end[r])
                }
                if (length(removed) > 10) message("   ... and ", length(removed)-10, " more.")
            }
            bed_file <- bed_file[keep, ]
        }
    }

    # Sort BED
    chr_levels <- fasta_chroms[fasta_chroms %in% unique(bed_file$chromosome)]
    bed_file$chromosome <- factor(bed_file$chromosome, levels = chr_levels)
    data.table::setorder(bed_file, chromosome, start, end)
    bed_file$chromosome <- as.character(bed_file$chromosome)
    bed_file$original_chrom <- NULL

    if (verbose) message("[INFO] BED file cleaned. Final rows: ", nrow(bed_file))
    if (nrow(bed_file) == 0) stop("[ERROR] No valid intervals remain.")

    # -------------------------------------------------------------------------
    # Compute counts
    # -------------------------------------------------------------------------
    counts <- tryCatch({
        ExomeDepth::getBamCounts(
            bed.frame = as.data.frame(bed_file),
            bam.files = bams,
            include.chr = include.chr,
            referenceFasta = fasta
        )
    }, error = function(e) {
        msg <- conditionMessage(e)
        row_num <- suppressWarnings(as.numeric(gsub(".*record ([0-9]+).*", "\\1", msg)))
        if (!is.na(row_num) && row_num <= nrow(bed_file)) {
            stop("[ERROR] Interval extraction failed for row ", row_num, ": ",
                 bed_file$chromosome[row_num], ":", bed_file$start[row_num], "-", bed_file$end[row_num],
                 "\nThis interval was not caught by pre‑validation. Check your FASTA.\n",
                 "Original error: ", msg)
        } else {
            stop(e)
        }
    })

    # -------------------------------------------------------------------------
    # Sample column mapping
    # -------------------------------------------------------------------------
    normalised_bams <- vapply(bams, function(b) make.names(basename(b)), character(1))
    sample_cols <- intersect(colnames(counts), normalised_bams)
    if (length(sample_cols) != length(bams)) {
        sample_cols <- intersect(colnames(counts), paste0("X", normalised_bams))
    }
    if (length(sample_cols) != length(bams)) {
        metadata_cols <- c("chromosome", "start", "end", "exon", "GC", "exon_number")
        sample_cols <- setdiff(colnames(counts), metadata_cols)
        sample_cols <- sample_cols[order(match(normalised_bams, sample_cols))]
        if (verbose) message("[INFO] Using alternative sample column detection and reordering")
    }
    name_mapping <- setNames(sample_names, sample_cols)
    for (old_name in sample_cols) {
        new_name <- name_mapping[old_name]
        colnames(counts)[colnames(counts) == old_name] <- new_name
        if (verbose) message("[INFO] Renamed column: ", old_name, " -> ", new_name)
    }

    keep_cols <- c("chromosome", "start", "end", "exon", sample_names)
    if ("exon_number" %in% colnames(counts)) keep_cols <- c(keep_cols, "exon_number")
    counts <- counts[, keep_cols, drop = FALSE]

    count_sample_cols <- setdiff(colnames(counts), c("chromosome", "start", "end", "exon", "exon_number"))
    if (!identical(sort(sample_names), sort(count_sample_cols))) {
        stop("[ERROR] Mismatch between sample_names list and coverage data frame columns.\n",
             "sample_names: ", paste(sample_names, collapse=", "), "\n",
             "coverage columns: ", paste(count_sample_cols, collapse=", "))
    }
    if (verbose) message("[INFO] Verified: sample_names match coverage columns.")

    if ("exon_number" %in% colnames(bed_file)) {
        exon_num <- bed_file[["exon_number"]]
        counts <- cbind(counts, exon_number = exon_num)
        exon_pos <- which(names(counts) == "exon")
        if (length(exon_pos) == 1) {
            counts <- counts[, c(names(counts)[1:exon_pos], "exon_number",
                                 names(counts)[(exon_pos+1):(ncol(counts)-1)])]
        }
    }

    dup_mask <- duplicated(counts[, c("chromosome", "start", "end", "exon")])
    if (any(dup_mask)) {
        message("[WARNING] Removed ", sum(dup_mask), " duplicate target region(s)")
        counts <- counts[!dup_mask, , drop = FALSE]
        bed_file <- bed_file[!dup_mask, ]
    }

    if (nrow(bed_file) != nrow(counts))
        stop("[ERROR] bed_file and counts row counts do not match (", nrow(bed_file), " vs ", nrow(counts), ")")

    bed_file <- as.data.frame(bed_file)
    counts <- as.data.frame(counts)

    save(counts, bams, bed_file, sample_names, fasta, name_mapping, file = data_out)
    if (verbose) message("[INFO] ", Sys.time(), " END bam coverage calculation")
    invisible(NULL)
}
