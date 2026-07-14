#' Resolve the reference FASTA used for coverage/GC-content calculation
#'
#' Supports three sources:
#' \itemize{
#'   \item \code{fasta_source = "file"} – use a user-supplied FASTA file (original behaviour).
#'   \item \code{fasta_source = "bsgenome"} – use a BSgenome package determined by \code{genome_version}
#'         (e.g., "hg19" → \code{BSgenome.Hsapiens.UCSC.hg19}).
#'   \item \code{fasta_source = "BSgenome.Hsapiens.UCSC.hg19"} – directly specify the installed BSgenome
#'         package name (no \code{genome_version} mapping required).
#' }
#'
#' @param fasta Character string or NULL. Path to a FASTA file (required when
#'   \code{fasta_source = "file"}).
#' @param fasta_source One of \code{"file"}, \code{"bsgenome"}, or a character string
#'   naming a Bioconductor BSgenome package (e.g. \code{"BSgenome.Hsapiens.UCSC.hg19"}).
#' @param genome_version Character string, e.g. \code{"hg19"} or \code{"hg38"}. Required
#'   only when \code{fasta_source = "bsgenome"} (ignored if a package name is given directly).
#' @param bed Character string or NULL. Path to the BED file; used only
#'   to determine which chromosomes are needed (for BSgenome sources).
#' @param cache_dir Character string or NULL. Where to cache the
#'   BSgenome-derived FASTA. Defaults to a per-user cache directory.
#' @param verbose Logical. Print progress messages.
#' @return Character string: path to a FASTA file, ready to pass to
#'   \code{Rsamtools::FaFile()} / \code{ExomeDepth::getBamCounts()}.
#' @export
resolve_reference_fasta <- function(fasta = NULL,
                                     fasta_source = c("file", "bsgenome"),
                                     genome_version = NULL,
                                     bed = NULL,
                                     cache_dir = NULL,
                                     verbose = TRUE) {
  # If fasta_source is a character string that is not "file" or "bsgenome",
  # treat it as a BSgenome package name.
  if (is.character(fasta_source) && length(fasta_source) == 1 &&
      !fasta_source %in% c("file", "bsgenome")) {
    pkg <- fasta_source
    fasta_source <- "bsgenome_pkg"   # internal marker
  } else {
    fasta_source <- match.arg(fasta_source)
  }

  if (fasta_source == "file") {
    stop_if_not_file(fasta, "[ERROR] FASTA file not found")
    return(fasta)
  }

  # ---- fasta_source is either "bsgenome" (use genome_version mapping) or "bsgenome_pkg" (direct) ----
  if (fasta_source == "bsgenome_pkg") {
    # pkg already assigned
  } else { # fasta_source == "bsgenome"
    if (is.null(genome_version)) {
      stop("[ERROR] genome_version is required when fasta_source = 'bsgenome'")
    }
    pkg <- switch(tolower(genome_version),
      hg19   = "BSgenome.Hsapiens.UCSC.hg19",
      grch37 = "BSgenome.Hsapiens.UCSC.hg19",
      hg38   = "BSgenome.Hsapiens.UCSC.hg38",
      grch38 = "BSgenome.Hsapiens.UCSC.hg38",
      stop("[ERROR] No BSgenome package mapping for genome_version '", genome_version,
           "'. Supported: hg19, hg38. Use fasta_source = 'file' for other builds.")
    )
  }

  # ---- Validate the package exists ----
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("[ERROR] Package '", pkg, "' is required for fasta_source = '", 
         if (fasta_source == "bsgenome") "bsgenome" else "bsgenome_pkg", 
         "' but is not installed.\n",
         "  Install with: BiocManager::install(\"", pkg, "\")\n",
         "  Or use fasta_source = 'file' with your own FASTA.")
  }

  bsgenome        <- getExportedValue(pkg, pkg)
  genome_seqnames <- GenomeInfoDb::seqnames(bsgenome)
  # Which chromosomes does the BED actually need? (cheap first-column peek)
  chroms_needed <- NULL
  if (!is.null(bed) && file.exists(bed)) {
    bed_chr_col <- tryCatch(data.table::fread(bed, header = FALSE, select = 1L)[[1]],
                             error = function(e) NULL)
    if (!is.null(bed_chr_col)) chroms_needed <- unique(trimws(as.character(bed_chr_col)))
  }

  map_one <- function(x) {
    if (x %in% genome_seqnames) return(x)
    alt <- if (grepl("^chr", x)) sub("^chr", "", x) else paste0("chr", x)
    if (alt %in% genome_seqnames) return(alt)
    NA_character_
  }
  if (!is.null(chroms_needed)) {
    mapped <- unique(stats::na.omit(vapply(chroms_needed, map_one, character(1))))
    if (!length(mapped)) {
      message("[WARNING] None of the BED chromosome names matched ", pkg,
              "; falling back to the full genome (this will be slower and larger).")
      mapped <- genome_seqnames
    }
  } else {
    mapped <- genome_seqnames
  }

  if (is.null(cache_dir)) {
    cache_dir <- tryCatch(tools::R_user_dir("ECHO", which = "cache"), error = function(e) tempdir())
  }
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  cache_fasta <- file.path(cache_dir, paste0("ECHO_ref_", genome_version, "_bsgenome.fa"))
  cache_fai   <- paste0(cache_fasta, ".fai")

  existing_chroms <- character(0)
  if (file.exists(cache_fasta) && file.exists(cache_fai)) {
    existing_chroms <- tryCatch(data.table::fread(cache_fai, header = FALSE, select = 1L)[[1]],
                                 error = function(e) character(0))
  }
  missing_chroms <- setdiff(mapped, existing_chroms)

  if (length(missing_chroms) > 0) {
    if (verbose) message("[INFO] Building reference FASTA from ", pkg, " for ", length(missing_chroms),
                          " sequence(s): ", paste(missing_chroms, collapse = ", "))
    new_seqs <- Biostrings::getSeq(bsgenome, names = missing_chroms)
    names(new_seqs) <- missing_chroms
    # Rebuilding the index from scratch immediately after every write means
    # the .fai can never end up stale relative to the file's own contents.
    Biostrings::writeXStringSet(new_seqs, filepath = cache_fasta, append = file.exists(cache_fasta))
    if (file.exists(cache_fai)) file.remove(cache_fai)
    Rsamtools::indexFa(cache_fasta)
  } else if (verbose) {
    message("[INFO] Reusing cached BSgenome-derived reference FASTA: ", cache_fasta)
  }

  cache_fasta
}

#' Extract BAM Coverage
#'
#' Reads BAM files and computes exon-level read counts using
#' \code{ExomeDepth::getBamCounts}.
#'
#' @param bamfiles Character string. File containing BAM file paths (one per line).
#' @param bamdir Character string. Directory containing BAM files.
#' @param bed Character string. Path to target BED file.
#' @param fasta Character string. Path to reference genome FASTA file.
#'   Required when \code{fasta_source = "file"} (the default); ignored
#'   when \code{fasta_source = "bsgenome"}.
#' @param fasta_source One of \code{"file"} (use \code{fasta} directly,
#'   the historical behaviour) or \code{"bsgenome"} (fetch sequence from
#'   an installed Bioconductor BSgenome package for \code{genome_version},
#'   no local FASTA file needed). Default \code{"file"}.
#' @param genome_version Character string, e.g. \code{"hg19"} or
#'   \code{"hg38"}. Required when \code{fasta_source = "bsgenome"}.
#' @param bsgenome_cache_dir Character string or NULL. Directory used to
#'   cache the BSgenome-derived FASTA across runs. Defaults to a per-user
#'   cache directory. Ignored when \code{fasta_source = "file"}.
#' @param rbams Character string. Optional TSV file with a \code{bam} column.
#' @param data_out Character string. Output \code{.RData} file path.
#' @param include.chr Logical. Passed to \code{ExomeDepth::getBamCounts}.
#' @param verbose Logical. Print progress messages.
#' @param sample_name_delim Character string. Regex delimiter to split filenames.
#' @param sample_name_keep Character string. Which parts to keep (e.g., "1", "1-2").
#' @param sample_name_collapse Character string. Separator for rejoining parts. Default NULL.
#' @param custom_sample_names Optional character vector.
#' @param bed_zero_based Logical. If TRUE, BED start coordinates are 0-based.
#' @param skip_invalid_intervals Logical. If TRUE, validate each interval against FASTA.
#' @param pad_terminal_exons Integer >= 0. Bases to extend the outward-facing
#'   edge of each gene's first and last exon (both edges, for a single-exon
#'   gene), to reduce the chance of a 0/low count right at a gene boundary
#'   where there's no neighbouring exon to carry the signal. Never creates
#'   an overlap with a neighbouring interval or crosses a contig boundary --
#'   padding is applied "if possible" and clamped short otherwise. See
#'   \code{\link{pad_gene_terminal_exons}}. Default \code{0} (disabled).
#' @return Invisibly \code{NULL}. Coverage object saved to disk.
#' @export
run_bam_coverage <- function(
  bamfiles = NULL,
  bamdir   = NULL,
  bed,
  fasta    = NULL,
  fasta_source       = "file",
  genome_version     = NULL,
  bsgenome_cache_dir = NULL,
  rbams    = NULL,
  data_out = "./ECHO_coverage.Rdata",
  include.chr = FALSE,
  verbose  = TRUE,
  sample_name_delim    = "\\.",
  sample_name_keep     = "1",
  sample_name_collapse = NULL,
  custom_sample_names  = NULL,
  bed_zero_based       = TRUE,
  skip_invalid_intervals = FALSE,
  pad_terminal_exons   = 0
) {
  if (verbose) message("[INFO] ", Sys.time(), " BEGIN bam coverage calculation")
  if (is.null(bamfiles) && is.null(bamdir)) stop("[ERROR] Either bamfiles or bamdir must be provided")
  stop_if_not_file(bed, "[ERROR] BED file not found")

  fasta <- resolve_reference_fasta(
    fasta          = fasta,
    fasta_source   = fasta_source,
    genome_version = genome_version,
    bed            = bed,
    cache_dir      = bsgenome_cache_dir,
    verbose        = verbose
  )

  fai_path    <- paste0(fasta, ".fai")
  needs_index <- !file.exists(fai_path)
  if (!needs_index && file.info(fasta)$mtime > file.info(fai_path)$mtime) {
    # A .fai older than the FASTA it indexes is a classic, silent cause of
    # "record N failed" errors deep inside Rsamtools/ExomeDepth: the byte
    # offsets it stores no longer match the file's actual contents (e.g.
    # the FASTA was edited, re-wrapped, or its line endings changed after
    # the index was built). Rebuild rather than trust a stale index.
    if (verbose) message("[WARNING] FASTA index (.fai) is older than the FASTA file itself - ",
                          "it may be stale. Rebuilding index.")
    needs_index <- TRUE
  }
  if (needs_index) {
    if (requireNamespace("Rsamtools", quietly = TRUE)) {
      message("[INFO] Creating FASTA index: ", fai_path)
      if (file.exists(fai_path)) file.remove(fai_path)
      Rsamtools::indexFa(fasta)
    } else {
      stop("[ERROR] FASTA index missing and Rsamtools not available. Run: samtools faidx ", fasta)
    }
  }

  fa           <- Rsamtools::FaFile(fasta)
  seq_info     <- Rsamtools::seqinfo(fa)
  fasta_chroms <- trimws(as.character(GenomeInfoDb::seqnames(seq_info)))
  fasta_lengths <- GenomeInfoDb::seqlengths(seq_info)
  names(fasta_lengths) <- fasta_chroms

  if (verbose) message("[INFO] FASTA contains ", length(fasta_chroms),
                        " sequence(s) (source: ", fasta_source, ")")

  bams <- character()
  if (!is.null(bamfiles) && file.exists(bamfiles)) bams <- c(bams, readLines(bamfiles))
  if (!is.null(bamdir)   && dir.exists(bamdir))    bams <- c(bams, list.files(bamdir, pattern = "\\.bam$", full.names = TRUE))
  if (!is.null(rbams)    && file.exists(rbams))    bams <- c(utils::read.csv(rbams, header = TRUE, sep = "\t")$bam, bams)
  bams <- unique(bams[nzchar(bams)])
  if (!length(bams)) stop("[ERROR] No BAM files found")

  missing_idx <- !vapply(bams, bam_has_index, logical(1))
  if (any(missing_idx)) stop("[ERROR] Missing BAM index for: ", paste(bams[missing_idx], collapse = ", "))

  # ── Sample name extraction ────────────────────────────────────────────────
  if (!is.null(custom_sample_names)) {
    if (length(custom_sample_names) != length(bams)) stop("[ERROR] custom_sample_names length mismatch")
    sample_names <- custom_sample_names
  } else {
    base_names <- tools::file_path_sans_ext(basename(bams))

    if (is.null(sample_name_collapse)) {
      sample_name_collapse <- if (nchar(sample_name_delim) == 1 && sample_name_delim != "\\")
        sample_name_delim else "."
    }

    keep_parts <- function(x, delim, keep_spec, collapse) {
      parts <- strsplit(x, delim)[[1]]
      if (grepl("-", keep_spec)) {
        idx   <- as.numeric(strsplit(keep_spec, "-")[[1]])
        parts <- parts[idx[1]:idx[2]]
      } else {
        idx   <- as.numeric(keep_spec)
        if (is.na(idx) || idx < 1L || idx > length(parts)) {
          warning("[WARNING] sample_name_keep index ", idx, " out of range for '",
                  x, "' (", length(parts), " parts). Using last part.")
          idx <- length(parts)
        }
        parts <- parts[idx]
      }
      paste(parts, collapse = collapse)
    }

    sample_names <- vapply(base_names, function(fname) {
      keep_parts(fname, sample_name_delim, sample_name_keep, sample_name_collapse)
    }, character(1))
    names(sample_names) <- NULL

    if (verbose) {
      message("[INFO] Extracted sample names using delimiter: ", sample_name_delim,
              ", keeping part(s): ", sample_name_keep,
              ", collapse: ", sample_name_collapse)
    }
  }
  if (any(duplicated(sample_names))) stop("[ERROR] Duplicate sample names after extraction")

  # ── BED loading ──────────────────────────────────────────────────────────
  bed_file <- data.table::fread(bed, header = FALSE)
  if (ncol(bed_file) < 4) stop("[ERROR] BED file must have at least 4 columns (chr, start, end, name)")

  core_cols <- c("chromosome", "start", "end", "gene")
  if (ncol(bed_file) >= 5) {
    bed_file <- bed_file[, 1:5, with = FALSE]
    colnames(bed_file) <- c(core_cols, "exon_number")
  } else {
    bed_file <- bed_file[, 1:4, with = FALSE]
    colnames(bed_file) <- core_cols
  }

  bed_file$chromosome <- trimws(as.character(bed_file$chromosome))

  if (bed_zero_based) {
    bed_file$start <- bed_file$start + 1L
    if (verbose) message("[INFO] Converted BED start coordinates from 0-based to 1-based")
  }

  map_chrom <- function(bed_chr) {
    if (bed_chr %in% fasta_chroms) return(bed_chr)
    with_chr    <- ifelse(grepl("^chr", bed_chr), bed_chr, paste0("chr", bed_chr))
    if (with_chr %in% fasta_chroms) return(with_chr)
    without_chr <- sub("^chr", "", bed_chr)
    if (without_chr %in% fasta_chroms) return(without_chr)
    return(NA_character_)
  }

  bed_file$original_chrom <- bed_file$chromosome
  bed_file$chromosome     <- vapply(bed_file$chromosome, map_chrom, character(1))
  removed_rows            <- is.na(bed_file$chromosome)
  if (any(removed_rows)) {
    if (verbose) message("[WARNING] Removed ", sum(removed_rows), " rows with chromosomes not in FASTA")
    bed_file <- bed_file[!removed_rows, ]
  }

  bed_file$start <- as.numeric(bed_file$start)
  bed_file$end   <- as.numeric(bed_file$end)
  valid_coords   <- which(bed_file$start > 0 & bed_file$end >= bed_file$start)
  if (length(valid_coords) < nrow(bed_file)) {
    message("[WARNING] Removing ", nrow(bed_file) - length(valid_coords), " rows with invalid coordinates")
    bed_file <- bed_file[valid_coords, ]
  }

  chr_len       <- fasta_lengths[bed_file$chromosome]
  out_of_bounds <- which(bed_file$end > chr_len | bed_file$start < 1)
  if (length(out_of_bounds) > 0) {
    message("[WARNING] Removing ", length(out_of_bounds), " rows outside chromosome boundaries")
    bed_file <- bed_file[-out_of_bounds, ]
  }

  # ── Terminal-exon padding (optional) ───────────────────────────────────────
  # Applied here -- after chromosomes are mapped/filtered and coordinates are
  # clean numerics, but before the reference-sequence extractability check
  # below -- so the (already boundary-aware) padded coordinates still get
  # validated against the FASTA/BSgenome sequence for free, and so every
  # downstream step (CNV calling, plots, VCF, report) sees the padded
  # windows via the bed_file saved in data_out, consistent with how they
  # already consume bed_file rather than re-reading the original BED path.
  if (nrow(bed_file) > 0 && !is.null(pad_terminal_exons) && !is.na(pad_terminal_exons) && pad_terminal_exons > 0) {
    bed_file <- pad_gene_terminal_exons(bed_file, padding = pad_terminal_exons,
                                         chr_lengths = fasta_lengths, verbose = verbose)
  }

  if (nrow(bed_file) > 0) {
    # Always test extractability against the reference, regardless of
    # skip_invalid_intervals. Previously this whole block was gated on
    # that flag, so when it was FALSE (as in some configs) a single bad
    # interval would sail through undetected here and only surface much
    # later as an opaque crash from deep inside ExomeDepth::getBamCounts()
    # (after that call had already started working through every BAM
    # file). Running the check unconditionally means bad intervals are
    # always caught here, cheaply, before the expensive step begins.
    # skip_invalid_intervals now only controls what happens with what's
    # found: TRUE drops the offending rows and continues; FALSE stops
    # with a precise, actionable report instead of letting ExomeDepth's
    # cryptic "record N failed" reach the user.
    if (verbose) message("[INFO] Validating ", nrow(bed_file), " interval(s) against the reference sequence...")
    keep       <- logical(nrow(bed_file))
    chunk_size <- 5000L
    for (chunk_start in seq(1L, nrow(bed_file), chunk_size)) {
      chunk_end   <- min(chunk_start + chunk_size - 1L, nrow(bed_file))
      chunk_idx   <- chunk_start:chunk_end
      gr <- GenomicRanges::GRanges(
        seqnames = bed_file$chromosome[chunk_idx],
        ranges   = IRanges::IRanges(start = bed_file$start[chunk_idx], end = bed_file$end[chunk_idx])
      )
      chunk_ok <- tryCatch({
        sq <- Rsamtools::scanFa(fa, param = gr)
        BiocGenerics::width(sq) > 0L
      }, error = function(e) {
        # Fall back to per-row testing only if the whole-chunk call errors
        # (e.g. one bad interval in the batch), so a single bad row doesn't
        # cost the vectorized speedup for every other chunk.
        vapply(chunk_idx, function(i) {
          tryCatch({
            sq_i <- Rsamtools::scanFa(fa, param = gr[match(i, chunk_idx)])
            length(sq_i) == 1L && BiocGenerics::width(sq_i)[1] > 0L
          }, error = function(e) FALSE)
        }, logical(1))
      })
      keep[chunk_idx] <- chunk_ok
    }

    if (!all(keep)) {
      removed  <- which(!keep)
      offender <- utils::head(removed, 10)
      detail   <- paste0("   Row ", offender, ": ", bed_file$chromosome[offender], ":",
                          bed_file$start[offender], "-", bed_file$end[offender], collapse = "\n")
      if (length(removed) > length(offender)) {
        detail <- paste0(detail, "\n   ... and ", length(removed) - length(offender), " more.")
      }

      if (skip_invalid_intervals) {
        message("[WARNING] Removing ", length(removed),
                " interval(s) that cannot be extracted from the reference sequence:\n", detail)
        bed_file <- bed_file[keep, ]
      } else {
        stop("[ERROR] ", length(removed), " interval(s) cannot be extracted from the reference sequence ",
             "(this usually means the FASTA file or its .fai index is stale, truncated, or simply ",
             "doesn't cover these coordinates):\n", detail,
             "\nOptions: set skip_invalid_intervals = TRUE to drop these rows automatically and continue; ",
             "delete the FASTA's .fai and let ECHO rebuild it (in case the index was stale); or set ",
             "fasta_source = 'bsgenome' to source sequence from an installed Bioconductor BSgenome ",
             "package instead of the local FASTA file.")
      }
    } else if (verbose) {
      message("[INFO] All intervals validated successfully against the reference sequence.")
    }
  }

  chr_levels      <- fasta_chroms[fasta_chroms %in% unique(bed_file$chromosome)]
  bed_file$chromosome <- factor(bed_file$chromosome, levels = chr_levels)
  data.table::setorder(bed_file, chromosome, start, end)
  bed_file$chromosome     <- as.character(bed_file$chromosome)
  bed_file$original_chrom <- NULL
  if (verbose) message("[INFO] BED file cleaned. Final rows: ", nrow(bed_file))
  if (nrow(bed_file) == 0) stop("[ERROR] No valid intervals remain.")

  counts <- tryCatch({
    suppressWarnings(
      ExomeDepth::getBamCounts(
        bed.frame      = as.data.frame(bed_file),
        bam.files      = bams,
        include.chr    = include.chr,
        referenceFasta = fasta
      )
    )
  }, error = function(e) {
    msg     <- conditionMessage(e)
    row_num <- suppressWarnings(as.numeric(gsub(".*record ([0-9]+).*", "\\1", msg)))
    if (!is.na(row_num) && row_num <= nrow(bed_file)) {
      # This is a backstop, not the primary defense: every interval already
      # passed the scanFa validation above, so reaching here means
      # ExomeDepth::getBamCounts() used a slightly different/larger window
      # internally than what was tested (rare, but possible).
      stop("[ERROR] Interval extraction failed for row ", row_num, ": ",
           bed_file$chromosome[row_num], ":", bed_file$start[row_num], "-", bed_file$end[row_num],
           " (unexpected: this interval passed the earlier reference-validation check; ExomeDepth's ",
           "internal window may extend slightly beyond it).")
    } else stop(e)
  })

  normalised_bams <- vapply(bams, function(b) make.names(basename(b)), character(1))
  sample_cols     <- intersect(colnames(counts), normalised_bams)
  if (length(sample_cols) != length(bams))
    sample_cols <- intersect(colnames(counts), paste0("X", normalised_bams))
  if (length(sample_cols) != length(bams)) {
    warning("[WARNING] Sample name matching failed; falling back to column order. Verify sample names.")
    metadata_cols <- c("chromosome", "start", "end", "exon", "GC", "exon_number")
    sample_cols   <- setdiff(colnames(counts), metadata_cols)
    sample_cols   <- sample_cols[order(match(normalised_bams, sample_cols))]
  }
  name_mapping <- setNames(sample_names, sample_cols)
  for (old_name in sample_cols) {
    colnames(counts)[colnames(counts) == old_name] <- name_mapping[old_name]
  }

  keep_cols <- c("chromosome", "start", "end", "exon", sample_names)
  if ("GC" %in% colnames(counts)) keep_cols <- c(keep_cols, "GC")
  if ("exon_number" %in% colnames(counts)) keep_cols <- c(keep_cols, "exon_number")
  counts <- counts[, keep_cols, drop = FALSE]

  if ("exon_number" %in% colnames(bed_file)) {
    counts$exon_number <- bed_file$exon_number
  }

  dup_mask <- duplicated(counts[, c("chromosome", "start", "end", "exon")])
  if (any(dup_mask)) {
    message("[WARNING] Removed ", sum(dup_mask), " duplicate target region(s)")
    counts   <- counts[!dup_mask, , drop = FALSE]
    bed_file <- bed_file[!dup_mask, ]
  }

  if (nrow(bed_file) != nrow(counts)) stop("[ERROR] bed_file and counts row counts do not match")
  bed_file <- as.data.frame(bed_file)
  counts   <- as.data.frame(counts)

  bed_file           <- assign_exon_numbers_per_gene(bed_file)
  counts$exon_number <- bed_file$exon_number
  if ("GC" %in% colnames(counts)) bed_file$GC <- counts$GC

  save(counts, bams, bed_file, sample_names, fasta, name_mapping, file = data_out)
  if (verbose) message("[INFO] ", Sys.time(), " END bam coverage calculation")
  invisible(NULL)
}