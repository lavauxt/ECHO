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

  direct_pkg <- is.character(fasta_source) && length(fasta_source) == 1 &&
    !fasta_source %in% c("file", "bsgenome")
  if (direct_pkg) {
    pkg <- fasta_source
  } else {
    fasta_source <- match.arg(fasta_source)
  }

  if (fasta_source == "file") {
    stop_if_not_file(fasta, "[ERROR] FASTA file not found")
    return(fasta)
  }

  if (!direct_pkg) {
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

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("[ERROR] Package '", pkg, "' is required for fasta_source = '",
         if (direct_pkg) pkg else "bsgenome",
         "' but is not installed.\n",
         "  Install with: BiocManager::install(\"", pkg, "\")\n",
         "  Or use fasta_source = 'file' with your own FASTA.")
  }

  bsgenome        <- getExportedValue(pkg, pkg)
  genome_seqnames <- GenomeInfoDb::seqnames(bsgenome)
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

  cache_tag <- if (direct_pkg) pkg else genome_version
  cache_fasta <- file.path(cache_dir, paste0("ECHO_ref_", cache_tag, "_bsgenome.fa"))
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
    # Append only to a cache whose index we could read (existing_chroms is
    # non-empty). A FASTA without a usable .fai is rewritten from scratch:
    # appending to it would duplicate sequence names and break indexFa().
    Biostrings::writeXStringSet(new_seqs, filepath = cache_fasta,
                                append = length(existing_chroms) > 0L)
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
#' @param rbams Character string. Optional path to external reference BAMs,
#'   either (a) a directory containing \code{.bam} files, or (b) a TSV file
#'   with a \code{bam} column listing reference BAM paths.
#' @param data_out Character string. Output \code{.RData} file path.
#' @param include.chr Logical. Passed to \code{ExomeDepth::getBamCounts}.
#' @param verbose Logical. Print progress messages.
#' @param sample_name_delim Character string. Regex delimiter to split filenames.
#' @param sample_name_keep Character string. Which parts to keep (e.g., "1", "1-2").
#' @param sample_name_collapse Character string. Separator for rejoining parts. Default NULL.
#' @param custom_sample_names Optional character vector, one name per analysis
#'   BAM. Reference BAMs use filename-derived names; REF_ is added only if a
#'   reference name collides with an analysis sample name.
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
  data_out = "./ECHO_coverage.RData",
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

  # Keep analysis BAMs separate from reference-only BAMs. External references
  # contribute counts to ExomeDepth reference selection but must never become
  # test samples or receive their own CNV calls.
  test_bams <- character()
  if (!is.null(bamfiles) && length(bamfiles) == 1L && file.exists(bamfiles)) {
    test_bams <- c(test_bams, trimws(readLines(bamfiles, warn = FALSE)))
  }
  if (!is.null(bamdir) && length(bamdir) == 1L && dir.exists(bamdir)) {
    test_bams <- c(test_bams, list.files(bamdir, pattern = "\\.bam$", full.names = TRUE))
  }
  test_bams <- unique(test_bams[!is.na(test_bams) & nzchar(test_bams)])

  reference_bams <- character()
  if (!is.null(rbams)) {
    if (length(rbams) == 1L && dir.exists(rbams)) {
      reference_bams <- list.files(rbams, pattern = "\\.bam$", full.names = TRUE)
      if (verbose) message("[INFO] Found ", length(reference_bams),
                           " external reference BAM(s) in directory: ", rbams)
    } else if (length(rbams) == 1L && file.exists(rbams)) {
      ref_table <- utils::read.delim(rbams, header = TRUE, stringsAsFactors = FALSE,
                                     check.names = FALSE)
      if (!"bam" %in% names(ref_table)) {
        stop("[ERROR] External reference BAM list must contain a 'bam' column: ", rbams)
      }
      reference_bams <- trimws(as.character(ref_table$bam))
    } else {
      warning("[WARNING] rbams path not found, skipping external reference BAMs: ", rbams)
    }
  }
  reference_bams <- unique(reference_bams[!is.na(reference_bams) & nzchar(reference_bams)])
  # A path supplied in both sources is an analysis BAM, not a reference-only BAM.
  reference_bams <- setdiff(reference_bams, test_bams)
  if (!length(test_bams)) stop("[ERROR] No analysis BAM files found (bamfiles/bamdir is empty)")

  bams <- c(test_bams, reference_bams)
  missing_idx <- !vapply(bams, bam_has_index, logical(1))
  if (any(missing_idx)) stop("[ERROR] Missing BAM index for: ", paste(bams[missing_idx], collapse = ", "))

  if (is.null(sample_name_collapse)) {
    sample_name_collapse <- if (nchar(sample_name_delim) == 1 && sample_name_delim != "\\")
      sample_name_delim else "."
  }
  derive_sample_names <- function(paths) {
    if (!length(paths)) return(character(0))
    base_names <- basename(paths)
    base_names <- sub("\\.(cram|bam)$", "", base_names, ignore.case = TRUE)
    base_names <- sub("\\.(sort|markdup|dedup|recal|fixmate|rg|merged)$", "",
                      base_names, ignore.case = TRUE)
    keep_parts <- function(x, delim, keep_spec, collapse) {
      parts <- strsplit(x, delim)[[1]]
      n <- length(parts)
      if (grepl("-", keep_spec)) {
        idx <- suppressWarnings(as.numeric(strsplit(keep_spec, "-")[[1]]))
        idx <- idx[!is.na(idx)]
        if (length(idx) != 2L) {
          warning("[WARNING] Malformed sample_name_keep range '", keep_spec,
                  "' for '", x, "'. Using last part.")
          return(parts[n])
        }
        lo <- max(1L, min(idx[1], n)); hi <- max(1L, min(idx[2], n))
        if (lo > hi) {
          warning("[WARNING] sample_name_keep range '", keep_spec,
                  "' exceeds filename '", x, "' (", n, " parts). Using last part.")
          return(parts[n])
        }
        parts <- parts[lo:hi]
      } else {
        idx <- suppressWarnings(as.numeric(keep_spec))
        if (is.na(idx) || idx < 1L || idx > n) {
          warning("[WARNING] sample_name_keep index ", idx, " out of range for '",
                  x, "' (", n, " parts). Using last part.")
          idx <- n
        }
        parts <- parts[idx]
      }
      paste(parts, collapse = collapse)
    }
    out <- vapply(base_names, keep_parts, character(1),
                  delim = sample_name_delim, keep_spec = sample_name_keep,
                  collapse = sample_name_collapse)
    unname(out)
  }

  if (!is.null(custom_sample_names)) {
    if (length(custom_sample_names) != length(test_bams)) {
      stop("[ERROR] custom_sample_names must contain one name per analysis BAM (",
           length(test_bams), "), not per external reference BAM.")
    }
    sample_names <- as.character(custom_sample_names)
  } else {
    sample_names <- derive_sample_names(test_bams)
  }
  reference_sample_names <- derive_sample_names(reference_bams)
  if (anyDuplicated(reference_sample_names)) {
    duplicate_refs <- unique(reference_sample_names[duplicated(reference_sample_names)])
    stop("[ERROR] External reference BAMs produce duplicate sample names: ",
         paste(duplicate_refs, collapse = ", "),
         ". Adjust sample_name_delim/sample_name_keep so each reference has a unique name.")
  }
  ref_name_collisions <- reference_sample_names %in% sample_names
  if (any(ref_name_collisions)) {
    warning("[WARNING] External reference sample names collide with analysis samples; ",
            "adding REF_ to: ", paste(reference_sample_names[ref_name_collisions], collapse = ", "),
            ". For mode X/Y, use these REF_-prefixed names in sample_table.")
    reference_sample_names[ref_name_collisions] <- paste0("REF_", reference_sample_names[ref_name_collisions])
  }
  all_bam_sample_names <- c(sample_names, reference_sample_names)
  if (anyNA(all_bam_sample_names) || any(!nzchar(all_bam_sample_names))) {
    stop("[ERROR] Sample names cannot be NA or empty.")
  }
  if (anyDuplicated(sample_names)) stop("[ERROR] Duplicate analysis sample names after extraction")
  if (anyDuplicated(all_bam_sample_names)) stop("[ERROR] Sample name collision between analysis and reference BAMs")
  reserved_names <- c("chromosome", "start", "end", "exon", "GC", "exon_number", "exon_label")
  if (any(all_bam_sample_names %in% reserved_names)) {
    stop("[ERROR] Sample name(s) clash with reserved column names (",
         paste(intersect(all_bam_sample_names, reserved_names), collapse = ", "),
         "); change sample_name_delim / sample_name_keep or use custom_sample_names.")
  }
  if (verbose) {
    message("[INFO] Identified ", length(sample_names), " analysis sample(s) and ",
            length(reference_sample_names), " external reference sample(s).")
    message("[INFO] Sample names extracted using delimiter: ", sample_name_delim,
            ", keeping part(s): ", sample_name_keep,
            ", collapse: ", sample_name_collapse)
  }

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

  if (nrow(bed_file) > 0 && !is.null(pad_terminal_exons) && !is.na(pad_terminal_exons) && pad_terminal_exons > 0) {
    bed_file <- pad_gene_terminal_exons(bed_file, padding = pad_terminal_exons,
                                         chr_lengths = fasta_lengths, verbose = verbose)
  }

  if (nrow(bed_file) > 0) {
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

  # ExomeDepth::getBamCounts() follows the BED convention: it expects 0-based
  # starts and adds 1 itself (IRanges(start = start + 1, end = end)). bed_file
  # holds 1-based starts at this point (needed for the FASTA checks and the
  # padding above), so hand ExomeDepth start - 1. Passing the 1-based start
  # made every target 1 bp too short and shifted every reported Start (and the
  # VCF POS) by one.
  bed_for_counts       <- as.data.frame(bed_file)
  bed_for_counts$start <- bed_for_counts$start - 1

  counts <- tryCatch({
    suppressWarnings(
      ExomeDepth::getBamCounts(
        bed.frame      = bed_for_counts,
        bam.files      = bams,
        include.chr    = include.chr,
        referenceFasta = fasta
      )
    )
  }, error = function(e) {
    msg     <- conditionMessage(e)
    row_num <- suppressWarnings(as.numeric(gsub(".*record ([0-9]+).*", "\\1", msg)))
    if (!is.na(row_num) && row_num <= nrow(bed_file)) {
      stop("[ERROR] Interval extraction failed for row ", row_num, ": ",
           bed_file$chromosome[row_num], ":", bed_file$start[row_num], "-", bed_file$end[row_num],
           " (unexpected: this interval passed the earlier reference-validation check; ExomeDepth's ",
           "internal window may extend slightly beyond it).")
    } else stop(e)
  })

  # One count column per BAM, named make.names(basename(bam)) by data.frame()
  # (prefixed with "X" when the name starts with a digit). Match each BAM to
  # its own column, so that analysis and reference names cannot drift from
  # the order of the BAMs/columns.
  normalised_bams <- unname(vapply(bams, function(b) make.names(basename(b)), character(1)))
  sample_cols <- vapply(normalised_bams, function(nb) {
    cand <- c(nb, paste0("X", nb))
    cand <- cand[cand %in% colnames(counts)]
    if (length(cand)) cand[1] else NA_character_
  }, character(1), USE.NAMES = FALSE)
  if (anyNA(sample_cols) || anyDuplicated(sample_cols)) {
    metadata_cols <- c("chromosome", "start", "end", "exon", "GC", "exon_number")
    positional    <- setdiff(colnames(counts), metadata_cols)
    if (length(positional) != length(bams)) {
      stop("[ERROR] Could not match the ", length(bams), " BAM file(s) to the ",
           length(positional), " count column(s) returned by ExomeDepth::getBamCounts().")
    }
    warning("[WARNING] Sample name matching failed; assuming getBamCounts() returned the BAM ",
            "columns in input order. Verify sample names.")
    sample_cols <- positional
  }
  name_mapping <- setNames(all_bam_sample_names, sample_cols)
  colnames(counts)[match(sample_cols, colnames(counts))] <- all_bam_sample_names

  keep_cols <- c("chromosome", "start", "end", "exon", all_bam_sample_names)
  if ("GC" %in% colnames(counts)) keep_cols <- c(keep_cols, "GC")
  if ("exon_number" %in% colnames(counts)) keep_cols <- c(keep_cols, "exon_number")
  counts <- counts[, keep_cols, drop = FALSE]

  if (nrow(bed_file) != nrow(counts)) stop("[ERROR] bed_file and counts row counts do not match")
  bed_file <- as.data.frame(bed_file)
  counts   <- as.data.frame(counts)

  # getBamCounts() re-sorts the targets by (chromosome, start + end), which is
  # not the (start, end) order used above when targets overlap or are nested.
  # Everything downstream is positional (counts row i <-> bed_file row i), so
  # align the annotation to the count matrix explicitly instead of trusting
  # that the two orders agree.
  strip_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
  interval_key <- function(chr, start, end) paste(strip_chr(chr), as.numeric(start), as.numeric(end), sep = "\r")
  occurrence   <- function(k) stats::ave(seq_along(k), k, FUN = seq_along)
  k_bed <- interval_key(bed_file$chromosome, bed_file$start, bed_file$end)
  k_cnt <- interval_key(counts$chromosome, counts$start, counts$end)
  if (!identical(k_bed, k_cnt)) {
    ord <- match(paste(k_cnt, occurrence(k_cnt)), paste(k_bed, occurrence(k_bed)))
    if (anyNA(ord)) {
      stop("[ERROR] The intervals returned by ExomeDepth::getBamCounts() do not correspond to the BED ",
           "rows (coordinates differ); cannot align the target annotation to the count matrix.")
    }
    if (verbose) message("[INFO] Re-ordering the target annotation to follow the count matrix ",
                         "(ExomeDepth sorts targets by start + end).")
    bed_file <- bed_file[ord, , drop = FALSE]
    rownames(bed_file) <- NULL
  }

  # Preserve user-curated exon numbering (col 5 of the BED) under a distinct
  # column name, because assign_exon_numbers_per_gene() overwrites `exon_number`
  # with a derived 1..N numbering (needed for terminal_only detection).
  if ("exon_number" %in% colnames(bed_file) &&
      any(!is.na(bed_file$exon_number))) {
    message("[INFO] Preserving user-supplied exon_number from BED column 5 ",
            "as `exon_label`; derived exon numbering stored as `exon_number`.")
    bed_file$exon_label <- bed_file$exon_number
  }

  bed_file <- assign_exon_numbers_per_gene(bed_file)
  counts$exon_number <- bed_file$exon_number
  if ("exon_label" %in% colnames(bed_file)) {
    counts$exon_label <- bed_file$exon_label
    bed_file$exon_label <- NULL
  }
  if ("GC" %in% colnames(counts)) bed_file$GC <- counts$GC

  save(counts, bams, bed_file, sample_names, reference_sample_names,
       fasta, name_mapping, file = data_out)
  if (verbose) message("[INFO] ", Sys.time(), " END bam coverage calculation")
  invisible(NULL)
}