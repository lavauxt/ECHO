
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
load_rdata  <- function(path, required = NULL) {
stop_if_not_file(path, paste0("[ERROR] RData file not found: ", path))
env  <- new.env()
load(path, envir = env)
if (!is.null(required)) {
missing  <- setdiff(required, ls(env))
if (length(missing)  > 0) {
stop("[ERROR] Missing variables in RData: ", paste(missing, collapse = ", "),
"\nAvailable objects: ", paste(ls(env), collapse = ", "))
}
}
as.list(env)
}
#' Extract gene symbols from structured BED target annotations
#'
#' Recognises RefSeq transcript prefixes, exon/UTR labels, score annotations,
#' and trailing genomic coordinates. Underscores and hyphens inside gene
#' symbols are preserved (for example, SCG5_GREM1 and BCR-ABL1).
#' Multiple comma-separated annotations are reduced to their unique symbols.
#' @param x Character vector of BED annotation names.
#' @param collapse Separator used when an interval maps to distinct symbols.
#' @return Character vector; NA where no structured annotation was recognised.
#' @noRd
extract_bed_gene_symbols <- function(x, collapse = ";") {
    if (is.null(x)) return(character(0))

    extract_one <- function(annotation) {
        annotation <- trimws(as.character(annotation))
        if (is.na(annotation) || !nzchar(annotation) || annotation == ".") {
            return(NA_character_)
        }

        has_coordinates <- grepl(
            "_(?:chr[^_]+|[0-9]+|X|Y|M|MT)_[0-9]+_[0-9]+$",
            annotation, perl = TRUE, ignore.case = TRUE
        )
        has_refseq <- grepl("^(?:NM|NR|XM|XR)_[0-9]+_", annotation,
                            perl = TRUE, ignore.case = TRUE)
        has_ensembl <- grepl("^ENST[0-9]+(?:\\.[0-9]+)?_", annotation,
                             perl = TRUE, ignore.case = TRUE)
        has_exon_or_utr <- grepl(
            "_(?:(?:5|3)?utr[0-9]*(?:[-_]ex[0-9]+)?|ex(?:on)?[_-]?[0-9]+)(?:$|[_(-])",
            annotation, perl = TRUE, ignore.case = TRUE
        )

        # Do not guess at custom formats. The caller can use its legacy,
        # explicitly configured field selection for those names.
        if (!(has_coordinates || has_refseq || has_ensembl || has_exon_or_utr)) {
            return(NA_character_)
        }

        # Remove the terminal coordinate triplet before symbol extraction.
        annotation <- sub(
            "_(?:chr[^_]+|[0-9]+|X|Y|M|MT)_[0-9]+_[0-9]+$",
            "", annotation, perl = TRUE, ignore.case = TRUE
        )
        # Scores such as (251/4972) are metadata, not part of the symbol.
        annotation <- gsub("\\s*\\([^()]*\\)", "", annotation, perl = TRUE)

        # Strip transcript accessions, not arbitrary underscore-delimited
        # prefixes. This is what protects fused symbols such as SCG5_GREM1.
        annotation <- sub("^(?:NM|NR|XM|XR)_[0-9]+_", "", annotation,
                          perl = TRUE, ignore.case = TRUE)
        annotation <- sub("^ENST[0-9]+(?:\\.[0-9]+)?_", "", annotation,
                          perl = TRUE, ignore.case = TRUE)

        # Remove structural exon/UTR suffixes while leaving the complete
        # preceding symbol (including internal underscores and hyphens).
        annotation <- sub(
            "_(?:(?:5|3)?utr[0-9]*(?:[-_]ex[0-9]+)?|ex(?:on)?[_-]?[0-9]+)(?:[_-].*)?$",
            "", annotation, perl = TRUE, ignore.case = TRUE
        )
        annotation <- trimws(annotation)
        if (!nzchar(annotation) || annotation == ".") return(NA_character_)
        annotation
    }

    vapply(as.character(x), function(value) {
        if (is.na(value) || !nzchar(trimws(value)) || trimws(value) == ".") {
            return(NA_character_)
        }
        annotations <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
        genes <- vapply(annotations, extract_one, character(1))
        genes <- genes[!is.na(genes) & nzchar(genes)]
        if (!length(genes)) return(NA_character_)
        paste(unique(genes), collapse = collapse)
    }, character(1), USE.NAMES = FALSE)
}

#' Check that a BAM file has an index on disk
#' @param bam Path to a BAM file.
#' @return Logical.
#' @noRd
bam_has_index <- function(bam) {
# Rsamtools accepts both "x.bam.bai" and "x.bai". The old second test looked
# for "x.bam.bam.bai", which never exists, so BAMs indexed as "x.bai" were
# wrongly rejected.
file.exists(paste0(bam, ".bai")) ||
  file.exists(sub("\\.bam$", ".bai", bam, ignore.case = TRUE))
}
#' Normalise a chromosome name to match a reference naming style (vectorised)
#' @param chr_vec Character vector of chromosome names.
#' @param ref_chromosomes Character vector of reference chromosome names.
#' @return Character vector with prefix aligned to reference style.
#' @noRd
normalize_chromosome_vec  <- function(chr_vec, ref_chromosomes) {
has_chr_ref  <- any(grepl("^chr", ref_chromosomes))
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
filter_chromosomes  <- function(df, include = NULL, exclude = NULL) {
if (is.null(df) || nrow(df) == 0) return(df)
chrom_col  <- intersect(c("chromosome", "Chromosome"), colnames(df))[1]
if (is.na(chrom_col)) {
stop("Data frame must contain a 'chromosome' or 'Chromosome' column.")
}
norm  <- function(x) {
unique(c(x, sub("^chr", "", x), paste0("chr", sub("^chr", "", x))))
}
if (!is.null(include)) {
df  <- df[df[[chrom_col]] %in% norm(include), ]
}
if (!is.null(exclude)) {
df  <- df[!df[[chrom_col]] %in% norm(exclude), ]
}
df
}

#' Quantile function of the Beta-Binomial distribution
#'
#' NOTE: uses the large-`size` approximation `qbeta(p, a, b) * size`.
#' For thin panels with small per-exon read counts the returned intervals
#' may deviate from the true beta-binomial quantiles. Use
#' \code{VGAM::qbetabinom.ab()} for exact quantiles if accuracy in that
#' regime matters.
#' @noRd
qbetabinom <- function(p, size, rho, prob) {
a <- prob * (1 - rho) / rho
b <- (1 - prob) * (1 - rho) / rho
qbeta(p, a, b) * size
}

#' Null coalescing operator
#' @export
`%||%`<- function(a, b) {
if (!is.null(a)) a else b
}
#' @import data.table
NULL

#' Column names of a BED-like data frame (lowercase or capitalised style)
#' @noRd
bed_cols <- function(df) {
  pick <- function(a, b) if (a %in% names(df)) a else b
  list(chrom = pick("chromosome", "Chr"), start = pick("start", "Start"),
       end = pick("end", "End"), gene = pick("gene", "Gene"))
}

#' Natural chromosome order: 1-22, X, Y, M (with and without "chr"), then the rest
#' @noRd
chrom_sort_levels <- function(chrs = character()) {
  base <- c(as.character(1:22), "X", "Y", "M")
  unique(c(paste0("chr", base), base, unique(as.character(chrs))))
}

#' Label given to retained off-target (HorsROI) targets
#' @noRd
OFF_TARGET_LABEL <- "HorsROI"

#' Which gene labels denote off-target (HorsROI) targets
#'
#' @param genes Character vector of gene labels.
#' @param pattern Optional extra regex (the run's \code{off_target_pattern}).
#' @return Logical vector.
#' @noRd
is_off_target_gene <- function(genes, pattern = "^HorsROI") {
  g   <- as.character(genes)
  hit <- !is.na(g) & toupper(g) == toupper(OFF_TARGET_LABEL)
  if (length(pattern) == 1L && !is.na(pattern) && nzchar(pattern)) {
    hit <- hit | (!is.na(g) & grepl(pattern, g, perl = TRUE))
  }
  hit
}

#' Column / row medians (matrixStats when installed)
#' @noRd
col_medians <- function(m) {
  if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::colMedians(m, na.rm = TRUE)
  else apply(m, 2, stats::median, na.rm = TRUE)
}
#' @noRd
row_medians <- function(m) {
  if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::rowMedians(m, na.rm = TRUE)
  else apply(m, 1, stats::median, na.rm = TRUE)
}

#' Genes covered by target index ranges
#'
#' @param genes Gene label of every target.
#' @param lo,hi Target index ranges (one pair per call).
#' @param drop_unknown Also drop \code{"Unknown"} labels.
#' @return One space-separated gene string per range (\code{"Unknown"} if none).
#' @noRd
genes_in_targets <- function(genes, lo, hi, drop_unknown = FALSE) {
  n <- length(genes)
  vapply(seq_along(lo), function(j) {
    a <- suppressWarnings(as.integer(lo[j])); b <- suppressWarnings(as.integer(hi[j]))
    if (is.na(a) || is.na(b) || b < a || a > n || b < 1L) return("Unknown")
    g <- unique(as.character(genes[seq.int(max(1L, a), min(b, n))]))
    g <- g[!is.na(g) & nzchar(g) & (!drop_unknown | g != "Unknown")]
    if (length(g)) paste(g, collapse = " ") else "Unknown"
  }, character(1))
}
#' Sanitize a string for use as a filename
#' @export
sanitize_filename <- function(name) {
gsub("[^[:alnum:]]", "_", name)
}
#' Reclassify off-target / filler BED intervals before exon numbering
#' @export
handle_off_target_regions  <- function(df, pattern = "^HorsROI",
handling = c("na", "remove", "merge"),
verbose = TRUE) {
handling  <- match.arg(handling)
if (is.null(pattern) || !nzchar(pattern)) return(df)
if (is.null(df) || nrow(df) == 0) return(df)
bc <- bed_cols(df)
chrom_col <- bc$chrom; start_col <- bc$start; end_col <- bc$end; gene_col <- bc$gene
if (!gene_col %in% names(df)) return(df)

if ("OriginalName" %in% names(df)) {
    name_vals <- as.character(df[["OriginalName"]])
    is_off <- !is.na(name_vals) & grepl(pattern, name_vals, perl = TRUE)
} else {
    gene_vals <- as.character(df[[gene_col]])
    is_off <- !is.na(gene_vals) & grepl(pattern, gene_vals, perl = TRUE)
}

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

# handling == "merge"
ord     <- order(df[[chrom_col]], df[[start_col]], df[[end_col]])
n       <- length(ord)
chrom_s <- as.character(df[[chrom_col]])[ord]
start_s <- df[[start_col]][ord]
end_s   <- df[[end_col]][ord]
gene_s  <- as.character(df[[gene_col]])[ord]
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
    }
}
df[[gene_col]][ord] <- new_gene_s

# <-- NEW: Drop the off-target marker from `OriginalName` for merged intervals
if ("OriginalName" %in% names(df) && any(merged_s)) {
    df[["OriginalName"]][ord[merged_s]] <- new_gene_s[merged_s]
}
    
if (verbose) {
    message(sprintf(
        "[INFO] handle_off_target_regions: merged %d/%d off-target interval(s) matching /%s/ into their nearest neighbouring gene (%d had no real gene on their chromosome to attach to).",
        n_merged, n_off, pattern, n_off - n_merged))
}
df
}

#' Assign sequential exon numbers within each chromosome/gene locus
#' @noRd
assign_exon_numbers_per_gene  <- function(bed_file) {
bc <- bed_cols(bed_file)
chrom_col <- bc$chrom; start_col <- bc$start; end_col <- bc$end; gene_col <- bc$gene
stopifnot(all(c(chrom_col, start_col, end_col, gene_col) %in% names(bed_file)))
n_in <- nrow(bed_file)
dt <- data.table::as.data.table(bed_file)
dt[, .orig_row := .I]
data.table::setnames(dt, c(chrom_col, start_col, end_col, gene_col),
                     c("chrom", "start", "end", "gene"))
dt[, start := as.numeric(start)]
dt[, end   := as.numeric(end)]
dt[, .key := paste(chrom, start, end, gene, sep = "\u241f")]
dup_rows <- duplicated(dt, by = ".key")
if (any(dup_rows)) {
    warning(sprintf(
        "Found %d duplicate row(s); duplicates share the exon_number of their first occurrence.",
        sum(dup_rows)), immediate. = TRUE)
}
dt_unique <- dt[!dup_rows]
dt_unique[, .chrom_fac := factor(chrom, levels = chrom_sort_levels(chrom))]
data.table::setorder(dt_unique, .chrom_fac, start, end)
dt_unique[, .chrom_fac := NULL]

# Retained off-target (HorsROI) intervals are not exons of a gene: no number.
valid_gene_rows <- !is.na(dt_unique$gene) & !dt_unique$gene %in% c("", ".") &
  !is_off_target_gene(dt_unique$gene, pattern = NULL)
dt_unique[, exon_number := NA_integer_]
if (any(valid_gene_rows)) {
    dt_unique[valid_gene_rows, exon_number := seq_len(.N), by = .(chrom, gene)]
}
dt[, exon_number := dt_unique$exon_number[match(.key, dt_unique$.key)]]
data.table::setorder(dt, .orig_row)
dt[, c(".orig_row", ".key") := NULL]
data.table::setnames(dt, c("chrom", "start", "end", "gene"),
                     c(chrom_col, start_col, end_col, gene_col))
out <- as.data.frame(dt)
stopifnot(nrow(out) == n_in)
out
}

#' Pad the outer edge of each gene's first and last exon
#' @export
pad_gene_terminal_exons  <- function(bed_file, padding = 0, chr_lengths = NULL, verbose = TRUE) {
if (is.null(padding) || length(padding) != 1 || is.na(padding) || padding  <= 0) {
return(bed_file)
}
padding  <- as.integer(round(padding))
bc <- bed_cols(bed_file)
chrom_col <- bc$chrom; start_col <- bc$start; end_col <- bc$end; gene_col <- bc$gene
stopifnot(all(c(chrom_col, start_col, end_col, gene_col) %in% names(bed_file)))
n_in <- nrow(bed_file)
if (n_in == 0) return(bed_file)
numbered <- assign_exon_numbers_per_gene(bed_file)
original_exon_number <- if ("exon_number" %in% names(bed_file)) bed_file$exon_number else NULL
dt <- data.table::as.data.table(numbered)
dt[, .orig_row := .I]
data.table::setnames(dt, c(chrom_col, start_col, end_col, gene_col),
                     c("chrom", "start", "end", "gene"))
dt[, is_first := exon_number == 1L]
dt[, is_last  := exon_number == max(exon_number), by = .(chrom, gene)]
# Unnumbered targets (no gene, HorsROI) are never padded.
no_gene <- is.na(dt$exon_number) | is.na(dt$gene) | dt$gene %in% c("", ".", "Unknown")
dt[no_gene, c("is_first", "is_last") := FALSE]
dt[, .chrom_fac := factor(chrom, levels = chrom_sort_levels(chrom))]
data.table::setorder(dt, .chrom_fac, start, end)
n        <- nrow(dt)
chrom_id <- as.integer(dt$.chrom_fac)
start_v  <- dt$start
end_v    <- dt$end
is_first_v <- dt$is_first
is_last_v  <- dt$is_last
right_extend <- integer(n)
left_extend  <- integer(n)
if (n > 1) {
    same_chr_pair <- chrom_id[-n] == chrom_id[-1]
    gap        <- pmax(start_v[-1] - end_v[-n] - 1L, 0L)
    want_left  <- ifelse(is_last_v[-n],  padding, 0L)
    want_right <- ifelse(is_first_v[-1], padding, 0L)
    demand     <- want_left + want_right
    grant_left  <- integer(n - 1L)
    grant_right <- integer(n - 1L)
    has_demand  <- same_chr_pair & demand > 0
    fits        <- has_demand & demand <= gap
    tight       <- has_demand & demand > gap
    grant_left[fits]  <- want_left[fits]
    grant_right[fits] <- want_right[fits]
    grant_left[tight]  <- as.integer(floor(gap[tight] * want_left[tight] / demand[tight]))
    grant_right[tight] <- gap[tight] - grant_left[tight]
    right_extend[-n] <- grant_left
    left_extend[-1]  <- grant_right
}
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
    avail <- ifelse(!is.na(chr_len_here), chr_len_here - end_v[no_next_last], 1e12)
    right_extend[no_next_last] <- as.integer(pmin(padding, pmax(avail, 0)))
}
n_padded_start <- sum(left_extend > 0L)
n_padded_end   <- sum(right_extend > 0L)
n_clamped      <- sum(is_first_v & left_extend  < padding) +
                   sum(is_last_v  & right_extend < padding)
dt[, start := start_v - left_extend]
dt[, end   := end_v   + right_extend]
data.table::setorder(dt, .orig_row)
dt[, c(".orig_row", ".chrom_fac", "is_first", "is_last", "exon_number") := NULL]
data.table::setnames(dt, c("chrom", "start", "end", "gene"),
                     c(chrom_col, start_col, end_col, gene_col))
out <- as.data.frame(dt)
if (!is.null(original_exon_number)) {
    out$exon_number <- original_exon_number
    out <- out[, names(bed_file), drop = FALSE]
}
stopifnot(nrow(out) == n_in)
if (verbose) {
    message(sprintf(
        "[INFO] pad_gene_terminal_exons: requested %d bp padding -- extended %d gene start(s) and %d gene end(s); %d side(s) received less than the full request.",
        padding, n_padded_start, n_padded_end, n_clamped))
}
out
}

#' Compute gap-inserted x-axis positions for CNV window plots
#'
#' A gap of \code{gap} units separates consecutive genes. When \code{bed_file}
#' has a logical \code{off_target} column, a run of consecutive off-target
#' (HorsROI) intervals forms one block, separated from its neighbours.
#' @export
compute_gene_gap_positions  <- function(bed_file, exon_range, gap = 1) {
n <- length(exon_range)
if (n == 0) {
return(data.frame(idx = integer(0), px = numeric(0),
gene_break = logical(0), gene_group = integer(0)))
}
gap  <- if (is.null(gap) || is.na(gap) || gap  < 0) 1 else gap
genes    <- as.character(bed_file[[bed_cols(bed_file)$gene]][exon_range])
genes[is.na(genes)] <- ""
off <- if ("off_target" %in% names(bed_file)) as.logical(bed_file$off_target[exon_range]) else rep(FALSE, n)
off[is.na(off)] <- FALSE
# Consecutive off-target intervals form one block, even with distinct labels.
both_off      <- off[-1] & off[-n]
gene_break    <- c(FALSE, (genes[-1] != genes[-n] & !both_off) | off[-1] != off[-n])
step          <- ifelse(gene_break, 1 + gap, 1)
step[1]       <- 0
px            <- cumsum(step) + 1
data.frame(idx = exon_range, px = px, gene_break = gene_break,
          gene_group = cumsum(gene_break) + 1L, stringsAsFactors = FALSE)
}

#' Convert global exon indices to within-gene indices
#' @export
add_within_gene_indices <- function(cnv_calls, bed_file) {
if (nrow(cnv_calls) == 0) return(cnv_calls)
# Guard against double-conversion: once start.p has been rewritten to
# within-gene indices, it no longer equals global_start.
if ("global_start" %in% names(cnv_calls) &&
    any(cnv_calls$global_start != as.integer(cnv_calls$start.p), na.rm = TRUE)) {
  stop("[ERROR] add_within_gene_indices(): 'start.p' has already been ",
       "converted to within-gene indices; refusing to convert again.")
}
bed_sorted <- assign_exon_numbers_per_gene(bed_file)
exon_numbers <- bed_sorted$exon_number
cnv_calls$global_start <- as.integer(cnv_calls$start.p)
cnv_calls$global_end   <- as.integer(cnv_calls$end.p)
cnv_calls$start.p <- exon_numbers[cnv_calls$global_start]
cnv_calls$end.p   <- exon_numbers[cnv_calls$global_end]
cnv_calls
}

#' Compute exon index for plotting
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
#' @export
load_config  <- function(yaml_path) {
if (!file.exists(yaml_path)) {
found_path  <- system.file(yaml_path, package = "ECHO")
if (found_path == "") {
stop("[ERROR] Config file not found at: ", yaml_path)
}
yaml_path  <- found_path
}
cfg  <- yaml::read_yaml(yaml_path)
required  <- c("input", "output")
if (!all(required %in% names(cfg))) {
stop("[ERROR] YAML must contain 'input' and 'output' sections")
}
return(cfg)
}

#' Plot PCA of Sample Coverage Profiles
#' @export
plot_coverage_pca <- function(counts, sample_names, output_pdf = NULL,
color_by = NULL, scale = TRUE) {
if (length(sample_names) < 3) {
warning("PCA needs at least 3 samples (got ", length(sample_names), "); skipped.")
return(invisible(NULL))
}
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
# The grouping column must exist BEFORE ggplot() is called: the plot keeps its
# own copy of the data, so a column added to pca_df afterwards is never seen
# ("object 'Group' not found") and the coloured PCA silently failed.
has_group <- !is.null(color_by) && length(color_by) == nrow(pca_df)
if (has_group) pca_df$Group <- as.factor(color_by)
p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, label = Sample))
if (has_group) {
    p <- p + ggplot2::aes(colour = Group)
    group_sizes <- table(pca_df$Group)
    if (sum(group_sizes >= 3) >= 1) {
        ell_df <- pca_df[pca_df$Group %in% names(group_sizes)[group_sizes >= 3], , drop = FALSE]
        p <- p + ggplot2::stat_ellipse(data = ell_df, ggplot2::aes(x = PC1, y = PC2, colour = Group),
                                       inherit.aes = FALSE, type = "norm", linetype = "dashed")
    }
}
p <- p +
    (if (has_group) ggplot2::geom_point(size = 3) else ggplot2::geom_point(size = 3, colour = "steelblue")) +
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
if (!has_group) {
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

#' Validate a BED file
#'
#' @param bed_path Path to the BED file.
#' @param verbose Print messages.
#' @param bed_zero_based Logical. With the default \code{TRUE} (UCSC BED), an
#'   interval is invalid if \code{end <= start} or \code{start < 0}. With
#'   \code{FALSE} (1-based, closed intervals) a one-base interval
#'   \code{start == end} is legitimate and only \code{end < start} or
#'   \code{start < 1} is invalid.
#' @export
validate_bed_regions  <- function(bed_path, verbose = TRUE, bed_zero_based = TRUE) {
bed  <- tryCatch(data.table::fread(bed_path, header = FALSE, sep = "\t"),
error = function(e) NULL)
if (is.null(bed) || nrow(bed) == 0) {
issues <- list(unreadable_file = bed_path)
if (verbose) message("[WARNING] Could not read a non-empty BED file for validation: ", bed_path)
return(invisible(issues))
}
if (ncol(bed)  < 3) {
issues <- list(insufficient_columns = ncol(bed))
if (verbose) message("[WARNING] BED file has fewer than 3 columns.")
return(invisible(issues))
}
data.table::setnames(bed, 1:3, c("chrom", "start", "end"))
if (ncol(bed)  >= 4) data.table::setnames(bed, 4, "name") else bed$name  <- NA_character_
issues <- list()
if (verbose) message("[INFO] Validating BED file: ", bed_path, " (", nrow(bed), " rows)")
start_num <- suppressWarnings(as.numeric(as.character(bed$start)))
end_num <- suppressWarnings(as.numeric(as.character(bed$end)))
bad_order <- if (isFALSE(bed_zero_based)) {
    end_num < start_num | start_num < 1
} else {
    end_num <= start_num | start_num < 0
}
bad_coords <- which(!is.finite(start_num) | !is.finite(end_num) |
                    start_num != floor(start_num) | end_num != floor(end_num) |
                    bad_order)
if (length(bad_coords)) {
    issues$invalid_coordinates <- as.data.frame(bed[bad_coords, ])
    if (verbose) message("[WARNING] ", length(bad_coords),
                         if (isFALSE(bed_zero_based))
                           " row(s) have non-integer coordinates, end < start, or a start below 1."
                         else
                           " row(s) have non-integer coordinates, end <= start, or a negative start.")
}
dup_rows <- which(duplicated(bed[, c("chrom", "start", "end")]))
if (length(dup_rows)) {
    issues$duplicate_regions <- as.data.frame(bed[dup_rows, ])
    if (verbose) message("[WARNING] ", length(dup_rows), " duplicate chrom/start/end row(s) found.")
}
if (verbose && length(issues) == 0) message("[INFO] BED validation found no issues.")
invisible(issues)
}

#' Flag Background-Exon Calibration Issues
#' @export
check_background_calibration  <- function(ratio, lo, hi, is_affected, min_n = 5) {
bg  <- !is_affected
n_bg  <- sum(bg, na.rm = TRUE)
if (n_bg == 0) {
return(list(n_background = 0L, n_outside = 0L, pct_outside = NA_real_, flag = FALSE))
}
outside  <- (ratio[bg]  < lo[bg]) | (ratio[bg]  > hi[bg])
outside[is.na(outside)]  <- FALSE
n_outside  <- sum(outside)
pct_outside  <- 100 * n_outside / n_bg
flag  <- n_bg  >= min_n  &&
stats::pbinom(n_outside - 1L, size = n_bg, prob = 0.05, lower.tail = FALSE)  < 0.05
list(n_background = n_bg, n_outside = n_outside, pct_outside = pct_outside, flag = flag)
}

#' Detect Outlier Samples
#' @export
detect_outlier_samples <- function(counts, pseudocount = 0.5, z_threshold = 3) {
counts <- as.matrix(counts)
log_counts   <- log2(counts + pseudocount)
sample_noise <- apply(log_counts, 2, stats::mad, na.rm = TRUE)
med  <- median(sample_noise, na.rm = TRUE)
mad0 <- stats::mad(sample_noise, na.rm = TRUE)
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
        stop("'gc' must have one value per row of count_matrix")
    gc_val       <- if (max(gc, na.rm = TRUE) > 1) gc / 100 else gc
    problematic  <- problematic | gc_val < gc_min | gc_val > gc_max |
        !is.finite(gc_val)
}
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

#' Flag CNV calls that affect only terminal exons
#' @export
flag_terminal_only_calls <- function(cnv_calls, bed_file) {
if (nrow(cnv_calls) == 0) return(cnv_calls)
if (!"exon_number" %in% colnames(bed_file)) {
bed_file <- assign_exon_numbers_per_gene(bed_file)
}
cnv_calls$terminal_only <- FALSE
chr_col <- bed_cols(bed_file)$chrom
for (i in seq_len(nrow(cnv_calls))) {
start_idx <- as.integer(cnv_calls$global_start[i])
end_idx   <- as.integer(cnv_calls$global_end[i])
if (is.na(start_idx) || is.na(end_idx)) next
affected_exons <- bed_file$exon_number[start_idx:end_idx]
if (length(affected_exons) == 0) next
genes <- unique(bed_file$gene[start_idx:end_idx])
if (length(genes) != 1) next
gene <- genes[1]
if (is.na(gene) || gene == "") next
gene_rows <- which(bed_file$gene == gene & bed_file[[chr_col]] == bed_file[[chr_col]][start_idx])
gene_exons <- bed_file$exon_number[gene_rows]
gene_exons <- gene_exons[!is.na(gene_exons)]
if (!length(gene_exons)) next
first_exon <- min(gene_exons)
last_exon  <- max(gene_exons)
# A call that covers EVERY exon of the gene is a whole-gene event, not a
# boundary artefact. Without this check every call in a gene of one or two
# exons was "terminal-only" (all its exons are first/last exons), so
# remove_terminal_only silently discarded whole-gene deletions/duplications.
covers_whole_gene <- all(gene_exons %in% affected_exons)
if (!covers_whole_gene && all(affected_exons %in% c(first_exon, last_exon))) {
  cnv_calls$terminal_only[i] <- TRUE
}
}
cnv_calls
}
#' Classify chromosome names as autosome / X / Y / other
#'
#' Matches exactly what \code{run_cnv_calling()} analyses: autosomes are
#' \code{1}-\code{22} (with or without a \code{chr} prefix), gonosomes are
#' \code{X} and \code{Y}; everything else (chrM, alt/decoy contigs...) is
#' \code{"other"} and is never analysed in any \code{modechrom}.
#' @param chr Character vector of chromosome names.
#' @return Character vector: \code{"autosome"}, \code{"X"}, \code{"Y"} or \code{"other"}.
#' @noRd
chromosome_class <- function(chr) {
  s   <- toupper(sub("^chr", "", trimws(as.character(chr)), ignore.case = TRUE))
  out <- rep("other", length(s))
  out[grepl("^([1-9]|1[0-9]|2[0-2])$", s)] <- "autosome"
  out[s == "X"] <- "X"
  out[s == "Y"] <- "Y"
  out
}

#' Which chromosome classes are in the target BED, and which this run analyses
#'
#' @param chromosomes Chromosome name of every target.
#' @param modechrom \code{"A"} (autosomes), \code{"X"} or \code{"Y"}; anything
#'   else is treated as \code{"A"}, like \code{run_cnv_calling()} does.
#' @return Data frame with columns \code{class}, \code{label}, \code{targets},
#'   \code{present} and \code{analysed}.
#' @noRd
chromosome_scope <- function(chromosomes, modechrom = "A") {
  modechrom <- toupper(as.character(modechrom %||% "A")[1])
  cls       <- chromosome_class(chromosomes)
  keys      <- c("autosome", "X", "Y", "other")
  analysed_class <- switch(modechrom, X = "X", Y = "Y", "autosome")
  scope <- data.frame(
    class   = keys,
    label   = c("Autosomes (chr1-22)", "chrX", "chrY", "Other contigs"),
    targets = vapply(keys, function(k) sum(cls == k), integer(1), USE.NAMES = FALSE),
    stringsAsFactors = FALSE
  )
  scope$present  <- scope$targets > 0L
  scope$analysed <- scope$class == analysed_class & scope$present
  scope
}

#' Warnings about the chromosome scope of a run
#'
#' ECHO calls CNVs on ONE chromosome class per run (\code{modechrom}): so every
#' run leaves the autosomes or the gonosomes (chrX/chrY) unanalysed, and a BED
#' file may also simply lack a class. Both situations are reported here so they
#' show up in the QC metrics table and in the HTML report.
#'
#' @param scope Output of \code{chromosome_scope()}.
#' @param modechrom The run's \code{modechrom}.
#' @return Data frame with columns \code{Type}, \code{Details}, \code{Severity}
#'   (zero rows when there is nothing to report).
#' @noRd
chromosome_scope_warnings <- function(scope, modechrom = "A") {
  modechrom <- toupper(as.character(modechrom %||% "A")[1])
  if (!modechrom %in% c("X", "Y")) modechrom <- "A"
  rows <- list()
  add  <- function(type, details, severity) {
    rows[[length(rows) + 1L]] <<- data.frame(Type = type, Details = details,
                                             Severity = severity, stringsAsFactors = FALSE)
  }
  get  <- function(k) scope[scope$class == k, , drop = FALSE]
  auto <- get("autosome"); x <- get("X"); y <- get("Y"); other <- get("other")

  # 1. Not in the BED at all.
  if (!auto$present) add("Missing Chromosome", "Autosomes (chr1-22) not present in BED file", "error")
  if (!x$present)    add("Missing Chromosome", "chrX not present in BED file", "warning")
  if (!y$present)    add("Missing Chromosome", "chrY not present in BED file", "warning")

  # 2. In the BED but not analysed by this run.
  mode_label <- switch(modechrom, A = "autosomes", X = "chrX", Y = "chrY")
  hint <- function(k) switch(k,
    autosome = "modechrom = 'A'",
    X = "modechrom = 'X' (needs a sample_table with gender)",
    Y = "modechrom = 'Y' (male samples only; needs a sample_table with gender)")
  skipped <- scope[scope$class %in% c("autosome", "X", "Y") & scope$present & !scope$analysed, ,
                   drop = FALSE]
  if ("autosome" %in% skipped$class) {
    add("Analysis scope",
        sprintf(paste0("Autosomes (chr1-22) NOT analysed in this run (%d targets; modechrom = '%s' ",
                       "restricts calling to %s): no autosomal CNV is reported. Re-run with %s to analyse them."),
                auto$targets, modechrom, mode_label, hint("autosome")),
        "warning")
  }
  gon <- skipped[skipped$class %in% c("X", "Y"), , drop = FALSE]
  if (nrow(gon)) {
    add("Analysis scope",
        sprintf(paste0("Gonosomes NOT analysed in this run (modechrom = '%s'): %s. No CNV can be detected ",
                       "on them; re-run with %s."),
                modechrom,
                paste(sprintf("%s (%d targets)", gon$label, gon$targets), collapse = " and "),
                paste(vapply(gon$class, hint, character(1)), collapse = " and ")),
        "warning")
  }
  if (other$present) {
    add("Analysis scope",
        sprintf(paste0("%d target(s) lie on non-canonical contigs (e.g. chrM, alt/decoy sequences); ",
                       "they are never analysed in any modechrom."), other$targets),
        "info")
  }
  if (!length(rows)) {
    return(data.frame(Type = character(), Details = character(), Severity = character(),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}
