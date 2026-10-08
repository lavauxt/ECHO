select_best_hits <- function(gr1, gr2, hits) {
  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  overlap_width <- IRanges::width(GenomicRanges::pintersect(gr1[qh], gr2[sh]))
  hit_data <- data.frame(q = qh, s = sh, width = overlap_width)
  hit_data <- hit_data[order(hit_data$q, -hit_data$width, hit_data$s), , drop = FALSE]
  hit_data <- hit_data[!duplicated(hit_data$q), , drop = FALSE]
  list(q = hit_data$q, s = hit_data$s)
}

annotate_bed_regen <- function(input_df, ref_df, bed_zero_based = TRUE) {
  input_start <- suppressWarnings(as.numeric(input_df$Start))
  input_end <- suppressWarnings(as.numeric(input_df$End))
  if (any(!is.finite(input_start)) || any(!is.finite(input_end)) ||
      any(input_start != floor(input_start)) || any(input_end != floor(input_end))) {
    stop("[ERROR] BED start/end coordinates must be numeric integers for REGEN.")
  }

  query_start <- if (bed_zero_based) input_start + 1 else input_start
  if (any(query_start < 1 | input_end < query_start)) {
    stop("[ERROR] BED contains invalid coordinates for the selected bed_zero_based setting.")
  }

  ref_chroms <- unique(ref_df$Chr)
  map_chromosome <- function(chr) {
    if (chr %in% ref_chroms) return(chr)
    alternative <- if (grepl("^chr", chr)) sub("^chr", "", chr) else paste0("chr", chr)
    if (alternative %in% ref_chroms) return(alternative)
    NA_character_
  }
  input_chroms <- vapply(as.character(input_df$Chr), map_chromosome, character(1))
  query_rows <- which(!is.na(input_chroms))

  query <- GenomicRanges::GRanges(
    seqnames = input_chroms[query_rows],
    ranges = IRanges::IRanges(start = query_start[query_rows], end = input_end[query_rows])
  )
  subject <- GenomicRanges::GRanges(
    seqnames = ref_df$Chr,
    ranges = IRanges::IRanges(start = ref_df$Start, end = ref_df$End)
  )
  hits <- GenomicRanges::findOverlaps(query, subject, ignore.strand = TRUE)
  best <- if (length(hits)) select_best_hits(query, subject, hits) else
    list(q = integer(), s = integer())

  result <- data.frame(
    chromosome = as.character(input_df$Chr),
    start = if (bed_zero_based) input_start else input_start - 1,
    end = input_end,
    gene = "Unknown",
    stringsAsFactors = FALSE
  )
  if (length(best$q)) result$gene[query_rows[best$q]] <- ref_df$Gene[best$s]

  result <- assign_exon_numbers_per_gene(result)
  unknown <- is.na(result$gene) | result$gene %in% c("", ".", "Unknown")
  result$exon_number[unknown] <- NA_integer_
  result
}
