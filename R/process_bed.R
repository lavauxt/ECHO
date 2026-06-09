#' Process BED File for CNV Analysis
#'
#' @param input_bed Path to input BED file.
#' @param output_bed Path to output processed BED file.
#' @param bed_process Mode: "STANDARD", "REGEN", or "NO".
#' @param bed_zero_based Logical: is the input BED 0‑based? Default TRUE.
#' @param unknown_gene Keep intervals without gene name? Default FALSE.
#' @param exon_sep Separator for splitting gene names (e.g., "_"). Default "_".
#' @param customexon Add a Custom.Exon column containing exon numbers. Default FALSE.
#' @param auto_exon_number If TRUE, assign sequential exon numbers per gene based on genomic order.
#'        If FALSE, try to extract from name; if none found, use NA. Default TRUE.
#' @param panel_files Character vector of BED file paths (one per line) or a single file containing paths.
#' @param genome_version "hg19" or "hg38" (REGEN mode).
#' @param gene_field_index 1‑based index of the gene name field after splitting.
#'        If NULL, auto‑detect (3rd for NM‑prefixed, 1st otherwise).
#' @param ... other parameters for compatibility.
#'
#' @export
process_bed_file <- function(input_bed, output_bed, bed_process = "STANDARD",
                             bed_zero_based = TRUE,
                             refseqgene = NULL, transcripts_file = NULL,
                             unknown_gene = FALSE, gene_list_restrict = NULL,
                             chr_list_restrict = NULL,
                             exon_sep = "_", kmer = NULL, customexon = FALSE,
                             auto_exon_number = TRUE,
                             list_genes = NULL, genes_file = NULL,
                             panel_files = NULL,
                             genome_version = NULL, txdb = NULL,
                             gene_field_index = NULL) {

  # ----------------------------
  # Helper: extract gene name and optional exon number
  # ----------------------------
  parse_bed_name <- function(name_vec, exon_sep = "_", gene_field_index = NULL, auto_exon = TRUE) {
    if (is.null(exon_sep) || exon_sep == "") exon_sep <- "_"
    
    gene_names <- character(length(name_vec))
    exon_numbers <- integer(length(name_vec))
    
    for (i in seq_along(name_vec)) {
      g <- as.character(name_vec[i])
      if (is.na(g) || g == "" || g == ".") {
        gene_names[i] <- NA_character_
        exon_numbers[i] <- NA_integer_
        next
      }
      
      # Remove any trailing comma-separated extras (keep first)
      g <- strsplit(g, ",")[[1]][1]
      
      # ----- Extract exon number from name (if present) -----
      exon_num <- NA_integer_
      # Pattern: _ex(\d+)_ , _ex(\d+)\( , ex(\d+)_ , ex(\d+)\(
      m <- regexpr("_ex([0-9]+)", g)
      if (m == -1) m <- regexpr("ex([0-9]+)", g)
      if (m != -1) {
        exon_str <- regmatches(g, m)
        exon_num <- as.integer(gsub("[^0-9]", "", exon_str))
      }
      # If auto_exon is FALSE but we didn't find a number, keep NA; we will fill later if auto_exon is TRUE.
      # However, we cannot fill here because we need per-gene ordering; that will be done after all parsing.
      exon_numbers[i] <- exon_num
      
      # ----- Extract gene name -----
      parts <- strsplit(g, exon_sep, fixed = TRUE)[[1]]
      if (!is.null(gene_field_index) && is.numeric(gene_field_index) && gene_field_index >= 1) {
        if (length(parts) >= gene_field_index) {
          gene <- parts[gene_field_index]
        } else {
          gene <- parts[1]
        }
      } else {
        # Auto-detect
        if (grepl("^(NM_|NR_|XM_)", g)) {
          if (length(parts) >= 3) {
            gene <- parts[3]
          } else {
            gene <- parts[1]
          }
        } else {
          gene <- parts[1]
        }
      }
      # Remove parentheses content
      gene <- gsub("\\(.*$", "", gene)
      gene_names[i] <- gene
    }
    list(gene = gene_names, exon = exon_numbers)
  }

  # ----------------------------
  # Helper: assign sequential exon numbers per gene based on genomic order
  # ----------------------------
  assign_sequential_exons <- function(df, start_col = "Start", gene_col = "Gene") {
    chrom_order <- c(paste0("chr", c(1:22, "X", "Y", "M")))
    if (!any(grepl("^chr", df$Chr))) chrom_order <- sub("^chr", "", chrom_order)

    df$.orig_order <- seq_len(nrow(df))
    df$.chr_order <- factor(df$Chr, levels = chrom_order)
    df$ExonNum <- NA_integer_

    valid_gene <- !is.na(df[[gene_col]]) & df[[gene_col]] != "" & df[[gene_col]] != "."
    gene_groups <- split(which(valid_gene), df[[gene_col]][valid_gene])

    for (idx in gene_groups) {
      idx_ordered <- idx[order(df$.chr_order[idx], df[[start_col]][idx], df$End[idx], na.last = TRUE)]
      df$ExonNum[idx_ordered] <- seq_along(idx_ordered)
    }

    df <- df[order(df$.chr_order, df[[start_col]], df$End, na.last = TRUE), ]
    df$.orig_order <- NULL
    df$.chr_order <- NULL
    df
  }

  # ----------------------------
  # Load input BED
  # ----------------------------
  input_df <- utils::read.table(input_bed, sep = "\t", header = FALSE,
                                stringsAsFactors = FALSE)
  if (ncol(input_df) < 3) stop("Input BED must have at least 3 columns.")
  colnames(input_df)[1:3] <- c("Chr", "Start", "End")
  if (ncol(input_df) >= 4) {
    colnames(input_df)[4] <- "OriginalName"
  } else {
    input_df$OriginalName <- NA_character_
  }
  
  # Parse gene names and try to extract exon numbers
  parsed <- parse_bed_name(input_df$OriginalName, exon_sep, gene_field_index, auto_exon_number)
  input_df$Gene <- parsed$gene
  input_df$ExtractedExon <- parsed$exon
  
  # If auto_exon_number is TRUE, we will assign sequential numbers later (overriding any extracted)
  # If FALSE, we keep extracted numbers (which may be NA)
  if (!auto_exon_number) {
    input_df$ExonNum <- input_df$ExtractedExon
  }

  # Helper for overlap selection
  select_best_hits <- function(gr1, gr2, hits) {
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    ov <- GenomicRanges::pintersect(gr1[qh], gr2[sh])
    w <- IRanges::width(ov)
    df <- data.frame(q = qh, s = sh, w = w)
    df <- df[order(df$q, -df$w), ]
    df <- df[!duplicated(df$q), ]
    list(q = df$q, s = df$s)
  }

  # ----------------------------
  # REGEN MODE (RefSeq Select)
  # ----------------------------
  if (bed_process == "REGEN") {
    message("[INFO] REGEN mode: using internal RefSeq Select-like DB")
    if (is.null(genome_version)) stop("REGEN mode requires genome_version = 'hg19' or 'hg38'")
    if (genome_version == "hg19") {
      if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE))
        stop("Install TxDb.Hsapiens.UCSC.hg19.knownGene")
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    } else if (genome_version == "hg38") {
      if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE))
        stop("Install TxDb.Hsapiens.UCSC.hg38.knownGene")
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    } else stop("genome_version must be 'hg19' or 'hg38'")

    ex <- GenomicFeatures::exons(txdb, columns = c("tx_name", "GENEID"))
    nm_str   <- S4Vectors::unstrsplit(ex$tx_name, sep = ",")
    gene_str <- S4Vectors::unstrsplit(ex$GENEID, sep = ",")
    ref_df <- data.frame(
      Chr = as.character(GenomicRanges::seqnames(ex)),
      Start = GenomicRanges::start(ex),
      End = GenomicRanges::end(ex),
      Strand = as.character(GenomicRanges::strand(ex)),
      Transcript = sub(",.*", "", nm_str),
      GeneID = sub(",.*", "", gene_str),
      stringsAsFactors = FALSE
    )
    ref_df <- ref_df[ref_df$GeneID != "", ]
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      ref_df$Gene <- suppressMessages(AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = ref_df$GeneID,
        keytype = "ENTREZID",
        column = "SYMBOL",
        multiVals = "first"
      ))
    } else ref_df$Gene <- ref_df$GeneID
    ref_df <- ref_df[!is.na(ref_df$Gene), ]
    # Clean gene names using same parser (though usually already clean)
    ref_parsed <- parse_bed_name(ref_df$Gene, exon_sep, gene_field_index, auto_exon_number)
    ref_df$Gene <- ref_parsed$gene
    if (any(grepl("^NM_", ref_df$Transcript))) ref_df <- ref_df[grepl("^NM_", ref_df$Transcript), ]
    
    split_ref <- split(ref_df, ref_df$Transcript)
    processed_list <- lapply(split_ref, function(sub_df) {
      if (sub_df$Strand[1] == "-") sub_df <- sub_df[order(sub_df$Start, decreasing = TRUE), ]
      else sub_df <- sub_df[order(sub_df$Start, decreasing = FALSE), ]
      sub_df$Exon <- seq_len(nrow(sub_df))
      sub_df
    })
    ref_df <- do.call(rbind, processed_list)
    rownames(ref_df) <- NULL

    ref_gr <- GenomicRanges::GRanges(
      seqnames = ref_df$Chr,
      ranges = IRanges::IRanges(start = ref_df$Start, end = ref_df$End),
      Gene = ref_df$Gene,
      Transcript = ref_df$Transcript,
      Exon = ref_df$Exon
    )

    input_start_1based <- if (bed_zero_based) input_df$Start + 1 else input_df$Start
    bed_gr <- GenomicRanges::GRanges(
      seqnames = input_df$Chr,
      ranges = IRanges::IRanges(start = input_start_1based, end = input_df$End)
    )
    hits <- GenomicRanges::findOverlaps(bed_gr, ref_gr)
    if (length(hits) == 0) stop("[ERROR] No overlaps between BED and RefSeq-like transcripts")
    
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)
    ov <- GenomicRanges::pintersect(bed_gr[qh], ref_gr[sh])
    w <- IRanges::width(ov)
    hit_df <- data.frame(q = qh, s = sh, w = w)
    hit_df <- hit_df[order(hit_df$q, -hit_df$w), ]
    hit_df <- hit_df[!duplicated(hit_df$q), ]
    qh <- hit_df$q
    sh <- hit_df$s

    output_start <- GenomicRanges::start(bed_gr[qh])
    if (bed_zero_based) output_start <- output_start - 1
    df <- data.frame(
      Chr = as.character(GenomicRanges::seqnames(bed_gr[qh])),
      Start = output_start,
      End = GenomicRanges::end(bed_gr[qh]),
      Gene = ref_gr$Gene[sh],
      Transcript = ref_gr$Transcript[sh],
      Exon = ref_gr$Exon[sh],
      stringsAsFactors = FALSE
    )
    df$Custom.Exon <- df$Exon
    # REGEN mode uses RefSeq exon numbers, no auto-renumbering needed
    # But if user wants to force auto_exon_number, we could still apply
    if (auto_exon_number && !customexon) {
      # Not applicable because Custom.Exon won't be written; but we could renumber anyway.
      # For simplicity, we skip.
    }

  } else if (bed_process == "STANDARD") {
    # ----------------------------
    # STANDARD MODE
    # ----------------------------
    panel_bed_paths <- NULL
    if (!is.null(list_genes) && file.exists(list_genes)) panel_bed_paths <- readLines(list_genes)
    else if (!is.null(genes_file) && file.exists(genes_file)) panel_bed_paths <- readLines(genes_file)
    else if (!is.null(panel_files)) {
      if (length(panel_files) == 1 && file.exists(panel_files[1])) panel_bed_paths <- readLines(panel_files[1])
      else panel_bed_paths <- panel_files
    }

    if (!is.null(panel_bed_paths) && length(panel_bed_paths) > 0) {
      all_panels <- data.frame()
      for (f in panel_bed_paths) {
        if (!file.exists(f)) next
        tmp <- utils::read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
        if (ncol(tmp) < 4) next
        tmp <- tmp[, 1:4]
        colnames(tmp) <- c("Chr", "Start", "End", "OriginalName")
        all_panels <- rbind(all_panels, tmp)
      }
      all_panels <- unique(all_panels)
      panel_parsed <- parse_bed_name(all_panels$OriginalName, exon_sep, gene_field_index, auto_exon_number)
      all_panels$Gene <- panel_parsed$gene
      all_panels$ExtractedExon <- panel_parsed$exon
      if (!auto_exon_number) all_panels$ExonNum <- all_panels$ExtractedExon else all_panels$ExonNum <- NA
      
      # Overlap with input BED intervals
      input_start_1based <- if (bed_zero_based) input_df$Start + 1 else input_df$Start
      bed_gr <- GenomicRanges::GRanges(
        seqnames = input_df$Chr,
        ranges = IRanges::IRanges(start = input_start_1based, end = input_df$End)
      )
      panel_gr <- GenomicRanges::GRanges(
        seqnames = all_panels$Chr,
        ranges = IRanges::IRanges(start = all_panels$Start + 1, end = all_panels$End),
        Gene = all_panels$Gene,
        ExonNum = all_panels$ExonNum
      )
      hits <- GenomicRanges::findOverlaps(bed_gr, panel_gr)
      sel <- select_best_hits(bed_gr, panel_gr, hits)
      qh <- sel$q; sh <- sel$s
      output_start <- GenomicRanges::start(bed_gr[qh])
      if (bed_zero_based) output_start <- output_start - 1
      df <- data.frame(
        Chr = as.character(GenomicRanges::seqnames(bed_gr[qh])),
        Start = output_start,
        End = GenomicRanges::end(bed_gr[qh]),
        Gene = panel_gr$Gene[sh],
        ExonNum = panel_gr$ExonNum[sh],
        stringsAsFactors = FALSE
      )
      # For panel mode, we might also want to auto-number exons if requested
      if (auto_exon_number) {
        # Need to reassign ExonNum sequentially per gene based on order in df (which is already overlapped)
        # But careful: df may not be sorted; we sort and assign
        df <- assign_sequential_exons(df, start_col = "Start", gene_col = "Gene")
      }
      df$Custom.Exon <- df$ExonNum
    } else {
      # No panel file: use original BED with parsed gene and exon numbers
      df <- data.frame(
        Chr = input_df$Chr,
        Start = if (bed_zero_based) input_df$Start else input_df$Start - 1,
        End = input_df$End,
        Gene = input_df$Gene,
        ExonNum = if (auto_exon_number) NA else input_df$ExtractedExon,
        stringsAsFactors = FALSE
      )
      if (auto_exon_number) {
        # Assign sequential exon numbers per gene based on genomic order
        df <- assign_sequential_exons(df, start_col = "Start", gene_col = "Gene")
      }
      df$Custom.Exon <- df$ExonNum
    }
  } else {
    stop("bed_process must be 'STANDARD', 'REGEN', or 'NO'")
  }

  # ----------------------------
  # POST-PROCESSING (common)
  # ----------------------------
  df$Gene <- as.character(df$Gene)
  # Remove any remaining commas (keep first)
  df$Gene <- vapply(strsplit(df$Gene, ",", fixed = TRUE), function(x) x[1], character(1))

  if (!unknown_gene) {
    df <- df[!is.na(df$Gene) & df$Gene != "" & df$Gene != ".", ]
  } else {
    df$Gene[is.na(df$Gene) | df$Gene == ""] <- "Unknown"
  }
  df <- unique(df)

  # Genomic sort (again, to ensure final order)
  chrom_base <- c(as.character(1:22), "X", "Y", "M")
  if (any(grepl("^chr", df$Chr))) chrom_levels <- paste0("chr", chrom_base) else chrom_levels <- chrom_base
  df$Chr <- factor(df$Chr, levels = chrom_levels)
  df <- df[order(df$Chr, df$Start), ]
  df$Chr <- as.character(df$Chr)

  # Output columns. The 5th column is always the per-gene exon number when available.
  if ("Custom.Exon" %in% colnames(df)) {
    out_df <- df[, c("Chr", "Start", "End", "Gene", "Custom.Exon")]
  } else if ("ExonNum" %in% colnames(df)) {
    out_df <- df[, c("Chr", "Start", "End", "Gene", "ExonNum")]
  } else {
    out_df <- df[, c("Chr", "Start", "End", "Gene")]
  }

  utils::write.table(out_df, file = output_bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  message("[INFO] BED written: ", output_bed)
  invisible(output_bed)
}