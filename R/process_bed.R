#' Process BED File for CNV Analysis (RefSeq Select version)
#'
#' REGEN mode now uses a curated RefSeq Select exon BED/GTF-style file
#' instead of TxDb inference (which is not clinically stable).
#'
#' @export
process_bed_file <- function(input_bed, output_bed, bed_process = "STANDARD",
                             refseqgene = NULL, transcripts_file = NULL,
                             unknown_gene = FALSE, gene_list_restrict = NULL,
                             chr_list_restrict = NULL,
                             exon_sep = NULL, kmer = NULL, customexon = FALSE,
                             list_genes = NULL, genes_file = NULL,
                             genome_version = NULL, txdb = NULL,
                             gene_field_index = 1) {

  # ----------------------------
  # Load input BED
  # ----------------------------
  input_df <- utils::read.table(input_bed, sep = "\t", header = FALSE,
                                stringsAsFactors = FALSE)

  if (ncol(input_df) < 3) {
    stop("Input BED must have at least 3 columns.")
  }

  colnames(input_df)[1:3] <- c("Chr", "Start", "End")

  if (ncol(input_df) >= 4) {
    colnames(input_df)[4] <- "Gene"
  } else {
    input_df$Gene <- NA_character_
  }

  input_df$Chr <- trimws(as.character(input_df$Chr))

  align_to_input_chr <- function(chr_vec) {
    normalize_chromosome_vec(as.character(chr_vec), input_df$Chr)
  }

  # ----------------------------
  # Helper: largest overlap selection
  # ----------------------------
  select_best_hits <- function(gr1, gr2, hits) {

    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)

    ov <- suppressWarnings(GenomicRanges::pintersect(gr1[qh], gr2[sh]))
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

    # -----------------------------
    # Load TxDb (hg38 default)
    # -----------------------------
    if (is.null(genome_version)) {
      stop("REGEN mode requires genome_version = 'hg19' or 'hg38'")
    }

    if (genome_version == "hg19") {
      if (!requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE))
        stop("Install TxDb.Hsapiens.UCSC.hg19.knownGene")
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    } else if (genome_version == "hg38") {
      if (!requireNamespace("TxDb.Hsapiens.UCSC.hg38.knownGene", quietly = TRUE))
        stop("Install TxDb.Hsapiens.UCSC.hg38.knownGene")
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    } else {
      stop("genome_version must be 'hg19' or 'hg38'")
    }

    # -----------------------------
    # Extract exon model
    # -----------------------------
    ex <- GenomicFeatures::exons(txdb, columns = c("tx_name", "GENEID"))

    # S4-Aware Unpacking: Rapidly collapse lists into strings, 
    # then take the first item to prevent data.frame mismatch crashes
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
    
    # Remove rows entirely missing a GeneID
    ref_df <- ref_df[ref_df$GeneID != "", ]

    # -----------------------------
    # Gene symbol mapping (safe)
    # -----------------------------
    if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      ref_df$Gene <- suppressMessages(AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = ref_df$GeneID,
        keytype = "ENTREZID",
        column = "SYMBOL",
        multiVals = "first"
      ))
    } else {
      ref_df$Gene <- ref_df$GeneID
    }

    # remove missing genes
    ref_df <- ref_df[!is.na(ref_df$Gene), ]
    ref_df$Chr <- align_to_input_chr(ref_df$Chr)

    # -----------------------------
    # RefSeq-like filtering
    # -----------------------------
    # keep NM_ transcripts if present
    if (any(grepl("^NM_", ref_df$Transcript))) {
      ref_df <- ref_df[grepl("^NM_", ref_df$Transcript), ]
    }

    # -----------------------------
    # Transcript ordering (strand-aware stable exon numbering)
    # -----------------------------
    split_ref <- split(ref_df, ref_df$Transcript)
    processed_list <- lapply(split_ref, function(sub_df) {
      if (sub_df$Strand[1] == "-") {
        # Negative strand transcripts transcribe backwards
        sub_df <- sub_df[order(sub_df$Start, decreasing = TRUE), ]
      } else {
        # Positive strand transcripts transcribe forwards
        sub_df <- sub_df[order(sub_df$Start, decreasing = FALSE), ]
      }
      sub_df$Exon <- seq_len(nrow(sub_df))
      return(sub_df)
    })
    
    ref_df <- do.call(rbind, processed_list)
    rownames(ref_df) <- NULL

    # -----------------------------
    # GRanges conversion
    # -----------------------------
    ref_gr <- GenomicRanges::GRanges(
      seqnames = ref_df$Chr,
      ranges = IRanges::IRanges(start = ref_df$Start,
                                end = ref_df$End),
      Gene = ref_df$Gene,
      Transcript = ref_df$Transcript,
      Exon = ref_df$Exon
    )

    bed_gr <- GenomicRanges::GRanges(
      seqnames = input_df$Chr,
      ranges = IRanges::IRanges(start = input_df$Start + 1,
                                end = input_df$End)
    )

    hits <- suppressWarnings(GenomicRanges::findOverlaps(bed_gr, ref_gr))

    if (length(hits) == 0) {
      stop("[ERROR] No overlaps between BED and RefSeq-like transcripts")
    }

    # -----------------------------
    # best-overlap selection
    # -----------------------------
    qh <- S4Vectors::queryHits(hits)
    sh <- S4Vectors::subjectHits(hits)

    ov <- suppressWarnings(GenomicRanges::pintersect(bed_gr[qh], ref_gr[sh]))
    w <- IRanges::width(ov)

    hit_df <- data.frame(q = qh, s = sh, w = w)
    hit_df <- hit_df[order(hit_df$q, -hit_df$w), ]
    hit_df <- hit_df[!duplicated(hit_df$q), ]

    qh <- hit_df$q
    sh <- hit_df$s

    # -----------------------------
    # FINAL OUTPUT (REGEN)
    # -----------------------------
    df <- data.frame(
      Chr = as.character(GenomicRanges::seqnames(bed_gr[qh])),
      Start = GenomicRanges::start(bed_gr[qh]) - 1,
      End = GenomicRanges::end(bed_gr[qh]),
      Gene = ref_gr$Gene[sh],
      Transcript = ref_gr$Transcript[sh],
      Exon = ref_gr$Exon[sh],
      stringsAsFactors = FALSE
    )

    df$Custom.Exon <- df$Exon

  } else if (bed_process == "STANDARD") {
  # ----------------------------
  # STANDARD MODE
  # ----------------------------

    panel_files <- NULL

    if (!is.null(list_genes) && file.exists(list_genes)) {
      panel_files <- readLines(list_genes)
    } else if (!is.null(genes_file) && file.exists(genes_file)) {
      panel_files <- readLines(genes_file)
    }

    if (!is.null(panel_files) && length(panel_files) > 0) {

      all_panels <- data.frame()

      for (f in panel_files) {
        if (!file.exists(f)) next

        tmp <- utils::read.table(f, sep = "\t", header = FALSE,
                                 stringsAsFactors = FALSE)

        if (ncol(tmp) < 4) next

        tmp <- tmp[, 1:4]
        colnames(tmp) <- c("Chr", "Start", "End", "Gene")

        all_panels <- rbind(all_panels, tmp)
      }

      all_panels <- unique(all_panels)
      all_panels$Chr <- align_to_input_chr(all_panels$Chr)

      bed_gr <- GenomicRanges::GRanges(
        seqnames = input_df$Chr,
        ranges = IRanges::IRanges(start = input_df$Start + 1,
                                  end = input_df$End)
      )

      panel_gr <- GenomicRanges::GRanges(
        seqnames = all_panels$Chr,
        ranges = IRanges::IRanges(start = all_panels$Start + 1,
                                  end = all_panels$End),
        Gene = all_panels$Gene
      )

      hits <- suppressWarnings(GenomicRanges::findOverlaps(bed_gr, panel_gr))

      sel <- select_best_hits(bed_gr, panel_gr, hits)
      qh <- sel$q
      sh <- sel$s

      df <- data.frame(
        Chr = as.character(GenomicRanges::seqnames(bed_gr[qh])),
        Start = GenomicRanges::start(bed_gr[qh]) - 1,
        End = GenomicRanges::end(bed_gr[qh]),
        Gene = panel_gr$Gene[sh],
        Custom.Exon = NA
      )

    } else {

      df <- input_df[, 1:4]
      colnames(df) <- c("Chr", "Start", "End", "Gene")
      df$Custom.Exon <- NA

    }
  }

  # ----------------------------
  # POST-PROCESSING
  # ----------------------------
  df$Gene <- as.character(df$Gene)

  # Extract pure Gene Symbol if overlapping regions yield comma-separated combinations
  df$Gene <- vapply(strsplit(df$Gene, ",", fixed = TRUE), function(x) {
    if (length(x) > 0) x[1] else NA_character_
  }, character(1))

  if (!unknown_gene) {
    df <- df[!is.na(df$Gene) & df$Gene != "" & df$Gene != ".", ]
  } else {
    df$Gene[is.na(df$Gene) | df$Gene == ""] <- "Unknown"
  }

  df <- unique(df)

  # genomic sort
  chrom_base <- c(as.character(1:22), "X", "Y", "M")
  chrom_levels <- if (any(grepl("^chr", df$Chr))) paste0("chr", chrom_base) else chrom_base
  df$Chr <- factor(df$Chr, levels = chrom_levels)

  df <- df[order(df$Chr, df$Start), ]
  df$Chr <- as.character(df$Chr)

  # output
  if ("Custom.Exon" %in% colnames(df)) {
    out_df <- df[, c("Chr", "Start", "End", "Gene", "Custom.Exon")]
  } else {
    out_df <- df[, c("Chr", "Start", "End", "Gene")]
  }

  utils::write.table(out_df, file = output_bed,
                     sep = "\t", row.names = FALSE,
                     col.names = FALSE, quote = FALSE)

  message("[INFO] BED written: ", output_bed)

  invisible(output_bed)
}