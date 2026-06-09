#' Process BED File for CNV Analysis
#'
#' @param input_bed Path to input BED file.
#' @param output_bed Path to output processed BED file.
#' @param bed_process Mode: "STANDARD", "REGEN", or "NO".
#' @param bed_zero_based Logical: is the input BED 0‑based? Default TRUE.
#' @param unknown_gene Keep intervals without gene name? Default FALSE.
#' @param exon_sep Character string or list/vector of characters to split gene names (e.g., "[_()]" or c("_", "(")). Default "_".
#' @param gene_name_collapse Character token used to rejoin kept parts. Default "_".
#' @param customexon Add a Custom.Exon column containing exon numbers. Default FALSE.
#' @param auto_exon_number If TRUE, assign sequential exon numbers per gene based on genomic order.
#'        If FALSE, check region_numbering_mode configuration. Default TRUE.
#' @param region_numbering_mode If auto_exon_number is FALSE, specifies numbering type:
#'        "bed_text" (extract from string pattern) or "file_order" (sequential per gene based on local file occurrence).
#' @param gene_name_keep Specifies parts to keep after splitting by exon_sep (e.g., "1", "1-2", "1,3").
#' @param panel_files Character vector of BED file paths (one per line) or a single file containing paths.
#' @param genome_version "hg19" or "hg38" (REGEN mode).
#' @param gene_field_index 1‑based index of the gene name field after splitting (legacy fallback for gene_name_keep).
#' @param ... other parameters for compatibility.
#'
#' @export
process_bed_file <- function(input_bed, output_bed, bed_process = "STANDARD",
                             bed_zero_based = TRUE,
                             refseqgene = NULL, transcripts_file = NULL,
                             unknown_gene = FALSE, gene_list_restrict = NULL,
                             chr_list_restrict = NULL,
                             exon_sep = "_", gene_name_collapse = "_", kmer = NULL, customexon = FALSE,
                             auto_exon_number = TRUE, region_numbering_mode = "bed_text",
                             gene_name_keep = NULL,
                             list_genes = NULL, genes_file = NULL,
                             panel_files = NULL,
                             genome_version = NULL, txdb = NULL,
                             gene_field_index = NULL) {

  # ----------------------------
  # Helper: Parse keep indices (supports "1", "1-2", "1,3")
  # ----------------------------
  parse_keep_indices <- function(keep_str, max_len) {
    if (is.null(keep_str) || keep_str == "") return(NULL)
    if (is.numeric(keep_str)) return(keep_str)
    keep_str <- as.character(keep_str)
    
    if (grepl("-", keep_str)) {
      parts <- as.numeric(strsplit(keep_str, "-")[[1]])
      if (length(parts) == 2 && !any(is.na(parts))) {
        return(seq(parts[1], min(parts[2], max_len)))
      }
    }
    if (grepl(",", keep_str)) {
      parts <- as.numeric(strsplit(keep_str, ",")[[1]])
      return(parts[!is.na(parts)])
    }
    idx <- as.numeric(keep_str)
    if (!is.na(idx)) return(idx)
    return(NULL)
  }

  # ----------------------------
  # Helper: extract gene name and optional exon number
  # ----------------------------
  parse_bed_name <- function(name_vec, exon_sep = "_", gene_field_index = NULL, 
                             gene_name_keep = NULL, auto_exon = TRUE, gene_name_collapse = "_") {
    if (is.null(exon_sep) || length(exon_sep) == 0 || any(exon_sep == "")) exon_sep <- "_"
    if (is.null(gene_name_collapse) || gene_name_collapse == "") gene_name_collapse <- "_"
    
    # Evaluate if multi-token regex pattern compilation or literal evaluation is needed
    if (length(exon_sep) > 1) {
      escaped_seps <- vapply(exon_sep, function(x) {
        gsub("([\\\\^\\$\\.\\|\\?\\*\\+\\(\\)\\[\\{\\]\\}])", "\\\\\\1", x)
      }, character(1))
      split_pat <- paste(escaped_seps, collapse = "|")
      use_fixed <- FALSE
    } else {
      if (grepl("[\\[\\]\\(\\)\\|\\.\\*\\+\\?\\^\\$]", exon_sep, perl = TRUE)) {
        split_pat <- exon_sep
        use_fixed <- FALSE
      } else {
        split_pat <- exon_sep
        use_fixed <- TRUE
      }
    }
    
    gene_names <- character(length(name_vec))
    exon_numbers <- integer(length(name_vec))
    
    for (i in seq_along(name_vec)) {
      g <- as.character(name_vec[i])
      if (is.na(g) || g == "" || g == ".") {
        gene_names[i] <- NA_character_
        exon_numbers[i] <- NA_integer_
        next
      }
      
      g <- strsplit(g, ",")[[1]][1]
      
      # Extract exon number
      exon_num <- NA_integer_
      m <- regexpr("_ex([0-9]+)", g)
      if (m == -1) m <- regexpr("ex([0-9]+)", g)
      if (m != -1) {
        exon_str <- regmatches(g, m)
        exon_num <- as.integer(gsub("[^0-9]", "", exon_str))
      }
      exon_numbers[i] <- exon_num
      
      # -----------------------------------------------------------------
      # GENERALIZED CLEANING (No Hardcoded Genes)
      # -----------------------------------------------------------------
      # 1. Strip parenthesis and anything inside them (e.g., "(zerze)")
      g_clean <- gsub("\\s*\\(.*?\\)", "", g)
      
      # 2. Split the string using your defined separators (e.g., "_")
      parts <- strsplit(g_clean, split = split_pat, fixed = use_fixed)[[1]]
      parts <- parts[parts != ""]  
      
      # 3. Apply the index rule (gene_name_keep)
      keep_idx <- NULL
      if (!is.null(gene_name_keep)) {
        keep_idx <- parse_keep_indices(gene_name_keep, length(parts))
      } else if (!is.null(gene_field_index)) {
        keep_idx <- parse_keep_indices(gene_field_index, length(parts))
      }
      
      if (!is.null(keep_idx)) {
        valid_idx <- keep_idx[keep_idx <= length(parts) & keep_idx >= 1]
        
        if (length(valid_idx) > 0) {
          # Standard case: NM_000251_MSH2 -> grabs MSH2
          gene <- paste(parts[valid_idx], collapse = gene_name_collapse)
        } else {
          # DYNAMIC FALLBACK: If requested index is larger than available parts
          # grab the LAST part available.
          gene <- parts[length(parts)]
        }
      } else {
        # Default behavior if no keep index is defined
        if (grepl("^(NM_|NR_|XM_)", g_clean)) {
          if (length(parts) >= 3) gene <- parts[3] else gene <- parts[length(parts)]
        } else {
          gene <- parts[1]
        }
      }
      
      # 4. AUTOMATED COORDINATE SAFEGUARD (No Hardcoding)
      # If the extracted string consists entirely of digits (e.g., "58149635"),
      # it represents a structural coordinate position, not a gene descriptor.
      # Automatically revert to the primary label component.
      if (!is.na(gene) && grepl("^[0-9]+$", gene)) {
        gene <- parts[1]
      }
      
      gene_names[i] <- gene
    }
    list(gene = gene_names, exon = exon_numbers)
  }

  # ----------------------------
  # Helper: assign sequential exons per gene based on genomic order
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
  # Helper: assign region index based purely on local file occurrence order
  # ----------------------------
  assign_file_order_exons <- function(df, gene_col = "Gene") {
    df$ExonNum <- NA_integer_
    valid_gene <- !is.na(df[[gene_col]]) & df[[gene_col]] != "" & df[[gene_col]] != "."
    
    gene_counts <- list()
    for (i in seq_len(nrow(df))) {
      if (valid_gene[i]) {
        g <- df[[gene_col]][i]
        if (is.null(gene_counts[[g]])) {
          gene_counts[[g]] <- 1
        } else {
          gene_counts[[g]] <- gene_counts[[g]] + 1
        }
        df$ExonNum[i] <- gene_counts[[g]]
      }
    }
    df
  }

  # ----------------------------
  # Load input BED
  # ----------------------------
  input_df <- utils::read.table(input_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  if (ncol(input_df) < 3) stop("Input BED must have at least 3 columns.")
  colnames(input_df)[1:3] <- c("Chr", "Start", "End")
  if (ncol(input_df) >= 4) {
    colnames(input_df)[4] <- "OriginalName"
  } else {
    input_df$OriginalName <- NA_character_
  }
  
  parsed <- parse_bed_name(input_df$OriginalName, exon_sep, gene_field_index, gene_name_keep, auto_exon_number, gene_name_collapse)
  input_df$Gene <- parsed$gene
  input_df$ExtractedExon <- parsed$exon
  
  if (!auto_exon_number && region_numbering_mode == "bed_text") {
    input_df$ExonNum <- input_df$ExtractedExon
  } else {
    input_df$ExonNum <- NA_integer_
  }

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

  # REGEN MODE
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
        org.Hs.eg.db, keys = ref_df$GeneID, keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"
      ))
    } else ref_df$Gene <- ref_df$GeneID
    ref_df <- ref_df[!is.na(ref_df$Gene), ]
    ref_parsed <- parse_bed_name(ref_df$Gene, exon_sep, gene_field_index, gene_name_keep, auto_exon_number, gene_name_collapse)
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
      seqnames = ref_df$Chr, ranges = IRanges::IRanges(start = ref_df$Start, end = ref_df$End),
      Gene = ref_df$Gene, Transcript = ref_df$Transcript, Exon = ref_df$Exon
    )

    input_start_1based <- if (bed_zero_based) input_df$Start + 1 else input_df$Start
    bed_gr <- GenomicRanges::GRanges(seqnames = input_df$Chr, ranges = IRanges::IRanges(start = input_start_1based, end = input_df$End))
    hits <- GenomicRanges::findOverlaps(bed_gr, ref_gr)
    if (length(hits) == 0) stop("[ERROR] No overlaps between BED and RefSeq-like transcripts")
    
    sel <- select_best_hits(bed_gr, ref_gr, hits)
    qh <- sel$q; sh <- sel$s

    output_start <- GenomicRanges::start(bed_gr[qh])
    if (bed_zero_based) output_start <- output_start - 1
    df <- data.frame(
      Chr = as.character(GenomicRanges::seqnames(bed_gr[qh])), Start = output_start, End = GenomicRanges::end(bed_gr[qh]),
      Gene = ref_gr$Gene[sh], Transcript = ref_gr$Transcript[sh], Exon = ref_gr$Exon[sh], stringsAsFactors = FALSE
    )
    df$Custom.Exon <- df$Exon

  } else if (bed_process == "STANDARD") {
    # STANDARD MODE
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
      panel_parsed <- parse_bed_name(all_panels$OriginalName, exon_sep, gene_field_index, gene_name_keep, auto_exon_number, gene_name_collapse)
      all_panels$Gene <- panel_parsed$gene
      all_panels$ExtractedExon <- panel_parsed$exon
      if (!auto_exon_number && region_numbering_mode == "bed_text") {
        all_panels$ExonNum <- all_panels$ExtractedExon
      } else {
        all_panels$ExonNum <- NA_integer_
      }
      
      input_start_1based <- if (bed_zero_based) input_df$Start + 1 else input_df$Start
      bed_gr <- GenomicRanges::GRanges(seqnames = input_df$Chr, ranges = IRanges::IRanges(start = input_start_1based, end = input_df$End))
      panel_gr <- GenomicRanges::GRanges(
        seqnames = all_panels$Chr, ranges = IRanges::IRanges(start = all_panels$Start + 1, end = all_panels$End),
        Gene = all_panels$Gene, ExonNum = all_panels$ExonNum
      )
      hits <- GenomicRanges::findOverlaps(bed_gr, panel_gr)
      sel <- select_best_hits(bed_gr, panel_gr, hits)
      qh <- sel$q; sh <- sel$s
      output_start <- GenomicRanges::start(bed_gr[qh])
      if (bed_zero_based) output_start <- output_start - 1
      df <- data.frame(
        Chr = as.character(GenomicRanges::seqnames(bed_gr[qh])), Start = output_start, End = GenomicRanges::end(bed_gr[qh]),
        Gene = panel_gr$Gene[sh], ExonNum = panel_gr$ExonNum[sh], stringsAsFactors = FALSE
      )

      if (auto_exon_number) {
        df <- assign_sequential_exons(df, start_col = "Start", gene_col = "Gene")
      } else if (region_numbering_mode == "file_order") {
        df <- assign_file_order_exons(df, gene_col = "Gene")
      }
      df$Custom.Exon <- df$ExonNum
    } else {
      df <- data.frame(
        Chr = input_df$Chr, Start = if (bed_zero_based) input_df$Start else input_df$Start - 1, End = input_df$End,
        Gene = input_df$Gene, ExonNum = input_df$ExonNum, stringsAsFactors = FALSE
      )
      if (auto_exon_number) {
        df <- assign_sequential_exons(df, start_col = "Start", gene_col = "Gene")
      } else if (region_numbering_mode == "file_order") {
        df <- assign_file_order_exons(df, gene_col = "Gene")
      }
      df$Custom.Exon <- df$ExonNum

    }
  } else {
    stop("bed_process must be 'STANDARD', 'REGEN', or 'NO'")
  }

  # Post-Processing
  df$Gene <- as.character(df$Gene)
  df$Gene <- vapply(strsplit(df$Gene, ",", fixed = TRUE), function(x) x[1], character(1))

  if (!unknown_gene) {
    df <- df[!is.na(df$Gene) & df$Gene != "" & df$Gene != ".", ]
  } else {
    df$Gene[is.na(df$Gene) | df$Gene == ""] <- "Unknown"
  }
  df <- unique(df)

  # Final Export Genomic Sorting
  chrom_base <- c(as.character(1:22), "X", "Y", "M")
  if (any(grepl("^chr", df$Chr))) chrom_levels <- paste0("chr", chrom_base) else chrom_levels <- chrom_base
  df$Chr <- factor(df$Chr, levels = chrom_levels)
  df <- df[order(df$Chr, df$Start), ]
  df$Chr <- as.character(df$Chr)

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