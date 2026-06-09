#' Process BED File for CNV Analysis
#'
#' @param input_bed Path to input BED file.
#' @param output_bed Path to output processed BED file.
#' @param bed_process Mode: "STANDARD", "REGEN", or "NO".
#' @param bed_zero_based Logical: is the input BED 0‑based? Default TRUE.
#' @param refseqgene Not used (legacy).
#' @param transcripts_file Not used.
#' @param unknown_gene Keep intervals without gene name? Default FALSE.
#' @param gene_list_restrict Not used.
#' @param chr_list_restrict Not used.
#' @param exon_sep Not used.
#' @param kmer Not used.
#' @param customexon Not used.
#' @param list_genes Deprecated: use `panel_files` instead.
#' @param genes_file Deprecated: use `panel_files` instead.
#' @param panel_files Character vector of BED file paths (one per line) or a single file containing paths.
#' @param genome_version "hg19" or "hg38" (REGEN mode).
#' @param txdb Not used.
#' @param gene_field_index Not used.
#'
#' @export
process_bed_file <- function(input_bed, output_bed, bed_process = "STANDARD",
                             bed_zero_based = TRUE,
                             refseqgene = NULL, transcripts_file = NULL,
                             unknown_gene = FALSE, gene_list_restrict = NULL,
                             chr_list_restrict = NULL,
                             exon_sep = NULL, kmer = NULL, customexon = FALSE,
                             list_genes = NULL, genes_file = NULL,
                             panel_files = NULL,
                             genome_version = NULL, txdb = NULL,
                             gene_field_index = 1) {

    # ----------------------------
    # Load input BED using data.table
    # ----------------------------
    input_df <- data.table::fread(input_bed, header = FALSE, stringsAsFactors = FALSE)
    if (ncol(input_df) < 3) stop("Input BED must have at least 3 columns.")
    colnames(input_df)[1:3] <- c("Chr", "Start", "End")
    if (ncol(input_df) >= 4) colnames(input_df)[4] <- "Gene" else input_df$Gene <- NA_character_
    input_df$Chr <- trimws(as.character(input_df$Chr))

    align_to_input_chr <- function(chr_vec) {
        normalize_chromosome_vec(as.character(chr_vec), input_df$Chr)
    }

    # Helper: select best overlap
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
                org.Hs.eg.db::org.Hs.eg.db,
                keys = ref_df$GeneID,
                keytype = "ENTREZID",
                column = "SYMBOL",
                multiVals = "first"
            ))
        } else ref_df$Gene <- ref_df$GeneID
        ref_df <- ref_df[!is.na(ref_df$Gene), ]
        ref_df$Chr <- align_to_input_chr(ref_df$Chr)
        if (any(grepl("^NM_", ref_df$Transcript))) ref_df <- ref_df[grepl("^NM_", ref_df$Transcript), ]

        # Strand-aware exon ordering
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

        # Convert input coordinates to 1‑based for GRanges (if BED is 0‑based)
        input_start_1based <- if (bed_zero_based) input_df$Start + 1 else input_df$Start
        bed_gr <- GenomicRanges::GRanges(
            seqnames = input_df$Chr,
            ranges = IRanges::IRanges(start = input_start_1based, end = input_df$End)
        )
        hits <- suppressWarnings(GenomicRanges::findOverlaps(bed_gr, ref_gr))
        if (length(hits) == 0) stop("[ERROR] No overlaps between BED and RefSeq-like transcripts")
        sel <- select_best_hits(bed_gr, ref_gr, hits)
        qh <- sel$q; sh <- sel$s

        # Build output (1‑based internally, but write 0‑based if input was 0‑based)
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

    } else if (bed_process == "STANDARD") {
        # ----------------------------
        # STANDARD MODE with optional panel BED list
        # ----------------------------
        panel_bed_paths <- NULL
        # Backward compatibility: list_genes / genes_file
        if (!is.null(list_genes) && file.exists(list_genes)) panel_bed_paths <- readLines(list_genes)
        else if (!is.null(genes_file) && file.exists(genes_file)) panel_bed_paths <- readLines(genes_file)
        if (!is.null(panel_files)) {
            if (length(panel_files) == 1 && file.exists(panel_files[1])) {
                panel_bed_paths <- readLines(panel_files[1])
            } else {
                panel_bed_paths <- panel_files
            }
        }

        if (!is.null(panel_bed_paths) && length(panel_bed_paths) > 0) {
            all_panels <- data.table::rbindlist(lapply(panel_bed_paths, function(f) {
                if (!file.exists(f)) {
                    warning("Panel file not found: ", f)
                    return(NULL)
                }
                tmp <- data.table::fread(f, header = FALSE, stringsAsFactors = FALSE)
                if (ncol(tmp) < 4) return(NULL)
                setnames(tmp, paste0("V", 1:4), c("Chr", "Start", "End", "Gene"))
                tmp[, Chr := trimws(as.character(Chr))]
                tmp
            }), fill = TRUE)
            if (nrow(all_panels) == 0) stop("No valid panel BED rows loaded.")
            all_panels <- unique(all_panels)
            all_panels$Chr <- align_to_input_chr(all_panels$Chr)

            input_start_1based <- if (bed_zero_based) input_df$Start + 1 else input_df$Start
            bed_gr <- GenomicRanges::GRanges(
                seqnames = input_df$Chr,
                ranges = IRanges::IRanges(start = input_start_1based, end = input_df$End)
            )
            panel_gr <- GenomicRanges::GRanges(
                seqnames = all_panels$Chr,
                ranges = IRanges::IRanges(start = all_panels$Start + 1, end = all_panels$End),
                Gene = all_panels$Gene
            )
            hits <- suppressWarnings(GenomicRanges::findOverlaps(bed_gr, panel_gr))
            sel <- select_best_hits(bed_gr, panel_gr, hits)
            qh <- sel$q; sh <- sel$s
            output_start <- GenomicRanges::start(bed_gr[qh])
            if (bed_zero_based) output_start <- output_start - 1
            df <- data.frame(
                Chr = as.character(GenomicRanges::seqnames(bed_gr[qh])),
                Start = output_start,
                End = GenomicRanges::end(bed_gr[qh]),
                Gene = panel_gr$Gene[sh],
                Custom.Exon = NA
            )
        } else {
            # No panel restriction: keep original BED but adjust coordinates
            df <- input_df[, 1:4]
            colnames(df) <- c("Chr", "Start", "End", "Gene")
            if (bed_zero_based) {
                # Already 0‑based, keep as is
                df$Custom.Exon <- NA
            } else {
                # Input is 1‑based; we need to output 0‑based (standard BED)
                df$Start <- df$Start - 1
                df$Custom.Exon <- NA
            }
        }
    } else {
        stop("bed_process must be 'STANDARD', 'REGEN', or 'NO'")
    }

    # ----------------------------
    # POST-PROCESSING (common to all modes)
    # ----------------------------
    df$Gene <- as.character(df$Gene)
    df$Gene <- vapply(strsplit(df$Gene, ",", fixed = TRUE), function(x) if (length(x) > 0) x[1] else NA_character_, character(1))

    if (!unknown_gene) {
        df <- df[!is.na(df$Gene) & df$Gene != "" & df$Gene != ".", ]
    } else {
        df$Gene[is.na(df$Gene) | df$Gene == ""] <- "Unknown"
    }
    df <- unique(df)

    # Genomic sort: convert to data.table for efficient ordering
    chrom_base <- c(as.character(1:22), "X", "Y", "M")
    chrom_levels <- if (any(grepl("^chr", df$Chr))) paste0("chr", chrom_base) else chrom_base
    df$Chr <- factor(df$Chr, levels = chrom_levels)
    
    dt <- data.table::as.data.table(df)
    data.table::setorder(dt, Chr, Start, End)
    dt[, Chr := as.character(Chr)]

    # Select output columns
    out_cols <- if ("Custom.Exon" %in% colnames(dt)) {
        c("Chr", "Start", "End", "Gene", "Custom.Exon")
    } else {
        c("Chr", "Start", "End", "Gene")
    }

    # Write output using data.table (correct syntax)
    data.table::fwrite(dt[, ..out_cols], file = output_bed, sep = "\t", col.names = FALSE, quote = FALSE)
    message("[INFO] BED written: ", output_bed)
    invisible(output_bed)
}