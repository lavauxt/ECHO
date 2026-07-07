#' Export CNV calls to VCF format
#'
#' Converts a data frame of CNV calls (from \code{run_cnv_calling}) into VCF.
#' When \code{sample_name} is NULL, a single multi-sample VCF is written.
#'
#' @param cnv_calls Data frame of CNV calls (must contain columns: Sample, Chromosome,
#'   Start, End, Type, Confidence, etc.)
#' @param output_vcf Character string. Output VCF file path.
#' @param sample_name Character string. If provided, write a single-sample VCF.
#'   If NULL, write a single multi-sample VCF.
#' @param source Character string. Source description for VCF header.
#' @param reference Character string. Reference genome assembly (e.g., "hg19").
#' @return Invisibly returns the written VCF path.
#' @export
export_cnvs_to_vcf <- function(cnv_calls, output_vcf, sample_name = NULL,
                               source = "ECHO", reference = "hg19") {
  if (is.null(cnv_calls) || nrow(cnv_calls) == 0) {
    message("[INFO] No CNVs to export.")
    return(invisible(NULL))
  }
  if (!is.null(sample_name)) {
    cnv_calls <- cnv_calls[cnv_calls$Sample == sample_name, , drop = FALSE]
    if (nrow(cnv_calls) == 0) {
      warning("No CNVs found for sample ", sample_name)
      return(invisible(NULL))
    }
  }

  safe_char <- function(x, default = ".") if (length(x) == 0 || is.na(x) || x == "NA") default else as.character(x)
  samples   <- sort(unique(as.character(cnv_calls$Sample)))

  header <- c(
    "##fileformat=VCFv4.2",
    paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
    paste0("##source=", source),
    paste0("##reference=", reference),
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant (DUP/DEL)">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
    '##INFO=<ID=CONFIDENCE,Number=1,Type=String,Description="Confidence level (HIGH/MEDIUM/LOW)">',
    '##INFO=<ID=REFS,Number=.,Type=String,Description="Reference samples used (space-separated)">',
    '##INFO=<ID=CORR,Number=1,Type=Float,Description="Correlation with reference">',
    '##INFO=<ID=NCOMP,Number=1,Type=Integer,Description="Number of reference samples">',
    '##INFO=<ID=BF,Number=1,Type=Float,Description="Bayes factor">',
    '##INFO=<ID=READS_RATIO,Number=1,Type=Float,Description="Observed/expected read ratio">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Estimated copy number state (1=DEL,2=normal,3=DUP)">',
    '##FORMAT=<ID=FR,Number=1,Type=Float,Description="Observed/expected fold change">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=DEL,Description="Deletion">',
    paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t",
           paste(samples, collapse = "\t"))
  )

  key_df <- unique(data.frame(
    Chromosome = as.character(cnv_calls$Chromosome),
    Start      = as.integer(cnv_calls$Start),
    End        = as.integer(cnv_calls$End),
    Type       = as.character(cnv_calls$Type),
    stringsAsFactors = FALSE
  ))

  chrom_order <- function(chr) {
    stripped <- sub("^chr", "", chr)
    num      <- suppressWarnings(as.integer(stripped))
    ifelse(is.na(num), Inf, num)   # X/Y/M sort after numbered chromosomes
  }
  key_df <- key_df[order(chrom_order(key_df$Chromosome),
                         key_df$Chromosome,   # tie-break for X, Y, M
                         key_df$Start,
                         key_df$End,
                         key_df$Type), , drop = FALSE]

  make_sample_field <- function(row_df, sv_type) {
    if (nrow(row_df) == 0) return("0/0:2:1.00")
    ratio <- suppressWarnings(as.numeric(row_df$Reads.ratio[1]))
    if (is.na(ratio)) ratio <- 1

    if (sv_type == "DEL" && ratio < 0.3) {
      cn <- 0L; gt <- "1/1"
    } else if (sv_type == "DEL") {
      cn <- 1L; gt <- "0/1"
    } else {
      cn <- 3L; gt <- "0/1"
    }
    paste0(gt, ":", cn, ":", format(round(ratio, 2), nsmall = 2, trim = TRUE))
  }

  records <- character(nrow(key_df))
  for (i in seq_len(nrow(key_df))) {
    key_row  <- key_df[i, ]
    chrom    <- as.character(key_row[["Chromosome"]])
    start    <- as.integer(key_row[["Start"]])
    end      <- as.integer(key_row[["End"]])
    cnv_type <- tolower(as.character(key_row[["Type"]]))
    is_dup   <- grepl("dup", cnv_type)
    alt      <- if (is_dup) "<DUP>" else "<DEL>"
    sv_type  <- if (is_dup) "DUP" else "DEL"
    sv_len   <- if (is_dup) end - start + 1L else -(end - start + 1L)

    row_hits    <- cnv_calls[
      as.character(cnv_calls$Chromosome) == chrom &
      as.integer(cnv_calls$Start)        == start &
      as.integer(cnv_calls$End)          == end   &
      tolower(as.character(cnv_calls$Type)) == cnv_type, , drop = FALSE]
    representative <- row_hits[1, , drop = FALSE]

    id <- if ("CNV.ID" %in% colnames(representative))
      paste0("CNV_", representative$`CNV.ID`)
    else
      paste0(gsub("^chr", "", chrom), "_", start, "_", end, "_", sv_type)

    info_items <- c(
      paste0("SVTYPE=", sv_type),
      paste0("SVLEN=", sv_len),
      paste0("END=", end)
    )
    conf <- safe_char(representative$Confidence)
    if (conf != ".") info_items <- c(info_items, paste0("CONFIDENCE=", conf))

    if ("Comparator.name" %in% colnames(representative)) {
      refs_val <- safe_char(representative$Comparator.name, "")
      if (nzchar(refs_val)) {

        refs_val <- gsub(",\\s*", ",", refs_val)
        info_items <- c(info_items, paste0("REFS=", refs_val))
      }
    }
    if ("Correlation" %in% colnames(representative)) {
      corr <- representative$Correlation
      if (!is.na(corr) && corr != "NA") info_items <- c(info_items, paste0("CORR=", corr))
    }
    if ("N.comp" %in% colnames(representative)) {
      ncomp <- representative$`N.comp`
      if (!is.na(ncomp) && ncomp != "NA") info_items <- c(info_items, paste0("NCOMP=", ncomp))
    }
    if ("BF" %in% colnames(representative)) {
      bf <- representative$BF
      if (!is.na(bf) && bf != "NA") info_items <- c(info_items, paste0("BF=", bf))
    }
    if ("Reads.ratio" %in% colnames(representative)) {
      ratio <- representative$`Reads.ratio`
      if (!is.na(ratio) && ratio != "NA")
        info_items <- c(info_items, paste0("READS_RATIO=", round(as.numeric(ratio), 2)))
    }

    format_field  <- "GT:CN:FR"
    sample_fields <- vapply(samples, function(s) {
      sample_row <- row_hits[row_hits$Sample == s, , drop = FALSE]
      make_sample_field(sample_row, sv_type)
    }, character(1))

    records[i] <- paste(chrom, start, id, "N", alt, ".", "PASS",
                        paste(info_items, collapse = ";"),
                        format_field,
                        paste(sample_fields, collapse = "\t"),
                        sep = "\t")
  }

  con <- file(output_vcf, "w")
  on.exit(close(con), add = TRUE)
  writeLines(header,  con)
  writeLines(records, con)
  message("[INFO] Wrote ", length(records), " CNV records to ", output_vcf)
  invisible(output_vcf)
}
