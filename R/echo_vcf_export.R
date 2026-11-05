#' Export CNV calls to VCF format
#'
#' Converts ECHO CNV calls to VCF 4.2. When `sample_name` is NULL, a single
#' multi-sample VCF is written. `sample_names` can be supplied to retain
#' analysis samples that have no CNV records in the multi-sample header.
#' Empty calls are supported: a valid header-only VCF is written instead of
#' silently returning without creating a file.
#'
#' @param cnv_calls Data frame of CNV calls. For non-empty input it must contain
#'   `Sample`, `Chromosome`, `Start`, `End`, `Type` and the fields used in INFO.
#' @param output_vcf Character string. Output VCF file path.
#' @param sample_name Optional single sample name. If provided, write a
#'   single-sample VCF, including a header-only VCF when that sample has no calls.
#' @param source Character string. Source description for the VCF header.
#' @param reference Character string. Reference genome assembly (e.g. `hg19`).
#' @param sample_names Optional character vector of analysis sample names. Used
#'   to include samples with no calls in a combined VCF header. Ignored when
#'   `sample_name` is supplied. If omitted, names are inferred from `cnv_calls`.
#' @return Invisibly returns the written VCF path.
#' @export
export_cnvs_to_vcf <- function(cnv_calls, output_vcf, sample_name = NULL,
                               source = "ECHO", reference = "hg19",
                               sample_names = NULL) {
  if (is.null(cnv_calls)) cnv_calls <- data.frame()
  cnv_calls <- as.data.frame(cnv_calls, stringsAsFactors = FALSE,
                             check.names = FALSE)

  required <- c("Sample", "Chromosome", "Start", "End", "Type")
  if (nrow(cnv_calls) > 0L) {
    missing <- setdiff(required, names(cnv_calls))
    if (length(missing)) {
      stop("[ERROR] CNV calls are missing required VCF column(s): ",
           paste(missing, collapse = ", "))
    }
  }

  if (!is.null(sample_name)) {
    if (length(sample_name) != 1L || is.na(sample_name) || !nzchar(sample_name)) {
      stop("[ERROR] sample_name must be one non-empty sample name or NULL.")
    }
    sample_name <- as.character(sample_name)
    if (nrow(cnv_calls) > 0L) {
      cnv_calls <- cnv_calls[as.character(cnv_calls$Sample) == sample_name,
                             , drop = FALSE]
    }
    vcf_samples <- sample_name
  } else if (!is.null(sample_names) && length(sample_names) > 0L) {
    vcf_samples <- unique(as.character(sample_names))
  } else if (nrow(cnv_calls) > 0L) {
    vcf_samples <- sort(unique(as.character(cnv_calls$Sample)))
  } else {
    vcf_samples <- character(0)
  }

  if (anyNA(vcf_samples) || any(!nzchar(vcf_samples)) ||
      any(grepl("[\\t\\r\\n]", vcf_samples))) {
    stop("[ERROR] VCF sample names must be non-empty and cannot contain tabs or newlines.")
  }
  if (anyDuplicated(vcf_samples)) {
    stop("[ERROR] VCF sample names must be unique.")
  }
  if (nrow(cnv_calls) > 0L) {
    calls_samples <- unique(as.character(cnv_calls$Sample))
    absent_samples <- setdiff(calls_samples, vcf_samples)
    if (length(absent_samples)) {
      stop("[ERROR] CNV calls contain sample(s) absent from the requested VCF header: ",
           paste(absent_samples, collapse = ", "))
    }
  }

  safe_char <- function(x, default = ".") {
    if (length(x) == 0L || is.na(x[1]) || as.character(x[1]) == "NA") {
      default
    } else {
      as.character(x[1])
    }
  }

  header <- c(
    "##fileformat=VCFv4.2",
    paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
    paste0("##source=", source),
    paste0("##reference=", reference),
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant (DUP/DEL)">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
    '##INFO=<ID=CONFIDENCE,Number=1,Type=String,Description="Confidence level (HIGH/MEDIUM/LOW)">',
    '##INFO=<ID=REFS,Number=.,Type=String,Description="Reference samples used (comma-separated)">',
    '##INFO=<ID=CORR,Number=1,Type=Float,Description="Correlation with reference">',
    '##INFO=<ID=NCOMP,Number=1,Type=Integer,Description="Number of reference samples">',
    '##INFO=<ID=BF,Number=1,Type=Float,Description="Bayes factor">',
    '##INFO=<ID=READS_RATIO,Number=1,Type=Float,Description="Observed/expected read ratio">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Estimated copy number state (1=DEL,2=normal,3=DUP)">',
    '##FORMAT=<ID=FR,Number=1,Type=Float,Description="Observed/expected fold change">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##ALT=<ID=DEL,Description="Deletion">'
  )

  fixed_columns <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  if (length(vcf_samples) > 0L) {
    header <- c(header, paste(c(fixed_columns, "FORMAT", vcf_samples), collapse = "\t"))
  } else {
    # VCFs without sample columns end at INFO; do not emit a FORMAT-only column.
    header <- c(header, paste(fixed_columns, collapse = "\t"))
  }

  key_df <- if (nrow(cnv_calls) > 0L) {
    unique(data.frame(
      Chromosome = as.character(cnv_calls$Chromosome),
      Start = as.integer(cnv_calls$Start),
      End = as.integer(cnv_calls$End),
      Type = as.character(cnv_calls$Type),
      stringsAsFactors = FALSE
    ))
  } else {
    data.frame(Chromosome = character(0), Start = integer(0), End = integer(0),
               Type = character(0), stringsAsFactors = FALSE)
  }

  chrom_order <- function(chr) {
    stripped <- sub("^chr", "", chr, ignore.case = TRUE)
    num <- suppressWarnings(as.integer(stripped))
    special <- match(toupper(stripped), c("X", "Y", "M", "MT"))
    ifelse(!is.na(num), num,
           ifelse(!is.na(special), 22 + special, Inf))
  }
  if (nrow(key_df) > 0L) {
    key_df <- key_df[order(chrom_order(key_df$Chromosome), key_df$Chromosome,
                           key_df$Start, key_df$End, key_df$Type), , drop = FALSE]
  }

  multi_sample <- length(vcf_samples) > 1L
  make_sample_field <- function(row_df, sv_type) {
    if (nrow(row_df) == 0L) {
      # Multi-sample VCF: absence of a call at this locus is homozygous reference.
      # Single-sample VCF: only called loci are emitted, so this is missing.
      return(if (multi_sample) "0/0:2:1.00" else "./.:.:.")
    }
    ratio <- suppressWarnings(as.numeric(row_df$Reads.ratio[1]))
    if (!length(ratio) || !is.finite(ratio)) ratio <- 1
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
    key_row <- key_df[i, , drop = FALSE]
    chrom <- as.character(key_row$Chromosome)
    start <- as.integer(key_row$Start)
    end <- as.integer(key_row$End)
    cnv_type <- tolower(as.character(key_row$Type))
    is_dup <- grepl("dup", cnv_type)
    alt <- if (is_dup) "<DUP>" else "<DEL>"
    sv_type <- if (is_dup) "DUP" else "DEL"
    sv_len <- if (is_dup) end - start + 1L else -(end - start + 1L)

    row_hits <- cnv_calls[
      as.character(cnv_calls$Chromosome) == chrom &
        as.integer(cnv_calls$Start) == start &
        as.integer(cnv_calls$End) == end &
        tolower(as.character(cnv_calls$Type)) == cnv_type,
      , drop = FALSE
    ]
    representative <- row_hits[1, , drop = FALSE]

    id <- if ("CNV.ID" %in% names(representative)) {
      paste0("CNV_", representative$`CNV.ID`)
    } else {
      paste0(gsub("^chr", "", chrom, ignore.case = TRUE), "_", start, "_", end, "_", sv_type)
    }
    info_items <- c(paste0("SVTYPE=", sv_type), paste0("SVLEN=", sv_len),
                    paste0("END=", end))

    conf <- safe_char(representative$Confidence)
    if (conf != ".") info_items <- c(info_items, paste0("CONFIDENCE=", conf))
    if ("Comparator.name" %in% names(representative)) {
      refs_val <- safe_char(representative$Comparator.name, "")
      if (nzchar(refs_val)) {
        refs_val <- gsub(",\\s*", ",", refs_val)
        info_items <- c(info_items, paste0("REFS=", refs_val))
      }
    }
    if ("Correlation" %in% names(representative)) {
      corr <- representative$Correlation[1]
      if (length(corr) && !is.na(corr) && corr != "NA")
        info_items <- c(info_items, paste0("CORR=", corr))
    }
    if ("N.comp" %in% names(representative)) {
      ncomp <- representative$`N.comp`[1]
      if (length(ncomp) && !is.na(ncomp) && ncomp != "NA")
        info_items <- c(info_items, paste0("NCOMP=", ncomp))
    }
    if ("BF" %in% names(representative)) {
      bf <- representative$BF[1]
      if (length(bf) && !is.na(bf) && bf != "NA")
        info_items <- c(info_items, paste0("BF=", bf))
    }
    if ("Reads.ratio" %in% names(representative)) {
      ratio <- representative$`Reads.ratio`[1]
      if (length(ratio) && !is.na(ratio) && ratio != "NA")
        info_items <- c(info_items, paste0("READS_RATIO=", round(as.numeric(ratio), 2)))
    }

    sample_fields <- vapply(vcf_samples, function(s) {
      sample_row <- row_hits[as.character(row_hits$Sample) == s, , drop = FALSE]
      make_sample_field(sample_row, sv_type)
    }, character(1))

    record_fields <- c(chrom, start, id, "N", alt, ".", "PASS",
                       paste(info_items, collapse = ";"))
    if (length(vcf_samples) > 0L) {
      record_fields <- c(record_fields, "GT:CN:FR", sample_fields)
    }
    records[i] <- paste(record_fields, collapse = "\t")
  }

  dir.create(dirname(output_vcf), showWarnings = FALSE, recursive = TRUE)
  con <- file(output_vcf, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(header, con, useBytes = TRUE)
  if (length(records)) writeLines(records, con, useBytes = TRUE)
  message("[INFO] Wrote ", length(records), " CNV record(s) to ", output_vcf)
  invisible(output_vcf)
}
