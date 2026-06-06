#' Export CNV calls to VCF format
#'
#' Converts a data frame of CNV calls (from \code{run_cnv_calling}) into a VCF file
#' with conformant structural variant records.
#'
#' @param cnv_calls Data frame of CNV calls (must contain columns: Sample, Chromosome, Start, End, Type, Confidence, etc.)
#' @param output_vcf Character string. Output VCF file path.
#' @param sample_name Character string. Single sample name to export. If NULL, all samples are written to separate VCFs (appending sample name).
#' @param source Character string. Source description for VCF header.
#' @param reference Character string. Reference genome assembly (e.g., "hg19").
#'
#' @return Invisibly returns a vector of the written VCF file path(s).
#' @export
export_cnvs_to_vcf <- function(cnv_calls, output_vcf, sample_name = NULL,
                               source = "ECHO", reference = "hg19") {
    if (nrow(cnv_calls) == 0) {
        message("[INFO] No CNVs to export.")
        return(invisible(NULL))
    }

    # Filter by sample if needed
    if (!is.null(sample_name)) {
        cnv_calls <- cnv_calls[cnv_calls$Sample == sample_name, ]
        if (nrow(cnv_calls) == 0) {
            warning("No CNVs found for sample ", sample_name)
            return(invisible(NULL))
        }
        output_file <- output_vcf
    } else {
        # Write each sample to a separate file
        samples <- unique(cnv_calls$Sample)
        out_files <- c()
        for (s in samples) {
            # Safely insert sample name before extension
            base <- tools::file_path_sans_ext(output_vcf)
            ext <- tools::file_ext(output_vcf)
            out_file <- paste0(base, "_", s, ".", ext)
            export_cnvs_to_vcf(cnv_calls[cnv_calls$Sample == s, ], out_file, sample_name = s, source = source, reference = reference)
            out_files <- c(out_files, out_file)
        }
        return(invisible(out_files))
    }

    # Prepare fully conformant standard VCF v4.2 metadata headers
    header <- c(
        "##fileformat=VCFv4.2",
        paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
        paste0("##source=", source),
        paste0("##reference=", reference),
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant (DUP/DEL)\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
        "##INFO=<ID=CONFIDENCE,Number=1,Type=String,Description=\"Confidence level (HIGH/MEDIUM/LOW)\">",
        "##INFO=<ID=REFS,Number=.,Type=String,Description=\"Reference samples used\">",
        "##INFO=<ID=CORR,Number=1,Type=Float,Description=\"Correlation with reference\">",
        "##INFO=<ID=NCOMP,Number=1,Type=Integer,Description=\"Number of reference samples\">",
        "##INFO=<ID=BF,Number=1,Type=Float,Description=\"Bayes factor\">",
        "##INFO=<ID=READS_RATIO,Number=1,Type=Float,Description=\"Observed/expected read ratio\">",
        "##ALT=<ID=DUP,Description=\"Duplication\">",
        "##ALT=<ID=DEL,Description=\"Deletion\">",
        paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    )

    # Helper to safely get non-NA character
    safe_char <- function(x, default = ".") {
        if (is.na(x) || x == "NA") default else as.character(x)
    }

    # Build VCF records
    records <- apply(cnv_calls, 1, function(row) {
        chrom <- as.character(row["Chromosome"])
        start <- as.integer(row["Start"])
        end   <- as.integer(row["End"])
        cnv_type <- tolower(as.character(row["Type"]))
        
        is_dup <- grepl("dup", cnv_type)
        alt <- ifelse(is_dup, "<DUP>", "<DEL>")
        sv_type <- ifelse(is_dup, "DUP", "DEL")
        
        sv_len <- end - start + 1
        if (!is_dup) sv_len <- -sv_len
        
        id <- paste0(row["Sample"], "_", row["CNV.ID"])
        qual <- "."
        filter <- "PASS"
        
        info_items <- c(
            paste0("SVTYPE=", sv_type),
            paste0("SVLEN=", sv_len),
            paste0("END=", end)
        )
        
        # Add optional fields only if they exist and are not NA/empty
        conf <- safe_char(row["Confidence"])
        if (conf != ".") info_items <- c(info_items, paste0("CONFIDENCE=", conf))
        
        # --- FIX: Keep commas as separators (VCF spec) ---
        refs_val <- safe_char(row["Comparator.name"], "")
        if (refs_val != "") {
            # No transformation needed – commas are valid in VCF INFO lists
            info_items <- c(info_items, paste0("REFS=", refs_val))
        }
        
        corr <- row["Correlation"]
        if (!is.na(corr) && corr != "NA") info_items <- c(info_items, paste0("CORR=", corr))
        
        ncomp <- row["N.comp"]
        if (!is.na(ncomp) && ncomp != "NA") info_items <- c(info_items, paste0("NCOMP=", ncomp))
        
        bf <- row["BF"]
        if (!is.na(bf) && bf != "NA") info_items <- c(info_items, paste0("BF=", bf))
        
        ratio <- row["Reads.ratio"]
        if (!is.na(ratio) && ratio != "NA") info_items <- c(info_items, paste0("READS_RATIO=", ratio))
        
        info_str <- paste(info_items, collapse = ";")
        
        paste(chrom, start, id, "N", alt, qual, filter, info_str, sep = "\t")
    })

    con <- file(output_file, "w")
    writeLines(header, con)
    writeLines(records, con)
    close(con)
    message("[INFO] Wrote ", nrow(cnv_calls), " CNV records to ", output_file)
    invisible(output_file)
}