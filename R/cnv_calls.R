#' Score CNV Confidence (HIGH/MEDIUM/LOW)
#'
#' @param cnv_calls Data frame of ExomeDepth CNV calls.
#' @param bed_file BED annotation data frame aligned with the coverage matrix rows.
#' @param current_cor Numeric correlation between test and reference coverage.
#' @param num_refs Integer number of reference samples used for the call.
#' @param score_high_corr Minimum correlation for HIGH confidence.
#' @param score_med_corr Minimum correlation for MEDIUM confidence.
#' @param score_high_refs Minimum number of reference samples for HIGH.
#' @param score_med_refs Minimum number of reference samples for MEDIUM.
#' @param score_low_ratio_low Lower bound of the borderline read-ratio interval (LOW).
#' @param score_low_ratio_high Upper bound of the borderline read-ratio interval (LOW).
#' @param score_high_ratio_low Ratio below this (with good metrics) gives HIGH.
#' @param score_high_ratio_high Ratio above this (with good metrics) gives HIGH.
#' @param score_med_ratio_low Ratio between this and low_ratio_low -> MEDIUM (if criteria met).
#' @param score_med_ratio_high Ratio between low_ratio_high and this -> MEDIUM.
#' @param score_low_confidence_genes Character vector of gene symbols always flagged LOW.
#'
#' @return The input \code{cnv_calls} data frame with confidence column added.
#' @export
score_cnv_confidence <- function(cnv_calls, bed_file, current_cor, num_refs,
                                 score_high_corr = 0.985,
                                 score_med_corr  = 0.95,
                                 score_high_refs = 3,
                                 score_med_refs  = 2,
                                 score_low_ratio_low   = 0.75,
                                 score_low_ratio_high  = 1.25,
                                 score_high_ratio_low  = 0.70,
                                 score_high_ratio_high = 1.30,
                                 score_med_ratio_low   = 0.60,
                                 score_med_ratio_high  = 1.40,
                                 score_low_confidence_genes = c("PMS2")) {
  if (nrow(cnv_calls) == 0) return(cnv_calls)

  cnv_calls$Gene <- vapply(seq_len(nrow(cnv_calls)), function(i) {
    start_idx <- max(1L, cnv_calls$start.p[i])
    end_idx   <- min(nrow(bed_file), cnv_calls$end.p[i])
    paste(unique(bed_file$gene[start_idx:end_idx]), collapse = " ")
  }, character(1))

  cnv_calls$Start.b <- NULL
  cnv_calls$End.b   <- NULL
  cnv_calls$Confidence <- "LOW"

  high_corr_ok  <- !is.na(current_cor) && current_cor >= score_high_corr
  high_refs_ok  <- num_refs >= score_high_refs
  ratio_high_ok <- cnv_calls$reads.ratio <= score_high_ratio_low |
                   cnv_calls$reads.ratio >= score_high_ratio_high

  med_corr_ok    <- !is.na(current_cor) && current_cor >= score_med_corr
  med_refs_ok    <- num_refs >= score_med_refs
  ratio_not_low  <- cnv_calls$reads.ratio <= score_low_ratio_low |
                    cnv_calls$reads.ratio >= score_low_ratio_high

  # FIX BUG-7: sapply() can silently return a list when the function produces
  # irregular output; vapply() enforces a logical(1) return type.
  gene_flagged <- vapply(cnv_calls$Gene, function(g_str) {
    any(trimws(unlist(strsplit(g_str, "[ ,]"))) %in% score_low_confidence_genes)
  }, logical(1))

  high_cond <- !gene_flagged & high_corr_ok & high_refs_ok & ratio_high_ok
  med_cond  <- !gene_flagged & !high_cond & med_corr_ok & med_refs_ok & ratio_not_low

  cnv_calls$Confidence[high_cond] <- "HIGH"
  cnv_calls$Confidence[med_cond]  <- "MEDIUM"

  return(cnv_calls)
}

#' Run CNV Calling using ExomeDepth
#'
#' @param rdata_file Character string. Output from \code{\link{run_bam_coverage}}.
#' @param output_file Character string. TSV output for CNV calls.
#' @param out_rdata Character string. RData output for plotting.
#' @param transition.probability Numeric. HMM transition probability.
#' @param expected.CNV.length Numeric. Expected CNV length in bp.
#' @param n.bins.reduced Integer. Number of bins to subsample for reference selection.
#' @param phi.bins Integer. Number of bins for over-dispersion parameter.
#' @param formula Character. Formula for the model.
#' @param data Data frame. Optional covariate data for the formula.
#' @param save_ed_objects Logical. If TRUE, save full ExomeDepth objects per sample.
#' @param modechrom Character. Mode for sex chromosome handling.
#' @param sample_table Optional character string. Path to sample table.
#' @param ... Additional parameters passed to \code{score_cnv_confidence}.
#'
#' @return Invisibly returns \code{NULL}. Writes TSV and RData outputs to disk.
#' @export
run_cnv_calling <- function(rdata_file,
                            output_file = "./CNV_calls.tsv",
                            out_rdata   = "./ECHO_summary.RData",
                            transition.probability = 1e-4,
                            expected.CNV.length    = 50000,
                            n.bins.reduced         = 10000,
                            phi.bins               = 1,
                            formula                = "cbind(test, reference) ~ 1",
                            data                   = NULL,
                            save_ed_objects        = FALSE,
                            modechrom              = "A",
                            sample_table           = NULL,
                            ...) {
  message("[INFO] ", Sys.time(), " BEGIN CNV calls")
  objs         <- load_rdata(rdata_file, required = c("counts", "sample_names", "bed_file"))
  counts       <- objs$counts
  sample_names <- objs$sample_names
  bed_file     <- objs$bed_file

  if (modechrom == "X") {
    message("[INFO] Mode X: restricting to chromosome X")
    chr_x    <- if (any(grepl("^chrX", bed_file$chromosome))) "chrX" else "X"
    keep_idx <- which(bed_file$chromosome == chr_x)
    if (length(keep_idx) == 0) stop("[ERROR] No chrX intervals found")
    bed_file <- bed_file[keep_idx, ]
    counts   <- counts[keep_idx, ]
  } else if (modechrom == "Y") {
    message("[INFO] Mode Y: restricting to chromosome Y")
    chr_y    <- if (any(grepl("^chrY", bed_file$chromosome))) "chrY" else "Y"
    keep_idx <- which(bed_file$chromosome == chr_y)
    if (length(keep_idx) == 0) stop("[ERROR] No chrY intervals found")
    bed_file <- bed_file[keep_idx, ]
    counts   <- counts[keep_idx, ]
  } else {
    message("[INFO] Mode A: autosomes only (excluding X and Y)")
    autosomes <- if (any(grepl("^chr", bed_file$chromosome))) paste0("chr", 1:22) else as.character(1:22)
    keep_idx  <- which(bed_file$chromosome %in% autosomes)
    if (length(keep_idx) == 0) stop("[ERROR] No autosomal intervals found")
    bed_file <- bed_file[keep_idx, ]
    counts   <- counts[keep_idx, ]
  }

  gender_map <- NULL
  if (modechrom %in% c("X", "Y")) {
    if (is.null(sample_table) || !file.exists(sample_table))
      stop("[ERROR] sample_table required for mode ", modechrom)
    gender_df <- utils::read.table(sample_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (!all(c("sample_name", "gender") %in% colnames(gender_df)))
      stop("[ERROR] sample_table must have 'sample_name' and 'gender'")
    gender_df$gender <- toupper(substr(gender_df$gender, 1, 1))
    gender_map <- setNames(gender_df$gender, gender_df$sample_name)
    missing    <- setdiff(sample_names, names(gender_map))
    if (length(missing) > 0) stop("[ERROR] Missing gender for samples: ", paste(missing, collapse = ", "))
  }

  if (modechrom == "Y") {
    male_samples <- names(gender_map[gender_map == "M"])
    sample_names <- intersect(sample_names, male_samples)
    if (length(sample_names) == 0) stop("[ERROR] No male samples for mode Y")
    message("[INFO] Mode Y: keeping only male samples (", length(sample_names), " samples)")
  }

  clean_chroms <- sub("^chr", "", counts$chromosome)
  results      <- list()
  ed_objects   <- list()

  for (sample in sample_names) {
    message("[INFO] Processing sample: ", sample)
    if (modechrom == "A") {
      ref_samples <- sample_names[sample_names != sample]
    } else {
      this_gender <- gender_map[[sample]]
      if (modechrom == "X") {
        ref_samples <- names(gender_map[gender_map == this_gender & names(gender_map) != sample])
      } else {
        ref_samples <- sample_names[sample_names != sample]
      }
    }

    if (length(ref_samples) == 0) {
      message("[WARNING] No reference samples available for ", sample, ". Skipping.")
      results[[sample]] <- list(calls = data.frame(), model = NULL, refs = character(0))
      next
    }

    test_counts <- counts[, sample, drop = TRUE]
    ref_matrix  <- as.matrix(counts[, ref_samples, drop = FALSE])

    best_refs <- suppressWarnings(tryCatch({
      ExomeDepth::select.reference.set(
        test.counts      = test_counts,
        reference.counts = ref_matrix,
        bin.length       = (counts$end - counts$start),
        n.bins.reduced   = n.bins.reduced,
        phi.bins         = phi.bins,
        formula          = formula,
        data             = data
      )$reference.choice
    }, error = function(e) {
      message("[WARNING] Reference selection failed for ", sample, ": ", e$message)
      character(0)
    }))

    if (length(best_refs) == 0) {
      message("[WARNING] Model selected 0 references for ", sample)
      results[[sample]] <- list(calls = data.frame(), model = NULL, refs = character(0))
      next
    }

    ref_counts <- rowSums(counts[, best_refs, drop = FALSE])
    suppressWarnings({
      ed <- methods::new("ExomeDepth", test = test_counts, reference = ref_counts, formula = formula)
      ed <- ExomeDepth::CallCNVs(
        x                      = ed,
        transition.probability = transition.probability,
        expected.CNV.length    = expected.CNV.length,
        chromosome             = clean_chroms,
        start                  = counts$start,
        end                    = counts$end,
        name                   = counts$exon
      )
    })

    raw_calls       <- ed@CNV.calls
    processed_calls <- data.frame()
    if (!is.null(raw_calls) && nrow(raw_calls) > 0) {
      current_cor     <- as.numeric(stats::cor(counts[, sample], ref_counts))
      num_refs        <- length(best_refs)
      processed_calls <- score_cnv_confidence(raw_calls, bed_file, current_cor, num_refs, ...)
      processed_calls <- add_within_gene_indices(processed_calls, bed_file)
      processed_calls <- cbind(Sample = sample, processed_calls)
      processed_calls$Correlation      <- current_cor
      processed_calls$N.comp           <- num_refs
      processed_calls$Comparator.name  <- paste(best_refs, collapse = ", ")

      names(processed_calls)[names(processed_calls) == "start.p"]       <- "Start.p"
      names(processed_calls)[names(processed_calls) == "end.p"]         <- "End.p"
      names(processed_calls)[names(processed_calls) == "type"]          <- "Type"
      names(processed_calls)[names(processed_calls) == "nexons"]        <- "Nexons"
      names(processed_calls)[names(processed_calls) == "start"]         <- "Start"
      names(processed_calls)[names(processed_calls) == "end"]           <- "End"
      names(processed_calls)[names(processed_calls) == "chromosome"]    <- "Chromosome"
      names(processed_calls)[names(processed_calls) == "bayesfactor"]   <- "BF"
      names(processed_calls)[names(processed_calls) == "reads.expected"] <- "Reads.expected"
      names(processed_calls)[names(processed_calls) == "reads.observed"] <- "Reads.observed"
      names(processed_calls)[names(processed_calls) == "reads.ratio"]   <- "Reads.ratio"
      processed_calls$Genomic.ID <- paste0(processed_calls$Chromosome, ":",
                                            processed_calls$Start, "-", processed_calls$End)
    }

    results[[sample]] <- list(
      calls = processed_calls,
      model = c(ed@phi,
                sum(counts[, sample]) / (sum(counts[, sample]) + sum(ref_counts))),
      refs  = best_refs
    )
    if (save_ed_objects) ed_objects[[sample]] <- ed
  }

  all_calls_list <- lapply(results, function(x) x$calls)
  all_calls_list <- Filter(Negate(is.null), all_calls_list)
  models         <- lapply(results, function(x) x$model)
  names(models)  <- sample_names
  refs           <- lapply(results, function(x) x$refs)
  names(refs)    <- sample_names

  cnv_calls <- if (length(all_calls_list) > 0) do.call(rbind, all_calls_list) else data.frame()
  if (nrow(cnv_calls) > 0) {
    cnv_calls  <- cbind(CNV.ID = seq_len(nrow(cnv_calls)), cnv_calls)
    exclude_cols <- c("global_start", "global_end")
    tsv_cols     <- setdiff(colnames(cnv_calls), exclude_cols)
    utils::write.table(cnv_calls[, tsv_cols], file = output_file,
                       row.names = FALSE, sep = "\t", quote = FALSE)
    message("[INFO] Saved CNV calls to: ", output_file)
  } else {
    message("[INFO] No CNVs detected across any samples.")
  }

  if (save_ed_objects) {
    save(counts, cnv_calls, refs, models, bed_file, ed_objects, gender_map, file = out_rdata)
  } else {
    save(counts, cnv_calls, refs, models, bed_file, gender_map, file = out_rdata)
  }
  message("[INFO] Saved summary RData to: ", out_rdata)
  message("[INFO] ", Sys.time(), " END CNV calls")
  invisible(NULL)
}
