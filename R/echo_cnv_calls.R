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
#' @param penalize_terminal_only Logical. If TRUE, any call flagged as `terminal_only`
#'   is forced to "LOW" confidence.
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
                                 score_low_confidence_genes = c("PMS2", "SMN1", "CYP2D6", "HBA1", "HBA2",
                                                                 "STRC", "CYP21A2", "GBA1", "CFTR"),
                                 penalize_terminal_only = FALSE) {
  if (nrow(cnv_calls) == 0) return(cnv_calls)
  
  if (!"terminal_only" %in% colnames(cnv_calls)) {
    cnv_calls <- flag_terminal_only_calls(cnv_calls, bed_file)
  }
  
  cnv_calls$Gene <- genes_in_targets(bed_file$gene, cnv_calls$start.p, cnv_calls$end.p)
  
  cnv_calls <- cnv_calls[, !tolower(names(cnv_calls)) %in% c("start.b", "end.b"), drop = FALSE]
  cnv_calls$Confidence <- "LOW"
  
  high_corr_ok  <- !is.na(current_cor) && current_cor >= score_high_corr
  high_refs_ok  <- num_refs >= score_high_refs
  ratio <- as.numeric(cnv_calls$reads.ratio)
  ratio_high_ok <- !is.na(ratio) &
    (ratio <= score_high_ratio_low | ratio >= score_high_ratio_high)
  
  med_corr_ok <- !is.na(current_cor) && current_cor >= score_med_corr
  med_refs_ok <- num_refs >= score_med_refs
  ratio_medium_ok <- !is.na(ratio) &
    ((ratio >= score_med_ratio_low & ratio <= score_low_ratio_low) |
     (ratio >= score_low_ratio_high & ratio <= score_med_ratio_high))
  
  gene_flagged <- vapply(cnv_calls$Gene, function(g_str) {
    any(trimws(unlist(strsplit(g_str, "[ ,;]"))) %in% score_low_confidence_genes)
  }, logical(1))
  
  high_cond <- !gene_flagged & high_corr_ok & high_refs_ok & ratio_high_ok
  med_cond  <- !gene_flagged & !high_cond & med_corr_ok & med_refs_ok & ratio_medium_ok
  
  cnv_calls$Confidence[high_cond] <- "HIGH"
  cnv_calls$Confidence[med_cond]  <- "MEDIUM"
  
  terminal_idx <- which(!is.na(cnv_calls$terminal_only) & cnv_calls$terminal_only)
  if (penalize_terminal_only && length(terminal_idx)) {
    cnv_calls$Confidence[terminal_idx] <- "LOW"
  }
  
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
#' @param sample_qc Logical. Exclude outlier samples (robust z-score on
#'   cross-target noise) before calling.
#' @param exon_qc Logical. Exclude problematic exons before calling.
#' @param qc_zscore Numeric. Robust z-score threshold for \code{sample_qc}.
#' @param outlier_sides \code{"upper"} (default; only unusually noisy samples are
#'   outliers, like CANOPE) or \code{"both"} (also unusually clean ones; the
#'   pre-1.x ECHO behaviour). See \code{\link{detect_outlier_samples}}.
#' @param exon_mad_quantile Numeric. MAD quantile threshold for \code{exon_qc}.
#' @param gc_extreme_filter Numeric[2] or \code{NULL}. Min/max GC fraction to retain.
#' @param min_exon_mean Numeric. Minimum mean count to retain a target.
#' @param remove_terminal_only Logical. If TRUE, remove calls that affect only
#'   first/last exons of a gene (after scoring). Default FALSE.
#' @param penalize_terminal_only Logical. If TRUE, force such calls to LOW
#'   confidence (instead of removing them). Default FALSE.
#' @param recurrence_action What to do with calls that recur in many samples of
#'   the cohort (a hallmark of a noisy target rather than a germline CNV):
#'   \code{"flag"} (default; adds \code{Carriers}/\code{Recurrent}),
#'   \code{"downgrade"} (also force LOW confidence), \code{"remove"} or
#'   \code{"off"}. See \code{\link{cnv_apply_recurrence}}.
#' @param recurrence_min_samples,recurrence_min_fraction,recurrence_min_overlap
#'   A call is recurrent when >= \code{min_samples} samples (and >=
#'   \code{min_fraction} of the analysed cohort) carry a same-type call with
#'   >= \code{min_overlap} reciprocal overlap.
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
                            sample_qc              = TRUE,
                            exon_qc                = TRUE,
                            qc_zscore              = 3,
                            outlier_sides          = "upper",
                            exon_mad_quantile      = 0.90,
                            gc_extreme_filter      = c(0.15, 0.85),
                            min_exon_mean          = 20,
                            remove_terminal_only   = FALSE,
                            penalize_terminal_only = FALSE,
                            recurrence_action      = "flag",
                            recurrence_min_samples = 3L,
                            recurrence_min_fraction = 0.25,
                            recurrence_min_overlap = 0.5,
                            ...) {
  message("[INFO] BEGIN CNV calls")
  md <- echo_mode_data(objs = load_rdata(rdata_file, required = c("counts", "sample_names", "bed_file")),
                       modechrom = modechrom, sample_table = sample_table)
  modechrom              <- md$modechrom
  counts                 <- md$counts
  bed_file               <- md$bed_file
  sample_names           <- md$sample_names
  reference_sample_names <- md$reference_sample_names
  candidate_sample_names <- md$candidate
  gender_map             <- md$gender_map
  n_targets_mode         <- nrow(counts)

  # ONE QC decision for both pipelines and for the QC-metrics step (see
  # cnv_apply_qc()): outlier samples among all candidate samples, then the
  # problematic targets among the samples that remain.
  gc_for_qc <- if ("GC" %in% colnames(counts)) as.numeric(counts$GC) else NULL
  if (isTRUE(exon_qc) && is.null(gc_for_qc) && !is.null(gc_extreme_filter))
    message("[INFO] exon_qc: no 'GC' column in the coverage RData (re-run ",
            "run_bam_coverage() to pick it up); skipping the GC-extreme check.")
  qc <- cnv_apply_qc(
    counts[, candidate_sample_names, drop = FALSE], chromosomes = counts$chromosome, gc = gc_for_qc,
    sample_qc = sample_qc, exon_qc = exon_qc, qc_zscore = qc_zscore, outlier_sides = outlier_sides,
    exon_mad_quantile = exon_mad_quantile, min_exon_mean = min_exon_mean, gc_range = gc_extreme_filter)

  excl <- qc$excluded_samples$sample
  if (length(excl) > 0) {
    message(sprintf("[WARNING] sample_qc: excluding %d outlier sample(s) from calling (used neither as test nor reference): %s",
                    length(excl), paste(sprintf("%s (z = %.2f)", excl,
                      qc$sample_qc$robust_z[match(excl, qc$sample_qc$sample)]), collapse = ", ")))
    sample_names           <- setdiff(sample_names, excl)
    reference_sample_names <- setdiff(reference_sample_names, excl)
    candidate_sample_names <- setdiff(candidate_sample_names, excl)
    if (length(sample_names) == 0)
      stop("[ERROR] sample_qc excluded every analysis sample; relax qc_zscore or set sample_qc = FALSE.")
  }

  removed_targets <- NULL
  if (isTRUE(exon_qc) && length(qc$keep_targets) < nrow(counts)) {
    drop_idx <- setdiff(seq_len(nrow(counts)), qc$keep_targets)
    removed_targets <- data.frame(
      chromosome = as.character(bed_file$chromosome[drop_idx]), start = bed_file$start[drop_idx],
      end = bed_file$end[drop_idx], gene = as.character(bed_file$gene[drop_idx]),
      label = target_labels(bed_file)[drop_idx], reason = qc$target_qc$reason[drop_idx],
      stringsAsFactors = FALSE)
    message(sprintf("[INFO] exon_qc: excluding %d/%d problematic target(s) before calling",
                    length(drop_idx), nrow(counts)))
    counts   <- counts[qc$keep_targets, , drop = FALSE]
    bed_file <- bed_file[qc$keep_targets, , drop = FALSE]
    rownames(counts)   <- NULL
    rownames(bed_file) <- NULL
  }
  qc_info <- list(
    params = list(sample_qc = sample_qc, exon_qc = exon_qc, qc_zscore = qc_zscore,
                  outlier_sides = outlier_sides, exon_mad_quantile = exon_mad_quantile,
                  min_exon_mean = min_exon_mean, gc_extreme_filter = gc_extreme_filter),
    sample_qc = qc$sample_qc, excluded_samples = qc$excluded_samples,
    removed_targets = removed_targets,
    n_targets = list(mode = n_targets_mode, analysed = nrow(counts)),
    candidate_samples = md$candidate)

  results      <- list()
  ed_objects   <- list()
  no_result    <- list(calls = data.frame(), model = NULL, refs = character(0))

  # Column mapping shared by the real TSV write and the header-only TSV below,
  # so empty and non-empty runs always emit the same schema. ExomeDepth's
  # CNV.calls columns (start, end, chromosome, id, type, BF, reads.expected,
  # reads.observed, reads.ratio, start.p, end.p, nexons) keep their name unless
  # renamed here; `id` is ExomeDepth's internal target id and is dropped.
  cnv_rename_map <- c(start.p = "Start.p", end.p = "End.p", type = "Type", nexons = "Nexons",
                      start = "Start", end = "End", chromosome = "Chromosome",
                      reads.expected = "Reads.expected", reads.observed = "Reads.observed",
                      reads.ratio = "Reads.ratio")
  cnv_exclude_cols <- c("global_start", "global_end", "terminal_only", "id")
  # ExomeDepth CNV.calls schema, in constructor order; the pipeline appends
  # Gene, Confidence, Correlation, N.comp, Comparator.name, Genomic.ID and
  # prepends CNV.ID / Sample before writing.
  ed_call_columns <- c("start", "end", "chromosome", "id", "type", "BF",
                       "reads.expected", "reads.observed", "reads.ratio",
                       "start.p", "end.p", "nexons")
  cnv_tsv_columns <- c(
    "CNV.ID", "Sample",
    unname(ifelse(ed_call_columns %in% names(cnv_rename_map),
                  cnv_rename_map[ed_call_columns], ed_call_columns)),
    "Gene", "Confidence", "Correlation", "N.comp", "Comparator.name", "Genomic.ID",
    "Carriers", "Recurrent"
  )
  cnv_tsv_columns <- setdiff(cnv_tsv_columns, cnv_exclude_cols)

  for (sample in sample_names) {
    message("[INFO] Processing sample: ", sample)
    if (modechrom == "A") {
      ref_samples <- candidate_sample_names[candidate_sample_names != sample]
    } else {
      this_gender <- gender_map[[sample]]
      if (modechrom == "X") {
        # Same-sex references, restricted to the samples that are actually in
        # the coverage matrix and survived sample_qc (the sample table may list
        # more samples than were sequenced, and outliers must be used neither
        # as test nor as reference).
        ref_samples <- intersect(
          names(gender_map[gender_map == this_gender & names(gender_map) != sample]),
          candidate_sample_names)
      } else {
        ref_samples <- candidate_sample_names[candidate_sample_names != sample]
      }
    }

    if (length(ref_samples) == 0) {
      message("[WARNING] No reference samples available for ", sample, ". Skipping.")
      results[[sample]] <- no_result
      next
    }

    test_counts <- counts[, sample, drop = TRUE]
    ref_matrix  <- as.matrix(counts[, ref_samples, drop = FALSE])

    best_refs <- suppressWarnings(tryCatch({
      ExomeDepth::select.reference.set(
        test.counts      = test_counts,
        reference.counts = ref_matrix,
        bin.length       = (counts$end - counts$start + 1L),
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
      results[[sample]] <- no_result
      next
    }

    ref_counts <- rowSums(counts[, best_refs, drop = FALSE])
    message(sprintf("[DEBUG] %s: %d reference(s) selected from %d candidate(s); correlation with the pooled reference %.4f; first references: %s",
                    sample, length(best_refs), length(ref_samples),
                    suppressWarnings(as.numeric(stats::cor(counts[, sample], ref_counts))),
                    paste(utils::head(best_refs, 5L), collapse = ", ")))
    suppressWarnings({
      ed <- methods::new("ExomeDepth", test = test_counts, reference = ref_counts,
                         formula = formula, data = data)
      ed <- ExomeDepth::CallCNVs(
        x                      = ed,
        transition.probability = transition.probability,
        expected.CNV.length    = expected.CNV.length,
        chromosome             = counts$chromosome,
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
      
      raw_calls$global_start <- as.integer(raw_calls$start.p)
      raw_calls$global_end   <- as.integer(raw_calls$end.p)
      
      raw_calls <- flag_terminal_only_calls(raw_calls, bed_file)
      
      processed_calls <- tryCatch({
        score_cnv_confidence(
          raw_calls, bed_file, current_cor, num_refs,
          penalize_terminal_only = penalize_terminal_only,
          ...
        )
      }, error = function(e) {
        message("[WARNING] Confidence scoring failed for ", sample, ": ", e$message)
        # Keep the schema identical to scored calls so that rbind() across
        # samples and the report/VCF code do not break on a missing column.
        raw_calls$Gene       <- "Unknown"
        raw_calls$Confidence <- "LOW"
        raw_calls
      })
      
      if (remove_terminal_only && any(processed_calls$terminal_only, na.rm = TRUE)) {
        processed_calls <- processed_calls[!processed_calls$terminal_only, ]
      }
      
      if (nrow(processed_calls) > 0) {
        processed_calls <- add_within_gene_indices(processed_calls, bed_file)
        processed_calls <- cbind(Sample = sample, processed_calls)
        processed_calls$Correlation      <- current_cor
        processed_calls$N.comp           <- num_refs
        processed_calls$Comparator.name  <- paste(best_refs, collapse = ", ")

        # cnv_rename_map / cnv_exclude_cols are shared with the empty-run
        # header below so the two TSV schemas cannot drift apart.
        hit <- names(processed_calls) %in% names(cnv_rename_map)
        names(processed_calls)[hit] <- cnv_rename_map[names(processed_calls)[hit]]
        processed_calls$Genomic.ID <- paste0(processed_calls$Chromosome, ":",
                                              processed_calls$Start, "-", processed_calls$End)
      } else {
        processed_calls <- data.frame()
      }
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
  all_calls_list <- Filter(function(x) !is.null(x) && nrow(x) > 0, all_calls_list)
  models         <- lapply(results, function(x) x$model)
  names(models)  <- names(results)        
  models         <- models[sample_names]
  refs           <- lapply(results, function(x) x$refs)
  names(refs)    <- names(results)
  refs           <- refs[sample_names]
  
  cnv_calls <- if (length(all_calls_list) > 0)
    as.data.frame(data.table::rbindlist(all_calls_list, fill = TRUE)) else data.frame()
  if (nrow(cnv_calls) > 0) {
    # Cohort-level recurrence: the same call in many samples is a noisy target,
    # not a germline CNV (applied once, on the whole cohort).
    cnv_calls <- cnv_apply_recurrence(
      cnv_calls, list(sample = "Sample", chr = "Chromosome", start = "Start", end = "End",
                      type = "Type", confidence = "Confidence"),
      n_samples = length(sample_names), action = recurrence_action,
      min_samples = recurrence_min_samples, min_fraction = recurrence_min_fraction,
      min_overlap = recurrence_min_overlap)
  }
  if (nrow(cnv_calls) > 0) {
    cnv_calls  <- cbind(CNV.ID = seq_len(nrow(cnv_calls)), cnv_calls)
    tsv_cols     <- setdiff(colnames(cnv_calls), cnv_exclude_cols)
    utils::write.table(cnv_calls[, tsv_cols, drop = FALSE], file = output_file,
                      row.names = FALSE, sep = "\t", quote = FALSE)
    message("[INFO] Saved CNV calls to: ", output_file)
  } else {
    # A header-only TSV represents a successful run with zero calls. Opening
    # the file for writing also replaces any calls left by an earlier run.
    # The header is derived from the same rename/exclude mapping used for real
    # calls (cnv_rename_map / cnv_exclude_cols above), so an empty run emits
    # exactly the columns -- and column order -- of a run with calls.
    empty_cnv_calls <- as.data.frame(
      setNames(lapply(cnv_tsv_columns, function(x) character(0)), cnv_tsv_columns),
      check.names = FALSE
    )
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(empty_cnv_calls, file = output_file,
                       row.names = FALSE, sep = "\t", quote = FALSE)
    message("[INFO] No CNVs detected across any samples; wrote header-only CNV TSV: ",
            output_file)
  }

  qc_info$recurrence <- list(action = recurrence_action, min_samples = recurrence_min_samples,
                             min_fraction = recurrence_min_fraction, min_overlap = recurrence_min_overlap)
  if (save_ed_objects) {
    save(counts, cnv_calls, refs, models, bed_file, ed_objects, gender_map,
         sample_names, reference_sample_names, qc_info, modechrom, file = out_rdata)
  } else {
    save(counts, cnv_calls, refs, models, bed_file, gender_map,
         sample_names, reference_sample_names, qc_info, modechrom, file = out_rdata)
  }
  message("[INFO] Saved summary RData to: ", out_rdata)
  message("[INFO] END CNV calls")
  invisible(NULL)
}