#' Run the Full ECHO Pipeline
#'
#' Runs coverage extraction, QC metrics, CNV calling, plot generation, and optional VCF export.
#'
#' @param config_path Optional character string. Path to a \code{config.yaml} file.
#' @param ... Optional named overrides (see details).
#' @param vcf_output Character string. If \code{NULL}, a default VCF path inside the output
#'   directory is used.  Set to \code{FALSE} to skip VCF export. Default \code{NULL}.
#' @param save_ed_objects Logical. Save full ExomeDepth objects? Default \code{FALSE}.
#' @param report Logical. Generate interactive HTML reports? Default \code{TRUE}.
#' @param plots Logical. Generate per-CNV PDF plots? Default \code{FALSE}. Can also be set
#'   in the YAML config under \code{plots} or \code{settings$plots}.
#' @param vcf_per_sample Logical. Write a separate VCF file for each sample? Default \code{FALSE}.
#' @param reuse_coverage Logical. Skip the BAM counting when the coverage RData of a previous run
#'   exists and the BAMs, the BED, the FASTA and every counting setting (padding, sample names,
#'   genome, BED processing) are unchanged (checked with a fingerprint stored next to the RData).
#'   Re-running with other calling / QC / report options is then fast. Default \code{FALSE}.
#' @param sample_name_delim Character string. Delimiter(s) to split the filename.
#' @param sample_name_keep Character string. Specifies which parts to keep after splitting.
#' @param sample_name_collapse Character string. Separator for rejoining parts. Default NULL.
#' @param custom_sample_names Optional character vector. Custom sample names.
#' @param log_file Character string. Path to log file.
#' @param gene_field_index Integer or \code{NULL} (default). Which field to keep after
#'   splitting by exon_sep; \code{NULL} lets \code{process_bed_file()} decide.
#' @param sample_table Optional character string. Path to a CSV/TSV with sample_name / gender.
#' @param ref_bams Optional character string. Path to TSV with external reference BAM list.
#' @param panel_files Optional character vector. BED file paths for target restriction.
#'
#' @return Invisibly returns \code{TRUE} on success.
#' @export
echo <- function(config_path = NULL, vcf_output = NULL, save_ed_objects = FALSE,
                 report = TRUE,
                 plots = NULL,  # NULL means "use config"
                 vcf_per_sample = FALSE,
                 reuse_coverage = FALSE,
                 sample_name_delim = "\\.",
                 sample_name_keep = "1",
                 sample_name_collapse = NULL,
                 custom_sample_names = NULL,
                 log_file = NULL,
                 gene_field_index = NULL,
                 sample_table = NULL,
                 ref_bams = NULL,
                 panel_files = NULL,
                 ...) {
  args <- list(...)

  if (!is.null(config_path)) {
    cfg <- load_config(config_path)
    if (!is.null(cfg$bed_preprocess)) {
      safe_keys <- setdiff(names(cfg$bed_preprocess), c("input", "output", "settings"))
      for (name in safe_keys) {
        cfg[[name]] <- cfg$bed_preprocess[[name]]
      }
    }
  } else if (length(args) > 0) {
    cfg <- list(
      input = list(
        bamdir             = args$bamdir %||% "./data",
        bamfiles           = args$bamfiles,
        bed                = args$bed,
        fasta              = args$fasta,
        fasta_source       = args$fasta_source %||% "file",
        bsgenome_cache_dir = args$bsgenome_cache_dir,
        rbams              = args$rbams
      ),
      output = list(
        dir    = args$outdir %||% "./result",
        prefix = args$prefix %||% "ECHO"
      ),
      settings = list(
        modechrom                  = args$modechrom %||% "A",
        min_corr                   = args$min_corr %||% 0.98,
        min_cov                    = args$min_cov %||% 100,
        min_total_reads            = args$min_total_reads %||% 300000,
        max_exon_cv                = args$max_exon_cv %||% 0.5,
        transition_probability     = args$transition_prob %||% args$transition_probability %||% 1e-4,
        expected_CNV_length        = args$expected_CNV_length %||% 50000,
        n_bins_reduced             = args$n_bins_reduced %||% 10000,
        phi_bins                   = args$phi_bins %||% 1,
        formula                    = args$formula %||% "cbind(test, reference) ~ 1",
        score_high_corr            = args$score_high_corr %||% 0.985,
        score_med_corr             = args$score_med_corr %||% 0.95,
        score_high_refs            = args$score_high_refs %||% 3,
        score_med_refs             = args$score_med_refs %||% 2,
        score_low_ratio_low        = args$score_low_ratio_low %||% 0.75,
        score_low_ratio_high       = args$score_low_ratio_high %||% 1.25,
        score_high_ratio_low       = args$score_high_ratio_low %||% 0.70,
        score_high_ratio_high      = args$score_high_ratio_high %||% 1.30,
        score_med_ratio_low        = args$score_med_ratio_low %||% 0.60,
        score_med_ratio_high       = args$score_med_ratio_high %||% 1.40,
        score_low_confidence_genes = args$score_low_confidence_genes %||%
          c("PMS2", "SMN1", "CYP2D6", "HBA1", "HBA2", "STRC", "CYP21A2", "GBA1", "CFTR"),
        sample_qc                  = args$sample_qc %||% TRUE,
        exon_qc                    = args$exon_qc %||% TRUE,
        qc_zscore                  = args$qc_zscore %||% 3,
        outlier_sides              = args$outlier_sides %||% "upper",
        recurrence_action          = args$recurrence_action %||% "flag",
        recurrence_min_samples     = args$recurrence_min_samples %||% 3L,
        recurrence_min_fraction    = args$recurrence_min_fraction %||% 0.25,
        recurrence_min_overlap     = args$recurrence_min_overlap %||% 0.5,
        exon_mad_quantile          = args$exon_mad_quantile %||% 0.90,
        gc_extreme_filter          = args$gc_extreme_filter %||% c(0.15, 0.85),
        min_exon_mean              = args$min_exon_mean %||% 20,
        pad_terminal_exons         = args$pad_terminal_exons %||% 0,
        remove_terminal_only       = args$remove_terminal_only %||% FALSE,
        penalize_terminal_only     = args$penalize_terminal_only %||% FALSE
      ),
      bed_process            = args$bed_process %||% "NO",
      refseqgene             = args$refseqgene,
      transcripts_file       = args$transcripts_file,
      unknown_gene           = args$unknown_gene %||% FALSE,
      gene_list_restrict     = args$gene_list_restrict,
      chr_list_restrict      = args$chr_list_restrict,
      exon_sep               = args$exon_sep,
      gene_name_keep         = args$gene_name_keep,
      region_numbering_mode  = args$region_numbering_mode %||% "bed_text",
      customexon             = args$customexon %||% FALSE,
      list_genes             = args$list_genes,
      genes_file             = args$genes_file,
      genome_version         = args$genome_version,
      bed_zero_based         = args$bed_zero_based %||% TRUE,
      skip_invalid_intervals = args$skip_invalid_intervals %||% TRUE,
      off_target_pattern     = args$off_target_pattern %||% "^HorsROI",
      off_target_handling    = args$off_target_handling %||% "na",
      keep_hors_roi          = args$keep_hors_roi %||% TRUE,
      sample_name_collapse   = args$sample_name_collapse %||% sample_name_collapse,
      panel_files            = args$panel_files %||% panel_files,
      auto_regen             = args$auto_regen %||% TRUE
    )
  } else {
    stop("[ERROR] Either config_path or pipeline parameters must be provided.")
  }

  # A formal argument overrides the YAML/variadic value and must be used by
  # coverage extraction as well as by the report.
  if (!is.null(ref_bams)) cfg$input$rbams <- ref_bams

  # Same plumbing for sample_table as for ref_bams above: a formal argument
  # overrides the YAML, and the YAML key is accepted at the top level or under
  # input: (mirroring ref_bams / input$rbams). Without this, a sample_table
  # set only in the YAML was silently ignored and mode X/Y runs failed later
  # with "[ERROR] sample_table required for mode X/Y".
  if (is.null(sample_table)) {
    sample_table <- cfg$sample_table %||% cfg$input$sample_table
  }

  # Fill in settings the config may omit. Without this a missing key reached
  # the callees as NULL (e.g. `modechrom == "A"` -> "argument is of length zero",
  # `min_corr = NULL` silently disabling the correlation QC).
  settings_defaults <- list(
    modechrom = "A", min_corr = 0.98, min_cov = 100, min_total_reads = 300000,
    max_exon_cv = 0.5, transition_probability = 1e-4, expected_CNV_length = 50000,
    n_bins_reduced = 10000, phi_bins = 1, formula = "cbind(test, reference) ~ 1"
  )
  if (is.null(cfg$settings)) cfg$settings <- list()
  for (nm in names(settings_defaults)) {
    if (is.null(cfg$settings[[nm]])) cfg$settings[[nm]] <- settings_defaults[[nm]]
  }
  cfg$settings$modechrom <- toupper(as.character(cfg$settings$modechrom)[1])
  if (!cfg$settings$modechrom %in% c("A", "X", "Y")) {
    stop("[ERROR] settings$modechrom must be 'A' (autosomes), 'X' or 'Y' (got '",
         cfg$settings$modechrom, "').")
  }
  cfg$output$prefix <- cfg$output$prefix %||% "ECHO"

  if (is.null(plots)) {
    plots <- cfg$plots %||% cfg$settings$plots %||% FALSE
  }

  if (!is.null(cfg$gene_field_index)) {
    gene_field_index <- cfg$gene_field_index
  } else if (!is.null(cfg$bed_preprocess$gene_field_index)) {
    gene_field_index <- cfg$bed_preprocess$gene_field_index
  }

  # file.path(getwd(), "/abs/dir") produced "<wd>//abs/dir": only prepend the
  # working directory to relative paths.
  out_dir <- cfg$output$dir %||% "./result"
  out_is_absolute <- grepl("^(/|~|[A-Za-z]:[/\\\\]|\\\\\\\\)", out_dir)
  cfg$output$dir <- normalizePath(
    if (out_is_absolute) path.expand(out_dir) else file.path(getwd(), out_dir),
    mustWork = FALSE)
  if (!dir.exists(cfg$output$dir)) {
    dir.create(cfg$output$dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (is.null(log_file)) {
    log_file <- file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_pipeline.log"))
  }
  log_dir <- dirname(log_file)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  # ONE logger for the whole run (shared with CANOPE): timestamped lines in the
  # file AND on the console, the file is truncated at start so the report only
  # shows this run, and lg$capture() routes everything the callees emit with
  # message() / warning() into the same file.
  lg <- cnv_logger(log_file, tag = "ECHO", level = cfg$settings$log_level %||% cfg$log_level %||% "INFO",
                    total_steps = NA_integer_)
  log_msg <- function(msg, type = "INFO") lg$log(msg, type)
  lg$header(config = if (is.null(config_path)) "command-line parameters" else config_path,
            packages = c("ExomeDepth", "data.table", "GenomicRanges"))
  log_msg(paste("Output directory:", basename(cfg$output$dir)))

  auto_regen <- cfg$auto_regen %||% TRUE
  if (auto_regen && !is.null(cfg$input$bed) && is.null(cfg$bed_process)) {
    cfg$bed_process <- "REGEN"
    log_msg("BED preprocessing auto-set to 'REGEN'")
  }
  if (is.null(cfg$genome_version) && !is.null(cfg$settings$genome_version)) {
    cfg$genome_version <- cfg$settings$genome_version
  }
  # config.yaml documents genome_version under `input:` (option B) as well.
  if (is.null(cfg$genome_version) && !is.null(cfg$input$genome_version)) {
    cfg$genome_version <- cfg$input$genome_version
  }
    if (is.null(cfg$genome_version)) {
    source <- cfg$input$fasta_source %||% "file"
    if (is.character(source) && length(source) == 1L &&
        grepl("^BSgenome\\.", source, ignore.case = TRUE)) {
      if (grepl("(^|\\.)(hg19|grch37)($|\\.)", source, ignore.case = TRUE)) {
        cfg$genome_version <- "hg19"
      } else if (grepl("(^|\\.)(hg38|grch38)($|\\.)", source, ignore.case = TRUE)) {
        cfg$genome_version <- "hg38"
      } else {
        stop("[ERROR] Cannot infer genome_version from BSgenome package '",
             source, "'. Set genome_version explicitly (supported: hg19, hg38) ",
             "or use fasta_source = 'file'.")
      }
    } else if (!is.null(cfg$input$fasta)) {
      fname <- basename(cfg$input$fasta)
      if (grepl("hg19|grch37", fname, ignore.case = TRUE)) {
        cfg$genome_version <- "hg19"
      } else if (grepl("hg38|grch38", fname, ignore.case = TRUE)) {
        cfg$genome_version <- "hg38"
      } else {
        stop("[ERROR] Cannot infer genome_version from FASTA filename '",
             fname, "'. Set genome_version explicitly in the config.")
      }
    } else {
      stop("[ERROR] genome_version could not be inferred (no FASTA path, ",
           "no BSgenome package, no genome_version setting).")
    }
  }

  if (is.character(cfg$genome_version) && length(cfg$genome_version) == 1L) {
    cfg$genome_version <- switch(tolower(cfg$genome_version),
                                 grch37 = "hg19", grch38 = "hg38", cfg$genome_version)
  }

  ok <- FALSE
  on.exit(lg$footer(ok), add = TRUE)

  tryCatch({
    lg$capture({
      log_msg(paste("Plots enabled:", plots))
      log_msg(paste("BED processing mode:", cfg$bed_process))
      log_msg(paste("Genome version:", cfg$genome_version))
      log_msg(paste("Terminal-exon padding:", paste0(cfg$settings$pad_terminal_exons %||% 0, " bp")))
      log_msg(paste("Remove terminal-only calls:", cfg$settings$remove_terminal_only %||% FALSE))
      log_msg(paste("Penalize terminal-only calls:", cfg$settings$penalize_terminal_only %||% FALSE))
      if (!is.null(cfg$bed_process) && cfg$bed_process != "NO") {
        log_msg(paste0("Off-target region handling: pattern = '", cfg$off_target_pattern %||% "^HorsROI",
                       "', handling = '", cfg$off_target_handling %||% "na", "'"))
        log_msg(paste("Keep HorsROI intervals:", cfg$keep_hors_roi %||% TRUE))
      }

      paths <- list(
        rdata   = file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_coverage.RData")),
        metrics = file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_QC_metrics.tsv")),
        cnvs    = file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_CNV_calls.tsv")),
        summary = file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_summary.RData")),
        plots   = file.path(cfg$output$dir, "Plots")
      )

      if (!is.null(args$rdata))   paths$rdata   <- args$rdata
      if (!is.null(args$metrics)) paths$metrics <- args$metrics
      if (!is.null(args$cnvs))    paths$cnvs    <- args$cnvs
      if (!is.null(args$summary)) paths$summary <- args$summary

      lapply(paths[!names(paths) %in% "plots"], function(p) {
        if (grepl("\\.[a-zA-Z]+$", p)) {
          dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
        } else {
          dir.create(p, showWarnings = FALSE, recursive = TRUE)
        }
      })

      if (is.null(vcf_output)) {
        vcf_output <- file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_CNVs.vcf"))
      } else if (identical(vcf_output, FALSE)) {
        vcf_output <- NULL
      }

      stop_if_not_file(cfg$input$bed, "[ERROR] BED file missing")

      fasta_source <- cfg$input$fasta_source %||% "file"
      if (identical(fasta_source, "file")) {
        stop_if_not_file(cfg$input$fasta, "[ERROR] FASTA file missing")
      } else {
        log_msg(paste0("Reference sequence source: fasta_source = '", fasta_source,
                        "' (genome_version = ", cfg$genome_version, "); no local FASTA required."))
      }

      bed_zero_based_run <- cfg$bed_zero_based %||% TRUE
      bed_issues <- validate_bed_regions(cfg$input$bed, verbose = TRUE,
                                         bed_zero_based = bed_zero_based_run)
      fatal_bed_issues <- setdiff(names(bed_issues), "duplicate_regions")
      if (length(fatal_bed_issues)) {
        stop("[ERROR] BED validation failed (", paste(fatal_bed_issues, collapse = ", "),
             "). Correct the BED file before continuing.")
      }

      if (!is.null(cfg$bed_process) && cfg$bed_process != "NO") {
        log_msg(paste("Preprocessing BED file using mode:", cfg$bed_process))
        processed_bed <- file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_targets.bed"))
        process_bed_file(
          input_bed             = cfg$input$bed,
          output_bed            = processed_bed,
          bed_process           = cfg$bed_process,
          bed_zero_based        = cfg$bed_zero_based %||% TRUE,
          refseqgene            = cfg$refseqgene %||% NULL,
          transcripts_file      = cfg$transcripts_file %||% NULL,
          unknown_gene          = cfg$unknown_gene %||% FALSE,
          gene_list_restrict    = cfg$gene_list_restrict %||% NULL,
          exon_sep              = cfg$exon_sep %||% NULL,
          gene_name_collapse    = cfg$gene_name_collapse %||% "_",
          customexon            = cfg$customexon %||% FALSE,
          auto_exon_number      = cfg$auto_exon_number %||% TRUE,
          region_numbering_mode = cfg$region_numbering_mode %||% "bed_text",
          gene_name_keep        = cfg$gene_name_keep %||% NULL,
          list_genes            = cfg$list_genes %||% NULL,
          genes_file            = cfg$genes_file %||% NULL,
          panel_files           = cfg$panel_files %||% NULL,
          genome_version        = cfg$genome_version,
          txdb                  = cfg$txdb %||% NULL,
          gene_field_index      = gene_field_index,
          off_target_pattern    = cfg$off_target_pattern %||% "^HorsROI",
          off_target_handling   = cfg$off_target_handling %||% "na",
          keep_hors_roi          = cfg$keep_hors_roi %||% TRUE
        )
        cfg$input$bed <- processed_bed
        # process_bed_file() always writes a 0-based BED (it subtracts 1 from
        # 1-based input), so the coverage step must read it as 0-based
        # whatever bed_zero_based said about the ORIGINAL file. Reusing the
        # original flag shifted every target by 1 bp when bed_zero_based was
        # false.
        bed_zero_based_run <- TRUE
        log_msg(paste("Processed BED saved to:", processed_bed))
      } else {
        log_msg("BED preprocessing skipped (bed_process = 'NO').")
      }

      steps_total <- 3L +
                     as.integer(cfg$settings$pca_plot %||% TRUE) +
                     as.integer(plots) +
                     as.integer(report) +
                     as.integer(!is.null(vcf_output) || vcf_per_sample)
      lg$set_total(steps_total)
      # a failed step stops the run (`on_error = "stop"`); optional outputs
      # (PCA, plots, report, VCF) only warn so the calls are never lost
      run_step <- function(step_name, expr, on_error = "stop") lg$step(step_name, expr, on_error = on_error)

      cov_fp <- NULL; cov_fp_file <- paste0(paths$rdata, ".fingerprint")
      if (isTRUE(reuse_coverage)) {
        cov_fp <- echo_coverage_fingerprint(
          files = c(echo_input_bams(cfg$input$bamfiles, cfg$input$bamdir, cfg$input$rbams),
                    cfg$input$bed, if (identical(fasta_source, "file")) cfg$input$fasta, cfg$input$rbams),
          settings = list(pad = cfg$settings$pad_terminal_exons %||% 0, zero_based = bed_zero_based_run,
                          fasta_source = fasta_source, genome = cfg$genome_version,
                          delim = sample_name_delim, keep = sample_name_keep,
                          collapse = cfg$sample_name_collapse %||% sample_name_collapse,
                          custom = custom_sample_names, skip_invalid = cfg$skip_invalid_intervals %||% TRUE,
                          bed_process = cfg$bed_process))
      }
      run_step("Extracting BAM coverage", {
        if (!is.null(cov_fp) && file.exists(paths$rdata) && file.exists(cov_fp_file) &&
            identical(readLines(cov_fp_file, n = 1L, warn = FALSE), cov_fp)) {
          log_msg(paste("Reusing", basename(paths$rdata), "(BAMs, BED, FASTA and counting settings unchanged since it was written)"))
        } else {
          if (!is.null(cov_fp) && file.exists(paths$rdata))
            log_msg("reuse_coverage: inputs or counting settings changed -> counting again")
        run_bam_coverage(
          bamfiles             = cfg$input$bamfiles,
          bamdir               = cfg$input$bamdir,
          bed                  = cfg$input$bed,
          fasta                = cfg$input$fasta,
          fasta_source         = fasta_source,
          genome_version       = cfg$genome_version,
          bsgenome_cache_dir   = cfg$input$bsgenome_cache_dir,
          rbams                = cfg$input$rbams,
          data_out             = paths$rdata,
          verbose              = TRUE,
          sample_name_delim    = sample_name_delim,
          sample_name_keep     = sample_name_keep,
          sample_name_collapse = cfg$sample_name_collapse %||% sample_name_collapse,
          custom_sample_names  = custom_sample_names,
          bed_zero_based       = bed_zero_based_run,
          skip_invalid_intervals = cfg$skip_invalid_intervals %||% TRUE,
          pad_terminal_exons   = cfg$settings$pad_terminal_exons %||% 0
        )
          if (!is.null(cov_fp)) writeLines(cov_fp, cov_fp_file)
        }
      })

      st <- cfg$settings
      if (st$pca_plot %||% TRUE) {
        run_step("Generating PCA plot", {
          objs         <- load_rdata(paths$rdata, required = c("counts", "sample_names"))
          group_info   <- NULL
          if (!is.null(sample_table) && file.exists(sample_table)) {
            sample_df  <- utils::read.table(sample_table, header = TRUE, sep = "\t")
            group_info <- sample_df$gender[match(objs$sample_names, sample_df$sample_name)]
          }
          plot_coverage_pca(objs$counts, objs$sample_names,
                            output_pdf = file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_PCA.pdf")),
                            color_by = group_info)
        }, on_error = "warn")
      }

      # the SAME QC arguments are given to the QC-metrics step and to the caller,
      # so the table in the report describes exactly what was excluded
      qc_args <- list(
        sample_table = sample_table,
        sample_qc = st$sample_qc %||% TRUE, exon_qc = st$exon_qc %||% TRUE,
        qc_zscore = st$qc_zscore %||% 3, outlier_sides = st$outlier_sides %||% "upper",
        exon_mad_quantile = st$exon_mad_quantile %||% 0.90,
        gc_extreme_filter = st$gc_extreme_filter %||% c(0.15, 0.85),
        min_exon_mean = st$min_exon_mean %||% 20)

      run_step("Running QC metrics", {
        do.call(run_qc_metrics, c(list(
          rdata_file = paths$rdata, output_file = paths$metrics,
          min_corr = st$min_corr, min_cov = st$min_cov, min_total_reads = st$min_total_reads,
          max_exon_cv = st$max_exon_cv, modechrom = st$modechrom), qc_args))
      }, on_error = "warn")

      run_step("Calling CNVs", {
        cnv_args <- c(list(
          rdata_file             = paths$rdata,
          output_file            = paths$cnvs,
          out_rdata              = paths$summary,
          transition.probability = st$transition_probability,
          expected.CNV.length    = st$expected_CNV_length,
          n.bins.reduced         = st$n_bins_reduced,
          phi.bins               = st$phi_bins,
          formula                = st$formula,
          data                   = NULL,
          save_ed_objects        = save_ed_objects,
          modechrom              = st$modechrom,
          remove_terminal_only   = st$remove_terminal_only %||% FALSE,
          penalize_terminal_only = st$penalize_terminal_only %||% FALSE,
          recurrence_action      = st$recurrence_action %||% "flag",
          recurrence_min_samples = st$recurrence_min_samples %||% 3L,
          recurrence_min_fraction = st$recurrence_min_fraction %||% 0.25,
          recurrence_min_overlap = st$recurrence_min_overlap %||% 0.5),
          qc_args, st[grepl("^score_", names(st))])
        do.call(run_cnv_calling, cnv_args)
      })

      if (plots) {
        run_step("Generating plots", {
          generate_plots(rdata_file = paths$summary, output_dir = paths$plots,
                         modechrom = st$modechrom, prefix = cfg$output$prefix,
                         gene_gap = st$plot_gene_gap %||% 1,
                         off_target_pattern = cfg$off_target_pattern %||% "^HorsROI")
        }, on_error = "warn")
      } else {
        log_msg("Plot generation skipped (plots = FALSE).")
      }

      if (report) {
        run_step("Generating HTML report", {
          generate_report(summary_rdata = paths$summary, qc_metrics_file = paths$metrics,
                          output_dir = cfg$output$dir, settings = st, config = cfg,
                          sample_table = sample_table, ref_bams = ref_bams %||% cfg$input$rbams,
                          log_file = log_file, pdf_output = st$pdf_output %||% FALSE,
                          prefix = cfg$output$prefix, coverage_rdata = paths$rdata)
        }, on_error = "warn")
      } else {
        log_msg("Report generation skipped (report = FALSE).")
      }

      if (!is.null(vcf_output) || vcf_per_sample) {
        run_step("Exporting CNVs to VCF", {
          if (!file.exists(paths$summary)) stop("Summary RData not found, cannot export the VCF")
          env <- new.env(parent = emptyenv())
          load(paths$summary, envir = env)
          calls <- env$cnv_calls %||% data.frame()
          samples_all <- unique(as.character(env$sample_names %||% character()))
          samples_all <- samples_all[!is.na(samples_all) & nzchar(samples_all)]
          if (!is.null(vcf_output)) {
            export_cnvs_to_vcf(calls, vcf_output, sample_name = NULL,
                               reference = cfg$genome_version, sample_names = samples_all)
          }
          if (vcf_per_sample) {
            # one VCF per analysed sample, header-only when it has no call; old
            # per-sample files of this prefix are removed first
            plots_root <- file.path(cfg$output$dir, "Plots")
            if (dir.exists(plots_root)) {
              old <- list.files(plots_root, pattern = "\\.vcf$", recursive = TRUE, full.names = TRUE)
              old <- old[startsWith(basename(old), paste0("ECHO_", cfg$output$prefix, "_"))]
              if (length(old)) unlink(old)
            }
            for (s in samples_all) {
              folder <- file.path(plots_root, sanitize_filename(s))
              dir.create(folder, recursive = TRUE, showWarnings = FALSE)
              export_cnvs_to_vcf(calls, file.path(folder, paste0("ECHO_", cfg$output$prefix, "_", s, ".vcf")),
                                 sample_name = s, reference = cfg$genome_version, sample_names = s)
            }
          }
        }, on_error = "warn")
      } else {
        log_msg("VCF export skipped.")
      }
    })
    ok <- TRUE
    invisible(TRUE)
  }, error = function(e) {
    lg$error(paste("FATAL:", conditionMessage(e)))
    stop(e)
  })
}
