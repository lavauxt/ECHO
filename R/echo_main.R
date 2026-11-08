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
#' @param sample_name_delim Character string. Delimiter(s) to split the filename.
#' @param sample_name_keep Character string. Specifies which parts to keep after splitting.
#' @param sample_name_collapse Character string. Separator for rejoining parts. Default NULL.
#' @param custom_sample_names Optional character vector. Custom sample names.
#' @param log_file Character string. Path to log file.
#' @param gene_field_index Integer. Which field to keep after splitting by exon_sep.
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
                 sample_name_delim = "\\.",
                 sample_name_keep = "1",
                 sample_name_collapse = NULL,
                 custom_sample_names = NULL,
                 log_file = NULL,
                 gene_field_index = 1,
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

  # log_msg() writes to the file itself; the message handler below captures
  # everything the pipeline functions emit with message() (QC exclusions, BED
  # warnings, ...) so that the log -- and the report built from it -- is
  # complete. `in_log_msg` stops log_msg()'s own message() being logged twice.
  in_log_msg <- FALSE
  log_msg <- function(msg, type = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    formatted <- paste0("[", type, "] ", timestamp, " ", msg)
    in_log_msg <<- TRUE
    on.exit(in_log_msg <<- FALSE, add = TRUE)
    message(formatted)
    cat(formatted, "\n", file = log_file, append = TRUE, sep = "")
  }

  log_msg("ECHO pipeline started")
  log_msg(paste("Configuration:", ifelse(is.null(config_path), "command-line parameters", config_path)))
  log_msg(paste("Output directory:", basename(cfg$output$dir)))
  log_msg(paste("Log file:", basename(log_file)))

  si <- utils::sessionInfo()
  log_msg(paste0("R version: ", si$R.version$version.string))
  log_msg(paste0("Platform: ", si$platform))
  log_msg("Packages:")
  for (pkg in names(si$otherPkgs)) {
    log_msg(paste0("   ", pkg, ": ", si$otherPkgs[[pkg]]$Version))
  }
  log_msg(" ")

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

  old_warn <- options(warn = 1)
  on.exit(options(warn = old_warn$warn), add = TRUE)

  tryCatch({
    withCallingHandlers({
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

      # The PCA block below runs outside run_step(); count it too, or the
      # "Step k/N" log counter under-reports N whenever PCA is enabled.
      steps_total   <- 3L +  
                       as.integer(cfg$settings$pca_plot %||% TRUE) +
                       as.integer(plots) +
                       as.integer(report) +
                       as.integer(!is.null(vcf_output) || vcf_per_sample)
      step_current  <- 0L
      run_step <- function(step_name, expr) {
        step_current <<- step_current + 1L
        log_msg(paste0("Step ", step_current, "/", steps_total, ": ", step_name, "..."))
        result <- expr
        log_msg(paste0("Step ", step_current, "/", steps_total, " completed."))
        invisible(result)
      }

      run_step("Extracting BAM coverage", {
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
      })

      if (cfg$settings$pca_plot %||% TRUE) {
        # PCA is counted in steps_total, so give it its own "Step k/N" log
        # lines. It used to run unlogged between steps 1 and 2, leaving the
        # counter one short (a 7-step run ended at "Step 6/7 completed").
        run_step("Generating PCA plot", {
          tryCatch({
            objs         <- load_rdata(paths$rdata, required = c("counts", "sample_names"))
            counts       <- objs$counts
            sample_names <- objs$sample_names
            pca_file     <- file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_PCA.pdf"))
            group_info   <- NULL
            if (!is.null(sample_table) && file.exists(sample_table)) {
              sample_df  <- read.table(sample_table, header = TRUE, sep = "\t")
              group_info <- sample_df$gender[match(sample_names, sample_df$sample_name)]
            }
            plot_coverage_pca(counts, sample_names, output_pdf = pca_file, color_by = group_info)
          }, error = function(e) log_msg(paste("PCA plot step failed:", conditionMessage(e)), "WARNING")
          )
        })
      }

      run_step("Running QC metrics", {
        tryCatch(
          run_qc_metrics(
            rdata_file      = paths$rdata,
            output_file     = paths$metrics,
            min_corr        = cfg$settings$min_corr,
            min_cov         = cfg$settings$min_cov,
            min_total_reads = cfg$settings$min_total_reads,
            max_exon_cv     = cfg$settings$max_exon_cv,
            modechrom       = cfg$settings$modechrom
          ),
          error = function(e) log_msg(paste("QC metrics step failed:", conditionMessage(e)), "WARNING")
        )
      })

      run_step("Calling CNVs", {
        cnv_args <- cfg$settings[grepl("^score_", names(cfg$settings))]
        cnv_args <- c(list(
          rdata_file             = paths$rdata,
          output_file            = paths$cnvs,
          out_rdata              = paths$summary,
          transition.probability = cfg$settings$transition_probability,
          expected.CNV.length    = cfg$settings$expected_CNV_length,
          n.bins.reduced         = cfg$settings$n_bins_reduced,
          phi.bins               = cfg$settings$phi_bins,
          formula                = cfg$settings$formula,
          data                   = NULL,
          save_ed_objects        = save_ed_objects,
          modechrom              = cfg$settings$modechrom,
          sample_table           = sample_table,
          sample_qc              = cfg$settings$sample_qc %||% TRUE,
          exon_qc                = cfg$settings$exon_qc %||% TRUE,
          qc_zscore              = cfg$settings$qc_zscore %||% 3,
          exon_mad_quantile      = cfg$settings$exon_mad_quantile %||% 0.90,
          gc_extreme_filter      = cfg$settings$gc_extreme_filter %||% c(0.15, 0.85),
          min_exon_mean          = cfg$settings$min_exon_mean %||% 20,
          remove_terminal_only   = cfg$settings$remove_terminal_only %||% FALSE,
          penalize_terminal_only = cfg$settings$penalize_terminal_only %||% FALSE
        ), cnv_args)
        do.call(run_cnv_calling, cnv_args)
      })

      if (plots) {
        run_step("Generating plots", {
          tryCatch(
            generate_plots(
              rdata_file = paths$summary,
              output_dir = paths$plots,
              modechrom  = cfg$settings$modechrom,
              prefix     = cfg$output$prefix,
              log_file   = NULL,  # its message() output is captured by the handler below
              gene_gap   = cfg$settings$plot_gene_gap %||% 1,
              off_target_pattern = cfg$off_target_pattern %||% "^HorsROI"
            ),
            error = function(e) log_msg(paste("Plot generation failed:", conditionMessage(e)), "WARNING")
          )
        })
      } else {
        log_msg("Plot generation skipped (plots = FALSE).")
      }

      if (report) {
        run_step("Generating HTML report", {
          tryCatch(
            generate_report(
              summary_rdata   = paths$summary,
              qc_metrics_file = paths$metrics,
              output_dir      = cfg$output$dir,
              settings        = cfg$settings,
              config          = cfg,
              sample_table    = sample_table,
              ref_bams        = ref_bams %||% cfg$input$rbams,
              log_file        = log_file,
              pdf_output      = cfg$settings$pdf_output %||% FALSE,
              prefix          = cfg$output$prefix,
              coverage_rdata  = paths$rdata
            ),
            error = function(e) log_msg(paste("Report generation failed:", conditionMessage(e)), "WARNING")
          )
        })
      } else {
        log_msg("Report generation skipped (report = FALSE).")
      }

      if (!is.null(vcf_output) || vcf_per_sample) {
        run_step("Exporting CNVs to VCF", {
          tryCatch({
          if (file.exists(paths$summary)) {
            summary_env <- new.env(parent = emptyenv())
            load(paths$summary, envir = summary_env)
            cnv_calls_local <- summary_env$cnv_calls
            if (is.null(cnv_calls_local)) cnv_calls_local <- data.frame()
            analysis_samples_local <- as.character(
              summary_env$sample_names %||% character(0)
            )
            analysis_samples_local <- unique(
              analysis_samples_local[!is.na(analysis_samples_local) &
                                       nzchar(analysis_samples_local)]
            )
            n_calls_local <- nrow(cnv_calls_local)

            if (!is.null(vcf_output)) {
              export_cnvs_to_vcf(
                cnv_calls_local, vcf_output, sample_name = NULL,
                reference = cfg$genome_version,
                sample_names = analysis_samples_local
              )
              if (n_calls_local > 0) {
                log_msg(paste("Combined VCF written to:", basename(vcf_output)))
              } else {
                log_msg(paste("Header-only combined VCF written (zero CNVs):",
                              basename(vcf_output)))
              }
            }

            if (vcf_per_sample) {
              # Remove old per-sample VCFs for this ECHO prefix first, including
              # files for samples excluded from the current run. Then regenerate
              # one VCF per current analysis sample, even if it has no records.
              plots_root <- file.path(cfg$output$dir, "Plots")
              if (dir.exists(plots_root)) {
                old_vcfs <- list.files(plots_root, pattern = "\\.vcf$",
                                       recursive = TRUE, full.names = TRUE)
                prefix_start <- paste0("ECHO_", cfg$output$prefix, "_")
                old_vcfs <- old_vcfs[startsWith(basename(old_vcfs), prefix_start)]
                if (length(old_vcfs) && !all(file.remove(old_vcfs))) {
                  stop("Could not remove one or more stale per-sample VCF files under: ",
                       plots_root)
                }
              }

              for (s in analysis_samples_local) {
                sample_folder <- file.path(plots_root, sanitize_filename(s))
                if (!dir.exists(sample_folder))
                  dir.create(sample_folder, recursive = TRUE, showWarnings = FALSE)
                vcf_file <- file.path(sample_folder,
                                      paste0("ECHO_", cfg$output$prefix, "_", s, ".vcf"))
                export_cnvs_to_vcf(
                  cnv_calls_local, vcf_file, sample_name = s,
                  reference = cfg$genome_version,
                  sample_names = s
                )
                if (n_calls_local > 0 && "Sample" %in% names(cnv_calls_local) &&
                    any(cnv_calls_local$Sample == s, na.rm = TRUE)) {
                  log_msg(paste("Per-sample VCF written:", basename(vcf_file)))
                } else {
                  log_msg(paste("Header-only per-sample VCF written (zero CNVs):",
                                basename(vcf_file)))
                }
              }
            }
            if (n_calls_local == 0) {
              log_msg("No CNV calls detected; valid header-only VCF output(s) were written.")
            }
          } else {
            log_msg("Summary RData not found – cannot export VCF.", "WARNING")
          }
          }, error = function(e) log_msg(paste("VCF export failed:", conditionMessage(e)), "WARNING")
        )
        })
      } else {
        log_msg("VCF export skipped.")
      }

      log_msg("Pipeline finished successfully!")
    }, message = function(m) {
      if (in_log_msg) return(invisible(NULL))
      txt <- sub("\n+$", "", conditionMessage(m))
      if (!nzchar(trimws(txt))) return(invisible(NULL))
      ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      line <- if (grepl("^\\[(INFO|WARNING|ERROR)\\] [0-9]{4}-[0-9]{2}-[0-9]{2} ", txt)) {
        txt
      } else if (grepl("^\\[(INFO|WARNING|ERROR)\\]", txt)) {
        paste0("[", sub("^\\[([A-Z]+)\\].*$", "\\1", txt), "] ", ts, " ",
               sub("^\\[[A-Z]+\\]\\s*", "", txt))
      } else {
        paste0("[INFO] ", ts, " ", txt)
      }
      cat(line, "\n", file = log_file, append = TRUE, sep = "")
    }, warning = function(w) {
      if (grepl("sequence levels not in the other", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
      log_msg(sub("^\\[WARNING\\]\\s*", "", conditionMessage(w)), "WARNING")
      invokeRestart("muffleWarning")
    })
    invisible(TRUE)
  }, error = function(e) {
    log_msg(paste("FATAL:", e$message), "ERROR")
    stop(e)
  })
}