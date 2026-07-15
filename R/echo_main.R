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

  # ---- Load or build configuration -----------------------------------------
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
        min_total_reads            = args$min_total_reads %||% 5e6,
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
      genome_version         = args$genome_version %||% "hg19",
      bed_zero_based         = args$bed_zero_based %||% TRUE,
      skip_invalid_intervals = args$skip_invalid_intervals %||% TRUE,
      off_target_pattern     = args$off_target_pattern %||% "^HorsROI",
      off_target_handling    = args$off_target_handling %||% "na",
      sample_name_collapse   = args$sample_name_collapse %||% sample_name_collapse,
      panel_files            = args$panel_files %||% panel_files
    )
  } else {
    stop("[ERROR] Either config_path or pipeline parameters must be provided.")
  }

  # ---- Handle plots argument ----------------------------------------------
  if (is.null(plots)) {
    plots <- cfg$plots %||% cfg$settings$plots %||% FALSE
  }

  # ---- Override gene_field_index from config if present --------------------
  if (!is.null(cfg$gene_field_index)) {
    gene_field_index <- cfg$gene_field_index
  } else if (!is.null(cfg$bed_preprocess$gene_field_index)) {
    gene_field_index <- cfg$bed_preprocess$gene_field_index
  }

  # ---- Set up output directory and log file --------------------------------
  cfg$output$dir <- normalizePath(file.path(getwd(), cfg$output$dir), mustWork = FALSE)
  if (!dir.exists(cfg$output$dir)) {
    dir.create(cfg$output$dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (is.null(log_file)) {
    log_file <- file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_pipeline.log"))
  }
  log_dir <- dirname(log_file)
  if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- Define log_msg early (available for all subsequent code) -----------
  log_msg <- function(msg, type = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    formatted <- paste0("[", type, "] ", timestamp, " ", msg)
    message(formatted)
    cat(formatted, "\n", file = log_file, append = TRUE)
  }

  # ---- Write session info to log ------------------------------------------
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

  # ---- Auto‑set bed_process to "REGEN" ONLY if not explicitly defined ----
  auto_regen <- cfg$auto_regen %||% TRUE
  if (auto_regen && !is.null(cfg$input$bed) && is.null(cfg$bed_process)) {
    cfg$bed_process <- "REGEN"
    log_msg("BED preprocessing auto‑set to 'REGEN' (using TxDb/org.Hs.eg.db)")
  }

  # ---- Infer genome_version if not set -------------------------------------
  if (is.null(cfg$genome_version)) {
    if (!is.null(cfg$settings$genome_version)) {
      cfg$genome_version <- cfg$settings$genome_version
    } else if (!is.null(cfg$input$fasta)) {
      fname <- basename(cfg$input$fasta)
      if (grepl("hg19", fname, ignore.case = TRUE)) {
        cfg$genome_version <- "hg19"
      } else if (grepl("hg38", fname, ignore.case = TRUE)) {
        cfg$genome_version <- "hg38"
      } else {
        cfg$genome_version <- "hg19"
        log_msg("Genome version not inferred; defaulting to 'hg19'.")
      }
    } else {
      cfg$genome_version <- "hg19"
    }
  }

  # ---- Save original warning option and set up error handling -------------
  old_warn <- options(warn = 1)
  on.exit(options(warn = old_warn$warn), add = TRUE)

  # ---- Main pipeline execution --------------------------------------------
  tryCatch({
    withCallingHandlers({
      log_msg(paste("Plots enabled:", plots))
      log_msg(paste("BED processing mode:", cfg$bed_process))
      log_msg(paste("Genome version:", cfg$genome_version))
      log_msg(paste("Terminal-exon padding:", paste0(cfg$settings$pad_terminal_exons %||% 0, " bp")))
      log_msg(paste("Remove terminal-only calls:", cfg$settings$remove_terminal_only %||% FALSE))
      log_msg(paste("Penalize terminal-only calls:", cfg$settings$penalize_terminal_only %||% FALSE))
      if (!is.null(cfg$bed_process) && cfg$bed_process != "NO" && !is.null(cfg$off_target_pattern %||% "^HorsROI")) {
        log_msg(paste0("Off-target region handling: pattern = '", cfg$off_target_pattern %||% "^HorsROI",
                       "', handling = '", cfg$off_target_handling %||% "na", "'"))
      }

      # Define output paths
      paths <- list(
        rdata   = file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_coverage.Rdata")),
        metrics = file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_QC_metrics.tsv")),
        cnvs    = file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_CNV_calls.tsv")),
        summary = file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_summary.RData")),
        plots   = file.path(cfg$output$dir, "Plots")
      )

      # Allow command‑line overrides for paths
      if (!is.null(args$rdata))   paths$rdata   <- args$rdata
      if (!is.null(args$metrics)) paths$metrics <- args$metrics
      if (!is.null(args$cnvs))    paths$cnvs    <- args$cnvs
      if (!is.null(args$summary)) paths$summary <- args$summary

      # Create directories for output files
      lapply(paths[!names(paths) %in% "plots"], function(p) {
        if (grepl("\\.[a-zA-Z]+$", p)) {
          dir.create(dirname(p), showWarnings = FALSE, recursive = TRUE)
        } else {
          dir.create(p, showWarnings = FALSE, recursive = TRUE)
        }
      })

      # VCF output path
      if (is.null(vcf_output)) {
        vcf_output <- file.path(cfg$output$dir, paste0("ECHO_", cfg$output$prefix, "_CNVs.vcf"))
      } else if (identical(vcf_output, FALSE)) {
        vcf_output <- NULL
      }

      # Validate BED file existence
      stop_if_not_file(cfg$input$bed, "[ERROR] BED file missing")

      # Reference sequence source
      fasta_source <- cfg$input$fasta_source %||% "file"
      if (identical(fasta_source, "file")) {
        stop_if_not_file(cfg$input$fasta, "[ERROR] FASTA file missing")
      } else {
        log_msg(paste0("Reference sequence source: fasta_source = '", fasta_source,
                        "' (genome_version = ", cfg$genome_version, "); no local FASTA required."))
      }

      # Non‑fatal BED validation
      tryCatch(
        validate_bed_regions(cfg$input$bed, verbose = TRUE),
        error = function(e) log_msg(paste("BED validation skipped:", conditionMessage(e)), "WARNING")
      )

      # ---- BED preprocessing (if not "NO") ---------------------------------
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
          txdb                  = NULL,
          gene_field_index      = gene_field_index,
          off_target_pattern    = cfg$off_target_pattern %||% "^HorsROI",
          off_target_handling   = cfg$off_target_handling %||% "na"
        )
        cfg$input$bed <- processed_bed
        log_msg(paste("Processed BED saved to:", processed_bed))
      } else {
        log_msg("BED preprocessing skipped (bed_process = 'NO').")
      }

      # ---- Step counter ---------------------------------------------------
      steps_total   <- 3L +  
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

      # ---- Step 1: BAM coverage extraction --------------------------------
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
          bed_zero_based       = cfg$bed_zero_based %||% TRUE,
          skip_invalid_intervals = cfg$skip_invalid_intervals %||% TRUE,
          pad_terminal_exons   = cfg$settings$pad_terminal_exons %||% 0
        )
      })

      # ---- PCA plot (optional) --------------------------------------------
      if (cfg$settings$pca_plot %||% TRUE) {
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
      }

      # ---- Step 2: QC metrics ---------------------------------------------
      run_step("Running QC metrics", {
        tryCatch(
          run_qc_metrics(
            rdata_file      = paths$rdata,
            output_file     = paths$metrics,
            min_corr        = cfg$settings$min_corr,
            min_cov         = cfg$settings$min_cov,
            min_total_reads = cfg$settings$min_total_reads,
            max_exon_cv     = cfg$settings$max_exon_cv
          ),
          error = function(e) log_msg(paste("QC metrics step failed:", conditionMessage(e)), "WARNING")
        )
      })

      # ---- Step 3: CNV calling --------------------------------------------
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

      # ---- Plot generation (optional) -------------------------------------
      if (plots) {
        run_step("Generating plots", {
          tryCatch(
            generate_plots(
              rdata_file = paths$summary,
              output_dir = paths$plots,
              modechrom  = cfg$settings$modechrom,
              prefix     = cfg$output$prefix,
              log_file   = log_file,
              gene_gap   = cfg$settings$plot_gene_gap %||% 1
            ),
            error = function(e) log_msg(paste("Plot generation failed:", conditionMessage(e)), "WARNING")
          )
        })
      } else {
        log_msg("Plot generation skipped (plots = FALSE).")
      }

      # ---- HTML report (optional) -----------------------------------------
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
              prefix          = cfg$output$prefix
            ),
            error = function(e) log_msg(paste("Report generation failed:", conditionMessage(e)), "WARNING")
          )
        })
      } else {
        log_msg("Report generation skipped (report = FALSE).")
      }

      # ---- VCF export (optional) ------------------------------------------
      if (!is.null(vcf_output) || vcf_per_sample) {
        run_step("Exporting CNVs to VCF", {
          tryCatch({
          if (file.exists(paths$summary)) {
            cnv_calls_local <- local({
              env <- new.env()
              load(paths$summary, envir = env)
              env$cnv_calls
            })
            if (!is.null(cnv_calls_local) && nrow(cnv_calls_local) > 0) {
              if (!is.null(vcf_output)) {
                export_cnvs_to_vcf(cnv_calls_local, vcf_output, sample_name = NULL)
                log_msg(paste("Combined VCF written to:", basename(vcf_output)))
              }
              if (vcf_per_sample) {
                samples <- unique(cnv_calls_local$Sample)
                for (s in samples) {
                  sample_folder <- file.path(cfg$output$dir, "Plots", sanitize_filename(s))
                  if (!dir.exists(sample_folder))
                    dir.create(sample_folder, recursive = TRUE, showWarnings = FALSE)
                  vcf_file <- file.path(sample_folder,
                                        paste0("ECHO_", cfg$output$prefix, "_", s, ".vcf"))
                  export_cnvs_to_vcf(cnv_calls_local, vcf_file, sample_name = s)
                  log_msg(paste("Per-sample VCF written to:", basename(vcf_file)))
                }
              }
            } else {
              log_msg("No CNV calls to export to VCF.")
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
    }, warning = function(w) {
      if (grepl("sequence levels not in the other", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
      log_msg(conditionMessage(w), "WARNING")
      invokeRestart("muffleWarning")
    })
    invisible(TRUE)
  }, error = function(e) {
    log_msg(paste("FATAL:", e$message), "ERROR")
    stop(e)
  })
}