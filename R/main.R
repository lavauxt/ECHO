#' Run the Full ECHO Pipeline
#'
#' Runs coverage extraction, QC metrics, CNV calling, plot generation, and optional VCF export.
#'
#' @param config_path Optional character string. Path to a \code{config.yaml} file.
#' @param ... Optional named overrides (see details).
#' @param vcf_output Character string. If \code{NULL}, no VCF is written. Otherwise, path to output VCF file.
#'   Default is \code{NULL} (legacy), but when using a config file it will be set to a default path inside the output directory.
#' @param save_ed_objects Logical. Save full ExomeDepth objects? Default \code{FALSE}.
#' @param report Logical. Generate interactive HTML reports? Default \code{TRUE}.
#' @param sample_name_delim Character string. Delimiter(s) to split the filename.
#'   Can be a single character (e.g., `"."`, `"_"`) or a regex like `"[._]"`. Default `"\\."` (dot only).
#' @param sample_name_keep Character string. Specifies which parts to keep after splitting.
#'   Examples: `"1"` (first part), `"1-2"` (first two parts), `"2"` (second part), `"2-3"` (parts 2 and 3).
#'   Default `"1"` (keep the first part before the first delimiter).
#' @param sample_name_collapse Character string. Separator for rejoining parts. Default NULL (uses first character of delim if single char else ".").
#' @param custom_sample_names Optional character vector. Custom sample names (must match number of BAMs).
#' @param log_file Character string. Path to log file. If \code{NULL}, a default log file is created inside the output directory.
#' @param gene_field_index Integer. Which field to keep after splitting by exon_sep in BED processing.
#'   Default 1 (first field). Use 3 for BEDs like "NM_006015_ARID1A_ex5...".
#' @param sample_table Optional character string. Path to a CSV/TSV file with columns
#'   \code{sample_name} and \code{gender} (M/F or male/female). Required for sex chromosome modes.
#' @param ref_bams Optional character string. Path to a TSV file with a \code{bam} column
#'   listing external reference BAMs (used for reporting only).
#' @param panel_files Optional character vector. In STANDARD mode, a vector of BED file paths
#'   (one per line) or a single file containing paths to panel BEDs used to restrict targets.
#'
#' @return Invisibly returns \code{TRUE} on success.
#' @export
echo <- function(config_path = NULL, vcf_output = NULL, save_ed_objects = FALSE,
                 report = TRUE,
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

    # -------------------------------------------------------------------------
    # 1. Load configuration
    # -------------------------------------------------------------------------
    if (!is.null(config_path)) {
        cfg <- load_config(config_path)
        
        # Merge BED preprocessing parameters from YAML if present
        if (!is.null(cfg$bed_preprocess)) {
            for (name in names(cfg$bed_preprocess)) {
                cfg[[name]] <- cfg$bed_preprocess[[name]]
            }
        }
        
    } else if (length(args) > 0) {
        cfg <- list(
            input = list(
                bamdir   = args$bamdir %||% "./data",
                bamfiles = args$bamfiles,
                bed      = args$bed,
                fasta    = args$fasta,
                rbams    = args$rbams
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
                score_low_confidence_genes = args$score_low_confidence_genes %||% c("PMS2", "SMN1", "CYP2D6", "HBA1", "HBA2", "STRC", "CYP21A2", "GBA1", "CFTR")
            ),
            # Added BED preprocessing parameters (optional)
            bed_process       = args$bed_process %||% "NO",
            refseqgene        = args$refseqgene,
            transcripts_file  = args$transcripts_file,
            unknown_gene      = args$unknown_gene %||% FALSE,
            gene_list_restrict = args$gene_list_restrict,
            chr_list_restrict = args$chr_list_restrict,
            exon_sep          = args$exon_sep,
            customexon        = args$customexon %||% FALSE,
            list_genes        = args$list_genes,
            genes_file        = args$genes_file,
            genome_version    = args$genome_version %||% "hg19",
            bed_zero_based    = args$bed_zero_based %||% TRUE,
            skip_invalid_intervals = args$skip_invalid_intervals %||% TRUE,
            sample_name_collapse = args$sample_name_collapse %||% sample_name_collapse,
            panel_files       = args$panel_files %||% panel_files
        )
    } else {
        stop("[ERROR] Either config_path or pipeline parameters must be provided.")
    }

    # Allow gene_field_index from config
    if (!is.null(cfg$gene_field_index)) {
        gene_field_index <- cfg$gene_field_index
    } else if (!is.null(cfg$bed_preprocess$gene_field_index)) {
        gene_field_index <- cfg$bed_preprocess$gene_field_index
    }

    # Convert output directory to absolute path and create it
    cfg$output$dir <- normalizePath(file.path(getwd(), cfg$output$dir), mustWork = FALSE)
    if (!dir.exists(cfg$output$dir)) {
        dir.create(cfg$output$dir, recursive = TRUE, showWarnings = FALSE)
    }

    # -------------------------------------------------------------------------
    # 2. Setup logging
    # -------------------------------------------------------------------------
    if (is.null(log_file)) {
        log_file <- file.path(cfg$output$dir, "pipeline.log")
    }
    log_dir <- dirname(log_file)
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

    log_con <- file(log_file, open = "wt")
    writeLines(paste0("# ECHO Pipeline Log - ", Sys.time()), log_con)
    writeLines(paste0("# Output directory: ", basename(cfg$output$dir)), log_con)
    writeLines("#", log_con)
    writeLines("## Session Info", log_con)
    si <- utils::sessionInfo()
    writeLines(paste0("R version: ", si$R.version$version.string), log_con)
    writeLines(paste0("Platform: ", si$platform), log_con)
    writeLines("Packages:", log_con)
    for (pkg in names(si$otherPkgs)) {
        writeLines(paste0("  ", pkg, ": ", si$otherPkgs[[pkg]]$Version), log_con)
    }
    writeLines("", log_con)
    close(log_con)

    log_msg <- function(msg, type = "INFO") {
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        formatted <- paste0("[", type, "] ", timestamp, " ", msg)
        message(formatted)
        cat(formatted, "\n", file = log_file, append = TRUE)
    }

    old_warn <- options(warn = 1)
    on.exit(options(warn = old_warn$warn), add = TRUE)

    tryCatch({
        withCallingHandlers({

            log_msg("ECHO pipeline started")
            log_msg(paste("Configuration:", ifelse(is.null(config_path), "command-line parameters", config_path)))
            log_msg(paste("Output directory:", basename(cfg$output$dir)))
            log_msg(paste("Log file:", basename(log_file)))

            # -----------------------------------------------------------------
            # 3. Define output file paths
            # -----------------------------------------------------------------
            paths <- list(
                rdata   = file.path(cfg$output$dir, "ECHO_coverage.Rdata"),
                metrics = file.path(cfg$output$dir, "QC_metrics.tsv"),
                cnvs    = file.path(cfg$output$dir, "CNV_calls.tsv"),
                summary = file.path(cfg$output$dir, "ECHO_summary.RData"),
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
                vcf_output <- file.path(cfg$output$dir, "CNVs.vcf")
            }

            stop_if_not_file(cfg$input$bed, "[ERROR] BED file missing")
            stop_if_not_file(cfg$input$fasta, "[ERROR] FASTA file missing")

            # -----------------------------------------------------------------
            # 4. Optional BED preprocessing
            # -----------------------------------------------------------------
            if (!is.null(cfg$bed_process) && cfg$bed_process != "NO") {
                log_msg("Preprocessing BED file using mode: ", cfg$bed_process)
                processed_bed <- file.path(cfg$output$dir, paste0(cfg$output$prefix, "_targets.bed"))
                process_bed_file(
                    input_bed = cfg$input$bed,
                    output_bed = processed_bed,
                    bed_process = cfg$bed_process,
                    bed_zero_based = cfg$bed_zero_based %||% TRUE,
                    refseqgene = cfg$refseqgene %||% NULL,
                    transcripts_file = cfg$transcripts_file %||% NULL,
                    unknown_gene = cfg$unknown_gene %||% FALSE,
                    gene_list_restrict = cfg$gene_list_restrict %||% NULL,
                    exon_sep = cfg$exon_sep %||% NULL,
                    customexon = cfg$customexon %||% FALSE,
                    list_genes = cfg$list_genes %||% NULL,
                    genes_file = cfg$genes_file %||% NULL,
                    panel_files = cfg$panel_files %||% NULL,
                    genome_version = cfg$genome_version %||% "hg19",
                    txdb = NULL,
                    gene_field_index = gene_field_index
                )
                cfg$input$bed <- processed_bed
                log_msg("Processed BED saved to: ", processed_bed)
            } else {
                log_msg("BED preprocessing skipped (bed_process = 'NO').")
            }

            # -----------------------------------------------------------------
            # 5. Run pipeline steps with warning capture
            # -----------------------------------------------------------------
            log_msg("Step 1/6: Extracting BAM coverage...")
            run_bam_coverage(
                bamfiles = cfg$input$bamfiles,
                bamdir   = cfg$input$bamdir,
                bed      = cfg$input$bed,
                fasta    = cfg$input$fasta,
                rbams    = cfg$input$rbams,
                data_out = paths$rdata,
                verbose = TRUE,
                sample_name_delim = sample_name_delim,
                sample_name_keep = sample_name_keep,
                sample_name_collapse = cfg$sample_name_collapse %||% sample_name_collapse,
                custom_sample_names = custom_sample_names,
                bed_zero_based = cfg$bed_zero_based %||% TRUE,
                skip_invalid_intervals = cfg$skip_invalid_intervals %||% TRUE
            )
            log_msg("Step 1/6 completed.")

            log_msg("Step 2/6: Running QC metrics...")
            run_qc_metrics(
                rdata_file      = paths$rdata,
                output_file     = paths$metrics,
                min_corr        = cfg$settings$min_corr,
                min_cov         = cfg$settings$min_cov,
                min_total_reads = cfg$settings$min_total_reads,
                max_exon_cv     = cfg$settings$max_exon_cv
            )
            log_msg("Step 2/6 completed.")

            log_msg("Step 3/6: Calling CNVs...")
            cnv_args <- cfg$settings[grepl("^score_", names(cfg$settings))]
            cnv_args <- c(list(
                rdata_file            = paths$rdata,
                output_file           = paths$cnvs,
                out_rdata             = paths$summary,
                transition.probability = cfg$settings$transition_probability,
                expected.CNV.length   = cfg$settings$expected_CNV_length,
                n.bins.reduced        = cfg$settings$n_bins_reduced,
                phi.bins              = cfg$settings$phi_bins,
                formula               = cfg$settings$formula,
                data                  = NULL,
                save_ed_objects       = save_ed_objects,
                modechrom             = cfg$settings$modechrom,
                sample_table          = sample_table
            ), cnv_args)

            do.call(run_cnv_calling, cnv_args)
            log_msg("Step 3/6 completed.")

            log_msg("Step 4/6: Generating plots...")
            generate_plots(
                rdata_file = paths$summary,
                output_dir = paths$plots,
                modechrom  = cfg$settings$modechrom,
                prefix     = cfg$output$prefix,
                log_file   = log_file
            )
            log_msg("Step 4/6 completed.")

            if (report) {
                log_msg("Step 5/6: Generating HTML report...")
                generate_report(summary_rdata = paths$summary,
                                qc_metrics_file = paths$metrics,
                                output_dir = cfg$output$dir,
                                settings = cfg$settings,
                                config = cfg,
                                sample_table = sample_table,
                                ref_bams = ref_bams %||% cfg$input$rbams,
                                log_file = log_file,
                                pdf_output = cfg$settings$pdf_output %||% FALSE)
                log_msg("Step 5/6 completed.")
            } else {
                log_msg("Step 5/6: Report generation skipped (report = FALSE).")
            }

            if (!is.null(vcf_output)) {
                log_msg("Step 6/6: Exporting CNVs to VCF...")
                if (file.exists(paths$summary)) {
                    cnv_calls_local <- local({
                        env <- new.env()
                        load(paths$summary, envir = env)
                        env$cnv_calls
                    })
                    if (exists("cnv_calls_local") && !is.null(cnv_calls_local) && nrow(cnv_calls_local) > 0) {
                        export_cnvs_to_vcf(cnv_calls_local, vcf_output, sample_name = NULL)
                        log_msg(paste("VCF written to:", basename(vcf_output)))
                    } else {
                        log_msg("No CNV calls to export to VCF.")
                    }
                } else {
                    log_msg("Summary RData not found – cannot export VCF.", "WARNING")
                }
                log_msg("Step 6/6 completed.")
            } else {
                log_msg("Step 6/6: VCF export skipped (vcf_output = NULL).")
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