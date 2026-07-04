#' Generate Interactive HTML Report
#'
#' @param summary_rdata Path to ECHO_summary.RData
#' @param qc_metrics_file Path to QC_Metrics.tsv
#' @param output_dir Output directory for the report (will be created if missing)
#' @param settings List of pipeline settings
#' @param config Original config list (optional)
#' @param global Logical. If TRUE, generate a single global report for all samples.
#'   If FALSE (default), generate one report per sample (legacy mode).
#' @param sample_name Optional single sample name (only used when global = FALSE).
#' @param sample_table Optional path to sample table (CSV/TSV) to display in report.
#' @param ref_bams Optional path to external reference BAM list (TSV with 'bam' column).
#' @param log_file Path to log file (optional).
#' @param pdf_output Logical. If TRUE, also generate a PDF version using pagedown.
#' @param prefix Character string. Prefix to use for output file name (e.g., "GOM").
#'   When provided, the global report will be named `ECHO_<prefix>_report.html`.
#'
#' @return Invisibly returns the output file path(s).
#' @export
generate_report <- function(summary_rdata, qc_metrics_file, output_dir, 
                            settings, config = NULL, global = TRUE, sample_name = NULL,
                            sample_table = NULL, ref_bams = NULL, log_file = NULL, 
                            pdf_output = FALSE, prefix = NULL) {
  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("Package 'rmarkdown' is required for report generation.")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) {
      stop("Failed to create output directory: ", output_dir)
    }
  }

  # Load data into local environment to avoid conflicts
  local_env <- new.env()
  load(summary_rdata, envir = local_env)
  cnv_calls <- local_env$cnv_calls
  counts    <- local_env$counts
  bed_file  <- local_env$bed_file
  models    <- local_env$models
  refs      <- local_env$refs
  gender_map <- local_env$gender_map
  sample_names <- local_env$sample_names   

  # Use fread to handle missing fields gracefully
  qc_metrics <- data.table::fread(qc_metrics_file, header = TRUE, fill = TRUE, data.table = FALSE)

  # Load sample table if provided and not already in RData
  sample_df <- NULL
  if (!is.null(sample_table) && file.exists(sample_table)) {
    sample_df <- utils::read.table(sample_table, header = TRUE, sep = "\t",
                                   stringsAsFactors = FALSE)
    if (!all(c("sample_name", "gender") %in% colnames(sample_df))) {
      warning("sample_table missing required columns 'sample_name' or 'gender'")
      sample_df <- NULL
    }
  } else if (!is.null(gender_map)) {
    # Convert gender_map to data frame
    sample_df <- data.frame(sample_name = names(gender_map),
                            gender = gender_map,
                            stringsAsFactors = FALSE)
  }

  # Load external reference BAM list if provided
  ref_bams_df <- NULL
  if (!is.null(ref_bams) && file.exists(ref_bams)) {
    ref_bams_df <- utils::read.table(ref_bams, header = TRUE, sep = "\t",
                                     stringsAsFactors = FALSE)
    if (!"bam" %in% colnames(ref_bams_df)) {
      warning("ref_bams file missing 'bam' column")
      ref_bams_df <- NULL
    }
  }

  # Captured immediately, before anything else in this function could change
  # it, so relative paths built by the caller (log_file, in particular)
  # keep meaning what the caller meant by them once rmarkdown::render()
  # changes the working directory during knitting — see knit_root_dir below.
  caller_wd <- getwd()

  if (global) {
    # BUG FIX: the template ships at inst/rmd/, so system.file() resolves
    # it under "rmd/..." once installed, not "rmarkdown/...". The old path
    # never matched anything post-install (system.file() silently returns
    # "" for a nonexistent subdirectory), so this always failed.
    template <- system.file("rmd/ECHO_global_report.Rmd", package = "ECHO")
    if (template == "") stop("Global report template not found. Reinstall package.")
    
    # Determine output file name
    if (!is.null(prefix) && nzchar(prefix)) {
      out_file <- file.path(output_dir, paste0("ECHO_", prefix, "_report.html"))
    } else {
      out_file <- file.path(output_dir, "ECHO_global_report.html")
    }
    # BUG FIX: rmarkdown::render() does not resolve a relative output_file
    # against the caller's working directory — it resolves it against the
    # directory containing the *input* Rmd file (template, above, which
    # system.file() resolves to somewhere inside the installed package,
    # essentially never the same place as output_dir). A perfectly valid,
    # already-created output_dir would make render() look for it relative
    # to the package's own installation folder and fail with "The
    # directory '...' does not exist" — confirmed by direct reproduction.
    # output_dir is already guaranteed to exist by this point, so
    # normalizePath() on it is safe; only the (not-yet-existing) filename
    # is appended afterward.
    out_file <- file.path(normalizePath(output_dir, mustWork = TRUE), basename(out_file))
    
    rmarkdown::render(template,
                      params = list(
                        cnv_calls = cnv_calls,
                        qc_metrics = qc_metrics,
                        counts = counts,
                        bed_file = bed_file,
                        models = models,
                        refs = refs,
                        settings = settings,
                        config = config,
                        sample_table = sample_df,
                        ref_bams = ref_bams_df,
                        log_file = log_file,
                        pdf_output = pdf_output
                      ),
                      output_file = out_file,
                      # BUG FIX: render() also evaluates the Rmd's own code
                      # chunks with the working directory set to the Rmd's
                      # directory by default, not the caller's — a second,
                      # quieter instance of the same problem. The "Pipeline
                      # Log" section reads params$log_file (built by the
                      # caller relative to *its own* working directory)
                      # behind a file.exists() guard, so this wasn't fatal,
                      # just a silently-empty log section with no error.
                      # Pinning knit_root_dir keeps every relative path
                      # inside the Rmd meaning what its caller intended.
                      knit_root_dir = caller_wd,
                      quiet = FALSE)
    message("[INFO] Global report written: ", out_file)
    return(invisible(out_file))
  } else {
    if (is.null(sample_name)) {
      samples <- unique(cnv_calls$Sample)
      if (length(samples) == 0) {
        message("[INFO] No CNV calls – skipping report.")
        return(invisible(NULL))
      }
    } else {
      samples <- sample_name
    }
    template <- system.file("rmd/ECHO_report.Rmd", package = "ECHO")
    if (template == "") stop("Report template not found. Reinstall package.")
    out_files <- c()
    for (s in samples) {
      sample_dir <- file.path(output_dir, sanitize_filename(s))
      dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
      out_file <- file.path(sample_dir, paste0("ECHO_report_", sanitize_filename(s), ".html"))
      out_file <- file.path(normalizePath(sample_dir, mustWork = TRUE), basename(out_file))
      rmarkdown::render(template,
                        params = list(
                          sample = s,
                          cnv_calls = cnv_calls,
                          qc_metrics = qc_metrics,
                          counts = counts,
                          bed_file = bed_file,
                          models = models,
                          refs = refs,
                          settings = settings,
                          config = config,
                          sample_table = sample_df,
                          ref_bams = ref_bams_df,
                          log_file = log_file,
                          pdf_output = pdf_output,
                          sample_names = sample_names,
                          gender_map = gender_map     
                        ),
                        output_file = out_file,
                        knit_root_dir = caller_wd,
                        quiet = FALSE)
      out_files <- c(out_files, out_file)
      message("[INFO] Report written: ", out_file)
    }
    invisible(out_files)
  }
}