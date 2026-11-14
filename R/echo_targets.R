# Opt-in {targets} layer for ECHO.
#
# BAM counting (the slow step) is skipped by echo(reuse_coverage = TRUE) whenever
# the BAMs, BED, FASTA and counting settings are unchanged, so with this layer
# a change of the calling / QC / report options re-runs the calling only, and a
# change of a BAM or of the BED re-counts. All the logic sits in plain
# functions (echo_targets_*); echo_targets() only wires them into a DAG.

#' ECHO as a {targets} pipeline (opt-in)
#'
#' Returns a list of targets for a \code{_targets.R} file:
#' \preformatted{
#' library(targets)
#' list(ECHO::echo_targets("config.yaml"))
#' }
#' Targets: \code{echo_config} (the YAML, tracked), \code{echo_bams} (BAM files,
#' tracked by size / mtime), \code{echo_run} (every file \code{echo()} wrote) and
#' \code{echo_calls} (the call table as a data frame, when one was written).
#'
#' Requires the \pkg{targets} package (>= 1.3.0), a suggested dependency.
#'
#' @param config_path Path to the YAML configuration (see \code{\link{echo}}).
#' @param ... Further arguments of \code{\link{echo}} (e.g. \code{report}, \code{plots}).
#' @return A list of \code{targets::tar_target()} objects.
#' @export
echo_targets <- function(config_path = "config.yaml", ...) {
  cnv_require_targets()
  dots <- list(...)
  list(
    targets::tar_target(echo_config, config_path, format = "file"),
    targets::tar_target(echo_bams, echo_targets_bams(echo_config), format = "file_fast"),
    targets::tar_target(echo_run, echo_targets_run(echo_config, echo_bams, !!!dots), format = "file"),
    targets::tar_target(echo_calls, echo_targets_calls(echo_run))
  )
}

#' The BAM files of the config (test and external reference BAMs)
#' @noRd
echo_targets_bams <- function(config_path) {
  cfg <- load_config(config_path)
  echo_input_bams(cfg$input$bamfiles, cfg$input$bamdir, cfg$input$rbams)
}

#' Run echo() reusing the coverage table; returns the files it wrote
#' @noRd
echo_targets_run <- function(config_path, bams = NULL, ...) {
  t0 <- Sys.time()
  args <- list(...)
  if (is.null(args$reuse_coverage)) args$reuse_coverage <- TRUE
  do.call(echo, c(list(config_path = config_path), args))
  out_dir <- load_config(config_path)$output$dir
  files <- cnv_files_since(out_dir, t0)
  cnv <- files[grepl("_CNV_calls\\.tsv$", files)]
  unique(c(cnv, files))
}

#' The call table written by an echo_targets_run() (NULL when none)
#' @noRd
echo_targets_calls <- function(files) {
  f <- files[grepl("_CNV_calls\\.tsv$", files)]
  if (!length(f)) return(NULL)
  utils::read.delim(f[1L], check.names = FALSE, stringsAsFactors = FALSE)
}
