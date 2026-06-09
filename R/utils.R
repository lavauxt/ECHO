#' Stop if value is NULL or empty
#'
#' @param val Object to validate.
#' @param msg Error message to display if validation fails.
#' @export
stop_if_missing <- function(val, msg) {
    if (is.null(val) || length(val) == 0) stop(msg)
}

#' Stop if file does not exist
#'
#' @param path Path to the file.
#' @param msg Error message to display if file is missing.
#' @export
stop_if_not_file <- function(path, msg) {
    if (is.null(path) || !file.exists(path)) stop(msg)
}

#' Load an RData file into an isolated environment
#'
#' @param path Path to the \code{.RData} file.
#' @param required Optional character vector of required object names.
#' @return A named list of objects from the file.
#' @noRd
load_rdata <- function(path, required = NULL) {
    stop_if_not_file(path, paste0("[ERROR] RData file not found: ", path))
    env <- new.env()
    load(path, envir = env)
    if (!is.null(required)) {
        missing <- setdiff(required, ls(env))
        if (length(missing) > 0) {
            stop("[ERROR] Missing variables in RData: ", paste(missing, collapse = ", "),
                 "\nAvailable objects: ", paste(ls(env), collapse = ", "))
        }
    }
    as.list(env)
}

#' Check that a BAM file has an index on disk
#'
#' @param bam Path to a BAM file.
#' @return Logical.
#' @noRd
bam_has_index <- function(bam) {
    file.exists(paste0(bam, ".bai")) || file.exists(paste0(bam, ".bam.bai"))
}

#' Normalise a chromosome name to match a reference naming style (vectorised)
#'
#' @param chr_vec Character vector of chromosome names.
#' @param ref_chromosomes Character vector of reference chromosome names (e.g. BED).
#' @return Character vector with prefix aligned to reference style.
#' @noRd
normalize_chromosome_vec <- function(chr_vec, ref_chromosomes) {
    has_chr_ref <- any(grepl("^chr", ref_chromosomes))
    if (has_chr_ref) {
        ifelse(grepl("^chr", chr_vec), chr_vec, paste0("chr", chr_vec))
    } else {
        ifelse(grepl("^chr", chr_vec), sub("^chr", "", chr_vec), chr_vec)
    }
}

#' Filter a data frame by chromosome column
#'
#' @param df Input data frame.
#' @param include Optional vector of chromosomes to include.
#' @param exclude Optional vector of chromosomes to exclude.
#'
#' @return Filtered data frame (empty if input is empty).
#' @export
filter_chromosomes <- function(df, include = NULL, exclude = NULL) {
    if (is.null(df) || nrow(df) == 0) return(df)
    chrom_col <- intersect(c("chromosome", "Chromosome"), colnames(df))[1]
    if (is.na(chrom_col)) {
        stop("Data frame must contain a 'chromosome' or 'Chromosome' column.")
    }
    norm <- function(x) {
        unique(c(x, sub("^chr", "", x), paste0("chr", sub("^chr", "", x))))
    }
    if (!is.null(include)) {
        df <- df[df[[chrom_col]] %in% norm(include), ]
    }
    if (!is.null(exclude)) {
        df <- df[!df[[chrom_col]] %in% norm(exclude), ]
    }
    df
}

#' Null coalescing operator
#'
#' Returns the first argument if not NULL, otherwise the second.
#'
#' @param a Primary value.
#' @param b Fallback value.
#' @export
`%||%` <- function(a, b) {
    if (!is.null(a)) a else b
}

#' @import data.table
NULL

#' Sanitize a string for use as a filename
#'
#' @param name Character string to sanitize.
#' @export
sanitize_filename <- function(name) {
    gsub("[^[:alnum:]]", "_", name)
}

#' Convert global exon indices to within‑gene indices (using start‑order)
#'
#' For a set of CNV calls (each with global start.p and end.p), this function
#' replaces those columns with per‑gene exon numbers. The original global indices
#' are preserved as `global_start` and `global_end`.
#'
#' @param cnv_calls Data frame of CNV calls (from ExomeDepth).
#' @param bed_file BED annotation data frame (must have columns `gene`, `start`).
#'
#' @return The same data frame with updated `start.p` / `end.p` (per‑gene) and
#'   added `global_start` / `global_end` columns.
#' @export
add_within_gene_indices <- function(cnv_calls, bed_file) {
    if (nrow(cnv_calls) == 0) return(cnv_calls)
    cnv_calls$global_start <- as.numeric(cnv_calls$start.p)
    cnv_calls$global_end   <- as.numeric(cnv_calls$end.p)
    # Order exons by start coordinate within each gene
    exon_in_gene <- ave(seq_len(nrow(bed_file)), bed_file$gene,
                        FUN = function(idx) order(bed_file$start[idx]))
    cnv_calls$start.p <- exon_in_gene[cnv_calls$global_start]
    cnv_calls$end.p   <- exon_in_gene[cnv_calls$global_end]
    cnv_calls
}

#' Compute exon index for plotting (consistent with add_within_gene_indices)
#'
#' @param bed_file BED annotation data frame.
#' @return Integer vector of per‑gene exon numbers, one per row.
#' @export
compute_exon_index <- function(bed_file) {
    ave(seq_len(nrow(bed_file)), bed_file$gene,
        FUN = function(idx) order(bed_file$start[idx]))
}

#' Load and validate YAML configuration
#'
#' @param yaml_path Character string. Path to a YAML file, or a filename resolved
#'   via \code{system.file()} when not found locally.
#'
#' @return Named list with \code{input}, \code{output}, and optional \code{settings}
#'   sections from the YAML file.
#' @export
load_config <- function(yaml_path) {
    if (!file.exists(yaml_path)) {
        found_path <- system.file(yaml_path, package = "ECHO")
        if (found_path == "") {
            stop("[ERROR] Config file not found at: ", yaml_path)
        }
        yaml_path <- found_path
    }
    cfg <- yaml::read_yaml(yaml_path)
    required <- c("input", "output")
    if (!all(required %in% names(cfg))) {
        stop("[ERROR] YAML must contain 'input' and 'output' sections")
    }
    return(cfg)
}

#' Initialize a log file
#'
#' @param log_file Path to log file.
#' @param session_info Logical; include session info.
#' @return Invisibly returns the log file path.
#' @noRd
init_log <- function(log_file, session_info = TRUE) {
    log_dir <- dirname(log_file)
    if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    log_con <- file(log_file, open = "wt")
    writeLines(paste0("# ECHO Pipeline Log - ", Sys.time()), log_con)
    writeLines(paste0("# Log file: ", log_file), log_con)
    writeLines("", log_con)
    if (session_info) {
        writeLines("## Session Info", log_con)
        si <- utils::sessionInfo()
        writeLines(paste0("R version: ", si$R.version$version.string), log_con)
        writeLines("Platform:", log_con)
        writeLines(paste0("  ", si$platform), log_con)
        writeLines("Packages:", log_con)
        for (pkg in names(si$otherPkgs)) {
            writeLines(paste0("  ", pkg, ": ", si$otherPkgs[[pkg]]$Version), log_con)
        }
        writeLines("", log_con)
    }
    close(log_con)
    invisible(log_file)
}

#' Write a message to both console and log file
#'
#' @param msg Character string.
#' @param log_file Path to log file.
#' @param type One of "INFO", "WARNING", "ERROR".
#' @noRd
log_message <- function(msg, log_file, type = "INFO") {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    formatted <- paste0("[", type, "] ", timestamp, " ", msg)
    message(formatted)
    cat(formatted, "\n", file = log_file, append = TRUE)
}

#' Capture warnings and redirect to log
#'
#' @param expr Expression to evaluate.
#' @param log_file Path to log file.
#' @noRd
with_warning_log <- function(expr, log_file) {
    withCallingHandlers(
        expr,
        warning = function(w) {
            log_message(conditionMessage(w), log_file, "WARNING")
            invokeRestart("muffleWarning")
        }
    )
}