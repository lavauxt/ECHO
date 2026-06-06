#' Stop if value is NULL or empty
#'
#' @param val Object to validate.
#' @param msg Error message to display if validation fails.
#'
#' @return Invisibly returns \code{NULL}, or stops with an error.
#' @export
stop_if_missing <- function(val, msg) {
    if (is.null(val) || length(val) == 0) stop(msg)
}

#' Stop if file does not exist
#'
#' @param path Path to the file.
#' @param msg Error message to display if file is missing.
#'
#' @return Invisibly returns \code{NULL}, or stops with an error.
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
    env <- new.env(parent = emptyenv())
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

#' Append metric failure entries to an existing metrics list
#'
#' @param m Existing metrics list.
#' @param samples Vector of sample names.
#' @param exons Vector of exon identifiers.
#' @param type Metric type label.
#' @param details Vector of detail messages.
#' @param genes Vector of gene names.
#'
#' @return A list containing updated metric fields.
#' @export
append_metric <- function(m, samples, exons, type, details, genes) {
    list(
        Sample  = c(m$Sample,  samples),
        Exon    = c(m$Exon,    exons),
        Type    = c(m$Type,    rep(type, length(samples))),
        Details = c(m$Details, details),
        Gene    = c(m$Gene,    genes)
    )
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
    # If input is empty, return empty data frame
    if (is.null(df) || nrow(df) == 0) {
        return(df)
    }
    
    chrom_col <- intersect(c("chromosome", "Chromosome"), colnames(df))[1]

    if (is.na(chrom_col)) {
        stop("Data frame must contain a 'chromosome' or 'Chromosome' column.")
    }

    norm <- function(x) {
        unique(c(
            x,
            sub("^chr", "", x),
            paste0("chr", sub("^chr", "", x))
        ))
    }

    if (!is.null(include)) {
        df <- df[df[[chrom_col]] %in% norm(include), ]
    }

    if (!is.null(exclude)) {
        df <- df[!df[[chrom_col]] %in% norm(exclude), ]
    }

    df
}

#' Convert wide coverage matrix to long-format data frame
#'
#' Base R replacement for reshape::melt.
#'
#' @param mat Coverage matrix.
#' @param exon_range Vector of exon indices.
#'
#' @return Long-format data frame with exon, sample, and coverage columns.
#' @export
melt_coverage <- function(mat, exon_range) {
    cols <- colnames(mat)

    data.frame(
        exon_idx = rep(exon_range, length(cols)),
        sample   = rep(cols, each = length(exon_range)),
        coverage = unlist(mat, use.names = FALSE),
        stringsAsFactors = FALSE
    )
}

#' Centralized environment validation
#'
#' @param config Named list with at least \code{bed} (path) and \code{outdir}
#'   (output directory).
#'
#' @return A named list of output paths (\code{rdata}, \code{metrics}, \code{plots}).
#' @export
check_env <- function(config) {
    stop_if_not_file(config$bed, "[ERROR] BED file missing")

    dir.create(
        config$outdir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    list(
        rdata   = file.path(config$outdir, "ECHO_coverage.Rdata"),
        metrics = file.path(config$outdir, "QC_metrics.tsv"),
        plots   = file.path(config$outdir, "Plots")
    )
}

#' Validate BED file schema
#'
#' @param bed_path Character string. Path to the BED file.
#'
#' @return Invisibly returns \code{NULL} if validation succeeds.
#' @export
validate_bed <- function(bed_path) {
    bed <- read.table(
        bed_path,
        header = FALSE,
        stringsAsFactors = FALSE
    )

    if (ncol(bed) < 4) {
        stop("BED file must have at least 4 columns (chrom, start, end, gene)")
    }

    if (!is.numeric(bed[[2]]) || !is.numeric(bed[[3]])) {
        stop("Coordinates in BED must be numeric")
    }

    message("[INFO] BED file validated successfully")
    invisible(NULL)
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

        found_path <- system.file(
            yaml_path,
            package = "ECHO"
        )

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

#' Initialize pipeline output directories (returns paths, does NOT create dirs – done by caller)
#'
#' Creates the output directory tree and returns standard file paths.
#'
#' @param output_dir Character string. Base output directory.
#'
#' @return Named list with elements \code{rdata}, \code{metrics}, \code{cnvs},
#'   \code{summary}, and \code{plots}.
#' @export
init_pipeline_dirs <- function(output_dir) {

    paths <- list(
        rdata   = file.path(output_dir, "ECHO_coverage.Rdata"),
        metrics = file.path(output_dir, "QC_metrics.tsv"),
        cnvs    = file.path(output_dir, "CNV_calls.tsv"),
        summary = file.path(output_dir, "ECHO_summary.RData"),
        plots   = file.path(output_dir, "Plots")
    )

    return(paths)
}

#' Null coalescing operator
#'
#' Returns the first argument if not NULL, otherwise the second.
#'
#' @param a Primary value.
#' @param b Fallback value.
#'
#' @return Either `a` or `b`.
#' @export
`%||%` <- function(a, b) {
    if (!is.null(a)) a else b
}

#' @import data.table
NULL

#' Sort a BED file
#'
#' Sorts a BED file by chromosome, start, and end coordinates.
#'
#' @param infile Input BED file path.
#' @param outfile Output BED file path.
#' @param genomic_order Logical; whether to enforce genomic chromosome ordering.
#'
#' @return Invisibly returns \code{NULL}. Writes the sorted BED file to disk.
#' @export
sort_bed <- function(infile, outfile, genomic_order = TRUE) {

    bed <- data.table::fread(infile, header = FALSE)

    if (genomic_order) {

        chr_levels <- c(
            paste0("chr", 1:22),
            "chrX",
            "chrY",
            "chrM"
        )

        bed[, V1 := factor(V1, levels = chr_levels)]
    }

    data.table::setorder(bed, V1, V2, V3)

    data.table::fwrite(
        bed,
        outfile,
        sep = "\t",
        col.names = FALSE
    )
    invisible(NULL)
}

#' Sanitize a string for use as a filename
#'
#' @param name Character string to sanitize.
#'
#' @return Sanitized filename string.
#' @export
sanitize_filename <- function(name) {
    gsub("[^[:alnum:]]", "_", name)
}

#' Convert global exon indices to within‑gene indices
#'
#' For a set of CNV calls (each with global start.p and end.p), this function
#' replaces those columns with per‑gene exon numbers. The original global indices
#' are preserved as `global_start` and `global_end`. The input `bed_file` is not modified.
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
    
    # Create a local copy of the exon-in-gene mapping without modifying input bed_file
    exon_in_gene <- ave(seq_len(nrow(bed_file)), bed_file$gene,
                        FUN = function(idx) order(bed_file$start[idx]))
    
    cnv_calls$start.p <- exon_in_gene[cnv_calls$global_start]
    cnv_calls$end.p   <- exon_in_gene[cnv_calls$global_end]
    
    cnv_calls
}


#' Initialize a log file
#'
#' @param log_file Path to log file.
#' @param session_info Logical; include session info.
#'
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