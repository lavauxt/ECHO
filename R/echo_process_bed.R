#' Process BED File for CNV Analysis
#' @export
process_bed_file  <- function(input_bed, output_bed, bed_process = "STANDARD",
bed_zero_based = TRUE,
refseqgene = NULL, transcripts_file = NULL,
unknown_gene = FALSE, gene_list_restrict = NULL,
chr_list_restrict = NULL,
exon_sep = " ", gene_name_collapse = " ", customexon = FALSE,
auto_exon_number = TRUE, region_numbering_mode = "bed_text",
gene_name_keep = NULL,
list_genes = NULL, genes_file = NULL,
panel_files = NULL,
genome_version = NULL, txdb = NULL,
gene_field_index = NULL,
off_target_pattern = "^HorsROI",
off_target_handling = c("na", "remove", "merge")) {
off_target_handling  <- match.arg(off_target_handling)
parse_keep_indices  <- function(keep_str, max_len) {
if (is.null(keep_str) || keep_str == "") return(NULL)
if (is.numeric(keep_str)) return(keep_str)
keep_str  <- as.character(keep_str)
if (grepl("-", keep_str)) {
parts  <- as.numeric(strsplit(keep_str, "-")[[1]])
if (length(parts) == 2 && !any(is.na(parts))) {
return(seq(parts[1], min(parts[2], max_len)))
}
}
if (grepl(",", keep_str)) {
parts  <- as.numeric(strsplit(keep_str, ",")[[1]])
return(parts[!is.na(parts)])
}
idx  <- as.numeric(keep_str)
if (!is.na(idx)) return(idx)
return(NULL)
}
parse_bed_name  <- function(name_vec, exon_sep = " ", gene_field_index = NULL,
gene_name_keep = NULL, auto_exon = TRUE, gene_name_collapse = " ") {
if (is.null(exon_sep) || length(exon_sep) == 0 || any(exon_sep == "")) exon_sep  <- " "
if (is.null(gene_name_collapse) || gene_name_collapse == "") gene_name_collapse  <- " "
if (length(exon_sep) > 1) {
   escaped_seps <- vapply(exon_sep, function(x) gsub("([\\\\^\\$\\.\\|\\?\\*\\+\\(\\)\\[\\{\\]\\}])", "\\\\\\1", x), character(1))
   split_pat <- paste(escaped_seps, collapse = "|")
   use_fixed <- FALSE
} else {
   if (grepl("[\\[\\]\\(\\)\\|\\.\\*\\+\\?\\^\\$]", exon_sep, perl = TRUE)) {
     split_pat <- exon_sep
     use_fixed <- FALSE
   } else {
     split_pat <- exon_sep
     use_fixed <- TRUE
   }
}
g_clean <- gsub("\\s*\\(.*?\\)", "", as.character(name_vec))
g_clean <- sub("^([^,]+),.*$", "\\1", g_clean)
exon_pattern <- "(?<![A-Za-z0-9])ex(?:on)?[_-]?([0-9]+)"
exon_all <- gregexpr(exon_pattern, g_clean, perl = TRUE)
exon_numbers <- vapply(seq_along(g_clean), function(i) {
    m <- exon_all[[i]]
    if (length(m) == 0 || m[1] == -1) return(NA_integer_)
    cap_start <- attr(m, "capture.start")[, 1]
    cap_len   <- attr(m, "capture.length")[, 1]
    last <- length(m)
    digits <- substr(g_clean[i], cap_start[last], cap_start[last] + cap_len[last] - 1)
    suppressWarnings(as.integer(digits))
}, integer(1))
parts_list <- strsplit(g_clean, split = split_pat, fixed = use_fixed)
gene_names <- vapply(seq_along(parts_list), function(i) {
  g_orig <- as.character(name_vec[i])
  if (is.na(g_orig) || g_orig == "" || g_orig == ".") return(NA_character_)
  parts <- parts_list[[i]]
  parts <- parts[parts != ""]
  keep_idx <- NULL
  if (!is.null(gene_name_keep)) {
    keep_str <- as.character(gene_name_keep)
    if (grepl("-", keep_str)) {
      p <- as.numeric(strsplit(keep_str, "-")[[1]])
      if (length(p) == 2 && !any(is.na(p))) keep_idx <- seq(p[1], min(p[2], length(parts)))
    } else if (grepl(",", keep_str)) {
      keep_idx <- as.numeric(strsplit(keep_str, ",")[[1]])
      keep_idx <- keep_idx[!is.na(keep_idx)]
    } else {
      idx <- as.numeric(keep_str)
      if (!is.na(idx)) keep_idx <- idx
    }
  } else if (!is.null(gene_field_index)) {
    idx <- as.numeric(gene_field_index)
    if (!is.na(idx)) keep_idx <- idx
  }
  if (!is.null(keep_idx)) {
    valid_idx <- keep_idx[keep_idx <= length(parts) & keep_idx >= 1]
    if (length(valid_idx) > 0) {
      gene <- paste(parts[valid_idx], collapse = gene_name_collapse)
    } else {
      gene <- parts[length(parts)]
    }
  } else {
    if (grepl("^(NM_|NR_|XM_)", g_orig)) {
      if (length(parts) >= 3) gene <- parts[3] else gene <- parts[length(parts)]
    } else {
      gene <- parts[1]
    }
  }
  if (!is.na(gene) && grepl("^[0-9]+$", gene) && length(parts) > 1) gene <- parts[1]
  return(gene)
}, character(1))
list(gene = gene_names, exon = exon_numbers)
}

select_best_hits <- function(gr1, gr2, hits) {
qh <- S4Vectors::queryHits(hits)
sh <- S4Vectors::subjectHits(hits)
ov <- GenomicRanges::pintersect(gr1[qh], gr2[sh])
w <- IRanges::width(ov)
df_hits <- data.frame(q = qh, s = sh, w = w)
df_hits <- df_hits[order(df_hits$q, -df_hits$w), ]
df_hits <- df_hits[!duplicated(df_hits$q), ]
list(q = df_hits$q, s = df_hits$s)
}

input_df  <- utils::read.table(input_bed, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
if (ncol(input_df)  < 3) stop("Input BED must have at least 3 columns.")
colnames(input_df)[1:3]  <- c("Chr", "Start", "End")
if (ncol(input_df)  >= 4) {
colnames(input_df)[4]  <- "OriginalName"
} else {
input_df$OriginalName  <- NA_character_
}
parsed <- parse_bed_name(input_df$OriginalName, exon_sep, gene_field_index, gene_name_keep, auto_exon_number, gene_name_collapse)
input_df$Gene <- parsed$gene
input_df$ExtractedExon <- parsed$exon
if (!auto_exon_number && region_numbering_mode == "bed_text") {
input_df$ExonNum <- input_df$ExtractedExon
} else {
input_df$ExonNum <- NA_integer_
}

if (bed_process == "REGEN") {
# ... (REGEN logic remains identical to your original) ...
stop("REGEN logic omitted for brevity, please keep your existing REGEN block")
} else if (bed_process == "STANDARD") {
panel_bed_paths <- NULL
if (!is.null(list_genes) && file.exists(list_genes)) panel_bed_paths <- readLines(list_genes)
else if (!is.null(genes_file) && file.exists(genes_file)) panel_bed_paths <- readLines(genes_file)
else if (!is.null(panel_files)) {
if (length(panel_files) == 1 && file.exists(panel_files[1])) panel_bed_paths <- readLines(panel_files[1])
else panel_bed_paths <- panel_files
}
if (!is.null(panel_bed_paths) && length(panel_bed_paths) > 0) {
   panel_list <- lapply(panel_bed_paths, function(f) {
     if (!file.exists(f)) return(NULL)
     tmp <- utils::read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
     if (ncol(tmp) < 4) return(NULL)
     tmp <- tmp[, 1:4]
     colnames(tmp) <- c("Chr", "Start", "End", "OriginalName")
     return(tmp)
   })
   all_panels <- data.table::rbindlist(Filter(Negate(is.null), panel_list), fill = TRUE)
   if (nrow(all_panels) > 0) {
     all_panels <- unique(as.data.frame(all_panels))
     panel_parsed <- parse_bed_name(all_panels$OriginalName, exon_sep, gene_field_index, gene_name_keep, auto_exon_number, gene_name_collapse)
     all_panels$Gene <- panel_parsed$gene
     all_panels$ExtractedExon <- panel_parsed$exon
     if (!auto_exon_number && region_numbering_mode == "bed_text") {
       all_panels$ExonNum <- all_panels$ExtractedExon
     } else {
       all_panels$ExonNum <- NA_integer_
     }
     input_start_1based <- if (bed_zero_based) input_df$Start + 1 else input_df$Start
     bed_gr <- GenomicRanges::GRanges(seqnames = input_df$Chr, ranges = IRanges::IRanges(start = input_start_1based, end = input_df$End))
     
     # FIX: Added OriginalName to GRanges metadata
     panel_gr <- GenomicRanges::GRanges(
       seqnames = all_panels$Chr, ranges = IRanges::IRanges(start = all_panels$Start + 1, end = all_panels$End),
       Gene = all_panels$Gene, ExonNum = all_panels$ExonNum, OriginalName = all_panels$OriginalName
     )
     hits <- GenomicRanges::findOverlaps(bed_gr, panel_gr)
     sel <- select_best_hits(bed_gr, panel_gr, hits)
     qh <- sel$q; sh <- sel$s
     output_start <- GenomicRanges::start(bed_gr[qh])
     if (bed_zero_based) output_start <- output_start - 1
     
     # FIX: Added OriginalName to df
     df <- data.frame(
       Chr = as.character(GenomicRanges::seqnames(bed_gr[qh])), Start = output_start, End = GenomicRanges::end(bed_gr[qh]),
       Gene = panel_gr$Gene[sh], ExonNum = panel_gr$ExonNum[sh], OriginalName = panel_gr$OriginalName[sh], stringsAsFactors = FALSE
     )
     
     names(df)[names(df) == "Chr"] <- "chromosome"
     names(df)[names(df) == "Start"] <- "start"
     names(df)[names(df) == "End"] <- "end"
     names(df)[names(df) == "Gene"] <- "gene"
     df <- handle_off_target_regions(df, pattern = off_target_pattern, handling = off_target_handling)
     df <- assign_exon_numbers_per_gene(df)
     if (!"exon_number" %in% names(df) || length(df$exon_number) != nrow(df)) {
       warning("assign_exon_numbers_per_gene() did not produce a valid exon_number column; falling back to sequential numbering.")
       df$exon_number <- seq_len(nrow(df))
     }
     names(df)[names(df) == "chromosome"] <- "Chr"
     names(df)[names(df) == "start"] <- "Start"
     names(df)[names(df) == "end"] <- "End"
     names(df)[names(df) == "gene"] <- "Gene"
     df$ExonNum <- df$exon_number
     df$exon_number <- NULL
     df$Custom.Exon <- df$ExonNum
   } else {
     # FIX: Added OriginalName to df fallback 1
     df <- data.frame(
       Chr = input_df$Chr, Start = if (bed_zero_based) input_df$Start else input_df$Start - 1, End = input_df$End,
       Gene = input_df$Gene, ExonNum = input_df$ExonNum, OriginalName = input_df$OriginalName, stringsAsFactors = FALSE
     )
     names(df)[names(df) == "Chr"] <- "chromosome"
     names(df)[names(df) == "Start"] <- "start"
     names(df)[names(df) == "End"] <- "end"
     names(df)[names(df) == "Gene"] <- "gene"
     df <- handle_off_target_regions(df, pattern = off_target_pattern, handling = off_target_handling)
     df <- assign_exon_numbers_per_gene(df)
     if (!"exon_number" %in% names(df) || length(df$exon_number) != nrow(df)) {
       warning("assign_exon_numbers_per_gene() did not produce a valid exon_number column; falling back to sequential numbering.")
       df$exon_number <- seq_len(nrow(df))
     }
     names(df)[names(df) == "chromosome"] <- "Chr"
     names(df)[names(df) == "start"] <- "Start"
     names(df)[names(df) == "end"] <- "End"
     names(df)[names(df) == "gene"] <- "Gene"
     df$ExonNum <- df$exon_number
     df$exon_number <- NULL
     df$Custom.Exon <- df$ExonNum
   }
} else {
  # FIX: Added OriginalName to df fallback 2
  df <- data.frame(
    Chr = input_df$Chr, Start = if (bed_zero_based) input_df$Start else input_df$Start - 1, End = input_df$End,
    Gene = input_df$Gene, ExonNum = input_df$ExonNum, OriginalName = input_df$OriginalName, stringsAsFactors = FALSE
  )
  names(df)[names(df) == "Chr"] <- "chromosome"
  names(df)[names(df) == "Start"] <- "start"
  names(df)[names(df) == "End"] <- "end"
  names(df)[names(df) == "Gene"] <- "gene"
  df <- handle_off_target_regions(df, pattern = off_target_pattern, handling = off_target_handling)
  df <- assign_exon_numbers_per_gene(df)
  if (!"exon_number" %in% names(df) || length(df$exon_number) != nrow(df)) {
    warning("assign_exon_numbers_per_gene() did not produce a valid exon_number column; falling back to sequential numbering.")
    df$exon_number <- seq_len(nrow(df))
  }
  names(df)[names(df) == "chromosome"] <- "Chr"
  names(df)[names(df) == "start"] <- "Start"
  names(df)[names(df) == "end"] <- "End"
  names(df)[names(df) == "gene"] <- "Gene"
  df$ExonNum <- df$exon_number
  df$exon_number <- NULL
  df$Custom.Exon <- df$ExonNum
}
} else {
stop("bed_process must be 'STANDARD', 'REGEN', or 'NO'")
}

df$Gene  <- as.character(df$Gene)
df$Gene  <- vapply(strsplit(df$Gene, ",", fixed = TRUE), function(x) x[1], character(1))
if (!unknown_gene) {
df  <- df[!is.na(df$Gene) & df$Gene != "" & df$Gene != ".", ]
} else {
df$Gene[is.na(df$Gene) | df$Gene == ""]  <- "Unknown"
}
df  <- unique(df)
chrom_base  <- c(as.character(1:22), "X", "Y", "M")
if (any(grepl("^chr", df$Chr))) {
chrom_levels  <- paste0("chr", chrom_base)
} else {
chrom_levels  <- chrom_base
}
unique_chrs  <- unique(df$Chr)
final_levels  <- c(chrom_levels, setdiff(unique_chrs, chrom_levels))
df$Chr  <- factor(df$Chr, levels = final_levels)
df <- df[order(df$Chr, df$Start), ]
df$Chr <- as.character(df$Chr)
if ("Custom.Exon" %in% colnames(df)) {
out_df  <- df[, c("Chr", "Start", "End", "Gene", "Custom.Exon")]
} else if ("ExonNum" %in% colnames(df)) {
out_df  <- df[, c("Chr", "Start", "End", "Gene", "ExonNum")]
} else {
out_df  <- df[, c("Chr", "Start", "End", "Gene")]
}
dir.create(dirname(output_bed), showWarnings = FALSE, recursive = TRUE)
utils::write.table(out_df, file = output_bed, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
message("[INFO] BED written: ", output_bed)
invisible(output_bed)
}