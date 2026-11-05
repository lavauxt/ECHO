# Small synthetic cohort used by the QC / report tests: 5 chr1 genes, 3 chrX
# genes, a handful of samples, one shallow sample and one DEL call.
make_synthetic_run <- function(dir, with_x = TRUE, n_samples = 8, seed = 1) {
  set.seed(seed)
  genes <- list(chr1 = paste0("GA", 1:5), chr2 = paste0("GB", 1:3))
  if (with_x) genes$chrX <- paste0("GX", 1:3)
  rows <- list()
  for (ch in names(genes)) {
    pos <- 1e6
    for (g in genes[[ch]]) {
      for (e in 1:5) {
        rows[[length(rows) + 1]] <- data.frame(chromosome = ch, start = pos, end = pos + 150,
                                               gene = g, stringsAsFactors = FALSE)
        pos <- pos + 1000
      }
      pos <- pos + 20000
    }
  }
  bed_file <- do.call(rbind, rows)
  bed_file$exon_number <- stats::ave(seq_len(nrow(bed_file)), bed_file$gene, FUN = seq_along)
  samples <- sprintf("T%02d", seq_len(n_samples))
  mu <- stats::rlnorm(nrow(bed_file), log(400), 0.3)
  m <- sapply(samples, function(s) stats::rnbinom(nrow(bed_file), mu = mu, size = 80))
  if ("T03" %in% colnames(m)) m[, "T03"] <- stats::rbinom(nrow(bed_file), m[, "T03"], 0.05)  # shallow sample
  del <- which(bed_file$gene == "GA2")[2:4]
  m[del, "T01"] <- stats::rbinom(3, m[del, "T01"], 0.5)
  counts <- cbind(bed_file[, c("chromosome", "start", "end")],
                  exon = paste0(bed_file$gene, "_", bed_file$exon_number),
                  GC = 0.5, as.data.frame(m), exon_number = bed_file$exon_number)
  sample_names <- samples
  cov_file <- file.path(dir, "coverage.RData")
  save(counts, bed_file, sample_names, file = cov_file)

  cnv_calls <- data.frame(
    CNV.ID = 1L, Sample = "T01", Type = "deletion", Nexons = 3L,
    Start = bed_file$start[del[1]], End = bed_file$end[del[3]], Chromosome = "chr1",
    BF = 20, Reads.expected = 1200, Reads.observed = 600, Reads.ratio = 0.5,
    Gene = "GA2", Confidence = "HIGH", Start.p = 2L, End.p = 4L,
    global_start = del[1], global_end = del[3], terminal_only = FALSE,
    Correlation = 0.99, N.comp = 5L, Comparator.name = "T02, T04, T05",
    Genomic.ID = paste0("chr1:", bed_file$start[del[1]], "-", bed_file$end[del[3]]),
    stringsAsFactors = FALSE)
  refs <- stats::setNames(lapply(samples, function(s) utils::head(setdiff(samples, c(s, "T03")), 5)), samples)
  models <- stats::setNames(lapply(samples, function(s) c(0.02, 0.5)), samples)
  gender_map <- stats::setNames(rep(c("F", "M"), length.out = n_samples), samples)
  summary_file <- file.path(dir, "summary.RData")
  save(counts, cnv_calls, refs, models, bed_file, gender_map, sample_names, file = summary_file)
  list(cov = cov_file, summary = summary_file, samples = samples)
}

report_json_of <- function(html_file) {
  txt <- paste(readLines(html_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  m <- regmatches(txt, regexpr("<script id=\"echo-data\" type=\"application/json\">.*?</script>", txt, perl = TRUE))
  json <- sub("</script>$", "", sub("^<script id=\"echo-data\" type=\"application/json\">", "", m))
  list(raw = json, data = jsonlite::fromJSON(gsub("<\\/", "</", json, fixed = TRUE), simplifyVector = FALSE))
}
