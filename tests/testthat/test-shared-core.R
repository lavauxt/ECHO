# Tests of the code shared (byte-identically) by ECHO and CANOPE.
# They run against whichever package they are shipped with.

core <- function(name) get(name, envir = asNamespace(pkg_name()), inherits = FALSE)
pkg_name <- function() {
  d <- testthat::test_path("..", "..", "DESCRIPTION")
  if (file.exists(d)) read.dcf(d, "Package")[1, 1] else "ECHO"
}

make_counts <- function(n_t = 120, n_s = 8, seed = 1) {
  set.seed(seed)
  base <- rpois(n_t, 200)
  m <- sapply(seq_len(n_s), function(i) rpois(n_t, base * runif(1, 0.8, 1.2)))
  colnames(m) <- paste0("S", seq_len(n_s))
  m
}

test_that("detect_outlier_samples flags only noisy samples by default", {
  m <- make_counts()
  m[, "S3"] <- rpois(nrow(m), 5) * sample(c(0, 1, 5), nrow(m), TRUE)   # very noisy
  r <- detect_outlier_samples(m, z_threshold = 3)
  expect_true(r$is_outlier[r$sample == "S3"])
  expect_equal(sum(r$is_outlier), 1L)
  expect_named(r, c("sample", "noise_score", "robust_z", "is_outlier"))
  # sides = "both" can only flag the same or more
  rb <- detect_outlier_samples(m, z_threshold = 3, sides = "both")
  expect_gte(sum(rb$is_outlier), sum(r$is_outlier))
})

test_that("cnv_apply_qc excludes outliers and keeps >= 1 target per chromosome", {
  m <- make_counts()
  m[, "S3"] <- rpois(nrow(m), 5) * sample(c(0, 1, 5), nrow(m), TRUE)
  q <- cnv_apply_qc(m, chromosomes = rep("chr1", 120), gc = rep(0.5, 120))
  expect_equal(q$excluded_samples$sample, "S3")
  expect_false("S3" %in% q$kept_samples)
  q0 <- cnv_apply_qc(m, chromosomes = rep("chr1", 120), sample_qc = FALSE, exon_qc = FALSE)
  expect_equal(nrow(q0$excluded_samples), 0L)
  expect_equal(q0$keep_targets, seq_len(120))
  # a chromosome whose targets are all problematic keeps its best one
  m2 <- make_counts()
  chr <- rep(c("chr1", "chr2"), each = 60)
  m2[61:120, ] <- 1L
  q2 <- cnv_apply_qc(m2, chromosomes = chr, min_exon_mean = 20)
  expect_equal(sum(chr[q2$keep_targets] == "chr2"), 1L)
  expect_true(all(q2$target_qc$reason[!q2$target_qc$problematic & chr == "chr2"] == ""))
})

test_that("cnv_qc_metrics has the unified schema and flags what it should", {
  m <- make_counts()
  tg <- data.frame(chromosome = "chr1", start = seq(100, by = 500, length.out = 120),
                   end = seq(300, by = 500, length.out = 120),
                   GENE = rep(c("A", "B", "C"), each = 40), stringsAsFactors = FALSE)
  m[5, ] <- 0L
  q  <- cnv_apply_qc(m, chromosomes = tg$chromosome)
  r  <- cnv_qc_metrics(m, tg, min_total_reads = 1e9, qc = q)
  expect_named(r, c("Sample", "Exon", "Type", "Details", "Gene", "Severity"))
  expect_true(all(r$Severity %in% c("error", "warning", "info")))
  expect_true("Whole sample" %in% r$Type)                       # low total reads
  expect_true(any(r$Type == "Whole exon" & r$Exon == "A:5"))
  small <- cnv_qc_metrics(m[, 1:2], tg, qc = NULL)
  expect_true("Cohort size" %in% small$Type)
})

test_that("chromosome scope reports what a run left out", {
  chrs <- c(rep("chr1", 5), rep("chrX", 3), "chrM")
  sc <- chromosome_scope(chrs, analysed_classes = "autosome")
  expect_equal(sc$targets[sc$class == "autosome"], 5L)
  expect_true(sc$analysed[sc$class == "autosome"])
  expect_false(sc$analysed[sc$class == "X"])
  w <- chromosome_scope_warnings(sc, hints = c(X = "modechrom = 'XX'"))
  expect_true(any(w$Type == "Analysis scope" & grepl("chrX", w$Details)))
  expect_true(any(w$Type == "Missing Chromosome" & grepl("chrY", w$Details)))
  expect_equal(chromosome_class(c("chr1", "23", "Y", "chrM", "chrUn_x")),
               c("autosome", "X", "Y", "other", "other"))
})

test_that("exon numbering is per gene, genomic, row-order independent, NA for off-target", {
  b <- data.frame(chromosome = c("chr1", "chr1", "chr1", "chr1"),
                  start = c(300, 100, 200, 400), end = c(350, 150, 250, 450),
                  GENE = c("G", "G", "G", NA), stringsAsFactors = FALSE)
  out <- assign_exon_numbers_per_gene(b)
  expect_equal(out$exon_number, c(3L, 1L, 2L, NA))
  expect_equal(nrow(out), 4L)
  lab <- target_labels(out)
  expect_equal(lab[1:3], c("G:3", "G:1", "G:2"))
  expect_match(lab[4], "^chr1:")
})

test_that("off-target handling: na / remove / merge", {
  b <- data.frame(chromosome = "chr1", start = c(100, 200, 300), end = c(150, 250, 350),
                  GENE = c("G1", "HorsROI_1", "G2"), stringsAsFactors = FALSE)
  expect_true(is.na(handle_off_target_regions(b, handling = "na", verbose = FALSE)$GENE[2]))
  expect_equal(nrow(handle_off_target_regions(b, handling = "remove", verbose = FALSE)), 2L)
  expect_equal(handle_off_target_regions(b, handling = "merge", verbose = FALSE)$GENE[2], "G1")
})

test_that("terminal-exon padding never overlaps a neighbour", {
  b <- data.frame(chromosome = "chr1", start = c(100, 205, 1000), end = c(200, 300, 1100),
                  GENE = c("A", "B", "C"), stringsAsFactors = FALSE)
  p <- pad_gene_terminal_exons(b, padding = 50, verbose = FALSE)
  expect_true(all(p$end[-nrow(p)] < p$start[-1]))
  expect_equal(nrow(p), nrow(b))
})

test_that("gene gap positions insert space only at gene boundaries", {
  b <- data.frame(chromosome = "chr1", start = 1:4, end = 2:5, GENE = c("A", "A", "B", "B"))
  g <- compute_gene_gap_positions(b, 1:4, gap = 1)
  expect_equal(g$px, c(1, 2, 4, 5))
  expect_equal(g$gene_group, c(1L, 1L, 2L, 2L))
})

test_that("terminal-only calls are flagged, whole-gene calls are not", {
  b <- data.frame(chromosome = "chr1", start = seq(100, by = 100, length.out = 6),
                  end = seq(150, by = 100, length.out = 6), GENE = "G", stringsAsFactors = FALSE)
  calls <- data.frame(global_start = c(1L, 1L, 3L), global_end = c(1L, 6L, 4L))
  out <- flag_terminal_only_calls(calls, b)
  expect_equal(out$terminal_only, c(TRUE, FALSE, FALSE))
})

test_that("recurrence counts same-type overlapping calls across samples", {
  r <- cnv_recurrence(chr = rep("chr1", 5), start = c(100, 100, 105, 100, 5000),
                      end = c(200, 200, 200, 200, 5100), type = rep("DEL", 5),
                      sample = paste0("S", 1:5), n_samples = 10, min_samples = 3, min_fraction = 0.25)
  expect_equal(r$n_carriers, c(4L, 4L, 4L, 4L, 1L))
  expect_equal(r$recurrent, c(TRUE, TRUE, TRUE, TRUE, FALSE))
  calls <- data.frame(SAMPLE = paste0("S", 1:5), CHR = "chr1", Start = c(100, 100, 105, 100, 5000),
                      End = c(200, 200, 200, 200, 5100), CNV = "DEL", Confidence = "HIGH",
                      stringsAsFactors = FALSE)
  cols <- list(sample = "SAMPLE", chr = "CHR", start = "Start", end = "End", type = "CNV", confidence = "Confidence")
  expect_equal(nrow(suppressMessages(cnv_apply_recurrence(calls, cols, 10, "remove"))), 1L)
  d <- suppressMessages(cnv_apply_recurrence(calls, cols, 10, "downgrade"))
  expect_equal(d$Confidence, c("LOW", "LOW", "LOW", "LOW", "HIGH"))
})

test_that("logger writes timestamped, levelled lines, captures messages and truncates", {
  f <- tempfile(fileext = ".log")
  writeLines("OLD CONTENT", f)
  lg <- cnv_logger(f, tag = "TST", console = FALSE, total_steps = 2)
  lg$header()
  lg$step("first", lg$info("inside"))
  lg$capture({ message("[INFO] hello from a function"); warning("careful"); message("plain") })
  r <- lg$step("second", stop("boom"), on_error = "warn")
  expect_null(r)
  lg$footer(ok = TRUE)
  x <- readLines(f)
  expect_false(any(grepl("OLD CONTENT", x)))
  expect_true(all(grepl("^\\[(DEBUG|INFO|WARNING|ERROR)\\] [0-9]{4}-[0-9]{2}-[0-9]{2} [0-9:]{8} ", x)))
  expect_true(any(grepl("^\\[WARNING\\] .* careful$", x)))
  expect_true(any(grepl("^\\[INFO\\] .* hello from a function$", x)))
  expect_true(any(grepl("Step 1/2", x)) && any(grepl("Step 2/2 failed", x)))
  expect_equal(unname(lg$counts()["WARNING"]), 2L)   # "careful" + failed step
  # append mode keeps previous content
  lg2 <- cnv_logger(f, console = FALSE, append = TRUE); lg2$info("again")
  expect_gt(length(readLines(f)), length(x))
})

test_that("vcf_write: LF endings, unpadded POS, header-only file is valid", {
  f <- tempfile(fileext = ".vcf")
  rec <- data.frame(CHROM = c("chr1", "chr10"), POS = c(5L, 1234567L), ID = c("a", "b"),
                    ALT = c("<DEL>", "<DUP>"), INFO = c("SVTYPE=DEL", "SVTYPE=DUP"), stringsAsFactors = FALSE)
  sf <- matrix(c("0/1:1", "./.:.", "./.:.", "0/1:3"), 2, dimnames = list(NULL, c("S1", "S2")))
  suppressMessages(vcf_write(f, c("##fileformat=VCFv4.2"), c("S1", "S2"), rec, "GT:CN", sf))
  raw <- readBin(f, "raw", file.size(f))
  expect_false(as.raw(13) %in% raw)
  x <- readLines(f)
  expect_equal(x[2], paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "S1", "S2"), collapse = "\t"))
  expect_equal(strsplit(x[3], "\t")[[1]][2], "5")
  g <- tempfile(fileext = ".vcf")
  suppressMessages(vcf_write(g, c("##fileformat=VCFv4.2"), "S1", rec[0, ], "GT", NULL))
  expect_equal(length(readLines(g)), 2L)
})

test_that("palette is complete and group colours are stable", {
  p <- cnv_palette()
  expect_true(all(c("del", "dup", "high", "medium", "low", "fail", "warn", "info", "dark") %in% names(p)))
  expect_false(identical(p$del, p$dup))
  expect_equal(unname(cnv_group_colours(c("F", "M", "unknown"))[3]), "#A6A6A0")
})

test_that("the report CSS tokens are the cnv_palette() colours (PDF, CSS and JS paint with one palette)", {
  p <- cnv_palette(); css <- report_css()
  tok <- regmatches(css, gregexpr("--c-[a-z]+:\\s*#[0-9A-Fa-f]{6}", css))[[1]]
  nm <- sub("^--c-([a-z]+):.*$", "\\1", tok); val <- toupper(sub("^.*:\\s*", "", tok))
  light <- !duplicated(nm)                                   # first definition = light theme
  for (i in which(light)) if (nm[i] %in% names(p)) expect_equal(val[i], toupper(p[[nm[i]]]), info = nm[i])
  for (i in which(!light)) if (nm[i] %in% names(p$dark)) expect_equal(val[i], toupper(p$dark[[nm[i]]]), info = paste("dark", nm[i]))
})

test_that("small utilities", {
  expect_equal(NULL %||% 3, 3)
  expect_equal(sanitize_filename("a b/c:d.e"), "a_b_c_d.e")
  expect_equal(normalize_chromosome_vec(c("1", "chr2"), "chrX"), c("chr1", "chr2"))
  expect_equal(normalize_chromosome_vec(c("1", "chr2"), "X"), c("1", "2"))
})
