test_that("chromosome_class recognises autosomes, gonosomes and other contigs", {
  expect_equal(ECHO:::chromosome_class(c("chr1", "22", "chrX", "Y", "chrM", "chr23", "chrUn_gl000")),
               c("autosome", "autosome", "X", "Y", "other", "other", "other"))
})

test_that("scope warnings flag the classes the selected modechrom does NOT analyse", {
  chr <- c(rep("chr1", 3), rep("chrX", 2), "chrY")
  w_a <- ECHO:::chromosome_scope_warnings(ECHO:::chromosome_scope(chr, "A"), "A")
  expect_true(any(w_a$Type == "Analysis scope" & grepl("^Gonosomes NOT analysed", w_a$Details)))
  expect_false(any(grepl("^Autosomes .*NOT analysed", w_a$Details)))
  expect_false("Missing Chromosome" %in% w_a$Type)

  w_x <- ECHO:::chromosome_scope_warnings(ECHO:::chromosome_scope(chr, "X"), "X")
  expect_true(any(grepl("^Autosomes .*NOT analysed", w_x$Details)))
  expect_true(any(grepl("chrY", w_x$Details) & grepl("^Gonosomes", w_x$Details)))

  w_m <- ECHO:::chromosome_scope_warnings(ECHO:::chromosome_scope(rep("chr1", 4), "A"), "A")
  expect_equal(sum(w_m$Type == "Missing Chromosome"), 2L)
  expect_false(any(w_m$Type == "Analysis scope"))
})

test_that("run_qc_metrics reports severities, scope and exact duplicates", {
  d <- withr::local_tempdir()
  run <- make_synthetic_run(d)
  out <- file.path(d, "qc.tsv")
  qc <- suppressMessages(run_qc_metrics(run$cov, min_corr = 0.5, min_cov = 100, min_total_reads = 20000,
                                        output_file = out, modechrom = "A"))
  expect_true(all(c("Sample", "Exon", "Type", "Details", "Gene", "Severity") %in% colnames(qc)))
  expect_true(all(qc$Severity %in% c("error", "warning", "info")))
  expect_true(any(qc$Sample == "T03" & qc$Type == "Whole sample"))
  expect_true(any(qc$Type == "Analysis scope" & grepl("Gonosomes", qc$Details)))
  expect_equal(nrow(data.table::fread(out)), nrow(qc))

  # a duplicated sample has r = 1: it must not be dropped from the correlation
  e <- new.env(); load(run$cov, envir = e)
  counts <- e$counts; counts$T02 <- counts$T01
  bed_file <- e$bed_file; sample_names <- e$sample_names
  dup <- file.path(d, "dup.RData"); save(counts, bed_file, sample_names, file = dup)
  qc2 <- suppressMessages(run_qc_metrics(dup, min_corr = 1.01, output_file = file.path(d, "qc2.tsv")))
  expect_true(all(c("T01", "T02") %in% qc2$Sample[grepl("Low correlation", qc2$Details)]))
  expect_false(any(grepl("NA", qc2$Details[qc2$Sample %in% c("T01", "T02")])))
})

test_that("small cohorts raise a Cohort size warning", {
  d <- withr::local_tempdir()
  run <- make_synthetic_run(d, n_samples = 2)
  qc <- suppressMessages(run_qc_metrics(run$cov, output_file = file.path(d, "qc.tsv")))
  expect_true("Cohort size" %in% qc$Type)
})
