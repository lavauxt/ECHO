# Inputs of the BAM-counting step and the fingerprint behind echo(reuse_coverage = TRUE).

touch <- function(path, txt = "x") { writeLines(txt, path); path }

test_that("echo_input_bams resolves test and reference BAMs like run_bam_coverage()", {
  d <- tempfile("bams"); dir.create(d)
  b <- c(touch(file.path(d, "a.bam")), touch(file.path(d, "b.bam")))
  expect_equal(sort(echo_input_bams(bamdir = d)), sort(b))
  expect_equal(echo_input_bams(bamfiles = b), b)
  expect_equal(echo_input_bams(bamfiles = b[1]), b[1])
  lst <- touch(file.path(d, "list.txt"), b)
  expect_equal(echo_input_bams(bamfiles = lst), b)
  rt <- file.path(d, "refs.tsv"); writeLines(c("bam", "r1.bam", b[2]), rt)
  out <- echo_input_bams(bamfiles = b[1], rbams = rt)
  expect_true(file.path(normalizePath(d), "r1.bam") %in% out)       # relative path anchored to the table
  expect_true(b[2] %in% out)                                         # absolute path kept
  expect_equal(echo_input_bams(), character())
})

test_that("the coverage fingerprint changes with an input file or a counting setting, not otherwise", {
  d <- tempfile("fp"); dir.create(d)
  bam <- touch(file.path(d, "a.bam")); bed <- touch(file.path(d, "t.bed"), "chr1\t1\t100")
  set <- list(pad = 0, genome = "hg38")
  f0 <- echo_coverage_fingerprint(c(bam, bed), set)
  expect_equal(echo_coverage_fingerprint(c(bam, bed), set), f0)
  expect_equal(echo_coverage_fingerprint(c(bed, bam, bam), set), f0)          # order / duplicates do not matter
  expect_false(identical(echo_coverage_fingerprint(c(bam, bed), list(pad = 50, genome = "hg38")), f0))
  Sys.sleep(1.1); writeLines("changed", bed); Sys.setFileTime(bed, Sys.time() + 5)
  expect_false(identical(echo_coverage_fingerprint(c(bam, bed), set), f0))
  expect_true(nzchar(echo_coverage_fingerprint(character(), list())))
})

test_that("the targets layer's helpers read the call table and list the BAMs", {
  d <- tempfile("tg"); dir.create(d)
  tb <- file.path(d, "ECHO_X_CNV_calls.tsv")
  utils::write.table(data.frame(Sample = "S1", Type = "deletion", Confidence = "HIGH"), tb, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_equal(echo_targets_calls(c(file.path(d, "ECHO_X.log"), tb))$Sample, "S1")
  expect_null(echo_targets_calls(file.path(d, "ECHO_X.log")))
  b <- touch(file.path(d, "a.bam"))
  cfg <- file.path(d, "config.yaml")
  yaml::write_yaml(list(input = list(bamfiles = b, bed = "x.bed"), output = list(dir = d, prefix = "X")), cfg)
  expect_equal(echo_targets_bams(cfg), b)
})

test_that("echo_targets() builds the DAG, or says what is missing", {
  if (requireNamespace("targets", quietly = TRUE)) {
    nm <- vapply(echo_targets("config.yaml"), function(t) t$settings$name, "")
    expect_equal(nm, c("echo_config", "echo_bams", "echo_run", "echo_calls"))
  } else {
    expect_error(echo_targets("config.yaml"), "targets")
  }
})
