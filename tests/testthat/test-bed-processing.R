test_that("annotate_bed_regen labels intervals by best overlap and numbers exons per transcript", {
  skip_if_not_installed("GenomicRanges")
  ref <- data.frame(Chr = "chr1", Start = c(100, 300, 500, 100), End = c(200, 400, 600, 200),
                    Strand = c("+", "+", "+", "-"), Transcript = c("NM_1", "NM_1", "NM_1", "NM_2"),
                    Gene = c("G1", "G1", "G1", "G2"), stringsAsFactors = FALSE)
  inp <- data.frame(Chr = "chr1", Start = c(99, 299, 700), End = c(200, 400, 800),
                    OriginalName = c("a", "b", "c"), stringsAsFactors = FALSE)
  out <- annotate_bed_regen(inp, ref, bed_zero_based = TRUE)
  expect_equal(nrow(out), 3L)
  expect_equal(out$gene[2:3], c("G1", "Unknown"))
  expect_equal(out$ExonNum[2], 2L)
  expect_true(is.na(out$ExonNum[3]))
  expect_named(out, c("chromosome", "start", "end", "gene", "ExonNum", "OriginalName"))
})

test_that("process_bed_file STANDARD works with a panel file (select_best_hits was undefined)", {
  skip_if_not_installed("GenomicRanges")
  d <- tempfile(); dir.create(d)
  writeLines(c("chr1\t99\t200\tx", "chr1\t299\t400\ty", "chr1\t999\t1100\tHorsROI_1"), file.path(d, "in.bed"))
  writeLines(c("chr1\t99\t200\tGENEA ex1", "chr1\t299\t400\tGENEA ex2"), file.path(d, "panel.bed"))
  out <- file.path(d, "out.bed")
  suppressMessages(process_bed_file(file.path(d, "in.bed"), out, bed_process = "STANDARD",
                                    panel_files = file.path(d, "panel.bed")))
  x <- read.table(out, sep = "\t", stringsAsFactors = FALSE)
  expect_equal(x$V4[1:2], c("GENEA", "GENEA"))
  expect_equal(x$V4[3], "HorsROI")
})
