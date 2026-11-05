test_that("bam_has_index accepts both .bam.bai and .bai", {
  d <- withr::local_tempdir()
  b1 <- file.path(d, "a.bam"); file.create(b1)
  expect_false(ECHO:::bam_has_index(b1))
  file.create(paste0(b1, ".bai")); expect_true(ECHO:::bam_has_index(b1))
  b2 <- file.path(d, "b.bam"); file.create(b2); file.create(file.path(d, "b.bai"))
  expect_true(ECHO:::bam_has_index(b2))
})

test_that("a call covering a whole short gene is not 'terminal only'", {
  bed <- data.frame(chromosome = "chr1", start = c(1, 100, 200, 300, 400), end = c(50, 150, 250, 350, 450),
                    gene = c("A", "A", "B", "B", "B"), stringsAsFactors = FALSE)
  calls <- data.frame(global_start = c(1L, 3L, 3L), global_end = c(2L, 3L, 5L))
  res <- flag_terminal_only_calls(calls, bed)
  expect_equal(res$terminal_only, c(FALSE, TRUE, FALSE))
})

test_that("coloured PCA no longer fails", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  counts <- as.data.frame(matrix(rpois(60 * 6, 100), ncol = 6, dimnames = list(NULL, paste0("S", 1:6))))
  out <- withr::local_tempfile(fileext = ".pdf")
  res <- suppressWarnings(suppressMessages(plot_coverage_pca(counts, paste0("S", 1:6), output_pdf = out,
                                            color_by = rep(c("F", "M"), 3))))
  expect_true(file.exists(out))
})

test_that("BED validation understands one-based intervals", {
  f <- withr::local_tempfile(fileext = ".bed")
  writeLines(c("chr1\t100\t100\tA", "chr1\t200\t300\tB"), f)
  expect_null(suppressMessages(validate_bed_regions(f, verbose = FALSE, bed_zero_based = FALSE))$invalid_coordinates)
  expect_false(is.null(suppressMessages(validate_bed_regions(f, verbose = FALSE, bed_zero_based = TRUE))$invalid_coordinates))
})
