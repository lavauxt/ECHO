# Per-chromosome depth and sex check (shared by ECHO and CANOPE).

make_xy_counts <- function() {
  set.seed(3)
  chr <- c(rep("chr1", 40), rep("chr2", 40), rep("chrX", 20), rep("chrY", 6))
  mu <- rep(200, length(chr))
  male <- c(M1 = TRUE, M2 = TRUE, F1 = FALSE, F2 = FALSE)
  m <- vapply(male, function(is_male) {
    f <- ifelse(chr == "chrX", if (is_male) 0.5 else 1, ifelse(chr == "chrY", if (is_male) 0.6 else 0.01, 1))
    rpois(length(chr), mu * f)
  }, numeric(length(chr)))
  list(counts = as.data.frame(m), chr = chr)
}

test_that("cnv_chrom_depth is relative to the autosomal median and keeps input order", {
  x <- make_xy_counts()
  cd <- cnv_chrom_depth(x$counts, x$chr)
  expect_equal(rownames(cd), c("chr1", "chr2", "chrX", "chrY"))
  expect_equal(colnames(cd), colnames(x$counts))
  expect_true(all(abs(cd[c("chr1", "chr2"), ] - 1) < 0.1))
  expect_true(all(abs(cd["chrX", c("M1", "M2")] - 0.5) < 0.1))
  expect_true(all(abs(cd["chrX", c("F1", "F2")] - 1) < 0.1))
  expect_null(cnv_chrom_depth(x$counts[1:10, ], x$chr))                      # length mismatch
  expect_null(cnv_chrom_depth(x$counts[x$chr == "chrX", ], x$chr[x$chr == "chrX"]))   # no autosomes
})

test_that("cnv_sex_check infers sex and flags annotation / mode mismatches", {
  x <- make_xy_counts(); cd <- cnv_chrom_depth(x$counts, x$chr)
  inf <- cnv_infer_sex(cd)
  expect_equal(inf$inferred_sex, c("M", "M", "F", "F"))
  expect_null(cnv_infer_sex(cd[rownames(cd) != "chrX", ]))
  st <- data.frame(sample_name = c("M1", "M2", "F1"), gender = c("male", "F", "XX"), stringsAsFactors = FALSE)
  lg <- cnv_logger(tempfile(fileext = ".log"), "INFO")
  sx <- cnv_sex_check(cd, st, expected_sex = "F", analysed = c("M1", "F1", "F2"), mode = "XX", lg = lg)
  expect_equal(sx$annotated_sex, c("M", "F", "F", NA))
  expect_equal(sx$mismatch, c(FALSE, TRUE, FALSE, FALSE))          # M2 annotated F but looks male
  expect_equal(sx$mode_mismatch, c(TRUE, FALSE, FALSE, FALSE))     # the annotation wins over the inference: only M1 is male in a female-mode run
  fl <- report_sex_flags(sx, analysed = c("M1", "F1", "F2"), modechrom = "XX")
  expect_equal(sort(fl$Sample), c("M1", "M2"))                      # M2 is not analysed: only its annotation mismatch shows
  expect_true(all(fl$Type == "Whole sample"))
  expect_null(report_sex_flags(NULL, character(), "A"))
})

test_that("report_chrom_depth keeps the natural chromosome order and the report shape", {
  x <- make_xy_counts(); cd <- cnv_chrom_depth(x$counts, rev(x$chr))
  rp <- report_chrom_depth(cd)
  expect_equal(as.character(rp$chromosomes), c("chr1", "chr2", "chrX", "chrY"))
  expect_equal(length(rp$z), 4L)
  expect_equal(length(rp$z[[1]]), ncol(x$counts))
  expect_null(report_chrom_depth(NULL))
})
