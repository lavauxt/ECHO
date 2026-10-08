test_that("REGEN keeps unmatched intervals and annotates largest overlaps", {
  input <- data.frame(
    Chr = c("chr1", "1", "chr9"),
    Start = c(100, 200, 0),
    End = c(160, 220, 25)
  )
  references <- data.frame(
    Chr = c("1", "1", "1"),
    Start = c(101, 120, 201),
    End = c(150, 155, 220),
    Gene = c("GENE_A", "GENE_B", "GENE_C"),
    Transcript = c("NM_A", "NM_B", "NM_C"),
    stringsAsFactors = FALSE
  )

  result <- annotate_bed_regen(input, references, bed_zero_based = TRUE)

  expect_equal(result$gene, c("GENE_A", "GENE_C", "Unknown"))
  expect_true(is.na(result$exon_number[3]))
  expect_equal(result$start, input$Start)
  expect_equal(result$end, input$End)
  expect_identical(result$chromosome, input$Chr)
})

test_that("REGEN handles BED and reference chromosome prefixes", {
  input <- data.frame(Chr = "1", Start = 100, End = 110)
  references <- data.frame(Chr = "chr1", Start = 101, End = 110, Gene = "GENE_A")

  result <- annotate_bed_regen(input, references, bed_zero_based = TRUE)

  expect_equal(result$gene, "GENE_A")
})

test_that("REGEN converts one-based input back to BED start coordinates", {
  input <- data.frame(Chr = "chr1", Start = 101, End = 110)
  references <- data.frame(Chr = "chr1", Start = 101, End = 110, Gene = "GENE_A")

  result <- annotate_bed_regen(input, references, bed_zero_based = FALSE)

  expect_equal(result$start, 100)
  expect_equal(result$end, 110)
  expect_equal(result$gene, "GENE_A")
})

test_that("equal overlap ties resolve deterministically", {
  input <- data.frame(Chr = "chr1", Start = 100, End = 110)
  references <- data.frame(
    Chr = c("chr1", "chr1"),
    Start = c(101, 101),
    End = c(110, 110),
    Gene = c("GENE_FIRST", "GENE_SECOND")
  )

  result <- annotate_bed_regen(input, references, bed_zero_based = TRUE)

  expect_equal(result$gene, "GENE_FIRST")
})

test_that("REGEN rejects invalid coordinates", {
  input <- data.frame(Chr = "chr1", Start = -1, End = 10)
  references <- data.frame(Chr = "chr1", Start = 1, End = 10, Gene = "GENE_A")

  expect_error(
    annotate_bed_regen(input, references, bed_zero_based = TRUE),
    "invalid coordinates"
  )
})

test_that("BED validation accepts files with no name column", {
  bed_file <- tempfile(fileext = ".bed")
  on.exit(unlink(bed_file))
  writeLines("chr1\t0\t10", bed_file)

  expect_no_error(validate_bed_regions(bed_file, verbose = FALSE))
})

test_that("REGEN reports a missing genome-specific TxDb", {
  skip_if(requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE))

  bed_file <- tempfile(fileext = ".bed")
  output_file <- tempfile(fileext = ".bed")
  on.exit(unlink(c(bed_file, output_file)))
  writeLines("chr1\t0\t10", bed_file)

  expect_error(
    process_bed_file(
      bed_file,
      output_file,
      bed_process = "REGEN",
      genome_version = "hg19"
    ),
    "TxDb.Hsapiens.UCSC.hg19.knownGene"
  )
})
