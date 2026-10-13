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

  result <- ECHO:::annotate_bed_regen(input, references, bed_zero_based = TRUE)

  expect_equal(result$gene, c("GENE_A", "GENE_C", "Unknown"))
  expect_true(is.na(result$exon_number[3]))
  expect_equal(result$start, input$Start)
  expect_equal(result$end, input$End)
  expect_identical(result$chromosome, input$Chr)
})

test_that("REGEN handles BED and reference chromosome prefixes", {
  input <- data.frame(Chr = "1", Start = 100, End = 110)
  references <- data.frame(Chr = "chr1", Start = 101, End = 110, Gene = "GENE_A")

  result <- ECHO:::annotate_bed_regen(input, references, bed_zero_based = TRUE)

  expect_equal(result$gene, "GENE_A")
})

test_that("REGEN converts one-based input back to BED start coordinates", {
  input <- data.frame(Chr = "chr1", Start = 101, End = 110)
  references <- data.frame(Chr = "chr1", Start = 101, End = 110, Gene = "GENE_A")

  result <- ECHO:::annotate_bed_regen(input, references, bed_zero_based = FALSE)

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

  result <- ECHO:::annotate_bed_regen(input, references, bed_zero_based = TRUE)

  expect_equal(result$gene, "GENE_FIRST")
})

test_that("empty overlap selection returns no selected hits", {
  query <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
  subject <- GenomicRanges::GRanges("chr2", IRanges::IRanges(1, 10))
  hits <- suppressWarnings(GenomicRanges::findOverlaps(query, subject))

  expect_identical(ECHO:::select_best_hits(query, subject, hits),
                   list(q = integer(), s = integer()))
})

test_that("REGEN rejects invalid coordinates", {
  input <- data.frame(Chr = "chr1", Start = -1, End = 10)
  references <- data.frame(Chr = "chr1", Start = 1, End = 10, Gene = "GENE_A")

  expect_error(
    ECHO:::annotate_bed_regen(input, references, bed_zero_based = TRUE),
    "invalid coordinates"
  )
})

test_that("REGEN preserves unmatched intervals through the exported API", {
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")

  bed_file <- tempfile(fileext = ".bed")
  output_file <- tempfile(fileext = ".bed")
  on.exit(unlink(c(bed_file, output_file)))
  writeLines("chrZZZ\t0\t10", bed_file)

  process_bed_file(
    bed_file,
    output_file,
    bed_process = "REGEN",
    genome_version = "hg19"
  )

  result <- utils::read.table(output_file, sep = "\t", header = FALSE,
                              stringsAsFactors = FALSE)
  expect_equal(nrow(result), 1L)
  expect_equal(result[[4]], "Unknown")
})

test_that("BED validation accepts files with no name column", {
  bed_file <- tempfile(fileext = ".bed")
  on.exit(unlink(bed_file))
  writeLines("chr1\t0\t10", bed_file)

  expect_no_error(validate_bed_regions(bed_file, verbose = FALSE))
})

test_that("BED validation identifies malformed and invalid coordinates", {
  bed_file <- tempfile(fileext = ".bed")
  on.exit(unlink(bed_file))
  writeLines(c(
    "chr1\t0\t10",
    "chr1\t5\t5",
    "chr2\t-1\t10",
    "chr3\tstart\t20",
    "chr4\t1.5\t20"
  ), bed_file)

  issues <- validate_bed_regions(bed_file, verbose = FALSE)

  expect_equal(nrow(issues$invalid_coordinates), 4L)
})

test_that("panel BED processing preserves input intervals without panel overlaps", {
  input_bed <- tempfile(fileext = ".bed")
  panel_bed <- tempfile(fileext = ".bed")
  output_bed <- tempfile(fileext = ".bed")
  on.exit(unlink(c(input_bed, panel_bed, output_bed)))
  writeLines(c("chr1\t0\t10\tGENE_A exon_1", "chr1\t30\t40\tGENE_B exon_1"), input_bed)
  writeLines("chr1\t0\t10\tGENE_A exon_1", panel_bed)

  process_bed_file(input_bed, output_bed, bed_process = "STANDARD",
                   panel_files = panel_bed, exon_sep = " ", unknown_gene = TRUE)
  result <- utils::read.table(output_bed, sep = "\t", header = FALSE,
                              stringsAsFactors = FALSE)

  expect_equal(nrow(result), 2L)
  expect_equal(result[[1]], c("chr1", "chr1"))
  expect_equal(result[[2]], c(0, 30))
  expect_equal(result[[3]], c(10, 40))
  expect_equal(result[[4]], c("GENE_A", "Unknown"))
})

test_that("panel BED processing preserves all intervals when there are no overlaps", {
  input_bed <- tempfile(fileext = ".bed")
  panel_bed <- tempfile(fileext = ".bed")
  output_bed <- tempfile(fileext = ".bed")
  on.exit(unlink(c(input_bed, panel_bed, output_bed)))
  writeLines(c("chr1\t0\t10\tGENE_A exon_1", "chr1\t30\t40\tGENE_B exon_1"), input_bed)
  writeLines("chr2\t0\t10\tOTHER exon_1", panel_bed)

  process_bed_file(input_bed, output_bed, bed_process = "STANDARD",
                   panel_files = panel_bed, exon_sep = " ")
  result <- utils::read.table(output_bed, sep = "\t", header = FALSE,
                              stringsAsFactors = FALSE)

  expect_equal(nrow(result), 2L)
  expect_equal(result[[2]], c(0, 30))
  expect_equal(result[[3]], c(10, 40))
  expect_equal(result[[4]], c("Unknown", "Unknown"))
})

test_that("BED processing preserves duplicate and unnamed target intervals", {
  input_bed <- tempfile(fileext = ".bed")
  output_bed <- tempfile(fileext = ".bed")
  on.exit(unlink(c(input_bed, output_bed)))
  writeLines(c("chr1\t0\t10", "chr1\t0\t10"), input_bed)

  suppressWarnings(process_bed_file(input_bed, output_bed, bed_process = "STANDARD",
                                    bed_zero_based = TRUE))
  result <- utils::read.table(output_bed, sep = "\t", header = FALSE,
                              stringsAsFactors = FALSE)

  expect_equal(nrow(result), 2L)
  expect_equal(result[[1]], c("chr1", "chr1"))
  expect_equal(result[[2]], c(0, 0))
  expect_equal(result[[3]], c(10, 10))
  expect_equal(result[[4]], c("Unknown", "Unknown"))
})

test_that("BED processing stops rather than dropping invalid intervals", {
  input_bed <- tempfile(fileext = ".bed")
  output_bed <- tempfile(fileext = ".bed")
  on.exit(unlink(c(input_bed, output_bed)))
  writeLines(c("chr1\t0\t10", "chr1\t5\t5"), input_bed)

  expect_error(process_bed_file(input_bed, output_bed), "invalid_coordinates")
  expect_false(file.exists(output_bed))
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
