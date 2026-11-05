settings <- list(modechrom = "A", min_corr = 0.5, min_total_reads = 20000)

build_run <- function(d, ...) {
  run <- make_synthetic_run(d, ...)
  qc  <- file.path(d, "qc.tsv")
  suppressMessages(run_qc_metrics(run$cov, min_corr = 0.5, min_total_reads = 20000,
                                  output_file = qc, modechrom = "A"))
  c(run, list(qc = qc))
}

test_that("the global report is one self-contained, multi-sample HTML file", {
  skip_if_not_installed("jsonlite")
  d <- withr::local_tempdir()
  run <- build_run(d)
  f <- suppressMessages(generate_report(run$summary, run$qc, file.path(d, "out"), settings,
                                        global = TRUE, prefix = "T", coverage_rdata = run$cov))
  expect_true(file.exists(f))
  expect_match(basename(f), "^ECHO_T_report\\.html$")
  html <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  expect_false(grepl("shiny", html, ignore.case = TRUE))
  expect_match(html, "id=\"echo-data\"", fixed = TRUE)
  expect_false(grepl("<script src=\"https://", html, fixed = TRUE))   # plotly is inlined

  p <- report_json_of(f)$data
  expect_equal(p$meta$modechrom, "A")
  expect_equal(length(p$samples$data$Sample), 8)
  expect_equal(p$cnvs$n, 1)
  expect_equal(unlist(p$cnvs$data$Type), "DEL")                      # ExomeDepth's "deletion" -> DEL
  expect_equal(unlist(p$samples$data$DEL)[match("T01", unlist(p$samples$data$Sample))], 1)
  expect_true("1" %in% names(p$windows))                              # coverage profile embedded
  expect_true(length(p$windows[["1"]]$ref_names) == 5)
  # the scope warning is in the alerts, whichever QC file was used
  titles <- vapply(p$alerts, function(a) a$title, "")
  expect_true(any(grepl("Gonosomes .* NOT analysed", titles)))
  expect_true(any(vapply(p$alerts, function(a) a$kind == "sample" && grepl("T03", a$text), TRUE)))
  # status of the shallow sample
  st <- unlist(p$samples$data$Status)[match("T03", unlist(p$samples$data$Sample))]
  expect_equal(st, "Warning")
})

test_that("missing chromosome classes and a missing QC file are reported", {
  d <- withr::local_tempdir()
  run <- build_run(d, with_x = FALSE)
  f <- suppressMessages(generate_report(run$summary, NULL, file.path(d, "out"), list(modechrom = "A")))
  p <- report_json_of(f)$data
  titles <- vapply(p$alerts, function(a) a$title, "")
  expect_true(any(grepl("missing from the BED", titles)))
  expect_true(any(grepl("No QC metrics file", titles)))
  expect_false(p$meta$qc_available)
})

test_that("mode X flags autosomes as not analysed", {
  d <- withr::local_tempdir()
  run <- build_run(d)
  f <- suppressMessages(generate_report(run$summary, run$qc, file.path(d, "out"), list(modechrom = "X")))
  p <- report_json_of(f)$data
  titles <- vapply(p$alerts, function(a) a$title, "")
  expect_true(any(grepl("Autosomes are NOT analysed", titles)))
})

test_that("per-sample reports only carry that sample's data", {
  d <- withr::local_tempdir()
  run <- build_run(d)
  files <- suppressMessages(generate_report(run$summary, run$qc, file.path(d, "out"), settings,
                                            global = FALSE, log_file = NULL))
  expect_length(files, 1)                                             # only T01 has a call
  p <- report_json_of(files[1])$data
  expect_equal(p$meta$focus, "T01")
  expect_equal(unlist(p$samples$data$Sample), "T01")
  expect_null(p$pca); expect_null(p$corr); expect_null(p$sample_table)
  # no other sample may appear, except as one of T01's own reference samples
  refs <- unlist(p$refs$T01)
  others <- setdiff(run$samples, c("T01", refs))
  raw <- report_json_of(files[1])$raw
  for (o in others) expect_false(grepl(paste0("\"", o, "\""), raw, fixed = TRUE), info = o)
})

test_that("JSON payload is inert inside the script element", {
  expect_false(grepl("</", ECHO:::report_json(list(x = "</script><!-- x")), fixed = TRUE))
})

test_that("the report works without any CNV call", {
  d <- withr::local_tempdir()
  run <- build_run(d)
  e <- new.env(); load(run$summary, envir = e)
  cnv_calls <- data.frame()
  counts <- e$counts; bed_file <- e$bed_file; refs <- e$refs; models <- e$models
  gender_map <- e$gender_map; sample_names <- e$sample_names
  s0 <- file.path(d, "s0.RData")
  save(counts, cnv_calls, refs, models, bed_file, gender_map, sample_names, file = s0)
  f <- suppressMessages(generate_report(s0, run$qc, file.path(d, "out"), settings))
  p <- report_json_of(f)$data
  expect_equal(p$cnvs$n, 0)
  expect_true(any(vapply(p$alerts, function(a) a$title == "No CNV called", TRUE)))
})
