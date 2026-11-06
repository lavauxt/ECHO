# Run from the ECHO package root with:
#   source("R/echo_utils.R")
#   source("test_gene_parser_regression.R")

if (!exists("extract_bed_gene_symbols", mode = "function")) {
    stop("Source R/echo_utils.R before running this regression script.")
}

bed_names <- c(
    "SCG5_GREM1(120/20883)_chr15_32990274_32990394",
    paste0(
        "NM_000059_BRCA2_ex11(251/4972)_chr13_32911039_32911290,",
        "NM_000059_BRCA2_ex11(317/4972)_chr13_32910347_32910698,",
        "NM_000059_BRCA2_ex11(335/4972)_chr13_32910699_32911034,",
        "NM_000059_BRCA2_ex11(4059/4972)_chr13_32911294_32915399"
    ),
    "NM_000059_BRCA2_ex12_chr13_32918659_32918824",
    "NM_000059_BRCA2_ex23,NM_000059_BRCA2_ex24_chr13_32953836_32954379",
    "NM_133509_RAD51B_5utr1-ex2_chr14_68290169_68290409",
    "NM_001144757_SCG5_ex3_chr15_32971912_32972143",
    "HorsROI_chr14_68288843_68289083",
    paste0(
        "NM_001144757_SCG5_ex6(225/625),",
        "SCG5_GREM1(39/20883)_chr15_32989215_32989359"
    ),
    "BCR-ABL1_ex1_chr9_100_200"
)

expected <- c(
    "SCG5_GREM1", "BRCA2", "BRCA2", "BRCA2", "RAD51B", "SCG5",
    "HorsROI", "SCG5;SCG5_GREM1", "BCR-ABL1"
)
observed <- extract_bed_gene_symbols(bed_names)
stopifnot(identical(observed, expected))

message("All BED gene-symbol parser regression cases passed.")
