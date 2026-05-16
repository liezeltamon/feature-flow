#!/usr/bin/env Rscript

renv::load()
setwd(here::here())
suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(tibble)
  source("code/feature_module_reduction/io_helpers.R")
  source("utils/test_normality.R")
})

# %% Functions

build_normality_summary <- function(results_df, alpha, correct, test) {
  tibble(
    n_features_total = nrow(results_df),
    n_features_tested = sum(!is.na(results_df$p_value)),
    n_features_untested = sum(is.na(results_df$p_value)),
    n_features_normal = sum(results_df$is_normal == 1, na.rm = TRUE),
    n_features_non_normal = sum(results_df$is_normal == 0, na.rm = TRUE),
    prop_features_normal = sum(results_df$is_normal == 1, na.rm = TRUE) / nrow(results_df),
    alpha = alpha,
    multiple_testing_correction = correct,
    test = test
  )
}

# %% Parameters

parser <- ArgumentParser(
  description = "Wrapper around utils/test_normality.R for pipeline-style normality testing."
)
parser$add_argument("--feature_path", required = TRUE, help = "Wide feature table with explicit sample ID column.")
parser$add_argument("--sample_key", required = TRUE, help = "Column name containing sample IDs.")
parser$add_argument("--out_dir", required = TRUE, help = "Output directory for CSV and PDF diagnostics.")
parser$add_argument("--alpha", type = "double", default = 0.05)
parser$add_argument("--correct", type = "character", default = "BH")
parser$add_argument("--test", type = "character", default = "auto")
args <- parser$parse_args()

# %% ----- MAIN -----

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

# A. Load feature table and validate the sample ID column.
feature_info <- read_feature_table_generic(
  args$feature_path,
  sample_key = args$sample_key
)
feature_df <- feature_info$feature_df
numeric_cols <- feature_info$numeric_cols
if (length(numeric_cols) == 0) {
  stop("No numeric feature columns found after excluding sample_key.")
}

# B. Test only numeric feature columns; non-numeric columns are ignored here.
results_df <- test_normality(
  x = feature_df %>% select(all_of(numeric_cols)),
  output_dir = args$out_dir,
  alpha = args$alpha,
  correct = args$correct,
  test = args$test,
  verbose = TRUE
)

# C. Save a one-row summary answering whether most features look normal.
summary_df <- build_normality_summary(
  results_df = results_df,
  alpha = args$alpha,
  correct = args$correct,
  test = args$test
)
write_csv(summary_df, file.path(args$out_dir, "normality_summary.csv"))

message(
  "Finished normality testing. ",
  summary_df$n_features_normal, " / ", summary_df$n_features_total,
  " features classified as normal."
)
