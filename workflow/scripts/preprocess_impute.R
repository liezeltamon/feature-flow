# %% Preprocess and impute (optional)
# Saving table again for now if not imputing because preprocessing might change input table but could considering just creating symlink

# %% Setup

renv::load()
setwd(here::here())
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(WGCNA)
  library(vdjremix)
  source("code/feature_module_reduction/io_helpers.R")
  source("code/feature_module_reduction/helpers.R")
  source("utils/benchmark_imputation_variant.R")
  source("utils/filter_features_by_imputation_error_variant.R")
  source("utils/filter_features_by_distribution_shift.R")
})

# %% Functions

build_preprocess_summary <- function(
    n_features_orig,
    n_features_after_preprocess,
    n_features_after_imputation_error_filter,
    n_features_after_distribution_shift_filter,
    n_features_final,
    had_missing_values_before_imputation,
    imputation_requested,
    imputation_applied,
    run_imputation_benchmark,
    used_benchmark_for_feature_filtering,
    imputation_method,
    distribution_shift_filter_enabled,
    distribution_shift_filter_applied,
    distribution_shift_group_key,
    distribution_shift_grouped
) {
  tibble(
    n_features_orig = n_features_orig,
    n_features_after_preprocess = n_features_after_preprocess,
    n_features_after_imputation_error_filter = n_features_after_imputation_error_filter,
    n_features_after_distribution_shift_filter = n_features_after_distribution_shift_filter,
    n_features_final = n_features_final,
    pct_features_after_preprocess = 100 * n_features_after_preprocess / n_features_orig,
    pct_features_after_imputation_error_filter = 100 * n_features_after_imputation_error_filter / n_features_orig,
    pct_features_after_distribution_shift_filter = 100 * n_features_after_distribution_shift_filter / n_features_orig,
    pct_features_final = 100 * n_features_final / n_features_orig,
    had_missing_values_before_imputation = had_missing_values_before_imputation,
    imputation_requested = imputation_requested,
    imputation_applied = imputation_applied,
    run_imputation_benchmark = run_imputation_benchmark,
    used_benchmark_for_feature_filtering = used_benchmark_for_feature_filtering,
    imputation_method = imputation_method,
    distribution_shift_filter_enabled = distribution_shift_filter_enabled,
    distribution_shift_filter_applied = distribution_shift_filter_applied,
    distribution_shift_group_key = if (is.null(distribution_shift_group_key)) {
      NA_character_
    } else {
      distribution_shift_group_key
    },
    distribution_shift_grouped = distribution_shift_grouped
  )
}

# %% Parameters

parser <- ArgumentParser(description = "Preprocess and impute feature table")
parser$add_argument("--feature_path", type = "character", required = TRUE)
parser$add_argument("--metadata_path", type = "character", default = NULL)
parser$add_argument("--sample_key", type = "character", required = TRUE)
parser$add_argument("--out_dir", type = "character", required = TRUE)
parser$add_argument("--missingness_threshold", type = "double", default = 0.2)
parser$add_argument("--min_unique_values", type = "integer", default = NULL)
parser$add_argument("--impute", type = "character", default = "TRUE")
parser$add_argument("--exclude_missing_samples", type = "character", default = "FALSE")
parser$add_argument("--run_imputation_benchmark", type = "character", default = "TRUE")
parser$add_argument("--use_benchmark_for_feature_filtering", type = "character", default = "TRUE")
parser$add_argument("--benchmark_methods", type = "character", default = "mean,missForest")
parser$add_argument("--imputation_method", type = "character", default = NULL)
parser$add_argument("--na_frequencies", type = "character", default = "0.05,0.10,0.15,0.20,0.30,0.40")
parser$add_argument("--error_threshold", type = "double", default = 0.40)
parser$add_argument("--error_metric", type = "character", default = "nrmse_mad")
parser$add_argument("--run_distribution_shift_filter", type = "character", default = "FALSE")
parser$add_argument("--distribution_shift_z_threshold", type = "double", default = 0.1)
parser$add_argument("--distribution_shift_ks_alpha", type = "double", default = 0.05)
parser$add_argument("--distribution_shift_group_key", type = "character", default = NULL)
args <- parser$parse_args()
list2env(args, envir = environment())

# Note that filter_features_by_imputation_error() need benchmark results to work

# vdjremix params
impute = parse_bool(impute, "impute")
exclude_missing_samples = parse_bool(exclude_missing_samples, "exclude_missing_samples")
run_imputation_benchmark = parse_bool(
  run_imputation_benchmark, "run_imputation_benchmark"
)
use_benchmark_for_feature_filtering = parse_bool(
  use_benchmark_for_feature_filtering, "use_benchmark_for_feature_filtering"
)
run_distribution_shift_filter = parse_bool(
  run_distribution_shift_filter, "run_distribution_shift_filter"
)
metadata_path = parse_optional_string(metadata_path)
distribution_shift_group_key = parse_optional_string(distribution_shift_group_key)
benchmark_methods = strsplit(benchmark_methods, ",")[[1]]
na_frequencies = as.numeric(strsplit(na_frequencies, ",")[[1]])

if (!run_imputation_benchmark && use_benchmark_for_feature_filtering) {
  stop(
    "use_benchmark_for_feature_filtering = TRUE requires ",
    "run_imputation_benchmark = TRUE."
  )
}

if (!is.null(distribution_shift_group_key) && is.null(metadata_path)) {
  stop("distribution_shift_group_key requires metadata_path.")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# %% ----- MAIN -----

feature_info <- read_feature_table_generic(feature_path, sample_key = sample_key)
df <- feature_info$feature_df

# Sample for testing
# set.seed(210)
# df <- df %>% slice_sample(n = 50)

mx <- feature_table_to_matrix(df, sample_key, table_label = "input feature table")

n_features_orig <- ncol(mx)

# A. Preprocess

if (is.null(min_unique_values)) {
  min_unique_values <- max(6, ceiling(0.05 * nrow(mx)))
}
message(paste0(
  "Using min_unique_values = ", min_unique_values
))

mat_filtered <- vdjremix::preprocess_features(
    mx,
    missingness_threshold = missingness_threshold,
    min_unique_values = min_unique_values
)

# Print number of features before and after filtering
message(paste0(
  "Number of features before filtering: ", ncol(mx),
  "; after filtering: ", ncol(mat_filtered)
))

n_features_after_preprocess <- ncol(mat_filtered)

# B. Impute

had_missing_values_before_imputation <- any(is.na(mat_filtered))
imputation_requested <- impute
imputation_applied <- FALSE
benchmark <- NULL
imputation_method_used <- NA_character_
mat_filtered2 <- mat_filtered
distribution_shift_filter_enabled <- run_distribution_shift_filter
distribution_shift_filter_applied <- FALSE
distribution_shift_grouped <- (
  run_distribution_shift_filter && !is.null(distribution_shift_group_key)
)
used_benchmark_for_feature_filtering <- FALSE
n_features_after_distribution_shift_filter <- ncol(mat_filtered2)

if (had_missing_values_before_imputation && impute) {

  imputation_applied <- TRUE
  # 1. Benchmark imputation
  if (run_imputation_benchmark) {
    message("Benchmarking imputation methods...")
    t0 <- proc.time()
    benchmark <- benchmark_imputation_variant( # vdjremix::benchmark_imputation
      feature_matrix = mat_filtered,
      na_frequencies = na_frequencies,
      methods = benchmark_methods
    )
    message("benchmark_imputation: ", elapsed_sec(t0), " sec")
    saveRDS(benchmark, file.path(out_dir, "benchmark.rds"))
  } else {
    saveRDS(NULL, file.path(out_dir, "benchmark.rds"))
  }

  # 2. Select best imputation method
  if (is.null(imputation_method)) {
    if (is.null(benchmark)) {
      stop("Either run imputation benchmarking or provide imputation_method.")
    }
    imputation_method <- vdjremix::select_imputation_method(benchmark)$best_method
  }
  imputation_method_used <- imputation_method
  writeLines(
    imputation_method, file.path(out_dir, "chosen_imputation_method.txt")
  )

  # 3. Filter features by imputation error
  if (use_benchmark_for_feature_filtering) {
    message("Filtering features by imputation error...")
    t0 <- proc.time()
    # filtered <- vdjremix::filter_features_by_imputation_error(
    #   feature_matrix = mat_filtered,
    #   benchmark_results = benchmark,
    #   method = imputation_method,
    #   rmse_threshold = error_threshold
    # )
    filtered <- filter_features_by_imputation_error_variant(
      feature_matrix = mat_filtered,
      benchmark_results = benchmark,
      method = imputation_method,
      error_threshold = error_threshold,
      error_metric = error_metric
    )
    message(
      "filter_features_by_imputation_error(): ", elapsed_sec(t0), " sec"
    )
    mat_filtered2 <- filtered$filtered_matrix
    used_benchmark_for_feature_filtering <- TRUE
  } else {
    mat_filtered2 <- mat_filtered
  }

  # 4. Impute features
  message("Imputing features with method: ", imputation_method)
  t0 <- proc.time()
  mat_imputed <- vdjremix::impute_features(
    feature_matrix = mat_filtered2,
    method = imputation_method
  )
  message("impute_features(): ", elapsed_sec(t0), " sec")

  mat_processed <- mat_imputed

  if (run_distribution_shift_filter) {
    message("Filtering features by distribution shift...")
    grouping_vector <- NULL
    if (distribution_shift_grouped) {
      message(
        "Using metadata column for distribution-shift groups: ",
        distribution_shift_group_key
      )
      metadata_df <- read_metadata_aligned(
        metadata_path,
        sample_key = sample_key,
        sample_ids = rownames(mat_filtered2),
        required_cols = distribution_shift_group_key
      )
      grouping_vector <- metadata_df[[distribution_shift_group_key]]
    }

    t0 <- proc.time()
    distribution_shift_filtered <- filter_features_by_distribution_shift(
      pre_imputation_matrix = mat_filtered2,
      post_imputation_matrix = mat_imputed,
      grouping_vector = grouping_vector,
      z_threshold = distribution_shift_z_threshold,
      ks_alpha = distribution_shift_ks_alpha
    )
    message(
      "filter_features_by_distribution_shift(): ",
      elapsed_sec(t0),
      " sec"
    )
    mat_processed <- distribution_shift_filtered$filtered_matrix
    distribution_shift_filter_applied <- TRUE

    write_csv(
      distribution_shift_filtered$summary,
      file.path(out_dir, "distribution_shift_summary.csv")
    )
    if (!is.null(distribution_shift_filtered$group_summary)) {
      write_csv(
        distribution_shift_filtered$group_summary,
        file.path(out_dir, "distribution_shift_group_summary.csv")
      )
    }
  }

  n_features_after_distribution_shift_filter <- ncol(mat_processed)

} else {
  message("No missing values found or imputation not requested; skipping imputation.")
  mat_processed <- mat_filtered
  if (had_missing_values_before_imputation && exclude_missing_samples) {
    n_before <- nrow(mat_processed)
    mat_processed <- mat_processed[complete.cases(mat_processed), ]
    message(
      n_before - nrow(mat_processed), " samples excluded due to missing features; ",
      nrow(mat_processed), " samples retained."
    )
  }
  saveRDS(NULL, file.path(out_dir, "benchmark.rds"))
  n_features_after_distribution_shift_filter <- ncol(mat_processed)
}

n_features_after_imputation_error_filter <- ncol(mat_filtered2)
n_features_final <- ncol(mat_processed)

# C. Save processed feature table

if (any(is.na(mat_processed))) {
  stop(
    "NA values found in processed feature matrix. ",
    "Set impute=TRUE or exclude_missing_samples=TRUE."
  )
}

df <- mat_processed %>%
  as.data.frame() %>%
  rownames_to_column(var = sample_key)

write_csv(
  df, file.path(out_dir, "bulk_x_features.csv")
)

# D. Save one-row preprocess run summary for downstream workflow decisions and
# quick retention checks across preprocessing/filtering stages.

summary_df <- build_preprocess_summary(
  n_features_orig = n_features_orig,
  n_features_after_preprocess = n_features_after_preprocess,
  n_features_after_imputation_error_filter = n_features_after_imputation_error_filter,
  n_features_after_distribution_shift_filter = n_features_after_distribution_shift_filter,
  n_features_final = n_features_final,
  had_missing_values_before_imputation = had_missing_values_before_imputation,
  imputation_requested = imputation_requested,
  imputation_applied = imputation_applied,
  run_imputation_benchmark = run_imputation_benchmark,
  used_benchmark_for_feature_filtering = used_benchmark_for_feature_filtering,
  imputation_method = imputation_method_used,
  distribution_shift_filter_enabled = distribution_shift_filter_enabled,
  distribution_shift_filter_applied = distribution_shift_filter_applied,
  distribution_shift_group_key = distribution_shift_group_key,
  distribution_shift_grouped = distribution_shift_grouped
)

write_csv(
  summary_df,
  file.path(out_dir, "preprocess_summary.csv")
)

# rm(list=ls()); gc()
