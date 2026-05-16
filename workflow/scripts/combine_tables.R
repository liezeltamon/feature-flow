#!/usr/bin/env Rscript

renv::load()
setwd(here::here())
suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(dplyr)
})
source("code/feature_module_reduction/io_helpers.R")

fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

preview_values <- function(x, n = 5) {
  paste(utils::head(unique(x), n), collapse = ", ")
}

require_column <- function(df, column_name, table_label, input_id) {
  if (!(column_name %in% colnames(df))) {
    fail(
      "Column '", column_name, "' not found in ", table_label,
      " for input '", input_id, "'."
    )
  }
}

assert_no_duplicates <- function(values, column_name, table_label, input_id) {
  if (anyDuplicated(values) > 0) {
    duplicated_values <- values[duplicated(values)]
    fail(
      "Duplicate values found in column '", column_name, "' for ",
      table_label, " in input '", input_id, "': ",
      preview_values(duplicated_values), "."
    )
  }
}

align_to_sample_ids <- function(df, target_ids, sample_key, table_label, input_id) {
  require_column(df, sample_key, table_label, input_id)
  assert_no_duplicates(df[[sample_key]], sample_key, table_label, input_id)

  current_ids <- df[[sample_key]]
  missing_ids <- setdiff(target_ids, current_ids)
  extra_ids <- setdiff(current_ids, target_ids)

  if (length(missing_ids) > 0 || length(extra_ids) > 0) {
    fail(
      table_label, " in input '", input_id,
      "' does not match the expected sample set. Missing: ",
      if (length(missing_ids) == 0) "none" else preview_values(missing_ids),
      ". Extra: ",
      if (length(extra_ids) == 0) "none" else preview_values(extra_ids),
      "."
    )
  }

  aligned_df <- df[match(target_ids, current_ids), , drop = FALSE]
  if (!identical(aligned_df[[sample_key]], target_ids)) {
    fail(
      "Failed to align ", table_label, " for input '", input_id,
      "' by sample_key '", sample_key, "'."
    )
  }

  aligned_df
}

validate_within_input <- function(feature_df, metadata_df, sample_key, input_id, mode) {
  require_column(feature_df, sample_key, "feature table", input_id)
  require_column(metadata_df, sample_key, "metadata table", input_id)
  assert_no_duplicates(feature_df[[sample_key]], sample_key, "feature table", input_id)
  assert_no_duplicates(metadata_df[[sample_key]], sample_key, "metadata table", input_id)

  if (mode == "row_bind") {
    if (!identical(feature_df[[sample_key]], metadata_df[[sample_key]])) {
      fail(
        "Feature and metadata sample order do not match for input '", input_id,
        "' in row_bind mode."
      )
    }
    return(list(feature_df = feature_df, metadata_df = metadata_df))
  }

  metadata_df <- align_to_sample_ids(
    metadata_df,
    feature_df[[sample_key]],
    sample_key,
    "metadata table",
    input_id
  )

  list(feature_df = feature_df, metadata_df = metadata_df)
}

parser <- ArgumentParser(
  description = "Combine feature and metadata tables across input subdirectories"
)
parser$add_argument("--mode", type = "character", default = "row_bind")
parser$add_argument("--input_dir", type = "character", required = TRUE)
parser$add_argument("--feature_filename", type = "character", default = "feature_table.csv")
parser$add_argument("--metadata_filename", type = "character", default = "sample_metadata.csv")
parser$add_argument("--sample_key", type = "character", required = TRUE)
parser$add_argument("--out_dir", type = "character", required = TRUE)
args <- parser$parse_args()

mode <- args$mode
input_dir <- args$input_dir
feature_filename <- args$feature_filename
metadata_filename <- args$metadata_filename
sample_key <- args$sample_key
out_dir <- args$out_dir

if (!(mode %in% c("row_bind", "column_bind"))) {
  fail("Unsupported mode: ", mode)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
input_dirs <- input_dirs[
  file.exists(file.path(input_dirs, feature_filename)) &
    file.exists(file.path(input_dirs, metadata_filename))
]

if (length(input_dirs) == 0) {
  fail("No valid input subdirectories found in: ", input_dir)
}

input_ids <- basename(input_dirs)
feature_df_list <- list()
metadata_df_list <- list()
manifest_df_list <- list()

for (i in seq_along(input_dirs)) {
  input_id <- input_ids[[i]]
  input_dir_i <- input_dirs[[i]]

  feature_info <- read_feature_table_generic(
    file.path(input_dir_i, feature_filename),
    sample_key = sample_key,
    table_label = paste0("feature table for input '", input_id, "'")
  )
  feature_df <- feature_info$feature_df
  metadata_df <- read_csv(
    file.path(input_dir_i, metadata_filename),
    show_col_types = FALSE
  )

  validated <- validate_within_input(
    feature_df = feature_df,
    metadata_df = metadata_df,
    sample_key = sample_key,
    input_id = input_id,
    mode = mode
  )

  feature_df <- validated$feature_df
  metadata_df <- validated$metadata_df

  feature_df_list[[input_id]] <- feature_df
  metadata_df_list[[input_id]] <- metadata_df
  manifest_df_list[[input_id]] <- tibble(
    input_id = input_id,
    n_samples = nrow(feature_df),
    n_feature_columns = length(feature_info$numeric_cols),
    n_non_numeric_feature_columns = length(feature_info$ignored_cols),
    combine_mode = mode
  )
}

if (mode == "row_bind") {
  feature_colnames <- lapply(feature_df_list, colnames)
  metadata_colnames <- lapply(metadata_df_list, colnames)

  if (!all(vapply(
    feature_colnames,
    identical,
    logical(1),
    y = feature_colnames[[1]]
  ))) {
    fail("Feature table columns do not match across inputs in row_bind mode.")
  }

  if (!all(vapply(
    metadata_colnames,
    identical,
    logical(1),
    y = metadata_colnames[[1]]
  ))) {
    fail("Metadata table columns do not match across inputs in row_bind mode.")
  }

  feature_df <- bind_rows(feature_df_list)
  metadata_df <- bind_rows(metadata_df_list)
  manifest_df <- bind_rows(manifest_df_list)

  if (anyDuplicated(feature_df[[sample_key]]) > 0) {
    fail("Combined feature table contains duplicated sample IDs in row_bind mode.")
  }
  if (anyDuplicated(metadata_df[[sample_key]]) > 0) {
    fail("Combined metadata table contains duplicated sample IDs in row_bind mode.")
  }
  if (!identical(feature_df[[sample_key]], metadata_df[[sample_key]])) {
    fail("Combined feature and metadata sample order do not match in row_bind mode.")
  }

  write_csv(feature_df, file.path(out_dir, "feature_table.csv"))
  write_csv(metadata_df, file.path(out_dir, "sample_metadata.csv"))
  write_csv(manifest_df, file.path(out_dir, "input_manifest.csv"))
  quit(save = "no")
}

base_input_id <- input_ids[[1]]
base_sample_ids <- feature_df_list[[base_input_id]][[sample_key]]
base_metadata_df <- align_to_sample_ids(
  metadata_df_list[[base_input_id]],
  base_sample_ids,
  sample_key,
  "metadata table",
  base_input_id
)

combined_numeric_parts <- list()
seen_numeric_cols <- character()

for (input_id in input_ids) {
  feature_df <- align_to_sample_ids(
    feature_df_list[[input_id]],
    base_sample_ids,
    sample_key,
    "feature table",
    input_id
  )
  metadata_df <- align_to_sample_ids(
    metadata_df_list[[input_id]],
    base_sample_ids,
    sample_key,
    "metadata table",
    input_id
  )

  if (!identical(colnames(metadata_df), colnames(base_metadata_df)) ||
      !identical(metadata_df, base_metadata_df)) {
    fail(
      "Metadata table for input '", input_id,
      "' does not match the reference metadata in column_bind mode."
    )
  }

  numeric_cols <- setdiff(colnames(feature_df), sample_key)
  duplicate_numeric_cols <- intersect(seen_numeric_cols, numeric_cols)
  if (length(duplicate_numeric_cols) > 0) {
    fail(
      "Duplicate numeric feature columns found across inputs in column_bind mode: ",
      preview_values(duplicate_numeric_cols), "."
    )
  }

  combined_numeric_parts[[input_id]] <- feature_df[, numeric_cols, drop = FALSE]
  seen_numeric_cols <- c(seen_numeric_cols, numeric_cols)
}

combined_feature_df <- tibble(!!sample_key := base_sample_ids)
for (input_id in input_ids) {
  numeric_part <- combined_numeric_parts[[input_id]]
  if (ncol(numeric_part) == 0) {
    next
  }
  combined_feature_df <- bind_cols(combined_feature_df, numeric_part)
}

if (!(sample_key %in% colnames(combined_feature_df))) {
  fail("Combined feature table is missing sample_key column '", sample_key, "'.")
}
if (anyDuplicated(combined_feature_df[[sample_key]]) > 0) {
  fail("Combined feature table contains duplicated sample IDs in column_bind mode.")
}
if (!identical(combined_feature_df[[sample_key]], base_metadata_df[[sample_key]])) {
  fail("Combined feature and metadata sample order do not match in column_bind mode.")
}

manifest_df <- bind_rows(manifest_df_list)

write_csv(combined_feature_df, file.path(out_dir, "feature_table.csv"))
write_csv(base_metadata_df, file.path(out_dir, "sample_metadata.csv"))
write_csv(manifest_df, file.path(out_dir, "input_manifest.csv"))
