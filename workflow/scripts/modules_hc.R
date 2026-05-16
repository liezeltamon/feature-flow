# %% Identify modules with using hierarchical clustering via vdjremix

# %% Setup

renv::load()
setwd(here::here())
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(WGCNA)
  library(vdjremix)
  source("utils/robust_correlation_variant.R")
  source("code/feature_module_reduction/io_helpers.R")
  source("code/feature_module_reduction/helpers.R")
})

# %% Parameters

parser <- ArgumentParser(description = "Identify modules using hierarchical clustering")
parser$add_argument("--src_id", type = "character", required = TRUE)
parser$add_argument("--feature_path", type = "character", required = TRUE)
parser$add_argument("--preprocess_summary_path", type = "character", required = TRUE)
parser$add_argument("--metadata_path", type = "character", default = NULL)
parser$add_argument("--out_dir", type = "character", required = TRUE)
parser$add_argument("--n_cores", type = "integer", default = 1)
parser$add_argument("--sample_key", type = "character", default = "filter.sample_id")
parser$add_argument("--subset_key", type = "character", default = NULL)
parser$add_argument("--subset_values", type = "character", default = NULL,
                    help = "Comma-separated values to keep (e.g. healthy,nonprogb)")
parser$add_argument("--subset_name", type = "character", default = "all")
parser$add_argument("--min_cluster_size", type = "integer", default = 10)
parser$add_argument("--n_repeats", type = "integer", default = 1)
parser$add_argument("--filter_features_by_subset_unique_values", type = "character", default = "FALSE")
parser$add_argument("--subset_min_unique_values", type = "integer", default = NULL)
args <- parser$parse_args()
list2env(args, envir = environment())

filter_features_by_subset_unique_values <- parse_bool(
  filter_features_by_subset_unique_values,
  "filter_features_by_subset_unique_values"
)

# Parse subset_values if provided
if (!is.null(subset_values)) {
  subset_values <- strsplit(subset_values, ",")[[1]]
}
src_path <- feature_path
out_dir <- file.path(out_dir, src_id, subset_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Using n_cores = ", n_cores)

# %% ----- MAIN -----

feature_info <- read_feature_table_generic(src_path, sample_key = sample_key)
df <- feature_info$feature_df
preprocess_summary <- read_csv(preprocess_summary_path, show_col_types = FALSE)
stopifnot(nrow(preprocess_summary) == 1)
imputation_applied <- preprocess_summary$imputation_applied[[1]]
configured_n_repeats <- n_repeats
metadata_df <- NULL
if (!is.null(metadata_path) && nzchar(metadata_path)) {
  required_metadata_cols <- NULL
  if (!is.null(subset_key) && !is.null(subset_values)) {
    required_metadata_cols <- subset_key
  }
  metadata_df <- read_metadata_aligned(
    metadata_path,
    sample_key = sample_key,
    sample_ids = df[[sample_key]],
    required_cols = required_metadata_cols
  )
}

# %% A. Optionally subset input matrix

if (!is.null(subset_key) && !is.null(subset_values)) {
  if (is.null(metadata_df)) {
    stop("metadata_path is required when subsetting")
  }
  if (!all(subset_values %in% metadata_df[[subset_key]])) {
    stop("Not all subset_values found in metadata column '", subset_key, "'")
  }
  keep_idx <- metadata_df[[subset_key]] %in% subset_values
  df <- df[keep_idx, , drop = FALSE]
  metadata_df <- metadata_df[keep_idx, , drop = FALSE]
  message("Subsetting on ", subset_key, " in: ", paste(subset_values, collapse = ", "))
}

# %% B. Build matrix

mx <- feature_table_to_matrix(df, sample_key, table_label = "processed feature table")

if (!all(is.finite(mx))) {
  stop("Non-finite values found in processed feature matrix")
}

# Filter features by unique values within subset if requested (to avoid issues with zero-variance features in correlation calculation)

n_features_before_subset_filter <- ncol(mx)
n_features_after_subset_filter <- n_features_before_subset_filter
subset_min_unique_values_used <- NA_integer_

if (filter_features_by_subset_unique_values) {
  if (is.null(subset_min_unique_values)) {
    subset_min_unique_values <- max(6, ceiling(0.05 * nrow(mx)))
  }
  subset_min_unique_values_used <- subset_min_unique_values

  message(
    "Filtering module-learning features by subset unique values ",
    "(min_unique_values=", subset_min_unique_values_used, ")"
  )

  mx <- vdjremix::preprocess_features(
    mx,
    missingness_threshold = Inf,
    min_unique_values = subset_min_unique_values_used
  )
  n_features_after_subset_filter <- ncol(mx)

  if (n_features_after_subset_filter <= 0) {
    stop("No features remain after subset unique-value filtering")
  }
}

write_csv(
  tibble(
    filter_enabled = filter_features_by_subset_unique_values,
    n_samples_used = nrow(mx),
    min_unique_values = subset_min_unique_values_used,
    n_features_before_subset_filter = n_features_before_subset_filter,
    n_features_after_subset_filter = n_features_after_subset_filter
  ),
  file.path(out_dir, "module_feature_filter_summary.csv")
)

# %% C. Correlate features

if (!isTRUE(imputation_applied)) {
  n_repeats <- 1L
}

subsampling_enabled <- n_repeats > 1
depth <- if (subsampling_enabled) vdjremix::get_subsample_depth(mx) else nrow(mx)
message(
  subset_name,
  " - Computing correlation matrix ",
  if (subsampling_enabled) "with" else "without",
  " subsampling",
  " (configured_n_repeats=", configured_n_repeats,
  ", effective_n_repeats=", n_repeats,
  ", imputation_applied=", imputation_applied, ")"
)

t0 <- proc.time()
# Get feature correlation matrix
# cor_mat <- vdjremix::robust_correlation(
#   feature_matrix = mx,
#   subsample_depth = depth,
#   n_repeats = n_repeats,
#   n_cores = n_cores
# )
cor_mat <- robust_correlation_variant(
  feature_matrix = mx,
  subsample_depth = depth,
  n_repeats = n_repeats,
  n_cores = n_cores,
  use_crossprod = TRUE
)
message(subset_name, " - Time taken to compute correlation matrix: ", elapsed_sec(t0), " sec")

if (!all(is.finite(cor_mat))) {
  stop(subset_name, ": NA values found in feature correlation matrix")
}

# %% D. Identify modules

t0 <- proc.time()
clustering <- vdjremix::identify_modules(
  correlation_matrix = cor_mat,
  min_cluster_size = min_cluster_size
)
message(subset_name, " - Time taken to identify modules: ", elapsed_sec(t0), " sec")
saveRDS(clustering, file = file.path(out_dir, "module_members.rds"))

# %% E. Compute eigengenes

eig <- vdjremix::compute_eigengenes(
  scaled_matrix = scale(mx),
  modules = clustering$modules
)
saveRDS(eig, file = file.path(out_dir, "module_eigengenes_hc.rds"))

# Save each module to a separate file m<index>.txt # including singletons if present (1 feature per file)

out_dir_tmp <- file.path(out_dir, "modules")
dir.create(out_dir_tmp, showWarnings = FALSE)
for (i in seq_along(clustering$modules)) {
  writeLines(
    clustering$modules[[i]], file.path(out_dir_tmp, paste0("m", i, ".txt"))
  )
}

if (length(eig$singletons) > 0) {
  for (i in seq_along(eig$singletons)) {
    writeLines(
      eig$singletons[[i]], file.path(out_dir_tmp, paste0("s", i, ".txt"))
    )
  }
  message("Also saved ", length(eig$singletons), " singletons")
}

# %% F. Save module information

# Save structured per-module member tables for downstream annotation

out_dir_tmp2 <- file.path(out_dir, "module_members")
dir.create(out_dir_tmp2, showWarnings = FALSE)

for (i in seq_along(clustering$modules)) {
  module_id <- paste(subset_name, src_id, "m", i, sep = "..")
  module_dir <- file.path(out_dir_tmp2, module_id)
  dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)

  w <- eig$loadings[[i]]
  p <- length(w)
  cutoff <- sqrt(1 / p)

  module_df <- tibble(
    feature = names(w),
    loading = as.numeric(w),
    important_loading = abs(as.numeric(w)) >= cutoff
  ) %>%
    arrange(desc(abs(loading)))

  write_csv(
    module_df,
    file.path(module_dir, "members.csv")
  )
}

# Save feature table with module eigengenes

eig_mx <- eig$eigengenes
# Assert eig_mx has identical rownames with mx and no non-finite values
if (!all(rownames(eig_mx) == rownames(mx))) {
  stop("Rownames of eigengenes matrix do not match original matrix")
}
if (!all(is.finite(eig_mx))) {
  stop("Non-finite values found in eigengenes matrix")
}

as.data.frame(eig_mx) %>%
  rownames_to_column(var = sample_key) %>%
  # Rename eigengene columns to have prefix <subset_name>__<src_id>__m<current column name>
  rename_with(
    ~ paste(subset_name, src_id, "m", ., sep = ".."),
    .cols = -all_of(sample_key)
  ) %>%
  write_csv(file.path(out_dir, "bulk_x_features.csv"))

message("Done for subset: ", subset_name)
