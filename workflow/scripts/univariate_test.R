# %% Compare all groups per feature fitting 1 model
# sbatch -J univariate_test --cpus-per-task 1 --mem 50G --output %x.log.out --error %x.log.err --wrap "Rscript univariate_test.R --feature_id proportions"

renv::load()
setwd(here::here())
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(lme4)
  library(emmeans)
  source(file.path("utils", "fit_lmm.R"))
})
source("code/feature_module_reduction/io_helpers.R")
source("code/feature_module_reduction/helpers.R")

parse_additional_contrasts <- function(additional_contrasts_json) {
  if (is.null(additional_contrasts_json) || additional_contrasts_json %in% c("", "{}")) {
    return(NULL)
  }

  additional_contrasts <- jsonlite::fromJSON(
    additional_contrasts_json,
    simplifyVector = FALSE
  )
  if (length(additional_contrasts) == 0) {
    return(NULL)
  }

  lapply(names(additional_contrasts), function(contrast_term) {
    contrast_list <- additional_contrasts[[contrast_term]]
    if (!is.list(contrast_list) || is.null(names(contrast_list))) {
      stop("additional_contrasts_json must be nested by contrast term.")
    }

    lapply(names(contrast_list), function(contrast_name) {
      weights <- contrast_list[[contrast_name]]
      if (!is.list(weights) || is.null(names(weights)) || any(names(weights) == "")) {
        stop("additional contrast '", contrast_name, "' must be a named group-to-weight map.")
      }
      weights <- suppressWarnings(as.numeric(unlist(weights, use.names = FALSE)))
      if (any(is.na(weights))) {
        stop("additional contrast '", contrast_name, "' contains non-numeric weights.")
      }
      setNames(weights, names(contrast_list[[contrast_name]]))
    }) %>%
      setNames(names(contrast_list))
  }) %>%
    setNames(names(additional_contrasts))
}

get_additional_contrast_ids <- function(additional_contrasts, contrast_terms) {
  if (is.null(additional_contrasts)) {
    return(character())
  }

  contrast_terms <- intersect(contrast_terms, names(additional_contrasts))
  unique(unlist(lapply(additional_contrasts[contrast_terms], names), use.names = FALSE))
}

parse_equivalence_deltas <- function(equivalence_deltas) {
  if (is.null(equivalence_deltas) || equivalence_deltas %in% c("", "none", "NULL")) {
    return(numeric())
  }
  deltas <- strsplit(equivalence_deltas, ",", fixed = TRUE)[[1]]
  deltas <- trimws(deltas)
  deltas <- deltas[deltas != ""]
  deltas <- suppressWarnings(as.numeric(deltas))
  if (any(is.na(deltas)) || any(deltas <= 0)) {
    stop("equivalence_deltas must be comma-separated positive numeric values.")
  }
  unique(deltas)
}

format_equivalence_delta_label <- function(delta) {
  paste0(
    "delta_",
    gsub(
      "\\.",
      "_",
      format(delta, scientific = FALSE, trim = TRUE)
    )
  )
}

# %% Parameters

# %% Parameters

arg_parser <- ArgumentParser(description = "Fit LMM for each feature and collect results for heatmap plotting")
arg_parser$add_argument("--feature_path", type = "character", required = TRUE, 
                        help = "Full path to bulk_x_features.csv file")
arg_parser$add_argument("--sample_metadata_path", type = "character", required = TRUE,
                        help = "Full path to sample metadata table")
arg_parser$add_argument("--out_dir", type = "character", required = TRUE,
                        help = "Output directory for results")
arg_parser$add_argument("--sample_key", type = "character", required = TRUE)
arg_parser$add_argument("--contrast", type = "character", nargs = "+", required = TRUE)
arg_parser$add_argument("--group_key_levels", type = "character", nargs = "+", required = TRUE,
                        help = "Ordered levels for group_key factor, with 1st element as reference (e.g. healthy nonprogb progb progf active inactive pss ra)")
arg_parser$add_argument("--random_effects", type = "character", nargs = "+", default = NULL)
arg_parser$add_argument("--fixed_effects", type = "character", nargs = "+", default = NULL)
arg_parser$add_argument("--method", type = "character", default = "lmer")
arg_parser$add_argument("--additional_contrasts_json", type = "character", default = "{}")
arg_parser$add_argument("--equivalence_deltas", type = "character", default = "",
                        help = "Comma-separated positive deltas for emmeans equivalence tests.")
arg_parser$add_argument("--compute_model_assumptions", type = "character", default = "FALSE")
arg_parser$add_argument("--missing_response_policy", type = "character", default = "fail",
                        help = "How to handle missing feature values: fail or drop_feature_sample")
arg_parser$add_argument("--min_samples_per_group", type = "integer", default = 3,
                        help = "Minimum samples required per group after feature-specific filtering")
args <- arg_parser$parse_args()
list2env(args, envir = environment())

group_key = contrast[1]
#group_key_levels = rev(c("healthy", "nonprogb", "progb", "progf", "active", "inactive", "pss", "ra"))
contrast_terms = contrast
compute_model_assumptions <- parse_bool(compute_model_assumptions, "compute_model_assumptions")
additional_contrasts <- parse_additional_contrasts(additional_contrasts_json)
equivalence_deltas <- parse_equivalence_deltas(equivalence_deltas)
equivalence_delta_labels <- vapply(
  equivalence_deltas,
  format_equivalence_delta_label,
  character(1)
)

allowed_missing_response_policies <- c("fail", "drop_feature_sample")
if (!(missing_response_policy %in% allowed_missing_response_policies)) {
  stop(
    "missing_response_policy must be one of: ",
    paste(allowed_missing_response_policies, collapse = ", ")
  )
}
if (min_samples_per_group < 1) {
  stop("min_samples_per_group must be at least 1")
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# %% ----- MAIN -----

group_key_levels <- rev(group_key_levels) # Reverse to have 1st element as reference for contrasts

if (!is.null(additional_contrasts) && group_key %in% names(additional_contrasts)) {
  for (contrast_name in names(additional_contrasts[[group_key]])) {
    weights <- additional_contrasts[[group_key]][[contrast_name]]
    weight_names <- names(weights)
    if (length(weights) != length(group_key_levels)) {
      stop(
        "additional_contrasts[['", group_key, "']][['", contrast_name,
        "']] has length ", length(weights),
        " but group_key_levels has ", length(group_key_levels), " levels."
      )
    }
    if (is.null(weight_names) || any(weight_names == "")) {
      stop(
        "additional_contrasts[['", group_key, "']][['", contrast_name,
        "']] must name every group level."
      )
    }
    missing_levels <- setdiff(group_key_levels, weight_names)
    extra_levels <- setdiff(weight_names, group_key_levels)
    if (length(missing_levels) > 0 || length(extra_levels) > 0) {
      stop(
        "Names for additional_contrasts[['", group_key, "']][['",
        contrast_name, "']] must match group_key_levels. Missing: ",
        paste(missing_levels, collapse = ", "), ". Extra: ",
        paste(extra_levels, collapse = ", "), "."
      )
    }
  }
}

# %% Load data

required_metadata_cols <- unique(c(group_key, random_effects, fixed_effects))
required_metadata_cols <- required_metadata_cols[!is.null(required_metadata_cols)]
required_metadata_cols <- required_metadata_cols[required_metadata_cols != ""]

feature_info <- read_feature_table_generic(feature_path, sample_key = sample_key)
feature_df_raw <- feature_info$feature_df
meta_df <- read_metadata_aligned(
  sample_metadata_path,
  sample_key = sample_key,
  sample_ids = feature_df_raw[[sample_key]],
  required_cols = required_metadata_cols,
  preserve_cols = random_effects
)

keep_idx <- meta_df[[group_key]] %in% group_key_levels
feature_df_raw <- feature_df_raw[keep_idx, , drop = FALSE]
meta_df <- meta_df[keep_idx, , drop = FALSE]

feature_df <- feature_df_raw %>%
  column_to_rownames(sample_key) %>%
  as.data.frame()
features <- colnames(feature_df)
if (length(features) == 0) {
  stop("No numeric feature columns found in processed feature table.")
}

meta_df <- meta_df %>%
  mutate(
    across(all_of(group_key), ~ factor(., levels = group_key_levels))
  ) %>%
  column_to_rownames(sample_key)
stopifnot(identical(rownames(feature_df), rownames(meta_df)))

response_col <- ".feature_response"

prepare_feature_dataset <- function(feature) {
  feature_values <- feature_df[[feature]]
  keep <- rep(TRUE, length(feature_values))
  n_samples_dropped_na <- 0L

  if (missing_response_policy == "drop_feature_sample") {
    feature_is_missing <- is.na(feature_values)
    n_samples_dropped_na <- sum(feature_is_missing)
    keep <- keep & !feature_is_missing

    if (length(required_metadata_cols) > 0) {
      metadata_is_complete <- stats::complete.cases(
        meta_df[, required_metadata_cols, drop = FALSE]
      )
      keep <- keep & metadata_is_complete
    }
  }

  feature_meta_df <- meta_df[keep, , drop = FALSE]
  feature_values <- feature_values[keep]
  n_samples_used <- length(feature_values)

  if (missing_response_policy == "drop_feature_sample") {
    group_counts <- table(
      factor(feature_meta_df[[group_key]], levels = group_key_levels)
    )
    low_groups <- names(group_counts)[group_counts < min_samples_per_group]
    if (length(low_groups) > 0) {
      return(list(
        dataset = NULL,
        n_samples_used = n_samples_used,
        n_samples_dropped_na = n_samples_dropped_na,
        error_message = paste0(
          "Fewer than ", min_samples_per_group,
          " samples after feature-specific filtering for group(s): ",
          paste(low_groups, collapse = ", ")
        )
      ))
    }
  }

  feature_sd <- stats::sd(feature_values, na.rm = TRUE)
  if (all(is.na(feature_values)) || !is.finite(feature_sd) || feature_sd == 0) {
    return(list(
      dataset = NULL,
      n_samples_used = n_samples_used,
      n_samples_dropped_na = n_samples_dropped_na,
      error_message = "Feature has zero or undefined variance after filtering."
    ))
  }

  feature_dataset <- feature_meta_df
  feature_dataset[[response_col]] <- as.numeric(scale(feature_values))

  list(
    dataset = feature_dataset,
    n_samples_used = n_samples_used,
    n_samples_dropped_na = n_samples_dropped_na,
    error_message = NA_character_
  )
}
 
# %% Fit models for each feature

out <- setNames(lapply(features, function(x) {
  prepared <- prepare_feature_dataset(x)
  if (!is.na(prepared$error_message)) {
    return(list(
      fit_output = NULL,
      n_samples_used = prepared$n_samples_used,
      n_samples_dropped_na = prepared$n_samples_dropped_na,
      error_message = prepared$error_message
    ))
  }

  tryCatch({
    list(
      fit_output = fit_lmm(
        prepared$dataset,
        response = response_col,
        contrast = contrast_terms,
        random_effects = random_effects,
        fixed_effects = fixed_effects,
        verbose = FALSE,
        save_models = FALSE,
        save_model_dir = ".",
        method = method,
        additional_contrasts = additional_contrasts,
        equivalence_deltas = equivalence_deltas
      ),
      n_samples_used = prepared$n_samples_used,
      n_samples_dropped_na = prepared$n_samples_dropped_na,
      error_message = NA_character_
    )
  }, error = function(e) {
    list(
      fit_output = NULL,
      n_samples_used = prepared$n_samples_used,
      n_samples_dropped_na = prepared$n_samples_dropped_na,
      error_message = conditionMessage(e)
    )
  })
}), features)

failed_features_df <- tibble(
  feature = features,
  method = method,
  error_message = sapply(out, function(x) x$error_message)
) %>%
  filter(!is.na(error_message))

# %% Collect results for heatmap plotting in univariate.R

# Get all possible pairs from group_key_levels
group_pairs <- t(combn(group_key_levels, 2)) %>%
  as.data.frame() %>%
  unite("contrast", V1, V2, sep = "_vs_") %>%
  pull(contrast)
additional_contrast_ids <- get_additional_contrast_ids(additional_contrasts, contrast_terms)
result_rows <- c(group_pairs, additional_contrast_ids)

# Initialise matrices to store estimates and p-values for heatmap plotting
mx <- matrix(NA_real_, ncol = length(features), nrow = length(result_rows))
rownames(mx) <- result_rows
colnames(mx) <- features
result <- list(effect_size = mx, p_value = mx)
for (delta_label in equivalence_delta_labels) {
  result[[paste0("equiv_p_value_", delta_label)]] <- mx
  result[[paste0("equiv_t_ratio_", delta_label)]] <- mx
}

collect_specs <- c(effect_size = "cohens.d", p_value = "coef.pvalue")
for (delta_label in equivalence_delta_labels) {
  collect_specs[paste0("equiv_p_value_", delta_label)] <- paste0(
    "equiv.pvalue.",
    delta_label
  )
  collect_specs[paste0("equiv_t_ratio_", delta_label)] <- paste0(
    "equiv.t.ratio.",
    delta_label
  )
}

for (i in seq_along(features)) {

  feature <- features[i]

  fit_output <- out[[feature]]$fit_output
  if (is.null(fit_output)) {
    next
  }
  model_df <- fit_output$model_df

  # Collect
  for (x in names(collect_specs)) {
    statistic_name <- collect_specs[[x]]
    mx <- result[[x]]
    collect_df <- model_df %>%
      filter(statistic == statistic_name) %>%
      mutate(
        comparison_id = if_else(
          is.na(contrast),
          var,
          paste(var, contrast, sep = "_vs_")
        )
      ) %>%
      filter(comparison_id %in% rownames(mx)) %>%
      distinct(comparison_id, .keep_all = TRUE)
    if (nrow(collect_df) > 0) {
      result[[x]][collect_df$comparison_id, feature] <- collect_df$value
    }
  }

}

# %% Get global p-value for each feature from out

global_p_value <- sapply(features, function(feature) {
  fit_output <- out[[feature]]$fit_output
  if (is.null(fit_output)) {
    return(NA_real_)
  }
  model_df <- fit_output$model_df
  # Get p-value for group_key from model_terms_test
  p_value <- model_df %>%
    filter(statistic == paste0(contrast_terms[1], ".Pr(>F)")) %>%
    pull(value)
  if (length(p_value) == 0) {
    return(NA_real_)
  }
  return(p_value[[1]])
})

# Add global p-value to result list
result$global_p_value <- global_p_value
names(result$global_p_value) <- features

# FDR adjust global p-values
result$global_p_value_adj <- safe_p_adjust(result$global_p_value, method = "BH")

# %% Save per-feature model assumption summaries

shapiro_pvalue <- setNames(rep(NA_real_, length(features)), features)
shapiro_padj <- setNames(rep(NA_real_, length(features)), features)
levene_pvalue <- setNames(rep(NA_real_, length(features)), features)
levene_padj <- setNames(rep(NA_real_, length(features)), features)
model_assumptions_met <- setNames(rep(NA, length(features)), features)

if (compute_model_assumptions) {
  shapiro_pvalue <- sapply(features, function(feature) {
    fit_output <- out[[feature]]$fit_output
    if (is.null(fit_output)) {
      return(NA_real_)
    }
    fit_output$assumption_tests$shapiro$p.value
  })
  names(shapiro_pvalue) <- features
  shapiro_padj <- safe_p_adjust(shapiro_pvalue, method = "BH")
  names(shapiro_padj) <- features

  levene_pvalue <- sapply(features, function(feature) {
    fit_output <- out[[feature]]$fit_output
    if (is.null(fit_output) || is.null(fit_output$assumption_tests$levene)) {
      return(NA_real_)
    }
    as.numeric(fit_output$assumption_tests$levene[1, "Pr(>F)"])
  })
  names(levene_pvalue) <- features

  non_na_levene <- !is.na(levene_pvalue)
  if (any(non_na_levene)) {
    levene_padj[non_na_levene] <- p.adjust(levene_pvalue[non_na_levene], method = "BH")
  }
  names(levene_padj) <- features

  model_assumptions_met <- shapiro_padj >= 0.05 & (is.na(levene_padj) | levene_padj >= 0.05)
  names(model_assumptions_met) <- features
}

result$shapiro_pvalue <- shapiro_pvalue
result$shapiro_padj <- shapiro_padj
result$levene_pvalue <- levene_pvalue
result$levene_padj <- levene_padj
result$model_assumptions_met <- model_assumptions_met
result$n_samples_used <- setNames(
  as.integer(sapply(out, function(x) x$n_samples_used)),
  features
)
result$n_samples_dropped_na <- setNames(
  as.integer(sapply(out, function(x) x$n_samples_dropped_na)),
  features
)

# %% Save

stopifnot(identical(names(result$global_p_value), colnames(result$effect_size)))
stopifnot(identical(names(result$shapiro_pvalue), colnames(result$effect_size)))
stopifnot(identical(names(result$shapiro_padj), colnames(result$effect_size)))
stopifnot(identical(names(result$levene_pvalue), colnames(result$effect_size)))
stopifnot(identical(names(result$levene_padj), colnames(result$effect_size)))
stopifnot(identical(names(result$model_assumptions_met), colnames(result$effect_size)))
stopifnot(identical(names(result$n_samples_used), colnames(result$effect_size)))
stopifnot(identical(names(result$n_samples_dropped_na), colnames(result$effect_size)))
matrix_result_names <- names(result)[vapply(result, is.matrix, logical(1))]
for (mx_name in matrix_result_names) {
  mx <- result[[mx_name]]
  stopifnot(all(!duplicated(rownames(mx))))
  stopifnot(all(!duplicated(colnames(mx))))
}

write_csv(failed_features_df, file.path(out_dir, "failed_features.csv"))

saveRDS(
  result, file.path(out_dir, "result_matrices.rds")
)
