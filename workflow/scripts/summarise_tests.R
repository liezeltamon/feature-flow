# Summarise univariate test outputs into one wide table for downstream filtering.

# Setup

renv::load()
setwd(here::here())
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
})
source("code/feature_module_reduction/helpers.R")

# Functions

adjust_p_value_matrix <- function(p_values, p_adjust_method = "fdr", adjust_by = "column") {
  if (adjust_by == "all") {
    adjusted <- safe_p_adjust(as.vector(p_values), method = p_adjust_method)
    return(matrix(
      adjusted,
      nrow = nrow(p_values),
      ncol = ncol(p_values),
      dimnames = dimnames(p_values)
    ))
  } else if (adjust_by == "row") {
    p_value_adj <- t(apply(p_values, 1, function(row) {
      safe_p_adjust(row, method = p_adjust_method)
    }))
    dimnames(p_value_adj) <- dimnames(p_values)
    return(p_value_adj)
  } else if (adjust_by == "column") {
    p_value_adj <- matrix(
      NA_real_,
      nrow = nrow(p_values),
      ncol = ncol(p_values),
      dimnames = dimnames(p_values)
    )
    for (j in seq_len(ncol(p_values))) {
      p_value_adj[, j] <- safe_p_adjust(p_values[, j], method = p_adjust_method)
    }
    return(p_value_adj)
  } else {
    stop("adjust_by must be 'all', 'row', or 'column'")
  }
}

process_results <- function(
  association_results, p_adjust_method = "fdr", adjust_by = "column"
) {
  association_results$p_value_adj <- adjust_p_value_matrix(
    association_results$p_value,
    p_adjust_method = p_adjust_method,
    adjust_by = adjust_by
  )

  equiv_pvalue_names <- grep(
    "^equiv_p_value_",
    names(association_results),
    value = TRUE
  )
  for (pvalue_name in equiv_pvalue_names) {
    delta_label <- sub("^equiv_p_value_", "", pvalue_name)
    association_results[[paste0("equiv_p_adj_", delta_label)]] <- adjust_p_value_matrix(
      association_results[[pvalue_name]],
      p_adjust_method = p_adjust_method,
      adjust_by = adjust_by
    )
  }

  association_results
}

build_wide_block <- function(mx, value_prefix) {
  as.data.frame(mx) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("feature") %>%
    setNames(c("feature", paste0(value_prefix, "...", rownames(mx))))
}

get_method_summary <- function(result_matrices, method_name) {
  result_matrices <- process_results(
    result_matrices,
    p_adjust_method = "fdr",
    adjust_by = "column"
  )

  effect_size_df <- build_wide_block(
    result_matrices$effect_size,
    paste0("effect_size_", method_name)
  )
  pvalue_df <- build_wide_block(
    result_matrices$p_value,
    paste0("pvalue_", method_name)
  )
  padj_df <- build_wide_block(
    result_matrices$p_value_adj,
    paste0("padj_", method_name)
  )

  equiv_delta_labels <- sub(
    "^equiv_p_value_",
    "",
    grep("^equiv_p_value_", names(result_matrices), value = TRUE)
  )
  equivalence_dfs <- unlist(
    lapply(equiv_delta_labels, function(delta_label) {
      list(
        build_wide_block(
          result_matrices[[paste0("equiv_p_value_", delta_label)]],
          paste0("equiv_pvalue_", delta_label, "_", method_name)
        ),
        build_wide_block(
          result_matrices[[paste0("equiv_p_adj_", delta_label)]],
          paste0("equiv_padj_", delta_label, "_", method_name)
        ),
        build_wide_block(
          result_matrices[[paste0("equiv_t_ratio_", delta_label)]],
          paste0("equiv_t_ratio_", delta_label, "_", method_name)
        )
      )
    }),
    recursive = FALSE
  )

  global_df <- tibble(
    feature = names(result_matrices$global_p_value),
    !!paste0("global_pvalue_", method_name) := as.numeric(result_matrices$global_p_value),
    !!paste0("global_padj_", method_name) := as.numeric(result_matrices$global_p_value_adj)
  )

  sample_count_df <- tibble(feature = colnames(result_matrices$effect_size))
  if ("n_samples_used" %in% names(result_matrices)) {
    sample_count_df[[paste0("n_samples_used_", method_name)]] <- as.integer(
      result_matrices$n_samples_used[sample_count_df$feature]
    )
  }
  if ("n_samples_dropped_na" %in% names(result_matrices)) {
    sample_count_df[[paste0("n_samples_dropped_na_", method_name)]] <- as.integer(
      result_matrices$n_samples_dropped_na[sample_count_df$feature]
    )
  }

  list(
    method_df = reduce(
      c(
        list(effect_size_df, pvalue_df, padj_df),
        equivalence_dfs,
        list(global_df, sample_count_df)
      ),
      full_join,
      by = "feature"
    ),
    result_matrices = result_matrices
  )
}

status_row <- function(
  source_id,
  status,
  reason,
  source_paths,
  n_features = NA_integer_
) {
  tibble(
    source = source_id,
    status = status,
    reason = reason,
    methods_found = paste(sort(unique(source_paths$method)), collapse = ","),
    result_paths = paste(sort(unique(source_paths$path)), collapse = ";"),
    n_features = as.integer(n_features)
  )
}

empty_status_df <- function() {
  tibble(
    source = character(),
    status = character(),
    reason = character(),
    methods_found = character(),
    result_paths = character(),
    n_features = integer()
  )
}

empty_summary_df <- function() {
  tibble(source = character(), feature = character())
}

read_normality_results <- function(module_normality_dir, available_only = FALSE) {
  module_normality_paths <- list.files(
    module_normality_dir,
    pattern = "^normality_test_results\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(module_normality_paths) == 0) {
    return(tibble(source = character(), feature = character(), is_normal = logical()))
  }

  normality_root <- normalizePath(module_normality_dir, mustWork = FALSE)
  normality_prefix <- paste0(normality_root, .Platform$file.sep)

  normality_df_list <- lapply(module_normality_paths, function(path) {
    tryCatch({
      path_norm <- normalizePath(path, mustWork = FALSE)
      if (startsWith(path_norm, normality_prefix)) {
        path_rel <- substring(path_norm, nchar(normality_prefix) + 1)
        dataset_id <- strsplit(path_rel, .Platform$file.sep, fixed = TRUE)[[1]][[1]]
      } else {
        path_parts <- strsplit(path, .Platform$file.sep, fixed = TRUE)[[1]]
        dataset_id <- path_parts[length(path_parts) - 2]
      }

      read_csv(path, show_col_types = FALSE) %>%
        transmute(
          source = dataset_id,
          feature,
          is_normal = case_when(
            is_normal == 1 ~ TRUE,
            is_normal == 0 ~ FALSE,
            TRUE ~ NA
          )
        )
    }, error = function(e) {
      if (!available_only) {
        stop(e)
      }
      warning(
        "Skipping unreadable normality file: ",
        path,
        " (",
        conditionMessage(e),
        ")"
      )
      NULL
    })
  })

  normality_df_list <- normality_df_list[
    !vapply(normality_df_list, is.null, logical(1))
  ]

  if (length(normality_df_list) == 0) {
    return(tibble(source = character(), feature = character(), is_normal = logical()))
  }

  bind_rows(normality_df_list)
}

build_source_summary <- function(
  source_id,
  path_df,
  required_methods,
  module_normality_df,
  out_dir,
  parametric_method,
  nonparametric_method,
  available_only = FALSE
) {
  all_source_paths <- path_df %>%
    filter(source == source_id)
  source_paths <- all_source_paths %>%
    filter(method %in% required_methods)

  missing_methods <- setdiff(required_methods, unique(source_paths$method))
  if (length(missing_methods) > 0) {
    reason <- paste("missing required method(s):", paste(missing_methods, collapse = ", "))
    if (available_only) {
      warning("Skipping source ", source_id, ": ", reason)
      return(list(
        summary = NULL,
        status = status_row(source_id, "skipped", reason, all_source_paths)
      ))
    }
    stop("Source ", source_id, " does not contain both required methods.")
  }

  make_summary <- function() {
    method_summaries <- setNames(
      lapply(required_methods, function(method_name) {
        method_path <- sort(source_paths$path[source_paths$method == method_name])[[1]]
        result_matrices <- readRDS(method_path)
        get_method_summary(result_matrices, method_name)
      }),
      required_methods
    )

    parametric_results <- method_summaries[[parametric_method]]$result_matrices
    nonparametric_results <- method_summaries[[nonparametric_method]]$result_matrices

    if (!identical(
      colnames(parametric_results$effect_size),
      colnames(nonparametric_results$effect_size)
    )) {
      stop("Parametric and nonparametric result matrices have different features.")
    }
    if (!identical(
      rownames(parametric_results$effect_size),
      rownames(nonparametric_results$effect_size)
    )) {
      stop("Parametric and nonparametric result matrices have different comparisons.")
    }

    assumption_df <- tibble(
      feature = names(parametric_results$model_assumptions_met),
      model_assumptions_met = as.logical(parametric_results$model_assumptions_met),
      shapiro_pvalue = as.numeric(parametric_results$shapiro_pvalue),
      shapiro_padj = as.numeric(parametric_results$shapiro_padj),
      levene_pvalue = as.numeric(parametric_results$levene_pvalue),
      levene_padj = as.numeric(parametric_results$levene_padj)
    )

    combined_df <- reduce(
      c(
        list(assumption_df),
        lapply(method_summaries, `[[`, "method_df")
      ),
      full_join,
      by = "feature"
    ) %>%
      mutate(source = source_id, .before = 1)

    comparisons <- rownames(parametric_results$effect_size)
    equiv_delta_labels <- sub(
      "^equiv_p_value_",
      "",
      grep("^equiv_p_value_", names(parametric_results), value = TRUE)
    )
    for (comparison in comparisons) {
      param_effect_col <- paste0("effect_size_", parametric_method, "...", comparison)
      param_pvalue_col <- paste0("pvalue_", parametric_method, "...", comparison)
      param_padj_col <- paste0("padj_", parametric_method, "...", comparison)

      nonparam_effect_col <- paste0("effect_size_", nonparametric_method, "...", comparison)
      nonparam_pvalue_col <- paste0("pvalue_", nonparametric_method, "...", comparison)
      nonparam_padj_col <- paste0("padj_", nonparametric_method, "...", comparison)

      combined_df[[paste0("effect_size_final...", comparison)]] <- ifelse(
        combined_df$model_assumptions_met,
        combined_df[[param_effect_col]],
        combined_df[[nonparam_effect_col]]
      )
      combined_df[[paste0("pvalue_final...", comparison)]] <- ifelse(
        combined_df$model_assumptions_met,
        combined_df[[param_pvalue_col]],
        combined_df[[nonparam_pvalue_col]]
      )
      combined_df[[paste0("padj_final...", comparison)]] <- ifelse(
        combined_df$model_assumptions_met,
        combined_df[[param_padj_col]],
        combined_df[[nonparam_padj_col]]
      )

      for (delta_label in equiv_delta_labels) {
        param_equiv_pvalue_col <- paste0(
          "equiv_pvalue_", delta_label, "_", parametric_method, "...", comparison
        )
        param_equiv_padj_col <- paste0(
          "equiv_padj_", delta_label, "_", parametric_method, "...", comparison
        )
        param_equiv_t_ratio_col <- paste0(
          "equiv_t_ratio_", delta_label, "_", parametric_method, "...", comparison
        )
        nonparam_equiv_pvalue_col <- paste0(
          "equiv_pvalue_", delta_label, "_", nonparametric_method, "...", comparison
        )
        nonparam_equiv_padj_col <- paste0(
          "equiv_padj_", delta_label, "_", nonparametric_method, "...", comparison
        )
        nonparam_equiv_t_ratio_col <- paste0(
          "equiv_t_ratio_", delta_label, "_", nonparametric_method, "...", comparison
        )

        required_equiv_cols <- c(
          param_equiv_pvalue_col,
          param_equiv_padj_col,
          param_equiv_t_ratio_col,
          nonparam_equiv_pvalue_col,
          nonparam_equiv_padj_col,
          nonparam_equiv_t_ratio_col
        )
        if (all(required_equiv_cols %in% colnames(combined_df))) {
          combined_df[[paste0("equiv_pvalue_", delta_label, "_final...", comparison)]] <- ifelse(
            combined_df$model_assumptions_met,
            combined_df[[param_equiv_pvalue_col]],
            combined_df[[nonparam_equiv_pvalue_col]]
          )
          combined_df[[paste0("equiv_padj_", delta_label, "_final...", comparison)]] <- ifelse(
            combined_df$model_assumptions_met,
            combined_df[[param_equiv_padj_col]],
            combined_df[[nonparam_equiv_padj_col]]
          )
          combined_df[[paste0("equiv_t_ratio_", delta_label, "_final...", comparison)]] <- ifelse(
            combined_df$model_assumptions_met,
            combined_df[[param_equiv_t_ratio_col]],
            combined_df[[nonparam_equiv_t_ratio_col]]
          )
        }
      }
    }

    param_global_p_col <- paste0("global_pvalue_", parametric_method)
    param_global_padj_col <- paste0("global_padj_", parametric_method)
    nonparam_global_p_col <- paste0("global_pvalue_", nonparametric_method)
    nonparam_global_padj_col <- paste0("global_padj_", nonparametric_method)

    combined_df$global_pvalue_final <- ifelse(
      combined_df$model_assumptions_met,
      combined_df[[param_global_p_col]],
      combined_df[[nonparam_global_p_col]]
    )
    combined_df$global_padj_final <- ifelse(
      combined_df$model_assumptions_met,
      combined_df[[param_global_padj_col]],
      combined_df[[nonparam_global_padj_col]]
    )

    source_summary_df <- combined_df %>%
      left_join(
        module_normality_df %>% filter(source == source_id),
        by = c("source", "feature")
      )

    source_out_dir <- file.path(out_dir, source_id)
    dir.create(source_out_dir, showWarnings = FALSE, recursive = TRUE)
    write_csv(source_summary_df, file.path(source_out_dir, "summary_df.csv"))

    source_summary_df
  }

  if (available_only) {
    error_reason <- NULL
    source_summary_df <- tryCatch(
      make_summary(),
      error = function(e) {
        error_reason <<- paste(
          "error while reading or summarising source:",
          conditionMessage(e)
        )
        warning("Skipping source ", source_id, ": ", error_reason)
        NULL
      }
    )

    if (is.null(source_summary_df)) {
      return(list(
        summary = NULL,
        status = status_row(source_id, "skipped", error_reason, all_source_paths)
      ))
    }
  } else {
    source_summary_df <- make_summary()
  }

  list(
    summary = source_summary_df,
    status = status_row(
      source_id,
      "included",
      "",
      all_source_paths,
      nrow(source_summary_df)
    )
  )
}

# Parameters

parser <- ArgumentParser(description = "Summarise univariate tests into one wide CSV")
parser$add_argument("--src_dir", type = "character", required = TRUE)
parser$add_argument("--module_normality_dir", type = "character", required = TRUE)
parser$add_argument("--out_dir", type = "character", required = TRUE)
parser$add_argument("--parametric_method", type = "character", default = "lmer")
parser$add_argument("--nonparametric_method", type = "character", default = "rank_transform")
parser$add_argument("--available_only", type = "character", default = "false")
args <- parser$parse_args()
list2env(args, envir = environment())
available_only <- parse_bool(available_only, "available_only")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----- MAIN -----

result_matrices_paths <- list.files(
  src_dir,
  pattern = "^result_matrices\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(result_matrices_paths) == 0) {
  if (available_only) {
    warning("No result_matrices.rds files found in: ", src_dir)
    write_csv(empty_summary_df(), file.path(out_dir, "summary_df.csv"))
    write_csv(empty_status_df(), file.path(out_dir, "summarise_status.csv"))
    quit(save = "no", status = 0)
  } else {
    stop("No result_matrices.rds files found in: ", src_dir)
  }
}

path_df <- tibble(path = result_matrices_paths) %>%
  mutate(
    method = basename(dirname(path)),
    source = basename(dirname(dirname(path)))
  )

required_methods <- c(parametric_method, nonparametric_method)
missing_methods <- setdiff(required_methods, unique(path_df$method))
if (!available_only && length(missing_methods) > 0) {
  stop("Missing required methods in src_dir: ", paste(missing_methods, collapse = ", "))
}

module_normality_df <- read_normality_results(
  module_normality_dir,
  available_only = available_only
)

source_results <- lapply(sort(unique(path_df$source)), function(source_id) {
  build_source_summary(
    source_id = source_id,
    path_df = path_df,
    required_methods = required_methods,
    module_normality_df = module_normality_df,
    out_dir = out_dir,
    parametric_method = parametric_method,
    nonparametric_method = nonparametric_method,
    available_only = available_only
  )
})

summary_df_list <- lapply(source_results, `[[`, "summary")
summary_df_list <- summary_df_list[
  !vapply(summary_df_list, is.null, logical(1))
]

summary_df <- if (length(summary_df_list) == 0) {
  empty_summary_df()
} else {
  bind_rows(summary_df_list)
}
write_csv(summary_df, file.path(out_dir, "summary_df.csv"))

if (available_only) {
  status_df <- bind_rows(lapply(source_results, `[[`, "status"))
  write_csv(status_df, file.path(out_dir, "summarise_status.csv"))
}
