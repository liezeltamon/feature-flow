# Identify appropriate RMSE cutoff for filtering features before imputation
# Use normalised RMSE especially when working with heterogeneous features

# Setup

renv::load()
setwd(here::here())
suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(foreach)
  library(doParallel)
  source("utils/filter_features_by_imputation_error_variant.R")
})

# Parameters

parser <- ArgumentParser(description = "Evaluate preprocessing and imputation filtering")
parser$add_argument("--feature_table_dir", type = "character", required = TRUE)
parser$add_argument("--feature_filename", type = "character", required = TRUE)
parser$add_argument("--benchmark_dir", type = "character", required = TRUE)
parser$add_argument("--error_metric", type = "character", default = "nrmse_mad")
parser$add_argument("--imputation_method", type = "character", default = NULL)
parser$add_argument("--n_jobs", type = "integer", default = 1)
parser$add_argument(
  "--thresholds",
  type = "character",
  default = "0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0"
)
args <- parser$parse_args()
list2env(args, globalenv())

thresholds <- as.numeric(strsplit(thresholds, ",")[[1]])
out_dir <- file.path(dirname(benchmark_dir), "evaluate_preprocessing", error_metric)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

DETAIL_DATASETS_PER_PAGE <- 12
OVERVIEW_LEGEND_MAX_DATASETS <- DETAIL_DATASETS_PER_PAGE

# Functions

resolve_feature_path <- function(dataset_id, feature_table_dir, feature_filename) {
  candidate_paths <- c(
    file.path(feature_table_dir, dataset_id, feature_filename),
    file.path(feature_table_dir, feature_filename)
  )
  existing_paths <- candidate_paths[file.exists(candidate_paths)]

  if (length(existing_paths) == 0) {
    stop("Feature table not found for dataset_id: ", dataset_id)
  }

  existing_paths[[1]]
}

write_placeholder_plot <- function(path, label, width = 12, height = 6) {
  grDevices::pdf(path, width = width, height = height)
  plot.new()
  text(0.5, 0.5, label)
  dev.off()
}

make_dataset_palette <- function(dataset_ids) {
  dataset_ids <- as.character(dataset_ids)
  setNames(scales::hue_pal()(length(dataset_ids)), dataset_ids)
}

chunk_dataset_ids <- function(dataset_ids, chunk_size = DETAIL_DATASETS_PER_PAGE) {
  split(dataset_ids, ceiling(seq_along(dataset_ids) / chunk_size))
}

page_subtitle <- function(dataset_ids, all_dataset_ids) {
  start_idx <- match(dataset_ids[[1]], all_dataset_ids)
  end_idx <- match(dataset_ids[[length(dataset_ids)]], all_dataset_ids)
  paste0("Datasets ", start_idx, "-", end_idx, " of ", length(all_dataset_ids))
}

line_plot_dimensions <- function(n_datasets, show_legend) {
  width <- if (show_legend) {
    min(18, max(10, 8 + ceiling(n_datasets / 4) * 1.2))
  } else {
    min(16, max(10, 8 + n_datasets * 0.08))
  }
  height <- if (show_legend) {
    min(12, max(6, 5 + ceiling(n_datasets / 6) * 0.7))
  } else {
    min(8, max(6, 5 + n_datasets * 0.02))
  }

  list(width = width, height = height)
}

error_distribution_dimensions <- function(n_datasets) {
  list(width = 12, height = max(6, 3 + n_datasets * 0.24))
}

write_plot_pages <- function(path, plot_list, width, height) {
  grDevices::pdf(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (plot_obj in plot_list) {
    print(plot_obj)
  }
}

build_retention_plot <- function(
  plot_df,
  dataset_ids,
  dataset_palette,
  y_col,
  y_label,
  title,
  subtitle = NULL,
  show_legend = TRUE
) {
  plot_df <- plot_df %>%
    filter(dataset_id %in% dataset_ids) %>%
    mutate(dataset_id = factor(dataset_id, levels = dataset_ids))

  p <- ggplot(
    plot_df,
    aes(x = stage, y = .data[[y_col]], color = dataset_id, group = dataset_id)
  ) +
    geom_line() +
    geom_point() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "",
      y = y_label,
      color = "Dataset"
    ) +
    scale_color_manual(
      values = dataset_palette[dataset_ids],
      limits = dataset_ids,
      breaks = dataset_ids
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 9)
    )

  if (identical(y_col, "pct_features")) {
    p <- p + scale_y_continuous(limits = c(0, 100))
  }

  if (show_legend) {
    p <- p +
      guides(color = guide_legend(ncol = if (length(dataset_ids) > 6) 2 else 1)) +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7)
      )
  } else {
    p <- p + theme(legend.position = "none")
  }

  p
}

build_threshold_plot <- function(
  plot_df,
  dataset_ids,
  dataset_palette,
  y_col,
  y_label,
  title,
  thresholds,
  subtitle = NULL,
  show_legend = TRUE
) {
  plot_df <- plot_df %>%
    filter(src_id %in% dataset_ids) %>%
    mutate(src_id = factor(src_id, levels = dataset_ids))

  p <- ggplot(
    plot_df,
    aes(x = threshold, y = .data[[y_col]], color = src_id, group = src_id)
  ) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = thresholds, linetype = "dashed", color = "grey80") +
    labs(
      title = title,
      subtitle = subtitle,
      x = paste0("Imputation error threshold (", error_metric, ")"),
      y = y_label,
      color = "Dataset"
    ) +
    scale_color_manual(
      values = dataset_palette[dataset_ids],
      limits = dataset_ids,
      breaks = dataset_ids
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 9)
    )

  if (identical(y_col, "perc_features_left")) {
    p <- p + scale_y_continuous(limits = c(0, 100))
  }

  if (show_legend) {
    p <- p +
      guides(color = guide_legend(ncol = if (length(dataset_ids) > 6) 2 else 1)) +
      theme(
        legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7)
      )
  } else {
    p <- p + theme(legend.position = "none")
  }

  p
}

make_retention_plot_df <- function(retention_summary_df) {
  retention_summary_df %>%
    transmute(
      dataset_id,
      Original = 100,
      `After preprocess` = pct_features_after_preprocess,
      `After imputation error filter` = pct_features_after_imputation_error_filter,
      `After distribution shift filter` = pct_features_after_distribution_shift_filter,
      Final = pct_features_final
    ) %>%
    pivot_longer(
      cols = -dataset_id,
      names_to = "stage",
      values_to = "pct_features"
    ) %>%
    mutate(
      stage = factor(
        stage,
        levels = c(
          "Original",
          "After preprocess",
          "After imputation error filter",
          "After distribution shift filter",
          "Final"
        )
      )
    )
}

make_retention_log10n_df <- function(retention_summary_df) {
  retention_summary_df %>%
    transmute(
      dataset_id,
      Original = n_features_orig,
      `After preprocess` = n_features_after_preprocess,
      `After imputation error filter` = n_features_after_imputation_error_filter,
      `After distribution shift filter` = n_features_after_distribution_shift_filter,
      Final = n_features_final
    ) %>%
    pivot_longer(
      cols = -dataset_id,
      names_to = "stage",
      values_to = "n_features"
    ) %>%
    mutate(
      stage = factor(
        stage,
        levels = c(
          "Original",
          "After preprocess",
          "After imputation error filter",
          "After distribution shift filter",
          "Final"
        )
      ),
      log10_n_features = log10(n_features + 1)
    )
}

# ----- MAIN -----

# A. Find benchmark and preprocessing summary files across datasets

benchmark_paths <- list.files(
  benchmark_dir,
  pattern = "benchmark.rds",
  recursive = TRUE,
  full.names = TRUE
)
summary_paths <- list.files(
  benchmark_dir,
  pattern = "preprocess_summary.csv",
  recursive = TRUE,
  full.names = TRUE
)

if (length(summary_paths) == 0) {
  stop("No preprocess_summary.csv files found in ", benchmark_dir)
}

# B. Combine feature-retention summaries across datasets

retention_summary_df <- lapply(summary_paths, function(path) {
  dataset_id <- basename(dirname(path))
  preprocess_summary_df <- read_csv(path, show_col_types = FALSE)
  stopifnot(nrow(preprocess_summary_df) == 1)

  preprocess_summary_df %>%
    mutate(dataset_id = dataset_id) %>%
    select(
      dataset_id,
      imputation_method,
      n_features_orig,
      n_features_after_preprocess,
      n_features_after_imputation_error_filter,
      n_features_after_distribution_shift_filter,
      n_features_final,
      pct_features_after_preprocess,
      pct_features_after_imputation_error_filter,
      pct_features_after_distribution_shift_filter,
      pct_features_final,
      distribution_shift_filter_applied
    )
}) %>%
  bind_rows() %>%
  arrange(dataset_id)

write_csv(
  retention_summary_df,
  file.path(out_dir, "feature_retention_summary.csv")
)

imputation_method_lookup <- setNames(
  retention_summary_df$imputation_method,
  retention_summary_df$dataset_id
)
dataset_ids <- retention_summary_df$dataset_id
dataset_palette <- make_dataset_palette(dataset_ids)
show_overview_legend <- length(dataset_ids) <= OVERVIEW_LEGEND_MAX_DATASETS
retention_dims <- line_plot_dimensions(length(dataset_ids), show_overview_legend)
retention_chunks <- chunk_dataset_ids(dataset_ids)

# C. Plot cross-dataset feature retention

retention_plot_df <- make_retention_plot_df(retention_summary_df)

p_retention <- build_retention_plot(
  retention_plot_df,
  dataset_ids = dataset_ids,
  dataset_palette = dataset_palette,
  y_col = "pct_features",
  y_label = "Percent of original features retained",
  title = "Feature retention across preprocessing and filtering",
  show_legend = show_overview_legend
)

ggsave(
  file.path(out_dir, "feature_retention_summary.pdf"),
  p_retention,
  width = retention_dims$width,
  height = retention_dims$height,
  limitsize = FALSE
)

write_plot_pages(
  file.path(out_dir, "feature_retention_summary_detailed.pdf"),
  lapply(retention_chunks, function(chunk_ids) {
    build_retention_plot(
      retention_plot_df,
      dataset_ids = chunk_ids,
      dataset_palette = dataset_palette,
      y_col = "pct_features",
      y_label = "Percent of original features retained",
      title = "Feature retention across preprocessing and filtering",
      subtitle = page_subtitle(chunk_ids, dataset_ids),
      show_legend = TRUE
    )
  }),
  width = 11,
  height = 8.5
)

retention_log10n_df <- make_retention_log10n_df(retention_summary_df)

p_retention_log10n <- build_retention_plot(
  retention_log10n_df,
  dataset_ids = dataset_ids,
  dataset_palette = dataset_palette,
  y_col = "log10_n_features",
  y_label = "log10(number of features retained + 1)",
  title = "Feature retention across preprocessing and filtering",
  show_legend = show_overview_legend
)

ggsave(
  file.path(out_dir, "feature_retention_summary_log10n.pdf"),
  p_retention_log10n,
  width = retention_dims$width,
  height = retention_dims$height,
  limitsize = FALSE
)

write_plot_pages(
  file.path(out_dir, "feature_retention_summary_log10n_detailed.pdf"),
  lapply(retention_chunks, function(chunk_ids) {
    build_retention_plot(
      retention_log10n_df,
      dataset_ids = chunk_ids,
      dataset_palette = dataset_palette,
      y_col = "log10_n_features",
      y_label = "log10(number of features retained + 1)",
      title = "Feature retention across preprocessing and filtering",
      subtitle = page_subtitle(chunk_ids, dataset_ids),
      show_legend = TRUE
    )
  }),
  width = 11,
  height = 8.5
)

# D. Write placeholder outputs if no benchmark results exist

if (length(benchmark_paths) == 0) {
  write_placeholder_plot(
    file.path(out_dir, "error_threshold.pdf"),
    "No benchmark.rds files found"
  )
  write_placeholder_plot(
    file.path(out_dir, "error_threshold_log10n.pdf"),
    "No benchmark.rds files found"
  )
  write_placeholder_plot(
    file.path(out_dir, "error_threshold_detailed.pdf"),
    "No benchmark.rds files found"
  )
  write_placeholder_plot(
    file.path(out_dir, "error_threshold_log10n_detailed.pdf"),
    "No benchmark.rds files found"
  )
  write_placeholder_plot(
    file.path(out_dir, "error_distribution.pdf"),
    "No benchmark.rds files found"
  )
  saveRDS(list(), file.path(out_dir, "error_metric_list.rds"))
  saveRDS(list(), file.path(out_dir, "features_remaining_after_imputation_list.rds"))
  quit(save = "no", status = 0)
}

# E. Recompute imputation-error summaries across datasets

src_ids <- basename(dirname(benchmark_paths))
names(benchmark_paths) <- src_ids

if (n_jobs > 1) {
  cl <- parallel::makeCluster(n_jobs)
  registerDoParallel(cl)
  `%op%` <- `%dopar%`
  print(paste0("Running with ", n_jobs, " cores."), quote = FALSE)
} else {
  `%op%` <- `%do%`
}

output_list <- foreach(
  src_id = src_ids,
  .inorder = FALSE,
  .packages = c("dplyr", "readr"),
  .export = c(
    "benchmark_paths",
    "error_metric",
    "feature_filename",
    "feature_table_dir",
    "filter_features_by_imputation_error_variant",
    "imputation_method",
    "imputation_method_lookup",
    "resolve_feature_path",
    "thresholds"
  )
) %op% {
  benchmark <- readRDS(benchmark_paths[[src_id]])

  if (is.null(benchmark)) {
    return(NULL)
  }

  feature_path <- resolve_feature_path(src_id, feature_table_dir, feature_filename)
  feature_mx <- readr::read_csv(feature_path, show_col_types = FALSE) %>%
    select(where(is.numeric))

  imputation_method_used <- imputation_method_lookup[[src_id]]
  if (is.null(imputation_method_used) || is.na(imputation_method_used) || imputation_method_used == "") {
    imputation_method_used <- imputation_method
  }
  if (is.null(imputation_method_used) || is.na(imputation_method_used) || imputation_method_used == "") {
    stop("No imputation method available for dataset: ", src_id)
  }

  feature_retention_by_threshold_persrc_list <- list()
  error_metrics <- NULL

  for (t in seq_along(thresholds)) {
    threshold <- thresholds[t]

    filtered <- filter_features_by_imputation_error_variant(
      feature_mx,
      benchmark,
      method = imputation_method_used,
      error_threshold = threshold,
      error_metric = error_metric
    )

    n_features_after_imputation_error_filter <- ncol(filtered$filtered_matrix)

    feature_retention_by_threshold_persrc_list[[paste0(src_id, ".", threshold, ".error")]] <- data.frame(
      src_id = src_id,
      threshold = threshold,
      stage = "after_imputation_error_filter",
      n_features_left = n_features_after_imputation_error_filter,
      perc_features_left = n_features_after_imputation_error_filter / ncol(feature_mx) * 100
    )

    if (t == 1) {
      error_metrics <- data.frame(
        src_id = src_id,
        error = filtered$summary$max_rmse
      )
    }
  }

  list(
    error_metrics = error_metrics,
    feature_retention_by_threshold = bind_rows(feature_retention_by_threshold_persrc_list)
  )
}

if (n_jobs > 1) {
  parallel::stopCluster(cl)
}

output_list <- Filter(Negate(is.null), output_list)

if (length(output_list) == 0) {
  write_placeholder_plot(
    file.path(out_dir, "error_threshold.pdf"),
    "No imputed datasets with benchmark results found"
  )
  write_placeholder_plot(
    file.path(out_dir, "error_threshold_log10n.pdf"),
    "No imputed datasets with benchmark results found"
  )
  write_placeholder_plot(
    file.path(out_dir, "error_threshold_detailed.pdf"),
    "No imputed datasets with benchmark results found"
  )
  write_placeholder_plot(
    file.path(out_dir, "error_threshold_log10n_detailed.pdf"),
    "No imputed datasets with benchmark results found"
  )
  write_placeholder_plot(
    file.path(out_dir, "error_distribution.pdf"),
    "No imputed datasets with benchmark results found"
  )
  saveRDS(list(), file.path(out_dir, "error_metric_list.rds"))
  saveRDS(list(), file.path(out_dir, "features_remaining_after_imputation_list.rds"))
  quit(save = "no", status = 0)
}

# F. Save aggregate benchmark summaries

error_metric_list <- unlist(
  lapply(output_list, `[`, "error_metrics"),
  recursive = FALSE
)
feature_retention_by_threshold_list <- unlist(
  lapply(output_list, `[`, "feature_retention_by_threshold"),
  recursive = FALSE
)

saveRDS(error_metric_list, file.path(out_dir, "error_metric_list.rds"))
saveRDS(
  feature_retention_by_threshold_list,
  file.path(out_dir, "features_remaining_after_imputation_list.rds")
)

# G. Plot retained features versus error threshold

p_df <- bind_rows(feature_retention_by_threshold_list) %>%
  mutate(log10_n_features_left = log10(n_features_left + 1))
threshold_dataset_ids <- dataset_ids[dataset_ids %in% unique(p_df$src_id)]
threshold_palette <- dataset_palette[threshold_dataset_ids]
threshold_chunks <- chunk_dataset_ids(threshold_dataset_ids)
threshold_show_overview_legend <- length(threshold_dataset_ids) <= OVERVIEW_LEGEND_MAX_DATASETS
threshold_dims <- line_plot_dimensions(length(threshold_dataset_ids), threshold_show_overview_legend)

p_error_threshold <- build_threshold_plot(
  p_df,
  dataset_ids = threshold_dataset_ids,
  dataset_palette = threshold_palette,
  y_col = "perc_features_left",
  y_label = "Percent of features retained",
  title = paste0(
    "Features remaining after filtering by imputation error (", error_metric, ")"
  ),
  thresholds = thresholds,
  show_legend = threshold_show_overview_legend
)

ggsave(
  file.path(out_dir, "error_threshold.pdf"),
  p_error_threshold,
  width = threshold_dims$width,
  height = threshold_dims$height,
  limitsize = FALSE
)

write_plot_pages(
  file.path(out_dir, "error_threshold_detailed.pdf"),
  lapply(threshold_chunks, function(chunk_ids) {
    build_threshold_plot(
      p_df,
      dataset_ids = chunk_ids,
      dataset_palette = threshold_palette,
      y_col = "perc_features_left",
      y_label = "Percent of features retained",
      title = paste0(
        "Features remaining after filtering by imputation error (", error_metric, ")"
      ),
      thresholds = thresholds,
      subtitle = page_subtitle(chunk_ids, threshold_dataset_ids),
      show_legend = TRUE
    )
  }),
  width = 11,
  height = 8.5
)

p_error_threshold_log10n <- build_threshold_plot(
  p_df,
  dataset_ids = threshold_dataset_ids,
  dataset_palette = threshold_palette,
  y_col = "log10_n_features_left",
  y_label = "log10(number of features left + 1)",
  title = paste0(
    "Features remaining after filtering by imputation error (", error_metric, ")"
  ),
  thresholds = thresholds,
  show_legend = threshold_show_overview_legend
)

ggsave(
  file.path(out_dir, "error_threshold_log10n.pdf"),
  p_error_threshold_log10n,
  width = threshold_dims$width,
  height = threshold_dims$height,
  limitsize = FALSE
)

write_plot_pages(
  file.path(out_dir, "error_threshold_log10n_detailed.pdf"),
  lapply(threshold_chunks, function(chunk_ids) {
    build_threshold_plot(
      p_df,
      dataset_ids = chunk_ids,
      dataset_palette = threshold_palette,
      y_col = "log10_n_features_left",
      y_label = "log10(number of features left + 1)",
      title = paste0(
        "Features remaining after filtering by imputation error (", error_metric, ")"
      ),
      thresholds = thresholds,
      subtitle = page_subtitle(chunk_ids, threshold_dataset_ids),
      show_legend = TRUE
    )
  }),
  width = 11,
  height = 8.5
)

# H. Plot feature-level imputation error distributions

p_df <- bind_rows(error_metric_list) %>%
  mutate(src_id = factor(src_id, levels = rev(threshold_dataset_ids)))
stopifnot(all(p_df$error > 0))

plot_thresholds <- thresholds[thresholds > 0]
error_dims <- error_distribution_dimensions(length(threshold_dataset_ids))
p <- ggplot(p_df, aes(y = src_id, x = log10(error))) +
  geom_boxplot() +
  geom_vline(xintercept = log10(plot_thresholds), linetype = "dashed", color = "grey80") +
  labs(
    title = "Distribution of Normalised RMSE for Features with NA Values",
    x = paste0("log10 imputation error (", error_metric, ")"),
    y = ""
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10),
    axis.text.y = element_text(size = if (length(threshold_dataset_ids) > 40) 5 else 7)
  )

ggsave(
  file.path(out_dir, "error_distribution.pdf"),
  p,
  width = error_dims$width,
  height = error_dims$height,
  limitsize = FALSE
)
