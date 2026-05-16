#!/usr/bin/env Rscript

# Plot features in modules identified by hierarchical clustering

renv::load()
setwd(here::here())
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(ggplot2)
  library(gridExtra)
  library(vdjremix)
})
source("code/feature_module_reduction/io_helpers.R")

# Functions

plot_module_loadings_variant <- function(
    loadings,
    module,
    feature_groupings = NULL,
    top_n = NULL,
    important_only = FALSE,
    top_n_per_group = NULL
) {
  stopifnot(is.list(loadings))
  stopifnot(is.list(feature_groupings) || is.null(feature_groupings))
  stopifnot(module %in% names(loadings))

  w <- loadings[[module]]
  df <- data.frame(
    feature = names(w),
    loading = as.numeric(w)
  )

  if (!is.null(feature_groupings)) {
    g <- feature_groupings[[module]]
    stopifnot(length(g) == nrow(df))
    df$group <- g
  }

  df <- df[order(abs(df$loading), decreasing = TRUE), ]

  p <- nrow(df)
  cutoff <- sqrt(1 / p)
  df$important <- abs(df$loading) >= cutoff

  if (!is.null(top_n_per_group) && !is.null(feature_groupings)) {
    df <- df %>%
      dplyr::group_by(group) %>%
      dplyr::slice_max(order_by = abs(loading), n = top_n_per_group) %>%
      dplyr::ungroup()
  } else if (important_only) {
    df <- df[df$important, ]
  } else if (!is.null(top_n)) {
    df <- df[1:min(top_n, nrow(df)), ]
  }

  df$feature <- factor(df$feature, levels = df$feature)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = feature, y = loading, colour = important)
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_hline(
      yintercept = c(-cutoff, cutoff),
      linetype = "dotted",
      colour = "red"
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Feature",
      y = "PC1",
      colour = "Important\ncontributor",
      title = paste0("m", module, ": feature loadings")
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )
    )

  if (!is.null(feature_groupings)) {
    p <- p + ggplot2::facet_grid(~group, scales = "free_x")
  }

  p
}

# Parameters

parser <- ArgumentParser(description = "Plot module loadings and feature distributions")
parser$add_argument("--plot_feature_path", type = "character", required = TRUE)
parser$add_argument("--module_feature_path", type = "character", default = NULL)
parser$add_argument("--eigengenes_path", type = "character", default = NULL)
parser$add_argument("--sample_metadata_path", type = "character", required = TRUE)
parser$add_argument("--module_annotations_path", type = "character", default = NULL)
parser$add_argument("--out_dir", type = "character", required = TRUE)
parser$add_argument("--sample_key", type = "character", required = TRUE)
parser$add_argument("--group_key", type = "character", required = TRUE)
parser$add_argument("--module_prefix", type = "character", required = TRUE)
parser$add_argument("--plot_group_order", type = "character", nargs = "+", required = TRUE)
args <- parser$parse_args()

plot_feature_path <- args$plot_feature_path
module_feature_path <- args$module_feature_path
eigengenes_path <- args$eigengenes_path
sample_metadata_path <- args$sample_metadata_path
module_annotations_path <- args$module_annotations_path
out_dir <- args$out_dir
sample_key <- args$sample_key
group_key <- args$group_key
module_prefix <- args$module_prefix
plot_group_order <- args$plot_group_order

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- MAIN -----

# B. Load plot feature table and (optionally) module feature table.

plot_feature_info <- read_feature_table_generic(
  plot_feature_path,
  sample_key = sample_key
)
plot_feature_df <- plot_feature_info$feature_df

if (is.null(module_feature_path)) {
  module_feature_path <- plot_feature_path
}

module_feature_info <- read_feature_table_generic(
  module_feature_path,
  sample_key = sample_key
)
module_feature_df <- module_feature_info$feature_df

# A. Load eigengene object or build a synthetic one (all features as one module).

if (!is.null(eigengenes_path)) {
  eig <- readRDS(eigengenes_path)
  p <- plot_variance_explained(
    eig$variance_explained, reorder = TRUE
  )
  ggsave(
    filename = file.path(out_dir, "variance_explained.svg"),
    plot = p, width = 6, height = 4
  )
} else {
  feature_names <- setdiff(colnames(plot_feature_df), sample_key)
  n <- length(feature_names)
  eig <- list(
    variance_explained = NULL,
    loadings = list("1" = setNames(rep(1 / sqrt(n), n), feature_names))
  )
}

plot_metadata_df <- read_metadata_aligned(
  sample_metadata_path,
  sample_key = sample_key,
  sample_ids = plot_feature_df[[sample_key]],
  required_cols = group_key
)
module_metadata_df <- read_metadata_aligned(
  sample_metadata_path,
  sample_key = sample_key,
  sample_ids = module_feature_df[[sample_key]],
  required_cols = group_key
)

plot_feature_df <- plot_feature_df %>%
  left_join(
    plot_metadata_df %>%
      select(all_of(c(sample_key, group_key))),
    by = sample_key
  ) %>%
  mutate(!!group_key := factor(.data[[group_key]], levels = plot_group_order))

module_feature_df <- module_feature_df %>%
  left_join(
    module_metadata_df %>%
      select(all_of(c(sample_key, group_key))),
    by = sample_key
  )

module_groups_present <- unique(as.character(module_feature_df[[group_key]]))
module_groups_present <- module_groups_present[!is.na(module_groups_present)]
module_group_order <- plot_group_order[plot_group_order %in% module_groups_present]
stopifnot(length(module_group_order) > 0)

module_feature_df <- module_feature_df %>%
  mutate(!!group_key := factor(.data[[group_key]], levels = module_group_order))

module_score_lookup <- c()
if (!is.null(module_annotations_path) && nzchar(module_annotations_path)) {
  module_annotations_df <- read_csv(module_annotations_path, show_col_types = FALSE)
  stopifnot(all(c("module_id", "renamed_module") %in% colnames(module_annotations_df)))
  module_score_lookup <- setNames(
    module_annotations_df$renamed_module,
    module_annotations_df$module_id
  )
}

# C. Build feature groupings for grouped loadings plots.

feature_groupings <- sapply(
  eig$loadings, USE.NAMES = TRUE, simplify = FALSE, function(x) {
    stringr::str_remove(names(x), ".*?__")
  }
)

# D. Plot each module using module scores from the reduced table and feature values from the plot table.

for (m in seq_along(eig$loadings)) {
  m_index <- names(eig$loadings)[m]

  message("Plotting features for module ", m_index)

  out_m_dir <- file.path(out_dir, paste0("m", m_index))
  dir.create(out_m_dir, recursive = TRUE, showWarnings = FALSE)

  p <- vdjremix::plot_module_loadings(
    eig$loadings, module = m_index
  )
  p_width <- 10 / 50 * nrow(p$data)
  ggsave(
    filename = file.path(out_m_dir, "loadings.svg"),
    plot = p, width = p_width, height = 10, limitsize = FALSE
  )

  p <- plot_module_loadings_variant(
    eig$loadings, module = m_index, feature_groupings = feature_groupings,
    top_n = 50
  ) +
    theme(
      strip.text = element_blank(),
      strip.background = element_blank()
    )
  p_width <- 15 / 50 * nrow(p$data)
  ggsave(
    filename = file.path(out_m_dir, "loadings_grouped.svg"),
    plot = p, width = p_width, height = 6, limitsize = FALSE
  )

  m_loadings_sorted <- sort(eig$loadings[[m_index]], decreasing = TRUE)
  features <- names(m_loadings_sorted)
  p_feature_list <- list()

  raw_module_score_col <- paste0(module_prefix, m_index)
  module_score_col <- module_score_lookup[[raw_module_score_col]]
  if (is.null(module_score_col) || is.na(module_score_col) || module_score_col == "") {
    module_score_col <- raw_module_score_col
  }

  if (module_score_col %in% colnames(module_feature_df)) {
    features_to_plot <- c(module_score_col, features)
  } else {
    features_to_plot <- features
  }

  for (feature in features_to_plot) {
    message("Plotting feature ", feature, " in module ", m_index)

    if (identical(feature, module_score_col)) {
      feature_plot_df <- module_feature_df %>%
        filter(
          !is.na(.data[[group_key]]),
          !is.na(.data[[feature]])
        ) %>%
        droplevels()
    } else {
      feature_plot_df <- plot_feature_df %>%
        filter(
          !is.na(.data[[group_key]]),
          !is.na(.data[[feature]])
        ) %>%
        droplevels()
    }

    p_feature_list[[feature]] <- ggplot(
      feature_plot_df,
      aes(x = .data[[group_key]], y = .data[[feature]])
    ) +
      geom_violin(fill = "lightblue", alpha = 0.5, trim = TRUE) +
      geom_boxplot(width = 0.2) +
      geom_point(alpha = 0.5, size = 1) +
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 3,
        size = 2,
        color = "red"
      ) +
      theme_bw() +
      labs(
        y = feature, x = group_key, title = feature
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 10)
      )
  }

  nrow <- 5
  ncol <- 5
  layout_mx <- matrix(
    1:(nrow * ncol), nrow = nrow, ncol = ncol, byrow = TRUE
  )
  ggsave(
    filename = file.path(out_m_dir, "features.pdf"),
    plot = marrangeGrob(
      p_feature_list,
      layout_matrix = layout_mx,
      top = paste0("m", m_index, " features (ordered by loadings ignoring groups)")
    ),
    width = 20, height = 15
  )

  message("Finished plotting features for module ", m_index)
}

message("Finished plotting features for all modules")
