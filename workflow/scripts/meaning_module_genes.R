#!/usr/bin/env Rscript

# %% Ordered gene-set meaning with online g:Profiler or offline mHG/GMT

setwd(here::here())
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
})
source("code/feature_module_reduction/helpers.R")

# Functions

run_ordered_gprofiler <- function(df, organism = "hsapiens") {
  stopifnot(is.data.frame(df))
  stopifnot(all(c("gene", "loading") %in% colnames(df)))
  stopifnot(!anyDuplicated(df$gene))

  ranked_genes <- df %>%
    arrange(desc(loading)) %>%
    pull(gene)

  gprofiler2::gost(
    query = ranked_genes,
    organism = organism,
    ordered_query = TRUE,
    significant = FALSE,
    exclude_iea = TRUE,
    measure_underrepresentation = FALSE,
    evcodes = TRUE,
    user_threshold = 0.05,
    correction_method = "fdr",
    sources = c("GO:BP")
  )
}

format_gmt_term_name <- function(term_id) {
  term_id %>%
    sub("^GOBP_", "", .) %>%
    gsub("_", " ", ., fixed = TRUE) %>%
    tolower()
}

run_mhg_gmt <- function(df, gmt_terms, min_hits = 2L, max_frac = 0.5) {
  stopifnot(is.data.frame(df))
  stopifnot(all(c("gene", "loading") %in% colnames(df)))
  stopifnot(!anyDuplicated(df$gene))

  ranked_genes <- df %>%
    arrange(desc(loading)) %>%
    pull(gene)

  term_results <- lapply(gmt_terms, function(term) {
    term_members <- unique(term$members)
    term_members <- term_members[!is.na(term_members) & term_members != ""]
    intersection_genes <- intersect(ranked_genes, term_members)

    if (length(term_members) == 0) {
      mhg_res <- data.frame(
        status = "empty_term",
        n_query = length(ranked_genes),
        n_term_total = 0L,
        n_term_in_query = 0L,
        n_term_missing = 0L,
        min_hits = min_hits,
        max_k = NA_integer_,
        n_cutoffs_tested = 0L,
        best_k = NA_integer_,
        hits_at_k = NA_integer_,
        expected_at_k = NA_real_,
        fold_enrichment = NA_real_,
        scan_pvalue = NA_real_,
        bonferroni_pvalue = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      mhg_res <- test_mhg_enrichment(
        query = ranked_genes,
        term = term_members,
        min_hits = min_hits,
        max_frac = max_frac
      )
    }

    tibble(
      term_id = term$term_id,
      term_name = format_gmt_term_name(term$term_id),
      term_size = length(term_members),
      intersection = paste(intersection_genes, collapse = ","),
      intersection_size = length(intersection_genes),
      raw_p_value = mhg_res$bonferroni_pvalue[[1]],
      p_value = NA_real_,
      significant = FALSE,
      parents = NA_character_,
      status = mhg_res$status[[1]],
      best_k = mhg_res$best_k[[1]],
      hits_at_k = mhg_res$hits_at_k[[1]],
      fold_enrichment = mhg_res$fold_enrichment[[1]],
      n_query = mhg_res$n_query[[1]],
      n_term_total = mhg_res$n_term_total[[1]],
      n_term_in_query = mhg_res$n_term_in_query[[1]],
      n_term_missing = mhg_res$n_term_missing[[1]],
      min_hits = mhg_res$min_hits[[1]],
      max_k = mhg_res$max_k[[1]],
      n_cutoffs_tested = mhg_res$n_cutoffs_tested[[1]],
      scan_pvalue = mhg_res$scan_pvalue[[1]]
    )
  })

  result_df <- bind_rows(term_results)
  valid <- result_df$status == "ok" & !is.na(result_df$raw_p_value)
  result_df$p_value[valid] <- p.adjust(result_df$raw_p_value[valid], method = "BH")
  result_df$significant <- !is.na(result_df$p_value) & result_df$p_value < 0.05
  result_df
}

sanitize_label <- function(x) {
  x %>%
    gsub(" ", "_", ., fixed = TRUE) %>%
    gsub("[^A-Za-z0-9_.-]", "_", .)
}

sanitize_contract_field <- function(x, default = "unknown") {
  x <- sanitize_label(x) %>%
    gsub("_+", "_", .) %>%
    gsub("^_+|_+$", "", .)

  if (is.na(x) || x == "") {
    return(default)
  }
  x
}

extract_enrichment_id <- function(feature, feature_id_separator) {
  if (is.na(feature_id_separator) || feature_id_separator == "") {
    stop("feature_id_separator must be a non-empty string")
  }

  vapply(
    strsplit(as.character(feature), feature_id_separator, fixed = TRUE),
    function(parts) parts[[1]],
    character(1)
  )
}

deduplicate_enrichment_ids <- function(df, duplicate_enrichment_id_policy) {
  if (duplicate_enrichment_id_policy != "max_abs_loading") {
    stop(
      "duplicate_enrichment_id_policy must be max_abs_loading; got: ",
      duplicate_enrichment_id_policy
    )
  }

  df %>%
    group_by(gene) %>%
    slice_max(order_by = abs(loading), n = 1, with_ties = FALSE) %>%
    ungroup()
}

add_top_gene_overlap <- function(df, top_genes) {
  if (nrow(df) == 0) {
    return(df %>% mutate(contains_top_n_genes = logical()))
  }

  df %>%
    mutate(
      contains_top_n_genes = vapply(intersection, function(x) {
        if (is.na(x) || x == "") {
          return(FALSE)
        }
        overlap_genes <- strsplit(x, ",")[[1]]
        any(top_genes %in% overlap_genes)
      }, logical(1))
    )
}

select_pathway <- function(result_df, top_genes) {
  empty_result <- list(
    selected = tibble(),
    primary = tibble(),
    relaxed = tibble(),
    status = "none"
  )

  if (is.null(result_df) || nrow(result_df) == 0) {
    return(empty_result)
  }

  result_df <- result_df %>%
    add_top_gene_overlap(top_genes)

  if ("status" %in% colnames(result_df)) {
    result_df <- result_df %>%
      filter(status == "ok")
  }

  primary_df <- result_df %>%
    filter(term_size >= 5 & term_size <= 250 & significant & intersection_size >= 2) %>%
    # Prefer terms covering more module genes; p-value breaks ties.
    arrange(desc(intersection_size), p_value)

  if (nrow(primary_df) > 0) {
    return(list(
      selected = slice_head(primary_df, n = 1),
      primary = primary_df,
      relaxed = tibble(),
      status = "significant"
    ))
  }

  relaxed_df <- result_df %>%
    filter(term_size >= 5 & term_size <= 500 & intersection_size >= 2) %>%
    arrange(desc(intersection_size), p_value)

  if (nrow(relaxed_df) > 0) {
    return(list(
      selected = slice_head(relaxed_df, n = 1),
      primary = primary_df,
      relaxed = relaxed_df,
      status = "suggestive"
    ))
  }

  empty_result
}

format_pathway_parents <- function(parents) {
  parent_text <- paste(parents, collapse = "__")
  if (is.na(parent_text) || parent_text == "") {
    return("none_parent")
  }

  parent_text %>%
    gsub(",", "__", ., fixed = TRUE) %>%
    sanitize_label()
}

parse_module_context <- function(module_id) {
  parts <- strsplit(module_id, "..", fixed = TRUE)[[1]]

  module_short <- sanitize_label(module_id)
  source_id <- "global__notapp.module"

  if (length(parts) >= 4 && parts[[length(parts) - 1]] == "m") {
    module_short <- paste0("m", parts[[length(parts)]])
    source_id <- parts[[length(parts) - 2]]
  }

  source_parts <- strsplit(source_id, "__", fixed = TRUE)[[1]]
  if (length(source_parts) >= 2) {
    cell_type <- paste(source_parts[-length(source_parts)], collapse = "__")
    level_feature_type <- source_parts[[length(source_parts)]]
  } else {
    cell_type <- source_id
    level_feature_type <- "notapp.module"
  }

  level_parts <- strsplit(level_feature_type, ".", fixed = TRUE)[[1]]
  if (length(level_parts) >= 2) {
    cell_level <- paste(level_parts[-length(level_parts)], collapse = ".")
    feature_type <- level_parts[[length(level_parts)]]
  } else {
    cell_level <- level_feature_type
    feature_type <- "module"
  }

  list(
    module_short = sanitize_contract_field(module_short),
    cell_type = sanitize_contract_field(cell_type),
    cell_level = sanitize_contract_field(cell_level),
    feature_type = sanitize_contract_field(feature_type)
  )
}

build_renamed_module <- function(module_id, top_genes_text, top_pathways_text) {
  module_context <- parse_module_context(module_id)
  feature_label <- paste(
    top_genes_text,
    top_pathways_text,
    module_context$module_short,
    sep = "_"
  ) %>%
    sanitize_contract_field()

  paste(
    feature_label,
    module_context$cell_type,
    module_context$cell_level,
    module_context$feature_type,
    sep = "__"
  )
}

# Parameters

parser <- ArgumentParser(description = "Add gene-set meaning to module columns")
parser$add_argument("--module_feature_path", type = "character", required = TRUE)
parser$add_argument("--modules_dir", type = "character", required = TRUE)
parser$add_argument("--out_dir", type = "character", required = TRUE)
parser$add_argument("--backend", type = "character", default = "gprofiler")
parser$add_argument("--gmt_path", type = "character", default = "")
parser$add_argument("--min_hits", type = "integer", default = 2)
parser$add_argument("--max_frac", type = "double", default = 0.5)
parser$add_argument("--organism", type = "character", default = "hsapiens")
parser$add_argument("--top_n", type = "integer", default = 2)
parser$add_argument("--use_only_important_loading", type = "character", default = "FALSE")
parser$add_argument("--feature_id_separator", type = "character", default = "__")
parser$add_argument(
  "--duplicate_enrichment_id_policy",
  type = "character",
  default = "max_abs_loading"
)
args <- parser$parse_args()

module_feature_path <- args$module_feature_path
modules_dir <- args$modules_dir
out_dir <- args$out_dir
backend <- args$backend
gmt_path <- args$gmt_path
min_hits <- args$min_hits
max_frac <- args$max_frac
organism <- args$organism
top_n <- args$top_n
feature_id_separator <- args$feature_id_separator
duplicate_enrichment_id_policy <- args$duplicate_enrichment_id_policy
use_only_important_loading <- parse_bool(
  args$use_only_important_loading,
  "use_only_important_loading"
)

if (!(backend %in% c("gprofiler", "mhg_gmt"))) {
  stop("backend must be one of gprofiler, mhg_gmt; got: ", backend)
}
if (is.na(top_n) || top_n < 1L) {
  stop("top_n must be a positive integer")
}
if (is.na(min_hits) || min_hits < 1L) {
  stop("min_hits must be a positive integer")
}
if (is.na(max_frac) || max_frac <= 0 || max_frac > 1) {
  stop("max_frac must be > 0 and <= 1")
}

gmt_terms <- NULL
if (backend == "mhg_gmt") {
  if (is.na(gmt_path) || gmt_path == "") {
    stop("gmt_path is required when backend is mhg_gmt")
  }
  source(file.path("utils", "test_mhg_enrichment.R"))
  gmt_terms <- read_gmt(gmt_path)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
module_out_dir <- file.path(out_dir, "modules")
dir.create(module_out_dir, recursive = TRUE, showWarnings = FALSE)

# ----- MAIN -----

# A. Load module table and discover per-module ranking CSVs.

module_feature_df <- read_csv(module_feature_path, show_col_types = FALSE)

module_csv_paths <- c(
  list.files(modules_dir, pattern = "\\.csv$", full.names = TRUE),
  list.files(modules_dir, pattern = "members\\.csv$", recursive = TRUE, full.names = TRUE)
)
module_csv_paths <- unique(module_csv_paths)

if (length(module_csv_paths) == 0) {
  stop("No module CSV files found in modules_dir: ", modules_dir)
}

module_ids <- sapply(module_csv_paths, function(x) {
  if (basename(x) == "members.csv") {
    basename(dirname(x))
  } else {
    tools::file_path_sans_ext(basename(x))
  }
}, USE.NAMES = FALSE)

stopifnot(!anyDuplicated(module_ids))

# B. Run pathway annotation per module and collect rename labels.

annotation_df_list <- list()

for (i in seq_along(module_csv_paths)) {
  module_id <- module_ids[[i]]
  module_path <- module_csv_paths[[i]]

  message("Annotating module: ", module_id)

  stopifnot(module_id %in% colnames(module_feature_df))

  # Keep the full feature ID as official; enrichment uses the left side only.
  genes_df <- read_csv(module_path, show_col_types = FALSE) %>%
    mutate(gene = extract_enrichment_id(feature, feature_id_separator)) %>%
    arrange(desc(loading))

  stopifnot(all(c("feature", "gene", "loading") %in% colnames(genes_df)))
  if (any(is.na(genes_df$gene) | genes_df$gene == "")) {
    stop("Blank enrichment gene IDs found for module: ", module_id)
  }

  if (use_only_important_loading) {
    stopifnot("important_loading" %in% colnames(genes_df))
    genes_df <- genes_df %>%
      filter(important_loading)
  }

  genes_df <- genes_df %>%
    deduplicate_enrichment_ids(duplicate_enrichment_id_policy) %>%
    arrange(desc(loading))

  stopifnot(nrow(genes_df) > 0)

  out_dir_i <- file.path(module_out_dir, module_id)
  dir.create(out_dir_i, recursive = TRUE, showWarnings = FALSE)

  write_csv(
    genes_df,
    file.path(out_dir_i, "ranked_genes.csv")
  )

  top_genes <- genes_df %>%
    slice_head(n = top_n) %>%
    pull(gene) %>%
    sort()

  pathway_result_df <- NULL
  if (backend == "gprofiler") {
    gost_res <- tryCatch(
      run_ordered_gprofiler(genes_df, organism = organism),
      error = function(e) {
        message("g:Profiler failed for ", module_id, ": ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(gost_res) && !is.null(gost_res$result)) {
      pathway_result_df <- gost_res$result
    }
  } else if (backend == "mhg_gmt") {
    pathway_result_df <- run_mhg_gmt(
      genes_df,
      gmt_terms = gmt_terms,
      min_hits = min_hits,
      max_frac = max_frac
    )
    write_csv(pathway_result_df, file.path(out_dir_i, "mhg_gmt.results.csv"))
  }

  result_filtered_df <- tibble()
  result_relaxed_df <- tibble()
  top_pathways_text <- "nonesig"
  top_pathway_parent_text <- "none_parent"
  top_pathway_status <- "none"
  top_pathway_p_value <- NA_real_
  top_pathway_intersection_size <- NA_integer_
  top_pathway_contains_top_genes <- NA

  if (!is.null(pathway_result_df)) {
    pathway_choice <- select_pathway(pathway_result_df, top_genes)
    result_filtered_df <- pathway_choice$primary
    result_relaxed_df <- pathway_choice$relaxed
    top_pathway_status <- pathway_choice$status

    if (nrow(result_filtered_df) > 0) {
      write_csv(
        result_filtered_df,
        file.path(out_dir_i, "ordered_ora.filtered_results.csv")
      )
    }

    if (nrow(result_relaxed_df) > 0) {
      write_csv(
        result_relaxed_df,
        file.path(out_dir_i, "ordered_ora.relaxed_results.csv")
      )
    }

    if (nrow(pathway_choice$selected) > 0) {
      top_pathway_row <- pathway_choice$selected

      top_pathways_text <- top_pathway_row %>%
        pull(term_name) %>%
        sanitize_label()

      top_pathway_parent_text <- top_pathway_row %>%
        pull(parents) %>%
        format_pathway_parents()

      top_pathway_p_value <- top_pathway_row %>%
        pull(p_value)

      top_pathway_intersection_size <- top_pathway_row %>%
        pull(intersection_size)

      top_pathway_contains_top_genes <- top_pathway_row %>%
        pull(contains_top_n_genes)
    }
  }

  top_genes_text <- paste(top_genes, collapse = "_") %>%
    sanitize_label()
  meaning_label <- paste(top_genes_text, top_pathways_text, sep = "_") %>%
    sanitize_label()

  writeLines(
    meaning_label,
    con = file.path(out_dir_i, "meaning.txt")
  )

  annotation_df_list[[module_id]] <- tibble(
    module_id = module_id,
    module_csv_path = module_path,
    module_meaning_backend = backend,
    module_meaning_gmt_path = ifelse(backend == "mhg_gmt", gmt_path, NA_character_),
    top_genes = top_genes_text,
    top_pathway = top_pathways_text,
    top_pathway_parents = top_pathway_parent_text,
    top_pathway_status = top_pathway_status,
    top_pathway_p_value = top_pathway_p_value,
    top_pathway_intersection_size = top_pathway_intersection_size,
    top_pathway_contains_top_genes = top_pathway_contains_top_genes,
    meaning_label = meaning_label,
    renamed_module = build_renamed_module(
      module_id,
      top_genes_text,
      top_pathways_text
    )
  )
}

annotation_df <- bind_rows(annotation_df_list)
write_csv(annotation_df, file.path(out_dir, "module_annotations.csv"))

# C. Rename module columns and save the interpreted module table.

rename_lookup <- annotation_df$module_id
names(rename_lookup) <- annotation_df$renamed_module

module_feature_df %>%
  rename(all_of(rename_lookup)) %>%
  write_csv(file.path(out_dir, "bulk_x_features.csv"))
