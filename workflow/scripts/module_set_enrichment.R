# %% Test ranked module members for enrichment of external feature sets.

# %% Setup

renv::load()
setwd(here::here())
suppressPackageStartupMessages({
  library(argparse)
  source("utils/test_mhg_enrichment.R")
  source("code/feature_module_reduction/helpers.R")
})

# %% Functions

extract_module_member_feature_ids <- function(feature, feature_id_separator) {
  if (is.na(feature_id_separator) || feature_id_separator == "") {
    stop("feature_id_separator must be a non-empty string")
  }
  vapply(
    strsplit(as.character(feature), feature_id_separator, fixed = TRUE),
    function(parts) parts[[1]],
    character(1)
  )
}

build_query_ids <- function(feature, feature_id_format, feature_id_separator) {
  if (feature_id_format == "module_member_feature_ids") {
    return(extract_module_member_feature_ids(feature, feature_id_separator))
  }
  if (feature_id_format == "full_feature_ids") {
    return(as.character(feature))
  }
  stop(
    "feature_id_format must be module_member_feature_ids or full_feature_ids; got: ",
    feature_id_format
  )
}

deduplicate_query_ids <- function(df, duplicate_enrichment_id_policy) {
  if (duplicate_enrichment_id_policy != "max_abs_loading") {
    stop(
      "duplicate_enrichment_id_policy must be max_abs_loading; got: ",
      duplicate_enrichment_id_policy
    )
  }
  df <- df[order(-df$abs_loading, df$query_id, df$feature), , drop = FALSE]
  df[!duplicated(df$query_id), , drop = FALSE]
}

empty_mhg_result <- function(query, term, status, min_hits, max_frac) {
  n_query <- length(query)
  n_term_total <- length(term)
  n_term_in_query <- length(intersect(query, term))
  data.frame(
    status = status,
    n_query = n_query,
    n_term_total = n_term_total,
    n_term_in_query = n_term_in_query,
    n_term_missing = n_term_total - n_term_in_query,
    min_hits = min_hits,
    max_k = if (n_query > 0) min(n_query, as.integer(floor(max_frac * n_query))) else NA_integer_,
    n_cutoffs_tested = 0L,
    best_k = NA_integer_,
    hits_at_k = NA_integer_,
    expected_at_k = NA_real_,
    fold_enrichment = NA_real_,
    scan_pvalue = NA_real_,
    bonferroni_pvalue = NA_real_,
    stringsAsFactors = FALSE
  )
}

read_module_query <- function(
    module_path,
    feature_id_format,
    feature_id_separator,
    duplicate_enrichment_id_policy,
    use_only_important_loading) {
  module_df <- read.csv(
    module_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c("feature", "loading")
  missing_cols <- setdiff(required_cols, colnames(module_df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in module file ", module_path, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (use_only_important_loading) {
    if (!"important_loading" %in% colnames(module_df)) {
      stop("Missing important_loading column in module file: ", module_path)
    }
    module_df <- module_df[
      module_df$important_loading %in% c(TRUE, "TRUE", "true", "1", 1),
      ,
      drop = FALSE
    ]
  }

  module_df$loading <- as.numeric(module_df$loading)
  if (any(is.na(module_df$loading))) {
    stop("NA loading values found in module file: ", module_path)
  }

  module_df$query_id <- build_query_ids(
    module_df$feature,
    feature_id_format,
    feature_id_separator
  )
  if (any(is.na(module_df$query_id) | module_df$query_id == "")) {
    stop("Blank query_id values found in module file: ", module_path)
  }

  module_df$abs_loading <- abs(module_df$loading)
  module_df <- module_df[order(-module_df$abs_loading, module_df$feature), , drop = FALSE]
  module_df <- deduplicate_query_ids(module_df, duplicate_enrichment_id_policy)
  module_df <- module_df[order(-module_df$abs_loading, module_df$feature), , drop = FALSE]
  module_df$rank <- seq_len(nrow(module_df))

  if (nrow(module_df) == 0) {
    stop("No module members remain after filtering: ", module_path)
  }

  module_df
}

# %% Parameters

parser <- ArgumentParser(description = "Test ranked module members for enrichment of external feature sets")
parser$add_argument("--modules_dir", type = "character", required = TRUE)
parser$add_argument("--gmt_path", type = "character", required = TRUE)
parser$add_argument("--out_dir", type = "character", required = TRUE)
parser$add_argument("--feature_id_format", type = "character", required = TRUE)
parser$add_argument("--feature_id_separator", type = "character", required = TRUE)
parser$add_argument("--duplicate_enrichment_id_policy", type = "character", required = TRUE)
parser$add_argument("--min_hits", type = "integer", required = TRUE)
parser$add_argument("--max_frac", type = "double", required = TRUE)
parser$add_argument("--use_only_important_loading", type = "character", required = TRUE)
args <- parser$parse_args()
list2env(args, envir = environment())

use_only_important_loading <- parse_bool(use_only_important_loading, "use_only_important_loading")

if (min_hits < 1L) {
  stop("min_hits must be a positive integer")
}
if (max_frac <= 0 || max_frac > 1) {
  stop("max_frac must be > 0 and <= 1")
}
if (!dir.exists(modules_dir)) {
  stop("modules_dir does not exist: ", modules_dir)
}

# %% ----- MAIN -----

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

module_paths <- list.files(
  modules_dir,
  pattern = "members\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
module_paths <- sort(unique(module_paths))
if (length(module_paths) == 0) {
  stop("No members.csv files found in modules_dir: ", modules_dir)
}

terms <- read_gmt(gmt_path)
term_ids <- vapply(terms, function(x) x$term_id, character(1))

# %% A. Run enrichment tests

result_list <- list()
row_i <- 1L

for (module_path in module_paths) {
  module_id <- basename(dirname(module_path))
  message("Testing module: ", module_id)

  query_df <- read_module_query(
    module_path = module_path,
    feature_id_format = feature_id_format,
    feature_id_separator = feature_id_separator,
    duplicate_enrichment_id_policy = duplicate_enrichment_id_policy,
    use_only_important_loading = use_only_important_loading
  )
  query <- query_df$query_id

  for (term in terms) {
    term_members <- term$members
    overlap_query_ids <- intersect(query, term_members)

    if (length(term_members) == 0) {
      mhg_res <- empty_mhg_result(
        query = query,
        term = term_members,
        status = "empty_term",
        min_hits = min_hits,
        max_frac = max_frac
      )
    } else {
      mhg_res <- test_mhg_enrichment(
        query = query,
        term = term_members,
        min_hits = min_hits,
        max_frac = max_frac
      )
    }

    result_list[[row_i]] <- cbind(
      data.frame(
        module_id = module_id,
        term_id = term$term_id,
        term_description = term$term_description,
        feature_id_format = feature_id_format,
        overlap_query_ids = paste(overlap_query_ids, collapse = ","),
        stringsAsFactors = FALSE
      ),
      mhg_res
    )
    row_i <- row_i + 1L
  }
}

# %% B. Adjust p-values and sort

results_df <- do.call(rbind, result_list)
results_df$padj <- NA_real_
valid_idx <- results_df$status == "ok" & !is.na(results_df$bonferroni_pvalue)
if (any(valid_idx)) {
  results_df$padj[valid_idx] <- p.adjust(results_df$bonferroni_pvalue[valid_idx], method = "BH")
}

results_df <- results_df[
  order(
    is.na(results_df$padj),
    results_df$padj,
    results_df$bonferroni_pvalue,
    results_df$module_id,
    results_df$term_id
  ),
  ,
  drop = FALSE
]

# %% C. Save results

write.csv(
  results_df,
  file = file.path(out_dir, "mhg_enrichment.csv"),
  row.names = FALSE,
  na = ""
)

message("Module set enrichment written")
message("  modules=", length(module_paths))
message("  terms=", length(term_ids))
message("  out_path=", normalizePath(file.path(out_dir, "mhg_enrichment.csv"), mustWork = FALSE))
