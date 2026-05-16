# Stop with a concise validation error.
fail_contract <- function(...) {
  stop(paste0(...), call. = FALSE)
}


# Read a CSV while keeping sample/key columns as character IDs.
read_csv_preserve_sample_key <- function(
    path,
    sample_key,
    show_col_types = FALSE,
    preserve_cols = character()
) {
  raw_df <- readr::read_csv(path, show_col_types = show_col_types)

  preserve_cols <- unique(c(sample_key, preserve_cols))
  preserve_cols <- preserve_cols[!is.na(preserve_cols) & preserve_cols != ""]
  preserve_cols <- intersect(preserve_cols, colnames(raw_df))

  if (length(preserve_cols) == 0) {
    return(raw_df)
  }

  specs <- stats::setNames(
    rep(list(readr::col_character()), length(preserve_cols)),
    preserve_cols
  )
  preserve_spec <- do.call(readr::cols_only, specs)

  preserve_df <- readr::read_csv(
    path, col_types = preserve_spec, na = "", show_col_types = show_col_types
  )

  for (col in preserve_cols) {
    raw_df[[col]] <- preserve_df[[col]]
  }
  raw_df
}


# Check that the sample key exists and uniquely identifies rows.
validate_sample_key_column <- function(
    df,
    sample_key,
    table_label = "feature table",
    require_first = TRUE
) {
  if (!(sample_key %in% colnames(df))) {
    fail_contract("Column '", sample_key, "' not found in ", table_label, ".")
  }

  if (require_first && (colnames(df)[[1]] != sample_key)) {
    fail_contract(
      "Column '", sample_key, "' must be the first column in ", table_label, "."
    )
  }

  sample_ids <- df[[sample_key]]
  if (any(is.na(sample_ids) | sample_ids == "")) {
    fail_contract(
      "Column '", sample_key, "' in ", table_label,
      " contains missing or empty values."
    )
  }

  if (anyDuplicated(sample_ids) > 0) {
    fail_contract(
      "Column '", sample_key, "' in ", table_label,
      " contains duplicated values."
    )
  }

  invisible(TRUE)
}


# Keep numeric feature columns and track ignored non-feature columns.
extract_numeric_feature_table <- function(df, sample_key, table_label = "feature table") {
  validate_sample_key_column(
    df, sample_key, table_label = table_label, require_first = TRUE
  )

  is_numeric_col <- vapply(df, is.numeric, logical(1))
  numeric_cols <- setdiff(names(is_numeric_col)[is_numeric_col], sample_key)
  ignored_cols <- setdiff(colnames(df), c(sample_key, numeric_cols))

  feature_df <- df[, c(sample_key, numeric_cols), drop = FALSE]

  list(
    feature_df = feature_df,
    numeric_cols = numeric_cols,
    ignored_cols = ignored_cols
  )
}


# Load a feature table in the standard pipeline format.
read_feature_table_generic <- function(path, sample_key, table_label = NULL) {
  if (is.null(table_label)) {
    table_label <- paste0("feature table '", path, "'")
  }

  raw_df <- read_csv_preserve_sample_key(
    path,
    sample_key = sample_key,
    show_col_types = FALSE
  )
  parsed <- extract_numeric_feature_table(
    raw_df,
    sample_key = sample_key,
    table_label = table_label
  )

  list(
    feature_df = parsed$feature_df,
    numeric_cols = parsed$numeric_cols,
    ignored_cols = parsed$ignored_cols
  )
}


# Convert a wide feature table to a samples x features matrix.
feature_table_to_matrix <- function(feature_df, sample_key, table_label = "feature table") {
  validate_sample_key_column(
    feature_df,
    sample_key,
    table_label = table_label,
    require_first = TRUE
  )

  numeric_cols <- setdiff(colnames(feature_df), sample_key)
  if (length(numeric_cols) == 0) {
    fail_contract("No numeric feature columns found in ", table_label, ".")
  }

  feature_df %>%
    as.data.frame() %>%
    tibble::column_to_rownames(var = sample_key) %>%
    as.matrix()
}


# Load metadata and align it to feature-table sample order.
read_metadata_aligned <- function(
    path,
    sample_key,
    sample_ids,
    required_cols = NULL,
    table_label = NULL,
    preserve_cols = character()
) {
  if (is.null(table_label)) {
    table_label <- paste0("metadata table '", path, "'")
  }

  metadata_df <- read_csv_preserve_sample_key(
    path,
    sample_key = sample_key,
    show_col_types = FALSE,
    preserve_cols = preserve_cols
  )
  validate_sample_key_column(
    metadata_df,
    sample_key,
    table_label = table_label,
    require_first = FALSE
  )

  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, colnames(metadata_df))
    if (length(missing_cols) > 0) {
      fail_contract(
        "Missing required columns in ", table_label, ": ",
        paste(missing_cols, collapse = ", "), "."
      )
    }
  }

  metadata_ids <- metadata_df[[sample_key]]
  missing_ids <- setdiff(sample_ids, metadata_ids)
  if (length(missing_ids) > 0) {
    fail_contract(
      table_label, " is missing sample IDs present in the feature table: ",
      paste(utils::head(missing_ids, 5), collapse = ", "),
      if (length(missing_ids) > 5) " ..." else "",
      "."
    )
  }

  aligned_df <- metadata_df[match(sample_ids, metadata_ids), , drop = FALSE]
  if (!identical(aligned_df[[sample_key]], sample_ids)) {
    fail_contract(
      "Failed to align ", table_label, " to feature-table sample order."
    )
  }

  aligned_df
}
