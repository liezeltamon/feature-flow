# Shared small utilities for the feature module reduction pipeline.

parse_bool <- function(x, arg_name) {
  if (is.logical(x)) return(x)
  x <- tolower(as.character(x))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(arg_name, " must be TRUE or FALSE")
}

parse_optional_string <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NULL)
  x
}

safe_p_adjust <- function(x, method = "BH") {
  out <- rep(NA_real_, length(x))
  non_na <- !is.na(x)
  if (any(non_na)) out[non_na] <- p.adjust(x[non_na], method = method)
  out
}

read_gmt <- function(gmt_path) {
  if (!file.exists(gmt_path)) {
    stop("GMT file does not exist: ", gmt_path)
  }

  lines <- readLines(gmt_path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) {
    stop("GMT file is empty: ", gmt_path)
  }

  terms <- lapply(seq_along(lines), function(i) {
    parts <- strsplit(lines[[i]], "\t", fixed = TRUE)[[1]]
    if (length(parts) < 2) {
      stop("Invalid GMT line ", i, ": expected at least term and description")
    }

    members <- character(0)
    if (length(parts) >= 3) {
      members <- unique(parts[seq.int(3L, length(parts))])
      members <- members[!is.na(members) & members != ""]
    }

    list(
      term_id = parts[[1]],
      term_description = parts[[2]],
      members = members
    )
  })

  term_ids <- vapply(terms, function(x) x$term_id, character(1))
  if (any(is.na(term_ids) | term_ids == "")) {
    stop("Blank term_id found in GMT file: ", gmt_path)
  }
  if (anyDuplicated(term_ids)) {
    dupes <- unique(term_ids[duplicated(term_ids)])
    stop(
      "Duplicated term_id values in GMT file. Examples: ",
      paste(utils::head(dupes, 5), collapse = ", ")
    )
  }

  terms
}

elapsed_sec <- function(t0) round((proc.time() - t0)[["elapsed"]], 2)
