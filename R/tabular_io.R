# Internal fast tabular, RDS, and Arrow/Parquet I/O helpers

arrow_is_available <- function() {
  requireNamespace("arrow", quietly = TRUE)
}

require_arrow_pkg <- function(context = "Parquet I/O") {
  if (!arrow_is_available()) {
    stop(
      "Package 'arrow' must be installed for ", context,
      ". Install it with install.packages('arrow').",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

is_parquet_path <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    return(FALSE)
  }
  if (grepl("\\.(parquet|pq)$", path, ignore.case = TRUE)) {
    return(TRUE)
  }
  if (dir.exists(path)) {
    files <- list.files(
      path,
      pattern = "(\\.(parquet|pq)$|^_common_metadata$|^_metadata$)",
      recursive = TRUE,
      all.files = TRUE,
      ignore.case = TRUE
    )
    return(length(files) > 0L)
  }
  FALSE
}

peek_delim_header <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else if (grepl("\\.bz2$", path, ignore.case = TRUE)) {
    bzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(con), add = TRUE)
  line1 <- readLines(con, n = 1L, warn = FALSE)
  if (!length(line1) || !nzchar(line1)) {
    return(character(0))
  }
  names(data.table::fread(text = paste0(line1, "\n"), nrows = 0L))
}

validate_split_filter <- function(split_filter) {
  if (is.null(split_filter)) {
    return(NULL)
  }
  checkmate::assert_list(split_filter, names = "unique", min.len = 1L)
  for (nm in names(split_filter)) {
    vals <- split_filter[[nm]]
    if (!is.atomic(vals) || length(vals) == 0L) {
      stop(
        "Each entry of 'split_filter' must be a non-empty atomic vector of allowed values (column '",
        nm, "').",
        call. = FALSE
      )
    }
  }
  split_filter
}

apply_split_filter_dt <- function(df, split_filter, cols = NULL) {
  if (is.null(split_filter) || length(split_filter) == 0L) {
    return(df)
  }
  use_cols <- intersect(names(split_filter), names(df))
  if (!is.null(cols)) {
    use_cols <- intersect(use_cols, cols)
  }
  if (!length(use_cols) || nrow(df) == 0L) {
    return(df)
  }
  keep <- rep(TRUE, nrow(df))
  for (nm in use_cols) {
    vals <- split_filter[[nm]]
    col_v <- df[[nm]]
    if (is.factor(col_v)) {
      col_v <- as.character(col_v)
      vals <- as.character(vals)
    }
    keep <- keep & !is.na(col_v) & (col_v %in% vals)
  }
  if (all(keep)) {
    return(df)
  }
  if (data.table::is.data.table(df)) {
    df[which(keep)]
  } else {
    df[which(keep), , drop = FALSE]
  }
}

arrow_write_parquet_impl <- function(df, path) {
  require_arrow_pkg("writing Parquet files")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(as.data.frame(df), sink = path)
  invisible(path)
}

arrow_dataset_names <- function(path) {
  require_arrow_pkg("reading Parquet schema")
  ds <- arrow::open_dataset(path)
  names(ds)
}

arrow_read_parquet_impl <- function(
  path,
  col_select = NULL,
  split_filter = NULL,
  as_data_table = TRUE
) {
  require_arrow_pkg("reading Parquet files")
  ds <- arrow::open_dataset(path)
  avail_cols <- names(ds)

  cols_use <- if (!is.null(col_select)) {
    hit <- intersect(unique(as.character(col_select)), avail_cols)
    if (length(hit) > 0L) hit else avail_cols
  } else {
    avail_cols
  }

  filter_cols <- if (!is.null(split_filter)) {
    intersect(names(split_filter), avail_cols)
  } else {
    character(0)
  }

  # Ensure columns referenced by pushed-down filters are available before projection
  scan_cols <- unique(c(cols_use, filter_cols))
  qry <- dplyr::select(ds, tidyselect::all_of(scan_cols))
  for (nm in filter_cols) {
    vals <- split_filter[[nm]]
    qry <- dplyr::filter(qry, .data[[nm]] %in% vals)
  }
  if (!identical(scan_cols, cols_use)) {
    qry <- dplyr::select(qry, tidyselect::all_of(cols_use))
  }
  out <- dplyr::collect(qry)
  if (isTRUE(as_data_table)) {
    data.table::as.data.table(out)
  } else {
    as.data.frame(out)
  }
}

medusa_write_table <- function(df, path) {
  checkmate::assert_string(path, min.chars = 1L)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (grepl("\\.(parquet|pq)$", path, ignore.case = TRUE)) {
    arrow_write_parquet_impl(df, path)
  } else if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    saveRDS(df, file = path)
  } else {
    data.table::fwrite(df, file = path, na = "NA")
  }
  invisible(path)
}

medusa_read_table <- function(
  path,
  col_select = NULL,
  split_filter = NULL,
  as_data_table = TRUE
) {
  checkmate::assert_string(path, min.chars = 1L)
  if (is_parquet_path(path)) {
    return(arrow_read_parquet_impl(
      path = path,
      col_select = col_select,
      split_filter = split_filter,
      as_data_table = as_data_table
    ))
  }

  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if (!is.data.frame(obj)) {
      stop("RDS file does not contain a data.frame/data.table: ", path, call. = FALSE)
    }
    if (!is.null(col_select)) {
      keep_c <- intersect(unique(as.character(col_select)), names(obj))
      if (length(keep_c) > 0L && length(keep_c) < ncol(obj)) {
        obj <- if (data.table::is.data.table(obj)) {
          obj[, keep_c, with = FALSE]
        } else {
          obj[, keep_c, drop = FALSE]
        }
      }
    }
    out <- if (isTRUE(as_data_table)) {
      if (data.table::is.data.table(obj)) data.table::copy(obj) else data.table::as.data.table(obj)
    } else {
      as.data.frame(obj)
    }
    return(apply_split_filter_dt(out, split_filter))
  }

  if (grepl("\\.(csv|tsv|dat|txt)(\\.(gz|bz2))?$", path, ignore.case = TRUE, perl = TRUE)) {
    if (grepl("\\.(gz|bz2)$", path, ignore.case = TRUE) &&
        requireNamespace("R.utils", quietly = TRUE)) {
      loadNamespace("R.utils")
    }
    select_arg <- NULL
    if (!is.null(col_select)) {
      avail_cols <- peek_delim_header(path)
      if (length(avail_cols) > 0L) {
        hit <- intersect(unique(as.character(col_select)), avail_cols)
        if (length(hit) > 0L && length(hit) < length(avail_cols)) {
          select_arg <- hit
        }
      }
    }
    out <- if (!is.null(select_arg)) {
      data.table::fread(path, select = select_arg, data.table = isTRUE(as_data_table))
    } else {
      data.table::fread(path, data.table = isTRUE(as_data_table))
    }
    return(apply_split_filter_dt(out, split_filter))
  }

  stop("Unable to sort out this data input: ", path, call. = FALSE)
}
