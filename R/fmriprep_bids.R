# Lightweight BIDS events and fMRIPrep confounds / subject ingestion helpers

strip_neuro_ext <- function(x) {
  out <- basename(as.character(x))
  out <- sub("\\.(nii\\.gz|nii|tsv\\.gz|tsv|csv\\.gz|csv|txt\\.gz|txt|1d\\.gz|1d|json)$", "", out, ignore.case = TRUE, perl = TRUE)
  out
}

parse_bids_entities_one <- function(path) {
  p_chr <- as.character(path)
  if (is.na(p_chr) || !nzchar(p_chr)) {
    return(list(entities = character(0L), suffix = NA_character_))
  }
  stem <- strip_neuro_ext(p_chr)
  parts <- unlist(strsplit(stem, "_", fixed = TRUE), use.names = FALSE)
  parts <- parts[nzchar(parts)]
  ent <- character(0L)
  suffix <- NA_character_
  if (length(parts) > 0L) {
    kv_mask <- grepl("^[[:alnum:]]+-[[:alnum:].-]+$", parts)
    for (tok in parts[kv_mask]) {
      dash <- regexpr("-", tok, fixed = TRUE)
      k <- tolower(substr(tok, 1L, dash - 1L))
      v <- substr(tok, dash + 1L, nchar(tok))
      ent[[k]] <- v
    }
    if (!kv_mask[length(parts)]) {
      suffix <- parts[length(parts)]
    }
  }
  dir_part <- dirname(p_chr)
  if (!identical(dir_part, ".") && nzchar(dir_part)) {
    dir_toks <- unlist(strsplit(dir_part, "[/\\\\]", perl = TRUE), use.names = FALSE)
    if (!("sub" %in% names(ent))) {
      sub_hits <- dir_toks[grepl("^sub-[[:alnum:].-]+$", dir_toks)]
      if (length(sub_hits) >= 1L) ent[["sub"]] <- sub("^sub-", "", sub_hits[length(sub_hits)])
    }
    if (!("ses" %in% names(ent))) {
      ses_hits <- dir_toks[grepl("^ses-[[:alnum:].-]+$", dir_toks)]
      if (length(ses_hits) >= 1L) ent[["ses"]] <- sub("^ses-", "", ses_hits[length(ses_hits)])
    }
  }
  list(entities = ent, suffix = suffix)
}

parse_bids_entity_table <- function(paths) {
  parsed <- lapply(paths, parse_bids_entities_one)
  get_ent <- function(key) {
    vapply(parsed, function(p) {
      if (key %in% names(p$entities)) unname(p$entities[[key]]) else NA_character_
    }, character(1L))
  }
  data.frame(
    file = as.character(paths),
    sub = get_ent("sub"),
    ses = get_ent("ses"),
    task = get_ent("task"),
    acq = get_ent("acq"),
    ce = get_ent("ce"),
    rec = get_ent("rec"),
    dir = get_ent("dir"),
    run = get_ent("run"),
    echo = get_ent("echo"),
    hemi = get_ent("hemi"),
    space = get_ent("space"),
    cohort = get_ent("cohort"),
    res = get_ent("res"),
    den = get_ent("den"),
    desc = get_ent("desc"),
    suffix = vapply(parsed, function(p) p$suffix, character(1L)),
    stringsAsFactors = FALSE
  )
}

extract_bold_run_context <- function(bold = NULL, run_names = NULL) {
  if (!is.null(run_names)) {
    checkmate::assert_character(run_names, min.len = 1L, any.missing = FALSE, min.chars = 1L, unique = TRUE)
  }
  if (is.null(bold)) {
    return(list(
      run_names = run_names,
      n_time = NULL,
      bold_files = NULL,
      tr = NULL
    ))
  }
  if (inherits(bold, "medusa_nifti_ts")) {
    b_runs <- names(bold$bold)
    orient <- if (!is.null(bold$orientation)) bold$orientation else "time_by_voxel"
    t_dim <- if (identical(orient, "voxel_by_time")) ncol else nrow
    n_t <- stats::setNames(vapply(bold$bold, t_dim, integer(1L)), b_runs)
    b_files <- if (!is.null(bold$run_info) && all(c("run", "file") %in% names(bold$run_info))) {
      stats::setNames(as.character(bold$run_info$file[match(b_runs, bold$run_info$run)]), b_runs)
    } else {
      NULL
    }
    b_tr <- if (is.numeric(bold$tr) && length(bold$tr) == 1L && is.finite(bold$tr) && bold$tr > 0) {
      as.double(bold$tr)
    } else {
      NULL
    }
    if (!is.null(run_names) && !identical(run_names, b_runs)) {
      stop("run_names does not match names(bold$bold).", call. = FALSE)
    }
    return(list(run_names = b_runs, n_time = n_t, bold_files = b_files, tr = b_tr))
  }
  if (is.matrix(bold)) {
    orient <- attr(bold, "orientation", exact = TRUE)
    n_t <- if (identical(orient, "voxel_by_time")) ncol(bold) else nrow(bold)
    r_nm <- if (!is.null(run_names)) {
      if (length(run_names) != 1L) stop("Single-matrix 'bold' requires length(run_names) == 1.", call. = FALSE)
      run_names
    } else {
      "run1"
    }
    geom <- attr(bold, "geometry", exact = TRUE)
    b_file <- if (is.list(geom) && is.character(geom$reference_file) && length(geom$reference_file) == 1L &&
                  !is.na(geom$reference_file)) {
      stats::setNames(geom$reference_file, r_nm)
    } else {
      NULL
    }
    attr_tr <- attr(bold, "tr", exact = TRUE)
    b_tr <- if (is.numeric(attr_tr) && length(attr_tr) == 1L && is.finite(attr_tr) && attr_tr > 0) {
      as.double(attr_tr)
    } else if (is.list(geom) && is.numeric(geom$header_tr) && length(geom$header_tr) == 1L &&
               is.finite(geom$header_tr) && geom$header_tr > 0) {
      as.double(geom$header_tr)
    } else {
      NULL
    }
    return(list(run_names = r_nm, n_time = stats::setNames(as.integer(n_t), r_nm), bold_files = b_file, tr = b_tr))
  }
  if (is.character(bold) && is.null(run_names)) {
    checkmate::assert_character(bold, min.len = 1L, any.missing = FALSE, min.chars = 1L, unique = TRUE)
    return(list(run_names = bold, n_time = NULL, bold_files = NULL, tr = NULL))
  }
  if (is.list(bold) && !is.data.frame(bold) && length(bold) >= 1L &&
      all(vapply(bold, is.matrix, logical(1L)))) {
    b_runs <- names(bold)
    if (is.null(b_runs) || any(!nzchar(b_runs)) || anyDuplicated(b_runs)) {
      if (!is.null(run_names)) {
        if (length(run_names) != length(bold)) {
          stop("Length of 'run_names' must match length(bold).", call. = FALSE)
        }
        b_runs <- run_names
      } else {
        b_runs <- paste0("run", seq_along(bold))
      }
    } else if (!is.null(run_names) && !identical(run_names, b_runs)) {
      stop("run_names does not match names(bold).", call. = FALSE)
    }
    n_t <- stats::setNames(vapply(bold, nrow, integer(1L)), b_runs)
    m1 <- bold[[1L]]
    attr_tr <- attr(m1, "tr", exact = TRUE)
    geom1 <- attr(m1, "geometry", exact = TRUE)
    b_tr <- if (is.numeric(attr_tr) && length(attr_tr) == 1L && is.finite(attr_tr) && attr_tr > 0) {
      as.double(attr_tr)
    } else if (is.list(geom1) && is.numeric(geom1$header_tr) && length(geom1$header_tr) == 1L &&
               is.finite(geom1$header_tr) && geom1$header_tr > 0) {
      as.double(geom1$header_tr)
    } else {
      NULL
    }
    return(list(run_names = b_runs, n_time = n_t, bold_files = NULL, tr = b_tr))
  }
  stop(
    "'bold' must be NULL, a 'medusa_nifti_ts' object, a numeric matrix, a named list of run matrices, or a character vector of run names.",
    call. = FALSE
  )
}

validate_bids_entity_compatibility <- function(in_tbl, bold_tbl, input_label) {
  for (field in c("sub", "ses", "task", "space", "cohort")) {
    b_vals <- unique(bold_tbl[[field]][!is.na(bold_tbl[[field]]) & nzchar(bold_tbl[[field]])])
    i_vals <- unique(in_tbl[[field]][!is.na(in_tbl[[field]]) & nzchar(in_tbl[[field]])])
    if (length(b_vals) > 0L && length(i_vals) > 0L && !all(i_vals %in% b_vals)) {
      stop(
        sprintf(
          "BIDS '%s-' entity mismatch between BOLD runs (%s) and %s (%s).",
          field,
          paste(b_vals, collapse = ", "),
          input_label,
          paste(i_vals, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

bids_run_signature <- function(tbl) {
  canon_r <- canonical_run_key(tbl$run)
  run_tok <- ifelse(!is.na(canon_r), canon_r, ifelse(is.na(tbl$run), "", tolower(tbl$run)))
  paste(
    ifelse(is.na(tbl$ses), "", tolower(tbl$ses)),
    ifelse(is.na(tbl$task), "", tolower(tbl$task)),
    ifelse(is.na(tbl$acq), "", tolower(tbl$acq)),
    ifelse(is.na(tbl$dir), "", tolower(tbl$dir)),
    run_tok,
    sep = "|"
  )
}

default_run_names_from_bids <- function(tbl) {
  n <- nrow(tbl)
  if (n == 0L) return(character(0L))
  has_multi_ses <- length(unique(tbl$ses[!is.na(tbl$ses) & nzchar(tbl$ses)])) > 1L
  has_run <- !is.na(tbl$run) & nzchar(tbl$run)
  if (all(has_run)) {
    r_clean <- vapply(tbl$run, function(r) {
      if (grepl("^[0-9]+$", r)) as.character(as.integer(r)) else r
    }, character(1L))
    cand <- if (has_multi_ses) {
      paste0("ses-", tbl$ses, "_run", r_clean)
    } else {
      paste0("run", r_clean)
    }
    if (!anyDuplicated(cand)) return(cand)
  }
  paste0("run", seq_len(n))
}

align_inputs_to_runs <- function(inputs, ctx, input_label = "inputs", is_paths = FALSE) {
  n_in <- length(inputs)
  in_names <- names(inputs)
  has_names <- !is.null(in_names) && any(nzchar(in_names))
  if (has_names && (any(!nzchar(in_names)) || anyDuplicated(in_names))) {
    stop("When ", input_label, " is named, all names must be nonempty and unique.", call. = FALSE)
  }

  in_tbl <- if (is_paths) parse_bids_entity_table(inputs) else NULL
  bold_tbl <- if (!is.null(ctx$bold_files)) parse_bids_entity_table(ctx$bold_files) else NULL
  if (!is.null(in_tbl) && !is.null(bold_tbl)) {
    validate_bids_entity_compatibility(in_tbl, bold_tbl, input_label)
  }

  target_runs <- ctx$run_names
  if (is.null(target_runs)) {
    resolved_runs <- if (has_names) {
      in_names
    } else if (!is.null(in_tbl)) {
      default_run_names_from_bids(in_tbl)
    } else {
      paste0("run", seq_len(n_in))
    }
    return(list(items = stats::setNames(inputs, resolved_runs), run_names = resolved_runs, bids_table = in_tbl))
  }

  n_target <- length(target_runs)
  if (n_in != n_target) {
    stop(
      sprintf(
        "Length of '%s' (%d) must match number of BOLD runs (%d: %s).",
        input_label, n_in, n_target, paste(target_runs, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  ord <- NULL
  if (has_names) {
    if (setequal(in_names, target_runs)) {
      ord <- match(target_runs, in_names)
    } else {
      c_in <- canonical_run_key(in_names)
      c_tg <- canonical_run_key(target_runs)
      if (!anyNA(c_in) && !anyNA(c_tg) && !anyDuplicated(c_in) && !anyDuplicated(c_tg) && setequal(c_in, c_tg)) {
        ord <- match(c_tg, c_in)
      } else {
        stop(
          sprintf(
            "Names of '%s' (%s) do not match BOLD run names (%s).",
            input_label, paste(in_names, collapse = ", "), paste(target_runs, collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }
  } else if (n_target == 1L) {
    ord <- 1L
  } else if (!is.null(in_tbl)) {
    if (!is.null(bold_tbl)) {
      sig_in <- bids_run_signature(in_tbl)
      sig_bd <- bids_run_signature(bold_tbl)
      if (!anyDuplicated(sig_in) && !anyDuplicated(sig_bd) && setequal(sig_in, sig_bd)) {
        ord <- match(sig_bd, sig_in)
      }
    }
    if (is.null(ord) && any(!is.na(in_tbl$run))) {
      c_in <- canonical_run_key(in_tbl$run)
      c_tg <- if (!is.null(bold_tbl) && any(!is.na(bold_tbl$run))) {
        canonical_run_key(bold_tbl$run)
      } else {
        canonical_run_key(target_runs)
      }
      if (!anyNA(c_in) && !anyNA(c_tg) && !anyDuplicated(c_in) && !anyDuplicated(c_tg) && setequal(c_in, c_tg)) {
        ord <- match(c_tg, c_in)
      } else {
        stop(
          sprintf(
            "BIDS run labels in '%s' (%s) do not match BOLD runs (%s); pass a named vector to specify run pairing explicitly.",
            input_label,
            paste(ifelse(is.na(in_tbl$run), "<none>", in_tbl$run), collapse = ", "),
            paste(target_runs, collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }
  }

  if (is.null(ord)) ord <- seq_len(n_target)
  reordered <- inputs[ord]
  names(reordered) <- target_runs
  if (!is.null(in_tbl)) {
    in_tbl <- in_tbl[ord, , drop = FALSE]
    rownames(in_tbl) <- NULL
  }
  list(items = reordered, run_names = target_runs, bids_table = in_tbl)
}

default_bids_na_strings <- function() {
  c("n/a", "N/A", "NA", ".", "NaN", "")
}

read_bids_delim <- function(path, na_strings = default_bids_na_strings()) {
  checkmate::assert_file_exists(path)
  if (R.utils::isGzipped(path, method = "extension")) {
    loadNamespace("R.utils")
  }
  is_tsv <- grepl("\\.tsv(\\.gz)?$", path, ignore.case = TRUE)
  dt <- if (is_tsv) {
    data.table::fread(path, sep = "\t", na.strings = na_strings, data.table = FALSE)
  } else {
    data.table::fread(path, na.strings = na_strings, data.table = FALSE)
  }
  as.data.frame(dt, stringsAsFactors = FALSE)
}

coerce_bids_numeric_col <- function(x, col_name, run_label) {
  if (is.logical(x) && all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  if (is.numeric(x)) {
    out <- as.double(x)
    out[is.nan(out)] <- NA_real_
    return(out)
  }
  stop(
    sprintf("Column '%s' in %s must be numeric (got %s).", col_name, run_label, class(x)[1L]),
    call. = FALSE
  )
}

parse_json_numeric_field <- function(json_path, key) {
  if (!is.character(json_path) || length(json_path) != 1L || is.na(json_path) || !file.exists(json_path)) {
    return(NA_real_)
  }
  txt <- paste(readLines(json_path, warn = FALSE), collapse = "\n")
  pat <- paste0('"', key, '"\\s*:\\s*([-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?)')
  m <- regexec(pat, txt, perl = TRUE)
  reg <- regmatches(txt, m)[[1L]]
  if (length(reg) >= 2L) as.numeric(reg[2L]) else NA_real_
}

inspect_braingnomes_provenance <- function(bold_path, use_jsonlite = TRUE) {
  empty <- list(
    is_braingnomes = FALSE,
    operations = character(0L),
    scrub_timepoints = FALSE,
    confound_regression = FALSE,
    apply_aroma = FALSE,
    temporal_filter = FALSE,
    scrub_interpolate = FALSE,
    spatial_mask = NA_character_,
    upstream_confounds = NA_character_,
    melodic_mix = NA_character_,
    aroma_metrics = NA_character_
  )
  if (!is.character(bold_path) || length(bold_path) != 1L || is.na(bold_path) || !nzchar(bold_path)) {
    return(empty)
  }
  json_path <- if (grepl("\\.json$", bold_path, ignore.case = TRUE)) {
    bold_path
  } else {
    sub("\\.nii(\\.gz)?$", ".json", bold_path, ignore.case = TRUE)
  }
  if (!file.exists(json_path)) return(empty)

  txt <- paste(readLines(json_path, warn = FALSE), collapse = "\n")
  if (!grepl('"BrainGnomes"\\s*:', txt, perl = TRUE)) return(empty)

  active_ops <- NULL
  spatial_mask <- NA_character_
  upstream_conf <- NA_character_
  melodic_mix <- NA_character_
  aroma_metrics <- NA_character_

  if (isTRUE(use_jsonlite) && requireNamespace("jsonlite", quietly = TRUE)) {
    parsed_json <- tryCatch(jsonlite::fromJSON(json_path, simplifyVector = FALSE), error = function(e) NULL)
    if (is.list(parsed_json) && is.list(parsed_json$BrainGnomes)) {
      bg <- parsed_json$BrainGnomes
      if (is.list(bg$Operations)) {
        op_keep <- vapply(bg$Operations, function(op) {
          is.list(op) && is.character(op$Name) && length(op$Name) == 1L &&
            is.character(op$Status) && length(op$Status) == 1L &&
            op$Status %in% c("executed", "reused")
        }, logical(1L))
        active_ops <- unique(vapply(bg$Operations[op_keep], function(op) op$Name, character(1L)))
      } else {
        active_ops <- character(0L)
      }
      if (is.list(bg$SourceFiles)) {
        for (sf in bg$SourceFiles) {
          if (is.list(sf) && is.character(sf$role) && is.character(sf$path) && length(sf$path) == 1L) {
            if (identical(sf$role, "spatial_mask") && is.na(spatial_mask)) spatial_mask <- sf$path
            if (identical(sf$role, "confounds") && is.na(upstream_conf)) upstream_conf <- sf$path
            if (identical(sf$role, "melodic_mix") && is.na(melodic_mix)) melodic_mix <- sf$path
            if (identical(sf$role, "aroma_metrics") && is.na(aroma_metrics)) aroma_metrics <- sf$path
          }
        }
      }
    }
  }

  if (is.null(active_ops)) {
    # Pure base-R fallback: strip nested SourceMetadata block before matching Operations/SourceFiles
    txt_self <- sub('(?s)"SourceMetadata"\\s*:.*$', "", txt, perl = TRUE)
    op_pat <- '\\{[[:space:]]*"Name"[[:space:]]*:[[:space:]]*"([^"]+)"[[:space:]]*,[[:space:]]*"Status"[[:space:]]*:[[:space:]]*"([^"]+)"'
    m_all <- gregexpr(op_pat, txt_self, perl = TRUE)
    chunks <- regmatches(txt_self, m_all)[[1L]]
    active_ops <- character(0L)
    for (ch in chunks) {
      sub_m <- regmatches(ch, regexec(op_pat, ch, perl = TRUE))[[1L]]
      if (length(sub_m) >= 3L && sub_m[3L] %in% c("executed", "reused")) {
        active_ops <- c(active_ops, sub_m[2L])
      }
    }
    active_ops <- unique(active_ops)

    extract_sf_path <- function(role_name) {
      sf_pat <- paste0('"role"\\s*:\\s*"', role_name, '"[^}]*?"path"\\s*:\\s*"([^"]+)"')
      sm <- regmatches(txt_self, regexec(sf_pat, txt_self, perl = TRUE))[[1L]]
      if (length(sm) >= 2L) return(sm[2L])
      sf_pat2 <- paste0('"path"\\s*:\\s*"([^"]+)"[^}]*?"role"\\s*:\\s*"', role_name, '"')
      sm2 <- regmatches(txt_self, regexec(sf_pat2, txt_self, perl = TRUE))[[1L]]
      if (length(sm2) >= 2L) return(sm2[2L])
      NA_character_
    }
    spatial_mask <- extract_sf_path("spatial_mask")
    upstream_conf <- extract_sf_path("confounds")
    melodic_mix <- extract_sf_path("melodic_mix")
    aroma_metrics <- extract_sf_path("aroma_metrics")
  }

  list(
    is_braingnomes = TRUE,
    operations = active_ops,
    scrub_timepoints = "scrub_timepoints" %in% active_ops,
    confound_regression = "confound_regression" %in% active_ops,
    apply_aroma = "apply_aroma" %in% active_ops,
    temporal_filter = "temporal_filter" %in% active_ops,
    scrub_interpolate = "scrub_interpolate" %in% active_ops,
    spatial_mask = spatial_mask,
    upstream_confounds = upstream_conf,
    melodic_mix = melodic_mix,
    aroma_metrics = aroma_metrics
  )
}

assert_no_scrubbed_bold_provenance <- function(bold_files) {
  if (is.null(bold_files) || !length(bold_files)) return(invisible(list()))
  prov_list <- lapply(bold_files, inspect_braingnomes_provenance)
  names(prov_list) <- names(bold_files)
  for (i in seq_along(bold_files)) {
    if (isTRUE(prov_list[[i]]$scrub_timepoints)) {
      stop(
        sprintf(
          "BOLD derivative '%s' records BrainGnomes 'scrub_timepoints' (physical volume deletion), which breaks MEDuSA's continuous acquisition clock. Use un-scrubbed postprocessed outputs with '*_censor.1D' / valid_rows censoring instead.",
          basename(bold_files[i])
        ),
        call. = FALSE
      )
    }
  }
  invisible(prov_list)
}

#' Read and normalize BIDS event files for MEDuSA HRF, deconvolution, and alignment workflows
#'
#' Reads one or more BIDS \code{*_events.tsv} files (or normalizes in-memory BIDS
#' \code{data.frame}s) into MEDuSA's event table contract: finite numeric
#' \code{onset} and non-negative \code{duration} (in seconds), non-empty
#' \code{condition} labels (mapped by default from BIDS \code{trial_type}),
#' run-matched \code{run} labels in \code{bold} run order, and unique within-run
#' \code{trial} indices. Standard BIDS missing tokens (\code{"n/a"}) are parsed
#' automatically without coercing numeric behavioral columns to \code{character}.
#'
#' @param events Character vector of BIDS event file paths (\code{.tsv},
#'   \code{.tsv.gz}, \code{.csv}, \code{.csv.gz}), a single \code{data.frame},
#'   or a list of per-run \code{data.frame}s.
#' @param bold Optional \code{"medusa_nifti_ts"} object, named list of run BOLD
#'   matrices, single BOLD matrix, or character vector of target run names used
#'   to validate BIDS \code{sub-}/\code{ses-}/\code{task-} entities and order
#'   event runs identically to BOLD runs.
#' @param run_names Optional explicit character vector of unique run names when
#'   \code{bold} is \code{NULL}.
#' @param condition_col Column name in the BIDS event table supplying condition
#'   labels when \code{condition} is not already present (or when different from
#'   \code{"condition"}). Default \code{"trial_type"}.
#' @param condition_map Optional named character vector recoding raw condition
#'   values (e.g., \code{c(win_fb = "feedback", loss_fb = "feedback")}).
#'   Unmatched condition labels are preserved unchanged.
#' @param default_condition Optional scalar string assigned to \code{condition}
#'   only when neither \code{"condition"} nor \code{condition_col} exists in a
#'   single-condition BIDS event table (e.g., \code{"event"}).
#' @param exclude_conditions Optional character vector of raw condition labels
#'   (prior to \code{condition_map}) to exclude from the returned table.
#' @param na_duration Replacement duration (in seconds, \code{>= 0}) used when
#'   \code{duration} is missing or \code{"n/a"} (default \code{0} for impulse
#'   events). Set \code{NULL} to error on missing durations.
#' @param drop_missing_onset Logical; if \code{TRUE}, rows with missing or
#'   non-finite \code{onset} are dropped instead of raising an error. Default
#'   \code{FALSE}.
#' @param drop_missing_condition Logical; if \code{TRUE}, rows with missing or
#'   empty \code{condition} labels are dropped instead of raising an error.
#'   Default \code{FALSE}.
#' @param time_shift Numeric scalar or run-ordered numeric vector added to
#'   \code{onset} (in seconds; default \code{0}).
#' @param trial_col Name of the within-run unique trial index column populated
#'   when not already present in \code{events}. Default \code{"trial"}.
#' @param run_col Name of the run identifier column populated on returned event
#'   tables. Default \code{"run"}.
#' @param sort_by_onset Logical; if \code{TRUE}, sorts events chronologically by
#'   \code{onset} within each run while recording original 1-based file row
#'   numbers in \code{bids_event_row}. Default \code{FALSE} preserves exact file
#'   row order for 1-to-1 behavioral joins.
#' @param na_strings Character vector of missing-value tokens passed to
#'   \code{data.table::fread}. Default \code{c("n/a", "N/A", "NA", ".", "NaN", "")}.
#' @param output Return container format: \code{"data.frame"} (default; single
#'   combined multi-run \code{data.frame} ready to pass directly to
#'   \code{\link{estimate_crossrun_hrf}}, \code{\link{decon_epoch_regression}},
#'   and \code{\link{align_medusa_ts}}) or \code{"list"} (run-named list of
#'   per-run \code{data.frame}s).
#'
#' @return Either a combined multi-run \code{data.frame} (\code{output = "data.frame"})
#'   or a run-named list of \code{data.frame}s (\code{output = "list"}).
#'
#' @seealso \code{\link{read_fmriprep_confounds}}, \code{\link{read_fmriprep_subject}},
#'   \code{\link{read_nifti_ts}}, \code{\link{estimate_crossrun_hrf}},
#'   \code{\link{decon_epoch_regression}}
#'
#' @examples
#' f1 <- tempfile("sub-01_task-mid_run-01_", fileext = "_events.tsv")
#' f2 <- tempfile("sub-01_task-mid_run-02_", fileext = "_events.tsv")
#' writeLines(c(
#'   "onset\tduration\ttrial_type\tresponse_time\trpe",
#'   "2.0\t0.5\tcue\tn/a\tn/a",
#'   "6.0\t0.8\tfeedback\t0.52\t0.75"
#' ), f1)
#' writeLines(c(
#'   "onset\tduration\ttrial_type\tresponse_time\trpe",
#'   "3.0\t0.5\tcue\tn/a\tn/a",
#'   "7.5\tn/a\tfeedback\t0.61\t-0.40"
#' ), f2)
#'
#' ev <- read_bids_events(c(f1, f2))
#' unlink(c(f1, f2))
#'
#' @export
read_bids_events <- function(
  events,
  bold = NULL,
  run_names = NULL,
  condition_col = "trial_type",
  condition_map = NULL,
  default_condition = NULL,
  exclude_conditions = NULL,
  na_duration = 0,
  drop_missing_onset = FALSE,
  drop_missing_condition = FALSE,
  time_shift = 0,
  trial_col = "trial",
  run_col = "run",
  sort_by_onset = FALSE,
  na_strings = c("n/a", "N/A", "NA", ".", "NaN", ""),
  output = c("data.frame", "list")
) {
  output <- match.arg(output)
  checkmate::assert_string(condition_col, min.chars = 1L)
  checkmate::assert_string(default_condition, min.chars = 1L, null.ok = TRUE)
  checkmate::assert_character(exclude_conditions, min.chars = 1L, any.missing = FALSE, null.ok = TRUE)
  if (!is.null(condition_map)) {
    checkmate::assert_character(condition_map, min.len = 1L, any.missing = FALSE, min.chars = 1L, names = "unique")
  }
  if (!is.null(na_duration)) {
    checkmate::assert_number(na_duration, lower = 0, finite = TRUE)
  }
  checkmate::assert_flag(drop_missing_onset)
  checkmate::assert_flag(drop_missing_condition)
  checkmate::assert_flag(sort_by_onset)
  checkmate::assert_string(trial_col, min.chars = 1L)
  checkmate::assert_string(run_col, min.chars = 1L)
  checkmate::assert_character(na_strings, any.missing = FALSE)

  ctx <- extract_bold_run_context(bold = bold, run_names = run_names)

  bids_tbl <- NULL
  if (is.character(events)) {
    checkmate::assert_character(events, min.len = 1L, any.missing = FALSE, min.chars = 1L)
    aligned <- align_inputs_to_runs(events, ctx, input_label = "events", is_paths = TRUE)
    resolved_runs <- aligned$run_names
    bids_tbl <- aligned$bids_table
    raw_list <- lapply(aligned$items, read_bids_delim, na_strings = na_strings)
    names(raw_list) <- resolved_runs
  } else if (is.data.frame(events)) {
    if (!is.null(ctx$run_names)) {
      raw_list <- split_multirun_events(events, ctx$run_names, run_col = run_col, allow_empty_runs = TRUE)
      resolved_runs <- ctx$run_names
    } else if (run_col %in% names(events)) {
      r_vals <- as.character(events[[run_col]])
      u_runs <- unique(r_vals[!is.na(r_vals) & nzchar(r_vals)])
      if (!length(u_runs)) u_runs <- "run1"
      raw_list <- split_multirun_events(events, u_runs, run_col = run_col, allow_empty_runs = TRUE)
      resolved_runs <- u_runs
    } else {
      resolved_runs <- "run1"
      raw_list <- list(run1 = as.data.frame(events, stringsAsFactors = FALSE))
    }
  } else if (is.list(events)) {
    checkmate::assert_list(events, types = "data.frame", min.len = 1L)
    aligned <- align_inputs_to_runs(events, ctx, input_label = "events", is_paths = FALSE)
    resolved_runs <- aligned$run_names
    raw_list <- lapply(aligned$items, as.data.frame, stringsAsFactors = FALSE)
  } else {
    stop("'events' must be a character vector of file paths, a data.frame, or a list of data.frames.", call. = FALSE)
  }

  n_runs <- length(resolved_runs)
  checkmate::assert_numeric(time_shift, any.missing = FALSE, finite = TRUE)
  if (length(time_shift) == 1L) {
    shifts <- stats::setNames(rep(as.numeric(time_shift), n_runs), resolved_runs)
  } else if (length(time_shift) == n_runs) {
    if (!is.null(names(time_shift)) && any(nzchar(names(time_shift)))) {
      if (!setequal(names(time_shift), resolved_runs)) {
        stop("Named 'time_shift' must match resolved run names.", call. = FALSE)
      }
      shifts <- as.numeric(time_shift[resolved_runs])
    } else {
      shifts <- as.numeric(time_shift)
    }
    names(shifts) <- resolved_runs
  } else {
    stop("'time_shift' must be scalar or have length equal to the number of runs.", call. = FALSE)
  }

  norm_list <- stats::setNames(vector("list", n_runs), resolved_runs)
  for (r in seq_len(n_runs)) {
    rname <- resolved_runs[r]
    df <- as.data.frame(raw_list[[r]], stringsAsFactors = FALSE)
    run_label <- sprintf("events[['%s']]", rname)

    if (!("onset" %in% names(df))) {
      stop("Missing required BIDS 'onset' column in ", run_label, ".", call. = FALSE)
    }

    # Resolve raw condition vector before row filtering
    use_explicit_col <- condition_col %in% names(df) &&
      (!("condition" %in% names(df)) || !identical(condition_col, "trial_type"))
    if (use_explicit_col) {
      cond_raw <- as.character(df[[condition_col]])
    } else if ("condition" %in% names(df)) {
      cond_raw <- as.character(df$condition)
    } else if (condition_col %in% names(df)) {
      cond_raw <- as.character(df[[condition_col]])
    } else if (!is.null(default_condition)) {
      cond_raw <- rep(default_condition, nrow(df))
    } else {
      stop(
        sprintf(
          "Neither 'condition' nor '%s' found in %s; specify 'condition_col' or 'default_condition'.",
          condition_col, run_label
        ),
        call. = FALSE
      )
    }

    df$.medusa_orig_row <- seq_len(nrow(df))

    # Optional condition exclusion prior to mapping
    if (!is.null(exclude_conditions) && nrow(df) > 0L) {
      keep_ex <- is.na(cond_raw) | !(cond_raw %in% exclude_conditions)
      df <- df[keep_ex, , drop = FALSE]
      cond_raw <- cond_raw[keep_ex]
    }

    onset_num <- coerce_bids_numeric_col(df$onset, "onset", run_label)
    bad_onset <- !is.finite(onset_num)
    if (any(bad_onset)) {
      if (isTRUE(drop_missing_onset)) {
        df <- df[!bad_onset, , drop = FALSE]
        onset_num <- onset_num[!bad_onset]
        cond_raw <- cond_raw[!bad_onset]
      } else {
        stop(
          sprintf(
            "Non-finite or missing 'onset' at row(s) %s in %s; set drop_missing_onset = TRUE to omit them.",
            paste(utils::head(df$.medusa_orig_row[bad_onset], 6L), collapse = ", "),
            run_label
          ),
          call. = FALSE
        )
      }
    }
    df$onset <- onset_num + shifts[r]

    # Condition mapping and validation
    if (!is.null(condition_map) && length(cond_raw) > 0L) {
      hit_map <- !is.na(cond_raw) & (cond_raw %in% names(condition_map))
      cond_raw[hit_map] <- unname(condition_map[cond_raw[hit_map]])
    }
    bad_cond <- is.na(cond_raw) | !nzchar(trimws(cond_raw))
    if (any(bad_cond)) {
      if (isTRUE(drop_missing_condition)) {
        df <- df[!bad_cond, , drop = FALSE]
        cond_raw <- cond_raw[!bad_cond]
      } else {
        stop(
          sprintf(
            "Missing or empty condition label at row(s) %s in %s; set drop_missing_condition = TRUE or exclude_conditions.",
            paste(utils::head(df$.medusa_orig_row[bad_cond], 6L), collapse = ", "),
            run_label
          ),
          call. = FALSE
        )
      }
    }
    df$condition <- as.character(cond_raw)

    # Duration normalization
    if ("duration" %in% names(df)) {
      dur_num <- coerce_bids_numeric_col(df$duration, "duration", run_label)
      if (anyNA(dur_num)) {
        if (is.null(na_duration)) {
          stop("Missing 'duration' values in ", run_label, " (na_duration = NULL).", call. = FALSE)
        }
        dur_num[is.na(dur_num)] <- na_duration
      }
      if (any(!is.finite(dur_num) | dur_num < 0)) {
        stop("Event 'duration' must be finite and >= 0 in ", run_label, ".", call. = FALSE)
      }
      df$duration <- dur_num
    } else {
      fill_dur <- if (is.null(na_duration)) 0 else as.numeric(na_duration)
      df$duration <- rep(fill_dur, nrow(df))
    }

    # Optional chronological sorting
    if (isTRUE(sort_by_onset) && nrow(df) > 1L && is.unsorted(df$onset)) {
      if (!("bids_event_row" %in% names(df))) {
        df$bids_event_row <- df$.medusa_orig_row
      }
      df <- df[order(df$onset, df$.medusa_orig_row), , drop = FALSE]
    }
    df$.medusa_orig_row <- NULL

    # Attach BIDS entities (id, session, task) when parsed from filenames and not already present
    if (!is.null(bids_tbl)) {
      if (!("id" %in% names(df)) && !is.na(bids_tbl$sub[r]) && nzchar(bids_tbl$sub[r])) {
        df$id <- rep(bids_tbl$sub[r], nrow(df))
      }
      if (!("session" %in% names(df)) && !is.na(bids_tbl$ses[r]) && nzchar(bids_tbl$ses[r])) {
        df$session <- rep(bids_tbl$ses[r], nrow(df))
      }
      if (!("task" %in% names(df)) && !is.na(bids_tbl$task[r]) && nzchar(bids_tbl$task[r])) {
        df$task <- rep(bids_tbl$task[r], nrow(df))
      }
    }

    df[[run_col]] <- rep(rname, nrow(df))
    if (!(trial_col %in% names(df))) {
      df[[trial_col]] <- seq_len(nrow(df))
    }

    # Reorder standard leading columns cleanly while preserving all extra behavioral columns
    lead_cols <- intersect(c("id", "session", "task", run_col, trial_col, "onset", "duration", "condition"), names(df))
    other_cols <- setdiff(names(df), lead_cols)
    df <- df[, c(lead_cols, other_cols), drop = FALSE]
    rownames(df) <- NULL
    norm_list[[r]] <- df
  }

  if (identical(output, "list")) {
    return(norm_list)
  }

  # Harmonize any all-NA logical columns across runs prior to stacking so rbindlist preserves types
  all_col_names <- unique(unlist(lapply(norm_list, names), use.names = FALSE))
  for (cname in all_col_names) {
    non_na_types <- unique(unlist(lapply(norm_list, function(d) {
      if (!(cname %in% names(d))) return(character(0L))
      v <- d[[cname]]
      if (is.logical(v) && all(is.na(v))) character(0L) else typeof(v)
    }), use.names = FALSE))
    if (length(non_na_types) == 1L && non_na_types %in% c("double", "integer", "character")) {
      for (r in seq_len(n_runs)) {
        if (cname %in% names(norm_list[[r]])) {
          v <- norm_list[[r]][[cname]]
          if (is.logical(v) && all(is.na(v))) {
            storage.mode(norm_list[[r]][[cname]]) <- non_na_types
          }
        }
      }
    }
  }

  has_empty_run <- any(vapply(norm_list, nrow, integer(1L)) == 0L)
  out_df <- as.data.frame(
    data.table::rbindlist(norm_list, use.names = TRUE, fill = TRUE),
    stringsAsFactors = FALSE
  )
  if (has_empty_run && (run_col %in% names(out_df))) {
    out_df[[run_col]] <- factor(as.character(out_df[[run_col]]), levels = resolved_runs)
  }
  rownames(out_df) <- NULL
  out_df
}

confound_strategy_columns <- function(strategy) {
  base_6p <- c("trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z")
  deriv_6p <- paste0(base_6p, "_derivative1")
  pow2_6p <- paste0(base_6p, "_power2")
  dpow2_6p <- paste0(base_6p, "_derivative1_power2")
  csf_wm_base <- c("csf", "white_matter")
  csf_wm_8p <- c(
    csf_wm_base,
    paste0(csf_wm_base, "_derivative1"),
    paste0(csf_wm_base, "_power2"),
    paste0(csf_wm_base, "_derivative1_power2")
  )
  gsr_4p <- c(
    "global_signal",
    "global_signal_derivative1",
    "global_signal_power2",
    "global_signal_derivative1_power2"
  )
  switch(
    tolower(strategy),
    "6p" = base_6p,
    "12p" = c(base_6p, deriv_6p),
    "24p" = c(base_6p, deriv_6p, pow2_6p, dpow2_6p),
    "27p" = c(base_6p, deriv_6p, pow2_6p, dpow2_6p, "csf", "white_matter", "global_signal"),
    "csf_wm" = csf_wm_base,
    "csf_wm_8p" = csf_wm_8p,
    "24p_csf_wm" = ,
    "32p" = c(base_6p, deriv_6p, pow2_6p, dpow2_6p, csf_wm_8p),
    "36p" = c(base_6p, deriv_6p, pow2_6p, dpow2_6p, csf_wm_8p, gsr_4p),
    "compcor" = c(base_6p, deriv_6p),
    "custom" = character(0L),
    "none" = character(0L),
    NULL
  )
}

expand_confound_column_specs <- function(specs, avail_names) {
  if (is.null(specs) || !length(specs)) return(character(0L))
  out <- character(0L)
  for (sp in specs) {
    sc <- confound_strategy_columns(sp)
    if (!is.null(sc) && tolower(sp) %in% c("6p", "12p", "24p", "27p", "csf_wm", "csf_wm_8p", "24p_csf_wm", "32p", "36p")) {
      out <- c(out, sc)
      next
    }
    if (grepl("^[^<]*<[^>]+>.*$", sp)) {
      pre <- sub("^([^<]*)<[^>]+>.*$", "\\1", sp)
      idx_str <- sub("^[^<]*<([^>]+)>.*$", "\\1", sp)
      post <- sub("^[^<]*<[^>]+>(.*)$", "\\1", sp)
      parts <- trimws(unlist(strsplit(idx_str, ",", fixed = TRUE), use.names = FALSE))
      if (!length(parts) || any(!nzchar(parts))) {
        stop(sprintf("Invalid column range specification '%s'.", sp), call. = FALSE)
      }
      idx_vals <- integer(0L)
      for (pt in parts) {
        if (grepl("^[0-9]+$", pt)) {
          idx_vals <- c(idx_vals, as.integer(pt))
        } else if (grepl("^[0-9]+-[0-9]+$", pt)) {
          rng <- as.integer(unlist(strsplit(pt, "-", fixed = TRUE), use.names = FALSE))
          idx_vals <- c(idx_vals, seq.int(rng[1L], rng[2L]))
        } else {
          stop(
            sprintf("Invalid column range specification '%s'; expected integer indices like 'a_comp_cor_<0-4>'.", sp),
            call. = FALSE
          )
        }
      }
      expanded_i <- vapply(idx_vals, function(k) {
        cand_pad <- paste0(pre, sprintf("%02d", k), post)
        cand_raw <- paste0(pre, k, post)
        if (!(cand_pad %in% avail_names) && (cand_raw %in% avail_names)) cand_raw else cand_pad
      }, character(1L))
      out <- c(out, expanded_i)
    } else {
      out <- c(out, sp)
    }
  }
  unique(out)
}

can_resolve_confound_column <- function(col_name, avail_names) {
  if (col_name %in% avail_names) return(TRUE)
  if (grepl("_derivative1_power2$", col_name)) {
    base_col <- sub("_derivative1_power2$", "", col_name)
    return((base_col %in% avail_names) || (paste0(base_col, "_derivative1") %in% avail_names))
  }
  if (grepl("_(derivative1|power2)$", col_name)) {
    base_col <- sub("_(derivative1|power2)$", "", col_name)
    return(base_col %in% avail_names)
  }
  FALSE
}

expand_fmriprep_column <- function(df, col_name, run_label) {
  if (col_name %in% names(df)) {
    return(coerce_bids_numeric_col(df[[col_name]], col_name, run_label))
  }
  # Synthesize standard Friston-24 / 36p derivative and quadratic expansions when base column exists
  if (grepl("_derivative1_power2$", col_name)) {
    base_col <- sub("_derivative1_power2$", "", col_name)
    deriv_col <- paste0(base_col, "_derivative1")
    if (deriv_col %in% names(df)) {
      d <- coerce_bids_numeric_col(df[[deriv_col]], deriv_col, run_label)
      return(d^2)
    }
    if (base_col %in% names(df)) {
      b <- coerce_bids_numeric_col(df[[base_col]], base_col, run_label)
      d <- if (length(b) == 0L) numeric(0L) else c(0, diff(b))
      return(d^2)
    }
  } else if (grepl("_derivative1$", col_name)) {
    base_col <- sub("_derivative1$", "", col_name)
    if (base_col %in% names(df)) {
      b <- coerce_bids_numeric_col(df[[base_col]], base_col, run_label)
      return(if (length(b) == 0L) numeric(0L) else c(0, diff(b)))
    }
  } else if (grepl("_power2$", col_name)) {
    base_col <- sub("_power2$", "", col_name)
    if (base_col %in% names(df)) {
      b <- coerce_bids_numeric_col(df[[base_col]], base_col, run_label)
      return(b^2)
    }
  }
  stop(sprintf("Required confound column '%s' not found in %s.", col_name, run_label), call. = FALSE)
}

design_iir_notch <- function(f0, bandwidth, fs) {
  Q <- f0 / bandwidth
  w0 <- 2 * pi * (f0 / fs)
  cw <- cos(w0)
  sw <- sin(w0)
  alpha <- sw / (2 * Q)
  b <- c(1, -2 * cw, 1) / (1 + alpha)
  a <- c(1, -2 * cw / (1 + alpha), (1 - alpha) / (1 + alpha))
  list(b = as.double(b), a = as.double(a))
}

design_butterworth_lowpass <- function(order_n, cutoff_hz, fs) {
  nyq <- fs / 2
  wn <- cutoff_hz / nyq
  v <- tan(pi * wn / 2)
  k <- seq_len(order_n)
  s_poles <- v * exp(1i * pi * (2 * k + order_n - 1) / (2 * order_n))
  z_poles <- (1 + s_poles) / (1 - s_poles)
  poly_a <- c(1 + 0i)
  for (zp in z_poles) {
    poly_a <- c(poly_a, 0i) - zp * c(0i, poly_a)
  }
  a <- Re(poly_a)
  b_unscaled <- choose(order_n, 0:order_n)
  b <- b_unscaled * (sum(a) / sum(b_unscaled))
  list(b = as.double(b / a[1L]), a = as.double(a / a[1L]))
}

lfilter_zi_r <- function(b, a) {
  k <- max(length(b), length(a))
  if (length(b) < k) b <- c(b, rep(0, k - length(b)))
  if (length(a) < k) a <- c(a, rep(0, k - length(a)))
  m <- k - 1L
  if (m <= 0L) return(numeric(0L))
  A_mat <- matrix(0, nrow = m, ncol = m)
  A_mat[, 1L] <- -a[2:k]
  if (m > 1L) {
    for (j in seq_len(m - 1L)) A_mat[j, j + 1L] <- 1
  }
  rhs <- b[2:k] - a[2:k] * b[1L]
  as.double(solve(diag(m) - A_mat, rhs))
}

lfilter_df2t_r <- function(x, b, a, zi = NULL) {
  k <- max(length(b), length(a))
  if (length(b) < k) b <- c(b, rep(0, k - length(b)))
  if (length(a) < k) a <- c(a, rep(0, k - length(a)))
  m <- k - 1L
  n <- length(x)
  if (n == 0L) return(numeric(0L))
  if (m == 0L) return(as.double(b[1L] * x))
  z <- if (is.null(zi)) numeric(m) else as.double(zi)
  y <- numeric(n)
  b1 <- b[1L]
  b_tail <- b[2:k]
  a_tail <- a[2:k]
  for (t in seq_len(n)) {
    xt <- x[t]
    yt <- b1 * xt + z[1L]
    y[t] <- yt
    if (m > 1L) {
      z[1:(m - 1L)] <- b_tail[1:(m - 1L)] * xt - a_tail[1:(m - 1L)] * yt + z[2:m]
    }
    z[m] <- b_tail[m] * xt - a_tail[m] * yt
  }
  y
}

filtfilt_1d <- function(x, b, a, padlen = NULL, padtype = c("constant", "odd", "even", "zero"), use_zi = TRUE) {
  padtype <- match.arg(padtype)
  x <- as.double(x)
  n <- length(x)
  b <- as.double(b) / a[1L]
  a <- as.double(a) / a[1L]
  k <- max(length(b), length(a))
  if (n <= 1L || k <= 1L) return(b[1L]^2 * x)
  if (is.null(padlen) || padlen < 0L) {
    padlen <- 3L * k
  }
  padlen <- as.integer(min(padlen, n - 1L))
  if (padlen > 0L) {
    left_pad <- switch(
      padtype,
      constant = rep(x[1L], padlen),
      odd = 2 * x[1L] - x[(padlen + 1L):2L],
      even = x[(padlen + 1L):2L],
      zero = numeric(padlen)
    )
    right_pad <- switch(
      padtype,
      constant = rep(x[n], padlen),
      odd = 2 * x[n] - x[(n - 1L):(n - padlen)],
      even = x[(n - 1L):(n - padlen)],
      zero = numeric(padlen)
    )
    x_ext <- c(left_pad, x, right_pad)
  } else {
    x_ext <- x
  }
  zi <- if (isTRUE(use_zi)) lfilter_zi_r(b, a) else NULL
  y_fwd <- lfilter_df2t_r(x_ext, b, a, zi = if (!is.null(zi)) zi * x_ext[1L] else NULL)
  y_rev <- rev(y_fwd)
  y_bwd <- lfilter_df2t_r(y_rev, b, a, zi = if (!is.null(zi)) zi * y_rev[1L] else NULL)
  y_out <- rev(y_bwd)
  if (padlen > 0L) y_out[(padlen + 1L):(padlen + n)] else y_out
}

infer_tr_from_confound_sidecar <- function(conf_path, bold_file = NULL) {
  if (is.character(bold_file) && length(bold_file) == 1L && !is.na(bold_file) && file.exists(bold_file)) {
    json_b <- sub("\\.nii(\\.gz)?$", ".json", bold_file, ignore.case = TRUE)
    j_tr <- parse_json_numeric_field(json_b, "RepetitionTime")
    if (is.finite(j_tr) && j_tr > 0) return(as.double(j_tr))
    geom <- tryCatch(resolve_nifti_spatial_geometry(bold_file, label = basename(bold_file)), error = function(e) NULL)
    if (is.list(geom) && is.finite(geom$header_tr) && geom$header_tr > 0) return(as.double(geom$header_tr))
  }
  if (!is.character(conf_path) || length(conf_path) != 1L || is.na(conf_path) || !file.exists(conf_path)) {
    return(NULL)
  }
  d <- dirname(conf_path)
  conf_row <- parse_bids_entity_table(conf_path)
  json_cands <- list.files(d, pattern = "_bold\\.json$", full.names = TRUE, ignore.case = TRUE)
  if (length(json_cands)) {
    j_tbl <- parse_bids_entity_table(json_cands)
    # Match on non-spatial entities (sub, ses, task, run, acq, dir) since confounds omit space-
    hit <- find_matching_bids_sidecar(conf_row, j_tbl, require_space = FALSE, allow_optional_space = FALSE)
    if (is.na(hit)) {
      # When multiple template spaces exist for the same run, RepetitionTime is identical across spaces
      same_run <- rep(TRUE, nrow(j_tbl))
      for (f in c("sub", "ses", "task", "acq", "ce", "rec", "dir", "echo")) {
        cv <- conf_row[[f]]
        jv <- j_tbl[[f]]
        if (!is.na(cv) && nzchar(cv)) {
          same_run <- same_run & !is.na(jv) & tolower(jv) == tolower(cv)
        } else {
          same_run <- same_run & (is.na(jv) | !nzchar(jv))
        }
      }
      cr <- canonical_run_key(conf_row$run)
      jr <- canonical_run_key(j_tbl$run)
      same_run <- if (!is.na(cr)) same_run & !is.na(jr) & jr == cr else same_run & (is.na(jr) | !nzchar(jr))
      if (any(same_run)) hit <- j_tbl$file[which(same_run)[1L]]
    }
    if (!is.na(hit)) {
      j_tr <- parse_json_numeric_field(hit, "RepetitionTime")
      if (is.finite(j_tr) && j_tr > 0) return(as.double(j_tr))
    }
  }
  nii_cands <- list.files(d, pattern = "_bold\\.nii(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  if (length(nii_cands)) {
    n_tbl <- parse_bids_entity_table(nii_cands)
    same_run <- rep(TRUE, nrow(n_tbl))
    for (f in c("sub", "ses", "task", "acq", "ce", "rec", "dir", "echo")) {
      cv <- conf_row[[f]]
      nv <- n_tbl[[f]]
      if (!is.na(cv) && nzchar(cv)) {
        same_run <- same_run & !is.na(nv) & tolower(nv) == tolower(cv)
      } else {
        same_run <- same_run & (is.na(nv) | !nzchar(nv))
      }
    }
    cr <- canonical_run_key(conf_row$run)
    nr <- canonical_run_key(n_tbl$run)
    same_run <- if (!is.na(cr)) same_run & !is.na(nr) & nr == cr else same_run & (is.na(nr) | !nzchar(nr))
    if (any(same_run)) {
      geom <- tryCatch(resolve_nifti_spatial_geometry(n_tbl$file[which(same_run)[1L]], label = "bold"), error = function(e) NULL)
      if (is.list(geom) && is.finite(geom$header_tr) && geom$header_tr > 0) return(as.double(geom$header_tr))
    }
  }
  NULL
}

filter_motion_parameters_for_fd <- function(
  df,
  tr,
  run_label,
  motion_bandstop_bpm = NULL,
  motion_lowpass_hz = NULL,
  motion_filter_order = NULL,
  motion_filter_pad = "constant"
) {
  mot_cols <- c("rot_x", "rot_y", "rot_z", "trans_x", "trans_y", "trans_z")
  miss <- setdiff(mot_cols, names(df))
  if (length(miss) > 0L) {
    stop(
      sprintf(
        "Cannot apply respiratory motion filter in %s: missing required rigid-body column(s): %s",
        run_label, paste(miss, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  n_vol <- nrow(df)
  if (n_vol <= 2L) {
    warning(
      sprintf("Skipping respiratory motion filter in %s because run has only %d volume(s).", run_label, n_vol),
      call. = FALSE
    )
    return(list(applied = FALSE, reason = "too_few_volumes", df = df))
  }

  fs <- 1 / tr
  nyquist <- fs / 2
  nyquist_bpm <- nyquist * 60

  if (!is.null(motion_bandstop_bpm)) {
    order_eff <- if (is.null(motion_filter_order)) 4L else as.integer(motion_filter_order)
    n_passes <- order_eff %/% 4L
    b_min <- motion_bandstop_bpm[1L]
    b_max <- motion_bandstop_bpm[2L]
    stop_hz <- c(b_min, b_max) / 60

    if (stop_hz[2L] < nyquist) {
      eff_hz <- stop_hz
      aliased <- FALSE
    } else {
      # Check whether both bounds lie strictly inside the same Nyquist zone (k * nyq, (k + 1) * nyq)
      near_mult <- abs(stop_hz / nyquist - round(stop_hz / nyquist)) < 1e-9
      zone_idx <- floor(stop_hz / nyquist)
      if (any(near_mult) || zone_idx[1L] != zone_idx[2L]) {
        warning(
          sprintf(
            "motion_bandstop_bpm (%g-%g BPM) straddles the Nyquist frequency (%.2f BPM for TR = %.4g s) in %s; skipping motion filter and using unfiltered FD.",
            b_min, b_max, nyquist_bpm, tr, run_label
          ),
          call. = FALSE
        )
        return(list(applied = FALSE, reason = "straddles_nyquist", df = df))
      }
      eff_hz <- sort(abs(stop_hz - floor((stop_hz + nyquist) / fs) * fs))
      if (eff_hz[1L] <= 0 || eff_hz[2L] >= nyquist || (eff_hz[2L] - eff_hz[1L]) <= 1e-9) {
        warning(
          sprintf(
            "Aliased motion_bandstop_bpm (%g-%g BPM) is invalid for Nyquist %.2f BPM (TR = %.4g s) in %s; skipping motion filter and using unfiltered FD.",
            b_min, b_max, nyquist_bpm, tr, run_label
          ),
          call. = FALSE
        )
        return(list(applied = FALSE, reason = "invalid_aliased_band", df = df))
      }
      aliased <- TRUE
      warning(
        sprintf(
          "motion_bandstop_bpm (%g-%g BPM) exceeds Nyquist (%.2f BPM for TR = %.4g s) in %s; folding to aliased stopband %.2f-%.2f BPM.",
          b_min, b_max, nyquist_bpm, tr, run_label, eff_hz[1L] * 60, eff_hz[2L] * 60
        ),
        call. = FALSE
      )
    }

    f0 <- mean(eff_hz)
    bw <- eff_hz[2L] - eff_hz[1L]
    coeffs <- design_iir_notch(f0 = f0, bandwidth = bw, fs = fs)
    ftype <- "notch"
  } else {
    order_eff <- if (is.null(motion_filter_order)) 4L else as.integer(motion_filter_order)
    order_1way <- order_eff %/% 2L
    n_passes <- 1L
    if (motion_lowpass_hz >= nyquist) {
      warning(
        sprintf(
          "motion_lowpass_hz (%.4g Hz) >= Nyquist frequency (%.4g Hz for TR = %.4g s) in %s; skipping motion filter and using unfiltered FD.",
          motion_lowpass_hz, nyquist, tr, run_label
        ),
        call. = FALSE
      )
      return(list(applied = FALSE, reason = "above_nyquist", df = df))
    }
    coeffs <- design_butterworth_lowpass(order_n = order_1way, cutoff_hz = motion_lowpass_hz, fs = fs)
    ftype <- "lowpass"
    aliased <- FALSE
    eff_hz <- c(0, motion_lowpass_hz)
  }

  df_out <- df
  for (cn in mot_cols) {
    v <- coerce_bids_numeric_col(df[[cn]], cn, run_label)
    if (any(!is.finite(v))) {
      fin_idx <- which(is.finite(v))
      if (length(fin_idx) == 0L) {
        v <- numeric(length(v))
      } else if (length(fin_idx) == 1L) {
        v[!is.finite(v)] <- v[fin_idx]
      } else {
        v <- stats::approx(x = fin_idx, y = v[fin_idx], xout = seq_along(v), method = "linear", rule = 2)$y
      }
    }
    vf <- v
    for (p in seq_len(n_passes)) {
      vf <- filtfilt_1d(vf, b = coeffs$b, a = coeffs$a, padtype = motion_filter_pad, use_zi = TRUE)
    }
    df_out[[cn]] <- vf
    d1 <- c(0, diff(vf))
    c_d1 <- paste0(cn, "_derivative1")
    c_p2 <- paste0(cn, "_power2")
    c_d1p2 <- paste0(cn, "_derivative1_power2")
    if (c_d1 %in% names(df_out)) df_out[[c_d1]] <- d1
    if (c_p2 %in% names(df_out)) df_out[[c_p2]] <- vf^2
    if (c_d1p2 %in% names(df_out)) df_out[[c_d1p2]] <- d1^2
  }

  list(
    applied = TRUE,
    type = ftype,
    order = order_eff,
    tr = tr,
    aliased = aliased,
    effective_hz = eff_hz,
    df = df_out
  )
}

compute_power_fd <- function(df, run_label, head_radius = 50, rot_units = c("rad", "deg")) {
  rot_units <- match.arg(rot_units)
  rot_cols <- c("rot_x", "rot_y", "rot_z")
  trans_cols <- c("trans_x", "trans_y", "trans_z")
  req <- c(rot_cols, trans_cols)
  if (!all(req %in% names(df))) return(NULL)
  n <- nrow(df)
  if (n == 0L) return(numeric(0L))
  mat <- do.call(cbind, lapply(req, function(cn) coerce_bids_numeric_col(df[[cn]], cn, run_label)))
  storage.mode(mat) <- "double"
  if (anyNA(mat)) mat[is.na(mat)] <- 0
  if (identical(rot_units, "deg")) {
    mat[, 1:3] <- mat[, 1:3, drop = FALSE] * (pi / 180)
  }
  if (n == 1L) return(0)
  diffs <- rbind(rep(0, 6L), diff(mat))
  fd <- rowSums(abs(diffs[, 1:3, drop = FALSE])) * head_radius +
    rowSums(abs(diffs[, 4:6, drop = FALSE]))
  fd[1L] <- 0
  as.double(fd)
}

parse_censor_one <- function(item, n_vol = NULL, run_label = "censor", na_strings = default_bids_na_strings()) {
  if (is.character(item) && length(item) == 1L) {
    checkmate::assert_file_exists(item)
    if (grepl("\\.(1d|txt)(\\.gz)?$", item, ignore.case = TRUE)) {
      con <- if (grepl("\\.gz$", item, ignore.case = TRUE)) gzfile(item, "rt") else item
      if (inherits(con, "connection")) on.exit(close(con), add = TRUE)
      lines <- readLines(con, warn = FALSE)
      lines <- trimws(sub("#.*$", "", lines))
      lines <- lines[nzchar(lines)]
      vals <- suppressWarnings(as.numeric(unlist(strsplit(lines, "\\s+"), use.names = FALSE)))
      if (anyNA(vals) || any(!(vals %in% c(0, 1)))) {
        stop(sprintf("Censor file '%s' must contain only 0 (censor) and 1 (keep) values.", item), call. = FALSE)
      }
      valid_vec <- (vals == 1)
    } else {
      c_df <- read_bids_delim(item, na_strings = na_strings)
      is_scrub <- grepl("_scrub\\.(tsv|csv)(\\.gz)?$", item, ignore.case = TRUE) ||
        (ncol(c_df) >= 1L && all(grepl("^(spike|motion_outlier|non_steady_state_outlier)", names(c_df))))
      if (is_scrub) {
        if (ncol(c_df) == 0L) {
          valid_vec <- rep(TRUE, nrow(c_df))
        } else {
          censor_any <- rep(FALSE, nrow(c_df))
          for (cn in names(c_df)) {
            vc <- coerce_bids_numeric_col(c_df[[cn]], cn, basename(item))
            censor_any <- censor_any | (is.finite(vc) & vc != 0)
          }
          valid_vec <- !censor_any
        }
      } else if (ncol(c_df) == 1L) {
        vc <- coerce_bids_numeric_col(c_df[[1L]], names(c_df)[1L], basename(item))
        if (anyNA(vc) || any(!(vc %in% c(0, 1)))) {
          stop(sprintf("Single-column censor file '%s' must contain only 0 and 1 values.", item), call. = FALSE)
        }
        valid_vec <- (vc == 1)
      } else {
        stop(
          sprintf(
            "Unrecognized censor table format in '%s'; supply a 1-column '*_censor.1D' (1=keep, 0=censor) or '*_scrub.tsv' spike indicator table.",
            item
          ),
          call. = FALSE
        )
      }
    }
  } else if (is.logical(item)) {
    if (anyNA(item)) stop(sprintf("%s cannot contain NA values.", run_label), call. = FALSE)
    valid_vec <- as.logical(item)
  } else if (is.numeric(item)) {
    if (anyNA(item) || any(!(item %in% c(0, 1)))) {
      stop(sprintf("Numeric %s must contain only 0 (censor) and 1 (keep) values.", run_label), call. = FALSE)
    }
    valid_vec <- (as.numeric(item) == 1)
  } else {
    stop(sprintf("Unsupported %s input of class '%s'.", run_label, class(item)[1L]), call. = FALSE)
  }

  if (!is.null(n_vol) && length(valid_vec) != n_vol) {
    stop(
      sprintf(
        "%s has length %d, but %s has %d volumes. Censor masks must match full acquisition volume counts.",
        run_label, length(valid_vec), "confounds/BOLD", n_vol
      ),
      call. = FALSE
    )
  }
  valid_vec
}

resolve_censor_inputs <- function(censor, resolved_runs, ctx, na_strings = default_bids_na_strings()) {
  if (is.null(censor)) return(NULL)
  ctx_censor <- list(run_names = resolved_runs, n_time = ctx$n_time, bold_files = ctx$bold_files)
  if (is.character(censor)) {
    checkmate::assert_character(censor, min.len = 1L, any.missing = FALSE, min.chars = 1L)
    al <- align_inputs_to_runs(censor, ctx_censor, input_label = "censor", is_paths = TRUE)
    return(al$items)
  }
  if (is.logical(censor) || is.numeric(censor)) {
    if (length(resolved_runs) != 1L) {
      stop("Multi-run confounds require a character vector of censor file paths or a list of per-run censor vectors.", call. = FALSE)
    }
    return(stats::setNames(list(censor), resolved_runs))
  }
  if (is.list(censor)) {
    al <- align_inputs_to_runs(censor, ctx_censor, input_label = "censor", is_paths = FALSE)
    return(al$items)
  }
  stop("'censor' must be NULL, a character vector of file paths, a logical/0-1 vector, or a list of per-run vectors.", call. = FALSE)
}

select_compcor_columns <- function(col_names, n_compcor, compcor_type) {
  if (is.null(n_compcor) || n_compcor <= 0L) return(character(0L))
  pick_top <- function(prefix, k) {
    hits <- grep(paste0("^", prefix, "_[0-9]+$"), col_names, value = TRUE)
    if (!length(hits)) return(character(0L))
    idx <- as.integer(sub(paste0("^", prefix, "_"), "", hits))
    hits <- hits[order(idx)]
    utils::head(hits, k)
  }
  if (identical(compcor_type, "anatomical")) {
    ac <- pick_top("a_comp_cor", n_compcor)
    if (length(ac)) return(ac)
    return(c(pick_top("c_comp_cor", n_compcor), pick_top("w_comp_cor", n_compcor)))
  }
  if (identical(compcor_type, "split")) {
    sc <- c(pick_top("c_comp_cor", n_compcor), pick_top("w_comp_cor", n_compcor))
    if (length(sc)) return(sc)
    return(pick_top("a_comp_cor", n_compcor))
  }
  pick_top("t_comp_cor", n_compcor)
}

dilate_censor_mask <- function(mask, censor_prev = 0L, censor_next = 0L) {
  n <- length(mask)
  if (n == 0L || !any(mask) || (censor_prev == 0L && censor_next == 0L)) return(mask)
  out <- mask
  idx <- which(mask)
  for (lag in seq.int(-censor_prev, censor_next)) {
    if (lag == 0L) next
    shifted <- idx + lag
    shifted <- shifted[shifted >= 1L & shifted <= n]
    out[shifted] <- TRUE
  }
  out
}

enforce_min_segment_length <- function(valid, min_len = 1L) {
  if (min_len <= 1L || !any(valid)) return(valid)
  runs <- rle(valid)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1L
  out <- valid
  short_valid <- which(runs$values & runs$lengths < min_len)
  for (k in short_valid) {
    out[starts[k]:ends[k]] <- FALSE
  }
  out
}

impute_nuisance_missing_only <- function(z_mat, valid_r, impute = "mean") {
  if (ncol(z_mat) == 0L || nrow(z_mat) == 0L) return(z_mat)
  for (j in seq_len(ncol(z_mat))) {
    vj <- z_mat[, j]
    if (any(!is.finite(vj))) {
      ok_valid <- valid_r & is.finite(vj)
      fill_val <- if (identical(impute, "mean") && any(ok_valid)) mean(vj[ok_valid]) else 0
      vj[!is.finite(vj)] <- fill_val
      z_mat[, j] <- vj
    }
  }
  z_mat
}

finalize_nuisance_matrix <- function(
  z_mat,
  valid_r,
  impute = "mean",
  impute_censored = TRUE,
  center = TRUE,
  scale = FALSE,
  drop_constant = TRUE
) {
  if (ncol(z_mat) == 0L) return(z_mat)
  keep_col <- rep(TRUE, ncol(z_mat))
  for (j in seq_len(ncol(z_mat))) {
    vj <- z_mat[, j]
    ok_valid <- valid_r & is.finite(vj)
    fill_val <- if (identical(impute, "mean") && any(ok_valid)) mean(vj[ok_valid]) else 0
    bad_valid <- valid_r & !is.finite(vj)
    if (any(bad_valid)) vj[bad_valid] <- fill_val
    if (isTRUE(impute_censored)) {
      if (any(!valid_r)) vj[!valid_r] <- fill_val
    } else {
      bad_cens <- !valid_r & !is.finite(vj)
      if (any(bad_cens)) vj[bad_cens] <- fill_val
    }
    if (isTRUE(center)) {
      mu_v <- if (any(valid_r)) mean(vj[valid_r]) else 0
      vj <- vj - mu_v
    }
    sd_v <- if (sum(valid_r) >= 2L) stats::sd(vj[valid_r]) else 0
    if (isTRUE(scale) && is.finite(sd_v) && sd_v > 1e-12) {
      vj <- vj / sd_v
    }
    if (isTRUE(drop_constant) && (!is.finite(sd_v) || sd_v <= 1e-12)) {
      keep_col[j] <- FALSE
    }
    z_mat[, j] <- vj
  }
  z_mat[, keep_col, drop = FALSE]
}

find_matching_aroma_sidecar <- function(ref_row, cand_tbl) {
  if (nrow(cand_tbl) == 0L) return(NA_character_)
  keep <- rep(TRUE, nrow(cand_tbl))
  for (field in c("sub", "ses", "task", "acq", "ce", "rec", "dir", "echo")) {
    b_val <- ref_row[[field]]
    c_val <- cand_tbl[[field]]
    if (!is.na(b_val) && nzchar(b_val)) {
      keep <- keep & !is.na(c_val) & tolower(c_val) == tolower(b_val)
    } else {
      keep <- keep & (is.na(c_val) | !nzchar(c_val))
    }
  }
  b_run <- canonical_run_key(ref_row$run)
  c_run <- canonical_run_key(cand_tbl$run)
  if (!is.na(b_run)) {
    keep <- keep & !is.na(c_run) & (c_run == b_run)
  } else if (!is.na(ref_row$run) && nzchar(ref_row$run)) {
    keep <- keep & !is.na(cand_tbl$run) & (tolower(cand_tbl$run) == tolower(ref_row$run))
  } else {
    keep <- keep & (is.na(cand_tbl$run) | !nzchar(cand_tbl$run))
  }
  idx <- which(keep)
  if (!length(idx)) return(NA_character_)
  if (length(idx) == 1L) return(cand_tbl$file[idx])

  sub_tbl <- cand_tbl[idx, , drop = FALSE]
  stems <- strip_neuro_ext(sub_tbl$file)
  if (length(unique(stems)) == 1L) return(sub_tbl$file[1L])

  # Score candidates: prefer exact space/res match when available and fMRIPost-AROMA lowercase desc
  score <- integer(nrow(sub_tbl))
  if (!is.na(ref_row$space) && nzchar(ref_row$space)) {
    score <- score + 4L * as.integer(!is.na(sub_tbl$space) & sub_tbl$space == ref_row$space)
  }
  if (!is.na(ref_row$res) && nzchar(ref_row$res)) {
    score <- score + 2L * as.integer(!is.na(sub_tbl$res) & sub_tbl$res == ref_row$res)
  }
  score <- score + as.integer(!is.na(sub_tbl$desc) & sub_tbl$desc %in% c("melodic", "aroma"))
  best <- which(score == max(score))
  best_files <- sub_tbl$file[best]
  if (length(unique(strip_neuro_ext(best_files))) == 1L) return(best_files[1L])
  stop(
    sprintf(
      "Multiple AROMA sidecar files matched run '%s': %s",
      basename(ref_row$file), paste(basename(best_files), collapse = ", ")
    ),
    call. = FALSE
  )
}

discover_aroma_classification_for_mixing <- function(mix_path) {
  if (!is.character(mix_path) || length(mix_path) != 1L || is.na(mix_path) || !file.exists(mix_path)) {
    return(NA_character_)
  }
  d <- dirname(mix_path)
  mix_row <- parse_bids_entity_table(mix_path)
  met_cands <- list.files(d, pattern = "_desc-aroma_metrics\\.tsv(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  if (length(met_cands)) {
    hit <- find_matching_aroma_sidecar(mix_row, parse_bids_entity_table(met_cands))
    if (!is.na(hit)) return(hit)
  }
  csv_cands <- list.files(d, pattern = "_AROMAnoiseICs\\.(csv|txt)(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  if (length(csv_cands)) {
    hit <- find_matching_aroma_sidecar(mix_row, parse_bids_entity_table(csv_cands))
    if (!is.na(hit)) return(hit)
  }
  NA_character_
}

read_melodic_mixing_one <- function(path, run_label = "mixing", na_strings = default_bids_na_strings()) {
  checkmate::assert_file_exists(path)
  if (R.utils::isGzipped(path, method = "extension")) {
    loadNamespace("R.utils")
  }
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else path
  lines_head <- local({
    if (inherits(con, "connection")) on.exit(close(con), add = TRUE)
    readLines(con, n = 10L, warn = FALSE)
  })
  lines_head <- trimws(lines_head)
  lines_head <- lines_head[nzchar(lines_head) & !grepl("^#", lines_head)]
  if (!length(lines_head)) {
    stop(sprintf("AROMA mixing file '%s' (%s) is empty.", path, run_label), call. = FALSE)
  }
  first_line <- lines_head[1L]
  sep_guess <- if (grepl("\t", first_line, fixed = TRUE)) {
    "\t"
  } else if (grepl(",", first_line, fixed = TRUE)) {
    ","
  } else {
    "\\s+"
  }
  toks <- if (identical(sep_guess, "\\s+")) {
    unlist(strsplit(first_line, "\\s+"), use.names = FALSE)
  } else {
    unlist(strsplit(first_line, sep_guess, fixed = TRUE), use.names = FALSE)
  }
  toks <- gsub('^["\']|["\']$', "", trimws(toks))
  toks <- toks[nzchar(toks)]
  has_header <- !(length(toks) > 0L && !anyNA(suppressWarnings(as.numeric(toks))))

  dt <- if (identical(sep_guess, "\t")) {
    data.table::fread(path, sep = "\t", header = has_header, na.strings = na_strings, data.table = FALSE)
  } else if (identical(sep_guess, ",")) {
    data.table::fread(path, sep = ",", header = has_header, na.strings = na_strings, data.table = FALSE)
  } else {
    data.table::fread(path, header = has_header, na.strings = na_strings, data.table = FALSE)
  }
  if (nrow(dt) < 3L || ncol(dt) < 1L) {
    stop(
      sprintf(
        "AROMA mixing matrix '%s' (%s) must have at least 3 rows and 1 component column (got %d x %d).",
        basename(path), run_label, nrow(dt), ncol(dt)
      ),
      call. = FALSE
    )
  }
  cols <- lapply(seq_len(ncol(dt)), function(j) {
    coerce_bids_numeric_col(dt[[j]], names(dt)[j], sprintf("%s (%s)", run_label, basename(path)))
  })
  mat <- do.call(cbind, cols)
  storage.mode(mat) <- "double"
  if (any(!is.finite(mat))) {
    stop(
      sprintf("AROMA mixing matrix '%s' (%s) contains missing or non-finite values.", basename(path), run_label),
      call. = FALSE
    )
  }
  if (!has_header || all(grepl("^V[0-9]+$", names(dt)))) {
    colnames(mat) <- sprintf("IC_%02d", seq_len(ncol(mat)))
  } else {
    colnames(mat) <- names(dt)
  }
  mat
}

normalize_mixing_matrix_one <- function(item, run_label = "mixing", na_strings = default_bids_na_strings()) {
  if (is.character(item) && length(item) == 1L) {
    return(list(matrix = read_melodic_mixing_one(item, run_label = run_label, na_strings = na_strings), file = item))
  }
  if (is.data.frame(item)) {
    cols <- lapply(seq_len(ncol(item)), function(j) coerce_bids_numeric_col(item[[j]], names(item)[j], run_label))
    mat <- if (length(cols)) do.call(cbind, cols) else matrix(numeric(0L), nrow = nrow(item), ncol = 0L)
    if (ncol(mat) > 0L) colnames(mat) <- names(item)
  } else if (is.matrix(item) && is.numeric(item)) {
    mat <- item
    storage.mode(mat) <- "double"
    if (is.null(colnames(mat)) && ncol(mat) > 0L) {
      colnames(mat) <- sprintf("IC_%02d", seq_len(ncol(mat)))
    }
  } else {
    stop(sprintf("Unsupported %s input of class '%s'; supply file path(s) or numeric T x K matrix/data.frame.", run_label, class(item)[1L]), call. = FALSE)
  }
  if (nrow(mat) < 3L || ncol(mat) < 1L) {
    stop(sprintf("%s must have at least 3 timepoints and 1 component column (got %d x %d).", run_label, nrow(mat), ncol(mat)), call. = FALSE)
  }
  if (any(!is.finite(mat))) {
    stop(sprintf("%s contains missing or non-finite values.", run_label), call. = FALSE)
  }
  list(matrix = mat, file = NA_character_)
}

validate_aroma_ic_indices <- function(idx, n_comp, run_label) {
  if (length(idx) == 0L) return(integer(0L))
  if (!is.numeric(idx) || any(!is.finite(idx)) || any(idx != floor(idx)) || any(idx < 1) || any(idx > n_comp)) {
    bad <- idx[!is.finite(idx) | idx != floor(idx) | idx < 1 | idx > n_comp]
    stop(
      sprintf(
        "AROMA rejected component indices in %s must be 1-based integers in 1..%d (invalid value(s): %s).",
        run_label, n_comp, paste(utils::head(bad, 6L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  sort(unique(as.integer(idx)))
}

parse_aroma_metrics_df <- function(tbl, n_comp, run_label) {
  nms_low <- tolower(names(tbl))
  cls_pos <- match("classification", nms_low)
  if (is.na(cls_pos)) {
    stop(sprintf("AROMA metrics table in %s must contain a 'classification' column.", run_label), call. = FALSE)
  }
  if (nrow(tbl) != n_comp) {
    stop(
      sprintf(
        "AROMA metrics table in %s has %d rows, but mixing matrix has %d components.",
        run_label, nrow(tbl), n_comp
      ),
      call. = FALSE
    )
  }
  cls_raw <- tolower(trimws(as.character(tbl[[cls_pos]])))
  valid_labels <- c("accepted", "rejected", "signal", "noise", "motion")
  if (anyNA(cls_raw) || any(!(cls_raw %in% valid_labels))) {
    bad_lbl <- unique(tbl[[cls_pos]][is.na(cls_raw) | !(cls_raw %in% valid_labels)])
    stop(
      sprintf(
        "Unrecognized AROMA classification value(s) in %s: %s (expected 'accepted' or 'rejected').",
        run_label, paste(utils::head(as.character(bad_lbl), 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  is_rej <- cls_raw %in% c("rejected", "noise", "motion")
  ic_pos <- match("ic", nms_low)
  if (!is.na(ic_pos) && is.numeric(tbl[[ic_pos]]) && !anyNA(tbl[[ic_pos]])) {
    ic_vals <- as.numeric(tbl[[ic_pos]])
    if (all(ic_vals == floor(ic_vals)) && setequal(as.integer(ic_vals), seq_len(n_comp))) {
      rej_idx <- sort(as.integer(ic_vals[is_rej]))
    } else {
      rej_idx <- which(is_rej)
    }
  } else {
    rej_idx <- which(is_rej)
  }
  out_df <- as.data.frame(tbl, stringsAsFactors = FALSE)
  if (!("IC" %in% names(out_df)) && is.na(ic_pos)) {
    out_df <- cbind(data.frame(IC = seq_len(n_comp), stringsAsFactors = FALSE), out_df)
  }
  list(noise_ics = as.integer(rej_idx), metrics = out_df)
}

parse_aroma_classification_one <- function(cls_item, n_comp, run_label = "classification", na_strings = default_bids_na_strings()) {
  if (is.character(cls_item) && length(cls_item) == 1L) {
    checkmate::assert_file_exists(cls_item)
    if (R.utils::isGzipped(cls_item, method = "extension")) {
      loadNamespace("R.utils")
    }
    con <- if (grepl("\\.gz$", cls_item, ignore.case = TRUE)) gzfile(cls_item, "rt") else cls_item
    lines <- local({
      if (inherits(con, "connection")) on.exit(close(con), add = TRUE)
      readLines(con, warn = FALSE)
    })
    clean_lines <- trimws(sub("#.*$", "", lines))
    clean_lines <- clean_lines[nzchar(clean_lines)]

    is_index_file <- grepl("AROMAnoiseICs|\\.(1d|txt)(\\.gz)?$", cls_item, ignore.case = TRUE) ||
      (length(clean_lines) == 0L) ||
      (length(clean_lines) >= 1L && !grepl("classification", clean_lines[1L], ignore.case = TRUE) &&
         !anyNA(suppressWarnings(as.numeric(unlist(strsplit(clean_lines[1L], "[,[:space:]]+"), use.names = FALSE)))))

    if (is_index_file) {
      if (!length(clean_lines)) {
        rej_idx <- integer(0L)
      } else {
        toks <- unlist(strsplit(paste(clean_lines, collapse = ","), "[,[:space:]]+"), use.names = FALSE)
        toks <- toks[nzchar(toks)]
        vals <- suppressWarnings(as.numeric(toks))
        rej_idx <- validate_aroma_ic_indices(vals, n_comp, sprintf("%s ('%s')", run_label, basename(cls_item)))
      }
      met_df <- data.frame(
        IC = seq_len(n_comp),
        classification = ifelse(seq_len(n_comp) %in% rej_idx, "rejected", "accepted"),
        stringsAsFactors = FALSE
      )
      return(list(noise_ics = rej_idx, metrics = met_df, file = cls_item))
    }

    tbl <- read_bids_delim(cls_item, na_strings = na_strings)
    parsed <- parse_aroma_metrics_df(tbl, n_comp = n_comp, run_label = sprintf("%s ('%s')", run_label, basename(cls_item)))
    return(list(noise_ics = parsed$noise_ics, metrics = parsed$metrics, file = cls_item))
  }

  if (is.data.frame(cls_item)) {
    parsed <- parse_aroma_metrics_df(cls_item, n_comp = n_comp, run_label = run_label)
    return(list(noise_ics = parsed$noise_ics, metrics = parsed$metrics, file = NA_character_))
  }

  if (is.logical(cls_item)) {
    if (length(cls_item) != n_comp || anyNA(cls_item)) {
      stop(
        sprintf("Logical AROMA %s must have length %d (matching mixing components) with no NA values.", run_label, n_comp),
        call. = FALSE
      )
    }
    rej_idx <- which(cls_item)
    met_df <- data.frame(
      IC = seq_len(n_comp),
      classification = ifelse(cls_item, "rejected", "accepted"),
      stringsAsFactors = FALSE
    )
    return(list(noise_ics = as.integer(rej_idx), metrics = met_df, file = NA_character_))
  }

  if (is.numeric(cls_item) || is.integer(cls_item)) {
    rej_idx <- validate_aroma_ic_indices(cls_item, n_comp = n_comp, run_label = run_label)
    met_df <- data.frame(
      IC = seq_len(n_comp),
      classification = ifelse(seq_len(n_comp) %in% rej_idx, "rejected", "accepted"),
      stringsAsFactors = FALSE
    )
    return(list(noise_ics = rej_idx, metrics = met_df, file = NA_character_))
  }

  stop(
    sprintf(
      "Unsupported AROMA %s of class '%s'; supply path(s) to '*_desc-aroma_metrics.tsv' / '*_AROMAnoiseICs.csv', a metrics data.frame, a logical mask, or 1-based integer indices.",
      run_label, class(cls_item)[1L]
    ),
    call. = FALSE
  )
}

aroma_regress_matrix <- function(
  Y,
  M,
  noise_ics,
  include_rows = NULL,
  mode = c("nonaggressive", "aggressive"),
  preserve_mean = TRUE,
  run_label = "run"
) {
  mode <- match.arg(mode)
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  n_t <- nrow(Y)
  if (nrow(M) != n_t) {
    stop(
      sprintf("AROMA mixing matrix in %s has %d rows, but signal matrix has %d rows.", run_label, nrow(M), n_t),
      call. = FALSE
    )
  }
  if (is.null(include_rows)) {
    include_idx <- seq_len(n_t)
  } else {
    if (length(include_rows) != n_t || anyNA(include_rows)) {
      stop(sprintf("valid_rows in %s must have length %d with no NA values.", run_label, n_t), call. = FALSE)
    }
    include_idx <- which(as.logical(include_rows))
  }
  n_fit <- length(include_idx)
  if (n_fit < 3L) {
    stop(sprintf("AROMA regression in %s requires at least 3 usable timepoints (got %d).", run_label, n_fit), call. = FALSE)
  }
  if (any(!is.finite(M))) {
    stop(sprintf("AROMA mixing matrix in %s contains non-finite values.", run_label), call. = FALSE)
  }
  if (ncol(Y) > 0L && any(!is.finite(Y[include_idx, , drop = FALSE]))) {
    stop(sprintf("Signal matrix in %s contains non-finite values on fitted timepoints.", run_label), call. = FALSE)
  }

  noise_ics <- validate_aroma_ic_indices(noise_ics, ncol(M), run_label)
  if (ncol(Y) == 0L) {
    return(list(
      Y = Y,
      applied = length(noise_ics) > 0L,
      reason = if (length(noise_ics) > 0L) NULL else "no_rejected_components",
      n_fit = n_fit,
      n_components = ncol(M),
      n_accepted = ncol(M) - length(noise_ics),
      n_rejected = length(noise_ics),
      mean_pct_var_removed = NA_real_,
      median_pct_var_removed = NA_real_
    ))
  }

  M_sub_raw <- M[include_idx, , drop = FALSE]
  col_rng <- apply(M_sub_raw, 2L, function(col) abs(max(col) - min(col)))
  keep <- which(is.finite(col_rng) & col_rng >= 1e-8)

  is_rej_kept <- keep %in% noise_ics
  rej_pos <- which(is_rej_kept)
  acc_pos <- which(!is_rej_kept)

  if (length(rej_pos) == 0L) {
    Y_out <- if (!isTRUE(preserve_mean)) {
      sweep(Y, 2L, colMeans(Y[include_idx, , drop = FALSE]), "-")
    } else {
      Y
    }
    dimnames(Y_out) <- dimnames(Y)
    return(list(
      Y = Y_out,
      applied = FALSE,
      reason = "no_rejected_components",
      n_fit = n_fit,
      n_components = ncol(M),
      n_accepted = ncol(M) - length(noise_ics),
      n_rejected = 0L,
      mean_pct_var_removed = 0,
      median_pct_var_removed = 0
    ))
  }

  if (identical(mode, "nonaggressive") && length(acc_pos) == 0L) {
    warning(
      sprintf(
        "All %d non-constant AROMA component(s) in %s are classified as rejected; non-aggressive and aggressive regression coincide.",
        length(keep), run_label
      ),
      call. = FALSE
    )
  }

  fit_pos <- if (identical(mode, "nonaggressive")) c(acc_pos, rej_pos) else rej_pos
  if (n_fit <= length(fit_pos) + 1L) {
    stop(
      sprintf(
        "Not enough fitted volumes (%d) in %s to estimate %d AROMA component(s) plus intercept.",
        n_fit, run_label, length(fit_pos)
      ),
      call. = FALSE
    )
  }

  M_use <- M[, keep, drop = FALSE]
  mu_M <- colMeans(M_use[include_idx, , drop = FALSE])
  M_cent <- sweep(M_use, 2L, mu_M, "-")

  mu_Y <- colMeans(Y[include_idx, , drop = FALSE])
  Y_sub <- sweep(Y[include_idx, , drop = FALSE], 2L, mu_Y, "-")

  # Place accepted signal components first so pivoted QR never drops signal in favor of collinear motion
  X_sub <- M_cent[include_idx, fit_pos, drop = FALSE]
  qr_X <- qr(X_sub, tol = 1e-7)
  B_fit <- qr.coef(qr_X, Y_sub)
  if (!is.matrix(B_fit)) {
    B_fit <- matrix(B_fit, nrow = ncol(X_sub), ncol = ncol(Y))
  }
  B_fit[!is.finite(B_fit)] <- 0

  rej_in_fit <- if (identical(mode, "nonaggressive")) {
    seq.int(length(acc_pos) + 1L, length(fit_pos))
  } else {
    seq_along(rej_pos)
  }
  B_rej <- B_fit[rej_in_fit, , drop = FALSE]
  Y_clean <- Y - M_cent[, rej_pos, drop = FALSE] %*% B_rej
  if (!isTRUE(preserve_mean)) {
    Y_clean <- sweep(Y_clean, 2L, mu_Y, "-")
  }
  dimnames(Y_clean) <- dimnames(Y)

  tss_before <- colSums(Y_sub^2)
  Y_clean_sub <- sweep(Y_clean[include_idx, , drop = FALSE], 2L, colMeans(Y_clean[include_idx, , drop = FALSE]), "-")
  tss_after <- colSums(Y_clean_sub^2)
  active_vox <- is.finite(tss_before) & tss_before > 1e-12
  pct_removed <- if (any(active_vox)) {
    100 * (tss_before[active_vox] - tss_after[active_vox]) / tss_before[active_vox]
  } else {
    numeric(0L)
  }

  list(
    Y = Y_clean,
    applied = TRUE,
    reason = NULL,
    n_fit = n_fit,
    n_components = ncol(M),
    n_accepted = ncol(M) - length(noise_ics),
    n_rejected = length(noise_ics),
    mean_pct_var_removed = if (length(pct_removed)) mean(pct_removed) else NA_real_,
    median_pct_var_removed = if (length(pct_removed)) stats::median(pct_removed) else NA_real_
  )
}

align_medusa_aroma_to_runs <- function(aroma, resolved_runs, n_time = NULL) {
  if (!inherits(aroma, "medusa_aroma")) {
    stop("'aroma' must be a 'medusa_aroma' object produced by read_fmriprep_aroma().", call. = FALSE)
  }
  a_runs <- names(aroma$mixing)
  if (length(a_runs) != length(resolved_runs)) {
    stop(
      sprintf(
        "AROMA object has %d run(s) (%s), which does not match target run count (%d: %s).",
        length(a_runs), paste(a_runs, collapse = ", "),
        length(resolved_runs), paste(resolved_runs, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (identical(a_runs, resolved_runs)) {
    ord <- seq_along(resolved_runs)
  } else if (setequal(a_runs, resolved_runs)) {
    ord <- match(resolved_runs, a_runs)
  } else {
    ca <- canonical_run_key(a_runs)
    ct <- canonical_run_key(resolved_runs)
    if (!anyNA(ca) && !anyNA(ct) && !anyDuplicated(ca) && !anyDuplicated(ct) && setequal(ca, ct)) {
      ord <- match(ct, ca)
    } else if (length(resolved_runs) == 1L) {
      ord <- 1L
    } else {
      stop(
        sprintf(
          "AROMA run names (%s) do not match target runs (%s).",
          paste(a_runs, collapse = ", "), paste(resolved_runs, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }
  out <- aroma
  out$mixing <- stats::setNames(aroma$mixing[ord], resolved_runs)
  out$noise_ics <- stats::setNames(aroma$noise_ics[ord], resolved_runs)
  out$signal_ics <- stats::setNames(aroma$signal_ics[ord], resolved_runs)
  out$metrics <- stats::setNames(aroma$metrics[ord], resolved_runs)
  if (!is.null(aroma$summary) && nrow(aroma$summary) == length(resolved_runs)) {
    out$summary <- aroma$summary[ord, , drop = FALSE]
    out$summary$run <- resolved_runs
    rownames(out$summary) <- NULL
  }
  if (!is.null(n_time)) {
    for (rn in resolved_runs) {
      if (nrow(out$mixing[[rn]]) != n_time[[rn]]) {
        stop(
          sprintf(
            "AROMA mixing matrix for '%s' has %d rows, but target data has %d volumes.",
            rn, nrow(out$mixing[[rn]]), n_time[[rn]]
          ),
          call. = FALSE
        )
      }
    }
  }
  out
}

#' Read fMRIPrep or postprocessed confound tables into clock-preserving MEDuSA nuisance and censoring inputs
#'
#' Converts fMRIPrep \code{*_desc-confounds_timeseries.tsv} or \code{BrainGnomes}
#' \code{*_desc-<desc>_confounds.tsv} files (or in-memory confound
#' \code{data.frame}s) into run-matched \code{valid_rows} logical masks and
#' finite \code{nuisance} design matrices ready for
#' \code{\link{estimate_crossrun_hrf}} and \code{\link{deconvolve_crossrun_hrf}}.
#' Enforces MEDuSA's clock-preservation invariant: censored volumes are marked
#' \code{FALSE} in \code{valid_rows} rather than deleted, initial \code{"n/a"}
#' derivative entries on volume 1 are zero-filled, and censored or missing
#' nuisance entries are safely imputed so \code{nuisance} is strictly finite on
#' all rows.
#'
#' @param confounds Character vector of confound TSV paths
#'   (\code{*_desc-confounds_timeseries.tsv}, \code{*_desc-<desc>_confounds.tsv},
#'   \code{.tsv.gz}, \code{.csv}), a single \code{data.frame}, or a list of
#'   per-run \code{data.frame}s.
#' @param bold Optional \code{"medusa_nifti_ts"} object, named list of run BOLD
#'   matrices, or single BOLD matrix used to match run ordering, validate BIDS
#'   \code{sub-}/\code{ses-}/\code{task-}/\code{space-}/\code{cohort-} entities,
#'   check \code{BrainGnomes} provenance sidecars against volume deletion, and
#'   verify exact volume counts (\code{nrow(confounds) == nrow(bold)}).
#' @param run_names Optional explicit character vector of unique run names when
#'   \code{bold} is \code{NULL}.
#' @param strategy Predefined nuisance regressor set:
#'   \itemize{
#'     \item \code{"24p"} (default): 6 realignment parameters (\code{trans_x},
#'       \code{trans_y}, \code{trans_z}, \code{rot_x}, \code{rot_y}, \code{rot_z}),
#'       their temporal derivatives, and quadratic terms (Friston 24-parameter model).
#'     \item \code{"12p"}: 6 motion parameters + 6 temporal derivatives.
#'     \item \code{"6p"}: 6 rigid-body motion parameters.
#'     \item \code{"csf_wm"}: mean \code{csf} and \code{white_matter} physiological
#'       signals (2p model; recommended alongside non-aggressive ICA-AROMA per
#'       Pruim et al., 2015 and Parkes et al., 2018).
#'     \item \code{"csf_wm_8p"}: \code{csf} and \code{white_matter} with their
#'       temporal derivatives and quadratic terms (8p physiological model).
#'     \item \code{"24p_csf_wm"}: \code{"24p"} + \code{csf} and \code{white_matter}
#'       with their derivatives and quadratic terms (32p physiological model without GSR).
#'     \item \code{"27p"}: \code{"24p"} + mean \code{csf}, \code{white_matter},
#'       and \code{global_signal} without physiological derivatives or quadratics.
#'     \item \code{"36p"}: \code{"24p_csf_wm"} + \code{global_signal} with its
#'       derivative and quadratic terms (Satterthwaite 36p model).
#'     \item \code{"compcor"}: \code{"12p"} motion + top \code{n_compcor}
#'       (default \code{5}) CompCor components + discrete cosine high-pass basis.
#'     \item \code{"custom"}: only columns requested via \code{columns},
#'       \code{column_pattern}, and/or \code{n_compcor}.
#'     \item \code{"none"}: empty 0-column nuisance matrices (censoring only).
#'   }
#' @param columns Optional character vector of explicit confound column names,
#'   shorthand sets (\code{"6p"}, \code{"12p"}, \code{"24p"}, \code{"27p"},
#'   \code{"csf_wm"}, \code{"csf_wm_8p"}, \code{"36p"}), or zero-padded index
#'   range specifications such as \code{"a_comp_cor_<0-4>"} or
#'   \code{"cosine<0,2-3>"} (in addition to \code{strategy}, or exact set when
#'   \code{strategy = "custom"}).
#' @param column_pattern Optional regular expression selecting additional columns
#'   present in each confound table (e.g., \code{"^aroma_motion_"}).
#' @param n_compcor Optional integer count of leading CompCor components to add
#'   (defaults to \code{5L} when \code{strategy = "compcor"} and \code{0L}
#'   otherwise).
#' @param compcor_type CompCor family when \code{n_compcor > 0}:
#'   \code{"anatomical"} (\code{a_comp_cor_*}, falling back to split
#'   \code{c_comp_cor_*} + \code{w_comp_cor_*}), \code{"split"}
#'   (\code{c_comp_cor_*} + \code{w_comp_cor_*}), or \code{"temporal"}
#'   (\code{t_comp_cor_*}).
#' @param include_cosines Whether to include fMRIPrep discrete cosine high-pass
#'   columns (\code{cosine00}, \code{cosine01}, ...): \code{"auto"} (default;
#'   included automatically whenever CompCor columns are selected, per fMRIPrep
#'   recommendations), \code{TRUE}, or \code{FALSE}.
#' @param censor Optional external censoring mask(s): character vector of
#'   \code{*_censor.1D} (\code{1} = keep, \code{0} = censor) or
#'   \code{*_scrub.tsv} (one-hot spike indicator) paths, a logical or \code{0/1}
#'   vector (single run), or a list of per-run logical/\code{0-1} vectors. When
#'   \code{censor} is supplied and \code{fd_threshold} is not explicitly set,
#'   \code{fd_threshold} defaults to \code{NULL} so the external censor mask is
#'   not overridden by a second FD cutoff.
#' @param fd_threshold Framewise displacement censoring cutoff in mm (default
#'   \code{0.5}; volumes with \code{FD > fd_threshold} are marked \code{FALSE}
#'   in \code{valid_rows}). Set \code{NULL} or \code{Inf} to disable FD censoring.
#' @param fd_col Column name supplying framewise displacement, or \code{"auto"}
#'   (default), which prefers \code{"framewise_displacement_filtered"} when
#'   present (e.g., respiratory-filtered \code{BrainGnomes} postprocessing
#'   outputs), then \code{"framewise_displacement"}, and falls back to computing
#'   Power (2012) FD from the 6 rigid-body motion parameters when missing or all
#'   \code{NA}.
#' @param head_radius Head radius in mm (default \code{50}) used when computing
#'   Power (2012) framewise displacement from the 6 rigid-body motion columns.
#' @param rot_units Rotation parameter units (\code{"rad"} default for fMRIPrep,
#'   or \code{"deg"}) used when computing Power (2012) FD from the 6 rigid-body
#'   motion columns.
#' @param tr Optional repetition time in seconds (scalar or per-run numeric
#'   vector), required only when \code{motion_bandstop_bpm} or
#'   \code{motion_lowpass_hz} is specified and TR cannot be inferred
#'   automatically from \code{bold} or a sibling \code{*_bold.json} /
#'   \code{*_bold.nii(.gz)} file.
#' @param motion_bandstop_bpm Optional length-2 numeric vector
#'   \code{c(min_bpm, max_bpm)} specifying a respiratory band-stop (notch) filter
#'   in breaths per minute (e.g., \code{c(12, 18)} for typical adult breathing
#'   at 0.2--0.3 Hz, or \code{c(18.58, 25.73)} per Fair et al., 2020) applied
#'   zero-phase via forward-backward biquad filtering to the 6 rigid-body motion
#'   parameters before recomputing Power (2012) FD. Suppresses fast-TR multiband
#'   respiratory B0 pseudomotion that would otherwise over-censor clean volumes.
#'   Mutually exclusive with \code{motion_lowpass_hz}. Default \code{NULL}.
#' @param motion_lowpass_hz Optional numeric scalar specifying a Butterworth
#'   low-pass cutoff frequency in Hz (e.g., \code{0.2} Hz; Gratton et al., 2020)
#'   applied zero-phase via forward-backward filtering to the 6 rigid-body
#'   motion parameters before recomputing Power (2012) FD. Mutually exclusive
#'   with \code{motion_bandstop_bpm}. Default \code{NULL}.
#' @param motion_filter_order Optional integer effective zero-phase filter order
#'   (default \code{NULL}, which uses \code{4L}). Must be a positive multiple of
#'   \code{4} for \code{motion_bandstop_bpm} (each forward-backward 2nd-order
#'   biquad pass contributes order 4) or an even integer \code{>= 2} for
#'   \code{motion_lowpass_hz} (a forward-backward \code{N}-pole Butterworth
#'   low-pass has effective order \code{2 * N}).
#' @param motion_filter_pad Boundary reflection padding mode for zero-phase
#'   forward-backward filtering: \code{"constant"} (default), \code{"odd"},
#'   \code{"even"}, or \code{"zero"}.
#' @param motion_filter_target Whether respiratory motion filtering updates only
#'   Power (2012) FD and censoring (\code{"fd_only"}, default; preserves
#'   unfiltered motion nuisance columns so respiratory B0 variance remaining in
#'   unfiltered BOLD is still modeled) or also replaces the 6 rigid-body motion
#'   columns and their Friston-24 derivatives/quadratics in \code{nuisance}
#'   (\code{"fd_and_nuisance"}).
#' @param dvars_threshold Optional DVARS censoring cutoff (e.g., \code{1.5} on
#'   \code{std_dvars}). Default \code{NULL} (disabled).
#' @param dvars_col Column name supplying DVARS when \code{dvars_threshold} is
#'   set. Default \code{"std_dvars"}.
#' @param include_non_steady_state Logical; if \code{TRUE} (default), marks any
#'   volume flagged in \code{non_steady_state_outlier*} columns as \code{FALSE}
#'   in \code{valid_rows}.
#' @param include_motion_outliers Logical; if \code{TRUE}, also marks volumes
#'   flagged in fMRIPrep's precomputed \code{motion_outlier*} columns as
#'   \code{FALSE} in \code{valid_rows}. Default \code{FALSE} when controlled via
#'   \code{fd_threshold}/\code{dvars_threshold}.
#' @param drop_initial_volumes Integer count of initial dummy volumes to mark
#'   \code{FALSE} in \code{valid_rows} without deleting rows (default \code{0L}).
#' @param censor_prev Integer count of preceding volumes to scrub before each
#'   FD/DVARS/motion-outlier spike (default \code{0L}).
#' @param censor_next Integer count of subsequent volumes to scrub after each
#'   FD/DVARS/motion-outlier spike (default \code{0L}).
#' @param min_segment_length Minimum contiguous uncensored segment length in
#'   volumes; shorter isolated islands between motion spikes are also marked
#'   \code{FALSE} (default \code{1L}).
#' @param impute Imputation rule for missing entries (beyond volume-1 derivative
#'   zeros) and censored rows: \code{"mean"} (usable-row column mean, default)
#'   or \code{"zero"}.
#' @param impute_censored Logical; if \code{TRUE} (default), replaces censored
#'   rows (\code{!valid_rows}) in \code{nuisance} with the imputed baseline value
#'   so motion spikes have zero leverage even if inspected or used unmasked.
#' @param center Logical; if \code{TRUE} (default), centers each nuisance column
#'   on \code{valid_rows}.
#' @param scale Logical; if \code{TRUE}, scales each nuisance column to unit SD
#'   on \code{valid_rows}. Default \code{FALSE}.
#' @param drop_constant Logical; if \code{TRUE} (default), removes columns with
#'   near-zero variance on \code{valid_rows} (such as all-zero outlier columns).
#' @param aroma Optional \code{"medusa_aroma"} object from
#'   \code{\link{read_fmriprep_aroma}}. When supplied, applies the matching
#'   ICA-AROMA projection (\code{aroma_mode}) to selected nuisance columns prior
#'   to final centering and censoring sanitization so subsequent joint nuisance
#'   regression in \code{\link{estimate_crossrun_hrf}} and
#'   \code{\link{deconvolve_crossrun_hrf}} cannot re-introduce AROMA motion
#'   variance into AROMA-denoised BOLD (Lindquist et al., 2019). Framewise
#'   displacement and \code{valid_rows} remain strictly derived from physical
#'   (or respiratory-filtered) head-motion parameters.
#' @param aroma_mode AROMA regression mode when \code{aroma} is supplied:
#'   \code{"nonaggressive"} (default; Pruim et al., 2015 partial regression) or
#'   \code{"aggressive"}.
#' @param aroma_fit_valid Logical; if \code{TRUE}, estimates AROMA component
#'   weights on uncensored \code{valid_rows} only rather than full MELODIC run
#'   timecourses (\code{FALSE}, default).
#' @param na_strings Character vector of missing-value tokens passed to
#'   \code{data.table::fread}. Default \code{c("n/a", "N/A", "NA", ".", "NaN", "")}.
#'
#' @return An S3 list of class \code{"medusa_confounds"} containing:
#'   \itemize{
#'     \item \code{nuisance}: Run-named list of finite numeric \code{T x P}
#'       matrices.
#'     \item \code{valid_rows}: Run-named list of length-\code{T} logical
#'       observation masks (\code{TRUE} = usable, \code{FALSE} = censored).
#'     \item \code{summary}: Per-run QC \code{data.frame} reporting volume
#'       counts, censored counts/percentages, outlier breakdowns, mean/max FD,
#'       whether respiratory motion filtering was applied, and retained nuisance
#'       column counts.
#'     \item \code{framewise_displacement}: Run-named list of numeric FD vectors
#'       used for QC and censoring (respiratory-filtered when
#'       \code{motion_bandstop_bpm} or \code{motion_lowpass_hz} is active).
#'     \item \code{framewise_displacement_raw}: Run-named list of unfiltered FD
#'       vectors when \code{motion_bandstop_bpm} or \code{motion_lowpass_hz} is
#'       requested (\code{NULL} otherwise).
#'     \item \code{variance_confounds}: Run-named list of numeric frame-level
#'       variance-envelope matrices (\code{framewise_displacement},
#'       \code{std_dvars}/\code{dvars}, and optional \code{rvt}) for
#'       \code{\link{decon_epoch_regression}}.
#'   }
#'   The \code{"medusa_confounds"} object can be passed directly as
#'   \code{nuisance = conf} to \code{\link{estimate_crossrun_hrf}} (automatically
#'   supplying both \code{nuisance} and \code{valid_rows}), as
#'   \code{variance_confounds = conf} to \code{\link{decon_epoch_regression}}, or as
#'   \code{valid = conf} to \code{\link{align_medusa_ts}}.
#'
#' @seealso \code{\link{read_bids_events}}, \code{\link{read_fmriprep_subject}},
#'   \code{\link{read_fmriprep_aroma}}, \code{\link{apply_aroma}},
#'   \code{\link{read_nifti_ts}}, \code{\link{estimate_crossrun_hrf}}
#'
#' @examples
#' f_conf <- tempfile("sub-01_task-mid_run-01_", fileext = "_desc-confounds_timeseries.tsv")
#' conf_df <- data.frame(
#'   trans_x = c(0.01, 0.02, 0.40, 0.03, 0.02),
#'   trans_y = c(0.00, -0.01, 0.25, 0.01, 0.00),
#'   trans_z = c(0.02, 0.03, 0.80, 0.04, 0.03),
#'   rot_x = c(0.001, 0.001, 0.01, 0.002, 0.001),
#'   rot_y = c(0.000, 0.001, 0.005, 0.001, 0.001),
#'   rot_z = c(-0.001, 0.000, 0.004, 0.000, -0.001),
#'   framewise_displacement = c("n/a", "0.08", "0.92", "0.12", "0.05"),
#'   non_steady_state_outlier00 = c(1L, 0L, 0L, 0L, 0L),
#'   stringsAsFactors = FALSE
#' )
#' utils::write.table(conf_df, f_conf, sep = "\t", row.names = FALSE, quote = FALSE)
#'
#' conf <- read_fmriprep_confounds(f_conf, strategy = "12p", fd_threshold = 0.5)
#' conf$valid_rows$run1
#' conf$summary
#' unlink(f_conf)
#'
#' @export
read_fmriprep_confounds <- function(
  confounds,
  bold = NULL,
  run_names = NULL,
  strategy = c("24p", "6p", "12p", "csf_wm", "csf_wm_8p", "24p_csf_wm", "27p", "36p", "compcor", "custom", "none"),
  columns = NULL,
  column_pattern = NULL,
  n_compcor = NULL,
  compcor_type = c("anatomical", "split", "temporal"),
  include_cosines = "auto",
  censor = NULL,
  fd_threshold = 0.5,
  fd_col = "auto",
  head_radius = 50,
  rot_units = c("rad", "deg"),
  tr = NULL,
  motion_bandstop_bpm = NULL,
  motion_lowpass_hz = NULL,
  motion_filter_order = NULL,
  motion_filter_pad = c("constant", "odd", "even", "zero"),
  motion_filter_target = c("fd_only", "fd_and_nuisance"),
  dvars_threshold = NULL,
  dvars_col = "std_dvars",
  include_non_steady_state = TRUE,
  include_motion_outliers = FALSE,
  drop_initial_volumes = 0L,
  censor_prev = 0L,
  censor_next = 0L,
  min_segment_length = 1L,
  impute = c("mean", "zero"),
  impute_censored = TRUE,
  center = TRUE,
  scale = FALSE,
  drop_constant = TRUE,
  aroma = NULL,
  aroma_mode = c("nonaggressive", "aggressive"),
  aroma_fit_valid = FALSE,
  na_strings = c("n/a", "N/A", "NA", ".", "NaN", "")
) {
  fd_threshold_supplied <- !missing(fd_threshold)
  strategy <- match.arg(strategy)
  compcor_type <- match.arg(compcor_type)
  rot_units <- match.arg(rot_units)
  motion_filter_pad <- match.arg(motion_filter_pad)
  motion_filter_target <- match.arg(motion_filter_target)
  impute <- match.arg(impute)
  aroma_mode <- match.arg(aroma_mode)
  checkmate::assert_flag(aroma_fit_valid)
  checkmate::assert_character(columns, min.chars = 1L, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_string(column_pattern, min.chars = 1L, null.ok = TRUE)
  if (is.null(n_compcor)) {
    n_compcor <- if (identical(strategy, "compcor")) 5L else 0L
  }
  checkmate::assert_integerish(n_compcor, lower = 0L, len = 1L, any.missing = FALSE)
  n_compcor <- as.integer(n_compcor)
  if (!identical(include_cosines, "auto")) {
    checkmate::assert_flag(include_cosines)
  }
  if (!is.null(censor) && !fd_threshold_supplied) {
    fd_threshold <- NULL
  }
  if (!is.null(fd_threshold) && !(is.numeric(fd_threshold) && length(fd_threshold) == 1L && is.infinite(fd_threshold) && fd_threshold > 0)) {
    checkmate::assert_number(fd_threshold, lower = 0, finite = TRUE)
  } else {
    fd_threshold <- NULL
  }
  checkmate::assert_string(fd_col, min.chars = 1L)
  checkmate::assert_number(head_radius, lower = 1, finite = TRUE)
  if (!is.null(tr)) {
    checkmate::assert_numeric(tr, lower = 0.01, finite = TRUE, min.len = 1L, any.missing = FALSE)
    tr_names <- names(tr)
    tr <- as.double(tr)
    names(tr) <- tr_names
  }
  if (!is.null(motion_bandstop_bpm) && !is.null(motion_lowpass_hz)) {
    stop("Specify at most one of 'motion_bandstop_bpm' or 'motion_lowpass_hz', not both.", call. = FALSE)
  }
  if (!is.null(motion_bandstop_bpm)) {
    checkmate::assert_numeric(motion_bandstop_bpm, lower = 1e-6, finite = TRUE, len = 2L, any.missing = FALSE)
    motion_bandstop_bpm <- as.double(motion_bandstop_bpm)
    if (abs(motion_bandstop_bpm[2L] - motion_bandstop_bpm[1L]) <= 1e-9) {
      stop("'motion_bandstop_bpm' lower and upper bounds must be distinct.", call. = FALSE)
    }
    if (motion_bandstop_bpm[2L] < motion_bandstop_bpm[1L]) {
      warning(
        sprintf(
          "'motion_bandstop_bpm' bounds (%g, %g) are reversed; reordering to c(%g, %g).",
          motion_bandstop_bpm[1L], motion_bandstop_bpm[2L],
          motion_bandstop_bpm[2L], motion_bandstop_bpm[1L]
        ),
        call. = FALSE
      )
      motion_bandstop_bpm <- rev(motion_bandstop_bpm)
    }
    if (!is.null(motion_filter_order)) {
      checkmate::assert_integerish(motion_filter_order, lower = 4L, len = 1L, any.missing = FALSE)
      motion_filter_order <- as.integer(motion_filter_order)
      if (motion_filter_order %% 4L != 0L) {
        stop(
          "For band-stop notch filtering ('motion_bandstop_bpm'), 'motion_filter_order' must be a positive multiple of 4 (e.g., 4 or 8).",
          call. = FALSE
        )
      }
    }
  } else if (!is.null(motion_lowpass_hz)) {
    checkmate::assert_number(motion_lowpass_hz, lower = 1e-4, finite = TRUE)
    motion_lowpass_hz <- as.double(motion_lowpass_hz)
    if (!is.null(motion_filter_order)) {
      checkmate::assert_integerish(motion_filter_order, lower = 2L, len = 1L, any.missing = FALSE)
      motion_filter_order <- as.integer(motion_filter_order)
      if (motion_filter_order %% 2L != 0L) {
        stop(
          "For zero-phase Butterworth low-pass filtering ('motion_lowpass_hz'), 'motion_filter_order' must be an even integer >= 2 (e.g., 2, 4, or 8).",
          call. = FALSE
        )
      }
    }
  } else if (!is.null(motion_filter_order)) {
    checkmate::assert_integerish(motion_filter_order, lower = 2L, len = 1L, any.missing = FALSE)
    motion_filter_order <- as.integer(motion_filter_order)
  }
  wants_motion_filter <- !is.null(motion_bandstop_bpm) || !is.null(motion_lowpass_hz)

  if (!is.null(dvars_threshold) && !(is.numeric(dvars_threshold) && length(dvars_threshold) == 1L && is.infinite(dvars_threshold) && dvars_threshold > 0)) {
    checkmate::assert_number(dvars_threshold, lower = 0, finite = TRUE)
  } else {
    dvars_threshold <- NULL
  }
  checkmate::assert_string(dvars_col, min.chars = 1L)
  checkmate::assert_flag(include_non_steady_state)
  checkmate::assert_flag(include_motion_outliers)
  checkmate::assert_integerish(drop_initial_volumes, lower = 0L, len = 1L, any.missing = FALSE)
  drop_initial_volumes <- as.integer(drop_initial_volumes)
  checkmate::assert_integerish(censor_prev, lower = 0L, len = 1L, any.missing = FALSE)
  censor_prev <- as.integer(censor_prev)
  checkmate::assert_integerish(censor_next, lower = 0L, len = 1L, any.missing = FALSE)
  censor_next <- as.integer(censor_next)
  checkmate::assert_integerish(min_segment_length, lower = 1L, len = 1L, any.missing = FALSE)
  min_segment_length <- as.integer(min_segment_length)
  checkmate::assert_flag(impute_censored)
  checkmate::assert_flag(center)
  checkmate::assert_flag(scale)
  checkmate::assert_flag(drop_constant)
  checkmate::assert_character(na_strings, any.missing = FALSE)

  ctx <- extract_bold_run_context(bold = bold, run_names = run_names)
  prov_list <- assert_no_scrubbed_bold_provenance(ctx$bold_files)
  if (length(prov_list) > 0L && any(vapply(prov_list, function(p) isTRUE(p$confound_regression), logical(1L)))) {
    wants_nuis <- !identical(strategy, "none") || !is.null(columns) || !is.null(column_pattern) || n_compcor > 0L
    if (wants_nuis) {
      warning(
        "BrainGnomes provenance indicates 'confound_regression' was already executed on BOLD input(s); including nuisance regressors again may re-introduce filtered artifact variance (set strategy = 'none' unless intentional).",
        call. = FALSE
      )
    }
  }
  src_files <- NULL

  if (is.character(confounds)) {
    checkmate::assert_character(confounds, min.len = 1L, any.missing = FALSE, min.chars = 1L)
    aligned <- align_inputs_to_runs(confounds, ctx, input_label = "confounds", is_paths = TRUE)
    resolved_runs <- aligned$run_names
    src_files <- as.character(aligned$items)
    raw_list <- lapply(aligned$items, read_bids_delim, na_strings = na_strings)
    names(raw_list) <- resolved_runs
  } else if (is.data.frame(confounds)) {
    if (!is.null(ctx$run_names) && length(ctx$run_names) > 1L) {
      stop("Multi-run 'bold' requires a character vector of confound paths or a list of per-run data.frames.", call. = FALSE)
    }
    resolved_runs <- if (!is.null(ctx$run_names)) ctx$run_names else "run1"
    raw_list <- stats::setNames(list(as.data.frame(confounds, stringsAsFactors = FALSE)), resolved_runs)
    src_files <- NA_character_
  } else if (is.list(confounds)) {
    checkmate::assert_list(confounds, types = "data.frame", min.len = 1L)
    aligned <- align_inputs_to_runs(confounds, ctx, input_label = "confounds", is_paths = FALSE)
    resolved_runs <- aligned$run_names
    raw_list <- lapply(aligned$items, as.data.frame, stringsAsFactors = FALSE)
    src_files <- rep(NA_character_, length(resolved_runs))
  } else {
    stop("'confounds' must be a character vector of file paths, a data.frame, or a list of data.frames.", call. = FALSE)
  }

  censor_items <- resolve_censor_inputs(censor, resolved_runs, ctx, na_strings = na_strings)
  aroma_aligned <- if (!is.null(aroma)) align_medusa_aroma_to_runs(aroma, resolved_runs) else NULL
  strategy_cols <- confound_strategy_columns(strategy)

  n_runs <- length(resolved_runs)
  nuis_list <- raw_nuis_list <- valid_list <- fd_list <- fd_raw_list <- vc_list <- stats::setNames(vector("list", n_runs), resolved_runs)
  sum_rows <- vector("list", n_runs)

  for (r in seq_len(n_runs)) {
    rname <- resolved_runs[r]
    df <- raw_list[[r]]
    n_vol <- nrow(df)
    run_label <- sprintf("confounds[['%s']]", rname)

    if (!is.null(ctx$n_time) && n_vol != ctx$n_time[[rname]]) {
      stop(
        sprintf(
          "%s has %d rows, but BOLD run '%s' has %d volumes. Do not delete volumes prior to MEDuSA; keep full acquisition grids and censor via valid_rows.",
          run_label, n_vol, rname, ctx$n_time[[rname]]
        ),
        call. = FALSE
      )
    }

    # 1. Censoring / valid_rows construction
    valid_ext <- if (!is.null(censor_items)) {
      parse_censor_one(
        censor_items[[r]],
        n_vol = n_vol,
        run_label = sprintf("censor[['%s']]", rname),
        na_strings = na_strings
      )
    } else {
      rep(TRUE, n_vol)
    }
    censor_ext <- !valid_ext

    censor_init <- rep(FALSE, n_vol)
    if (drop_initial_volumes > 0L && n_vol > 0L) {
      censor_init[seq_len(min(drop_initial_volumes, n_vol))] <- TRUE
    }

    censor_nss <- rep(FALSE, n_vol)
    if (isTRUE(include_non_steady_state) && n_vol > 0L) {
      nss_cols <- grep("^non_steady_state_outlier", names(df), value = TRUE)
      for (nc in nss_cols) {
        v_nss <- coerce_bids_numeric_col(df[[nc]], nc, run_label)
        censor_nss <- censor_nss | (is.finite(v_nss) & v_nss != 0)
      }
    }

    fd_vec <- rep(NA_real_, n_vol)
    fd_col_resolved <- if (identical(fd_col, "auto")) {
      if ("framewise_displacement_filtered" %in% names(df)) {
        "framewise_displacement_filtered"
      } else {
        "framewise_displacement"
      }
    } else {
      fd_col
    }
    allow_fd_fallback <- fd_col %in% c("auto", "framewise_displacement", "framewise_displacement_filtered")

    if (fd_col_resolved %in% names(df)) {
      fd_cand <- coerce_bids_numeric_col(df[[fd_col_resolved]], fd_col_resolved, run_label)
      non_row1_all_na <- (n_vol == 0L) || all(is.na(if (n_vol >= 2L) fd_cand[-1L] else fd_cand))
      if (non_row1_all_na && allow_fd_fallback) {
        fd_fb <- compute_power_fd(df, run_label, head_radius = head_radius, rot_units = rot_units)
        if (!is.null(fd_fb)) {
          fd_vec <- fd_fb
        } else {
          fd_vec <- fd_cand
          if (n_vol >= 1L && is.na(fd_vec[1L])) fd_vec[1L] <- 0
        }
      } else {
        fd_vec <- fd_cand
        if (n_vol >= 1L && is.na(fd_vec[1L])) fd_vec[1L] <- 0
      }
    } else if (allow_fd_fallback) {
      fd_fb <- compute_power_fd(df, run_label, head_radius = head_radius, rot_units = rot_units)
      if (!is.null(fd_fb)) {
        fd_vec <- fd_fb
      } else if (!is.null(fd_threshold) && !wants_motion_filter) {
        stop(
          sprintf(
            "Neither '%s' nor the 6 rigid-body motion columns found in %s; set fd_threshold = NULL to disable FD censoring or specify fd_col.",
            fd_col_resolved, run_label
          ),
          call. = FALSE
        )
      }
    } else if (!is.null(fd_threshold) && !wants_motion_filter) {
      stop(
        sprintf(
          "Column '%s' not found in %s; set fd_threshold = NULL to disable FD censoring or specify fd_col.",
          fd_col, run_label
        ),
        call. = FALSE
      )
    }

    # Optional respiratory notch / low-pass filtering of the 6 rigid-body motion parameters prior to Power FD
    mf_applied <- FALSE
    df_nuis <- df
    fd_raw_r <- fd_vec
    if (wants_motion_filter) {
      has_unfilt_fd <- FALSE
      if ("framewise_displacement" %in% names(df)) {
        fd_unfilt <- coerce_bids_numeric_col(df[["framewise_displacement"]], "framewise_displacement", run_label)
        if (n_vol >= 1L && is.na(fd_unfilt[1L])) fd_unfilt[1L] <- 0
        if (!all(is.na(if (n_vol >= 2L) fd_unfilt[-1L] else fd_unfilt))) {
          fd_raw_r <- fd_unfilt
          has_unfilt_fd <- TRUE
        }
      }
      if (!has_unfilt_fd) {
        fd_fb_raw <- compute_power_fd(df, run_label, head_radius = head_radius, rot_units = rot_units)
        if (!is.null(fd_fb_raw)) {
          fd_raw_r <- fd_fb_raw
          has_unfilt_fd <- TRUE
        }
      }
      if ("framewise_displacement_filtered" %in% names(df) &&
          fd_col %in% c("auto", "framewise_displacement_filtered")) {
        warning(
          sprintf(
            "%s already contains 'framewise_displacement_filtered'; skipping duplicate respiratory motion filter and using existing filtered FD.",
            run_label
          ),
          call. = FALSE
        )
      } else {
        if (!(fd_col %in% c("auto", "framewise_displacement", "framewise_displacement_filtered"))) {
          warning(
            sprintf(
              "Respiratory motion filter requested in %s; recomputing Power (2012) FD from filtered motion parameters instead of fd_col = '%s'.",
              run_label, fd_col
            ),
            call. = FALSE
          )
        }
        bold_file_r <- if (!is.null(ctx$bold_files) && rname %in% names(ctx$bold_files)) ctx$bold_files[[rname]] else NULL
        tr_r <- if (!is.null(tr)) {
          if (length(tr) == 1L) {
            tr[1L]
          } else if (!is.null(names(tr)) && rname %in% names(tr)) {
            tr[[rname]]
          } else if (length(tr) == n_runs) {
            tr[r]
          } else {
            stop(
              sprintf("Length of 'tr' (%d) must be 1 or match number of runs (%d).", length(tr), n_runs),
              call. = FALSE
            )
          }
        } else if (!is.null(ctx$tr)) {
          ctx$tr
        } else {
          infer_tr_from_confound_sidecar(src_files[r], bold_file = bold_file_r)
        }
        if (is.null(tr_r) || !is.finite(tr_r) || tr_r <= 0) {
          stop(
            sprintf(
              "Sampling interval 'tr' (in seconds) is required when 'motion_bandstop_bpm' or 'motion_lowpass_hz' is set (could not infer TR for %s); supply 'tr' or pass a 'medusa_nifti_ts' object to 'bold'.",
              run_label
            ),
            call. = FALSE
          )
        }
        mf_res <- filter_motion_parameters_for_fd(
          df = df,
          tr = tr_r,
          run_label = run_label,
          motion_bandstop_bpm = motion_bandstop_bpm,
          motion_lowpass_hz = motion_lowpass_hz,
          motion_filter_order = motion_filter_order,
          motion_filter_pad = motion_filter_pad
        )
        if (isTRUE(mf_res$applied)) {
          fd_vec <- compute_power_fd(mf_res$df, run_label, head_radius = head_radius, rot_units = rot_units)
          mf_applied <- TRUE
          if (identical(motion_filter_target, "fd_and_nuisance")) {
            df_nuis <- mf_res$df
            if ("framewise_displacement" %in% names(df_nuis)) {
              df_nuis[["framewise_displacement"]] <- fd_vec
            }
          }
        } else {
          non_row1_na <- (n_vol > 0L) && all(is.na(if (n_vol >= 2L) fd_vec[-1L] else fd_vec))
          if (non_row1_na && has_unfilt_fd) {
            fd_vec <- fd_raw_r
          } else if (non_row1_na && !is.null(fd_threshold)) {
            stop(
              sprintf(
                "Column '%s' and the 6 rigid-body motion columns are unavailable for fallback FD censoring in %s; set fd_threshold = NULL to disable FD censoring.",
                fd_col, run_label
              ),
              call. = FALSE
            )
          }
        }
      }
    }
    fd_raw_list[[r]] <- fd_raw_r

    censor_fd <- if (!is.null(fd_threshold) && n_vol > 0L) {
      !is.finite(fd_vec) | (fd_vec > fd_threshold)
    } else {
      rep(FALSE, n_vol)
    }

    censor_dvars <- rep(FALSE, n_vol)
    if (!is.null(dvars_threshold) && n_vol > 0L) {
      if (!(dvars_col %in% names(df))) {
        stop(
          sprintf(
            "Column '%s' not found in %s; set dvars_threshold = NULL or specify dvars_col.",
            dvars_col, run_label
          ),
          call. = FALSE
        )
      }
      dv_vec <- coerce_bids_numeric_col(df[[dvars_col]], dvars_col, run_label)
      if (n_vol >= 1L && is.na(dv_vec[1L])) dv_vec[1L] <- 0
      censor_dvars <- !is.finite(dv_vec) | (dv_vec > dvars_threshold)
    }

    censor_mo <- rep(FALSE, n_vol)
    if (isTRUE(include_motion_outliers) && n_vol > 0L) {
      mo_cols <- grep("^motion_outlier", names(df), value = TRUE)
      for (mc in mo_cols) {
        v_mo <- coerce_bids_numeric_col(df[[mc]], mc, run_label)
        censor_mo <- censor_mo | (is.finite(v_mo) & v_mo != 0)
      }
      if (mf_applied && any(censor_mo & !censor_fd)) {
        warning(
          sprintf(
            "include_motion_outliers = TRUE in %s censored %d additional volume(s) from precomputed 'motion_outlier*' columns (typically derived from unfiltered FD/DVARS); set include_motion_outliers = FALSE to rely on respiratory-filtered FD.",
            run_label, sum(censor_mo & !censor_fd)
          ),
          call. = FALSE
        )
      }
    }

    spike_censor <- dilate_censor_mask(censor_fd | censor_dvars | censor_mo, censor_prev, censor_next)
    valid_r <- valid_ext & !(censor_init | censor_nss | spike_censor)
    valid_r <- enforce_min_segment_length(valid_r, min_segment_length)
    valid_list[[r]] <- valid_r
    fd_list[[r]] <- fd_vec

    # Extract non-negative/physiological variance envelopes (FD, std_dvars/dvars, RVT) for Stage-3 covariance regression
    vc_cols_r <- list()
    if (n_vol > 0L && any(is.finite(fd_vec))) {
      fd_vc <- as.numeric(fd_vec)
      if (!is.finite(fd_vc[1L])) fd_vc[1L] <- 0.0
      bad_fd <- !is.finite(fd_vc)
      if (any(bad_fd)) {
        fin_fd <- fd_vc[!bad_fd]
        fd_vc[bad_fd] <- if (length(fin_fd)) stats::median(fin_fd) else 0.0
      }
      fd_vc <- pmax(0.0, fd_vc)
      if (n_vol >= 2L && is.finite(stats::sd(fd_vc)) && stats::sd(fd_vc) > 1e-10) {
        vc_cols_r[["framewise_displacement"]] <- fd_vc
      }
    }
    dv_cand <- if (!is.null(dvars_col) && nzchar(dvars_col) && (dvars_col %in% names(df))) {
      dvars_col
    } else {
      intersect(c("std_dvars", "dvars"), names(df))[1L]
    }
    if (!is.na(dv_cand) && n_vol > 0L) {
      dv_vc <- tryCatch(coerce_bids_numeric_col(df[[dv_cand]], dv_cand, run_label), error = function(e) NULL)
      if (!is.null(dv_vc) && any(is.finite(dv_vc))) {
        fin_dv <- dv_vc[is.finite(dv_vc)]
        med_dv <- if (length(fin_dv)) stats::median(fin_dv) else 0.0
        dv_vc[!is.finite(dv_vc)] <- med_dv
        dv_vc <- pmax(0.0, dv_vc)
        if (n_vol >= 2L && is.finite(stats::sd(dv_vc)) && stats::sd(dv_vc) > 1e-10) {
          vc_cols_r[[dv_cand]] <- dv_vc
        }
      }
    }
    for (phys_cn in intersect(c("rvt", "respiratory_volume_per_time"), names(df))) {
      if (n_vol > 0L) {
        ph_vc <- tryCatch(coerce_bids_numeric_col(df[[phys_cn]], phys_cn, run_label), error = function(e) NULL)
        if (!is.null(ph_vc) && any(is.finite(ph_vc))) {
          fin_ph <- ph_vc[is.finite(ph_vc)]
          ph_vc[!is.finite(ph_vc)] <- if (length(fin_ph)) stats::median(fin_ph) else 0.0
          if (n_vol >= 2L && is.finite(stats::sd(ph_vc)) && stats::sd(ph_vc) > 1e-10) {
            vc_cols_r[[phys_cn]] <- ph_vc
          }
        }
      }
    }
    vc_list[[r]] <- if (length(vc_cols_r) > 0L) {
      do.call(cbind, vc_cols_r)
    } else {
      matrix(numeric(0L), nrow = n_vol, ncol = 0L)
    }

    # 2. Nuisance regressor selection & finite sanitization
    if (identical(strategy, "none") && is.null(columns) && is.null(column_pattern) && n_compcor == 0L) {
      z_mat <- matrix(numeric(0L), nrow = n_vol, ncol = 0L)
    } else {
      cc_cols <- select_compcor_columns(names(df_nuis), n_compcor, compcor_type)
      if (identical(strategy, "compcor") && n_compcor > 0L && length(cc_cols) == 0L) {
        stop(sprintf("No CompCor columns matching compcor_type = '%s' found in %s.", compcor_type, run_label), call. = FALSE)
      }
      pat_cols <- if (!is.null(column_pattern)) grep(column_pattern, names(df_nuis), value = TRUE) else character(0L)
      exp_columns <- expand_confound_column_specs(columns, names(df_nuis))
      if (length(exp_columns)) {
        resolvable <- vapply(exp_columns, can_resolve_confound_column, logical(1L), avail_names = names(df_nuis))
        if (any(!resolvable)) {
          stop(
            sprintf(
              "Requested column(s) not found in %s: %s",
              run_label, paste(exp_columns[!resolvable], collapse = ", ")
            ),
            call. = FALSE
          )
        }
      }
      req_cols <- unique(c(strategy_cols, cc_cols, exp_columns, pat_cols))
      use_cos <- if (identical(include_cosines, "auto")) {
        any(grepl("_comp_cor_", req_cols))
      } else {
        isTRUE(include_cosines)
      }
      if (use_cos) {
        cos_cols <- grep("^cosine[0-9]+$", names(df_nuis), value = TRUE)
        if (length(cos_cols)) {
          cos_idx <- as.integer(sub("^cosine", "", cos_cols))
          cos_cols <- cos_cols[order(cos_idx)]
          req_cols <- unique(c(req_cols, cos_cols))
        }
      }

      if (length(req_cols) == 0L) {
        z_mat <- matrix(numeric(0L), nrow = n_vol, ncol = 0L)
        raw_nuis_list[[r]] <- z_mat
        aroma_cleaned_r <- !is.null(aroma_aligned) && length(aroma_aligned$noise_ics[[r]]) > 0L
      } else {
        col_vecs <- lapply(req_cols, function(cn) expand_fmriprep_column(df_nuis, cn, run_label))
        z_mat <- do.call(cbind, col_vecs)
        colnames(z_mat) <- req_cols
        storage.mode(z_mat) <- "double"

        is_deriv_like <- grepl("(_derivative1|^framewise_displacement(_filtered)?$|^dvars$|^std_dvars$|^rmsd$)", req_cols)
        if (n_vol >= 1L && any(is_deriv_like)) {
          d_idx <- which(is_deriv_like)
          bad_row1 <- d_idx[!is.finite(z_mat[1L, d_idx])]
          if (length(bad_row1)) z_mat[1L, bad_row1] <- 0
        }

        z_raw_r <- impute_nuisance_missing_only(z_mat, valid_r = valid_r, impute = impute)
        raw_nuis_list[[r]] <- z_raw_r

        if (!is.null(aroma_aligned)) {
          ar_r <- aroma_regress_matrix(
            Y = z_raw_r,
            M = aroma_aligned$mixing[[r]],
            noise_ics = aroma_aligned$noise_ics[[r]],
            include_rows = if (isTRUE(aroma_fit_valid)) valid_r else NULL,
            mode = aroma_mode,
            preserve_mean = TRUE,
            run_label = run_label
          )
          z_for_final <- ar_r$Y
          aroma_cleaned_r <- isTRUE(ar_r$applied)
        } else {
          z_for_final <- z_raw_r
          aroma_cleaned_r <- FALSE
        }

        z_mat <- finalize_nuisance_matrix(
          z_mat = z_for_final,
          valid_r = valid_r,
          impute = impute,
          impute_censored = impute_censored,
          center = center,
          scale = scale,
          drop_constant = drop_constant
        )
      }
    }
    if (is.null(raw_nuis_list[[r]])) {
      raw_nuis_list[[r]] <- z_mat
      aroma_cleaned_r <- !is.null(aroma_aligned) && length(aroma_aligned$noise_ics[[r]]) > 0L
    }
    nuis_list[[r]] <- z_mat

    fd_fin <- fd_vec[is.finite(fd_vec)]
    fd_step <- if (n_vol >= 2L) fd_vec[-1L][is.finite(fd_vec[-1L])] else fd_fin
    sum_rows[[r]] <- data.frame(
      run = rname,
      file = src_files[r],
      n_volumes = as.integer(n_vol),
      n_valid = as.integer(sum(valid_r)),
      n_censored = as.integer(sum(!valid_r)),
      pct_censored = if (n_vol > 0L) 100 * sum(!valid_r) / n_vol else 0,
      n_fd_outliers = as.integer(sum(censor_fd)),
      n_dvars_outliers = as.integer(sum(censor_dvars)),
      n_non_steady_state = as.integer(sum(censor_nss)),
      n_external_censor = as.integer(sum(censor_ext)),
      mean_fd = if (length(fd_step)) mean(fd_step) else NA_real_,
      max_fd = if (length(fd_fin)) max(fd_fin) else NA_real_,
      motion_filtered = mf_applied,
      aroma_cleaned = aroma_cleaned_r,
      n_nuisance_cols = as.integer(ncol(z_mat)),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, sum_rows)
  rownames(summary_df) <- NULL

  motion_filter_meta <- if (wants_motion_filter) {
    list(
      type = if (!is.null(motion_bandstop_bpm)) "notch" else "lowpass",
      bandstop_bpm = motion_bandstop_bpm,
      lowpass_hz = motion_lowpass_hz,
      order = if (is.null(motion_filter_order)) 4L else motion_filter_order,
      pad = motion_filter_pad,
      target = motion_filter_target
    )
  } else {
    NULL
  }

  aroma_meta <- if (!is.null(aroma_aligned)) {
    list(
      applied = any(summary_df$aroma_cleaned),
      mode = aroma_mode,
      fit_valid_rows = aroma_fit_valid
    )
  } else {
    NULL
  }

  conf_out <- structure(
    list(
      nuisance = nuis_list,
      valid_rows = valid_list,
      summary = summary_df,
      framewise_displacement = fd_list,
      framewise_displacement_raw = if (wants_motion_filter) fd_raw_list else NULL,
      variance_confounds = vc_list,
      settings = list(
        strategy = strategy,
        fd_threshold = fd_threshold,
        fd_col = fd_col,
        dvars_threshold = dvars_threshold,
        has_external_censor = !is.null(censor_items),
        include_non_steady_state = include_non_steady_state,
        include_motion_outliers = include_motion_outliers,
        drop_initial_volumes = drop_initial_volumes,
        censor_prev = censor_prev,
        censor_next = censor_next,
        min_segment_length = min_segment_length,
        impute = impute,
        impute_censored = impute_censored,
        center = center,
        scale = scale,
        drop_constant = drop_constant,
        motion_filter = motion_filter_meta,
        aroma = aroma_meta
      )
    ),
    class = "medusa_confounds"
  )
  attr(conf_out, "raw_nuisance") <- raw_nuis_list
  conf_out
}

#' @export
print.medusa_confounds <- function(x, ...) {
  n_runs <- length(x$valid_rows)
  fd_lbl <- if (!is.null(x$settings$fd_threshold)) {
    sprintf("%.3g mm", x$settings$fd_threshold)
  } else if (isTRUE(x$settings$has_external_censor)) {
    "external censor"
  } else {
    "none"
  }
  mf <- x$settings$motion_filter
  mf_lbl <- if (!is.null(mf)) {
    if (identical(mf$type, "notch")) {
      sprintf(" | fd_filter: notch %g-%g BPM", mf$bandstop_bpm[1L], mf$bandstop_bpm[2L])
    } else {
      sprintf(" | fd_filter: lowpass %g Hz", mf$lowpass_hz)
    }
  } else {
    ""
  }
  ar <- x$settings$aroma
  ar_lbl <- if (is.list(ar) && isTRUE(ar$applied)) {
    sprintf(" | aroma: %s", ar$mode)
  } else {
    ""
  }
  cat(sprintf(
    "<medusa_confounds> %d run%s | strategy: %s | fd_threshold: %s%s%s\n",
    n_runs,
    if (n_runs == 1L) "" else "s",
    x$settings$strategy,
    fd_lbl,
    mf_lbl,
    ar_lbl
  ))
  if (!is.null(x$summary) && nrow(x$summary) > 0L) {
    disp <- x$summary[, c("run", "n_volumes", "n_valid", "n_censored", "pct_censored", "mean_fd", "max_fd", "n_nuisance_cols"), drop = FALSE]
    disp$pct_censored <- round(disp$pct_censored, 1L)
    disp$mean_fd <- round(disp$mean_fd, 3L)
    disp$max_fd <- round(disp$max_fd, 3L)
    print(disp, row.names = FALSE)
  }
  invisible(x)
}

#' Read ICA-AROMA mixing matrices and component classifications for MEDuSA denoising
#'
#' Ingests \code{fMRIPost-AROMA} (\code{*_desc-melodic_mixing.tsv} +
#' \code{*_desc-aroma_metrics.tsv}), legacy fMRIPrep \code{<= 23.0}
#' (\code{*_desc-MELODIC_mixing.tsv} + \code{*_AROMAnoiseICs.csv}), or in-memory
#' MELODIC mixing matrices and component classifications into a run-aligned
#' \code{"medusa_aroma"} bundle for \code{\link{apply_aroma}},
#' \code{\link{read_fmriprep_confounds}}, and \code{\link{read_fmriprep_subject}}.
#'
#' Both headerless MELODIC mixing TSVs (standard FSL/\code{fMRIPost-AROMA}
#' format) and headered TSV/CSV tables are detected automatically. BIDS file
#' inputs are matched to \code{bold} runs using acquisition entities
#' (\code{sub-}, \code{ses-}, \code{task-}, \code{acq-}, \code{dir-},
#' \code{run-}) without requiring template \code{space-} or \code{res-} tags to
#' match, since \code{fMRIPost-AROMA} derives run-level MELODIC timecourses in
#' \code{MNI152NLin6Asym} (\code{res-2}) regardless of the analysis output space.
#'
#' @param mixing Character vector of MELODIC mixing file paths
#'   (\code{*_desc-melodic_mixing.tsv}, \code{*_desc-MELODIC_mixing.tsv},
#'   \code{.tsv.gz}, \code{.csv}), a single numeric \code{T x K} matrix or
#'   \code{data.frame} (single run), or a list of per-run \code{T x K} matrices /
#'   \code{data.frame}s.
#' @param classification Component classification specification identifying
#'   rejected motion components:
#'   \itemize{
#'     \item \code{NULL} (default when \code{mixing} is file paths): auto-discovers
#'       sibling \code{*_desc-aroma_metrics.tsv(.gz)} or
#'       \code{*_AROMAnoiseICs.csv} files in the same directory as each mixing
#'       file.
#'     \item Character vector of paths to \code{*_desc-aroma_metrics.tsv}
#'       (\code{classification == "rejected"}) or \code{*_AROMAnoiseICs.csv}
#'       (comma-separated 1-based component indices; empty file indicates zero
#'       rejected components).
#'     \item Single 1-based integer vector (e.g., \code{c(2L, 5L)}), logical
#'       mask of length \code{K} (\code{TRUE} = rejected), or metrics
#'       \code{data.frame} with a \code{classification} column (for single-run
#'       inputs), or a per-run list of any of these.
#'   }
#' @param bold Optional \code{"medusa_nifti_ts"} object, named list of run BOLD
#'   matrices, single BOLD matrix, or character vector of run names used to
#'   validate acquisition entities, align run ordering, and verify volume counts.
#' @param run_names Optional explicit character vector of unique run names when
#'   \code{bold} is \code{NULL}.
#' @param na_strings Character vector of missing-value tokens passed to
#'   \code{data.table::fread}. Default \code{c("n/a", "N/A", "NA", ".", "NaN", "")}.
#'
#' @return An S3 list of class \code{"medusa_aroma"} containing:
#'   \itemize{
#'     \item \code{mixing}: Run-named list of finite \code{T x K} numeric
#'       MELODIC mixing matrices.
#'     \item \code{noise_ics}: Run-named list of 1-based integer vectors of
#'       rejected motion component indices.
#'     \item \code{signal_ics}: Run-named list of 1-based integer vectors of
#'       accepted signal component indices.
#'     \item \code{metrics}: Run-named list of per-component classification
#'       \code{data.frame}s.
#'     \item \code{summary}: Per-run \code{data.frame} summarizing volume counts
#'       and accepted/rejected component counts.
#'   }
#'
#' @seealso \code{\link{apply_aroma}}, \code{\link{read_fmriprep_subject}},
#'   \code{\link{read_fmriprep_confounds}}
#'
#' @examples
#' f_mix <- tempfile("sub-01_task-mid_run-01_", fileext = "_desc-melodic_mixing.tsv")
#' f_met <- sub("_desc-melodic_mixing\\.tsv$", "_desc-aroma_metrics.tsv", f_mix)
#' mix_mat <- matrix(rnorm(30 * 4), nrow = 30, ncol = 4)
#' utils::write.table(mix_mat, f_mix, sep = "\t", row.names = FALSE, col.names = FALSE)
#' met_df <- data.frame(
#'   IC = 1:4,
#'   classification = c("accepted", "rejected", "accepted", "rejected")
#' )
#' utils::write.table(met_df, f_met, sep = "\t", row.names = FALSE, quote = FALSE)
#'
#' ar <- read_fmriprep_aroma(f_mix)
#' ar
#' unlink(c(f_mix, f_met))
#'
#' @export
read_fmriprep_aroma <- function(
  mixing,
  classification = NULL,
  bold = NULL,
  run_names = NULL,
  na_strings = c("n/a", "N/A", "NA", ".", "NaN", "")
) {
  checkmate::assert_character(na_strings, any.missing = FALSE)
  ctx <- extract_bold_run_context(bold = bold, run_names = run_names)

  validate_aroma_entities <- function(in_tbl, bold_tbl, input_label) {
    for (field in c("sub", "ses", "task")) {
      b_vals <- unique(bold_tbl[[field]][!is.na(bold_tbl[[field]]) & nzchar(bold_tbl[[field]])])
      i_vals <- unique(in_tbl[[field]][!is.na(in_tbl[[field]]) & nzchar(in_tbl[[field]])])
      if (length(b_vals) > 0L && length(i_vals) > 0L && !all(i_vals %in% b_vals)) {
        stop(
          sprintf(
            "BIDS '%s-' entity mismatch between BOLD runs (%s) and %s (%s).",
            field, paste(b_vals, collapse = ", "), input_label, paste(i_vals, collapse = ", ")
          ),
          call. = FALSE
        )
      }
    }
  }

  ctx_no_space <- ctx
  ctx_no_space$bold_files <- NULL

  if (is.character(mixing)) {
    checkmate::assert_character(mixing, min.len = 1L, any.missing = FALSE, min.chars = 1L)
    if (!is.null(ctx$bold_files)) {
      validate_aroma_entities(parse_bids_entity_table(mixing), parse_bids_entity_table(ctx$bold_files), "AROMA mixing")
      aligned_mix <- align_inputs_to_runs(
        mixing,
        list(run_names = ctx$run_names, n_time = ctx$n_time, bold_files = NULL),
        input_label = "mixing",
        is_paths = TRUE
      )
      # Re-check run-signature ordering against bold_files when unnamed multi-run paths are passed
      if ((is.null(names(mixing)) || !any(nzchar(names(mixing)))) && length(mixing) > 1L) {
        m_tbl <- parse_bids_entity_table(mixing)
        b_tbl <- parse_bids_entity_table(ctx$bold_files)
        sig_m <- bids_run_signature(m_tbl)
        sig_b <- bids_run_signature(b_tbl)
        if (!anyDuplicated(sig_m) && !anyDuplicated(sig_b) && setequal(sig_m, sig_b)) {
          ord_m <- match(sig_b, sig_m)
          aligned_mix$items <- stats::setNames(mixing[ord_m], ctx$run_names)
        }
      }
    } else {
      aligned_mix <- align_inputs_to_runs(mixing, ctx_no_space, input_label = "mixing", is_paths = TRUE)
    }
    resolved_runs <- aligned_mix$run_names
    mix_items <- as.list(aligned_mix$items)
  } else if (is.matrix(mixing) || is.data.frame(mixing)) {
    if (!is.null(ctx$run_names) && length(ctx$run_names) > 1L) {
      stop("Multi-run 'bold' requires a character vector of mixing paths or a list of per-run mixing matrices.", call. = FALSE)
    }
    resolved_runs <- if (!is.null(ctx$run_names)) ctx$run_names else "run1"
    mix_items <- stats::setNames(list(mixing), resolved_runs)
  } else if (is.list(mixing)) {
    checkmate::assert_list(mixing, min.len = 1L)
    aligned_mix <- align_inputs_to_runs(mixing, ctx_no_space, input_label = "mixing", is_paths = FALSE)
    resolved_runs <- aligned_mix$run_names
    mix_items <- aligned_mix$items
  } else {
    stop("'mixing' must be a character vector of file paths, a numeric matrix/data.frame, or a list of matrices.", call. = FALSE)
  }

  n_runs <- length(resolved_runs)
  if (is.null(classification)) {
    if (!is.character(mixing)) {
      stop("'classification' must be supplied when 'mixing' is an in-memory matrix or list.", call. = FALSE)
    }
    auto_cls <- vapply(mix_items, discover_aroma_classification_for_mixing, character(1L))
    if (anyNA(auto_cls)) {
      miss_r <- resolved_runs[is.na(auto_cls)]
      stop(
        sprintf(
          "Could not auto-discover sibling '*_desc-aroma_metrics.tsv' or '*_AROMAnoiseICs.csv' for run(s): %s; supply 'classification' explicitly.",
          paste(miss_r, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    cls_items <- as.list(stats::setNames(auto_cls, resolved_runs))
  } else if (is.character(classification)) {
    checkmate::assert_character(classification, min.len = 1L, any.missing = FALSE, min.chars = 1L)
    if (!is.null(ctx$bold_files)) {
      validate_aroma_entities(parse_bids_entity_table(classification), parse_bids_entity_table(ctx$bold_files), "AROMA classification")
    }
    ctx_cls <- list(run_names = resolved_runs, n_time = ctx$n_time, bold_files = NULL)
    al_cls <- align_inputs_to_runs(classification, ctx_cls, input_label = "classification", is_paths = TRUE)
    if ((is.null(names(classification)) || !any(nzchar(names(classification)))) && n_runs > 1L && !is.null(ctx$bold_files)) {
      c_tbl <- parse_bids_entity_table(classification)
      b_tbl <- parse_bids_entity_table(ctx$bold_files)
      sig_c <- bids_run_signature(c_tbl)
      sig_b <- bids_run_signature(b_tbl)
      if (!anyDuplicated(sig_c) && !anyDuplicated(sig_b) && setequal(sig_c, sig_b)) {
        al_cls$items <- stats::setNames(classification[match(sig_b, sig_c)], resolved_runs)
      }
    }
    cls_items <- as.list(al_cls$items)
  } else if (is.data.frame(classification) || is.logical(classification) || is.numeric(classification) || is.integer(classification)) {
    if (n_runs != 1L) {
      stop("Multi-run AROMA inputs require a character vector of classification file paths or a per-run list.", call. = FALSE)
    }
    cls_items <- stats::setNames(list(classification), resolved_runs)
  } else if (is.list(classification)) {
    al_cls <- align_inputs_to_runs(
      classification,
      list(run_names = resolved_runs, n_time = ctx$n_time, bold_files = NULL),
      input_label = "classification",
      is_paths = FALSE
    )
    cls_items <- al_cls$items
  } else {
    stop(
      "'classification' must be NULL, file path(s), a metrics data.frame, a logical mask, a 1-based integer vector, or a per-run list.",
      call. = FALSE
    )
  }

  mix_list <- noise_list <- signal_list <- metrics_list <- stats::setNames(vector("list", n_runs), resolved_runs)
  sum_rows <- vector("list", n_runs)

  for (r in seq_len(n_runs)) {
    rname <- resolved_runs[r]
    m_norm <- normalize_mixing_matrix_one(mix_items[[r]], run_label = sprintf("mixing[['%s']]", rname), na_strings = na_strings)
    m_mat <- m_norm$matrix
    if (!is.null(ctx$n_time) && nrow(m_mat) != ctx$n_time[[rname]]) {
      stop(
        sprintf(
          "AROMA mixing[['%s']] has %d rows, but BOLD run '%s' has %d volumes.",
          rname, nrow(m_mat), rname, ctx$n_time[[rname]]
        ),
        call. = FALSE
      )
    }
    n_comp <- ncol(m_mat)
    c_norm <- parse_aroma_classification_one(
      cls_items[[r]],
      n_comp = n_comp,
      run_label = sprintf("classification[['%s']]", rname),
      na_strings = na_strings
    )
    rej <- c_norm$noise_ics
    acc <- setdiff(seq_len(n_comp), rej)

    mix_list[[r]] <- m_mat
    noise_list[[r]] <- rej
    signal_list[[r]] <- acc
    metrics_list[[r]] <- c_norm$metrics
    sum_rows[[r]] <- data.frame(
      run = rname,
      mixing_file = m_norm$file,
      classification_file = c_norm$file,
      n_volumes = as.integer(nrow(m_mat)),
      n_components = as.integer(n_comp),
      n_accepted = as.integer(length(acc)),
      n_rejected = as.integer(length(rej)),
      pct_rejected = 100 * length(rej) / n_comp,
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, sum_rows)
  rownames(summary_df) <- NULL

  structure(
    list(
      mixing = mix_list,
      noise_ics = noise_list,
      signal_ics = signal_list,
      metrics = metrics_list,
      summary = summary_df
    ),
    class = "medusa_aroma"
  )
}

#' @export
print.medusa_aroma <- function(x, ...) {
  n_runs <- length(x$mixing)
  tot_rej <- if (!is.null(x$summary)) sum(x$summary$n_rejected) else 0L
  tot_comp <- if (!is.null(x$summary)) sum(x$summary$n_components) else 0L
  cat(sprintf(
    "<medusa_aroma> %d run%s | %d/%d total component(s) classified as rejected\n",
    n_runs,
    if (n_runs == 1L) "" else "s",
    tot_rej,
    tot_comp
  ))
  if (!is.null(x$summary) && nrow(x$summary) > 0L) {
    disp <- x$summary[, c("run", "n_volumes", "n_components", "n_accepted", "n_rejected", "pct_rejected"), drop = FALSE]
    disp$pct_rejected <- round(disp$pct_rejected, 1L)
    print(disp, row.names = FALSE)
  }
  invisible(x)
}

apply_aroma_to_confounds_obj <- function(
  conf_obj,
  aroma_obj,
  mode = c("nonaggressive", "aggressive"),
  fit_valid_rows = FALSE,
  force = FALSE
) {
  mode <- match.arg(mode)
  if (!isTRUE(force) && is.list(conf_obj$settings$aroma) && isTRUE(conf_obj$settings$aroma$applied)) {
    warning("Confound nuisance matrices were already AROMA-cleaned; skipping duplicate AROMA projection (set force = TRUE to override).", call. = FALSE)
    return(conf_obj)
  }
  runs <- names(conf_obj$nuisance)
  n_t <- stats::setNames(vapply(conf_obj$valid_rows, length, integer(1L)), runs)
  aroma_al <- align_medusa_aroma_to_runs(aroma_obj, runs, n_time = n_t)
  raw_nuis <- attr(conf_obj, "raw_nuisance", exact = TRUE)

  imp <- if (!is.null(conf_obj$settings$impute)) conf_obj$settings$impute else "mean"
  imp_cens <- if (!is.null(conf_obj$settings$impute_censored)) conf_obj$settings$impute_censored else TRUE
  ctr <- if (!is.null(conf_obj$settings$center)) conf_obj$settings$center else TRUE
  scl <- if (!is.null(conf_obj$settings$scale)) conf_obj$settings$scale else FALSE
  drp <- if (!is.null(conf_obj$settings$drop_constant)) conf_obj$settings$drop_constant else TRUE

  cleaned_flags <- logical(length(runs))
  for (i in seq_along(runs)) {
    rn <- runs[i]
    valid_r <- conf_obj$valid_rows[[rn]]
    z_in <- if (is.list(raw_nuis) && !is.null(raw_nuis[[rn]])) raw_nuis[[rn]] else conf_obj$nuisance[[rn]]
    if (ncol(z_in) == 0L) {
      cleaned_flags[i] <- length(aroma_al$noise_ics[[rn]]) > 0L
      next
    }
    ar_r <- aroma_regress_matrix(
      Y = z_in,
      M = aroma_al$mixing[[rn]],
      noise_ics = aroma_al$noise_ics[[rn]],
      include_rows = if (isTRUE(fit_valid_rows)) valid_r else NULL,
      mode = mode,
      preserve_mean = TRUE,
      run_label = sprintf("confounds[['%s']]", rn)
    )
    z_clean <- finalize_nuisance_matrix(
      z_mat = ar_r$Y,
      valid_r = valid_r,
      impute = imp,
      impute_censored = imp_cens,
      center = ctr,
      scale = scl,
      drop_constant = drp
    )
    conf_obj$nuisance[[rn]] <- z_clean
    cleaned_flags[i] <- isTRUE(ar_r$applied)
    if (!is.null(conf_obj$summary) && nrow(conf_obj$summary) >= i) {
      conf_obj$summary$n_nuisance_cols[i] <- as.integer(ncol(z_clean))
      conf_obj$summary$aroma_cleaned[i] <- cleaned_flags[i]
    }
  }
  conf_obj$settings$aroma <- list(
    applied = any(cleaned_flags),
    mode = mode,
    fit_valid_rows = fit_valid_rows
  )
  conf_obj
}

#' Apply non-aggressive (or aggressive) ICA-AROMA denoising to MEDuSA time-series or confound objects
#'
#' Removes ICA-AROMA motion/artifact components in memory using vectorized
#' QR least-squares projection across all voxels simultaneously, without writing
#' intermediate 4D NIfTIs to disk or breaking MEDuSA's continuous acquisition
#' clock.
#'
#' Under \code{mode = "nonaggressive"} (the default; Pruim et al., 2015;
#' matching FSL \code{fsl_regfilt} and \code{BrainGnomes}), all accepted and
#' rejected MELODIC components are fit jointly with an intercept
#' (\eqn{Y = \mathbf{1}\mu^\top + M_{\mathrm{acc}} B_{\mathrm{acc}} + M_{\mathrm{rej}} B_{\mathrm{rej}} + E}),
#' and \emph{only} the partial contribution of rejected components
#' (\eqn{\tilde{M}_{\mathrm{rej}} \hat{B}_{\mathrm{rej}}}) is subtracted. Shared
#' variance between task-evoked neural components (\eqn{M_{\mathrm{acc}}}) and
#' task-correlated motion components (\eqn{M_{\mathrm{rej}}}) is therefore
#' retained rather than attenuated, and temporal degrees of freedom in Step 1
#' (\code{\link{estimate_crossrun_hrf}}) and Step 2
#' (\code{\link{deconvolve_crossrun_hrf}}) are not consumed by dozens of ICA
#' columns. When \code{clean_nuisance = TRUE}, attached nuisance regressors
#' (\code{$nuisance} / \code{$confounds}) receive the identical AROMA projection
#' so downstream joint nuisance regression cannot re-introduce AROMA motion
#' variance (Lindquist et al., 2019).
#'
#' @param x Target object to denoise: a \code{"medusa_nifti_ts"} object (updates
#'   \code{$bold} and, when \code{clean_nuisance = TRUE}, attached
#'   \code{$nuisance}/\code{$confounds}), a \code{"medusa_confounds"} object
#'   (updates \code{$nuisance}), a run-named list of \code{T x V} numeric BOLD
#'   matrices, or a single \code{T x V} numeric matrix.
#' @param aroma Optional \code{"medusa_aroma"} object from
#'   \code{\link{read_fmriprep_aroma}}. Required unless \code{mixing} is supplied
#'   directly.
#' @param mixing Optional MELODIC mixing file path(s) or \code{T x K} matrix /
#'   list of matrices passed to \code{\link{read_fmriprep_aroma}} when
#'   \code{aroma} is \code{NULL}.
#' @param classification Optional component classification path(s), metrics
#'   \code{data.frame}, logical mask, or 1-based integer vector of rejected ICs
#'   passed to \code{\link{read_fmriprep_aroma}} when \code{aroma} is
#'   \code{NULL}.
#' @param mode Denoising projection mode: \code{"nonaggressive"} (default;
#'   Pruim et al., 2015 partial regression protecting accepted signal ICs) or
#'   \code{"aggressive"} (full orthogonal projection of rejected ICs only).
#' @param clean_nuisance Logical; if \code{TRUE} (default) and \code{x} is a
#'   \code{"medusa_nifti_ts"} object with \code{$confounds} or \code{$nuisance}
#'   attached, applies the identical AROMA projection to nuisance columns so
#'   downstream joint regression cannot re-introduce removed motion variance.
#' @param fit_valid_rows Logical; if \code{TRUE}, fits AROMA component weights
#'   using only uncensored \code{valid_rows} while projecting across all rows.
#'   Default \code{FALSE} fits across full continuous MELODIC run series
#'   (matching FSL \code{fsl_regfilt}, Pruim et al., 2015, and
#'   \code{BrainGnomes}).
#' @param valid_rows Optional logical vector or run-named list of logical masks
#'   used when \code{fit_valid_rows = TRUE} (defaults to \code{x$valid_rows}
#'   when present on \code{x}).
#' @param preserve_mean Logical; if \code{TRUE} (default), preserves each
#'   voxel's temporal mean across fitted timepoints. Ignored for
#'   \code{"medusa_confounds"} inputs (which follow \code{x$settings$center}).
#' @param force Logical; if \code{TRUE}, re-applies AROMA even when metadata or
#'   \code{BrainGnomes} provenance indicates AROMA was already executed. Default
#'   \code{FALSE}.
#'
#' @return Object of the same class and structure as \code{x} with denoised
#'   signal matrices and AROMA diagnostic metadata attached (\code{x$aroma} on
#'   \code{"medusa_nifti_ts"}, \code{x$settings$aroma} on
#'   \code{"medusa_confounds"}, or \code{attr(x, "aroma")} on plain matrices /
#'   lists).
#'
#' @seealso \code{\link{read_fmriprep_aroma}}, \code{\link{read_fmriprep_subject}},
#'   \code{\link{read_fmriprep_confounds}}, \code{\link{estimate_crossrun_hrf}}
#'
#' @examples
#' set.seed(1)
#' s_task <- sin(seq(0, 6 * pi, length.out = 80))
#' m_motion <- 0.5 * s_task + rnorm(80, sd = 0.6)
#' mix <- cbind(IC_01 = s_task, IC_02 = m_motion, IC_03 = rnorm(80))
#' Y <- 500 + cbind(v1 = 2.0 * s_task + 1.5 * m_motion, v2 = 1.0 * s_task - 0.8 * m_motion) +
#'   matrix(rnorm(160, sd = 0.05), 80, 2)
#'
#' ar <- read_fmriprep_aroma(mix, classification = 2L)
#' Y_clean <- apply_aroma(Y, ar, mode = "nonaggressive")
#' attr(Y_clean, "aroma")$summary
#'
#' @export
apply_aroma <- function(
  x,
  aroma = NULL,
  mixing = NULL,
  classification = NULL,
  mode = c("nonaggressive", "aggressive"),
  clean_nuisance = TRUE,
  fit_valid_rows = FALSE,
  valid_rows = NULL,
  preserve_mean = TRUE,
  force = FALSE
) {
  mode <- match.arg(mode)
  checkmate::assert_flag(clean_nuisance)
  checkmate::assert_flag(fit_valid_rows)
  checkmate::assert_flag(preserve_mean)
  checkmate::assert_flag(force)

  if (is.null(aroma)) {
    if (is.null(mixing)) {
      stop("Supply either 'aroma' (a 'medusa_aroma' object) or 'mixing' (plus optional 'classification').", call. = FALSE)
    }
    bold_ref <- if (inherits(x, "medusa_confounds")) names(x$nuisance) else x
    aroma <- read_fmriprep_aroma(mixing = mixing, classification = classification, bold = bold_ref)
  } else if (!inherits(aroma, "medusa_aroma")) {
    stop("'aroma' must be a 'medusa_aroma' object produced by read_fmriprep_aroma().", call. = FALSE)
  }

  # 1. Dispatch on medusa_confounds
  if (inherits(x, "medusa_confounds")) {
    return(apply_aroma_to_confounds_obj(
      conf_obj = x,
      aroma_obj = aroma,
      mode = mode,
      fit_valid_rows = fit_valid_rows,
      force = force
    ))
  }

  resolve_valid_for_runs <- function(vr, runs, n_t_map) {
    if (!isTRUE(fit_valid_rows)) {
      return(stats::setNames(vector("list", length(runs)), runs))
    }
    if (is.null(vr)) {
      stop("fit_valid_rows = TRUE requires 'valid_rows' (or x$valid_rows on a 'medusa_nifti_ts' object).", call. = FALSE)
    }
    if (inherits(vr, "medusa_confounds")) vr <- vr$valid_rows
    if (is.logical(vr) && !is.list(vr)) {
      if (length(runs) != 1L) stop("Multi-run AROMA fitting requires a per-run list for 'valid_rows'.", call. = FALSE)
      vr <- stats::setNames(list(vr), runs)
    }
    al <- align_inputs_to_runs(
      vr,
      list(run_names = runs, n_time = n_t_map, bold_files = NULL),
      input_label = "valid_rows",
      is_paths = FALSE
    )
    al$items
  }

  # 2. Dispatch on medusa_nifti_ts
  if (inherits(x, "medusa_nifti_ts")) {
    if (!isTRUE(force)) {
      if (is.list(x$aroma) && isTRUE(x$aroma$applied)) {
        warning("BOLD data in 'x' was already denoised with apply_aroma(); skipping duplicate AROMA regression (set force = TRUE to override).", call. = FALSE)
        return(x)
      }
      prov <- if (is.list(x$bids) && is.list(x$bids$postproc_provenance)) {
        x$bids$postproc_provenance
      } else if (!is.null(x$run_info) && "file" %in% names(x$run_info)) {
        lapply(x$run_info$file, inspect_braingnomes_provenance)
      } else {
        list()
      }
      if (length(prov) > 0L && any(vapply(prov, function(p) isTRUE(p$apply_aroma), logical(1L)))) {
        warning("BrainGnomes provenance indicates 'apply_aroma' was already executed on BOLD input(s); skipping duplicate AROMA regression (set force = TRUE to override).", call. = FALSE)
        return(x)
      }
    }

    runs <- names(x$bold)
    vox_by_time <- identical(x$orientation, "voxel_by_time")
    n_t <- stats::setNames(vapply(x$bold, if (vox_by_time) ncol else nrow, integer(1L)), runs)
    aroma_al <- align_medusa_aroma_to_runs(aroma, runs, n_time = n_t)
    vr_use <- if (!is.null(valid_rows)) valid_rows else x$valid_rows
    vr_list <- resolve_valid_for_runs(vr_use, runs, n_t)

    sum_rows <- vector("list", length(runs))
    any_applied <- FALSE
    for (i in seq_along(runs)) {
      rn <- runs[i]
      mat_in <- if (vox_by_time) t(x$bold[[rn]]) else x$bold[[rn]]
      res_r <- aroma_regress_matrix(
        Y = mat_in,
        M = aroma_al$mixing[[rn]],
        noise_ics = aroma_al$noise_ics[[rn]],
        include_rows = vr_list[[rn]],
        mode = mode,
        preserve_mean = preserve_mean,
        run_label = sprintf("bold[['%s']]", rn)
      )
      x$bold[[rn]] <- if (vox_by_time) t(res_r$Y) else res_r$Y
      if (isTRUE(res_r$applied)) any_applied <- TRUE
      sum_rows[[i]] <- data.frame(
        run = rn,
        mode = mode,
        n_volumes = as.integer(n_t[[rn]]),
        n_fit_volumes = as.integer(res_r$n_fit),
        n_components = as.integer(res_r$n_components),
        n_accepted = as.integer(res_r$n_accepted),
        n_rejected = as.integer(res_r$n_rejected),
        applied = isTRUE(res_r$applied),
        mean_pct_var_removed = res_r$mean_pct_var_removed,
        median_pct_var_removed = res_r$median_pct_var_removed,
        stringsAsFactors = FALSE
      )
    }

    nuis_cleaned <- FALSE
    if (isTRUE(clean_nuisance)) {
      if (!is.null(x$confounds) && inherits(x$confounds, "medusa_confounds")) {
        x$confounds <- apply_aroma_to_confounds_obj(
          x$confounds,
          aroma_al,
          mode = mode,
          fit_valid_rows = fit_valid_rows,
          force = force
        )
        x$nuisance <- x$confounds$nuisance
        nuis_cleaned <- isTRUE(x$confounds$settings$aroma$applied)
      } else if (!is.null(x$nuisance) && is.list(x$nuisance)) {
        for (i in seq_along(runs)) {
          rn <- runs[i]
          z_r <- x$nuisance[[rn]]
          if (is.matrix(z_r) && ncol(z_r) > 0L) {
            ar_z <- aroma_regress_matrix(
              Y = z_r,
              M = aroma_al$mixing[[rn]],
              noise_ics = aroma_al$noise_ics[[rn]],
              include_rows = vr_list[[rn]],
              mode = mode,
              preserve_mean = TRUE,
              run_label = sprintf("nuisance[['%s']]", rn)
            )
            v_mask <- if (!is.null(x$valid_rows) && !is.null(x$valid_rows[[rn]])) x$valid_rows[[rn]] else rep(TRUE, nrow(z_r))
            x$nuisance[[rn]] <- finalize_nuisance_matrix(
              ar_z$Y,
              valid_r = v_mask,
              impute = "mean",
              impute_censored = TRUE,
              center = TRUE,
              scale = FALSE,
              drop_constant = TRUE
            )
            if (isTRUE(ar_z$applied)) nuis_cleaned <- TRUE
          }
        }
      }
    }

    sum_df <- do.call(rbind, sum_rows)
    rownames(sum_df) <- NULL
    x$aroma <- list(
      applied = any_applied,
      mode = mode,
      clean_nuisance = clean_nuisance,
      nuisance_cleaned = nuis_cleaned,
      fit_valid_rows = fit_valid_rows,
      preserve_mean = preserve_mean,
      summary = sum_df,
      components = aroma_al
    )
    return(x)
  }

  # 3. Dispatch on single matrix or list of run matrices
  if (!isTRUE(force) && is.list(attr(x, "aroma", exact = TRUE)) && isTRUE(attr(x, "aroma", exact = TRUE)$applied)) {
    warning("Input matrix/list was already denoised with apply_aroma(); skipping duplicate AROMA regression (set force = TRUE to override).", call. = FALSE)
    return(x)
  }

  single_mat <- is.matrix(x) && is.numeric(x)
  if (single_mat) {
    vox_by_time <- identical(attr(x, "orientation", exact = TRUE), "voxel_by_time")
    runs <- "run1"
    mat_list <- list(run1 = x)
  } else if (is.list(x) && !is.data.frame(x) && length(x) >= 1L && all(vapply(x, function(m) is.matrix(m) && is.numeric(m), logical(1L)))) {
    vox_by_time <- FALSE
    runs <- names(x)
    if (is.null(runs) || any(!nzchar(runs))) runs <- paste0("run", seq_along(x))
    mat_list <- stats::setNames(x, runs)
  } else {
    stop("'x' must be a 'medusa_nifti_ts' object, a 'medusa_confounds' object, a numeric T x V matrix, or a named list of numeric matrices.", call. = FALSE)
  }

  n_t <- stats::setNames(vapply(mat_list, if (vox_by_time) ncol else nrow, integer(1L)), runs)
  aroma_al <- align_medusa_aroma_to_runs(aroma, runs, n_time = n_t)
  vr_list <- resolve_valid_for_runs(valid_rows, runs, n_t)

  out_list <- stats::setNames(vector("list", length(runs)), runs)
  sum_rows <- vector("list", length(runs))
  any_applied <- FALSE
  for (i in seq_along(runs)) {
    rn <- runs[i]
    m_orig <- mat_list[[rn]]
    m_in <- if (vox_by_time) t(m_orig) else m_orig
    res_r <- aroma_regress_matrix(
      Y = m_in,
      M = aroma_al$mixing[[rn]],
      noise_ics = aroma_al$noise_ics[[rn]],
      include_rows = vr_list[[rn]],
      mode = mode,
      preserve_mean = preserve_mean,
      run_label = sprintf("bold[['%s']]", rn)
    )
    m_clean <- if (vox_by_time) t(res_r$Y) else res_r$Y
    # Preserve custom matrix attributes (e.g. read_nifti_ts(..., drop = TRUE) geometry/voxel_map/tr)
    orig_attrs <- attributes(m_orig)
    if (!is.null(orig_attrs)) {
      extra_nms <- setdiff(names(orig_attrs), c("dim", "dimnames"))
      for (an in extra_nms) attr(m_clean, an) <- orig_attrs[[an]]
    }
    out_list[[rn]] <- m_clean
    if (isTRUE(res_r$applied)) any_applied <- TRUE
    sum_rows[[i]] <- data.frame(
      run = rn,
      mode = mode,
      n_volumes = as.integer(n_t[[rn]]),
      n_fit_volumes = as.integer(res_r$n_fit),
      n_components = as.integer(res_r$n_components),
      n_accepted = as.integer(res_r$n_accepted),
      n_rejected = as.integer(res_r$n_rejected),
      applied = isTRUE(res_r$applied),
      mean_pct_var_removed = res_r$mean_pct_var_removed,
      median_pct_var_removed = res_r$median_pct_var_removed,
      stringsAsFactors = FALSE
    )
  }

  sum_df <- do.call(rbind, sum_rows)
  rownames(sum_df) <- NULL
  aroma_attr <- list(
    applied = any_applied,
    mode = mode,
    fit_valid_rows = fit_valid_rows,
    preserve_mean = preserve_mean,
    summary = sum_df
  )
  if (single_mat) {
    out_mat <- out_list[[1L]]
    attr(out_mat, "aroma") <- aroma_attr
    return(out_mat)
  }
  attr(out_list, "aroma") <- aroma_attr
  out_list
}

resolve_bids_subject_dir <- function(root_dir, sub_id, label = "fmriprep_dir") {
  checkmate::assert_directory_exists(root_dir)
  sub_tag <- paste0("sub-", sub_id)
  norm_root <- normalizePath(root_dir, winslash = "/", mustWork = TRUE)
  if (identical(basename(norm_root), sub_tag)) {
    return(norm_root)
  }
  for (cand in c(
    file.path(root_dir, sub_tag),
    file.path(root_dir, "fmriprep", sub_tag),
    file.path(root_dir, "postproc", sub_tag)
  )) {
    if (dir.exists(cand)) return(normalizePath(cand, winslash = "/", mustWork = TRUE))
  }
  stop(sprintf("Could not locate subject directory '%s' inside %s ('%s').", sub_tag, label, root_dir), call. = FALSE)
}

find_matching_bids_sidecar <- function(bold_row, candidates_tbl, require_space = FALSE, allow_optional_space = TRUE) {
  if (nrow(candidates_tbl) == 0L) return(NA_character_)
  keep <- rep(TRUE, nrow(candidates_tbl))
  for (field in c("sub", "ses", "task", "acq", "ce", "rec", "dir", "echo")) {
    b_val <- bold_row[[field]]
    c_val <- candidates_tbl[[field]]
    if (!is.na(b_val) && nzchar(b_val)) {
      keep <- keep & !is.na(c_val) & tolower(c_val) == tolower(b_val)
    } else {
      keep <- keep & (is.na(c_val) | !nzchar(c_val))
    }
  }
  # Canonical run match (run-01 <-> run-1, or both NA)
  b_run <- canonical_run_key(bold_row$run)
  c_run <- canonical_run_key(candidates_tbl$run)
  if (!is.na(b_run)) {
    keep <- keep & !is.na(c_run) & (c_run == b_run)
  } else if (!is.na(bold_row$run) && nzchar(bold_row$run)) {
    keep <- keep & !is.na(candidates_tbl$run) & (tolower(candidates_tbl$run) == tolower(bold_row$run))
  } else {
    keep <- keep & (is.na(candidates_tbl$run) | !nzchar(candidates_tbl$run))
  }

  spatial_fields <- c("space", "cohort", "res", "den", "hemi")
  if (isTRUE(require_space)) {
    for (field in spatial_fields) {
      if (!(field %in% names(candidates_tbl))) next
      b_val <- bold_row[[field]]
      c_val <- candidates_tbl[[field]]
      if (!is.na(b_val) && nzchar(b_val)) {
        keep <- keep & !is.na(c_val) & (c_val == b_val)
      } else {
        keep <- keep & (is.na(c_val) | !nzchar(c_val))
      }
    }
  } else if (isTRUE(allow_optional_space)) {
    # Sidecars with explicit space/cohort/res/den/hemi (e.g. BrainGnomes postproc confounds/censor)
    # must match BOLD when present, while space-agnostic sidecars (e.g. fMRIPrep confounds/BIDS events) match any space
    for (field in spatial_fields) {
      if (!(field %in% names(candidates_tbl))) next
      b_val <- bold_row[[field]]
      c_val <- candidates_tbl[[field]]
      has_c <- !is.na(c_val) & nzchar(c_val)
      if (!is.na(b_val) && nzchar(b_val)) {
        keep <- keep & (!has_c | (c_val == b_val))
      } else {
        keep <- keep & !has_c
      }
    }
  }

  hits <- candidates_tbl$file[which(keep)]
  if (length(hits) == 0L) return(NA_character_)
  if (length(hits) > 1L) {
    # Prefer .nii.gz / .tsv over uncompressed/compressed duplicates of the same stem
    stems <- strip_neuro_ext(hits)
    if (length(unique(stems)) == 1L) return(hits[1L])
    stop(sprintf("Multiple sidecar files matched BOLD run '%s': %s", basename(bold_row$file), paste(basename(hits), collapse = ", ")), call. = FALSE)
  }
  hits[1L]
}

find_matching_confounds_sidecar <- function(bold_row, ctbl, upstream_conf = NA_character_) {
  b_desc <- if (!is.na(bold_row$desc) && nzchar(bold_row$desc)) bold_row$desc else NA_character_
  is_std_fprep <- (is.na(ctbl$desc) & tolower(ctbl$suffix) == "confounds") |
    (!is.na(ctbl$desc) & tolower(ctbl$desc) == "confounds")

  if (nrow(ctbl) > 0L) {
    if (!is.na(b_desc) && !identical(tolower(b_desc), "preproc")) {
      # Postprocessed BOLD: prefer desc-<b_desc>_confounds.tsv first, then standard fMRIPrep confounds
      exact_desc <- !is.na(ctbl$desc) & (ctbl$desc == b_desc) & tolower(ctbl$suffix) == "confounds"
      if (any(exact_desc)) {
        hit <- find_matching_bids_sidecar(bold_row, ctbl[exact_desc, , drop = FALSE], require_space = FALSE)
        if (!is.na(hit)) return(hit)
      }
      if (any(is_std_fprep)) {
        hit <- find_matching_bids_sidecar(bold_row, ctbl[is_std_fprep, , drop = FALSE], require_space = FALSE)
        if (!is.na(hit)) return(hit)
      }
    } else {
      # Raw fMRIPrep BOLD: prefer standard desc-confounds_timeseries.tsv over postprocessed confounds
      if (any(is_std_fprep)) {
        hit <- find_matching_bids_sidecar(bold_row, ctbl[is_std_fprep, , drop = FALSE], require_space = FALSE)
        if (!is.na(hit)) return(hit)
      }
      hit <- find_matching_bids_sidecar(bold_row, ctbl, require_space = FALSE)
      if (!is.na(hit)) return(hit)
    }
  }

  if (is.character(upstream_conf) && length(upstream_conf) == 1L && !is.na(upstream_conf) && file.exists(upstream_conf)) {
    return(upstream_conf)
  }
  NA_character_
}

find_matching_censor_sidecar <- function(bold_row, censor_tbl) {
  if (nrow(censor_tbl) == 0L) return(NA_character_)
  b_desc <- if (!is.na(bold_row$desc) && nzchar(bold_row$desc)) bold_row$desc else NA_character_
  if (!is.na(b_desc) && !identical(tolower(b_desc), "preproc")) {
    exact_desc <- !is.na(censor_tbl$desc) & (censor_tbl$desc == b_desc)
    if (any(exact_desc)) {
      hit <- find_matching_bids_sidecar(bold_row, censor_tbl[exact_desc, , drop = FALSE], require_space = FALSE)
      if (!is.na(hit)) return(hit)
    }
    undesc <- is.na(censor_tbl$desc) | !nzchar(censor_tbl$desc)
    if (any(undesc)) {
      return(find_matching_bids_sidecar(bold_row, censor_tbl[undesc, , drop = FALSE], require_space = FALSE))
    }
    return(NA_character_)
  }
  allowed <- is.na(censor_tbl$desc) | !nzchar(censor_tbl$desc) |
    (!is.na(b_desc) & tolower(censor_tbl$desc) == tolower(b_desc))
  if (!any(allowed)) return(NA_character_)
  find_matching_bids_sidecar(bold_row, censor_tbl[allowed, , drop = FALSE], require_space = FALSE)
}

find_inherited_bids_events <- function(bold_row, bids_dir, sub_dir_bids = NULL) {
  # BIDS Inheritance Principle fallback: subject/session task-level or root task-level *_events.tsv
  search_dirs <- character(0L)
  if (!is.null(sub_dir_bids) && dir.exists(sub_dir_bids)) {
    if (!is.na(bold_row$ses) && nzchar(bold_row$ses)) {
      ses_d <- file.path(sub_dir_bids, paste0("ses-", bold_row$ses), "func")
      if (dir.exists(ses_d)) search_dirs <- c(search_dirs, ses_d)
      ses_root <- file.path(sub_dir_bids, paste0("ses-", bold_row$ses))
      if (dir.exists(ses_root)) search_dirs <- c(search_dirs, ses_root)
    }
    func_d <- file.path(sub_dir_bids, "func")
    if (dir.exists(func_d)) search_dirs <- c(search_dirs, func_d)
    search_dirs <- c(search_dirs, sub_dir_bids)
  }
  if (!is.null(bids_dir) && dir.exists(bids_dir)) {
    search_dirs <- c(search_dirs, bids_dir)
  }
  search_dirs <- unique(search_dirs)
  for (d in search_dirs) {
    ev_files <- list.files(d, pattern = "_events\\.tsv(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
    if (!length(ev_files)) next
    etbl <- parse_bids_entity_table(ev_files)
    #Exact run match first
    exact_hit <- find_matching_bids_sidecar(bold_row, etbl, require_space = FALSE)
    if (!is.na(exact_hit)) return(exact_hit)
    # Inherited run-agnostic task-level events file
    inh_keep <- is.na(etbl$run) &
      (is.na(etbl$sub) | etbl$sub == bold_row$sub) &
      (is.na(etbl$ses) | (!is.na(bold_row$ses) && etbl$ses == bold_row$ses)) &
      (!is.na(etbl$task) & !is.na(bold_row$task) & tolower(etbl$task) == tolower(bold_row$task))
    if (sum(inh_keep) == 1L) return(etbl$file[which(inh_keep)])
  }
  NA_character_
}

combine_fmriprep_masks <- function(mask_files, mode = c("intersection", "union", "first"), threshold = 0) {
  mode <- match.arg(mode)
  if (length(mask_files) == 1L || identical(mode, "first")) {
    return(mask_files[1L])
  }
  ref_img <- RNifti::readNifti(mask_files[1L])
  ref_geom <- resolve_nifti_spatial_geometry(mask_files[1L], label = basename(mask_files[1L]))
  acc <- is.finite(ref_img) & (as.numeric(ref_img) > threshold)
  for (i in seq.int(2L, length(mask_files))) {
    mi_info <- read_3d_roi_values(mask_files[i], label = basename(mask_files[i]))
    validate_nifti_spatial_match(ref_geom, mi_info$geometry, basename(mask_files[1L]), basename(mask_files[i]))
    mi_on <- is.finite(mi_info$values) & (mi_info$values > threshold)
    acc <- if (identical(mode, "intersection")) (acc & mi_on) else (acc | mi_on)
  }
  arr <- array(as.integer(acc), dim = ref_geom$dim)
  RNifti::updateNifti(arr, template = ref_img)
}

collect_subject_func_dirs <- function(sub_dir_path) {
  non_func_modalities <- c("anat", "dwi", "fmap", "perf", "meg", "eeg", "ieeg", "beh", "figures", "log", "logs", "work", "scratch")
  f_dirs <- character(0L)
  direct_f <- file.path(sub_dir_path, "func")
  if (dir.exists(direct_f)) f_dirs <- c(f_dirs, direct_f)
  ch_dirs <- list.dirs(sub_dir_path, full.names = TRUE, recursive = FALSE)
  ses_d <- ch_dirs[grepl("^ses-", basename(ch_dirs))]
  for (sd in ses_d) {
    fd_p <- file.path(sd, "func")
    if (dir.exists(fd_p)) f_dirs <- c(f_dirs, fd_p)
  }
  flat_d <- c(
    sub_dir_path,
    ch_dirs[!(tolower(basename(ch_dirs)) %in% non_func_modalities) & !grepl("\\.provenance$", basename(ch_dirs), ignore.case = TRUE)]
  )
  unique(c(f_dirs, flat_d))
}

discover_aroma_sidecars_for_subject <- function(btbl, sub_id, sub_dir, func_dirs, prov_list = list(), aroma_dir = NULL) {
  sub_tag <- paste0("sub-", sub_id)
  dir_groups <- list()

  if (!is.null(aroma_dir)) {
    checkmate::assert_directory_exists(aroma_dir)
    norm_ad <- normalizePath(aroma_dir, winslash = "/", mustWork = TRUE)
    ad_sub <- NULL
    if (identical(basename(norm_ad), sub_tag)) {
      ad_sub <- norm_ad
    } else {
      for (cand in c(
        file.path(norm_ad, sub_tag),
        file.path(norm_ad, "fmripost_aroma", sub_tag),
        file.path(norm_ad, "fmripost-aroma", sub_tag),
        file.path(norm_ad, "aroma", sub_tag)
      )) {
        if (dir.exists(cand)) {
          ad_sub <- normalizePath(cand, winslash = "/", mustWork = TRUE)
          break
        }
      }
    }
    dir_groups[[length(dir_groups) + 1L]] <- if (!is.null(ad_sub)) collect_subject_func_dirs(ad_sub) else norm_ad
  }

  dir_groups[[length(dir_groups) + 1L]] <- func_dirs

  if (is.null(aroma_dir)) {
    p1 <- dirname(sub_dir)
    p2 <- dirname(p1)
    sib_dirs <- character(0L)
    for (root_p in unique(c(p2, p1))) {
      for (dname in c("fmripost_aroma", "fmripost-aroma", "aroma", "fmriprep")) {
        cand_s <- file.path(root_p, dname, sub_tag)
        if (dir.exists(cand_s)) {
          cand_norm <- normalizePath(cand_s, winslash = "/", mustWork = TRUE)
          if (!identical(cand_norm, sub_dir)) {
            sib_dirs <- c(sib_dirs, collect_subject_func_dirs(cand_norm))
          }
        }
      }
    }
    if (length(sib_dirs)) {
      dir_groups[[length(dir_groups) + 1L]] <- unique(sib_dirs)
    }
  }

  n_runs <- nrow(btbl)
  mix_hits <- rep(NA_character_, n_runs)
  cls_hits <- rep(NA_character_, n_runs)

  for (grp in dir_groups) {
    if (!length(grp)) next
    mix_cands <- unique(unlist(lapply(grp, function(d) {
      list.files(d, pattern = "_(desc-melodic_mixing|melodic_mixing)\\.tsv(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
    }), use.names = FALSE))
    met_cands <- unique(unlist(lapply(grp, function(d) {
      list.files(d, pattern = "_desc-aroma_metrics\\.tsv(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
    }), use.names = FALSE))
    csv_cands <- unique(unlist(lapply(grp, function(d) {
      list.files(d, pattern = "_AROMAnoiseICs\\.(csv|txt)(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
    }), use.names = FALSE))

    mix_tbl <- if (length(mix_cands)) parse_bids_entity_table(mix_cands) else parse_bids_entity_table(character(0L))
    met_tbl <- if (length(met_cands)) parse_bids_entity_table(met_cands) else parse_bids_entity_table(character(0L))
    csv_tbl <- if (length(csv_cands)) parse_bids_entity_table(csv_cands) else parse_bids_entity_table(character(0L))

    for (i in seq_len(n_runs)) {
      brow <- btbl[i, , drop = FALSE]
      if (is.na(mix_hits[i]) && nrow(mix_tbl) > 0L) {
        mix_hits[i] <- find_matching_aroma_sidecar(brow, mix_tbl)
      }
      if (is.na(cls_hits[i])) {
        if (nrow(met_tbl) > 0L) {
          cls_hits[i] <- find_matching_aroma_sidecar(brow, met_tbl)
        }
        if (is.na(cls_hits[i]) && nrow(csv_tbl) > 0L) {
          cls_hits[i] <- find_matching_aroma_sidecar(brow, csv_tbl)
        }
      }
    }
    if (all(!is.na(mix_hits)) && all(!is.na(cls_hits))) break
  }

  for (i in seq_len(n_runs)) {
    if (length(prov_list) >= i) {
      pm <- prov_list[[i]]$melodic_mix
      pa <- prov_list[[i]]$aroma_metrics
      if (is.na(mix_hits[i]) && is.character(pm) && length(pm) == 1L && !is.na(pm) && file.exists(pm)) {
        mix_hits[i] <- pm
      }
      if (is.na(cls_hits[i]) && is.character(pa) && length(pa) == 1L && !is.na(pa) && file.exists(pa)) {
        cls_hits[i] <- pa
      }
    }
    if (is.na(cls_hits[i]) && !is.na(mix_hits[i])) {
      cls_hits[i] <- discover_aroma_classification_for_mixing(mix_hits[i])
    }
  }

  list(mixing_files = mix_hits, classification_files = cls_hits)
}

#' Discover and ingest an fMRIPrep or postprocessed participant into a unified MEDuSA time-series bundle
#'
#' Locates a participant's preprocessed or \code{BrainGnomes}-postprocessed BOLD
#' runs (\code{*_desc-preproc_bold.nii.gz} or \code{*_desc-<desc>_bold.nii.gz}),
#' brain or template masks (\code{*_desc-brain_mask.nii.gz},
#' \code{*_templatemask.nii.gz}), confound tables
#' (\code{*_desc-confounds_timeseries.tsv}, \code{*_desc-<desc>_confounds.tsv}),
#' optional \code{fMRIPost-AROMA} / legacy ICA-AROMA mixing and classification
#' sidecars, optional \code{*_censor.1D} masks, and optional BIDS event tables
#' (\code{*_events.tsv}), returning a \code{"medusa_nifti_ts"} object populated
#' with C++-extracted \code{T x V} \code{$bold} matrices, \code{$voxel_map},
#' \code{$parcel_ids}, \code{$tr}, \code{$nuisance}, \code{$valid_rows},
#' \code{$confounds}, and \code{$events} in locked run order, ready to pass
#' directly into \code{\link{estimate_crossrun_hrf}} and
#' \code{\link{deconvolve_crossrun_hrf}}.
#'
#' @param fmriprep_dir Path to the fMRIPrep or \code{BrainGnomes} postprocessing
#'   derivatives root directory (or directly to \code{sub-<id>}).
#' @param sub Participant identifier (with or without \code{"sub-"} prefix, e.g.,
#'   \code{"01"} or \code{"sub-01"}).
#' @param bids_dir Optional path to the raw BIDS root directory containing
#'   \code{*_events.tsv} files (when \code{events = "auto"}).
#' @param ses Optional session filter (with or without \code{"ses-"} prefix).
#' @param task Optional BIDS task label (with or without \code{"task-"} prefix).
#'   Required when multiple tasks exist for the participant.
#' @param runs Optional vector of run identifiers to include (e.g.,
#'   \code{1:3} or \code{c("01", "02")}).
#' @param space Optional template/output space filter (e.g.,
#'   \code{"MNI152NLin2009cAsym"}). Auto-selected when only one space is present;
#'   required when multiple spaces exist.
#' @param cohort Optional TemplateFlow \code{cohort-*} entity filter (e.g.,
#'   \code{"1"} or \code{"2"} for \code{MNIPediatricAsym}). Auto-selected when
#'   unique; required when multiple cohorts exist.
#' @param res Optional BIDS resolution entity (\code{res-*}). Auto-selected when
#'   unique.
#' @param desc BOLD derivative \code{desc-*} entity (default \code{"preproc"} for
#'   fMRIPrep; use e.g. \code{"postproc"} for \code{BrainGnomes} postprocessed
#'   outputs, or \code{NULL} to match un-described BOLD files).
#' @param mask Either \code{"auto"} (default; discovers per-run
#'   \code{*_desc-brain_mask.nii.gz} or \code{*_templatemask.nii.gz} files and
#'   combines them via \code{mask_mode}, falling back to anatomical mask,
#'   provenance spatial mask, or \code{atlas}), an explicit 3D NIfTI path /
#'   \code{RNifti} image, or \code{NULL}.
#' @param mask_mode Rule for combining per-run brain masks when
#'   \code{mask = "auto"}: \code{"intersection"} (default; retains voxels in-mask
#'   across every run), \code{"union"}, or \code{"first"}.
#' @param atlas Optional 3D parcellation NIfTI path or \code{RNifti} image on the
#'   BOLD spatial grid.
#' @param atlas_subset Optional integer/numeric vector of parcel labels to retain.
#' @param mask_threshold Numeric threshold (\code{>}) applied to \code{mask}.
#'   Default \code{0}.
#' @param atlas_threshold Numeric threshold (\code{>}) applied to \code{atlas}.
#'   Default \code{0}.
#' @param tr Optional repetition time in seconds (inferred from NIfTI headers or
#'   \code{*_bold.json} \code{RepetitionTime} when \code{NULL}).
#' @param events Either \code{"auto"} (default; discovers matching
#'   \code{*_events.tsv} in \code{bids_dir} / \code{fmriprep_dir}), explicit
#'   event file paths / \code{data.frame}s passed to \code{\link{read_bids_events}},
#'   or \code{NULL} to skip event ingestion (e.g., resting-state).
#' @param condition_col,condition_map,default_condition,exclude_conditions,na_duration
#'   Event normalization options forwarded to \code{\link{read_bids_events}}.
#' @param aroma ICA-AROMA denoising policy applied in memory during ingestion:
#'   \itemize{
#'     \item \code{"none"} (default): no AROMA projection applied.
#'     \item \code{"nonaggressive"}: discovers \code{fMRIPost-AROMA}
#'       (\code{*_desc-melodic_mixing.tsv} + \code{*_desc-aroma_metrics.tsv}) or
#'       legacy fMRIPrep \code{<= 23.0} (\code{*_desc-MELODIC_mixing.tsv} +
#'       \code{*_AROMAnoiseICs.csv}) sidecars and applies Pruim et al. (2015)
#'       non-aggressive partial regression to \code{$bold} (and selected
#'       \code{$nuisance} columns when \code{aroma_clean_nuisance = TRUE}).
#'     \item \code{"auto"}: applies \code{"nonaggressive"} AROMA when matching
#'       AROMA sidecars exist for all runs and skips silently otherwise.
#'     \item \code{"aggressive"}: applies full orthogonal removal of rejected
#'       AROMA components without protecting accepted signal components.
#'     \item A pre-constructed \code{"medusa_aroma"} object from
#'       \code{\link{read_fmriprep_aroma}} (applied non-aggressively).
#'   }
#'   When \code{BrainGnomes} provenance indicates \code{apply_aroma} was already
#'   executed upstream on the BOLD derivative, duplicate AROMA regression on
#'   \code{$bold} is skipped automatically.
#' @param aroma_dir Optional path to an \code{fMRIPost-AROMA} derivatives root
#'   (or subject directory). When \code{NULL}, \code{read_fmriprep_subject()}
#'   searches \code{fmriprep_dir}, \code{BrainGnomes} provenance
#'   \code{SourceFiles}, and sibling \code{fmripost_aroma} /
#'   \code{fmripost-aroma} / \code{aroma} derivative directories automatically.
#' @param aroma_clean_nuisance Logical; if \code{TRUE} (default), applies the
#'   identical AROMA projection operator to selected \code{$nuisance} columns so
#'   downstream joint regression in \code{\link{estimate_crossrun_hrf}} and
#'   \code{\link{deconvolve_crossrun_hrf}} cannot re-introduce removed AROMA
#'   motion variance (Lindquist et al., 2019).
#' @param aroma_fit_valid Logical; if \code{TRUE}, estimates AROMA component
#'   weights using uncensored \code{valid_rows} only rather than full MELODIC
#'   timecourses (\code{FALSE}, default).
#' @param strategy,columns,column_pattern,n_compcor,compcor_type,include_cosines
#'   Nuisance selection options forwarded to \code{\link{read_fmriprep_confounds}}
#'   (including \code{"csf_wm"} and \code{"csf_wm_8p"} for pairing physiological
#'   regressors with \code{aroma = "nonaggressive"}). When \code{BrainGnomes}
#'   provenance indicates \code{confound_regression} was already applied to the
#'   BOLD derivative and \code{strategy} is not explicitly specified,
#'   \code{strategy} defaults automatically to \code{"none"} to avoid
#'   double-regressing unfiltered confounds.
#' @param censor Either \code{"auto"} (default; automatically uses matching
#'   \code{*_censor.1D} sidecars when present for postprocessed outputs), explicit
#'   censor file paths / vectors forwarded to \code{\link{read_fmriprep_confounds}},
#'   or \code{NULL} to ignore \code{*_censor.1D} sidecars.
#' @param fd_threshold,fd_col,head_radius,rot_units,motion_bandstop_bpm,motion_lowpass_hz,motion_filter_order,motion_filter_pad,motion_filter_target,dvars_threshold,dvars_col,include_non_steady_state,include_motion_outliers,drop_initial_volumes,censor_prev,censor_next,min_segment_length
#'   Clock-preserving censoring and optional respiratory motion filtering
#'   options forwarded to \code{\link{read_fmriprep_confounds}}.
#' @param id_style Voxel ID convention forwarded to \code{\link{read_nifti_ts}}
#'   (\code{"ijk"} or \code{"vnum"}).
#' @param ncpus Number of parallel workers forwarded to \code{\link{read_nifti_ts}}.
#'
#' @return An S3 object of class \code{"medusa_nifti_ts"} enriched with
#'   \code{$nuisance}, \code{$valid_rows}, \code{$confounds}
#'   (\code{"medusa_confounds"}), \code{$aroma} (when AROMA is applied),
#'   \code{$events} (when discovered/supplied), and \code{$bids} provenance
#'   metadata.
#'
#' @seealso \code{\link{read_nifti_ts}}, \code{\link{read_fmriprep_confounds}},
#'   \code{\link{read_fmriprep_aroma}}, \code{\link{apply_aroma}},
#'   \code{\link{read_bids_events}}, \code{\link{estimate_crossrun_hrf}},
#'   \code{\link{deconvolve_crossrun_hrf}}
#'
#' @examples
#' root <- tempfile("medusa_bids_demo_")
#' fprep_func <- file.path(root, "derivatives", "fmriprep", "sub-01", "func")
#' bids_func <- file.path(root, "bids", "sub-01", "func")
#' dir.create(fprep_func, recursive = TRUE)
#' dir.create(bids_func, recursive = TRUE)
#'
#' mask_img <- RNifti::asNifti(array(1L, dim = c(2L, 2L, 2L)))
#' RNifti::pixdim(mask_img) <- c(2, 2, 2)
#' for (r in 1:2) {
#'   stem <- sprintf("sub-01_task-mid_run-0%d", r)
#'   bold_arr <- array(rnorm(2L * 2L * 2L * 10L), dim = c(2L, 2L, 2L, 10L))
#'   bold_img <- RNifti::asNifti(bold_arr)
#'   RNifti::pixdim(bold_img) <- c(2, 2, 2, 1.0)
#'   RNifti::writeNifti(
#'     bold_img,
#'     file.path(fprep_func, paste0(stem, "_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"))
#'   )
#'   RNifti::writeNifti(
#'     mask_img,
#'     file.path(fprep_func, paste0(stem, "_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz"))
#'   )
#'   conf_df <- data.frame(
#'     trans_x = rnorm(10, sd = 0.02), trans_y = rnorm(10, sd = 0.02),
#'     trans_z = rnorm(10, sd = 0.02), rot_x = rnorm(10, sd = 0.002),
#'     rot_y = rnorm(10, sd = 0.002), rot_z = rnorm(10, sd = 0.002),
#'     framewise_displacement = c("n/a", rep("0.08", 9))
#'   )
#'   utils::write.table(
#'     conf_df,
#'     file.path(fprep_func, paste0(stem, "_desc-confounds_timeseries.tsv")),
#'     sep = "\t", row.names = FALSE, quote = FALSE
#'   )
#'   ev_df <- data.frame(onset = c(2, 6), duration = c(0.5, 0.5), trial_type = c("cue", "feedback"))
#'   utils::write.table(
#'     ev_df,
#'     file.path(bids_func, paste0(stem, "_events.tsv")),
#'     sep = "\t", row.names = FALSE, quote = FALSE
#'   )
#' }
#'
#' sub01 <- read_fmriprep_subject(
#'   fmriprep_dir = file.path(root, "derivatives", "fmriprep"),
#'   bids_dir = file.path(root, "bids"),
#'   sub = "01",
#'   task = "mid",
#'   strategy = "6p"
#' )
#' sub01
#' unlink(root, recursive = TRUE)
#'
#' @export
read_fmriprep_subject <- function(
  fmriprep_dir,
  sub,
  bids_dir = NULL,
  ses = NULL,
  task = NULL,
  runs = NULL,
  space = NULL,
  cohort = NULL,
  res = NULL,
  desc = "preproc",
  mask = "auto",
  mask_mode = c("intersection", "union", "first"),
  atlas = NULL,
  atlas_subset = NULL,
  mask_threshold = 0,
  atlas_threshold = 0,
  tr = NULL,
  events = "auto",
  condition_col = "trial_type",
  condition_map = NULL,
  default_condition = NULL,
  exclude_conditions = NULL,
  na_duration = 0,
  aroma = "none",
  aroma_dir = NULL,
  aroma_clean_nuisance = TRUE,
  aroma_fit_valid = FALSE,
  strategy = c("24p", "6p", "12p", "csf_wm", "csf_wm_8p", "24p_csf_wm", "27p", "36p", "compcor", "custom", "none"),
  columns = NULL,
  column_pattern = NULL,
  n_compcor = NULL,
  compcor_type = c("anatomical", "split", "temporal"),
  include_cosines = "auto",
  censor = "auto",
  fd_threshold = 0.5,
  fd_col = "auto",
  head_radius = 50,
  rot_units = c("rad", "deg"),
  motion_bandstop_bpm = NULL,
  motion_lowpass_hz = NULL,
  motion_filter_order = NULL,
  motion_filter_pad = c("constant", "odd", "even", "zero"),
  motion_filter_target = c("fd_only", "fd_and_nuisance"),
  dvars_threshold = NULL,
  dvars_col = "std_dvars",
  include_non_steady_state = TRUE,
  include_motion_outliers = FALSE,
  drop_initial_volumes = 0L,
  censor_prev = 0L,
  censor_next = 0L,
  min_segment_length = 1L,
  id_style = c("ijk", "vnum"),
  ncpus = 1L
) {
  strategy_supplied <- !missing(strategy)
  fd_threshold_supplied <- !missing(fd_threshold)
  mask_mode <- match.arg(mask_mode)
  strategy <- match.arg(strategy)
  compcor_type <- match.arg(compcor_type)
  rot_units <- match.arg(rot_units)
  motion_filter_pad <- match.arg(motion_filter_pad)
  motion_filter_target <- match.arg(motion_filter_target)
  id_style <- match.arg(id_style)
  checkmate::assert_flag(aroma_clean_nuisance)
  checkmate::assert_flag(aroma_fit_valid)
  if (!is.null(aroma_dir)) checkmate::assert_directory_exists(aroma_dir)
  if (inherits(aroma, "medusa_aroma")) {
    aroma_req <- "object"
    aroma_mode_eff <- "nonaggressive"
  } else if (is.character(aroma) && length(aroma) == 1L && !is.na(aroma)) {
    aroma_req <- match.arg(aroma, c("none", "nonaggressive", "auto", "aggressive"))
    aroma_mode_eff <- if (identical(aroma_req, "aggressive")) "aggressive" else "nonaggressive"
  } else {
    stop("'aroma' must be one of 'none', 'nonaggressive', 'auto', 'aggressive', or a 'medusa_aroma' object.", call. = FALSE)
  }
  checkmate::assert_string(sub, min.chars = 1L)
  sub_id <- sub("^sub-", "", trimws(sub))
  sub_dir <- resolve_bids_subject_dir(fmriprep_dir, sub_id, label = "fmriprep_dir")
  if (!is.null(bids_dir)) checkmate::assert_directory_exists(bids_dir)

  ses_clean <- if (!is.null(ses)) sub("^ses-", "", trimws(as.character(ses))) else NULL
  task_clean <- if (!is.null(task)) sub("^task-", "", trimws(as.character(task))) else NULL
  space_clean <- if (!is.null(space)) sub("^space-", "", trimws(as.character(space))) else NULL
  cohort_clean <- if (!is.null(cohort)) sub("^cohort-", "", trimws(as.character(cohort))) else NULL
  res_clean <- if (!is.null(res)) sub("^res-", "", trimws(as.character(res))) else NULL
  desc_clean <- if (!is.null(desc)) sub("^desc-", "", trimws(as.character(desc))) else NULL

  # Locate functional directories under sub_dir (standard BIDS func/ + flat BrainGnomes postproc layout)
  func_dirs <- character(0L)
  direct_func <- file.path(sub_dir, "func")
  if (dir.exists(direct_func)) func_dirs <- c(func_dirs, direct_func)
  child_dirs <- list.dirs(sub_dir, full.names = TRUE, recursive = FALSE)
  ses_dirs <- child_dirs[grepl("^ses-", basename(child_dirs))]
  for (sd in ses_dirs) {
    fd_path <- file.path(sd, "func")
    if (dir.exists(fd_path)) func_dirs <- c(func_dirs, fd_path)
  }

  # Fallback for BrainGnomes postproc directories that write directly into sub-XX/, sub-XX/ses-YY/, or stream folders
  non_func_modalities <- c("anat", "dwi", "fmap", "perf", "meg", "eeg", "ieeg", "beh", "figures", "log", "logs", "work", "scratch")
  flat_dirs <- c(
    sub_dir,
    child_dirs[!(tolower(basename(child_dirs)) %in% non_func_modalities) & !grepl("\\.provenance$", basename(child_dirs), ignore.case = TRUE)]
  )
  if (!length(func_dirs)) {
    func_dirs <- unique(flat_dirs)
  } else {
    has_bold_in_func <- any(vapply(func_dirs, function(d) {
      length(list.files(d, pattern = "_bold\\.nii(\\.gz)?$", ignore.case = TRUE)) > 0L
    }, logical(1L)))
    if (!has_bold_in_func) func_dirs <- unique(c(func_dirs, flat_dirs))
  }

  bold_candidates <- unique(unlist(lapply(func_dirs, function(d) {
    list.files(d, pattern = "_bold\\.nii(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE))
  if (!length(bold_candidates)) {
    stop(sprintf("No '*_bold.nii(.gz)' files found under '%s'.", sub_dir), call. = FALSE)
  }

  btbl <- parse_bids_entity_table(bold_candidates)
  if (!is.null(desc_clean)) {
    btbl <- btbl[!is.na(btbl$desc) & btbl$desc == desc_clean, , drop = FALSE]
  }
  if (!is.null(ses_clean)) {
    btbl <- btbl[!is.na(btbl$ses) & btbl$ses %in% ses_clean, , drop = FALSE]
  }
  if (!is.null(task_clean)) {
    btbl <- btbl[!is.na(btbl$task) & tolower(btbl$task) %in% tolower(task_clean), , drop = FALSE]
  } else {
    avail_tasks <- unique(btbl$task[!is.na(btbl$task) & nzchar(btbl$task)])
    if (length(avail_tasks) > 1L) {
      stop(
        sprintf(
          "Multiple BIDS tasks found for sub-%s (%s); specify 'task' explicitly.",
          sub_id, paste(avail_tasks, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(space_clean)) {
    btbl <- btbl[!is.na(btbl$space) & btbl$space == space_clean, , drop = FALSE]
  } else {
    avail_spaces <- unique(btbl$space[!is.na(btbl$space) & nzchar(btbl$space)])
    if (length(avail_spaces) > 1L) {
      stop(
        sprintf(
          "Multiple output spaces found for sub-%s (%s); specify 'space' explicitly.",
          sub_id, paste(avail_spaces, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(cohort_clean)) {
    btbl <- btbl[!is.na(btbl$cohort) & btbl$cohort == cohort_clean, , drop = FALSE]
  } else {
    avail_cohorts <- unique(btbl$cohort[!is.na(btbl$cohort) & nzchar(btbl$cohort)])
    if (length(avail_cohorts) > 1L) {
      stop(
        sprintf(
          "Multiple template cohorts found for sub-%s (%s); specify 'cohort' explicitly.",
          sub_id, paste(avail_cohorts, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(res_clean)) {
    btbl <- btbl[!is.na(btbl$res) & btbl$res == res_clean, , drop = FALSE]
  } else {
    avail_res <- unique(btbl$res[!is.na(btbl$res) & nzchar(btbl$res)])
    if (length(avail_res) > 1L) {
      stop(
        sprintf(
          "Multiple spatial resolutions found for sub-%s (%s); specify 'res' explicitly.",
          sub_id, paste(avail_res, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  avail_hemi <- unique(btbl$hemi[!is.na(btbl$hemi) & nzchar(btbl$hemi)])
  if (length(avail_hemi) > 1L) {
    stop(
      sprintf(
        "Multiple hemisphere ('hemi-') entities found for sub-%s (%s); ingest one hemisphere at a time.",
        sub_id, paste(avail_hemi, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.null(runs)) {
    req_canon <- canonical_run_key(runs)
    b_canon <- canonical_run_key(btbl$run)
    keep_r <- (!is.na(b_canon) & b_canon %in% req_canon) |
      (!is.na(btbl$run) & tolower(btbl$run) %in% tolower(as.character(runs)))
    btbl <- btbl[keep_r, , drop = FALSE]
  }

  if (nrow(btbl) == 0L) {
    stop(
      sprintf(
        "No BOLD runs matched filters for sub-%s (ses=%s, task=%s, space=%s, cohort=%s, desc=%s).",
        sub_id,
        if (is.null(ses_clean)) "*" else paste(ses_clean, collapse = ","),
        if (is.null(task_clean)) "*" else paste(task_clean, collapse = ","),
        if (is.null(space_clean)) "*" else space_clean,
        if (is.null(cohort_clean)) "*" else cohort_clean,
        if (is.null(desc_clean)) "*" else desc_clean
      ),
      call. = FALSE
    )
  }

  # Deduplicate uncompressed/compressed stem collisions & sort deterministically
  btbl <- btbl[!duplicated(strip_neuro_ext(btbl$file)), , drop = FALSE]
  run_num <- suppressWarnings(as.numeric(canonical_run_key(btbl$run)))
  ord_b <- order(
    ifelse(is.na(btbl$ses), "", btbl$ses),
    ifelse(is.na(run_num), Inf, run_num),
    ifelse(is.na(btbl$run), "", btbl$run),
    btbl$file
  )
  btbl <- btbl[ord_b, , drop = FALSE]
  rownames(btbl) <- NULL
  run_names <- default_run_names_from_bids(btbl)
  named_bold_files <- stats::setNames(btbl$file, run_names)

  # Fail fast before reading NIfTIs if BrainGnomes provenance records physical volume deletion (scrub_timepoints)
  prov_list <- assert_no_scrubbed_bold_provenance(named_bold_files)
  bg_conf_regressed <- length(prov_list) > 0L && all(vapply(prov_list, function(p) isTRUE(p$confound_regression), logical(1L)))
  if (bg_conf_regressed && !strategy_supplied && is.null(columns) && is.null(column_pattern) && is.null(n_compcor)) {
    strategy <- "none"
  }

  # Resolve brain/template mask when mask = "auto"
  mask_files_used <- character(0L)
  mask_resolved <- mask
  if (identical(mask, "auto")) {
    mask_cands <- unique(unlist(lapply(func_dirs, function(d) {
      list.files(d, pattern = "_(mask|brainmask|templatemask)\\.nii(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
    }), use.names = FALSE))
    mtbl <- if (length(mask_cands)) parse_bids_entity_table(mask_cands) else parse_bids_entity_table(character(0L))
    if (nrow(mtbl) > 0L) {
      # Exclude internal non-brain masks such as BrainGnomes desc-intensityReferenceCore_mask.nii.gz
      valid_mask_row <- tolower(mtbl$suffix) %in% c("brainmask", "templatemask") |
        (tolower(mtbl$suffix) == "mask" & (is.na(mtbl$desc) | tolower(mtbl$desc) %in% c("brain", "brainmask")))
      mtbl <- mtbl[valid_mask_row, , drop = FALSE]
    }
    mtbl_brain <- if (nrow(mtbl) > 0L) {
      mtbl[tolower(mtbl$suffix) != "templatemask", , drop = FALSE]
    } else {
      mtbl
    }
    mtbl_tmpl <- if (nrow(mtbl) > 0L) {
      mtbl[tolower(mtbl$suffix) == "templatemask", , drop = FALSE]
    } else {
      mtbl
    }

    matched_masks <- vapply(seq_len(nrow(btbl)), function(i) {
      h <- find_matching_bids_sidecar(btbl[i, , drop = FALSE], mtbl_brain, require_space = TRUE)
      if (is.na(h) && nrow(mtbl_tmpl) > 0L) {
        h <- find_matching_bids_sidecar(btbl[i, , drop = FALSE], mtbl_tmpl, require_space = TRUE)
      }
      h
    }, character(1L))

    if (all(!is.na(matched_masks))) {
      mask_files_used <- matched_masks
      mask_resolved <- combine_fmriprep_masks(matched_masks, mode = mask_mode, threshold = mask_threshold)
    } else {
      # Fallback: check anat/*_desc-brain_mask.nii.gz or BrainGnomes provenance spatial_mask on matching grid
      anat_dir <- file.path(sub_dir, "anat")
      anat_masks <- if (dir.exists(anat_dir)) {
        list.files(anat_dir, pattern = "_desc-brain_mask\\.nii(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
      } else {
        character(0L)
      }
      prov_masks <- unique(vapply(prov_list, function(p) p$spatial_mask, character(1L)))
      prov_masks <- prov_masks[!is.na(prov_masks) & nzchar(prov_masks) & file.exists(prov_masks)]
      grid_cands <- unique(c(anat_masks, prov_masks))

      bold_geom1 <- resolve_nifti_spatial_geometry(btbl$file[1L], label = basename(btbl$file[1L]))
      grid_hit <- NA_character_
      for (am in grid_cands) {
        am_geom <- resolve_nifti_spatial_geometry(am, label = basename(am))
        same_grid <- identical(bold_geom1$dim, am_geom$dim) &&
          all(abs(bold_geom1$spacing - am_geom$spacing) <= 1e-4) &&
          all(abs(bold_geom1$affine - am_geom$affine) <= 1e-4)
        if (same_grid) {
          grid_hit <- am
          break
        }
      }
      if (!is.na(grid_hit)) {
        mask_files_used <- grid_hit
        mask_resolved <- grid_hit
      } else if (!is.null(atlas)) {
        mask_resolved <- NULL
      } else {
        stop(
          sprintf(
            "mask = 'auto' could not find matching '*_desc-brain_mask.nii(.gz)' or '*_templatemask.nii(.gz)' files for sub-%s; supply 'mask', 'atlas', or set mask = NULL explicitly.",
            sub_id
          ),
          call. = FALSE
        )
      }
    }
  }

  # Optional JSON sidecar RepetitionTime fallback if NIfTI header TR is absent
  tr_use <- tr
  if (is.null(tr_use)) {
    json_cand <- sub("\\.nii(\\.gz)?$", ".json", btbl$file[1L], ignore.case = TRUE)
    json_tr <- parse_json_numeric_field(json_cand, "RepetitionTime")
    hdr_geom1 <- resolve_nifti_spatial_geometry(btbl$file[1L], label = basename(btbl$file[1L]))
    if (!is.finite(hdr_geom1$header_tr) && is.finite(json_tr) && json_tr > 0) {
      tr_use <- json_tr
    }
  }

  ts <- read_nifti_ts(
    niftis = named_bold_files,
    mask = mask_resolved,
    atlas = atlas,
    atlas_subset = atlas_subset,
    mask_threshold = mask_threshold,
    atlas_threshold = atlas_threshold,
    tr = tr_use,
    run_names = run_names,
    id_style = id_style,
    ncpus = ncpus
  )

  # Optional external censor.1D discovery when censor = "auto"
  censor_files_used <- character(0L)
  censor_resolved <- censor
  if (identical(censor, "auto")) {
    censor_cands <- unique(unlist(lapply(func_dirs, function(d) {
      list.files(d, pattern = "_censor\\.1d(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
    }), use.names = FALSE))
    censor_tbl <- if (length(censor_cands)) parse_bids_entity_table(censor_cands) else parse_bids_entity_table(character(0L))
    censor_hits <- vapply(seq_len(nrow(btbl)), function(i) {
      find_matching_censor_sidecar(btbl[i, , drop = FALSE], censor_tbl)
    }, character(1L))
    if (all(!is.na(censor_hits))) {
      censor_files_used <- stats::setNames(censor_hits, run_names)
      censor_resolved <- censor_files_used
    } else {
      censor_resolved <- NULL
    }
  } else if (is.character(censor)) {
    censor_files_used <- censor
  }

  fd_threshold_use <- if (!is.null(censor_resolved) && !fd_threshold_supplied) NULL else fd_threshold

  # Locate fMRIPrep or BrainGnomes postprocessed confound TSVs per run
  conf_cands <- unique(unlist(lapply(func_dirs, function(d) {
    list.files(
      d,
      pattern = "_(desc-confounds_(timeseries|regressors)|confounds)\\.tsv(\\.gz)?$",
      full.names = TRUE,
      ignore.case = TRUE
    )
  }), use.names = FALSE))
  ctbl <- if (length(conf_cands)) parse_bids_entity_table(conf_cands) else parse_bids_entity_table(character(0L))
  conf_files <- vapply(seq_len(nrow(btbl)), function(i) {
    up_conf <- if (length(prov_list) >= i) prov_list[[i]]$upstream_confounds else NA_character_
    find_matching_confounds_sidecar(btbl[i, , drop = FALSE], ctbl, upstream_conf = up_conf)
  }, character(1L))

  can_build_without_conf_tsv <- identical(strategy, "none") &&
    is.null(columns) && is.null(column_pattern) &&
    (is.null(n_compcor) || n_compcor == 0L) &&
    is.null(fd_threshold_use) && is.null(dvars_threshold) &&
    !isTRUE(include_motion_outliers) &&
    is.null(motion_bandstop_bpm) && is.null(motion_lowpass_hz)

  if (anyNA(conf_files)) {
    if (can_build_without_conf_tsv) {
      conf_input <- stats::setNames(lapply(run_names, function(rn) {
        data.frame(matrix(nrow = nrow(ts$bold[[rn]]), ncol = 0L))
      }), run_names)
      conf_files <- character(0L)
    } else {
      miss_r <- run_names[is.na(conf_files)]
      stop(
        sprintf(
          "Missing confounds TSV for run(s): %s under '%s'.",
          paste(miss_r, collapse = ", "), sub_dir
        ),
        call. = FALSE
      )
    }
  } else {
    names(conf_files) <- run_names
    conf_input <- conf_files
  }

  # Optional ICA-AROMA sidecar discovery and BrainGnomes provenance coordination
  bg_aroma_applied <- length(prov_list) > 0L && any(vapply(prov_list, function(p) isTRUE(p$apply_aroma), logical(1L)))
  wants_nuis_cols <- !identical(strategy, "none") || !is.null(columns) || !is.null(column_pattern) || (!is.null(n_compcor) && n_compcor > 0L)
  is_raw_fprep_conf <- length(conf_files) > 0L &&
    any(grepl("_desc-confounds_(timeseries|regressors)\\.tsv", basename(conf_files), ignore.case = TRUE))

  aroma_obj <- NULL
  aroma_mix_used <- character(0L)
  aroma_cls_used <- character(0L)
  apply_aroma_to_bold <- FALSE
  apply_aroma_to_conf <- FALSE

  if (bg_aroma_applied) {
    if (aroma_req %in% c("nonaggressive", "aggressive", "object")) {
      warning(
        "BrainGnomes provenance indicates 'apply_aroma' was already executed on BOLD input(s); skipping duplicate AROMA regression on BOLD.",
        call. = FALSE
      )
    }
    if (is_raw_fprep_conf && wants_nuis_cols) {
      if (identical(aroma_req, "none")) {
        warning(
          "BrainGnomes provenance indicates 'apply_aroma' was executed on BOLD, while confounds resolved to un-denoised fMRIPrep 'desc-confounds_timeseries.tsv'; pass aroma = 'auto' (or use postprocessed confounds / strategy = 'none') so nuisance regression does not re-introduce AROMA motion variance.",
          call. = FALSE
        )
      } else if (isTRUE(aroma_clean_nuisance)) {
        if (identical(aroma_req, "object")) {
          aroma_obj <- align_medusa_aroma_to_runs(
            aroma, run_names,
            n_time = stats::setNames(vapply(ts$bold, nrow, integer(1L)), run_names)
          )
          apply_aroma_to_conf <- TRUE
        } else {
          disc_ar <- discover_aroma_sidecars_for_subject(btbl, sub_id, sub_dir, func_dirs, prov_list = prov_list, aroma_dir = aroma_dir)
          if (all(!is.na(disc_ar$mixing_files)) && all(!is.na(disc_ar$classification_files))) {
            aroma_mix_used <- stats::setNames(disc_ar$mixing_files, run_names)
            aroma_cls_used <- stats::setNames(disc_ar$classification_files, run_names)
            aroma_obj <- read_fmriprep_aroma(mixing = aroma_mix_used, classification = aroma_cls_used, bold = ts)
            apply_aroma_to_conf <- TRUE
          }
        }
      }
    }
  } else if (!identical(aroma_req, "none")) {
    if (identical(aroma_req, "object")) {
      aroma_obj <- align_medusa_aroma_to_runs(
        aroma, run_names,
        n_time = stats::setNames(vapply(ts$bold, nrow, integer(1L)), run_names)
      )
      apply_aroma_to_bold <- TRUE
      apply_aroma_to_conf <- isTRUE(aroma_clean_nuisance)
    } else {
      disc_ar <- discover_aroma_sidecars_for_subject(btbl, sub_id, sub_dir, func_dirs, prov_list = prov_list, aroma_dir = aroma_dir)
      if (all(!is.na(disc_ar$mixing_files)) && all(!is.na(disc_ar$classification_files))) {
        aroma_mix_used <- stats::setNames(disc_ar$mixing_files, run_names)
        aroma_cls_used <- stats::setNames(disc_ar$classification_files, run_names)
        aroma_obj <- read_fmriprep_aroma(mixing = aroma_mix_used, classification = aroma_cls_used, bold = ts)
        apply_aroma_to_bold <- TRUE
        apply_aroma_to_conf <- isTRUE(aroma_clean_nuisance)
      } else if (!identical(aroma_req, "auto")) {
        miss_mix <- run_names[is.na(disc_ar$mixing_files)]
        miss_cls <- run_names[is.na(disc_ar$classification_files)]
        has_legacy_aroma_cols <- length(conf_files) > 0L && any(vapply(conf_files, function(cf) {
          hdr <- tryCatch(names(data.table::fread(cf, nrows = 0L, data.table = FALSE)), error = function(e) character(0L))
          any(grepl("^aroma_motion_", hdr))
        }, logical(1L)))
        extra_hint <- if (has_legacy_aroma_cols) {
          " Note: 'aroma_motion_*' columns in confounds TSV contain only rejected components; non-aggressive AROMA requires the full MELODIC mixing matrix ('*_desc-MELODIC_mixing.tsv' + '*_AROMAnoiseICs.csv' or fMRIPost-AROMA '*_desc-melodic_mixing.tsv' + '*_desc-aroma_metrics.tsv')."
        } else {
          ""
        }
        stop(
          sprintf(
            "Could not locate AROMA sidecars for sub-%s (%s)%s; supply 'aroma_dir' pointing to fMRIPost-AROMA derivatives or set aroma = 'none'.%s",
            sub_id,
            paste(c(
              if (length(miss_mix)) sprintf("missing '*_desc-melodic_mixing.tsv' for %s", paste(miss_mix, collapse = ", ")) else NULL,
              if (length(miss_cls)) sprintf("missing '*_desc-aroma_metrics.tsv' / '*_AROMAnoiseICs.csv' for %s", paste(miss_cls, collapse = ", ")) else NULL
            ), collapse = "; "),
            if (!is.null(aroma_dir)) sprintf(" in aroma_dir ('%s')", aroma_dir) else "",
            extra_hint
          ),
          call. = FALSE
        )
      }
    }
  }

  conf_obj <- read_fmriprep_confounds(
    confounds = conf_input,
    bold = ts,
    strategy = strategy,
    columns = columns,
    column_pattern = column_pattern,
    n_compcor = n_compcor,
    compcor_type = compcor_type,
    include_cosines = include_cosines,
    censor = censor_resolved,
    fd_threshold = fd_threshold_use,
    fd_col = fd_col,
    head_radius = head_radius,
    rot_units = rot_units,
    tr = tr,
    motion_bandstop_bpm = motion_bandstop_bpm,
    motion_lowpass_hz = motion_lowpass_hz,
    motion_filter_order = motion_filter_order,
    motion_filter_pad = motion_filter_pad,
    motion_filter_target = motion_filter_target,
    dvars_threshold = dvars_threshold,
    dvars_col = dvars_col,
    include_non_steady_state = include_non_steady_state,
    include_motion_outliers = include_motion_outliers,
    drop_initial_volumes = drop_initial_volumes,
    censor_prev = censor_prev,
    censor_next = censor_next,
    min_segment_length = min_segment_length,
    aroma = if (apply_aroma_to_conf) aroma_obj else NULL,
    aroma_mode = aroma_mode_eff,
    aroma_fit_valid = aroma_fit_valid
  )

  if (apply_aroma_to_bold && !is.null(aroma_obj)) {
    ts <- apply_aroma(
      x = ts,
      aroma = aroma_obj,
      mode = aroma_mode_eff,
      clean_nuisance = FALSE,
      fit_valid_rows = aroma_fit_valid,
      valid_rows = conf_obj$valid_rows
    )
    ts$aroma$clean_nuisance <- aroma_clean_nuisance
    ts$aroma$nuisance_cleaned <- isTRUE(conf_obj$settings$aroma$applied)
  }

  # Optional BIDS events discovery / ingestion
  ev_out <- NULL
  ev_files_used <- character(0L)
  if (identical(events, "auto")) {
    sub_dir_bids <- if (!is.null(bids_dir)) {
      cand_b <- file.path(bids_dir, paste0("sub-", sub_id))
      if (dir.exists(cand_b)) cand_b else NULL
    } else {
      NULL
    }
    ev_hits <- vapply(seq_len(nrow(btbl)), function(i) {
      h <- NA_character_
      if (!is.null(bids_dir)) {
        h <- find_inherited_bids_events(btbl[i, , drop = FALSE], bids_dir = bids_dir, sub_dir_bids = sub_dir_bids)
      }
      if (is.na(h)) {
        h <- find_inherited_bids_events(btbl[i, , drop = FALSE], bids_dir = fmriprep_dir, sub_dir_bids = sub_dir)
      }
      h
    }, character(1L))
    if (all(!is.na(ev_hits))) {
      ev_files_used <- stats::setNames(ev_hits, run_names)
      ev_out <- read_bids_events(
        events = ev_files_used,
        bold = ts,
        condition_col = condition_col,
        condition_map = condition_map,
        default_condition = default_condition,
        exclude_conditions = exclude_conditions,
        na_duration = na_duration
      )
    } else if (!is.null(bids_dir)) {
      miss_ev <- run_names[is.na(ev_hits)]
      stop(
        sprintf(
          "Could not locate BIDS '*_events.tsv' in bids_dir ('%s') for run(s): %s; pass events = NULL for task-free data.",
          bids_dir, paste(miss_ev, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  } else if (!is.null(events)) {
    ev_out <- read_bids_events(
      events = events,
      bold = ts,
      condition_col = condition_col,
      condition_map = condition_map,
      default_condition = default_condition,
      exclude_conditions = exclude_conditions,
      na_duration = na_duration
    )
    if (is.character(events)) ev_files_used <- events
  }

  ts$nuisance <- conf_obj$nuisance
  ts$valid_rows <- conf_obj$valid_rows
  ts$confounds <- conf_obj
  ts$events <- ev_out
  ts$bids <- list(
    sub = sub_id,
    ses = unique(btbl$ses[!is.na(btbl$ses)]),
    task = unique(btbl$task[!is.na(btbl$task)]),
    space = unique(btbl$space[!is.na(btbl$space)]),
    cohort = unique(btbl$cohort[!is.na(btbl$cohort)]),
    desc = unique(btbl$desc[!is.na(btbl$desc)]),
    bold_files = named_bold_files,
    confound_files = conf_files,
    censor_files = censor_files_used,
    mask_files = mask_files_used,
    event_files = ev_files_used,
    aroma_mixing_files = aroma_mix_used,
    aroma_classification_files = aroma_cls_used,
    postproc_provenance = prov_list
  )
  ts
}
