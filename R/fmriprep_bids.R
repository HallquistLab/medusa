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
      bold_files = NULL
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
    if (!is.null(run_names) && !identical(run_names, b_runs)) {
      stop("run_names does not match names(bold$bold).", call. = FALSE)
    }
    return(list(run_names = b_runs, n_time = n_t, bold_files = b_files))
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
    return(list(run_names = r_nm, n_time = stats::setNames(as.integer(n_t), r_nm), bold_files = b_file))
  }
  if (is.character(bold) && is.null(run_names)) {
    checkmate::assert_character(bold, min.len = 1L, any.missing = FALSE, min.chars = 1L, unique = TRUE)
    return(list(run_names = bold, n_time = NULL, bold_files = NULL))
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
    return(list(run_names = b_runs, n_time = n_t, bold_files = NULL))
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
    temporal_filter = FALSE,
    scrub_interpolate = FALSE,
    spatial_mask = NA_character_,
    upstream_confounds = NA_character_
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
  }

  list(
    is_braingnomes = TRUE,
    operations = active_ops,
    scrub_timepoints = "scrub_timepoints" %in% active_ops,
    confound_regression = "confound_regression" %in% active_ops,
    temporal_filter = "temporal_filter" %in% active_ops,
    scrub_interpolate = "scrub_interpolate" %in% active_ops,
    spatial_mask = spatial_mask,
    upstream_confounds = upstream_conf
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
    if (!is.null(sc) && tolower(sp) %in% c("6p", "12p", "24p", "27p", "24p_csf_wm", "32p", "36p")) {
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
#'   \code{"36p"}), or zero-padded index range specifications such as
#'   \code{"a_comp_cor_<0-4>"} or \code{"cosine<0,2-3>"} (in addition to
#'   \code{strategy}, or exact set when \code{strategy = "custom"}).
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
#'   or \code{"deg"}) used only when computing fallback Power FD.
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
#'       and retained nuisance column counts.
#'     \item \code{framewise_displacement}: Run-named list of numeric FD vectors.
#'   }
#'   The \code{"medusa_confounds"} object can be passed directly as
#'   \code{nuisance = conf} to \code{\link{estimate_crossrun_hrf}} (automatically
#'   supplying both \code{nuisance} and \code{valid_rows}) or as
#'   \code{valid = conf} to \code{\link{align_medusa_ts}}.
#'
#' @seealso \code{\link{read_bids_events}}, \code{\link{read_fmriprep_subject}},
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
  strategy = c("24p", "6p", "12p", "24p_csf_wm", "27p", "36p", "compcor", "custom", "none"),
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
  na_strings = c("n/a", "N/A", "NA", ".", "NaN", "")
) {
  fd_threshold_supplied <- !missing(fd_threshold)
  strategy <- match.arg(strategy)
  compcor_type <- match.arg(compcor_type)
  rot_units <- match.arg(rot_units)
  impute <- match.arg(impute)
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
  strategy_cols <- confound_strategy_columns(strategy)

  n_runs <- length(resolved_runs)
  nuis_list <- valid_list <- fd_list <- stats::setNames(vector("list", n_runs), resolved_runs)
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
      } else if (!is.null(fd_threshold)) {
        stop(
          sprintf(
            "Neither '%s' nor the 6 rigid-body motion columns found in %s; set fd_threshold = NULL to disable FD censoring or specify fd_col.",
            fd_col_resolved, run_label
          ),
          call. = FALSE
        )
      }
    } else if (!is.null(fd_threshold)) {
      stop(
        sprintf(
          "Column '%s' not found in %s; set fd_threshold = NULL to disable FD censoring or specify fd_col.",
          fd_col, run_label
        ),
        call. = FALSE
      )
    }

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
    }

    spike_censor <- dilate_censor_mask(censor_fd | censor_dvars | censor_mo, censor_prev, censor_next)
    valid_r <- valid_ext & !(censor_init | censor_nss | spike_censor)
    valid_r <- enforce_min_segment_length(valid_r, min_segment_length)
    valid_list[[r]] <- valid_r
    fd_list[[r]] <- fd_vec

    # 2. Nuisance regressor selection & finite sanitization
    if (identical(strategy, "none") && is.null(columns) && is.null(column_pattern) && n_compcor == 0L) {
      z_mat <- matrix(numeric(0L), nrow = n_vol, ncol = 0L)
    } else {
      cc_cols <- select_compcor_columns(names(df), n_compcor, compcor_type)
      if (identical(strategy, "compcor") && n_compcor > 0L && length(cc_cols) == 0L) {
        stop(sprintf("No CompCor columns matching compcor_type = '%s' found in %s.", compcor_type, run_label), call. = FALSE)
      }
      pat_cols <- if (!is.null(column_pattern)) grep(column_pattern, names(df), value = TRUE) else character(0L)
      exp_columns <- expand_confound_column_specs(columns, names(df))
      if (length(exp_columns)) {
        resolvable <- vapply(exp_columns, can_resolve_confound_column, logical(1L), avail_names = names(df))
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
        cos_cols <- grep("^cosine[0-9]+$", names(df), value = TRUE)
        if (length(cos_cols)) {
          cos_idx <- as.integer(sub("^cosine", "", cos_cols))
          cos_cols <- cos_cols[order(cos_idx)]
          req_cols <- unique(c(req_cols, cos_cols))
        }
      }

      if (length(req_cols) == 0L) {
        z_mat <- matrix(numeric(0L), nrow = n_vol, ncol = 0L)
      } else {
        col_vecs <- lapply(req_cols, function(cn) expand_fmriprep_column(df, cn, run_label))
        z_mat <- do.call(cbind, col_vecs)
        colnames(z_mat) <- req_cols
        storage.mode(z_mat) <- "double"

        is_deriv_like <- grepl("(_derivative1|^framewise_displacement(_filtered)?$|^dvars$|^std_dvars$|^rmsd$)", req_cols)
        if (n_vol >= 1L && any(is_deriv_like)) {
          d_idx <- which(is_deriv_like)
          bad_row1 <- d_idx[!is.finite(z_mat[1L, d_idx])]
          if (length(bad_row1)) z_mat[1L, bad_row1] <- 0
        }

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
        z_mat <- z_mat[, keep_col, drop = FALSE]
      }
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
      n_nuisance_cols = as.integer(ncol(z_mat)),
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, sum_rows)
  rownames(summary_df) <- NULL

  structure(
    list(
      nuisance = nuis_list,
      valid_rows = valid_list,
      summary = summary_df,
      framewise_displacement = fd_list,
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
        min_segment_length = min_segment_length
      )
    ),
    class = "medusa_confounds"
  )
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
  cat(sprintf(
    "<medusa_confounds> %d run%s | strategy: %s | fd_threshold: %s\n",
    n_runs,
    if (n_runs == 1L) "" else "s",
    x$settings$strategy,
    fd_lbl
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

#' Discover and ingest an fMRIPrep or postprocessed participant into a unified MEDuSA time-series bundle
#'
#' Locates a participant's preprocessed or \code{BrainGnomes}-postprocessed BOLD
#' runs (\code{*_desc-preproc_bold.nii.gz} or \code{*_desc-<desc>_bold.nii.gz}),
#' brain or template masks (\code{*_desc-brain_mask.nii.gz},
#' \code{*_templatemask.nii.gz}), confound tables
#' (\code{*_desc-confounds_timeseries.tsv}, \code{*_desc-<desc>_confounds.tsv}),
#' optional \code{*_censor.1D} masks, and optional BIDS event tables
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
#' @param strategy,columns,column_pattern,n_compcor,compcor_type,include_cosines
#'   Nuisance selection options forwarded to \code{\link{read_fmriprep_confounds}}.
#'   When \code{BrainGnomes} provenance indicates \code{confound_regression} was
#'   already applied to the BOLD derivative and \code{strategy} is not explicitly
#'   specified, \code{strategy} defaults automatically to \code{"none"} to avoid
#'   double-regressing unfiltered confounds.
#' @param censor Either \code{"auto"} (default; automatically uses matching
#'   \code{*_censor.1D} sidecars when present for postprocessed outputs), explicit
#'   censor file paths / vectors forwarded to \code{\link{read_fmriprep_confounds}},
#'   or \code{NULL} to ignore \code{*_censor.1D} sidecars.
#' @param fd_threshold,fd_col,dvars_threshold,dvars_col,include_non_steady_state,include_motion_outliers,drop_initial_volumes,censor_prev,censor_next,min_segment_length
#'   Clock-preserving censoring options forwarded to
#'   \code{\link{read_fmriprep_confounds}}.
#' @param id_style Voxel ID convention forwarded to \code{\link{read_nifti_ts}}
#'   (\code{"ijk"} or \code{"vnum"}).
#' @param ncpus Number of parallel workers forwarded to \code{\link{read_nifti_ts}}.
#'
#' @return An S3 object of class \code{"medusa_nifti_ts"} enriched with
#'   \code{$nuisance}, \code{$valid_rows}, \code{$confounds}
#'   (\code{"medusa_confounds"}), \code{$events} (when discovered/supplied), and
#'   \code{$bids} provenance metadata.
#'
#' @seealso \code{\link{read_nifti_ts}}, \code{\link{read_fmriprep_confounds}},
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
  strategy = c("24p", "6p", "12p", "24p_csf_wm", "27p", "36p", "compcor", "custom", "none"),
  columns = NULL,
  column_pattern = NULL,
  n_compcor = NULL,
  compcor_type = c("anatomical", "split", "temporal"),
  include_cosines = "auto",
  censor = "auto",
  fd_threshold = 0.5,
  fd_col = "auto",
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
  id_style <- match.arg(id_style)
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
    !isTRUE(include_motion_outliers)

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
    dvars_threshold = dvars_threshold,
    dvars_col = dvars_col,
    include_non_steady_state = include_non_steady_state,
    include_motion_outliers = include_motion_outliers,
    drop_initial_volumes = drop_initial_volumes,
    censor_prev = censor_prev,
    censor_next = censor_next,
    min_segment_length = min_segment_length
  )

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
    postproc_provenance = prov_list
  )
  ts
}
