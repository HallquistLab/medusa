# Fast C++-backed matrix-native event alignment and epoch interpolation

normalize_matrix_orientation <- function(mat, orientation = NULL, label = "ts") {
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop(label, " must be a numeric matrix.", call. = FALSE)
  }
  orient <- if (!is.null(orientation)) {
    match.arg(orientation, c("time_by_voxel", "voxel_by_time"))
  } else if (!is.null(attr(mat, "orientation", exact = TRUE))) {
    match.arg(attr(mat, "orientation", exact = TRUE), c("time_by_voxel", "voxel_by_time"))
  } else {
    "time_by_voxel"
  }
  if (identical(orient, "voxel_by_time")) {
    out <- t(mat)
    for (at in c("tr", "time_offset", "observation_mask", "voxel_map", "parcel_ids", "time")) {
      val <- attr(mat, at, exact = TRUE)
      if (!is.null(val)) attr(out, at) <- val
    }
    return(out)
  }
  mat
}

build_unit_metadata_df <- function(n_units, col_names = NULL, unit_data = NULL, parcel_ids = NULL) {
  if (!is.null(unit_data)) {
    checkmate::assert_data_frame(unit_data)
    if (nrow(unit_data) != n_units) {
      stop("unit_data must have one row per spatial column (", n_units, " columns, got ",
           nrow(unit_data), " rows).", call. = FALSE)
    }
    udf <- as.data.frame(unit_data, stringsAsFactors = FALSE)
  } else {
    udf <- data.frame(
      voxel = if (!is.null(col_names) && length(col_names) == n_units) {
        as.character(col_names)
      } else {
        sprintf("v%d", seq_len(n_units))
      },
      vnum = seq_len(n_units),
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(parcel_ids)) {
    if (length(parcel_ids) != n_units) {
      stop("parcel_ids must have length equal to the number of spatial columns (",
           n_units, ", got ", length(parcel_ids), ").", call. = FALSE)
    }
    if (!is.null(names(parcel_ids)) && !is.null(col_names) && length(col_names) == n_units) {
      if (!setequal(names(parcel_ids), col_names)) {
        stop("Named parcel_ids do not match matrix column names.", call. = FALSE)
      }
      parcel_vals <- unname(parcel_ids[col_names])
    } else {
      parcel_vals <- unname(parcel_ids)
    }
    udf$parcel <- parcel_vals
    if (!("atlas_value" %in% names(udf))) {
      udf$atlas_value <- parcel_vals
    }
  } else if ("atlas_value" %in% names(udf) && !("parcel" %in% names(udf))) {
    udf$parcel <- udf$atlas_value
  } else if ("parcel" %in% names(udf) && !("atlas_value" %in% names(udf))) {
    udf$atlas_value <- udf$parcel
  }

  udf
}

align_matrix_single_run <- function(
  signals,
  raw_time,
  n_time_raw,
  n_units,
  unit_df,
  events,
  event,
  trial_col = "trial",
  time_before = -3.0,
  time_after = 3.0,
  collide_before = NULL,
  collide_after = NULL,
  pad_before = -1.5,
  pad_after = 1.5,
  output_resolution = NULL,
  tr = NULL,
  group_by = "trial",
  time_audit = FALSE,
  overlap_policy = NULL,
  raw_valid = NULL,
  neighbor_events = NULL,
  max_gap = NULL,
  readout_kernel = "linear",
  readout_width = NULL,
  coerce_factor_trial = FALSE
) {
  overlap_policy <- resolve_overlap_policy(overlap_policy, collide_before, collide_after)
  if (is.null(neighbor_events)) neighbor_events <- event
  checkmate::assert_character(neighbor_events, min.len = 1L, unique = TRUE)
  checkmate::assert_number(time_before, upper = 0, finite = TRUE)
  checkmate::assert_number(time_after, lower = 0, finite = TRUE)
  checkmate::assert_number(pad_before, upper = 0, finite = TRUE)
  checkmate::assert_number(pad_after, lower = 0, finite = TRUE)
  checkmate::assert_data_frame(events)

  if (n_time_raw < 1L || n_units < 1L) {
    stop("No time series samples supplied for event alignment")
  }
  checkmate::assert_string(event)
  checkmate::assert_string(trial_col)
  checkmate::assert_subset(c(event, collide_before, collide_after, neighbor_events), names(events))
  if (nrow(events) == 0L) stop("No trial events supplied for alignment")
  if (!(trial_col %in% names(events))) {
    stop("Event data must contain trial identifier column '", trial_col, "'")
  }
  if (anyNA(events[[trial_col]]) || anyDuplicated(events[[trial_col]])) {
    stop("Event data must contain unique, nonmissing trial identifiers")
  }
  for (column in unique(c(event, collide_before, collide_after, neighbor_events))) {
    checkmate::assert_numeric(events[[column]], all.missing = TRUE)
  }

  if (is.null(output_resolution)) {
    if (is.null(tr) || !is.finite(tr) || tr <= 0) {
      stop("Specify 'output_resolution' or a positive 'tr'.")
    }
    output_resolution <- tr
  }
  tout <- epoch_output_grid(time_before, time_after, output_resolution)
  if (is.null(max_gap)) {
    if (is.null(tr) || !is.finite(tr) || tr <= 0) {
      stop("Specify 'max_gap' or a positive 'tr'.")
    }
    max_gap <- tr
  }
  checkmate::assert_number(max_gap, lower = .Machine$double.eps, finite = TRUE)
  quad <- epoch_kernel_quadrature(readout_kernel, readout_width, max_gap)

  sig_names <- names(signals)
  if (is.null(sig_names) || any(!nzchar(sig_names)) || anyDuplicated(sig_names)) {
    stop("Signal matrices must have unique, nonempty names.")
  }
  for (s in seq_along(signals)) {
    checkmate::assert_numeric(signals[[s]], all.missing = TRUE)
    if (!is.double(signals[[s]])) {
      signals[[s]] <- as.double(signals[[s]])
    }
  }

  avail_group_cols <- unique(c(trial_col, names(unit_df)))
  checkmate::assert_subset(group_by, avail_group_cols)

  # Normalize clock, signals, and optional censor mask
  checkmate::assert_numeric(raw_time, len = n_time_raw, all.missing = TRUE)
  fin_mask <- is.finite(raw_time)

  cens_is_matrix <- FALSE
  if (!is.null(raw_valid)) {
    if (!is.logical(raw_valid)) stop("valid_col must name a logical column")
    vlen <- length(raw_valid)
    if (vlen == n_time_raw) {
      cens_is_matrix <- FALSE
    } else if (vlen == as.numeric(n_time_raw) * as.numeric(n_units)) {
      cens_is_matrix <- TRUE
    } else {
      stop("valid mask must have length equal to nrow(ts) or nrow(ts) * ncol(ts).")
    }
  }

  if (all(fin_mask) && !is.unsorted(raw_time, strictly = TRUE)) {
    clock <- as.double(raw_time)
    censored <- if (is.null(raw_valid)) logical(0L) else (is.na(raw_valid) | !raw_valid)
  } else if (!any(fin_mask)) {
    clock <- numeric(0L)
    signals <- lapply(signals, function(x) numeric(0L))
    censored <- logical(0L)
    cens_is_matrix <- FALSE
  } else {
    fin_idx <- which(fin_mask)
    fin_time <- as.double(raw_time[fin_idx])
    ord <- order(fin_time)
    clock <- fin_time[ord]
    if (anyDuplicated(clock)) {
      stop("Duplicate source times within an epoch unit")
    }
    keep_rows <- fin_idx[ord]
    signals <- lapply(signals, function(x) {
      dim(x) <- c(n_time_raw, n_units)
      as.double(x[keep_rows, , drop = FALSE])
    })
    if (is.null(raw_valid)) {
      censored <- logical(0L)
    } else if (!cens_is_matrix) {
      rv <- raw_valid[keep_rows]
      censored <- is.na(rv) | !rv
    } else {
      dim(raw_valid) <- c(n_time_raw, n_units)
      rv <- raw_valid[keep_rows, , drop = FALSE]
      censored <- as.logical(is.na(rv) | !rv)
    }
  }

  run_df <- data.table::as.data.table(data.table::copy(events))
  if (isTRUE(coerce_factor_trial) && is.factor(run_df[[trial_col]])) {
    run_df[[trial_col]] <- as.character(run_df[[trial_col]])
  }
  data.table::setorderv(run_df, trial_col)
  n_trials <- nrow(run_df)

  annotations <- epoch_event_annotations(
    run_df, event, trial_col, neighbor_events,
    time_before, time_after,
    if (length(clock)) clock[1L] else NA_real_,
    if (length(clock)) clock[length(clock)] else NA_real_,
    overlap_policy
  )

  onsets <- as.double(run_df[[event]])
  evt_before <- rep(-Inf, n_trials)
  evt_after <- rep(Inf, n_trials)
  if (identical(overlap_policy, "truncate")) {
    if (length(collide_before)) {
      bvals <- unlist(run_df[, collide_before, with = FALSE], use.names = FALSE)
      bvals <- sort(unique(as.double(bvals[is.finite(bvals)])))
      if (length(bvals)) {
        fin_o <- which(is.finite(onsets))
        idx_b <- findInterval(onsets[fin_o], bvals, left.open = TRUE)
        has_b <- idx_b > 0L
        evt_before[fin_o[has_b]] <- bvals[idx_b[has_b]]
      }
    }
    if (length(collide_after)) {
      avals <- unlist(run_df[, collide_after, with = FALSE], use.names = FALSE)
      avals <- sort(unique(as.double(avals[is.finite(avals)])))
      if (length(avals)) {
        fin_o <- which(is.finite(onsets))
        idx_a <- findInterval(onsets[fin_o], avals)
        has_a <- idx_a < length(avals)
        evt_after[fin_o[has_a]] <- avals[idx_a[has_a] + 1L]
      }
    }
  }

  trial_grouped <- trial_col %in% group_by
  spatial_cols <- setdiff(group_by, trial_col)
  if (length(spatial_cols) == 0L) {
    n_groups <- 1L
    unit_group <- rep(0L, n_units)
    group_dt <- NULL
  } else {
    u_dt <- data.table::as.data.table(unit_df[, spatial_cols, drop = FALSE])
    group_dt <- unique(u_dt)
    data.table::setorderv(group_dt, spatial_cols)
    group_dt[, .gid := seq_len(.N) - 1L]
    u_dt[, .u_ord := seq_len(.N)]
    u_matched <- merge(u_dt, group_dt, by = spatial_cols, sort = FALSE)
    data.table::setorderv(u_matched, ".u_ord")
    unit_group <- as.integer(u_matched$.gid)
    group_dt[, .gid := NULL]
    n_groups <- nrow(group_dt)
  }

  spatial_first <- trial_grouped && length(spatial_cols) > 0L && !identical(group_by[1L], trial_col)

  cpp_res <- align_interpolate_matrix_cpp(
    signals = signals,
    clock = clock,
    censored = censored,
    cens_is_matrix = cens_is_matrix,
    unit_group = unit_group,
    n_groups = n_groups,
    onsets = onsets,
    evt_before = evt_before,
    evt_after = evt_after,
    overlap_n = annotations$align_window_overlap_n,
    overlap_policy = overlap_policy,
    time_before = time_before,
    time_after = time_after,
    pad_before = pad_before,
    pad_after = pad_after,
    tout = tout,
    max_gap = max_gap,
    readout_kernel = quad$kernel,
    quad_offsets = quad$offsets,
    quad_weights = quad$weights,
    trial_grouped = trial_grouped,
    spatial_first = spatial_first
  )

  n_lags <- length(tout)
  out_cols <- list()

  if (trial_grouped) {
    t_rep <- if (spatial_first) {
      rep(rep(seq_len(n_trials), each = n_lags), times = n_groups)
    } else {
      rep(seq_len(n_trials), each = n_groups * n_lags)
    }
    g_rep <- if (length(spatial_cols) > 0L) {
      if (spatial_first) {
        rep(seq_len(n_groups), each = n_trials * n_lags)
      } else {
        rep(rep(seq_len(n_groups), each = n_lags), times = n_trials)
      }
    } else {
      NULL
    }
    for (col in group_by) {
      if (identical(col, trial_col)) {
        out_cols[[col]] <- run_df[[trial_col]][t_rep]
      } else {
        out_cols[[col]] <- group_dt[[col]][g_rep]
      }
    }
  } else {
    for (col in spatial_cols) {
      out_cols[[col]] <- rep(group_dt[[col]], each = n_lags)
    }
  }

  out_cols[["evt_time"]] <- rep(tout, times = if (trial_grouped) n_trials * n_groups else n_groups)

  stat_names <- c(
    "mean", "median", "sd", "n", "n_total", "n_censored",
    "n_missing", "n_gap", "n_outside", "source_ids", "source_max_epochs"
  )
  for (st in stat_names) {
    for (s_idx in seq_along(sig_names)) {
      out_cols[[paste0(sig_names[s_idx], "_", st)]] <- cpp_res$signals[[s_idx]][[st]]
    }
  }

  if (trial_grouped) {
    out_cols[["epoch_status"]] <- cpp_res$epoch_status
    for (acol in names(annotations)) {
      out_cols[[acol]] <- annotations[[acol]][t_rep]
    }
    if (isTRUE(time_audit)) {
      eo <- onsets[t_rep]
      out_cols[["evt_onset"]] <- eo
      out_cols[["collide_before"]] <- evt_before[t_rep]
      out_cols[["collide_after"]] <- evt_after[t_rep]
      out_cols[["time"]] <- ifelse(is.finite(eo), out_cols[["evt_time"]] + eo, NA_real_)
    }
  }

  data.table::setDT(out_cols)
  if (trial_grouped) {
    if (length(group_by) >= 3L && !(match(trial_col, group_by) %in% c(1L, length(group_by)))) {
      data.table::setorderv(out_cols, c(group_by, "evt_time"))
    }
  } else {
    data.table::setkeyv(out_cols, c(group_by, "evt_time"))
  }

  data.table::setattr(
    out_cols,
    "source_clock",
    data.frame(source_id = seq_along(clock), time = clock)
  )
  out_cols
}

extract_fmri_ts_matrix_view <- function(fmri_obj, valid_col = NULL, group_by = NULL) {
  ts_dt <- fmri_obj$ts_data
  vm <- fmri_obj$vm
  vmvec <- fmri_obj$get_vmvec()
  if (is.null(ts_dt) || nrow(ts_dt) == 0L) {
    stop("No time series samples supplied for event alignment")
  }

  ts_orig_cols <- ifelse(names(ts_dt) %in% names(vmvec), unname(vmvec[names(ts_dt)]), names(ts_dt))
  all_orig_cols <- unique(c(vm$key, ts_orig_cols))
  if (any(grepl("^(align_|\\.align_)", all_orig_cols))) {
    stop("Input contains reserved alignment columns")
  }

  checkmate::assert_string(valid_col, null.ok = TRUE)
  raw_valid <- NULL
  if (!is.null(valid_col)) {
    checkmate::assert_subset(valid_col, all_orig_cols)
    if (valid_col %in% vm$key) {
      return(NULL)
    }
    int_valid <- if (valid_col %in% vmvec) names(vmvec)[match(valid_col, vmvec)] else valid_col
    raw_valid <- ts_dt[[int_valid]]
    if (!is.logical(raw_valid)) stop("valid_col must name a logical column")
  }

  spatial_req <- setdiff(group_by, vm$trial)
  if (length(spatial_req) > 0L && !all(spatial_req %in% vm$key)) {
    return(NULL)
  }

  ts_keys <- fmri_obj$ts_keys
  if (is.null(ts_keys) || length(ts_keys) == 0L) {
    n_units <- 1L
    n_time_raw <- nrow(ts_dt)
    unit_df <- data.frame(.unit_id = 1L)
  } else {
    ends <- sort(unique(unlist(lapply(ts_keys, function(r) cumsum(r$lengths)), use.names = FALSE)))
    starts <- c(1L, head(ends, -1L) + 1L)
    u_lens <- ends - starts + 1L
    n_time_raw <- u_lens[1L]
    if (!all(u_lens == n_time_raw)) {
      return(NULL)
    }
    n_units <- length(u_lens)
    if (!check_shared_clock_cpp(ts_dt$time, n_time_raw, n_units)) {
      return(NULL)
    }
    unit_cols <- vector("list", length(ts_keys))
    names(unit_cols) <- vm$key
    for (kk in seq_along(ts_keys)) {
      rk <- ts_keys[[kk]]
      r_ends <- cumsum(rk$lengths)
      run_idx <- findInterval(starts - 1L, r_ends) + 1L
      unit_cols[[kk]] <- rk$values[run_idx]
    }
    unit_df <- as.data.frame(unit_cols, stringsAsFactors = FALSE, check.names = FALSE)
  }

  vvars <- if (length(vm$value) > 1L) paste0("value", seq_along(vm$value)) else "value"
  signals <- stats::setNames(vector("list", length(vvars)), vm$value)
  for (vv in seq_along(vvars)) {
    signals[[vv]] <- ts_dt[[vvars[vv]]]
  }

  raw_time <- ts_dt$time[seq_len(n_time_raw)]
  list(
    signals = signals,
    raw_time = raw_time,
    n_time_raw = n_time_raw,
    n_units = n_units,
    unit_df = unit_df,
    raw_valid = raw_valid,
    events = fmri_obj$get_events(orig_names = TRUE),
    trial_col = vm$trial,
    tr = fmri_obj$tr
  )
}

canonical_run_key <- function(x) {
  s <- tolower(trimws(as.character(x)))
  is_simple <- !is.na(s) & grepl("^(run[-_]?)?[0-9]+[[:alpha:]]*$", s, perl = TRUE)
  out <- rep(NA_character_, length(s))
  if (any(is_simple)) {
    stripped <- sub("^run[-_]?", "", s[is_simple], perl = TRUE)
    stripped <- sub("^0+(?=[0-9])", "", stripped, perl = TRUE)
    out[is_simple] <- stripped
  }
  out
}

split_multirun_events <- function(events, run_names, run_col = "run", allow_empty_runs = FALSE) {
  n_runs <- length(run_names)
  if (is.data.frame(events)) {
    if (n_runs == 1L && !(run_col %in% names(events))) {
      return(stats::setNames(list(events), run_names))
    }
    if (!(run_col %in% names(events))) {
      stop("Multi-run event alignment requires either a list of per-run event tables or a '",
           run_col, "' column in 'events'.", call. = FALSE)
    }
    ev_col <- events[[run_col]]
    ev_runs <- as.character(ev_col)
    canon_runs <- canonical_run_key(run_names)
    canon_ev <- canonical_run_key(ev_runs)
    ev_levels <- if (is.factor(ev_col)) levels(ev_col) else character(0L)
    canon_levels <- canonical_run_key(ev_levels)
    out <- stats::setNames(vector("list", n_runs), run_names)
    for (r in seq_len(n_runs)) {
      rname <- run_names[r]
      hit <- !is.na(ev_runs) & (ev_runs == rname)
      can_use_canon <- !is.na(canon_runs[r]) && sum(canon_runs == canon_runs[r], na.rm = TRUE) == 1L
      if (!any(hit) && can_use_canon) {
        hit <- !is.na(canon_ev) & (canon_ev == canon_runs[r])
      }
      if (!any(hit)) {
        in_factor_levels <- (rname %in% ev_levels) ||
          (can_use_canon && any(!is.na(canon_levels) & canon_levels == canon_runs[r]))
        if (isTRUE(allow_empty_runs) && in_factor_levels) {
          out[[r]] <- events[0L, , drop = FALSE]
          next
        }
        stop("No rows in 'events' matched run '", rname, "' in column '", run_col, "'.", call. = FALSE)
      }
      out[[r]] <- events[which(hit), , drop = FALSE]
    }
    return(out)
  }

  if (is.list(events)) {
    if (length(events) != n_runs) {
      stop("Length of 'events' list (", length(events), ") must match number of runs (",
           n_runs, ").", call. = FALSE)
    }
    if (!is.null(names(events)) && all(nzchar(names(events)))) {
      if (!setequal(names(events), run_names)) {
        stop("Named 'events' list must match run names: ",
             paste(run_names, collapse = ", "), ".", call. = FALSE)
      }
      events <- events[run_names]
    }
    for (r in seq_len(n_runs)) {
      checkmate::assert_data_frame(events[[r]])
    }
    return(stats::setNames(events, run_names))
  }

  stop("'events' must be a data.frame or a list of per-run data.frames.", call. = FALSE)
}

resolve_run_valid_mask <- function(valid, r_idx, r_name, mat) {
  if (is.null(valid)) return(NULL)
  if (is.character(valid) && length(valid) == 1L && identical(valid, "observation_mask")) {
    om <- attr(mat, "observation_mask", exact = TRUE)
    if (is.null(om)) {
      stop("valid = 'observation_mask' requested, but matrix has no 'observation_mask' attribute.", call. = FALSE)
    }
    return(om)
  }
  if (is.list(valid)) {
    if (!is.null(names(valid)) && r_name %in% names(valid)) {
      return(valid[[r_name]])
    }
    return(valid[[r_idx]])
  }
  valid
}

resolve_run_time_vec <- function(time, r_idx, r_name, n_t, tr, time_offset, mat) {
  if (!is.null(time)) {
    if (is.list(time)) {
      tvec <- if (!is.null(names(time)) && r_name %in% names(time)) time[[r_name]] else time[[r_idx]]
      return(as.double(tvec))
    }
    return(as.double(time))
  }
  mat_time <- attr(mat, "time", exact = TRUE)
  if (!is.null(mat_time)) {
    return(as.double(mat_time))
  }
  tr_r <- if (!is.null(tr)) tr else attr(mat, "tr", exact = TRUE)
  if (is.null(tr_r) || !is.finite(tr_r) || tr_r <= 0) {
    stop("Specify a positive 'tr' or explicit acquisition 'time' vector.", call. = FALSE)
  }
  off_r <- if (!is.null(time_offset)) {
    if (length(time_offset) > 1L) {
      if (!is.null(names(time_offset)) && r_name %in% names(time_offset)) {
        time_offset[[r_name]]
      } else {
        time_offset[r_idx]
      }
    } else {
      time_offset
    }
  } else if (!is.null(attr(mat, "time_offset", exact = TRUE))) {
    as.double(attr(mat, "time_offset", exact = TRUE))
  } else {
    0.0
  }
  (seq_len(n_t) - 1.0) * as.double(tr_r) + as.double(off_r)
}

#' Event-lock and interpolate MEDuSA time-series matrices onto an epoch grid
#'
#' @description
#' Matrix-native counterpart to [event_lock_ts()] and [get_medusa_interpolated_ts()]
#' that aligns and interpolates directly on `time x voxels` (`T x V`) numeric
#' matrices without expanding intermediate long tables or constructing `fmri_ts`
#' R6 objects.
#'
#' Preserves every audit, censoring, gap-fencing, collision/overlap, and kernel
#' quadrature guarantee of [get_medusa_interpolated_ts()], while running in C++
#' across single-run matrices, `"medusa_nifti_ts"` bundles from [read_nifti_ts()],
#' deconvolution outputs from [deconvolve_crossrun_hrf()], multi-run matrix lists,
#' multivariate signal lists, and legacy `fmri_ts` objects.
#'
#' @param ts Time-series input: a `T x V` numeric `matrix`, a `"medusa_nifti_ts"`
#'   object from [read_nifti_ts()], a deconvolution result list from
#'   [deconvolve_crossrun_hrf()], a named `list` of `T x V` matrices (multi-run or
#'   multivariate signals; see `list_mode`), or a legacy `fmri_ts` object.
#' @param events Trial-level `data.frame` (or named `list` of per-run `data.frame`s
#'   for multi-run inputs). Optional when `ts` is an `fmri_ts` object or a
#'   `deconvolve_crossrun_hrf()` fit whose stored events include `trial_col` and `event`.
#' @param event Name of the target event onset column in `events`.
#' @param time_before Earliest epoch lag relative to `event` in seconds (`<= 0`). Default `-3.0`.
#' @param time_after Latest epoch lag relative to `event` in seconds (`>= 0`). Default `3.0`.
#' @param collide_before Optional character vector of column names in `events`
#'   truncating earlier observations (`time > evt_before`) when `overlap_policy = "truncate"`.
#' @param collide_after Optional character vector of column names in `events`
#'   truncating later observations (`time < evt_after`) when `overlap_policy = "truncate"`.
#' @param pad_before Extra window padding (in seconds, `<= 0`) retained before
#'   `time_before` to supply interpolation brackets at the left boundary. Default `-1.5`.
#' @param pad_after Extra window padding (in seconds, `>= 0`) retained after
#'   `time_after` to supply interpolation brackets at the right boundary. Default `1.5`.
#' @param output_resolution Epoch output step in seconds. Defaults to `tr`.
#' @param tr Optional positive scalar repetition time in seconds. Inferred from
#'   `ts` (`ts$tr` or `attr(ts, "tr")`) when available.
#' @param time Optional numeric vector of length `T` (or per-run list) giving
#'   acquisition timestamps in seconds. Defaults to `(0:(T - 1)) * tr + time_offset`.
#' @param time_offset Optional scalar (or per-run numeric vector) start-time offset
#'   in seconds when `time` is not supplied. Defaults to `attr(mat, "time_offset")` or `0`.
#' @param parcel_ids Optional vector of length `V` assigning each spatial column
#'   to a parcel/ROI (populates both `parcel` and `atlas_value` for `group_by`).
#' @param unit_data Optional `V`-row `data.frame` of spatial unit metadata (e.g.,
#'   `ts$voxel_map`) whose columns may be referenced in `group_by`.
#' @param group_by Character vector of grouping variables, or `NULL`. For
#'   `align_medusa_ts()`, `missing(group_by)` defaults to `c(trial_col, "parcel")`
#'   when `parcel_ids` (or `ts$parcel_ids`) is present and `trial_col` otherwise;
#'   passing explicit `group_by = NULL` aggregates across trials and units.
#' @param value_name Signal name prefix used for single-signal matrix inputs
#'   (e.g., `"value"` produces `value_mean`, `value_sd`, ...). Default `"value"`.
#' @param trial_col Name of the unique within-run trial identifier column in
#'   `events`. Default `"trial"`.
#' @param run_col Name of the run identifier column used when `events` is a
#'   combined multi-run `data.frame` and attached to multi-run outputs. Default `"run"`.
#' @param valid Optional logical validity mask (`TRUE` = usable, `FALSE`/`NA` =
#'   censored): a length-`T` vector, a `T x V` logical matrix, a per-run list,
#'   or `"observation_mask"` to use `attr(mat, "observation_mask")`. For `fmri_ts`
#'   inputs, a single column name string is also accepted (matching `valid_col`).
#' @param overlap_policy One of `"allow"`, `"truncate"`, or `"exclude"`. `NULL`
#'   uses `"truncate"` when collision columns are supplied and `"allow"` otherwise.
#' @param neighbor_events Event columns used for `align_*` nearest-event
#'   annotations. Defaults to `event`.
#' @param max_gap Maximum allowed interpolation bracket duration in seconds
#'   (defaults to `tr`).
#' @param readout_kernel Readout kernel at each output lag: `"linear"` (default),
#'   `"boxcar"`, `"triangle"`, or `"gaussian"`.
#' @param readout_width Full window width in seconds for non-linear
#'   `readout_kernel` modes (defaults to `1.5 * max_gap`).
#' @param time_audit Logical; if `TRUE`, appends `evt_onset`, `collide_before`,
#'   `collide_after`, and absolute `time` when `trial_col` is in `group_by`.
#' @param list_mode Interpretation when `ts` is a plain `list` of matrices:
#'   `"auto"` (default), `"runs"` (multi-run list), or `"signals"` (multivariate
#'   signals for one run).
#' @param orientation Optional matrix orientation override (`"time_by_voxel"` or
#'   `"voxel_by_time"`).
#' @param ncpus Integer number of parallel workers used across runs when aligning
#'   multi-run inputs. Default `1L`.
#'
#' @return A `data.table` of aligned, interpolated epoch summaries and `align_*` /
#'   `epoch_status` diagnostics with a `source_clock` `data.frame` attribute.
#'
#' @details
#' `align_medusa_ts()` performs descriptive event-locked interpolation and
#' overlap/censoring auditing without assuming a forward generative design. When
#' `ts` is a [deconvolve_crossrun_hrf()] result and the goal is to remove
#' non-target event bleed, invert target-to-target cross-talk, model variable
#' durations, or propagate analytical per-trial standard errors (`se_demixed`),
#' use [decon_epoch_regression()] instead.
#'
#' @seealso [read_nifti_ts()], [deconvolve_crossrun_hrf()],
#'   [decon_epoch_regression()], [get_medusa_interpolated_ts()],
#'   [evt_align_decon_files()], [get_medusa_compression_score()]
#'
#' @examples
#' set.seed(1)
#' Y <- matrix(rnorm(60 * 8), nrow = 60, ncol = 8)
#' parcels <- rep(c("ROI_A", "ROI_B"), each = 4)
#' events <- data.frame(trial = 1:4, onset = c(8.5, 20.2, 33.7, 46.1))
#' aligned <- align_medusa_ts(
#'   Y, events = events, event = "onset", tr = 1.0,
#'   parcel_ids = parcels, time_before = -2, time_after = 4, output_resolution = 1.0
#' )
#' head(aligned)
#'
#' @export
align_medusa_ts <- function(
  ts,
  events = NULL,
  event = NULL,
  time_before = -3.0,
  time_after = 3.0,
  collide_before = NULL,
  collide_after = NULL,
  pad_before = -1.5,
  pad_after = 1.5,
  output_resolution = NULL,
  tr = NULL,
  time = NULL,
  time_offset = NULL,
  parcel_ids = NULL,
  unit_data = NULL,
  group_by = "trial",
  value_name = "value",
  trial_col = "trial",
  run_col = "run",
  valid = NULL,
  overlap_policy = NULL,
  neighbor_events = NULL,
  max_gap = NULL,
  readout_kernel = "linear",
  readout_width = NULL,
  time_audit = FALSE,
  list_mode = c("auto", "runs", "signals"),
  orientation = NULL,
  ncpus = 1L
) {
  group_by_missing <- missing(group_by)
  list_mode <- match.arg(list_mode)
  checkmate::assert_string(value_name, min.chars = 1L)
  checkmate::assert_string(trial_col, min.chars = 1L)
  checkmate::assert_string(run_col, min.chars = 1L)
  checkmate::assert_integerish(ncpus, lower = 1L, len = 1L, any.missing = FALSE)
  ncpus <- as.integer(ncpus)

  # 1. Legacy fmri_ts R6 input
  if (inherits(ts, "fmri_ts")) {
    valid_col <- if (is.character(valid) && length(valid) == 1L) valid else NULL
    if (!is.null(valid) && is.null(valid_col)) {
      stop("For fmri_ts objects, 'valid' must be NULL or a column name string.", call. = FALSE)
    }
    return(get_medusa_interpolated_ts(
      fmri_obj = ts,
      event = event,
      time_before = time_before,
      time_after = time_after,
      collide_before = collide_before,
      collide_after = collide_after,
      pad_before = pad_before,
      pad_after = pad_after,
      output_resolution = output_resolution,
      group_by = if (group_by_missing) ts$vm$trial else group_by,
      time_audit = time_audit,
      overlap_policy = overlap_policy,
      valid_col = valid_col,
      neighbor_events = neighbor_events,
      max_gap = max_gap,
      readout_kernel = readout_kernel,
      readout_width = readout_width
    ))
  }

  # 2. Unpack structured MEDuSA containers ("medusa_nifti_ts" or deconvolve_crossrun_hrf output)
  is_multi_run_container <- FALSE
  if (inherits(ts, "medusa_nifti_ts")) {
    if (is.null(tr) && is.finite(ts$tr)) tr <- ts$tr
    if (is.null(unit_data) && !is.null(ts$voxel_map)) unit_data <- ts$voxel_map
    if (is.null(parcel_ids) && !is.null(ts$parcel_ids)) parcel_ids <- ts$parcel_ids
    if (is.null(orientation) && !is.null(ts$orientation)) orientation <- ts$orientation
    ts <- ts$bold
    is_multi_run_container <- TRUE
  } else if (inherits(ts, "crossrun_decon") ||
             (is.list(ts) && !is.data.frame(ts) && all(c("activity", "hrf_fit") %in% names(ts)))) {
    if (is.null(tr) && !is.null(ts$hrf_fit$settings$tr)) tr <- ts$hrf_fit$settings$tr
    if (is.null(events) && !is.null(ts$hrf_fit$specification$events)) {
      events <- ts$hrf_fit$specification$events
    }
    if (is.null(parcel_ids) && !is.null(ts$hrf_fit$specification$spatial_groups)) {
      parcel_ids <- ts$hrf_fit$specification$spatial_groups
    }
    if (is.null(unit_data)) {
      vm_dec <- if (!is.null(ts$voxel_map)) ts$voxel_map else ts$hrf_fit$voxel_map
      if (!is.null(vm_dec)) unit_data <- vm_dec
    }
    ts <- ts$activity
    is_multi_run_container <- TRUE
  }

  # Default group_by when omitted on matrix-native calls
  if (group_by_missing) {
    if (!is.null(parcel_ids) || (is.matrix(ts) && !is.null(attr(ts, "parcel_ids", exact = TRUE))) ||
        (!is.null(unit_data) && "parcel" %in% names(unit_data))) {
      group_by <- c(trial_col, "parcel")
    } else if (!is.null(unit_data) && "atlas_value" %in% names(unit_data)) {
      group_by <- c(trial_col, "atlas_value")
    } else {
      group_by <- trial_col
    }
  }

  # 3. Single numeric matrix
  if (is.matrix(ts)) {
    mat <- normalize_matrix_orientation(ts, orientation = orientation, label = "ts")
    if (is.null(events)) stop("'events' must be supplied when aligning a matrix.", call. = FALSE)
    if (is.null(unit_data) && !is.null(attr(mat, "voxel_map", exact = TRUE))) {
      unit_data <- attr(mat, "voxel_map", exact = TRUE)
    }
    if (is.null(parcel_ids) && !is.null(attr(mat, "parcel_ids", exact = TRUE))) {
      parcel_ids <- attr(mat, "parcel_ids", exact = TRUE)
    }
    tr_eff <- if (!is.null(tr)) tr else attr(mat, "tr", exact = TRUE)
    udf <- build_unit_metadata_df(ncol(mat), colnames(mat), unit_data = unit_data, parcel_ids = parcel_ids)
    t_vec <- resolve_run_time_vec(time, 1L, "run1", nrow(mat), tr_eff, time_offset, mat)
    v_mask <- resolve_run_valid_mask(valid, 1L, "run1", mat)
    sigs <- stats::setNames(list(mat), value_name)

    return(align_matrix_single_run(
      signals = sigs,
      raw_time = t_vec,
      n_time_raw = nrow(mat),
      n_units = ncol(mat),
      unit_df = udf,
      events = events,
      event = event,
      trial_col = trial_col,
      time_before = time_before,
      time_after = time_after,
      collide_before = collide_before,
      collide_after = collide_after,
      pad_before = pad_before,
      pad_after = pad_after,
      output_resolution = output_resolution,
      tr = tr_eff,
      group_by = group_by,
      time_audit = time_audit,
      overlap_policy = overlap_policy,
      raw_valid = v_mask,
      neighbor_events = neighbor_events,
      max_gap = max_gap,
      readout_kernel = readout_kernel,
      readout_width = readout_width
    ))
  }

  # 4. List of matrices (either multi-run or multivariate single-run)
  if (is.list(ts) && !is.data.frame(ts) && length(ts) >= 1L) {
    mats <- lapply(seq_along(ts), function(i) {
      normalize_matrix_orientation(ts[[i]], orientation = orientation, label = paste0("ts[[", i, "]]"))
    })
    names(mats) <- names(ts)

    mode_resolved <- if (is_multi_run_container) {
      "runs"
    } else if (!identical(list_mode, "auto")) {
      list_mode
    } else {
      same_dims <- length(unique(lapply(mats, dim))) == 1L
      ev_is_run_list <- is.list(events) && !is.data.frame(events)
      ev_multi_runs <- is.data.frame(events) && (run_col %in% names(events)) &&
        length(unique(events[[run_col]])) > 1L
      if (ev_is_run_list || ev_multi_runs || !same_dims) {
        "runs"
      } else if (length(mats) == 1L && (is.null(names(mats)) || !nzchar(names(mats)[1L]))) {
        "runs"
      } else {
        "signals"
      }
    }

    if (identical(mode_resolved, "signals")) {
      if (is.null(events) || !is.data.frame(events)) {
        stop("Multivariate signal alignment requires a single-run 'events' data.frame.", call. = FALSE)
      }
      if (is.null(names(mats)) || any(!nzchar(names(mats))) || anyDuplicated(names(mats))) {
        stop("Multivariate signal list 'ts' must have unique, nonempty signal names.", call. = FALSE)
      }
      dims <- lapply(mats, dim)
      if (length(unique(dims)) != 1L) {
        stop("All signal matrices in multivariate 'ts' must share identical T x V dimensions.", call. = FALSE)
      }
      ref_mat <- mats[[1L]]
      if (is.null(unit_data) && !is.null(attr(ref_mat, "voxel_map", exact = TRUE))) {
        unit_data <- attr(ref_mat, "voxel_map", exact = TRUE)
      }
      if (is.null(parcel_ids) && !is.null(attr(ref_mat, "parcel_ids", exact = TRUE))) {
        parcel_ids <- attr(ref_mat, "parcel_ids", exact = TRUE)
      }
      tr_eff <- if (!is.null(tr)) tr else attr(ref_mat, "tr", exact = TRUE)
      udf <- build_unit_metadata_df(ncol(ref_mat), colnames(ref_mat), unit_data = unit_data, parcel_ids = parcel_ids)
      t_vec <- resolve_run_time_vec(time, 1L, "run1", nrow(ref_mat), tr_eff, time_offset, ref_mat)
      v_mask <- resolve_run_valid_mask(valid, 1L, "run1", ref_mat)

      return(align_matrix_single_run(
        signals = mats,
        raw_time = t_vec,
        n_time_raw = nrow(ref_mat),
        n_units = ncol(ref_mat),
        unit_df = udf,
        events = events,
        event = event,
        trial_col = trial_col,
        time_before = time_before,
        time_after = time_after,
        collide_before = collide_before,
        collide_after = collide_after,
        pad_before = pad_before,
        pad_after = pad_after,
        output_resolution = output_resolution,
        tr = tr_eff,
        group_by = group_by,
        time_audit = time_audit,
        overlap_policy = overlap_policy,
        raw_valid = v_mask,
        neighbor_events = neighbor_events,
        max_gap = max_gap,
        readout_kernel = readout_kernel,
        readout_width = readout_width
      ))
    }

    # Multi-run alignment
    if (is.null(events)) {
      stop("'events' must be supplied for multi-run alignment.", call. = FALSE)
    }
    run_names <- names(mats)
    if (is.null(run_names) || any(!nzchar(run_names))) {
      run_names <- paste0("run", seq_along(mats))
      names(mats) <- run_names
    } else if (anyDuplicated(run_names)) {
      stop("Multi-run matrix list must have unique run names.", call. = FALSE)
    }

    ev_by_run <- split_multirun_events(events, run_names, run_col = run_col)
    group_by_inner <- setdiff(group_by, run_col)

    run_worker <- function(r_idx) {
      r_name <- run_names[r_idx]
      mat_r <- mats[[r_idx]]
      ev_r <- ev_by_run[[r_idx]]
      ud_r <- if (!is.null(unit_data)) unit_data else attr(mat_r, "voxel_map", exact = TRUE)
      pid_r <- if (!is.null(parcel_ids)) parcel_ids else attr(mat_r, "parcel_ids", exact = TRUE)
      tr_r <- if (!is.null(tr)) tr else attr(mat_r, "tr", exact = TRUE)
      udf_r <- build_unit_metadata_df(ncol(mat_r), colnames(mat_r), unit_data = ud_r, parcel_ids = pid_r)
      t_vec_r <- resolve_run_time_vec(time, r_idx, r_name, nrow(mat_r), tr_r, time_offset, mat_r)
      v_mask_r <- resolve_run_valid_mask(valid, r_idx, r_name, mat_r)

      res_r <- align_matrix_single_run(
        signals = stats::setNames(list(mat_r), value_name),
        raw_time = t_vec_r,
        n_time_raw = nrow(mat_r),
        n_units = ncol(mat_r),
        unit_df = udf_r,
        events = ev_r,
        event = event,
        trial_col = trial_col,
        time_before = time_before,
        time_after = time_after,
        collide_before = collide_before,
        collide_after = collide_after,
        pad_before = pad_before,
        pad_after = pad_after,
        output_resolution = output_resolution,
        tr = tr_r,
        group_by = group_by_inner,
        time_audit = time_audit,
        overlap_policy = overlap_policy,
        raw_valid = v_mask_r,
        neighbor_events = neighbor_events,
        max_gap = max_gap,
        readout_kernel = readout_kernel,
        readout_width = readout_width
      )

      run_val <- if (run_col %in% names(ev_r) && length(unique(ev_r[[run_col]])) == 1L) {
        ev_r[[run_col]][1L]
      } else {
        r_name
      }
      clk_r <- attr(res_r, "source_clock")
      clk_r[[run_col]] <- rep(run_val, nrow(clk_r))
      res_r[, (run_col) := rep(run_val, .N)]
      data.table::setcolorder(res_r, c(run_col, setdiff(names(res_r), run_col)))
      list(ts = res_r, clock = clk_r)
    }

    parts <- medusa_parallel_lapply(seq_along(mats), run_worker, ncpus = ncpus)
    combined <- data.table::rbindlist(lapply(parts, "[[", "ts"))
    combined_clock <- as.data.frame(data.table::rbindlist(lapply(parts, "[[", "clock")))
    data.table::setattr(combined, "source_clock", combined_clock)
    return(combined)
  }

  stop("Unsupported 'ts' input type. Supply a numeric matrix, 'medusa_nifti_ts', ",
       "deconvolution fit, list of matrices, or 'fmri_ts' object.", call. = FALSE)
}

compression_score_matrix_single_run <- function(
  signals,
  raw_time,
  n_time_raw,
  n_units,
  unit_df,
  events,
  event,
  trial_col = "trial",
  time_before = -3,
  time_after = 3,
  collide_before = NULL,
  collide_after = NULL,
  group_by = NULL,
  pexp_target = 0.9,
  scale_columns = TRUE,
  coerce_factor_trial = FALSE
) {
  checkmate::assert_string(event)
  checkmate::assert_number(time_before, upper = 0, finite = TRUE)
  checkmate::assert_number(time_after, lower = 0, finite = TRUE)
  checkmate::assert_character(collide_before, null.ok = TRUE, unique = TRUE)
  checkmate::assert_character(collide_after, null.ok = TRUE, unique = TRUE)
  checkmate::assert_data_frame(events)
  checkmate::assert_subset(c(event, collide_before, collide_after), names(events))
  if (nrow(events) == 0L) stop("No trial events supplied for alignment")
  if (anyNA(events[[trial_col]]) || anyDuplicated(events[[trial_col]])) {
    stop("Event data must contain unique, nonmissing trial identifiers")
  }

  overlap_policy <- resolve_overlap_policy(NULL, collide_before, collide_after)
  fin_mask <- is.finite(raw_time)
  fin_idx <- which(fin_mask)
  fin_time <- as.double(raw_time[fin_idx])
  ord <- order(fin_time)
  clock <- fin_time[ord]
  keep_rows <- fin_idx[ord]

  sig_mats <- lapply(signals, function(x) {
    x <- as.double(x)
    dim(x) <- c(n_time_raw, n_units)
    if (length(keep_rows) == n_time_raw && identical(ord, seq_len(n_time_raw))) {
      x
    } else {
      x[keep_rows, , drop = FALSE]
    }
  })
  sig_names <- names(sig_mats)

  run_df <- data.table::as.data.table(data.table::copy(events))
  if (isTRUE(coerce_factor_trial) && is.factor(run_df[[trial_col]])) {
    run_df[[trial_col]] <- as.character(run_df[[trial_col]])
  }
  data.table::setorderv(run_df, trial_col)
  n_trials <- nrow(run_df)
  onsets <- as.double(run_df[[event]])

  evt_before <- rep(-Inf, n_trials)
  evt_after <- rep(Inf, n_trials)
  if (identical(overlap_policy, "truncate")) {
    if (length(collide_before)) {
      bvals <- sort(unique(as.double(unlist(run_df[, collide_before, with = FALSE], use.names = FALSE))))
      bvals <- bvals[is.finite(bvals)]
      if (length(bvals)) {
        fin_o <- which(is.finite(onsets))
        idx_b <- findInterval(onsets[fin_o], bvals, left.open = TRUE)
        has_b <- idx_b > 0L
        evt_before[fin_o[has_b]] <- bvals[idx_b[has_b]]
      }
    }
    if (length(collide_after)) {
      avals <- sort(unique(as.double(unlist(run_df[, collide_after, with = FALSE], use.names = FALSE))))
      avals <- avals[is.finite(avals)]
      if (length(avals)) {
        fin_o <- which(is.finite(onsets))
        idx_a <- findInterval(onsets[fin_o], avals)
        has_a <- idx_a < length(avals)
        evt_after[fin_o[has_a]] <- avals[idx_a[has_a] + 1L]
      }
    }
  }

  avail_group_cols <- unique(c(trial_col, names(unit_df)))
  checkmate::assert_subset(group_by, avail_group_cols)
  spatial_cols <- setdiff(group_by, trial_col)

  if (length(spatial_cols) == 0L) {
    group_dt <- NULL
    group_cols_list <- list(seq_len(n_units))
    n_groups <- 1L
  } else {
    u_dt <- data.table::as.data.table(unit_df[, spatial_cols, drop = FALSE])
    group_dt <- unique(u_dt)
    data.table::setorderv(group_dt, spatial_cols)
    group_dt[, .gid := seq_len(.N)]
    u_dt[, .u_ord := seq_len(.N)]
    u_matched <- merge(u_dt, group_dt, by = spatial_cols, sort = FALSE)
    data.table::setorderv(u_matched, ".u_ord")
    group_cols_list <- split(seq_len(n_units), u_matched$.gid)
    group_dt[, .gid := NULL]
    n_groups <- nrow(group_dt)
  }

  # Trial window index slices (pad_before = 0, pad_after = 0 in get_medusa_compression_score)
  win_slices <- vector("list", n_trials)
  for (t in seq_len(n_trials)) {
    onset <- onsets[t]
    if (!is.finite(onset) || length(clock) == 0L) {
      win_slices[[t]] <- integer(0L)
    } else {
      rt <- clock - onset
      win_slices[[t]] <- which(
        rt >= time_before & rt <= time_after & clock > evt_before[t] & clock < evt_after[t]
      )
    }
  }

  n_combos <- n_trials * n_groups
  t_idx_vec <- rep(seq_len(n_trials), each = n_groups)
  g_idx_vec <- rep(seq_len(n_groups), times = n_trials)

  out_cols <- list()
  for (col in group_by) {
    if (identical(col, trial_col)) {
      out_cols[[col]] <- run_df[[trial_col]][t_idx_vec]
    } else {
      out_cols[[col]] <- group_dt[[col]][g_idx_vec]
    }
  }

  metric_keys <- c(
    as.vector(rbind(paste0("compress_", pexp_target), paste0("acompress_", pexp_target))),
    "compress_nt"
  )
  # Match dcast column ordering: metric-major, signal-minor
  sig_metrics <- vector("list", length(sig_names))
  for (s_idx in seq_along(sig_names)) {
    ymat <- sig_mats[[s_idx]]
    res_list <- vector("list", n_combos)
    for (k in seq_len(n_combos)) {
      rows_k <- win_slices[[t_idx_vec[k]]]
      cols_k <- group_cols_list[[g_idx_vec[k]]]
      submat <- if (length(rows_k) == 0L) {
        matrix(NA_real_, nrow = 1L, ncol = length(cols_k))
      } else {
        ymat[rows_k, cols_k, drop = FALSE]
      }
      res_list[[k]] <- compress_mts_pca(submat, pexp_target = pexp_target, scale_columns = scale_columns)
    }
    sig_metrics[[s_idx]] <- data.table::rbindlist(res_list)
  }

  for (mk in metric_keys) {
    for (s_idx in seq_along(sig_names)) {
      out_cols[[paste0(sig_names[s_idx], "_", mk)]] <- sig_metrics[[s_idx]][[mk]]
    }
  }

  data.table::setDT(out_cols)
  data.table::setkeyv(out_cols, group_by)
  out_cols
}
