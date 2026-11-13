#' function to get interpolated event locked data
#' @param fmri_obj An fmri_ts object containing time series data that are optionally keyed on one or more
#'   grouping variables (e.g., ROI). The fmri_obj must also contain trial-level data (in the $event_data field)
#'   for event-locking to proceed.
#' @param event Name of column in \code{fmri_obj$event_data} that identifies onset for event of interest
#' @param time_before How many seconds before the \code{event} do we want data. Default: -3
#' @param time_after How many seconds after the \code{event} do we want data. Default: 3
#' @param collide_before An optional vector of column names in \code{trial_df} that set a boundary on the earliest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_before}.
#' @param collide_after An optional vector of column names in \code{trial_df} that set a boundary on the latest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_after}.
#' @param pad_before Number of seconds to include in the epoch time window before the event of interest. Interpolation
#'   spans the window from \code{time_before} to \code{time_after}, but padding includes data points at the
#'   boundary that can help to have sufficient data to interpolate early and late times within the epoch.
#' @param pad_after Number of seconds to include in the epoch time window after the event of interest.
#' @param time_audit If TRUE, additional columns will be added to the output showing how alignment is calculated vis-a-vis event timing
#' @param overlap_policy One of \code{"allow"}, \code{"truncate"}, or \code{"exclude"}.
#'   NULL uses truncation when collision columns are supplied and allows overlap
#'   otherwise. Explicit allow retains full windows. Exclude removes whole epochs
#'   whose closed requested windows overlap another target epoch, including ties.
#' @param valid_col Optional logical source-data column. FALSE or NA censors that
#'   row for all value columns without deleting its acquisition time.
#' @param neighbor_events Event columns used for nearest-event annotations;
#'   defaults to the target event column. Does not itself censor observations.
#' @details Nonfinite target onsets and unsupported windows retain placeholder
#'   rows with missing time and values and an \code{epoch_status} column.
#'   Nonfinite collision candidates are ignored. Nearest-event boundaries remain
#'   strict and do not guarantee disjoint trial windows. For regular shared-clock
#'   \code{time x voxels} matrices or \code{"medusa_nifti_ts"} inputs, prefer
#'   \code{\link{align_medusa_ts}} to align and interpolate directly in C++
#'   without expanding trial-by-voxel \code{fmri_ts} tables.
#'
#' @seealso \code{\link{align_medusa_ts}}, \code{\link{get_medusa_interpolated_ts}}
#'
#' @importFrom dplyr filter pull
#' @importFrom data.table rbindlist
#' @importFrom checkmate assert_numeric assert_data_frame assert_string assert_subset
#'
#' @export
event_lock_ts <- function(fmri_obj, event=NULL, time_before=-3, time_after=3,
                          collide_before=NULL, collide_after=NULL, pad_before=-1, pad_after=1, time_audit=FALSE,
                          overlap_policy=NULL, valid_col=NULL, neighbor_events=NULL) {

  ts_data <- fmri_obj$get_ts(orig_names=TRUE) #so that vm pass-through to returned object works
  run_df <- fmri_obj$get_events(orig_names = TRUE)
  vm <- fmri_obj$vm
  overlap_policy <- resolve_overlap_policy(overlap_policy, collide_before, collide_after)
  checkmate::assert_string(valid_col, null.ok = TRUE)
  if (is.null(neighbor_events)) neighbor_events <- event
  checkmate::assert_character(neighbor_events, min.len = 1L, unique = TRUE)

  checkmate::assert_number(time_before, upper = 0, finite = TRUE)
  checkmate::assert_number(time_after, lower = 0, finite = TRUE)
  checkmate::assert_number(pad_before, upper = 0, finite = TRUE)
  checkmate::assert_number(pad_after, lower = 0, finite = TRUE)
  checkmate::assert_data_frame(ts_data)
  checkmate::assert_data_frame(run_df)
  if (nrow(ts_data) == 0L) stop("No time series samples supplied for event alignment")
  checkmate::assert_string(event)
  checkmate::assert_subset(c(event, collide_before, collide_after, neighbor_events), names(run_df))
  if (nrow(run_df) == 0L) stop("No trial events supplied for alignment")
  if (anyNA(run_df[[vm$trial]]) || anyDuplicated(run_df[[vm$trial]])) {
    stop("Event data must contain unique, nonmissing trial identifiers")
  }
  for (column in unique(c(event, collide_before, collide_after, neighbor_events))) {
    checkmate::assert_numeric(run_df[[column]], all.missing = TRUE)
  }
  run_df <- run_df[order(run_df[[vm$trial]]), ]
  trials <- run_df[[vm$trial]]
  if (any(grepl("^(align_|\\.align_)", names(ts_data)))) {
    stop("Input contains reserved alignment columns")
  }
  clock <- sort(unique(ts_data[[vm$time]][is.finite(ts_data[[vm$time]])]))
  ts_data$.align_source_id <- match(ts_data[[vm$time]], clock)
  ts_data$.align_censored <- FALSE
  if (!is.null(valid_col)) {
    checkmate::assert_subset(valid_col, names(ts_data))
    if (!is.logical(ts_data[[valid_col]])) stop("valid_col must name a logical column")
    ts_data$.align_censored <- is.na(ts_data[[valid_col]]) | !ts_data[[valid_col]]
    for (column in vm$value) ts_data[[column]][ts_data$.align_censored] <- NA_real_
  }
  annotations <- epoch_event_annotations(run_df, event, vm$trial, neighbor_events,
    time_before, time_after, if (length(clock)) min(clock) else NA_real_,
    if (length(clock)) max(clock) else NA_real_, overlap_policy)
  finite_events <- function(columns) {
    values <- unlist(run_df[, columns, with = FALSE], use.names = FALSE)
    values[is.finite(values)]
  }
  before_events <- finite_events(collide_before)
  after_events <- finite_events(collide_after)
  tlist <- vector("list", length(trials))

  # A placeholder per original unit preserves missing/fully censored epochs.
  # Its missing time and values cannot be mistaken for observed samples.
  units <- if (length(vm$key)) {
    ts_data[!duplicated(ts_data[, vm$key, with = FALSE]), ]
  } else ts_data[seq_len(min(1L, nrow(ts_data))), ]
  observed_times <- ts_data[[vm$time]]
  observed_times <- observed_times[is.finite(observed_times)]

  for (t in seq_along(trials)) {
    evt_onset <- run_df[[event]][t]
    evt_before <- -Inf
    evt_after <- Inf
    status <- "ok"
    if (!is.finite(evt_onset)) {
      status <- "missing_onset"
      trial_ts <- units[0L, ]
    } else {
      past <- before_events[before_events < evt_onset]
      future <- after_events[after_events > evt_onset]
      if (overlap_policy == "truncate") {
        if (length(past)) evt_before <- max(past)
        if (length(future)) evt_after <- min(future)
      }
      relative_time <- ts_data[[vm$time]] - evt_onset
      keep <- is.finite(relative_time) &
        relative_time >= time_before + pad_before & relative_time <= time_after + pad_after &
        ts_data[[vm$time]] > evt_before & ts_data[[vm$time]] < evt_after
      if (overlap_policy == "exclude" && annotations$align_window_overlap_n[t] > 0L) {
        keep[] <- FALSE
        status <- "overlap_excluded"
      }
      trial_ts <- data.table::copy(ts_data[which(keep), ])
      trial_ts$evt_time <- relative_time[which(keep)]
      if (nrow(trial_ts) == 0L && status == "ok") {
        status <- if (!length(observed_times)) "no_data" else if (
          evt_onset + time_after + pad_after < min(observed_times) ||
          evt_onset + time_before + pad_before > max(observed_times)
        ) "outside_run" else "censored"
      }
    }
    if (nrow(trial_ts) == 0L) {
      trial_ts <- data.table::copy(units)
      trial_ts[[vm$time]] <- rep(NA_real_, nrow(trial_ts))
      for (column in vm$value) trial_ts[[column]] <- rep(NA_real_, nrow(trial_ts))
      trial_ts$evt_time <- rep(NA_real_, nrow(trial_ts))
      trial_ts$.align_source_id <- NA_integer_
      trial_ts$.align_censored <- FALSE
    }
    trial_ts[[vm$trial]] <- rep(trials[t], nrow(trial_ts))
    trial_ts$epoch_status <- rep(status, nrow(trial_ts))
    trial_ts$.align_epoch_id <- rep(t, nrow(trial_ts))
    for (column in names(annotations)) trial_ts[[column]] <- rep(annotations[[column]][t], nrow(trial_ts))
    if (isTRUE(time_audit)) {
      trial_ts$evt_onset <- rep(evt_onset, nrow(trial_ts))
      trial_ts$collide_before <- rep(evt_before, nrow(trial_ts))
      trial_ts$collide_after <- rep(evt_after, nrow(trial_ts))
    }
    tlist[[t]] <- trial_ts
  }

  names(tlist) <- as.character(trials)
  tdf <- data.table::rbindlist(tlist)
  res <- fmri_ts$new(ts_data = tdf, event_data = run_df, vm = vm, tr = fmri_obj$tr)
  res$replace_vm(time = "evt_time")
  res$add_keys(vm$trial)
  attr(res, "source_clock") <- data.frame(source_id = seq_along(clock), time = clock)

  return(res)

}

#' front-end function for taking a list of windowed time series by mask value, interpolating them onto a time grid,
#' and (optionally) averaging across voxels/units within a value to derive the mean interpolated time series
#'
#' @param fmri_obj an fmri_ts object containing a single run of data with corresponding events
#' @param event the event to which the fmri time series should be aligned (column in \code{fmri_obj$event_data})
#' @param time_before How many seconds before the \code{event} do we want data. Default: -3
#' @param time_after How many seconds after the \code{event} do we want data. Default: 3
#' @param collide_before An optional vector of column names in \code{trial_df} that set a boundary on the earliest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_before}.
#' @param collide_after An optional vector of column names in \code{trial_df} that set a boundary on the latest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_after}.
#' @param pad_before Number of seconds to include in the epoch time window before the event of interest. Interpolation
#'   spans the window from \code{time_before} to \code{time_after}, but padding includes data points at the
#'   boundary that can help provide sufficient data to interpolate early and late times within the epoch.
#' @param pad_after Number of seconds to include in the epoch time window after the event of interest.
#' @param output_resolution Time step, in seconds, for the interpolated data. Defaults to the repetition time
#'   stored in \code{fmri_obj}.
#' @param group_by return interpolated time series for each combination of group_by variables. Default is
#'   to provide one interpolated time series per trial.
#' @param time_audit If TRUE, additional columns will be added to the output showing how alignment is calculated vis-a-vis event timing
#' @param overlap_policy One of \code{"allow"}, \code{"truncate"}, or \code{"exclude"}.
#'   NULL uses truncation when collision columns are supplied and allows overlap
#'   otherwise. Explicit allow retains full windows. Exclude removes whole epochs
#'   whose closed requested windows overlap another target epoch, including ties.
#' @param valid_col Optional logical source-data column. FALSE or NA censors that
#'   row for all value columns without deleting its acquisition time.
#' @param neighbor_events Event columns used for nearest-event annotations;
#'   defaults to the target event column. Does not itself censor observations.
#' @param max_gap Maximum interpolation bracket duration in seconds, defaulting
#'   to the acquisition TR. Explicit invalid rows always block interpolation;
#'   this tolerance only permits wider brackets across deleted rows.
#' @param readout_kernel Readout kernel applied at each output grid time:
#'   \code{"linear"} (default two-point bracket interpolation), \code{"boxcar"},
#'   \code{"triangle"}, or \code{"gaussian"} local window integration within each
#'   contiguous uncensored segment.
#' @param readout_width Full window width (in seconds) for non-linear
#'   \code{readout_kernel} modes. Defaults to \code{1.5 * max_gap}.
#' @param events Optional trial-level \code{data.frame} (or list of per-run
#'   \code{data.frame}s) used when \code{fmri_obj} is a matrix, \code{"medusa_nifti_ts"},
#'   or deconvolution fit rather than an \code{fmri_ts} object.
#' @param tr Optional repetition time in seconds when \code{fmri_obj} is a matrix
#'   without a \code{tr} attribute.
#' @param time Optional numeric vector of acquisition timestamps (length \code{nrow(fmri_obj)})
#'   for matrix inputs.
#' @param time_offset Optional scalar start-time offset in seconds for matrix inputs.
#' @param parcel_ids Optional vector assigning matrix columns to parcels/ROIs.
#' @param unit_data Optional \code{data.frame} of unit metadata (one row per matrix column).
#' @param value_name Signal name prefix for single-signal matrix inputs. Default \code{"value"}.
#' @param trial_col Trial identifier column name for matrix inputs. Default \code{"trial"}.
#' @param use_legacy Logical; if \code{TRUE}, forces the legacy \code{event_lock_ts()} +
#'   \code{interpolate_fmri_epochs()} path on \code{fmri_ts} inputs (primarily for testing).
#' @details The interpolation grid includes zero and remains within the requested
#'   interval, including zero-width and one-sided windows. When trial is a grouping
#'   variable, missing and unsupported epochs retain missing estimates and an
#'   \code{epoch_status} column: \code{"ok"}, \code{"partial"},
#'   \code{"insufficient_data"}, \code{"missing_onset"}, \code{"outside_run"},
#'   \code{"censored"}, or \code{"no_data"}. These describe timing and available
#'   output support, not physiological isolation from other trials.
#'   Interpolation occurs within each unit before aggregation, preserving explicit
#'   missing rows and acquisition gaps. Off-grid median/SD now summarize the
#'   interpolated units. Each value has contributor/reason counts, source IDs, and
#'   the maximum distinct-epoch reuse count. Source IDs index observed acquisition
#'   times within this run; reuse is computed within the returned grouping.
#'   The \code{source_clock} attribute maps IDs to absolute acquisition times.
#'   Epoch annotations beginning with \code{align_} are always included when trial
#'   is a grouping variable. They describe timing, not physiological attribution.
#'   Whole-window exclusion adds \code{"overlap_excluded"} to epoch statuses.
#'   When spatial units share a common acquisition clock, \code{fmri_ts} inputs
#'   are aligned and interpolated via a zero-copy C++ matrix view, automatically
#'   falling back to the ragged per-unit \code{data.table} engine when unit
#'   clocks differ (e.g., slice-specific timestamps).
#'
#' @seealso \code{\link{align_medusa_ts}}, \code{\link{evt_align_decon_files}},
#'   \code{\link{decon_epoch_regression}}
#'
#' @author Michael Hallquist
#' @export
get_medusa_interpolated_ts <- function(fmri_obj, event=NULL, time_before=-3.0, time_after=3.0,
                                       collide_before=NULL, collide_after=NULL,
                                       pad_before=-1.5, pad_after=1.5, output_resolution=NULL,
                                       group_by="trial", time_audit=FALSE, overlap_policy=NULL,
                                       valid_col=NULL, neighbor_events=NULL, max_gap=NULL,
                                       readout_kernel="linear", readout_width=NULL,
                                       events=NULL, tr=NULL, time=NULL, time_offset=NULL,
                                       parcel_ids=NULL, unit_data=NULL, value_name="value",
                                       trial_col="trial", use_legacy=FALSE) {

  if (!inherits(fmri_obj, "fmri_ts")) {
    return(align_medusa_ts(
      ts = fmri_obj,
      events = events,
      event = event,
      time_before = time_before,
      time_after = time_after,
      collide_before = collide_before,
      collide_after = collide_after,
      pad_before = pad_before,
      pad_after = pad_after,
      output_resolution = output_resolution,
      tr = tr,
      time = time,
      time_offset = time_offset,
      parcel_ids = parcel_ids,
      unit_data = unit_data,
      group_by = group_by,
      value_name = value_name,
      trial_col = trial_col,
      valid = valid_col,
      overlap_policy = overlap_policy,
      neighbor_events = neighbor_events,
      max_gap = max_gap,
      readout_kernel = readout_kernel,
      readout_width = readout_width,
      time_audit = time_audit
    ))
  }

  #make default output_resolution equal to fmri TR
  if (is.null(output_resolution)) { output_resolution <- fmri_obj$tr }

  if (!isTRUE(use_legacy)) {
    view <- extract_fmri_ts_matrix_view(fmri_obj, valid_col = valid_col, group_by = group_by)
    if (!is.null(view)) {
      return(align_matrix_single_run(
        signals = view$signals,
        raw_time = view$raw_time,
        n_time_raw = view$n_time_raw,
        n_units = view$n_units,
        unit_df = view$unit_df,
        events = view$events,
        event = event,
        trial_col = view$trial_col,
        time_before = time_before,
        time_after = time_after,
        collide_before = collide_before,
        collide_after = collide_after,
        pad_before = pad_before,
        pad_after = pad_after,
        output_resolution = output_resolution,
        tr = view$tr,
        group_by = group_by,
        time_audit = time_audit,
        overlap_policy = overlap_policy,
        raw_valid = view$raw_valid,
        neighbor_events = neighbor_events,
        max_gap = max_gap,
        readout_kernel = readout_kernel,
        readout_width = readout_width,
        coerce_factor_trial = TRUE
      ))
    }
  }

  # talign is an fmri_ts object keyed by trial (and other keying variables)
  talign <- event_lock_ts(fmri_obj,
    event = event, time_before = time_before, time_after = time_after,
    pad_before = pad_before, pad_after = pad_after, collide_before = collide_before, collide_after = collide_after,
    time_audit = time_audit, overlap_policy = overlap_policy,
    valid_col = valid_col, neighbor_events = neighbor_events
  )

  #need to interpolate by key variables
  interpolated_epochs <- interpolate_fmri_epochs(talign, time_before=time_before, time_after=time_after,
    output_resolution=output_resolution, group_by=group_by, time_audit = time_audit, max_gap = max_gap,
    readout_kernel = readout_kernel, readout_width = readout_width)

  return(interpolated_epochs)
}

#' Compute event-locked multivariate PCA compression scores across spatial units
#'
#' @param fmri_obj an fmri_ts object, numeric \code{T x V} matrix, or single-run \code{"medusa_nifti_ts"} object
#' @param event the event to which the fmri time series should be aligned (column in \code{fmri_obj$event_data} or \code{events})
#' @param time_before How many seconds before the \code{event} do we want data. Default: -3
#' @param time_after How many seconds after the \code{event} do we want data. Default: 3
#' @param collide_before An optional vector of column names in \code{trial_df} that set a boundary on the earliest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_before}.
#' @param collide_after An optional vector of column names in \code{trial_df} that set a boundary on the latest
#'   time point used in interpolation. This effectively truncates the time series for interpolation to a smaller window
#'   than specified by \code{time_after}.
#' @param group_by return compression scores for each combination of group_by variables.
#' @param events Optional trial-level \code{data.frame} for matrix inputs.
#' @param tr Optional repetition time in seconds for matrix inputs.
#' @param time Optional numeric vector of acquisition timestamps for matrix inputs.
#' @param time_offset Optional scalar start-time offset in seconds for matrix inputs.
#' @param parcel_ids Optional vector assigning matrix columns to parcels/ROIs.
#' @param unit_data Optional \code{data.frame} of unit metadata for matrix inputs.
#' @param value_name Signal name prefix for single-signal matrix inputs. Default \code{"value"}.
#' @param trial_col Trial identifier column name for matrix inputs. Default \code{"trial"}.
#' @param pexp_target Target cumulative variance proportion(s) passed to \code{compress_mts_pca()}. Default \code{0.9}.
#' @param scale_columns Logical; standardize columns prior to SVD. Default \code{TRUE}.
#' @param use_legacy Logical; force legacy \code{fmri_ts} path when \code{TRUE}.
#'
#' @importFrom checkmate assert_string assert_numeric assert_character
#' @importFrom data.table dcast melt
#' @author Michael Hallquist
#' @export
get_medusa_compression_score <- function(fmri_obj, event=NULL, time_before=-3, time_after=3,
                                         collide_before=NULL, collide_after=NULL, group_by=NULL,
                                         events=NULL, tr=NULL, time=NULL, time_offset=NULL,
                                         parcel_ids=NULL, unit_data=NULL, value_name="value",
                                         trial_col="trial", pexp_target=0.9, scale_columns=TRUE,
                                         use_legacy=FALSE) {

  #input validation
  checkmate::assert_string(event)
  checkmate::assert_numeric(time_before, upper=0, max.len=1)
  checkmate::assert_numeric(time_after, lower=0, max.len=1)
  checkmate::assert_character(collide_before, null.ok=TRUE, unique=TRUE)
  checkmate::assert_character(collide_after, null.ok=TRUE, unique=TRUE)

  if (!inherits(fmri_obj, "fmri_ts")) {
    group_by_missing <- missing(group_by)
    if (inherits(fmri_obj, "medusa_nifti_ts")) {
      if (is.null(tr) && is.finite(fmri_obj$tr)) tr <- fmri_obj$tr
      if (is.null(unit_data) && !is.null(fmri_obj$voxel_map)) unit_data <- fmri_obj$voxel_map
      if (is.null(parcel_ids) && !is.null(fmri_obj$parcel_ids)) parcel_ids <- fmri_obj$parcel_ids
      if (length(fmri_obj$bold) != 1L) {
        stop("Pass a single-run matrix or single-run 'medusa_nifti_ts' to get_medusa_compression_score().", call. = FALSE)
      }
      mat <- normalize_matrix_orientation(fmri_obj$bold[[1L]], orientation = fmri_obj$orientation)
    } else {
      mat <- normalize_matrix_orientation(fmri_obj, label = "fmri_obj")
    }
    if (is.null(events)) stop("'events' must be supplied when computing compression scores on a matrix.", call. = FALSE)
    if (is.null(unit_data) && !is.null(attr(mat, "voxel_map", exact = TRUE))) {
      unit_data <- attr(mat, "voxel_map", exact = TRUE)
    }
    if (is.null(parcel_ids) && !is.null(attr(mat, "parcel_ids", exact = TRUE))) {
      parcel_ids <- attr(mat, "parcel_ids", exact = TRUE)
    }
    tr_eff <- if (!is.null(tr)) tr else attr(mat, "tr", exact = TRUE)
    udf <- build_unit_metadata_df(ncol(mat), colnames(mat), unit_data = unit_data, parcel_ids = parcel_ids)
    t_vec <- resolve_run_time_vec(time, 1L, "run1", nrow(mat), tr_eff, time_offset, mat)
    if (group_by_missing || is.null(group_by)) {
      has_p <- !is.null(parcel_ids) || !is.null(attr(mat, "parcel_ids", exact = TRUE))
      group_by <- if (has_p) c(trial_col, "parcel") else trial_col
    }
    return(compression_score_matrix_single_run(
      signals = stats::setNames(list(mat), value_name),
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
      group_by = group_by,
      pexp_target = pexp_target,
      scale_columns = scale_columns
    ))
  }

  if (!isTRUE(use_legacy) && !is.null(group_by) && (fmri_obj$vm$trial %in% group_by)) {
    view <- extract_fmri_ts_matrix_view(fmri_obj, group_by = group_by)
    if (!is.null(view)) {
      return(compression_score_matrix_single_run(
        signals = view$signals,
        raw_time = view$raw_time,
        n_time_raw = view$n_time_raw,
        n_units = view$n_units,
        unit_df = view$unit_df,
        events = view$events,
        event = event,
        trial_col = view$trial_col,
        time_before = time_before,
        time_after = time_after,
        collide_before = collide_before,
        collide_after = collide_after,
        group_by = group_by,
        pexp_target = pexp_target,
        scale_columns = scale_columns,
        coerce_factor_trial = TRUE
      ))
    }
  }

  #talign is an aligned fmri_ts object keyed by trial (and other keying variables)
  talign <- event_lock_ts(fmri_obj, event=event, time_before=time_before, time_after=time_after,
    pad_before=0, pad_after=0, collide_before=collide_before, collide_after=collide_after)

  interp_dt <- talign$get_ts(orig_names=TRUE) %>%
    data.table::melt(measure.vars=talign$vm$value, variable.name=".vkey", values.name="value")

  #small wrapper function to cast a given matrix into 
  calculate_compression <- function(dt, keys=NULL, ...) {
    #convert to wide format with multiple signals (e.g., from many voxels) as columns
    ff <- as.formula(paste0("evt_time ~ ", paste(keys, collapse=" + ")))
    thismat <- data.table::dcast(dt, ff, value.var="value")
    thismat[, evt_time := NULL] #omit time column itself from matrix
    thismat <- as.matrix(thismat)
    compress_mts_pca(thismat, ...)
  }

  #find the key variables that are not in the group_by
  other_keys <- setdiff(talign$vm$key, group_by)

  #so, on the columns of the dcast, it's all keys that are not in the group_by... making them simultaneous columns
  comp_calc <- interp_dt[, calculate_compression(.SD, keys=other_keys, pexp_target=pexp_target, scale_columns=scale_columns), by=c(".vkey", group_by)]

  #reshape multiple signals to wide format again
  output_dt <- data.table::dcast(comp_calc, formula = ... ~ .vkey, value.var=grep("a*compress_.*", names(comp_calc), value=TRUE), sep="|")

  #make names more readable: like decon_compress_0.9, not compress_0.9_decon
  setnames(output_dt, names(output_dt), sub("(a*compress_[^|]+)\\|(.+)", "\\2_\\1", names(output_dt), perl=TRUE))

  return(output_dt)

}

#' Interpolate an aligned fmri_ts object onto a consistent time grid with respect to a target event
#'
#' @param a_obj an fmri_ts object that has data aligned to an event of interest
#' @param evt_time column name in \code{a_obj} event data to which time series should be aligned and interpolated
#' @param time_before The earliest time point (in seconds) relative to the event to be output in interpolation
#' @param time_after The latest time point (in seconds) relative to the event to be output in interpolation
#' @param output_resolution The timestep (in seconds) used for interpolation
#' @param group_by a character vector of keying variables used for aggregation of data prior to interpolation
#' @param time_audit If TRUE, additional columns will be added to the output showing how alignment is calculated vis-a-vis event timing
#'
#' @param max_gap Maximum source bracket duration (default acquisition TR).
#'
#' @importFrom data.table dcast melt setnames
#' @importFrom dplyr summarize group_by across
#' @author Michael Hallquist
#' @keywords internal
interpolate_fmri_epochs <- function(a_obj, evt_time="evt_time", time_before=-3, time_after=3,
                                    output_resolution=1.0, group_by=NULL, time_audit=FALSE,
                                    max_gap=NULL, readout_kernel="linear", readout_width=NULL) {
  stopifnot(inherits(a_obj, "fmri_ts"))
  tout <- epoch_output_grid(time_before, time_after, output_resolution)
  if (is.null(max_gap)) max_gap <- a_obj$tr
  checkmate::assert_number(max_gap, lower = .Machine$double.eps, finite = TRUE)
  quad <- epoch_kernel_quadrature(readout_kernel, readout_width, max_gap)
  input <- a_obj$get_ts(orig_names = TRUE)
  for (column in a_obj$vm$value) {
    checkmate::assert_numeric(input[[column]], all.missing = TRUE)
    input[[column]] <- as.numeric(input[[column]])
  }
  interp_dt <- data.table::melt(input,
    measure.vars = a_obj$vm$value, variable.name = ".vkey", value.name = "value")
  checkmate::assert_subset(group_by, names(interp_dt))
  unit_keys <- setdiff(a_obj$vm$key, group_by)
  time_col <- evt_time
  trial_col <- a_obj$vm$trial
  trial_grouped <- trial_col %in% group_by
  audit_columns <- c("epoch_status", grep("^align_", names(interp_dt), value = TRUE),
    if (time_audit) c("evt_onset", "collide_before", "collide_after"))
  audit_columns <- intersect(audit_columns, names(interp_dt))
  if (trial_grouped) {
    audit_dt <- unique(interp_dt[, c(trial_col, audit_columns), with = FALSE])
    if (anyDuplicated(audit_dt[[trial_col]])) stop("Inconsistent epoch annotations within trial")
  }
  # The group worker interpolates each unit before summarizing, using shared
  # matrix brackets on regular grids and a bounded-memory ragged-grid fallback.
  interp_dt <- interp_dt[, interpolate_epoch_group(.SD, unit_keys, time_col, tout,
    max_gap, trial_grouped, readout_kernel = quad$kernel, readout_width = quad$width), by = c(".vkey", group_by)]
  interp_dt[, .result_id := seq_len(.N)]
  source_groups <- c(".vkey", setdiff(group_by, trial_col))
  usage <- interp_dt[, list(.member = as.character(unlist(.usage, use.names = FALSE))),
    by = c(source_groups, ".result_id")]
  interp_dt[, value_source_max_epochs := 0L]
  if (nrow(usage)) {
    usage[, .epoch := sub(":.*$", "", .member)]
    usage[, .source := sub("^[^:]+:", "", .member)]
    counts <- unique(usage[, c(source_groups, ".epoch", ".source"), with = FALSE])[
      , list(.epochs = data.table::uniqueN(.epoch)), by = c(source_groups, ".source")]
    usage <- merge(usage, counts, by = c(source_groups, ".source"), sort = FALSE)
    row_counts <- usage[, list(.epochs = max(.epochs)), by = ".result_id"]
    interp_dt$value_source_max_epochs[row_counts$.result_id] <- row_counts$.epochs
  }
  interp_dt[, c(".usage", ".result_id") := NULL]
  stats_columns <- setdiff(names(interp_dt), c(".vkey", group_by, "evt_time"))
  if (anyDuplicated(interp_dt[, c(".vkey", group_by, "evt_time"), with = FALSE])) {
    stop("Duplicate interpolation keys; refusing to aggregate signal values as counts")
  }
  output_dt <- data.table::dcast(interp_dt, formula = ... ~ .vkey,
    value.var = stats_columns, sep = "|")
  data.table::setnames(output_dt, names(output_dt),
    sub("value_([^|]+)\\|(.+)", "\\2_\\1", names(output_dt), perl = TRUE))
  if (trial_grouped) {
    output_dt <- merge(output_dt, audit_dt, by = trial_col, all.x = TRUE, sort = FALSE)
  }
  if (trial_grouped && nrow(output_dt)) {
    mean_columns <- intersect(paste0(a_obj$vm$value, "_mean"), names(output_dt))
    support <- output_dt[, list(.support_status = {
      usable <- rowSums(is.finite(as.matrix(.SD))) > 0L
      if (!any(usable)) "insufficient_data" else if (all(usable)) "ok" else "partial"
    }), by = group_by, .SDcols = mean_columns]
    output_dt <- merge(output_dt, support, by = group_by, all.x = TRUE, sort = FALSE)
    output_dt[epoch_status == "ok", epoch_status := .support_status]
    output_dt[, .support_status := NULL]
  }
  if (isTRUE(time_audit) && trial_grouped) {
    output_dt[, time := ifelse(is.finite(evt_onset), evt_time + evt_onset, NA_real_)]
  }
  attr(output_dt, "source_clock") <- attr(a_obj, "source_clock")
  output_dt
}


#' Core function for computing multivariate time series compression scores by principal components
#'
#' @param mts multivariate time series structured time x signals
#' @param pexp_target Proportion of variance explained by principal components. This can be a vector,
#'   in which case compression is calculated at different thresholds.
#' @param scale_columns whether to z-score the time series prior to eigendecomposition (recommended)
#'
#' @details This function accepts a time x signals (e.g., voxels) time series matrix. It computes the
#'   singular value decomposition (SVD) and then examines how many eigenvectors are needed to explain at least
#'   \code{pexp_target} proportion of variance.
#'
#' Compression scores are normalized 0 -- 1.0 by the equation: 1 - (n_components / n_timeseries). Thus,
#'   if 6 components explain 91% of the variance in a 12-time series matrix, then the compression score is 0.5.
#'
#' Given that the number of eigenvectors is an integer, a linear approximation to the exact proportion of variance
#'   explained is also calculated. For example, if 3 components explain 84% of variance, while 4 components explain
#'   92%, we 'overshoot' our target of 90% slightly with 4 components. We can then calculated how much compression would
#'   get us to precisely 90% by linear interpolation.
#'
#' @return a list containing compression estimates of the matrix. For each pexp_target value, two values are included,
#'  one representing the compression calculated using integer
#' @importFrom checkmate assert_matrix assert_numeric assert_logical
#' @author Michael Hallquist
#' @export
compress_mts_pca <- function(mts, pexp_target=0.9, scale_columns=TRUE) {
  checkmate::assert_matrix(mts)
  checkmate::assert_numeric(pexp_target, lower=1e-2, upper=1.0, unique=TRUE)
  checkmate::assert_logical(scale_columns, max.len=1L)
  
  #orig <- mts
  mts <- mts[, !apply(mts, 2, function(x) all(is.na(x))), drop=FALSE] #drop columns/time series that are all NA (edge of brain)
  mts <- na.omit(mts) #now drop rows that have NAs (event censoring)

  #screen for lousy time series (no variation)
  sds <- apply(mts, 2, sd)

  bad <- !is.finite(sds) | sds < 1e-3 #very small or nonfinite SD
  mts <- mts[, !bad, drop=FALSE]

  if (nrow(mts) <= 3 || ncol(mts) <= 1) {
    warning("mts matrix cannot be compressed. Dims: ", paste(dim(mts), collapse=", "))
    compress_ret <- list() #setup dummy list
    for (pp in pexp_target) {
      compress_ret[[paste0("compress_", pp)]] <- NA_real_
      compress_ret[[paste0("acompress_", pp)]] <- NA_real_
    }

    compress_ret[["compress_nt"]] <- NA_integer_
  } else {
    #if (isTRUE(scale_columns)) { mts <- scale(mts) } #z-score columns of matrix before compression (standard practice)
    #for some reason, scale is blowing up svd (infinite/missing values -- just use manual z-scoring)
    if (isTRUE(scale_columns)) { 
      mts <- apply(mts, 2, function(x) { (x - mean(x)) / sd(x) })
    }

    dd <- svd(mts)
    pexp <- cumsum(dd$d^2)/sum(dd$d^2) #proportion of variance explained

    compress_ret <- list()
    for (pp in pexp_target) {
      nexceed <- min(which(pexp > pp))

      #linear approximation of rank, allowing non-integer values
      napprox <- approx(y=1:length(dd$d), x=pexp, xout=pp, ties="ordered")$y
      compress_ret[[paste0("compress_", pp)]] <- 1-(nexceed/length(dd$d))
      compress_ret[[paste0("acompress_", pp)]] <- 1-(napprox/length(dd$d))
    }

    compress_ret[["compress_nt"]] <- length(dd$d) #number of timepoints

  }

  return(compress_ret)
}


#' Align a set of deconvolved time series files from \code{voxelwise_deconvolution} to a set of events in a task-related design
#' @param atlas_files A character vector containing filenames of atlases that used for aggregating voxel-level deconvolved time series
#' @param decon_dir The directory containing the deconvolved data from \code{voxelwise_deconvolution}. Should be the same as
#'   \code{out_dir} from \code{voxelwise_deconvolution}.
#' @param trial_df Trial-level data used for event alignment.
#' @param alignments A list of alignment specifications, each of which will be computed by this function.
#'   See \code{get_medusa_interpolated_ts} for help with arguments, and see Details.
#'   Settings include overlap_policy, valid_col, neighbor_events, and max_gap.
#' @param nbins For atlases with continuous values, how many bins should be used to discretize the mask values, leading to aggregation
#'   by bin.
#' @param aggregate_by The column name in the individual deconvolved files used for averaging repeated units (e.g., voxels) into single
#'   event-aligned time series. Most commonly, this is "atlas_value", which will lead to averaging of voxels within each parcel in the mask.
#' @param atlas_subset Optional subset of atlas values to retain.
#' @param overwrite if TRUE, overwrite existing output files
#' @param tr The repetition time of the scan in seconds
#' @param ncpus The number of cores to use for each parallel job
#' @param mem_per_cpu The amount of memory (RAM) allocated for each CPU job. Default is "4g" (4 gigabytes)
#' @param walltime The time requested for each event alignment job. Default: 1:00:00 (1 hour).
#' @param scheduler The job scheduler to be used in alignment. Options are "slurm" (default), "torque", or "local"
#' @param time_audit If TRUE, additional columns will be added to the output showing how alignment is calculated vis-a-vis event timing
#' @param file_manifest Optional data.frame with \code{file}, \code{id},
#'   \code{run}, and optional \code{session} columns. Paths must identify the
#'   deconvolution files. Overrides per-file metadata and legacy filename parsing.
#' @param fail_on_error Stop if any file fails alignment (default TRUE). Set FALSE
#'   explicitly to allow partial output with a warning and a file status manifest.
#' @param input_format Format of deconvolution files to discover inside each atlas's
#'   \code{deconvolved/} directory: \code{"auto"} (default, discovers
#'   \code{_deconvolved.csv.gz}, \code{_deconvolved.rds}, \code{_deconvolved.nii.gz},
#'   and \code{_deconvolved.parquet}, deduplicating stems if multiple formats coexist),
#'   or an explicit format (\code{"csv.gz"}, \code{"rds"}, \code{"nifti"}, \code{"parquet"}).
#' @param out_format Output serialization format for each aligned panel:
#'   \code{"csv.gz"} (default, written via \code{data.table::fwrite}),
#'   \code{"parquet"} (columnar Parquet via \code{arrow::write_parquet}, directly
#'   consumable by \code{\link{mixed_by}} with column projection and filter pushdown),
#'   or \code{"rds"}.
#'
#' @details
#'   This function will create a separate R batch job on the scheduler for each combination of atlas files and alignments. Each of these
#'   jobs can be allocated a number of cores, which will be used to loop over the deconvolved time series files for event alignment. The
#'   number of cores requested for each alignment job is specified using \code{ncpus}, which defaults to 8.
#'
#'   The alignments argument is a list containing information about how to event-align the deconvolved time series. Elements include:
#'
#' \itemize{
#'   \item \code{$evt_col}: (Required) The name of the column in \code{trial_df} containing the event of interest.
#'   \item \code{$time_before}: (Required) The number of seconds before the event of interest to include in the interpolation window.
#'   \item \code{$time_after}: (Required) The number of seconds after the event of interest to include in the interpolation window.
#'   \item \code{$pad_before}: (Optional) The number of seconds before the earliest event to include in the interpolation window
#'      to avoid missing values at the edge of the interpolation window. Default: -1.5
#'   \item \code{$pad_after}: (Optional) The number of seconds after the latest event to include in the interpolation window
#'      to avoid missing values at the edge of the interpolation window. Default: 1.5
#' }
#' 
#' @examples 
#'
#' \dontrun{
#'   atlas_dir <- file.path(
#'     "/proj/mnhallqlab/projects/clock_analysis/fmri",
#'     "ph_da_striatum/masks"
#'   )
#'   atlas_files <- file.path(
#'     atlas_dir,
#'     c(
#'       "bilateral_striatum_tight_7Networks_2.3mm.nii.gz",
#'       "pauli_combined_integermask_2.3mm.nii.gz"
#'     )
#'   )
#'
#'   # Outputs of voxelwise_deconvolution() for these atlases.
#'   decon_dir <- "/proj/mnhallqlab/users/michael/sceptic_decon"
#'   trial_df <- get_trial_data(
#'     repo_directory = "/proj/mnhallqlab/projects/clock_analysis",
#'     dataset = "mmclock_fmri",
#'     groupfixed = TRUE
#'   )
#'
#'   alignments <- list(
#'     clock_long = list(
#'       evt_col = "clock_onset",
#'       time_before = -5,
#'       time_after = 10,
#'       collide_before = "iti_onset", # censor data if we bump into the end of the prior trial
#'       collide_after = "clock_onset" # censor data it we hit the next trial
#'     ),
#'     rt_long = list(
#'       evt_col = "rt_time",
#'       time_before = -4,
#'       time_after = 7,
#'       collide_before = "iti_onset",
#'       collide_after = "clock_onset"
#'     )
#'   )
#'
#'   # run all atlases and alignments as separate slurm jobs (2 hours each, 8 cpus) 
#'   run_decon_alignment(atlas_files, decon_dir, trial_df, alignments, 
#'     overwrite = TRUE, tr = 1.0, ncpus = 8, walltime = "2:00:00", scheduler = "slurm")
#'
#' }
#' @export
run_decon_alignment <- function(atlas_files, decon_dir, trial_df, alignments = list(), nbins = 12, aggregate_by = "atlas_value", 
    atlas_subset=NULL, overwrite = FALSE, tr = NULL, ncpus = 8, mem_per_cpu = "4g", walltime = "1:00:00", scheduler = "slurm", time_audit=FALSE,
    file_manifest = NULL, fail_on_error = TRUE,
    input_format = c("auto", "csv.gz", "rds", "nifti", "parquet"),
    out_format = c("csv.gz", "parquet", "rds")) {
  checkmate::assert_file_exists(atlas_files)
  checkmate::assert_directory_exists(decon_dir)
  checkmate::assert_data_frame(trial_df)
  checkmate::assert_list(alignments, names = "named") # require named list input
  checkmate::assert_integerish(nbins, lower = 1, upper = 1e4, len = 1L)
  checkmate::assert_string(aggregate_by)
  checkmate::assert_numeric(atlas_subset, null.ok = TRUE)
  checkmate::assert_logical(overwrite, len = 1L)
  checkmate::assert_number(tr, lower = 0.01, upper = 1000)
  checkmate::assert_integerish(ncpus, lower = 1, upper = 1e3, len = 1L)
  walltime <- validate_dhms(walltime) # validate and correct ambiguities
  checkmate::assert_string(scheduler)
  checkmate::assert_subset(scheduler, c("local", "slurm", "torque"))
  checkmate::assert_flag(fail_on_error)
  input_format <- match.arg(input_format)
  out_format <- match.arg(out_format)
  if (identical(input_format, "parquet") || identical(out_format, "parquet")) {
    require_arrow_pkg("run_decon_alignment() with Parquet files")
  }
  out_ext <- switch(
    out_format,
    "csv.gz" = ".csv.gz",
    "parquet" = ".parquet",
    "rds" = ".rds"
  )

  alignments <- validate_alignments(alignments, tr, trial_df)

  # loop over atlases and alignments
  for (af in atlas_files) {
    if (!checkmate::test_file_exists(af)) {
      message("Cannot find input atlas file: {af}. Skipping alignment for this atlas.")
      next
    }

    aname <- file_sans_ext(basename(af)) # atlas name

    mask <- RNifti::readNifti(af)
    mi <- which(mask > 0, arr.ind = TRUE)
    maskvals <- unique(mask[mi])

    if (isTRUE(checkmate::test_integerish(maskvals))) {
      continuous <- FALSE
      message("Atlas: ", aname, " has integer values. Using unique mask values for aggregation in event-aligning.")
      atlas_cuts <- NULL
    } else {
      continuous <- TRUE
      message("Atlas: ", aname, "has continuous values. Dividing into ", nbins, "bins")
      atlas_cuts <- seq(min(mask[mi]) - 1e-5, max(mask[mi]) + 1e-5, length.out = nbins + 1)
    }

    adir <- file.path(decon_dir, aname, "deconvolved")
    if (!dir.exists(adir)) {
      message(glue("Cannot find expected deconvolved output directory: {adir}. Have you run voxelwise_deconvolution yet?"))
      next
    }

    d_files <- discover_decon_files(adir, input_format = input_format)
    if (length(d_files) == 0L) {
      message(glue("No deconvolved files found in directory: {adir}"))
      next
    }

    atlas_manifest <- validate_decon_manifest(file_manifest, d_files)

    # loop over event alignments within this atlas
    for (ee in names(alignments)) {
      out_file <- file.path(decon_dir, aname, glue("{aname}_{ee}_decon_aligned{out_ext}"))
      this_alignment <- alignments[[ee]]
      this_alignment$aname <- aname # populate atlas name so that it can be added in downstream alignment function

      if (file.exists(out_file) && isFALSE(overwrite)) {
        message("Output file already exists: ", out_file)
        next
      }

      # mask <- list(continuous = continuous, atlas_cuts = atlas_cuts, nifti = mask)

      d_batch <- create_decon_batch_job(
        job_name = glue("evtalign_{aname}_{ee}"), n_cpus = ncpus, mem_per_cpu = mem_per_cpu,
        wall_time = walltime, scheduler = scheduler,
        # pass relevant vars to the batch
        input_objects = named_list(d_files, trial_df, this_alignment, aggregate_by, tr, atlas_cuts, atlas_subset, out_file, ncpus, time_audit, atlas_manifest, fail_on_error),
        r_packages = "medusa",
        r_code = c(
          "evt_align_decon_files(d_files, trial_df, this_alignment, aggregate_by, tr, atlas_cuts, atlas_subset, out_file, ncpus, time_audit, file_manifest = atlas_manifest, fail_on_error = fail_on_error)"
        )
      )

      d_batch$submit()
    }
  }
}

#' Align a set of deconvolved time series files to an event of interest function.
#' @details
#'   This is intended to be used internally by \code{run_decon_alignment}, which accepts a set of mask/atlas files
#'   and alignments, then processes these in parallel. Long \code{.csv}/\code{.csv.gz} and \code{.parquet} tables
#'   sharing a run clock are pivoted directly into \code{time x voxels} matrices in C++ without intermediate
#'   \code{fmri_ts} allocation (automatically falling back to \code{fmri_ts} for ragged per-unit clocks); compact
#'   \code{.rds} matrix bundles (\code{list(decon = Y, coords = ..., time = ...)}) and 4D \code{.nii}/\code{.nii.gz}
#'   deconvolved images (paired with \code{.coords.rds} sidecars) are also supported directly.
#' @param d_files A vector of deconvolved \code{.csv}/\code{.csv.gz}, \code{.rds}, \code{.nii}/\code{.nii.gz}, or
#'   \code{.parquet} files (created by \code{voxelwise_deconvolution}).
#' @param trial_df The trial-level data.frame containing id and run for each subject represented in \code{d_files}.
#' @param alignment A list containing alignment details passed to \code{\link{align_medusa_ts}} /
#'   \code{\link{get_medusa_interpolated_ts}}, including \code{evt_col}, \code{time_before}, \code{time_after},
#'   \code{pad_before}, \code{pad_after}, \code{output_resolution}, \code{collide_before}, \code{collide_after},
#'   \code{overlap_policy}, \code{valid_col}, \code{neighbor_events}, \code{max_gap}, \code{readout_kernel},
#'   and \code{readout_width}.
#' @param aggregate_by The column name in the individual deconvolved files used for averaging repeated units (e.g., voxels) into single
#'   event-aligned time series. Most commonly, this is "atlas_value", which will lead to averaging of voxels within each parcel in the mask.
#' @param tr The repetition time of the sequence in seconds.
#' @param atlas_cuts For a continuous-valued atlas (e.g., containing a gradient of interest), a vector of cut points for binning values.
#'   These cuts are applied to the \code{aggregate_by} column, commonly "atlas_value"
#' @param atlas_subset An optional numeric vector containing values of \code{aggregate_by} 
#' @param out_file Optional output file path for the event-aligned panel (\code{.csv}, \code{.csv.gz}, \code{.parquet},
#'   or \code{.rds}). If NULL, nothing is written to disk (the event-aligned data are always returned).
#' @param ncpus Number of local worker processes to use.
#' @param time_audit If TRUE, additional columns will be added to the output showing how alignment is calculated vis-a-vis event timing
#' @param file_manifest Optional data.frame with \code{file}, \code{id},
#'   \code{run}, and optional \code{session} columns. Overrides per-file metadata
#'   and legacy filename parsing; ambiguous identities are errors.
#' @param fail_on_error Stop if any file fails (default TRUE). FALSE permits
#'   partial output with a warning. A status row is retained for every file.
#' @return Aligned data with a \code{file_status} data.frame attribute. If
#'   \code{out_file} is supplied, the manifest is saved to
#'   \code{paste0(out_file, ".files.csv")}, including on failure. With
#'   \code{fail_on_error = TRUE}, a \code{decon_alignment_error} condition includes
#'   the manifest and aligned output is not written if any file fails. Successful
#'   results also carry a \code{source_clock} table, saved to
#'   \code{paste0(out_file, ".sources.csv")} when output is requested. The
#'   \code{alignment_file} column scopes acquisition-time IDs to their input file.
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach registerDoSEQ
#' @importFrom dplyr select filter
#' @importFrom tidyselect all_of
#' @importFrom rlang sym
#' @export
evt_align_decon_files <- function(d_files, trial_df, alignment = list(), aggregate_by="atlas_value", 
    tr = NULL, atlas_cuts = NULL, atlas_subset = NULL, out_file = NULL, ncpus = 8, time_audit=FALSE, file_manifest = NULL, fail_on_error = TRUE) {
  checkmate::assert_character(d_files, min.len = 1L, unique = TRUE)
  for (df_path in d_files) {
    if (dir.exists(df_path)) {
      if (!is_parquet_path(df_path)) {
        stop("Directory does not contain Parquet files: ", df_path, call. = FALSE)
      }
    } else {
      checkmate::assert_file_exists(df_path)
    }
  }
  checkmate::assert_data_frame(trial_df)
  checkmate::assert_list(alignment)
  checkmate::assert_number(tr, lower = 0.01, upper = 1000)
  checkmate::assert_numeric(atlas_cuts, null.ok = TRUE)
  checkmate::assert_numeric(atlas_subset, null.ok = TRUE)
  checkmate::assert_string(out_file, null.ok = TRUE)
  if (!is.null(out_file) && grepl("\\.(parquet|pq)$", out_file, ignore.case = TRUE)) {
    require_arrow_pkg("evt_align_decon_files(out_file = '*.parquet')")
  }
  checkmate::assert_integerish(ncpus, lower=1L, len=1L)

  checkmate::assert_flag(fail_on_error)
  file_manifest <- validate_decon_manifest(file_manifest, d_files)
  alignment <- validate_alignments(alignment, tr, trial_df)

  align_one_file <- function(fname) {
    identity <- NULL
    tryCatch({
      is_rds <- grepl("\\.rds$", fname, ignore.case = TRUE)
      is_nifti <- !is_rds && grepl("\\.nii(\\.gz)?$", fname, ignore.case = TRUE)
      is_pq <- !is_rds && !is_nifti && is_parquet_path(fname)
      rds_bundle <- NULL
      if (is_rds) {
        obj_in <- readRDS(fname)
        if (is.data.frame(obj_in)) {
          d <- as.data.frame(obj_in)
        } else if (is.matrix(obj_in)) {
          rds_bundle <- list(decon = obj_in, coords = attr(obj_in, "voxel_map", exact = TRUE),
                             time = attr(obj_in, "time", exact = TRUE),
                             valid = attr(obj_in, "observation_mask", exact = TRUE),
                             meta = data.frame())
        } else if (is.list(obj_in) && "decon" %in% names(obj_in) && is.matrix(obj_in$decon)) {
          rds_bundle <- list(decon = obj_in$decon, coords = obj_in$coords,
                             time = obj_in$time, valid = obj_in$valid,
                             meta = if (is.data.frame(obj_in$metadata)) obj_in$metadata else data.frame())
        } else {
          stop("Unsupported RDS object in ", fname)
        }
      } else if (is_nifti) {
        coords_path <- paste0(fname, ".coords.rds")
        if (!file.exists(coords_path)) {
          stop("Missing coordinate sidecar required to align deconvolved NIfTI: ", coords_path)
        }
        cmeta <- readRDS(coords_path)
        if (!is.list(cmeta) || is.null(cmeta$voxel_indices) || is.null(cmeta$coords)) {
          stop("Invalid coordinate sidecar in ", coords_path)
        }
        mat_nii <- extract_nifti_voxel_ts(fname, cmeta$voxel_indices, TRUE)
        if (any(!is.finite(mat_nii))) {
          mat_nii[!is.finite(mat_nii)] <- NA_real_
        }
        rds_bundle <- list(
          decon = mat_nii,
          coords = cmeta$coords,
          time = cmeta$time,
          valid = cmeta$valid,
          meta = if (is.data.frame(cmeta$metadata)) cmeta$metadata else data.frame()
        )
      } else if (is_pq) {
        avail_cols <- arrow_dataset_names(fname)
        drop_cand <- setdiff(
          c("volume", "x", "y", "z", "atlas_name", "decon_status"),
          c(aggregate_by, alignment$valid_col)
        )
        cols_read <- setdiff(avail_cols, drop_cand)
        d <- medusa_read_table(fname, col_select = cols_read, as_data_table = FALSE)
      } else {
        if (grepl("\\.gz$", fname, ignore.case = TRUE) && requireNamespace("R.utils", quietly = TRUE)) {
          loadNamespace("R.utils")
        }
        avail_cols <- peek_delim_header(fname)
        drop_cand <- setdiff(
          c("volume", "x", "y", "z", "atlas_name", "decon_status"),
          c(aggregate_by, alignment$valid_col)
        )
        cols_read <- setdiff(avail_cols, drop_cand)
        if (length(cols_read) > 0L && length(cols_read) < length(avail_cols)) {
          d <- data.table::fread(fname, select = cols_read, data.table = FALSE)
        } else {
          d <- data.table::fread(fname, data.table = FALSE)
        }
      }

      rk_eff <- if (is.null(alignment$readout_kernel)) "linear" else alignment$readout_kernel
      grp_eff <- unique(c(aggregate_by, "trial"))

      if (!is.null(rds_bundle)) {
        mat <- rds_bundle$decon
        coords <- rds_bundle$coords
        if (is.null(coords) || !(aggregate_by %in% names(coords))) {
          stop("Missing required signal columns: ", aggregate_by)
        }
        identity <- resolve_decon_identity(fname, rds_bundle$meta, file_manifest)
        subj_df <- match_decon_trials(identity, as.data.frame(trial_df))
        identity <- as.list(subj_df[1L, intersect(c("id", "run", "session"), names(subj_df)), drop = FALSE])
        if (length(mat) == 0L || !any(is.finite(mat))) stop("No finite deconvolution estimates in file")
        agg_vec <- if (!is.null(atlas_cuts)) {
          as.character(cut(coords[[aggregate_by]], atlas_cuts))
        } else {
          as.numeric(coords[[aggregate_by]])
        }
        keep_cols <- if (!is.null(atlas_subset)) which(agg_vec %in% atlas_subset) else seq_len(ncol(mat))
        if (length(keep_cols) == 0L) stop("No signals remain after atlas selection")
        mat <- mat[, keep_cols, drop = FALSE]
        udf <- data.frame(vnum = seq_along(keep_cols), stringsAsFactors = FALSE)
        udf[[aggregate_by]] <- agg_vec[keep_cols]
        t_vec <- if (!is.null(rds_bundle$time)) as.double(rds_bundle$time) else (seq_len(nrow(mat)) - 1.0) * tr
        v_mask <- if (!is.null(alignment$valid_col)) rds_bundle$valid else NULL
        subj_align <- align_matrix_single_run(
          signals = list(decon = mat),
          raw_time = t_vec,
          n_time_raw = nrow(mat),
          n_units = ncol(mat),
          unit_df = udf,
          events = subj_df,
          event = alignment$evt_col,
          time_before = alignment$time_before,
          time_after = alignment$time_after,
          collide_before = alignment$collide_before,
          collide_after = alignment$collide_after,
          pad_before = alignment$pad_before,
          pad_after = alignment$pad_after,
          output_resolution = alignment$output_resolution,
          tr = tr,
          group_by = grp_eff,
          time_audit = time_audit,
          overlap_policy = alignment$overlap_policy,
          raw_valid = v_mask,
          neighbor_events = alignment$neighbor_events,
          max_gap = alignment$max_gap,
          readout_kernel = rk_eff,
          readout_width = alignment$readout_width,
          coerce_factor_trial = TRUE
        )
      } else {
        required <- unique(c("decon", "time", "vnum", aggregate_by))
        if (!all(required %in% names(d))) {
          stop("Missing required signal columns: ", paste(setdiff(required, names(d)), collapse = ", "))
        }
        identity <- resolve_decon_identity(fname, d, file_manifest)
        subj_df <- match_decon_trials(identity, as.data.frame(trial_df))
        # Use the matched trial identities, retaining their original types.
        identity <- as.list(subj_df[1L, intersect(c("id", "run", "session"), names(subj_df)), drop = FALSE])
        if (nrow(d) == 0L || !any(is.finite(d$decon))) stop("No finite deconvolution estimates in file")
        if (!is.null(atlas_cuts)) {
          d[[aggregate_by]] <- as.character(cut(d[[aggregate_by]], atlas_cuts))
        } else {
          d[[aggregate_by]] <- as.numeric(d[[aggregate_by]])
        }
        if (!is.null(atlas_subset)) d <- d[d[[aggregate_by]] %in% atlas_subset, , drop = FALSE]
        if (nrow(d) == 0L) stop("No signals remain after atlas selection")

        d_keep_names <- setdiff(names(d), c("atlas_name", "x", "y", "z"))
        if (any(grepl("^(align_|\\.align_)", d_keep_names))) {
          stop("Input contains reserved alignment columns")
        }
        valid_vec <- NULL
        if (!is.null(alignment$valid_col)) {
          checkmate::assert_subset(alignment$valid_col, d_keep_names)
          valid_vec <- d[[alignment$valid_col]]
          if (!is.logical(valid_vec)) stop("valid_col must name a logical column")
        }

        grp_vals <- sort(unique(d[[aggregate_by]]), na.last = TRUE)
        grp_idx <- match(d[[aggregate_by]], grp_vals) - 1L
        vnum_num <- if (is.numeric(d$vnum)) as.double(d$vnum) else as.double(match(d$vnum, unique(d$vnum)))
        checkmate::assert_numeric(d$time, all.missing = TRUE)
        checkmate::assert_numeric(d$decon, all.missing = TRUE)

        piv <- fast_pivot_decon_cpp(
          as.double(d$time), vnum_num, as.integer(grp_idx), as.double(d$decon), valid_vec
        )

        if (isTRUE(piv$rectangular)) {
          n_u_piv <- ncol(piv$y)
          udf <- data.frame(vnum = seq_len(n_u_piv), stringsAsFactors = FALSE)
          udf[[aggregate_by]] <- grp_vals[piv$unit_group + 1L]
          if (length(piv$clock) == 0L) {
            y_in <- matrix(NA_real_, nrow = 1L, ncol = n_u_piv)
            t_in <- NA_real_
            nt_in <- 1L
            v_in <- NULL
          } else {
            y_in <- piv$y
            t_in <- piv$clock
            nt_in <- length(piv$clock)
            v_in <- if (length(piv$censored) > 0L) !piv$censored else NULL
          }
          subj_align <- align_matrix_single_run(
            signals = list(decon = y_in),
            raw_time = t_in,
            n_time_raw = nt_in,
            n_units = n_u_piv,
            unit_df = udf,
            events = subj_df,
            event = alignment$evt_col,
            time_before = alignment$time_before,
            time_after = alignment$time_after,
            collide_before = alignment$collide_before,
            collide_after = alignment$collide_after,
            pad_before = alignment$pad_before,
            pad_after = alignment$pad_after,
            output_resolution = alignment$output_resolution,
            tr = tr,
            group_by = grp_eff,
            time_audit = time_audit,
            overlap_policy = alignment$overlap_policy,
            raw_valid = v_in,
            neighbor_events = alignment$neighbor_events,
            max_gap = alignment$max_gap,
            readout_kernel = rk_eff,
            readout_width = alignment$readout_width,
            coerce_factor_trial = TRUE
          )
        } else {
          d <- d[d_keep_names]
          tsobj <- fmri_ts$new(ts_data = d, event_data = subj_df, tr = tr,
            vm = list(value = "decon", key = unique(c("vnum", aggregate_by))))
          subj_align <- get_medusa_interpolated_ts(tsobj, event = alignment$evt_col,
            time_before = alignment$time_before, time_after = alignment$time_after,
            collide_before = alignment$collide_before, collide_after = alignment$collide_after,
            pad_before = alignment$pad_before, pad_after = alignment$pad_after,
            output_resolution = alignment$output_resolution,
            group_by = grp_eff, time_audit = time_audit,
            overlap_policy = alignment$overlap_policy, valid_col = alignment$valid_col,
            neighbor_events = alignment$neighbor_events, max_gap = alignment$max_gap,
            readout_kernel = rk_eff, readout_width = alignment$readout_width,
            use_legacy = TRUE)
        }
      }

      clock <- attr(subj_align, "source_clock")
      file_key <- normalizePath(fname, mustWork = TRUE)
      clock$alignment_file <- rep(file_key, nrow(clock))
      data.table::set(subj_align, j = "alignment_file", value = file_key)
      for (column in names(identity)) {
        data.table::set(subj_align, j = column, value = identity[[column]])
        clock[[column]] <- rep(identity[[column]], nrow(clock))
      }
      list(ts = subj_align, source_clock = clock, identity = identity, status = "ok", message = "")
    }, error = function(err) list(ts = NULL, identity = identity,
      status = "error", message = conditionMessage(err)))
  }

  n_workers <- min(as.integer(ncpus), length(d_files))
  if (n_workers <= 1L) {
    elist <- lapply(d_files, align_one_file)
  } else if (identical(.Platform$OS.type, "unix") && !nzchar(Sys.getenv("RSTUDIO"))) {
    elist <- parallel::mclapply(d_files, align_one_file, mc.cores = n_workers, mc.preschedule = TRUE)
  } else {
    clusterobj <- parallel::makeCluster(n_workers)
    on.exit(try(parallel::stopCluster(clusterobj), silent = TRUE), add = TRUE)
    lib_paths <- .libPaths()
    parallel::clusterCall(clusterobj, function(lp) .libPaths(lp), lib_paths)
    elist <- parallel::parLapply(clusterobj, d_files, align_one_file)
  }

  # Even an unexpected worker-level error must have a manifest row.
  elist <- lapply(elist, function(x) {
    if (inherits(x, "try-error")) {
      cond <- attr(x, "condition")
      list(ts = NULL, identity = NULL, status = "error",
           message = if (!is.null(cond)) conditionMessage(cond) else as.character(x))
    } else if (inherits(x, "error")) {
      list(ts = NULL, identity = NULL, status = "error", message = conditionMessage(x))
    } else {
      x
    }
  })
  file_status <- dplyr::bind_rows(lapply(seq_along(d_files), function(i) {
    entry <- elist[[i]]
    identity <- entry$identity
    data.frame(file = d_files[i], id = if (is.null(identity$id)) NA_character_ else as.character(identity$id),
      run = if (is.null(identity$run)) NA_character_ else as.character(identity$run),
      session = if (is.null(identity$session)) NA_character_ else as.character(identity$session),
      status = entry$status, message = entry$message,
      n_rows = if (is.null(entry$ts)) 0L else nrow(entry$ts))
  }))
  if (!is.null(out_file)) medusa_write_table(file_status, paste0(out_file, ".files.csv"))
  failed <- file_status$status != "ok"
  if (any(failed)) {
    message <- paste0(sum(failed), " deconvolution file(s) failed alignment: ",
      paste(paste0(basename(file_status$file[failed]), ": ", file_status$message[failed]), collapse = "; "))
    if (fail_on_error) {
      stop(structure(list(message = message, call = NULL, file_status = file_status),
        class = c("decon_alignment_error", "error", "condition")))
    }
    warning(message, call. = FALSE)
  }
  all_e <- data.table::as.data.table(dplyr::bind_rows(lapply(elist, "[[", "ts")))
  if (nrow(all_e) && !is.null(alignment$aname)) all_e[, atlas := alignment$aname]
  attr(all_e, "file_status") <- file_status
  source_clock <- dplyr::bind_rows(lapply(elist, "[[", "source_clock"))
  if (!ncol(source_clock)) source_clock <- data.frame(source_id = integer(), time = numeric(), alignment_file = character())
  attr(all_e, "source_clock") <- source_clock
  if (!is.null(out_file)) {
    message("Writing output: ", out_file)
    medusa_write_table(all_e, out_file)
    medusa_write_table(source_clock, paste0(out_file, ".sources.csv"))
  }

  return(all_e)
}

discover_decon_files <- function(
  adir,
  input_format = c("auto", "csv.gz", "rds", "nifti", "parquet")
) {
  input_format <- match.arg(input_format)
  pat <- switch(
    input_format,
    "auto" = "_deconvolved\\.(rds|parquet|pq|csv\\.gz|nii\\.gz|nii)$",
    "csv.gz" = "_deconvolved\\.csv\\.gz$",
    "rds" = "_deconvolved\\.rds$",
    "nifti" = "_deconvolved\\.nii(\\.gz)?$",
    "parquet" = "_deconvolved\\.(parquet|pq)$"
  )
  files <- list.files(
    path = adir,
    pattern = pat,
    full.names = TRUE,
    all.files = TRUE,
    ignore.case = TRUE
  )
  if (length(files) <= 1L || !identical(input_format, "auto")) {
    return(files)
  }
  bnames <- basename(files)
  stems <- sub(
    "_deconvolved\\.(rds|parquet|pq|csv\\.gz|nii\\.gz|nii)$",
    "",
    bnames,
    ignore.case = TRUE
  )
  if (!anyDuplicated(stems)) {
    return(files)
  }
  ext_rank <- rep(5L, length(files))
  ext_rank[grepl("_deconvolved\\.rds$", bnames, ignore.case = TRUE)] <- 1L
  ext_rank[grepl("_deconvolved\\.(parquet|pq)$", bnames, ignore.case = TRUE)] <- 2L
  ext_rank[grepl("_deconvolved\\.csv\\.gz$", bnames, ignore.case = TRUE)] <- 3L
  ext_rank[grepl("_deconvolved\\.nii(\\.gz)?$", bnames, ignore.case = TRUE)] <- 4L
  ord <- order(stems, ext_rank)
  keep <- ord[!duplicated(stems[ord])]
  keep <- sort(keep)
  message(
    "Multiple deconvolution file formats found for the same run stem in ", adir,
    "; selecting one file per stem (priority: rds > parquet > csv.gz > nifti)."
  )
  files[keep]
}


validate_alignments <- function(alignments, tr, trial_df) {
  checkmate::assert_list(alignments, min.len = 1L)
  checkmate::assert_data_frame(trial_df)

  # if a single alignment is passed in, wrap it as a single-element list
  single_align <- FALSE
  if (!is.list(alignments[[1]])) {
    alignments <- list(a=alignments)
    single_align <- TRUE
  }
  

  # validate padding and alignment settings, input defaults
  for (aa in seq_along(alignments)) {
    if (!"pad_before" %in% names(alignments[[aa]])) {
      message(glue("For alignment {names(alignments)[aa]}, defaulting pad_before to -1.5"))
      alignments[[aa]]$pad_before <- -1.5
    } else {
      checkmate::assert_number(alignments[[aa]]$pad_before, upper = 0)
    }

    if (!"pad_after" %in% names(alignments[[aa]])) {
      message(glue("For alignment {names(alignments)[aa]}, defaulting pad_after to 1.5"))
      alignments[[aa]]$pad_after <- 1.5
    } else {
      checkmate::assert_number(alignments[[aa]]$pad_after, lower = 0)
    }

    if (!"output_resolution" %in% names(alignments[[aa]])) {
      message(glue("For alignment {names(alignments)[aa]}, defaulting output_resolution (sampling rate) of {tr}"))
      alignments[[aa]]$output_resolution <- tr
    } else {
      checkmate::assert_number(alignments[[aa]]$output_resolution, lower = 0.01, upper = 1e2)
    }

    has_ingredients <- c("pad_before", "pad_after", "evt_col", "time_before", "time_after", "output_resolution") %in% names(alignments[[aa]])
    if (!all(has_ingredients)) {
      stop(glue("For alignment {names(alignments)[aa]}, missing required fields: {paste(names(alignments[[aa]]), collapse=', ')}"))
    }

    checkmate::assert_string(alignments[[aa]]$evt_col)
    alignments[[aa]]$overlap_policy <- resolve_overlap_policy(alignments[[aa]]$overlap_policy,
      alignments[[aa]]$collide_before, alignments[[aa]]$collide_after)
    checkmate::assert_string(alignments[[aa]]$valid_col, null.ok = TRUE)
    checkmate::assert_character(alignments[[aa]]$neighbor_events, min.len = 1L,
      unique = TRUE, null.ok = TRUE)
    checkmate::assert_subset(alignments[[aa]]$neighbor_events, names(trial_df))
    checkmate::assert_number(alignments[[aa]]$max_gap, lower = .Machine$double.eps,
      finite = TRUE, null.ok = TRUE)
    if (!is.null(alignments[[aa]]$readout_kernel)) {
      alignments[[aa]]$readout_kernel <- match.arg(alignments[[aa]]$readout_kernel,
        c("linear", "boxcar", "triangle", "gaussian"))
    }
    checkmate::assert_number(alignments[[aa]]$readout_width, lower = .Machine$double.eps,
      finite = TRUE, null.ok = TRUE)
    checkmate::assert_number(alignments[[aa]]$time_before, upper = 0, finite = TRUE)
    checkmate::assert_number(alignments[[aa]]$time_after, lower = 0, finite = TRUE)
    checkmate::assert_subset(c(alignments[[aa]]$evt_col, alignments[[aa]]$collide_before,
      alignments[[aa]]$collide_after), names(trial_df))
  }

  if (single_align) alignments <- alignments[[1L]]
  return(alignments)

}
