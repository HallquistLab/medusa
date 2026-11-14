make_ergonomics_nifti_fixture <- function(tmp_dir, n_time = 64L, tr = 1.25, seed = 2026L) {
  set.seed(seed)
  dim3 <- c(5L, 4L, 3L)
  n_spatial <- prod(dim3)

  affine <- matrix(
    c(
      -2.0, 0.0, 0.0, 40.0,
       0.0, 2.0, 0.0, -50.0,
       0.0, 0.0, 2.5, -20.0,
       0.0, 0.0, 0.0, 1.0
    ),
    nrow = 4L, byrow = TRUE
  )
  attr(affine, "code") <- 4L

  write_nii <- function(arr, path, img_tr = NULL, datatype = "float") {
    nii <- RNifti::asNifti(arr)
    RNifti::pixdim(nii) <- if (length(dim(arr)) == 4L) {
      c(2, 2, 2.5, if (is.null(img_tr)) tr else img_tr)
    } else {
      c(2, 2, 2.5)
    }
    RNifti::sform(nii) <- affine
    RNifti::writeNifti(nii, path, datatype = datatype)
    path
  }

  mask_arr <- array(0L, dim = dim3)
  mask_arr[2:4, 1:3, 1:2] <- 1L

  atlas_arr <- array(0L, dim = dim3)
  atlas_arr[2:3, 1:2, 1] <- 10L # 4 voxels in parcel 10
  atlas_arr[4, 1:3, 2] <- 20L   # 3 voxels in parcel 20
  vox_lin <- which(mask_arr > 0L & atlas_arr > 0L)
  n_vox <- length(vox_lin)

  # Build combined trial_df with interleaved cue (rpe = NA) and feedback rows across 2 runs
  make_run_events <- function(run_num, offset_shift = 0) {
    n_trials <- 6L
    cue_onsets <- seq(4, by = 9.5, length.out = n_trials) + offset_shift
    fb_onsets <- cue_onsets + 2.5
    rpe <- round(seq(-1.5, 1.5, length.out = n_trials) + rnorm(n_trials, sd = 0.15), 3)
    rt <- round(runif(n_trials, 0.45, 1.15), 3)
    is_rew <- rpe > 0
    valence <- factor(ifelse(rpe >= 0, "pos", "neg"), levels = c("neg", "pos"))
    rows <- vector("list", 2L * n_trials)
    for (i in seq_len(n_trials)) {
      rows[[2L * i - 1L]] <- data.frame(
        id = "sub-01",
        session = "ses-01",
        run = run_num,
        trial = i,
        condition = "cue",
        onset = cue_onsets[i],
        duration = 0,
        rpe = NA_real_,
        abs_rpe = NA_real_,
        rt = NA_real_,
        is_rewarded = NA,
        valence = factor(NA, levels = c("neg", "pos")),
        stringsAsFactors = FALSE
      )
      rows[[2L * i]] <- data.frame(
        id = "sub-01",
        session = "ses-01",
        run = run_num,
        trial = i,
        condition = "feedback",
        onset = fb_onsets[i],
        duration = 0,
        rpe = rpe[i],
        abs_rpe = abs(rpe[i]),
        rt = rt[i],
        is_rewarded = is_rew[i],
        valence = valence[i],
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }

  trial_df <- rbind(make_run_events(1L, 0.0), make_run_events(2L, 0.4))
  rownames(trial_df) <- NULL

  # Generate task-driven 4D BOLD arrays so HRF fitting & epoch regression see realistic signal
  k_lib <- hrf_library(tr, peak_delays = c(5, 6.5), dispersions = 1, undershoot_ratios = 6)
  sim_run_array <- function(ev_run, seed_offset) {
    set.seed(seed + seed_offset)
    arr <- array(rnorm(n_spatial * n_time, mean = 100, sd = 0.4), dim = c(dim3, n_time))
    t_grid <- (seq_len(n_time) - 1L) * tr
    for (j in seq_along(vox_lin)) {
      parcel_j <- atlas_arr[vox_lin[j]]
      h <- k_lib[, if (parcel_j == 10L) 1L else 2L]
      sig <- numeric(n_time)
      for (e in seq_len(nrow(ev_run))) {
        amp <- if (ev_run$condition[e] == "cue") {
          0.8
        } else {
          1.2 + 0.9 * ev_run$rpe[e]
        }
        dt <- t_grid - ev_run$onset[e]
        active <- which(dt >= 0 & dt <= (length(h) - 1L) * tr)
        if (length(active) > 0L) {
          sig[active] <- sig[active] + amp * stats::approx(
            x = (seq_along(h) - 1L) * tr, y = h, xout = dt[active], rule = 2
          )$y
        }
      }
      ijk <- arrayInd(vox_lin[j], dim3)
      arr[ijk[1L], ijk[2L], ijk[3L], ] <- 100 + sig + stats::arima.sim(
        model = list(ar = 0.25), n = n_time, sd = 0.25
      )
    }
    arr
  }

  arr1 <- sim_run_array(trial_df[trial_df$run == 1L, ], 11L)
  arr2 <- sim_run_array(trial_df[trial_df$run == 2L, ], 22L)

  f_run1 <- write_nii(arr1, file.path(tmp_dir, "sub-01_run-01_bold.nii.gz"), img_tr = tr)
  f_run2 <- write_nii(arr2, file.path(tmp_dir, "sub-01_run-02_bold.nii.gz"), img_tr = tr)
  f_mask <- write_nii(mask_arr, file.path(tmp_dir, "mask.nii.gz"), datatype = "int16")
  f_atlas <- write_nii(atlas_arr, file.path(tmp_dir, "atlas.nii.gz"), datatype = "int16")

  list(
    f_runs = c("run-01" = f_run1, "run-02" = f_run2),
    f_mask = f_mask,
    f_atlas = f_atlas,
    trial_df = trial_df,
    tr = tr,
    affine = affine,
    dim3 = dim3,
    n_vox = n_vox
  )
}

test_that("Proposal 3.5: medusa_nifti_ts flows directly through HRF estimation, deconvolution, NIfTI export, backprojection, and alignment", {
  tmp_dir <- tempfile("ergo_35_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  fx <- make_ergonomics_nifti_fixture(tmp_dir)
  ts <- read_nifti_ts(fx$f_runs, mask = fx$f_mask, atlas = fx$f_atlas)

  # Direct ergonomic call: ts + combined data.frame events (integer run 1/2 -> BIDS "run-01"/"run-02")
  hrf_direct <- estimate_crossrun_hrf(
    ts,
    events = fx$trial_df,
    pooling = "parcel",
    min_parcel_voxels = 2L
  )
  expect_s3_class(hrf_direct, "crossrun_hrf")
  expect_output(print(hrf_direct), "<crossrun_hrf> 2 runs \\| 7 voxels")
  expect_equal(hrf_direct$settings$tr, fx$tr, tolerance = 1e-6)
  expect_equal(
    hrf_direct$specification$spatial_groups,
    stats::setNames(as.character(ts$parcel_ids), names(ts$parcel_ids))
  )
  expect_null(hrf_direct$specification$spatial_coords)
  expect_equal(hrf_direct$voxel_map, ts$voxel_map)
  expect_equal(hrf_direct$geometry$dim, fx$dim3)

  # Manual unpacked reference call must produce bit-for-bit identical numerical results
  manual_events <- list(
    "run-01" = fx$trial_df[fx$trial_df$run == 1L, ],
    "run-02" = fx$trial_df[fx$trial_df$run == 2L, ]
  )
  hrf_manual <- estimate_crossrun_hrf(
    bold = ts$bold,
    events = manual_events,
    tr = ts$tr,
    kernel_grid = hrf_library(ts$tr, preset = "medusa"),
    spatial_groups = ts$parcel_ids,
    pooling = "parcel",
    min_parcel_voxels = 2L
  )
  expect_equal(hrf_direct$kernels, hrf_manual$kernels)
  expect_equal(hrf_direct$selections, hrf_manual$selections)
  expect_equal(hrf_direct$noise, hrf_manual$noise)

  # Deconvolution directly on `ts` propagates geometry & voxel_map and matches manual BOLD list
  decon_direct <- deconvolve_crossrun_hrf(
    ts, hrf_direct, lambda = 0.2, scaling = "noise_sd", penalty = "sobolev"
  )
  decon_manual <- deconvolve_crossrun_hrf(
    ts$bold, hrf_manual, lambda = 0.2, scaling = "noise_sd", penalty = "sobolev"
  )
  expect_s3_class(decon_direct, "crossrun_decon")
  expect_output(print(decon_direct), "<crossrun_decon> 2 runs \\| 7 voxels")
  expect_output(print(decon_direct), "NIfTI geometry attached")
  expect_equal(decon_direct$scaling, decon_manual$scaling)
  expect_equal(decon_direct$noise_innovations, decon_manual$noise_innovations)
  for (r in names(ts$bold)) {
    expect_equal(
      as.matrix(decon_direct$activity[[r]]),
      as.matrix(decon_manual$activity[[r]]),
      ignore_attr = TRUE
    )
    expect_equal(attr(decon_direct$activity[[r]], "voxel_map", exact = TRUE), ts$voxel_map)
    expect_equal(attr(decon_direct$activity[[r]], "geometry", exact = TRUE)$dim, fx$dim3)
  }

  # Zero-reference write_nifti_ts(decon, ...) round-trips 4D activity cleanly
  out_r1 <- file.path(tmp_dir, "decon_r1.nii.gz")
  out_r2 <- file.path(tmp_dir, "decon_r2.nii.gz")
  write_nifti_ts(decon_direct, c(out_r1, out_r2), datatype = "double")
  rt_decon <- read_nifti_ts(c("run-01" = out_r1, "run-02" = out_r2), mask = fx$f_mask, atlas = fx$f_atlas)
  expect_equal(as.numeric(rt_decon$bold[["run-01"]]), as.numeric(decon_direct$activity[["run-01"]]), tolerance = 1e-10)
  expect_equal(as.numeric(rt_decon$bold[["run-02"]]), as.numeric(decon_direct$activity[["run-02"]]), tolerance = 1e-10)
  expect_equal(rt_decon$tr, fx$tr, tolerance = 1e-6)
  expect_equal(rt_decon$geometry$affine, fx$affine, ignore_attr = TRUE, tolerance = 1e-6)

  # Single-run activity matrix and 3D summary vector with `reference = decon_direct` / `hrf_direct`
  out_single <- file.path(tmp_dir, "decon_r1_single.nii.gz")
  write_nifti_ts(decon_direct$activity[["run-01"]], out_single, datatype = "double")
  rt_single <- read_nifti_ts(out_single, mask = fx$f_mask, atlas = fx$f_atlas)
  expect_equal(as.numeric(rt_single$bold[[1L]]), as.numeric(decon_direct$activity[["run-01"]]), tolerance = 1e-10)

  stat_vec <- colMeans(decon_direct$activity[["run-01"]])
  out_stat <- file.path(tmp_dir, "stat_3d.nii.gz")
  write_nifti_ts(stat_vec, out_stat, reference = hrf_direct, datatype = "double")
  rt_stat <- read_nifti_ts(out_stat, mask = fx$f_mask, atlas = fx$f_atlas)
  expect_equal(as.numeric(rt_stat$bold[[1L]]), as.numeric(stat_vec), tolerance = 1e-10)

  # diagnose_crossrun_recovery accepts combined data.frame events directly
  diag_df <- diagnose_crossrun_recovery(decon_direct, events = fx$trial_df, lags = c(0, 1.25))
  diag_list <- diagnose_crossrun_recovery(decon_direct, events = manual_events, lags = c(0, 1.25))
  expect_equal(diag_df$events, diag_list$events)

  # align_medusa_ts(decon_direct, ...) defaults to parcel aggregation and preserves voxel_map for voxel-level calls
  fb_events <- fx$trial_df[fx$trial_df$condition == "feedback", ]
  aligned_parcel <- align_medusa_ts(
    decon_direct, events = fb_events, event = "onset", time_before = -1.25, time_after = 2.5
  )
  expect_equal(sort(unique(aligned_parcel$parcel)), c("10", "20"))
  aligned_vox <- align_medusa_ts(
    decon_direct, events = fb_events, event = "onset", time_before = -1.25, time_after = 2.5,
    group_by = c("trial", "voxel", "x", "y", "z", "atlas_value")
  )
  expect_true(all(c("voxel", "x", "y", "z", "atlas_value") %in% names(aligned_vox)))
  expect_equal(sort(unique(aligned_vox$atlas_value)), c(10, 20))

  # backproject_medusa accepts `brain_mask = decon_direct` (and `ts`) directly
  coef_df <- expand.grid(roi = c(10L, 20L), time = c(0, 1.25), term = "rpe", stringsAsFactors = FALSE)
  coef_df$estimate <- c(0.45, -0.30, 0.80, -0.15)
  bp_dir_file <- file.path(tmp_dir, "bp_file")
  bp_dir_obj <- file.path(tmp_dir, "bp_obj")
  f_bp_file <- backproject_medusa(coef_df, brain_mask = fx$f_atlas, plot_cols = "estimate", output_dir = bp_dir_file, ncpus = 1L)
  f_bp_obj <- backproject_medusa(coef_df, brain_mask = decon_direct, plot_cols = "estimate", output_dir = bp_dir_obj, ncpus = 1L)
  expect_equal(as.numeric(RNifti::readNifti(f_bp_obj)), as.numeric(RNifti::readNifti(f_bp_file)), tolerance = 1e-6)

  # In-memory `nifti_header` fallback succeeds even after original NIfTI files are deleted from disk
  unlink(c(fx$f_runs, fx$f_mask, fx$f_atlas))
  out_after_unlink <- c(file.path(tmp_dir, "unlinked_r1.nii.gz"), file.path(tmp_dir, "unlinked_r2.nii.gz"))
  write_nifti_ts(decon_direct, out_after_unlink, datatype = "double")
  expect_true(all(file.exists(out_after_unlink)))
  bp_dir_unlinked <- file.path(tmp_dir, "bp_unlinked")
  f_bp_unlinked <- backproject_medusa(coef_df, brain_mask = decon_direct, plot_cols = "estimate", output_dir = bp_dir_unlinked, ncpus = 1L)
  expect_equal(as.numeric(RNifti::readNifti(f_bp_unlinked)), as.numeric(RNifti::readNifti(f_bp_file)), tolerance = 1e-6)
})

test_that("Proposal 3.5 edge cases: voxel_by_time, drop = TRUE, spatial_fwhm auto-coords, data.frame nuisance, and multi-ID/session guards", {
  tmp_dir <- tempfile("ergo_35_edge_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  fx <- make_ergonomics_nifti_fixture(tmp_dir)
  k_small <- hrf_library(fx$tr, peak_delays = c(5, 6.5), dispersions = 1, undershoot_ratios = 6)

  # 1. orientation = "voxel_by_time" matches "time_by_voxel"
  ts_tv <- read_nifti_ts(fx$f_runs, mask = fx$f_mask, atlas = fx$f_atlas, orientation = "time_by_voxel")
  ts_vt <- read_nifti_ts(fx$f_runs, mask = fx$f_mask, atlas = fx$f_atlas, orientation = "voxel_by_time")
  fit_tv <- estimate_crossrun_hrf(ts_tv, fx$trial_df, kernel_grid = k_small, pooling = "parcel", min_parcel_voxels = 2L)
  fit_vt <- estimate_crossrun_hrf(ts_vt, fx$trial_df, kernel_grid = k_small, pooling = "parcel", min_parcel_voxels = 2L)
  expect_equal(fit_vt$selections, fit_tv$selections)
  dec_vt <- deconvolve_crossrun_hrf(ts_vt, fit_vt, lambda = 0.25)
  dec_tv <- deconvolve_crossrun_hrf(ts_tv, fit_tv, lambda = 0.25)
  expect_equal(dec_vt$activity, dec_tv$activity)

  # 2. spatial_fwhm auto-populates spatial_coords from ts$voxel_map[, c("x", "y", "z")]
  fit_fwhm <- estimate_crossrun_hrf(
    ts_tv, fx$trial_df, kernel_grid = k_small, spatial_fwhm = 4.0, spatial_strength = 0.05
  )
  expect_equal(
    fit_fwhm$specification$spatial_coords,
    as.matrix(ts_tv$voxel_map[, c("x", "y", "z")]),
    ignore_attr = TRUE
  )

  # 3. drop = TRUE single-run matrix + data.frame nuisance + single-run logical valid_rows
  mat_r1 <- read_nifti_ts(fx$f_runs[1L], mask = fx$f_mask, atlas = fx$f_atlas, drop = TRUE)
  ev_r1 <- fx$trial_df[fx$trial_df$run == 1L, ]
  ev_r1$run <- "run1"
  nuis_df <- data.frame(motion1 = sin(seq_len(nrow(mat_r1)) / 5), motion2 = cos(seq_len(nrow(mat_r1)) / 7))
  valid_vec <- rep(TRUE, nrow(mat_r1))
  valid_vec[c(1L, 25L)] <- FALSE
  fit_single <- estimate_crossrun_hrf(
    mat_r1, events = ev_r1, kernel_grid = k_small, nuisance = nuis_df, valid_rows = valid_vec
  )
  expect_equal(fit_single$run_ids, "run1")
  dec_single <- deconvolve_crossrun_hrf(mat_r1, fit_single, lambda = 0.3)
  out_single_drop <- file.path(tmp_dir, "drop_single.nii.gz")
  write_nifti_ts(dec_single, out_single_drop)
  expect_true(file.exists(out_single_drop))

  # 4. Multi-participant unfiltered events table raises clear error
  bad_events <- fx$trial_df
  bad_events$id[bad_events$run == 2L] <- "sub-02"
  expect_error(
    estimate_crossrun_hrf(ts_tv, events = bad_events, kernel_grid = k_small),
    "multiple participant IDs"
  )

  # 5. Multi-session labels ("ses-01_run-01" vs "ses-02_run-01") do not collide under canonical_run_key()
  ts_ses <- ts_tv
  names(ts_ses$bold) <- c("ses-01_run-01", "ses-02_run-01")
  ev_ses <- fx$trial_df
  ev_ses$run <- ifelse(ev_ses$run == 1L, "ses-01_run-01", "ses-02_run-01")
  fit_ses <- estimate_crossrun_hrf(ts_ses, events = ev_ses, kernel_grid = k_small)
  expect_equal(fit_ses$run_ids, c("ses-01_run-01", "ses-02_run-01"))
})

test_that("Proposal 3.6: decon_epoch_regression supports column-name and formula covariates, column-name trial_durations, and combined events", {
  tmp_dir <- tempfile("ergo_36_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  fx <- make_ergonomics_nifti_fixture(tmp_dir)
  ts <- read_nifti_ts(fx$f_runs, mask = fx$f_mask, atlas = fx$f_atlas)
  k_small <- hrf_library(fx$tr, peak_delays = c(5, 6.5), dispersions = 1, undershoot_ratios = 6)

  hrf_fit <- estimate_crossrun_hrf(
    ts, events = fx$trial_df, kernel_grid = k_small, pooling = "parcel", min_parcel_voxels = 2L
  )
  decon <- deconvolve_crossrun_hrf(ts, hrf_fit, lambda = 0.2, scaling = "noise_sd", penalty = "sobolev")

  # Reference manual per-run covariates sliced to target_condition == "feedback"
  manual_cov_list <- list(
    "run-01" = as.matrix(fx$trial_df[fx$trial_df$run == 1L & fx$trial_df$condition == "feedback", c("rpe", "abs_rpe")]),
    "run-02" = as.matrix(fx$trial_df[fx$trial_df$run == 2L & fx$trial_df$condition == "feedback", c("rpe", "abs_rpe")])
  )
  lags_grid <- c(0, 1.25, 2.5)

  reg_manual <- decon_epoch_regression(
    decon,
    covariates = manual_cov_list,
    target_condition = "feedback",
    lags = lags_grid,
    restore_scale = TRUE,
    aggregate_by_parcel = TRUE,
    whiten_panel = TRUE,
    return_panel = TRUE
  )

  # 1. Column-name covariates extracted automatically from fit$specification$events (with cue NAs ignored!)
  reg_colnames <- decon_epoch_regression(
    decon,
    covariates = c("rpe", "abs_rpe"),
    target_condition = "feedback",
    lags = lags_grid,
    restore_scale = TRUE,
    aggregate_by_parcel = TRUE,
    whiten_panel = TRUE,
    return_panel = TRUE
  )
  expect_equal(reg_colnames$coefficients, reg_manual$coefficients)
  expect_equal(reg_colnames$summary, reg_manual$summary)
  expect_equal(reg_colnames$panel, reg_manual$panel)
  expect_true(all(c("id", "session", "rpe", "abs_rpe", "rpe_whitened", "abs_rpe_whitened") %in% names(reg_colnames$panel)))
  expect_equal(unique(reg_colnames$panel$id), "sub-01")
  expect_equal(unique(reg_colnames$panel$session), "ses-01")

  # 2. Combined multi-run data.frame `events` passed directly to decon_epoch_regression() with newly added column
  trial_df_extra <- fx$trial_df
  trial_df_extra$rpe_sq <- trial_df_extra$rpe^2
  reg_df_events <- decon_epoch_regression(
    decon,
    covariates = c("rpe", "rpe_sq"),
    events = trial_df_extra,
    target_condition = "feedback",
    lags = c(0, 1.25),
    restore_scale = TRUE,
    return_panel = TRUE
  )
  expect_equal(sort(unique(reg_df_events$coefficients$term)), c("rpe", "rpe_sq"))
  expect_true("rpe_sq" %in% names(reg_df_events$panel))

  # 3. One-sided formula covariates (`~ rpe + abs_rpe` and factor interaction `~ rpe * valence`)
  reg_formula <- decon_epoch_regression(
    decon,
    covariates = ~ rpe + abs_rpe,
    target_condition = "feedback",
    lags = lags_grid,
    restore_scale = TRUE,
    aggregate_by_parcel = TRUE,
    whiten_panel = TRUE,
    return_panel = TRUE
  )
  expect_equal(reg_formula$coefficients, reg_manual$coefficients)
  expect_equal(reg_formula$panel, reg_manual$panel)

  reg_interact <- decon_epoch_regression(
    decon,
    covariates = ~ rpe * valence,
    target_condition = "feedback",
    lags = c(0, 1.25)
  )
  expect_true(all(c("rpe", "valencepos", "rpe:valencepos") %in% unique(reg_interact$coefficients$term)))

  # 4. Logical column name in `covariates` coerces cleanly to 0/1 predictor
  reg_logical <- decon_epoch_regression(
    decon,
    covariates = c("rpe", "is_rewarded"),
    target_condition = "feedback",
    lags = c(0, 1.25)
  )
  expect_equal(sort(unique(reg_logical$coefficients$term)), c("is_rewarded", "rpe"))

  # 5. Column-name `trial_durations = "rt"` auto-selects `source_model = "duration"` and matches explicit numeric vector
  fb_rts <- fx$trial_df$rt[fx$trial_df$condition == "feedback"]
  reg_dur_col <- decon_epoch_regression(
    decon,
    covariates = "rpe",
    target_condition = "feedback",
    lags = c(0, 1.25),
    trial_durations = "rt",
    restore_scale = TRUE
  )
  reg_dur_vec <- decon_epoch_regression(
    decon,
    covariates = "rpe",
    target_condition = "feedback",
    lags = c(0, 1.25),
    trial_durations = fb_rts,
    restore_scale = TRUE
  )
  expect_equal(reg_dur_col$settings$source_model, "duration")
  expect_equal(reg_dur_col$coefficients, reg_dur_vec$coefficients)

  # 6. Clear validation errors on missing, non-numeric, or reserved covariate column names
  expect_error(
    decon_epoch_regression(decon, covariates = c("rpe", "missing_col"), target_condition = "feedback"),
    "missing_col"
  )
  expect_error(
    decon_epoch_regression(decon, covariates = "condition", target_condition = "feedback"),
    "must be numeric or logical"
  )
  expect_error(
    decon_epoch_regression(decon, covariates = "trial", target_condition = "feedback"),
    "reserved panel output columns"
  )
})
