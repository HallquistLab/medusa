make_synthetic_bids_fmriprep_tree <- function(root_dir, n_time = 64L, tr = 1.25, seed = 2026L, use_json_tr = FALSE) {
  set.seed(seed)
  dim3 <- c(5L, 4L, 3L)
  n_spatial <- prod(dim3)

  fprep_dir <- file.path(root_dir, "derivatives", "fmriprep")
  fprep_func <- file.path(fprep_dir, "sub-01", "ses-01", "func")
  bids_dir <- file.path(root_dir, "bids")
  bids_func <- file.path(bids_dir, "sub-01", "ses-01", "func")
  dir.create(fprep_func, recursive = TRUE)
  dir.create(bids_func, recursive = TRUE)

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

  write_nii <- function(arr, path, img_tr = tr, datatype = "float") {
    nii <- RNifti::asNifti(arr)
    RNifti::pixdim(nii) <- if (length(dim(arr)) == 4L) {
      c(2, 2, 2.5, img_tr)
    } else {
      c(2, 2, 2.5)
    }
    RNifti::sform(nii) <- affine
    RNifti::writeNifti(nii, path, datatype = datatype)
    path
  }

  # Mask 1 & Mask 2 differ by 1 extra voxel in run 2 to test intersection vs union
  mask1 <- array(0L, dim = dim3)
  mask1[2:4, 1:3, 1:2] <- 1L
  mask2 <- mask1
  mask2[1, 1, 1] <- 1L

  atlas_arr <- array(0L, dim = dim3)
  atlas_arr[2:3, 1:2, 1] <- 10L # 4 voxels in parcel 10
  atlas_arr[4, 1:3, 2] <- 20L   # 3 voxels in parcel 20
  atlas_path <- write_nii(atlas_arr, file.path(root_dir, "parcellation_atlas.nii.gz"), datatype = "int16")

  vox_lin <- which(mask1 > 0L & atlas_arr > 0L)
  k_lib <- hrf_library(tr, peak_delays = c(5, 6.5), dispersions = 1, undershoot_ratios = 6)

  bold_files <- mask_files <- conf_files <- event_files <- character(2L)
  spike_vol <- c(max(5L, as.integer(floor(n_time * 0.35))), max(6L, as.integer(floor(n_time * 0.55))))

  for (r in 1:2) {
    set.seed(seed + r * 10L)
    stem_bids <- sprintf("sub-01_ses-01_task-mid_run-0%d", r)
    stem_space <- paste0(stem_bids, "_space-MNI152NLin2009cAsym")

    # Events TSV with BIDS "n/a" tokens on cue rows and impulse feedback durations
    n_trials <- 6L
    cue_onsets <- seq(5.0, by = 9.0, length.out = n_trials) + (r - 1L) * 0.4
    fb_onsets <- cue_onsets + 2.5
    rpe_vals <- round(seq(-1.4, 1.4, length.out = n_trials), 3)
    rt_vals <- round(seq(0.45, 0.95, length.out = n_trials), 3)

    ev_rows <- vector("list", 2L * n_trials)
    for (i in seq_len(n_trials)) {
      ev_rows[[2L * i - 1L]] <- data.frame(
        onset = cue_onsets[i],
        duration = "0.5",
        trial_type = "cue",
        response_time = "n/a",
        rpe = "n/a",
        stringsAsFactors = FALSE
      )
      ev_rows[[2L * i]] <- data.frame(
        onset = fb_onsets[i],
        duration = "n/a",
        trial_type = if (rpe_vals[i] >= 0) "win_feedback" else "loss_feedback",
        response_time = sprintf("%.3f", rt_vals[i]),
        rpe = sprintf("%.3f", rpe_vals[i]),
        stringsAsFactors = FALSE
      )
    }
    ev_df <- do.call(rbind, ev_rows)
    ev_path <- file.path(bids_func, paste0(stem_bids, "_events.tsv"))
    utils::write.table(ev_df, ev_path, sep = "\t", row.names = FALSE, quote = FALSE)
    event_files[r] <- ev_path

    # Confounds TSV with row-1 "n/a" derivatives, non-steady-state row 1, and 1 FD spike
    tx <- rnorm(n_time, sd = 0.03)
    ty <- rnorm(n_time, sd = 0.03)
    tz <- rnorm(n_time, sd = 0.04)
    rx <- rnorm(n_time, sd = 0.002)
    ry <- rnorm(n_time, sd = 0.002)
    rz <- rnorm(n_time, sd = 0.002)
    tx[spike_vol[r]] <- 1.2
    fd_num <- stats::runif(n_time, 0.04, 0.18)
    fd_num[spike_vol[r]] <- 0.95
    fd_chr <- c("n/a", sprintf("%.4f", fd_num[-1L]))
    dvars_chr <- c("n/a", sprintf("%.4f", stats::runif(n_time - 1L, 0.8, 1.2)))
    nss <- integer(n_time)
    nss[1L] <- 1L

    mk_deriv_chr <- function(v) c("n/a", sprintf("%.6f", diff(v)))
    mk_dpow2_chr <- function(v) c("n/a", sprintf("%.6f", diff(v)^2))

    conf_df <- data.frame(
      csf = rnorm(n_time, mean = 200, sd = 1.5),
      white_matter = rnorm(n_time, mean = 300, sd = 1.2),
      global_signal = rnorm(n_time, mean = 250, sd = 1.8),
      std_dvars = dvars_chr,
      dvars = c("n/a", sprintf("%.3f", stats::runif(n_time - 1L, 15, 25))),
      framewise_displacement = fd_chr,
      rmsd = c("n/a", sprintf("%.4f", stats::runif(n_time - 1L, 0.02, 0.12))),
      a_comp_cor_00 = rnorm(n_time),
      a_comp_cor_01 = rnorm(n_time),
      a_comp_cor_02 = rnorm(n_time),
      c_comp_cor_00 = rnorm(n_time),
      w_comp_cor_00 = rnorm(n_time),
      t_comp_cor_00 = rnorm(n_time),
      cosine00 = cos(pi * (seq_len(n_time) - 0.5) / n_time),
      cosine01 = cos(2 * pi * (seq_len(n_time) - 0.5) / n_time),
      trans_x = tx,
      trans_x_derivative1 = mk_deriv_chr(tx),
      trans_x_power2 = tx^2,
      trans_x_derivative1_power2 = mk_dpow2_chr(tx),
      trans_y = ty,
      trans_y_derivative1 = mk_deriv_chr(ty),
      trans_y_power2 = ty^2,
      trans_y_derivative1_power2 = mk_dpow2_chr(ty),
      trans_z = tz,
      trans_z_derivative1 = mk_deriv_chr(tz),
      trans_z_power2 = tz^2,
      trans_z_derivative1_power2 = mk_dpow2_chr(tz),
      rot_x = rx,
      rot_x_derivative1 = mk_deriv_chr(rx),
      rot_x_power2 = rx^2,
      rot_x_derivative1_power2 = mk_dpow2_chr(rx),
      rot_y = ry,
      rot_y_derivative1 = mk_deriv_chr(ry),
      rot_y_power2 = ry^2,
      rot_y_derivative1_power2 = mk_dpow2_chr(ry),
      rot_z = rz,
      rot_z_derivative1 = mk_deriv_chr(rz),
      rot_z_power2 = rz^2,
      rot_z_derivative1_power2 = mk_dpow2_chr(rz),
      non_steady_state_outlier00 = nss,
      motion_outlier00 = as.integer(seq_len(n_time) == spike_vol[r]),
      aroma_motion_01 = rnorm(n_time),
      constant_zero_col = rep(0, n_time),
      stringsAsFactors = FALSE
    )
    conf_path <- file.path(fprep_func, paste0(stem_bids, "_desc-confounds_timeseries.tsv"))
    utils::write.table(conf_df, conf_path, sep = "\t", row.names = FALSE, quote = FALSE)
    conf_files[r] <- conf_path

    # 4D BOLD NIfTI with parcel-specific HRF signal + motion artifact at spike_vol[r]
    arr <- array(rnorm(n_spatial * n_time, mean = 100, sd = 0.35), dim = c(dim3, n_time))
    t_grid <- (seq_len(n_time) - 1L) * tr
    for (j in seq_along(vox_lin)) {
      parcel_j <- atlas_arr[vox_lin[j]]
      h <- k_lib[, if (parcel_j == 10L) 1L else 2L]
      sig <- numeric(n_time)
      for (i in seq_len(n_trials)) {
        dt_c <- t_grid - cue_onsets[i]
        act_c <- which(dt_c >= 0 & dt_c <= (length(h) - 1L) * tr)
        if (length(act_c)) {
          sig[act_c] <- sig[act_c] + 0.8 * stats::approx((seq_along(h) - 1L) * tr, h, xout = dt_c[act_c], rule = 2)$y
        }
        dt_f <- t_grid - fb_onsets[i]
        act_f <- which(dt_f >= 0 & dt_f <= (length(h) - 1L) * tr)
        if (length(act_f)) {
          sig[act_f] <- sig[act_f] + (1.2 + 0.9 * rpe_vals[i]) *
            stats::approx((seq_along(h) - 1L) * tr, h, xout = dt_f[act_f], rule = 2)$y
        }
      }
      sig[spike_vol[r]] <- sig[spike_vol[r]] + 25.0 # huge motion spike on censored frame
      ijk <- arrayInd(vox_lin[j], dim3)
      arr[ijk[1L], ijk[2L], ijk[3L], ] <- arr[ijk[1L], ijk[2L], ijk[3L], ] + sig + 2.0 * tx
    }

    bold_path <- file.path(fprep_func, paste0(stem_space, "_desc-preproc_bold.nii.gz"))
    write_nii(arr, bold_path, img_tr = if (use_json_tr) 0 else tr)
    if (use_json_tr) {
      json_path <- sub("\\.nii\\.gz$", ".json", bold_path)
      writeLines(sprintf('{\n  "RepetitionTime": %.4f,\n  "TaskName": "mid"\n}', tr), json_path)
    }
    bold_files[r] <- bold_path

    mask_path <- file.path(fprep_func, paste0(stem_space, "_desc-brain_mask.nii.gz"))
    write_nii(if (r == 1L) mask1 else mask2, mask_path, datatype = "uint8")
    mask_files[r] <- mask_path
  }

  list(
    fmriprep_dir = fprep_dir,
    bids_dir = bids_dir,
    atlas_path = atlas_path,
    bold_files = bold_files,
    mask_files = mask_files,
    conf_files = conf_files,
    event_files = event_files,
    spike_vol = spike_vol,
    tr = tr,
    n_time = n_time
  )
}

test_that("read_bids_events normalizes BIDS TSVs, preserves numeric types across n/a, and aligns run order", {
  tmp_dir <- tempfile("bids_ev_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  f1 <- file.path(tmp_dir, "sub-01_ses-02_task-mid_run-01_events.tsv")
  f2 <- file.path(tmp_dir, "sub-01_ses-02_task-mid_run-02_events.tsv")

  # Note: rpe is ALL "n/a" in run 1 (fread reads as logical NA) and numeric in run 2!
  writeLines(c(
    "onset\tduration\ttrial_type\tresponse_time\trpe",
    "6.0\t0.8\twin_fb\t0.52\tn/a",
    "2.0\tn/a\tcue\tn/a\tn/a",
    "10.0\t0.5\tfixation\tn/a\tn/a"
  ), f1)
  writeLines(c(
    "onset\tduration\ttrial_type\tresponse_time\trpe",
    "3.0\t0.5\tcue\tn/a\tn/a",
    "7.5\tn/a\tloss_fb\t0.61\t-0.85"
  ), f2)

  # Supply files in reverse order (run-02 then run-01) against target BOLD run names run1, run2
  ev <- read_bids_events(
    events = c(f2, f1),
    bold = c("run1", "run2"),
    condition_map = c(win_fb = "feedback", loss_fb = "feedback"),
    exclude_conditions = "fixation",
    na_duration = 0,
    time_shift = c(run1 = 0.5, run2 = -0.5),
    sort_by_onset = TRUE
  )

  expect_s3_class(ev, "data.frame")
  expect_identical(ev$run, c("run1", "run1", "run2", "run2"))
  expect_identical(ev$id, rep("01", 4L))
  expect_identical(ev$session, rep("02", 4L))
  expect_identical(ev$task, rep("mid", 4L))
  expect_identical(ev$trial, c(1L, 2L, 1L, 2L))
  # Chronological sorting within run1 + original file row tracking + per-run time_shift
  expect_equal(ev$onset, c(2.5, 6.5, 2.5, 7.0))
  expect_identical(ev$bids_event_row[ev$run == "run1"], c(2L, 1L))
  expect_equal(ev$duration, c(0.0, 0.8, 0.5, 0.0))
  expect_identical(ev$condition, c("cue", "feedback", "cue", "feedback"))
  # Numeric behavioral columns must remain double, never character or logical
  expect_type(ev$response_time, "double")
  expect_type(ev$rpe, "double")
  expect_equal(ev$rpe, c(NA_real_, NA_real_, NA_real_, -0.85))

  # na_duration = NULL errors on "n/a" duration
  expect_error(read_bids_events(f1, na_duration = NULL), "Missing 'duration'")

  # Single-condition table without trial_type uses default_condition
  f_single <- file.path(tmp_dir, "sub-01_task-motor_run-01_events.tsv")
  writeLines(c("onset\tduration", "1.0\t0.2", "5.0\t0.2"), f_single)
  expect_error(read_bids_events(f_single), "Neither 'condition' nor 'trial_type'")
  ev_single <- read_bids_events(f_single, default_condition = "tap")
  expect_identical(ev_single$condition, c("tap", "tap"))

  # Missing onset / condition handling
  f_bad <- file.path(tmp_dir, "sub-01_task-mid_run-01_bad_events.tsv")
  writeLines(c("onset\tduration\ttrial_type", "n/a\t0.5\tcue", "4.0\t0.5\tn/a", "8.0\t0.5\tcue"), f_bad)
  expect_error(read_bids_events(f_bad), "Non-finite or missing 'onset'")
  expect_error(read_bids_events(f_bad, drop_missing_onset = TRUE), "Missing or empty condition")
  ev_clean <- read_bids_events(f_bad, drop_missing_onset = TRUE, drop_missing_condition = TRUE)
  expect_equal(nrow(ev_clean), 1L)
  expect_equal(ev_clean$onset, 8.0)

  # Empty run after exclusion preserves factor levels for multi-run downstream splitting
  ev_empty_r1 <- read_bids_events(
    c(f1, f2),
    bold = c("run1", "run2"),
    exclude_conditions = c("win_fb", "cue", "fixation")
  )
  expect_s3_class(ev_empty_r1$run, "factor")
  expect_identical(levels(ev_empty_r1$run), c("run1", "run2"))
  split_ev <- medusa:::split_multirun_events(ev_empty_r1, c("run1", "run2"), allow_empty_runs = TRUE)
  expect_equal(nrow(split_ev$run1), 0L)
  expect_equal(nrow(split_ev$run2), 1L)
})

test_that("read_bids_events and read_fmriprep_confounds enforce BIDS entity and volume-count safety checks", {
  fix <- make_synthetic_bids_fmriprep_tree(tempfile("bids_safety_"), n_time = 20L)
  on.exit(unlink(dirname(fix$bids_dir), recursive = TRUE), add = TRUE)

  ts <- read_nifti_ts(
    niftis = c(run1 = fix$bold_files[1L], run2 = fix$bold_files[2L]),
    mask = fix$mask_files[1L],
    tr = fix$tr
  )

  # Wrong subject in events path -> explicit BIDS sub- mismatch error
  wrong_sub_ev <- file.path(dirname(fix$bids_dir), "sub-02_ses-01_task-mid_run-01_events.tsv")
  file.copy(fix$event_files[1L], wrong_sub_ev)
  expect_error(
    read_bids_events(c(wrong_sub_ev, fix$event_files[2L]), bold = ts),
    "BIDS 'sub-' entity mismatch"
  )

  # Wrong task in confounds path -> explicit BIDS task- mismatch error
  wrong_task_conf <- file.path(dirname(fix$bids_dir), "sub-01_ses-01_task-rest_run-01_desc-confounds_timeseries.tsv")
  file.copy(fix$conf_files[1L], wrong_task_conf)
  expect_error(
    read_fmriprep_confounds(c(wrong_task_conf, fix$conf_files[2L]), bold = ts),
    "BIDS 'task-' entity mismatch"
  )

  # Volume count mismatch (e.g., if user pre-deleted volumes!) -> explicit clock-preservation error
  trunc_df <- utils::read.delim(fix$conf_files[1L], stringsAsFactors = FALSE)[1:18, ]
  trunc_conf <- file.path(dirname(fix$bids_dir), "sub-01_ses-01_task-mid_run-01_desc-confounds_timeseries.tsv")
  utils::write.table(trunc_df, trunc_conf, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(
    read_fmriprep_confounds(c(trunc_conf, fix$conf_files[2L]), bold = ts),
    "Do not delete volumes prior to MEDuSA"
  )
})

test_that("read_fmriprep_confounds builds finite nuisance strategies and clock-preserving censor masks", {
  fix <- make_synthetic_bids_fmriprep_tree(tempfile("conf_strat_"), n_time = 30L)
  on.exit(unlink(dirname(fix$bids_dir), recursive = TRUE), add = TRUE)

  # 1. Default 24p strategy with reverse-ordered files matched to run1, run2
  conf_24p <- read_fmriprep_confounds(
    c(fix$conf_files[2L], fix$conf_files[1L]),
    run_names = c("run1", "run2")
  )
  # Without BOLD file entities, unnamed BIDS run-02, run-01 reorders to run1, run2!
  expect_s3_class(conf_24p, "medusa_confounds")
  expect_identical(names(conf_24p$nuisance), c("run1", "run2"))
  expect_equal(nrow(conf_24p$nuisance$run1), 30L)
  expect_equal(ncol(conf_24p$nuisance$run1), 24L)
  expect_true(all(is.finite(conf_24p$nuisance$run1)))
  expect_true(all(is.finite(conf_24p$nuisance$run2)))

  # Check exact censoring on run 1: non-steady-state vol 1 + FD spike at spike_vol[1] (22)
  expect_false(conf_24p$valid_rows$run1[1L])
  expect_false(conf_24p$valid_rows$run1[fix$spike_vol[1L]])
  expect_equal(sum(!conf_24p$valid_rows$run1), 2L)
  # Censored rows have zero leverage (0 after usable-row mean imputation + centering)
  expect_equal(unname(conf_24p$nuisance$run1[fix$spike_vol[1L], ]), rep(0, 24L), tolerance = 1e-12)
  expect_equal(unname(colMeans(conf_24p$nuisance$run1[conf_24p$valid_rows$run1, ])), rep(0, 24L), tolerance = 1e-12)

  # 2. Spike dilation (censor_prev = 1, censor_next = 2) dilates FD spike ONLY, NOT non-steady-state row 1!
  conf_dil <- read_fmriprep_confounds(
    fix$conf_files[1L],
    strategy = "6p",
    fd_threshold = 0.5,
    censor_prev = 1L,
    censor_next = 2L
  )
  # Row 1 is non-steady-state, so row 2 & 3 MUST remain TRUE (not dilated from row 1!)
  expect_false(conf_dil$valid_rows$run1[1L])
  expect_true(conf_dil$valid_rows$run1[2L])
  expect_true(conf_dil$valid_rows$run1[3L])
  # Spike at sp1 dilates to (sp1 - 1):(sp1 + 2)
  sp1 <- fix$spike_vol[1L]
  expect_false(any(conf_dil$valid_rows$run1[(sp1 - 1L):(sp1 + 2L)]))
  expect_true(conf_dil$valid_rows$run1[sp1 - 2L])
  expect_true(conf_dil$valid_rows$run1[sp1 + 3L])

  # 3. min_segment_length scrubs short isolated usable islands between censored frames
  # Spike at sp1 is FALSE; if drop_initial_volumes = sp1 - 3, rows (sp1 - 2):(sp1 - 1) form a 2-TR island
  conf_seg <- read_fmriprep_confounds(
    fix$conf_files[1L],
    strategy = "6p",
    drop_initial_volumes = sp1 - 3L,
    min_segment_length = 4L
  )
  expect_false(any(conf_seg$valid_rows$run1[seq_len(sp1)]))
  expect_true(all(conf_seg$valid_rows$run1[(sp1 + 1L):30L]))

  # 4. Synthesized 24p & 36p expansions when only base columns exist + row-1 derivative finiteness when row 1 is usable
  minimal_df <- data.frame(
    trans_x = seq(0.1, 1.0, length.out = 10L),
    trans_y = rnorm(10L, sd = 0.05),
    trans_z = rnorm(10L, sd = 0.05),
    rot_x = rnorm(10L, sd = 0.005),
    rot_y = rnorm(10L, sd = 0.005),
    rot_z = rnorm(10L, sd = 0.005),
    csf = rnorm(10L),
    white_matter = rnorm(10L),
    global_signal = rnorm(10L),
    framewise_displacement = c(NA_real_, rep(0.1, 9L))
  )
  conf_36p <- read_fmriprep_confounds(minimal_df, strategy = "36p", include_non_steady_state = FALSE)
  expect_true(conf_36p$valid_rows$run1[1L]) # row 1 is usable!
  expect_equal(ncol(conf_36p$nuisance$run1), 36L)
  expect_true(all(is.finite(conf_36p$nuisance$run1)))

  # 5. CompCor + auto cosine selection + constant-column dropping in custom mode
  conf_cc <- read_fmriprep_confounds(fix$conf_files[1L], strategy = "compcor", n_compcor = 2L)
  expect_true(all(c("a_comp_cor_00", "a_comp_cor_01", "cosine00", "cosine01") %in% colnames(conf_cc$nuisance$run1)))
  expect_false("a_comp_cor_02" %in% colnames(conf_cc$nuisance$run1))

  conf_split <- read_fmriprep_confounds(fix$conf_files[1L], strategy = "compcor", n_compcor = 1L, compcor_type = "split")
  expect_true(all(c("c_comp_cor_00", "w_comp_cor_00", "cosine00") %in% colnames(conf_split$nuisance$run1)))

  conf_custom <- read_fmriprep_confounds(
    fix$conf_files[1L],
    strategy = "custom",
    columns = c("csf", "constant_zero_col"),
    column_pattern = "^aroma_motion_"
  )
  expect_identical(colnames(conf_custom$nuisance$run1), c("csf", "aroma_motion_01"))

  conf_none <- read_fmriprep_confounds(fix$conf_files[1L], strategy = "none")
  expect_equal(ncol(conf_none$nuisance$run1), 0L)

  out_txt <- paste(capture.output(print(conf_24p)), collapse = "\n")
  expect_match(out_txt, "<medusa_confounds> 2 runs")
})

test_that("read_fmriprep_subject runs end-to-end across HRF estimation, deconvolution, epoch regression, and alignment", {
  fix <- make_synthetic_bids_fmriprep_tree(
    tempfile("fprep_e2e_"),
    n_time = 64L,
    tr = 1.25,
    use_json_tr = TRUE # NIfTI header pixdim[4] = 0 -> must recover TR = 1.25 from *_bold.json!
  )
  on.exit(unlink(dirname(fix$bids_dir), recursive = TRUE), add = TRUE)

  cond_map <- c(win_feedback = "feedback", loss_feedback = "feedback")

  # 1. One-call fMRIPrep + BIDS subject discovery
  sub01 <- read_fmriprep_subject(
    fmriprep_dir = fix$fmriprep_dir,
    bids_dir = fix$bids_dir,
    sub = "01",
    ses = "01",
    task = "mid",
    atlas = fix$atlas_path,
    condition_map = cond_map,
    strategy = "24p",
    fd_threshold = 0.5
  )

  expect_s3_class(sub01, "medusa_nifti_ts")
  expect_equal(sub01$tr, 1.25)
  expect_equal(nrow(sub01$voxel_map), 7L) # 4 in parcel 10 + 3 in parcel 20
  expect_s3_class(sub01$confounds, "medusa_confounds")
  expect_equal(nrow(sub01$events), 24L)
  expect_match(paste(capture.output(print(sub01)), collapse = "\n"), "Confounds attached")
  expect_match(paste(capture.output(print(sub01)), collapse = "\n"), "BIDS events attached")

  # Verify mask_mode = "intersection" vs "union" when atlas is omitted
  sub_inter <- read_fmriprep_subject(fix$fmriprep_dir, sub = "sub-01", events = NULL, mask_mode = "intersection")
  sub_union <- read_fmriprep_subject(fix$fmriprep_dir, sub = "sub-01", events = NULL, mask_mode = "union")
  expect_equal(nrow(sub_inter$voxel_map), 18L)
  expect_equal(nrow(sub_union$voxel_map), 19L)

  # 2. End-to-end MEDuSA pipeline directly from `sub01` with zero manual unpacking
  k_grid <- hrf_library(1.25, peak_delays = c(5, 6.5), dispersions = 1, undershoot_ratios = 6)
  hrf_fit <- estimate_crossrun_hrf(sub01, kernel_grid = k_grid, pooling = "parcel", min_parcel_voxels = 2L)
  decon <- deconvolve_crossrun_hrf(sub01, hrf_fit, scaling = "noise_sd", penalty = "sobolev")
  fb_fit <- decon_epoch_regression(
    decon,
    covariates = "rpe",
    target_condition = "feedback",
    lags = seq(-1, 4, by = 1),
    aggregate_by_parcel = TRUE,
    whiten_panel = TRUE,
    return_panel = TRUE
  )
  expect_true(all(is.finite(fb_fit$coefficients$beta_gls)))
  expect_true("rpe_whitened" %in% names(fb_fit$panel))

  # 3. Exact numerical equivalence against modular helpers & explicit list arguments
  ts_plain <- read_nifti_ts(
    niftis = c(run1 = fix$bold_files[1L], run2 = fix$bold_files[2L]),
    mask = sub01$geometry$mask,
    atlas = fix$atlas_path,
    tr = 1.25
  )
  ev_mod <- read_bids_events(fix$event_files, bold = ts_plain, condition_map = cond_map)
  conf_mod <- read_fmriprep_confounds(fix$conf_files, bold = ts_plain, strategy = "24p", fd_threshold = 0.5)

  hrf_mod <- estimate_crossrun_hrf(
    ts_plain,
    events = ev_mod,
    kernel_grid = k_grid,
    nuisance = conf_mod,
    pooling = "parcel",
    min_parcel_voxels = 2L
  )
  hrf_explicit <- estimate_crossrun_hrf(
    ts_plain$bold,
    events = split(ev_mod, ev_mod$run),
    tr = 1.25,
    kernel_grid = k_grid,
    nuisance = conf_mod$nuisance,
    valid_rows = conf_mod$valid_rows,
    spatial_groups = ts_plain$parcel_ids,
    pooling = "parcel",
    min_parcel_voxels = 2L
  )

  expect_equal(hrf_fit$selected_kernels, hrf_mod$selected_kernels)
  expect_equal(hrf_fit$selected_kernels, hrf_explicit$selected_kernels)
  expect_equal(hrf_fit$noise, hrf_mod$noise)

  # 4. Direct raw-BOLD matrix alignment uses sub01$valid_rows automatically to fence the motion spike!
  # Spike in run 1 is at volume 22 -> t = 21 * 1.25 = 26.25 s (near trial 3 feedback at 25.5 s)
  al_cens <- align_medusa_ts(sub01, event = "onset", time_before = 0, time_after = 2.5, output_resolution = 1.25)
  al_mod <- align_medusa_ts(ts_plain, events = ev_mod, event = "onset", valid = conf_mod,
                            time_before = 0, time_after = 2.5, output_resolution = 1.25)
  al_uncens <- align_medusa_ts(sub01, event = "onset", valid = NULL,
                               time_before = 0, time_after = 2.5, output_resolution = 1.25)

  expect_equal(al_cens$value_mean, al_mod$value_mean)
  expect_gt(sum(al_cens$value_n_censored), 0L)
  expect_equal(sum(al_uncens$value_n_censored), 0L)
})

test_that("read_fmriprep_subject supports BIDS Inheritance Principle and disambiguates multi-task / multi-space dirs", {
  fix <- make_synthetic_bids_fmriprep_tree(tempfile("fprep_inh_"), n_time = 16L)
  on.exit(unlink(dirname(fix$bids_dir), recursive = TRUE), add = TRUE)

  # Remove run-level event files and place a single inherited task-mid_events.tsv at BIDS root
  inherited_ev <- file.path(fix$bids_dir, "task-mid_events.tsv")
  file.copy(fix$event_files[1L], inherited_ev)
  unlink(fix$event_files)

  sub_inh <- read_fmriprep_subject(
    fmriprep_dir = fix$fmriprep_dir,
    bids_dir = fix$bids_dir,
    sub = "01",
    task = "mid",
    strategy = "6p"
  )
  expect_equal(nrow(sub_inh$events), 24L)
  expect_identical(sort(unique(sub_inh$events$run)), c("run1", "run2"))

  # Add a second output space file -> omitting 'space' must prompt explicit disambiguation
  t1w_bold <- sub("space-MNI152NLin2009cAsym", "space-T1w", fix$bold_files[1L])
  file.copy(fix$bold_files[1L], t1w_bold)
  expect_error(
    read_fmriprep_subject(fix$fmriprep_dir, sub = "01", task = "mid", events = NULL),
    "Multiple output spaces found"
  )
  sub_mni <- read_fmriprep_subject(
    fix$fmriprep_dir, sub = "01", task = "mid", space = "MNI152NLin2009cAsym", events = NULL, strategy = "6p"
  )
  expect_equal(length(sub_mni$bold), 2L)

  # Add a second task file -> omitting 'task' must prompt explicit disambiguation
  rest_bold <- sub("task-mid", "task-rest", fix$bold_files[1L])
  file.copy(fix$bold_files[1L], rest_bold)
  expect_error(
    read_fmriprep_subject(fix$fmriprep_dir, sub = "01", space = "MNI152NLin2009cAsym", events = NULL),
    "Multiple BIDS tasks found"
  )
})

test_that("Tier A confound enhancements: '.'/'NaN' NA tokens, Power FD fallback, 27p strategy, and <0-k> range specs", {
  tmp_dir <- tempfile("tierA_conf_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # 1. '.' and 'NaN' missing tokens in events and confounds TSVs
  ev_file <- file.path(tmp_dir, "sub-01_task-mid_run-01_events.tsv")
  writeLines(c(
    "onset\tduration\ttrial_type\tresponse_time\trpe",
    "2.0\t.\tcue\t.\tNaN",
    "6.0\t0.5\tfeedback\t0.55\t0.90"
  ), ev_file)
  ev <- read_bids_events(ev_file, na_duration = 0)
  expect_equal(ev$duration, c(0, 0.5))
  expect_equal(ev$response_time, c(NA_real_, 0.55))
  expect_equal(ev$rpe, c(NA_real_, 0.90))

  # 2. Power (2012) FD fallback when framewise_displacement is absent or all NA/NaN
  n_t <- 8L
  mot_df <- data.frame(
    trans_x = c(0.0, 0.1, 0.1, 0.5, 0.5, 0.5, 0.5, 0.5),
    trans_y = rep(0.0, n_t),
    trans_z = c(0.0, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2),
    rot_x = c(0.0, 0.0, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004), # diff = 0.004 rad -> * 50 = 0.20 mm at vol 3
    rot_y = rep(0.0, n_t),
    rot_z = rep(0.0, n_t),
    csf = seq_len(n_t) * 0.1,
    white_matter = cos(seq_len(n_t)),
    global_signal = sin(seq_len(n_t)),
    a_comp_cor_00 = rnorm(n_t),
    a_comp_cor_01 = rnorm(n_t),
    a_comp_cor_02 = rnorm(n_t),
    cosine00 = cos(pi * (seq_len(n_t) - 0.5) / n_t),
    cosine01 = cos(2 * pi * (seq_len(n_t) - 0.5) / n_t)
  )
  # Expected Power FD: vol1=0, vol2=0.1, vol3=0.004*50=0.2, vol4=0.4+0.2=0.6 (spike > 0.5!), vol5..8=0
  conf_fb <- read_fmriprep_confounds(mot_df, strategy = "27p", fd_threshold = 0.5, drop_constant = FALSE)
  expect_equal(ncol(conf_fb$nuisance$run1), 27L)
  expect_true(all(c("csf", "white_matter", "global_signal") %in% colnames(conf_fb$nuisance$run1)))
  expect_false(any(c("csf_derivative1", "white_matter_power2", "global_signal_derivative1") %in% colnames(conf_fb$nuisance$run1)))
  expect_equal(conf_fb$framewise_displacement$run1[1:5], c(0, 0.1, 0.2, 0.6, 0.0), tolerance = 1e-12)
  expect_identical(conf_fb$valid_rows$run1, c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE))

  # Also falls back when framewise_displacement column exists in TSV with '.' and 'NaN'
  mot_tsv <- file.path(tmp_dir, "sub-01_task-mid_run-01_confounds.tsv")
  mot_df_na <- mot_df
  mot_df_na$framewise_displacement <- c("n/a", ".", "NaN", "NA", ".", ".", ".", ".")
  utils::write.table(mot_df_na, mot_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  conf_fb2 <- read_fmriprep_confounds(mot_tsv, strategy = "6p", fd_threshold = 0.5)
  expect_equal(conf_fb2$framewise_displacement$run1, conf_fb$framewise_displacement$run1, tolerance = 1e-12)

  # Degree conversion in Power FD fallback + explicit custom missing fd_col error
  mot_deg <- mot_df
  mot_deg$rot_x <- mot_df$rot_x * (180 / pi)
  conf_deg <- read_fmriprep_confounds(mot_deg, strategy = "6p", rot_units = "deg", fd_threshold = 0.5)
  expect_equal(conf_deg$framewise_displacement$run1, conf_fb$framewise_displacement$run1, tolerance = 1e-12)
  expect_error(read_fmriprep_confounds(mot_df, fd_col = "nonexistent_fd"), "Column 'nonexistent_fd' not found")

  # 3. Shorthand + <0-2> zero-padded range expansion in 'columns'
  conf_rng <- read_fmriprep_confounds(
    mot_df,
    strategy = "custom",
    columns = c("12p", "a_comp_cor_<0-1>", "cosine<0,1>"),
    include_cosines = FALSE
  )
  expect_true(all(c("trans_x", "trans_x_derivative1", "a_comp_cor_00", "a_comp_cor_01", "cosine00", "cosine01") %in%
                    colnames(conf_rng$nuisance$run1)))
  expect_false("a_comp_cor_02" %in% colnames(conf_rng$nuisance$run1))
  expect_error(
    read_fmriprep_confounds(mot_df, strategy = "custom", columns = "a_comp_cor_<bad>"),
    "Invalid column range specification"
  )
})

test_that("BrainGnomes postprocessing handoff: censor.1D, scrub.tsv, filtered FD, cohort/templatemask, and provenance guards", {
  fix <- make_synthetic_bids_fmriprep_tree(tempfile("bg_postproc_"), n_time = 24L, tr = 1.0)
  on.exit(unlink(dirname(fix$bids_dir), recursive = TRUE), add = TRUE)

  # 1. External *_censor.1D and *_scrub.tsv + framewise_displacement_filtered in read_fmriprep_confounds
  c_df <- utils::read.delim(fix$conf_files[1L], stringsAsFactors = FALSE)
  # Simulate respiratory pseudomotion: raw FD = 0.72 at vol 10, but filtered FD = 0.12
  c_df$framewise_displacement[10L] <- "0.7200"
  c_df$framewise_displacement_filtered <- c_df$framewise_displacement
  c_df$framewise_displacement_filtered[10L] <- "0.1200"
  bg_conf_tsv <- file.path(dirname(fix$bids_dir), "sub-01_ses-01_task-mid_run-01_desc-postproc_confounds.tsv")
  utils::write.table(c_df, bg_conf_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

  # fd_col = "auto" picks framewise_displacement_filtered -> vol 10 stays TRUE!
  conf_filt <- read_fmriprep_confounds(bg_conf_tsv, strategy = "6p", fd_threshold = 0.5)
  expect_true(conf_filt$valid_rows$run1[10L])
  expect_false(conf_filt$valid_rows$run1[fix$spike_vol[1L]])

  # *_censor.1D with comment header: censor volumes 5 and 6 only (leave raw FD spike uncensored unless fd_threshold supplied)
  censor_vec <- rep(1L, 24L)
  censor_vec[c(5L, 6L)] <- 0L
  censor_1d <- file.path(dirname(fix$bids_dir), "sub-01_ses-01_task-mid_run-01_desc-postproc_censor.1D")
  writeLines(c("# BrainGnomes censor file (1=keep, 0=censor)", as.character(censor_vec)), censor_1d)

  conf_c1d <- read_fmriprep_confounds(bg_conf_tsv, censor = censor_1d, strategy = "6p")
  expect_false(conf_c1d$valid_rows$run1[1L]) # non_steady_state_outlier00
  expect_false(any(conf_c1d$valid_rows$run1[5:6]))
  expect_true(conf_c1d$valid_rows$run1[fix$spike_vol[1L]]) # default fd_threshold disabled when censor supplied
  expect_equal(conf_c1d$summary$n_external_censor, 2L)

  # Explicit fd_threshold + censor intersects both
  conf_both <- read_fmriprep_confounds(bg_conf_tsv, censor = censor_1d, fd_threshold = 0.5, strategy = "6p")
  expect_false(any(conf_both$valid_rows$run1[c(5L, 6L, fix$spike_vol[1L])]))

  # *_scrub.tsv one-hot spike indicator matrix
  scrub_df <- data.frame(
    spike_1 = as.integer(seq_len(24L) == 7L),
    spike_2 = as.integer(seq_len(24L) == 8L)
  )
  scrub_tsv <- file.path(dirname(fix$bids_dir), "sub-01_ses-01_task-mid_run-01_desc-postproc_scrub.tsv")
  utils::write.table(scrub_df, scrub_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  conf_scrub <- read_fmriprep_confounds(bg_conf_tsv, censor = scrub_tsv, strategy = "none", include_non_steady_state = FALSE)
  expect_identical(which(!conf_scrub$valid_rows$run1), c(7L, 8L))

  # 2. Full BrainGnomes flat postproc/sub-01/ layout + cohort + templatemask + intensityReferenceCore exclusion + provenance
  postproc_root <- file.path(dirname(fix$bids_dir), "postproc")
  postproc_sub <- file.path(postproc_root, "sub-01")
  dir.create(postproc_sub, recursive = TRUE)

  # Create a tiny decoy 1-voxel intensityReferenceCore_mask that MUST be ignored in favor of 18-voxel templatemask!
  core_arr <- array(0L, dim = c(5L, 4L, 3L))
  core_arr[2L, 1L, 1L] <- 1L
  core_nii <- RNifti::updateNifti(core_arr, template = RNifti::readNifti(fix$mask_files[1L]))

  for (r in 1:2) {
    stem_pp <- sprintf("sub-01_ses-01_task-mid_run-0%d_space-MNIPediatricAsym_cohort-2", r)
    pp_bold <- file.path(postproc_sub, paste0(stem_pp, "_desc-postproc_bold.nii.gz"))
    file.copy(fix$bold_files[r], pp_bold)

    # Decoy cohort-1 BOLD in run 1 to test cohort disambiguation!
    if (r == 1L) {
      coh1_bold <- file.path(postproc_sub, sprintf("sub-01_ses-01_task-mid_run-01_space-MNIPediatricAsym_cohort-1_desc-postproc_bold.nii.gz"))
      file.copy(fix$bold_files[r], coh1_bold)
    }

    # TemplateFlow resampled mask (*_templatemask.nii.gz) + decoy intensityReferenceCore_mask
    tmpl_mask <- file.path(postproc_sub, paste0(stem_pp, "_desc-preproc_templatemask.nii.gz"))
    file.copy(fix$mask_files[1L], tmpl_mask)
    decoy_core <- file.path(postproc_sub, paste0(stem_pp, "_desc-intensityReferenceCore_mask.nii.gz"))
    RNifti::writeNifti(core_nii, decoy_core)

    # Companion *_censor.1D and *_desc-postproc_confounds.tsv
    c_vec <- rep(1L, 24L)
    c_vec[fix$spike_vol[r]] <- 0L
    writeLines(as.character(c_vec), file.path(postproc_sub, paste0(stem_pp, "_desc-postproc_censor.1D")))
    file.copy(fix$conf_files[r], file.path(postproc_sub, paste0(stem_pp, "_desc-postproc_confounds.tsv")))

    # BrainGnomes derivative provenance JSON sidecar (confound_regression executed, scrub_timepoints skipped)
    prov_json <- file.path(postproc_sub, paste0(stem_pp, "_desc-postproc_bold.json"))
    writeLines(c(
      "{",
      '  "Description": "BrainGnomes postprocessing derivative",',
      '  "RepetitionTime": 1.0,',
      '  "BrainGnomes": {',
      '    "SchemaVersion": "brain-gnomes-derivative-v1",',
      '    "Stage": "postprocessing",',
      '    "Operations": [',
      '      { "Name": "spatial_smooth", "Status": "executed" },',
      '      { "Name": "confound_regression", "Status": "executed" },',
      '      { "Name": "scrub_timepoints", "Status": "skipped" }',
      "    ],",
      '    "SourceFiles": [',
      sprintf('      { "role": "confounds", "path": "%s" }', fix$conf_files[r]),
      "    ],",
      '    "SourceMetadata": {',
      '      "BrainGnomes": {',
      '        "Operations": [ { "Name": "scrub_timepoints", "Status": "executed" } ]',
      "      }",
      "    }",
      "  }",
      "}"
    ), prov_json)
  }

  # Verify both jsonlite and pure-R provenance parsers ignore nested SourceMetadata and detect confound_regression
  p_jl <- medusa:::inspect_braingnomes_provenance(
    file.path(postproc_sub, "sub-01_ses-01_task-mid_run-01_space-MNIPediatricAsym_cohort-2_desc-postproc_bold.nii.gz"),
    use_jsonlite = TRUE
  )
  p_rx <- medusa:::inspect_braingnomes_provenance(
    file.path(postproc_sub, "sub-01_ses-01_task-mid_run-01_space-MNIPediatricAsym_cohort-2_desc-postproc_bold.nii.gz"),
    use_jsonlite = FALSE
  )
  expect_identical(p_jl$confound_regression, TRUE)
  expect_identical(p_rx$confound_regression, TRUE)
  expect_identical(p_jl$scrub_timepoints, FALSE)
  expect_identical(p_rx$scrub_timepoints, FALSE)
  expect_identical(p_rx$upstream_confounds, fix$conf_files[1L])

  # Multi-cohort ambiguity check without 'cohort'
  expect_error(
    read_fmriprep_subject(postproc_root, sub = "01", desc = "postproc", events = NULL),
    "Multiple template cohorts found"
  )

  # Ingest cohort = "2": auto-uses templatemask (18 voxels, NOT 1-voxel intensityReferenceCore!),
  # auto-uses *_censor.1D, and auto-defaults strategy = "none" because confound_regression already ran!
  sub_pp <- read_fmriprep_subject(
    fmriprep_dir = postproc_root,
    bids_dir = fix$bids_dir,
    sub = "01",
    cohort = "2",
    desc = "postproc"
  )
  expect_equal(nrow(sub_pp$voxel_map), 18L)
  expect_equal(ncol(sub_pp$nuisance$run1), 0L)
  expect_identical(sub_pp$confounds$settings$strategy, "none")
  expect_false(sub_pp$valid_rows$run1[fix$spike_vol[1L]])
  expect_equal(length(sub_pp$bids$censor_files), 2L)

  # Explicit strategy = "6p" on already-regressed BOLD emits informative orthogonalization warning
  expect_warning(
    sub_pp_warn <- read_fmriprep_subject(
      fmriprep_dir = postproc_root,
      sub = "01",
      cohort = "2",
      desc = "postproc",
      strategy = "6p",
      events = NULL
    ),
    "confound_regression"
  )
  expect_equal(ncol(sub_pp_warn$nuisance$run1), 6L)

  # Even if *_desc-postproc_confounds.tsv is absent in postproc_sub, upstream/censor-only ingestion succeeds!
  unlink(list.files(postproc_sub, pattern = "_confounds\\.tsv$", full.names = TRUE))
  sub_pp_noconf <- read_fmriprep_subject(
    fmriprep_dir = postproc_root,
    sub = "01",
    cohort = "2",
    desc = "postproc",
    events = NULL
  )
  expect_equal(ncol(sub_pp_noconf$nuisance$run1), 0L)
  expect_false(sub_pp_noconf$valid_rows$run1[fix$spike_vol[1L]])

  # Flip scrub_timepoints Status to "executed" in run 1 JSON -> MUST error immediately!
  r1_json <- file.path(postproc_sub, "sub-01_ses-01_task-mid_run-01_space-MNIPediatricAsym_cohort-2_desc-postproc_bold.json")
  j_lines <- readLines(r1_json)
  j_lines <- sub('"Name": "scrub_timepoints", "Status": "skipped"', '"Name": "scrub_timepoints", "Status": "executed"', j_lines, fixed = TRUE)
  writeLines(j_lines, r1_json)
  expect_error(
    read_fmriprep_subject(postproc_root, sub = "01", cohort = "2", desc = "postproc", events = NULL),
    "scrub_timepoints"
  )
})

test_that("Respiratory motion filtering for FD suppresses fast-TR breathing pseudomotion while preserving head-jerk censoring", {
  # Simulate fast-TR multiband acquisition (TR = 0.8 s, N = 120 volumes, Nyquist = 0.625 Hz = 37.5 BPM)
  # with in-phase 15 BPM (0.25 Hz) breathing oscillation in trans_y (0.30 mm) and rot_x (0.003 rad = 0.15 mm on 50 mm sphere)
  # plus a true sharp head-jerk step (+0.85 mm in trans_z) starting at volume 60.
  tr_fast <- 0.8
  n_vol <- 120L
  t_sec <- (seq_len(n_vol) - 1L) * tr_fast

  tx <- 0.01 * sin(2 * pi * 0.03 * t_sec)
  ty <- 0.30 * sin(2 * pi * 0.25 * t_sec)
  tz <- ifelse(seq_len(n_vol) >= 60L, 0.85, 0.0)
  rx <- 0.003 * sin(2 * pi * 0.25 * t_sec)
  ry <- 0.0002 * sin(2 * pi * 0.04 * t_sec)
  rz <- -0.0002 * cos(2 * pi * 0.02 * t_sec)

  mk_d1 <- function(v) c(NA_real_, diff(v))
  mk_dp2 <- function(v) c(NA_real_, diff(v)^2)

  mot_mat <- cbind(rx, ry, rz, tx, ty, tz)
  d_mot <- rbind(rep(0, 6L), diff(mot_mat))
  raw_fd <- rowSums(abs(d_mot[, 1:3, drop = FALSE])) * 50 + rowSums(abs(d_mot[, 4:6, drop = FALSE]))
  raw_fd_col <- c(NA_real_, raw_fd[-1L])

  resp_df <- data.frame(
    trans_x = tx, trans_x_derivative1 = mk_d1(tx), trans_x_power2 = tx^2, trans_x_derivative1_power2 = mk_dp2(tx),
    trans_y = ty, trans_y_derivative1 = mk_d1(ty), trans_y_power2 = ty^2, trans_y_derivative1_power2 = mk_dp2(ty),
    trans_z = tz, trans_z_derivative1 = mk_d1(tz), trans_z_power2 = tz^2, trans_z_derivative1_power2 = mk_dp2(tz),
    rot_x = rx, rot_x_derivative1 = mk_d1(rx), rot_x_power2 = rx^2, rot_x_derivative1_power2 = mk_dp2(rx),
    rot_y = ry, rot_y_derivative1 = mk_d1(ry), rot_y_power2 = ry^2, rot_y_derivative1_power2 = mk_dp2(ry),
    rot_z = rz, rot_z_derivative1 = mk_d1(rz), rot_z_power2 = rz^2, rot_z_derivative1_power2 = mk_dp2(rz),
    framewise_displacement = raw_fd_col
  )

  # 1. Unfiltered FD at 0.30 mm severely over-censors clean breathing frames
  conf_raw <- read_fmriprep_confounds(
    resp_df,
    strategy = "24p",
    fd_threshold = 0.30,
    include_non_steady_state = FALSE,
    impute_censored = FALSE,
    center = FALSE
  )
  expect_null(conf_raw$framewise_displacement_raw)
  expect_false(conf_raw$summary$motion_filtered)
  expect_gt(sum(!conf_raw$valid_rows$run1), 40L)
  expect_gt(conf_raw$summary$mean_fd, 0.25)

  # 2. Band-stop notch filter (12-18 BPM) suppresses 15 BPM breathing while still censoring the true jerk at vol 60!
  conf_notch <- read_fmriprep_confounds(
    resp_df,
    strategy = "24p",
    fd_threshold = 0.30,
    tr = tr_fast,
    motion_bandstop_bpm = c(12, 18),
    motion_filter_target = "fd_only",
    include_non_steady_state = FALSE,
    impute_censored = FALSE,
    center = FALSE
  )
  expect_true(conf_notch$summary$motion_filtered)
  expect_lt(conf_notch$summary$mean_fd, 0.05)
  expect_false(conf_notch$valid_rows$run1[60L]) # true step jerk at vol 60 is still censored!
  expect_true(all(conf_notch$valid_rows$run1[15:50])) # steady-state breathing volumes stay usable!
  expect_true(all(conf_notch$valid_rows$run1[70:110]))
  expect_equal(conf_notch$framewise_displacement_raw$run1, conf_raw$framewise_displacement$run1, tolerance = 1e-12)
  expect_match(paste(capture.output(print(conf_notch)), collapse = "\n"), "fd_filter: notch 12-18 BPM")

  # Default motion_filter_target = "fd_only" leaves all 24 Friston nuisance columns untouched!
  expect_equal(conf_notch$nuisance$run1, conf_raw$nuisance$run1, tolerance = 1e-12)

  # 3. motion_filter_target = "fd_and_nuisance" updates all 24 Friston motion columns (base + deriv + quadratics)
  conf_notch_nuis <- read_fmriprep_confounds(
    resp_df,
    strategy = "24p",
    fd_threshold = 0.30,
    tr = tr_fast,
    motion_bandstop_bpm = c(12, 18),
    motion_filter_target = "fd_and_nuisance",
    include_non_steady_state = FALSE,
    impute_censored = FALSE,
    center = FALSE
  )
  expect_identical(conf_notch_nuis$valid_rows$run1, conf_notch$valid_rows$run1)
  expect_equal(conf_notch_nuis$framewise_displacement$run1, conf_notch$framewise_displacement$run1, tolerance = 1e-12)
  expect_lt(stats::sd(conf_notch_nuis$nuisance$run1[20:50, "trans_y"]), 0.01)
  expect_lt(stats::sd(conf_notch_nuis$nuisance$run1[20:50, "trans_y_derivative1"]), 0.01)
  expect_lt(stats::sd(conf_notch_nuis$nuisance$run1[20:50, "trans_y_power2"]), 0.001)
  expect_lt(stats::sd(conf_notch_nuis$nuisance$run1[20:50, "trans_y_derivative1_power2"]), 0.001)

  # Ensure column_pattern on a 6-DOF-only table does not pick up surprise synthesized columns under fd_and_nuisance
  df_6dof <- resp_df[, c("trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z")]
  conf_pat <- read_fmriprep_confounds(
    df_6dof,
    strategy = "custom",
    column_pattern = "^trans_",
    tr = tr_fast,
    motion_bandstop_bpm = c(12, 18),
    motion_filter_target = "fd_and_nuisance",
    drop_constant = FALSE
  )
  expect_identical(colnames(conf_pat$nuisance$run1), c("trans_x", "trans_y", "trans_z"))

  # 4. Butterworth low-pass filter (0.20 Hz) also suppresses 0.25 Hz breathing while catching the step jerk at vol 60
  conf_lp <- read_fmriprep_confounds(
    resp_df,
    strategy = "6p",
    fd_threshold = 0.30,
    tr = tr_fast,
    motion_lowpass_hz = 0.20,
    include_non_steady_state = FALSE
  )
  expect_true(conf_lp$summary$motion_filtered)
  expect_lt(conf_lp$summary$mean_fd, 0.12)
  expect_false(conf_lp$valid_rows$run1[60L])
  expect_true(all(conf_lp$valid_rows$run1[15:50]))
  expect_match(paste(capture.output(print(conf_lp)), collapse = "\n"), "fd_filter: lowpass 0.2 Hz")

  # Order 8 cascaded notch + alternate padding modes + reversed bandstop bounds warning
  conf_ord8 <- read_fmriprep_confounds(
    resp_df, strategy = "none", fd_threshold = 0.30, tr = tr_fast,
    motion_bandstop_bpm = c(12, 18), motion_filter_order = 8L, motion_filter_pad = "odd",
    include_non_steady_state = FALSE
  )
  expect_true(conf_ord8$summary$motion_filtered)
  expect_warning(
    conf_rev <- read_fmriprep_confounds(
      resp_df, strategy = "none", fd_threshold = 0.30, tr = tr_fast,
      motion_bandstop_bpm = c(18, 12), include_non_steady_state = FALSE
    ),
    "reversed"
  )
  expect_equal(conf_rev$framewise_displacement$run1, conf_notch$framewise_displacement$run1, tolerance = 1e-12)
})

test_that("Respiratory motion filter handles TR inference, Nyquist edge cases, validation guards, and DSP parity", {
  # 1. Validation & Nyquist edge cases
  n_v <- 80L
  t2 <- (seq_len(n_v) - 1L) * 2.0 # TR = 2.0 s -> fs = 0.5 Hz, Nyquist = 0.25 Hz (15 BPM)
  # Signal at 0.35 Hz (21 BPM) aliases under TR = 2.0 s to |0.35 - 0.50| = 0.15 Hz (9 BPM)
  df_slow <- data.frame(
    trans_x = rep(0, n_v),
    trans_y = 0.25 * sin(2 * pi * 0.35 * t2),
    trans_z = rep(0, n_v),
    rot_x = rep(0, n_v),
    rot_y = rep(0, n_v),
    rot_z = rep(0, n_v)
  )

  # Missing TR on in-memory data.frame -> clear error
  expect_error(
    read_fmriprep_confounds(df_slow, motion_bandstop_bpm = c(12, 18)),
    "Sampling interval 'tr'"
  )
  # Mutual exclusivity of notch and lowpass
  expect_error(
    read_fmriprep_confounds(df_slow, tr = 0.8, motion_bandstop_bpm = c(12, 18), motion_lowpass_hz = 0.2),
    "at most one"
  )
  # Distinct bounds & valid filter orders
  expect_error(
    read_fmriprep_confounds(df_slow, tr = 0.8, motion_bandstop_bpm = c(15, 15)),
    "distinct"
  )
  expect_error(
    read_fmriprep_confounds(df_slow, tr = 0.8, motion_bandstop_bpm = c(12, 18), motion_filter_order = 6L),
    "multiple of 4"
  )
  expect_error(
    read_fmriprep_confounds(df_slow, tr = 0.8, motion_lowpass_hz = 0.2, motion_filter_order = 3L),
    "even integer"
  )

  # Straddling Nyquist (12-18 BPM across 15 BPM Nyquist at TR = 2.0 s) -> warns and falls back to unfiltered FD
  raw_slow <- read_fmriprep_confounds(df_slow, strategy = "none", fd_threshold = 0.30)
  expect_warning(
    strad <- read_fmriprep_confounds(df_slow, strategy = "none", fd_threshold = 0.30, tr = 2.0, motion_bandstop_bpm = c(12, 18)),
    "straddles the Nyquist frequency"
  )
  expect_false(strad$summary$motion_filtered)
  expect_equal(strad$framewise_displacement$run1, raw_slow$framewise_displacement$run1, tolerance = 1e-12)

  # Second Nyquist zone (18-24 BPM strictly inside (15, 30) BPM at TR = 2.0 s) -> folds to 6-12 BPM and attenuates 21 BPM!
  expect_warning(
    folded <- read_fmriprep_confounds(df_slow, strategy = "none", fd_threshold = 0.30, tr = 2.0, motion_bandstop_bpm = c(18, 24)),
    "folding to aliased stopband 6.00-12.00 BPM"
  )
  expect_true(folded$summary$motion_filtered)
  expect_lt(folded$summary$mean_fd, 0.05 * raw_slow$summary$mean_fd)

  # Low-pass cutoff >= Nyquist (0.30 Hz >= 0.25 Hz at TR = 2.0 s) -> warns and falls back to unfiltered FD
  expect_warning(
    lp_hi <- read_fmriprep_confounds(df_slow, strategy = "none", fd_threshold = 0.30, tr = 2.0, motion_lowpass_hz = 0.30),
    ">= Nyquist frequency"
  )
  expect_false(lp_hi$summary$motion_filtered)
  expect_equal(lp_hi$framewise_displacement$run1, raw_slow$framewise_displacement$run1, tolerance = 1e-12)

  # Pre-existing framewise_displacement_filtered column -> skips duplicate filtering with warning
  # and still records true unfiltered Power FD in framewise_displacement_raw!
  df_already <- df_slow
  df_already$framewise_displacement_filtered <- rep(0.05, n_v)
  expect_warning(
    dup_filt <- read_fmriprep_confounds(df_already, strategy = "none", tr = 0.8, motion_bandstop_bpm = c(12, 18)),
    "skipping duplicate respiratory motion filter"
  )
  expect_false(dup_filt$summary$motion_filtered)
  expect_equal(dup_filt$framewise_displacement$run1, rep(0.05, n_v))
  expect_equal(dup_filt$framewise_displacement_raw$run1, raw_slow$framewise_displacement$run1, tolerance = 1e-12)

  # Custom missing fd_col + Nyquist-straddling filter falls back to unfiltered Power FD instead of NA-censoring all rows
  expect_warning(
    expect_warning(
      fb_custom <- read_fmriprep_confounds(
        df_slow, strategy = "none", fd_col = "missing_custom_fd", fd_threshold = 0.60,
        tr = 2.0, motion_bandstop_bpm = c(12, 18)
      ),
      "recomputing Power"
    ),
    "straddles the Nyquist frequency"
  )
  expect_true(all(fb_custom$valid_rows$run1))
  expect_equal(fb_custom$framewise_displacement$run1, raw_slow$framewise_displacement$run1, tolerance = 1e-12)

  # Isolated NA in a drifted rigid-body column interpolates smoothly instead of injecting a step spike to 0
  t_fast <- (seq_len(n_v) - 1L) * 0.8
  df_drift_na <- data.frame(
    trans_x = rep(0, n_v),
    trans_y = 3.0 + 0.25 * sin(2 * pi * 0.25 * t_fast),
    trans_z = rep(0, n_v),
    rot_x = rep(0, n_v),
    rot_y = rep(0, n_v),
    rot_z = rep(0, n_v)
  )
  df_drift_na$trans_y[40L] <- NA_real_
  cf_drift_na <- read_fmriprep_confounds(
    df_drift_na, strategy = "none", fd_threshold = 0.30, tr = 0.8,
    motion_bandstop_bpm = c(12, 18), include_non_steady_state = FALSE
  )
  expect_true(all(cf_drift_na$valid_rows$run1[15:65]))

  # Explicit include_motion_outliers = TRUE warns when raw motion_outlier* re-censors filtered breathing frames
  df_with_mo <- df_drift_na
  df_with_mo$trans_y[40L] <- 3.0 + 0.25 * sin(2 * pi * 0.25 * t_fast[40L])
  df_with_mo$motion_outlier00 <- as.integer(seq_len(n_v) == 30L)
  expect_warning(
    cf_mo <- read_fmriprep_confounds(
      df_with_mo, strategy = "none", fd_threshold = 0.30, tr = 0.8,
      motion_bandstop_bpm = c(12, 18), include_motion_outliers = TRUE
    ),
    "include_motion_outliers = TRUE"
  )
  expect_false(cf_mo$valid_rows$run1[30L])

  # Automatic TR inference from a named list of BOLD matrices carrying attr(, "tr") and per-run TR vectors
  m_r1 <- matrix(0, nrow = n_v, ncol = 2L)
  attr(m_r1, "tr") <- 0.8
  cf_list_tr <- read_fmriprep_confounds(
    list(run1 = df_with_mo[, 1:6], run2 = df_with_mo[, 1:6]),
    bold = list(run1 = m_r1, run2 = m_r1),
    strategy = "none",
    motion_bandstop_bpm = c(12, 18)
  )
  expect_true(all(cf_list_tr$summary$motion_filtered))
  cf_vec_tr <- read_fmriprep_confounds(
    list(run1 = df_with_mo[, 1:6], run2 = df_with_mo[, 1:6]),
    tr = c(run1 = 0.8, run2 = 0.8),
    strategy = "none",
    motion_bandstop_bpm = c(12, 18)
  )
  expect_true(all(cf_vec_tr$summary$motion_filtered))

  # 2. Automatic TR inference from sibling *_bold.json and from read_fmriprep_subject(..., motion_bandstop_bpm = ...)
  fix <- make_synthetic_bids_fmriprep_tree(tempfile("fd_filt_subj_"), n_time = 40L, tr = 0.8, use_json_tr = TRUE)
  on.exit(unlink(dirname(fix$bids_dir), recursive = TRUE), add = TRUE)

  # Standalone confound file path infers RepetitionTime = 0.8 from sibling *_bold.json without 'bold' or 'tr'!
  conf_sidecar_tr <- read_fmriprep_confounds(fix$conf_files[1L], strategy = "6p", motion_bandstop_bpm = c(12, 18))
  expect_true(conf_sidecar_tr$summary$motion_filtered)

  sub_filt <- read_fmriprep_subject(
    fmriprep_dir = fix$fmriprep_dir,
    bids_dir = fix$bids_dir,
    sub = "01",
    task = "mid",
    strategy = "24p",
    motion_bandstop_bpm = c(12, 18),
    fd_threshold = 0.5
  )
  expect_true(all(sub_filt$confounds$summary$motion_filtered))
  expect_equal(length(sub_filt$confounds$framewise_displacement_raw), 2L)

  # 3. Exact numerical parity against signal::butter and BrainGnomes::filtfilt_cpp when installed
  if (requireNamespace("signal", quietly = TRUE)) {
    for (ord in 1:4) {
      bw_medusa <- medusa:::design_butterworth_lowpass(ord, cutoff_hz = 0.2, fs = 1.25)
      bw_sig <- signal::butter(ord, 0.2 / 0.625, type = "low")
      expect_equal(bw_medusa$b, as.double(bw_sig$b), tolerance = 1e-13)
      expect_equal(bw_medusa$a, as.double(bw_sig$a), tolerance = 1e-13)
    }
  }
  if (requireNamespace("BrainGnomes", quietly = TRUE)) {
    bg_ns <- asNamespace("BrainGnomes")
    if (exists("filtfilt_cpp", envir = bg_ns, inherits = FALSE)) {
      bg_filtfilt <- get("filtfilt_cpp", envir = bg_ns)
      set.seed(99L)
      x_test <- cumsum(rnorm(150L, sd = 0.05)) + 0.2 * sin(2 * pi * 0.25 * seq(0, by = 0.8, length.out = 150L))
      notch_c <- medusa:::design_iir_notch(f0 = 0.25, bandwidth = 0.10, fs = 1.25)
      for (pt in c("constant", "odd", "even", "zero")) {
        y_med <- medusa:::filtfilt_1d(x_test, b = notch_c$b, a = notch_c$a, padtype = pt, use_zi = TRUE)
        y_bg <- as.double(bg_filtfilt(x_test, b = notch_c$b, a = notch_c$a, padtype = pt, use_zi = TRUE))
        expect_equal(y_med, y_bg, tolerance = 1e-13)
      }
    }
  }
})
