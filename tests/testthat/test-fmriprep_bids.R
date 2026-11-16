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

test_that("read_fmriprep_aroma and apply_aroma support fMRIPost-AROMA, legacy AROMA, and BrainGnomes parity", {
  fix <- make_synthetic_bids_fmriprep_tree(tempfile("aroma_core_"), n_time = 40L, tr = 1.0)
  on.exit(unlink(dirname(fix$bids_dir), recursive = TRUE), add = TRUE)

  ts <- read_nifti_ts(
    niftis = c(run1 = fix$bold_files[1L], run2 = fix$bold_files[2L]),
    mask = fix$mask_files[1L],
    tr = fix$tr
  )

  # 1. fMRIPost-AROMA sidecars: headerless res-2 space-MNI152NLin6Asym mixing + space-less aroma_metrics
  #    matched cleanly against space-MNI152NLin2009cAsym BOLD even when passed in reversed run order!
  aroma_dir <- file.path(dirname(fix$bids_dir), "aroma_files")
  dir.create(aroma_dir)

  set.seed(20261116L)
  mix_files <- met_files <- character(2L)
  for (r in 1:2) {
    t_seq <- seq_len(40L)
    s_acc <- sin(2 * pi * t_seq / 10) + rnorm(40L, sd = 0.15)
    m_rej1 <- 0.55 * s_acc + rnorm(40L, sd = 0.5) + 1.5 # nonzero mean + task-correlated motion!
    m_acc2 <- cos(2 * pi * t_seq / 14) + rnorm(40L, sd = 0.2)
    m_rej2 <- rnorm(40L, sd = 0.7) - 0.8
    mix_r <- cbind(s_acc, m_rej1, m_acc2, m_rej2)

    mix_files[r] <- file.path(
      aroma_dir,
      sprintf("sub-01_ses-01_task-mid_run-0%d_space-MNI152NLin6Asym_res-2_desc-melodic_mixing.tsv", r)
    )
    met_files[r] <- file.path(
      aroma_dir,
      sprintf("sub-01_ses-01_task-mid_run-0%d_desc-aroma_metrics.tsv", r)
    )
    utils::write.table(mix_r, mix_files[r], sep = "\t", row.names = FALSE, col.names = FALSE)
    met_df <- data.frame(
      IC = 1:4,
      classification = if (r == 1L) c("accepted", "rejected", "accepted", "rejected") else c("accepted", "rejected", "accepted", "accepted"),
      max_rp_corr = c(0.12, 0.84, 0.09, 0.71),
      stringsAsFactors = FALSE
    )
    utils::write.table(met_df, met_files[r], sep = "\t", row.names = FALSE, quote = FALSE)
  }

  # Auto-discovers sibling *_desc-aroma_metrics.tsv when classification = NULL, and reorders reversed run-02, run-01!
  ar_auto <- read_fmriprep_aroma(mixing = rev(mix_files), bold = ts)
  expect_s3_class(ar_auto, "medusa_aroma")
  expect_identical(names(ar_auto$mixing), c("run1", "run2"))
  expect_identical(ar_auto$noise_ics$run1, c(2L, 4L))
  expect_identical(ar_auto$signal_ics$run1, c(1L, 3L))
  expect_identical(ar_auto$noise_ics$run2, 2L)
  expect_identical(colnames(ar_auto$mixing$run1), c("IC_01", "IC_02", "IC_03", "IC_04"))
  expect_match(paste(capture.output(print(ar_auto)), collapse = "\n"), "3/8 total component")

  # Headered mixing TSV + legacy *_AROMAnoiseICs.csv (including empty 0-rejected CSV!)
  leg_dir <- file.path(dirname(fix$bids_dir), "legacy_aroma")
  dir.create(leg_dir)
  leg_mix1 <- file.path(leg_dir, "sub-01_ses-01_task-mid_run-01_desc-MELODIC_mixing.tsv")
  leg_csv1 <- file.path(leg_dir, "sub-01_ses-01_task-mid_run-01_AROMAnoiseICs.csv")
  leg_mix2 <- file.path(leg_dir, "sub-01_ses-01_task-mid_run-02_desc-MELODIC_mixing.tsv")
  leg_csv2 <- file.path(leg_dir, "sub-01_ses-01_task-mid_run-02_AROMAnoiseICs.csv")

  hdr_df1 <- as.data.frame(ar_auto$mixing$run1)
  names(hdr_df1) <- paste0("melodic_", 1:4)
  utils::write.table(hdr_df1, leg_mix1, sep = "\t", row.names = FALSE, quote = FALSE)
  writeLines("2, 4", leg_csv1)
  utils::write.table(ar_auto$mixing$run2, leg_mix2, sep = "\t", row.names = FALSE, col.names = FALSE)
  writeLines("", leg_csv2) # 0 rejected ICs in run 2!

  ar_leg <- read_fmriprep_aroma(c(leg_mix1, leg_mix2), bold = ts)
  expect_identical(colnames(ar_leg$mixing$run1), paste0("melodic_", 1:4))
  expect_identical(ar_leg$noise_ics$run1, c(2L, 4L))
  expect_identical(ar_leg$noise_ics$run2, integer(0L))

  # Applying AROMA when run2 has 0 rejected ICs leaves run2 bit-for-bit unchanged
  ts_leg <- apply_aroma(ts, ar_leg)
  expect_true(ts_leg$aroma$summary$applied[1L])
  expect_false(ts_leg$aroma$summary$applied[2L])
  expect_equal(ts_leg$bold$run2, ts$bold$run2, tolerance = 0)
  expect_match(paste(capture.output(print(ts_leg)), collapse = "\n"), "AROMA denoising: nonaggressive")

  # Duplicate apply_aroma guard warns and skips unless force = TRUE
  expect_warning(ts_dup <- apply_aroma(ts_leg, ar_leg), "already denoised")
  expect_equal(ts_dup$bold$run1, ts_leg$bold$run1, tolerance = 0)

  # Validation guards: out-of-range 0-based index, wrong volume count, and insufficient df
  expect_error(read_fmriprep_aroma(ar_auto$mixing$run1, classification = c(0L, 2L)), "1-based integers")
  expect_error(read_fmriprep_aroma(ar_auto$mixing$run1, classification = 9L), "1-based integers")
  expect_error(read_fmriprep_aroma(ar_auto$mixing$run1[1:10, ], classification = 2L, bold = ts$bold$run1), "has 10 rows")
  expect_error(
    apply_aroma(matrix(rnorm(12L), nrow = 6L, ncol = 2L), mixing = matrix(rnorm(30L), nrow = 6L, ncol = 5L), classification = 1L),
    "Not enough fitted volumes"
  )

  # 2. Exact numerical parity against lm.fit partial regression and BrainGnomes:::lmfit_residuals_mat
  Y_test <- ts$bold$run1
  M_test <- cbind(ar_auto$mixing$run1, const_ic = 5.0) # include constant IC column!
  valid_sub <- rep(TRUE, 40L)
  valid_sub[c(1L, 12L, 28L)] <- FALSE
  rej_test <- c(2L, 4L, 5L)

  lm_ref_aroma <- function(Y, M, rej_idx, inc, nonagg = TRUE) {
    idx <- which(inc)
    rng <- apply(M[idx, , drop = FALSE], 2L, function(v) abs(max(v) - min(v)))
    keep <- which(rng >= 1e-8)
    M_k <- M[, keep, drop = FALSE]
    rej_k <- which(keep %in% rej_idx)
    X_all <- cbind(1, M_k)
    fit_cols <- if (nonagg) seq_len(ncol(X_all)) else c(1L, 1L + rej_k)
    reg_cols <- c(1L, 1L + rej_k)
    pred_pos <- match(reg_cols, fit_cols)
    out <- Y
    for (v in seq_len(ncol(Y))) {
      y <- Y[, v]
      cf <- stats::lm.fit(X_all[idx, fit_cols, drop = FALSE], y[idx])$coefficients
      cf[is.na(cf)] <- 0
      res <- y - as.vector(X_all[, reg_cols, drop = FALSE] %*% cf[pred_pos])
      out[, v] <- res - mean(res[idx]) + mean(y[idx])
    }
    out
  }

  y_na_all <- apply_aroma(Y_test, mixing = M_test, classification = rej_test, mode = "nonaggressive")
  expect_equal(as.matrix(y_na_all), lm_ref_aroma(Y_test, M_test, rej_test, rep(TRUE, 40L), TRUE),
               tolerance = 1e-12, ignore_attr = TRUE)

  y_na_val <- apply_aroma(Y_test, mixing = M_test, classification = rej_test, mode = "nonaggressive",
                          fit_valid_rows = TRUE, valid_rows = valid_sub)
  expect_equal(as.matrix(y_na_val), lm_ref_aroma(Y_test, M_test, rej_test, valid_sub, TRUE),
               tolerance = 1e-12, ignore_attr = TRUE)

  y_ag_all <- apply_aroma(Y_test, mixing = M_test, classification = rej_test, mode = "aggressive")
  expect_equal(as.matrix(y_ag_all), lm_ref_aroma(Y_test, M_test, rej_test, rep(TRUE, 40L), FALSE),
               tolerance = 1e-12, ignore_attr = TRUE)

  # Voxel-by-time medusa_nifti_ts orientation + drop = TRUE single-matrix attributes
  ts_vt <- read_nifti_ts(c(run1 = fix$bold_files[1L]), mask = fix$mask_files[1L], tr = fix$tr, orientation = "voxel_by_time")
  ts_vt_clean <- apply_aroma(ts_vt, mixing = M_test, classification = rej_test)
  expect_equal(t(ts_vt_clean$bold$run1), as.matrix(y_na_all), tolerance = 1e-12, ignore_attr = TRUE)

  mat_drop <- read_nifti_ts(fix$bold_files[1L], mask = fix$mask_files[1L], tr = fix$tr, drop = TRUE)
  mat_drop_clean <- apply_aroma(mat_drop, mixing = M_test, classification = rej_test)
  expect_identical(attr(mat_drop_clean, "voxel_map"), attr(mat_drop, "voxel_map"))
  expect_identical(attr(mat_drop_clean, "tr"), attr(mat_drop, "tr"))

  if (requireNamespace("BrainGnomes", quietly = TRUE)) {
    bg_ns <- asNamespace("BrainGnomes")
    if (exists("lmfit_residuals_mat", envir = bg_ns, inherits = FALSE)) {
      bg_lm <- get("lmfit_residuals_mat", envir = bg_ns)
      # Center MELODIC columns on valid_sub so older installed BrainGnomes builds (pre-offset-recenter fix) also match
      M_bg <- sweep(M_test, 2L, colMeans(M_test[valid_sub, , drop = FALSE]), "-")
      M_bg[, 5L] <- 5.0
      y_bg <- suppressWarnings(bg_lm(
        Y = Y_test, X = M_bg, include_rows = valid_sub,
        add_intercept = TRUE, preserve_mean = TRUE,
        regress_cols = rej_test, exclusive = FALSE
      ))
      expect_equal(as.matrix(y_na_val), y_bg, tolerance = 1e-12, ignore_attr = TRUE)
    }
  }
  expect_equal(colMeans(y_na_val[valid_sub, , drop = FALSE]), colMeans(Y_test[valid_sub, , drop = FALSE]), tolerance = 1e-12)
})

test_that("read_fmriprep_confounds and read_fmriprep_subject coordinate AROMA with nuisance regressors and BrainGnomes provenance", {
  fix <- make_synthetic_bids_fmriprep_tree(tempfile("aroma_subj_"), n_time = 36L, tr = 1.0)
  on.exit(unlink(dirname(fix$bids_dir), recursive = TRUE), add = TRUE)

  # Populate sibling derivatives/fmripost-aroma/sub-01/ses-01/func/ to test automatic sibling discovery!
  deriv_root <- dirname(fix$fmriprep_dir)
  aroma_func <- file.path(deriv_root, "fmripost-aroma", "sub-01", "ses-01", "func")
  dir.create(aroma_func, recursive = TRUE)

  set.seed(77L)
  for (r in 1:2) {
    c_df <- utils::read.delim(fix$conf_files[r], stringsAsFactors = FALSE)
    s_sig <- sin(2 * pi * seq_len(36L) / 9)
    m_mot <- 0.5 * s_sig + rnorm(36L, sd = 0.6)
    mix_r <- cbind(s_sig, m_mot, rnorm(36L))

    # Inject m_mot into csf and make trans_x purely proportional to m_mot (should drop after AROMA cleaning!)
    c_df$csf <- 800 + 2.5 * m_mot + rnorm(36L, sd = 0.2)
    c_df$white_matter <- 650 - 1.2 * m_mot + rnorm(36L, sd = 0.2)
    c_df$trans_x <- 0.35 * m_mot
    utils::write.table(c_df, fix$conf_files[r], sep = "\t", row.names = FALSE, quote = FALSE)

    stem_r <- sprintf("sub-01_ses-01_task-mid_run-0%d", r)
    utils::write.table(
      mix_r,
      file.path(aroma_func, paste0(stem_r, "_res-2_desc-melodic_mixing.tsv")),
      sep = "\t", row.names = FALSE, col.names = FALSE
    )
    utils::write.table(
      data.frame(IC = 1:3, classification = c("accepted", "rejected", "accepted")),
      file.path(aroma_func, paste0(stem_r, "_desc-aroma_metrics.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }

  # 1. Verify csf_wm / csf_wm_8p strategies and exact commutativity of confound AROMA cleaning under censoring
  ts_raw <- read_nifti_ts(c(run1 = fix$bold_files[1L], run2 = fix$bold_files[2L]), mask = fix$mask_files[1L], tr = fix$tr)
  ar_obj <- read_fmriprep_aroma(
    list.files(aroma_func, pattern = "_desc-melodic_mixing\\.tsv$", full.names = TRUE),
    bold = ts_raw
  )

  conf_pre <- read_fmriprep_confounds(fix$conf_files, bold = ts_raw, strategy = "csf_wm", fd_threshold = 0.5)
  expect_identical(colnames(conf_pre$nuisance$run1), c("csf", "white_matter"))
  conf_8p <- read_fmriprep_confounds(fix$conf_files, bold = ts_raw, strategy = "csf_wm_8p", fd_threshold = 0.5)
  expect_equal(ncol(conf_8p$nuisance$run1), 8L)

  conf_direct <- read_fmriprep_confounds(fix$conf_files, bold = ts_raw, strategy = "csf_wm", fd_threshold = 0.5, aroma = ar_obj)
  conf_post <- apply_aroma(conf_pre, ar_obj)
  expect_equal(conf_direct$nuisance, conf_post$nuisance, tolerance = 1e-13)
  expect_identical(conf_direct$valid_rows, conf_pre$valid_rows)
  expect_equal(conf_direct$framewise_displacement, conf_pre$framewise_displacement, tolerance = 0)
  expect_match(paste(capture.output(print(conf_direct)), collapse = "\n"), "aroma: nonaggressive")

  # trans_x was purely collinear with rejected IC 2 -> automatically dropped as zero-variance after AROMA cleaning!
  conf_6p_ar <- read_fmriprep_confounds(fix$conf_files, bold = ts_raw, strategy = "6p", fd_threshold = NULL,
                                        include_non_steady_state = FALSE, aroma = ar_obj)
  expect_false("trans_x" %in% colnames(conf_6p_ar$nuisance$run1))

  # 2. Turnkey read_fmriprep_subject auto-discovers sibling derivatives/fmripost-aroma without aroma_dir!
  sub_ar <- read_fmriprep_subject(
    fmriprep_dir = fix$fmriprep_dir,
    bids_dir = fix$bids_dir,
    sub = "01",
    task = "mid",
    aroma = "nonaggressive",
    strategy = "csf_wm"
  )
  expect_true(sub_ar$aroma$applied)
  expect_true(sub_ar$aroma$nuisance_cleaned)
  expect_equal(length(sub_ar$bids$aroma_mixing_files), 2L)
  expect_equal(sub_ar$nuisance, conf_direct$nuisance, tolerance = 1e-13)

  # Remove sibling fmripost-aroma: aroma = "auto" succeeds as graceful no-op, while aroma = "nonaggressive" errors clearly
  unlink(file.path(deriv_root, "fmripost-aroma"), recursive = TRUE)
  sub_auto_noop <- read_fmriprep_subject(fix$fmriprep_dir, bids_dir = fix$bids_dir, sub = "01", task = "mid", aroma = "auto", strategy = "csf_wm")
  expect_null(sub_auto_noop$aroma)

  # Add legacy aroma_motion_01 column to confounds TSV -> error includes helpful hint about full MELODIC mixing matrix
  c1 <- utils::read.delim(fix$conf_files[1L], stringsAsFactors = FALSE)
  c1$aroma_motion_01 <- rnorm(36L)
  utils::write.table(c1, fix$conf_files[1L], sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(
    read_fmriprep_subject(fix$fmriprep_dir, bids_dir = fix$bids_dir, sub = "01", task = "mid", aroma = "nonaggressive"),
    "aroma_motion_"
  )

  # 3. BrainGnomes provenance guard when 'apply_aroma' was already executed upstream
  pp_dir <- file.path(deriv_root, "postproc", "sub-01")
  dir.create(pp_dir, recursive = TRUE)
  aroma_ext <- file.path(deriv_root, "custom_aroma")
  dir.create(aroma_ext, recursive = TRUE)

  for (r in 1:2) {
    stem_pp <- sprintf("sub-01_ses-01_task-mid_run-0%d_space-MNI152NLin2009cAsym", r)
    pp_bold <- file.path(pp_dir, paste0(stem_pp, "_desc-postproc_bold.nii.gz"))
    file.copy(fix$bold_files[r], pp_bold)
    file.copy(fix$mask_files[r], file.path(pp_dir, paste0(stem_pp, "_desc-brain_mask.nii.gz")))

    mix_p <- file.path(aroma_ext, sprintf("sub-01_ses-01_task-mid_run-0%d_desc-melodic_mixing.tsv", r))
    met_p <- file.path(aroma_ext, sprintf("sub-01_ses-01_task-mid_run-0%d_desc-aroma_metrics.tsv", r))
    utils::write.table(ar_obj$mixing[[r]], mix_p, sep = "\t", row.names = FALSE, col.names = FALSE)
    utils::write.table(ar_obj$metrics[[r]], met_p, sep = "\t", row.names = FALSE, quote = FALSE)

    writeLines(c(
      "{",
      '  "BrainGnomes": {',
      '    "Operations": [ { "Name": "apply_aroma", "Status": "executed" } ],',
      '    "SourceFiles": [',
      sprintf('      { "role": "confounds", "path": "%s" },', fix$conf_files[r]),
      sprintf('      { "role": "melodic_mix", "path": "%s" },', mix_p),
      sprintf('      { "role": "aroma_metrics", "path": "%s" }', met_p),
      "    ]",
      "  }",
      "}"
    ), sub("\\.nii\\.gz$", ".json", pp_bold))
  }

  # Calling aroma = "none" with raw upstream confounds + strategy = "csf_wm" warns about potential artifact re-introduction
  expect_warning(
    read_fmriprep_subject(file.path(deriv_root, "postproc"), sub = "01", desc = "postproc",
                          aroma = "none", strategy = "csf_wm", events = NULL),
    "re-introduce AROMA motion variance"
  )

  # Calling aroma = "auto" skips duplicate BOLD AROMA while cleaning the raw upstream confounds via provenance SourceFiles!
  sub_bg_auto <- read_fmriprep_subject(
    file.path(deriv_root, "postproc"), sub = "01", desc = "postproc",
    aroma = "auto", strategy = "csf_wm", events = NULL
  )
  expect_null(sub_bg_auto$aroma) # BOLD was NOT double-regressed
  expect_equal(sub_bg_auto$bold, ts_raw$bold, tolerance = 0)
  expect_true(isTRUE(sub_bg_auto$confounds$settings$aroma$applied)) # raw confounds WERE cleaned!

  # Calling aroma = "nonaggressive" explicitly warns that BOLD already had apply_aroma executed
  expect_warning(
    read_fmriprep_subject(file.path(deriv_root, "postproc"), sub = "01", desc = "postproc",
                          aroma = "nonaggressive", strategy = "csf_wm", events = NULL),
    "skipping duplicate AROMA regression on BOLD"
  )
})

test_that("End-to-end MEDuSA Steps 1-3 with non-aggressive AROMA preserves task-correlated neural signal and blocks artifact re-introduction", {
  set.seed(20261117L)
  root <- tempfile("medusa_aroma_e2e_")
  fprep_func <- file.path(root, "derivatives", "fmriprep", "sub-01", "func")
  aroma_func <- file.path(root, "derivatives", "fmripost-aroma", "sub-01", "func")
  bids_func <- file.path(root, "bids", "sub-01", "func")
  dir.create(fprep_func, recursive = TRUE)
  dir.create(aroma_func, recursive = TRUE)
  dir.create(bids_func, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  tr <- 1.0
  n_time <- 96L
  dim3 <- c(3L, 2L, 2L) # 12 spatial voxels, 6 in-mask parcel voxels
  mask_arr <- array(0L, dim = dim3)
  mask_arr[, , 1L] <- 1L
  atlas_arr <- mask_arr

  affine_e2e <- diag(c(2, 2, 2, 1))
  attr(affine_e2e, "code") <- 4L
  write_nii_local <- function(arr, path, img_tr = tr) {
    nii <- RNifti::asNifti(arr)
    RNifti::pixdim(nii) <- if (length(dim(arr)) == 4L) c(2, 2, 2, img_tr) else c(2, 2, 2)
    RNifti::sform(nii) <- affine_e2e
    RNifti::writeNifti(nii, path)
    path
  }
  atlas_path <- write_nii_local(atlas_arr, file.path(root, "atlas.nii.gz"))
  mask_template <- write_nii_local(mask_arr, file.path(root, "mask.nii.gz"))

  k_lib <- hrf_library(tr, peak_delays = c(5.0, 6.5), dispersions = 1.0, undershoot_ratios = 6)
  h_true <- as.numeric(k_lib[, 1L])
  true_beta <- 0.60

  for (r in 1:3) {
    onsets <- seq(6, 83, by = 7) # 12 feedback events per run
    rpe <- round(as.numeric(scale(sin(seq_along(onsets) * 1.3 + r))), 3L)
    ev_df <- data.frame(
      onset = onsets,
      duration = rep(0, length(onsets)),
      trial_type = rep("feedback", length(onsets)),
      rpe = rpe
    )

    # True neural stick train: mean feedback response (1.0) + RPE modulation (0.60 * rpe)
    neural <- numeric(n_time)
    idx_on <- as.integer(round(onsets / tr)) + 1L
    neural[idx_on] <- 1.0 + true_beta * rpe
    conv_full <- stats::convolve(neural, rev(h_true), type = "open")
    s_task <- conv_full[seq_len(n_time)]

    # Task-correlated motion component (r ~ 0.65 with task BOLD!) + independent motion component
    m_rej1 <- 0.85 * as.numeric(scale(s_task)) + rnorm(n_time, sd = 0.65)
    m_rej2 <- sin(2 * pi * seq_len(n_time) / 5.5) + rnorm(n_time, sd = 0.5)
    s_bg <- cos(2 * pi * seq_len(n_time) / 18)

    # True physiological CSF/WM noise (shared with BOLD) + contamination from motion ICs
    u_csf <- rnorm(n_time, sd = 0.35)
    u_wm <- rnorm(n_time, sd = 0.35)
    csf_raw <- 500 + 0.8 * m_rej1 + u_csf
    wm_raw <- 400 + 0.8 * m_rej2 + u_wm

    # Voxel BOLD = baseline + true task signal + background neural IC + heavy motion ICs + physiological + low white noise
    vox_ts <- matrix(0, nrow = n_time, ncol = 6L)
    mot_w1 <- c(0.9, 1.1, 0.8, 1.2, 1.0, 0.85)
    mot_w2 <- c(-0.6, 0.7, -0.5, 0.6, -0.7, 0.55)
    for (v in 1:6) {
      vox_ts[, v] <- 100 + s_task + 0.03 * s_bg + mot_w1[v] * m_rej1 + mot_w2[v] * m_rej2 +
        0.8 * u_csf + 0.8 * u_wm + rnorm(n_time, sd = 0.03)
    }

    arr4d <- array(0, dim = c(dim3, n_time))
    for (t_i in seq_len(n_time)) arr4d[, , 1L, t_i] <- matrix(vox_ts[t_i, ], nrow = 3L, ncol = 2L)

    stem <- sprintf("sub-01_task-mid_run-0%d", r)
    write_nii_local(arr4d, file.path(fprep_func, paste0(stem, "_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz")))
    file.copy(mask_template, file.path(fprep_func, paste0(stem, "_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz")))
    utils::write.table(ev_df, file.path(bids_func, paste0(stem, "_events.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

    conf_df <- data.frame(
      trans_x = rnorm(n_time, sd = 0.01), trans_y = rnorm(n_time, sd = 0.01), trans_z = rnorm(n_time, sd = 0.01),
      rot_x = rnorm(n_time, sd = 0.001), rot_y = rnorm(n_time, sd = 0.001), rot_z = rnorm(n_time, sd = 0.001),
      csf = csf_raw, white_matter = wm_raw,
      framewise_displacement = c("n/a", rep("0.08", n_time - 1L)),
      stringsAsFactors = FALSE
    )
    utils::write.table(conf_df, file.path(fprep_func, paste0(stem, "_desc-confounds_timeseries.tsv")),
                       sep = "\t", row.names = FALSE, quote = FALSE)

    # MELODIC mixing matrix: IC 1 = task neural, IC 2 = background neural, IC 3 & 4 = rejected motion ICs
    mix_mat <- cbind(s_task, s_bg, m_rej1, m_rej2)
    utils::write.table(mix_mat, file.path(aroma_func, paste0(stem, "_res-2_desc-melodic_mixing.tsv")),
                       sep = "\t", row.names = FALSE, col.names = FALSE)
    utils::write.table(
      data.frame(IC = 1:4, classification = c("accepted", "accepted", "rejected", "rejected")),
      file.path(aroma_func, paste0(stem, "_desc-aroma_metrics.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }

  run_pipeline <- function(aroma_opt, clean_nuis = TRUE, strat = "csf_wm") {
    sub_ts <- read_fmriprep_subject(
      fmriprep_dir = file.path(root, "derivatives", "fmriprep"),
      bids_dir = file.path(root, "bids"),
      sub = "01",
      task = "mid",
      atlas = atlas_path,
      aroma = aroma_opt,
      aroma_clean_nuisance = clean_nuis,
      strategy = strat,
      fd_threshold = 0.5
    )
    h_fit <- estimate_crossrun_hrf(
      sub_ts,
      kernel_grid = k_lib,
      pooling = "parcel",
      amplitude_model = "trial",
      noise_model = "ar1",
      evidence = "calibrated"
    )
    dec <- deconvolve_crossrun_hrf(
      sub_ts,
      h_fit,
      lambda = 0.05,
      scaling = "noise_sd",
      nuisance_mode = if (identical(strat, "none")) "ignore" else "joint"
    )
    ep <- decon_epoch_regression(
      dec,
      covariates = "rpe",
      target_condition = "feedback",
      lags = 0,
      aggregate_by_parcel = TRUE,
      restore_scale = TRUE
    )
    cf0 <- subset(ep$coefficients, term == "rpe" & abs(lag) < 1e-9)
    list(
      noise_sd = median(unlist(lapply(dec$scaling, function(s) s$noise_sd))),
      beta = cf0$beta_corrected[1L],
      se = cf0$se_weighted[1L]
    )
  }

  res_raw <- run_pipeline("none", strat = "none")
  res_agg <- run_pipeline("aggressive", clean_nuis = TRUE, strat = "csf_wm")
  res_unclean_nuis <- run_pipeline("nonaggressive", clean_nuis = FALSE, strat = "csf_wm")
  res_nonagg <- run_pipeline("nonaggressive", clean_nuis = TRUE, strat = "csf_wm")

  # 1. Non-aggressive AROMA + cleaned CSF/WM dramatically lowers Step-2 observation noise SD vs un-denoised BOLD
  expect_lt(res_nonagg$noise_sd, 0.15)
  expect_gt(res_raw$noise_sd, 3.0 * res_nonagg$noise_sd)

  # 2. Aggressive AROMA strips task-correlated motion shared with IC 1, severely attenuating the lag-0 RPE slope!
  expect_lt(res_agg$beta, 0.65 * true_beta)

  # 3. Non-aggressive AROMA preserves the true lag-0 RPE slope accurately with much smaller SE than un-denoised BOLD
  expect_lt(abs(res_nonagg$beta - true_beta), 0.10)
  expect_gt(res_nonagg$beta, 1.4 * res_agg$beta)
  expect_lt(res_nonagg$se, 0.5 * res_raw$se)

  # 4. Lindquist et al. (2019) check: failing to AROMA-clean CSF/WM re-introduces m_rej1 during joint regression!
  expect_gt(res_unclean_nuis$noise_sd, 2.0 * res_nonagg$noise_sd)
})
