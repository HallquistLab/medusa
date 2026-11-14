test_that("matrix-native alignment matches legacy fmri_ts oracle across groupings, kernels, gaps, and censors", {
  set.seed(20261113)
  # Include a deleted acquisition timepoint (gap between t=14 and t=18)
  times <- c(0:14, 18:40)
  n_t <- length(times)
  n_v <- 12L
  parcels <- rep(c(10, 20, 30), each = 4L)

  y_mat <- matrix(rnorm(n_t * n_v), nrow = n_t, ncol = n_v)
  colnames(y_mat) <- sprintf("v%d", seq_len(n_v))
  # Inject explicit NA values in specific voxels/times
  y_mat[5, 2] <- NA_real_
  y_mat[16, 1:4] <- NA_real_ # all units in parcel 10 missing at one sample
  y_mat[, 12] <- NA_real_    # completely missing voxel in parcel 30

  valid_vec <- rep(TRUE, n_t)
  valid_vec[9] <- FALSE
  valid_vec[20] <- NA

  # Construct matching long data.frame for fmri_ts
  long_df <- expand.grid(idx = seq_len(n_t), vnum = seq_len(n_v))
  long_df$time <- times[long_df$idx]
  long_df$atlas_value <- parcels[long_df$vnum]
  long_df$decon <- y_mat[cbind(long_df$idx, long_df$vnum)]
  long_df$usable <- valid_vec[long_df$idx]
  long_df$idx <- NULL

  ev_df <- data.frame(
    id = "sub01",
    run = 1L,
    trial = c(4L, 1L, 3L, 2L, 5L, 6L),
    onset = c(15.5, 4.3, 9.2, 8.6, NA_real_, 95.0),
    prev_bound = c(12.0, NA_real_, 8.8, Inf, 1.0, 90.0),
    next_bound = c(22.0, 8.5, 13.5, -Inf, 5.0, 99.0)
  )

  f_obj <- fmri_ts$new(
    ts_data = long_df,
    event_data = ev_df,
    tr = 1.0,
    vm = list(value = "decon", key = c("vnum", "atlas_value"))
  )

  unit_df <- data.frame(vnum = seq_len(n_v), atlas_value = parcels)

  for (gb in list("trial", c("trial", "atlas_value"), c("atlas_value", "trial"), "atlas_value", NULL)) {
    for (pol in c("allow", "truncate", "exclude")) {
      for (rk in c("linear", "boxcar", "triangle", "gaussian")) {
        for (mg in c(1.0, 5.0)) {
          oracle <- get_medusa_interpolated_ts(
            f_obj, event = "onset",
            time_before = -2.0, time_after = 3.0,
            collide_before = "prev_bound", collide_after = "next_bound",
            pad_before = -1.5, pad_after = 1.5,
            output_resolution = 0.5,
            group_by = gb,
            time_audit = TRUE,
            overlap_policy = pol,
            valid_col = "usable",
            neighbor_events = c("onset", "prev_bound"),
            max_gap = mg,
            readout_kernel = rk,
            use_legacy = TRUE
          )

          fast_fmri <- get_medusa_interpolated_ts(
            f_obj, event = "onset",
            time_before = -2.0, time_after = 3.0,
            collide_before = "prev_bound", collide_after = "next_bound",
            pad_before = -1.5, pad_after = 1.5,
            output_resolution = 0.5,
            group_by = gb,
            time_audit = TRUE,
            overlap_policy = pol,
            valid_col = "usable",
            neighbor_events = c("onset", "prev_bound"),
            max_gap = mg,
            readout_kernel = rk,
            use_legacy = FALSE
          )

          mat_native <- align_medusa_ts(
            y_mat,
            events = ev_df,
            event = "onset",
            time = times,
            tr = 1.0,
            unit_data = unit_df,
            group_by = gb,
            value_name = "decon",
            valid = valid_vec,
            time_before = -2.0, time_after = 3.0,
            collide_before = "prev_bound", collide_after = "next_bound",
            pad_before = -1.5, pad_after = 1.5,
            output_resolution = 0.5,
            time_audit = TRUE,
            overlap_policy = pol,
            neighbor_events = c("onset", "prev_bound"),
            max_gap = mg,
            readout_kernel = rk
          )

          expect_equal(fast_fmri, oracle, tolerance = 1e-12)
          expect_equal(attr(fast_fmri, "source_clock"), attr(oracle, "source_clock"))
          expect_equal(mat_native, oracle, tolerance = 1e-12)
          expect_equal(attr(mat_native, "source_clock"), attr(oracle, "source_clock"))
        }
      }
    }
  }
})

test_that("multivariate signals and voxel-by-time censor matrices match legacy oracle", {
  set.seed(20261114)
  n_t <- 25L
  n_v <- 6L
  times <- (seq_len(n_t) - 1) * 1.5
  sig_a <- matrix(rnorm(n_t * n_v), nrow = n_t, ncol = n_v)
  sig_b <- matrix(rnorm(n_t * n_v), nrow = n_t, ncol = n_v)
  sig_b[, 1:3] <- NA_real_

  valid_mat <- matrix(sample(c(TRUE, FALSE), n_t * n_v, replace = TRUE, prob = c(0.85, 0.15)),
                      nrow = n_t, ncol = n_v)
  parcels <- rep(c("ROI_1", "ROI_2"), each = 3L)

  long_df <- expand.grid(idx = seq_len(n_t), vnum = seq_len(n_v))
  long_df$time <- times[long_df$idx]
  long_df$region <- parcels[long_df$vnum]
  long_df$decon <- sig_a[cbind(long_df$idx, long_df$vnum)]
  long_df$bold <- sig_b[cbind(long_df$idx, long_df$vnum)]
  long_df$ok <- valid_mat[cbind(long_df$idx, long_df$vnum)]
  long_df$idx <- NULL

  ev_df <- data.frame(
    id = "S1", run = 1L,
    trial = c("t2", "t1", "t3"),
    onset = c(14.2, 6.1, 27.8)
  )

  f_obj <- fmri_ts$new(
    ts_data = long_df,
    event_data = ev_df,
    tr = 1.5,
    vm = list(value = c("decon", "bold"), key = c("vnum", "region"))
  )

  oracle <- get_medusa_interpolated_ts(
    f_obj, event = "onset", time_before = -3, time_after = 3,
    output_resolution = 0.75, group_by = c("region", "trial"),
    valid_col = "ok", time_audit = TRUE, use_legacy = TRUE
  )

  fast_fmri <- get_medusa_interpolated_ts(
    f_obj, event = "onset", time_before = -3, time_after = 3,
    output_resolution = 0.75, group_by = c("region", "trial"),
    valid_col = "ok", time_audit = TRUE, use_legacy = FALSE
  )

  mat_multi <- align_medusa_ts(
    list(decon = sig_a, bold = sig_b),
    events = ev_df,
    event = "onset",
    time = times,
    tr = 1.5,
    unit_data = data.frame(vnum = seq_len(n_v), region = parcels, stringsAsFactors = FALSE),
    group_by = c("region", "trial"),
    valid = valid_mat,
    time_before = -3, time_after = 3,
    output_resolution = 0.75,
    time_audit = TRUE,
    list_mode = "signals"
  )

  expect_equal(fast_fmri, oracle, tolerance = 1e-12)
  expect_equal(mat_multi, oracle, tolerance = 1e-12)
})

test_that("ragged per-unit clocks in fmri_ts fall back cleanly and preserve unit-specific gaps", {
  # Unit 1 has t = 0..6; Unit 2 is missing row t = 3 (deleted row only in Unit 2)
  d1 <- data.frame(time = 0:6, unit = "u1", value = as.numeric(0:6))
  d2 <- data.frame(time = c(0:2, 4:6), unit = "u2", value = as.numeric(c(0:2, 4:6)) * 2)
  d_ragged <- rbind(d1, d2)
  ev <- data.frame(id = "A", run = 1L, trial = 1L, onset = 3.0)
  f_ragged <- fmri_ts$new(d_ragged, ev, tr = 1.0, vm = list(key = "unit"))

  expect_null(medusa:::extract_fmri_ts_matrix_view(f_ragged))
  expect_error(f_ragged$to_matrix(), "common acquisition clock")

  # At max_gap = 1, unit 2 has a gap at evt_time = 0; at max_gap = 2, unit 2 bridges [2, 4] -> 6
  r_strict <- get_medusa_interpolated_ts(f_ragged, event = "onset", time_before = 0, time_after = 0, max_gap = 1)
  expect_equal(r_strict$value_n, 1L)
  expect_equal(r_strict$value_n_gap, 1L)
  expect_equal(r_strict$value_mean, 3.0)

  r_bridged <- get_medusa_interpolated_ts(f_ragged, event = "onset", time_before = 0, time_after = 0, max_gap = 2)
  expect_equal(r_bridged$value_n, 2L)
  expect_equal(r_bridged$value_mean, mean(c(3.0, 6.0)))
})

test_that("read_nifti_ts and deconvolve_crossrun_hrf integrate directly with align_medusa_ts", {
  set.seed(20261115)
  dim3 <- c(3L, 3L, 2L)
  n_t <- 40L
  tr <- 1.0

  atlas_arr <- array(0L, dim = dim3)
  atlas_arr[1:2, 1:2, 1] <- 1L
  atlas_arr[2:3, 2:3, 2] <- 2L
  f_atlas <- tempfile(fileext = ".nii.gz")
  nii_atlas <- RNifti::asNifti(atlas_arr)
  RNifti::pixdim(nii_atlas) <- c(2, 2, 2)
  RNifti::writeNifti(nii_atlas, f_atlas)

  f_r1 <- tempfile(fileext = ".nii.gz")
  f_r2 <- tempfile(fileext = ".nii.gz")
  arr1 <- array(rnorm(prod(dim3) * n_t), dim = c(dim3, n_t))
  arr2 <- array(rnorm(prod(dim3) * n_t), dim = c(dim3, n_t))
  nii1 <- RNifti::asNifti(arr1)
  nii2 <- RNifti::asNifti(arr2)
  RNifti::pixdim(nii1) <- c(2, 2, 2, tr)
  RNifti::pixdim(nii2) <- c(2, 2, 2, tr)
  RNifti::writeNifti(nii1, f_r1)
  RNifti::writeNifti(nii2, f_r2)
  on.exit(unlink(c(f_atlas, f_r1, f_r2)), add = TRUE)

  ts_bundle <- read_nifti_ts(c(run1 = f_r1, run2 = f_r2), atlas = f_atlas)
  ev_combined <- data.frame(
    id = "sub01",
    run = c("run1", "run1", "run2", "run2"),
    trial = c(1L, 2L, 1L, 2L),
    onset = c(8.5, 22.0, 10.25, 26.75),
    duration = 1.0
  )

  aligned_multi <- align_medusa_ts(
    ts_bundle,
    events = ev_combined,
    event = "onset",
    time_before = -2,
    time_after = 3,
    output_resolution = 1.0,
    ncpus = 2L
  )

  # Compare against per-run alignment
  a_r1 <- align_medusa_ts(
    ts_bundle$bold$run1,
    events = ev_combined[ev_combined$run == "run1", ],
    event = "onset",
    tr = tr,
    parcel_ids = ts_bundle$parcel_ids,
    time_before = -2,
    time_after = 3,
    output_resolution = 1.0
  )
  expect_equal(
    aligned_multi[run == "run1", !"run"],
    a_r1,
    ignore_attr = TRUE,
    tolerance = 1e-12
  )
  expect_equal(sort(unique(aligned_multi$parcel)), c(1, 2))
  expect_equal(sort(unique(attr(aligned_multi, "source_clock")$run)), c("run1", "run2"))

  # Cross-run HRF + deconvolution fit with observation_mask and time_offset
  ev_hrf <- ev_combined
  ev_hrf$condition <- "cue"
  ev_list <- list(
    run1 = ev_hrf[ev_hrf$run == "run1", ],
    run2 = ev_hrf[ev_hrf$run == "run2", ]
  )
  valid_rows <- list(
    run1 = seq_len(n_t) != 9L,
    run2 = rep(TRUE, n_t)
  )
  hrf_fit <- estimate_crossrun_hrf(
    ts_bundle$bold,
    events = ev_list,
    tr = tr,
    time_offset = c(run1 = 0.0, run2 = 2.0),
    valid_rows = valid_rows,
    spatial_groups = ts_bundle$parcel_ids,
    kernel_grid = hrf_library(tr, peak_delays = c(5, 6), dispersions = 1, undershoot_ratios = 6)
  )
  dec <- deconvolve_crossrun_hrf(ts_bundle$bold, hrf_fit, lambda = 0.5)

  aligned_dec <- align_medusa_ts(
    dec,
    event = "onset",
    parcel_ids = ts_bundle$parcel_ids,
    valid = "observation_mask",
    time_before = -1,
    time_after = 2,
    output_resolution = 1.0,
    time_audit = TRUE
  )
  # Frame 9 in run1 corresponds to t = 8.0s, bracketing run1 trial 1 onset = 8.5s at evt_time = 0
  r1_t1_zero <- aligned_dec[run == "run1" & trial == 1L & evt_time == 0]
  expect_true(all(is.na(r1_t1_zero$value_mean)))
  expect_true(all(r1_t1_zero$value_n_censored > 0L))
  # Run 2 source clock starts at time_offset = 2.0s
  clk_dec <- attr(aligned_dec, "source_clock")
  expect_equal(min(clk_dec$time[clk_dec$run == "run2"]), 2.0)
})

test_that("matrix-native PCA compression scores and RDS/CSV file alignment match legacy outputs", {
  set.seed(20261116)
  n_t <- 50L
  n_v <- 10L
  parcels <- rep(1:2, each = 5L)
  t_seq <- (seq_len(n_t) - 1) * 1.0
  y_mat <- outer(sin(t_seq / 3), seq(0.5, 1.5, length.out = n_v)) +
    matrix(rnorm(n_t * n_v, sd = 0.2), nrow = n_t, ncol = n_v)

  long_df <- expand.grid(idx = seq_len(n_t), vnum = seq_len(n_v))
  long_df$time <- t_seq[long_df$idx]
  long_df$atlas_value <- parcels[long_df$vnum]
  long_df$decon <- y_mat[cbind(long_df$idx, long_df$vnum)]
  long_df$idx <- NULL

  ev_df <- data.frame(
    id = "A",
    run = 1L,
    trial = 1:4,
    onset = c(6, 17, 29, 41),
    bound = c(3, 14, 25, 48)
  )

  f_obj <- fmri_ts$new(
    long_df, ev_df, tr = 1.0,
    vm = list(value = "decon", key = c("vnum", "atlas_value"))
  )

  comp_legacy <- get_medusa_compression_score(
    f_obj, event = "onset", time_before = -2, time_after = 4,
    collide_before = "bound", collide_after = "bound",
    group_by = c("trial", "atlas_value"), use_legacy = TRUE
  )
  comp_fast <- get_medusa_compression_score(
    f_obj, event = "onset", time_before = -2, time_after = 4,
    collide_before = "bound", collide_after = "bound",
    group_by = c("trial", "atlas_value"), use_legacy = FALSE
  )
  comp_mat <- get_medusa_compression_score(
    y_mat, events = ev_df, event = "onset", tr = 1.0,
    unit_data = data.frame(vnum = seq_len(n_v), atlas_value = parcels),
    value_name = "decon",
    time_before = -2, time_after = 4,
    collide_before = "bound", collide_after = "bound",
    group_by = c("trial", "atlas_value")
  )
  expect_equal(comp_fast, comp_legacy, tolerance = 1e-12)
  expect_equal(comp_mat, comp_legacy, tolerance = 1e-12)

  # Test evt_align_decon_files with both .csv.gz (mat2df layout) and .rds matrix bundle
  tmp_dir <- tempfile("decon_align_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  csv_path <- file.path(tmp_dir, "subA_run1_deconvolved.csv.gz")
  rds_path <- file.path(tmp_dir, "subA_run1_deconvolved.rds")
  coords <- data.frame(vnum = seq_len(n_v), atlas_value = parcels)
  m_df <- medusa:::mat2df(unname(t(y_mat)), value_name = "decon", varnames = c("vnum", "volume"))
  m_df$time <- (m_df$volume - 1) * 1.0
  m_df <- merge(m_df, coords, by = "vnum", sort = FALSE)
  data.table::fwrite(m_df, csv_path)
  saveRDS(list(decon = y_mat, coords = coords, time = t_seq), rds_path)

  al_spec <- list(
    evt_col = "onset", time_before = -2, time_after = 3,
    pad_before = -1, pad_after = 1, output_resolution = 0.5,
    collide_before = "bound", collide_after = "bound",
    readout_kernel = "triangle"
  )
  res_csv <- evt_align_decon_files(csv_path, ev_df, al_spec, tr = 1.0, ncpus = 1L, time_audit = TRUE)
  res_rds <- evt_align_decon_files(rds_path, ev_df, al_spec, tr = 1.0, ncpus = 1L, time_audit = TRUE)

  oracle_file <- get_medusa_interpolated_ts(
    f_obj, event = "onset", time_before = -2, time_after = 3,
    pad_before = -1, pad_after = 1, output_resolution = 0.5,
    collide_before = "bound", collide_after = "bound",
    group_by = c("atlas_value", "trial"), readout_kernel = "triangle",
    time_audit = TRUE, use_legacy = TRUE
  )

  shared_cols <- names(oracle_file)
  expect_equal(res_csv[, ..shared_cols], oracle_file, ignore_attr = TRUE, tolerance = 1e-12)
  expect_equal(res_rds[, ..shared_cols], oracle_file, ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("matrix-native alignment substantially outperforms legacy fmri_ts workflow", {
  set.seed(20261117)
  n_t <- 200L
  n_v <- 800L
  n_tr <- 30L
  parcels <- rep(1:10, length.out = n_v)
  y_mat <- matrix(rnorm(n_t * n_v), nrow = n_t, ncol = n_v)
  times <- (seq_len(n_t) - 1) * 1.0
  ev_df <- data.frame(
    id = "sub01", run = 1L, trial = seq_len(n_tr),
    onset = seq(6, by = 6, length.out = n_tr) + runif(n_tr, -0.4, 0.4)
  )

  long_df <- expand.grid(idx = seq_len(n_t), vnum = seq_len(n_v))
  long_df$time <- times[long_df$idx]
  long_df$atlas_value <- parcels[long_df$vnum]
  long_df$decon <- y_mat[cbind(long_df$idx, long_df$vnum)]
  long_df$idx <- NULL

  t_legacy <- system.time({
    f_obj <- fmri_ts$new(long_df, ev_df, tr = 1.0, vm = list(value = "decon", key = c("vnum", "atlas_value")))
    out_legacy <- get_medusa_interpolated_ts(
      f_obj, event = "onset", time_before = -2, time_after = 4,
      group_by = c("atlas_value", "trial"), use_legacy = TRUE
    )
  })[["elapsed"]]

  t_matrix <- system.time({
    out_matrix <- align_medusa_ts(
      y_mat, events = ev_df, event = "onset", tr = 1.0,
      unit_data = data.frame(vnum = seq_len(n_v), atlas_value = parcels),
      group_by = c("atlas_value", "trial"), value_name = "decon",
      time_before = -2, time_after = 4
    )
  })[["elapsed"]]

  expect_equal(out_matrix, out_legacy, tolerance = 1e-12)
  expect_lt(t_matrix, t_legacy / 5)
})
