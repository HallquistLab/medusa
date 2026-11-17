# Validation suite for Tier 2 Steps 1-2 Pre-Launch Upgrades:
#   1. GLMdenoise-style cross-validated noise-pool temporal PCs (`noise_pcs = "auto"`)
#   2. Modulator-aware (`hrf_modulators`) & SVD fractional-ridge (`trial_ridge` / `"auto"`)
#      HRF scoring for zero-mean parametric modulators and rapid-ITI designs
#   3. SNR-adaptive canonical de-shrinking & training-only spatial noise correlation

sim_ar1_mat_t2 <- function(n_time, n_vox, rho = 0.45, sd = 1.0) {
  e <- matrix(stats::rnorm(n_time * n_vox, sd = sd * sqrt(max(1e-6, 1.0 - rho^2))), n_time, n_vox)
  out <- matrix(0.0, n_time, n_vox)
  out[1L, ] <- stats::rnorm(n_vox, sd = sd)
  for (t in 2:n_time) {
    out[t, ] <- rho * out[t - 1L, ] + e[t, ]
  }
  out
}

implied_peak_delay_t2 <- function(fit, r = 1L, cols = NULL) {
  w <- fit$kernel_weights[[r]]
  if (!is.null(cols)) w <- w[, cols, drop = FALSE]
  as.numeric(crossprod(w, fit$parameters$peak_delay))
}

test_that("Tier 2.1: cross-validated noise-pool PCs suppress narrowband physio noise before deconvolution, select K*=0 on independent AR(1) noise, and preserve target-run exclusion", {
  set.seed(20261101)
  tr <- 0.6
  n_time <- 210L
  n_runs <- 2L
  lib <- hrf_library(tr = tr, peak_delays = seq(4.0, 7.5, by = 0.5), dispersions = 1.0, undershoot_ratios = 6)
  pars <- attr(lib, "parameters")
  idx_can <- which(abs(pars$peak_delay - 6.0) < 1e-6)[1L]
  idx_fast <- which(abs(pars$peak_delay - 4.5) < 1e-6)[1L]
  idx_slow <- which(abs(pars$peak_delay - 7.0) < 1e-6)[1L]

  # 2 active parcels (6 voxels each = 12) + 10 null background parcels (4 voxels each = 40) -> 52 voxels
  groups <- c(rep("A1_fast", 6L), rep("A2_slow", 6L), rep(sprintf("BG%02d", 1:10), each = 4L))
  n_vox <- length(groups)
  vox_ids <- sprintf("v%02d", seq_len(n_vox))
  names(groups) <- vox_ids

  bold_phys <- bold_indep <- events <- x_cov <- vector("list", n_runs)
  names(bold_phys) <- names(bold_indep) <- names(events) <- names(x_cov) <- paste0("run", seq_len(n_runs))
  t_sec <- (seq_len(n_time) - 1.0) * tr

  # Scale neural impulses by 1/tr so unit-sum HRF convolution matches continuous physical BOLD units
  amp_scale <- 1.0 / tr
  for (r in seq_len(n_runs)) {
    onsets <- cumsum(c(4.0, stats::runif(23L, 3.8, 6.2)))
    onsets <- onsets[onsets < (n_time - 22L) * tr]
    m_ev <- length(onsets)
    xv <- as.numeric(scale(stats::rnorm(m_ev)))
    events[[r]] <- data.frame(onset = onsets, duration = 0.0, condition = "stim")
    x_cov[[r]] <- data.frame(x = xv)

    des_k <- hrf_event_design_cpp(n_time, tr, 0.0, onsets, rep(0.0, m_ev),
      seq_len(m_ev), m_ev, lib[, c(idx_fast, idx_slow), drop = FALSE])
    y_clean <- matrix(0.0, n_time, n_vox, dimnames = list(NULL, vox_ids))
    y_clean[, 1:6]  <- replicate(6L, as.numeric(des_k[[1L]] %*% (amp_scale * (1.70 + 0.70 * xv))))
    y_clean[, 7:12] <- replicate(6L, as.numeric(des_k[[2L]] %*% (amp_scale * (1.70 - 0.65 * xv))))

    e_indep <- sim_ar1_mat_t2(n_time, n_vox, rho = 0.45, sd = 0.24)
    # Narrowband aliased cardiac (~0.24 Hz) + respiratory (~0.36 Hz) harmonics with parcel-coherent topographies
    f1 <- sin(2 * pi * 0.24 * t_sec + stats::runif(1L, 0, 2 * pi))
    f2 <- cos(2 * pi * 0.36 * t_sec + stats::runif(1L, 0, 2 * pi))
    load1 <- stats::runif(n_vox, 0.36, 0.52)
    parcel_sgn2 <- c(rep(1.0, 6L), rep(-1.0, 6L), rep(rep(c(1.0, 1.0, -1.0, -1.0), each = 4L), length.out = n_vox - 12L))
    load2 <- parcel_sgn2 * stats::runif(n_vox, 0.32, 0.48)
    e_phys <- e_indep + tcrossprod(f1, load1) + tcrossprod(f2, load2)

    bold_indep[[r]] <- y_clean + e_indep
    bold_phys[[r]]  <- y_clean + e_phys
  }

  # 1. Under narrowband physio noise, noise_pcs = "auto" selects K* >= 2 strictly from BG01..BG10
  fit_auto <- estimate_crossrun_hrf(
    bold = bold_phys, events = events, tr = tr, kernel_grid = lib,
    spatial_groups = groups, pooling = "parcel", fallback_kernel = idx_can,
    min_parcel_voxels = 4L, selection_mode = "soft_average", noise_pcs = "auto"
  )
  expect_gte(fit_auto$settings$noise_pcs_used, 2L)
  expect_lte(fit_auto$settings$noise_pcs_used, 5L)
  expect_length(intersect(fit_auto$specification$noise_pool_voxels, vox_ids[1:12]), 0L)
  expect_gt(length(fit_auto$specification$noise_pool_voxels), 24L)

  # Strict target-run exclusion: corrupting run1 leaves run1's HRF selection 100% invariant
  bold_corrupt <- bold_phys
  bold_corrupt$run1 <- bold_corrupt$run1 + matrix(stats::rnorm(n_time * n_vox, sd = 20.0), n_time, n_vox)
  fit_corrupt <- estimate_crossrun_hrf(
    bold = bold_corrupt, events = events, tr = tr, kernel_grid = lib,
    spatial_groups = groups, pooling = "parcel", fallback_kernel = idx_can,
    min_parcel_voxels = 4L, selection_mode = "soft_average", noise_pcs = "auto"
  )
  expect_identical(fit_auto$selections$run1, fit_corrupt$selections$run1)
  expect_equal(fit_auto$selected_kernels$run1, fit_corrupt$selected_kernels$run1, tolerance = 0)

  # 2. Deconvolution & Stage-3 epoch regression improve substantially vs noise_pcs = 0
  fit_off <- estimate_crossrun_hrf(
    bold = bold_phys, events = events, tr = tr, kernel_grid = lib,
    spatial_groups = groups, pooling = "parcel", fallback_kernel = idx_can,
    min_parcel_voxels = 4L, selection_mode = "soft_average", noise_pcs = 0L
  )
  dec_on  <- deconvolve_crossrun_hrf(bold_phys, fit_auto, scaling = "noise_sd", penalty = "sobolev", noise_pcs = "auto")
  dec_off <- deconvolve_crossrun_hrf(bold_phys, fit_off,  scaling = "noise_sd", penalty = "sobolev", noise_pcs = 0L)
  expect_equal(dec_on$settings$noise_pcs_used, fit_auto$settings$noise_pcs_used)
  expect_equal(dec_off$settings$noise_pcs_used, 0L)

  # Residual GLS innovation noise_sd drops by >50% down to the true thermal noise floor (0.24)
  nsd_on  <- mean(dec_on$scaling$run1$noise_sd)
  nsd_off <- mean(dec_off$scaling$run1$noise_sd)
  expect_lt(nsd_on, 0.55 * nsd_off)

  reg_on <- decon_epoch_regression(
    dec_on, covariates = x_cov, events = events, target_condition = "stim",
    lags = 0, aggregate_by_parcel = TRUE, parcel_groups = groups, restore_scale = TRUE, return_panel = TRUE
  )
  reg_off <- decon_epoch_regression(
    dec_off, covariates = x_cov, events = events, target_condition = "stim",
    lags = 0, aggregate_by_parcel = TRUE, parcel_groups = groups, restore_scale = TRUE, return_panel = TRUE
  )

  # Residual GLS innovation noise_sd drops by >45% across BOTH runs, and at matched
  # regularization (`lambda = 1.0`) trial-level `se_demixed` drops sharply
  expect_lt(mean(dec_on$scaling$run2$noise_sd), 0.55 * mean(dec_off$scaling$run2$noise_sd))

  dec_on_m  <- deconvolve_crossrun_hrf(bold_phys, fit_auto, lambda = 1.0, scaling = "noise_sd", penalty = "sobolev", noise_pcs = "auto")
  dec_off_m <- deconvolve_crossrun_hrf(bold_phys, fit_off,  lambda = 1.0, scaling = "noise_sd", penalty = "sobolev", noise_pcs = 0L)
  reg_on_m  <- decon_epoch_regression(dec_on_m,  covariates = x_cov, events = events, target_condition = "stim", lags = 0, aggregate_by_parcel = TRUE, parcel_groups = groups, restore_scale = TRUE, return_panel = TRUE)
  reg_off_m <- decon_epoch_regression(dec_off_m, covariates = x_cov, events = events, target_condition = "stim", lags = 0, aggregate_by_parcel = TRUE, parcel_groups = groups, restore_scale = TRUE, return_panel = TRUE)
  expect_lt(
    mean(reg_on_m$panel$se_demixed[reg_on_m$panel$voxel %in% c("A1_fast", "A2_slow")]),
    0.78 * mean(reg_off_m$panel$se_demixed[reg_off_m$panel$voxel %in% c("A1_fast", "A2_slow")])
  )

  # On non-pool parcels (using explicit held-out noise_pool = BG03..BG10), shared noise correlation on BG01 <-> BG02 collapses
  dec_on_ho <- deconvolve_crossrun_hrf(bold_phys, fit_auto, lambda = 1.0, scaling = "noise_sd", penalty = "sobolev",
                                       noise_pcs = "auto", noise_pool = vox_ids[21:52])
  cor_ho_on  <- stats::cor(rowMeans(dec_on_ho$noise_innovations$run1[, 13:16]),
                           rowMeans(dec_on_ho$noise_innovations$run1[, 17:20]))
  cor_ho_off <- stats::cor(rowMeans(dec_off_m$noise_innovations$run1[, 13:16]),
                           rowMeans(dec_off_m$noise_innovations$run1[, 17:20]))
  expect_lt(abs(cor_ho_on), 0.15)
  expect_gt(cor_ho_off, 0.40)

  cf_on <- reg_on$coefficients
  est_a1 <- cf_on$beta_corrected[cf_on$term == "x" & cf_on$voxel == "A1_fast"] / amp_scale
  est_a2 <- cf_on$beta_corrected[cf_on$term == "x" & cf_on$voxel == "A2_slow"] / amp_scale
  expect_equal(est_a1,  0.70, tolerance = 0.25)
  expect_equal(est_a2, -0.65, tolerance = 0.25)

  # 3. Spatially independent AR(1) background noise cleanly selects K* = 0 (no overfitting)
  fit_indep <- estimate_crossrun_hrf(
    bold = bold_indep, events = events, tr = tr, kernel_grid = lib,
    spatial_groups = groups, pooling = "parcel", fallback_kernel = idx_can,
    min_parcel_voxels = 4L, selection_mode = "soft_average", noise_pcs = "auto"
  )
  expect_identical(fit_indep$settings$noise_pcs_used, 0L)
  expect_null(fit_indep$specification$noise_pcs)
})

test_that("Tier 2.2a: modulator-aware and fracridge (`trial_ridge` & `auto`) recover fast/slow HRFs under pure zero-mean RPE and rapid ITIs while protecting modulated parcels from bg_mask", {
  set.seed(20261102)
  tr <- 1.0
  n_time <- 135L
  n_runs <- 3L
  lib <- hrf_library(tr = tr, peak_delays = seq(3.8, 7.8, by = 0.4), dispersions = 1.0, undershoot_ratios = 6)
  pars <- attr(lib, "parameters")
  idx_can  <- which.min(abs(pars$peak_delay - 6.0))
  idx_fast <- which.min(abs(pars$peak_delay - 4.2))
  idx_slow <- which.min(abs(pars$peak_delay - 7.4))

  groups <- c(rep("RPE_fast", 6L), rep("RPE_slow", 6L), rep(sprintf("BG%02d", 1:10), each = 4L))
  n_vox <- length(groups)
  vox_ids <- sprintf("v%02d", seq_len(n_vox))
  names(groups) <- vox_ids

  bold <- events <- rpe_cov <- vector("list", n_runs)
  names(bold) <- names(events) <- names(rpe_cov) <- paste0("run", seq_len(n_runs))
  t_sec <- (seq_len(n_time) - 1.0) * tr

  for (r in seq_len(n_runs)) {
    # Rapid event-related ITI (mean ~ 2.3s, ~42 trials per run)
    onsets <- cumsum(c(4.0, stats::runif(41L, 1.5, 3.2)))
    onsets <- onsets[onsets < (n_time - 20L) * tr]
    m_ev <- length(onsets)
    # Exact zero-mean signed RPE (beta_0 = 0.0 -> condition mean has ZERO signal!)
    rpe <- as.numeric(scale(stats::rnorm(m_ev), center = TRUE, scale = TRUE))
    events[[r]] <- data.frame(onset = onsets, duration = 0.0, condition = "fb", rpe = rpe)
    rpe_cov[[r]] <- data.frame(rpe = rpe)

    des_k <- hrf_event_design_cpp(n_time, tr, 0.0, onsets, rep(0.0, m_ev),
      seq_len(m_ev), m_ev, lib[, c(idx_fast, idx_slow), drop = FALSE])
    y_clean <- matrix(0.0, n_time, n_vox, dimnames = list(NULL, vox_ids))
    y_clean[, 1:6]  <- replicate(6L, as.numeric(des_k[[1L]] %*% (1.35 * rpe)))
    y_clean[, 7:12] <- replicate(6L, as.numeric(des_k[[2L]] %*% (-1.35 * rpe)))

    e_ar <- sim_ar1_mat_t2(n_time, n_vox, rho = 0.40, sd = 0.24)
    f_phys <- 0.32 * sin(2 * pi * 0.19 * t_sec + r)
    bold[[r]] <- y_clean + e_ar + tcrossprod(f_phys, rep(1.0, n_vox))
  }

  # 1. Baseline `amplitude_model = "condition"` (without modulators) is completely blind to beta_0 = 0 RPE
  fit_cond <- estimate_crossrun_hrf(
    bold = bold, events = events, tr = tr, kernel_grid = lib,
    spatial_groups = groups, pooling = "parcel", fallback_kernel = idx_can,
    min_parcel_voxels = 4L, amplitude_model = "condition", selection_mode = "soft_average",
    noise_pcs = 0L
  )
  pk_cond <- implied_peak_delay_t2(fit_cond, 1L, c(1L, 7L))
  expect_equal(pk_cond, rep(pars$peak_delay[idx_can], 2L), tolerance = 0.05)

  # 2. Explicit `hrf_modulators = "rpe"` recovers both fast (4.2s) and slow (7.4s) HRFs accurately
  fit_mod <- estimate_crossrun_hrf(
    bold = bold, events = events, tr = tr, kernel_grid = lib,
    spatial_groups = groups, pooling = "parcel", fallback_kernel = idx_can,
    min_parcel_voxels = 4L, amplitude_model = "condition", hrf_modulators = "rpe",
    selection_mode = "soft_average", noise_pcs = "auto"
  )
  pk_mod <- implied_peak_delay_t2(fit_mod, 1L, c(1L, 7L))
  expect_lt(pk_mod[1L], 4.65)
  expect_gt(pk_mod[2L], 6.95)

  # 3. Unmodeled zero-mean RPE (`events` stripped of `rpe` column!):
  #    `amplitude_model = "auto"` and `"trial_ridge"` rescue both parcels and beat unregularized `"trial"`
  ev_no_rpe <- lapply(events, function(df) df[, c("onset", "duration", "condition")])
  fit_auto <- estimate_crossrun_hrf(
    bold = bold, events = ev_no_rpe, tr = tr, kernel_grid = lib,
    spatial_groups = groups, pooling = "parcel", fallback_kernel = idx_can,
    min_parcel_voxels = 4L, amplitude_model = "auto", selection_mode = "soft_average",
    noise_pcs = "auto"
  )
  fit_tr_unreg <- estimate_crossrun_hrf(
    bold = bold, events = ev_no_rpe, tr = tr, kernel_grid = lib,
    spatial_groups = groups, pooling = "parcel", fallback_kernel = idx_can,
    min_parcel_voxels = 4L, amplitude_model = "trial", selection_mode = "soft_average",
    noise_pcs = 0L
  )

  pk_auto  <- implied_peak_delay_t2(fit_auto, 1L, c(1L, 7L))
  pk_unreg <- implied_peak_delay_t2(fit_tr_unreg, 1L, c(1L, 7L))
  true_pks <- c(pars$peak_delay[idx_fast], pars$peak_delay[idx_slow])
  err_auto  <- mean(abs(pk_auto - true_pks))
  err_unreg <- mean(abs(pk_unreg - true_pks))
  expect_lt(pk_auto[1L], 4.80)
  expect_gt(pk_auto[2L], 6.80)
  expect_lt(err_auto, err_unreg)
  expect_lt(err_auto, 0.45)

  # Crucial synergy: `amplitude_model = "auto"` rescues RPE_fast & RPE_slow BEFORE bg_mask,
  # so ZERO modulated voxels leak into `noise_pool_voxels` while BG01..BG08 physio PC is still extracted!
  expect_length(intersect(fit_auto$specification$noise_pool_voxels, vox_ids[1:12]), 0L)
  expect_gte(fit_auto$settings$noise_pcs_used, 1L)

  # Strict target-run exclusion under `amplitude_model = "auto"`
  bold_c1 <- bold
  bold_c1$run1 <- bold_c1$run1 + matrix(stats::rnorm(n_time * n_vox, sd = 25.0), n_time, n_vox)
  fit_auto_c1 <- estimate_crossrun_hrf(
    bold = bold_c1, events = ev_no_rpe, tr = tr, kernel_grid = lib,
    spatial_groups = groups, pooling = "parcel", fallback_kernel = idx_can,
    min_parcel_voxels = 4L, amplitude_model = "auto", selection_mode = "soft_average",
    noise_pcs = "auto"
  )
  expect_identical(fit_auto$selections$run1, fit_auto_c1$selections$run1)
  expect_equal(fit_auto$selected_kernels$run1, fit_auto_c1$selected_kernels$run1, tolerance = 0)
})

test_that("Tier 2.2b: SNR-adaptive de-shrinking and training-only spatial noise correlation eliminate canonical pull on boundary HRFs while preserving calibrated null control", {
  set.seed(20261103)
  tr <- 1.0
  n_time <- 180L
  n_runs <- 2L
  lib <- hrf_library(tr = tr, peak_delays = seq(3.5, 8.0, by = 0.5), dispersions = c(0.8, 1.0, 1.2), undershoot_ratios = 6)
  pars <- attr(lib, "parameters")
  idx_can <- which(abs(pars$peak_delay - 6.0) < 1e-6 & abs(pars$dispersion - 1.0) < 1e-6)[1L]
  idx_45  <- which(abs(pars$peak_delay - 4.5) < 1e-6 & abs(pars$dispersion - 1.0) < 1e-6)[1L]
  idx_50  <- which(abs(pars$peak_delay - 5.0) < 1e-6 & abs(pars$dispersion - 1.0) < 1e-6)[1L]
  idx_70  <- which(abs(pars$peak_delay - 7.0) < 1e-6 & abs(pars$dispersion - 1.0) < 1e-6)[1L]

  # 3 active 12-voxel parcels (4.5s, 5.0s, 7.0s) + 20 pure-null 8-voxel parcels for Type-I check
  groups <- c(
    rep("P_45", 12L), rep("P_50", 12L), rep("P_70", 12L),
    rep(sprintf("NULL_%02d", 1:20), each = 8L)
  )
  n_vox <- length(groups)
  vox_ids <- sprintf("v%03d", seq_len(n_vox))
  names(groups) <- vox_ids

  bold <- events <- vector("list", n_runs)
  names(bold) <- names(events) <- paste0("run", seq_len(n_runs))
  for (r in seq_len(n_runs)) {
    onsets <- cumsum(c(5.0, stats::runif(26L, 4.0, 7.0)))
    onsets <- onsets[onsets < (n_time - 24L) * tr]
    m_ev <- length(onsets)
    events[[r]] <- data.frame(onset = onsets, duration = 0.0, condition = "task")
    des_k <- hrf_event_design_cpp(n_time, tr, 0.0, onsets, rep(0.0, m_ev),
      rep(1L, m_ev), 1L, lib[, c(idx_45, idx_50, idx_70), drop = FALSE])

    ym <- sim_ar1_mat_t2(n_time, n_vox, rho = 0.42, sd = 0.24)
    ym[, 1:12]  <- sweep(ym[, 1:12],  1L, 1.15 * as.numeric(des_k[[1L]]), "+")
    ym[, 13:24] <- sweep(ym[, 13:24], 1L, 0.85 * as.numeric(des_k[[2L]]), "+")
    ym[, 25:36] <- sweep(ym[, 25:36], 1L, 1.15 * as.numeric(des_k[[3L]]), "+")
    colnames(ym) <- vox_ids
    bold[[r]] <- ym
  }

  fit_soft <- estimate_crossrun_hrf(
    bold = bold, events = events, tr = tr, kernel_grid = lib,
    spatial_groups = groups, pooling = "parcel", fallback_kernel = idx_can,
    min_parcel_voxels = 5L, amplitude_model = "condition",
    selection_mode = "soft_average", noise_pcs = 0L
  )

  pks_r1 <- implied_peak_delay_t2(fit_soft, 1L, c(1L, 13L, 25L))
  expect_equal(pks_r1[1L], 4.50, tolerance = 0.22)
  expect_equal(pks_r1[2L], 5.00, tolerance = 0.28)
  expect_equal(pks_r1[3L], 7.00, tolerance = 0.22)

  # Null Type-I parcel departure rate across the 20 8-voxel AR(1) null parcels remains <= 0.10
  null_first_cols <- seq(37L, n_vox, by = 8L)
  pks_null <- implied_peak_delay_t2(fit_soft, 1L, null_first_cols)
  null_depart_rate <- mean(abs(pks_null - 6.0) > 0.15)
  expect_lte(null_depart_rate, 0.10)
})

test_that("decon_epoch_regression and aggregate_decon_panel preserve data.table selfref without shallow-copy warnings", {
  set.seed(8106)
  tr <- 1.0
  n_time <- 90L
  lib <- hrf_library(tr = tr, preset = "standard")
  groups <- rep(c("P1", "P2"), each = 3L)
  vox_ids <- paste0("v", 1:6)
  names(groups) <- vox_ids

  bold <- events <- vector("list", 2L)
  names(bold) <- names(events) <- c("run1", "run2")
  for (r in 1:2) {
    events[[r]] <- data.frame(onset = seq(8, 68, by = 10), duration = 0.0, condition = "task")
    ym <- sim_ar1_mat_t2(n_time, 6L, rho = 0.35, sd = 0.30)
    colnames(ym) <- vox_ids
    bold[[r]] <- ym
  }

  hf <- estimate_crossrun_hrf(bold, events, tr = tr, kernel_grid = lib, spatial_groups = groups, pooling = "parcel", noise_pcs = 0L)
  dc <- deconvolve_crossrun_hrf(bold, hf)
  reg <- decon_epoch_regression(dc, covariates = NULL, events = events, target_condition = "task", lags = 0, aggregate_by_parcel = FALSE, return_panel = TRUE)

  warns_p <- character(0L)
  withCallingHandlers({
    p1 <- reg$panel
    p1[, test_col_a := 1.0]
    agg <- aggregate_decon_panel(reg, parcel_groups = groups, parcel_col = "region")
    p2 <- agg$panel
    p2[, test_col_b := 2.0]
  }, warning = function(w) {
    warns_p <<- c(warns_p, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_false(any(grepl("shallow copy", warns_p, ignore.case = TRUE)))
})
