test_that("Track 2a: SNR-stratified lambda_mode = 'auto'/'stratified' enforces monotone tiers, quiets background, and collapses on null", {
  set.seed(20261114L)
  tr <- 1.0
  n_time <- 160L
  n_trials <- 28L
  lib <- hrf_library(tr, peak_delays = c(5.0, 6.0), dispersions = 1.0, undershoot_ratios = 6)
  h_true <- lib[, 2L, drop = FALSE]

  # 6 x 6 x 2 lattice (72 voxels): 12 high-SNR focal, 12 moderate periphery, 48 null background
  grid <- expand.grid(i = 1:6, j = 1:6, k = 1:2)
  V <- nrow(grid)
  vox_ids <- paste0("v", seq_len(V))
  vm <- data.frame(
    voxel = seq_len(V),
    i = grid$i, j = grid$j, k = grid$k,
    x = 3.0 * grid$i, y = 3.0 * grid$j, z = 3.0 * grid$k
  )
  geom <- list(dim = c(6L, 6L, 2L), spacing = c(3.0, 3.0, 3.0))

  focal_idx <- which(grid$i <= 2L & grid$j <= 3L)       # 12 voxels
  periph_idx <- which(grid$i == 3L & grid$j <= 6L)      # 12 voxels
  bg_idx <- setdiff(seq_len(V), c(focal_idx, periph_idx)) # 48 voxels

  sim_ar1 <- function(n, rho, sd_inn = 1.0) {
    as.numeric(stats::arima.sim(model = list(ar = rho), n = n, sd = sd_inn))
  }

  bold <- bold_null <- events <- vector("list", 2L)
  names(bold) <- names(bold_null) <- names(events) <- c("run1", "run2")
  neural_true_r1 <- matrix(0.0, n_time, length(focal_idx))

  for (r in 1:2) {
    onsets <- cumsum(c(6.0, stats::runif(n_trials - 1L, 3.8, 6.2)))
    amp <- 1.0 + 0.65 * as.numeric(scale(stats::rnorm(n_trials)))
    events[[r]] <- data.frame(onset = onsets, duration = 0.0, condition = "task")

    # True sparse neural impulse train + biphasic post-peak dip on focal voxels
    s_base <- numeric(n_time)
    t_grid <- (seq_len(n_time) - 1.0) * tr
    for (ti in seq_len(n_trials)) {
      k0 <- round(onsets[ti] / tr) + 1L
      if (k0 >= 1L && k0 <= n_time) s_base[k0] <- s_base[k0] + amp[ti]
      if (k0 + 2L <= n_time) s_base[k0 + 2L] <- s_base[k0 + 2L] - 0.45 * amp[ti]
    }
    y_sig <- stats::convolve(s_base, rev(as.numeric(h_true)), type = "open")[seq_len(n_time)]

    Y_r <- matrix(0.0, n_time, V, dimnames = list(NULL, vox_ids))
    Y_null_r <- matrix(0.0, n_time, V, dimnames = list(NULL, vox_ids))
    for (v in seq_len(V)) {
      eps <- sim_ar1(n_time, rho = 0.35, sd_inn = 1.0)
      Y_null_r[, v] <- eps
      if (v %in% focal_idx) {
        Y_r[, v] <- 7.5 * y_sig + eps
      } else if (v %in% periph_idx) {
        Y_r[, v] <- 3.2 * y_sig + eps
      } else {
        Y_r[, v] <- eps
      }
    }
    if (r == 1L) {
      for (fi in seq_along(focal_idx)) neural_true_r1[, fi] <- 7.5 * s_base
    }
    attr(Y_r, "voxel_map") <- vm
    attr(Y_r, "geometry") <- geom
    attr(Y_null_r, "voxel_map") <- vm
    attr(Y_null_r, "geometry") <- geom
    bold[[r]] <- Y_r
    bold_null[[r]] <- Y_null_r
  }

  fit_act <- estimate_crossrun_hrf(
    bold, events, lib, tr = tr, pooling = "voxel", noise_model = "ar1"
  )
  dec_auto <- deconvolve_crossrun_hrf(
    bold, fit_act, lambda = NULL, lambda_mode = "auto",
    scaling = "noise_sd", penalty = "sobolev"
  )
  dec_glob <- deconvolve_crossrun_hrf(
    bold, fit_act, lambda = NULL, lambda_mode = "global",
    scaling = "noise_sd", penalty = "sobolev"
  )

  expect_equal(dec_auto$settings$lambda_mode_used, "stratified")
  expect_true(is.data.frame(dec_auto$settings$lambda_strata))
  expect_gte(nrow(dec_auto$settings$lambda_strata), 2L)

  # Monotone non-increasing lambda across ascending task-informativeness strata
  lam_tiers <- dec_auto$settings$lambda_strata$lambda
  expect_true(all(diff(lam_tiers) <= 1e-12))
  expect_gt(lam_tiers[1L], lam_tiers[length(lam_tiers)])

  # Focal voxels receive strictly less damping than null background voxels
  lam_v <- dec_auto$scaling$run1$lambda
  expect_equal(length(lam_v), V)
  expect_lt(stats::median(lam_v[focal_idx]), stats::median(lam_v[bg_idx]))

  # Background variance is reduced relative to global lambda while focal neural correlation is preserved/sharpened
  bg_var_auto <- mean(apply(dec_auto$activity$run1[, bg_idx, drop = FALSE], 2L, stats::var))
  bg_var_glob <- mean(apply(dec_glob$activity$run1[, bg_idx, drop = FALSE], 2L, stats::var))
  foc_cor_auto <- mean(vapply(seq_along(focal_idx), function(i) {
    stats::cor(dec_auto$activity$run1[, focal_idx[i]], neural_true_r1[, i])
  }, numeric(1L)))
  foc_cor_glob <- mean(vapply(seq_along(focal_idx), function(i) {
    stats::cor(dec_glob$activity$run1[, focal_idx[i]], neural_true_r1[, i])
  }, numeric(1L)))

  expect_lt(bg_var_auto, bg_var_glob)
  expect_gte(foc_cor_auto, foc_cor_glob - 0.015)

  # Pure-null map automatically collapses to 1 global stratum under lambda_mode = 'auto'
  fit_null <- estimate_crossrun_hrf(
    bold_null, events, lib, tr = tr, pooling = "voxel", noise_model = "ar1"
  )
  dec_null <- deconvolve_crossrun_hrf(
    bold_null, fit_null, lambda = NULL, lambda_mode = "auto",
    scaling = "noise_sd", penalty = "sobolev"
  )
  expect_equal(dec_null$settings$lambda_mode_used, "global")
  expect_equal(length(unique(dec_null$scaling$run1$lambda)), 1L)
})

test_that("Track 2a: heterogeneous per-voxel lambda propagates cleanly through diagnose_crossrun_recovery and decon_epoch_regression", {
  set.seed(20261115L)
  tr <- 1.0
  n_time <- 180L
  n_trials <- 32L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  V <- 8L
  vox_ids <- paste0("v", seq_len(V))
  beta_true <- 1.25

  bold <- events <- covs <- vector("list", 2L)
  names(bold) <- names(events) <- names(covs) <- c("run1", "run2")
  for (r in 1:2) {
    onsets <- cumsum(c(6.0, stats::runif(n_trials - 1L, 4.2, 6.0)))
    x <- as.numeric(scale(stats::rnorm(n_trials)))
    events[[r]] <- data.frame(onset = onsets, duration = 0.0, condition = "cue")
    covs[[r]] <- data.frame(x = x)
    G <- hrf_event_design_cpp(
      n_time, tr, 0.0, onsets, rep(0.0, n_trials), seq_len(n_trials), n_trials, lib
    )[[1L]]
    y_clean <- as.numeric(G %*% (1.5 + beta_true * x))
    Y_r <- matrix(0.0, n_time, V, dimnames = list(NULL, vox_ids))
    for (v in seq_len(V)) {
      Y_r[, v] <- y_clean + as.numeric(stats::arima.sim(list(ar = 0.25), n = n_time, sd = 0.45))
    }
    bold[[r]] <- Y_r
  }

  fit_h <- estimate_crossrun_hrf(bold, events, lib, tr = tr, pooling = "voxel", noise_model = "ar1")
  # 15x spread in lambda across voxels within the same run and parcel
  lam_vec <- c(rep(0.08, 4L), rep(1.20, 4L))
  dec_het <- deconvolve_crossrun_hrf(
    bold, fit_h, lambda = lam_vec, scaling = "noise_sd", penalty = "sobolev"
  )
  expect_equal(dec_het$scaling$run1$lambda, lam_vec)

  diag_res <- diagnose_crossrun_recovery(dec_het, events = events)
  expect_equal(nrow(diag_res$voxels), 2L * V)
  expect_true("lambda" %in% names(diag_res$summary))
  expect_true(all(is.finite(diag_res$summary$mean_gain)))
  # Low-lambda strata have higher raw amplitude gain than heavily damped strata
  expect_gt(
    mean(diag_res$summary$mean_gain[diag_res$summary$lambda < 0.2]),
    mean(diag_res$summary$mean_gain[diag_res$summary$lambda > 0.5])
  )

  # Step 3 voxelwise: exact operator matching restores unbiased slope in BOTH low- and high-lambda strata
  reg_vox <- decon_epoch_regression(
    dec_het, covariates = covs, events = events, target_condition = "cue",
    lags = 0, restore_scale = TRUE, tau2_shrinkage = "none"
  )
  b_low <- mean(reg_vox$coefficients$beta_gls[1:4])
  b_high <- mean(reg_vox$coefficients$beta_gls[5:8])
  expect_lt(abs(b_low - beta_true), 0.15)
  expect_lt(abs(b_high - beta_true), 0.15)
  expect_lt(abs(b_low - b_high), 0.10)

  # Step 3 parcel pre-aggregation mixing heterogeneous lambda sub-strata
  p_groups <- stats::setNames(c(rep("P_mix", 6L), rep("P_hi", 2L)), vox_ids)
  reg_par <- decon_epoch_regression(
    dec_het, covariates = covs, events = events, target_condition = "cue",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, parcel_groups = p_groups
  )
  expect_equal(nrow(reg_par$coefficients), 2L)
  expect_true(all(abs(reg_par$coefficients$beta_gls - beta_true) < 0.18))
})

test_that("Track 2b: bilateral Graph-EB shrinks tau2_bio accurately, respects sulcal/parcel borders and hotspots, and preserves marginal covariate invariance", {
  set.seed(20261116L)
  tr <- 1.0
  n_time <- 175L
  n_trials <- 30L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)

  # Two cortical banks separated by a 1-voxel CSF gap (i = 4 omitted from mask):
  # Bank A (i = 1..3, j = 1..5, k = 1..2): 30 voxels, true biological variance tau2 = 1.40
  # Bank B (i = 5..7, j = 1..5, k = 1..2): 30 voxels, pure null tau2 = 0.0 except 1 hotspot voxel (tau2 = 8.0)
  gA <- expand.grid(i = 1:3, j = 1:5, k = 1:2)
  gB <- expand.grid(i = 5:7, j = 1:5, k = 1:2)
  grid <- rbind(gA, gB)
  V <- nrow(grid)
  vox_ids <- paste0("v", seq_len(V))
  vm <- data.frame(
    voxel = seq_len(V),
    i = grid$i, j = grid$j, k = grid$k,
    x = 3.0 * grid$i, y = 3.0 * grid$j, z = 3.0 * grid$k
  )
  geom <- list(dim = c(7L, 5L, 2L), spacing = c(3.0, 3.0, 3.0))

  bankA_idx <- 1:30
  bankB_idx <- 31:60
  hotspot_v <- 45L
  bankB_null <- setdiff(bankB_idx, hotspot_v)

  tau2_true <- numeric(V)
  tau2_true[bankA_idx] <- 2.00
  tau2_true[hotspot_v] <- 16.00

  beta_true <- numeric(V)
  beta_true[bankA_idx] <- 0.85

  bold <- bold_unstruct <- events <- covs <- vector("list", 2L)
  names(bold) <- names(bold_unstruct) <- names(events) <- names(covs) <- c("run1", "run2")

  for (r in 1:2) {
    onsets <- cumsum(c(6.0, stats::runif(n_trials - 1L, 4.0, 6.2)))
    x <- as.numeric(scale(stats::rnorm(n_trials)))
    events[[r]] <- data.frame(onset = onsets, duration = 0.0, condition = "stim")
    covs[[r]] <- data.frame(x = x)
    G <- hrf_event_design_cpp(
      n_time, tr, 0.0, onsets, rep(0.0, n_trials), seq_len(n_trials), n_trials, lib
    )[[1L]]

    Y_r <- matrix(0.0, n_time, V, dimnames = list(NULL, vox_ids))
    # Shared bank-level component + independent trial biological variability
    bio_shared_A <- stats::rnorm(n_trials, sd = sqrt(0.70))
    for (v in seq_len(V)) {
      bio_v <- if (v %in% bankA_idx) {
        bio_shared_A + stats::rnorm(n_trials, sd = sqrt(max(0.0, tau2_true[v] - 0.70)))
      } else if (v == hotspot_v) {
        stats::rnorm(n_trials, sd = sqrt(tau2_true[v]))
      } else {
        numeric(n_trials)
      }
      amps_v <- (if (v %in% bankA_idx || v == hotspot_v) 2.5 else 0.0) + beta_true[v] * x + bio_v
      eps_v <- as.numeric(stats::arima.sim(list(ar = 0.30), n = n_time, sd = 0.45))
      Y_r[, v] <- as.numeric(G %*% amps_v) + eps_v
    }
    bold_unstruct[[r]] <- Y_r
    attr(Y_r, "voxel_map") <- vm
    attr(Y_r, "geometry") <- geom
    bold[[r]] <- Y_r
  }

  # Even with Step 1 default pooling = 'voxel', NIfTI/voxel_map triggers 'graph_eb' in Step 3 'auto'
  fit_map <- estimate_crossrun_hrf(bold, events, lib, tr = tr, pooling = "voxel", noise_model = "ar1")
  dec_map <- deconvolve_crossrun_hrf(bold, fit_map, lambda = 0.4, scaling = "noise_sd", penalty = "sobolev")

  reg_auto <- decon_epoch_regression(
    dec_map, covariates = covs, events = events, target_condition = "stim",
    lags = 0, restore_scale = TRUE, tau2_shrinkage = "auto"
  )
  reg_none <- decon_epoch_regression(
    dec_map, covariates = covs, events = events, target_condition = "stim",
    lags = 0, restore_scale = TRUE, tau2_shrinkage = "none"
  )
  reg_glob <- decon_epoch_regression(
    dec_map, covariates = covs, events = events, target_condition = "stim",
    lags = 0, restore_scale = TRUE, tau2_shrinkage = "eb"
  )

  expect_equal(reg_auto$settings$tau2_shrinkage_used, "graph_eb")
  expect_equal(reg_none$settings$tau2_shrinkage_used, "none")
  expect_equal(reg_glob$settings$tau2_shrinkage_used, "eb")

  tau_ge <- reg_auto$coefficients$tau2_bio
  tau_no <- reg_none$coefficients$tau2_bio
  tau_gl <- reg_glob$coefficients$tau2_bio

  # 1. Substantial RMSE reduction in biological trial variance across the map
  rmse_ge <- sqrt(mean((tau_ge - tau2_true)^2))
  rmse_no <- sqrt(mean((tau_no - tau2_true)^2))
  expect_lt(rmse_ge, 0.75 * rmse_no)

  # 2. Eliminates false zero-truncation inside true-signal Bank A
  expect_equal(sum(tau_ge[bankA_idx] == 0), 0L)

  # 3. Geodesic CSF gap + bilateral stopping prevents Bank A from bleeding into null Bank B
  #    (whereas non-spatial global EB pulls null Bank B upward toward global mean!)
  expect_lt(stats::median(tau_ge[bankB_null]), 0.15)
  expect_lt(stats::median(tau_ge[bankB_null]), stats::median(tau_gl[bankB_null]))

  # 4. Efron-Morris limited translation protects the isolated outlier hotspot in Bank B
  expect_gt(tau_ge[hotspot_v], 3.0)

  # 5. Strict covariate invariance of tau2_mode = 'marginal' under Graph-EB and Global-EB
  reg_nocov <- decon_epoch_regression(
    dec_map, covariates = NULL, events = events, target_condition = "stim",
    lags = 0, restore_scale = TRUE, tau2_shrinkage = "graph_eb", return_panel = TRUE
  )
  reg_marg <- decon_epoch_regression(
    dec_map, covariates = covs, events = events, target_condition = "stim",
    lags = 0, restore_scale = TRUE, tau2_mode = "marginal",
    tau2_shrinkage = "graph_eb", return_panel = TRUE
  )
  expect_equal(reg_nocov$coefficients$tau2_bio, reg_marg$coefficients$tau2_bio, tolerance = 1e-12)
  expect_equal(reg_nocov$panel$precision_weight, reg_marg$panel$precision_weight, tolerance = 1e-12)

  # 6. Explicit spatial_graph with parcel_boundary = 'hard' isolates touching parcels completely
  vm_touch <- vm
  vm_touch$i[bankB_idx] <- vm_touch$i[bankB_idx] - 1L
  vm_touch$x[bankB_idx] <- 3.0 * vm_touch$i[bankB_idx]
  grp_touch <- c(rep("ParcelA", length(bankA_idx)), rep("ParcelB", length(bankB_idx)))
  sg_hard <- build_spatial_graph(
    list(voxel_map = vm_touch, geometry = list(dim = c(6L, 5L, 2L), spacing = c(3.0, 3.0, 3.0))),
    spatial_groups = grp_touch,
    parcel_boundary = "hard",
    voxel_ids = vox_ids
  )
  reg_hard <- decon_epoch_regression(
    dec_map, covariates = covs, events = events, target_condition = "stim",
    lags = 0, restore_scale = TRUE, tau2_shrinkage = "graph_eb", spatial_graph = sg_hard
  )
  # With no CSF gap, 'hard' parcel boundary yields the exact same isolated Bank A/B estimates as the CSF gap!
  expect_equal(reg_hard$coefficients$tau2_bio, reg_auto$coefficients$tau2_bio, tolerance = 1e-12)

  # 7. Unstructured matrix without coordinates defaults to 'none' in 'auto' and errors cleanly on 'graph_eb'
  fit_un <- estimate_crossrun_hrf(bold_unstruct, events, lib, tr = tr, pooling = "voxel", noise_model = "ar1")
  dec_un <- deconvolve_crossrun_hrf(bold_unstruct, fit_un, lambda = 0.4, scaling = "noise_sd", penalty = "sobolev")
  reg_un_auto <- decon_epoch_regression(dec_un, covariates = covs, events = events, lags = 0)
  reg_un_none <- decon_epoch_regression(dec_un, covariates = covs, events = events, lags = 0, tau2_shrinkage = "none")
  expect_equal(reg_un_auto$settings$tau2_shrinkage_used, "none")
  expect_equal(reg_un_auto$coefficients$tau2_bio, reg_un_none$coefficients$tau2_bio, tolerance = 1e-14)
  expect_equal(reg_un_auto$coefficients$se_gls, reg_un_none$coefficients$se_gls, tolerance = 1e-14)
  expect_error(
    decon_epoch_regression(dec_un, covariates = covs, events = events, lags = 0, tau2_shrinkage = "graph_eb"),
    "requires spatial coordinates"
  )
})

test_that("Track 2b: Isserlis spatial-correlation accounting maintains calibrated Type I error on spatially smooth 3D AR(2) null maps", {
  set.seed(20261117L)
  tr <- 1.0
  n_time <- 150L
  n_trials <- 26L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)

  # 6 x 6 x 3 volumetric patch (108 voxels) with 6 mm FWHM Gaussian spatial smoothing of AR(2) noise
  grid <- expand.grid(i = 1:6, j = 1:6, k = 1:3)
  V <- nrow(grid)
  vox_ids <- paste0("v", seq_len(V))
  coords <- 3.0 * as.matrix(grid)
  vm <- data.frame(
    voxel = seq_len(V),
    i = grid$i, j = grid$j, k = grid$k,
    x = coords[, 1L], y = coords[, 2L], z = coords[, 3L]
  )
  geom <- list(dim = c(6L, 6L, 3L), spacing = c(3.0, 3.0, 3.0))

  # Spatial smoothing matrix for 6 mm FWHM (sigma = 2.55 mm)
  Dmat <- as.matrix(stats::dist(coords))
  sigma_sm <- 6.0 / sqrt(8.0 * log(2.0))
  K_sm <- exp(-(Dmat^2) / (2.0 * sigma_sm^2))
  K_sm <- sweep(K_sm, 1L, sqrt(rowSums(K_sm^2)), "/")

  p_vals_ge <- numeric(0L)
  max_t_ge <- numeric(4L)
  max_t_no <- numeric(4L)

  for (rep_i in 1:4) {
    bold <- events <- covs <- vector("list", 2L)
    names(bold) <- names(events) <- names(covs) <- c("run1", "run2")
    for (r in 1:2) {
      onsets <- cumsum(c(6.0, stats::runif(n_trials - 1L, 4.0, 6.0)))
      x <- as.numeric(scale(stats::rnorm(n_trials)))
      events[[r]] <- data.frame(onset = onsets, duration = 0.0, condition = "cue")
      covs[[r]] <- data.frame(x = x)

      E_ind <- matrix(0.0, n_time, V)
      for (v in seq_len(V)) {
        E_ind[, v] <- as.numeric(stats::arima.sim(list(ar = c(0.40, -0.12)), n = n_time, sd = 1.0))
      }
      Y_r <- E_ind %*% t(K_sm)
      colnames(Y_r) <- vox_ids
      attr(Y_r, "voxel_map") <- vm
      attr(Y_r, "geometry") <- geom
      bold[[r]] <- Y_r
    }

    fit_h <- estimate_crossrun_hrf(bold, events, lib, tr = tr, pooling = "graph", noise_model = "ar1")
    dec_h <- deconvolve_crossrun_hrf(bold, fit_h, lambda = NULL, scaling = "noise_sd", penalty = "sobolev")
    reg_ge <- decon_epoch_regression(
      dec_h, covariates = covs, events = events, lags = 0,
      weighting = "gls", restore_scale = TRUE, tau2_shrinkage = "graph_eb"
    )
    reg_no <- decon_epoch_regression(
      dec_h, covariates = covs, events = events, lags = 0,
      weighting = "gls", restore_scale = TRUE, tau2_shrinkage = "none"
    )
    p_vals_ge <- c(p_vals_ge, reg_ge$coefficients$p_value)
    max_t_ge[rep_i] <- max(abs(reg_ge$coefficients$t_gls))
    max_t_no[rep_i] <- max(abs(reg_no$coefficients$t_gls))
  }

  fpr_ge <- mean(p_vals_ge < 0.05)
  expect_gte(fpr_ge, 0.015)
  expect_lte(fpr_ge, 0.085)
  # Scale moderation suppresses extreme spurious |t| spikes from single-voxel variance underestimation
  expect_lte(mean(max_t_ge), mean(max_t_no) + 0.15)
})

test_that("settings$lambda_by_voxel exports cleanly via write_nifti_ts and parcel-mode graph_eb auto-builds centroid graph", {
  set.seed(20261118L)
  tmp_dir <- withr::local_tempdir()
  tr <- 1.0
  n_time <- 140L
  n_trials <- 24L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)

  # 8 parcels (3 voxels each = 24 voxels): parcels 1..4 on active gyral cluster, parcels 5..8 distant null cluster
  grid <- expand.grid(i = 1:6, j = 1:4, k = 1L)
  V <- nrow(grid)
  vox_ids <- paste0("v", seq_len(V))
  parcel_num <- ((grid$j - 1L) * 2L) + ifelse(grid$i <= 3L, 1L, 2L) # 8 parcels of 3 voxels each
  # Separate odd vs even parcels across a wide 40 mm gap so centroids form two distinct 4-parcel clusters
  x_mm <- ifelse(grid$i <= 3L, 3.0 * grid$i, 40.0 + 3.0 * grid$i)
  y_mm <- 3.0 * grid$j
  z_mm <- rep(0.0, V)

  mask_arr <- array(0L, dim = c(6L, 4L, 1L))
  for (v in seq_len(V)) mask_arr[grid$i[v], grid$j[v], 1L] <- parcel_num[v]
  f_atlas <- file.path(tmp_dir, "atlas.nii.gz")
  RNifti::writeNifti(RNifti::asNifti(mask_arr, pixdim = c(1, 3, 3, 3)), f_atlas)

  vm <- data.frame(
    voxel = vox_ids,
    linear_index = seq_len(V),
    i = grid$i, j = grid$j, k = grid$k,
    x = x_mm, y = y_mm, z = z_mm,
    atlas_value = parcel_num,
    stringsAsFactors = FALSE
  )
  geom <- medusa:::resolve_nifti_spatial_geometry(f_atlas, label = "atlas")

  active_parcels <- c(1L, 3L, 5L, 7L) # x <= 9 mm cluster
  bold <- events <- covs <- vector("list", 2L)
  names(bold) <- names(events) <- names(covs) <- c("run1", "run2")
  for (r in 1:2) {
    onsets <- cumsum(c(6.0, stats::runif(n_trials - 1L, 4.0, 5.8)))
    x_cov <- as.numeric(scale(stats::rnorm(n_trials)))
    events[[r]] <- data.frame(onset = onsets, duration = 0.0, condition = "stim")
    covs[[r]] <- data.frame(x = x_cov)
    G <- hrf_event_design_cpp(
      n_time, tr, 0.0, onsets, rep(0.0, n_trials), seq_len(n_trials), n_trials, lib
    )[[1L]]
    Y_r <- matrix(0.0, n_time, V, dimnames = list(NULL, vox_ids))
    for (p_id in 1:8) {
      bio_p <- if (p_id %in% active_parcels) stats::rnorm(n_trials, sd = sqrt(1.8)) else numeric(n_trials)
      amp_p <- (if (p_id %in% active_parcels) 3.0 + 0.9 * x_cov else 0.0) + bio_p
      for (v in which(parcel_num == p_id)) {
        Y_r[, v] <- as.numeric(G %*% amp_p) + as.numeric(stats::arima.sim(list(ar = 0.25), n = n_time, sd = 0.5))
      }
    }
    attr(Y_r, "voxel_map") <- vm
    attr(Y_r, "geometry") <- geom
    bold[[r]] <- Y_r
  }

  p_labels <- stats::setNames(as.character(parcel_num), vox_ids)
  fit <- estimate_crossrun_hrf(bold, events, lib, tr = tr, spatial_groups = p_labels, pooling = "parcel", noise_model = "ar1")
  dec <- deconvolve_crossrun_hrf(bold, fit, lambda = NULL, lambda_mode = "stratified", scaling = "noise_sd", penalty = "sobolev")

  # Polish #4: settings$lambda_by_voxel is named by voxel_ids, matches scaling[[1]]$lambda, and writes 3D NIfTI in 1 line
  expect_identical(names(dec$settings$lambda_by_voxel), vox_ids)
  expect_equal(unname(dec$settings$lambda_by_voxel), dec$scaling$run1$lambda)
  f_lam <- file.path(tmp_dir, "log10_lambda.nii.gz")
  write_nifti_ts(log10(dec$settings$lambda_by_voxel), f_lam, reference = dec, datatype = "double")
  expect_equal(as.numeric(RNifti::readNifti(f_lam)), log10(unname(dec$settings$lambda_by_voxel)), tolerance = 1e-10)

  # Polish #5: parcel-mode tau2_shrinkage = 'auto' stays 'none', while 'graph_eb' auto-builds parcel-centroid graph
  reg_p_auto <- decon_epoch_regression(
    dec, covariates = covs, events = events, lags = 0,
    restore_scale = TRUE, parcel_groups = p_labels, tau2_shrinkage = "auto"
  )
  reg_p_ge <- decon_epoch_regression(
    dec, covariates = covs, events = events, lags = 0,
    restore_scale = TRUE, parcel_groups = p_labels, tau2_shrinkage = "graph_eb"
  )
  expect_equal(reg_p_auto$settings$tau2_shrinkage_used, "none")
  expect_equal(reg_p_ge$settings$tau2_shrinkage_used, "graph_eb")
  expect_equal(nrow(reg_p_ge$coefficients), 8L)

  act_rows <- match(as.character(active_parcels), reg_p_ge$coefficients$parcel)
  null_rows <- setdiff(seq_len(8L), act_rows)
  expect_gt(stats::median(reg_p_ge$coefficients$tau2_bio[act_rows]), 0.6)
  expect_lt(stats::median(reg_p_ge$coefficients$tau2_bio[null_rows]), 0.25)
})
