test_that("Tier 1.1: highpass DCT drift + preproject nuisance_mode annihilates drift and recovers PPI", {
  set.seed(20261101)
  tr <- 1.0
  n_time <- 180L
  n_trials <- 26L
  lib <- hrf_library(tr, peak_delays = c(5, 6, 7), dispersions = c(0.8, 1.0, 1.2), undershoot_ratios = 6)
  h_can <- lib[, which(abs(attr(lib, "parameters")$peak_delay - 6) < 1e-6 &
                       abs(attr(lib, "parameters")$dispersion - 1) < 1e-6)[1L], drop = FALSE]

  # 1. Operator-level exact drift & nuisance annihilation and C++ vs R agreement
  rows <- seq_len(n_time)
  u_nuis <- cbind(mot1 = rnorm(n_time), mot2 = rnorm(n_time))
  z_pre <- medusa:::hrf_build_decon_nuisance_matrix(
    n_obs = n_time, tr = tr, rows = rows, nuisance_r = u_nuis,
    nuisance_mode = "preproject", highpass_sec = 128.0
  )
  expect_gt(ncol(z_pre), 2L + ncol(u_nuis)) # includes intercept, trend, DCT(<1/128Hz), and mot1/mot2

  op_cpp <- medusa:::hrf_inverse_operator(
    n_time, h_can[, 1L], rows, 0.35, z_pre,
    settings = list(
      lambda = 1.0, penalty = "sobolev", ridge_floor = 1e-6,
      backend = "cpp", tail_padding = 0L, nuisance_mode = "preproject"
    )
  )
  op_r <- medusa:::hrf_inverse_operator(
    n_time, h_can[, 1L], rows, 0.35, z_pre,
    settings = list(
      lambda = 1.0, penalty = "sobolev", ridge_floor = 1e-6,
      backend = "r", tail_padding = 0L, nuisance_mode = "preproject"
    )
  )
  expect_equal(op_cpp, op_r, tolerance = 1e-10)
  expect_lt(max(abs(op_cpp %*% z_pre)), 1e-11)

  # 2. Modern vs legacy default preservation
  ev_simple <- list(
    run1 = data.frame(onset = seq(8, 150, length.out = 20), duration = 0, condition = "stim"),
    run2 = data.frame(onset = seq(8, 150, length.out = 20), duration = 0, condition = "stim")
  )
  y_simple <- list(
    run1 = matrix(rnorm(n_time * 4), n_time, 4, dimnames = list(NULL, paste0("v", 1:4))),
    run2 = matrix(rnorm(n_time * 4), n_time, 4, dimnames = list(NULL, paste0("v", 1:4)))
  )
  fit_mod <- estimate_crossrun_hrf(y_simple, ev_simple, tr = tr, kernel_grid = lib)
  dec_mod <- deconvolve_crossrun_hrf(y_simple, fit_mod, lambda = 1.0)
  expect_equal(fit_mod$settings$highpass, 128.0)
  expect_equal(dec_mod$settings$nuisance_mode, "preproject")

  fit_leg <- estimate_crossrun_hrf(y_simple, ev_simple, tr = tr, kernel_grid = lib,
                                   noise_model = "ols", evidence = "legacy")
  dec_leg <- deconvolve_crossrun_hrf(y_simple, fit_leg, lambda = 1.0)
  expect_true(is.infinite(fit_leg$settings$highpass))
  expect_equal(dec_leg$settings$nuisance_mode, "ignore")

  # 3. PPI slope recovery under strong low-frequency scanner drift
  n_subj <- 8L
  v_groups <- stats::setNames(rep(c("R_seed", "R_targ"), each = 3L), paste0("v", 1:6))
  t_idx <- seq(0, n_time - 1L)
  dct_low <- cbind(
    cos(pi * (2 * t_idx + 1) * 1 / (2 * n_time)),
    cos(pi * (2 * t_idx + 1) * 2 / (2 * n_time)),
    cos(pi * (2 * t_idx + 1) * 3 / (2 * n_time))
  )
  dct_low <- scale(dct_low)

  fits_pre <- vector("list", n_subj)
  fits_ign <- vector("list", n_subj)
  names(fits_pre) <- names(fits_ign) <- sprintf("sub%02d", seq_len(n_subj))

  for (s in seq_len(n_subj)) {
    bold_s <- ev_s <- cov_s <- list()
    for (r in 1:2) {
      rn <- paste0("run", r)
      ons <- cumsum(c(6.0, runif(n_trials - 1L, 4.2, 6.8)))
      xv <- as.numeric(scale(rnorm(n_trials)))
      eta_s <- rnorm(n_trials, 0, 0.70)
      eta_t <- (0.40 + 0.38 * xv) * eta_s + rnorm(n_trials, 0, 0.35)
      G_tr <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons, rep(0, n_trials), seq_len(n_trials), n_trials, h_can)[[1L]]
      s_clean <- as.numeric(G_tr %*% (0.8 + 0.3 * xv + eta_s))
      t_clean <- as.numeric(G_tr %*% (0.8 + 0.3 * xv + eta_t))
      drift_s <- as.numeric(dct_low %*% rnorm(3, 0, 0.65))
      drift_t <- as.numeric(dct_low %*% rnorm(3, 0, 0.65))
      ym <- matrix(0, n_time, 6L, dimnames = list(NULL, names(v_groups)))
      for (k in 1:3) ym[, k] <- s_clean + drift_s + as.numeric(stats::filter(rnorm(n_time, 0, 0.16), 0.35, method = "recursive"))
      for (k in 4:6) ym[, k] <- t_clean + drift_t + as.numeric(stats::filter(rnorm(n_time, 0, 0.16), 0.35, method = "recursive"))
      bold_s[[rn]] <- ym
      ev_s[[rn]] <- data.frame(onset = ons, duration = 0, condition = "fb", trial = seq_len(n_trials))
      cov_s[[rn]] <- data.frame(x = xv)
    }
    hf_pre <- estimate_crossrun_hrf(bold_s, ev_s, tr = tr, kernel_grid = lib, spatial_groups = v_groups, pooling = "parcel")
    dc_pre <- deconvolve_crossrun_hrf(bold_s, hf_pre, scaling = "noise_sd", penalty = "sobolev", nuisance_mode = "preproject")
    er_pre <- decon_epoch_regression(dc_pre, covariates = cov_s, target_condition = "fb", lags = 0,
                                     restore_scale = TRUE, aggregate_by_parcel = TRUE, parcel_groups = v_groups,
                                     whiten_panel = FALSE, return_panel = TRUE)
    er_pre$panel$id <- names(fits_pre)[s]
    fits_pre[[s]] <- er_pre

    hf_ign <- estimate_crossrun_hrf(bold_s, ev_s, tr = tr, kernel_grid = lib, spatial_groups = v_groups, pooling = "parcel", highpass = Inf)
    dc_ign <- deconvolve_crossrun_hrf(bold_s, hf_ign, scaling = "noise_sd", penalty = "sobolev", nuisance_mode = "ignore")
    er_ign <- decon_epoch_regression(dc_ign, covariates = cov_s, target_condition = "fb", lags = 0,
                                     restore_scale = TRUE, aggregate_by_parcel = TRUE, parcel_groups = v_groups,
                                     whiten_panel = FALSE, return_panel = TRUE)
    er_ign$panel$id <- names(fits_ign)[s]
    fits_ign[[s]] <- er_ign
  }

  ppi_pre <- decon_connectivity_repair(fits_pre, "R_seed", "R_targ", 0, 0, test_ppi = TRUE, method = "mmoms")$group
  ppi_ign <- decon_connectivity_repair(fits_ign, "R_seed", "R_targ", 0, 0, test_ppi = TRUE, method = "mmoms")$group
  b_pre <- ppi_pre[term == "seed_c:x_c"]$estimate[1L]
  b_ign <- ppi_ign[term == "seed_c:x_c"]$estimate[1L]
  expect_gt(b_pre, 0.24)
  expect_gt(abs(b_pre) - abs(b_ign), 0.10)
})

test_that("Tier 1.6.1 & 1.6.2: whole-brain (>4096 voxels) rwls screening and safe SVD fallback succeed", {
  set.seed(20261102)
  tr <- 1.0
  n_time <- 90L
  n_vox <- 4300L
  lib <- hrf_library(tr, peak_delays = c(5, 6, 7), dispersions = c(1.0), undershoot_ratios = 6)
  h_can <- lib[, 2L, drop = FALSE]
  vox_ids <- paste0("v", seq_len(n_vox))

  bold <- events <- list()
  spike_frames <- c(28L, 62L)
  for (r in 1:2) {
    rn <- paste0("run", r)
    ons <- seq(8, 76, by = 8.5)
    G <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons, rep(0, length(ons)), rep(1L, length(ons)), 1L, h_can)[[1L]]
    ym <- matrix(rnorm(n_time * n_vox, 0, 0.25), n_time, n_vox, dimnames = list(NULL, vox_ids))
    # Place task signal + massive motion artifact strictly in voxels > 4096
    act_cols <- 4150L:4300L
    ym[, act_cols] <- ym[, act_cols] + matrix(2.2 * as.numeric(G), n_time, length(act_cols))
    ym[spike_frames, ] <- ym[spike_frames, ] + 4.5
    bold[[rn]] <- ym
    events[[rn]] <- data.frame(onset = ons, duration = 0, condition = "task")
  }

  fit <- estimate_crossrun_hrf(bold, events, tr = tr, kernel_grid = lib, min_parcel_voxels = 3L)
  expect_true("selected_adjusted_partial_r2" %in% names(fit$selections[[1L]]))
  dec <- deconvolve_crossrun_hrf(bold, fit, lambda = 1.0, rwls = "woodbury")
  expect_lt(min(dec$frame_weights[[1L]][spike_frames]), 0.45)

  # Safe SVD helper handles clustered singular values and extreme observation weights cleanly
  M_deg <- diag(c(rep(1.0, 8L), rep(1e-15, 4L)))
  sv_safe <- medusa:::hrf_safe_svd(M_deg, nu = 12L, nv = 12L)
  expect_equal(length(sv_safe$d), 12L)
  expect_true(all(is.finite(sv_safe$u)) && all(is.finite(sv_safe$v)))

  w_ext <- rep(1.0, n_time)
  w_ext[c(15L, 16L, 45L)] <- 1e-4
  op_w <- medusa:::hrf_inverse_operator_cpp(
    n_time, as.numeric(h_can), seq_len(n_time), 0.4,
    matrix(0, n_time, 0L), 0.8, "sobolev", 1e-6, 0L,
    NULL, NULL, 1.0, w_ext, "woodbury"
  )
  expect_true(all(is.finite(op_w)))
})

test_that("Tier 1.3: variable-duration source model avoids collateral self-subtraction in contiguous WM trials", {
  set.seed(20261103)
  tr <- 1.0
  dt_fine <- 0.05
  n_trials <- 28L
  lib <- hrf_library(tr, peak_delays = c(5, 6, 7), dispersions = c(0.8, 1.0, 1.2), undershoot_ratios = 6)
  h_can <- lib[, which(abs(attr(lib, "parameters")$peak_delay - 6) < 1e-6 &
                       abs(attr(lib, "parameters")$dispersion - 1) < 1e-6)[1L], drop = FALSE]

  # Contiguous encode (0.8s) -> variable delay (2.5..8.0s) -> probe (1.0s) with 0s inter-phase gap
  bold <- events <- covs <- durs <- list()
  true_slope <- 0.45
  for (r in 1:2) {
    rn <- paste0("run", r)
    del_d <- runif(n_trials, 2.5, 8.0)
    itis <- runif(n_trials, 3.0, 6.0)
    ons_enc <- numeric(n_trials)
    ons_del <- numeric(n_trials)
    ons_prb <- numeric(n_trials)
    t_cur <- 6.0
    for (i in seq_len(n_trials)) {
      ons_enc[i] <- t_cur
      ons_del[i] <- t_cur + 0.8
      ons_prb[i] <- ons_del[i] + del_d[i]
      t_cur <- ons_prb[i] + 1.0 + itis[i]
    }
    n_time <- as.integer(ceiling(t_cur + 18.0))
    load_x <- as.numeric(scale(sample(rep(c(-1.5, -0.5, 0.5, 1.5), length.out = n_trials))))
    amp_del <- 0.55 + true_slope * load_x + rnorm(n_trials, 0, 0.12)

    ev_df <- rbind(
      data.frame(onset = ons_enc, duration = 0.8, condition = "encode"),
      data.frame(onset = ons_del, duration = del_d, condition = "delay"),
      data.frame(onset = ons_prb, duration = 1.0, condition = "probe")
    )
    ev_df <- ev_df[order(ev_df$onset), ]
    rownames(ev_df) <- NULL

    G_enc <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons_enc, rep(0.8, n_trials), rep(1L, n_trials), 1L, h_can)[[1L]]
    G_del <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons_del, del_d, seq_len(n_trials), n_trials, h_can)[[1L]]
    G_prb <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons_prb, rep(1.0, n_trials), rep(1L, n_trials), 1L, h_can)[[1L]]
    sig_clean <- as.numeric(0.70 * G_enc + G_del %*% amp_del + 0.65 * G_prb)

    ym <- matrix(0, n_time, 4L, dimnames = list(NULL, paste0("v", 1:4)))
    for (v in 1:4) {
      ym[, v] <- sig_clean + as.numeric(stats::filter(rnorm(n_time, 0, 0.14), 0.30, method = "recursive"))
    }
    bold[[rn]] <- ym
    events[[rn]] <- ev_df
    covs[[rn]] <- data.frame(load_z = load_x)
    durs[[rn]] <- del_d
  }

  v_grp <- stats::setNames(rep("dlPFC", 4L), paste0("v", 1:4))
  hf <- estimate_crossrun_hrf(bold, events, tr = tr, kernel_grid = lib, spatial_groups = v_grp, pooling = "parcel")
  dc <- deconvolve_crossrun_hrf(bold, hf, scaling = "noise_sd", penalty = "sobolev")
  er <- decon_epoch_regression(
    dc, covariates = covs, events = events, target_condition = "delay", lags = 0,
    source_model = "duration", trial_durations = durs, duration_scale = "rate",
    restore_scale = TRUE, aggregate_by_parcel = TRUE, parcel_groups = v_grp,
    return_panel = TRUE
  )
  cf_del <- er$coefficients[er$coefficients$term == "load_z", ]
  expect_equal(cf_del$beta_gls[1L], true_slope, tolerance = 0.10)
  emp_sd_resid <- stats::sd(stats::resid(stats::lm(decon_demixed ~ load_z, data = er$panel)))
  med_se_demix <- stats::median(er$panel$se_demixed)
  expect_lt(med_se_demix, 1.35 * emp_sd_resid)
})

test_that("Tier 1.4: other_event_model = 'auto' suppresses trial-varying cue confound while keeping constant-cue efficiency", {
  set.seed(20261104)
  tr <- 1.0
  n_time <- 185L
  n_trials <- 28L
  lib <- hrf_library(tr, peak_delays = c(6.0), dispersions = c(1.0), undershoot_ratios = 6)
  h_can <- lib[, 1L, drop = FALSE]
  v_grp <- stats::setNames(rep("R_test", 4L), paste0("v", 1:4))

  sim_cue_fb <- function(cue_mod_sd) {
    bold <- events <- covs <- list()
    for (r in 1:2) {
      rn <- paste0("run", r)
      ons_cue <- cumsum(c(6.0, runif(n_trials - 1L, 4.8, 6.5)))
      isi <- runif(n_trials, 1.5, 2.4)
      ons_fb <- ons_cue + isi
      xv <- as.numeric(scale(rnorm(n_trials)))
      cue_amp <- if (cue_mod_sd > 0) {
        1.2 + cue_mod_sd * (0.70 * xv + sqrt(1 - 0.70^2) * rnorm(n_trials))
      } else {
        rep(1.2, n_trials)
      }
      fb_amp <- if (cue_mod_sd > 0) rnorm(n_trials, 0.5, 0.15) else (0.5 + 0.45 * xv + rnorm(n_trials, 0, 0.15))

      G_cue <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons_cue, rep(0, n_trials), seq_len(n_trials), n_trials, h_can)[[1L]]
      G_fb <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons_fb, rep(0, n_trials), seq_len(n_trials), n_trials, h_can)[[1L]]
      sig <- as.numeric(G_cue %*% cue_amp + G_fb %*% fb_amp)
      ym <- matrix(0, n_time, 4L, dimnames = list(NULL, names(v_grp)))
      for (v in 1:4) ym[, v] <- sig + as.numeric(stats::filter(rnorm(n_time, 0, 0.14), 0.30, method = "recursive"))

      ev_df <- rbind(
        data.frame(onset = ons_cue, duration = 0, condition = "cue"),
        data.frame(onset = ons_fb, duration = 0, condition = "feedback")
      )
      ev_df <- ev_df[order(ev_df$onset), ]
      bold[[rn]] <- ym
      events[[rn]] <- ev_df
      covs[[rn]] <- data.frame(x = xv)
    }
    list(bold = bold, events = events, covs = covs)
  }

  # 1. Confounded trial-varying cue (true feedback x slope = 0)
  d_conf <- sim_cue_fb(cue_mod_sd = 0.85)
  hf1 <- estimate_crossrun_hrf(d_conf$bold, d_conf$events, tr = tr, kernel_grid = lib, spatial_groups = v_grp, pooling = "parcel")
  dc1 <- deconvolve_crossrun_hrf(d_conf$bold, hf1, lambda = 0.8, scaling = "noise_sd", penalty = "sobolev")
  er_cond <- decon_epoch_regression(dc1, covariates = d_conf$covs, target_condition = "feedback", lags = 0,
                                    restore_scale = TRUE, aggregate_by_parcel = TRUE, parcel_groups = v_grp,
                                    other_event_model = "condition")
  er_auto <- decon_epoch_regression(dc1, covariates = d_conf$covs, target_condition = "feedback", lags = 0,
                                    restore_scale = TRUE, aggregate_by_parcel = TRUE, parcel_groups = v_grp,
                                    other_event_model = "auto")
  b_cond <- er_cond$coefficients$beta_gls[er_cond$coefficients$term == "x"]
  b_auto <- er_auto$coefficients$beta_gls[er_auto$coefficients$term == "x"]
  expect_gt(b_cond, 0.18)
  expect_lt(abs(b_auto), 0.12)

  # 2. Constant cue (no cue modulation; auto preserves condition-mean efficiency)
  d_const <- sim_cue_fb(cue_mod_sd = 0.0)
  hf2 <- estimate_crossrun_hrf(d_const$bold, d_const$events, tr = tr, kernel_grid = lib, spatial_groups = v_grp, pooling = "parcel")
  dc2 <- deconvolve_crossrun_hrf(d_const$bold, hf2, lambda = 0.8, scaling = "noise_sd", penalty = "sobolev")
  er_auto2 <- decon_epoch_regression(dc2, covariates = d_const$covs, target_condition = "feedback", lags = 0,
                                     restore_scale = TRUE, aggregate_by_parcel = TRUE, parcel_groups = v_grp,
                                     other_event_model = "auto")
  er_trial2 <- decon_epoch_regression(dc2, covariates = d_const$covs, target_condition = "feedback", lags = 0,
                                      restore_scale = TRUE, aggregate_by_parcel = TRUE, parcel_groups = v_grp,
                                      other_event_model = "trial")
  se_auto2 <- er_auto2$coefficients$se_gls[er_auto2$coefficients$term == "x"]
  se_trial2 <- er_trial2$coefficients$se_gls[er_trial2$coefficients$term == "x"]
  expect_lt(se_auto2, se_trial2)
})

test_that("Tier 1.5: R_post within-trial deblurring + multilag_control isolate directed Seed(0s)->Target(4s) connectivity", {
  set.seed(20261105)
  tr <- 1.0
  n_time <- 260L
  n_trials <- 28L
  n_subj <- 20L
  lags_grid <- c(0, 2, 4)
  lib <- hrf_library(tr, peak_delays = c(5, 6, 7), dispersions = c(1.0), undershoot_ratios = 6)
  h_seed <- lib[, which(abs(attr(lib, "parameters")$peak_delay - 5) < 1e-6)[1L], drop = FALSE]
  h_targ <- lib[, which(abs(attr(lib, "parameters")$peak_delay - 6) < 1e-6)[1L], drop = FALSE]
  v_grp <- stats::setNames(rep(c("R_seed", "R_targ"), each = 3L), paste0("v", 1:6))

  fits_sub <- vector("list", n_subj)
  names(fits_sub) <- sprintf("sub%02d", seq_len(n_subj))

  for (s in seq_len(n_subj)) {
    bold_s <- ev_s <- cov_s <- list()
    for (r in 1:2) {
      rn <- paste0("run", r)
      ons_cue <- cumsum(c(6.0, runif(n_trials - 1L, 5.8, 8.8)))
      ons <- ons_cue + runif(n_trials, 1.5, 2.5)
      xv <- as.numeric(scale(rnorm(n_trials)))
      eta_s0 <- rnorm(n_trials, 0, 0.70)
      eta_s4 <- rnorm(n_trials, 0, 0.32)
      eta_t0 <- rnorm(n_trials, 0, 0.45)
      eta_t4 <- 0.55 * eta_s0 + rnorm(n_trials, 0, 0.35)

      G_sc <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons_cue, rep(0, n_trials), rep(1L, n_trials), 1L, h_seed)[[1L]]
      G_tc <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons_cue, rep(0, n_trials), rep(1L, n_trials), 1L, h_targ)[[1L]]
      G_s0 <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons, rep(0, n_trials), seq_len(n_trials), n_trials, h_seed)[[1L]]
      G_s4 <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons + 4.0, rep(0, n_trials), seq_len(n_trials), n_trials, h_seed)[[1L]]
      G_t0 <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons, rep(0, n_trials), seq_len(n_trials), n_trials, h_targ)[[1L]]
      G_t4 <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons + 4.0, rep(0, n_trials), seq_len(n_trials), n_trials, h_targ)[[1L]]

      sig_s <- as.numeric(1.2 * G_sc + G_s0 %*% (0.8 + 0.4 * xv + eta_s0) + G_s4 %*% eta_s4)
      sig_t <- as.numeric(0.9 * G_tc + G_t0 %*% (0.7 + eta_t0) + G_t4 %*% (0.4 + eta_t4))
      ym <- matrix(0, n_time, 6L, dimnames = list(NULL, names(v_grp)))
      for (k in 1:3) ym[, k] <- sig_s + as.numeric(stats::filter(rnorm(n_time, 0, 0.16), 0.35, method = "recursive"))
      for (k in 4:6) ym[, k] <- sig_t + as.numeric(stats::filter(rnorm(n_time, 0, 0.16), 0.35, method = "recursive"))
      bold_s[[rn]] <- ym
      ev_df <- rbind(
        data.frame(onset = ons_cue, duration = 0, condition = "cue"),
        data.frame(onset = ons, duration = 0, condition = "fb")
      )
      ev_df <- ev_df[order(ev_df$onset), ]
      ev_s[[rn]] <- ev_df
      cov_s[[rn]] <- data.frame(x = xv)
    }
    hf <- estimate_crossrun_hrf(bold_s, ev_s, tr = tr, kernel_grid = lib, spatial_groups = v_grp, pooling = "parcel", min_parcel_voxels = 3L)
    dc <- deconvolve_crossrun_hrf(bold_s, hf, scaling = "noise_sd", penalty = "sobolev")
    er <- decon_epoch_regression(
      dc, covariates = cov_s, target_condition = "fb", lags = lags_grid,
      demix_within_trial = TRUE, restore_scale = TRUE,
      aggregate_by_parcel = TRUE, parcel_groups = v_grp,
      whiten_panel = FALSE, return_panel = TRUE
    )
    er$panel$id <- names(fits_sub)[s]
    fits_sub[[s]] <- er
  }

  eval_pair <- function(sl, tl) {
    decon_connectivity_repair(
      fits_sub, "R_seed", "R_targ", seed_lag = sl, target_lag = tl,
      covariates = "x", multilag_control = TRUE, method = "mmoms"
    )$group[term == "seed_c"][1L]
  }

  g_04 <- eval_pair(0, 4)
  g_00 <- eval_pair(0, 0)
  g_02 <- eval_pair(0, 2)
  g_24 <- eval_pair(2, 4)
  g_40 <- eval_pair(4, 0)
  g_44 <- eval_pair(4, 4)

  expect_gt(g_04$estimate, 0.25)
  expect_lt(g_04$p.value, 0.01)
  for (g_null in list(g_00, g_02, g_24, g_40, g_44)) {
    expect_lt(abs(g_null$estimate), 0.18)
  }
})

test_that("Tier 1.6.3-1.6.5: poet_covreg j_rwls/marginal cov_list and medusa_between_subject_repair built-in fallback work", {
  set.seed(20261106)
  # 1. poet_covreg with moderate rwls weights (0.05 <= fw < 0.85) and tau2_mode = "marginal"
  tr <- 1.0
  n_time <- 110L
  lib <- hrf_library(tr, peak_delays = c(6.0), dispersions = c(1.0), undershoot_ratios = 6)
  h_can <- lib[, 1L, drop = FALSE]
  v_grp <- stats::setNames(rep(c("P1", "P2"), each = 3L), paste0("v", 1:6))

  bold <- events <- covs <- list()
  for (r in 1:2) {
    rn <- paste0("run", r)
    ons <- seq(8, 92, by = 6.0)
    nt <- length(ons)
    xv <- as.numeric(scale(rnorm(nt)))
    G <- medusa:::hrf_event_design_cpp(n_time, tr, 0, ons, rep(0, nt), seq_len(nt), nt, h_can)[[1L]]
    ym <- matrix(rnorm(n_time * 6L, 0, 0.25), n_time, 6L, dimnames = list(NULL, names(v_grp)))
    ym[, 1:3] <- ym[, 1:3] + as.numeric(G %*% (0.8 + 0.4 * xv))
    # Inject moderate transient jolt so rwls produces 0.05 <= fw < 0.85
    ym[c(25L, 65L), ] <- ym[c(25L, 65L), ] + rnorm(12L, 0, 1.8)
    bold[[rn]] <- ym
    events[[rn]] <- data.frame(onset = ons, duration = 0, condition = "fb")
    covs[[rn]] <- data.frame(x = xv)
  }
  hf <- estimate_crossrun_hrf(bold, events, tr = tr, kernel_grid = lib, spatial_groups = v_grp, pooling = "parcel")
  dc <- deconvolve_crossrun_hrf(bold, hf, lambda = 1.0, rwls = "woodbury")
  # Force sustained moderate frame weights (>= 12 frames across runs) to deterministically exercise j_rwls
  dc$frame_weights[[1L]][22L:29L] <- seq(0.22, 0.65, length.out = 8L)
  dc$frame_weights[[2L]][55L:62L] <- seq(0.25, 0.70, length.out = 8L)
  er_marg <- decon_epoch_regression(
    dc, covariates = covs, target_condition = "fb", lags = 0,
    tau2_mode = "marginal", aggregate_by_parcel = TRUE, parcel_groups = v_grp,
    shared_error_method = "poet_covreg", return_panel = TRUE
  )
  expect_false(is.null(er_marg$shared_error_model))
  expect_true(all(is.finite(er_marg$panel$se_demixed)))

  # 2. Built-in medusa_between_subject_repair fallback disattenuates Level-2 slopes
  n_s <- 55L
  n_tr <- 24L
  s_ids <- sprintf("s%03d", seq_len(n_s))
  b_vs_true <- rnorm(n_s, 0, 0.30)
  b_mp_true <- 1.25 * b_vs_true + rnorm(n_s, 0, 0.18)
  z_behav <- 1.40 * b_vs_true + rnorm(n_s, 0, 0.22)
  subj_df <- data.frame(id = s_ids, z_behav = z_behav, stringsAsFactors = FALSE)

  rows_l1 <- lapply(seq_len(n_s), function(i) {
    xv <- as.numeric(scale(rnorm(n_tr)))
    w_i <- runif(n_tr, 0.6, 1.4)
    w_i <- w_i / mean(w_i)
    y_vs <- (0.60 + b_vs_true[i]) * xv + rnorm(n_tr, 0, 0.75 / sqrt(w_i))
    y_mp <- (0.50 + b_mp_true[i]) * xv + rnorm(n_tr, 0, 0.70 / sqrt(w_i))
    data.frame(id = s_ids[i], rpe_whitened = xv, y_vs = y_vs, y_mp = y_mp, pw = w_i, stringsAsFactors = FALSE)
  })
  df_l1 <- do.call(rbind, rows_l1)
  m_vs <- suppressWarnings(lme4::lmer(y_vs ~ rpe_whitened + (1 + rpe_whitened | id), data = df_l1, weights = pw))
  m_mp <- suppressWarnings(lme4::lmer(y_mp ~ rpe_whitened + (1 + rpe_whitened | id), data = df_l1, weights = pw))

  rep_beh <- medusa_between_subject_repair(
    formula = z_behav ~ re(rpe_whitened, model = "vs"),
    models = list(vs = m_vs), data = subj_df, id = "id", regularization = "step_up"
  )
  cf_beh <- rep_beh$coefficients[grep("vs", rep_beh$coefficients$term), , drop = FALSE]
  expect_equal(nrow(cf_beh), 1L)
  expect_true(all(c("term", "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high") %in% names(cf_beh)))
  expect_equal(cf_beh$estimate[1L], 1.40, tolerance = 0.25)

  # Compare against attenuated naive OLS subject slopes
  ols_slopes <- vapply(split(df_l1, df_l1$id)[s_ids], function(d) stats::coef(stats::lm(y_vs ~ rpe_whitened, data = d, weights = pw))[[2L]], numeric(1L))
  b_naive <- stats::coef(stats::lm(z_behav ~ ols_slopes))[[2L]]
  expect_gt(cf_beh$estimate[1L], b_naive + 0.15)

  # Brain-to-brain random-slope structural regression (re on both LHS and RHS)
  rep_b2b <- medusa_between_subject_repair(
    formula = re(rpe_whitened, model = "mpfc") ~ re(rpe_whitened, model = "vs"),
    models = list(vs = m_vs, mpfc = m_mp), data = subj_df, id = "id", regularization = "step_up"
  )
  cf_b2b <- rep_b2b$coefficients[grep("vs", rep_b2b$coefficients$term), , drop = FALSE]
  expect_equal(cf_b2b$estimate[1L], 1.25, tolerance = 0.25)
})
