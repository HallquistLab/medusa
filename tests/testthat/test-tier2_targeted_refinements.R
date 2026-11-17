test_that("Tier 2B-1: decon_connectivity_repair match_resolution eliminates PPI product attenuation and K_S != K_T resolution bias", {
  set.seed(20261117L)
  tr <- 1.0
  n_time <- 215L
  n_trials <- 40L
  lib <- hrf_library(tr, peak_delays = c(4.5, 6.0, 7.5), dispersions = 1.0, undershoot_ratios = 6)
  pars <- attr(lib, "parameters")
  idx_fast <- which(abs(pars$peak_delay - 4.5) < 1e-6)[1L]
  idx_can <- which(abs(pars$peak_delay - 6.0) < 1e-6)[1L]
  hrf_fast <- lib[, idx_fast, drop = FALSE]
  hrf_can <- lib[, idx_can, drop = FALSE]

  reg_names <- c("R_seed", "R_targ", "R_null_s", "R_null_t", paste0("R_bg", 1:5))
  groups <- stats::setNames(rep(reg_names, each = 3L), paste0("v", seq_len(27L)))

  # 1A & 1C: Rapid/moderate ITI differential-HRF PPI (Seed 4.5s, Target 6.0s, Seed main x drive = 0.50)
  #          plus pervasive shared-noise null regions (R_null_s, R_null_t across 9 parcels)
  n_subj <- 18L
  true_bS <- 0.45
  true_bP <- 0.35
  fits_list <- vector("list", n_subj)
  names(fits_list) <- sprintf("sub%02d", seq_len(n_subj))

  for (s in seq_len(n_subj)) {
    bold_s <- ev_s <- cov_s <- list()
    for (r in 1:2) {
      itis <- runif(n_trials - 1L, 1.9, 3.5)
      t_ev <- cumsum(c(8.0, itis))
      xv <- as.numeric(scale(stats::rnorm(n_trials)))
      eta_s <- stats::rnorm(n_trials, 0, 0.70)
      eta_t <- (true_bS + true_bP * xv) * eta_s + stats::rnorm(n_trials, 0, 0.25)

      amp_seed <- 0.85 + 0.50 * xv + eta_s
      amp_targ <- 0.80 + 0.40 * xv + eta_t
      amp_ns <- 0.70 + stats::rnorm(n_trials, 0, 0.55)
      amp_nt <- 0.70 + stats::rnorm(n_trials, 0, 0.55)

      des_f <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials), seq_len(n_trials), n_trials, hrf_fast)[[1L]]
      des_c <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials), seq_len(n_trials), n_trials, hrf_can)[[1L]]

      sh_noise <- as.numeric(stats::arima.sim(list(ar = 0.35), n = n_time, sd = 0.28))
      bg_blks <- lapply(1:5, function(k) {
        sig_bg <- as.numeric(des_c %*% rep(0.55, n_trials))
        replicate(3L, 100 + sig_bg + sh_noise + stats::rnorm(n_time, sd = 0.15))
      })
      ym <- do.call(cbind, c(
        list(
          replicate(3L, 100 + as.numeric(des_f %*% amp_seed) + stats::rnorm(n_time, sd = 0.14)),
          replicate(3L, 100 + as.numeric(des_c %*% amp_targ) + stats::rnorm(n_time, sd = 0.14)),
          replicate(3L, 100 + as.numeric(des_f %*% amp_ns) + sh_noise + stats::rnorm(n_time, sd = 0.15)),
          replicate(3L, 100 + as.numeric(des_c %*% amp_nt) + sh_noise + stats::rnorm(n_time, sd = 0.15))
        ),
        bg_blks
      ))
      colnames(ym) <- names(groups)
      rn <- paste0("run", r)
      bold_s[[rn]] <- ym
      ev_s[[rn]] <- data.frame(onset = t_ev, duration = 0.0, condition = "fb")
      cov_s[[rn]] <- data.frame(x = xv)
    }
    hf <- estimate_crossrun_hrf(
      bold_s, ev_s, tr, lib, spatial_groups = groups, pooling = "parcel",
      min_parcel_voxels = 3L, noise_model = "ar1", noise_pcs = 0L, evidence = "calibrated"
    )
    dc <- deconvolve_crossrun_hrf(bold_s, hf, lambda = 0.8, scaling = "noise_sd", penalty = "sobolev", noise_pcs = 0L)
    rg <- decon_epoch_regression(
      dc, cov_s, events = ev_s, target_condition = "fb", lags = 0,
      demix_ridge = 0.05, restore_scale = TRUE, parcel_groups = groups,
      bold = bold_s, whiten_panel = FALSE, return_panel = TRUE
    )
    data.table::set(rg$panel, j = "id", value = sprintf("sub%02d", s))
    fits_list[[s]] <- rg
  }

  panel_stacked <- data.table::rbindlist(lapply(fits_list, `[[`, "panel"))
  expect_true("trial_res_vec" %in% names(panel_stacked))

  conn_unmatched <- decon_connectivity_repair(
    fits_list, "R_seed", "R_targ", covariates = "x", test_ppi = TRUE,
    method = "mmoms", match_resolution = FALSE, region_col = "parcel"
  )
  conn_matched_op <- decon_connectivity_repair(
    fits_list, "R_seed", "R_targ", covariates = "x", test_ppi = TRUE,
    method = "mmoms", match_resolution = TRUE, return_augmented_panel = TRUE, region_col = "parcel"
  )
  conn_matched_st <- decon_connectivity_repair(
    panel_stacked, "R_seed", "R_targ", covariates = "x", test_ppi = TRUE,
    method = "mmoms", match_resolution = TRUE, region_col = "parcel"
  )
  csv_tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(csv_tmp), add = TRUE)
  medusa:::medusa_write_table(panel_stacked, csv_tmp)
  conn_matched_csv <- decon_connectivity_repair(
    csv_tmp, "R_seed", "R_targ", covariates = "x", test_ppi = TRUE,
    method = "mmoms", match_resolution = TRUE, region_col = "parcel"
  )

  ppi_unm <- conn_unmatched$group[term == "seed_c:x_c"]$estimate
  ppi_op <- conn_matched_op$group[term == "seed_c:x_c"]$estimate
  ppi_st <- conn_matched_st$group[term == "seed_c:x_c"]$estimate
  ppi_csv <- conn_matched_csv$group[term == "seed_c:x_c"]$estimate

  # Resolution matching substantially reduces PPI product attenuation and agrees across operator, stacked, and CSV inputs
  expect_gt(ppi_op, ppi_unm + 0.04)
  expect_lt(abs(ppi_op - true_bP), 0.70 * abs(ppi_unm - true_bP))
  expect_equal(ppi_op, ppi_st, tolerance = 0.02)
  expect_equal(ppi_st, ppi_csv, tolerance = 1e-5)

  # Moment-preserving REPAIR augmented columns still reproduce subject M-MOMS slopes under match_resolution = TRUE
  p_aug1 <- conn_matched_op$panel_repair[id == "sub01"]
  lm_aug1 <- stats::lm(targ_c ~ seed_c_repair + x_c_repair + `seed_c:x_c_repair`, data = p_aug1, weights = repair_weight)
  expect_equal(unname(stats::coef(lm_aug1)[-1L]), conn_matched_op$subjects[id == "sub01"]$estimate, tolerance = 1e-5)

  # Shared-error null remains calibrated under match_resolution = TRUE
  conn_sh_null <- decon_connectivity_repair(
    fits_list, "R_null_s", "R_null_t", test_ppi = FALSE,
    method = "mmoms", shared_error = "auto", match_resolution = TRUE, region_col = "parcel"
  )
  conn_sh_zero <- decon_connectivity_repair(
    fits_list, "R_null_s", "R_null_t", test_ppi = FALSE,
    method = "mmoms", shared_error = "zero", match_resolution = TRUE, region_col = "parcel"
  )
  expect_lt(abs(conn_sh_null$group[term == "seed_c"]$estimate), 0.50 * abs(conn_sh_zero$group[term == "seed_c"]$estimate))
  expect_gt(conn_sh_null$group[term == "seed_c"]$p.value, 0.05)

  # 1B: Asymmetric raw resolution ||K_S||_F != ||K_T||_F (heavily smoothed Seed vs sharp Target, demix_target_trials = FALSE)
  set.seed(20261121L)
  fits_asym <- vector("list", 16L)
  groups_2 <- c(v1 = "Seed_smooth", v2 = "Targ_sharp")
  for (s in 1:16) {
    bold_2 <- ev_2 <- list()
    for (r in 1:2) {
      t_ev <- cumsum(c(8.0, runif(29L, 2.2, 4.2)))
      eta_s <- stats::rnorm(30L, 0, 0.75)
      eta_t <- 0.50 * eta_s + stats::rnorm(30L, 0, 0.30)
      des_c <- hrf_event_design_cpp(160L, tr, 0, t_ev, rep(0, 30L), 1:30, 30L, hrf_can)[[1L]]
      ym2 <- cbind(
        v1 = 100 + as.numeric(des_c %*% eta_s) + stats::rnorm(160L, sd = 0.14),
        v2 = 100 + as.numeric(des_c %*% eta_t) + stats::rnorm(160L, sd = 0.14)
      )
      bold_2[[paste0("run", r)]] <- ym2
      ev_2[[paste0("run", r)]] <- data.frame(onset = t_ev, duration = 0, condition = "fb")
    }
    hf2 <- estimate_crossrun_hrf(bold_2, ev_2, tr, lib[, idx_can, drop = FALSE], noise_model = "ols", evidence = "legacy", noise_pcs = 0L)
    dc2 <- deconvolve_crossrun_hrf(
      bold_2, hf2, lambda = c(5.0, 0.20),
      scaling = "noise_sd", penalty = "sobolev", noise_pcs = 0L
    )
    rg2 <- decon_epoch_regression(
      dc2, covariates = NULL, events = ev_2, target_condition = "fb", lags = 0,
      demix_target_trials = FALSE, restore_scale = TRUE, parcel_groups = groups_2, return_panel = TRUE
    )
    data.table::set(rg2$panel, j = "id", value = sprintf("sub%02d", s))
    fits_asym[[s]] <- rg2
  }
  panel_asym <- data.table::rbindlist(lapply(fits_asym, `[[`, "panel"))
  asym_unm <- decon_connectivity_repair(panel_asym, "Seed_smooth", "Targ_sharp", method = "mmoms", shared_error = "zero", match_resolution = FALSE, region_col = "parcel")$group$estimate
  asym_mat_st <- decon_connectivity_repair(panel_asym, "Seed_smooth", "Targ_sharp", method = "mmoms", shared_error = "zero", match_resolution = TRUE, region_col = "parcel")$group$estimate
  asym_mat_op <- decon_connectivity_repair(fits_asym, "Seed_smooth", "Targ_sharp", method = "mmoms", shared_error = "zero", match_resolution = TRUE, region_col = "parcel")$group$estimate
  # Unmatched stacked panel overshoots substantially (> 0.65 vs true 0.50); match_resolution eliminates overshoot across stacked and operator inputs
  expect_gt(asym_unm, 0.65)
  expect_lt(abs(asym_mat_st - 0.50), 0.35 * abs(asym_unm - 0.50))
  expect_lt(abs(asym_mat_op - 0.50), 0.065)
})

test_that("Tier 2B-2: demix_ridge = 'auto' sharpens high-SNR trial/capacity recovery and rapid-ITI slopes while remaining covariate-agnostic", {
  set.seed(20261118L)
  tr <- 1.0
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_c <- lib[, 1L, drop = FALSE]

  # 2A: Study 4B-style individual working-memory capacity K_s recovery + rapid-ITI alternating-trial slope
  n_subj <- 16L
  true_Ks <- seq(2.2, 5.2, length.out = n_subj)
  est_K_auto <- est_K_fix <- numeric(n_subj)
  cor_trial_auto <- cor_trial_fix <- numeric(n_subj)
  alt_beta_auto <- alt_beta_fix <- numeric(n_subj)

  for (s in seq_len(n_subj)) {
    bold_s <- ev_s <- cov_s <- list()
    th_true_all <- list()
    for (r in 1:2) {
      n_tr <- 30L
      # Rapid event spacing (2.0-3.2s) where fixed demix_ridge = 0.05 over-smooths adjacent trials
      t_ev <- cumsum(c(6.0, runif(n_tr - 1L, 2.0, 3.2)))
      n_t <- as.integer(ceiling(max(t_ev) + 24.0))
      loads <- sample(rep(1:6, length.out = n_tr))
      x_alt <- as.numeric(scale(stats::arima.sim(list(ar = -0.75), n = n_tr)))
      cap_sig <- pmin(loads, true_Ks[s])
      th_cap <- 0.55 * cap_sig + stats::rnorm(n_tr, 0, 0.12)
      th_alt <- 0.80 + 0.65 * x_alt + stats::rnorm(n_tr, 0, 0.15)
      th_true_all[[r]] <- list(cap = th_cap - mean(th_cap), loads = loads)

      des_r <- hrf_event_design_cpp(n_t, tr, 0, t_ev, rep(0, n_tr), seq_len(n_tr), n_tr, hrf_c)[[1L]]
      ym <- cbind(
        v_cap1 = 100 + as.numeric(des_r %*% th_cap) + stats::rnorm(n_t, sd = 0.025),
        v_cap2 = 100 + as.numeric(des_r %*% th_cap) + stats::rnorm(n_t, sd = 0.025),
        v_alt1 = 100 + as.numeric(des_r %*% th_alt) + stats::rnorm(n_t, sd = 0.04),
        v_null = 100 + stats::rnorm(n_t, sd = 0.25)
      )
      rn <- paste0("run", r)
      bold_s[[rn]] <- ym
      ev_s[[rn]] <- data.frame(onset = t_ev, duration = 0, condition = "delay", load = loads)
      cov_s[[rn]] <- data.frame(x_alt = x_alt)
    }
    hf <- estimate_crossrun_hrf(bold_s, ev_s, tr, lib, noise_model = "ols", evidence = "legacy", noise_pcs = 0L)
    dc <- deconvolve_crossrun_hrf(bold_s, hf, lambda = 0.5, scaling = "noise_sd", penalty = "sobolev", noise_pcs = 0L)

    rg_fix <- decon_epoch_regression(
      dc, cov_s, events = ev_s, target_condition = "delay", lags = 0,
      demix_ridge = 0.05, restore_scale = TRUE, return_panel = TRUE
    )
    rg_auto <- decon_epoch_regression(
      dc, cov_s, events = ev_s, target_condition = "delay", lags = 0,
      demix_ridge = "auto", restore_scale = TRUE, return_panel = TRUE
    )
    if (s == 1L) {
      # Verify covariate agnosticism: demix_ridge = 'auto' never peeks at covariates
      rg_nocov <- decon_epoch_regression(
        dc, covariates = NULL, events = ev_s, target_condition = "delay", lags = 0,
        demix_ridge = "auto", restore_scale = TRUE, return_panel = TRUE
      )
      expect_equal(rg_auto$panel$decon_demixed, rg_nocov$panel$decon_demixed, tolerance = 1e-12)
      t_null <- rg_auto$coefficients$t_stat[rg_auto$coefficients$voxel == "v_null"]
      expect_true(is.finite(t_null) && abs(t_null) < 3.0)
      # Pure-noise unit receives MORE conservative ridge (> 0.05), reducing trial error SE
      expect_lt(
        mean(rg_auto$panel[voxel == "v_null"]$se_demixed),
        mean(rg_fix$panel[voxel == "v_null"]$se_demixed)
      )
    }

    alt_beta_fix[s] <- stats::coef(stats::lm(
      decon_demixed ~ factor(run) + x_alt,
      data = rg_fix$panel[voxel == "v_alt1"]
    ))[["x_alt"]]
    alt_beta_auto[s] <- stats::coef(stats::lm(
      decon_demixed ~ factor(run) + x_alt,
      data = rg_auto$panel[voxel == "v_alt1"]
    ))[["x_alt"]]

    th_vec <- c(th_true_all[[1L]]$cap, th_true_all[[2L]]$cap)
    ld_vec <- c(th_true_all[[1L]]$loads, th_true_all[[2L]]$loads)
    z_fix <- rg_fix$panel[voxel == "v_cap1"]$decon_demixed
    z_fix_c <- ave(z_fix, rg_fix$panel[voxel == "v_cap1"]$run, FUN = function(u) u - mean(u))
    z_auto <- rg_auto$panel[voxel == "v_cap1"]$decon_demixed
    z_auto_c <- ave(z_auto, rg_auto$panel[voxel == "v_cap1"]$run, FUN = function(u) u - mean(u))

    cor_trial_fix[s] <- stats::cor(z_fix_c, th_vec)
    cor_trial_auto[s] <- stats::cor(z_auto_c, th_vec)

    fit_k <- function(yv) {
      k_grid <- seq(1.5, 5.8, by = 0.05)
      sse <- vapply(k_grid, function(kc) {
        xc <- pmin(ld_vec, kc) - mean(pmin(ld_vec, kc))
        b <- sum(xc * yv) / sum(xc * xc)
        sum((yv - b * xc)^2)
      }, numeric(1L))
      k_grid[which.min(sse)]
    }
    est_K_fix[s] <- fit_k(z_fix_c)
    est_K_auto[s] <- fit_k(z_auto_c)
  }

  # Single-trial correlation, individual K_s MSE, and rapid alternating-trial slope all improve over fixed 0.05
  expect_gt(mean(cor_trial_auto), mean(cor_trial_fix) + 0.008)
  expect_lt(mean((est_K_auto - true_Ks)^2), mean((est_K_fix - true_Ks)^2))
  expect_lt(abs(mean(alt_beta_auto) - 0.65), 0.80 * abs(mean(alt_beta_fix) - 0.65))
})

test_that("Tier 2B-3: single-run (R = 1) 16s guard-buffered LOBO CV controls false HRF departures while recovering fast/slow HRFs", {
  set.seed(20261119L)
  tr <- 1.0
  n_time <- 360L
  n_trials <- 60L
  lib <- hrf_library(tr, peak_delays = c(4.5, 6.0, 7.5), dispersions = 1.0, undershoot_ratios = 6)
  pars <- attr(lib, "parameters")
  idx_fast <- which(abs(pars$peak_delay - 4.5) < 1e-6)[1L]
  idx_can <- which(abs(pars$peak_delay - 6.0) < 1e-6)[1L]
  idx_slow <- which(abs(pars$peak_delay - 7.5) < 1e-6)[1L]

  groups <- stats::setNames(rep(c("R_fast", "R_can", "R_slow", "R_null"), each = 4L), paste0("v", 1:16))
  t_sec <- (0:(n_time - 1L)) * tr

  n_reps <- 30L
  can_false_depart <- integer(n_reps)
  null_false_task <- integer(n_reps)
  exact_order <- integer(n_reps)

  for (m in seq_len(n_reps)) {
    t_ev <- cumsum(c(8.0, runif(n_trials - 1L, 3.8, 6.8)))
    t_ev <- 8.0 + (t_ev - 8.0) * min(1.0, (n_time - 28.0) / max(max(t_ev) - 8.0, 1.0))
    amps <- 1.40 + stats::rnorm(n_trials, 0, 0.45)

    des_f <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials), seq_len(n_trials), n_trials, lib[, idx_fast, drop = FALSE])[[1L]]
    des_c <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials), seq_len(n_trials), n_trials, lib[, idx_can, drop = FALSE])[[1L]]
    des_s <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials), seq_len(n_trials), n_trials, lib[, idx_slow, drop = FALSE])[[1L]]

    mk_ar_drift <- function(sd_e) {
      drift <- 0.25 * sin(2 * pi * t_sec / 180.0 + runif(1L, 0, 2 * pi))
      drift + as.numeric(stats::arima.sim(list(ar = 0.45), n = n_time, sd = sd_e))
    }

    ym <- cbind(
      vapply(1:4, function(i) 100 + as.numeric(des_f %*% amps) + mk_ar_drift(0.16), numeric(n_time)),
      vapply(1:4, function(i) 100 + as.numeric(des_c %*% amps) + mk_ar_drift(0.16), numeric(n_time)),
      vapply(1:4, function(i) 100 + as.numeric(des_s %*% amps) + mk_ar_drift(0.16), numeric(n_time)),
      vapply(1:4, function(i) 100 + mk_ar_drift(0.20), numeric(n_time))
    )
    colnames(ym) <- names(groups)
    ev1 <- list(run1 = data.frame(onset = t_ev, duration = 0.0, condition = "fb"))

    fit1 <- estimate_crossrun_hrf(
      list(run1 = ym), ev1, tr, lib,
      spatial_groups = groups, pooling = "parcel", min_parcel_voxels = 3L,
      fallback_kernel = idx_can, noise_model = "ar1", noise_pcs = 0L,
      evidence = "calibrated", selection_mode = "soft_average"
    )
    sel <- fit1$selections$run1
    can_false_depart[m] <- as.integer(sel$kernel_index[5L] != idx_can)
    null_false_task[m] <- as.integer(sel$status[13L] != "fallback_prespecified")

    kw <- fit1$kernel_weights$run1
    pk_fast <- sum(pars$peak_delay * kw[, 1L])
    pk_can <- sum(pars$peak_delay * kw[, 5L])
    pk_slow <- sum(pars$peak_delay * kw[, 9L])
    exact_order[m] <- as.integer(pk_fast < pk_can - 0.5 && pk_can < pk_slow - 0.5)
  }

  # Canonical departure FPR and null-parcel FPR stay within nominal binomial bounds while Fast < Can < Slow order is recovered reliably
  expect_lte(mean(can_false_depart), 0.10)
  expect_lte(mean(null_false_task), 0.14)
  expect_gte(mean(exact_order), 0.90)
})

test_that("Tier 2B-4: mixed_by conditional tau2_bio reweighting eliminates marginal-panel slope attenuation", {
  set.seed(20261120L)
  tr <- 1.0
  n_time <- 160L
  n_trials <- 26L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_c <- lib[, 1L, drop = FALSE]

  n_subj <- 24L
  # Heterogeneous positive slopes: high-slope participants have much larger marginal Var(theta) = beta_s^2 Var(x) + tau_resid^2,
  # which systematically downweights them when Step 3 ran with covariates = NULL unless mixed_by reweights conditionally.
  subj_betas <- seq(0.10, 1.30, length.out = n_subj)
  true_pop_beta <- mean(subj_betas) # 0.70

  panels_marg <- panels_cond <- vector("list", n_subj)

  for (s in seq_len(n_subj)) {
    bold_s <- ev_s <- cov_s <- list()
    sid <- sprintf("sub%02d", s)
    beh_runs <- list()
    for (r in 1:2) {
      t_ev <- cumsum(c(8.0, runif(n_trials - 1L, 3.8, 6.2)))
      xv <- as.numeric(scale(stats::rnorm(n_trials)))
      theta <- 0.50 + subj_betas[s] * xv + stats::rnorm(n_trials, 0, 0.12)
      des_r <- hrf_event_design_cpp(n_time, tr, 0, t_ev, rep(0, n_trials), seq_len(n_trials), n_trials, hrf_c)[[1L]]
      ym <- cbind(v1 = 100 + as.numeric(des_r %*% theta) + stats::rnorm(n_time, sd = 0.16))
      rn <- paste0("run", r)
      bold_s[[rn]] <- ym
      ev_s[[rn]] <- data.frame(onset = t_ev, duration = 0.0, condition = "fb")
      cov_s[[rn]] <- data.frame(x = xv)
      beh_runs[[r]] <- data.table::data.table(id = sid, run = rn, trial = seq_len(n_trials), x = xv)
    }
    hf <- estimate_crossrun_hrf(bold_s, ev_s, tr, lib, noise_model = "ols", evidence = "legacy", noise_pcs = 0L)
    dc <- deconvolve_crossrun_hrf(bold_s, hf, lambda = 0.5, scaling = "noise_sd", penalty = "sobolev", noise_pcs = 0L)

    rg_m <- decon_epoch_regression(
      dc, covariates = NULL, events = ev_s, target_condition = "fb",
      lags = 0, restore_scale = TRUE, return_panel = TRUE
    )
    rg_c <- decon_epoch_regression(
      dc, covariates = cov_s, events = ev_s, target_condition = "fb",
      lags = 0, tau2_mode = "residual", restore_scale = TRUE, return_panel = TRUE
    )
    pm <- data.table::copy(rg_m$panel)
    pm[, id := sid]
    pm <- merge(pm, data.table::rbindlist(beh_runs), by = c("id", "run", "trial"))
    pc <- data.table::copy(rg_c$panel)
    pc[, id := sid]
    panels_marg[[s]] <- pm
    panels_cond[[s]] <- pc
  }

  dt_marg <- data.table::rbindlist(panels_marg)
  dt_cond <- data.table::rbindlist(panels_cond)

  # Confirm marginal panel severely downweights the highest-slope participants relative to lowest-slope participants
  w_low_marg <- mean(dt_marg[id == "sub01"]$precision_weight_cohort)
  w_high_marg <- mean(dt_marg[id == sprintf("sub%02d", n_subj)]$precision_weight_cohort)
  expect_gt(w_low_marg, 3.0 * w_high_marg)

  mb_marg_unfixed <- suppressMessages(suppressWarnings(mixed_by(
    dt_marg,
    outcomes = "decon_demixed",
    rhs_model_formulae = list(m1 = ~ x + (1 | id)),
    split_on = c("voxel", "evt_time"),
    weights = "precision_weight_cohort",
    reweight_conditional_tau2 = FALSE,
    calculate = "parameter_estimates_reml"
  )))
  mb_marg_auto <- suppressMessages(suppressWarnings(mixed_by(
    dt_marg,
    outcomes = "decon_demixed",
    rhs_model_formulae = list(m1 = ~ x + (1 | id)),
    split_on = c("voxel", "evt_time"),
    weights = "precision_weight_cohort",
    reweight_conditional_tau2 = "auto",
    calculate = "parameter_estimates_reml"
  )))
  mb_cond_ref <- suppressMessages(suppressWarnings(mixed_by(
    dt_cond,
    outcomes = "decon_demixed",
    rhs_model_formulae = list(m1 = ~ x + (1 | id)),
    split_on = c("voxel", "evt_time"),
    weights = "precision_weight_cohort",
    calculate = "parameter_estimates_reml"
  )))

  b_unfixed <- mb_marg_unfixed$coef_df_reml[term == "x"]$estimate
  b_auto <- mb_marg_auto$coef_df_reml[term == "x"]$estimate
  b_ref <- mb_cond_ref$coef_df_reml[term == "x"]$estimate

  # Unfixed marginal weights attenuate the group slope severely (< 0.55 vs true 0.70);
  # conditional reweighting eliminates the attenuation and matches Step-3 conditional weights!
  expect_lt(b_unfixed, true_pop_beta - 0.12)
  expect_lt(abs(b_auto - true_pop_beta), 0.06)
  expect_equal(b_auto, b_ref, tolerance = 0.025)

  # Also verify random-slope formula (1 + x | id) works cleanly and preserves caller data ownership
  dt_marg_snapshot <- data.table::copy(dt_marg)
  mb_rs <- suppressMessages(suppressWarnings(mixed_by(
    dt_marg,
    outcomes = "decon_demixed",
    rhs_model_formulae = list(m1 = ~ x + (1 + x | id)),
    split_on = c("voxel", "evt_time"),
    weights = "precision_weight_cohort",
    calculate = "parameter_estimates_reml"
  )))
  expect_equal(dt_marg$precision_weight_cohort, dt_marg_snapshot$precision_weight_cohort, tolerance = 1e-15)
  expect_lt(abs(mb_rs$coef_df_reml[term == "x"]$estimate - true_pop_beta), 0.06)
})
