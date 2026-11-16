test_that("v3 (#1, #2, #6): bold-free handoff, dynamic shared-error covreg, and Operator vs Stacked parity", {
  set.seed(20261120L)
  tr <- 1.0
  n_time <- 220L
  n_trials <- 26L
  lib <- hrf_library(tr, peak_delays = c(5.0, 6.0), dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 2L, drop = FALSE]
  groups <- stats::setNames(rep(c("R_seed", "R_targ", "R_null"), each = 3L), paste0("v", 1:9))

  n_subj <- 16L
  fits_null <- vector("list", n_subj)
  panels_null <- vector("list", n_subj)
  panels_true <- vector("list", n_subj)

  for (s in seq_len(n_subj)) {
    bold_n <- bold_t <- ev_s <- cov_s <- list()
    for (r in 1:2) {
      # Jittered ITIs (3.2s to 10.5s) so task envelope has genuine quiet baselines vs bursts
      itis <- pmin(10.5, 3.2 + stats::rexp(n_trials - 1L, rate = 1 / 3.2))
      t_ev <- cumsum(c(6.0, itis))
      t_ev <- 6.0 + (t_ev - 6.0) * min(1.0, (n_time - 24.0) / max(max(t_ev) - 6.0, 1.0))
      xv <- as.numeric(scale(stats::rnorm(n_trials)))
      des <- hrf_event_design_cpp(
        n_time, tr, 0.0, t_ev, rep(0.0, n_trials), seq_len(n_trials), n_trials, hrf_k
      )[[1L]]
      des_pk <- des / max(hrf_k)

      # Dynamic covariate-modulated shared physiological burst (zero true neural PPI in bold_n!)
      e_t0 <- as.numeric(scale(des_pk %*% rep(1.0, n_trials)))
      e_x0 <- as.numeric(scale(des_pk %*% xv))
      v_sh <- pmax(0.06, 0.26 + 0.045 * e_t0 + 0.060 * e_x0)
      e_sh <- as.numeric(stats::filter(sqrt(v_sh) * stats::rnorm(n_time), filter = 0.35, method = "recursive"))

      eta_s <- stats::rnorm(n_trials, 0, 0.85)
      eta_t_null <- stats::rnorm(n_trials, 0, 0.85)
      eta_t_true <- (0.45 + 0.40 * xv) * eta_s + stats::rnorm(n_trials, 0, 0.28)
      eta_0 <- stats::rnorm(n_trials, 0, 0.65)

      mk_mat <- function(eta_targ_v) {
        ym <- cbind(
          replicate(3L, 100 + as.numeric(des_pk %*% (0.8 + 0.30 * xv + eta_s)) + e_sh + stats::rnorm(n_time, sd = 0.25)),
          replicate(3L, 100 + as.numeric(des_pk %*% (0.8 + 0.30 * xv + eta_targ_v)) + e_sh + stats::rnorm(n_time, sd = 0.25)),
          replicate(3L, 100 + as.numeric(des_pk %*% (0.6 + eta_0)) + stats::rnorm(n_time, sd = 0.25))
        )
        colnames(ym) <- names(groups)
        ym
      }
      rn <- paste0("run", r)
      bold_n[[rn]] <- mk_mat(eta_t_null)
      bold_t[[rn]] <- mk_mat(eta_t_true)
      ev_s[[rn]] <- data.frame(onset = t_ev, duration = 0.0, condition = "fb")
      cov_s[[rn]] <- data.frame(x = xv)
    }

    fh_n <- estimate_crossrun_hrf(
      bold_n, ev_s, tr, lib, spatial_groups = groups,
      pooling = "parcel", amplitude_model = "trial", noise_model = "ar1"
    )
    dc_n <- deconvolve_crossrun_hrf(bold_n, fh_n, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev")
    expect_true(is.list(dc_n$bold_obs))
    expect_equal(ncol(dc_n$bold_obs[[1L]]), 9L)

    # Step 3 runs full v3 engine WITHOUT passing `bold`
    rg_n <- decon_epoch_regression(
      dc_n, cov_s, events = ev_s, target_condition = "fb",
      lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE
    )
    expect_equal(rg_n$shared_error_model$version, 3L)
    expect_equal(rg_n$shared_error_model$mode, "exact")
    expect_true(all(c("shared_noise_vec", "trial_cor_vec") %in% names(rg_n$panel)))

    sid <- sprintf("sub%02d", s)
    rg_n$panel$id <- sid
    fits_null[[s]] <- rg_n
    panels_null[[s]] <- data.table::copy(rg_n$panel)

    fh_t <- estimate_crossrun_hrf(
      bold_t, ev_s, tr, lib, spatial_groups = groups,
      pooling = "parcel", amplitude_model = "trial", noise_model = "ar1"
    )
    dc_t <- deconvolve_crossrun_hrf(bold_t, fh_t, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev")
    rg_t <- decon_epoch_regression(
      dc_t, cov_s, events = ev_s, target_condition = "fb",
      lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE
    )
    rg_t$panel$id <- sid
    panels_true[[s]] <- rg_t$panel
  }
  names(fits_null) <- sprintf("sub%02d", seq_len(n_subj))
  dt_null_stacked <- data.table::rbindlist(panels_null)
  dt_true_stacked <- data.table::rbindlist(panels_true)
  expect_null(attr(dt_null_stacked, "shared_error_model"))

  # Operator mode vs bare post-rbindlist Stacked mode agree closely and suppress dynamic shared-noise bias
  res_op <- decon_connectivity_repair(fits_null, "R_seed", "R_targ", test_ppi = TRUE, method = "mmoms", shared_error = "auto")
  res_st <- decon_connectivity_repair(dt_null_stacked, "R_seed", "R_targ", test_ppi = TRUE, method = "mmoms", shared_error = "auto")
  res_zero <- decon_connectivity_repair(dt_null_stacked, "R_seed", "R_targ", test_ppi = TRUE, method = "mmoms", shared_error = "zero")

  expect_lt(max(abs(res_op$subjects$estimate - res_st$subjects$estimate)), 0.04)
  expect_gt(res_zero$group[term == "seed_c"]$estimate, 0.18)
  expect_lt(abs(res_st$group[term == "seed_c"]$estimate), 0.14)
  expect_lt(abs(res_st$group[term == "seed_c:x_c"]$estimate), 0.16)
  expect_gt(res_st$group[term == "seed_c:x_c"]$p.value, 0.05)

  # True neural connectivity & PPI recovery
  res_true <- decon_connectivity_repair(dt_true_stacked, "R_seed", "R_targ", test_ppi = TRUE, method = "mmoms", shared_error = "auto")
  expect_lt(abs(res_true$group[term == "seed_c"]$estimate - 0.45), 0.15)
  expect_lt(abs(res_true$group[term == "seed_c:x_c"]$estimate - 0.40), 0.16)
  expect_lt(res_true$group[term == "seed_c:x_c"]$p.value, 0.01)
})

test_that("v3 (#3 / 3A): variance_confounds Schur-complement invariance to null columns and Step-2 rWLS burst handling", {
  set.seed(20261121L)
  tr <- 1.0
  n_time <- 180L
  n_trials <- 30L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]
  groups <- stats::setNames(rep(c("R_seed", "R_targ"), each = 3L), paste0("v", 1:6))

  bold_s <- ev_s <- cov_s <- vc_1 <- vc_6 <- list()
  for (r in 1:2) {
    t_ev <- cumsum(c(6.0, stats::runif(n_trials - 1L, 4.0, 6.2)))
    xv <- as.numeric(scale(stats::rnorm(n_trials)))
    des <- hrf_event_design_cpp(
      n_time, tr, 0.0, t_ev, rep(0.0, n_trials), seq_len(n_trials), n_trials, hrf_k
    )[[1L]]
    # Physical non-negative RVT envelope correlated with positive x trials + 2 motion spikes
    rvt_raw <- pmax(0.0, 0.5 + 0.85 * as.numeric(scale(des %*% pmax(0.0, xv))) + 0.35 * stats::rnorm(n_time))
    e_rvt <- (0.15 + 0.35 * rvt_raw) * stats::rnorm(n_time)
    spike_frames <- c(45L, 125L)
    e_rvt[spike_frames] <- e_rvt[spike_frames] + c(4.2, -3.8)

    eta_s <- stats::rnorm(n_trials, 0, 0.65)
    eta_t <- (0.40 + 0.35 * xv) * eta_s + stats::rnorm(n_trials, 0, 0.25)
    ym <- cbind(
      replicate(3L, 100 + as.numeric(des %*% (0.7 + eta_s)) + e_rvt + stats::rnorm(n_time, sd = 0.20)),
      replicate(3L, 100 + as.numeric(des %*% (0.7 + eta_t)) + e_rvt + stats::rnorm(n_time, sd = 0.20))
    )
    colnames(ym) <- names(groups)
    rn <- paste0("run", r)
    bold_s[[rn]] <- ym
    ev_s[[rn]] <- data.frame(onset = t_ev, duration = 0.0, condition = "fb")
    cov_s[[rn]] <- data.frame(x = xv)
    null_5 <- matrix(abs(stats::rnorm(n_time * 5L)), nrow = n_time, ncol = 5L)
    vc_1[[rn]] <- cbind(rvt = rvt_raw)
    vc_6[[rn]] <- cbind(rvt = rvt_raw, null_5)
  }

  fh <- estimate_crossrun_hrf(bold_s, ev_s, tr, lib, spatial_groups = groups, pooling = "parcel", noise_model = "ar1")
  dc <- deconvolve_crossrun_hrf(bold_s, fh, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev", rwls = "auto")

  rg_base <- decon_epoch_regression(
    dc, cov_s, events = ev_s, target_condition = "fb",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE
  )
  rg_vc1 <- decon_epoch_regression(
    dc, cov_s, events = ev_s, target_condition = "fb",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE,
    variance_confounds = vc_1
  )
  rg_vc6 <- decon_epoch_regression(
    dc, cov_s, events = ev_s, target_condition = "fb",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE,
    variance_confounds = vc_6
  )

  # Pure-null variance_confounds test on clean data: Schur-complement restoration keeps estimates invariant
  null_only_6 <- lapply(1:2, function(r) matrix(stats::rnorm(n_time * 6L), nrow = n_time, ncol = 6L))
  names(null_only_6) <- c("run1", "run2")
  rg_null6 <- decon_epoch_regression(
    dc, cov_s, events = ev_s, target_condition = "fb",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE,
    variance_confounds = null_only_6
  )

  c_base <- decon_connectivity_repair(rg_base, "R_seed", "R_targ", test_ppi = TRUE)
  c_vc1 <- decon_connectivity_repair(rg_vc1, "R_seed", "R_targ", test_ppi = TRUE)
  c_vc6 <- decon_connectivity_repair(rg_vc6, "R_seed", "R_targ", test_ppi = TRUE)
  c_n6 <- decon_connectivity_repair(rg_null6, "R_seed", "R_targ", test_ppi = TRUE)

  expect_true(all(is.finite(c_vc6$group$estimate)))
  expect_true(all(rg_vc6$shared_error_model$sig0_rob <= rg_vc6$shared_error_model$sig0_sd + 1e-12))
  # 6 pure-null variance confounds produce near-zero drift relative to no variance_confounds
  expect_lt(max(abs(c_n6$group$estimate - c_base$group$estimate)), 0.05)
  # Adding 5 noisy null columns alongside true RVT stays close to 1-column RVT
  expect_lt(max(abs(c_vc6$group$estimate - c_vc1$group$estimate)), 0.08)
})

test_that("v3 (#4 / 4A & #5 / 5A): duration/multi-lag M_s protection and single-subject N = 1 rapid-ITI cross-trial SEs", {
  set.seed(20261122L)
  tr <- 1.0
  n_time <- 210L
  n_trials <- 42L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]
  groups <- stats::setNames(rep(c("R_seed", "R_targ"), each = 3L), paste0("v", 1:6))

  # 1. Single-subject N = 1 rapid-ITI (1.8-3.4s) Monte Carlo calibration check (28 independent N = 1 subjects)
  n_rep <- 28L
  est_op <- se_ct_op <- se_noct_op <- cov_ct_op <- cov_ct_st <- df_ct_op <- numeric(n_rep)
  beta_true <- 0.35

  for (rep_i in seq_len(n_rep)) {
    bold_1 <- ev_1 <- cov_1 <- list()
    for (r in 1:2) {
      # Rapid event-related ITIs (1.9 to 3.3s) + variable durations (0.8 to 2.8s) + unmodeled AR(2)/drift noise
      t_ev <- cumsum(c(5.0, stats::runif(n_trials - 1L, 1.9, 3.3)))
      dur_v <- stats::runif(n_trials, 0.8, 2.8)
      xv <- as.numeric(scale(stats::rnorm(n_trials)))
      eta_s <- stats::rnorm(n_trials, 0, 0.65)
      eta_t <- beta_true * eta_s + stats::rnorm(n_trials, 0, 0.35)
      des <- hrf_event_design_cpp(
        n_time, tr, 0.0, t_ev, dur_v, seq_len(n_trials), n_trials, hrf_k
      )[[1L]]
      t_sec <- (seq_len(n_time) - 1.0) * tr
      mk_noise <- function() {
        as.numeric(stats::arima.sim(list(ar = c(0.55, 0.15)), n = n_time, sd = 0.25)) +
          0.14 * sin(2 * pi * 0.03 * t_sec + stats::runif(1, 0, 2 * pi))
      }
      e_sh <- mk_noise()
      ym <- cbind(
        replicate(3L, 100 + as.numeric(des %*% (0.6 + eta_s)) + 0.6 * e_sh + 0.8 * mk_noise()),
        replicate(3L, 100 + as.numeric(des %*% (0.6 + eta_t)) + 0.6 * e_sh + 0.8 * mk_noise())
      )
      colnames(ym) <- names(groups)
      rn <- paste0("run", r)
      bold_1[[rn]] <- ym
      ev_1[[rn]] <- data.frame(onset = t_ev, duration = dur_v, condition = "fb")
      cov_1[[rn]] <- data.frame(x = xv)
    }
    fh <- estimate_crossrun_hrf(bold_1, ev_1, tr, lib, spatial_groups = groups, pooling = "parcel", noise_model = "ar1")
    dc <- deconvolve_crossrun_hrf(bold_1, fh, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev")
    rg <- decon_epoch_regression(
      dc, cov_1, events = ev_1, target_condition = "fb",
      lags = c(0, 2), source_model = "duration", trial_durations = lapply(ev_1, `[[`, "duration"),
      restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE
    )
    expect_true("lag_noise_vec" %in% names(rg$panel))

    # N = 1 Operator with vs without crosstrial_se, and N = 1 bare Stacked data.table
    fit_ct <- decon_connectivity_repair(rg, "R_seed", "R_targ", crosstrial_se = TRUE)
    fit_no <- decon_connectivity_repair(rg, "R_seed", "R_targ", crosstrial_se = FALSE)
    dt_bare <- data.table::rbindlist(list(rg$panel))
    fit_st <- decon_connectivity_repair(dt_bare, "R_seed", "R_targ", crosstrial_se = TRUE)

    g_ct <- fit_ct$group[term == "seed_c"]
    g_no <- fit_no$group[term == "seed_c"]
    g_st <- fit_st$group[term == "seed_c"]

    est_op[rep_i] <- g_ct$estimate
    se_ct_op[rep_i] <- g_ct$std.error
    se_noct_op[rep_i] <- g_no$std.error
    df_ct_op[rep_i] <- g_ct$df
    cov_ct_op[rep_i] <- as.numeric(g_ct$conf.low <= beta_true && beta_true <= g_ct$conf.high)
    cov_ct_st[rep_i] <- as.numeric(g_st$conf.low <= beta_true && beta_true <= g_st$conf.high)
  }

  emp_sd <- stats::sd(est_op)
  expect_true(all(is.finite(df_ct_op) & df_ct_op < 2L * n_trials))
  # Cross-trial sandwich SE corrects downward bias of diagonal-only SE under rapid ITIs + slow noise
  expect_gt(mean(se_ct_op), mean(se_noct_op))
  expect_gt(mean(se_ct_op) / emp_sd, 0.80)
  expect_gte(mean(cov_ct_op), 0.82)
  expect_gte(mean(cov_ct_st), 0.80)
})

test_that("v3 (#7 & #8): P = 120 multi-lag Block-PAF, V = 600 voxel-safe lowrank_voxel mode, and aggregate_decon_panel()", {
  set.seed(20261123L)
  tr <- 1.0
  n_time <- 140L
  n_trials <- 24L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]

  # 1. High-P parcel scale (P = 120 parcels x 2 lags = 240 > old dim_al = 96 cutoff!)
  P_hi <- 120L
  p_ids <- sprintf("P%03d", seq_len(P_hi))
  bold_hi <- ev_hi <- cov_hi <- list()
  for (r in 1:2) {
    t_ev <- cumsum(c(6.0, stats::runif(n_trials - 1L, 4.2, 6.2)))
    xv <- as.numeric(scale(stats::rnorm(n_trials)))
    des <- hrf_event_design_cpp(n_time, tr, 0.0, t_ev, rep(0.0, n_trials), seq_len(n_trials), n_trials, hrf_k)[[1L]] / max(hrf_k)
    f_glob <- stats::rnorm(n_time, sd = 0.28)
    ym <- matrix(100 + stats::rnorm(n_time * P_hi, sd = 0.25), nrow = n_time, ncol = P_hi, dimnames = list(NULL, p_ids))
    # Shared global noise across first 40 parcels
    ym[, 1:40] <- ym[, 1:40] + f_glob
    eta_1 <- stats::rnorm(n_trials, 0, 0.65)
    eta_2 <- 0.50 * eta_1 + stats::rnorm(n_trials, 0, 0.25)
    ym[, 1L] <- ym[, 1L] + as.numeric(des %*% (0.8 + eta_1))
    ym[, 2L] <- ym[, 2L] + as.numeric(des %*% (0.8 + eta_2))
    rn <- paste0("run", r)
    bold_hi[[rn]] <- ym
    ev_hi[[rn]] <- data.frame(onset = t_ev, duration = 0.0, condition = "fb")
    cov_hi[[rn]] <- data.frame(x = xv)
  }
  fh_hi <- estimate_crossrun_hrf(bold_hi, ev_hi, tr, lib, fallback_kernel = 1L, pooling = "voxel", noise_model = "ar1")
  dc_hi <- deconvolve_crossrun_hrf(bold_hi, fh_hi, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev")
  rg_hi <- decon_epoch_regression(
    dc_hi, cov_hi, events = ev_hi, target_condition = "fb",
    lags = c(0, 1), restore_scale = TRUE, return_panel = TRUE
  )
  expect_equal(rg_hi$shared_error_model$mode, "exact")
  # Non-zero Block-PAF shared_noise_vec embeddings despite 40 * 2 = 80 / P = 120 * 2 = 240 > 96
  p1_vecs <- rg_hi$panel[voxel == "P001" & evt_time == 0]$shared_noise_vec
  expect_gt(mean(vapply(p1_vecs, function(u) sum(u^2), numeric(1L))), 1e-4)

  # 2. Voxel-safe lowrank_voxel mode (V = 600 > 512: zero V x V matrix allocated!) + aggregate_decon_panel()
  V_vox <- 600L
  N_par <- 60L # 10 voxels per parcel
  vox_names <- sprintf("v%04d", seq_len(V_vox))
  vox_to_par <- stats::setNames(rep(sprintf("G%02d", seq_len(N_par)), each = 10L), vox_names)

  bold_v <- ev_v <- cov_v <- list()
  for (r in 1:2) {
    t_ev <- cumsum(c(6.0, stats::runif(n_trials - 1L, 4.2, 6.2)))
    xv <- as.numeric(scale(stats::rnorm(n_trials)))
    des <- hrf_event_design_cpp(n_time, tr, 0.0, t_ev, rep(0.0, n_trials), seq_len(n_trials), n_trials, hrf_k)[[1L]] / max(hrf_k)
    f1 <- stats::rnorm(n_time, sd = 0.32)
    ym <- matrix(100 + stats::rnorm(n_time * V_vox, sd = 0.30), nrow = n_time, ncol = V_vox, dimnames = list(NULL, vox_names))
    # Pervasive global factor across first 300 voxels (parcels G01..G30)
    ym[, 1:300] <- ym[, 1:300] + f1
    eta_s <- stats::rnorm(n_trials, 0, 0.70)
    eta_t <- 0.55 * eta_s + stats::rnorm(n_trials, 0, 0.25)
    for (vv in 1:10) ym[, vv] <- ym[, vv] + as.numeric(des %*% (0.7 + eta_s))
    for (vv in 11:20) ym[, vv] <- ym[, vv] + as.numeric(des %*% (0.7 + eta_t))
    rn <- paste0("run", r)
    bold_v[[rn]] <- ym
    ev_v[[rn]] <- data.frame(onset = t_ev, duration = 0.0, condition = "fb")
    cov_v[[rn]] <- data.frame(x = xv)
  }

  fh_v <- estimate_crossrun_hrf(bold_v, ev_v, tr, lib, fallback_kernel = 1L, pooling = "voxel", noise_model = "ar1")
  dc_v <- deconvolve_crossrun_hrf(bold_v, fh_v, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev")

  # Unaggregated voxel fit at V = 600 > 512 triggers lowrank_voxel without any V x V allocation
  rg_vox <- decon_epoch_regression(
    dc_v, cov_v, events = ev_v, target_condition = "fb",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = FALSE, return_panel = TRUE
  )
  expect_null(rg_vox$unit_error_cor)
  expect_equal(rg_vox$shared_error_model$mode, "lowrank_voxel")
  expect_gte(rg_vox$shared_error_model$poet_K, 1L)

  # Direct voxel-to-voxel connectivity works in lowrank_voxel mode (both Operator & Stacked)
  c_vox_op <- decon_connectivity_repair(rg_vox, "v0001", "v0011", method = "mmoms")
  c_vox_st <- decon_connectivity_repair(data.table::rbindlist(list(rg_vox$panel)), "v0001", "v0011", method = "mmoms")
  expect_true(all(is.finite(c_vox_op$group$estimate)))
  expect_true(all(is.finite(c_vox_st$group$estimate)))
  expect_lt(abs(c_vox_op$group$estimate[1L] - c_vox_st$group$estimate[1L]), 0.14)

  # Post-hoc voxel -> parcel aggregation via aggregate_decon_panel() matches pre-aggregated parcel fit!
  rg_pre <- decon_epoch_regression(
    dc_v, cov_v, events = ev_v, target_condition = "fb",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, parcel_groups = vox_to_par, return_panel = TRUE
  )
  rg_post_fit <- aggregate_decon_panel(rg_vox, parcel_groups = vox_to_par)
  dt_post_bare <- aggregate_decon_panel(data.table::rbindlist(list(rg_vox$panel)), parcel_groups = vox_to_par)

  expect_true(is.matrix(rg_post_fit$unit_error_cor))
  expect_equal(dim(rg_post_fit$unit_error_cor), c(N_par, N_par))
  expect_gt(rg_post_fit$unit_error_cor["G01", "G02"], 0.25)

  c_pre <- decon_connectivity_repair(rg_pre, "G01", "G02", method = "mmoms")
  c_post_fit <- decon_connectivity_repair(rg_post_fit, "G01", "G02", method = "mmoms")
  c_post_bare <- decon_connectivity_repair(dt_post_bare, "G01", "G02", method = "mmoms")
  c_post_zero <- decon_connectivity_repair(dt_post_bare, "G01", "G03", method = "mmoms", shared_error = "zero")
  c_post_auto <- decon_connectivity_repair(dt_post_bare, "G01", "G03", method = "mmoms", shared_error = "auto")

  expect_lt(abs(c_post_fit$group[term == "seed_c"]$estimate - c_pre$group[term == "seed_c"]$estimate), 0.08)
  expect_lt(abs(c_post_bare$group[term == "seed_c"]$estimate - c_pre$group[term == "seed_c"]$estimate), 0.10)
  # Null pair G01-G03 shares f1 noise: shared_error = "auto" suppresses false positive relative to "zero"
  expect_lt(abs(c_post_auto$group[term == "seed_c"]$estimate), abs(c_post_zero$group[term == "seed_c"]$estimate))
})

test_that("v3 (#9): decon_connectivity_repair API parity (augmented panel, supplied err_xy, decon_raw gain, modulator)", {
  set.seed(20261124L)
  tr <- 1.0
  n_time <- 160L
  n_trials <- 28L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]
  groups <- stats::setNames(rep(c("R_seed", "R_targ"), each = 3L), paste0("v", 1:6))

  panels <- vector("list", 8L)
  for (s in 1:8) {
    bold_s <- ev_s <- cov_s <- list()
    for (r in 1:2) {
      t_ev <- cumsum(c(6.0, stats::runif(n_trials - 1L, 4.2, 6.5)))
      rpe <- as.numeric(scale(stats::rnorm(n_trials)))
      abs_rpe <- as.numeric(scale(abs(rpe) + 0.3 * stats::rnorm(n_trials)))
      eta_s <- stats::rnorm(n_trials, 0, 0.70)
      # True PPI is driven by `abs_rpe` (+0.45), NOT `rpe` (0.0)!
      eta_t <- (0.40 + 0.45 * abs_rpe) * eta_s + stats::rnorm(n_trials, 0, 0.25)
      des <- hrf_event_design_cpp(n_time, tr, 0.0, t_ev, rep(0.0, n_trials), seq_len(n_trials), n_trials, hrf_k)[[1L]] / max(hrf_k)
      e_sh <- stats::rnorm(n_time, sd = 0.25)
      ym <- cbind(
        replicate(3L, 100 + as.numeric(des %*% (0.8 + 0.2 * rpe + eta_s)) + e_sh + stats::rnorm(n_time, sd = 0.18)),
        replicate(3L, 100 + as.numeric(des %*% (0.8 + 0.2 * rpe + eta_t)) + e_sh + stats::rnorm(n_time, sd = 0.18))
      )
      colnames(ym) <- names(groups)
      rn <- paste0("run", r)
      bold_s[[rn]] <- ym
      ev_s[[rn]] <- data.frame(onset = t_ev, duration = 0.0, condition = "fb")
      cov_s[[rn]] <- data.frame(rpe = rpe, abs_rpe = abs_rpe)
    }
    # Asymmetric regularization across seed vs target voxels to test raw_gain_by_lag_unit equalization
    lam_het <- c(rep(0.15, 3L), rep(1.20, 3L))
    fh <- estimate_crossrun_hrf(bold_s, ev_s, tr, lib, spatial_groups = groups, pooling = "parcel", noise_model = "ar1")
    dc <- deconvolve_crossrun_hrf(bold_s, fh, lambda = lam_het, scaling = "noise_sd", penalty = "sobolev")
    rg <- decon_epoch_regression(
      dc, cov_s, events = ev_s, target_condition = "fb",
      lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE
    )
    rg$panel$id <- sprintf("sub%02d", s)
    rg$panel$custom_err_xy <- 0.04
    panels[[s]] <- rg
  }
  names(panels) <- sprintf("sub%02d", 1:8)

  # 1. Explicit `modulator = "abs_rpe"` with multi-covariate vector & `return_augmented_panel = TRUE`
  res_mod <- decon_connectivity_repair(
    panels, "R_seed", "R_targ",
    covariates = c("rpe", "abs_rpe"), modulator = "abs_rpe",
    test_ppi = TRUE, method = "mmoms", return_augmented_panel = TRUE
  )
  expect_true("seed_c:abs_rpe_c" %in% res_mod$group$term)
  expect_gt(res_mod$group[term == "seed_c:abs_rpe_c"]$estimate, 0.25)
  expect_lt(res_mod$group[term == "seed_c:abs_rpe_c"]$p.value, 0.05)

  # Augmented panel WLS reproduces subject M-MOMS slopes even with dynamic shared error + multi-covariates
  p1_aug <- res_mod$panel_repair[id == "sub01"]
  lm_p1 <- stats::lm(
    targ_c ~ seed_c_repair + rpe_c_repair + abs_rpe_c_repair + `seed_c:abs_rpe_c_repair`,
    data = p1_aug, weights = repair_weight
  )
  expect_equal(unname(stats::coef(lm_p1)[-1L]), res_mod$subjects[id == "sub01"]$estimate, tolerance = 1e-5)

  # 2. `shared_error = "supplied"` with `err_xy_col`
  dt_sup <- data.table::rbindlist(lapply(panels, `[[`, "panel"))
  res_sup <- decon_connectivity_repair(
    dt_sup, "R_seed", "R_targ",
    covariates = c("rpe", "abs_rpe"), modulator = "abs_rpe",
    test_ppi = TRUE, method = "mmoms",
    shared_error = "supplied", err_xy_col = "custom_err_xy"
  )
  expect_true(all(is.finite(res_sup$group$estimate)))
  expect_gt(mean(res_sup$subjects$r_shared_error), 0.05)

  # 3. `y_col = "decon_raw"` raw-gain equalization on heterogeneous-lambda single-subject panel
  res_raw <- decon_connectivity_repair(
    panels[[1L]], "R_seed", "R_targ",
    covariates = c("rpe", "abs_rpe"), y_col = "decon_raw", method = "mmoms"
  )
  expect_true(all(is.finite(res_raw$group$estimate)))
})

test_that("v3 (#6 & #7): whole-brain whitened_sketch handoff (n_vox > 16384) without bold_obs or V x V matrices", {
  set.seed(20261125L)
  tr <- 1.0
  n_time <- 72L
  n_trials <- 10L
  V_wb <- 16400L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]
  vox_ids <- sprintf("w%05d", seq_len(V_wb))
  par_map <- stats::setNames(rep(sprintf("WB%02d", 1:40), each = 410L), vox_ids)

  bold_wb <- ev_wb <- cov_wb <- list()
  for (r in 1:2) {
    t_ev <- seq(6.0, by = 5.6, length.out = n_trials)
    xv <- as.numeric(scale(stats::rnorm(n_trials)))
    des <- hrf_event_design_cpp(n_time, tr, 0.0, t_ev, rep(0.0, n_trials), seq_len(n_trials), n_trials, hrf_k)[[1L]] / max(hrf_k) / max(hrf_k)
    f_glob <- stats::rnorm(n_time, sd = 0.35)
    ym <- matrix(100 + stats::rnorm(n_time * V_wb, sd = 0.35), nrow = n_time, ncol = V_wb, dimnames = list(NULL, vox_ids))
    ym[, 1:4100] <- ym[, 1:4100] + f_glob
    eta_1 <- stats::rnorm(n_trials, 0, 0.70)
    eta_2 <- 0.50 * eta_1 + stats::rnorm(n_trials, 0, 0.25)
    ym[, 1:410] <- sweep(ym[, 1:410, drop = FALSE], 1L, as.numeric(des %*% (0.7 + eta_1)), "+")
    ym[, 411:820] <- sweep(ym[, 411:820, drop = FALSE], 1L, as.numeric(des %*% (0.7 + eta_2)), "+")
    rn <- paste0("run", r)
    bold_wb[[rn]] <- ym
    ev_wb[[rn]] <- data.frame(onset = t_ev, duration = 0.0, condition = "fb")
    cov_wb[[rn]] <- data.frame(x = xv)
  }

  fh_wb <- estimate_crossrun_hrf(bold_wb, ev_wb, tr, lib, fallback_kernel = 1L, pooling = "parcel", spatial_groups = par_map, noise_model = "ols")
  dc_wb <- deconvolve_crossrun_hrf(bold_wb, fh_wb, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev")

  # At V = 16,400 > 16,384, bold_obs is omitted and compact whitened_sketch is cached instead
  expect_null(dc_wb$bold_obs)
  expect_true(is.list(dc_wb$whitened_sketch))
  expect_equal(ncol(dc_wb$whitened_sketch[[1L]]$B_sk), V_wb)

  # Pre-aggregated parcel fit reconstructs parcel residuals from whitened_sketch without `bold`
  rg_wb_par <- decon_epoch_regression(
    dc_wb, cov_wb, events = ev_wb, target_condition = "fb",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, parcel_groups = par_map, return_panel = TRUE
  )
  expect_equal(rg_wb_par$shared_error_model$mode, "exact")
  expect_gt(rg_wb_par$unit_error_cor["WB01", "WB02"], 0.30)

  # Unaggregated whole-brain voxel fit runs in lowrank_voxel mode from whitened_sketch without `bold`
  rg_wb_vox <- decon_epoch_regression(
    dc_wb, cov_wb, events = ev_wb, target_condition = "fb",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = FALSE, return_panel = FALSE, return_epochs = TRUE
  )
  expect_equal(rg_wb_vox$shared_error_model$mode, "lowrank_voxel")
  expect_gte(rg_wb_vox$shared_error_model$poet_K, 1L)
})

test_that("CSV / .csv.gz round-trip restores panel vector list-columns across medusa_read_table, fread, connectivity, and aggregation", {
  set.seed(20261126L)
  tr <- 1.0
  n_time <- 150L
  n_trials <- 24L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]
  groups <- stats::setNames(rep(c("R_seed", "R_targ", "R_aux"), each = 2L), paste0("v", 1:6))

  subj_panels <- vector("list", 3L)
  csv_paths <- character(3L)
  for (s in 1:3) {
    bold_s <- ev_s <- cov_s <- list()
    for (r in 1:2) {
      t_ev <- cumsum(c(5.5, stats::runif(n_trials - 1L, 3.6, 5.8)))
      xv <- as.numeric(scale(stats::rnorm(n_trials)))
      des <- hrf_event_design_cpp(n_time, tr, 0.0, t_ev, rep(0.0, n_trials), seq_len(n_trials), n_trials, hrf_k)[[1L]] / max(hrf_k)
      e_sh <- as.numeric(stats::filter(stats::rnorm(n_time, sd = 0.28), filter = 0.35, method = "recursive"))
      eta_s <- stats::rnorm(n_trials, 0, 0.70)
      eta_t <- (0.45 + 0.35 * xv) * eta_s + stats::rnorm(n_trials, 0, 0.25)
      ym <- cbind(
        replicate(2L, 100 + as.numeric(des %*% (0.8 + eta_s)) + e_sh + stats::rnorm(n_time, sd = 0.20)),
        replicate(2L, 100 + as.numeric(des %*% (0.8 + eta_t)) + e_sh + stats::rnorm(n_time, sd = 0.20)),
        replicate(2L, 100 + stats::rnorm(n_time, sd = 0.25))
      )
      colnames(ym) <- names(groups)
      rn <- paste0("run", r)
      bold_s[[rn]] <- ym
      ev_s[[rn]] <- data.frame(onset = t_ev, duration = 0.0, condition = "fb")
      cov_s[[rn]] <- data.frame(x = xv)
    }
    fh <- estimate_crossrun_hrf(bold_s, ev_s, tr, lib, spatial_groups = groups, pooling = "parcel", noise_model = "ar1")
    dc <- deconvolve_crossrun_hrf(bold_s, fh, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev")
    rg <- decon_epoch_regression(
      dc, cov_s, events = ev_s, target_condition = "fb",
      lags = c(0, 1), restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE
    )
    rg$panel$id <- sprintf("sub%02d", s)
    subj_panels[[s]] <- data.table::copy(rg$panel)
    f_sub <- tempfile(sprintf("panel_sub%02d_", s), fileext = ".csv")
    medusa_write_table(rg$panel, f_sub)
    csv_paths[s] <- f_sub
  }
  on.exit(unlink(csv_paths), add = TRUE)

  dt_mem <- data.table::rbindlist(subj_panels)
  expect_true(all(c("shared_noise_vec", "trial_cor_vec", "lag_noise_vec") %in% names(dt_mem)))

  # 1. medusa_write_table / medusa_read_table (.csv and .csv.gz)
  f_csv <- tempfile("panel_all_", fileext = ".csv")
  f_gz <- tempfile("panel_all_", fileext = ".csv.gz")
  on.exit(unlink(c(f_csv, f_gz)), add = TRUE)
  medusa_write_table(dt_mem, f_csv)
  medusa_write_table(dt_mem, f_gz)

  dt_csv <- medusa_read_table(f_csv)
  dt_gz <- medusa_read_table(f_gz)
  for (cn in c("shared_noise_vec", "trial_cor_vec", "lag_noise_vec")) {
    expect_true(is.list(dt_csv[[cn]]))
    expect_true(is.list(dt_gz[[cn]]))
    expect_equal(lengths(dt_csv[[cn]]), lengths(dt_mem[[cn]]))
    expect_equal(unlist(dt_csv[[cn]]), unlist(dt_mem[[cn]]), tolerance = 1e-6)
    expect_equal(unlist(dt_gz[[cn]]), unlist(dt_mem[[cn]]), tolerance = 1e-6)
  }

  # 2. Raw data.table::fwrite / fread round-trip (bypassing medusa_read_table)
  f_raw <- tempfile("panel_raw_", fileext = ".csv")
  on.exit(unlink(f_raw), add = TRUE)
  data.table::fwrite(dt_mem, f_raw)
  dt_raw_fread <- data.table::fread(f_raw)
  expect_true(is.character(dt_raw_fread$trial_cor_vec))

  res_mem <- decon_connectivity_repair(
    dt_mem, "R_seed", "R_targ", seed_lag = 0, target_lag = 1,
    control_target_seed_lag = TRUE, test_ppi = TRUE
  )
  res_csv <- decon_connectivity_repair(
    dt_csv, "R_seed", "R_targ", seed_lag = 0, target_lag = 1,
    control_target_seed_lag = TRUE, test_ppi = TRUE
  )
  res_fread <- decon_connectivity_repair(
    dt_raw_fread, "R_seed", "R_targ", seed_lag = 0, target_lag = 1,
    control_target_seed_lag = TRUE, test_ppi = TRUE
  )
  res_path_gz <- decon_connectivity_repair(
    f_gz, "R_seed", "R_targ", seed_lag = 0, target_lag = 1,
    control_target_seed_lag = TRUE, test_ppi = TRUE
  )
  res_path_vec <- decon_connectivity_repair(
    csv_paths, "R_seed", "R_targ", seed_lag = 0, target_lag = 1,
    control_target_seed_lag = TRUE, test_ppi = TRUE
  )

  expect_equal(res_csv$group$estimate, res_mem$group$estimate, tolerance = 1e-5)
  expect_equal(res_csv$group$std.error, res_mem$group$std.error, tolerance = 1e-5)
  expect_equal(res_fread$group$estimate, res_mem$group$estimate, tolerance = 1e-5)
  expect_equal(res_path_gz$group$estimate, res_mem$group$estimate, tolerance = 1e-5)
  expect_equal(res_path_vec$group$estimate, res_mem$group$estimate, tolerance = 1e-5)

  # 3. Length-1 numeric & empty NA list-column edge cases via raw fread + aggregate_decon_panel()
  dt_vox_edge <- data.table::data.table(
    id = rep("sub01", 4L),
    run = rep("run1", 4L),
    trial = rep(1:2, each = 2L),
    evt_time = 0,
    voxel = rep(c("v1", "v2"), times = 2L),
    decon_raw = c(1.1, 1.3, -0.4, -0.2),
    decon_demixed = c(1.2, 1.4, -0.5, -0.3),
    se_raw = rep(0.4, 4L),
    se_demixed = rep(0.5, 4L),
    error_ac1 = rep(0.1, 4L),
    tau2_bio = rep(0.3, 4L),
    shared_noise_vec = list(0.25, 0.25, 0.25, 0.25),
    trial_cor_vec = replicate(4L, numeric(0L), simplify = FALSE)
  )
  f_edge <- tempfile("panel_edge_", fileext = ".csv")
  on.exit(unlink(f_edge), add = TRUE)
  data.table::fwrite(dt_vox_edge, f_edge)
  dt_edge_fread <- data.table::fread(f_edge)
  # fread parses length-1 numbers as atomic numeric and empty strings as logical NA
  expect_true(is.numeric(dt_edge_fread$shared_noise_vec))
  expect_true(is.logical(dt_edge_fread$trial_cor_vec))

  pg_edge <- c(v1 = "P1", v2 = "P1")
  agg_mem <- aggregate_decon_panel(dt_vox_edge, parcel_groups = pg_edge)
  agg_fread <- aggregate_decon_panel(dt_edge_fread, parcel_groups = pg_edge)
  agg_path <- aggregate_decon_panel(f_edge, parcel_groups = pg_edge)
  expect_equal(agg_fread$decon_demixed, agg_mem$decon_demixed, tolerance = 1e-6)
  expect_equal(agg_fread$se_demixed, agg_mem$se_demixed, tolerance = 1e-6)
  expect_equal(unlist(agg_fread$shared_noise_vec), unlist(agg_mem$shared_noise_vec), tolerance = 1e-6)
  expect_equal(agg_path$se_demixed, agg_mem$se_demixed, tolerance = 1e-6)
})

test_that("variance_confounds accepts medusa_confounds, medusa_nifti_ts, and stacked data.frames", {
  set.seed(20261127L)
  tr <- 1.0
  n_time <- 140L
  n_trials <- 24L
  lib <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_k <- lib[, 1L, drop = FALSE]
  groups <- stats::setNames(rep(c("R_seed", "R_targ"), each = 3L), paste0("v", 1:6))

  bold_s <- ev_s <- cov_s <- tsv_files <- list()
  for (r in 1:2) {
    rn <- paste0("run", r)
    t_ev <- cumsum(c(6.0, stats::runif(n_trials - 1L, 3.8, 5.8)))
    xv <- as.numeric(scale(stats::rnorm(n_trials)))
    des <- hrf_event_design_cpp(n_time, tr, 0.0, t_ev, rep(0.0, n_trials), seq_len(n_trials), n_trials, hrf_k)[[1L]] / max(hrf_k)

    fd_true <- pmax(0.0, 0.08 + 0.12 * as.numeric(scale(des %*% pmax(0.0, xv))) + 0.04 * stats::rnorm(n_time))
    fd_true[50L] <- 0.85 # censored motion outlier > 0.5 mm
    dvars_true <- pmax(0.4, 1.0 + 0.35 * scale(fd_true)[, 1L] + 0.08 * stats::rnorm(n_time))
    e_mot <- (0.12 + 0.80 * fd_true) * stats::rnorm(n_time)

    eta_s <- stats::rnorm(n_trials, 0, 0.65)
    eta_t <- 0.40 * eta_s + stats::rnorm(n_trials, 0, 0.25)
    ym <- cbind(
      replicate(3L, 100 + as.numeric(des %*% (0.7 + eta_s)) + e_mot + stats::rnorm(n_time, sd = 0.20)),
      replicate(3L, 100 + as.numeric(des %*% (0.7 + eta_t)) + e_mot + stats::rnorm(n_time, sd = 0.20))
    )
    colnames(ym) <- names(groups)
    bold_s[[rn]] <- ym
    ev_s[[rn]] <- data.frame(onset = t_ev, duration = 0.0, condition = "fb")
    cov_s[[rn]] <- data.frame(x = xv)

    conf_df <- data.frame(
      trans_x = stats::rnorm(n_time, 0, 0.02),
      trans_y = stats::rnorm(n_time, 0, 0.02),
      trans_z = stats::rnorm(n_time, 0, 0.02),
      rot_x = stats::rnorm(n_time, 0, 0.002),
      rot_y = stats::rnorm(n_time, 0, 0.002),
      rot_z = stats::rnorm(n_time, 0, 0.002),
      framewise_displacement = c("n/a", sprintf("%.5f", fd_true[-1L])),
      std_dvars = c("n/a", sprintf("%.5f", dvars_true[-1L])),
      stringsAsFactors = FALSE
    )
    f_tsv <- tempfile(sprintf("sub-01_task-mid_run-0%d_", r), fileext = "_desc-confounds_timeseries.tsv")
    utils::write.table(conf_df, file = f_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
    tsv_files[[rn]] <- f_tsv
  }
  on.exit(unlink(unlist(tsv_files)), add = TRUE)

  conf <- read_fmriprep_confounds(
    unlist(tsv_files),
    bold = bold_s,
    strategy = "6p",
    fd_threshold = 0.5
  )
  expect_true(is.list(conf$variance_confounds))
  expect_identical(names(conf$variance_confounds), c("run1", "run2"))
  expect_identical(colnames(conf$variance_confounds$run1), c("framewise_displacement", "std_dvars"))
  expect_equal(nrow(conf$variance_confounds$run1), n_time)
  expect_equal(unname(conf$variance_confounds$run1[1L, "framewise_displacement"]), 0.0)
  expect_gt(unname(conf$variance_confounds$run1[1L, "std_dvars"]), 0.5)
  # Frame 50 was censored via fd_threshold = 0.5, so valid_rows has < n_time TRUE entries
  expect_false(conf$valid_rows$run1[50L])

  fh <- estimate_crossrun_hrf(
    bold_s, ev_s, tr, lib,
    nuisance = conf, spatial_groups = groups, pooling = "parcel", noise_model = "ar1"
  )
  dc <- deconvolve_crossrun_hrf(bold_s, fh, lambda = 0.6, scaling = "noise_sd", penalty = "sobolev")

  rg_conf <- decon_epoch_regression(
    dc, cov_s, events = ev_s, target_condition = "fb",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE,
    variance_confounds = conf
  )
  rg_manual <- decon_epoch_regression(
    dc, cov_s, events = ev_s, target_condition = "fb",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE,
    variance_confounds = conf$variance_confounds
  )
  expect_equal(rg_conf$panel$se_demixed, rg_manual$panel$se_demixed, tolerance = 1e-12)
  expect_equal(rg_conf$shared_error_model$B_dyn_gated, rg_manual$shared_error_model$B_dyn_gated, tolerance = 1e-12)

  # Reversed run order in medusa_confounds aligns by run name automatically
  conf_rev <- conf
  conf_rev$variance_confounds <- conf$variance_confounds[c("run2", "run1")]
  rg_rev <- decon_epoch_regression(
    dc, cov_s, events = ev_s, target_condition = "fb",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE,
    variance_confounds = conf_rev
  )
  expect_equal(rg_rev$panel$se_demixed, rg_manual$panel$se_demixed, tolerance = 1e-12)

  # medusa_nifti_ts wrapper with $confounds works directly
  ts_mock <- structure(list(bold = bold_s, confounds = conf), class = "medusa_nifti_ts")
  rg_ts <- decon_epoch_regression(
    dc, cov_s, events = ev_s, target_condition = "fb",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE,
    variance_confounds = ts_mock
  )
  expect_equal(rg_ts$panel$se_demixed, rg_manual$panel$se_demixed, tolerance = 1e-12)

  # Stacked data.frame with `run` column works directly
  vc_df_stacked <- rbind(
     cbind(data.frame(run = "run1", stringsAsFactors = FALSE), as.data.frame(conf$variance_confounds$run1)),
     cbind(data.frame(run = "run2", stringsAsFactors = FALSE), as.data.frame(conf$variance_confounds$run2))
  )
  rg_df <- decon_epoch_regression(
    dc, cov_s, events = ev_s, target_condition = "fb",
    lags = 0, restore_scale = TRUE, aggregate_by_parcel = TRUE, return_panel = TRUE,
    variance_confounds = vc_df_stacked
  )
  expect_equal(rg_df$panel$se_demixed, rg_manual$panel$se_demixed, tolerance = 1e-12)

  # Mismatched run names raise a clear error immediately
  conf_bad <- conf
  names(conf_bad$variance_confounds) <- c("runA", "runB")
  expect_error(
    decon_epoch_regression(dc, cov_s, events = ev_s, target_condition = "fb", lags = 0, variance_confounds = conf_bad),
    "Run mismatch in 'variance_confounds'"
  )
})
