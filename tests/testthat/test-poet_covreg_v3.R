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

test_that("factored operators, bounded pair batching, compact V_obs_spec, and Kronecker-free Ci preserve exact numerical parity", {
  set.seed(20261208L)
  m <- 96L
  rows <- sort(sample(seq_len(112L), m))
  L_max <- 8L
  d_mat <- abs(outer(rows, rows, "-"))
  adj_idx1 <- which(abs(rows[2:m] - rows[1:(m - 1L)]) == 1L)
  lag_pairs <- lapply(seq_len(L_max), function(d) which(d_mat == d, arr.ind = TRUE))
  f_obs <- cbind(scale(stats::rnorm(m)), scale(stats::rnorm(m)), pmax(0.0, stats::rnorm(m)))

  ar_list <- list(c(0.38), c(0.44, 0.14))
  Q_list <- lapply(c(14L, 18L), function(rk) qr.Q(qr(matrix(stats::rnorm(m * rk), m, rk))))

  # 1. poet_tr_MM_factored_cpp: exact parity on both direct (S <= 8) and blocked dsyrk (S = 40) paths
  Q_big <- lapply(seq_len(40L), function(s) qr.Q(qr(matrix(stats::rnorm(m * (24L + (s %% 5L))), m, 24L + (s %% 5L)))))
  M_big_flat <- vapply(Q_big, function(Qs) as.numeric(diag(m) - tcrossprod(Qs)), numeric(m * m))
  tr_ref_big <- crossprod(M_big_flat)
  tr_fac_small <- poet_tr_MM_factored_cpp(Q_big[1:6], m)
  tr_fac_big <- poet_tr_MM_factored_cpp(Q_big, m)
  expect_lt(max(abs(tr_fac_small - tr_ref_big[1:6, 1:6])), 1e-12)
  expect_lt(max(abs(tr_fac_big - tr_ref_big)), 1e-12)

  # 2. poet_build_pair_X_factored_cpp vs dense poet_build_pair_X_cpp across none, innovation, and woodbury rWLS
  fw_eff <- rep(1.0, m)
  fw_eff[c(12L, 47L, 81L)] <- c(0.008, 0.015, 0.022)
  for (rw_mode in c("none", "innovation", "woodbury")) {
    has_rw <- rw_mode != "none"
    e_rw <- NULL
    if (rw_mode == "woodbury") {
      act <- which(fw_eff < 1.0 - 1e-6)
      e_rw <- matrix(0.0, m, length(act))
      e_rw[cbind(act, seq_along(act))] <- sqrt(1.0 / fw_eff[act] - 1.0)
    }

    dense_M <- dense_G <- dense_H <- fac_ops <- vector("list", 2L)
    for (s in 1:2) {
      Qs <- Q_list[[s]]
      Ms <- diag(m) - tcrossprod(Qs)
      Ws_dense <- hrf_whitening_matrix(rows, ar_list[[s]], unit_marginal = TRUE)
      cs <- sqrt(mean(rowSums(Ws_dense^2)))
      Wt0_dense <- Ws_dense / cs
      Ws_csr <- .scale_csr_rows(.build_whitening_csr(rows, ar_list[[s]], unit_marginal = TRUE), 1.0 / cs)

      if (!has_rw) {
        Gs <- Ms
        Hs <- Ms %*% Wt0_dense
        Wb_csr <- Ws_csr
        VHQ <- poet_csr_transpose_mult_cpp(Wb_csr$row_ptr, Wb_csr$col_idx, Wb_csr$val, m, m, Qs)
        fac_ops[[s]] <- list(
          m = m, Q = Qs, VH_Q = VHQ,
          W_row_ptr = Wb_csr$row_ptr, W_col_idx = Wb_csr$col_idx, W_val = Wb_csr$val,
          has_woodbury = FALSE
        )
      } else if (rw_mode == "innovation") {
        sw <- sqrt(fw_eff)
        Gs <- t(sweep(Ms, 1L, sw, "*"))
        Hs <- Ms %*% sweep(Wt0_dense, 1L, sw, "*")
        Wb_csr <- .scale_csr_rows(Ws_csr, sw)
        VHQ <- poet_csr_transpose_mult_cpp(Wb_csr$row_ptr, Wb_csr$col_idx, Wb_csr$val, m, m, Qs)
        fac_ops[[s]] <- list(
          m = m, Q = Qs, VH_Q = VHQ,
          W_row_ptr = Wb_csr$row_ptr, W_col_idx = Wb_csr$col_idx, W_val = Wb_csr$val,
          has_woodbury = FALSE
        )
      } else {
        we_w <- hrf_apply_whitening(e_rw, rows = rows, ar_vec = ar_list[[s]], unit_marginal = TRUE)
        sv_w <- svd(we_w, nv = 0L)
        uw <- sv_w$u
        fw_fac <- 1.0 / sqrt(1.0 + sv_w$d^2) - 1.0
        app_w <- function(A) A + uw %*% sweep(crossprod(uw, A), 1L, fw_fac, "*")
        Gs <- t(app_w(Ms))
        Hs <- Ms %*% app_w(Wt0_dense)
        VHQ <- poet_csr_transpose_mult_cpp(Ws_csr$row_ptr, Ws_csr$col_idx, Ws_csr$val, m, m, Qs)
        MUw <- uw - Qs %*% crossprod(Qs, uw)
        Pw <- sweep(
          poet_csr_transpose_mult_cpp(Ws_csr$row_ptr, Ws_csr$col_idx, Ws_csr$val, m, m, uw),
          2L, fw_fac, "*"
        )
        fac_ops[[s]] <- list(
          m = m, Q = Qs, VH_Q = VHQ,
          W_row_ptr = Ws_csr$row_ptr, W_col_idx = Ws_csr$col_idx, W_val = Ws_csr$val,
          has_woodbury = TRUE, u_we = uw, MU_we = MUw, fac_we = fw_fac, P_we = Pw
        )
      }
      dense_M[[s]] <- Ms; dense_G[[s]] <- Gs; dense_H[[s]] <- Hs
    }

    f_inn <- if (rw_mode == "innovation") sweep(f_obs, 1L, fw_eff, "*") else f_obs
    for (pair in list(c(1L, 1L), c(1L, 2L), c(2L, 2L))) {
      sa <- pair[1L]; sb <- pair[2L]
      X_den <- poet_build_pair_X_cpp(
        dense_M[[sa]], dense_M[[sb]], dense_G[[sa]], dense_G[[sb]], dense_H[[sa]], dense_H[[sb]],
        adj_idx1, lag_pairs, f_obs, L_max
      )
      X_fac <- poet_build_pair_X_factored_cpp(
        fac_ops[[sa]], fac_ops[[sb]], adj_idx1, lag_pairs, f_inn, f_obs, L_max
      )
      expect_lt(max(abs(X_fac - X_den)), 1e-13)
    }

    # 3. Stage D batch splitting is bitwise invariant across arbitrary chunk sizes (with HAC + non-HAC + VC)
    P_b <- 6L
    E_r <- list(matrix(stats::rnorm(m * P_b), m, P_b), matrix(stats::rnorm(m * P_b), m, P_b))
    X11 <- poet_build_pair_X_factored_cpp(fac_ops[[1L]], fac_ops[[1L]], adj_idx1, lag_pairs, f_inn, f_obs, L_max)
    X12 <- poet_build_pair_X_factored_cpp(fac_ops[[1L]], fac_ops[[2L]], adj_idx1, lag_pairs, f_inn, f_obs, L_max)
    slot_XtX <- list(2.0 * crossprod(X11), 2.0 * crossprod(X12))
    slot_X_runs <- list(list(X11, X11), list(X12, X12))
    L1 <- L_max + 1L; M_dyn <- ncol(f_obs); B_tot <- L1 + 2L * M_dyn
    A_tie <- matrix(0.0, B_tot, L1 + M_dyn)
    A_tie[1:L1, 1:L1] <- diag(L1)
    for (j in seq_len(M_dyn)) { A_tie[L1 + j, L1 + j] <- 0.25; A_tie[L1 + M_dyn + j, L1 + j] <- 0.75 }
    R_p <- stats::cov2cor(crossprod(matrix(stats::rnorm(P_b * P_b), P_b, P_b)) + diag(P_b))
    s0 <- rep(1.1, P_b)
    off_ij <- which(row(R_p) < col(R_p), arr.ind = TRUE)
    n_op <- nrow(off_ij)
    p_slot <- rep(c(1L, 2L), length.out = n_op)
    p_hac <- rep(c(TRUE, FALSE, FALSE), length.out = n_op)
    p_rhs <- rep(FALSE, n_op)
    rhs_z <- matrix(0.0, B_tot, n_op)

    fit_one <- poet_fit_pair_batch_cpp(
      E_r, list(adj_idx1, adj_idx1), list(lag_pairs, lag_pairs),
      slot_XtX, slot_X_runs, A_tie, L_max, L1, M_dyn, 2L, 3L, s0, R_p,
      off_ij[, 1L], off_ij[, 2L], p_slot, p_hac, p_rhs, rhs_z
    )
    chunks <- split(seq_len(n_op), ceiling(seq_len(n_op) / 4L))
    b_cat <- do.call(cbind, lapply(chunks, function(ix) {
      poet_fit_pair_batch_cpp(
        E_r, list(adj_idx1, adj_idx1), list(lag_pairs, lag_pairs),
        slot_XtX, slot_X_runs, A_tie, L_max, L1, M_dyn, 2L, 3L, s0, R_p,
        off_ij[ix, 1L], off_ij[ix, 2L], p_slot[ix], p_hac[ix], p_rhs[ix], rhs_z[, ix, drop = FALSE]
      )$b_dyn
    }))
    v_cat <- do.call(cbind, lapply(chunks, function(ix) {
      poet_fit_pair_batch_cpp(
        E_r, list(adj_idx1, adj_idx1), list(lag_pairs, lag_pairs),
        slot_XtX, slot_X_runs, A_tie, L_max, L1, M_dyn, 2L, 3L, s0, R_p,
        off_ij[ix, 1L], off_ij[ix, 2L], p_slot[ix], p_hac[ix], p_rhs[ix], rhs_z[, ix, drop = FALSE]
      )$v_dyn
    }))
    expect_identical(b_cat, fit_one$b_dyn)
    expect_identical(v_cat, fit_one$v_dyn)
  }

  # 4. poet_build_trial_Ci_cpp matches explicit Kronecker trial covariance construction to < 1e-13
  nC <- 7L; n_sC <- 4L; n_lags <- 3L
  map_s <- c(1L, 2L, 1L, 3L, 4L, 2L, 4L)
  U_ti <- matrix(stats::rnorm(m * n_sC * n_lags), m, n_sC * n_lags)
  V_ti <- matrix(stats::rnorm(m * n_sC * n_lags), m, n_sC * n_lags)
  G_sub <- lapply(seq_len(L_max + 1L), function(l) { A <- matrix(stats::rnorm(nC * nC), nC, nC); 0.5 * (A + t(A)) })
  B_sub <- lapply(1:3, function(j) { A <- matrix(stats::rnorm(nC * nC), nC, nC); 0.5 * (A + t(A)) })
  gm_v <- c(0.8, 1.0, 0.8)
  rsc <- seq(0.85, 1.25, length.out = nC)
  s2_m <- matrix(abs(stats::rnorm(nC * n_lags)) + 0.5, nC, n_lags)

  Ci_cpp <- poet_build_trial_Ci_cpp(U_ti, V_ti, map_s, G_sub, B_sub, lag_pairs, f_obs, gm_v, rsc, s2_m, n_lags, n_sC, L_max)

  ones_l <- matrix(1.0, n_lags, n_lags)
  U_exp <- matrix(0.0, m, nC * n_lags); V_exp <- matrix(0.0, m, nC * n_lags)
  for (la in seq_len(n_lags)) {
    for (ic in seq_len(nC)) {
      U_exp[, (la - 1L) * nC + ic] <- U_ti[, (la - 1L) * n_sC + map_s[ic]]
      V_exp[, (la - 1L) * nC + ic] <- V_ti[, (la - 1L) * n_sC + map_s[ic]]
    }
  }
  Ci_ref <- kronecker(ones_l, G_sub[[1L]]) * crossprod(U_exp)
  for (d in seq_len(L_max)) {
    pd <- lag_pairs[[d]]
    Ci_ref <- Ci_ref + kronecker(ones_l, G_sub[[d + 1L]]) * crossprod(U_exp[pd[, 1L], , drop = FALSE], U_exp[pd[, 2L], , drop = FALSE])
  }
  for (j in 1:3) {
    fj <- f_obs[, j]; gm <- gm_v[j]
    Ci_ref <- Ci_ref + kronecker(ones_l, B_sub[[j]]) * ((1.0 - gm) * crossprod(U_exp, fj * U_exp) + gm * crossprod(V_exp, fj * V_exp))
  }
  Ci_ref <- Ci_ref * tcrossprod(rep(rsc, times = n_lags))
  for (la in seq_len(n_lags)) diag(Ci_ref)[(la - 1L) * nC + seq_len(nC)] <- s2_m[, la]
  expect_lt(max(abs(Ci_cpp - Ci_ref)), 1e-13)

  # 5. Compact V_obs_spec + .sem_get_V_obs exact reconstruction under multi-stratum parcel_mix & Woodbury rWLS
  tr <- 1.0; n_t <- 130L; n_tr <- 20L
  lib <- hrf_library(tr, peak_delays = c(5.0, 6.0), dispersions = 1.0, undershoot_ratios = 6)
  p_grp <- stats::setNames(rep(c("P1", "P2", "P3"), each = 2L), paste0("v", 1:6))
  bold_m <- ev_m <- cov_m <- valid_m <- list()
  for (r in 1:2) {
    rn <- paste0("run", r)
    tev <- cumsum(c(6.0, stats::runif(n_tr - 1L, 4.0, 5.8)))
    xv <- as.numeric(scale(stats::rnorm(n_tr)))
    des <- hrf_event_design_cpp(n_t, tr, 0.0, tev, rep(0.0, n_tr), seq_len(n_tr), n_tr, lib[, 1L, drop = FALSE])[[1L]] / max(lib[, 1L])
    ym <- matrix(100 + stats::rnorm(n_t * 6L, sd = 0.25), n_t, 6L, dimnames = list(NULL, names(p_grp)))
    f_sh <- stats::rnorm(n_t, sd = 0.30)
    ym[, 1:4] <- ym[, 1:4] + f_sh + as.numeric(des %*% (0.7 + 0.3 * xv))
    ym[c(35L, 90L), ] <- ym[c(35L, 90L), ] + c(8.0, -7.5)
    kl <- rep(TRUE, n_t); kl[c(20L, 21L, 70L)] <- FALSE
    valid_m[[rn]] <- kl; bold_m[[rn]] <- ym
    ev_m[[rn]] <- data.frame(onset = tev, duration = 0.0, condition = "fb")
    cov_m[[rn]] <- data.frame(x = xv)
  }
  fh_m <- estimate_crossrun_hrf(bold_m, ev_m, tr, lib, valid_rows = valid_m, pooling = "voxel", noise_model = "ar1")
  dc_m <- deconvolve_crossrun_hrf(bold_m, fh_m, lambda = seq(0.3, 0.9, length.out = 6L), scaling = "noise_sd", penalty = "sobolev", rwls = "auto")
  rg_m <- decon_epoch_regression(
    dc_m, cov_m, events = ev_m, target_condition = "fb",
    lags = c(0, 1), restore_scale = TRUE, aggregate_by_parcel = TRUE, parcel_groups = p_grp, return_panel = TRUE
  )
  rs1 <- rg_m$shared_error_model$runs[[1L]]
  expect_null(rs1$V_obs_strata)
  expect_true(is.list(rs1$V_obs_spec))
  V_full_1 <- .sem_get_V_obs(rs1, 1L, 1L)
  sub_cols <- c(7L, 2L, 14L, 2L)
  expect_lt(max(abs(.sem_get_V_obs(rs1, 1L, 1L, sub_cols) - V_full_1[, sub_cols, drop = FALSE])), 1e-14)

  # Materializing V_obs_strata yields identical decon_connectivity_repair output (< 1e-13)
  rg_dense_v <- rg_m
  for (ri in seq_along(rg_dense_v$shared_error_model$runs)) {
    rs_i <- rg_dense_v$shared_error_model$runs[[ri]]
    S_i <- length(rs_i$U_strata)
    L_i <- length(rs_i$U_strata[[1L]])
    rs_i$V_obs_strata <- lapply(seq_len(S_i), function(si) lapply(seq_len(L_i), function(li) .sem_get_V_obs(rs_i, si, li)))
    rs_i$V_obs_spec <- NULL
    rg_dense_v$shared_error_model$runs[[ri]] <- rs_i
  }
  c_spec <- decon_connectivity_repair(rg_m, "P1", "P2", seed_lag = 0, target_lag = 1, test_ppi = TRUE)
  c_dens <- decon_connectivity_repair(rg_dense_v, "P1", "P2", seed_lag = 0, target_lag = 1, test_ppi = TRUE)
  expect_equal(c_spec$group$estimate, c_dens$group$estimate, tolerance = 1e-13)
  expect_equal(c_spec$group$std.error, c_dens$group$std.error, tolerance = 1e-13)

  # 6. keep_sem_operators = FALSE preserves panel & stacked connectivity identically while releasing U_strata
  rg_no_ops <- decon_epoch_regression(
    dc_m, cov_m, events = ev_m, target_condition = "fb",
    lags = c(0, 1), restore_scale = TRUE, aggregate_by_parcel = TRUE, parcel_groups = p_grp,
    return_panel = TRUE, keep_sem_operators = FALSE
  )
  expect_null(rg_no_ops$shared_error_model$runs[[1L]]$U_strata)
  expect_null(rg_no_ops$shared_error_model$runs[[1L]]$V_obs_spec)
  expect_equal(rg_no_ops$panel$se_demixed, rg_m$panel$se_demixed, tolerance = 1e-14)
  expect_equal(unlist(rg_no_ops$panel$trial_cor_vec), unlist(rg_m$panel$trial_cor_vec), tolerance = 1e-14)
  expect_equal(unlist(rg_no_ops$panel$shared_noise_vec), unlist(rg_m$panel$shared_noise_vec), tolerance = 1e-12)
  c_st_m <- decon_connectivity_repair(data.table::rbindlist(list(rg_m$panel)), "P1", "P2", seed_lag = 0, target_lag = 1, test_ppi = TRUE)
  c_st_no <- decon_connectivity_repair(data.table::rbindlist(list(rg_no_ops$panel)), "P1", "P2", seed_lag = 0, target_lag = 1, test_ppi = TRUE)
  expect_equal(c_st_no$group$estimate, c_st_m$group$estimate, tolerance = 1e-12)
  expect_equal(c_st_no$group$std.error, c_st_m$group$std.error, tolerance = 1e-12)

  # 7. poet_cols_to_list_cpp matches column splitting into plain numeric vectors
  M_split <- matrix(stats::rnorm(5L * 11L), 5L, 11L)
  lst_cpp <- poet_cols_to_list_cpp(M_split)
  expect_length(lst_cpp, 11L)
  expect_null(dim(lst_cpp[[1L]]))
  expect_identical(do.call(cbind, lst_cpp), M_split)

  # 8. shallow_panel_dt strips heavy attributes without mutating caller data.table
  sh_dt <- shallow_panel_dt(rg_m$panel)
  expect_null(attr(sh_dt, "shared_error_model"))
  expect_false(is.null(attr(rg_m$panel, "shared_error_model")))
  orig_names <- copy(names(rg_m$panel))
  sh_dt[, .tmp_col := 1L]
  expect_identical(names(rg_m$panel), orig_names)
  c_panel_attr <- decon_connectivity_repair(rg_m$panel, "P1", "P2", seed_lag = 0, target_lag = 1, test_ppi = TRUE)
  expect_equal(c_panel_attr$group$estimate, c_spec$group$estimate, tolerance = 1e-13)
  expect_equal(c_panel_attr$group$std.error, c_spec$group$std.error, tolerance = 1e-13)
  expect_identical(names(rg_m$panel), orig_names)
})
