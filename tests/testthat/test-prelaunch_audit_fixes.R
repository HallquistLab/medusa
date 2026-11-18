test_that("Step 1 single-run amplitude_model='auto' handles 150 <= n_obs < 240 without crash and rescues zero-mean RPE", {
  set.seed(20261117L)
  tr <- 1.0
  n_obs <- 190L
  lib <- hrf_library(tr = tr, peak_delays = c(4.0, 6.0, 8.0), dispersions = 1.0, undershoot_ratios = 6)
  onsets <- cumsum(c(6.0, stats::runif(31L, 3.6, 6.2)))
  onsets <- onsets[onsets < 168.0]
  m_ev <- length(onsets)
  ev <- list(run1 = data.frame(onset = onsets, duration = 0, condition = "fb"))

  des_fast <- hrf_event_design_cpp(n_obs, tr, 0, onsets, rep(0, m_ev), seq_len(m_ev), m_ev, lib[, 1L, drop = FALSE])[[1L]]
  rpe <- as.numeric(scale(stats::rnorm(m_ev)))
  y_sig <- as.numeric(des_fast %*% (3.5 * rpe))
  Y <- cbind(
    v1 = 100 + y_sig + stats::rnorm(n_obs, sd = 0.12),
    v2 = 100 + y_sig + stats::rnorm(n_obs, sd = 0.12),
    v3 = 100 + y_sig + stats::rnorm(n_obs, sd = 0.12),
    v4 = 100 + y_sig + stats::rnorm(n_obs, sd = 0.12),
    v5 = 100 + y_sig + stats::rnorm(n_obs, sd = 0.12),
    v6 = 100 + stats::rnorm(n_obs, sd = 0.35)
  )
  grp <- c(rep("P_RPE", 5L), "P_Null")

  fit_auto <- expect_no_error(estimate_crossrun_hrf(
    bold = list(run1 = Y),
    events = ev,
    tr = tr,
    kernel_grid = lib,
    spatial_groups = grp,
    pooling = "parcel",
    min_parcel_voxels = 3L,
    amplitude_model = "auto",
    fallback_kernel = 2L,
    noise_pcs = 0L
  ))
  expect_equal(fit_auto$selections$run1$kernel_index[1L], 1L)
  expect_equal(fit_auto$selections$run1$kernel_index[6L], 2L)
})

test_that("Step 1 single-run LOBO CV preserves global condition order under hrf_modulators", {
  set.seed(20261118L)
  tr <- 1.0
  n_obs <- 210L
  lib <- hrf_library(tr = tr, peak_delays = c(4.0, 6.0, 8.0), dispersions = 1.0, undershoot_ratios = 6)
  onsets <- cumsum(c(6.0, stats::runif(35L, 3.6, 6.2)))
  onsets <- onsets[onsets < 186.0]
  m_ev <- length(onsets)
  # Deliberately flip first-appearance condition order across early vs middle/late sub-blocks
  conds <- rep(c("cue", "fb"), length.out = m_ev)
  conds[1L] <- "cue"
  conds[ceiling(m_ev / 3):m_ev] <- rep(c("fb", "cue"), length.out = m_ev - ceiling(m_ev / 3) + 1L)
  mod_v <- as.numeric(scale(stats::rnorm(m_ev)))
  ev <- list(run1 = data.frame(onset = onsets, duration = 0, condition = conds, mod = mod_v))

  des_late <- hrf_event_design_cpp(n_obs, tr, 0, onsets, rep(0, m_ev), seq_len(m_ev), m_ev, lib[, 3L, drop = FALSE])[[1L]]
  amps <- ifelse(conds == "cue", 2.4, -1.0) + ifelse(conds == "fb", 2.2 * mod_v, 0.0)
  y_sig <- as.numeric(des_late %*% amps)
  Y <- cbind(
    v1 = 100 + y_sig + stats::rnorm(n_obs, sd = 0.14),
    v2 = 100 + y_sig + stats::rnorm(n_obs, sd = 0.14),
    v3 = 100 + y_sig + stats::rnorm(n_obs, sd = 0.14),
    v4 = 100 + y_sig + stats::rnorm(n_obs, sd = 0.14),
    v5 = 100 + y_sig + stats::rnorm(n_obs, sd = 0.14)
  )

  fit_mod <- expect_no_error(estimate_crossrun_hrf(
    bold = list(run1 = Y),
    events = ev,
    tr = tr,
    kernel_grid = lib,
    spatial_groups = rep("P1", 5L),
    pooling = "parcel",
    min_parcel_voxels = 3L,
    amplitude_model = "condition",
    hrf_modulators = "mod",
    fallback_kernel = 2L,
    noise_pcs = 0L
  ))
  expect_equal(fit_mod$selections$run1$kernel_index[1L], 3L)
})

test_that("hrf_identify_noise_pool tolerates NA values on valid_rows-masked scrubbed frames", {
  set.seed(20261119L)
  n_t <- 120L
  n_vox <- 28L
  grp <- rep(paste0("BG", 1:7), each = 4L)
  Y1 <- matrix(100 + stats::rnorm(n_t * n_vox, sd = 0.3), n_t, n_vox, dimnames = list(NULL, paste0("v", seq_len(n_vox))))
  Y2 <- matrix(100 + stats::rnorm(n_t * n_vox, sd = 0.3), n_t, n_vox, dimnames = list(NULL, paste0("v", seq_len(n_vox))))
  vr1 <- rep(TRUE, n_t); vr1[c(15L, 16L, 60L)] <- FALSE
  vr2 <- rep(TRUE, n_t); vr2[c(30L, 31L)] <- FALSE
  Y1[!vr1, ] <- NA_real_
  Y2[!vr2, ] <- NA_real_

  pool <- medusa:::hrf_identify_noise_pool(
    bold = list(run1 = Y1, run2 = Y2),
    selections = NULL,
    spatial_groups = grp,
    noise_pool = rep(TRUE, n_vox),
    mode = "auto",
    valid_rows = list(run1 = vr1, run2 = vr2)
  )
  expect_equal(length(pool$pool_idx), n_vox)
})

test_that("multi-run estimate_crossrun_hrf projects DCT drift before scoring so silent parcels enter noise_pcs='auto'", {
  set.seed(20261122L)
  tr <- 1.0
  n_t <- 240L
  lib <- hrf_library(tr = tr, peak_delays = c(4.5, 6.0, 7.5), dispersions = 1.0, undershoot_ratios = 6)
  grp <- c(rep("Act1", 4L), rep(paste0("BG", 1:10), each = 3L))
  n_vox <- length(grp)
  t_idx <- seq(0, n_t - 1L)
  dct2 <- as.numeric(scale(cos(pi * (2 * t_idx + 1) * 2 / (2 * n_t))))
  t_sec <- t_idx * tr
  bold_l <- ev_l <- list()
  for (r in 1:2) {
    onsets <- cumsum(c(6.0, stats::runif(29L, 4.5, 7.5)))
    onsets <- onsets[onsets < 215.0]
    ev_r <- data.frame(onset = onsets, duration = 0, condition = "fb")
    des_f <- hrf_event_design_cpp(n_t, tr, 0, onsets, rep(0, nrow(ev_r)), 1:nrow(ev_r), nrow(ev_r), lib[, 1L, drop = FALSE])[[1L]]
    y_act <- as.numeric(des_f %*% (2.8 + stats::rnorm(nrow(ev_r), sd = 0.4)))
    phys <- 0.35 * sin(2 * pi * 0.21 * t_sec + stats::runif(1L, 0, 2 * pi))
    ym <- matrix(100.0, n_t, n_vox, dimnames = list(NULL, paste0("v", seq_len(n_vox))))
    for (j in seq_len(n_vox)) {
      ym[, j] <- ym[, j] + (if (j <= 4L) y_act else 0.0) +
        0.55 * stats::rnorm(1L, 1.0, 0.2) * dct2 +
        stats::runif(1L, 0.8, 1.2) * phys +
        stats::rnorm(n_t, sd = 0.16)
    }
    rn <- paste0("run", r)
    bold_l[[rn]] <- ym
    ev_l[[rn]] <- ev_r
  }
  fit_h <- estimate_crossrun_hrf(bold_l, ev_l, tr, lib, spatial_groups = grp, pooling = "parcel",
                                 min_parcel_voxels = 3L, fallback_kernel = 2L, highpass = 128.0, noise_pcs = "auto")
  expect_gte(fit_h$settings$noise_pcs_used, 1L)
  expect_gte(length(fit_h$specification$noise_pool_voxels), 24L)
})

test_that("encode/decode trial_res_band preserves exact lag offsets across censored mid-run trials", {
  tr_nums <- c(1L, 3L, 4L, 7L, 8L)
  m <- length(tr_nums)
  Q_m <- diag(m) - matrix(1 / m, m, m)
  B_raw <- diag(c(0.95, 0.90, 0.92, 0.88, 0.93))
  B_raw[2L, 3L] <- B_raw[3L, 2L] <- -0.18 # trial 3 <-> trial 4 (lag 1)
  B_raw[1L, 2L] <- B_raw[2L, 1L] <- -0.06 # trial 1 <-> trial 3 (lag 2)
  B_raw[3L, 4L] <- B_raw[4L, 3L] <- -0.04 # trial 4 <-> trial 7 (lag 3)
  M_conn <- Q_m %*% B_raw %*% Q_m

  enc8 <- medusa:::.encode_trial_res_band(M_conn, max_lag = 8L, tr_nums = tr_nums)
  M_dec8 <- medusa:::.decode_trial_res_band(lapply(seq_len(m), function(j) enc8[, j]), tr_nums = tr_nums)
  expect_equal(M_dec8, M_conn, tolerance = 1e-12)

  enc4 <- medusa:::.encode_trial_res_band(M_conn, max_lag = 4L, tr_nums = tr_nums)
  M_dec4 <- medusa:::.decode_trial_res_band(lapply(seq_len(m), function(j) enc4[, j]), tr_nums = tr_nums)
  expect_equal(M_dec4, M_conn, tolerance = 0.02)
})

test_that("demix_ridge='auto' stratifies parcels in aggregate_by_parcel=TRUE and avoids POET explosion", {
  set.seed(20261120L)
  tr <- 1.0
  lib <- hrf_library(tr = tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)
  hrf_c <- lib[, 1L, drop = FALSE]
  bold_l <- ev_l <- list()
  for (r in 1:2) {
    t_ev <- cumsum(c(6.0, stats::runif(27L, 2.0, 3.2)))
    n_t <- as.integer(ceiling(max(t_ev) + 22.0))
    des_r <- hrf_event_design_cpp(n_t, tr, 0, t_ev, rep(0, 28L), 1:28, 28L, hrf_c)[[1L]]
    th <- stats::rnorm(28L, sd = 1.2)
    ym <- cbind(
      v1 = 100 + as.numeric(des_r %*% th) + stats::rnorm(n_t, sd = 0.04),
      v2 = 100 + as.numeric(des_r %*% th) + stats::rnorm(n_t, sd = 0.04),
      v3 = 100 + stats::rnorm(n_t, sd = 0.35),
      v4 = 100 + stats::rnorm(n_t, sd = 0.35)
    )
    rn <- paste0("run", r)
    bold_l[[rn]] <- ym
    ev_l[[rn]] <- data.frame(onset = t_ev, duration = 0, condition = "fb")
  }
  grp <- c("P_Act", "P_Act", "P_Null", "P_Null")
  hf <- estimate_crossrun_hrf(bold_l, ev_l, tr, lib, spatial_groups = grp, pooling = "parcel",
                              noise_model = "ols", evidence = "legacy", noise_pcs = 0L)
  dc <- deconvolve_crossrun_hrf(bold_l, hf, lambda = 0.5, scaling = "noise_sd", penalty = "sobolev", noise_pcs = 0L)

  rg_par_fix <- decon_epoch_regression(dc, covariates = NULL, events = ev_l, target_condition = "fb", lags = 0,
                                       demix_ridge = 0.05, aggregate_by_parcel = TRUE, parcel_groups = grp, return_panel = TRUE)
  rg_par_auto <- decon_epoch_regression(dc, covariates = NULL, events = ev_l, target_condition = "fb", lags = 0,
                                        demix_ridge = "auto", aggregate_by_parcel = TRUE, parcel_groups = grp, return_panel = TRUE)
  expect_lt(
    mean(rg_par_auto$panel[parcel == "P_Null"]$se_demixed),
    mean(rg_par_fix$panel[parcel == "P_Null"]$se_demixed)
  )

  rg_vox_auto <- expect_no_error(decon_epoch_regression(
    dc, covariates = NULL, events = ev_l, target_condition = "fb", lags = 0,
    demix_ridge = "auto", aggregate_by_parcel = FALSE, whiten_panel = TRUE,
    shared_error_method = "poet_covreg", return_panel = TRUE
  ))
  agg_par <- aggregate_decon_panel(rg_vox_auto, parcel_groups = grp)
  expect_true(all(c("precision_weight_cohort", "precision_weight_op", "se_whitened") %in% names(agg_par$panel)))
})

test_that("mixed_by reweight_conditional_tau2='auto' skips decon_eb, centers whitened runs, and preserves Step-3 conditional weights", {
  set.seed(20261121L)
  n_subj <- 14L
  n_trials <- 24L
  rows_list <- list()
  for (s in seq_len(n_subj)) {
    sid <- sprintf("sub%02d", s)
    b_s <- 0.55 + 0.20 * ((s - mean(seq_len(n_subj))) / 4)
    se_s <- if (s %% 2L == 0L) 0.20 else 0.55
    tau2_true <- 0.16
    for (r in 1:2) {
      rn <- paste0("run", r)
      mu_run <- if (r == 1L) -2.5 else 2.5
      xv <- as.numeric(scale(stats::rnorm(n_trials)))
      eps <- stats::rnorm(n_trials, 0, sqrt(tau2_true + se_s^2))
      y_w <- mu_run + b_s * xv + eps
      rel <- tau2_true / (tau2_true + se_s^2)
      y_eb_w <- mu_run + rel * (y_w - mean(y_w))
      w_cond <- 1 / (tau2_true + se_s^2)
      rows_list[[length(rows_list) + 1L]] <- data.table::data.table(
        id = sid, run = rn, trial = seq_len(n_trials), voxel = "V1", evt_time = 0,
        x = xv, x_whitened = xv,
        decon_demixed = y_w, decon_demixed_whitened = y_w,
        decon_eb_whitened = y_eb_w,
        se_demixed = rep(se_s * 1.35, n_trials),
        se_whitened = rep(se_s, n_trials),
        tau2_bio = rep(tau2_true, n_trials),
        precision_weight_cohort = rep(w_cond, n_trials)
      )
    }
  }
  dt <- data.table::rbindlist(rows_list)

  # 1. EB outcome is NOT double-shrunk under 'auto' (matches FALSE to machine precision)
  mb_eb_auto <- suppressMessages(suppressWarnings(mixed_by(
    dt, outcomes = "decon_eb_whitened", rhs_model_formulae = list(m = ~ x_whitened + (1 | id)),
    split_on = c("voxel", "evt_time"), weights = "precision_weight_cohort",
    reweight_conditional_tau2 = "auto", calculate = "parameter_estimates_reml"
  )))
  mb_eb_false <- suppressMessages(suppressWarnings(mixed_by(
    dt, outcomes = "decon_eb_whitened", rhs_model_formulae = list(m = ~ x_whitened + (1 | id)),
    split_on = c("voxel", "evt_time"), weights = "precision_weight_cohort",
    reweight_conditional_tau2 = FALSE, calculate = "parameter_estimates_reml"
  )))
  expect_equal(
    mb_eb_auto$coef_df_reml[term == "x_whitened"]$std.error,
    mb_eb_false$coef_df_reml[term == "x_whitened"]$std.error,
    tolerance = 1e-12
  )

  # 2. Already-conditional Step-3 weights are preserved under 'auto' even after rbindlist()
  mb_cond_auto <- suppressMessages(suppressWarnings(mixed_by(
    dt, outcomes = "decon_demixed_whitened", rhs_model_formulae = list(m = ~ x_whitened + (1 | id)),
    split_on = c("voxel", "evt_time"), weights = "precision_weight_cohort",
    reweight_conditional_tau2 = "auto", calculate = "parameter_estimates_reml"
  )))
  mb_cond_false <- suppressMessages(suppressWarnings(mixed_by(
    dt, outcomes = "decon_demixed_whitened", rhs_model_formulae = list(m = ~ x_whitened + (1 | id)),
    split_on = c("voxel", "evt_time"), weights = "precision_weight_cohort",
    reweight_conditional_tau2 = FALSE, calculate = "parameter_estimates_reml"
  )))
  expect_equal(
    mb_cond_auto$coef_df_reml[term == "x_whitened"]$estimate,
    mb_cond_false$coef_df_reml[term == "x_whitened"]$estimate,
    tolerance = 1e-12
  )

  # 3. File-backed CSV marginal panel reweights cleanly and centers run means on decon_demixed_whitened
  dt_marg_file <- data.table::copy(dt)
  dt_marg_file[, x_whitened := NULL]
  dt_marg_file[, tau2_bio := mean((decon_demixed_whitened - ave(decon_demixed_whitened, run))^2), by = id]
  dt_marg_file[, precision_weight_cohort := 1 / (tau2_bio + se_whitened^2)]
  f_csv <- tempfile(fileext = ".csv")
  on.exit(unlink(f_csv), add = TRUE)
  data.table::fwrite(dt_marg_file, f_csv)
  mb_file <- suppressMessages(suppressWarnings(mixed_by(
    f_csv, outcomes = "decon_demixed_whitened", rhs_model_formulae = list(m = ~ x + (1 | id)),
    split_on = c("voxel", "evt_time"), weights = "precision_weight_cohort",
    reweight_conditional_tau2 = "auto", calculate = "parameter_estimates_reml"
  )))
  expect_equal(mb_file$coef_df_reml[term == "x"]$estimate, mb_cond_false$coef_df_reml[term == "x_whitened"]$estimate, tolerance = 0.04)
})
