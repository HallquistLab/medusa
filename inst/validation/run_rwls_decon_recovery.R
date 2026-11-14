#!/usr/bin/env Rscript
# Pillar 14: Adaptive rWLS Temporal Frame Weighting in MEDuSA Step 2 Deconvolution
#
# Evaluates `deconvolve_crossrun_hrf(..., rwls = "auto")` against `"none"`, `"innovation"`,
# and `"woodbury"` across 7 realistic fMRI noise/artifact/confounding regimes through the
# full Step 1 -> Step 2 -> Step 3 -> Step 4 MEDuSA pipeline:
#   R1_Clean_Spatial       : Clean spatially correlated runs (rho_spatial = 0.55, AR(3) cascade)
#   R2_SubScrub_Glitches   : Sub-scrubbing observation glitches (1.5% TRs, 4.8-7.2 sigma, incl. task peaks)
#   R3_Scrub_SpinHistory   : Hard-scrubbed motion gaps + 3-TR exponential spin-history recovery transients
#   R4_Physio_AR_Bursts    : Sustained 5-TR recursive AR(3) physiological innovation bursts
#   R5_HeavyNeural_Mismatch: Rapid-ITI 4.5x heavy-tailed neural bursts + 1.1s HRF latency mismatch (guard check)
#   R6_SmallROI_Subcortical: Small subcortical parcel (V = 8 voxels) with sub-scrubbing glitches
#   R7_Connectivity_MMOMS  : Stage-4 REPAIR (M-MOMS) true PPI (P1->P2) + null coupling (P1->P3) under shared glitches

Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(medusa))
}
suppressPackageStartupMessages(library(data.table))

out_csv_dir <- if (dir.exists("inst/validation")) {
  "inst/validation/rwls_results"
} else {
  file.path(getwd(), "rwls_results")
}
dir.create(out_csv_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
n_reps <- if (length(args) >= 1L) as.integer(args[1L]) else 100L
n_cores <- if (n_reps <= 2L) 1L else min(8L, max(1L, parallel::detectCores() - 2L))

hrf_ar_cascade_info_fn <- getFromNamespace("hrf_ar_cascade_info", "medusa")
hrf_event_design_cpp_fn <- getFromNamespace("hrf_event_design_cpp", "medusa")

simulate_ar_series <- function(n, phi, sd_u, innov = NULL) {
  u <- if (is.null(innov)) stats::rnorm(n, sd = sd_u) else innov
  p <- length(phi)
  if (p == 0L) return(u)
  burn <- stats::rnorm(40L, sd = sd_u)
  u_all <- c(burn, u)
  y_all <- numeric(length(u_all))
  for (t in seq_along(u_all)) {
    acc <- u_all[t]
    for (k in seq_len(min(p, t - 1L))) {
      acc <- acc + phi[k] * y_all[t - k]
    }
    y_all[t] <- acc
  }
  y_all[41L:(40L + n)]
}

build_regime_dataset <- function(seed, regime) {
  set.seed(seed)
  tr <- 0.8
  n_time <- 220L
  n_runs <- if (regime == "R7_Connectivity_MMOMS") 3L else 2L
  n_vox <- if (regime == "R6_SmallROI_Subcortical") 8L else 24L

  ar_casc <- c(0.56, 0.14, -0.20)
  info_ar <- hrf_ar_cascade_info_fn(ar_casc, max_lag = 16L)
  phi <- info_ar$phi
  sd_u <- sqrt(info_ar$sig_u2)

  lib <- hrf_library(tr, peak_delays = c(4.5, 5.5, 6.5), dispersions = c(0.9, 1.1), undershoot_ratios = 6)
  h_true_std <- lib[, 3L, drop = FALSE] # peak 5.5s
  h_true_late <- hrf_library(tr, peak_delays = 6.6, dispersions = 1.15, undershoot_ratios = 6)[, 1L, drop = FALSE]

  if (regime == "R6_SmallROI_Subcortical") {
    parcels <- stats::setNames(rep("P1", n_vox), paste0("v", seq_len(n_vox)))
  } else if (regime == "R7_Connectivity_MMOMS") {
    parcels <- stats::setNames(rep(c("P1", "P2", "P3"), each = 8L), paste0("v", seq_len(n_vox)))
  } else {
    parcels <- stats::setNames(rep(c("P1", "P2"), each = n_vox %/% 2L), paste0("v", seq_len(n_vox)))
  }
  coords <- cbind(seq_len(n_vox), rep(1:4, length.out = n_vox), 1)
  rownames(coords) <- names(parcels)

  gains_vox <- ifelse(parcels == "P1", 1.0, seq(0.85, 1.15, length.out = n_vox))
  beta_true <- 0.35
  psi_main_true <- 0.45
  psi_ppi_true <- 0.35

  bold_list <- events_list <- valid_list <- cov_list <- a_true_list <- amp_p1_list <- glitch_trs_list <- list()

  for (r in seq_len(n_runs)) {
    rname <- paste0("run", r)
    if (regime == "R5_HeavyNeural_Mismatch") {
      itis <- stats::runif(34L, min = 2.4, max = 4.4)
      onsets <- cumsum(c(6.0, itis))
      onsets <- onsets[onsets < (n_time * tr - 18.0)]
    } else {
      itis <- stats::runif(24L, min = 4.8, max = 7.6)
      onsets <- cumsum(c(7.0, itis))
      onsets <- onsets[onsets < (n_time * tr - 18.0)]
    }
    m_tr <- length(onsets)
    ev_df <- data.frame(onset = onsets, duration = 0.0, condition = "stim")
    x_mod <- as.numeric(scale(stats::rnorm(m_tr)))
    x_null <- as.numeric(scale(stats::rnorm(m_tr)))

    if (regime == "R5_HeavyNeural_Mismatch") {
      burst_tr <- sample.int(m_tr, max(2L, round(0.15 * m_tr)))
      eps_tr <- stats::rnorm(m_tr, sd = 0.25)
      eps_tr[burst_tr] <- sample(c(-1, 1), length(burst_tr), replace = TRUE) *
        stats::runif(length(burst_tr), 1.8, 2.8)
      amp_p1 <- 0.70 + beta_true * x_mod + eps_tr
      amp_p2 <- 0.55 + beta_true * x_mod + stats::rnorm(m_tr, sd = 0.22)
      amp_p3 <- amp_p2
    } else if (regime == "R7_Connectivity_MMOMS") {
      b_seed <- stats::rnorm(m_tr, sd = 0.60)
      amp_p1 <- 1.35 + beta_true * x_mod + b_seed
      s_c <- b_seed - mean(b_seed)
      amp_p2 <- 0.55 + 0.20 * x_mod + psi_main_true * s_c + psi_ppi_true * (s_c * x_mod) + stats::rnorm(m_tr, sd = 0.25)
      amp_p3 <- 0.55 + 0.20 * x_mod + stats::rnorm(m_tr, sd = 0.25)
    } else {
      amp_p1 <- 0.65 + beta_true * x_mod + stats::rnorm(m_tr, sd = 0.22)
      amp_p2 <- 0.55 + beta_true * x_mod + stats::rnorm(m_tr, sd = 0.22)
      amp_p3 <- amp_p2
    }

    des_std <- hrf_event_design_cpp_fn(
      n_time, tr, 0.0, onsets, rep(0.0, m_tr), seq_len(m_tr), m_tr, h_true_std
    )[[1L]]
    des_late <- hrf_event_design_cpp_fn(
      n_time, tr, 0.0, onsets, rep(0.0, m_tr), seq_len(m_tr), m_tr, h_true_late
    )[[1L]]
    des_neural <- hrf_event_design_cpp_fn(
      n_time, tr, 0.0, onsets, rep(0.0, m_tr), seq_len(m_tr), m_tr, matrix(c(1.0, 0.0), ncol = 1L)
    )[[1L]]

    noise_scale <- if (regime == "R7_Connectivity_MMOMS") 0.35 else 1.0
    g_innov <- stats::rnorm(n_time, sd = noise_scale * sd_u)
    if (regime == "R4_Physio_AR_Bursts") {
      for (bs in c(52L, 138L)) {
        g_innov[bs:(bs + 4L)] <- g_innov[bs:(bs + 4L)] + stats::rnorm(5L, sd = 4.8 * sd_u)
      }
    }
    g_ar <- simulate_ar_series(n_time, phi, noise_scale * sd_u, innov = g_innov)

    Y_r <- matrix(0.0, n_time, n_vox, dimnames = list(NULL, names(parcels)))
    A_r <- matrix(0.0, n_time, n_vox, dimnames = list(NULL, names(parcels)))
    rho_sp <- if (regime == "R1_Clean_Spatial") 0.55 else 0.35
    for (v in seq_len(n_vox)) {
      p_id <- parcels[v]
      base_amp <- if (p_id == "P1") amp_p1 else if (p_id == "P2") amp_p2 else amp_p3
      amp_v <- gains_vox[v] * base_amp
      use_late <- (regime == "R5_HeavyNeural_Mismatch") && (v %% 2L == 0L)
      sig_v <- as.numeric((if (use_late) des_late else des_std) %*% amp_v)
      A_r[, v] <- as.numeric(des_neural %*% amp_v)

      v_innov <- stats::rnorm(n_time, sd = noise_scale * sd_u)
      if (regime == "R4_Physio_AR_Bursts") {
        for (bs in c(52L, 138L)) {
          v_innov[bs:(bs + 4L)] <- v_innov[bs:(bs + 4L)] + stats::rnorm(5L, sd = 2.5 * sd_u)
        }
      }
      e_idio <- simulate_ar_series(n_time, phi, noise_scale * sd_u, innov = v_innov)
      eps_v <- sqrt(rho_sp) * g_ar + sqrt(1.0 - rho_sp) * e_idio
      Y_r[, v] <- 100.0 + seq(-0.6, 0.6, length.out = n_time) + sig_v + eps_v
    }

    valid_r <- rep(TRUE, n_time)
    sp_trs <- integer(0)
    if (regime %in% c("R2_SubScrub_Glitches", "R6_SmallROI_Subcortical", "R7_Connectivity_MMOMS")) {
      peak_trs <- pmin(n_time - 4L, pmax(5L, round(onsets[c(3L, 9L, 16L)] / tr) + 6L))
      rand_trs <- c(33L, 108L, 172L, 173L, 196L)
      sp_trs <- sort(unique(c(peak_trs, peak_trs[2L] + 1L, rand_trs)))
      for (st in sp_trs) {
        g_amp <- sample(c(-1, 1), 1L) * stats::runif(1L, 4.8, 7.2) * noise_scale
        vox_load <- stats::rnorm(n_vox, mean = 0.85, sd = 0.35)
        Y_r[st, ] <- Y_r[st, ] + g_amp * vox_load + stats::rnorm(n_vox, sd = 1.2 * noise_scale)
      }
    } else if (regime == "R3_Scrub_SpinHistory") {
      gap1 <- 60:63
      gap2 <- 142:145
      valid_r[c(gap1, gap2)] <- FALSE
      Y_r[c(gap1, gap2), ] <- Y_r[c(gap1, gap2), ] + 25.0
      for (g_end in c(max(gap1), max(gap2))) {
        rec_trs <- (g_end + 1L):(g_end + 3L)
        sp_trs <- c(sp_trs, rec_trs)
        dec_env <- c(6.2, 4.4, 2.8) * sample(c(-1, 1), 1L)
        vox_spin <- stats::rnorm(n_vox, mean = 0.9, sd = 0.3)
        for (ii in seq_along(rec_trs)) {
          Y_r[rec_trs[ii], ] <- Y_r[rec_trs[ii], ] + dec_env[ii] * vox_spin + stats::rnorm(n_vox, sd = 1.0)
        }
      }
    } else if (regime == "R4_Physio_AR_Bursts") {
      sp_trs <- c(52:56, 138:142)
    }

    bold_list[[rname]] <- Y_r
    events_list[[rname]] <- ev_df
    valid_list[[rname]] <- valid_r
    cov_list[[rname]] <- data.frame(mod = x_mod, null_cov = x_null)
    a_true_list[[rname]] <- A_r
    amp_p1_list[[rname]] <- amp_p1 - mean(amp_p1)
    glitch_trs_list[[rname]] <- sp_trs
  }

  list(
    tr = tr,
    n_time = n_time,
    n_runs = n_runs,
    n_vox = n_vox,
    lib = lib,
    parcels = parcels,
    coords = coords,
    gains_vox = gains_vox,
    beta_true = beta_true,
    psi_main_true = psi_main_true,
    psi_ppi_true = psi_ppi_true,
    bold_list = bold_list,
    events_list = events_list,
    valid_list = valid_list,
    cov_list = cov_list,
    a_true_list = a_true_list,
    amp_p1_list = amp_p1_list,
    glitch_trs_list = glitch_trs_list
  )
}

run_one_rep <- function(seed, regime) {
  ds <- build_regime_dataset(seed, regime)
  fit1 <- estimate_crossrun_hrf(
    bold = ds$bold_list,
    events = ds$events_list,
    tr = ds$tr,
    kernel_grid = if (regime == "R5_HeavyNeural_Mismatch") ds$lib[, 1:2, drop = FALSE] else ds$lib,
    spatial_groups = ds$parcels,
    spatial_coords = ds$coords,
    valid_rows = ds$valid_list,
    pooling = "parcel",
    selection_mode = "soft_average"
  )

  modes <- c("none", "innovation", "woodbury", "auto")
  amp_p1_true <- unlist(ds$amp_p1_list, use.names = FALSE)

  rows_out <- lapply(modes, function(md) {
    dec <- deconvolve_crossrun_hrf(
      bold = ds$bold_list,
      hrf_fit = fit1,
      lambda = NULL,
      scaling = "noise_sd",
      penalty = "sobolev",
      nuisance_mode = "joint",
      rwls = md
    )

    fw_all <- unlist(dec$frame_weights, use.names = FALSE)
    any_down <- any(fw_all < 0.99)
    frac_down <- mean(fw_all < 0.85)
    spaces <- vapply(dec$frame_weights, function(w) {
      if (!any(abs(as.numeric(w) - 1.0) > 1e-6)) return("none")
      sp <- attr(w, "weight_space")
      if (is.null(sp)) "woodbury" else as.character(sp)
    }, character(1L))

    ep_par <- decon_epoch_regression(
      deconvolution = dec,
      covariates = ds$cov_list,
      events = ds$events_list,
      target_condition = "stim",
      lags = 0.0,
      parcel_groups = ds$parcels,
      aggregate_by_parcel = TRUE,
      restore_scale = TRUE,
      return_panel = TRUE
    )
    cf_par <- ep_par$coefficients
    row_p1_mod <- cf_par[cf_par$voxel == "P1" & cf_par$term == "mod", , drop = FALSE]
    row_p1_nul <- cf_par[cf_par$voxel == "P1" & cf_par$term == "null_cov", , drop = FALSE]

    pan_p1 <- ep_par$panel[ep_par$panel$parcel == "P1", , drop = FALSE]
    pan_p1 <- pan_p1[order(pan_p1$run, pan_p1$trial), , drop = FALSE]
    trial_cor <- if (nrow(pan_p1) == length(amp_p1_true)) {
      suppressWarnings(stats::cor(pan_p1$decon_demixed, amp_p1_true, use = "complete.obs"))
    } else NA_real_

    b_p1 <- row_p1_mod$beta_corrected[1L]
    se_p1 <- row_p1_mod$se_corrected[1L]
    df_p1 <- row_p1_mod$df[1L]
    tcrit_p1 <- stats::qt(0.975, df = max(df_p1, 1))
    cov95_p1 <- as.numeric(abs(b_p1 - ds$beta_true) <= tcrit_p1 * se_p1)
    null_rej_p1 <- as.numeric(row_p1_nul$p_value[1L] < 0.05)

    ep_vox <- decon_epoch_regression(
      deconvolution = dec,
      covariates = ds$cov_list,
      events = ds$events_list,
      target_condition = "stim",
      lags = 0.0,
      aggregate_by_parcel = FALSE,
      restore_scale = TRUE,
      return_panel = FALSE
    )
    cf_vox <- ep_vox$coefficients
    sub_vm <- cf_vox[cf_vox$term == "mod", , drop = FALSE]
    ord_v <- match(names(ds$parcels), sub_vm$voxel)
    b_vox <- sub_vm$beta_corrected[ord_v]
    true_vox <- if (regime == "R7_Connectivity_MMOMS") {
      ifelse(ds$parcels == "P1", ds$beta_true, 0.20 * ds$gains_vox)
    } else {
      ds$beta_true * ds$gains_vox
    }
    vox_mse <- mean((b_vox - true_vox)^2, na.rm = TRUE)
    sub_vn <- cf_vox[cf_vox$term == "null_cov", , drop = FALSE]
    null_rej_vox <- mean(sub_vn$p_value < 0.05, na.rm = TRUE)

    main12_est <- ppi12_est <- main13_est <- ppi13_est <- NA_real_
    main13_rej <- ppi13_rej <- r_shared_13 <- NA_real_
    if (regime == "R7_Connectivity_MMOMS") {
      c12 <- decon_connectivity_repair(
        panel_dt = ep_par, seed_region = "P1", target_region = "P2",
        covariates = "mod", test_ppi = TRUE, use_whitened = FALSE,
        orthogonalize_seed_x = FALSE, shared_error = "auto",
        method = "mmoms", region_col = "parcel"
      )$subjects
      c13 <- decon_connectivity_repair(
        panel_dt = ep_par, seed_region = "P1", target_region = "P3",
        covariates = "mod", test_ppi = TRUE, use_whitened = FALSE,
        orthogonalize_seed_x = FALSE, shared_error = "auto",
        method = "mmoms", region_col = "parcel"
      )$subjects
      main12_est <- c12$estimate[c12$term == "seed_c"][1L]
      ppi12_est  <- c12$estimate[c12$term == "seed_c:x_c"][1L]
      main13_est <- c13$estimate[c13$term == "seed_c"][1L]
      main13_rej <- as.numeric(c13$p.value[c13$term == "seed_c"][1L] < 0.05)
      ppi13_est  <- c13$estimate[c13$term == "seed_c:x_c"][1L]
      ppi13_rej  <- as.numeric(c13$p.value[c13$term == "seed_c:x_c"][1L] < 0.05)
      r_shared_13 <- c13$r_shared_error[1L]
    }

    data.table(
      regime = regime,
      mode = md,
      seed = seed,
      lambda = dec$settings$lambda,
      trial_cor = trial_cor,
      any_down = as.numeric(any_down),
      frac_down = frac_down,
      frac_woodbury = mean(spaces == "woodbury"),
      frac_innovation = mean(spaces == "innovation"),
      b_p1 = b_p1,
      se_p1 = se_p1,
      sq_err_p1 = (b_p1 - ds$beta_true)^2,
      cov95_p1 = cov95_p1,
      vox_mse = vox_mse,
      null_rej_p1 = null_rej_p1,
      null_rej_vox = null_rej_vox,
      main12_est = main12_est,
      main12_sq = (main12_est - ds$psi_main_true)^2,
      ppi12_est = ppi12_est,
      ppi12_sq = (ppi12_est - ds$psi_ppi_true)^2,
      main13_est = main13_est,
      main13_sq = main13_est^2,
      main13_rej = main13_rej,
      ppi13_est = ppi13_est,
      ppi13_sq = ppi13_est^2,
      ppi13_rej = ppi13_rej,
      r_shared_13 = r_shared_13
    )
  })
  rbindlist(rows_out)
}

# Run forked mclapply BEFORE any parent-process LAPACK SVD calls
# to avoid macOS Accelerate Grand Central Dispatch fork hazards.
cat("=================================================================\n")
cat(sprintf("Part A: 7-Regime x 4-Mode Monte Carlo Study (%d reps, %d cores)\n", n_reps, n_cores))
cat("=================================================================\n")

regimes <- c(
  "R1_Clean_Spatial",
  "R2_SubScrub_Glitches",
  "R3_Scrub_SpinHistory",
  "R4_Physio_AR_Bursts",
  "R5_HeavyNeural_Mismatch",
  "R6_SmallROI_Subcortical",
  "R7_Connectivity_MMOMS"
)

all_res <- vector("list", length(regimes))
names(all_res) <- regimes
for (rg in regimes) {
  t0 <- proc.time()[3L]
  seeds <- 40000L + match(rg, regimes) * 1000L + seq_len(n_reps)
  res_rg <- if (n_cores > 1L) {
    parallel::mclapply(seeds, run_one_rep, regime = rg, mc.cores = n_cores)
  } else {
    lapply(seeds, run_one_rep, regime = rg)
  }
  all_res[[rg]] <- rbindlist(res_rg)
  cat(sprintf("  Finished %-25s in %.1fs (%d rows)\n", rg, proc.time()[3L] - t0, nrow(all_res[[rg]])))
}
raw_dt <- rbindlist(all_res)

cat("\n=================================================================\n")
cat("Part B: Generating Illustrative Single-Run Traces & Frame Weights\n")
cat("=================================================================\n")

extract_trace_df <- function(seed, regime, panel_label) {
  ds <- build_regime_dataset(seed, regime)
  fit1 <- estimate_crossrun_hrf(
    bold = ds$bold_list, events = ds$events_list, tr = ds$tr,
    kernel_grid = if (regime == "R5_HeavyNeural_Mismatch") ds$lib[, 1:2, drop = FALSE] else ds$lib,
    spatial_groups = ds$parcels, spatial_coords = ds$coords, valid_rows = ds$valid_list,
    pooling = "parcel", selection_mode = "soft_average"
  )
  dec_none <- deconvolve_crossrun_hrf(
    ds$bold_list, fit1, lambda = NULL, scaling = "noise_sd",
    penalty = "sobolev", nuisance_mode = "joint", rwls = "none"
  )
  dec_auto <- deconvolve_crossrun_hrf(
    ds$bold_list, fit1, lambda = NULL, scaling = "noise_sd",
    penalty = "sobolev", nuisance_mode = "joint", rwls = "auto"
  )
  p1_idx <- which(ds$parcels == "P1")
  bold_p1 <- rowMeans(ds$bold_list$run1[, p1_idx, drop = FALSE])
  bold_p1 <- bold_p1 - mean(bold_p1)
  true_p1 <- rowMeans(ds$a_true_list$run1[, p1_idx, drop = FALSE])
  hat_none <- rowMeans(dec_none$activity$run1[, p1_idx, drop = FALSE])
  hat_auto <- rowMeans(dec_auto$activity$run1[, p1_idx, drop = FALSE])
  fw_r1 <- as.numeric(dec_auto$frame_weights$run1)
  w_space <- if (!any(abs(fw_r1 - 1.0) > 1e-6)) {
    "none"
  } else if (identical(attr(dec_auto$frame_weights$run1, "weight_space"), "innovation")) {
    "innovation"
  } else {
    "woodbury"
  }
  onset_trs <- pmin(ds$n_time, pmax(1L, round(ds$events_list$run1$onset / ds$tr) + 1L))

  data.table(
    panel = panel_label,
    regime = regime,
    tr_idx = seq_len(ds$n_time),
    time_s = (seq_len(ds$n_time) - 1L) * ds$tr,
    bold_centered = bold_p1,
    true_neural = as.numeric(scale(true_p1)),
    decon_none = as.numeric(scale(hat_none)),
    decon_auto = as.numeric(scale(hat_auto)),
    frame_weight = fw_r1,
    selected_space = w_space,
    lambda_none = dec_none$settings$lambda,
    lambda_auto = dec_auto$settings$lambda,
    is_artifact_tr = seq_len(ds$n_time) %in% ds$glitch_trs_list$run1,
    is_trial_onset = seq_len(ds$n_time) %in% onset_trs
  )
}

trace_dt <- rbind(
  extract_trace_df(42007L, "R2_SubScrub_Glitches", "Sub-Scrubbing Motion Glitches (5-7 sigma)"),
  extract_trace_df(45007L, "R5_HeavyNeural_Mismatch", "4.5x Heavy Neural Bursts + 1.1s HRF Mismatch")
)
fwrite(trace_dt, file.path(out_csv_dir, "pillar14_rwls_illustrative_traces.csv"))

# Summarize Step 2-3 metrics across all 7 regimes x 4 modes
mode_levels <- c("none", "innovation", "woodbury", "auto")
raw_dt[, mode := factor(mode, levels = mode_levels)]
raw_dt[, regime := factor(regime, levels = regimes)]

sum_dt <- raw_dt[, .(
  n_reps = .N,
  trig_rate = mean(any_down, na.rm = TRUE),
  down_frame_frac = mean(frac_down, na.rm = TRUE),
  woodbury_share = mean(frac_woodbury, na.rm = TRUE),
  innovation_share = mean(frac_innovation, na.rm = TRUE),
  median_lambda = stats::median(lambda, na.rm = TRUE),
  trial_cor = mean(trial_cor, na.rm = TRUE),
  parcel_mean_est = mean(b_p1, na.rm = TRUE),
  parcel_rmse = sqrt(mean(sq_err_p1, na.rm = TRUE)),
  voxel_rmse = sqrt(mean(vox_mse, na.rm = TRUE)),
  se_cal_ratio = mean(se_p1, na.rm = TRUE) / pmax(1e-8, stats::sd(b_p1, na.rm = TRUE)),
  coverage_95 = mean(cov95_p1, na.rm = TRUE),
  null_type1_parcel = mean(null_rej_p1, na.rm = TRUE),
  null_type1_voxel = mean(null_rej_vox, na.rm = TRUE)
), by = .(regime, mode)]

sum_dt[, `:=`(
  parcel_rmse_chg_pct = 100.0 * (parcel_rmse / parcel_rmse[mode == "none"] - 1.0),
  voxel_rmse_chg_pct = 100.0 * (voxel_rmse / voxel_rmse[mode == "none"] - 1.0)
), by = regime]

# Summarize Stage-4 REPAIR M-MOMS Connectivity on R7_Connectivity_MMOMS
conn_sub <- raw_dt[regime == "R7_Connectivity_MMOMS"]
conn_dt <- conn_sub[, .(
  n_reps = .N,
  true_main_est = mean(main12_est, na.rm = TRUE),
  true_main_rmse = sqrt(mean(main12_sq, na.rm = TRUE)),
  true_ppi_est = mean(ppi12_est, na.rm = TRUE),
  true_ppi_rmse = sqrt(mean(ppi12_sq, na.rm = TRUE)),
  null_main_bias = mean(main13_est, na.rm = TRUE),
  null_main_rmse = sqrt(mean(main13_sq, na.rm = TRUE)),
  null_main_fp05 = mean(main13_rej, na.rm = TRUE),
  null_ppi_bias = mean(ppi13_est, na.rm = TRUE),
  null_ppi_rmse = sqrt(mean(ppi13_sq, na.rm = TRUE)),
  null_ppi_fp05 = mean(ppi13_rej, na.rm = TRUE),
  mean_r_shared = mean(r_shared_13, na.rm = TRUE)
), by = mode]
conn_dt[, `:=`(
  true_main_rmse_chg_pct = 100.0 * (true_main_rmse / true_main_rmse[mode == "none"] - 1.0),
  true_ppi_rmse_chg_pct = 100.0 * (true_ppi_rmse / true_ppi_rmse[mode == "none"] - 1.0),
  null_main_rmse_chg_pct = 100.0 * (null_main_rmse / null_main_rmse[mode == "none"] - 1.0)
)]

fwrite(sum_dt, file.path(out_csv_dir, "pillar14_rwls_regime_summary.csv"))
fwrite(conn_dt, file.path(out_csv_dir, "pillar14_rwls_connectivity_summary.csv"))
if (nzchar(Sys.getenv("MEDUSA_RWLS_RDS_PATH", unset = ""))) {
  saveRDS(
    list(
      regime_summary = sum_dt,
      connectivity_summary = conn_dt,
      illustrative_traces = trace_dt,
      replications = raw_dt
    ),
    Sys.getenv("MEDUSA_RWLS_RDS_PATH")
  )
}

cat("\n--- Step 2-3 Summary Across All 7 Regimes ---\n")
print(sum_dt[, .(regime, mode, trig_rate, woodbury_share, innovation_share,
                 median_lambda, trial_cor, parcel_rmse, parcel_rmse_chg_pct,
                 se_cal_ratio, coverage_95, null_type1_parcel)], digits = 3)

cat("\n--- Stage-4 M-MOMS Connectivity Summary (R7_Connectivity_MMOMS) ---\n")
print(conn_dt, digits = 3)
