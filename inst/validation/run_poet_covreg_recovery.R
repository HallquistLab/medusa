#!/usr/bin/env Rscript
# Pillar 15: Subspace-Debiased Adaptive-POET + Multi-Lag Toeplitz & General-W_s
#            Covariance-Regression Shared-Error Engine in MEDuSA Steps 3-6
#
# Self-contained validation benchmark evaluating `decon_epoch_regression(..., shared_error_method = "poet_covreg", variance_confounds = ...)`
# and `decon_connectivity_repair(..., shared_error = "auto", crosstrial_se = TRUE)` + `aggregate_decon_panel()` across:
#   A1_DynamicRVT_Plus4NullVC_Null : Dynamic task + RVT shared-noise null with 4 pure-null confound columns (LTZ + Schur gate)
#   A2_MotionBurst_AutoRWLS_Null   : Sub-scrubbing cross-parcel motion bursts with Step-2 `rwls = "auto"`
#   A3_GlobalFactor_TruePPI        : True baseline connectivity (0.45) & PPI (0.35) under pervasive global factor noise
#   B1_N1_RapidITI_SharedNull      : Single-subject (N = 1) rapid-ITI cross-trial sandwich SE & Satterthwaite df calibration
#   C1-C3 Scalability & Parity     : P = 120 multi-lag C++ batch, V = 600 `lowrank_voxel` + `aggregate_decon_panel()`, V = 16,500 `whitened_sketch`

Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

if (dir.exists("/tmp/medusa_lib") && !("/tmp/medusa_lib" %in% .libPaths())) {
  .libPaths(c("/tmp/medusa_lib", .libPaths()))
}
if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
  suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(medusa))
}
suppressPackageStartupMessages(library(data.table))

out_csv_dir <- if (dir.exists("inst/validation")) {
  "inst/validation/poet_covreg_results"
} else {
  file.path(getwd(), "poet_covreg_results")
}
dir.create(out_csv_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
n_reps_a <- if (length(args) >= 1L) as.integer(args[1L]) else 24L
n_reps_b <- if (length(args) >= 2L) as.integer(args[2L]) else 120L
n_cores <- if (n_reps_a <= 2L) 1L else min(6L, max(1L, parallel::detectCores() - 2L))

medusa_write_table_fn <- getFromNamespace("medusa_write_table", "medusa")
medusa_read_table_fn <- getFromNamespace("medusa_read_table", "medusa")

conv_open_trim <- function(s, h, n) as.numeric(stats::convolve(s, rev(h), type = "open"))[seq_len(n)]
build_imp_series <- function(onsets, weights, n, tr, h_pk) {
  s <- numeric(n)
  g <- (0:(n - 1L)) * tr
  for (j in seq_along(onsets)) s <- s + weights[j] * pmax(0, 1 - abs(g / tr - onsets[j] / tr))
  conv_open_trim(s, h_pk, n)
}

pack_rep_row <- function(dt_g, scen, term_nm, meth_nm, rep_i, true_v, n_subj) {
  r <- as.data.table(dt_g)[term == term_nm][1L]
  df_v <- if ("df" %in% names(r) && is.finite(r$df) && r$df > 0) as.numeric(r$df) else (if (n_subj > 1L) n_subj - 1.0 else Inf)
  t_crit <- if (is.finite(df_v)) stats::qt(0.975, df = df_v) else 1.959964
  cl <- if ("conf.low" %in% names(r) && is.finite(r$conf.low)) as.numeric(r$conf.low) else as.numeric(r$estimate) - t_crit * as.numeric(r$std.error)
  ch <- if ("conf.high" %in% names(r) && is.finite(r$conf.high)) as.numeric(r$conf.high) else as.numeric(r$estimate) + t_crit * as.numeric(r$std.error)
  data.table(
    scenario = scen, term = term_nm, method = meth_nm, rep = rep_i, n_subj = n_subj,
    true_val = true_v, estimate = as.numeric(r$estimate), std.error = as.numeric(r$std.error),
    df = df_v, p.value = as.numeric(r$p.value), conf.low = cl, conf.high = ch
  )
}

tr_v <- 1.0
lib_v <- hrf_library(tr_v, peak_delays = c(5.0, 6.0), dispersions = 1.0, undershoot_ratios = 6)
hrf_pk <- as.numeric(lib_v[, 2L]) / max(lib_v[, 2L])
parcels_27 <- stats::setNames(rep(c("Seed", "Target", paste0("Ref", 1:7)), each = 3L), paste0("v", seq_len(27L)))

run_one_a <- function(rep_i) {
  set.seed(151000L + rep_i)
  n_subj <- 18L
  out_a <- list()
  for (sc_code in c("A1_DynamicRVT_Plus4NullVC_Null", "A2_MotionBurst_AutoRWLS_Null", "A3_GlobalFactor_TruePPI")) {
    is_true <- identical(sc_code, "A3_GlobalFactor_TruePPI")
    is_jolt <- identical(sc_code, "A2_MotionBurst_AutoRWLS_Null")
    tv_main <- if (is_true) 0.45 else 0.00
    tv_ppi  <- if (is_true) 0.35 else 0.00
    fits_v3 <- fits_leg <- stats::setNames(vector("list", n_subj), sprintf("sub%02d", seq_len(n_subj)))
    for (s in seq_len(n_subj)) {
      bold_s <- ev_s <- cov_s <- vc_s <- list()
      for (r in 1:2) {
        n_tr <- 26L; n_tm <- 200L
        itis <- pmin(10.5, 3.6 + stats::rexp(n_tr - 1L, rate = 1 / 3.0))
        t_ev <- cumsum(c(6.0, itis))
        t_ev <- 6.0 + (t_ev - 6.0) * min(1.0, (n_tm - 24.0) / max(max(t_ev) - 6.0, 1.0))
        xv <- as.numeric(scale(stats::rnorm(n_tr)))
        e_t0 <- as.numeric(scale(build_imp_series(t_ev, rep(1.0, n_tr), n_tm, tr_v, hrf_pk)))
        e_x0 <- as.numeric(scale(build_imp_series(t_ev, xv, n_tm, tr_v, hrf_pk)))
        if (is_jolt) {
          jolt_idx <- sample(seq(18L, n_tm - 12L), 4L)
          v_sh <- rep(0.18, n_tm)
          for (ji in jolt_idx) v_sh[pmin(n_tm, ji:(ji + 1L))] <- 1.65
          e_sh <- as.numeric(stats::filter(sqrt(v_sh) * stats::rnorm(n_tm), filter = 0.30, method = "recursive"))
          vc_r <- NULL
        } else {
          t_g <- (0:(n_tm - 1L)) * tr_v
          rvt <- pmax(0.0, 0.5 + 0.35 * as.numeric(stats::arima.sim(list(ar = 0.80, ma = 0.25), n = n_tm)))
          for (ii in seq_along(t_ev)) {
            dt_e <- t_g - t_ev[ii]
            hit <- dt_e >= 6.0 & dt_e <= 15.0
            if (any(hit) && xv[ii] > 0.5) {
              rvt[hit] <- rvt[hit] + 1.8 * xv[ii] * exp(-0.5 * ((dt_e[hit] - 10.5) / 2.5)^2)
            }
          }
          e_sh <- (0.20 + (if (!is_true) 0.24 * rvt else 0.08 * e_t0)) * stats::rnorm(n_tm)
          vc_r <- cbind(rvt = rvt, matrix(stats::rnorm(n_tm * 4L), nrow = n_tm, ncol = 4L,
                                          dimnames = list(NULL, paste0("null_vc", 1:4))))
        }
        eta_s <- stats::rnorm(n_tr, 0, 0.85)
        eta_t <- if (is_true) (tv_main + tv_ppi * xv) * eta_s + stats::rnorm(n_tr, 0, 0.35) else stats::rnorm(n_tr, 0, 0.85)
        sig_s <- build_imp_series(t_ev, 0.7 + eta_s, n_tm, tr_v, hrf_pk)
        sig_t <- build_imp_series(t_ev, 0.7 + eta_t, n_tm, tr_v, hrf_pk)
        ref_mats <- lapply(1:7, function(k) {
          sig_k <- build_imp_series(t_ev, 0.6 + stats::rnorm(n_tr, 0, 0.65), n_tm, tr_v, hrf_pk)
          replicate(3L, 100 + sig_k + (if (k <= 2) 0.55 else 0.15) * e_sh + stats::rnorm(n_tm, sd = 0.22))
        })
        ym <- do.call(cbind, c(list(
          replicate(3L, 100 + sig_s + e_sh + stats::rnorm(n_tm, sd = 0.22)),
          replicate(3L, 100 + sig_t + e_sh + stats::rnorm(n_tm, sd = 0.22))
        ), ref_mats))
        colnames(ym) <- names(parcels_27)
        rn <- paste0("run", r)
        bold_s[[rn]] <- ym
        ev_s[[rn]] <- data.frame(onset = t_ev, duration = 0.0, condition = "fb")
        cov_s[[rn]] <- data.frame(x = xv)
        if (!is.null(vc_r)) vc_s[[rn]] <- vc_r
      }
      hf <- estimate_crossrun_hrf(bold_s, ev_s, tr_v, lib_v, spatial_groups = parcels_27,
                                  pooling = "parcel", amplitude_model = "trial", noise_model = "ar1")
      dc <- deconvolve_crossrun_hrf(bold_s, hf, lambda = NULL, rwls = "auto", scaling = "noise_sd", penalty = "sobolev")
      sid <- sprintf("sub%02d", s)
      rg_v3 <- decon_epoch_regression(dc, cov_s, events = ev_s, target_condition = "fb", lags = 0,
                                      parcel_groups = parcels_27, restore_scale = TRUE,
                                      variance_confounds = if (length(vc_s)) vc_s else NULL,
                                      shared_error_method = "poet_covreg", return_panel = TRUE)
      data.table::set(rg_v3$panel, j = "id", value = sid)
      fits_v3[[s]] <- rg_v3
      rg_lg <- decon_epoch_regression(dc, cov_s, events = ev_s, target_condition = "fb", lags = 0,
                                      parcel_groups = parcels_27, restore_scale = TRUE,
                                      shared_error_method = "legacy", return_panel = TRUE)
      data.table::set(rg_lg$panel, j = "id", value = sid)
      fits_leg[[s]] <- rg_lg
    }
    dt_stk <- rbindlist(lapply(fits_v3, `[[`, "panel"))
    attr(dt_stk, "shared_error_model") <- NULL
    csv_f <- tempfile(fileext = ".csv.gz")
    medusa_write_table_fn(dt_stk, csv_f)
    dt_csv <- medusa_read_table_fn(csv_f)
    unlink(csv_f)
    meth_map <- list(
      "POET + CovReg (Operator)" = decon_connectivity_repair(fits_v3, "Seed", "Target", covariates = "x", test_ppi = TRUE, shared_error = "auto", method = "mmoms", region_col = "parcel")$group,
      "POET + CovReg (CSV.gz)"   = decon_connectivity_repair(dt_csv,  "Seed", "Target", covariates = "x", test_ppi = TRUE, shared_error = "auto", method = "mmoms", region_col = "parcel")$group,
      "Legacy Stationary (auto)" = decon_connectivity_repair(fits_leg,"Seed", "Target", covariates = "x", test_ppi = TRUE, shared_error = "auto", method = "mmoms", region_col = "parcel")$group,
      "Zero Shared Error"        = decon_connectivity_repair(fits_v3, "Seed", "Target", covariates = "x", test_ppi = TRUE, shared_error = "zero", method = "mmoms", region_col = "parcel")$group
    )
    for (mn in names(meth_map)) {
      out_a[[length(out_a) + 1L]] <- pack_rep_row(meth_map[[mn]], sc_code, "seed_c",     mn, rep_i, tv_main, n_subj)
      out_a[[length(out_a) + 1L]] <- pack_rep_row(meth_map[[mn]], sc_code, "seed_c:x_c", mn, rep_i, tv_ppi,  n_subj)
    }
  }
  rbindlist(out_a)
}

res_a <- if (n_cores > 1L) parallel::mclapply(seq_len(n_reps_a), run_one_a, mc.cores = n_cores) else lapply(seq_len(n_reps_a), run_one_a)
dt_all <- rbindlist(res_a)
sum_dt <- dt_all[, .(
  n_reps = .N, true_val = first(true_val), mean_est = mean(estimate),
  bias = mean(estimate - true_val), rmse = sqrt(mean((estimate - true_val)^2)),
  coverage_95 = mean(conf.low <= true_val & conf.high >= true_val),
  reject_05 = mean(p.value <= 0.05)
), by = .(scenario, term, method)]

fwrite(sum_dt, file.path(out_csv_dir, "poet_covreg_recovery_summary.csv"))
print(sum_dt)
