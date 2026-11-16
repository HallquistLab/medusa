library(medusa)
library(data.table)
probe <- medusa::medusa_fit_mmom(rnorm(20), matrix(rnorm(20), ncol = 1),
                                 array(.01, dim = c(1, 1, 20)))
cat("M-MOM backend:", if (is.null(probe$repair_fit)) "medusa native" else "repair", "\n")
set.seed(6112026)
sims <- 120L
n_subject <- 16L
m <- 60L
noise_sd <- .7
# Unregularized moment sensitivity in the observed trial design. This is a
# diagnostic calculation with a stipulated error cross-moment, not a wrapper mode.
sensitivity_slopes <- function(z, shared) {
  z <- copy(z)
  z[, c("seed", "target") := lapply(.SD, function(v) v - mean(v)), by = run,
    .SDcols = c("seed", "target")]
  X <- cbind(1, z$seed, z$x, z$seed * z$x)
  v <- noise_sd^2 * (1 - 1 / m)
  A <- matrix(0, 4, 4)
  A[2, 2] <- A[4, 4] <- nrow(z) * v
  H <- crossprod(X) - A
  if (min(eigen(H, symmetric = TRUE, only.values = TRUE)$values) <= 0) return(NA_real_)
  rhs <- crossprod(X, z$target)
  rhs[4L] <- rhs[4L] - nrow(z) * shared * noise_sd^2 * (1 - 2 / m)
  solve(H, rhs)[4L]
}
rows <- list()
for (scenario in c("Genuine PPI", "Shared-error null")) {
  beta <- if (scenario == "Genuine PPI") .35 else 0
  corr <- if (scenario == "Shared-error null") .65 else 0
  for (i in seq_len(sims)) {
    z <- CJ(id = seq_len(n_subject), run = 1:2, trial = seq_len(m))
    z[, x := sample(rep(c(-1, 1), each = m / 2)), by = .(id, run)]
    z[, latent_seed := rnorm(.N)]
    z[, latent_target := (.4 * (beta != 0) + (beta + rep(rnorm(n_subject, sd = .11),
      each = 2 * m)) * x) * latent_seed + .2 * x + rnorm(.N, sd = .65)]
    z[, e := rnorm(.N)]
    z[, `:=`(seed = latent_seed + noise_sd * e,
             target = latent_target + noise_sd * (corr * x * e + sqrt(1 - corr^2) * rnorm(.N)))]
    panel <- rbindlist(list(z[, .(id, run, trial, event_index = trial, evt_time = 0,
      parcel = "Seed", decon_demixed = seed, se_demixed = noise_sd, x)],
      z[, .(id, run, trial, event_index = trial, evt_time = 0,
      parcel = "Target", decon_demixed = target, se_demixed = noise_sd, x)]))
    for (method in c("naive_unweighted", "mmoms")) {
      fit <- decon_connectivity_repair(panel, "Seed", "Target", covariates = "x",
        test_ppi = TRUE, orthogonalize_seed_x = FALSE, method = method,
        use_whitened = FALSE, region_col = "parcel")
      f <- fit$group[term == "seed_c:x_c"]
      rows[[length(rows)+1L]] <- f[, .(scenario, iteration = i, truth = beta,
        method = if (method == "mmoms") "M-MOMS" else "Naive",
        estimate, p.value, covers = conf.low <= beta & conf.high >= beta)]
    }
    slopes <- z[, .(slope = sensitivity_slopes(.SD, corr)), by = id]$slope
    t <- t.test(slopes, mu = 0)
    rows[[length(rows)+1L]] <- data.table(scenario, iteration = i, truth = beta,
      method = "Known-error moments", estimate = mean(slopes, na.rm = TRUE),
      p.value = t$p.value, covers = t$conf.int[1] <= beta & t$conf.int[2] >= beta)
  }
}
out <- rbindlist(rows)
stopifnot(nrow(out) == 2 * sims * 3, all(is.finite(out$estimate)))
summary <- out[, .(truth = truth[1], n = .N, mean = mean(estimate),
  bias = mean(estimate - truth), rmse = sqrt(mean((estimate - truth)^2)),
  reject = mean(p.value <= .05), coverage = mean(covers)), by = .(scenario, method)]
fwrite(summary, "ppi_summary.csv")
fwrite(out, "ppi_draws.csv")
print(summary)

# Part 2: End-to-end POET + covariance-regression shared-error & N = 1 cross-trial SE benchmark
set.seed(20261116L)
tr_b <- 1.0
n_time_b <- 200L
n_trials_b <- 26L
lib_b <- hrf_library(tr_b, preset = "medusa")
hrf_b <- as.numeric(hrf_library(tr_b, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6)[, 1L])
hrf_pk <- hrf_b / max(hrf_b)
reg_names_b <- c("Seed", "Target", paste0("Ref", 1:7))
parcels_b <- setNames(rep(reg_names_b, each = 3L), paste0("v", seq_len(27L)))
n_cores_b <- min(6L, max(1L, parallel::detectCores() - 2L))

conv_open_trim <- function(s, h, n) as.numeric(stats::convolve(s, rev(h), type = "open"))[seq_len(n)]
build_imp_des <- function(onsets, w, n, tr, h) {
  s <- numeric(n)
  g <- (0:(n - 1L)) * tr
  for (j in seq_along(onsets)) {
    s <- s + w[j] * pmax(0, 1 - abs(g / tr - onsets[j] / tr))
  }
  conv_open_trim(s, h, n)
}

run_poet_cohort <- function(rep_i) {
  set.seed(20261100L + rep_i)
  n_subj_poet <- 16L
  fits_poet <- fits_leg <- vector("list", n_subj_poet)
  names(fits_poet) <- names(fits_leg) <- sprintf("sub%02d", seq_len(n_subj_poet))
  for (s in seq_len(n_subj_poet)) {
    bold_s <- ev_s <- cov_s <- vc_s <- list()
    for (r in 1:2) {
      itis <- pmin(10.5, 3.6 + stats::rexp(n_trials_b - 1L, rate = 1 / 3.0))
      t_ev <- cumsum(c(6.0, itis))
      t_ev <- 6.0 + (t_ev - 6.0) * min(1.0, (n_time_b - 24.0) / max(max(t_ev) - 6.0, 1.0))
      xv <- as.numeric(scale(stats::rnorm(n_trials_b)))
      t_g <- (0:(n_time_b - 1L)) * tr_b
      e_x0 <- as.numeric(scale(build_imp_des(t_ev, xv, n_time_b, tr_b, hrf_pk)))
      rvt <- pmax(0.0, as.numeric(stats::arima.sim(list(ar = 0.82, ma = 0.25), n = n_time_b)))
      for (ii in seq_along(t_ev)) {
        dt_e <- t_g - t_ev[ii]
        hit <- dt_e >= 6.0 & dt_e <= 15.0
        if (any(hit) && xv[ii] > 0.6) {
          rvt[hit] <- rvt[hit] + 1.8 * xv[ii] * exp(-0.5 * ((dt_e[hit] - 10.5) / 2.5)^2)
        }
      }
      # Pervasive physiological factor across parcels + task/RVT-modulated burst on Seed & Target
      f_glob <- as.numeric(stats::arima.sim(list(ar = 0.45), n = n_time_b, sd = sqrt(pmax(0.02, 0.05 + 0.045 * e_x0))))
      sh_rvt <- 0.25 * rvt * stats::rnorm(n_time_b)

      eta_seed <- stats::rnorm(n_trials_b, 0, 0.85)
      eta_targ <- stats::rnorm(n_trials_b, 0, 0.85) # Zero true coupling & zero true PPI
      sig_seed <- build_imp_des(t_ev, 0.8 + eta_seed, n_time_b, tr_b, hrf_pk)
      sig_targ <- build_imp_des(t_ev, 0.8 + eta_targ, n_time_b, tr_b, hrf_pk)

      ref_blocks <- lapply(1:7, function(k) {
        sig_k <- build_imp_des(t_ev, 0.6 + stats::rnorm(n_trials_b, 0, 0.65), n_time_b, tr_b, hrf_pk)
        replicate(3L, 100 + sig_k + 0.85 * f_glob + stats::rnorm(n_time_b, sd = 0.20))
      })
      ym <- do.call(cbind, c(
        list(
          replicate(3L, 100 + sig_seed + f_glob + sh_rvt + stats::rnorm(n_time_b, sd = 0.20)),
          replicate(3L, 100 + sig_targ + f_glob + sh_rvt + stats::rnorm(n_time_b, sd = 0.20))
        ),
        ref_blocks
      ))
      colnames(ym) <- names(parcels_b)
      rn <- paste0("run", r)
      bold_s[[rn]] <- ym
      ev_s[[rn]] <- data.frame(onset = t_ev, duration = 0.0, condition = "fb")
      cov_s[[rn]] <- data.frame(x = xv)
      vc_s[[rn]] <- cbind(rvt = rvt, null_vc1 = stats::rnorm(n_time_b), null_vc2 = stats::rnorm(n_time_b))
    }
    hf <- estimate_crossrun_hrf(bold_s, ev_s, tr_b, lib_b, spatial_groups = parcels_b,
                                pooling = "parcel", amplitude_model = "condition", noise_model = "ar1",
                                evidence = "calibrated", selection_mode = "soft_average")
    dc <- deconvolve_crossrun_hrf(bold_s, hf, lambda = NULL, rwls = "auto", scaling = "noise_sd", penalty = "sobolev")
    sid <- sprintf("sub%02d", s)
    rg_p <- decon_epoch_regression(dc, cov_s, events = ev_s, target_condition = "fb",
                                   lags = 0, parcel_groups = parcels_b, restore_scale = TRUE,
                                   variance_confounds = vc_s, shared_error_method = "poet_covreg", return_panel = TRUE)
    data.table::set(rg_p$panel, j = "id", value = sid)
    fits_poet[[s]] <- rg_p
    rg_l <- decon_epoch_regression(dc, cov_s, events = ev_s, target_condition = "fb",
                                   lags = 0, parcel_groups = parcels_b, restore_scale = TRUE,
                                   shared_error_method = "legacy", return_panel = TRUE)
    data.table::set(rg_l$panel, j = "id", value = sid)
    fits_leg[[s]] <- rg_l
  }
  csv_tmp <- tempfile(fileext = ".csv")
  medusa:::medusa_write_table(rbindlist(lapply(fits_poet, `[[`, "panel")), csv_tmp)
  dt_csv <- medusa:::medusa_read_table(csv_tmp)
  unlink(csv_tmp)

  m_list <- list(
    "POET + CovReg (Operator)" = decon_connectivity_repair(fits_poet, "Seed", "Target", covariates = "x", test_ppi = TRUE, shared_error = "auto", method = "mmoms", region_col = "parcel"),
    "POET + CovReg (CSV Panel)" = decon_connectivity_repair(dt_csv, "Seed", "Target", covariates = "x", test_ppi = TRUE, shared_error = "auto", method = "mmoms", region_col = "parcel"),
    "Legacy Stationary (auto)" = decon_connectivity_repair(fits_leg, "Seed", "Target", covariates = "x", test_ppi = TRUE, shared_error = "auto", method = "mmoms", region_col = "parcel"),
    "Zero Shared Error" = decon_connectivity_repair(fits_poet, "Seed", "Target", covariates = "x", test_ppi = TRUE, shared_error = "zero", method = "mmoms", region_col = "parcel")
  )
  rows_r <- list()
  for (mn in names(m_list)) {
    for (tm in c("seed_c", "seed_c:x_c")) {
      g <- m_list[[mn]]$group[term == tm][1L]
      rows_r[[length(rows_r) + 1L]] <- data.table(
        regime = "Dynamic Shared-Noise Null (S = 16)", term = tm, method = mn,
        truth = 0.0, estimate = g$estimate, se = g$std.error, df = g$df,
        reject = g$p.value <= 0.05, covers = g$conf.low <= 0 & g$conf.high >= 0
      )
    }
  }
  rbindlist(rows_r)
}

run_n1_rep <- function(rep_i) {
  set.seed(20262100L + rep_i)
  bold_1 <- ev_1 <- cov_1 <- list()
  for (r in 1:2) {
    itis <- pmin(5.0, 1.3 + stats::rexp(31L, rate = 1 / 0.9))
    t_ev <- cumsum(c(5.0, itis))
    n_t1 <- max(120L, as.integer(ceiling(max(t_ev) + 22)))
    xv <- as.numeric(scale(stats::rnorm(32L)))
    f_glob <- as.numeric(stats::arima.sim(list(ar = 0.45), n = n_t1, sd = 0.24))
    sig_s <- build_imp_des(t_ev, 0.8 + stats::rnorm(32L, 0, 0.85), n_t1, tr_b, hrf_pk)
    sig_t <- build_imp_des(t_ev, 0.8 + stats::rnorm(32L, 0, 0.85), n_t1, tr_b, hrf_pk)
    ref_b1 <- lapply(1:7, function(k) {
      sig_k <- build_imp_des(t_ev, 0.6 + stats::rnorm(32L, 0, 0.65), n_t1, tr_b, hrf_pk)
      replicate(3L, 100 + sig_k + 0.85 * f_glob + stats::rnorm(n_t1, sd = 0.20))
    })
    ym <- do.call(cbind, c(
      list(
        replicate(3L, 100 + sig_s + f_glob + stats::rnorm(n_t1, sd = 0.20)),
        replicate(3L, 100 + sig_t + f_glob + stats::rnorm(n_t1, sd = 0.20))
      ),
      ref_b1
    ))
    colnames(ym) <- names(parcels_b)
    rn <- paste0("run", r)
    bold_1[[rn]] <- ym
    ev_1[[rn]] <- data.frame(onset = t_ev, duration = 0.0, condition = "fb")
    cov_1[[rn]] <- data.frame(x = xv)
  }
  hf1 <- estimate_crossrun_hrf(bold_1, ev_1, tr_b, lib_b, spatial_groups = parcels_b,
                               pooling = "parcel", amplitude_model = "condition", noise_model = "ar1",
                               evidence = "calibrated", selection_mode = "soft_average")
  dc1 <- deconvolve_crossrun_hrf(bold_1, hf1, lambda = NULL, rwls = "auto", scaling = "noise_sd", penalty = "sobolev")
  rg1 <- decon_epoch_regression(dc1, cov_1, events = ev_1, target_condition = "fb",
                                lags = 0, parcel_groups = parcels_b, restore_scale = TRUE,
                                shared_error_method = "poet_covreg", return_panel = TRUE)
  data.table::set(rg1$panel, j = "id", value = "sub01")
  f_ct <- decon_connectivity_repair(list(sub01 = rg1), "Seed", "Target", covariates = "x", test_ppi = TRUE,
                                    shared_error = "auto", crosstrial_se = TRUE, region_col = "parcel")$group[term == "seed_c"][1L]
  f_dg <- decon_connectivity_repair(list(sub01 = rg1), "Seed", "Target", covariates = "x", test_ppi = TRUE,
                                    shared_error = "auto", crosstrial_se = FALSE, region_col = "parcel")$group[term == "seed_c"][1L]
  rbindlist(lapply(list(list(nm = "N=1 Cross-Trial SE (crosstrial_se=TRUE)", g = f_ct),
                        list(nm = "N=1 Diagonal SE (crosstrial_se=FALSE)", g = f_dg)), function(pair) {
    data.table(
      regime = "Single-Subject Rapid-ITI Null (N = 1)", term = "seed_c", method = pair$nm,
      truth = 0.0, estimate = pair$g$estimate, se = pair$g$std.error, df = pair$g$df,
      reject = pair$g$p.value <= 0.05, covers = pair$g$conf.low <= 0 & pair$g$conf.high >= 0
    )
  }))
}

poet_rows <- c(
  parallel::mclapply(seq_len(24L), run_poet_cohort, mc.cores = n_cores_b),
  parallel::mclapply(seq_len(120L), run_n1_rep, mc.cores = n_cores_b)
)

poet_dt <- rbindlist(poet_rows)
poet_sum <- poet_dt[, .(
  n = .N, truth = truth[1L], mean = mean(estimate), bias = mean(estimate - truth),
  rmse = sqrt(mean((estimate - truth)^2)), mean_se = mean(se),
  mean_df = mean(ifelse(is.finite(df), df, 999)),
  reject = mean(reject), coverage = mean(covers)
), by = .(regime, term, method)]
fwrite(poet_sum, "poet_covreg_step6_summary.csv")
print(poet_sum)
