# Head-to-head benchmark comparing Parametric Canonical GLM, Condition-Modulated FIR,
# Least-Squares-Separate (LSS) Single-Trial Beta-Series (Mumford et al., 2012),
# and MEDuSA across three experimental regimes.

library(medusa)
library(data.table)
library(ggplot2)
library(patchwork)

set.seed(20261115)

tr <- 1.0
n_time <- 260L
n_subj <- 16L
n_vox <- 6L
run_ids <- c("run1", "run2", "run3")
vox_ids <- paste0("v", seq_len(n_vox))
parcel_ids <- setNames(rep("Focal_Parcel", n_vox), vox_ids)

# Candidate manifold: column 1 is canonical (6.0 s, disp 1.0, ratio 6)
candidates <- hrf_library(tr = tr, preset = "medusa")
pars <- attr(candidates, "parameters")
h_canon <- candidates[, 1L]
idx_fast <- which.min(abs(pars$peak_delay - 3.0) + abs(pars$dispersion - 0.8) + ifelse(is.infinite(pars$undershoot_ratio), 0, 0.5))
h_fast <- candidates[, idx_fast]

# Continuous event column helper (fractional bin occupancy / linear impulse split)
source_col <- function(onset, duration, n, tr) {
  left <- (0:(n - 1L)) * tr
  if (duration > 0) {
    pmax(0, pmin(left + tr, onset + duration) - pmax(left, onset)) / tr
  } else {
    pmax(0, 1 - abs(left / tr - onset / tr))
  }
}

conv_open <- function(x, h) {
  as.numeric(stats::convolve(x, rev(h), type = "open"))[seq_along(x)]
}

sim_ar1 <- function(n, rho, sd_marg) {
  as.numeric(stats::arima.sim(list(ar = rho), n = n, sd = sd_marg * sqrt(1 - rho^2)))
}

# AR(1) prewhitened OLS helper for GLM, LSS, and FIR comparators
ar1_whiten_ols <- function(X, Y) {
  # Y: T x V matrix; X: T x P design
  n <- nrow(Y)
  qr0 <- qr(X)
  res0 <- qr.resid(qr0, Y)
  rho_hat <- mean(vapply(seq_len(ncol(Y)), function(v) {
    r <- res0[, v]
    sum(r[-1L] * r[-n]) / max(sum(r^2), 1e-12)
  }, numeric(1)))
  rho_hat <- max(-0.75, min(0.85, rho_hat))
  # Exact AR(1) Cholesky prewhitening
  Xw <- X
  Yw <- Y
  s1 <- sqrt(1 - rho_hat^2)
  Xw[1L, ] <- s1 * X[1L, ]
  Yw[1L, ] <- s1 * Y[1L, ]
  Xw[-1L, ] <- X[-1L, ] - rho_hat * X[-n, ]
  Yw[-1L, ] <- Y[-1L, ] - rho_hat * Y[-n, ]
  coef_mat <- qr.coef(qr(Xw), Yw)
  rowMeans(coef_mat, na.rm = TRUE)
}

# Fit Parametric Canonical GLM on one subject
fit_parametric_glm <- function(sub_data, h_assumed = h_canon) {
  run_betas <- vapply(run_ids, function(r) {
    Y <- sub_data$bold[[r]]
    ev <- sub_data$events[[r]]
    n <- nrow(Y)
    t_grid <- seq(-1, 1, length.out = n)
    fb_rows <- which(ev$condition == "feedback")
    rpe_c <- sub_data$covariates[[r]]$rpe - mean(sub_data$covariates[[r]]$rpe)

    s_fb_main <- rowSums(vapply(fb_rows, function(i) source_col(ev$onset[i], 0, n, tr), numeric(n)))
    s_fb_rpe  <- as.numeric(vapply(fb_rows, function(i) source_col(ev$onset[i], 0, n, tr), numeric(n)) %*% rpe_c)

    cols <- list(
      intercept = rep(1, n),
      trend     = t_grid,
      fb_main   = conv_open(s_fb_main, h_assumed),
      fb_rpe    = conv_open(s_fb_rpe, h_assumed)
    )
    if (any(ev$condition == "cue")) {
      cue_rows <- which(ev$condition == "cue")
      s_cue <- rowSums(vapply(cue_rows, function(i) source_col(ev$onset[i], 0, n, tr), numeric(n)))
      cols$cue_main <- conv_open(s_cue, h_assumed)
    }
    if (any(ev$condition == "decision")) {
      dec_rows <- which(ev$condition == "decision")
      s_dec <- rowSums(vapply(dec_rows, function(i) source_col(ev$onset[i], ev$duration[i], n, tr), numeric(n)))
      cols$dec_main <- conv_open(s_dec, h_assumed)
    }
    X <- do.call(cbind, cols)
    b <- ar1_whiten_ols(X, Y)
    unname(b["fb_rpe"])
  }, numeric(1))
  mean(run_betas)
}

# Fit Mumford et al. (2012) LSS Single-Trial Beta-Series on one subject
fit_lss_subject <- function(sub_data, id_lbl, h_assumed = h_canon) {
  rbindlist(lapply(run_ids, function(r) {
    Y <- sub_data$bold[[r]]
    ev <- sub_data$events[[r]]
    n <- nrow(Y)
    t_grid <- seq(-1, 1, length.out = n)
    fb_rows <- which(ev$condition == "feedback")
    m <- length(fb_rows)
    fb_cols <- vapply(fb_rows, function(i) conv_open(source_col(ev$onset[i], 0, n, tr), h_assumed), numeric(n))

    nuis <- cbind(intercept = rep(1, n), trend = t_grid)
    if (any(ev$condition == "cue")) {
      cue_rows <- which(ev$condition == "cue")
      s_cue <- rowSums(vapply(cue_rows, function(i) source_col(ev$onset[i], 0, n, tr), numeric(n)))
      nuis <- cbind(nuis, cue_main = conv_open(s_cue, h_assumed))
    }
    if (any(ev$condition == "decision")) {
      dec_rows <- which(ev$condition == "decision")
      s_dec <- rowSums(vapply(dec_rows, function(i) source_col(ev$onset[i], ev$duration[i], n, tr), numeric(n)))
      nuis <- cbind(nuis, dec_main = conv_open(s_dec, h_assumed))
    }

    # Estimate AR(1) once from full-condition GLM on this run, then solve all m LSS models
    X_init <- cbind(nuis, fb_all = rowSums(fb_cols))
    res0 <- qr.resid(qr(X_init), Y)
    rho_hat <- max(-0.75, min(0.85, mean(vapply(seq_len(ncol(Y)), function(v) {
      rv <- res0[, v]; sum(rv[-1L] * rv[-n]) / max(sum(rv^2), 1e-12)
    }, numeric(1)))))
    whiten_mat <- function(M) {
      Mw <- M
      Mw[1L, ] <- sqrt(1 - rho_hat^2) * M[1L, ]
      Mw[-1L, ] <- M[-1L, ] - rho_hat * M[-n, ]
      Mw
    }
    Yw <- whiten_mat(Y)
    nuis_w <- whiten_mat(nuis)
    fb_w <- whiten_mat(fb_cols)
    fb_sum_w <- rowSums(fb_w)

    betas_i <- vapply(seq_len(m), function(i) {
      X_i <- cbind(target = fb_w[, i], other_fb = fb_sum_w - fb_w[, i], nuis_w)
      mean(qr.coef(qr(X_i), Yw)["target", ], na.rm = TRUE)
    }, numeric(1))

    rpe_vec <- sub_data$covariates[[r]]$rpe
    data.table(
      id = id_lbl,
      run = r,
      trial = seq_len(m),
      beta_lss = betas_i,
      rpe = rpe_vec,
      rpe_within = rpe_vec - mean(rpe_vec)
    )
  }))
}

# Fit Condition-Modulated FIR (-2 to +14 s) on one subject
fir_lags <- -2:14
fit_fir_subject <- function(sub_data) {
  run_curves <- vapply(run_ids, function(r) {
    Y <- sub_data$bold[[r]]
    ev <- sub_data$events[[r]]
    n <- nrow(Y)
    t_grid <- seq(-1, 1, length.out = n)
    fb_rows <- which(ev$condition == "feedback")
    rpe_c <- sub_data$covariates[[r]]$rpe - mean(sub_data$covariates[[r]]$rpe)

    B_main <- matrix(0, n, length(fir_lags))
    B_mod  <- matrix(0, n, length(fir_lags))
    for (k in seq_along(fir_lags)) {
      lag_k <- fir_lags[k]
      for (idx_t in seq_along(fb_rows)) {
        t_ev <- ev$onset[fb_rows[idx_t]] + lag_k
        if (t_ev >= 0 && t_ev <= (n - 1L) * tr) {
          sc <- source_col(t_ev, 0, n, tr)
          B_main[, k] <- B_main[, k] + sc
          B_mod[, k]  <- B_mod[, k] + sc * rpe_c[idx_t]
        }
      }
    }
    colnames(B_main) <- paste0("fir_m_", seq_along(fir_lags))
    colnames(B_mod)  <- paste0("fir_rpe_", seq_along(fir_lags))

    nuis <- cbind(intercept = rep(1, n), trend = t_grid)
    if (any(ev$condition == "cue")) {
      cue_rows <- which(ev$condition == "cue")
      s_cue <- rowSums(vapply(cue_rows, function(i) source_col(ev$onset[i], 0, n, tr), numeric(n)))
      nuis <- cbind(nuis, cue_main = conv_open(s_cue, h_canon))
    }
    if (any(ev$condition == "decision")) {
      dec_rows <- which(ev$condition == "decision")
      s_dec <- rowSums(vapply(dec_rows, function(i) source_col(ev$onset[i], ev$duration[i], n, tr), numeric(n)))
      nuis <- cbind(nuis, dec_main = conv_open(s_dec, h_canon))
    }
    X <- cbind(nuis, B_main, B_mod)
    # Small ridge stabilization in AR(1)-whitened space for overlapping multi-lag sticks
    res0 <- qr.resid(qr(cbind(nuis, rowSums(B_main))), Y)
    rho_hat <- max(-0.75, min(0.85, mean(vapply(seq_len(ncol(Y)), function(v) {
      rv <- res0[, v]; sum(rv[-1L] * rv[-n]) / max(sum(rv^2), 1e-12)
    }, numeric(1)))))
    Xw <- X; Yw <- Y
    Xw[1L, ] <- sqrt(1 - rho_hat^2) * X[1L, ]
    Yw[1L, ] <- sqrt(1 - rho_hat^2) * Y[1L, ]
    Xw[-1L, ] <- X[-1L, ] - rho_hat * X[-n, ]
    Yw[-1L, ] <- Y[-1L, ] - rho_hat * Y[-n, ]
    P_pen <- diag(c(rep(1e-6, ncol(nuis)), rep(0.05, 2L * length(fir_lags))))
    b <- solve(crossprod(Xw) + P_pen, crossprod(Xw, rowMeans(Yw)))
    as.numeric(b[colnames(B_mod), 1L])
  }, numeric(length(fir_lags)))
  rowMeans(run_curves)
}

# Simulate one subject under one of the three regimes
simulate_regime_subject <- function(regime, beta_sub) {
  bold <- events <- covariates <- list()
  for (r in run_ids) {
    if (regime == "Regime 1: Slow single-phase (canonical HRF)") {
      m <- 20L
      fb_onsets <- cumsum(c(8.0, runif(m - 1L, 7.5, 10.5)))
      fb_onsets <- fb_onsets[fb_onsets < (n_time - 22L)]
      m <- length(fb_onsets)
      rpe <- as.numeric(scale(rnorm(m)))
      ev <- data.frame(onset = fb_onsets, duration = 0, condition = "feedback", trial_id = seq_len(m))
      S_neural <- vapply(seq_len(m), function(i) source_col(ev$onset[i], 0, n_time, tr) * (0.75 + beta_sub * rpe[i]), numeric(n_time))
      s_true <- rowSums(S_neural)
      h_true <- h_canon
    } else if (regime == "Regime 2: Rapid RL cue-RT-feedback (fast 3.0 s HRF)") {
      m <- 18L
      cue <- cumsum(c(8.0, runif(m - 1L, 11.0, 14.5)))
      cue <- cue[cue < (n_time - 24L)]
      m <- length(cue)
      val <- as.numeric(scale(rnorm(m)))
      outcome <- rbinom(m, 1L, plogis(0.6 * val))
      # Signed TD prediction error correlates negatively with cue expected value
      rpe <- as.numeric(scale(outcome - plogis(0.6 * val) + rnorm(m, sd = 0.35)))
      rt <- pmin(2.2, pmax(0.55, 1.25 - 0.20 * abs(val) + rnorm(m, sd = 0.22)))
      dec <- cue + runif(m, 0.6, 1.1)
      fb <- dec + rt + runif(m, 1.0, 2.5)

      ev <- rbind(
        data.frame(onset = cue, duration = 0, condition = "cue", trial_id = seq_len(m)),
        data.frame(onset = dec, duration = rt, condition = "decision", trial_id = seq_len(m)),
        data.frame(onset = fb,  duration = 0, condition = "feedback", trial_id = seq_len(m))
      )
      ev <- ev[order(ev$onset), ]
      rownames(ev) <- NULL

      s_true <- numeric(n_time)
      for (i in seq_len(nrow(ev))) {
        tid <- ev$trial_id[i]
        amp <- if (ev$condition[i] == "cue") {
          1.00 + 0.35 * val[tid]
        } else if (ev$condition[i] == "decision") {
          0.55
        } else {
          0.85 + beta_sub * rpe[tid]
        }
        s_true <- s_true + source_col(ev$onset[i], ev$duration[i], n_time, tr) * amp
      }
      h_true <- h_fast
    } else {
      # Regime 3: Biphasic neural trajectory (+0.50 phasic onset burst at 0 s, -0.35 regulatory dip across 2..4 s)
      m <- 20L
      fb_onsets <- cumsum(c(8.0, runif(m - 1L, 8.0, 11.0)))
      fb_onsets <- fb_onsets[fb_onsets < (n_time - 24L)]
      m <- length(fb_onsets)
      rpe <- as.numeric(scale(rnorm(m)))
      ev <- data.frame(onset = fb_onsets, duration = 0, condition = "feedback", trial_id = seq_len(m))
      s_true <- numeric(n_time)
      for (i in seq_len(m)) {
        # Early phasic burst (0..1.0 s, +0.50) + late regulatory dip (1.6..3.6 s, -0.45)
        s_true <- s_true +
          source_col(ev$onset[i], 1.0, n_time, tr) * (0.75 + beta_sub * rpe[i]) +
          source_col(ev$onset[i] + 1.6, 2.0, n_time, tr) * (-0.45 * rpe[i])
      }
      h_true <- h_canon
    }

    bold_clean <- conv_open(s_true, h_true)
    # Shared + independent AR(1) voxel noise
    shared_err <- sim_ar1(n_time, rho = 0.55, sd_marg = 0.025)
    Y <- vapply(seq_len(n_vox), function(v) {
      100 + bold_clean + shared_err + sim_ar1(n_time, rho = 0.55, sd_marg = 0.038)
    }, numeric(n_time))
    colnames(Y) <- vox_ids
    bold[[r]] <- Y
    events[[r]] <- ev
    fb_tid <- ev$trial_id[ev$condition == "feedback"]
    covariates[[r]] <- data.frame(rpe = rpe[fb_tid])
  }
  list(bold = bold, events = events, covariates = covariates)
}

regimes <- c(
  "Regime 1: Slow single-phase (canonical HRF)",
  "Regime 2: Rapid RL cue-RT-feedback (fast 3.0 s HRF)",
  "Regime 3: Biphasic burst (+0.50) then late dip (-0.45)"
)

summary_rows <- list()
regime3_curves <- list()

for (reg in regimes) {
  cat("Running:", reg, "\n")
  subj_ids <- sprintf("S%02d", seq_len(n_subj))
  beta_subj <- 0.50 + rnorm(n_subj, sd = 0.10)

  glm_slopes <- numeric(n_subj)
  lss_panels <- vector("list", n_subj)
  fir_mat    <- matrix(NA_real_, nrow = length(fir_lags), ncol = n_subj)
  medusa_eps <- vector("list", n_subj)
  names(medusa_eps) <- subj_ids

  demix_r_val <- 0.05
  lags_eval <- if (grepl("^Regime 3", reg)) -2:6 else -1:3

  for (s_idx in seq_len(n_subj)) {
    sid <- subj_ids[s_idx]
    sub_d <- simulate_regime_subject(reg, beta_subj[s_idx])

    glm_slopes[s_idx] <- fit_parametric_glm(sub_d)
    lss_panels[[s_idx]] <- fit_lss_subject(sub_d, sid)
    fir_mat[, s_idx] <- fit_fir_subject(sub_d)

    hfit <- estimate_crossrun_hrf(
      bold = sub_d$bold, events = sub_d$events, tr = tr,
      kernel_grid = candidates, spatial_groups = parcel_ids,
      pooling = "parcel", amplitude_model = "trial",
      noise_model = "ar1", evidence = "calibrated",
      selection_mode = "soft_average"
    )
    dec <- deconvolve_crossrun_hrf(
      bold = sub_d$bold, hrf_fit = hfit,
      scaling = "noise_sd", penalty = "sobolev",
      lambda_mode = "global", nuisance_mode = "ignore",
      observation_mask = "fit"
    )
    ep <- decon_epoch_regression(
      deconvolution = dec, covariates = sub_d$covariates,
      events = sub_d$events, target_condition = "feedback",
      lags = lags_eval, parcel_groups = parcel_ids,
      demix_other_events = TRUE, demix_target_trials = TRUE,
      demix_within_trial = FALSE, demix_ridge = demix_r_val,
      restore_scale = TRUE, tau2_mode = "marginal",
      return_panel = TRUE
    )
    ep$panel[, id := sid]
    medusa_eps[[sid]] <- ep
  }

  # 1. Parametric GLM group test
  tt_glm <- stats::t.test(glm_slopes)
  # 2. LSS mixed model
  lss_dt <- rbindlist(lss_panels)
  lss_fit <- suppressMessages(lmerTest::lmer(beta_lss ~ rpe_within + (1 + rpe_within || id) + (1 | id:run), data = lss_dt))
  lss_sum <- summary(lss_fit)$coefficients["rpe_within", ]
  lss_ci  <- lss_sum["Estimate"] + c(-1, 1) * stats::qt(0.975, df = lss_sum["df"]) * lss_sum["Std. Error"]

  # 3. FIR peak-BOLD window (+4 to +6 s, scaled by positive-lobe sum of canonical peak ~ 0.22 -> report raw BOLD peak & normalized)
  # Note: Under positive-lobe HRF normalization, peak height of h_canon is max(h_canon); dividing FIR peak window [4..6 s]
  # by mean(h_canon[fir_lags %in% 4:6 + 3]) places FIR peak amplitude on the same unit-impulse scale as GLM/LSS/MEDuSA!
  canon_win_gain <- mean(h_canon[which(((0:(length(h_canon) - 1L)) * tr) %in% 4:6)])
  fir_win_subj <- colMeans(fir_mat[fir_lags %in% 4:6, , drop = FALSE]) / canon_win_gain
  tt_fir <- stats::t.test(fir_win_subj)

  # 4. MEDuSA planned onset mixed_by (lag = 0 s)
  med_panel <- rbindlist(lapply(medusa_eps, `[[`, "panel"))
  med_panel[, rpe_within := rpe - mean(rpe), by = .(id, run, parcel, evt_time)]
  med_lag0 <- med_panel[evt_time == 0]
  med_mb <- mixed_by(
    med_lag0,
    model_formulae = list(fb = decon_demixed ~ rpe_within + (1 + rpe_within || id) + (1 | id:run)),
    weights = "precision_weight_cohort",
    calculate = "parameter_estimates_reml",
    padjust_by = NULL
  )
  mb_row <- subset(med_mb$coef_df_reml, term == "rpe_within")

  summary_rows[[reg]] <- rbind(
    data.table(regime = reg, method = "Parametric GLM (Canonical)", estimate = unname(tt_glm$estimate),
               se = tt_glm$stderr, conf.low = tt_glm$conf.int[1L], conf.high = tt_glm$conf.int[2L], p.value = tt_glm$p.value, truth = 0.50),
    data.table(regime = reg, method = "LSS Beta-Series (Canonical)", estimate = unname(lss_sum["Estimate"]),
               se = unname(lss_sum["Std. Error"]), conf.low = lss_ci[1L], conf.high = lss_ci[2L], p.value = unname(lss_sum["Pr(>|t|)"]), truth = 0.50),
    data.table(regime = reg, method = "FIR (4-6 s BOLD Peak Window)", estimate = unname(tt_fir$estimate),
               se = tt_fir$stderr, conf.low = tt_fir$conf.int[1L], conf.high = tt_fir$conf.int[2L], p.value = tt_fir$p.value, truth = 0.50),
    data.table(regime = reg, method = "MEDuSA (Lag 0 s Onset)", estimate = mb_row$estimate,
               se = mb_row$std.error, conf.low = mb_row$conf.low, conf.high = mb_row$conf.high, p.value = mb_row$p.value, truth = 0.50)
  )

  if (grepl("^Regime 3", reg)) {
    traj <- medusa_multilag_infer(
      medusa_eps,
      covariates = "rpe",
      trajectory_mode = "auto",
      contrast_source = "unsmoothed",
      contrasts = list(late_minus_early = list(early = c(0, 1), late = c(2, 3, 4))),
      tr = tr,
      n_perm = 2000L,
      n_boot_latency = 0L,
      seed = 20261115L
    )
    fir_mean <- rowMeans(fir_mat) / max(h_canon)
    fir_se   <- apply(fir_mat / max(h_canon), 1L, stats::sd) / sqrt(n_subj)
    regime3_curves <- list(
      medusa_coef     = traj$coefficients,
      medusa_contrast = traj$contrasts,
      medusa_latency  = traj$latency,
      fir_curve       = data.table(
        lag = fir_lags,
        estimate = fir_mean,
        conf.low = fir_mean - stats::qt(0.975, n_subj - 1L) * fir_se,
        conf.high = fir_mean + stats::qt(0.975, n_subj - 1L) * fir_se
      ),
      glm_scalar      = unname(tt_glm$estimate),
      lss_scalar      = unname(lss_sum["Estimate"])
    )
  }
}

comp_dt <- rbindlist(summary_rows)
print(comp_dt)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
out_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
} else if (dir.exists("assets")) {
  "assets"
} else {
  "."
}
fwrite(comp_dt, file.path(out_dir, "glm_fir_lss_medusa_summary.csv"))
saveRDS(list(summary = comp_dt, regime3 = regime3_curves), file.path(out_dir, "glm_fir_lss_medusa.rds"))

# Build publication figure
method_levels <- c(
  "Parametric GLM (Canonical)",
  "LSS Beta-Series (Canonical)",
  "FIR (4-6 s BOLD Peak Window)",
  "MEDuSA (Lag 0 s Onset)"
)
comp_dt[, method := factor(method, levels = rev(method_levels))]
comp_dt[, regime_short := factor(
  regime,
  levels = regimes,
  labels = c(
    "A1. Slow Single-Phase\n(Canonical 6 s HRF)",
    "A2. Multi-Phase RL Trial\n(Fast 3.0 s HRF + Cue/RT)",
    "A3. Biphasic Trajectory\n(+0.50 Early, -0.45 Late Dip)"
  )
)]

pal <- c(
  "Parametric GLM (Canonical)"  = "#6A6A6A",
  "LSS Beta-Series (Canonical)" = "#B35806",
  "FIR (4-6 s BOLD Peak Window)" = "#7B52AB",
  "MEDuSA (Lag 0 s Onset)"      = "#1B6CA8"
)

p_top <- ggplot(comp_dt, aes(x = estimate, y = method, color = method)) +
  geom_vline(xintercept = 0.50, linetype = "dashed", color = "#176E66", linewidth = 0.8) +
  geom_vline(xintercept = 0.00, linetype = "dotted", color = "grey65") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0.22, linewidth = 0.9) +
  geom_point(size = 2.9) +
  facet_wrap(~regime_short, nrow = 1) +
  scale_color_manual(values = pal, guide = "none") +
  labs(
    title = "A. Onset RPE Slope Recovery Across Experimental Regimes (True Onset Slope = +0.50, Dashed Line)",
    x = "Estimated Population Feedback RPE Slope (95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 10.5),
    plot.title = element_text(face = "bold", size = 12.5),
    plot.title.position = "plot"
  )

# Bottom panel: Regime 3 biphasic trajectory (MEDuSA neural time vs FIR BOLD time vs scalar LSS)
med_c <- copy(regime3_curves$medusa_coef)
truth_c <- data.table(
  lag = seq(-2, 6, by = 0.1),
  truth = ifelse(seq(-2, 6, by = 0.1) >= 0.0 & seq(-2, 6, by = 0.1) <= 1.0, 0.50,
                 ifelse(seq(-2, 6, by = 0.1) >= 1.6 & seq(-2, 6, by = 0.1) <= 3.6, -0.45, 0.0))
)
fir_c <- copy(regime3_curves$fir_curve)[lag >= -2 & lag <= 12]

p_b1 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey55") +
  geom_hline(yintercept = regime3_curves$lss_scalar, linetype = "dotdash", color = "#B35806", linewidth = 0.8) +
  geom_ribbon(data = med_c, aes(x = lag, ymin = conf.low.sim, ymax = conf.high.sim), fill = "#CFE4F2", alpha = 0.75) +
  geom_step(data = truth_c, aes(x = lag, y = truth, color = "Generating neural truth"), linewidth = 0.9, linetype = "dashed") +
  geom_line(data = med_c, aes(x = lag, y = estimate, color = "MEDuSA (neural event time)"), linewidth = 1.05) +
  geom_point(data = med_c, aes(x = lag, y = estimate, color = "MEDuSA (neural event time)"), size = 2.1) +
  annotate("text", x = 5.9, y = regime3_curves$lss_scalar + 0.13, label = "Scalar LSS / GLM (phases cancel)", color = "#B35806", size = 3.1, hjust = 1) +
  scale_color_manual(values = c("MEDuSA (neural event time)" = "#1B6CA8", "Generating neural truth" = "#176E66")) +
  scale_x_continuous(breaks = -2:6) +
  labs(
    title = "B. Regime 3 Trajectory on Neural Time (MEDuSA)",
    x = "Lag from feedback onset (neural s)",
    y = "RPE slope",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 12), plot.title.position = "plot")

p_b2 <- ggplot(fir_c, aes(x = lag, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey55") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#E4DCF2", alpha = 0.75) +
  geom_line(color = "#7B52AB", linewidth = 1.05) +
  geom_point(color = "#7B52AB", size = 2.0) +
  scale_x_continuous(breaks = seq(-2, 12, by = 2)) +
  labs(
    title = "C. Same Signal in Condition-Modulated FIR",
    x = "Time from feedback onset (delayed BOLD s)",
    y = "Peak-scaled FIR RPE wave"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 12), plot.title.position = "plot")

p_full <- p_top / (p_b1 | p_b2) + plot_layout(heights = c(1.05, 1.05))
ggsave(file.path(out_dir, "glm-fir-lss-medusa.png"), p_full, width = 11.2, height = 7.4, dpi = 175, bg = "white")
cat("Saved comparison benchmark to", out_dir, "\n")
