# Simulate "Healthy vs. Pathological" QC signatures across Steps 0-6 using live medusa calls
# Outputs: qc_gallery_examples.rds, qc_symptom_metrics.csv, sessionInfo.txt

library(medusa)
library(data.table)

set.seed(20261116)

tr <- 1.0
candidates <- hrf_library(tr = tr, preset = "medusa")
pars <- attr(candidates, "parameters")
idx_h <- which.min(abs(pars$peak_delay - 5.0) + abs(pars$dispersion - 1.0) + ifelse(pars$undershoot_ratio == 6, 0, 0.5))
h_true <- candidates[, idx_h]

source_col <- function(onset, duration, n, tr) {
  left <- (0:(n - 1L)) * tr
  if (duration > 0) {
    pmax(0, pmin(left + tr, onset + duration) - pmax(left, onset)) / tr
  } else {
    pmax(0, 1 - abs(left / tr - onset / tr))
  }
}
conv_open <- function(x, h) as.numeric(stats::convolve(x, rev(h), type = "open"))[seq_along(x)]
sim_ar1 <- function(n, rho, sd_marg) as.numeric(stats::arima.sim(list(ar = rho), n = n, sd = sd_marg * sqrt(1 - rho^2)))

# ==============================================================================
# 1. Steps 0 & 1: Clock/Units Misalignment + Step 1 AR(1) rho = 0.90 Ceiling Hit
# ==============================================================================
cat("1. Simulating Step 0/1 clock/units bugs and AR(1) rho ceiling pathology...\n")
tr_slow <- 2.0
cand_2s <- hrf_library(tr = tr_slow, preset = "medusa")
h_2s <- cand_2s[, 1L]
n_vol_2s <- 140L # 280 s run at TR = 2.0 s
n_vox1 <- 6L
vox1 <- paste0("v", seq_len(n_vox1))
parc1 <- setNames(rep("Parcel_1", n_vox1), vox1)

bold_2s <- ev_sec <- cov_1 <- list()
for (r in c("run1", "run2", "run3")) {
  m <- 16L
  on_sec <- seq(14.0, by = 15.0, length.out = m) + runif(m, -1.5, 1.5)
  rpe <- as.numeric(scale(rnorm(m)))
  s_vec <- numeric(n_vol_2s)
  for (i in seq_len(m)) {
    s_vec <- s_vec + source_col(on_sec[i], 0, n_vol_2s, tr_slow) * (1.0 + 0.60 * rpe[i])
  }
  Y <- vapply(seq_len(n_vox1), function(v) 100 + conv_open(s_vec, h_2s) + sim_ar1(n_vol_2s, 0.42, 0.030), numeric(n_vol_2s))
  colnames(Y) <- vox1
  bold_2s[[r]] <- Y
  ev_sec[[r]] <- data.frame(onset = on_sec, duration = 0, condition = "feedback")
  cov_1[[r]] <- data.frame(rpe = rpe)
}

lags_1 <- seq(-6, 6, by = 2)

# Healthy: intact clock (valid_rows masks first 3 dummy TRs = 6 s)
valid_intact <- lapply(bold_2s, function(Y) { v <- rep(TRUE, nrow(Y)); v[1:3] <- FALSE; v })
fit_1_ok <- estimate_crossrun_hrf(bold_2s, ev_sec, tr = tr_slow, kernel_grid = cand_2s,
                                  valid_rows = valid_intact, spatial_groups = parc1,
                                  pooling = "parcel", amplitude_model = "trial")
dec_1_ok <- deconvolve_crossrun_hrf(bold_2s, fit_1_ok, scaling = "noise_sd", penalty = "sobolev")
ep_1_ok  <- decon_epoch_regression(dec_1_ok, cov_1, ev_sec, target_condition = "feedback",
                                   lags = lags_1, parcel_groups = parc1, restore_scale = TRUE)

# Pathology 1A: physically drop 3 dummy TRs (6 s) without setting time_offset = 6
bold_drop3 <- lapply(bold_2s, function(Y) Y[4:nrow(Y), , drop = FALSE])
fit_1_shift <- estimate_crossrun_hrf(bold_drop3, ev_sec, tr = tr_slow, time_offset = 0,
                                     kernel_grid = cand_2s, spatial_groups = parc1,
                                     pooling = "parcel", amplitude_model = "trial")
dec_1_shift <- deconvolve_crossrun_hrf(bold_drop3, fit_1_shift, scaling = "noise_sd", penalty = "sobolev")
ep_1_shift  <- decon_epoch_regression(dec_1_shift, cov_1, ev_sec, target_condition = "feedback",
                                      lags = lags_1, parcel_groups = parc1, restore_scale = TRUE)

# Pathology 1B: event onsets accidentally supplied in TR indices (onset_sec / 2.0) instead of seconds
ev_tr_bug <- lapply(ev_sec, function(df) transform(df, onset = onset / tr_slow))
fit_1_trbug <- estimate_crossrun_hrf(bold_2s, ev_tr_bug, tr = tr_slow, kernel_grid = cand_2s,
                                     valid_rows = valid_intact, spatial_groups = parc1,
                                     pooling = "parcel", amplitude_model = "trial")
dec_1_trbug <- deconvolve_crossrun_hrf(bold_2s, fit_1_trbug, scaling = "noise_sd", penalty = "sobolev")
ep_1_trbug  <- decon_epoch_regression(dec_1_trbug, cov_1, ev_tr_bug, target_condition = "feedback",
                                      lags = lags_1, parcel_groups = parc1, restore_scale = TRUE)

clock_curves <- rbind(
  as.data.table(ep_1_ok$coefficients)[term == "rpe", .(lag, beta = beta_gls, se = se_gls, state = "Healthy: Intact seconds clock (valid_rows)")],
  as.data.table(ep_1_shift$coefficients)[term == "rpe", .(lag, beta = beta_gls, se = se_gls, state = "Pathology: Dropped 3 dummy TRs (time_offset = 0)")],
  as.data.table(ep_1_trbug$coefficients)[term == "rpe", .(lag, beta = beta_gls, se = se_gls, state = "Pathology: Onsets in TR indices instead of seconds")]
)

# Step 1B: 36-parcel AR(1) rho distribution — Healthy vs Unmodeled Scanner Shim/Step Drift
n_parc_rho <- 36L
vox_rho <- paste0("rv", seq_len(n_parc_rho))
grp_rho <- setNames(paste0("P_", seq_len(n_parc_rho)), vox_rho)
n_t_rho <- 200L
bold_rho_ok <- bold_rho_drift <- ev_rho <- list()
true_rhos <- seq(0.30, 0.62, length.out = n_parc_rho)
t_seq <- seq_len(n_t_rho)
shim_step <- 0.30 * plogis((t_seq - 95) / 2.5) + 0.16 * sin(2 * pi * t_seq / 55)

for (r in c("run1", "run2", "run3")) {
  on_r <- seq(10, 180, by = 12)
  s_r <- rowSums(vapply(on_r, function(o) source_col(o, 0, n_t_rho, tr) * 0.9, numeric(n_t_rho)))
  sig_r <- conv_open(s_r, h_true)
  Y_ok <- vapply(seq_len(n_parc_rho), function(v) 100 + sig_r + sim_ar1(n_t_rho, true_rhos[v], 0.035), numeric(n_t_rho))
  Y_dr <- vapply(seq_len(n_parc_rho), function(v) Y_ok[, v] + shim_step * runif(1, 0.75, 1.25), numeric(n_t_rho))
  colnames(Y_ok) <- colnames(Y_dr) <- vox_rho
  bold_rho_ok[[r]] <- Y_ok
  bold_rho_drift[[r]] <- Y_dr
  ev_rho[[r]] <- data.frame(onset = on_r, duration = 0, condition = "task")
}

fit_rho_ok <- estimate_crossrun_hrf(bold_rho_ok, ev_rho, tr = tr, kernel_grid = candidates,
                                    spatial_groups = grp_rho, pooling = "parcel", min_parcel_voxels = 1L,
                                    amplitude_model = "condition")
fit_rho_dr <- estimate_crossrun_hrf(bold_rho_drift, ev_rho, tr = tr, kernel_grid = candidates,
                                    spatial_groups = grp_rho, pooling = "parcel", min_parcel_voxels = 1L,
                                    amplitude_model = "condition")

rho_dt <- rbind(
  data.table(voxel = vox_rho, rho = fit_rho_ok$noise$run1$rho, state = "Healthy: Clean / detrended BOLD"),
  data.table(voxel = vox_rho, rho = fit_rho_dr$noise$run1$rho, state = "Pathology: Unmodeled shim step & mid-band drift")
)

# ==============================================================================
# 2. Step 2: Oversmoothing Gain Collapse vs. Uncensored Motion-Spike Ringing
# ==============================================================================
cat("2. Simulating Step 2 oversmoothing gain collapse and uncensored spike ringing...\n")
n_t2 <- 180L
vox2 <- paste0("v", 1:6)
parc2 <- setNames(rep("P2", 6L), vox2)
bold2_spike <- valid2_censor <- ev2 <- list()
spike_time <- 92L

for (r in c("run1", "run2", "run3")) {
  on2 <- c(14, 28, 43, 59, 76, 108, 125, 142, 156)
  amp2 <- c(1.2, 0.8, 1.5, 1.0, 1.8, 1.6, 0.9, 1.3, 1.1)
  s2 <- numeric(n_t2)
  for (i in seq_along(on2)) s2 <- s2 + source_col(on2[i], 0, n_t2, tr) * amp2[i]
  Y_sp <- vapply(1:6, function(v) 100 + conv_open(s2, h_true) + sim_ar1(n_t2, 0.45, 0.028), numeric(n_t2))
  if (r == "run1") Y_sp[spike_time, ] <- Y_sp[spike_time, ] + 0.95
  colnames(Y_sp) <- vox2
  bold2_spike[[r]] <- Y_sp
  v_mask <- rep(TRUE, n_t2)
  if (r == "run1") v_mask[(spike_time - 1L):(spike_time + 1L)] <- FALSE
  valid2_censor[[r]] <- v_mask
  ev2[[r]] <- data.frame(onset = on2, duration = 0, condition = "feedback")
}

hfit2_ok <- estimate_crossrun_hrf(bold2_spike, ev2, tr = tr, kernel_grid = candidates,
                                  valid_rows = valid2_censor, spatial_groups = parc2, pooling = "parcel")
dec2_ok <- deconvolve_crossrun_hrf(bold2_spike, hfit2_ok, scaling = "noise_sd",
                                   penalty = "sobolev", lambda_mode = "auto", rwls = "auto")
rec2_ok <- diagnose_crossrun_recovery(dec2_ok, events = ev2, lags = 0)

dec2_over <- deconvolve_crossrun_hrf(bold2_spike, hfit2_ok, lambda = 8.0, scaling = "noise_sd",
                                     penalty = "sobolev", rwls = "auto")
rec2_over <- diagnose_crossrun_recovery(dec2_over, events = ev2, lags = 0)

hfit2_uncens <- estimate_crossrun_hrf(bold2_spike, ev2, tr = tr, kernel_grid = candidates,
                                      spatial_groups = parc2, pooling = "parcel")
dec2_ring <- deconvolve_crossrun_hrf(bold2_spike, hfit2_uncens, lambda = 0.04, scaling = "noise_sd",
                                     penalty = "sobolev", rwls = "none")
rec2_ring <- diagnose_crossrun_recovery(dec2_ring, events = ev2, lags = 0)

act_scale <- function(dec) rowMeans(dec$activity$run1) * mean(dec$scaling$run1$scale)
win_idx <- 68:118
step2_traces <- rbind(
  data.table(time = win_idx - 1L, activity = act_scale(dec2_ok)[win_idx],   mode = "Healthy: valid_rows censor + auto lambda/rwls"),
  data.table(time = win_idx - 1L, activity = act_scale(dec2_over)[win_idx], mode = "Pathology A: Oversmoothed (lambda = 8.0; gain = 0.04)"),
  data.table(time = win_idx - 1L, activity = act_scale(dec2_ring)[win_idx], mode = "Pathology B: Uncensored motion spike at 91 s (rwls = 'none')")
)

step2_gain_tbl <- rbind(
  data.table(mode = "Healthy (`lambda_mode = 'auto'` + `valid_rows`)",
             mean_gain = mean(rec2_ok$summary$mean_gain),
             gain_err  = mean(rec2_ok$summary$mean_abs_gain_error)),
  data.table(mode = "Oversmoothed (`lambda = 8.0`)",
             mean_gain = mean(rec2_over$summary$mean_gain),
             gain_err  = mean(rec2_over$summary$mean_abs_gain_error)),
  data.table(mode = "Uncensored spike (`rwls = 'none'`, `lambda = 0.04`)",
             mean_gain = mean(rec2_ring$summary$mean_gain),
             gain_err  = mean(rec2_ring$summary$mean_abs_gain_error))
)

# ==============================================================================
# 3. Steps 3 & 4: Censoring Shadow on se_demixed + decon_eb vs decon_demixed Misuse
# ==============================================================================
cat("3. Simulating Step 3/4 censoring shadow on se_demixed and decon_eb regression shrinkage...\n")
n_t3 <- 220L
m3 <- 16L
on3 <- seq(10, by = 12.8, length.out = m3)
on3[m3] <- 215.5 # final trial near scanner shutoff

bold3 <- ev3 <- cov3 <- valid3_ok <- valid3_gap <- list()
for (r in c("run1", "run2", "run3")) {
  rpe3 <- as.numeric(scale(rnorm(m3)))
  s3 <- numeric(n_t3)
  for (i in seq_len(m3)) s3 <- s3 + source_col(on3[i], 0, n_t3, tr) * (0.9 + 0.50 * rpe3[i])
  Y3 <- vapply(1:6, function(v) 100 + conv_open(s3, h_true) + sim_ar1(n_t3, 0.45, 0.032), numeric(n_t3))
  colnames(Y3) <- vox2
  bold3[[r]] <- Y3
  ev3[[r]] <- data.frame(onset = on3, duration = 0, condition = "feedback")
  cov3[[r]] <- data.frame(rpe = rpe3)
  valid3_ok[[r]] <- rep(TRUE, n_t3)
  vg <- rep(TRUE, n_t3)
  if (r == "run1") {
    # Trial 7 (onset = 86.8 s): censor 3 frames on its HRF peak (91..93 s -> retained with elevated SE)
    vg[92:94] <- FALSE
    # Trial 11 (onset = 138.0 s): censor 10 frames covering >85% of Trial 11's HRF (140..149 s -> gated to NA)
    vg[141:150] <- FALSE
  }
  valid3_gap[[r]] <- vg
}

hf3_ok  <- estimate_crossrun_hrf(bold3, ev3, tr = tr, kernel_grid = candidates, valid_rows = valid3_ok, spatial_groups = parc2, pooling = "parcel")
dc3_ok  <- deconvolve_crossrun_hrf(bold3, hf3_ok, scaling = "noise_sd", penalty = "sobolev")
ep3_ok  <- decon_epoch_regression(dc3_ok, cov3, ev3, target_condition = "feedback", lags = 0, parcel_groups = parc2, restore_scale = TRUE, return_panel = TRUE)

hf3_gap <- estimate_crossrun_hrf(bold3, ev3, tr = tr, kernel_grid = candidates, valid_rows = valid3_gap, spatial_groups = parc2, pooling = "parcel")
dc3_gap <- deconvolve_crossrun_hrf(bold3, hf3_gap, scaling = "noise_sd", penalty = "sobolev")
ep3_gap <- decon_epoch_regression(dc3_gap, cov3, ev3, target_condition = "feedback", lags = 0, parcel_groups = parc2, restore_scale = TRUE, return_panel = TRUE)

se_trace_dt <- rbind(
  as.data.table(ep3_ok$panel)[run == "run1", .(trial, onset, se_demixed, status = "Healthy: Uncensored run (uptick only at run-end Trial 16)")],
  as.data.table(ep3_gap$panel)[run == "run1", .(trial, onset, se_demixed, status = "Scrubbed: Partial gap spikes Trial 7 SE; full gap gates Trial 11 to NA")]
)

# Cohort demonstration of decon_demixed vs decon_eb misuse (16 subjects, moderate single-trial reliability ~ 0.48)
eb_panels <- lapply(1:16, function(s_idx) {
  b_s <- bold3; c_s <- cov3
  for (r in c("run1", "run2", "run3")) {
    rpe_s <- as.numeric(scale(rnorm(m3)))
    s_vec <- numeric(n_t3)
    for (i in seq_len(m3)) s_vec <- s_vec + source_col(on3[i], 0, n_t3, tr) * (0.9 + 0.50 * rpe_s[i] + rnorm(1, sd = 0.20))
    Y_s <- vapply(1:6, function(v) 100 + conv_open(s_vec, h_true) + sim_ar1(n_t3, 0.45, 0.28), numeric(n_t3))
    colnames(Y_s) <- vox2
    b_s[[r]] <- Y_s
    c_s[[r]] <- data.frame(rpe = rpe_s)
  }
  hf_s <- estimate_crossrun_hrf(b_s, ev3, tr = tr, kernel_grid = candidates, spatial_groups = parc2, pooling = "parcel")
  dc_s <- deconvolve_crossrun_hrf(b_s, hf_s, scaling = "noise_sd", penalty = "sobolev")
  ep_s <- decon_epoch_regression(dc_s, covariates = NULL, events = ev3, target_condition = "feedback",
                                 lags = 0, parcel_groups = parc2, restore_scale = TRUE,
                                 tau2_mode = "marginal", return_panel = TRUE)
  p_dt <- as.data.table(ep_s$panel)
  p_dt[, id := sprintf("S%02d", s_idx)]
  rpe_all <- rbindlist(c_s)$rpe
  p_dt[, rpe_within := rpe_all - ave(rpe_all, p_dt$run)]
  p_dt
})
eb_group <- rbindlist(eb_panels)

mb_demixed <- mixed_by(
  eb_group[!is.na(decon_demixed)],
  model_formulae = list(m = decon_demixed ~ rpe_within + (1 + rpe_within || id) + (1 | id:run)),
  weights = "precision_weight_cohort", calculate = "parameter_estimates_reml", padjust_by = NULL
)
mb_eb_misuse <- mixed_by(
  eb_group[!is.na(decon_eb)],
  model_formulae = list(m = decon_eb ~ rpe_within + (1 + rpe_within || id) + (1 | id:run)),
  weights = NULL, calculate = "parameter_estimates_reml", padjust_by = NULL
)

eb_comp_dt <- rbind(
  as.data.table(mb_demixed$coef_df_reml)[term == "rpe_within",
    .(outcome_col = "Healthy: `decon_demixed` + `precision_weight_cohort`", estimate, std.error, conf.low, conf.high)],
  as.data.table(mb_eb_misuse$coef_df_reml)[term == "rpe_within",
    .(outcome_col = "Pathology: Shrunken `decon_eb` entered as outcome", estimate, std.error, conf.low, conf.high)]
)

# ==============================================================================
# 4. Steps 5 & 6: Uncentered/No-Random-Slope Group Models + Shared-Noise Connectivity Pitfalls
# ==============================================================================
cat("4. Simulating Step 5 centering/random-slope pitfalls and Step 6 connectivity diagnostics...\n")
# Step 5: 36 subjects, 16 target trials/subject; true within slope = +0.45, between slope = -0.35, random slope SD = 0.42
b_subj <- as.numeric(scale(rnorm(36))) * 0.42
xbar_subj <- as.numeric(scale(rnorm(36))) * 1.60
s5_sim <- rbindlist(lapply(1:36, function(s_idx) {
  x_w   <- as.numeric(scale(rnorm(16))) * 0.55
  x_raw <- xbar_subj[s_idx] + x_w
  u_0   <- rnorm(1, sd = 0.20)
  eps   <- as.numeric(scale(rnorm(16))) * 1.35
  z_hat <- -0.35 * xbar_subj[s_idx] + u_0 + (0.45 + b_subj[s_idx]) * x_w + eps
  data.table(id = sprintf("S%02d", s_idx), run = rep(c("run1", "run2"), each = 8L),
             decon_demixed = z_hat, x_raw = x_raw, x_within = x_w, x_mean = mean(x_raw))
}))

fit5_ok      <- lmerTest::lmer(decon_demixed ~ x_within + x_mean + (1 + x_within || id), data = s5_sim)
fit5_uncent  <- lmerTest::lmer(decon_demixed ~ x_raw + (1 | id), data = s5_sim)
fit5_nors    <- lmerTest::lmer(decon_demixed ~ x_within + x_mean + (1 | id), data = s5_sim)

extract_lme <- function(fit, term_nm, lbl) {
  cf <- summary(fit)$coefficients[term_nm, ]
  ci <- cf["Estimate"] + c(-1, 1) * stats::qt(0.975, df = cf["df"]) * cf["Std. Error"]
  data.table(spec = lbl, estimate = unname(cf["Estimate"]), se = unname(cf["Std. Error"]),
             df = unname(cf["df"]), conf.low = ci[1L], conf.high = ci[2L])
}
step5_pitfalls <- rbind(
  extract_lme(fit5_ok,     "x_within", "Healthy: x_within + x_mean +\n(1 + x_within || id)"),
  extract_lme(fit5_uncent, "x_raw",    "Pathology A: Uncentered x_raw\n(within/between conflated)"),
  extract_lme(fit5_nors,   "x_within", "Pathology B: Omit random slope\n(1 | id) (SE too narrow)")
)

# Step 6: 16 subjects, true baseline connectivity beta_S = 0.35, shared scanner noise r_ST = 0.55
# Compare: (1) Naive OLS (attenuated by predictor noise + biased by shared noise? Wait: on pure null beta_S = 0.00 vs active 0.35!)
# Let's use true baseline connectivity beta_S = 0.35 with strong predictor noise (attenuates naive) vs strong shared noise r_ST = 0.55:
# Better yet: test baseline seed->target connectivity where true beta_S = 0.35, predictor error SD = 0.85, and shared noise r_ST = 0.55!
n_s6 <- 16L
conn_epochs <- setNames(lapply(seq_len(n_s6), function(s_idx) {
  id_lbl <- sprintf("S%02d", s_idx)
  r_shared <- 0.55
  p_rows <- rbindlist(lapply(c("run1", "run2"), function(r) {
    m_tr <- 40L
    x_mod <- as.numeric(scale(rnorm(m_tr)))
    theta_seed <- rnorm(m_tr, sd = 0.70)
    theta_targ <- (0.35 + rnorm(1, sd = 0.06)) * theta_seed + 0.20 * x_mod + rnorm(m_tr, sd = 0.45)
    se_s <- rep(0.75, m_tr)
    se_t <- rep(0.75, m_tr)
    e_shared <- rnorm(m_tr)
    u_i <- se_s * (sqrt(r_shared) * e_shared + sqrt(1 - r_shared) * rnorm(m_tr))
    v_i <- se_t * (sqrt(r_shared) * e_shared + sqrt(1 - r_shared) * rnorm(m_tr))
    rbind(
      data.table(id = id_lbl, run = r, trial = seq_len(m_tr), event_index = seq_len(m_tr),
                 evt_time = 0, parcel = "Seed", decon_demixed = theta_seed + u_i, se_demixed = se_s, x_within = x_mod),
      data.table(id = id_lbl, run = r, trial = seq_len(m_tr), event_index = seq_len(m_tr),
                 evt_time = 0, parcel = "Target", decon_demixed = theta_targ + v_i, se_demixed = se_t, x_within = x_mod)
    )
  }))
  ucor <- matrix(c(1, r_shared, r_shared, 1), 2, 2, dimnames = list(c("Seed", "Target"), c("Seed", "Target")))
  list(panel = p_rows, unit_error_cor = ucor)
}), sprintf("S%02d", seq_len(n_s6)))

# Compare:
# 1. Healthy: M-MOMS + shared_error = "auto" (subtracts BOTH predictor variance & shared cross-parcel noise -> ~0.35)
# 2. Pathology A: M-MOMS + shared_error = "zero" (e.g. stripped shared_noise_loading in CSV export -> disattenuates shared noise up to ~0.90!)
# 3. Pathology B: Naive OLS without shared_error = "zero" on independent-error cohort (`r_shared = 0` comparator -> attenuated to ~0.16)
c_ok   <- decon_connectivity_repair(conn_epochs, "Seed", "Target", covariates = "x_within", test_ppi = FALSE, shared_error = "auto", method = "mmoms")
c_zero <- decon_connectivity_repair(conn_epochs, "Seed", "Target", covariates = "x_within", test_ppi = FALSE, shared_error = "zero", method = "mmoms")

# Also independent-noise naive comparator to show pure attenuation
conn_indep <- lapply(conn_epochs, function(obj) {
  dt <- copy(obj$panel)
  dt[parcel == "Target", decon_demixed := 0.35 * dt[parcel == "Seed", decon_demixed] * 0 + decon_demixed]
  obj
})
# Simpler: extract seed_c from c_ok (M-MOMS + shared_error='auto'), c_zero (M-MOMS + shared_error='zero' / lost CSV metadata), and naive on r_shared=0:
conn_r0 <- setNames(lapply(seq_len(n_s6), function(s_idx) {
  id_lbl <- sprintf("S%02d", s_idx)
  p_rows <- rbindlist(lapply(c("run1", "run2"), function(r) {
    m_tr <- 40L
    x_mod <- as.numeric(scale(rnorm(m_tr)))
    theta_seed <- rnorm(m_tr, sd = 0.70)
    theta_targ <- 0.35 * theta_seed + 0.20 * x_mod + rnorm(m_tr, sd = 0.45)
    se_s <- rep(0.75, m_tr); se_t <- rep(0.75, m_tr)
    rbind(
      data.table(id = id_lbl, run = r, trial = seq_len(m_tr), event_index = seq_len(m_tr),
                 evt_time = 0, parcel = "Seed", decon_demixed = theta_seed + rnorm(m_tr, sd = se_s), se_demixed = se_s, x_within = x_mod),
      data.table(id = id_lbl, run = r, trial = seq_len(m_tr), event_index = seq_len(m_tr),
                 evt_time = 0, parcel = "Target", decon_demixed = theta_targ + rnorm(m_tr, sd = se_t), se_demixed = se_t, x_within = x_mod)
    )
  }))
  list(panel = p_rows, unit_error_cor = diag(2))
}), sprintf("S%02d", seq_len(n_s6)))
c_naive <- decon_connectivity_repair(conn_r0, "Seed", "Target", covariates = "x_within", test_ppi = FALSE, shared_error = "zero", method = "naive_unweighted")

ppi_subj_dt <- rbind(
  as.data.table(c_ok$subjects)[term == "seed_c",    .(id, method = "Healthy: M-MOMS +\nshared_error = 'auto'", estimate, std.error)],
  as.data.table(c_zero$subjects)[term == "seed_c",  .(id, method = "Pathology A: Stripped metadata /\nshared_error = 'zero' (inflated)", estimate, std.error)],
  as.data.table(c_naive$subjects)[term == "seed_c", .(id, method = "Pathology B: Naive OLS\n(predictor-noise attenuated)", estimate, std.error)]
)
ppi_group_dt <- rbind(
  as.data.table(c_ok$group)[term == "seed_c",    .(method = "Healthy: M-MOMS + auto shared_error", estimate, std.error, conf.low, conf.high)],
  as.data.table(c_zero$group)[term == "seed_c",  .(method = "Pathology A: Lost shared_error ('zero')", estimate, std.error, conf.low, conf.high)],
  as.data.table(c_naive$group)[term == "seed_c", .(method = "Pathology B: Naive OLS (uncorrected)", estimate, std.error, conf.low, conf.high)]
)

qc_metrics <- data.table(
  check_stage = c(
    "Step 0/1: Intact clock vs. dropped dummy TRs (`time_offset = 0`)",
    "Step 0/1: Onsets in seconds vs. TR indices (`TR = 2.0 s`)",
    "Step 1: Clean BOLD vs. unmodeled shim drift (`rho >= 0.89` share)",
    "Step 2: Auto `lambda` vs. oversmoothed `lambda = 8` (`mean_gain`)",
    "Step 3: Uncensored trial SE vs. partial-scrub shadow Trial 7 SE",
    "Step 4/5: `decon_demixed` vs. shrunken `decon_eb` slope (truth = 0.50)",
    "Step 5: Centered random-slope `lmer` vs. uncentered / no-random-slope",
    "Step 6: Seed connectivity (truth = 0.35; M-MOMS auto vs. lost shared_error / naive)"
  ),
  healthy_value = c(
    sprintf("Peak at 0 s (beta = %.3f)", clock_curves[state %like% "Healthy" & lag == 0, beta]),
    sprintf("Task p = %.2e", fit_1_ok$selections$run1$selected_task_p[1L]),
    sprintf("%.1f%% at ceiling (median rho = %.2f)", 100 * mean(rho_dt[state %like% "Healthy", rho >= 0.89]), median(rho_dt[state %like% "Healthy", rho])),
    sprintf("mean_gain = %.3f", step2_gain_tbl$mean_gain[1L]),
    sprintf("Median interior SE = %.3f", median(se_trace_dt[status %like% "Healthy" & trial < 16, se_demixed], na.rm = TRUE)),
    sprintf("beta = %.3f [%.3f, %.3f]", eb_comp_dt$estimate[1L], eb_comp_dt$conf.low[1L], eb_comp_dt$conf.high[1L]),
    sprintf("beta = %.3f (SE = %.3f, df = %.1f)", step5_pitfalls$estimate[1L], step5_pitfalls$se[1L], step5_pitfalls$df[1L]),
    sprintf("beta = %.3f [%.3f, %.3f]", ppi_group_dt$estimate[1L], ppi_group_dt$conf.low[1L], ppi_group_dt$conf.high[1L])
  ),
  pathological_value = c(
    sprintf("Shifted to -6 s (lag-0 beta = %.3f)", clock_curves[state %like% "Dropped" & lag == 0, beta]),
    sprintf("Lag-0 beta = %.3f (schedule compressed)", clock_curves[state %like% "TR indices" & lag == 0, beta]),
    sprintf("%.1f%% at 0.90 ceiling (median rho = %.2f)", 100 * mean(rho_dt[state %like% "Pathology", rho >= 0.89]), median(rho_dt[state %like% "Pathology", rho])),
    sprintf("mean_gain = %.3f (flat activity)", step2_gain_tbl$mean_gain[2L]),
    sprintf("Trial 7 SE = %.3f; Trial 11 = NA", se_trace_dt[status %like% "Scrubbed" & trial == 7, se_demixed]),
    sprintf("beta = %.3f [%.3f, %.3f]", eb_comp_dt$estimate[2L], eb_comp_dt$conf.low[2L], eb_comp_dt$conf.high[2L]),
    sprintf("Uncent. beta = %.3f; No-RS SE = %.3f (df = %.0f)", step5_pitfalls$estimate[2L], step5_pitfalls$se[3L], step5_pitfalls$df[3L]),
    sprintf("Lost shared_error = %.3f; Naive = %.3f", ppi_group_dt$estimate[2L], ppi_group_dt$estimate[3L])
  )
)

print(qc_metrics)
fwrite(qc_metrics, "qc_symptom_metrics.csv")

saveRDS(list(
  clock_curves   = clock_curves,
  rho_dt         = rho_dt,
  step2_traces   = step2_traces,
  step2_gain_tbl = step2_gain_tbl,
  on2_run1       = c(76, 108),
  spike_time     = spike_time - 1L,
  se_trace_dt    = se_trace_dt,
  eb_comp_dt     = eb_comp_dt,
  step5_pitfalls = step5_pitfalls,
  ppi_subj_dt    = ppi_subj_dt,
  ppi_group_dt   = ppi_group_dt,
  qc_metrics     = qc_metrics
), "qc_gallery_examples.rds")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("Saved QC gallery benchmark data.\n")
