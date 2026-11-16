# Step 0 simulations:
# 1. Experimental design & inter-phase jitter identifiability benchmarks
# 2. Live on-disk BIDS + fMRIPrep ingestion via read_fmriprep_subject(),
#    dummy-scan clock alignment checks, missed/timeout trial handling, and
#    non-aggressive ICA-AROMA + matched nuisance projection (fMRIPost-AROMA).

library(medusa)
library(data.table)
library(RNifti)

set.seed(20261100)

tr <- 1.0
candidates <- hrf_library(tr = tr, preset = "medusa")
pars <- attr(candidates, "parameters")
idx_true_hrf <- which.min(abs(pars$peak_delay - 4.5) + abs(pars$dispersion - 0.95) + ifelse(pars$undershoot_ratio == 6, 0, 0.5))
h_true <- candidates[, idx_true_hrf]
t_hrf <- (0:(nrow(candidates) - 1L)) * tr
sp_true <- stats::spline(t_hrf, h_true, n = 400L)
true_peak_sec <- sp_true$x[which.max(sp_true$y)]

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
# Part 1: Experimental Design Feasibility (ISI Jitter & Run Structure)
# ==============================================================================
cat("Part 1: Simulating experimental design & jitter benchmarks...\n")

jitter_regimes <- c(
  "Fixed 0.8 s (Unjittered rapid)",
  "Fixed 1.8 s (Unjittered short)",
  "Uniform jitter 1.0-4.5 s",
  "Exponential jitter 1.0-6.0 s"
)

sample_isi <- function(m, regime) {
  switch(
    regime,
    "Fixed 0.8 s (Unjittered rapid)" = rep(0.8, m),
    "Fixed 1.8 s (Unjittered short)" = rep(1.8, m),
    "Uniform jitter 1.0-4.5 s"       = runif(m, 1.0, 4.5),
    "Exponential jitter 1.0-6.0 s"   = pmin(6.0, 1.0 + rexp(m, rate = 1 / 1.6))
  )
}

n_rep_design <- 20L
vox_design <- c(paste0("Target_", 1:6), paste0("Cue_Only_", 1:6))
parc_design <- setNames(rep(c("Target_Parcel", "CueOnly_Parcel"), each = 6L), vox_design)

run_one_design_rep <- function(jitter_lbl, n_runs = 3L, trials_per_run = 20L, noise_sd = 0.035, noisy_first_run = FALSE) {
  bold <- events <- covariates <- list()
  rnames <- paste0("run", seq_len(n_runs))
  for (r in rnames) {
    sd_r <- if (noisy_first_run && r == "run1") noise_sd * 2.1 else noise_sd
    isi_cf <- sample_isi(trials_per_run, jitter_lbl)
    iti    <- runif(trials_per_run, 10.5, 14.5)
    cue_on <- numeric(trials_per_run)
    fb_on  <- numeric(trials_per_run)
    t_cur  <- 8.0
    for (i in seq_len(trials_per_run)) {
      cue_on[i] <- t_cur
      fb_on[i]  <- t_cur + isi_cf[i]
      t_cur     <- fb_on[i] + iti[i]
    }
    n_t <- as.integer(ceiling(max(fb_on) + 26))
    rpe <- as.numeric(scale(rnorm(trials_per_run)))
    # Cue covaries moderately with predictor (r ~ 0.50)
    cue_mod <- as.numeric(scale(0.50 * rpe + sqrt(1 - 0.50^2) * rnorm(trials_per_run)))

    ev <- rbind(
      data.frame(onset = cue_on, duration = 0, condition = "cue", trial_id = seq_len(trials_per_run)),
      data.frame(onset = fb_on,  duration = 0, condition = "feedback", trial_id = seq_len(trials_per_run))
    )
    ev <- ev[order(ev$onset), ]
    rownames(ev) <- NULL

    s_target <- s_cueonly <- numeric(n_t)
    for (i in seq_len(nrow(ev))) {
      tid <- ev$trial_id[i]
      sc <- source_col(ev$onset[i], 0, n_t, tr)
      if (ev$condition[i] == "cue") {
        s_target  <- s_target  + sc * (1.2 + 0.45 * cue_mod[tid])
        s_cueonly <- s_cueonly + sc * (1.4 + 0.55 * cue_mod[tid])
      } else {
        s_target  <- s_target  + sc * (0.9 + 0.60 * rpe[tid])
        # s_cueonly has ZERO true feedback response
      }
    }

    Y <- cbind(
      vapply(1:6, function(j) 100 + conv_open(s_target, h_true)  + sim_ar1(n_t, 0.45, sd_r), numeric(n_t)),
      vapply(1:6, function(j) 100 + conv_open(s_cueonly, h_true) + sim_ar1(n_t, 0.45, sd_r), numeric(n_t))
    )
    colnames(Y) <- vox_design
    bold[[r]] <- Y
    events[[r]] <- ev
    fb_tid <- ev$trial_id[ev$condition == "feedback"]
    covariates[[r]] <- data.frame(rpe = rpe[fb_tid])
  }

  hfit <- estimate_crossrun_hrf(
    bold = bold, events = events, tr = tr, kernel_grid = candidates,
    spatial_groups = parc_design, pooling = "parcel",
    amplitude_model = "condition", noise_model = "ar1",
    evidence = "calibrated", selection_mode = "soft_average"
  )
  est_peaks <- vapply(rnames, function(r) {
    k_wave <- hfit$selected_kernels[[r]][, 1L]
    sp <- stats::spline(t_hrf, k_wave, n = 400L)
    sp$x[which.max(sp$y)]
  }, numeric(1))
  hrf_mae <- mean(abs(est_peaks - true_peak_sec))

  dec <- deconvolve_crossrun_hrf(
    bold = bold, hrf_fit = hfit, scaling = "noise_sd",
    penalty = "sobolev", lambda_mode = "global",
    nuisance_mode = "ignore", observation_mask = "fit"
  )
  rec <- diagnose_crossrun_recovery(dec, events = events, lags = 0)
  mean_gain <- mean(rec$summary$mean_gain, na.rm = TRUE)

  ep <- decon_epoch_regression(
    deconvolution = dec, covariates = covariates, events = events,
    target_condition = "feedback", lags = 0, parcel_groups = parc_design,
    demix_other_events = TRUE, demix_target_trials = TRUE,
    demix_ridge = 0.02, restore_scale = TRUE, return_panel = FALSE
  )
  cf <- as.data.table(ep$coefficients)
  b_targ  <- cf[parcel == "Target_Parcel" & term == "rpe", beta_gls]
  se_targ <- cf[parcel == "Target_Parcel" & term == "rpe", se_gls]
  b_cue   <- cf[parcel == "CueOnly_Parcel" & term == "rpe", beta_gls]

  data.table(
    hrf_mae_sec   = hrf_mae,
    mean_gain     = mean_gain,
    beta_target   = b_targ,
    se_target     = se_targ,
    beta_cuebleed = b_cue
  )
}

jitter_summary <- rbindlist(lapply(jitter_regimes, function(jr) {
  reps <- rbindlist(lapply(seq_len(n_rep_design), function(i) run_one_design_rep(jr, n_runs = 3L, trials_per_run = 20L, noise_sd = 0.032)))
  data.table(
    dimension     = "Cue-Feedback ISI Schedule (3 runs x 20 trials)",
    setting       = jr,
    hrf_mae_sec   = mean(reps$hrf_mae_sec),
    mean_gain     = mean(reps$mean_gain),
    target_mean   = mean(reps$beta_target),
    target_rmse   = sqrt(mean((reps$beta_target - 0.60)^2)),
    mean_slope_se = mean(reps$se_target),
    cue_bleed_abs = mean(abs(reps$beta_cuebleed))
  )
}))

run_configs <- list(
  list(lbl = "1 run x 60 trials (Sub-block split)", n_runs = 1L, tpr = 60L),
  list(lbl = "2 runs x 30 trials (1 training run/fold)", n_runs = 2L, tpr = 30L),
  list(lbl = "3 runs x 20 trials (Recommended >=3 runs)", n_runs = 3L, tpr = 20L),
  list(lbl = "4 runs x 15 trials (4-fold cross-fit)", n_runs = 4L, tpr = 15L)
)

run_summary <- rbindlist(lapply(run_configs, function(cfg) {
  reps <- rbindlist(lapply(seq_len(n_rep_design), function(i) {
    run_one_design_rep("Uniform jitter 1.0-4.5 s", n_runs = cfg$n_runs, trials_per_run = cfg$tpr, noise_sd = 0.055, noisy_first_run = TRUE)
  }))
  data.table(
    dimension     = "Run Structure (60 total trials, moderate noise)",
    setting       = cfg$lbl,
    hrf_mae_sec   = mean(reps$hrf_mae_sec),
    mean_gain     = mean(reps$mean_gain),
    target_mean   = mean(reps$beta_target),
    target_rmse   = sqrt(mean((reps$beta_target - 0.60)^2)),
    mean_slope_se = mean(reps$se_target),
    cue_bleed_abs = mean(abs(reps$beta_cuebleed))
  )
}))

design_benchmarks <- rbind(jitter_summary, run_summary)
print(design_benchmarks)
fwrite(design_benchmarks, "design_benchmarks.csv")

# ==============================================================================
# Part 2: Live BIDS + fMRIPrep Subject Ingestion, Clock Alignment & Missed Trials
# ==============================================================================
cat("\nPart 2: Constructing mini BIDS + fMRIPrep dataset on disk and testing read_fmriprep_subject()...\n")

bids_root <- tempfile("medusa_step0_bids_")
fprep_dir <- file.path(bids_root, "derivatives", "fmriprep")
bids_dir  <- file.path(bids_root, "bids")
fprep_func <- file.path(fprep_dir, "sub-01", "func")
bids_func  <- file.path(bids_dir, "sub-01", "func")
dir.create(fprep_func, recursive = TRUE)
dir.create(bids_func, recursive = TRUE)

dims3 <- c(3L, 2L, 2L)
mask_arr <- array(1L, dim = dims3)
atlas_arr <- array(c(rep(1L, 6L), rep(2L, 6L)), dim = dims3)
mask_nii <- RNifti::asNifti(mask_arr); RNifti::pixdim(mask_nii) <- c(2, 2, 2)
atlas_path <- file.path(bids_root, "toy_atlas_2parcels.nii.gz")
atlas_nii <- RNifti::asNifti(atlas_arr); RNifti::pixdim(atlas_nii) <- c(2, 2, 2)
RNifti::writeNifti(atlas_nii, atlas_path)

n_vols <- 240L
n_dummy <- 4L # first 4 TRs (0..3 s) are non-steady-state dummy volumes

rl_model_rows <- list()

for (r in 1:3) {
  stem <- sprintf("sub-01_task-learn_run-0%d", r)
  # 18 total trials per run: trials 6, 12, 16 are missed/timeout trials
  n_tr_run <- 18L
  missed_idx <- c(6L, 12L, 16L)
  cue_on <- seq(10.0, by = 11.8, length.out = n_tr_run)
  isi_cf <- runif(n_tr_run, 1.5, 3.6)
  fb_on  <- cue_on + isi_cf

  rpe_valid <- as.numeric(scale(rnorm(n_tr_run - length(missed_idx))))
  rpe_all <- rep(NA_real_, n_tr_run)
  rpe_all[-missed_idx] <- rpe_valid
  rt_all <- ifelse(seq_len(n_tr_run) %in% missed_idx, NA_real_, round(runif(n_tr_run, 0.45, 1.15), 3))

  ev_cue <- data.frame(
    onset         = round(cue_on, 3),
    duration      = rep(0.5, n_tr_run),
    trial_type    = "cue",
    trial_number  = seq_len(n_tr_run),
    response_time = rep("n/a", n_tr_run),
    rpe           = rep("n/a", n_tr_run),
    stringsAsFactors = FALSE
  )
  ev_fb <- data.frame(
    onset         = round(fb_on, 3),
    duration      = rep("n/a", n_tr_run), # impulse coded as n/a in BIDS
    trial_type    = ifelse(seq_len(n_tr_run) %in% missed_idx, "missed_feedback", "feedback"),
    trial_number  = seq_len(n_tr_run),
    response_time = ifelse(is.na(rt_all), "n/a", sprintf("%.3f", rt_all)),
    rpe           = ifelse(is.na(rpe_all), "n/a", sprintf("%.4f", rpe_all)),
    stringsAsFactors = FALSE
  )
  ev_bids <- rbind(ev_cue, ev_fb)
  ev_bids <- ev_bids[order(as.numeric(ev_bids$onset)), ]
  rownames(ev_bids) <- NULL

  valid_fb_rows <- which(ev_bids$trial_type == "feedback")
  rl_model_rows[[r]] <- data.table(
    id           = "sub-01",
    run          = paste0("run", r),
    event_index  = valid_fb_rows,
    trial_number = ev_bids$trial_number[valid_fb_rows],
    rpe_rl       = as.numeric(ev_bids$rpe[valid_fb_rows])
  )

  # Simulate true neural activity across ALL events (missed/timeout trials evoke a strong 2.5-unit visual/salience burst)
  s_p1 <- s_p2 <- numeric(n_vols)
  for (i in seq_len(nrow(ev_bids))) {
    on_i <- as.numeric(ev_bids$onset[i])
    dur_i <- if (ev_bids$duration[i] == "n/a") 0 else as.numeric(ev_bids$duration[i])
    tt_i <- ev_bids$trial_type[i]
    sc <- source_col(on_i, dur_i, n_vols, tr)
    if (tt_i == "cue") {
      s_p1 <- s_p1 + 0.75 * sc
      s_p2 <- s_p2 + 0.60 * sc
    } else if (tt_i == "feedback") {
      rpe_i <- as.numeric(ev_bids$rpe[i])
      s_p1 <- s_p1 + (0.85 + 0.65 * rpe_i) * sc
      s_p2 <- s_p2 + 0.50 * sc
    } else if (tt_i == "missed_feedback") {
      s_p1 <- s_p1 + 2.60 * sc
      s_p2 <- s_p2 + 2.00 * sc
    }
  }

  trans_z <- cumsum(rnorm(n_vols, sd = 0.012))
  spike_trs <- c(52L, 162L)
  trans_z[spike_trs] <- trans_z[spike_trs] + c(0.68, -0.74)
  fd_vec <- c(NA_real_, abs(diff(trans_z)) + runif(n_vols - 1L, 0.03, 0.11))

  arr4d <- array(0, dim = c(dims3, n_vols))
  for (v in 1:12) {
    s_v <- if (v <= 6L) s_p1 else s_p2
    y_v <- 100 + conv_open(s_v, h_true) + 0.04 * scale(trans_z)[, 1L] + sim_ar1(n_vols, 0.45, 0.030)
    y_v[1:n_dummy] <- y_v[1:n_dummy] + c(2.4, 1.5, 0.8, 0.35)
    y_v[spike_trs] <- y_v[spike_trs] + c(0.85, -0.90)
    idx3 <- arrayInd(v, dims3)
    arr4d[idx3[1], idx3[2], idx3[3], ] <- y_v
  }

  bold_nii <- RNifti::asNifti(arr4d)
  RNifti::pixdim(bold_nii) <- c(2, 2, 2, tr)
  RNifti::writeNifti(bold_nii, file.path(fprep_func, paste0(stem, "_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz")))
  RNifti::writeNifti(mask_nii, file.path(fprep_func, paste0(stem, "_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz")))
  writeLines('{"RepetitionTime": 1.0, "TaskName": "learn"}', file.path(bids_func, paste0(stem, "_bold.json")))

  conf_df <- data.frame(
    trans_x = round(rnorm(n_vols, sd = 0.015), 5),
    trans_y = round(rnorm(n_vols, sd = 0.015), 5),
    trans_z = round(trans_z, 5),
    rot_x   = round(rnorm(n_vols, sd = 0.0015), 6),
    rot_y   = round(rnorm(n_vols, sd = 0.0015), 6),
    rot_z   = round(rnorm(n_vols, sd = 0.0015), 6),
    framewise_displacement = ifelse(is.na(fd_vec), "n/a", sprintf("%.4f", fd_vec)),
    non_steady_state_outlier00 = as.integer(seq_len(n_vols) == 1L),
    non_steady_state_outlier01 = as.integer(seq_len(n_vols) == 2L),
    non_steady_state_outlier02 = as.integer(seq_len(n_vols) == 3L),
    non_steady_state_outlier03 = as.integer(seq_len(n_vols) == 4L),
    stringsAsFactors = FALSE
  )
  utils::write.table(conf_df, file.path(fprep_func, paste0(stem, "_desc-confounds_timeseries.tsv")),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(ev_bids, file.path(bids_func, paste0(stem, "_events.tsv")),
                     sep = "\t", row.names = FALSE, quote = FALSE)
}

sub01 <- read_fmriprep_subject(
  fmriprep_dir             = fprep_dir,
  bids_dir                 = bids_dir,
  sub                      = "01",
  task                     = "learn",
  atlas                    = atlas_path,
  strategy                 = "6p",
  fd_threshold             = 0.5,
  include_non_steady_state = TRUE
)
qc_summary <- as.data.table(sub01$confounds$summary)
qc_summary[, file := basename(file)]
print(qc_summary)
fwrite(qc_summary, "confounds_qc_summary.csv")

# ==============================================================================
# Demonstration A: Clock Alignment (Intact valid_rows vs time_offset vs Broken Dummy Drop)
# ==============================================================================
lags_clock <- -5:4

fit_intact <- estimate_crossrun_hrf(
  sub01, kernel_grid = candidates, pooling = "parcel",
  amplitude_model = "trial", selection_mode = "soft_average"
)
dec_intact <- deconvolve_crossrun_hrf(
  sub01, fit_intact, scaling = "noise_sd", penalty = "sobolev",
  nuisance_mode = "joint", observation_mask = "fit"
)
ep_intact <- decon_epoch_regression(
  dec_intact, covariates = "rpe", target_condition = "feedback",
  lags = lags_clock, aggregate_by_parcel = TRUE,
  restore_scale = TRUE, return_panel = TRUE
)

bold_trunc  <- lapply(sub01$bold, function(Y) Y[(n_dummy + 1L):nrow(Y), , drop = FALSE])
nuis_trunc  <- lapply(sub01$nuisance, function(Z) Z[(n_dummy + 1L):nrow(Z), , drop = FALSE])
valid_trunc <- lapply(sub01$valid_rows, function(v) v[(n_dummy + 1L):length(v)])
fit_offset <- estimate_crossrun_hrf(
  bold = bold_trunc, events = sub01$events, tr = tr, time_offset = n_dummy * tr,
  kernel_grid = candidates, nuisance = nuis_trunc, valid_rows = valid_trunc,
  spatial_groups = sub01$parcel_ids, pooling = "parcel",
  amplitude_model = "trial", selection_mode = "soft_average"
)
dec_offset <- deconvolve_crossrun_hrf(
  bold_trunc, fit_offset, scaling = "noise_sd", penalty = "sobolev",
  nuisance_mode = "joint", observation_mask = "fit"
)
ep_offset <- decon_epoch_regression(
  dec_offset, covariates = "rpe", events = sub01$events, target_condition = "feedback",
  lags = lags_clock, parcel_groups = sub01$parcel_ids,
  restore_scale = TRUE, return_panel = FALSE
)

fit_broken <- estimate_crossrun_hrf(
  bold = bold_trunc, events = sub01$events, tr = tr, time_offset = 0,
  kernel_grid = candidates, nuisance = nuis_trunc, valid_rows = valid_trunc,
  spatial_groups = sub01$parcel_ids, pooling = "parcel",
  amplitude_model = "trial", selection_mode = "soft_average"
)
dec_broken <- deconvolve_crossrun_hrf(
  bold_trunc, fit_broken, scaling = "noise_sd", penalty = "sobolev",
  nuisance_mode = "joint", observation_mask = "fit"
)
ep_broken <- decon_epoch_regression(
  dec_broken, covariates = "rpe", events = sub01$events, target_condition = "feedback",
  lags = lags_clock, parcel_groups = sub01$parcel_ids,
  restore_scale = TRUE, return_panel = FALSE
)

peak_of <- function(hfit) {
  sp <- stats::spline(t_hrf, hfit$selected_kernels[[1L]][, 1L], n = 400L)
  sp$x[which.max(sp$y)]
}

cf_intact <- as.data.table(ep_intact$coefficients)
cf_offset <- as.data.table(ep_offset$coefficients)
cf_broken <- as.data.table(ep_broken$coefficients)
p1_lbl <- as.character(unique(cf_intact$parcel)[1L])

extract_clock_row <- function(lbl, hfit, cf_dt) {
  sub_cf <- cf_dt[as.character(parcel) == p1_lbl & term == "rpe"]
  data.table(
    handling          = lbl,
    hrf_peak_sec      = peak_of(hfit),
    lag0_beta         = sub_cf[lag == 0, beta_gls],
    lag0_se           = sub_cf[lag == 0, se_gls],
    apparent_peak_lag = sub_cf[which.max(beta_gls), lag]
  )
}

clock_tbl <- rbind(
  extract_clock_row("Intact clock + valid_rows mask (Recommended)", fit_intact, cf_intact),
  extract_clock_row("Deleted 4 dummy TRs + explicit time_offset = 4 s", fit_offset, cf_offset),
  extract_clock_row("Deleted 4 dummy TRs without time_offset (Clock shifted)", fit_broken, cf_broken)
)
print(clock_tbl)
fwrite(clock_tbl, "clock_alignment_metrics.csv")

# ==============================================================================
# Demonstration B: Missed/Timeout Trial Handling & RL Table Joins
# ==============================================================================
ev_dropped_missed <- subset(sub01$events, condition != "missed_feedback")
fit_drop_missed <- estimate_crossrun_hrf(
  bold = sub01$bold, events = ev_dropped_missed, tr = tr,
  kernel_grid = candidates, nuisance = sub01$confounds,
  spatial_groups = sub01$parcel_ids, pooling = "parcel",
  amplitude_model = "trial", selection_mode = "soft_average"
)
dec_drop_missed <- deconvolve_crossrun_hrf(
  sub01$bold, fit_drop_missed, scaling = "noise_sd", penalty = "sobolev",
  nuisance_mode = "joint", observation_mask = "fit"
)
ep_drop_missed <- decon_epoch_regression(
  dec_drop_missed, covariates = "rpe", events = ev_dropped_missed,
  target_condition = "feedback", lags = 0, parcel_groups = sub01$parcel_ids,
  restore_scale = TRUE, return_panel = TRUE
)
cf_drop <- as.data.table(ep_drop_missed$coefficients)

rl_dt <- rbindlist(rl_model_rows)
panel_lag0 <- copy(as.data.table(ep_intact$panel))[as.character(parcel) == p1_lbl & evt_time == 0]
joined_safe <- merge(panel_lag0, rl_dt, by = c("run", "event_index"))
stopifnot(nrow(joined_safe) == 45L, max(abs(joined_safe$rpe - joined_safe$rpe_rl)) < 1e-6)

# Fragile join: matching target-row counter `trial` (1..15) to task `trial_number` (1..18, skipping 6, 12, 16)
bad_join <- merge(panel_lag0, rl_dt, by.x = c("run", "trial"), by.y = c("run", "trial_number"))
fit_bad_join <- stats::lm(decon_demixed ~ rpe_rl, data = bad_join, weights = precision_weight)

missed_tbl <- rbind(
  data.table(
    strategy     = "Keep 'missed_feedback' in events + event_index join (Recommended)",
    est_noise_sd = mean(dec_intact$scaling[[1L]]$noise_sd[1:6]),
    est_ar1_rho  = mean(fit_intact$noise[[1L]]$rho[1:6]),
    lag0_beta    = cf_intact[as.character(parcel) == p1_lbl & lag == 0 & term == "rpe", beta_gls],
    lag0_se      = cf_intact[as.character(parcel) == p1_lbl & lag == 0 & term == "rpe", se_gls],
    matched_rows = nrow(joined_safe)
  ),
  data.table(
    strategy     = "Delete missed trials from events.tsv (Timeout BOLD unmodeled)",
    est_noise_sd = mean(dec_drop_missed$scaling[[1L]]$noise_sd[1:6]),
    est_ar1_rho  = mean(fit_drop_missed$noise[[1L]]$rho[1:6]),
    lag0_beta    = cf_drop[as.character(parcel) == p1_lbl & lag == 0 & term == "rpe", beta_gls],
    lag0_se      = cf_drop[as.character(parcel) == p1_lbl & lag == 0 & term == "rpe", se_gls],
    matched_rows = nrow(as.data.table(ep_drop_missed$panel)[as.character(parcel) == p1_lbl & evt_time == 0])
  ),
  data.table(
    strategy     = "Misaligned join on target counter `trial` == `trial_number`",
    est_noise_sd = mean(dec_intact$scaling[[1L]]$noise_sd[1:6]),
    est_ar1_rho  = mean(fit_intact$noise[[1L]]$rho[1:6]),
    lag0_beta    = unname(stats::coef(fit_bad_join)["rpe_rl"]),
    lag0_se      = summary(fit_bad_join)$coefficients["rpe_rl", "Std. Error"],
    matched_rows = nrow(bad_join)
  )
)
print(missed_tbl)
fwrite(missed_tbl, "missed_trial_metrics.csv")

ev_r1 <- as.data.table(sub01$events)[run == unique(sub01$events$run)[1L]]
run1_missed_time <- as.numeric(ev_r1[condition == "missed_feedback", onset])
run1_valid_time  <- as.numeric(ev_r1[condition == "feedback", onset])
act_intact_p1 <- rowMeans(dec_intact$activity[[1L]][, 1:6]) * mean(dec_intact$scaling[[1L]]$scale[1:6])
act_drop_p1   <- rowMeans(dec_drop_missed$activity[[1L]][, 1:6]) * mean(dec_drop_missed$scaling[[1L]]$scale[1:6])

# ==============================================================================
# Demonstration C: Non-Aggressive ICA-AROMA & Matched Nuisance Projection
# ==============================================================================
cat("\nPart 3: Simulating fMRIPost-AROMA non-aggressive denoising & matched confound projection...\n")
set.seed(20261116)

aroma_root      <- tempfile("medusa_step0_aroma_")
fprep_aroma_dir <- file.path(aroma_root, "derivatives", "fmriprep")
post_aroma_dir  <- file.path(aroma_root, "derivatives", "fmripost-aroma")
bids_aroma_dir  <- file.path(aroma_root, "bids")
dir.create(file.path(fprep_aroma_dir, "sub-01", "func"), recursive = TRUE)
dir.create(file.path(post_aroma_dir,  "sub-01", "func"), recursive = TRUE)
dir.create(file.path(bids_aroma_dir,  "sub-01", "func"), recursive = TRUE)

atlas_aroma_path <- file.path(aroma_root, "atlas_2parcels.nii.gz")
RNifti::writeNifti(atlas_nii, atlas_aroma_path)

n_vols_a <- 240L
n_tr_a   <- 20L
mot_ic_list <- list()

for (r in 1:3) {
  stem_a <- sprintf("sub-01_task-learn_run-0%d", r)
  cue_on <- seq(8.0, by = 10.8, length.out = n_tr_a)
  fb_on  <- cue_on + runif(n_tr_a, 1.4, 3.6)
  rpe_r  <- as.numeric(scale(rnorm(n_tr_a)))

  ev_a <- rbind(
    data.frame(onset = round(cue_on, 3), duration = 0.5, trial_type = "cue", rpe = "n/a", stringsAsFactors = FALSE),
    data.frame(onset = round(fb_on, 3),  duration = "n/a", trial_type = "feedback", rpe = sprintf("%.4f", rpe_r), stringsAsFactors = FALSE)
  )
  ev_a <- ev_a[order(as.numeric(ev_a$onset)), ]
  rownames(ev_a) <- NULL

  s_targ <- s_ctrl <- s_rpe_only <- numeric(n_vols_a)
  for (i in seq_len(nrow(ev_a))) {
    on_i  <- as.numeric(ev_a$onset[i])
    dur_i <- if (ev_a$duration[i] == "n/a") 0 else as.numeric(ev_a$duration[i])
    sc    <- source_col(on_i, dur_i, n_vols_a, tr)
    if (ev_a$trial_type[i] == "cue") {
      s_targ <- s_targ + 0.80 * sc
      s_ctrl <- s_ctrl + 0.90 * sc
    } else {
      rpe_i <- as.numeric(ev_a$rpe[i])
      s_targ <- s_targ + (0.90 + 0.65 * rpe_i) * sc
      s_ctrl <- s_ctrl + 0.55 * sc
      s_rpe_only <- s_rpe_only + (0.65 * rpe_i) * sc
    }
  }

  x_targ <- conv_open(s_targ, h_true)
  x_ctrl <- conv_open(s_ctrl, h_true)
  x_rpe  <- as.numeric(scale(conv_open(s_rpe_only, h_true)))

  # 2 accepted neural ICs + 6 rejected motion/edge ICs
  ic1_acc <- as.numeric(scale(x_targ + rnorm(n_vols_a, sd = 0.010)))
  ic2_acc <- as.numeric(scale(x_ctrl + rnorm(n_vols_a, sd = 0.010)))
  # IC_3: task-correlated head motion (feedback flinch / button press, r ~ 0.62 with convolved RPE)
  ic3_mot <- as.numeric(scale(0.62 * x_rpe + sqrt(1 - 0.62^2) * sim_ar1(n_vols_a, 0.30, 1.0)))
  ic4_mot <- as.numeric(scale(sim_ar1(n_vols_a, 0.35, 1.0)))
  ic5_mot <- as.numeric(scale(sim_ar1(n_vols_a, 0.25, 1.0)))
  ic6_mot <- as.numeric(scale(sim_ar1(n_vols_a, 0.40, 1.0)))
  ic7_mot <- as.numeric(scale(rnorm(n_vols_a)))
  ic8_mot <- as.numeric(scale(rnorm(n_vols_a)))
  mot_ic_list[[r]] <- cbind(ic3_mot, ic4_mot, ic5_mot, ic6_mot, ic7_mot, ic8_mot)

  # Physiological CSF / WM noise shared with gray matter AND contaminated by motion ICs in raw confounds
  u_csf <- sim_ar1(n_vols_a, 0.50, 0.080)
  u_wm  <- sim_ar1(n_vols_a, 0.55, 0.070)

  arr4d_a <- array(0, dim = c(dims3, n_vols_a))
  for (v in 1:12) {
    sig_v <- if (v <= 6L) x_targ else x_ctrl
    mot_v <- (if (v <= 6L) 0.085 else 0.060) * ic3_mot +
             0.110 * ic4_mot - 0.095 * ic5_mot + 0.080 * ic6_mot +
             0.055 * ic7_mot - 0.050 * ic8_mot
    phys_v <- 1.0 * u_csf + 0.9 * u_wm
    y_v <- 100 + sig_v + mot_v + phys_v + sim_ar1(n_vols_a, 0.45, 0.030)
    idx3 <- arrayInd(v, dims3)
    arr4d_a[idx3[1], idx3[2], idx3[3], ] <- y_v
  }

  bold_nii_a <- RNifti::asNifti(arr4d_a); RNifti::pixdim(bold_nii_a) <- c(2, 2, 2, tr)
  RNifti::writeNifti(bold_nii_a, file.path(fprep_aroma_dir, "sub-01", "func",
    paste0(stem_a, "_space-MNI152NLin6Asym_res-2_desc-preproc_bold.nii.gz")))
  RNifti::writeNifti(mask_nii, file.path(fprep_aroma_dir, "sub-01", "func",
    paste0(stem_a, "_space-MNI152NLin6Asym_res-2_desc-brain_mask.nii.gz")))
  writeLines('{"RepetitionTime": 1.0, "TaskName": "learn"}',
    file.path(bids_aroma_dir, "sub-01", "func", paste0(stem_a, "_bold.json")))
  utils::write.table(ev_a, file.path(bids_aroma_dir, "sub-01", "func", paste0(stem_a, "_events.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE)

  csf_raw <- u_csf + 0.085 * ic4_mot - 0.065 * ic5_mot + 0.050 * ic7_mot
  wm_raw  <- u_wm  - 0.070 * ic4_mot + 0.080 * ic6_mot - 0.050 * ic8_mot
  fd_a    <- c(NA_real_, round(runif(n_vols_a - 1L, 0.05, 0.24), 4))
  conf_a  <- data.frame(
    csf = round(csf_raw, 6),
    white_matter = round(wm_raw, 6),
    trans_x = round(0.08 * ic3_mot, 6),
    trans_y = round(0.10 * ic4_mot, 6),
    trans_z = round(0.09 * ic5_mot, 6),
    rot_x   = round(0.002 * ic6_mot, 6),
    rot_y   = round(0.002 * ic7_mot, 6),
    rot_z   = round(0.002 * ic8_mot, 6),
    framewise_displacement = ifelse(is.na(fd_a), "n/a", sprintf("%.4f", fd_a)),
    stringsAsFactors = FALSE
  )
  utils::write.table(conf_a, file.path(fprep_aroma_dir, "sub-01", "func",
    paste0(stem_a, "_desc-confounds_timeseries.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)

  # Sibling fMRIPost-AROMA outputs (omits space- entity, matching fMRIPost-AROMA >= 0.1)
  mix_mat <- round(cbind(ic1_acc, ic2_acc, ic3_mot, ic4_mot, ic5_mot, ic6_mot, ic7_mot, ic8_mot), 6)
  utils::write.table(mix_mat, file.path(post_aroma_dir, "sub-01", "func",
    paste0(stem_a, "_res-2_desc-melodic_mixing.tsv")), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  metrics_df <- data.frame(
    IC = paste0("IC_", 1:8),
    classification = c("accepted", "accepted", rep("rejected", 6L)),
    rationale = c("", "", rep("edge_fraction;high_frequency_content", 6L)),
    stringsAsFactors = FALSE
  )
  utils::write.table(metrics_df, file.path(post_aroma_dir, "sub-01", "func",
    paste0(stem_a, "_res-2_desc-aroma_metrics.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
}

run_aroma_regime <- function(lbl, aroma_mode, clean_nuis) {
  sub_a <- read_fmriprep_subject(
    fmriprep_dir         = fprep_aroma_dir,
    bids_dir             = bids_aroma_dir,
    sub                  = "01",
    task                 = "learn",
    atlas                = atlas_aroma_path,
    strategy             = "csf_wm",
    aroma                = aroma_mode,
    aroma_clean_nuisance = clean_nuis,
    fd_threshold         = 0.5
  )
  hfit <- estimate_crossrun_hrf(
    sub_a, kernel_grid = candidates, pooling = "parcel",
    amplitude_model = "trial", selection_mode = "soft_average"
  )
  dec <- deconvolve_crossrun_hrf(
    sub_a, hfit, scaling = "noise_sd", penalty = "sobolev",
    nuisance_mode = "joint", observation_mask = "fit"
  )
  ep <- decon_epoch_regression(
    dec, covariates = "rpe", target_condition = "feedback",
    lags = -2:4, aggregate_by_parcel = TRUE,
    restore_scale = TRUE, return_panel = FALSE
  )
  mot_r2 <- mean(vapply(seq_along(sub_a$bold), function(r_idx) {
    Y_p1 <- rowMeans(sub_a$bold[[r_idx]][, 1:6, drop = FALSE])
    Z_r  <- cbind(1, seq_along(Y_p1), sub_a$nuisance[[r_idx]])
    res_post <- qr.resid(qr(Z_r), Y_p1)
    max(0, summary(stats::lm(res_post ~ mot_ic_list[[r_idx]][, 2:6, drop = FALSE]))$adj.r.squared)
  }, numeric(1)))

  cf <- as.data.table(ep$coefficients)
  p1 <- as.character(unique(cf$parcel)[1L])
  sub_cf <- cf[as.character(parcel) == p1 & term == "rpe"]
  list(
    metrics = data.table(
      regime          = lbl,
      est_noise_sd    = mean(vapply(dec$scaling, function(s) mean(s$noise_sd[1:6]), numeric(1))),
      resid_motion_r2 = mot_r2,
      lag0_beta       = sub_cf[lag == 0, beta_gls],
      lag0_se         = sub_cf[lag == 0, se_gls],
      lag0_bias       = sub_cf[lag == 0, beta_gls] - 0.65
    ),
    profile = sub_cf[, .(lag, beta = beta_gls, se = se_gls, regime = lbl)]
  )
}

r_none   <- run_aroma_regime("No AROMA (csf_wm only)", "none", FALSE)
r_aggr   <- run_aroma_regime("Aggressive AROMA (aroma = 'aggressive')", "aggressive", TRUE)
r_uncl   <- run_aroma_regime("Non-aggressive AROMA + uncleaned csf_wm", "nonaggressive", FALSE)
r_nonagg <- run_aroma_regime("Non-aggressive AROMA + matched csf_wm (Recommended)", "nonaggressive", TRUE)

aroma_tbl      <- rbind(r_none$metrics, r_aggr$metrics, r_uncl$metrics, r_nonagg$metrics)
aroma_profiles <- rbind(r_none$profile, r_aggr$profile, r_uncl$profile, r_nonagg$profile)
print(aroma_tbl)
fwrite(aroma_tbl, "aroma_denoising_metrics.csv")
unlink(aroma_root, recursive = TRUE)

saveRDS(list(
  true_peak_sec = true_peak_sec,
  t_hrf = t_hrf,
  hrf_curves = data.table(
    time = rep(t_hrf, 3L),
    value = c(h_true, fit_intact$selected_kernels[[1L]][, 1L], fit_broken$selected_kernels[[1L]][, 1L]),
    curve = rep(c("Generating HRF", "Intact clock (`valid_rows`)", "Deleted 4 dummy TRs (no `time_offset`)"), each = length(t_hrf))
  ),
  clock_profiles = rbind(
    cf_intact[as.character(parcel) == p1_lbl & term == "rpe", .(lag, beta = beta_gls, se = se_gls, mode = "Intact clock (`valid_rows`)")],
    cf_offset[as.character(parcel) == p1_lbl & term == "rpe", .(lag, beta = beta_gls, se = se_gls, mode = "Truncated + `time_offset = 4`")],
    cf_broken[as.character(parcel) == p1_lbl & term == "rpe", .(lag, beta = beta_gls, se = se_gls, mode = "Deleted 4 dummy TRs (no `time_offset`)")]
  ),
  missed_trace = data.table(
    time = 0:(n_vols - 1L),
    intact = act_intact_p1,
    dropped = act_drop_p1
  ),
  run1_missed_time = run1_missed_time,
  run1_valid_time  = run1_valid_time,
  qc_summary = qc_summary,
  clock_tbl = clock_tbl,
  missed_tbl = missed_tbl,
  aroma_tbl = aroma_tbl,
  aroma_profiles = aroma_profiles,
  design_benchmarks = design_benchmarks
), "step0_examples.rds")

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
unlink(bids_root, recursive = TRUE)
cat("Step 0 simulation & BIDS ingestion checks complete.\n")
