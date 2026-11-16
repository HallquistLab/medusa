#!/usr/bin/env Rscript
# ==============================================================================
# Benchmark: Sparse Bilateral HRF & Noise Graph (pooling = "graph")
#
# Validates:
#   Scenario 1: Shifted person-specific functional/neurovascular boundaries
#               (Kong et al., 2021, Cerebral Cortex, PMID 33942058), 2 mm
#               out-of-mask sulcal CSF gap (geodesic vs Euclidean barrier),
#               and high-rho / high-variance pial draining-vein border rim.
#   Scenario 2: Spatially smooth 3D AR(3) null calibration (FWHM = 6 mm) for
#               both condition- and trial-amplitude HRF screening + Step 3
#               epoch regression Type I error control.
#   Scenario 3: Whole-brain 3D volumetric scalability up to V = 50,000 voxels
#               (geodesic CSR graph construction, bilateral HRF/noise pooling,
#               48-prototype quantization, and GLS deconvolution).
# ==============================================================================

Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

if (!("--installed" %in% commandArgs(trailingOnly = TRUE)) &&
    file.exists("DESCRIPTION") &&
    requireNamespace("pkgload", quietly = TRUE)) {
  tryCatch(
    suppressPackageStartupMessages(pkgload::load_all(".", quiet = TRUE)),
    error = function(e) suppressPackageStartupMessages(library(medusa))
  )
} else {
  suppressPackageStartupMessages(library(medusa))
}
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
smoke_mode <- "--smoke" %in% args

out_dir <- if (smoke_mode) {
  file.path(tempdir(), "sparse_graph_results_smoke")
} else if (dir.exists("inst/validation")) {
  "inst/validation/sparse_graph_results"
} else {
  file.path(getwd(), "sparse_graph_results")
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

hrf_event_design_cpp_fn <- getFromNamespace("hrf_event_design_cpp", "medusa")

generate_subject_schedule <- function(
  n_runs = 2L,
  n_trials = 30L,
  tr = 1.0,
  iti_mean = 3.5,
  iti_min = 1.2,
  iti_max = 9.0,
  cue_fb_min = 1.5,
  cue_fb_max = 2.6,
  rt_cor = 0.0,
  confound_cor = 0.0,
  x_ar1 = 0.0,
  isi_x_cor = 0.0,
  run_x_sds = NULL,
  x_run_mean_sd = 0.0,
  x_subj_mean = 0.0,
  x2_cor = 0.0,
  include_x2 = FALSE,
  tail_buffer_s = 18.0
) {
  if (is.null(run_x_sds)) run_x_sds <- rep(1.0, n_runs)
  runs_list <- vector("list", n_runs)
  names(runs_list) <- paste0("run", seq_len(n_runs))

  # Compute a common n_time across runs sufficient to hold all trials + tail
  expected_trial_span <- (cue_fb_min + cue_fb_max) / 2 + iti_mean + 1.5
  n_time <- as.integer(ceiling((10.0 + n_trials * expected_trial_span + tail_buffer_s) / tr))
  n_time <- max(n_time, as.integer(ceiling((12.0 + n_trials * (cue_fb_max + iti_min + 1.5) + tail_buffer_s) / tr)))

  for (r in seq_len(n_runs)) {
    # Generate target parametric modulator x (optionally with trial-level AR(1) autocorrelation)
    x_innov <- stats::rnorm(n_trials)
    x_raw <- numeric(n_trials)
    x_raw[1L] <- x_innov[1L]
    if (n_trials > 1L && abs(x_ar1) > 1e-6) {
      s_x <- sqrt(max(0, 1 - x_ar1^2))
      for (i in 2:n_trials) {
        x_raw[i] <- x_ar1 * x_raw[i - 1L] + s_x * x_innov[i]
      }
    } else {
      x_raw <- x_innov
    }
    x_within <- as.numeric(scale(x_raw)) * run_x_sds[min(r, length(run_x_sds))]
    x_run_shift <- if (abs(x_run_mean_sd) > 1e-8) stats::rnorm(1L, 0, x_run_mean_sd) else 0.0
    x <- x_within + x_subj_mean + x_run_shift
    x_std <- as.numeric(scale(x_within))

    # Optional correlated distractor covariate x2
    x2 <- NULL
    if (isTRUE(include_x2) || abs(x2_cor) > 1e-8) {
      x2_orth <- as.numeric(scale(stats::resid(stats::lm(stats::rnorm(n_trials) ~ x_std))))
      x2 <- x2_cor * x_std + sqrt(max(0, 1 - x2_cor^2)) * x2_orth
    }

    # Draw ITIs from truncated exponential distribution
    raw_iti <- stats::rexp(n_trials, rate = 1 / max(0.5, iti_mean - iti_min)) + iti_min
    itis <- pmin(iti_max, raw_iti)

    # Cue-to-feedback ISI (optionally correlated with x so shorter ISI co-occurs with larger x)
    if (abs(isi_x_cor) > 1e-6) {
      u_orth <- as.numeric(scale(stats::resid(stats::lm(stats::rnorm(n_trials) ~ x_std))))
      isi_z <- isi_x_cor * x_std + sqrt(max(0, 1 - isi_x_cor^2)) * u_orth
      isi_u <- stats::pnorm(isi_z)
      cue_to_fb <- cue_fb_min + (cue_fb_max - cue_fb_min) * isi_u
    } else {
      cue_to_fb <- stats::runif(n_trials, min = cue_fb_min, max = cue_fb_max)
    }

    # Rescale ITIs if total schedule exceeds run window
    max_usable_s <- (n_time * tr) - tail_buffer_s - 6.0
    total_s <- sum(cue_to_fb) + sum(itis)
    if (total_s > max_usable_s) {
      scale_f <- max(0.4, (max_usable_s - sum(cue_to_fb)) / sum(itis))
      itis <- pmax(iti_min * 0.75, itis * scale_f)
    }

    t_cue <- numeric(n_trials)
    t_fb <- numeric(n_trials)
    cursor <- 6.0 + stats::runif(1L, 0.0, 0.8)
    for (i in seq_len(n_trials)) {
      t_cue[i] <- cursor
      t_fb[i] <- cursor + cue_to_fb[i]
      cursor <- t_fb[i] + itis[i]
    }

    # Generate reaction time / choice duration (optionally correlated with x)
    z_orth <- as.numeric(scale(stats::resid(stats::lm(stats::rnorm(n_trials) ~ x))))
    rt_z <- rt_cor * x_std + sqrt(max(0, 1 - rt_cor^2)) * z_orth
    rt <- pmin(2.8, pmax(0.45, 1.45 + 0.55 * rt_z))

    # Preceding cue / surprise / trial confounds
    c_orth <- as.numeric(scale(stats::resid(stats::lm(stats::rnorm(n_trials) ~ x))))
    cue_mod <- confound_cor * x_std + sqrt(max(0, 1 - confound_cor^2)) * c_orth

    trial_df <- data.frame(
      run = paste0("run", r),
      trial = seq_len(n_trials),
      onset_cue = t_cue,
      onset_fb = t_fb,
      iti = itis,
      cue_fb_isi = cue_to_fb,
      rt = rt,
      x = x,
      x_within = x_within,
      x_run_mean = x_subj_mean + x_run_shift,
      cue_mod = cue_mod,
      stringsAsFactors = FALSE
    )
    if (!is.null(x2)) {
      trial_df$x2 <- x2
    }

    events_df <- rbind(
      data.frame(onset = t_cue, duration = 0, condition = "cue", trial = seq_len(n_trials)),
      data.frame(onset = t_fb, duration = 0, condition = "feedback", trial = seq_len(n_trials))
    )
    events_df <- events_df[order(events_df$onset), ]
    rownames(events_df) <- NULL

    runs_list[[r]] <- list(
      trials = trial_df,
      events = events_df,
      n_time = n_time
    )
  }

  list(runs = runs_list, n_time = n_time, tr = tr)
}

# Continuous peak time of a discrete HRF kernel via quadratic interpolation
hrf_continuous_peak_time <- function(h, tr = 1.0) {
  idx <- which.max(h)
  if (idx <= 1L || idx >= length(h)) return((idx - 1L) * tr)
  y1 <- h[idx - 1L]; y2 <- h[idx]; y3 <- h[idx + 1L]
  denom <- (y1 - 2 * y2 + y3)
  shift <- if (abs(denom) > 1e-12) 0.5 * (y1 - y3) / denom else 0.0
  shift <- max(-0.5, min(0.5, shift))
  ((idx - 1L) + shift) * tr
}

# Build unit-L2 spatial smoothing operator over 3D coordinates within a bank
make_spatial_smoother <- function(coords_mat, fwhm_mm = 5.0) {
  nv <- nrow(coords_mat)
  if (nv <= 1L || fwhm_mm <= 0) return(diag(nv))
  sigma_sp <- fwhm_mm / sqrt(8 * log(2))
  dmat <- as.matrix(stats::dist(coords_mat))
  k_sp <- exp(-0.5 * (dmat / sigma_sp)^2)
  sweep(k_sp, 1L, sqrt(rowSums(k_sp^2)), "/")
}

# Simulate stationary AR(p) timecourses from unit-variance spatially correlated innovations
sim_arp_from_innov <- function(Z_mat, ar_coefs = c(0.38), sd_vec = 1.0) {
  n_t <- nrow(Z_mat)
  n_v <- ncol(Z_mat)
  if (length(sd_vec) == 1L) sd_vec <- rep(sd_vec, n_v)
  E <- matrix(0.0, n_t, n_v)
  p_ord <- length(ar_coefs)
  # Empirical burn-in so marginal stationary SD matches sd_vec accurately
  n_burn <- 40L
  Z_full <- rbind(matrix(stats::rnorm(n_burn * n_v), n_burn, n_v), Z_mat)
  X_full <- matrix(0.0, n_burn + n_t, n_v)
  for (tt in seq_len(n_burn + n_t)) {
    acc <- Z_full[tt, ]
    for (lag_k in seq_len(p_ord)) {
      if (tt > lag_k) acc <- acc + ar_coefs[lag_k] * X_full[tt - lag_k, ]
    }
    X_full[tt, ] <- acc
  }
  # Theoretical AR(1) stationary SD when p = 1; empirical standardization fallback for AR(p)
  if (p_ord == 1L) {
    ar_sd <- 1.0 / sqrt(max(1e-6, 1.0 - ar_coefs[1L]^2))
  } else {
    # Compute impulse response L2 norm for exact stationary SD of AR(p)
    ir <- numeric(80L)
    ir[1L] <- 1.0
    for (tt in 2:80L) {
      s <- 0.0
      for (lag_k in seq_len(p_ord)) if (tt > lag_k) s <- s + ar_coefs[lag_k] * ir[tt - lag_k]
      ir[tt] <- s
    }
    ar_sd <- sqrt(sum(ir^2))
  }
  E <- sweep(X_full[(n_burn + 1L):(n_burn + n_t), , drop = FALSE] / ar_sd, 2L, sd_vec, "*")
  E
}

cat("======================================================================\n")
cat("SCENARIO 1: Person-Specific Shifted Boundaries + Sulcal Gap + Pial Vein\n")
cat("======================================================================\n")

set.seed(20261114L)
tr <- 1.0
n_subj_s1 <- if (smoke_mode) 2L else 24L
n_runs_s1 <- 3L
n_trials_s1 <- 26L
true_beta_s1 <- 0.65

# Candidate HRF library (15 shapes: 5 peak delays x 3 dispersions)
lib_s1 <- hrf_library(
  tr = tr,
  peak_delays = c(4.0, 5.0, 6.0, 7.0, 8.0),
  dispersions = c(0.9, 1.0, 1.1),
  undershoot_ratios = 6.0
)
pars_s1 <- attr(lib_s1, "parameters")
idx_can_s1 <- which(abs(pars_s1$peak_delay - 6.0) < 1e-6 & abs(pars_s1$dispersion - 1.0) < 1e-6)[1L]
n_kernel_s1 <- nrow(lib_s1)

# 3D Cortical Patch Geometry (2 mm isotropic voxels):
#   - Bank 1 (i = 1..12, j = 1..5, k = 1..2): 120 voxels
#       * Fixed Group Atlas splits at i <= 6 ("Atlas_P1") vs i >= 7 ("Atlas_P2")
#       * Pial-vein border rim at j == 5 (24 voxels: rho = 0.78, sd = 0.75)
#       * GM core at j <= 4 (96 voxels: rho = 0.38, sd = 0.26)
#   - Out-of-mask Sulcal CSF Gap at i == 13 (2 mm gap)
#   - Bank 2 (i = 14..16, j = 1..5, k = 1..2): 30 voxels ("Atlas_P3", pure null)
bank1_df <- expand.grid(i = 1:12, j = 1:5, k = 1:2)
bank2_df <- expand.grid(i = 14:16, j = 1:5, k = 1:2)
vm_s1 <- rbind(bank1_df, bank2_df)
vm_s1$x <- vm_s1$i * 2.0
vm_s1$y <- vm_s1$j * 2.0
vm_s1$z <- vm_s1$k * 2.0
vm_s1$voxel <- sprintf("v%d", seq_len(nrow(vm_s1)))
nvox_s1 <- nrow(vm_s1)
n_b1 <- nrow(bank1_df)
coords_s1 <- as.matrix(vm_s1[, c("x", "y", "z")])
rownames(coords_s1) <- vm_s1$voxel

atlas_groups_s1 <- stats::setNames(
  ifelse(vm_s1$i <= 6L, "Atlas_P1", ifelse(vm_s1$i <= 12L, "Atlas_P2", "Atlas_P3")),
  vm_s1$voxel
)
is_bank1 <- vm_s1$i <= 12L
is_vein  <- is_bank1 & (vm_s1$j == 5L)
is_gm1   <- is_bank1 & (vm_s1$j <= 4L)
is_sulc_border_b2 <- (vm_s1$i == 14L) # opposite sulcal bank border (4 mm Euclidean from i = 12)
is_bank2 <- vm_s1$i >= 14L

# Prebuild spatial graphs on this cortical mask
geom_s1 <- list(voxel_map = vm_s1, geometry = list(spacing = c(2.0, 2.0, 2.0)))
g_geo_free <- build_spatial_graph(geom_s1, spatial_fwhm = 5.5, metric = "geodesic",
                                  parcel_boundary = "ignore")
g_geo_soft <- build_spatial_graph(geom_s1, spatial_fwhm = 5.5, metric = "geodesic",
                                  spatial_groups = atlas_groups_s1, parcel_boundary = "soft")
g_euc_free <- build_spatial_graph(geom_s1, spatial_fwhm = 5.5, metric = "euclidean",
                                  parcel_boundary = "ignore")

K_sp_b1 <- make_spatial_smoother(coords_s1[1:n_b1, , drop = FALSE], fwhm_mm = 5.0)
K_sp_b2 <- make_spatial_smoother(coords_s1[(n_b1 + 1L):nvox_s1, , drop = FALSE], fwhm_mm = 5.0)

# Canonical single-column library for fixed-HRF baseline
lib_can_only <- lib_s1[, idx_can_s1, drop = FALSE]

shifts_pool <- c(-3L, -2L, 2L, 3L) # +/- 4 to 6 mm individual boundary shifts
s1_rows <- list()

t0_s1 <- proc.time()[["elapsed"]]
for (s in seq_len(n_subj_s1)) {
  shift_s <- shifts_pool[((s - 1L) %% length(shifts_pool)) + 1L]
  true_split_i <- 6L + shift_s # Atlas split is 6; true individual split is 3, 4, 8, or 9

  # Continuous off-grid subject-specific HRFs for Fast Territory A and Slow Territory B
  pk_fast <- min(5.0, max(3.9, 4.35 + stats::rnorm(1L, 0, 0.16)))
  dp_fast <- min(1.02, max(0.85, 0.92 + stats::rnorm(1L, 0, 0.03)))
  pk_slow <- min(8.1, max(6.8, 7.45 + stats::rnorm(1L, 0, 0.16)))
  dp_slow <- min(1.18, max(0.98, 1.08 + stats::rnorm(1L, 0, 0.03)))

  h_fast <- hrf_library(tr, peak_delays = pk_fast, dispersions = dp_fast, undershoot_ratios = 6.0)[, 1L]
  h_slow <- hrf_library(tr, peak_delays = pk_slow, dispersions = dp_slow, undershoot_ratios = 6.0)[, 1L]
  h_can  <- lib_s1[, idx_can_s1]

  true_k_mat <- matrix(0.0, n_kernel_s1, nvox_s1)
  true_peak_vec <- numeric(nvox_s1)
  for (v in seq_len(nvox_s1)) {
    if (!is_bank1[v]) {
      true_k_mat[, v] <- h_can
      true_peak_vec[v] <- hrf_continuous_peak_time(h_can, tr)
    } else if (vm_s1$i[v] <= true_split_i) {
      true_k_mat[, v] <- h_fast
      true_peak_vec[v] <- hrf_continuous_peak_time(h_fast, tr)
    } else {
      true_k_mat[, v] <- h_slow
      true_peak_vec[v] <- hrf_continuous_peak_time(h_slow, tr)
    }
  }

  is_shifted_gm <- is_gm1 & (
    (vm_s1$i <= 6L & vm_s1$i > true_split_i) |
    (vm_s1$i > 6L  & vm_s1$i <= true_split_i)
  )
  is_core_gm <- is_gm1 & !is_shifted_gm

  true_rho_vec <- ifelse(is_vein, 0.78, 0.38)
  true_sd_vec  <- ifelse(is_vein, 0.75, 0.26)

  sched <- generate_subject_schedule(n_runs = n_runs_s1, n_trials = n_trials_s1, tr = tr)
  n_t <- sched$n_time

  bold_list <- events_list <- x_list <- neural_true_list <- vector("list", n_runs_s1)
  names(bold_list) <- names(events_list) <- names(x_list) <- names(neural_true_list) <- paste0("run", seq_len(n_runs_s1))

  h_pair <- cbind(fast = h_fast, slow = h_slow)
  for (r in seq_len(n_runs_s1)) {
    r_info <- sched$runs[[r]]
    tr_df <- r_info$trials
    events_list[[r]] <- r_info$events
    x_list[[r]] <- tr_df$x

    amp_cue <- 0.55 + 0.15 * stats::rnorm(n_trials_s1)
    amp_fb  <- 0.90 + true_beta_s1 * tr_df$x + stats::rnorm(n_trials_s1, 0, 0.20)

    # Underlying discrete TR neural stick function for Step 2 correlation validation
    u_true <- numeric(n_t)
    idx_c <- pmin(n_t, pmax(1L, round(tr_df$onset_cue / tr) + 1L))
    idx_f <- pmin(n_t, pmax(1L, round(tr_df$onset_fb / tr) + 1L))
    for (ti in seq_len(n_trials_s1)) {
      u_true[idx_c[ti]] <- u_true[idx_c[ti]] + amp_cue[ti]
      u_true[idx_f[ti]] <- u_true[idx_f[ti]] + amp_fb[ti]
    }
    neural_true_list[[r]] <- u_true

    des_cue <- hrf_event_design_cpp_fn(n_t, tr, 0, tr_df$onset_cue, rep(0, n_trials_s1),
                                             seq_len(n_trials_s1), n_trials_s1, h_pair)
    des_fb  <- hrf_event_design_cpp_fn(n_t, tr, 0, tr_df$onset_fb, rep(0, n_trials_s1),
                                             seq_len(n_trials_s1), n_trials_s1, h_pair)
    sig_fast <- as.numeric(des_cue[[1L]] %*% amp_cue + des_fb[[1L]] %*% amp_fb)
    sig_slow <- as.numeric(des_cue[[2L]] %*% amp_cue + des_fb[[2L]] %*% amp_fb)

    # Spatially correlated innovations within Bank 1 and Bank 2 separately (CSF gap in between)
    Z_b1 <- matrix(stats::rnorm(n_t * n_b1), n_t, n_b1) %*% t(K_sp_b1)
    Z_b2 <- matrix(stats::rnorm(n_t * (nvox_s1 - n_b1)), n_t, nvox_s1 - n_b1) %*% t(K_sp_b2)
    Z_all <- cbind(Z_b1, Z_b2)

    E_gm   <- sim_arp_from_innov(Z_all, ar_coefs = 0.38, sd_vec = true_sd_vec)
    E_vein <- sim_arp_from_innov(Z_all[, is_vein, drop = FALSE], ar_coefs = 0.78, sd_vec = true_sd_vec[is_vein])
    E_all  <- E_gm
    E_all[, is_vein] <- E_vein

    Y_r <- matrix(100.0, n_t, nvox_s1, dimnames = list(NULL, vm_s1$voxel))
    for (v in seq_len(nvox_s1)) {
      if (is_gm1[v]) {
        gain_v <- if (is_shifted_gm[v]) 0.85 else 1.00 # slightly weaker SNR right on shifted border
        sig_v <- if (vm_s1$i[v] <= true_split_i) sig_fast else sig_slow
        Y_r[, v] <- 100.0 + gain_v * sig_v + E_all[, v]
      } else if (is_vein[v]) {
        # Pial vein has weak vascular spillover + heavy AR(1) physiological noise
        sig_v <- if (vm_s1$i[v] <= true_split_i) sig_fast else sig_slow
        Y_r[, v] <- 100.0 + 0.35 * sig_v + E_all[, v]
      } else {
        # Opposite gyral bank across CSF gap: pure null noise
        Y_r[, v] <- 100.0 + E_all[, v]
      }
    }
    bold_list[[r]] <- Y_r
  }

  true_groups_s <- stats::setNames(
    ifelse(vm_s1$i <= true_split_i, "True_Fast", ifelse(vm_s1$i <= 12L, "True_Slow", "True_Null")),
    vm_s1$voxel
  )

  # Fit all 8 pipelines (including person-specific Oracle parcellation ceiling)
  fits <- list(
    "1. Canonical Fixed (6.0s)" = estimate_crossrun_hrf(
      bold_list, events_list, tr, lib_can_only, fallback_kernel = 1L,
      pooling = "voxel", amplitude_model = "trial", selection_mode = "hard"
    ),
    "2. Voxel Soft-Average" = estimate_crossrun_hrf(
      bold_list, events_list, tr, lib_s1, fallback_kernel = idx_can_s1,
      pooling = "voxel", amplitude_model = "trial", selection_mode = "soft_average"
    ),
    "3. Group Parcel Soft-Average" = estimate_crossrun_hrf(
      bold_list, events_list, tr, lib_s1, spatial_groups = atlas_groups_s1,
      fallback_kernel = idx_can_s1, pooling = "parcel",
      amplitude_model = "trial", selection_mode = "soft_average"
    ),
    "4. Voxel + Euclidean Smooth" = estimate_crossrun_hrf(
      bold_list, events_list, tr, lib_s1, spatial_coords = coords_s1,
      spatial_fwhm = 5.5, spatial_strength = 0.25,
      fallback_kernel = idx_can_s1, pooling = "voxel",
      amplitude_model = "trial", selection_mode = "soft_average"
    ),
    "5. Graph (Euclidean, No Barrier)" = estimate_crossrun_hrf(
      bold_list, events_list, tr, lib_s1, spatial_graph = g_euc_free,
      fallback_kernel = idx_can_s1, pooling = "graph",
      amplitude_model = "trial", selection_mode = "soft_average"
    ),
    "6. Graph (Geodesic, Atlas-Free)" = estimate_crossrun_hrf(
      bold_list, events_list, tr, lib_s1, spatial_graph = g_geo_free,
      fallback_kernel = idx_can_s1, pooling = "graph",
      amplitude_model = "trial", selection_mode = "soft_average"
    ),
    "7. Graph (Geodesic + Soft Atlas)" = estimate_crossrun_hrf(
      bold_list, events_list, tr, lib_s1, spatial_groups = atlas_groups_s1,
      spatial_graph = g_geo_soft, parcel_boundary = "soft", parcel_cross_weight = 0.35,
      fallback_kernel = idx_can_s1, pooling = "graph",
      amplitude_model = "trial", selection_mode = "soft_average"
    ),
    "8. Oracle True-Territory Parcel" = estimate_crossrun_hrf(
      bold_list, events_list, tr, lib_s1, spatial_groups = true_groups_s,
      fallback_kernel = idx_can_s1, pooling = "parcel",
      amplitude_model = "trial", selection_mode = "soft_average"
    )
  )

  for (m_name in names(fits)) {
    fit_m <- fits[[m_name]]

    # 1. HRF Waveform & Peak-Latency Accuracy across folds
    rmse_gm_runs <- rmse_shift_runs <- rmse_core_runs <- rmse_sulc_runs <- numeric(n_runs_s1)
    mae_pk_shift_runs <- sulc_dep_runs <- rho_rmse_runs <- rho_vein_rmse_runs <- numeric(n_runs_s1)

    for (r in seq_len(n_runs_s1)) {
      km <- fit_m$selected_kernels[[r]]
      rmse_v <- sqrt(colMeans((km - true_k_mat)^2))
      rmse_gm_runs[r]    <- mean(rmse_v[is_gm1])
      rmse_shift_runs[r] <- mean(rmse_v[is_shifted_gm])
      rmse_core_runs[r]  <- mean(rmse_v[is_core_gm])
      rmse_sulc_runs[r]  <- mean(rmse_v[is_sulc_border_b2])

      est_pk <- apply(km[, is_shifted_gm, drop = FALSE], 2L, hrf_continuous_peak_time, tr = tr)
      mae_pk_shift_runs[r] <- mean(abs(est_pk - true_peak_vec[is_shifted_gm]))

      if (ncol(fit_m$kernels) == 1L) {
        sulc_dep_runs[r] <- 0.0
      } else {
        sulc_dep_runs[r] <- mean(fit_m$selections[[r]]$kernel_index[is_sulc_border_b2] != idx_can_s1)
      }

      rho_est <- fit_m$noise[[r]]$rho
      rho_rmse_runs[r] <- sqrt(mean((rho_est[is_bank1] - true_rho_vec[is_bank1])^2, na.rm = TRUE))
      rho_vein_rmse_runs[r] <- sqrt(mean((rho_est[is_vein] - true_rho_vec[is_vein])^2, na.rm = TRUE))
    }

    # 2. Step 2 Deconvolved Neural Timecourse Recovery
    dec_m <- deconvolve_crossrun_hrf(bold_list, fit_m, lambda = 0.8, scaling = "noise_sd", penalty = "sobolev")
    cor_decon_shift <- numeric(n_runs_s1)
    cor_decon_gm    <- numeric(n_runs_s1)
    for (r in seq_len(n_runs_s1)) {
      u_tr <- neural_true_list[[r]]
      # Smooth true stick slightly (1-TR triangle) for fair discrete comparison
      u_sm <- as.numeric(stats::filter(u_tr, c(0.25, 0.50, 0.25), sides = 2L))
      u_sm[is.na(u_sm)] <- u_tr[is.na(u_sm)]
      cors_r <- as.numeric(stats::cor(dec_m$activity[[r]], u_sm))
      cor_decon_shift[r] <- mean(cors_r[is_shifted_gm], na.rm = TRUE)
      cor_decon_gm[r]    <- mean(cors_r[is_gm1], na.rm = TRUE)
    }

    # 3. Step 3 Trial-Level Modulator Slope Recovery (lag = 0 on feedback)
    reg_m <- decon_epoch_regression(
      dec_m, covariates = x_list, events = events_list,
      target_condition = "feedback", lags = 0,
      demix_other_events = TRUE, demix_target_trials = TRUE,
      restore_scale = TRUE
    )
    coef_0 <- reg_m$coefficients[reg_m$coefficients$lag == 0 & reg_m$coefficients$term == "x1", , drop = FALSE]
    idx_map <- match(vm_s1$voxel, coef_0$voxel)
    beta_v <- coef_0$beta_corrected[idx_map]
    r_dem_v <- coef_0$r_demixed[idx_map]

    s1_rows[[length(s1_rows) + 1L]] <- data.table(
      subj = s,
      shift_vox = shift_s,
      method = m_name,
      hrf_rmse_gm = mean(rmse_gm_runs),
      hrf_rmse_shifted = mean(rmse_shift_runs),
      hrf_rmse_core = mean(rmse_core_runs),
      hrf_rmse_sulcal_b2 = mean(rmse_sulc_runs),
      peak_mae_shifted = mean(mae_pk_shift_runs),
      sulcal_bleed_rate = mean(sulc_dep_runs),
      noise_rho_rmse_all = mean(rho_rmse_runs),
      noise_rho_rmse_vein = mean(rho_vein_rmse_runs),
      decon_cor_gm = mean(cor_decon_gm),
      decon_cor_shifted = mean(cor_decon_shift),
      beta_mean_shifted = mean(beta_v[is_shifted_gm], na.rm = TRUE),
      beta_rmse_shifted = sqrt(mean((beta_v[is_shifted_gm] - true_beta_s1)^2, na.rm = TRUE)),
      r_trial_shifted = mean(r_dem_v[is_shifted_gm], na.rm = TRUE),
      beta_mean_sulcal_b2 = mean(abs(beta_v[is_sulc_border_b2]), na.rm = TRUE)
    )
  }
  if (s %% 6L == 0L) {
    cat(sprintf("  Completed subject %2d / %2d (%.1fs elapsed)\n",
                s, n_subj_s1, proc.time()[["elapsed"]] - t0_s1))
  }
}

s1_dt <- rbindlist(s1_rows)
s1_summary <- s1_dt[, .(
  hrf_rmse_gm = mean(hrf_rmse_gm),
  hrf_rmse_shifted = mean(hrf_rmse_shifted),
  hrf_rmse_core = mean(hrf_rmse_core),
  peak_mae_shifted_s = mean(peak_mae_shifted),
  sulcal_bleed_rate = mean(sulcal_bleed_rate),
  noise_rho_rmse_vein = mean(noise_rho_rmse_vein),
  decon_cor_shifted = mean(decon_cor_shifted),
  beta_shifted_mean = mean(beta_mean_shifted),
  beta_shifted_rmse = mean(beta_rmse_shifted),
  r_trial_shifted = mean(r_trial_shifted),
  cohort_t_shifted = mean(beta_mean_shifted) / (stats::sd(beta_mean_shifted) / sqrt(.N))
), by = .(method)]

cat("\n--- SCENARIO 1 SUMMARY (N = 24 Subjects, True Modulator Slope = 0.65) ---\n")
print(s1_summary, digits = 3, nrows = 50)
fwrite(s1_dt, file.path(out_dir, "scenario1_subject_metrics.csv"))
fwrite(s1_summary, file.path(out_dir, "scenario1_summary.csv"))


cat("\n======================================================================\n")
cat("SCENARIO 2: Spatially Smooth 3D AR(3) Null Calibration (FWHM = 6 mm)\n")
cat("======================================================================\n")

set.seed(20261115L)
n_subj_s2 <- if (smoke_mode) 2L else 24L
n_runs_s2 <- 3L
n_trials_s2 <- 24L
lib_s2 <- hrf_library(tr = 1.0, preset = "medusa")
idx_can_s2 <- 1L

vm_s2 <- expand.grid(i = 1:6, j = 1:6, k = 1:3) # 108 voxels at 2 mm spacing
vm_s2$x <- vm_s2$i * 2.0
vm_s2$y <- vm_s2$j * 2.0
vm_s2$z <- vm_s2$k * 2.0
vm_s2$voxel <- sprintf("v%d", seq_len(nrow(vm_s2)))
nvox_s2 <- nrow(vm_s2)
atlas_s2 <- stats::setNames(ifelse(vm_s2$i <= 3L, "Null_P1", "Null_P2"), vm_s2$voxel)
K_sp_s2 <- make_spatial_smoother(as.matrix(vm_s2[, c("x", "y", "z")]), fwhm_mm = 6.0)
geom_s2 <- list(voxel_map = vm_s2, geometry = list(spacing = c(2, 2, 2)))
g_s2_free <- build_spatial_graph(geom_s2, spatial_fwhm = 5.5, metric = "geodesic")
g_s2_soft <- build_spatial_graph(geom_s2, spatial_fwhm = 5.5, metric = "geodesic",
                                 spatial_groups = atlas_s2, parcel_boundary = "soft")

s2_rows <- list()
for (s in seq_len(n_subj_s2)) {
  sched2 <- generate_subject_schedule(n_runs = n_runs_s2, n_trials = n_trials_s2, tr = 1.0)
  n_t2 <- sched2$n_time
  bold2 <- events2 <- x2 <- vector("list", n_runs_s2)
  names(bold2) <- names(events2) <- names(x2) <- paste0("run", seq_len(n_runs_s2))
  for (r in seq_len(n_runs_s2)) {
    events2[[r]] <- sched2$runs[[r]]$events
    x2[[r]] <- sched2$runs[[r]]$trials$x
    Z_null <- matrix(stats::rnorm(n_t2 * nvox_s2), n_t2, nvox_s2) %*% t(K_sp_s2)
    E_null <- sim_arp_from_innov(Z_null, ar_coefs = c(0.45, -0.15, 0.08), sd_vec = 0.35)
    bold2[[r]] <- 100.0 + E_null
    colnames(bold2[[r]]) <- vm_s2$voxel
  }

  for (amp_mod in c("condition", "trial")) {
    null_fits <- list(
      "Voxel Soft-Average" = estimate_crossrun_hrf(
        bold2, events2, 1.0, lib_s2, fallback_kernel = idx_can_s2,
        pooling = "voxel", amplitude_model = amp_mod, selection_mode = "soft_average"
      ),
      "Group Parcel Soft-Average" = estimate_crossrun_hrf(
        bold2, events2, 1.0, lib_s2, spatial_groups = atlas_s2, fallback_kernel = idx_can_s2,
        pooling = "parcel", amplitude_model = amp_mod, selection_mode = "soft_average"
      ),
      "Graph (Geodesic, Atlas-Free)" = estimate_crossrun_hrf(
        bold2, events2, 1.0, lib_s2, spatial_graph = g_s2_free, fallback_kernel = idx_can_s2,
        pooling = "graph", amplitude_model = amp_mod, selection_mode = "soft_average"
      ),
      "Graph (Geodesic + Soft Atlas)" = estimate_crossrun_hrf(
        bold2, events2, 1.0, lib_s2, spatial_groups = atlas_s2, spatial_graph = g_s2_soft,
        parcel_boundary = "soft", fallback_kernel = idx_can_s2,
        pooling = "graph", amplitude_model = amp_mod, selection_mode = "soft_average"
      )
    )

    for (m_nm in names(null_fits)) {
      fm <- null_fits[[m_nm]]
      dep_r <- mean(vapply(seq_len(n_runs_s2), function(r) {
        mean(fm$selections[[r]]$kernel_index != idx_can_s2)
      }, numeric(1)))
      non_fb_w <- mean(vapply(seq_len(n_runs_s2), function(r) {
        mean(1.0 - fm$kernel_weights[[r]][idx_can_s2, ])
      }, numeric(1)))

      # Evaluate Step 3 epoch regression null Type I error on the trial branch
      vox_fpr <- NA_real_
      mean_beta_null <- NA_real_
      if (amp_mod == "trial") {
        dec_null <- deconvolve_crossrun_hrf(bold2, fm, lambda = 0.8, scaling = "noise_sd", penalty = "sobolev")
        reg_null <- decon_epoch_regression(
          dec_null, covariates = x2, events = events2,
          target_condition = "feedback", lags = 0,
          demix_other_events = TRUE, demix_target_trials = TRUE,
          restore_scale = TRUE
        )
        c0 <- reg_null$coefficients[reg_null$coefficients$lag == 0 & reg_null$coefficients$term == "x1", , drop = FALSE]
        vox_fpr <- mean(c0$p_value < 0.05, na.rm = TRUE)
        mean_beta_null <- mean(c0$beta_corrected, na.rm = TRUE)
      }

      s2_rows[[length(s2_rows) + 1L]] <- data.table(
        subj = s,
        amplitude_model = amp_mod,
        method = m_nm,
        false_departure_rate = dep_r,
        mean_non_fallback_weight = non_fb_w,
        voxel_step3_fpr = vox_fpr,
        subj_mean_beta_null = mean_beta_null
      )
    }
  }
}

s2_dt <- rbindlist(s2_rows)
s2_summary <- s2_dt[, .(
  false_departure_rate = mean(false_departure_rate),
  max_subj_departure = max(false_departure_rate),
  mean_non_fallback_weight = mean(mean_non_fallback_weight),
  voxel_step3_fpr = mean(voxel_step3_fpr, na.rm = TRUE),
  cohort_group_p = if (all(is.na(subj_mean_beta_null))) NA_real_ else stats::t.test(subj_mean_beta_null, mu = 0)$p.value
), by = .(amplitude_model, method)]

cat("\n--- SCENARIO 2 SUMMARY (3D Smooth AR(3) Pure Null, FWHM = 6 mm) ---\n")
print(s2_summary, digits = 3, nrows = 50)
fwrite(s2_dt, file.path(out_dir, "scenario2_null_subject_metrics.csv"))
fwrite(s2_summary, file.path(out_dir, "scenario2_null_summary.csv"))


cat("\n======================================================================\n")
cat("SCENARIO 3: Whole-Brain 3D Volumetric Scalability (Up to V = 50,000)\n")
cat("======================================================================\n")

set.seed(20261116L)
vox_targets <- if (smoke_mode) c(300L, 600L) else c(2000L, 10000L, 25000L, 50000L)
lib_s3 <- hrf_library(tr = 1.0, peak_delays = c(4.0, 5.0, 6.0, 7.0, 8.0),
                      dispersions = c(0.9, 1.0, 1.1), undershoot_ratios = 6.0)
idx_can_s3 <- which(abs(attr(lib_s3, "parameters")$peak_delay - 6.0) < 1e-6 &
                    abs(attr(lib_s3, "parameters")$dispersion - 1.0) < 1e-6)[1L]
n_t3 <- 130L
onsets_s3 <- seq(6, 118, by = 8)
ev_s3 <- list(
  run1 = data.frame(onset = onsets_s3, duration = 0, condition = "stim"),
  run2 = data.frame(onset = onsets_s3 + 1.0, duration = 0, condition = "stim")
)
des_s3_r1 <- hrf_event_design_cpp_fn(n_t3, 1.0, 0, ev_s3$run1$onset, rep(0, length(onsets_s3)),
                                           rep(1L, length(onsets_s3)), 1L, lib_s3)
des_s3_r2 <- hrf_event_design_cpp_fn(n_t3, 1.0, 0, ev_s3$run2$onset, rep(0, length(onsets_s3)),
                                           rep(1L, length(onsets_s3)), 1L, lib_s3)

s3_rows <- list()
for (V_targ in vox_targets) {
  # Construct 3D cortical shell Lattice of exact size V_targ
  nx <- ceiling(V_targ^(1 / 3) * 1.25)
  ny <- ceiling(V_targ^(1 / 3) * 1.25)
  nz <- ceiling(V_targ / (nx * ny)) + 1L
  grid_full <- expand.grid(i = seq_len(nx), j = seq_len(ny), k = seq_len(nz))
  vm_v <- grid_full[seq_len(V_targ), , drop = FALSE]
  vm_v$x <- vm_v$i * 2.0
  vm_v$y <- vm_v$j * 2.0
  vm_v$z <- vm_v$k * 2.0
  vm_v$voxel <- sprintf("v%d", seq_len(V_targ))

  # Simulate 2 runs with spatially varying Fast / Canonical / Slow territories
  terr <- ifelse(vm_v$i <= floor(nx / 3), 2L, ifelse(vm_v$i <= floor(2 * nx / 3), idx_can_s3, 5L))
  y1 <- matrix(stats::rnorm(n_t3 * V_targ, sd = 0.35), n_t3, V_targ, dimnames = list(NULL, vm_v$voxel))
  y2 <- matrix(stats::rnorm(n_t3 * V_targ, sd = 0.35), n_t3, V_targ, dimnames = list(NULL, vm_v$voxel))
  for (k_id in unique(terr)) {
    cols <- which(terr == k_id)
    s1_vec <- as.numeric(des_s3_r1[[k_id]] %*% matrix(1.1, 1L, 1L))
    s2_vec <- as.numeric(des_s3_r2[[k_id]] %*% matrix(1.1, 1L, 1L))
    y1[, cols] <- sweep(y1[, cols, drop = FALSE], 1L, 100.0 + s1_vec, "+")
    y2[, cols] <- sweep(y2[, cols, drop = FALSE], 1L, 100.0 + s2_vec, "+")
  }
  bold_v <- list(run1 = y1, run2 = y2)

  gc( FALSE )
  t0 <- proc.time()[["elapsed"]]
  g_v <- build_spatial_graph(
    list(voxel_map = vm_v, geometry = list(spacing = c(2, 2, 2))),
    spatial_fwhm = 5.5, max_neighbors = 27L, metric = "geodesic"
  )
  t_graph <- proc.time()[["elapsed"]] - t0

  t1 <- proc.time()[["elapsed"]]
  fit_v <- estimate_crossrun_hrf(
    bold_v, ev_s3, tr = 1.0, kernel_grid = lib_s3, fallback_kernel = idx_can_s3,
    pooling = "graph", spatial_graph = g_v,
    amplitude_model = "condition", selection_mode = "soft_average",
    max_hrf_strata = 48L, max_noise_strata = 10L
  )
  t_hrf <- proc.time()[["elapsed"]] - t1

  t2 <- proc.time()[["elapsed"]]
  dec_v <- deconvolve_crossrun_hrf(bold_v, fit_v, lambda = 0.8, scaling = "noise_sd", penalty = "sobolev")
  t_decon <- proc.time()[["elapsed"]] - t2

  graph_mb <- as.numeric(utils::object.size(g_v)) / (1024^2)
  dense_equiv_mb <- (as.numeric(V_targ)^2 * 8.0) / (1024^2)
  n_proto_r1 <- length(unique(fit_v$hrf_strata$run1))

  cat(sprintf("  V = %6d | graph: %5.2fs (%.1f MB vs %.0f MB dense) | HRF fit: %5.2fs (%2d strata) | decon: %5.2fs\n",
              V_targ, t_graph, graph_mb, dense_equiv_mb, t_hrf, n_proto_r1, t_decon))

  s3_rows[[length(s3_rows) + 1L]] <- data.table(
    n_voxels = V_targ,
    mean_degree = mean(diff(g_v$row_ptr)),
    graph_build_sec = t_graph,
    hrf_fit_sec = t_hrf,
    decon_sec = t_decon,
    total_sec = t_graph + t_hrf + t_decon,
    sparse_graph_mb = graph_mb,
    dense_matrix_equiv_mb = dense_equiv_mb,
    n_hrf_prototypes = n_proto_r1
  )
  rm(bold_v, y1, y2, fit_v, dec_v, g_v)
  gc(FALSE)
}

s3_dt <- rbindlist(s3_rows)
cat("\n--- SCENARIO 3 SUMMARY (Whole-Brain Scalability) ---\n")
print(s3_dt, digits = 3)
fwrite(s3_dt, file.path(out_dir, "scenario3_scalability_summary.csv"))

# Generate diagnostic publication-style figures
p1_plot_dt <- melt(
  s1_summary,
  id.vars = "method",
  measure.vars = c("hrf_rmse_shifted", "peak_mae_shifted_s", "beta_shifted_rmse", "sulcal_bleed_rate"),
  variable.name = "metric",
  value.name = "value"
)
metric_labels <- c(
  hrf_rmse_shifted = "Shifted-Border HRF RMSE",
  peak_mae_shifted_s = "Shifted-Border Peak MAE (s)",
  beta_shifted_rmse = "Step 3 Modulator Slope RMSE",
  sulcal_bleed_rate = "Cross-Sulcal Bleed Rate"
)
p1_plot_dt[, metric_label := factor(metric_labels[as.character(metric)], levels = unname(metric_labels))]

fig1 <- ggplot(p1_plot_dt, aes(x = reorder(method, -as.integer(factor(method))), y = value, fill = method)) +
  geom_col(width = 0.72, show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ metric_label, scales = "free_x", ncol = 2) +
  labs(
    title = "Sparse Bilateral Geodesic Graph vs. Parcel & Voxel Pooling",
    subtitle = "N = 24 simulated participants with +/-4-6 mm shifted boundaries, 2 mm CSF sulcal gap, and pial-vein border",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "#f0f4f8"), plot.title = element_text(face = "bold"))
ggsave(file.path(out_dir, "fig1_scenario1_accuracy.png"), fig1, width = 11, height = 6.5, dpi = 180)

s3_long <- melt(
  s3_dt,
  id.vars = "n_voxels",
  measure.vars = c("graph_build_sec", "hrf_fit_sec", "decon_sec"),
  variable.name = "stage",
  value.name = "seconds"
)
stage_labels <- c(
  graph_build_sec = "1. Geodesic CSR Graph Build",
  hrf_fit_sec = "2. Bilateral HRF + Noise Fit",
  decon_sec = "3. Prototype-Quantized Deconvolution"
)
s3_long[, stage_label := factor(stage_labels[as.character(stage)], levels = unname(stage_labels))]

fig2 <- ggplot(s3_long, aes(x = n_voxels, y = seconds, color = stage_label)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.6) +
  labs(
    title = "Linear O(V) Whole-Brain Scaling of Sparse Bilateral HRF Graph",
    subtitle = "3D cortical mask (2 mm voxels, 26-connected geodesic CSR, <= 48 HRF prototypes)",
    x = "Number of In-Mask Cortical Voxels (V)",
    y = "Wall-Clock Runtime (seconds)",
    color = "Pipeline Stage"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave(file.path(out_dir, "fig2_scenario3_scalability.png"), fig2, width = 8.5, height = 5.2, dpi = 180)

saveRDS(
  list(
    scenario1_dt = s1_dt,
    scenario1_summary = s1_summary,
    scenario2_dt = s2_dt,
    scenario2_summary = s2_summary,
    scenario3_scalability = s3_dt
  ),
  file = file.path(out_dir, "sparse_graph_hrf_benchmark.rds")
)
cat(sprintf("\nSaved all benchmark tables, RDS, and figures to %s\n", out_dir))
