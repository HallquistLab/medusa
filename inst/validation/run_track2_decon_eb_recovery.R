#!/usr/bin/env Rscript
# ==============================================================================
# Benchmark: Track 2 — Hyperparameter Sharing in Deconvolution & Epoch Regression
#
# Validates:
#   Experiment 1 (Track 2a): SNR-Stratified Adaptive lambda_g vs Global Adaptive
#                lambda vs Naive GCV on a heterogeneous 16x16 cortical slice
#                (10% focal biphasic core + 15% moderate periphery + 75% silent
#                spatially smooth AR(2) background; N = 20 subjects).
#   Experiment 2 (Track 2b): Bilateral Graph Empirical-Bayes tau2_bio & Sandwich
#                Scale Moderation across a folded 3D sulcus (CSF gap) + isolated
#                Efron-Morris hotspots under rapid jittered ISI (N = 24 subjects),
#                testing both Step-1 pooling = "voxel" (default NIfTI) and
#                Step-1 pooling = "graph".
#   Experiment 3 (Null Calibration): Spatially smooth 3D AR(3) null volumes
#                (6x6x6, FWHM = 6 mm, N = 32 subjects), verifying subject-level
#                and group-level Type I error calibration and suppression of
#                spurious single-voxel variance-underestimation |t| spikes.
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
  file.path(tempdir(), "track2_eb_results_smoke")
} else if (dir.exists("inst/validation")) {
  "inst/validation/track2_eb_results"
} else {
  file.path(getwd(), "track2_eb_results")
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

hrf_event_design_cpp_fn <- getFromNamespace("hrf_event_design_cpp", "medusa")

make_spatial_smoother <- function(coords_mat, fwhm_mm = 6.0) {
  nv <- nrow(coords_mat)
  if (nv <= 1L || fwhm_mm <= 0) return(diag(nv))
  sigma_sp <- fwhm_mm / sqrt(8.0 * log(2.0))
  dmat <- as.matrix(stats::dist(coords_mat))
  k_sp <- exp(-0.5 * (dmat / sigma_sp)^2)
  sweep(k_sp, 1L, sqrt(rowSums(k_sp^2)), "/")
}

sim_arp_from_innov <- function(Z_mat, ar_coefs = c(0.42, -0.12), sd_vec = 1.0) {
  n_t <- nrow(Z_mat)
  n_v <- ncol(Z_mat)
  if (length(sd_vec) == 1L) sd_vec <- rep(sd_vec, n_v)
  p_ord <- length(ar_coefs)
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
  ir <- numeric(80L)
  ir[1L] <- 1.0
  for (tt in 2:80L) {
    s <- 0.0
    for (lag_k in seq_len(p_ord)) if (tt > lag_k) s <- s + ar_coefs[lag_k] * ir[tt - lag_k]
    ir[tt] <- s
  }
  ar_sd <- sqrt(sum(ir^2))
  sweep(X_full[(n_burn + 1L):(n_burn + n_t), , drop = FALSE] / ar_sd, 2L, sd_vec, "*")
}

cat("======================================================================\n")
cat("EXPERIMENT 1: SNR-Stratified Deconvolution Regularization (Track 2a)\n")
cat("======================================================================\n")

set.seed(20261120L)
tr1 <- 1.0
n_time1 <- 190L
n_trials1 <- 34L
n_subj1 <- if (smoke_mode) 2L else 20L

# 16 x 16 x 1 slice (256 voxels, 3 mm isotropic)
grid1 <- expand.grid(i = 1:16, j = 1:16, k = 1L)
V1 <- nrow(grid1)
vox_ids1 <- paste0("v", seq_len(V1))
coords1 <- 3.0 * as.matrix(grid1)
vm1 <- data.frame(
  voxel = seq_len(V1),
  i = grid1$i, j = grid1$j, k = grid1$k,
  x = coords1[, 1L], y = coords1[, 2L], z = coords1[, 3L]
)
geom1 <- list(dim = c(16L, 16L, 1L), spacing = c(3.0, 3.0, 3.0))
K_sm1 <- make_spatial_smoother(coords1, fwhm_mm = 5.0)

# 25 focal core voxels (5x5), 39 moderate periphery ring, 192 silent background (75%)
dist_center1 <- sqrt((grid1$i - 5.0)^2 + (grid1$j - 5.0)^2)
focal_idx1 <- which(grid1$i >= 3L & grid1$i <= 7L & grid1$j >= 3L & grid1$j <= 7L) # 25
periph_idx1 <- setdiff(which(dist_center1 <= 4.6), focal_idx1)                      # 39
bg_idx1 <- setdiff(seq_len(V1), c(focal_idx1, periph_idx1))                         # 192

lib1 <- hrf_library(tr1, peak_delays = c(5.0, 6.0, 7.0), dispersions = 1.0, undershoot_ratios = 6)
h_true1 <- lib1[, 2L, drop = FALSE]

beta0_true1 <- numeric(V1)
beta2_true1 <- numeric(V1)
beta0_true1[focal_idx1] <- 1.10
beta2_true1[focal_idx1] <- -0.45
beta0_true1[periph_idx1] <- 0.45
beta2_true1[periph_idx1] <- -0.15

exp1_rows <- vector("list", n_subj1 * 3L)
rk1 <- 0L
t0_exp1 <- proc.time()[3L]

for (s in seq_len(n_subj1)) {
  bold_s <- events_s <- covs_s <- vector("list", 2L)
  names(bold_s) <- names(events_s) <- names(covs_s) <- c("run1", "run2")
  neural_focal_r1 <- matrix(0.0, n_time1, length(focal_idx1))

  for (r in 1:2) {
    onsets <- cumsum(c(6.0, stats::runif(n_trials1 - 1L, 3.4, 6.0)))
    x <- as.numeric(scale(stats::rnorm(n_trials1)))
    events_s[[r]] <- data.frame(onset = onsets, duration = 0.0, condition = "stim")
    covs_s[[r]] <- data.frame(x = x)

    Z_sp <- matrix(stats::rnorm(n_time1 * V1), n_time1, V1) %*% t(K_sm1)
    E_r <- sim_arp_from_innov(Z_sp, ar_coefs = c(0.40, -0.10), sd_vec = 0.55)
    Y_r <- E_r

    for (vi in c(focal_idx1, periph_idx1)) {
      s_v <- numeric(n_time1)
      base_amp <- if (vi %in% focal_idx1) 3.2 else 1.3
      b0 <- beta0_true1[vi]
      b2 <- beta2_true1[vi]
      bio_jit <- stats::rnorm(n_trials1, sd = if (vi %in% focal_idx1) 0.50 else 0.25)
      for (ti in seq_len(n_trials1)) {
        k0 <- round(onsets[ti] / tr1) + 1L
        if (k0 >= 1L && k0 <= n_time1) {
          s_v[k0] <- s_v[k0] + (base_amp + b0 * x[ti] + bio_jit[ti])
        }
        if (k0 + 2L <= n_time1) {
          s_v[k0 + 2L] <- s_v[k0 + 2L] + (-0.35 * base_amp + b2 * x[ti])
        }
      }
      if (r == 1L && vi %in% focal_idx1) {
        neural_focal_r1[, match(vi, focal_idx1)] <- s_v
      }
      y_clean <- stats::convolve(s_v, rev(as.numeric(h_true1)), type = "open")[seq_len(n_time1)]
      Y_r[, vi] <- Y_r[, vi] + y_clean
    }
    colnames(Y_r) <- vox_ids1
    attr(Y_r, "voxel_map") <- vm1
    attr(Y_r, "geometry") <- geom1
    bold_s[[r]] <- Y_r
  }

  fit_s <- estimate_crossrun_hrf(bold_s, events_s, lib1, tr = tr1, pooling = "voxel", noise_model = "ar1")

  # 1. Naive whole-map unweighted observation GCV minimum
  y_sc_r1 <- scale(bold_s$run1)
  naive_tun <- deconvolve_reglin(
    BOLDobs = y_sc_r1, kernel = as.numeric(h_true1), lambda = NULL,
    penalty = "sobolev", tune_by = "global", prewhiten_gcv = TRUE,
    prewhiten_rho = stats::median(fit_s$noise$run1$rho), return_diagnostics = TRUE
  )
  lam_naive <- as.numeric(naive_tun$lambda)[1L]
  dec_naive <- deconvolve_crossrun_hrf(bold_s, fit_s, lambda = lam_naive, scaling = "noise_sd", penalty = "sobolev")

  # 2. Global active-weighted adaptive lambda
  dec_glob <- deconvolve_crossrun_hrf(bold_s, fit_s, lambda = NULL, lambda_mode = "global", scaling = "noise_sd", penalty = "sobolev")

  # 3. SNR-stratified adaptive lambda_g (Track 2a)
  dec_strat <- deconvolve_crossrun_hrf(bold_s, fit_s, lambda = NULL, lambda_mode = "auto", scaling = "noise_sd", penalty = "sobolev")

  methods_1 <- list(
    naive_gcv = dec_naive,
    global_adaptive = dec_glob,
    stratified_auto = dec_strat
  )

  for (m_name in names(methods_1)) {
    dec_m <- methods_1[[m_name]]
    bg_sd <- mean(apply(dec_m$activity$run1[, bg_idx1, drop = FALSE], 2L, stats::sd))
    foc_cor <- mean(vapply(seq_along(focal_idx1), function(ii) {
      stats::cor(dec_m$activity$run1[, focal_idx1[ii]], neural_focal_r1[, ii])
    }, numeric(1L)))

    reg_m <- decon_epoch_regression(
      dec_m, covariates = covs_s, events = events_s, target_condition = "stim",
      lags = c(0, 2), restore_scale = TRUE, tau2_shrinkage = "none"
    )
    cf0 <- reg_m$coefficients[reg_m$coefficients$lag == 0, ]
    cf2 <- reg_m$coefficients[reg_m$coefficients$lag == 2, ]
    contrast_est <- cf0$beta_gls[focal_idx1] - cf2$beta_gls[focal_idx1]
    contrast_true <- beta0_true1[focal_idx1] - beta2_true1[focal_idx1]

    lam_sc <- dec_m$scaling$run1$lambda
    rk1 <- rk1 + 1L
    exp1_rows[[rk1]] <- data.table(
      subject = s,
      method = m_name,
      lambda_focal = stats::median(lam_sc[focal_idx1]),
      lambda_bg = stats::median(lam_sc[bg_idx1]),
      n_strata = length(unique(lam_sc)),
      bg_decon_sd = bg_sd,
      focal_neural_cor = foc_cor,
      focal_contrast_mean = mean(contrast_est),
      focal_contrast_rmse = sqrt(mean((contrast_est - contrast_true)^2)),
      focal_beta0_rmse = sqrt(mean((cf0$beta_gls[focal_idx1] - beta0_true1[focal_idx1])^2)),
      bg_beta0_rmse = sqrt(mean((cf0$beta_gls[bg_idx1] - 0.0)^2)),
      whole_map_beta0_rmse = sqrt(mean((cf0$beta_gls - beta0_true1)^2))
    )
  }
}

dt_exp1 <- rbindlist(exp1_rows)
sum_exp1 <- dt_exp1[, .(
  median_lambda_focal = stats::median(lambda_focal),
  median_lambda_bg = stats::median(lambda_bg),
  mean_n_strata = mean(n_strata),
  bg_decon_sd = mean(bg_decon_sd),
  focal_neural_cor = mean(focal_neural_cor),
  focal_contrast_mean = mean(focal_contrast_mean),
  focal_contrast_rmse = mean(focal_contrast_rmse),
  focal_beta0_rmse = mean(focal_beta0_rmse),
  bg_beta0_rmse = mean(bg_beta0_rmse),
  whole_map_beta0_rmse = mean(whole_map_beta0_rmse)
), by = method]

fwrite(dt_exp1, file.path(out_dir, "exp1_snr_stratified_lambda_by_subject.csv"))
fwrite(sum_exp1, file.path(out_dir, "exp1_snr_stratified_lambda_summary.csv"))
cat(sprintf("Experiment 1 completed in %.1f s:\n", proc.time()[3L] - t0_exp1))
print(sum_exp1)

cat("\n======================================================================\n")
cat("EXPERIMENT 2: Bilateral Graph-EB tau2_bio & Scale Moderation (Track 2b)\n")
cat("======================================================================\n")

set.seed(20261121L)
tr2 <- 1.0
n_time2 <- 185L
n_trials2 <- 32L
n_subj2 <- if (smoke_mode) 2L else 24L

# Two 3D cortical banks separated by a 3 mm CSF sulcal gap (i = 4 omitted):
# Bank A (i = 1..3, j = 1..10, k = 1..2): 60 voxels
#   - Medial patch (j <= 5, 30 voxels): high biological trial variance tau2 = 2.60, beta = 0.85
#   - Lateral patch (j > 5, 30 voxels): moderate biological trial variance tau2 = 0.75, beta = 0.55
# Bank B (i = 5..7, j = 1..10, k = 1..2): 60 voxels across sulcal gap
#   - 58 pure-null voxels (tau2 = 0.0, beta = 0.0)
#   - 2 isolated Efron-Morris hotspot voxels (tau2 = 14.0, beta = 1.25)
gA2 <- expand.grid(i = 1:3, j = 1:10, k = 1:2)
gB2 <- expand.grid(i = 5:7, j = 1:10, k = 1:2)
grid2 <- rbind(gA2, gB2)
V2 <- nrow(grid2)
vox_ids2 <- paste0("v", seq_len(V2))
coords2 <- 3.0 * as.matrix(grid2)
vm2 <- data.frame(
  voxel = seq_len(V2),
  i = grid2$i, j = grid2$j, k = grid2$k,
  x = coords2[, 1L], y = coords2[, 2L], z = coords2[, 3L]
)
geom2 <- list(dim = c(7L, 10L, 2L), spacing = c(3.0, 3.0, 3.0))

bankA_med <- which(grid2$i <= 3L & grid2$j <= 5L)   # 30
bankA_lat <- which(grid2$i <= 3L & grid2$j > 5L)    # 30
bankA_all <- c(bankA_med, bankA_lat)                # 60
bankB_all <- which(grid2$i >= 5L)                   # 60
hotspot_idx2 <- bankB_all[c(15L, 45L)]              # 2 isolated hotspots
bankB_null2 <- setdiff(bankB_all, hotspot_idx2)     # 58 pure null

tau2_true2 <- numeric(V2)
tau2_true2[bankA_med] <- 2.60
tau2_true2[bankA_lat] <- 0.75
tau2_true2[hotspot_idx2] <- 14.00

beta_true2 <- numeric(V2)
beta_true2[bankA_med] <- 0.85
beta_true2[bankA_lat] <- 0.55
beta_true2[hotspot_idx2] <- 1.25

# Smooth noise separately within each bank (sulcal CSF boundary in true anatomy)
K_A2 <- make_spatial_smoother(coords2[bankA_all, , drop = FALSE], fwhm_mm = 5.5)
K_B2 <- make_spatial_smoother(coords2[bankB_all, , drop = FALSE], fwhm_mm = 5.5)

lib2 <- hrf_library(tr2, peak_delays = c(5.0, 6.0, 7.0), dispersions = 1.0, undershoot_ratios = 6)
true_k_idx2 <- c(rep(1L, 30L), rep(2L, 30L), rep(3L, 60L))

pipelines2 <- c("voxel_none", "voxel_global_eb", "voxel_graph_eb", "graph_graph_eb")
beta_store2 <- lapply(pipelines2, function(nm) matrix(NA_real_, n_subj2, V2))
names(beta_store2) <- pipelines2
se_store2 <- lapply(pipelines2, function(nm) matrix(NA_real_, n_subj2, V2))
names(se_store2) <- pipelines2
df_store2 <- lapply(pipelines2, function(nm) matrix(NA_real_, n_subj2, V2))
names(df_store2) <- pipelines2
tau_store2 <- lapply(pipelines2, function(nm) matrix(NA_real_, n_subj2, V2))
names(tau_store2) <- pipelines2

t0_exp2 <- proc.time()[3L]
for (s in seq_len(n_subj2)) {
  bold_s <- events_s <- covs_s <- vector("list", 2L)
  names(bold_s) <- names(events_s) <- names(covs_s) <- c("run1", "run2")

  for (r in 1:2) {
    # Highly variable ISI (2.2s to 7.2s) creates strong trial-to-trial heteroscedasticity in se_demixed
    isis <- ifelse(stats::runif(n_trials2 - 1L) < 0.45,
                   stats::runif(n_trials2 - 1L, 2.2, 3.2),
                   stats::runif(n_trials2 - 1L, 4.8, 7.2))
    onsets <- cumsum(c(6.0, isis))
    x <- as.numeric(scale(stats::rnorm(n_trials2)))
    events_s[[r]] <- data.frame(onset = onsets, duration = 0.0, condition = "cue")
    covs_s[[r]] <- data.frame(x = x)

    G_by_k <- lapply(1:3, function(ki) {
      hrf_event_design_cpp_fn(
        n_time2, tr2, 0.0, onsets, rep(0.0, n_trials2), seq_len(n_trials2), n_trials2,
        lib2[, ki, drop = FALSE]
      )[[1L]]
    })

    ZA <- matrix(stats::rnorm(n_time2 * length(bankA_all)), n_time2, length(bankA_all)) %*% t(K_A2)
    ZB <- matrix(stats::rnorm(n_time2 * length(bankB_all)), n_time2, length(bankB_all)) %*% t(K_B2)
    Z_all <- cbind(ZA, ZB)
    # Heterogeneous observation noise SD across voxels (0.35 to 0.65)
    sd_vox <- 0.35 + 0.30 * ((grid2$j - 1.0) / 9.0)
    Y_r <- sim_arp_from_innov(Z_all, ar_coefs = c(0.38, -0.08), sd_vec = sd_vox)

    bio_shared_med <- stats::rnorm(n_trials2, sd = sqrt(0.80))
    bio_shared_lat <- stats::rnorm(n_trials2, sd = sqrt(0.25))
    for (vi in seq_len(V2)) {
      if (tau2_true2[vi] > 0 || beta_true2[vi] != 0) {
        bio_v <- if (vi %in% bankA_med) {
          bio_shared_med + stats::rnorm(n_trials2, sd = sqrt(tau2_true2[vi] - 0.80))
        } else if (vi %in% bankA_lat) {
          bio_shared_lat + stats::rnorm(n_trials2, sd = sqrt(tau2_true2[vi] - 0.25))
        } else {
          stats::rnorm(n_trials2, sd = sqrt(tau2_true2[vi]))
        }
        amps_v <- 2.2 + beta_true2[vi] * x + bio_v
        Y_r[, vi] <- Y_r[, vi] + as.numeric(G_by_k[[true_k_idx2[vi]]] %*% amps_v)
      }
    }
    colnames(Y_r) <- vox_ids2
    attr(Y_r, "voxel_map") <- vm2
    attr(Y_r, "geometry") <- geom2
    bold_s[[r]] <- Y_r
  }

  # Step 1: default NIfTI voxelwise vs sparse bilateral graph
  fit_vox <- estimate_crossrun_hrf(bold_s, events_s, lib2, tr = tr2, pooling = "voxel", noise_model = "ar1")
  dec_vox <- deconvolve_crossrun_hrf(bold_s, fit_vox, lambda = NULL, lambda_mode = "auto", scaling = "noise_sd", penalty = "sobolev")

  fit_grp <- estimate_crossrun_hrf(bold_s, events_s, lib2, tr = tr2, pooling = "graph", noise_model = "ar1")
  dec_grp <- deconvolve_crossrun_hrf(bold_s, fit_grp, lambda = NULL, lambda_mode = "auto", scaling = "noise_sd", penalty = "sobolev")

  r_v_none <- decon_epoch_regression(dec_vox, covs_s, events = events_s, lags = 0, weighting = "gls", restore_scale = TRUE, tau2_shrinkage = "none")
  r_v_glob <- decon_epoch_regression(dec_vox, covs_s, events = events_s, lags = 0, weighting = "gls", restore_scale = TRUE, tau2_shrinkage = "eb")
  r_v_ge   <- decon_epoch_regression(dec_vox, covs_s, events = events_s, lags = 0, weighting = "gls", restore_scale = TRUE, tau2_shrinkage = "auto")
  r_g_ge   <- decon_epoch_regression(dec_grp, covs_s, events = events_s, lags = 0, weighting = "gls", restore_scale = TRUE, tau2_shrinkage = "auto")

  res_list <- list(
    voxel_none = r_v_none,
    voxel_global_eb = r_v_glob,
    voxel_graph_eb = r_v_ge,
    graph_graph_eb = r_g_ge
  )
  for (p_nm in pipelines2) {
    cf <- res_list[[p_nm]]$coefficients
    beta_store2[[p_nm]][s, ] <- cf$beta_gls
    se_store2[[p_nm]][s, ] <- cf$se_gls
    df_store2[[p_nm]][s, ] <- cf$df
    tau_store2[[p_nm]][s, ] <- cf$tau2_bio
  }
}

exp2_sum_rows <- lapply(pipelines2, function(p_nm) {
  B <- beta_store2[[p_nm]]
  S <- se_store2[[p_nm]]
  D <- df_store2[[p_nm]]
  Tmat <- tau_store2[[p_nm]]

  true_B_mat <- matrix(beta_true2, nrow = n_subj2, ncol = V2, byrow = TRUE)
  true_T_mat <- matrix(tau2_true2, nrow = n_subj2, ncol = V2, byrow = TRUE)

  t_crit <- stats::qt(0.975, df = D)
  cov_mat <- abs(B - true_B_mat) <= t_crit * S

  # Group-level random-effects 1-sample t-stat across N = 24 subjects at each voxel
  grp_mean <- colMeans(B)
  grp_se <- apply(B, 2L, stats::sd) / sqrt(n_subj2)
  grp_t <- grp_mean / grp_se
  grp_p <- 2.0 * stats::pt(abs(grp_t), df = n_subj2 - 1L, lower.tail = FALSE)

  emp_sd_vox <- apply(B, 2L, stats::sd)
  mean_se_vox <- colMeans(S)

  data.table(
    pipeline = p_nm,
    tau2_rmse_all = sqrt(mean((Tmat - true_T_mat)^2)),
    tau2_rmse_bankA = sqrt(mean((Tmat[, bankA_all] - true_T_mat[, bankA_all])^2)),
    tau2_zero_rate_bankA = mean(Tmat[, bankA_all] == 0),
    tau2_bleed_bankB_null = mean(Tmat[, bankB_null2]),
    tau2_hotspot_mean = mean(Tmat[, hotspot_idx2]),
    beta_gls_rmse_bankA = sqrt(mean((B[, bankA_all] - true_B_mat[, bankA_all])^2)),
    beta_gls_rmse_all = sqrt(mean((B - true_B_mat)^2)),
    se_to_emp_sd_ratio = stats::median(mean_se_vox / emp_sd_vox),
    ci95_coverage_bankA = mean(cov_mat[, bankA_all]),
    subj_fpr_bankB_null = mean(2.0 * stats::pt(abs(B[, bankB_null2] / S[, bankB_null2]), df = D[, bankB_null2], lower.tail = FALSE) < 0.05),
    group_mean_t_bankA = mean(grp_t[bankA_all]),
    group_power_bankA_p001 = mean(grp_p[bankA_all] < 0.001),
    group_fpr_bankB_null = mean(grp_p[bankB_null2] < 0.05)
  )
})
dt_exp2_sum <- rbindlist(exp2_sum_rows)
fwrite(dt_exp2_sum, file.path(out_dir, "exp2_graph_eb_sulcal_summary.csv"))
cat(sprintf("Experiment 2 completed in %.1f s:\n", proc.time()[3L] - t0_exp2))
print(dt_exp2_sum)

cat("\n======================================================================\n")
cat("EXPERIMENT 3: Spatially Smooth 3D AR(3) Null Calibration (FWHM = 6 mm)\n")
cat("======================================================================\n")

set.seed(20261122L)
tr3 <- 1.0
n_time3 <- 165L
n_trials3 <- 28L
n_subj3 <- if (smoke_mode) 3L else 32L

grid3 <- expand.grid(i = 1:6, j = 1:6, k = 1:6)
V3 <- nrow(grid3) # 216 voxels
vox_ids3 <- paste0("v", seq_len(V3))
coords3 <- 3.0 * as.matrix(grid3)
vm3 <- data.frame(
  voxel = seq_len(V3),
  i = grid3$i, j = grid3$j, k = grid3$k,
  x = coords3[, 1L], y = coords3[, 2L], z = coords3[, 3L]
)
geom3 <- list(dim = c(6L, 6L, 6L), spacing = c(3.0, 3.0, 3.0))
K_sm3 <- make_spatial_smoother(coords3, fwhm_mm = 6.0)
lib3 <- hrf_library(tr3, peak_delays = c(5.0, 6.0, 7.0), dispersions = 1.0, undershoot_ratios = 6)

modes3 <- c("none", "eb", "graph_eb")
p_store3 <- lapply(modes3, function(m) matrix(NA_real_, n_subj3, V3))
names(p_store3) <- modes3
t_store3 <- lapply(modes3, function(m) matrix(NA_real_, n_subj3, V3))
names(t_store3) <- modes3
b_store3 <- lapply(modes3, function(m) matrix(NA_real_, n_subj3, V3))
names(b_store3) <- modes3

t0_exp3 <- proc.time()[3L]
for (s in seq_len(n_subj3)) {
  bold_s <- events_s <- covs_s <- vector("list", 2L)
  names(bold_s) <- names(events_s) <- names(covs_s) <- c("run1", "run2")
  for (r in 1:2) {
    onsets <- cumsum(c(6.0, stats::runif(n_trials3 - 1L, 3.8, 6.2)))
    x <- as.numeric(scale(stats::rnorm(n_trials3)))
    events_s[[r]] <- data.frame(onset = onsets, duration = 0.0, condition = "cue")
    covs_s[[r]] <- data.frame(x = x)
    Z_sp <- matrix(stats::rnorm(n_time3 * V3), n_time3, V3) %*% t(K_sm3)
    Y_r <- sim_arp_from_innov(Z_sp, ar_coefs = c(0.45, -0.15, 0.06), sd_vec = 1.0)
    colnames(Y_r) <- vox_ids3
    attr(Y_r, "voxel_map") <- vm3
    attr(Y_r, "geometry") <- geom3
    bold_s[[r]] <- Y_r
  }

  fit_s <- estimate_crossrun_hrf(bold_s, events_s, lib3, tr = tr3, pooling = "graph", noise_model = "ar1")
  dec_s <- deconvolve_crossrun_hrf(bold_s, fit_s, lambda = NULL, lambda_mode = "auto", scaling = "noise_sd", penalty = "sobolev")

  for (m_nm in modes3) {
    reg_m <- decon_epoch_regression(
      dec_s, covs_s, events = events_s, lags = 0,
      weighting = "gls", restore_scale = TRUE, tau2_shrinkage = m_nm
    )
    cf <- reg_m$coefficients
    p_store3[[m_nm]][s, ] <- cf$p_value
    t_store3[[m_nm]][s, ] <- cf$t_gls
    b_store3[[m_nm]][s, ] <- cf$beta_gls
  }
}

exp3_sum_rows <- lapply(modes3, function(m_nm) {
  P <- p_store3[[m_nm]]
  Tmat <- t_store3[[m_nm]]
  B <- b_store3[[m_nm]]

  grp_t <- colMeans(B) / (apply(B, 2L, stats::sd) / sqrt(n_subj3))
  grp_p <- 2.0 * stats::pt(abs(grp_t), df = n_subj3 - 1L, lower.tail = FALSE)

  data.table(
    tau2_shrinkage = m_nm,
    subj_fpr_05 = mean(P < 0.05),
    subj_fpr_01 = mean(P < 0.01),
    mean_subj_max_abs_t = mean(apply(abs(Tmat), 1L, max)),
    q99_abs_t = as.numeric(stats::quantile(abs(Tmat), 0.99)),
    group_fpr_05 = mean(grp_p < 0.05)
  )
})
dt_exp3_sum <- rbindlist(exp3_sum_rows)
fwrite(dt_exp3_sum, file.path(out_dir, "exp3_smooth_null_calibration.csv"))
cat(sprintf("Experiment 3 completed in %.1f s:\n", proc.time()[3L] - t0_exp3))
print(dt_exp3_sum)

# Diagnostic Figure: Track 2 Summary Panel
p1_df <- melt(
  dt_exp1[, .(subject, method, `Background Decon SD` = bg_decon_sd, `Focal Contrast RMSE` = focal_contrast_rmse)],
  id.vars = c("subject", "method"), variable.name = "metric", value.name = "value"
)
fig1 <- ggplot(p1_df, aes(x = method, y = value, fill = method)) +
  geom_boxplot(alpha = 0.85, width = 0.55, outlier.size = 1.2) +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(
    title = "Track 2a: SNR-Stratified Deconvolution Regularization (16x16 Slice, N = 20)",
    x = "Lambda Selection Mode", y = "Metric Value"
  ) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "fig1_track2a_snr_stratified_lambda.png"), fig1, width = 8.5, height = 4.2, dpi = 160)

tau_prof_df <- rbindlist(lapply(pipelines2, function(p_nm) {
  data.table(
    voxel_index = seq_len(V2),
    pipeline = p_nm,
    mean_tau2_hat = colMeans(tau_store2[[p_nm]]),
    tau2_true = tau2_true2
  )
}))
fig2 <- ggplot(tau_prof_df, aes(x = voxel_index)) +
  geom_line(aes(y = tau2_true), color = "black", linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(y = mean_tau2_hat, color = pipeline), linewidth = 0.75) +
  facet_wrap(~pipeline, ncol = 2L) +
  theme_bw(base_size = 11) +
  labs(
    title = "Track 2b: Spatial Recovery of Biological Trial Variance Across Sulcal Gap & Hotspots",
    subtitle = "Dashed black = ground truth (Bank A medial=2.6, lateral=0.75 | CSF gap | Bank B null=0 + 2 hotspots=14.0)",
    x = "Voxel Index (1..60 Bank A, 61..120 Bank B)", y = "Biological Trial Variance (tau2_bio)"
  ) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "fig2_track2b_tau2_spatial_profile.png"), fig2, width = 9.5, height = 5.5, dpi = 160)

cat("\nSaved all Track 2 benchmark tables and figures to:", out_dir, "\n")
