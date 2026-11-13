# Fixed-input demonstration of the public parcel panel and its working covariance on live library(medusa).
library(medusa)
set.seed(1104026)

n <- 300L
tr <- 1
rho <- 0.4
sigma <- c(0.013, 0.018, 0.022, 0.032)
spatial_rho <- c("Negative" = -0.2, "Independent" = 0, "Shared" = 0.65, "Nearly common" = 0.95)
groups <- setNames(rep(names(spatial_rho), each = 4L), sprintf("v%02d", seq_len(16L)))
K <- hrf_library(tr, peak_delays = 6, dispersions = 1)
V <- toeplitz(rho^(0:(n - 1L)))
Ltime <- t(chol(V))
bold <- events <- covs <- truth <- list()

event_design <- function(ev, kernels) {
  medusa:::hrf_event_design_cpp(n, tr, 0, ev$onset, ev$duration, seq_len(nrow(ev)), nrow(ev), kernels)
}

for (r in c("run1", "run2")) {
  onsets <- 14 + cumsum(c(0, sample(7:12, 26, replace = TRUE)))
  ev <- data.frame(onset = onsets, condition = "feedback", duration = 0)
  x <- as.numeric(scale(arima.sim(list(ar = 0.65), n = length(onsets))))
  a <- 0.7 + 0.18 * x + rnorm(length(onsets), sd = 0.09)
  Fmat <- event_design(ev, K)[[1L]]
  y <- 100 + as.numeric(Fmat %*% a)
  Y <- vapply(seq_along(groups), function(i) y, numeric(n))
  colnames(Y) <- names(groups)
  for (p in names(spatial_rho)) {
    R <- (1 - spatial_rho[p]) * diag(4) + spatial_rho[p]
    noise <- Ltime %*% matrix(rnorm(n * 4L), n) %*% chol(R * outer(sigma, sigma))
    Y[, groups == p] <- Y[, groups == p] + noise
  }
  bold[[r]] <- Y
  events[[r]] <- ev
  covs[[r]] <- data.frame(value = x, value_squared = x^2)
  truth[[r]] <- a
}

fit <- estimate_crossrun_hrf(
  bold, events, tr, K, noise_rho = rho,
  amplitude_model = "trial", selection_mode = "hard"
)
dec <- deconvolve_crossrun_hrf(
  bold, fit, lambda = 0.03, penalty = "sobolev",
  scaling = "noise_sd", nuisance_mode = "ignore", rwls = "none"
)
for (r in names(bold)) {
  # Supply true marginal SD for exact analytical SE propagation check
  dec$scaling[[r]]$noise_sd <- rep(sigma, length(spatial_rho))
}

epochs <- function(dc = dec, cv = covs, ...) {
  args <- list(
    deconvolution = dc, covariates = cv, target_condition = "feedback",
    lags = 0, parcel_groups = groups, source_model = "impulse",
    restore_scale = TRUE, tau2_mode = "marginal", return_panel = TRUE
  )
  do.call(decon_epoch_regression, modifyList(args, list(...)))
}

out <- epochs(whiten_panel = TRUE)
resid_out <- epochs(tau2_mode = "residual")
bad <- dec
bad$scaling$run2$status[5L] <- "invalid"
bad$activity$run2[, 5L] <- NA_real_
drop_out <- epochs(bad)
unaggregated <- decon_epoch_regression(
  dec, covariates = NULL, target_condition = "feedback", lags = 0,
  source_model = "impulse", restore_scale = TRUE, return_panel = TRUE
)

panel <- as.data.frame(out$panel)
panel_resid <- as.data.frame(resid_out$panel)
panel$truth <- unlist(truth, use.names = FALSE)[
  match(panel$run, names(truth)) * nrow(events$run1) - nrow(events$run1) + panel$trial
]
stopifnot(
  !anyDuplicated(panel[c("run", "event_index", "evt_time", "parcel")]),
  all(is.finite(panel$se_demixed)), all(panel$se_demixed > 0),
  all(is.finite(panel$precision_weight_cohort)), all(panel$precision_weight_cohort > 0)
)

h <- as.numeric(K[, 1L])
D <- medusa:::hrf_inverse_operator(n, h, seq_len(n), rho, matrix(numeric(), n, 0L), dec$settings)
Qobs <- diag(n) - matrix(1 / n, n, n)
metrics <- weights <- transform_check <- list()

for (r in names(events)) {
  r_idx <- match(r, names(events))
  ev <- events[[r]]
  m <- nrow(ev)
  Q <- diag(m) - matrix(1 / m, m, m)
  S <- medusa:::build_epoch_readout_matrix(ev$onset, 0, n, tr, 0, "triangle", NULL)
  G <- event_design(ev, K)[[1L]]
  Gcent <- Qobs %*% G
  DGcent <- D %*% Gcent
  sol <- medusa:::solve_duration_demix_operator(S, D, G, Gcent, "noise_sd", 0.05, 0.05)
  M_act <- (sweep(Q, 1L, sol$eff_gain_vec, "/") + diag(m) - Q) %*% sol$W %*% S
  B <- M_act %*% D
  Bcal <- Q %*% B %*% Qobs / sqrt(1 - 1 / m)
  unit_cov <- Bcal %*% V %*% t(Bcal)
  M_post <- Q %*% (M_act %*% DGcent) %*% Q
  Z_vox_r <- dec$noise_innovations[[r_idx]]

  for (p in names(spatial_rho)) {
    in_p <- groups == p
    u <- out$parcel_weights[[p]]
    y <- panel[panel$run == r & panel$parcel == p, ]
    yall <- panel[panel$parcel == p, ]
    a <- sweep(dec$activity[[r]][, in_p], 2L, dec$scaling[[r]]$scale[in_p], "*")
    d2 <- diff(a, differences = 2L)
    sd_ind <- sqrt(sum((u * sigma)^2))
    sd_cap <- sum(u * sigma)
    sd_true <- sqrt((1 - spatial_rho[p]) * sd_ind^2 + spatial_rho[p] * sd_cap^2)

    # Naive activity second-difference proxy (contaminated by shared task curvature)
    sd_naive_d2 <- min(sd_cap, max(sd_ind, sd_cap * mad(d2 %*% u) / sum(u * apply(d2, 2L, mad))))

    # Live MEDuSA task-free whitened innovation spatial scale
    zeta_g <- as.numeric(Z_vox_r[, in_p, drop = FALSE] %*% (u * sigma))
    sd_medusa <- sqrt(min(sd_cap^2, max(0.05 * sd_ind^2, sum(zeta_g^2))))

    Omega <- sd_medusa^2 * unit_cov
    Omega_cent <- (1 - 1 / m) * Omega
    se_err <- max(abs(y$se_demixed - sqrt(diag(Omega))))
    signal_err <- max(abs(y$decon_demixed - as.numeric(B %*% Qobs %*% bold[[r]][, in_p] %*% u)))
    stopifnot(se_err < 1e-8, signal_err < 1e-7, abs(mean(yall$precision_weight) - 1) < 1e-10)

    s2med <- median(yall$se_demixed^2)
    tau <- max(y$tau2_bio[1L], 0.10 * s2med)
    rel <- tau / (tau + y$se_demixed^2)
    mu <- mean(y$decon_demixed)
    stopifnot(
      max(abs(y$reliability_trial - rel)) < 1e-10,
      max(abs(y$decon_eb - (mu + rel * (y$decon_demixed - mu)))) < 1e-10
    )

    tau_zca <- max(tau, 0.15 * s2med)
    Stot <- Omega_cent + tau_zca * tcrossprod(M_post)
    Phi <- medusa:::compute_zca_whitening_matrix(Stot, eig_floor = 0.08, offdiag_scale = 0.90)
    xf_y <- function(x) mean(x) + as.numeric(Phi %*% (x - mean(x)))
    xf_x <- function(x) mean(x) + as.numeric(Phi %*% (M_post %*% (x - mean(x))))
    stopifnot(
      max(abs(y$value_whitened - xf_x(y$value))) < 1e-9,
      max(abs(y$value_squared_whitened - xf_x(y$value_squared))) < 1e-9,
      max(abs(y$decon_demixed_whitened - xf_y(y$decon_demixed))) < 1e-9,
      max(abs(y$se_whitened - sqrt(diag(Phi %*% Omega %*% t(Phi))))) < 1e-9
    )

    rms_cor <- function(Sigma) sqrt(mean(cov2cor(Sigma)[row(Sigma) != col(Sigma)]^2))
    metrics[[paste(r, p)]] <- data.frame(
      run = r, parcel = p, rho = spatial_rho[p],
      independence = sd_ind, oracle = sd_true, medusa = sd_medusa,
      naive_d2 = sd_naive_d2, ceiling = sd_cap,
      se_check = se_err, signal_check = signal_err, gamma = 0.90,
      rms_before = rms_cor(Stot), rms_after = rms_cor(Phi %*% Stot %*% t(Phi))
    )

    if (p == "Independent" && r == "run2") {
      keep <- 2:4
      u_run <- u[keep] / sum(u[keep])
      dropped <- subset(as.data.frame(drop_out$panel), parcel == p & run == "run2")
      expected <- B %*% Qobs %*% bold[[r]][, which(in_p)[keep]] %*% u_run
      stopifnot(max(abs(dropped$decon_demixed - expected)) < 1e-7)
    }
    if (p == "Independent") {
      weights[[r]] <- data.frame(
        run = r, voxel = paste0("Voxel ", 1:4), sd = sigma,
        nominal = u, realized = if (r == "run2") c(0, u_run) else u
      )
    }
    if (p == "Nearly common" && r == "run1") {
      transform_check <- list(
        cor_before = cov2cor(Stot),
        cor_after = cov2cor(Phi %*% Stot %*% t(Phi)),
        panel = y
      )
    }
  }
}

metrics <- do.call(rbind, metrics)
weights <- do.call(rbind, weights)
voxel_panel <- as.data.frame(unaggregated$panel)
voxel_panel$parcel <- unname(groups[voxel_panel$voxel])

# Mixed-stratum demonstration: 2 voxels with HRF peaks at 4 s and 8 s.
# Live MEDuSA accumulates exact per-stratum operators (\sum_g u_g D_g F_g) rather than a single averaged HRF.
Kmix <- hrf_library(tr, peak_delays = c(4, 8), dispersions = 1)
fwd <- event_design(data.frame(onset = 110, duration = 0), Kmix)
Dmix <- lapply(seq_len(ncol(Kmix)), function(i) {
  medusa:::hrf_inverse_operator(n, as.numeric(Kmix[, i]), seq_len(n), rho, matrix(numeric(), n, 0L), dec$settings)
})
hrep <- rowMeans(Kmix)
hrep <- hrep / sum(hrep)
Drep <- medusa:::hrf_inverse_operator(n, hrep, seq_len(n), rho, matrix(numeric(), n, 0L), dec$settings)
a_real <- (Dmix[[1]] %*% Qobs %*% fwd[[1]] + Dmix[[2]] %*% Qobs %*% fwd[[2]]) / 2
a_rep <- Drep %*% Qobs %*% ((fwd[[1]] + fwd[[2]]) / 2)
lag <- seq_len(n) - 1L - 110
mixed <- data.frame(lag = lag, actual = as.numeric(a_real), representative = as.numeric(a_rep))
mixed_ratio <- sqrt(sum((a_real - a_rep)^2) / sum(a_real^2))

write.csv(panel, "example_panel.csv", row.names = FALSE)
write.csv(metrics, "parcel_noise_check.csv", row.names = FALSE)
write.csv(weights, "voxel_weights.csv", row.names = FALSE)
saveRDS(
  list(
    panel = panel, panel_resid = panel_resid, voxel_panel = voxel_panel,
    metrics = metrics, weights = weights, mixed = mixed, mixed_ratio = mixed_ratio,
    transform_check = transform_check
  ),
  "step4_examples.rds"
)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
print(metrics, digits = 4, row.names = FALSE)
print(weights, digits = 3, row.names = FALSE)
cat("Relative mixed-HRF response difference (exact multi-stratum vs single averaged HRF):", mixed_ratio, "\n")
cat("Verified panel grain, innovation-based spatial SEs, dropped-voxel renormalization, precision normalization, EB and paired ZCA.\n")
