# Illustrative controlled examples generated with live library(medusa).
library(medusa)
set.seed(1102026)

tr <- 1
n <- 180L
tt <- (0:(n - 1L)) * tr
K <- hrf_library(tr, peak_delays = 6, dispersions = 1, undershoot_ratios = 6)
h <- as.numeric(K[, 1L])
L <- length(h)
Hfull <- medusa:::decon_convolution_matrix(n, h)
H <- Hfull[, (L - 1L) + 1:n, drop = FALSE]

ids <- sprintf("v%02d", 1:8)
groups <- setNames(rep(c("Region A", "Region B"), each = 4L), ids)
sigma <- rep(c(0.025, 0.045, 0.075, 0.12), 2L)
rho <- 0.5

bold <- events <- nuisance <- truth <- noises <- list()
for (r in 1:2) {
  rn <- paste0("run", r)
  ev <- data.frame(
    onset = c(20, 24, 40, 44, 65, 69, 86, 90, 110, 114, 135, 139, 155, 159, 174) + if (r == 2) 0.25 else 0,
    condition = "event",
    duration = rep(c(2, 3, 2), 5L)
  )
  neural <- medusa:::hrf_event_design_cpp(n, tr, 0, ev$onset, ev$duration, seq_len(nrow(ev)), nrow(ev), matrix(c(1, 0), ncol = 1L))[[1L]]
  a <- 1 + 0.35 * rnorm(nrow(ev))
  a[c(4, 9, 12)] <- -0.7
  s <- as.numeric(neural %*% a)
  z <- cbind(
    drift1 = cos(pi * (2 * (0:(n - 1L)) + 1) / (2 * n)),
    drift2 = cos(pi * (2 * (0:(n - 1L)) + 1) * 2 / (2 * n))
  )
  Y <- E <- matrix(0, n, 8L, dimnames = list(NULL, ids))
  S <- Y
  for (v in 1:8) {
    gain <- if (v <= 4L) 1 else 0.75
    e <- as.numeric(arima.sim(list(ar = rho), n = n, n.start = 600)) * sqrt(1 - rho^2) * sigma[v]
    S[, v] <- gain * s
    E[, v] <- e
    Y[, v] <- as.numeric(H %*% S[, v]) + 0.12 * z[, 1] + 0.05 * z[, 2] + e
  }
  bold[[rn]] <- Y
  events[[rn]] <- ev
  nuisance[[rn]] <- z
  truth[[rn]] <- S
  noises[[rn]] <- E
}

message("Construct known-HRF/noise Step 1 handoff")
fit <- estimate_crossrun_hrf(
  bold, events, tr, K, nuisance = nuisance, spatial_groups = groups,
  pooling = "parcel", amplitude_model = "trial", noise_rho = rho,
  return_scores = TRUE, selection_mode = "hard"
)

runfit <- function(lambda = NULL, penalty = "sobolev", scaling = "demean",
                   nuisance_mode = "joint", hrf_fit = fit, data = bold,
                   tail_padding = 0L, observation_mask = NULL, rwls = "none") {
  deconvolve_crossrun_hrf(
    data, hrf_fit, lambda = lambda, penalty = penalty, scaling = scaling,
    nuisance_mode = nuisance_mode, backend = "cpp", tail_padding = tail_padding,
    observation_mask = observation_mask, rwls = rwls
  )
}

message("Tune each penalty")
penalties <- setNames(
  lapply(c("ridge", "diff1", "diff2", "matern1", "sobolev"), function(p) runfit(penalty = p)),
  c("ridge", "diff1", "diff2", "matern1", "sobolev")
)
auto <- penalties$sobolev
fixed <- setNames(lapply(c(0.0001, 0.03, 3), function(l) runfit(lambda = l)), c("weak", "moderate", "strong"))

message("Scaling, mask, nuisance, boundary, and HRF sensitivity")
scales <- setNames(
  lapply(c("none", "demean", "zscore", "noise_sd"), function(sc) runfit(lambda = auto$settings$lambda, scaling = sc)),
  c("none", "demean", "zscore", "noise_sd")
)
ignore <- runfit(lambda = auto$settings$lambda, nuisance_mode = "ignore")
padded <- runfit(lambda = auto$settings$lambda, tail_padding = L)
Kbad <- hrf_library(tr, peak_delays = 8, dispersions = 1)
badfit <- estimate_crossrun_hrf(
  bold, events, tr, Kbad, nuisance = nuisance, spatial_groups = groups,
  pooling = "parcel", amplitude_model = "trial", noise_rho = rho, return_scores = TRUE
)
badhrf <- runfit(lambda = auto$settings$lambda, hrf_fit = badfit)

artifact <- bold
artifact$run1[71:81, ] <- artifact$run1[71:81, ] + 1.2
maskfit <- fit
maskfit$specification$valid_rows$run1[71:81] <- FALSE
masked <- runfit(lambda = auto$settings$lambda, hrf_fit = maskfit, data = artifact)
allrows <- runfit(lambda = auto$settings$lambda, hrf_fit = maskfit, data = artifact, observation_mask = "all")

message("Compute frozen inverse and exact dense verification")
rows <- 1:n
z <- cbind(1, seq(-1, 1, length.out = n), nuisance$run1)
operator <- function(result) medusa:::hrf_inverse_operator(n, h, rows, rho, z, result$settings)
ops <- lapply(c(fixed, list(automatic = auto)), operator)
V <- toeplitz(rho^(0:(n - 1L)))
yc <- bold$run1[, 1] - mean(bold$run1[, 1])
op <- ops$automatic
stopifnot(max(abs(op %*% yc - auto$activity$run1[, 1])) < 1e-8)

restore_err <- sapply(c("zscore", "noise_sd"), function(sc) {
  max(abs(sweep(scales[[sc]]$activity$run1, 2, scales[[sc]]$scaling$run1$scale, "*") - scales$demean$activity$run1))
})
stopifnot(max(restore_err) < 1e-8)

rfit <- deconvolve_crossrun_hrf(
  bold, fit, lambda = auto$settings$lambda, penalty = "sobolev",
  scaling = "demean", nuisance_mode = "joint", backend = "r", rwls = "none"
)
backend_err <- max(abs(rfit$activity$run1 - auto$activity$run1))
stopifnot(backend_err < 1e-7)

diag_fit <- diagnose_crossrun_recovery(auto, return_matrices = TRUE)
resolution <- lapply(ops, function(D) D %*% H)
impulse_curves <- lapply(resolution, function(R) R[, 91])
noise_sd_readout <- sapply(ops, function(D) sigma[1] * sqrt(as.numeric(D[91, , drop = FALSE] %*% V %*% t(D[91, , drop = FALSE]))))

# --- Event-informed continuous GLS vs. blind deconvolution ---
message("Event-informed continuous GLS vs blind deconvolution comparison")
set.seed(20261113)
burst_bold <- burst_truth <- list()
for (r in 1:2) {
  rn <- paste0("run", r)
  ev <- events[[rn]]
  neural <- medusa:::hrf_event_design_cpp(n, tr, 0, ev$onset, ev$duration, seq_len(nrow(ev)), nrow(ev), matrix(c(1, 0), ncol = 1L))[[1L]]
  a_b <- rep(0.9, nrow(ev))
  a_b[c(3, 8, 13)] <- -0.6
  a_b[5] <- 3.8 # Genuine heavy-tailed neural burst at onset ~65 s
  s_b <- as.numeric(neural %*% a_b)
  Yb <- Sb <- matrix(0, n, 8L, dimnames = list(NULL, ids))
  for (v in 1:8) {
    # Voxels 1-6 are task-active with varying neural gain; voxels 7-8 are pure noise
    gain <- if (v <= 3L) 1.0 else if (v <= 6L) 0.5 else 0.0
    sig_v <- 0.035
    e <- as.numeric(arima.sim(list(ar = rho), n = n, n.start = 600)) * sqrt(1 - rho^2) * sig_v
    Sb[, v] <- gain * s_b
    Yb[, v] <- as.numeric(H %*% Sb[, v]) + e
  }
  # Inject an isolated sub-scrubbing non-task scanner glitch at t = 124 s (index 125) in run 1
  if (r == 1L) {
    Yb[125L, ] <- Yb[125L, ] + 0.55
  }
  burst_bold[[rn]] <- Yb
  burst_truth[[rn]] <- Sb
}

fit_ev <- estimate_crossrun_hrf(
  burst_bold, events, tr, K, amplitude_model = "trial",
  noise_model = "ar1", evidence = "calibrated"
)
# 1. Event-informed continuous GLS (recommended: noise_sd + task-weighted GCV + rwls = "auto")
dec_event_auto <- deconvolve_crossrun_hrf(
  burst_bold, fit_ev, lambda = NULL, scaling = "noise_sd",
  penalty = "sobolev", nuisance_mode = "ignore", rwls = "auto"
)
# 2. Same continuous GLS without adaptive glitch weighting (rwls = "none")
dec_event_norwls <- deconvolve_crossrun_hrf(
  burst_bold, fit_ev, lambda = dec_event_auto$settings$lambda, scaling = "noise_sd",
  penalty = "sobolev", nuisance_mode = "ignore", rwls = "none"
)
# 3. Blind linear deconvolution (deconvolve_reglin: zscore + raw unweighted GCV across all voxels)
dec_blind_lin <- deconvolve_reglin(
  burst_bold$run1, kernel = h, lambda = NULL, penalty = "sobolev",
  scaling = "zscore", prewhiten_gcv = FALSE
)
# 4. Blind nonlinear Bush & Cisler (2011) sigmoid deconvolution on z-scored BOLD
y1_z <- as.numeric(scale(burst_bold$run1[, 1L]))
dec_blind_bush <- as.numeric(deconvolve_nlreg(matrix(y1_z, ncol = 1L), kernel = h))

# Scale each estimated voxel-1 series onto physical input units for fair visual comparison
s_true_v1 <- burst_truth$run1[, 1L]
s_ev_auto <- dec_event_auto$activity$run1[, 1L] * dec_event_auto$scaling$run1$scale[1L]
s_ev_norwls <- dec_event_norwls$activity$run1[, 1L] * dec_event_norwls$scaling$run1$scale[1L]
s_blind_lin <- dec_blind_lin[, 1L] * sd(burst_bold$run1[, 1L])
# Center and match baseline SD for Bush sigmoid encoding so waveform compression is visible
s_bush_centered <- (dec_blind_bush - median(dec_blind_bush)) * (sd(s_true_v1) / max(sd(dec_blind_bush), 1e-8))

fw_run1 <- dec_event_auto$frame_weights$run1
if (is.null(fw_run1)) fw_run1 <- rep(1.0, n)

blind_vs_event <- list(
  time = tt,
  truth = s_true_v1,
  event_auto = s_ev_auto,
  event_norwls = s_ev_norwls,
  blind_lin = s_blind_lin,
  blind_bush = s_bush_centered,
  frame_weights = as.numeric(fw_run1),
  glitch_time = tt[125L],
  burst_time = events$run1$onset[5L],
  scale_table = data.frame(
    voxel = ids,
    regime = c(rep("High task gain (1.0)", 3), rep("Moderate task gain (0.5)", 3), rep("Null voxel (0.0)", 2)),
    true_noise_sd = 0.035,
    blind_zscore_sd = apply(burst_bold$run1, 2, sd),
    event_gls_noise_sd = dec_event_auto$scaling$run1$scale
  ),
  summary_metrics = data.frame(
    method = c(
      "Event-informed GLS (rwls = 'auto')",
      "Event-informed GLS (rwls = 'none')",
      "Blind linear GCV (zscore)",
      "Blind nonlinear sigmoid (Bush 2011)"
    ),
    cor_with_truth = c(
      cor(s_ev_auto, s_true_v1),
      cor(s_ev_norwls, s_true_v1),
      cor(s_blind_lin, s_true_v1),
      cor(s_bush_centered, s_true_v1)
    ),
    rmse_physical = c(
      sqrt(mean((s_ev_auto - s_true_v1)^2)),
      sqrt(mean((s_ev_norwls - s_true_v1)^2)),
      sqrt(mean((s_blind_lin - s_true_v1)^2)),
      sqrt(mean((s_bush_centered - s_true_v1)^2))
    ),
    burst_peak_recovery = c(
      max(s_ev_auto[64:69]) / max(s_true_v1[64:69]),
      max(s_ev_norwls[64:69]) / max(s_true_v1[64:69]),
      max(s_blind_lin[64:69]) / max(s_true_v1[64:69]),
      max(s_bush_centered[64:69]) / max(s_true_v1[64:69])
    ),
    glitch_ringing_abs_max = c(
      max(abs(s_ev_auto[118:131])),
      max(abs(s_ev_norwls[118:131])),
      max(abs(s_blind_lin[118:131])),
      max(abs(s_bush_centered[118:131]))
    )
  )
)

metrics <- do.call(rbind, lapply(names(penalties), function(p) {
  x <- penalties[[p]]
  s <- x$activity$run1[, 1]
  data.frame(
    penalty = p, lambda = x$settings$lambda,
    correlation = cor(s, truth$run1[, 1]),
    RMSE = sqrt(mean((s - truth$run1[, 1])^2)),
    recon_RMSE = sqrt(mean((as.numeric(H %*% s) - as.numeric(H %*% truth$run1[, 1]))^2))
  )
}))
write.csv(metrics, "penalty_metrics.csv", row.names = FALSE)
write.csv(auto$settings$lambda_tuning, "lambda_tuning.csv", row.names = FALSE)
write.csv(diag_fit$summary, "resolution_summary.csv", row.names = FALSE)
write.csv(scales$noise_sd$scaling$run1, "scaling_noise_sd.csv", row.names = FALSE)
write.csv(blind_vs_event$summary_metrics, "blind_vs_event_metrics.csv", row.names = FALSE)

saveRDS(
  list(
    tr = tr, n = n, tt = tt, K = K, h = h, L = L, H = H, Hfull = Hfull,
    ids = ids, groups = groups, sigma = sigma, rho = rho,
    bold = bold, events = events, nuisance = nuisance, truth = truth, noises = noises,
    fit = fit, auto = auto, penalties = penalties, fixed = fixed, scales = scales,
    ignore = ignore, padded = padded, badhrf = badhrf, artifact = artifact,
    masked = masked, allrows = allrows, ops = ops, V = V, diag_fit = diag_fit,
    resolution = resolution, impulse_curves = impulse_curves, noise_sd_readout = noise_sd_readout,
    metrics = metrics, restore_err = restore_err, backend_err = backend_err,
    blind_vs_event = blind_vs_event
  ),
  "step2_examples.rds"
)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
writeLines(
  c(
    paste("Public inverse max error", max(abs(op %*% yc - auto$activity$run1[, 1]))),
    paste("Scaling restoration max errors", paste(restore_err, collapse = ", ")),
    paste("Backend agreement error", backend_err)
  ),
  "verification.txt"
)
message("Saved Step 2 examples")
