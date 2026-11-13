# Controlled examples generated with live library(medusa).
library(medusa)
set.seed(1103026)

tr <- 1
n <- 360L
tt <- 0:(n - 1L)
rho <- 0.4
K <- hrf_library(tr, peak_delays = 6, dispersions = 1)
h <- as.numeric(K[, 1L])
L <- length(h)
H <- medusa:::decon_convolution_matrix(n, h)[, (L - 1L) + 1:n, drop = FALSE]
ids <- c("Cue_only", "Target", "Noisy_target")
sigma <- c(0.015, 0.015, 0.04)

bold <- stress_bold <- events <- covs <- amp_target <- rt_list <- list()
for (r in 1:2) {
  rn <- paste0("run", r)
  on <- 10 + cumsum(c(0, sample(13:16, 19, TRUE)))
  rt <- sample(1:3, 20, TRUE)
  x <- as.numeric(scale(rt + 0.4 * rnorm(20)))
  feedback <- on + rt
  ev <- rbind(
    data.frame(onset = on, condition = "cue", duration = 0),
    data.frame(onset = feedback, condition = "feedback", duration = 0)
  )
  ev <- ev[order(ev$onset), ]
  rownames(ev) <- NULL
  X <- medusa:::hrf_event_design_cpp(n, tr, 0, ev$onset, ev$duration, seq_len(nrow(ev)), nrow(ev), K)[[1L]]
  ci <- which(ev$condition == "cue")
  fi <- which(ev$condition == "feedback")
  Y <- Ys <- matrix(0, n, 3L, dimnames = list(NULL, ids))
  a_mat <- matrix(0, 20L, 3L, dimnames = list(NULL, ids))
  for (v in 1:3) {
    a <- rep(0, nrow(ev))
    a[ci] <- 1.8
    a[fi] <- if (v == 1L) 0 else 1 + 0.7 * x
    as <- a
    as[ci] <- 1.8 + 0.9 * x
    e <- as.numeric(arima.sim(list(ar = rho), n = n, n.start = 600)) * sqrt(1 - rho^2) * sigma[v]
    Y[, v] <- 100 + as.numeric(X %*% a) + e
    Ys[, v] <- 100 + as.numeric(X %*% as) + e
    a_mat[, v] <- a[fi]
  }
  bold[[rn]] <- Y
  stress_bold[[rn]] <- Ys
  events[[rn]] <- ev
  covs[[rn]] <- data.frame(x = x)
  amp_target[[rn]] <- a_mat
  rt_list[[rn]] <- rt
}

make_decon <- function(y, ev = events) {
  f <- estimate_crossrun_hrf(y, ev, tr, K, noise_rho = rho, amplitude_model = "trial")
  dc <- deconvolve_crossrun_hrf(y, f, lambda = 0.03, penalty = "sobolev", scaling = "noise_sd", nuisance_mode = "ignore", rwls = "none")
  # Supply true marginal SD to isolate propagation algebra from scale estimation uncertainty
  for (r in names(y)) dc$scaling[[r]]$noise_sd <- sigma
  dc
}

message("Deconvolve main and trial-varying-cue examples")
dec <- make_decon(bold)
dec_stress <- make_decon(stress_bold)
lags <- -2:4
epoch <- function(dc = dec, cv = covs, ...) {
  args <- list(
    deconvolution = dc, covariates = cv, events = events, target_condition = "feedback",
    lags = lags, restore_scale = TRUE, tau2_mode = "marginal", return_panel = TRUE,
    return_epochs = TRUE, weighting = "none"
  )
  do.call(decon_epoch_regression, modifyList(args, list(...)))
}

message("Epoch policies")
raw <- epoch(demix_other_events = FALSE, demix_target_trials = FALSE)
target_only <- epoch(demix_other_events = FALSE)
default <- epoch()
joint <- epoch(demix_ridge = 0.02)
within <- epoch(demix_within_trial = TRUE)
agnostic <- epoch(cv = NULL, whiten_panel = TRUE)
stress_default <- epoch(dc = dec_stress)
stress_joint <- epoch(dc = dec_stress, demix_ridge = 0.02)

# Reconstruct the default single-lag composite operator using exact matrix algebra
r <- "run1"
ev <- events[[r]]
ti <- which(ev$condition == "feedback")
oi <- which(ev$condition != "feedback")
m <- length(ti)
D <- medusa:::hrf_inverse_operator(n, h, 1:n, rho, matrix(numeric(), n, 0L), dec$settings)
S <- medusa:::build_epoch_readout_matrix(ev$onset[ti], 0, n, tr, 0, "triangle", NULL)
Sall <- medusa:::build_epoch_readout_matrix(ev$onset, 0, n, tr, 0, "triangle", NULL)
Gunc <- medusa:::hrf_event_design_cpp(n, tr, 0, ev$onset, ev$duration, seq_len(nrow(ev)), nrow(ev), K)[[1L]]
G <- sweep(Gunc, 2L, colMeans(Gunc), "-")
DG <- D %*% G
other <- medusa:::build_other_demix_operator(Sall, DG, ev, oi, ti, 0.05)
Rother <- S %*% DG[, oi, drop = FALSE]
Sclean <- S - Rother %*% other$L_other
sol <- medusa:::solve_duration_demix_operator(S, D, Gunc[, ti, drop = FALSE], G[, ti, drop = FALSE], dec$settings$scaling, 0.05, 0.05)
Q <- diag(m) - matrix(1 / m, m, m)
A_gain <- sweep(Q, 1L, sol$eff_gain_vec, "/") + (diag(m) - Q)
M_act <- A_gain %*% sol$W %*% Sclean
B <- M_act %*% D
Braw <- S %*% D
Btarget <- A_gain %*% sol$W %*% S %*% D
Cobs <- diag(n) - matrix(1 / n, n, n)

y1 <- bold$run1[, 1L] - mean(bold$run1[, 1L])
pub <- subset(as.data.frame(default$panel), run == "run1" & voxel == ids[1L] & evt_time == 0)
public_error <- max(abs(as.numeric(B %*% y1) - pub$decon_demixed))
stopifnot(public_error < 1e-7)

V <- toeplitz(rho^(0:(n - 1L)))
Bcal <- Q %*% B %*% Cobs / sqrt(1 - 1 / m)
Omega <- sigma[1L]^2 * Bcal %*% V %*% t(Bcal)
se_error <- max(abs(sqrt(diag(Omega)) - pub$se_demixed))
stopifnot(se_error < 1e-8)

fact_cov <- sigma[1L]^2 * medusa:::compute_ar1_operator_cov(B, 1:n, rho, center_obs = TRUE)
fact_cov <- Q %*% fact_cov %*% Q / (1 - 1 / m)
cov_error <- max(abs(Omega - fact_cov))
stopifnot(cov_error < 1e-9)

set.seed(1103027)
Z <- matrix(rnorm(n * 2500L), n, 2500L)
E <- sigma[1L] * t(chol(V)) %*% Z
noise_trials <- Bcal %*% E
se_mc <- data.frame(trial = 1:m, predicted = sqrt(diag(Omega)), empirical = apply(noise_trials, 1L, sd))

# Regularized ZCA in the covariate-agnostic panel (including M_post biological covariance)
pag_all <- subset(as.data.frame(agnostic$panel), voxel == ids[1L] & evt_time == 0)
pag <- subset(pag_all, run == "run1")
smed <- median(pag_all$se_demixed^2)
tau_eff <- max(pag$tau2_bio[1L], 0.10 * smed)
tau_zca <- max(tau_eff, 0.15 * smed)
Omega_cent <- (1 - 1 / m) * Omega
M_post <- Q %*% (M_act %*% DG[, ti, drop = FALSE]) %*% Q
Sigma_tot <- Omega_cent + tau_zca * tcrossprod(M_post)
Phi <- medusa:::compute_zca_whitening_matrix(Sigma_tot, eig_floor = 0.08, offdiag_scale = 0.90)
zerr <- max(abs(mean(pag$decon_demixed) + as.numeric(Phi %*% (pag$decon_demixed - mean(pag$decon_demixed))) - pag$decon_demixed_whitened))
stopifnot(zerr < 1e-7)

# Duration example: constant source height, variable time-on-task; same activity, three estimands
message("Variable-duration source example")
set.seed(1103028)
bd <- ed <- cxd <- durs <- list()
true_rate <- 1.8
for (r in 1:2) {
  rn <- paste0("run", r)
  on <- seq(12, by = 12, length.out = 24L)
  du <- runif(24L, 0.4, 2.8)
  ev_d <- data.frame(onset = on, condition = "decision", duration = du)
  Xd <- medusa:::hrf_event_design_cpp(n, tr, 0, on, du, 1:24, 24L, K)[[1L]]
  yd <- 100 + as.numeric(Xd %*% rep(true_rate, 24L)) + as.numeric(arima.sim(list(ar = rho), n = n, n.start = 600)) * 0.008 * sqrt(1 - rho^2)
  bd[[rn]] <- matrix(yd, ncol = 1L, dimnames = list(NULL, "Duration"))
  ed[[rn]] <- ev_d
  cxd[[rn]] <- data.frame(duration = du)
  durs[[rn]] <- du
}
fd <- estimate_crossrun_hrf(bd, ed, tr, K, noise_rho = rho, amplitude_model = "trial")
dd <- deconvolve_crossrun_hrf(bd, fd, lambda = 0.03, penalty = "sobolev", scaling = "noise_sd", nuisance_mode = "ignore", rwls = "none")
for (r in names(dd$scaling)) dd$scaling[[r]]$noise_sd <- 0.008
durfit <- function(model, scale = "rate") {
  decon_epoch_regression(
    dd, cxd, events = ed, target_condition = "decision", lags = 0,
    source_model = model, duration_scale = scale, restore_scale = TRUE,
    weighting = "none", tau2_mode = "marginal", return_panel = TRUE, return_epochs = TRUE
  )
}
duration_impulse <- durfit("impulse")
duration_rate <- durfit("duration")
duration_integral <- durfit("duration", "integral")

readouts <- sapply(c("linear", "boxcar", "triangle", "gaussian"), function(k) {
  as.numeric(medusa:::build_epoch_readout_matrix(20.35, 0, 45, 1, 0, k, if (k == "linear") NULL else 2.5))
})
stopifnot(max(abs(colSums(readouts) - 1)) < 1e-12)

response_raw <- Braw %*% G[, ti, drop = FALSE]
response_target <- Btarget %*% G[, ti, drop = FALSE]
response_default <- B %*% G[, ti, drop = FALSE]
response_cue <- B %*% G[, oi, drop = FALSE]

metrics <- do.call(rbind, lapply(c("raw", "target_only", "default", "joint", "stress_default", "stress_joint"), function(k) {
  res <- get(k)
  cf <- res$coefficients
  cf <- cf[cf$lag == 0, c("voxel", "beta_raw", "beta_corrected", "mean_se_trial")]
  cf$policy <- k
  cf
}))
write.csv(metrics, "policy_metrics.csv", row.names = FALSE)
write.csv(se_mc, "conditional_se_check.csv", row.names = FALSE)
write.csv(as.data.frame(default$panel), "example_panel.csv", row.names = FALSE)
write.csv(
  data.frame(
    model = c("impulse", "duration_rate", "duration_integral"),
    beta = sapply(list(duration_impulse, duration_rate, duration_integral), function(x) x$coefficients$beta_corrected[1L])
  ),
  "duration_metrics.csv", row.names = FALSE
)
saveRDS(
  list(
    tr = tr, n = n, tt = tt, rho = rho, K = K, h = h, H = H, ids = ids, sigma = sigma,
    bold = bold, stress_bold = stress_bold, events = events, covs = covs,
    amp_target = amp_target, rt_list = rt_list, dec = dec, raw = raw,
    target_only = target_only, default = default, joint = joint, within = within,
    agnostic = agnostic, stress_default = stress_default, stress_joint = stress_joint,
    D = D, S = S, Sall = Sall, G = G, Gunc = Gunc, B = B, Braw = Braw, Bcal = Bcal,
    Q = Q, V = V, Omega = Omega, se_mc = se_mc, noise_trials = noise_trials,
    Sigma_tot = Sigma_tot, Phi = Phi, duration_impulse = duration_impulse,
    duration_rate = duration_rate, duration_integral = duration_integral,
    dd = dd, durs = durs, true_rate = true_rate, readouts = readouts,
    response_raw = response_raw, response_target = response_target,
    response_default = response_default, response_cue = response_cue,
    metrics = metrics, public_error = public_error, se_error = se_error,
    cov_error = cov_error, zerr = zerr
  ),
  "step3_examples.rds"
)
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
writeLines(
  c(
    paste("Public readout reproduction error", public_error),
    paste("Public SE reproduction error", se_error),
    paste("Independent dense AR covariance error", cov_error),
    paste("Public regularized-ZCA reproduction error", zerr),
    paste("Maximum readout weight-sum error", max(abs(colSums(readouts) - 1)))
  ),
  "verification.txt"
)
message("Saved Step 3 examples")
