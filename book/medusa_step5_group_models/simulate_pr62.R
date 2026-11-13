# Exercise joint multi-lag trajectory inference using live library(medusa).
library(medusa)
library(data.table)

panel <- fread("synthetic_panel.csv")
lags <- sort(unique(panel$evt_time))
L <- length(lags)
C <- diag(L)
R <- toeplitz(0.68^(seq_len(L) - 1L))
window <- list(late_minus_early = list(early = c(0, 1), late = c(3, 4)))
args <- list(
  data = panel, covariates = c("value_within", "entropy_within"),
  lag_resolution = C, lag_error_cor = R, use_whitened = FALSE,
  contrasts = window, n_perm = 4000L, n_boot_latency = 400L, seed = 6205L
)
result <- do.call(medusa_multilag_infer, c(args, list(trajectory_mode = "pspline")))
stopifnot(
  nrow(result$coefficients) == 3L * 2L * L,
  all(result$coefficients$n_subjects == 24L),
  all(result$coefficients$df == 23L),
  all(result$coefficients$p_maxT >= 1 / 4001),
  all(result$coefficients$conf.low.sim <= result$coefficients$conf.low),
  all(result$coefficients$conf.high.sim >= result$coefficients$conf.high)
)

# Audit covariance and contrast propagation (default contrast_source = "unsmoothed").
w <- ifelse(lags %in% 0:1, -0.5, ifelse(lags %in% 3:4, 0.5, 0))
for (p in unique(panel$parcel)) {
  for (tm in args$covariates) {
    c <- result$covariance[[p]][[tm]]
    out <- result$coefficients[parcel == p & term == tm]
    con <- result$contrasts[parcel == p & term == tm]
    plug_var <- diag(c$M_op %*% c$Sigma_raw %*% t(c$M_op))
    stopifnot(
      all(diag(c$Sigma_multilag) + 1e-12 >= plug_var),
      max(abs(c$M_op %*% out$estimate - out$estimate_multilag)) < 1e-9,
      abs(sum(w * out$estimate) - con$estimate) < 1e-9,
      abs(sqrt(drop(crossprod(w, c$Sigma_deblur %*% w))) - con$std.error) < 1e-9,
      isTRUE(all.equal(out$p_fdr_by, p.adjust(out$p.value, "BY"))),
      length(unique(out$omnibus_p)) == 1L
    )
  }
}

# Supply the same subject slopes and their SEs via the coefficient-table route.
coefs <- result$subjects[
  parcel == "Phasic" & term == "value_within",
  .(parcel, id, lag, term, estimate = estimate_raw, std.error = std.error_raw)
]
check <- medusa_multilag_infer(
  coefs, lag_resolution = C, lag_error_cor = R,
  trajectory_mode = "none", n_perm = 200, n_boot_latency = 0
)
observed <- result$coefficients[parcel == "Phasic" & term == "value_within"]
stopifnot(
  max(abs(check$coefficients$estimate - observed$estimate)) < 1e-12,
  all(check$coefficients$estimate == check$coefficients$estimate_multilag)
)

# Verify the named per-subject result route using the same panel and operators.
list_panel <- copy(panel[parcel == "Phasic"])
subject_inputs <- split(list_panel, by = "id", keep.by = FALSE)
subject_inputs <- lapply(subject_inputs, function(subject_panel) {
  list(
    panel = subject_panel, lag_resolution = list(Phasic = C),
    lag_error_cor = list(Phasic = R), settings = list(demix_within_trial = FALSE)
  )
})
list_check <- medusa_multilag_infer(
  subject_inputs, covariates = c("value_within", "entropy_within"), trajectory_mode = "none",
  use_whitened = FALSE, n_perm = 200L, n_boot_latency = 0L, seed = 6207L
)
list_obs <- list_check$coefficients[parcel == "Phasic" & term == "value_within"]
stopifnot(
  max(abs(list_obs$estimate - observed$estimate)) < 1e-12,
  max(abs(list_check$covariance$Phasic$value_within$Sigma_raw -
            result$covariance$Phasic$value_within$Sigma_raw)) < 1e-12
)

# Deterministic known-blur coefficient curve exercises spectral deblurring.
C_blur <- exp(-0.5 * (outer(lags, lags, "-") / 0.64)^2)
theta <- 0.27 * exp(-0.5 * ((lags - 1) / 0.85)^2)
blurred <- data.table(lag = lags, estimate = as.numeric(C_blur %*% theta), std.error = 0.025, df = 23)
blur_fit <- medusa_multilag_infer(
  blurred, lag_resolution = C_blur, lag_error_cor = R, trajectory_mode = "deblur",
  n_perm = 400, n_boot_latency = 0, seed = 6206L
)
stopifnot(max(abs(blur_fit$coefficients$estimate_multilag - blurred$estimate)) > 1e-5)

saveRDS(list(result = result, blur = blur_fit$coefficients, truth = theta), "pr62_examples.rds")
fwrite(result$coefficients, "pr62_coefficients.csv")
fwrite(result$contrasts, "pr62_contrasts.csv")
fwrite(result$latency, "pr62_latency.csv")
cat(
  "Verified panel, named subject-result list, and coefficient-table routes, Marra-Wood smoothing covariance,\n",
  "unsmoothed default contrasts, per-trajectory BY, and simultaneous bands.\n"
)
print(result$latency[, .(parcel, term, peak_lag, peak_lag_grid, omnibus_p, min_p_maxT, min_p_tfce)])
print(result$contrasts[, .(parcel, term, estimate, std.error, p.value)])
