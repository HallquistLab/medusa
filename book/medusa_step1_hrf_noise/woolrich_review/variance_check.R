source("compare_noise.R", echo = FALSE)
# Frozen inverse/readout: isolate covariance and marginal scale estimation on live medusa.
Vtrue <- toeplitz(as.numeric(ARMAacf(ar = d$phi, lag.max = d$n - 1L)))
res <- list()
for (r in names(d$bold)) {
  for (g in names(d$true_idx)) {
    ix <- which(d$groups == g)
    y <- d$bold[[r]][, ix, drop = FALSE]
    ev <- d$events[[r]]
    z <- cbind(1, seq(-1, 1, length.out = d$n), d$nuisance[[r]])
    h <- d$soft$selected_kernels[[r]][, ix[1]]
    pp <- read.csv(paste0("parameters_", r, "_", g, ".csv"))
    ar0 <- unlist(pp[2, c("rho", "psi1", "psi2")], use.names = FALSE)
    D <- medusa:::hrf_inverse_operator(
      d$n, h, 1:d$n, ar0, z,
      list(backend = "r", lambda = 0.01, penalty = "diff2", ridge_floor = 1e-8, tail_padding = 0L)
    )
    S <- matrix(0, nrow(ev), d$n)
    for (i in seq_len(nrow(ev))) {
      a_t <- ev$onset[i] / d$tr
      l <- floor(a_t)
      f <- a_t - l
      S[i, l + 1L] <- 1 - f
      S[i, l + 2L] <- f
    }
    B <- S %*% D
    B <- B - rowMeans(B)
    omega_true <- B %*% Vtrue %*% t(B)
    qtrue <- diag(omega_true)
    for (j in 1:2) {
      a <- unlist(pp[j, c("rho", "psi1", "psi2")], use.names = FALSE)
      vhat <- toeplitz(as.numeric(ARMAacf(ar = ar_to_phi(a), lag.max = d$n - 1L)))
      omega_hat <- B %*% vhat %*% t(B)
      qhat <- diag(omega_hat)
      s_d2 <- medusa:::estimate_decon_noise_sd(y, noise_rho = a[1], noise_psi1 = a[2], noise_psi2 = a[3])
      hmean <- if (j == 1L) d$K[, 1L] else h
      dm_j <- medusa:::hrf_build_noise_mean_design(d$n, ev, d$tr, 0, hmean, z, rep(TRUE, d$n))
      gls_sc <- medusa:::hrf_estimate_gls_noise_sd(
        y_rows = y, x_rows = dm_j$design, rows = seq_len(d$n),
        ar_vec = a, fallback_d2_sd = s_d2
      )
      ratio_d2 <- sqrt(outer(qhat / qtrue, (s_d2 / d$sd_vox[ix])^2))
      ratio_gls <- sqrt(outer(qhat / qtrue, (gls_sc / d$sd_vox[ix])^2))
      cov_err <- function(scales) {
        median(vapply(scales, function(sc) norm(omega_hat * sc^2 - omega_true, "F") / norm(omega_true, "F"), numeric(1)))
      }
      res[[length(res) + 1L]] <- data.frame(
        run = r, parcel = g, method = pp$pass[j],
        median_d2_SD_ratio = median(s_d2 / d$sd_vox[ix]),
        median_gls_SD_ratio = median(gls_sc / d$sd_vox[ix]),
        median_d2_SE_ratio = median(ratio_d2),
        median_gls_SE_ratio = median(ratio_gls),
        median_known_scale_SE_ratio = median(sqrt(qhat / qtrue)),
        median_full_cov_error_d2 = cov_err(s_d2 / d$sd_vox[ix]),
        median_full_cov_error_gls = cov_err(gls_sc / d$sd_vox[ix])
      )
    }
  }
}
out <- do.call(rbind, res)
write.csv(out, "variance_check.csv", row.names = FALSE)
print(aggregate(. ~ parcel + method, out[, setdiff(names(out), "run")], mean))
