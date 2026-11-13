# Diagnostic comparison on the Step 1 demonstration using live library(medusa).
library(medusa)
d <- readRDS("../step1_examples.rds")

ar_to_phi <- function(a) {
  rho <- if (length(a) >= 1L && is.finite(a[1])) a[1] else 0
  psi <- numeric(4L)
  if (length(a) >= 2L) {
    k <- min(4L, length(a) - 1L)
    psi[seq_len(k)] <- ifelse(is.finite(a[1L + seq_len(k)]), a[1L + seq_len(k)], 0)
  }
  poly <- stats::convolve(c(1, -rho), rev(c(1, -psi)), type = "open")
  p <- -poly[-1L]
  idx <- which(abs(p) > 1e-12)
  if (!length(idx)) numeric(0L) else p[seq_len(max(idx))]
}
ac <- function(a) {
  p <- ar_to_phi(a)
  if (!length(p) || all(abs(p) < 1e-12)) return(c(1, rep(0, 24)))
  as.numeric(ARMAacf(ar = p, lag.max = 24))
}
truth <- as.numeric(ARMAacf(ar = d$phi, lag.max = 24))

residuals_for <- function(y, ev, h, z) {
  dm <- medusa:::hrf_build_noise_mean_design(d$n, ev, d$tr, 0, h, z, rep(TRUE, d$n))
  des <- dm$design
  s <- svd(sweep(des, 2, pmax(sqrt(colSums(des^2)), 1e-12), "/"), nv = 0)
  rank <- sum(s$d > 1e-8 * s$d[1])
  q <- s$u[, seq_len(rank), drop = FALSE]
  res <- y - q %*% crossprod(q, y)
  attr(res, "q") <- q
  attr(res, "des") <- des
  res
}

ar_cols <- intersect(c("rho", "psi1", "psi2", "psi3", "psi4"), names(d$soft$noise$run1))
out <- curves <- list()
for (r in names(d$bold)) {
  for (g in names(d$true_idx)) {
    message(r, " ", g)
    ix <- which(d$groups == g)
    y <- d$bold[[r]][, ix, drop = FALSE]
    ev <- d$events[[r]]
    z <- cbind(intercept = 1, trend = seq(-1, 1, length.out = d$n), d$nuisance[[r]])
    hsel <- d$soft$selected_kernels[[r]][, ix[1]]
    htrue <- d$K[, d$true_idx[g]]
    fit_noise <- function(h) {
      medusa:::hrf_estimate_noise(
        y, ev, d$tr, 0, matrix(h, ncol = 1L), 1L, z, rep(TRUE, d$n),
        setNames(rep(g, length(ix)), colnames(y)), seq(-0.8, 0.9, 0.1), 2048L
      )
    }
    can_fit <- d$canonical_only$noise[[r]][ix[1], ]
    current_fit <- d$soft$noise[[r]][ix[1], ]
    oracle <- fit_noise(htrue)[1, ]
    methods <- list(
      Canonical_HRF = ac(unlist(can_fit[, ar_cols])),
      Current_MEDuSA = ac(unlist(current_fit[, ar_cols])),
      Oracle_true_HRF = ac(unlist(oracle[, ar_cols]))
    )
    ecan <- residuals_for(y, ev, d$K[, 1], z)
    esel <- residuals_for(y, ev, hsel, z)
    for (M in c(10, 15, 28)) {
      tk <- function(e) {
        a <- sapply(0:(M - 1), function(l) {
          mean(colSums(e[1:(nrow(e) - l), , drop = FALSE] * e[(1 + l):nrow(e), , drop = FALSE]) /
            (nrow(e) - l) / colMeans(e^2))
        })
        a <- 0.5 * (1 + cos(pi * (0:(M - 1)) / M)) * a
        c(a, rep(0, max(0, 25 - length(a))))[1:25]
      }
      methods[[paste0("Tukey", M, "_canonical")]] <- tk(ecan)
      methods[[paste0("Tukey", M, "_selected")]] <- tk(esel)
    }
    for (m in names(methods)) {
      a <- methods[[m]]
      out[[length(out) + 1L]] <- data.frame(
        run = r, parcel = g, method = m,
        acf_RMSE_1_24 = sqrt(mean((a[-1] - truth[-1])^2))
      )
      curves[[length(curves) + 1L]] <- data.frame(
        run = r, parcel = g, method = m, lag = 0:24, estimated = a, true = truth
      )
    }
    write.csv(
      data.frame(
        run = r, parcel = g,
        pass = c("Canonical_HRF", "Current_MEDuSA", "Oracle_true_HRF"),
        rbind(
          can_fit[, c("rho", "psi1", "psi2", "method")],
          current_fit[, c("rho", "psi1", "psi2", "method")],
          oracle[, c("rho", "psi1", "psi2", "method")]
        )
      ),
      paste0("parameters_", r, "_", g, ".csv"),
      row.names = FALSE
    )
  }
}
ans <- do.call(rbind, out)
write.csv(ans, "diagnostic_comparison.csv", row.names = FALSE)
write.csv(do.call(rbind, curves), "diagnostic_acfs.csv", row.names = FALSE)
summary <- aggregate(acf_RMSE_1_24 ~ parcel + method, ans, mean)
write.csv(summary, "diagnostic_summary.csv", row.names = FALSE)
print(summary)
writeLines(capture.output(sessionInfo()), "comparison_session.txt")
