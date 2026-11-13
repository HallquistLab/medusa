library(medusa)
library(data.table)
probe <- medusa::medusa_fit_mmom(rnorm(20), matrix(rnorm(20), ncol = 1),
                                 array(.01, dim = c(1, 1, 20)))
cat("M-MOM backend:", if (is.null(probe$repair_fit)) "medusa native" else "repair", "\n")
set.seed(6112026)
sims <- 120L
n_subject <- 16L
m <- 60L
noise_sd <- .7
# Unregularized moment sensitivity in the observed trial design. This is a
# diagnostic calculation with a stipulated error cross-moment, not a wrapper mode.
sensitivity_slopes <- function(z, shared) {
  z <- copy(z)
  z[, c("seed", "target") := lapply(.SD, function(v) v - mean(v)), by = run,
    .SDcols = c("seed", "target")]
  X <- cbind(1, z$seed, z$x, z$seed * z$x)
  v <- noise_sd^2 * (1 - 1 / m)
  A <- matrix(0, 4, 4)
  A[2, 2] <- A[4, 4] <- nrow(z) * v
  H <- crossprod(X) - A
  if (min(eigen(H, symmetric = TRUE, only.values = TRUE)$values) <= 0) return(NA_real_)
  rhs <- crossprod(X, z$target)
  rhs[4L] <- rhs[4L] - nrow(z) * shared * noise_sd^2 * (1 - 2 / m)
  solve(H, rhs)[4L]
}
rows <- list()
for (scenario in c("Genuine PPI", "Shared-error null")) {
  beta <- if (scenario == "Genuine PPI") .35 else 0
  corr <- if (scenario == "Shared-error null") .65 else 0
  for (i in seq_len(sims)) {
    z <- CJ(id = seq_len(n_subject), run = 1:2, trial = seq_len(m))
    z[, x := sample(rep(c(-1, 1), each = m / 2)), by = .(id, run)]
    z[, latent_seed := rnorm(.N)]
    z[, latent_target := (.4 * (beta != 0) + (beta + rep(rnorm(n_subject, sd = .11),
      each = 2 * m)) * x) * latent_seed + .2 * x + rnorm(.N, sd = .65)]
    z[, e := rnorm(.N)]
    z[, `:=`(seed = latent_seed + noise_sd * e,
             target = latent_target + noise_sd * (corr * x * e + sqrt(1 - corr^2) * rnorm(.N)))]
    panel <- rbindlist(list(z[, .(id, run, trial, event_index = trial, evt_time = 0,
      parcel = "Seed", decon_demixed = seed, se_demixed = noise_sd, x)],
      z[, .(id, run, trial, event_index = trial, evt_time = 0,
      parcel = "Target", decon_demixed = target, se_demixed = noise_sd, x)]))
    for (method in c("naive_unweighted", "mmoms")) {
      fit <- decon_connectivity_repair(panel, "Seed", "Target", covariates = "x",
        test_ppi = TRUE, orthogonalize_seed_x = FALSE, method = method,
        use_whitened = FALSE, region_col = "parcel")
      f <- fit$group[term == "seed_c:x_c"]
      rows[[length(rows)+1L]] <- f[, .(scenario, iteration = i, truth = beta,
        method = if (method == "mmoms") "M-MOMS" else "Naive",
        estimate, p.value, covers = conf.low <= beta & conf.high >= beta)]
    }
    slopes <- z[, .(slope = sensitivity_slopes(.SD, corr)), by = id]$slope
    t <- t.test(slopes, mu = 0)
    rows[[length(rows)+1L]] <- data.table(scenario, iteration = i, truth = beta,
      method = "Known-error moments", estimate = mean(slopes, na.rm = TRUE),
      p.value = t$p.value, covers = t$conf.int[1] <= beta & t$conf.int[2] >= beta)
  }
}
out <- rbindlist(rows)
stopifnot(nrow(out) == 2 * sims * 3, all(is.finite(out$estimate)))
summary <- out[, .(truth = truth[1], n = .N, mean = mean(estimate),
  bias = mean(estimate - truth), rmse = sqrt(mean((estimate - truth)^2)),
  reject = mean(p.value <= .05), coverage = mean(covers)), by = .(scenario, method)]
fwrite(summary, "ppi_summary.csv")
fwrite(out, "ppi_draws.csv")
print(summary)
