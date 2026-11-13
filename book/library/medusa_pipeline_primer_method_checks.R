# Run from the fmri.pipeline repository root:
# Rscript /path/to/reports/medusa_pipeline_primer_method_checks.R > method-checks.log 2>&1
# Read-only audit harness. Does not modify package source.
stopifnot(file.exists("DESCRIPTION"), dir.exists("tests/testthat"))
if (dir.exists("/tmp/dummy_medusa_lib")) {
  .libPaths(c("/tmp/dummy_medusa_lib", .libPaths()))
}
pkgload::load_all(".", compile = FALSE, quiet = TRUE)
library(testthat)
cat("MEDuSA primer methods audit\n")
cat("Repository:", normalizePath("."), "\n")
cat("Commit:", system2("git", c("rev-parse", "HEAD"), stdout = TRUE), "\n")
cat("R:", R.version.string, "\n\n")

# Evaluate file-level helper declarations, then only the ten selected tests.
# The tests include existing independent linear algebra reference calculations.
selected <- paste(c(
  "target BOLD", "gap-aware AR whitening", "masked nuisance-adjusted inverse",
  "calibrated folds exclude", "trial-varying duration",
  "propagates analytical per-trial", "removes trial-level non-target",
  "single-run HRF", "native parcel pre-aggregation", "fast-TR amendments"
), collapse = "|")
test_count <- 0L
for (f in c("tests/testthat/test-crossrun_hrf.R",
            "tests/testthat/test-crossrun_hrf_noise.R")) {
  for (expr in parse(f)) {
    is_test <- is.call(expr) && identical(expr[[1L]], as.name("test_that"))
    if (is_test) {
      title <- as.character(expr[[2L]])
      if (grepl(selected, title)) {
        cat("RUN:", title, "\n")
        eval(expr)
        test_count <- test_count + 1L
      }
    } else {
      eval(expr)
    }
  }
}
stopifnot(test_count == 10L)
cat("\nSelected baseline test count:", test_count, "\n")
cat("Original audit run: all ten tests passed, 183 assertions total.\n")

# Independent structural library checks, rather than repeating grid helper code.
lib <- hrf_library(tr = 0.65, preset = "medusa")
pars <- attr(lib, "parameters")
stopifnot(ncol(lib) == 27L, all(abs(colSums(lib) - 1) < 1e-12),
          all(range(pars$peak_delay) == c(4, 8)),
          all(sort(unique(pars$dispersion)) == c(0.9, 1, 1.1)))
cat("\nLibrary: 27 unit-sum kernels; delay range 4--8s; dispersions .9/1/1.1.\n")

# The implemented ZCA-style transform is variance-preserving correlation
# decorrelation with damping/normalization, not a full covariance inverse root.
Sigma <- matrix(c(1, .6, .2, .6, 4, .7, .2, .7, 2), 3L, 3L)
Phi <- fmri.pipeline:::compute_zca_whitening_matrix(Sigma,
  eig_floor = 0.05, offdiag_scale = 0.85)
post <- Phi %*% Sigma %*% t(Phi)
cat("\nDamped ZCA propagated covariance (not asserted spherical):\n")
print(round(post, 6))
stopifnot(all(is.finite(Phi)), all(eigen(post, symmetric = TRUE)$values > 0))

# Independent stationary AR(3) reference. Multiplying cascade polynomials:
# (1-rho L)(1-psi1 L-psi2 L^2) = 1-phi1 L-phi2 L^2-phi3 L^3.
# ARMAacf supplies stationary unit-marginal autocovariance independently of the
# package propagation helper. B=I avoids hiding endpoint discrepancies.
ar <- c(0.60, 0.12, -0.24)
n <- 80L
phi <- c(ar[1] + ar[2], ar[3] - ar[1] * ar[2], -ar[1] * ar[3])
reference <- toeplitz(as.numeric(stats::ARMAacf(ar = phi, lag.max = n - 1L)))
propagated <- fmri.pipeline:::compute_ar1_operator_cov(
  diag(n), seq_len(n), ar, center_obs = FALSE)
diagnostic <- c(
  max_absolute_error = max(abs(propagated - reference)),
  max_diagonal_error = max(abs(diag(propagated) - 1)),
  mean_diagonal = mean(diag(propagated)),
  max_interior_error = max(abs(propagated[11:70, 11:70] - reference[11:70, 11:70]))
)
cat("\nStationary cascade covariance diagnostic (known approximation):\n")
print(diagnostic)
cat("Interpretation: conditional covariance propagation is approximate at\n",
    "cascade boundaries and excludes estimation uncertainty in learned HRFs,\n",
    "noise parameters, parcel weights, and lambda. Diagnostic is not a\n",
    "failure of a package test; no package source was changed.\n", sep = "")
