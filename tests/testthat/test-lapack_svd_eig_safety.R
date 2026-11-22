test_that("medusa_svd_cpp and hrf_safe_svd match base::svd without calling La.svd", {
  hrf_safe_svd <- getFromNamespace("hrf_safe_svd", "medusa")

  set.seed(20261122)
  u0 <- qr.Q(qr(matrix(rnorm(36 * 10), 36, 10)))
  v0 <- qr.Q(qr(matrix(rnorm(10 * 10), 10, 10)))
  # Tightly clustered singular values that stress divide-and-conquer deflation
  d0 <- c(rep(1.0 + 1e-13 * seq_len(5), 2))

  cases <- list(
    tall = matrix(rnorm(20 * 6), 20, 6),
    wide = matrix(rnorm(6 * 20), 6, 20),
    square = matrix(rnorm(8 * 8), 8, 8),
    clustered = u0 %*% diag(d0) %*% t(v0),
    zero_mat = matrix(0.0, 9, 4)
  )

  la_calls <- 0L
  trace(base::La.svd, tracer = function() { la_calls <<- la_calls + 1L }, print = FALSE)
  on.exit(try(untrace(base::La.svd), silent = TRUE), add = TRUE)

  for (X in cases) {
    nr <- nrow(X)
    nc <- ncol(X)
    rm <- min(nr, nc)
    for (cb in list(c(rm, rm), c(rm, 0L), c(0L, rm), c(0L, 0L), c(nr, nc), c(2L, 2L))) {
      nu <- cb[1L]
      nv <- cb[2L]
      got <- hrf_safe_svd(X, nu = nu, nv = nv)
      expect_null(dim(got$d))
      expect_length(got$d, rm)
      if (nu > 0L) {
        expect_identical(dim(got$u), c(nr, nu))
        expect_lt(max(abs(crossprod(got$u) - diag(nu))), 1e-11)
      } else {
        expect_null(got$u)
      }
      if (nv > 0L) {
        expect_identical(dim(got$v), c(nc, nv))
        expect_lt(max(abs(crossprod(got$v) - diag(nv))), 1e-11)
      } else {
        expect_null(got$v)
      }
      if (nu == rm && nv == rm && rm > 0L) {
        expect_lt(max(abs(got$u %*% (got$d * t(got$v)) - X)), 1e-10)
      }
    }
  }

  expect_equal(la_calls, 0L)
  expect_error(hrf_safe_svd(matrix(c(1, NA_real_, 2, 3), 2, 2)), "infinite or missing")
  z0 <- hrf_safe_svd(matrix(0.0, 0, 3))
  expect_length(z0$d, 0L)
})

test_that("cross-run graph HRF, rwls Woodbury whitening, and POET epoch regression never invoke La.svd", {
  set.seed(1024)
  n_t <- 80L
  n_v <- 8L
  tr <- 1.0
  coords <- cbind(x = seq_len(n_v) * 2.0, y = 0.0, z = 0.0)
  groups <- rep(1:2, each = 4L)

  mk_events <- function(seed_offset) {
    set.seed(1024 + seed_offset)
    onsets <- cumsum(c(4.0, stats::runif(11L, 2.2, 5.5)))
    data.frame(
      onset = onsets,
      duration = rep(0.5, length(onsets)),
      condition = rep("task", length(onsets)),
      mod = as.numeric(scale(rnorm(length(onsets))))
    )
  }
  events <- list(run1 = mk_events(1L), run2 = mk_events(2L))

  lib <- hrf_library(tr = tr, preset = "standard")
  bold <- lapply(1:2, function(r) {
    Y <- matrix(stats::rnorm(n_t * n_v, sd = 0.5), n_t, n_v)
    colnames(Y) <- paste0("v", seq_len(n_v))
    # Inject outlier spikes to trigger Step-1/2/3 rwls = "auto" Woodbury path
    Y[c(15L, 42L), ] <- Y[c(15L, 42L), ] + 6.0
    Y
  })
  names(bold) <- c("run1", "run2")

  la_calls <- 0L
  trace(base::La.svd, tracer = function() { la_calls <<- la_calls + 1L }, print = FALSE)
  on.exit(try(untrace(base::La.svd), silent = TRUE), add = TRUE)

  hrf_fit <- estimate_crossrun_hrf(
    bold = bold,
    events = events,
    tr = tr,
    kernel_grid = lib,
    pooling = "graph",
    spatial_coords = coords,
    spatial_groups = groups
  )
  dec <- deconvolve_crossrun_hrf(
    bold = bold,
    hrf_fit = hrf_fit,
    rwls = "auto"
  )
  covs <- lapply(events, function(ev) data.frame(mod = ev$mod))
  reg <- decon_epoch_regression(
    deconvolution = dec,
    covariates = covs,
    events = events,
    target_condition = "task",
    lags = 0,
    tau2_shrinkage = "graph_eb",
    aggregate_by_parcel = FALSE
  )

  expect_equal(la_calls, 0L)
  expect_true(is.data.frame(reg$coefficients))
  expect_equal(nrow(reg$coefficients), n_v * 2L)
})
