test_that("default pooling stays 'voxel' for atlas-free NIfTI while pooling = 'graph' auto-builds geodesic graph", {
  tmp_dir <- withr::local_tempdir()
  dims <- c(4L, 3L, 2L, 80L)
  nvox <- prod(dims[1:3])
  tr <- 1.0
  lib <- hrf_library(tr, peak_delays = c(4.5, 6.0, 7.5), dispersions = 1.0, undershoot_ratios = 6)
  pars <- attr(lib, "parameters")
  idx_can <- which(pars$peak_delay == 6.0)[1L]
  idx_fast <- which(pars$peak_delay == 4.5)[1L]
  idx_slow <- which(pars$peak_delay == 7.5)[1L]

  mask_arr <- array(1L, dim = dims[1:3])
  f_mask <- file.path(tmp_dir, "mask.nii.gz")
  RNifti::writeNifti(RNifti::asNifti(mask_arr, pixdim = c(1, 2, 2, 2)), f_mask)

  set.seed(8101)
  ev_all <- list()
  f_runs <- character(2L)
  for (r in 1:2) {
    onsets <- seq(6, 66, by = 10)
    ev_r <- data.frame(run = paste0("run", r), onset = onsets, duration = 0, condition = "stim")
    ev_all[[r]] <- ev_r
    des <- medusa:::hrf_event_design_cpp(dims[4L], tr, 0, onsets, rep(0, length(onsets)),
                                         rep(1L, length(onsets)), 1L, lib)
    arr <- array(0.0, dim = dims)
    for (k in seq_len(dims[3L])) {
      for (j in seq_len(dims[2L])) {
        for (i in seq_len(dims[1L])) {
          true_k <- if (i <= 2L) idx_fast else idx_slow
          sig <- as.numeric(des[[true_k]] %*% matrix(1.2, 1L, 1L))
          arr[i, j, k, ] <- 100 + sig + as.numeric(arima.sim(list(ar = 0.35), n = dims[4L])) * 0.25
        }
      }
    }
    hdr <- RNifti::asNifti(arr, pixdim = c(1, 2, 2, 2, tr))
    f_runs[r] <- file.path(tmp_dir, sprintf("run%d.nii.gz", r))
    RNifti::writeNifti(hdr, f_runs[r])
  }
  events_df <- do.call(rbind, ev_all)

  ts <- read_nifti_ts(f_runs, mask = f_mask, run_names = c("run1", "run2"))
  expect_null(ts$parcel_ids)

  # User Choice 1: default pooling for NIfTI without atlas MUST remain "voxel"
  fit_default <- estimate_crossrun_hrf(ts, events = events_df, kernel_grid = lib, fallback_kernel = idx_can)
  expect_identical(fit_default$settings$pooling, "voxel")
  expect_null(fit_default$specification$spatial_graph)

  # Explicit opt-in to pooling = "graph" without an atlas uses 26-connected lattice geodesic graph
  fit_graph <- estimate_crossrun_hrf(
    ts, events = events_df, kernel_grid = lib, fallback_kernel = idx_can,
    pooling = "graph", selection_mode = "soft_average"
  )
  expect_identical(fit_graph$settings$pooling, "graph")
  expect_null(fit_graph$specification$spatial_groups)
  expect_s3_class(fit_graph$specification$spatial_graph, "medusa_spatial_graph")
  expect_identical(fit_graph$specification$spatial_graph$metric, "lattice_geodesic")

  dec_graph <- deconvolve_crossrun_hrf(ts, fit_graph, lambda = 0.25)
  expect_s3_class(dec_graph, "crossrun_decon")
  expect_equal(dim(dec_graph$activity$run1), c(dims[4L], nvox))
})

test_that("26-connected in-mask geodesic graph blocks cross-sulcal CSF mask gaps", {
  # Two gyral banks (i = 1..3 and i = 5..7) separated by an out-of-mask CSF slice at i = 4 (2 mm gap)
  bank_a <- expand.grid(i = 1:3, j = 1:4, k = 1:2)
  bank_b <- expand.grid(i = 5:7, j = 1:4, k = 1:2)
  vm <- rbind(bank_a, bank_b)
  vm$x <- vm$i * 2.0
  vm$y <- vm$j * 2.0
  vm$z <- vm$k * 2.0
  vm$voxel <- sprintf("v%d", seq_len(nrow(vm)))

  g_geo <- build_spatial_graph(
    list(voxel_map = vm, geometry = list(spacing = c(2, 2, 2))),
    spatial_fwhm = 6.0, metric = "geodesic"
  )
  g_euc <- build_spatial_graph(
    list(voxel_map = vm, geometry = list(spacing = c(2, 2, 2))),
    spatial_fwhm = 6.0, metric = "euclidean"
  )

  n_a <- nrow(bank_a)
  cross_edges_geo <- 0L
  cross_edges_euc <- 0L
  for (u in seq_len(n_a)) {
    nb_geo <- g_geo$col_idx[g_geo$row_ptr[u]:(g_geo$row_ptr[u + 1L] - 1L)]
    nb_euc <- g_euc$col_idx[g_euc$row_ptr[u]:(g_euc$row_ptr[u + 1L] - 1L)]
    cross_edges_geo <- cross_edges_geo + sum(nb_geo > n_a)
    cross_edges_euc <- cross_edges_euc + sum(nb_euc > n_a)
  }
  expect_equal(cross_edges_geo, 0L)
  expect_gt(cross_edges_euc, 0L)

  # Simulate Fast (4s) signal on Bank A and pure noise on Bank B:
  # geodesic graph must have zero cross-sulcal bleed onto Bank B border (i = 5)
  set.seed(8102)
  tr <- 1.0
  n_obs <- 140L
  lib <- hrf_library(tr, peak_delays = c(4.0, 6.0, 8.0), dispersions = 1.0, undershoot_ratios = 6)
  pars <- attr(lib, "parameters")
  idx_can <- which(pars$peak_delay == 6.0)[1L]
  idx_fast <- which(pars$peak_delay == 4.0)[1L]

  bold <- events <- list()
  for (r in 1:3) {
    onsets <- seq(6, 126, by = 8)
    ev <- data.frame(onset = onsets, duration = 0, condition = "stim")
    des <- medusa:::hrf_event_design_cpp(n_obs, tr, 0, onsets, rep(0, length(onsets)),
                                         seq_along(onsets), length(onsets), lib)
    amps <- rnorm(length(onsets), mean = 1.1, sd = 0.25)
    ym <- matrix(0.0, n_obs, nrow(vm), dimnames = list(NULL, vm$voxel))
    for (v in seq_len(nrow(vm))) {
      sig <- if (v <= n_a) as.numeric(des[[idx_fast]] %*% amps) else rep(0.0, n_obs)
      ym[, v] <- 100 + sig + rnorm(n_obs, sd = 0.35)
    }
    bold[[paste0("run", r)]] <- ym
    events[[paste0("run", r)]] <- ev
  }

  fit_geo <- estimate_crossrun_hrf(
    bold, events, tr, lib, fallback_kernel = idx_can,
    pooling = "graph", spatial_graph = g_geo,
    amplitude_model = "trial", selection_mode = "soft_average"
  )
  # Bank A recovers Fast HRF; Bank B stays at canonical fallback
  expect_true(mean(fit_geo$selections$run1$kernel_index[1:n_a] == idx_fast) >= 0.95)
  expect_true(all(fit_geo$selections$run1$kernel_index[(n_a + 1L):nrow(vm)] == idx_can))
})

test_that("soft parcel boundary recovers shifted subject boundaries while preserving true edges", {
  set.seed(8103)
  tr <- 1.0
  n_obs <- 180L
  n_trials <- 22L
  lib <- hrf_library(tr, peak_delays = c(4.0, 5.0, 6.0, 7.0, 8.0), dispersions = 1.0, undershoot_ratios = 6)
  pars <- attr(lib, "parameters")
  idx_can <- which(pars$peak_delay == 6.0)[1L]
  idx_fast <- which(pars$peak_delay == 4.0)[1L]
  idx_slow <- which(pars$peak_delay == 8.0)[1L]

  vm <- expand.grid(i = 1:10, j = 1:4, k = 1:2)
  vm$x <- vm$i * 2.0
  vm$y <- vm$j * 2.0
  vm$z <- vm$k * 2.0
  vm$voxel <- sprintf("v%d", seq_len(nrow(vm)))
  nvox <- nrow(vm)

  # Group atlas boundary is at i <= 5 ("P1") vs i >= 6 ("P2"),
  # whereas this participant's true functional/vascular boundary is at i <= 7 (Fast 4s) vs i >= 8 (Slow 8s).
  atlas_groups <- stats::setNames(ifelse(vm$i <= 5L, "P1", "P2"), vm$voxel)
  true_idx <- ifelse(vm$i <= 7L, idx_fast, idx_slow)
  shifted_strip <- which(vm$i %in% c(6L, 7L))
  true_p2_core <- which(vm$i >= 8L)

  bold <- events <- list()
  for (r in 1:3) {
    onsets <- cumsum(c(6, runif(n_trials - 1L, 4.5, 9.5)))
    ev <- data.frame(onset = onsets, duration = 0, condition = "stim")
    des <- medusa:::hrf_event_design_cpp(n_obs, tr, 0, onsets, rep(0, n_trials),
                                         seq_len(n_trials), n_trials, lib)
    amps <- 1.1 + 0.55 * rnorm(n_trials)
    ym <- matrix(0.0, n_obs, nvox, dimnames = list(NULL, vm$voxel))
    for (v in seq_len(nvox)) {
      sd_v <- if (vm$i[v] %in% c(6L, 7L)) 0.35 else 0.25
      sig <- as.numeric(des[[true_idx[v]]] %*% amps)
      ym[, v] <- 100 + sig + as.numeric(arima.sim(list(ar = 0.35), n = n_obs)) * sd_v
    }
    bold[[paste0("run", r)]] <- ym
    events[[paste0("run", r)]] <- ev
  }

  g_soft <- build_spatial_graph(list(voxel_map = vm, geometry = list(spacing = c(2, 2, 2))),
                                spatial_fwhm = 5.5, spatial_groups = atlas_groups, parcel_boundary = "soft")
  fit_parcel <- estimate_crossrun_hrf(
    bold, events, tr, lib, spatial_groups = atlas_groups,
    pooling = "parcel", fallback_kernel = idx_can,
    amplitude_model = "trial", selection_mode = "soft_average"
  )
  fit_graph_soft <- estimate_crossrun_hrf(
    bold, events, tr, lib, spatial_groups = atlas_groups,
    pooling = "graph", spatial_graph = g_soft, parcel_boundary = "soft",
    fallback_kernel = idx_can, amplitude_model = "trial", selection_mode = "soft_average"
  )

  true_k_mat <- lib[, true_idx, drop = FALSE]
  rmse_parcel_shift <- sqrt(mean((fit_parcel$selected_kernels$run1[, shifted_strip] - true_k_mat[, shifted_strip])^2))
  rmse_soft_shift <- sqrt(mean((fit_graph_soft$selected_kernels$run1[, shifted_strip] - true_k_mat[, shifted_strip])^2))
  rmse_soft_p2 <- sqrt(mean((fit_graph_soft$selected_kernels$run1[, true_p2_core] - true_k_mat[, true_p2_core])^2))

  expect_lt(rmse_soft_shift, 0.50 * rmse_parcel_shift)
  expect_lt(rmse_soft_p2, 0.03)
})

test_that("pooling = 'graph' enforces strict target-run exclusion and separates pial-vein noise borders", {
  set.seed(8104)
  tr <- 1.0
  n_obs <- 160L
  lib <- hrf_library(tr, peak_delays = c(4.5, 6.0, 7.5), dispersions = 1.0, undershoot_ratios = 6)
  pars <- attr(lib, "parameters")
  idx_can <- which(pars$peak_delay == 6.0)[1L]
  idx_fast <- which(pars$peak_delay == 4.5)[1L]

  vm <- expand.grid(i = 1:6, j = 1:5, k = 1:2)
  vm$x <- vm$i * 2.0
  vm$y <- vm$j * 2.0
  vm$z <- vm$k * 2.0
  vm$voxel <- sprintf("v%d", seq_len(nrow(vm)))
  nvox <- nrow(vm)
  is_vein <- vm$j == 5L # 20% high-rho, high-variance pial-vein border

  bold <- events <- list()
  for (r in 1:3) {
    onsets <- seq(8, 140, by = 9)
    ev <- data.frame(onset = onsets, duration = 0, condition = "stim")
    des <- medusa:::hrf_event_design_cpp(n_obs, tr, 0, onsets, rep(0, length(onsets)),
                                         seq_along(onsets), length(onsets), lib)
    amps <- rnorm(length(onsets), mean = 1.1, sd = 0.35)
    ym <- matrix(0.0, n_obs, nvox, dimnames = list(NULL, vm$voxel))
    for (v in seq_len(nvox)) {
      rho_v <- if (is_vein[v]) 0.80 else 0.25
      sd_v <- if (is_vein[v]) 2.20 else 0.35
      sig <- if (is_vein[v]) rep(0.0, n_obs) else as.numeric(des[[idx_fast]] %*% amps)
      ym[, v] <- 100 + sig + as.numeric(arima.sim(list(ar = rho_v), n = n_obs)) * sd_v
    }
    bold[[paste0("run", r)]] <- ym
    events[[paste0("run", r)]] <- ev
  }

  g_lat <- build_spatial_graph(list(voxel_map = vm, geometry = list(spacing = c(2, 2, 2))), spatial_fwhm = 5.0)
  fit_base <- estimate_crossrun_hrf(
    bold, events, tr, lib, fallback_kernel = idx_can,
    pooling = "graph", spatial_graph = g_lat,
    amplitude_model = "trial", selection_mode = "soft_average"
  )

  # Bilateral noise pooling separates pial-vein rho (~0.80) from GM rho (~0.25)
  expect_gt(mean(fit_base$noise$run1$rho[is_vein]), mean(fit_base$noise$run1$rho[!is_vein]) + 0.30)

  # Perturb target run1 heavily and verify zero leakage into run1's HRF selection or weights
  bold_pert <- bold
  bold_pert$run1 <- bold_pert$run1 * -5.0 + matrix(rnorm(n_obs * nvox, sd = 25.0), n_obs, nvox,
                                                   dimnames = dimnames(bold$run1))
  fit_pert <- estimate_crossrun_hrf(
    bold_pert, events, tr, lib, fallback_kernel = idx_can,
    pooling = "graph", spatial_graph = g_lat,
    amplitude_model = "trial", selection_mode = "soft_average"
  )

  expect_equal(fit_pert$kernel_weights$run1, fit_base$kernel_weights$run1, tolerance = 0)
  expect_equal(fit_pert$selected_kernels$run1, fit_base$selected_kernels$run1, tolerance = 0)
  expect_identical(fit_pert$selections$run1$kernel_index, fit_base$selections$run1$kernel_index)
  expect_identical(fit_pert$selections$run1$status, fit_base$selections$run1$status)
})

test_that("pooling = 'graph' maintains calibrated Type I control on spatially smooth AR(2) null volumes", {
  set.seed(8105)
  tr <- 1.0
  n_obs <- 140L
  lib <- hrf_library(tr, preset = "medusa")
  idx_can <- 1L

  vm <- expand.grid(i = 1:6, j = 1:6, k = 1:3)
  vm$x <- vm$i * 2.0
  vm$y <- vm$j * 2.0
  vm$z <- vm$k * 2.0
  vm$voxel <- sprintf("v%d", seq_len(nrow(vm)))
  nvox <- nrow(vm)
  coords <- as.matrix(vm[, c("x", "y", "z")])
  dmat <- as.matrix(stats::dist(coords))
  sigma_null <- 6.0 / sqrt(8 * log(2))
  l_sp <- chol(exp(-(dmat^2) / (2 * sigma_null^2)) + diag(1e-6, nvox))

  bold <- events <- list()
  for (r in 1:3) {
    onsets <- seq(8, 128, by = 8)
    events[[paste0("run", r)]] <- data.frame(onset = onsets, duration = 0, condition = "stim")
    innov <- matrix(rnorm(n_obs * nvox), n_obs, nvox) %*% l_sp
    ym <- matrix(0.0, n_obs, nvox, dimnames = list(NULL, vm$voxel))
    for (v in seq_len(nvox)) {
      ym[, v] <- 100 + as.numeric(stats::filter(innov[, v], c(0.45, -0.12), method = "recursive"))
    }
    bold[[paste0("run", r)]] <- ym
  }

  g_lat <- build_spatial_graph(list(voxel_map = vm, geometry = list(spacing = c(2, 2, 2))), spatial_fwhm = 5.0)
  fit_null <- estimate_crossrun_hrf(
    bold, events, tr, lib, fallback_kernel = idx_can,
    pooling = "graph", spatial_graph = g_lat,
    amplitude_model = "trial", selection_mode = "soft_average"
  )

  dep_rate <- mean(fit_null$selections$run1$kernel_index != idx_can)
  non_fb_weight <- mean(1.0 - fit_null$kernel_weights$run1[idx_can, ])
  expect_lte(dep_rate, 0.08)
  expect_lte(non_fb_weight, 0.06)
})

test_that("single-run sub-block mode, custom edge graph, prototype quantization, and decon_epoch_regression work", {
  set.seed(8106)
  tr <- 1.0
  n_obs <- 200L
  lib <- hrf_library(tr, preset = "extended")
  idx_can <- which(attr(lib, "parameters")$peak_delay == 6 & attr(lib, "parameters")$dispersion == 1)[1L]
  idx_fast <- which(attr(lib, "parameters")$peak_delay == 4 & attr(lib, "parameters")$dispersion == 1)[1L]

  nvox <- 36L
  edges <- data.frame(
    from = seq_len(nvox - 1L),
    to = seq.int(2L, nvox),
    dist = rep(2.0, nvox - 1L)
  )
  vox_ids <- sprintf("v%d", seq_len(nvox))
  g_chain <- build_spatial_graph(edges, n_voxels = nvox, voxel_ids = vox_ids, spatial_fwhm = 5.0)
  expect_identical(g_chain$metric, "custom")

  onsets <- cumsum(c(8, runif(21L, 5.5, 9.5)))
  mod_x <- as.numeric(scale(rnorm(length(onsets))))
  ev1 <- data.frame(onset = onsets, duration = 0, condition = "stim")
  des <- medusa:::hrf_event_design_cpp(n_obs, tr, 0, onsets, rep(0, length(onsets)),
                                       seq_along(onsets), length(onsets), lib)
  amps <- 1.0 + 0.65 * mod_x
  y1 <- matrix(0.0, n_obs, nvox, dimnames = list(NULL, vox_ids))
  for (v in seq_len(nvox)) {
    k_use <- if (v <= 18L) idx_fast else idx_can
    y1[, v] <- 100 + as.numeric(des[[k_use]] %*% amps) + rnorm(n_obs, sd = 0.22 + 0.004 * v)
  }

  fit_1r <- estimate_crossrun_hrf(
    list(run1 = y1), list(run1 = ev1), tr, lib, fallback_kernel = idx_can,
    pooling = "graph", spatial_graph = g_chain,
    amplitude_model = "trial", selection_mode = "soft_average",
    max_hrf_strata = 6L
  )
  expect_lte(length(unique(fit_1r$hrf_strata$run1)), 6L)

  dec_1r <- deconvolve_crossrun_hrf(list(run1 = y1), fit_1r)
  reg_1r <- decon_epoch_regression(
    dec_1r,
    covariates = list(run1 = mod_x),
    events = list(run1 = ev1),
    target_condition = "stim",
    lags = c(-1, 0, 1),
    restore_scale = TRUE
  )
  b0 <- reg_1r$summary$mean_beta_corrected[reg_1r$summary$lag == 0]
  r0 <- reg_1r$summary$mean_r_demixed[reg_1r$summary$lag == 0]
  expect_gt(b0, 0.35)
  expect_gt(r0, 0.35)
})

test_that("hrf_safe_svd fallback matches base::svd singular values, subspace, and OLS projection", {
  set.seed(8107)
  x <- cbind(matrix(rnorm(80 * 12), 80, 12), 0.0)
  x[, 13] <- x[, 1] + 1e-13 * rnorm(80)
  y <- matrix(rnorm(80 * 3), 80, 3)

  sv_ref <- base::svd(x)
  # Directly exercise the non-dgesdd polar-polished Gram eigendecomposition branch
  eg <- eigen(crossprod(x), symmetric = TRUE)
  d_all <- sqrt(pmax(eg$values, 0.0))
  v_all <- eg$vectors
  pos <- which(d_all > 1e-12 * d_all[1L])
  u_raw <- x %*% sweep(v_all[, pos, drop = FALSE], 2L, d_all[pos], "/")
  s_u <- eigen(crossprod(u_raw), symmetric = TRUE)
  inv_sqrt <- s_u$vectors %*% sweep(t(s_u$vectors), 1L, 1.0 / sqrt(pmax(s_u$values, 1e-14)), "*")
  u_pol <- u_raw %*% inv_sqrt

  rk <- sum(sv_ref$d > 1e-8 * sv_ref$d[1L])
  expect_equal(d_all[seq_len(rk)], sv_ref$d[seq_len(rk)], tolerance = 1e-10)
  expect_equal(crossprod(u_pol[, seq_len(rk)]), diag(rk), tolerance = 1e-12)
  proj_ref <- sv_ref$u[, seq_len(rk)] %*% crossprod(sv_ref$u[, seq_len(rk)], y)
  proj_pol <- u_pol[, seq_len(rk)] %*% crossprod(u_pol[, seq_len(rk)], y)
  expect_equal(proj_pol, proj_ref, tolerance = 1e-10)
  expect_equal(medusa:::hrf_safe_svd(x)$d[seq_len(rk)], sv_ref$d[seq_len(rk)], tolerance = 1e-12)
})

test_that("build_spatial_graph converts sparse Matrix distance graphs cleanly", {
  skip_if_not_installed("Matrix")
  sm <- Matrix::sparseMatrix(
    i = c(1L, 2L, 2L, 3L),
    j = c(2L, 1L, 3L, 2L),
    x = c(2.0, 2.0, 2.5, 2.5),
    dims = c(4L, 4L),
    dimnames = list(paste0("v", 1:4), paste0("v", 1:4))
  )
  g <- build_spatial_graph(sm, spatial_fwhm = 4.0, max_neighbors = 4L, max_radius = 6.0)
  expect_s3_class(g, "medusa_spatial_graph")
  expect_equal(g$n_voxels, 4L)
  expect_equal(g$voxel_ids, paste0("v", 1:4))
  expect_equal(diff(g$row_ptr), c(2L, 3L, 2L, 1L))
})

test_that("build_spatial_graph supports stats::dist, dense Matrix, square distance matrices, and x_type disambiguation", {
  coords <- as.matrix(expand.grid(x = c(0, 2, 4), y = c(0, 2), z = 0))
  rownames(coords) <- paste0("n", seq_len(nrow(coords)))

  g_coords <- build_spatial_graph(coords, max_neighbors = 8L)
  d_obj <- stats::dist(coords)
  d_mat <- as.matrix(d_obj)

  # 1. stats::dist and base R square distance matrix match coordinate graph (including auto FWHM/radius)
  g_dist <- build_spatial_graph(d_obj, max_neighbors = 8L)
  g_dmat <- build_spatial_graph(d_mat, max_neighbors = 8L)
  expect_identical(g_dist$voxel_ids, rownames(coords))
  expect_equal(g_dist$spatial_fwhm, g_coords$spatial_fwhm, tolerance = 1e-10)
  expect_equal(g_dist$max_radius, g_coords$max_radius, tolerance = 1e-10)
  expect_identical(g_dist$row_ptr, g_coords$row_ptr)
  expect_identical(g_dist$col_idx, g_coords$col_idx)
  expect_equal(g_dist$dist, g_coords$dist, tolerance = 1e-10)
  expect_identical(g_dmat$row_ptr, g_dist$row_ptr)
  expect_identical(g_dmat$col_idx, g_dist$col_idx)
  expect_equal(g_dmat$dist, g_dist$dist, tolerance = 1e-10)

  # 2. Inf entries in geodesic distance matrices act as disconnected barriers
  d_inf <- d_mat
  d_inf[1:3, 4:6] <- Inf
  d_inf[4:6, 1:3] <- Inf
  g_inf <- build_spatial_graph(d_inf, spatial_fwhm = 6.0, max_radius = 10.0)
  nbrs_1 <- g_inf$col_idx[g_inf$row_ptr[1L]:(g_inf$row_ptr[2L] - 1L)]
  expect_true(all(nbrs_1 <= 3L))

  # 3. Dense Matrix (dgeMatrix / dsyMatrix) and symmetric sparse dsCMatrix work without summary() crash
  if (requireNamespace("Matrix", quietly = TRUE)) {
    dm_dense <- Matrix::Matrix(d_mat, sparse = FALSE)
    dm_sym <- Matrix::forceSymmetric(dm_dense)
    g_dense_M <- build_spatial_graph(dm_dense, max_neighbors = 8L)
    g_sym_M <- build_spatial_graph(dm_sym, max_neighbors = 8L)
    expect_identical(g_dense_M$row_ptr, g_dist$row_ptr)
    expect_identical(g_dense_M$col_idx, g_dist$col_idx)
    expect_equal(g_dense_M$dist, g_dist$dist, tolerance = 1e-10)
    expect_identical(g_sym_M$col_idx, g_dist$col_idx)

    sm_sym <- Matrix::forceSymmetric(Matrix::sparseMatrix(
      i = c(1L, 2L), j = c(2L, 3L), x = c(2.0, 2.5), dims = c(4L, 4L)
    ))
    g_sm_sym <- build_spatial_graph(sm_sym)
    expect_equal(diff(g_sm_sym$row_ptr), c(2L, 3L, 2L, 1L))
  }

  # 4. Small 3 x 3 matrices: x_type = 'auto' keeps 3D coordinates, while x_type = 'distance' treats as distances
  c3 <- matrix(c(0, 0, 0, 2, 0, 0, 0, 2, 0), nrow = 3L, byrow = TRUE)
  expect_identical(build_spatial_graph(c3)$metric, "euclidean")
  d3 <- as.matrix(stats::dist(c3))
  g_d3 <- build_spatial_graph(d3, x_type = "distance")
  expect_identical(g_d3$metric, "custom")
  expect_equal(g_d3$n_voxels, 3L)

  # Negative off-diagonal distances fail cleanly
  d_bad <- d3
  d_bad[1L, 2L] <- -1.0
  expect_error(build_spatial_graph(d_bad, x_type = "distance"), "nonnegative")
})
