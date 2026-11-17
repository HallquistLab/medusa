# Helpers for the cross-run interface. Covariance parameters may use their own
# run's observations; shape selection and its evidence use OTHER runs only.

# Drop-in SVD wrapper guarding against LAPACK divide-and-conquer `dgesdd`
# secular-equation non-convergence (`error code 1`) on macOS Accelerate / OpenBLAS
# when multi-trial convolved design matrices have near-degenerate singular values.
hrf_safe_svd <- function(x, nu = min(nrow(x), ncol(x)), nv = min(nrow(x), ncol(x))) {
  res <- tryCatch(base::svd(x, nu = nu, nv = nv), error = function(e) NULL)
  if (!is.null(res)) return(res)

  nr <- nrow(x)
  nc <- ncol(x)
  r_min <- min(nr, nc)
  if (r_min == 0L) {
    return(list(
      d = numeric(0L),
      u = if (nu > 0L) matrix(0.0, nr, 0L) else NULL,
      v = if (nv > 0L) matrix(0.0, nc, 0L) else NULL
    ))
  }
  if (nr >= nc) {
    eg <- eigen(crossprod(x), symmetric = TRUE)
    d_all <- sqrt(pmax(eg$values, 0.0))
    v_all <- eg$vectors
    u_out <- NULL
    if (nu > 0L) {
      nu_use <- min(as.integer(nu), r_min)
      d_max <- max(d_all[1L], .Machine$double.eps)
      pos <- which(d_all > 1e-12 * d_max)
      u_mat <- matrix(0.0, nr, nu_use)
      if (length(pos)) {
        k_p <- min(length(pos), nu_use)
        u_raw <- x %*% sweep(v_all[, seq_len(k_p), drop = FALSE], 2L, d_all[seq_len(k_p)], "/")
        s_u <- eigen(crossprod(u_raw), symmetric = TRUE)
        inv_sqrt <- s_u$vectors %*% sweep(t(s_u$vectors), 1L, 1.0 / sqrt(pmax(s_u$values, 1e-14)), "*")
        u_mat[, seq_len(k_p)] <- u_raw %*% inv_sqrt
      }
      u_out <- u_mat
    }
    v_out <- if (nv > 0L) v_all[, seq_len(min(as.integer(nv), nc)), drop = FALSE] else NULL
    list(d = d_all[seq_len(r_min)], u = u_out, v = v_out)
  } else {
    tr_svd <- Recall(t(x), nu = nv, nv = nu)
    list(d = tr_svd$d, u = tr_svd$v, v = tr_svd$u)
  }
}

# Convert cascaded AR parameters c(rho1, psi_1, ..., psi_q) representing
# (1 - rho1 L)(1 - sum_{j=1}^q psi_j L^j) = 1 - sum_{k=1}^p phi_k L^k
# into stationary causal AR(p) coefficients, unit-marginal ACF, and exact
# innovation-to-marginal variance ratio sigma_u^2 / sigma_e^2.
hrf_ar_cascade_info <- function(ar_vec, max_lag = 16L) {
  max_lag <- max(1L, as.integer(max_lag))
  if (!length(ar_vec)) {
    return(list(phi = numeric(0), acf = c(1.0, rep(0.0, max_lag)), sig_u2 = 1.0))
  }
  rho1 <- if (is.finite(ar_vec[1L])) min(max(ar_vec[1L], -0.95), 0.95) else 0.0
  psi <- if (length(ar_vec) >= 2L) as.numeric(ar_vec[-1L]) else numeric(0)
  psi[!is.finite(psi)] <- 0.0
  while (length(psi) > 0L && abs(psi[length(psi)]) <= 1e-10) {
    psi <- psi[-length(psi)]
  }
  if (!length(psi)) {
    if (abs(rho1) <= 1e-10) {
      return(list(phi = numeric(0), acf = c(1.0, rep(0.0, max_lag)), sig_u2 = 1.0))
    }
    acf_v <- rho1^(0:max_lag)
    return(list(phi = rho1, acf = acf_v, sig_u2 = max(1e-6, 1.0 - rho1 * rho1)))
  }
  # Enforce strict causal stationarity of the residual psi polynomial
  roots_psi <- tryCatch(Mod(polyroot(c(1.0, -psi))), error = function(e) rep(0.0, length(psi)))
  if (any(roots_psi <= 1.02)) {
    rad <- 0.97 / max(1e-6, max(1.0 / pmax(roots_psi, 1e-6)))
    psi <- psi * (rad^seq_along(psi))
  }
  q <- length(psi)
  phi <- numeric(q + 1L)
  phi[1L] <- rho1 + psi[1L]
  if (q >= 2L) {
    for (k in 2L:q) {
      phi[k] <- psi[k] - rho1 * psi[k - 1L]
    }
  }
  phi[q + 1L] <- -rho1 * psi[q]
  while (length(phi) > 0L && abs(phi[length(phi)]) <= 1e-12) {
    phi <- phi[-length(phi)]
  }
  p <- length(phi)
  if (p == 0L) {
    return(list(phi = numeric(0), acf = c(1.0, rep(0.0, max_lag)), sig_u2 = 1.0))
  }
  roots_phi <- tryCatch(Mod(polyroot(c(1.0, -phi))), error = function(e) rep(0.0, p))
  if (any(roots_phi <= 1.015)) {
    rad <- 0.975 / max(1e-6, max(1.0 / pmax(roots_phi, 1e-6)))
    phi <- phi * (rad^seq_len(p))
  }
  lag_need <- max(max_lag, p)
  acf_v <- as.numeric(stats::ARMAacf(ar = phi, lag.max = lag_need))
  sig_u2 <- max(1e-6, 1.0 - sum(phi * acf_v[2L:(p + 1L)]))
  list(phi = phi, acf = acf_v[seq_len(max_lag + 1L)], sig_u2 = sig_u2)
}

# Build the nuisance + free-trial + condition-derivative + DF-adaptive DCT (< 1/128 Hz)
# mean design matrix for Step 1 noise estimation and GLS residual scale estimation.
hrf_build_dct_drift_cols <- function(n_obs, tr, n_mask, highpass_sec = 128.0) {
  if (!is.finite(highpass_sec) || highpass_sec <= 0 ||
      n_obs * tr < 0.75 * highpass_sec || n_mask < 40L) {
    return(matrix(numeric(), n_obs, 0L))
  }
  cap_cos <- max(6L, min(24L, as.integer(floor(0.15 * n_mask))))
  n_cos <- max(2L, min(cap_cos, as.integer(floor(2.0 * n_obs * tr / highpass_sec))))
  t_idx <- seq(0L, n_obs - 1L)
  vapply(seq_len(n_cos), function(k) {
    cos(pi * (2.0 * t_idx + 1.0) * k / (2.0 * n_obs))
  }, numeric(n_obs))
}

hrf_combine_run_nuisance <- function(nuisance_r = NULL, noise_pcs_r = NULL, n_obs = NULL) {
  has_nuis <- !is.null(nuisance_r) && is.matrix(nuisance_r) && ncol(nuisance_r) > 0L
  has_pcs  <- !is.null(noise_pcs_r) && is.matrix(noise_pcs_r) && ncol(noise_pcs_r) > 0L
  if (has_nuis && has_pcs) {
    cbind(nuisance_r, noise_pcs_r)
  } else if (has_nuis) {
    nuisance_r
  } else if (has_pcs) {
    noise_pcs_r
  } else if (!is.null(nuisance_r) && is.matrix(nuisance_r)) {
    nuisance_r
  } else {
    matrix(numeric(), if (is.null(n_obs)) 0L else as.integer(n_obs), 0L)
  }
}

# Identify task-unrelated background cortex ("noise pool") from frozen out-of-sample
# cross-run task informativeness (`selections`), or validate a user-supplied `noise_pool`
# (logical mask, integer indices, or voxel IDs). Protects active parcels against edge
# PSF bleed and requires adequate spatial breadth in automatic mode so noise-pool PCs
# never induce small-pool Global-Signal-Regression (GSR) negative correlations.
hrf_identify_noise_pool <- function(bold, selections, spatial_groups = NULL,
                                    spatial_graph = NULL, spatial_coords = NULL,
                                    spatial_fwhm = NULL, parcel_cross_weight = 0.35,
                                    voxel_map = NULL, geometry = NULL,
                                    noise_pool = NULL, mode = c("auto", "fixed"),
                                    min_auto_voxels = 24L, min_auto_parcels = 8L) {
  mode <- match.arg(mode)
  n_runs <- length(bold)
  y1 <- bold[[1L]]
  n_vox <- ncol(y1)
  voxel_ids <- colnames(y1)
  if (is.null(voxel_ids)) voxel_ids <- paste0("v", seq_len(n_vox))

  empty_pool <- list(
    pool_idx = integer(0L),
    pool_voxels = character(0L),
    eval_idx = integer(0L),
    act_idx = integer(0L),
    explicit_pool = !is.null(noise_pool)
  )

  # Ensure candidate columns have finite, non-constant observations across runs
  ok_all <- rep(TRUE, n_vox)
  for (r in seq_len(n_runs)) {
    yr <- bold[[r]]
    st_r <- decon_signal_status(yr)
    ok_all <- ok_all & (st_r == "ok")
  }
  if (sum(ok_all) < 4L) return(empty_pool)

  # User-supplied explicit noise pool (e.g., WM/CSF or external background voxels)
  if (!is.null(noise_pool)) {
    if (is.logical(noise_pool)) {
      if (length(noise_pool) != n_vox || anyNA(noise_pool)) {
        stop("'noise_pool' logical mask must have length ncol(bold[[1]]) with no NAs.", call. = FALSE)
      }
      cand_bg <- which(noise_pool & ok_all)
    } else if (is.numeric(noise_pool)) {
      idx_in <- unique(as.integer(noise_pool))
      if (anyNA(idx_in) || any(idx_in < 1L | idx_in > n_vox)) {
        stop("'noise_pool' integer indices must lie within 1..ncol(bold[[1]]).", call. = FALSE)
      }
      cand_bg <- idx_in[ok_all[idx_in]]
    } else if (is.character(noise_pool)) {
      if (!all(noise_pool %in% voxel_ids)) {
        stop("'noise_pool' character IDs must match column names of 'bold'.", call. = FALSE)
      }
      idx_in <- unique(match(noise_pool, voxel_ids))
      cand_bg <- idx_in[ok_all[idx_in]]
    } else {
      stop("'noise_pool' must be NULL, a logical mask, integer indices, or character voxel IDs.", call. = FALSE)
    }
    if (length(cand_bg) < 4L) return(empty_pool)
    act_idx <- setdiff(which(ok_all), cand_bg)
    # Stratified 75% / 25% split of cand_bg for bi-cross-validation when mode == "auto"
    if (length(cand_bg) >= 8L) {
      ev_pos <- seq(4L, length(cand_bg), by = 4L)
      eval_idx <- cand_bg[ev_pos]
      pool_sub <- cand_bg[-ev_pos]
    } else {
      eval_idx <- act_idx
      pool_sub <- cand_bg
    }
    return(list(
      pool_idx = cand_bg,
      pool_cv_idx = pool_sub,
      pool_voxels = voxel_ids[cand_bg],
      eval_idx = eval_idx,
      act_idx = act_idx,
      explicit_pool = TRUE
    ))
  }

  # Automatic noise-pool construction from cross-run selections (Step 2 bg_mask foundation)
  if (is.null(selections) || !length(selections)) return(empty_pool)
  if (mode == "auto" && n_vox < min_auto_voxels) return(empty_pool)

  r2_mat <- matrix(0.0, nrow = n_vox, ncol = length(selections))
  r2_sel_mat <- matrix(0.0, nrow = n_vox, ncol = length(selections))
  tp_mat <- matrix(1.0, nrow = n_vox, ncol = length(selections))
  sup_any <- rep(FALSE, n_vox)
  for (r in seq_along(selections)) {
    sel_r <- selections[[r]]
    if (is.data.frame(sel_r) && nrow(sel_r) == n_vox) {
      if ("status" %in% names(sel_r)) {
        st_r <- as.character(sel_r$status)
        sup_any <- sup_any | (!is.na(st_r) & nzchar(st_r) & st_r != "fallback_prespecified")
      }
      if ("best_adjusted_partial_r2" %in% names(sel_r)) {
        v_r2 <- as.numeric(sel_r$best_adjusted_partial_r2)
        v_r2[!is.finite(v_r2)] <- 0.0
        r2_mat[, r] <- pmax(0.0, v_r2)
        r2_sel_mat[, r] <- r2_mat[, r]
      }
      if ("selected_adjusted_partial_r2" %in% names(sel_r)) {
        v_sr2 <- as.numeric(sel_r$selected_adjusted_partial_r2)
        v_sr2[!is.finite(v_sr2)] <- 0.0
        r2_sel_mat[, r] <- pmax(0.0, v_sr2)
      }
      if ("selected_task_p" %in% names(sel_r)) {
        v_tp <- as.numeric(sel_r$selected_task_p)
        v_tp[!is.finite(v_tp)] <- 1.0
        tp_mat[, r] <- pmin(1.0, pmax(0.0, v_tp))
      }
    }
  }

  r2_oos <- rowMeans(0.65 * r2_sel_mat + 0.35 * r2_mat, na.rm = TRUE)
  r2_oos[!is.finite(r2_oos)] <- 0.0
  r2_max <- apply(r2_mat, 1L, max, na.rm = TRUE)
  r2_max[!is.finite(r2_max)] <- 0.0

  sg_use <- NULL
  if (inherits(spatial_graph, "medusa_spatial_graph") && spatial_graph$n_voxels == n_vox) {
    sg_use <- spatial_graph
  } else if (!is.null(voxel_map) && is.data.frame(voxel_map) && nrow(voxel_map) == n_vox) {
    sg_use <- tryCatch(
      build_spatial_graph(
        list(voxel_map = voxel_map, geometry = geometry),
        spatial_groups = spatial_groups,
        parcel_boundary = "soft",
        voxel_ids = voxel_ids
      ),
      error = function(e) NULL
    )
  } else if (!is.null(spatial_coords)) {
    sg_use <- tryCatch(
      build_spatial_graph(
        spatial_coords,
        spatial_fwhm = spatial_fwhm,
        spatial_groups = spatial_groups,
        parcel_boundary = "soft",
        voxel_ids = voxel_ids
      ),
      error = function(e) NULL
    )
  }

  if (!is.null(sg_use)) {
    pos_oos <- r2_oos[r2_oos > 0]
    rng_bw <- if (length(pos_oos) >= 4L) {
      max(0.025, 1.5 * stats::mad(pos_oos, constant = 1.4826))
    } else {
      0.035
    }
    g_ids_bg <- hrf_encode_group_ids(spatial_groups, n_vox)
    r2_stab <- as.numeric(smooth_scalar_csr_cpp(
      sg_use$row_ptr, sg_use$col_idx, sg_use$dist, sg_use$sigma,
      r2_oos, rng_bw, g_ids_bg, as.numeric(parcel_cross_weight)
    ))
    r2_stab[!is.finite(r2_stab)] <- r2_oos[!is.finite(r2_stab)]
  } else if (!is.null(spatial_groups) && length(spatial_groups) == n_vox) {
    g_chr <- as.character(spatial_groups)
    g_means <- tapply(r2_oos, g_chr, mean, na.rm = TRUE)
    r2_stab <- 0.50 * r2_oos + 0.50 * as.numeric(g_means[g_chr])
  } else {
    r2_stab <- r2_oos
  }

  p_comb_v <- stats::pchisq(
    -2.0 * rowSums(log(pmin(pmax(tp_mat, 1e-16), 1.0))),
    df = 2.0 * ncol(tp_mat),
    lower.tail = FALSE
  )
  p_comb_v[!is.finite(p_comb_v)] <- 1.0
  q_comb_v <- stats::p.adjust(p_comb_v, method = "BH")

  bg_mask <- ok_all & !sup_any & (q_comb_v > 0.20) & (p_comb_v > 0.12) & (r2_stab < 0.030) & (r2_max < 0.055)

  # Protect entire parcels that show parcel-pooled task support so edge PSF voxels
  # of active parcels never leak task signal into the background noise pool.
  if (!is.null(spatial_groups) && length(spatial_groups) == n_vox) {
    g_chr <- as.character(spatial_groups)
    for (gg in unique(g_chr)) {
      idx_g <- which(g_chr == gg)
      p_g <- hrf_cauchy_p(p_comb_v[idx_g])
      r2_g <- mean(r2_oos[idx_g], na.rm = TRUE)
      if (any(sup_any[idx_g]) || (is.finite(p_g) && p_g <= 0.10) || (is.finite(r2_g) && r2_g >= 0.020)) {
        bg_mask[idx_g] <- FALSE
      }
    }
  }

  cand_bg <- which(bg_mask)
  act_idx <- which(ok_all & !bg_mask & (sup_any | (q_comb_v < 0.15) | (r2_stab >= 0.020)))

  if (mode == "auto") {
    if (length(cand_bg) < min_auto_voxels) return(empty_pool)
    if (!is.null(spatial_groups) && length(spatial_groups) == n_vox) {
      n_bg_parcels <- length(unique(as.character(spatial_groups[cand_bg])))
      if (n_bg_parcels < min_auto_parcels) return(empty_pool)
    }
  } else {
    if (length(cand_bg) < 6L) return(empty_pool)
  }

  # Disjoint 75% / 25% spatial split of cand_bg (parcel-respecting when parcellated)
  # for Owen-Perry bi-cross-validation of K*
  if (!is.null(spatial_groups) && length(spatial_groups) == n_vox) {
    g_bg <- unique(as.character(spatial_groups[cand_bg]))
    if (length(g_bg) >= 4L) {
      ev_g <- g_bg[seq(4L, length(g_bg), by = 4L)]
      eval_idx <- cand_bg[as.character(spatial_groups[cand_bg]) %in% ev_g]
      pool_sub <- setdiff(cand_bg, eval_idx)
    } else {
      ev_pos <- seq(4L, length(cand_bg), by = 4L)
      eval_idx <- cand_bg[ev_pos]
      pool_sub <- cand_bg[-ev_pos]
    }
  } else {
    ev_pos <- seq(4L, length(cand_bg), by = 4L)
    eval_idx <- cand_bg[ev_pos]
    pool_sub <- cand_bg[-ev_pos]
  }

  list(
    pool_idx = cand_bg,
    pool_cv_idx = pool_sub,
    pool_voxels = voxel_ids[cand_bg],
    eval_idx = eval_idx,
    act_idx = act_idx,
    explicit_pool = FALSE
  )
}

# GLMdenoise-style cross-validated noise-pool temporal PC extraction:
#   1. Residualize noise-pool voxels against intercept, linear trend, supplied confounds,
#      and DCT (< 1/128 Hz) drift per run.
#   2. Prewhiten by median background AR(p) parameters and apply the Gavish-Donoho (2014)
#      optimal Marchenko-Pastur singular-value gate so spatially independent AR(p) noise
#      cleanly selects K* = 0.
#   3. Select K* in {0, ..., K_spiked} via df-adjusted AR-whitened residual variance
#      reduction on disjoint held-out background voxels (Owen-Perry bi-CV with Kay et al.
#      2013 95%-of-peak elbow rule) plus a Beta task-subspace orthogonality guard.
hrf_extract_cv_noise_pcs <- function(bold, events, tr, offsets, kernels, fallback_kernel,
                                     nuisance, valid_rows, selections, noise = NULL,
                                     spatial_groups = NULL, spatial_graph = NULL,
                                     spatial_coords = NULL, spatial_fwhm = NULL,
                                     parcel_cross_weight = 0.35,
                                     voxel_map = NULL, geometry = NULL,
                                     highpass_sec = 128.0,
                                     noise_pcs = "auto", max_noise_pcs = 10L,
                                     noise_pool = NULL) {
  n_runs <- length(bold)
  empty_res <- list(noise_pcs = NULL, k_selected = 0L, pool_voxels = character(0L), cv_mse = NULL)

  if (is.null(noise_pcs) || isFALSE(noise_pcs) || identical(noise_pcs, "none")) {
    return(empty_res)
  }
  is_auto <- isTRUE(noise_pcs) || identical(noise_pcs, "auto")
  if (!is_auto) {
    checkmate::assert_integerish(noise_pcs, lower = 0L, upper = 50L, len = 1L, any.missing = FALSE)
    k_fixed <- as.integer(noise_pcs)
    if (k_fixed == 0L) return(empty_res)
    max_k <- k_fixed
  } else {
    checkmate::assert_integerish(max_noise_pcs, lower = 0L, upper = 50L, len = 1L, any.missing = FALSE)
    max_k <- as.integer(max_noise_pcs)
    if (max_k == 0L) return(empty_res)
  }

  pool_info <- hrf_identify_noise_pool(
    bold = bold,
    selections = selections,
    spatial_groups = spatial_groups,
    spatial_graph = spatial_graph,
    spatial_coords = spatial_coords,
    spatial_fwhm = spatial_fwhm,
    parcel_cross_weight = parcel_cross_weight,
    voxel_map = voxel_map,
    geometry = geometry,
    noise_pool = noise_pool,
    mode = if (is_auto) "auto" else "fixed"
  )
  pool_idx <- pool_info$pool_idx
  if (length(pool_idx) < 4L) return(empty_res)

  pool_cv_idx <- if (!is.null(pool_info$pool_cv_idx) && length(pool_info$pool_cv_idx) >= 4L) {
    pool_info$pool_cv_idx
  } else {
    pool_idx
  }
  eval_idx <- pool_info$eval_idx

  # Helper: standardize columns robustly and mildly winsorize >5.5 sigma spikes so a single
  # glitch TR (already handled by rWLS) cannot masquerade as a persistent physiological PC.
  prep_bg_block <- function(ym_res) {
    sc <- apply(ym_res, 2L, stats::mad, center = 0.0, constant = 1.4826)
    bad <- !is.finite(sc) | sc <= 1e-8
    if (any(bad)) {
      sc2 <- sqrt(colMeans(ym_res^2))
      sc[bad] <- ifelse(is.finite(sc2[bad]) & sc2[bad] > 1e-8, sc2[bad], 1.0)
    }
    ys <- sweep(ym_res, 2L, sc, "/")
    pmin(pmax(ys, -5.5), 5.5)
  }

  get_bg_ar_vec <- function(r, idx_sub) {
    if (is.null(noise) || length(noise) < r || is.null(noise[[r]])) return(0.0)
    nr <- noise[[r]]
    if (!("rho" %in% names(nr))) return(0.0)
    val <- stats::median(nr[["rho"]][idx_sub], na.rm = TRUE)
    if (!is.finite(val)) return(0.0)
    # Prewhiten only broadband first-order AR(1) autocorrelation (never resonant psi1..psi4
    # notch poles, which would notch out narrowband physiological harmonics before SVD!).
    min(0.80, max(-0.80, as.numeric(val)))
  }

  run_cache <- vector("list", n_runs)
  k_cap_runs <- max_k
  for (r in seq_len(n_runs)) {
    yr <- bold[[r]]
    n_obs <- nrow(yr)
    mask_r <- if (!is.null(valid_rows) && length(valid_rows) >= r) {
      !is.na(valid_rows[[r]]) & valid_rows[[r]]
    } else {
      rep(TRUE, n_obs)
    }
    rows_r <- which(mask_r)
    m_r <- length(rows_r)
    if (m_r < 24L) return(empty_res)

    nuis_r <- if (!is.null(nuisance) && length(nuisance) >= r && !is.null(nuisance[[r]])) {
      nuisance[[r]]
    } else {
      matrix(numeric(), n_obs, 0L)
    }
    z_base <- cbind(intercept = 1.0, trend = seq(-1.0, 1.0, length.out = n_obs), nuis_r)
    dct_r <- hrf_build_dct_drift_cols(n_obs, tr, m_r, highpass_sec = highpass_sec)
    if (ncol(dct_r) > 0L) z_base <- cbind(z_base, dct_r)
    z_obs <- z_base[rows_r, , drop = FALSE]

    ar_bg <- get_bg_ar_vec(r, pool_idx)
    wz <- hrf_apply_whitening(z_obs, rows_r, ar_bg, unit_marginal = FALSE)
    nrm_wz <- sqrt(colSums(wz^2)); nrm_wz[nrm_wz == 0] <- 1.0
    sv_wz <- hrf_safe_svd(sweep(wz, 2L, nrm_wz, "/"), nv = 0L)
    rk_wz <- if (length(sv_wz$d) && sv_wz$d[1L] > 0) sum(sv_wz$d > 1e-8 * sv_wz$d[1L]) else 0L
    qz_w <- if (rk_wz > 0L) sv_wz$u[, seq_len(rk_wz), drop = FALSE] else matrix(0.0, m_r, 0L)
    df_r <- m_r - rk_wz
    if (df_r < 16L) return(empty_res)

    # Unwhitened Z orthonormal basis for mapping spatial noise loadings back to observation time
    nrm_zo <- sqrt(colSums(z_obs^2)); nrm_zo[nrm_zo == 0] <- 1.0
    sv_zo <- hrf_safe_svd(sweep(z_obs, 2L, nrm_zo, "/"), nv = 0L)
    rk_zo <- if (length(sv_zo$d) && sv_zo$d[1L] > 0) sum(sv_zo$d > 1e-8 * sv_zo$d[1L]) else 0L
    qz_o <- if (rk_zo > 0L) sv_zo$u[, seq_len(rk_zo), drop = FALSE] else matrix(0.0, m_r, 0L)

    k_cap_r <- min(
      max_k,
      length(pool_cv_idx) - 1L,
      max(0L, as.integer(floor(0.15 * df_r)))
    )
    k_cap_runs <- min(k_cap_runs, k_cap_r)

    # Build AR-whitened task subspace Q_task_w for task-orthogonality safeguard
    q_task_w <- matrix(0.0, m_r, 0L)
    ev_r <- if (!is.null(events) && length(events) >= r) events[[r]] else NULL
    if (!is.null(ev_r) && is.data.frame(ev_r) && nrow(ev_r) >= 1L && !is.null(kernels)) {
      lbl_r <- as.character(ev_r$condition)
      u_lbl <- unique(lbl_r)
      dur_r <- if (is.null(ev_r$duration)) rep(0.0, nrow(ev_r)) else as.numeric(ev_r$duration)
      off_r <- if (length(offsets) >= r) offsets[r] else 0.0
      k_fb <- kernels[, min(max(1L, as.integer(fallback_kernel)), ncol(kernels)), drop = FALSE]
      x_cond <- hrf_event_design_cpp(n_obs, tr, off_r, ev_r$onset, dur_r,
        match(lbl_r, u_lbl), length(u_lbl), k_fb)[[1L]][rows_r, , drop = FALSE]
      wx_c <- hrf_apply_whitening(x_cond, rows_r, ar_bg, unit_marginal = FALSE)
      if (rk_wz > 0L) wx_c <- wx_c - qz_w %*% crossprod(qz_w, wx_c)
      nx_c <- sqrt(colSums(wx_c^2)); pos_xc <- nx_c > 1e-8
      if (any(pos_xc)) {
        sv_xc <- hrf_safe_svd(sweep(wx_c[, pos_xc, drop = FALSE], 2L, nx_c[pos_xc], "/"), nv = 0L)
        rk_xc <- if (length(sv_xc$d) && sv_xc$d[1L] > 0) sum(sv_xc$d > 1e-6 * sv_xc$d[1L]) else 0L
        if (rk_xc > 0L) q_task_w <- sv_xc$u[, seq_len(rk_xc), drop = FALSE]
      }
    }

    # Subsample at most 1024 evenly strided pool voxels for fast O(m^2) SVD at whole-brain scale
    sub_full <- if (length(pool_idx) > 1024L) {
      pool_idx[unique(as.integer(round(seq(1, length(pool_idx), length.out = 1024L))))]
    } else {
      pool_idx
    }
    sub_cv <- if (length(pool_cv_idx) > 1024L) {
      pool_cv_idx[unique(as.integer(round(seq(1, length(pool_cv_idx), length.out = 1024L))))]
    } else {
      pool_cv_idx
    }

    eval_white_block <- function(yr_mat, rws, ar_v, qzw_mat, cols_use) {
      yb <- yr_mat[rws, cols_use, drop = FALSE]
      wb <- hrf_apply_whitening(yb, rws, ar_v, unit_marginal = FALSE)
      if (ncol(qzw_mat) > 0L) wb <- wb - qzw_mat %*% crossprod(qzw_mat, wb)
      prep_bg_block(wb)
    }

    w_cv <- eval_white_block(yr, rows_r, ar_bg, qz_w, sub_cv)
    sv_cv <- hrf_safe_svd(w_cv, nu = min(m_r, ncol(w_cv)), nv = min(max_k, ncol(w_cv)))
    w_eval <- if (length(eval_idx) >= 1L) eval_white_block(yr, rows_r, ar_bg, qz_w, eval_idx) else NULL

    run_cache[[r]] <- list(
      n_obs = n_obs,
      rows_r = rows_r,
      m_r = m_r,
      df_r = df_r,
      ar_bg = ar_bg,
      qz_w = qz_w,
      qz_o = qz_o,
      q_task_w = q_task_w,
      sub_full = sub_full,
      sv_cv = sv_cv,
      w_eval = w_eval
    )
  }

  if (k_cap_runs <= 0L) return(empty_res)

  if (!is_auto) {
    k_star <- k_cap_runs
    cv_mse_vec <- NULL
  } else {
    # Step A: Gavish-Donoho (2014) optimal singular-value threshold in AR-whitened space
    k_gd_runs <- integer(n_runs)
    for (r in seq_len(n_runs)) {
      rc <- run_cache[[r]]
      d_all <- rc$sv_cv$d
      n_nz <- min(rc$df_r, ncol(rc$sv_cv$v))
      if (length(d_all) >= 4L) {
        d_use <- d_all[seq_len(min(length(d_all), rc$df_r))]
        beta_mp <- min(rc$df_r, length(d_use)) / max(rc$df_r, ncol(rc$w_eval) + length(d_use), 1.0)
        # Exact Gavish-Donoho aspect ratio of the m_eff x V_cv matrix:
        v_cv_n <- nrow(rc$sv_cv$v)
        beta_mp <- min(rc$df_r, v_cv_n) / max(rc$df_r, v_cv_n)
        omega_mp <- 0.56 * (beta_mp^3) - 0.95 * (beta_mp^2) + 1.82 * beta_mp + 1.43
        gd_thr <- 1.02 * omega_mp * stats::median(d_use)
        k_gd_runs[r] <- min(k_cap_runs, sum(d_use > gd_thr))
      }
    }
    # Require factor spike support in at least half of runs (or the single run)
    k_spiked <- as.integer(floor(stats::median(k_gd_runs)))
    if (k_spiked <= 0L) return(empty_res)

    # Step B: Task-subspace orthogonality safeguard across runs (Beta null upper tail)
    k_safe <- 0L
    for (kk in seq_len(k_spiked)) {
      p_align_runs <- numeric(0L)
      for (r in seq_len(n_runs)) {
        rc <- run_cache[[r]]
        qt <- rc$q_task_w
        if (ncol(qt) > 0L && ncol(rc$sv_cv$u) >= kk && (rc$df_r - ncol(qt)) >= 6L) {
          uk <- rc$sv_cv$u[, kk]
          eta_k <- min(0.9999, max(0.0, sum(crossprod(qt, uk)^2)))
          p_b <- stats::pbeta(
            eta_k,
            shape1 = 0.5 * ncol(qt),
            shape2 = 0.5 * (rc$df_r - ncol(qt)),
            lower.tail = FALSE
          )
          p_align_runs <- c(p_align_runs, p_b)
        }
      }
      if (length(p_align_runs) > 0L) {
        p_fish <- stats::pchisq(
          -2.0 * sum(log(pmax(p_align_runs, 1e-16))),
          df = 2.0 * length(p_align_runs),
          lower.tail = FALSE
        )
        # Stop before any PC that aligns anomalously strongly with the task design
        if (is.finite(p_fish) && p_fish < 0.005) break
      }
      k_safe <- kk
    }
    if (k_safe <= 0L) return(empty_res)

    # Step C: Disjoint spatial bi-cross-validation MSE(k) on held-out background voxels
    has_eval <- all(vapply(run_cache, function(rc) !is.null(rc$w_eval) && ncol(rc$w_eval) >= 2L, logical(1L)))
    if (!has_eval) {
      k_star <- k_safe
      cv_mse_vec <- NULL
    } else {
      cv_mse_vec <- numeric(k_safe + 1L)
      for (kk in 0:k_safe) {
        mse_r <- numeric(n_runs)
        for (r in seq_len(n_runs)) {
          rc <- run_cache[[r]]
          we <- rc$w_eval
          if (kk > 0L) {
            uk <- rc$sv_cv$u[, seq_len(kk), drop = FALSE]
            we <- we - uk %*% crossprod(uk, we)
          }
          mse_r[r] <- sum(we^2) / (as.numeric(rc$df_r - kk) * ncol(we))
        }
        cv_mse_vec[kk + 1L] <- mean(mse_r)
      }
      mse_0 <- cv_mse_vec[1L]
      mse_min <- min(cv_mse_vec)
      rel_gain <- if (is.finite(mse_0) && mse_0 > 0) (mse_0 - mse_min) / mse_0 else 0.0
      if (!is.finite(rel_gain) || rel_gain < 0.015) {
        return(empty_res)
      }
      # Kay et al. (2013) GLMdenoise 95%-of-peak variance-reduction elbow rule:
      gains <- mse_0 - cv_mse_vec
      peak_gain <- max(gains)
      k_star <- which(gains >= 0.95 * peak_gain)[1L] - 1L
      if (!is.finite(k_star) || k_star <= 0L) return(empty_res)
    }
  }

  # Extract final observation-space noise PCs per run from full noise pool
  # (index `bold[[r]]` explicitly per run `r` so R `for`-loop scope never aliases across runs!).
  # In automatic mode, filter out:
  #   (a) temporally localized framewise spikes (`top3_share > 0.22`), leaving transient
  #       motion glitches to `rwls` frame downweighting rather than wasting nuisance DF;
  #   (b) temporally unstructured white-noise spatial factors (`max_acf < 0.28` and
  #       `spec_crest < 11.0` when `noise_pool` is not explicitly supplied), leaving
  #       i.i.d. spatial covariance to Stage-3 `poet_covreg` while targeting structured
  #       narrowband physiological harmonics and correlated global vascular waves.
  u_kept_runs <- vector("list", n_runs)
  for (r in seq_len(n_runs)) {
    rc <- run_cache[[r]]
    yr_r <- bold[[r]][rc$rows_r, rc$sub_full, drop = FALSE]
    wb_r <- hrf_apply_whitening(yr_r, rc$rows_r, rc$ar_bg, unit_marginal = FALSE)
    if (ncol(rc$qz_w) > 0L) wb_r <- wb_r - rc$qz_w %*% crossprod(rc$qz_w, wb_r)
    w_full <- prep_bg_block(wb_r)

    yb_o <- if (ncol(rc$qz_o) > 0L) yr_r - rc$qz_o %*% crossprod(rc$qz_o, yr_r) else yr_r
    y_obs_full <- prep_bg_block(yb_o)

    sv_f <- hrf_safe_svd(w_full, nu = 0L, nv = k_star)
    v_sp <- sv_f$v[, seq_len(k_star), drop = FALSE]
    u_obs <- y_obs_full %*% v_sp

    if (is_auto) {
      keep_cols <- logical(ncol(u_obs))
      for (jj in seq_len(ncol(u_obs))) {
        uj <- u_obs[, jj]
        nrm_u <- sqrt(sum(uj^2))
        if (!is.finite(nrm_u) || nrm_u <= 1e-12) next
        uj <- uj / nrm_u
        mu_n <- length(uj)
        top3_share <- sum(utils::head(sort(uj^2, decreasing = TRUE), min(3L, mu_n)))
        if (top3_share > 0.22) next
        if (is.null(pool_info$explicit_pool) || !isTRUE(pool_info$explicit_pool)) {
          max_lg <- min(6L, max(1L, mu_n - 6L))
          max_acf <- max(abs(vapply(seq_len(max_lg), function(lg) {
            sum(uj[seq_len(mu_n - lg)] * uj[(lg + 1L):mu_n])
          }, numeric(1L))))
          nf_half <- as.integer(floor(mu_n / 2L)) + 1L
          p_sp <- if (nf_half >= 4L) Mod(stats::fft(uj))[2:nf_half]^2 else 0.0
          spec_crest <- if (length(p_sp) >= 3L) max(p_sp) / max(stats::median(p_sp), 1e-12) else 0.0
          if (max_acf < 0.28 && spec_crest < 11.0) next
        }
        keep_cols[jj] <- TRUE
      }
      u_obs <- u_obs[, keep_cols, drop = FALSE]
    }
    u_kept_runs[[r]] <- u_obs
  }

  k_final <- min(vapply(u_kept_runs, ncol, integer(1L)))
  if (!is.finite(k_final) || k_final <= 0L) return(empty_res)

  noise_pcs_out <- vector("list", n_runs)
  names(noise_pcs_out) <- names(bold)
  for (r in seq_len(n_runs)) {
    rc <- run_cache[[r]]
    u_use <- u_kept_runs[[r]][, seq_len(k_final), drop = FALSE]
    sv_uo <- hrf_safe_svd(u_use, nu = k_final, nv = 0L)
    pc_mat <- matrix(0.0, nrow = rc$n_obs, ncol = k_final)
    pc_mat[rc$rows_r, ] <- sv_uo$u[, seq_len(k_final), drop = FALSE] * sqrt(rc$m_r)
    colnames(pc_mat) <- paste0("noise_pc", seq_len(k_final))
    noise_pcs_out[[r]] <- pc_mat
  }

  list(
    noise_pcs = noise_pcs_out,
    k_selected = as.integer(k_final),
    pool_voxels = pool_info$pool_voxels,
    cv_mse = cv_mse_vec
  )
}

hrf_build_decon_nuisance_matrix <- function(n_obs, tr, rows, nuisance_r = NULL,
                                            nuisance_mode = "ignore", highpass_sec = 128.0,
                                            drift_active = TRUE) {
  if (identical(nuisance_mode, "ignore")) {
    return(matrix(numeric(), n_obs, 0L))
  }
  nuis_use <- if (is.null(nuisance_r)) matrix(numeric(), n_obs, 0L) else nuisance_r
  if (identical(nuisance_mode, "preproject") && isFALSE(drift_active) && ncol(nuis_use) == 0L) {
    return(matrix(numeric(), n_obs, 0L))
  }
  z_base <- cbind(intercept = 1.0, trend = seq(-1.0, 1.0, length.out = n_obs), nuis_use)
  if (identical(nuisance_mode, "preproject") && !isFALSE(drift_active)) {
    dct_cols <- hrf_build_dct_drift_cols(n_obs, tr, length(rows), highpass_sec = highpass_sec)
    if (ncol(dct_cols) > 0L) z_base <- cbind(z_base, dct_cols)
  }
  z_base
}

hrf_detect_lowfreq_drift <- function(bold, hrf_fit, observation_mask = "fit", highpass_sec = 128.0) {
  n_runs <- length(bold)
  tr <- as.numeric(hrf_fit$settings$tr)
  fb_idx <- if (!is.null(hrf_fit$settings$fallback_kernel)) as.integer(hrf_fit$settings$fallback_kernel) else 1L
  k_fb <- hrf_fit$kernels[, min(max(1L, fb_idx), ncol(hrf_fit$kernels)), drop = FALSE]
  for (r in seq_len(n_runs)) {
    y <- bold[[r]]
    n <- nrow(y)
    mask <- if (identical(observation_mask, "fit")) hrf_fit$specification$valid_rows[[r]] else rep(TRUE, n)
    rows <- which(!is.na(mask) & mask)
    m <- length(rows)
    if (m < 28L) next
    dct_r <- hrf_build_dct_drift_cols(n, tr, m, highpass_sec = highpass_sec)
    z_dr <- cbind(trend = seq(-1.0, 1.0, length.out = n), dct_r)[rows, , drop = FALSE]
    z_dr <- sweep(z_dr, 2L, colMeans(z_dr), "-")
    if (ncol(z_dr) == 0L || m < ncol(z_dr) + 20L) next

    ev_r <- if (!is.null(hrf_fit$specification$events) && length(hrf_fit$specification$events) >= r) {
      hrf_fit$specification$events[[r]]
    } else NULL
    off_r <- if (length(hrf_fit$settings$time_offset) >= r) hrf_fit$settings$time_offset[r] else 0.0
    x0 <- matrix(1.0, nrow = m, ncol = 1L)
    if (!is.null(ev_r) && nrow(ev_r) > 0L) {
      dur_r <- if (is.null(ev_r$duration)) rep(0.0, nrow(ev_r)) else as.numeric(ev_r$duration)
      cond_r <- if (is.null(ev_r$condition)) rep(1L, nrow(ev_r)) else as.integer(factor(ev_r$condition))
      n_c <- max(1L, max(cond_r, na.rm = TRUE))
      g_c <- hrf_event_design_cpp(n, tr, off_r, ev_r$onset, dur_r, cond_r, n_c, k_fb)[[1L]][rows, , drop = FALSE]
      x0 <- cbind(x0, g_c)
    }

    y_obs <- y[rows, , drop = FALSE]
    ok_v <- which(colSums(!is.finite(y_obs)) == 0L &
      apply(y_obs, 2L, function(col) stats::sd(col) > .Machine$double.eps))
    if (!length(ok_v)) next
    if (length(ok_v) > 64L) {
      ok_v <- ok_v[unique(round(seq(1, length(ok_v), length.out = 64L)))]
    }
    y_sub <- y_obs[, ok_v, drop = FALSE]
    rho_r <- if (!is.null(hrf_fit$noise) && length(hrf_fit$noise) >= r && !is.null(hrf_fit$noise[[r]])) {
      stats::median(hrf_fit$noise[[r]]$rho[ok_v], na.rm = TRUE)
    } else 0.0
    if (!is.finite(rho_r)) rho_r <- 0.0
    rho_det <- min(max(rho_r, 0.0), 0.25)

    wb <- hrf_apply_whitening(list(wx = x0, wz = z_dr, wy = y_sub), rows = rows, ar_vec = rho_det, unit_marginal = TRUE)
    nx <- sqrt(colSums(wb$wx^2)); nx[nx == 0] <- 1.0
    sv0 <- hrf_safe_svd(sweep(wb$wx, 2L, nx, "/"), nv = 0L)
    rk0 <- if (length(sv0$d) && sv0$d[1L] > 0) sum(sv0$d > 1e-8 * sv0$d[1L]) else 0L
    df_tot <- m - rk0
    if (df_tot < ncol(z_dr) + 12L) next
    wz_res <- if (rk0 > 0L) wb$wz - sv0$u[, seq_len(rk0), drop = FALSE] %*% crossprod(sv0$u[, seq_len(rk0), drop = FALSE], wb$wz) else wb$wz
    wy_res <- if (rk0 > 0L) wb$wy - sv0$u[, seq_len(rk0), drop = FALSE] %*% crossprod(sv0$u[, seq_len(rk0), drop = FALSE], wb$wy) else wb$wy

    nz <- sqrt(colSums(wz_res^2)); keep_z <- nz > 1e-8
    if (!any(keep_z)) next
    sv_dr <- hrf_safe_svd(sweep(wz_res[, keep_z, drop = FALSE], 2L, nz[keep_z], "/"), nv = 0L)
    rk_dr <- if (length(sv_dr$d) && sv_dr$d[1L] > 0) sum(sv_dr$d > 1e-6 * sv_dr$d[1L]) else 0L
    df2 <- df_tot - rk_dr
    if (rk_dr < 1L || df2 < 12L) next
    Q_dr <- sv_dr$u[, seq_len(rk_dr), drop = FALSE]
    ss_dr <- colSums(crossprod(Q_dr, wy_res)^2)
    ss_err <- pmax(colSums(wy_res^2) - ss_dr, 1e-12)
    f_v <- (ss_dr / rk_dr) / (ss_err / df2)
    p_v <- stats::pf(f_v, df1 = rk_dr, df2 = df2, lower.tail = FALSE)
    p_run <- hrf_cauchy_p(p_v)
    if (is.finite(p_run) && p_run < 0.005) {
      return(TRUE)
    }
  }
  FALSE
}

hrf_build_noise_mean_design <- function(n_obs, events, tr, offset, kernel_col, z, mask,
                                        highpass_sec = 128.0) {
  n_mask <- sum(mask)
  rank_fn <- function(a) {
    if (ncol(a) == 0L || n_mask == 0L) return(0L)
    a <- a[mask, , drop = FALSE]
    norm <- sqrt(colSums(a^2)); norm[norm == 0] <- 1
    s <- hrf_safe_svd(sweep(a, 2L, norm, "/"), nu = 0, nv = 0)$d
    if (length(s) && s[1L] > 0) sum(s > 1e-8 * s[1L]) else 0L
  }
  if (is.null(events) || nrow(events) == 0L) {
    return(list(
      design = z,
      df_nominal = n_mask - rank_fn(z),
      trial_design = matrix(0, n_obs, 0L),
      used_trial_design = TRUE
    ))
  }
  durations <- if (is.null(events$duration)) rep(0, nrow(events)) else events$duration
  k_mat <- matrix(as.numeric(kernel_col), ncol = 1L)
  x <- hrf_event_design_cpp(n_obs, tr, offset, events$onset, durations,
    seq_len(nrow(events)), nrow(events), k_mat)[[1L]]
  cond_levels <- unique(events$condition)
  cond_ids <- match(events$condition, cond_levels)
  x_cond <- if (length(cond_levels) > 0L) {
    hrf_event_design_cpp(n_obs, tr, offset, events$onset, durations,
      cond_ids, length(cond_levels), k_mat)[[1L]]
  } else {
    matrix(0, n_obs, 0L)
  }
  derivative <- function(a) {
    if (ncol(a) == 0L || nrow(a) < 2L) return(a)
    rbind(a[2:nrow(a), , drop = FALSE], 0) - rbind(0, a[1:(nrow(a) - 1L), , drop = FALSE])
  }
  dct_cols <- hrf_build_dct_drift_cols(n_obs, tr, n_mask, highpass_sec = highpass_sec)
  z_dct <- if (ncol(dct_cols) > 0L) cbind(z, dct_cols) else z
  aug_trial_cond <- cbind(z_dct, x, derivative(x_cond), derivative(derivative(x_cond)))
  aug_trial <- cbind(z_dct, x)
  aug_cond <- cbind(z_dct, x_cond, derivative(x_cond), derivative(derivative(x_cond)))
  min_df_full <- max(36L, as.integer(0.68 * n_mask))
  used_trial <- FALSE
  mean_design <- if (n_mask - rank_fn(aug_trial_cond) >= min_df_full) {
    used_trial <- TRUE
    aug_trial_cond
  } else if (n_mask - rank_fn(aug_trial) >= min_df_full) {
    used_trial <- TRUE
    aug_trial
  } else if (n_mask - rank_fn(cbind(z, x)) >= 20L) {
    used_trial <- TRUE
    cbind(z, x)
  } else if (n_mask - rank_fn(aug_cond) >= 20L) {
    aug_cond
  } else if (n_mask - rank_fn(cbind(z, x_cond)) >= 12L) {
    cbind(z, x_cond)
  } else {
    z
  }
  list(
    design = mean_design,
    df_nominal = n_mask - rank_fn(mean_design),
    trial_design = x,
    used_trial_design = used_trial
  )
}

# Coherent, spike-robust leverage-standardized GLS/REML residual scale estimator
# (Woolrich et al. 2001 Eq. 4 / Seber 1977, robustified against un-scrubbed motion spikes).
# Under N(0, sigma^2 Sigma_cor), unit-marginal whitened residuals e = (I - Q Q') W_unit y
# satisfy Var(e_t) = sigma^2 (1 - h_t), where h_t = ||Q[t, ]||_2^2 is the GLS leverage.
# Leverage-standardized innovations z_t = e_t / sqrt(1 - h_t) have exact variance sigma^2;
# Winsorizing at c = 3.0 * s_mad with exact Gaussian consistency beta(3.0) = E[min(Z^2, 9)]
# retains >99% REML efficiency on clean Gaussian data while bounding motion spikes.
hrf_estimate_gls_noise_sd <- function(y_rows, x_rows, rows, ar_vec, fallback_d2_sd = NULL,
                                      return_innovations = FALSE, trial_rows = NULL,
                                      return_raw_innovations = FALSE, obs_weights = NULL) {
  y_rows <- as.matrix(y_rows)
  m <- length(rows)
  p_vox <- ncol(y_rows)
  want_list <- isTRUE(return_innovations) || isTRUE(return_raw_innovations)
  pack_out <- function(sd_vec, z_u = NULL, z_raw = NULL, r_obs = NULL) {
    if (!want_list) return(sd_vec)
    out_l <- list(sd = sd_vec)
    if (isTRUE(return_innovations)) out_l$z_unit <- z_u
    if (isTRUE(return_raw_innovations)) {
      out_l$z_raw <- z_raw
      out_l$r_obs <- r_obs
    }
    out_l
  }
  make_d2_unit <- function(mat) {
    z_fb <- matrix(0.0, m, p_vox)
    if (m < 3L || p_vox == 0L) return(z_fb)
    ok_c <- which(colSums(!is.finite(mat)) == 0L)
    if (!length(ok_c)) return(z_fb)
    tri_idx <- which((rows[3:m] - rows[1:(m - 2L)]) == 2L)
    if (length(tri_idx) >= 2L) {
      d2_c <- mat[tri_idx + 2L, ok_c, drop = FALSE] -
        2.0 * mat[tri_idx + 1L, ok_c, drop = FALSE] +
        mat[tri_idx, ok_c, drop = FALSE]
      out_rows <- tri_idx + 2L
    } else {
      d2_c <- mat[3:m, ok_c, drop = FALSE] - 2.0 * mat[2:(m - 1L), ok_c, drop = FALSE] + mat[1:(m - 2L), ok_c, drop = FALSE]
      out_rows <- 3:m
    }
    d2_c <- sweep(d2_c, 2L, colMeans(d2_c), "-")
    nrm_c <- sqrt(colSums(d2_c^2))
    good_c <- is.finite(nrm_c) & nrm_c > .Machine$double.eps
    if (any(good_c)) {
      z_fb[out_rows, ok_c[good_c]] <- sweep(d2_c[, good_c, drop = FALSE], 2L, nrm_c[good_c], "/")
    }
    z_fb
  }
  if (m < 4L || p_vox == 0L) {
    sd_na <- rep(NA_real_, p_vox)
    return(pack_out(sd_na, if (isTRUE(return_innovations)) matrix(0.0, m, p_vox) else NULL))
  }
  if (is.null(fallback_d2_sd)) {
    r1 <- if (length(ar_vec) >= 1L) ar_vec[1L] else 0.0
    p1 <- if (length(ar_vec) >= 2L) ar_vec[2L] else 0.0
    p2 <- if (length(ar_vec) >= 3L) ar_vec[3L] else 0.0
    p3 <- if (length(ar_vec) >= 4L) ar_vec[4L] else 0.0
    p4 <- if (length(ar_vec) >= 5L) ar_vec[5L] else 0.0
    fallback_d2_sd <- estimate_decon_noise_sd(
      y_rows, noise_rho = r1, noise_psi1 = p1, noise_psi2 = p2, noise_psi3 = p3, noise_psi4 = p4,
      rows = rows
    )
  }
  if (is.null(x_rows) || m < 20L) {
    sd_fb <- as.numeric(fallback_d2_sd)
    return(pack_out(sd_fb, if (isTRUE(return_innovations)) make_d2_unit(y_rows) else NULL))
  }

  x_rows <- as.matrix(x_rows)
  ok_cols <- which(colSums(!is.finite(y_rows)) == 0L)
  want_trial_proj <- (isTRUE(return_innovations) || isTRUE(return_raw_innovations)) &&
    !is.null(trial_rows) && ncol(trial_rows) > 0L
  y_ok <- if (length(ok_cols)) y_rows[, ok_cols, drop = FALSE] else NULL
  has_obs_w <- !is.null(obs_weights) && length(obs_weights) == m &&
    any(abs(as.numeric(obs_weights) - 1.0) > 1e-6)
  w_sp <- if (has_obs_w && identical(attr(obs_weights, "weight_space"), "innovation")) {
    "innovation"
  } else {
    "woodbury"
  }
  e_act_mat <- NULL
  if (has_obs_w && w_sp == "woodbury") {
    w_cl <- pmin(pmax(as.numeric(obs_weights), 1e-4), 1.0)
    act_idx <- which(w_cl < 1.0 - 1e-6)
    if (length(act_idx) > 0L) {
      e_act_mat <- matrix(0.0, m, length(act_idx))
      d_act <- sqrt(pmax(1.0 / w_cl[act_idx] - 1.0, 0.0))
      for (jj in seq_along(act_idx)) {
        e_act_mat[act_idx[jj], jj] <- d_act[jj]
      }
    }
  }
  w_bundle <- hrf_apply_whitening(
    list(
      wx = x_rows,
      wy = y_ok,
      wt = if (want_trial_proj) as.matrix(trial_rows) else NULL,
      we = e_act_mat
    ),
    rows = rows, ar_vec = ar_vec, unit_marginal = TRUE
  )
  if (has_obs_w) {
    if (w_sp == "innovation") {
      sw <- sqrt(pmax(as.numeric(obs_weights), 1e-8))
      w_bundle$wx <- sweep(w_bundle$wx, 1L, sw, "*")
      if (!is.null(w_bundle$wy)) w_bundle$wy <- sweep(w_bundle$wy, 1L, sw, "*")
      if (!is.null(w_bundle$wt)) w_bundle$wt <- sweep(w_bundle$wt, 1L, sw, "*")
    } else if (!is.null(w_bundle$we) && ncol(w_bundle$we) > 0L) {
      sv_we <- hrf_safe_svd(w_bundle$we, nv = 0L)
      pos_we <- sv_we$d > 1e-8
      if (any(pos_we)) {
        u_we <- sv_we$u[, pos_we, drop = FALSE]
        fac_we <- 1.0 / sqrt(1.0 + sv_we$d[pos_we]^2) - 1.0
        apply_wood <- function(mat_in) {
          if (is.null(mat_in)) return(NULL)
          mat_in + u_we %*% sweep(crossprod(u_we, mat_in), 1L, fac_we, "*")
        }
        w_bundle$wx <- apply_wood(w_bundle$wx)
        w_bundle$wy <- apply_wood(w_bundle$wy)
        w_bundle$wt <- apply_wood(w_bundle$wt)
      }
    }
  }
  wx <- w_bundle$wx
  nrm <- sqrt(colSums(wx^2)); nrm[nrm == 0] <- 1.0
  sv <- hrf_safe_svd(sweep(wx, 2L, nrm, "/"), nv = if (isTRUE(return_raw_innovations)) ncol(wx) else 0L)
  rk <- if (length(sv$d) && sv$d[1L] > 0) sum(sv$d > 1e-8 * sv$d[1L]) else 0L
  df_res <- m - rk
  if (df_res < 16L) {
    sd_fb <- as.numeric(fallback_d2_sd)
    return(pack_out(sd_fb, if (isTRUE(return_innovations)) make_d2_unit(y_rows) else NULL))
  }

  q <- if (rk > 0L) sv$u[, seq_len(rk), drop = FALSE] else matrix(0.0, m, 0L)
  out_sd <- rep(NA_real_, p_vox)
  z_unit <- if (isTRUE(return_innovations)) matrix(0.0, m, p_vox) else NULL
  if (!length(ok_cols)) {
    return(pack_out(out_sd, z_unit))
  }

  wy <- w_bundle$wy
  qtwy <- if (rk > 0L) crossprod(q, wy) else NULL
  e_mat <- if (rk > 0L) wy - q %*% qtwy else wy
  lev <- if (rk > 0L) rowSums(q^2) else rep(0.0, m)
  keep_t <- (1.0 - lev) > 0.10
  use_lev <- sum(keep_t) >= 16L
  inv_sd_t <- if (use_lev) 1.0 / sqrt(1.0 - lev[keep_t]) else NULL
  c_win <- 3.0
  beta_win <- (2 * stats::pnorm(c_win) - 1) - 2 * c_win * stats::dnorm(c_win) +
    2 * (c_win^2) * (1 - stats::pnorm(c_win))

  # For spatial/cross-parcel innovation covariance (`return_innovations = TRUE`) and
  # rWLS frame variance estimation (`return_raw_innovations = TRUE`), if `x_rows`
  # fell back to a condition-mean design (dense rapid ITI), project out the top passband
  # singular modes of `trial_rows` so trial-by-trial neural fluctuations never enter innovations.
  e_innov <- e_mat
  lev_innov <- lev
  keep_innov <- keep_t
  inv_sd_innov <- inv_sd_t
  res_obs_ok <- NULL
  if (isTRUE(return_raw_innovations)) {
    if (rk > 0L) {
      coef_x <- sweep(sv$v[, seq_len(rk), drop = FALSE], 2L, sv$d[seq_len(rk)], "/") %*% qtwy
      coef_x <- sweep(coef_x, 1L, nrm, "/")
      res_obs_ok <- y_ok - x_rows %*% coef_x
    } else {
      res_obs_ok <- y_ok
    }
  }
  if (want_trial_proj && (df_res > 24L)) {
    wt_raw <- w_bundle$wt
    wt <- if (rk > 0L) wt_raw - q %*% crossprod(q, wt_raw) else wt_raw
    nrm_t <- sqrt(colSums(wt^2))
    pos_t <- nrm_t > 1e-8
    if (any(pos_t)) {
      sv_t <- hrf_safe_svd(sweep(wt[, pos_t, drop = FALSE], 2L, nrm_t[pos_t], "/"),
                           nv = if (isTRUE(return_raw_innovations)) sum(pos_t) else 0L)
      rk_max <- max(0L, df_res - 24L)
      rk_t <- if (length(sv_t$d) && sv_t$d[1L] > 0) min(rk_max, sum(sv_t$d > 0.06 * sv_t$d[1L])) else 0L
      if (rk_t > 0L) {
        qt <- sv_t$u[, seq_len(rk_t), drop = FALSE]
        qte <- crossprod(qt, e_innov)
        e_innov <- e_innov - qt %*% qte
        lev_tot <- lev + rowSums(qt^2)
        lev_innov <- lev_tot
        keep_innov <- (1.0 - lev_tot) > 0.10
        if (sum(keep_innov) >= 16L) {
          inv_sd_innov <- 1.0 / sqrt(1.0 - lev_tot[keep_innov])
        } else {
          keep_innov <- keep_t
        }
        if (isTRUE(return_raw_innovations) && !is.null(res_obs_ok)) {
          # Map trial-mode projection back to observation space so res_obs_ok == W_unit^{-1} e_innov
          t_mat_pos <- as.matrix(trial_rows)[, pos_t, drop = FALSE]
          if (rk > 0L) {
            qtwt <- crossprod(q, wt_raw[, pos_t, drop = FALSE])
            coef_tx <- sweep(sv$v[, seq_len(rk), drop = FALSE], 2L, sv$d[seq_len(rk)], "/") %*% qtwt
            coef_tx <- sweep(coef_tx, 1L, nrm, "/")
            t_perp_obs <- t_mat_pos - x_rows %*% coef_tx
          } else {
            t_perp_obs <- t_mat_pos
          }
          coef_t <- sweep(sv_t$v[, seq_len(rk_t), drop = FALSE], 2L, sv_t$d[seq_len(rk_t)], "/") %*% qte
          coef_t <- sweep(coef_t, 1L, nrm_t[pos_t], "/")
          res_obs_ok <- res_obs_ok - t_perp_obs %*% coef_t
        }
      }
    }
  }
  use_lev_innov <- sum(keep_innov) >= 16L

  for (idx_k in seq_along(ok_cols)) {
    j <- ok_cols[idx_k]
    ej <- e_mat[, idx_k]
    s_gls <- sqrt(sum(ej^2) / df_res)
    if (!is.finite(s_gls) || s_gls <= .Machine$double.eps) {
      out_sd[j] <- fallback_d2_sd[j]
      next
    }
    s_est <- s_gls
    if (use_lev) {
      zj <- ej[keep_t] * inv_sd_t
      s0 <- stats::mad(zj, center = 0)
      if (is.finite(s0) && s0 > .Machine$double.eps) {
        u <- pmin(abs(zj), c_win * s0)
        s_est <- sqrt(mean(u^2) / beta_win)
      }
    }
    if (isTRUE(return_innovations)) {
      ej_in <- e_innov[, idx_k]
      if (use_lev_innov) {
        zi <- ej_in[keep_innov] * inv_sd_innov
        si0 <- stats::mad(zi, center = 0)
        ui <- if (is.finite(si0) && si0 > .Machine$double.eps) {
          sign(zi) * pmin(abs(zi), c_win * si0)
        } else zi
        ui <- ui - mean(ui)
        nrm_u <- sqrt(sum(ui^2))
        if (is.finite(nrm_u) && nrm_u > .Machine$double.eps) {
          z_unit[keep_innov, j] <- ui / nrm_u
        }
      } else {
        ui <- ej_in - mean(ej_in)
        nrm_u <- sqrt(sum(ui^2))
        if (is.finite(nrm_u) && nrm_u > .Machine$double.eps) {
          z_unit[, j] <- ui / nrm_u
        }
      }
    }
    # Guard against severe omitted-event misspecification inflating GLS residuals far above Delta^2 y
    d2_j <- fallback_d2_sd[j]
    if (is.finite(d2_j) && d2_j > .Machine$double.eps && s_est > 1.35 * d2_j) {
      s_est <- 1.25 * d2_j
    }
    out_sd[j] <- if (is.finite(s_est) && s_est > .Machine$double.eps) s_est else d2_j
  }

  z_raw_out <- r_obs_out <- NULL
  if (isTRUE(return_raw_innovations)) {
    z_raw_out <- r_obs_out <- matrix(NA_real_, m, p_vox)
    sd_ok <- pmax(out_sd[ok_cols], .Machine$double.eps)
    inv_lev_all <- 1.0 / sqrt(pmax(1.0 - lev_innov, 0.10))
    z_raw_out[, ok_cols] <- sweep(sweep(e_innov, 1L, inv_lev_all, "*"), 2L, sd_ok, "/")
    r_obs_out[, ok_cols] <- sweep(sweep(res_obs_ok, 1L, inv_lev_all, "*"), 2L, sd_ok, "/")
  }
  pack_out(out_sd, z_unit, z_raw_out, r_obs_out)
}

hrf_estimate_noise <- function(y, events, tr, offset, kernels, fallback, z, mask,
    groups, rho_grid, batch_size, fixed_rho = NULL, selected_kernels = NULL,
    spatial_coords = NULL, spatial_fwhm = NULL,
    spatial_graph = NULL, graph_groups = NULL, graph_cross_weight = 1.0,
    max_noise_strata = 3L, kernel_strata = NULL, return_v_quad = FALSE,
    return_group_cor = FALSE) {
  ids <- colnames(y)
  n_vox <- ncol(y)
  labels <- if (is.null(groups)) rep("run", n_vox) else unname(as.character(groups))
  n_obs <- nrow(y)
  rows <- which(mask)
  n_mask <- length(rows)

  attach_v_quad <- function(out_df, z_mat, df_res_est) {
    if (is.null(z_mat) || nrow(z_mat) < 4L || df_res_est <= 2.0) {
      return(out_df)
    }
    ok_z <- colSums(z_mat^2) > 0.5 & is.finite(out_df$rho) & is.finite(out_df$noise_sd)
    if (isTRUE(return_group_cor) && !is.null(groups)) {
      z_sp_dz2 <- if (nrow(z_mat) >= 16L) {
        dz2 <- z_mat[3:nrow(z_mat), , drop = FALSE] -
          2.0 * z_mat[2:(nrow(z_mat) - 1L), , drop = FALSE] +
          z_mat[1:(nrow(z_mat) - 2L), , drop = FALSE]
        nrm_dz2 <- sqrt(colSums(dz2 * dz2))
        nrm_dz2[nrm_dz2 <= .Machine$double.eps] <- 1.0
        sweep(dz2, 2L, nrm_dz2, "/")
      } else {
        NULL
      }
      u_grp <- unique(labels)
      g_cor <- stats::setNames(rep(NA_real_, length(u_grp)), u_grp)
      for (gg in u_grp) {
        ig <- which(labels == gg & ok_z)
        vg <- length(ig)
        if (vg >= 2L) {
          rs_g <- rowSums(z_mat[, ig, drop = FALSE])
          rc_g <- (sum(rs_g * rs_g) - vg) / (as.numeric(vg) * (vg - 1.0))
          if (!is.null(z_sp_dz2)) {
            rs_d <- rowSums(z_sp_dz2[, ig, drop = FALSE])
            rc_d <- (sum(rs_d * rs_d) - vg) / (as.numeric(vg) * (vg - 1.0))
            if (is.finite(rc_d)) rc_g <- min(rc_g, rc_d)
          }
          if (is.finite(rc_g)) g_cor[gg] <- min(0.95, max(0.0, rc_g))
        }
      }
      attr(out_df, "group_spatial_cor") <- g_cor
    }
    if (!isTRUE(return_v_quad) || is.null(spatial_graph)) {
      return(out_df)
    }
    lsd_v <- log(pmax(out_df$noise_sd, 1e-12))
    sd_l <- if (sum(ok_z) > 2L) max(0.50, stats::sd(lsd_v[ok_z], na.rm = TRUE)) else 0.50
    g_int <- if (!is.null(graph_groups) && length(graph_groups) == n_vox) {
      hrf_encode_group_ids(graph_groups, n_vox)
    } else {
      integer(0L)
    }
    vq <- csr_quad_cor_cpp(
      Z_unit = z_mat,
      d_res = as.numeric(df_res_est),
      row_ptr = spatial_graph$row_ptr,
      col_idx = spatial_graph$col_idx,
      dist = spatial_graph$dist,
      group_ids = g_int,
      rho_vec = as.numeric(out_df$rho),
      lsd_vec = as.numeric(lsd_v),
      usable_col = ok_z,
      sigma_sp = as.numeric(spatial_graph$sigma),
      cross_weight = as.numeric(graph_cross_weight),
      rho_scale = 0.25,
      lsd_scale = as.numeric(sd_l)
    )
    attr(out_df, "v_quad") <- as.numeric(vq)
    out_df
  }

  need_z_innov <- (isTRUE(return_v_quad) && !is.null(spatial_graph) && n_mask >= 16L) ||
    (isTRUE(return_group_cor) && !is.null(groups) && n_mask >= 12L)

  if (!is.null(fixed_rho)) {
    r_fix <- unname(fixed_rho)
    nsd_fix <- rep(NA_real_, n_vox)
    z_fix <- if (need_z_innov) matrix(0.0, n_mask, n_vox) else NULL
    df_fix <- max(1.0, n_mask - ncol(z))
    if (n_mask >= 4L) {
      if (is.null(selected_kernels)) {
        md_info <- hrf_build_noise_mean_design(n_obs, events, tr, offset, kernels[, fallback], z, mask)
        x_r <- if (md_info$df_nominal >= 16L) md_info$design[rows, , drop = FALSE] else NULL
        df_fix <- if (md_info$df_nominal >= 16L) md_info$df_nominal else max(1.0, n_mask - ncol(z))
        gls_f <- hrf_estimate_gls_noise_sd(y[rows, , drop = FALSE], x_r, rows, r_fix,
          return_innovations = !is.null(z_fix))
        nsd_fix <- if (is.list(gls_f)) gls_f$sd else gls_f
        if (is.list(gls_f) && !is.null(gls_f$z_unit)) z_fix <- gls_f$z_unit
      } else {
        k_sigs <- if (!is.null(kernel_strata) && length(kernel_strata) == n_vox) {
          paste0("p", kernel_strata)
        } else {
          hrf_fast_col_keys(selected_kernels, digits = 4L)
        }
        df_vec <- numeric(0L)
        for (sg in unique(k_sigs)) {
          idx_s <- which(k_sigs == sg)
          md_info <- hrf_build_noise_mean_design(n_obs, events, tr, offset, selected_kernels[, idx_s[1L]], z, mask)
          x_r <- if (md_info$df_nominal >= 16L) md_info$design[rows, , drop = FALSE] else NULL
          df_vec <- c(df_vec, if (md_info$df_nominal >= 16L) md_info$df_nominal else max(1.0, n_mask - ncol(z)))
          gls_f <- hrf_estimate_gls_noise_sd(y[rows, idx_s, drop = FALSE], x_r, rows, r_fix,
            return_innovations = !is.null(z_fix))
          nsd_fix[idx_s] <- if (is.list(gls_f)) gls_f$sd else gls_f
          if (is.list(gls_f) && !is.null(gls_f$z_unit)) z_fix[, idx_s] <- gls_f$z_unit
        }
        if (length(df_vec)) df_fix <- mean(df_vec)
      }
    }
    out_fix <- data.frame(voxel = ids, group = labels,
      rho = rep(r_fix, n_vox), psi1 = 0, psi2 = 0, psi3 = 0, psi4 = 0,
      noise_sd = nsd_fix, method = "supplied", n_voxels = NA_integer_)
    return(attach_v_quad(out_fix, z_fix, df_fix))
  }
  if (n_mask < 12L) return(data.frame(voxel = ids, group = labels,
    rho = 0, psi1 = 0, psi2 = 0, psi3 = 0, psi4 = 0, noise_sd = NA_real_,
    method = "insufficient_noise_df", n_voxels = 0L))

  # Partition voxels by frozen HRF kernel stratum (either fallback in Pass 1 or
  # cross-run selected_kernels in Pass 2) so non-canonical HRFs do not leak task
  # signal into noise REML profiles or GLS residual variance.
  if (is.null(selected_kernels)) {
    k_sigs <- rep("fb", n_vox)
    k_list <- list(fb = kernels[, fallback])
  } else {
    k_sigs <- if (!is.null(kernel_strata) && length(kernel_strata) == n_vox) {
      paste0("p", kernel_strata)
    } else {
      hrf_fast_col_keys(selected_kernels, digits = 4L)
    }
    u_sigs <- unique(k_sigs)
    k_list <- lapply(u_sigs, function(sg) selected_kernels[, which(k_sigs == sg)[1L]])
    names(k_list) <- u_sigs
  }
  md_list <- lapply(k_list, function(k_col) {
    hrf_build_noise_mean_design(n_obs, events, tr, offset, k_col, z, mask)
  })
  min_df_any <- min(vapply(md_list, `[[`, integer(1), "df_nominal"))
  if (min_df_any < 10L) return(data.frame(voxel = ids, group = labels,
    rho = 0, psi1 = 0, psi2 = 0, psi3 = 0, psi4 = 0, noise_sd = NA_real_,
    method = "insufficient_noise_df", n_voxels = 0L))

  des_list <- lapply(md_list, `[[`, "design")
  profiles <- matrix(NA_real_, length(rho_grid), n_vox)
  for (sg in names(des_list)) {
    idx_s <- which(k_sigs == sg)
    profiles[, idx_s] <- hrf_noise_profile_cpp(
      y[, idx_s, drop = FALSE], des_list[[sg]], rows, rho_grid, as.integer(batch_size)
    )
  }

  out <- data.frame(voxel = ids, group = labels, rho = 0, psi1 = 0, psi2 = 0,
    psi3 = 0, psi4 = 0, noise_sd = NA_real_,
    method = "insufficient_noise_voxels", n_voxels = 0L)

  q_use <- if (n_mask >= 96L && tr < 0.95) 4L else 2L
  lag_idx <- lapply(seq_len(q_use), function(L) {
    if (n_mask > L) which((rows[(1L + L):n_mask] - rows[1L:(n_mask - L)]) == L) else integer(0)
  })

  interp_rho <- function(prof_vec) {
    best <- which.min(prof_vec)
    rho_hat <- rho_grid[best]
    if (best > 1L && best < length(rho_grid) && length(rho_grid) >= 5L) {
      x1 <- rho_grid[best - 1L]; x2 <- rho_grid[best]; x3 <- rho_grid[best + 1L]
      y1 <- prof_vec[best - 1L]; y2 <- prof_vec[best]; y3 <- prof_vec[best + 1L]
      d1 <- x2 - x1; d2 <- x3 - x2
      if (d1 > 0 && d2 > 0) {
        g1 <- (y2 - y1) / d1; g2 <- (y3 - y2) / d2
        if (is.finite(g1) && is.finite(g2) && (g2 - g1) > 1e-12) {
          rho_interp <- x2 - (d1 * g2 + d2 * g1) / (2 * (g2 - g1))
          rho_hat <- round(min(max(rho_interp, x1), x3), 4L)
        }
      }
    }
    list(rho = rho_hat, best = best)
  }

  estimate_psi_weighted <- function(vox_idx, vox_weights, rho_val) {
    psi_zero <- rep(0.0, 4L)
    if (n_mask < 48L || length(lag_idx[[1L]]) < 24L || length(lag_idx[[2L]]) < 24L) {
      return(psi_zero)
    }
    q_eff <- if (q_use == 4L && length(lag_idx[[3L]]) >= 24L && length(lag_idx[[4L]]) >= 24L) 4L else 2L
    r_lags <- numeric(q_eff)
    w_tot <- 0.0
    df_acc <- 0.0
    u_chunks <- list()
    w_chunks <- list()
    vox_sigs <- k_sigs[vox_idx]
    u_sigs_psi <- unique(vox_sigs)
    in_psi_list <- vector("list", 2L * length(u_sigs_psi))
    for (si in seq_along(u_sigs_psi)) {
      sg <- u_sigs_psi[si]
      sub_v <- vox_idx[vox_sigs == sg]
      in_psi_list[[2L * si - 1L]] <- des_list[[sg]][rows, , drop = FALSE]
      in_psi_list[[2L * si]] <- y[rows, sub_v, drop = FALSE]
    }
    w_psi_list <- hrf_apply_whitening(in_psi_list, rows, rho_val, unit_marginal = FALSE)
    for (si in seq_along(u_sigs_psi)) {
      sg <- u_sigs_psi[si]
      sub_pos <- which(vox_sigs == sg)
      sub_w <- vox_weights[sub_pos]
      wz <- w_psi_list[[2L * si - 1L]]
      norm_z <- sqrt(colSums(wz^2)); norm_z[norm_z == 0] <- 1
      sv <- hrf_safe_svd(sweep(wz, 2L, norm_z, "/"), nv = 0)
      rk <- if (length(sv$d) && sv$d[1L] > 0) sum(sv$d > 1e-8 * sv$d[1L]) else 0L
      df_res <- n_mask - rk
      if (df_res < 20L) next
      q_des <- if (rk > 0L) sv$u[, seq_len(rk), drop = FALSE] else matrix(0, n_mask, 0L)
      wy <- w_psi_list[[2L * si]]
      res1 <- if (rk > 0L) wy - q_des %*% crossprod(q_des, wy) else wy
      v_norms <- sqrt(colSums(res1^2))
      ok_v <- is.finite(v_norms) & v_norms > 1e-10
      if (!any(ok_v)) next
      res1_ok <- res1[, ok_v, drop = FALSE]
      res1_n <- sweep(res1_ok, 2L, v_norms[ok_v], "/")
      sw_ok <- sub_w[ok_v]
      sw_sum <- sum(sw_ok)
      if (sw_sum <= 0) next
      m_diag <- if (rk > 0L) pmax(0.0, 1.0 - rowSums(q_des^2)) else rep(1.0, n_mask)
      for (L in seq_len(q_eff)) {
        li <- lag_idx[[L]]
        # Normalize each lag-L product by the local energy of the participating gap-free
        # index pairs (li, li + L) so scattered motion/censor gaps never attenuate emp_L.
        d_exp_L <- sqrt(max(sum(m_diag[li + L]) * sum(m_diag[li]), 1e-12))
        bL <- if (rk > 0L) (-sum(q_des[li + L, , drop = FALSE] * q_des[li, , drop = FALSE])) / d_exp_L else 0.0
        r_hi <- res1_ok[li + L, , drop = FALSE]
        r_lo <- res1_ok[li, , drop = FALSE]
        den_L <- sqrt(pmax(colSums(r_hi^2) * colSums(r_lo^2), 1e-24))
        emp_L <- colSums(r_hi * r_lo) / den_L
        r_lags[L] <- r_lags[L] + sum(sw_ok * (emp_L - bL))
      }
      u_chunks[[length(u_chunks) + 1L]] <- res1_n
      w_chunks[[length(w_chunks) + 1L]] <- sw_ok
      w_tot <- w_tot + sw_sum
      df_acc <- df_acc + sw_sum * df_res
    }
    if (w_tot <= 0) return(psi_zero)
    r_lags <- pmin(pmax(r_lags / w_tot, -0.80), 0.80)
    df_eff <- max(20.0, df_acc / w_tot)
    pair_frac_12 <- min(length(lag_idx[[1L]]), length(lag_idx[[2L]])) / max(1.0, n_mask - 2.0)
    df_eff_12 <- max(16.0, df_eff * min(1.0, pair_frac_12))

    # Exact Isserlis fourth-moment variance factor sum_{i,j} w_i w_j (u_i' u_j)^2
    # accounting for both parcel size and cross-voxel spatial noise correlation.
    u_all <- if (length(u_chunks) == 1L) u_chunks[[1L]] else do.call(cbind, u_chunks)
    w_norm_u <- unlist(w_chunks, use.names = FALSE) / w_tot
    u_w <- sweep(u_all, 2L, sqrt(pmax(w_norm_u, 0)), "*")
    sp_var_factor <- if (ncol(u_w) <= nrow(u_w)) {
      g_mat <- crossprod(u_w)
      sum(g_mat * g_mat)
    } else {
      m_mat <- tcrossprod(u_w)
      sum(m_mat * m_mat)
    }
    sigma_k2_12 <- max(1e-6, sp_var_factor / df_eff_12)
    tau12 <- max(0.060, sqrt(5.9915 * sigma_k2_12))

    if (q_eff == 2L) {
      r1 <- r_lags[1L]; r2 <- r_lags[2L]
      k1 <- r1
      k2 <- min(max((r2 - r1 * r1) / max(1.0 - r1 * r1, 0.10), -0.80), 0.80)
      k_norm2 <- k1 * k1 + k2 * k2
      if (!is.finite(k_norm2) || k_norm2 <= tau12 * tau12) return(psi_zero)
      shrink_eb <- max(0.0, 1.0 - (tau12 * tau12) / k_norm2)
      k1_s <- min(max(shrink_eb * k1, -0.70), 0.70)
      k2_s <- min(max(shrink_eb * k2, -0.70), 0.70)
      psi2 <- round(k2_s, 4L)
      psi1 <- round(k1_s * (1.0 - psi2), 4L)
      return(c(psi1, psi2, 0.0, 0.0))
    }

    # Step-up Levinson-Durbin up to q_eff = 4 with blockwise Empirical-Bayes PACF shrinkage
    pacf_k <- numeric(q_eff)
    phi_ld <- numeric(0)
    pred_v <- 1.0
    for (m_ord in seq_len(q_eff)) {
      num <- if (m_ord == 1L) r_lags[1L] else r_lags[m_ord] - sum(phi_ld * rev(r_lags[1L:(m_ord - 1L)]))
      km <- min(max(num / max(pred_v, 0.12), -0.75), 0.75)
      pacf_k[m_ord] <- km
      phi_ld <- if (m_ord == 1L) km else c(phi_ld - km * rev(phi_ld), km)
      pred_v <- pred_v * (1.0 - km * km)
    }
    pair_frac_34 <- min(length(lag_idx[[3L]]), length(lag_idx[[4L]])) / max(1.0, n_mask - 4.0)
    df_eff_34 <- max(16.0, df_eff * min(1.0, pair_frac_34))
    sigma_k2_34 <- max(1e-6, sp_var_factor / df_eff_34)
    tau34 <- max(0.090, sqrt(9.2103 * sigma_k2_34))
    n12 <- sum(pacf_k[1:2]^2)
    s12 <- if (is.finite(n12) && n12 > tau12^2) max(0.0, 1.0 - (tau12^2) / n12) else 0.0
    pacf_k[1:2] <- pmin(pmax(s12 * pacf_k[1:2], -0.70), 0.70)
    n34 <- sum(pacf_k[3:4]^2)
    s34 <- if (is.finite(n34) && n34 > tau34^2) max(0.0, 1.0 - (tau34^2) / n34) else 0.0
    pacf_k[3:4] <- pmin(pmax(s34 * pacf_k[3:4], -0.55), 0.55)
    if (all(abs(pacf_k) <= 1e-10)) return(psi_zero)
    if (all(abs(pacf_k[3:4]) <= 1e-10)) {
      psi2 <- round(pacf_k[2L], 4L)
      psi1 <- round(pacf_k[1L] * (1.0 - psi2), 4L)
      return(c(psi1, psi2, 0.0, 0.0))
    }
    psi_cur <- numeric(0)
    for (m_ord in seq_len(4L)) {
      km <- pacf_k[m_ord]
      psi_cur <- if (m_ord == 1L) km else c(psi_cur - km * rev(psi_cur), km)
    }
    round(psi_cur, 4L)
  }

  # Damped 2nd-pass alternating REML: pre-whiten demeaned y and design by AR(q)(psi_0),
  # re-profile rho, and refresh psi so short-lag physiological/MA structure does not pull rho.
  refine_ar_fit <- function(vox_idx, vox_weights) {
    w_norm <- vox_weights / sum(vox_weights)
    prof_0 <- as.numeric(profiles[, vox_idx, drop = FALSE] %*% w_norm)
    ip0 <- interp_rho(prof_0)
    rho_hat <- ip0$rho
    best_idx <- ip0$best
    psi_hat <- estimate_psi_weighted(vox_idx, w_norm, rho_hat)
    if (n_mask >= 64L && any(abs(psi_hat) > 1e-4)) {
      prof_2 <- matrix(NA_real_, length(rho_grid), length(vox_idx))
      vox_sigs <- k_sigs[vox_idx]
      u_sigs_ref <- unique(vox_sigs)
      in_ref_list <- vector("list", 2L * length(u_sigs_ref))
      for (si in seq_along(u_sigs_ref)) {
        sg <- u_sigs_ref[si]
        sub_v <- vox_idx[vox_sigs == sg]
        y_rows_sub <- y[rows, sub_v, drop = FALSE]
        in_ref_list[[2L * si - 1L]] <- sweep(y_rows_sub, 2L, colMeans(y_rows_sub), "-")
        in_ref_list[[2L * si]] <- des_list[[sg]][rows, , drop = FALSE]
      }
      w_ref_list <- hrf_apply_whitening(in_ref_list, rows, c(0.0, psi_hat), unit_marginal = TRUE)
      for (si in seq_along(u_sigs_ref)) {
        sg <- u_sigs_ref[si]
        sub_pos <- which(vox_sigs == sg)
        sub_v <- vox_idx[sub_pos]
        y_f <- y[, sub_v, drop = FALSE]
        d_f <- des_list[[sg]]
        y_f[rows, ] <- w_ref_list[[2L * si - 1L]]
        d_f[rows, ] <- w_ref_list[[2L * si]]
        prof_2[, sub_pos] <- hrf_noise_profile_cpp(
          y_f, d_f, rows, rho_grid, as.integer(batch_size)
        )
      }
      ok2 <- colSums(is.finite(prof_2)) == length(rho_grid)
      if (any(ok2)) {
        sw2 <- w_norm[ok2] / sum(w_norm[ok2])
        ip1 <- interp_rho(as.numeric(prof_2[, ok2, drop = FALSE] %*% sw2))
        rho_hat <- round(0.75 * ip1$rho + 0.25 * rho_hat, 4L)
        best_idx <- ip1$best
        psi_hat <- estimate_psi_weighted(vox_idx, w_norm, rho_hat)
      }
    }
    list(rho = rho_hat, psi = psi_hat, boundary = best_idx %in% c(1L, length(rho_grid)))
  }

  if (is.null(spatial_graph) && !is.null(spatial_coords) && n_vox > 2500L) {
    spatial_graph <- build_spatial_graph(
      spatial_coords,
      spatial_fwhm = spatial_fwhm,
      spatial_groups = labels,
      parcel_boundary = "hard",
      metric = "euclidean"
    )
    graph_groups <- labels
    graph_cross_weight <- 0.0
  }
  spatial_w <- if (is.null(spatial_graph) && !is.null(spatial_coords)) {
    hrf_compute_spatial_weights(spatial_coords, labels, spatial_fwhm)
  } else NULL

  # Deterministic 1D quantization into at most K clusters (never touches .Random.seed)
  cluster_1d_det <- function(x, max_k = 3L, include_tails = FALSE) {
    ux <- sort(unique(round(x, 2L)))
    k_cap <- as.integer(max_k)
    if (isTRUE(include_tails) && length(x) > 1L) {
      rng_x <- max(x) - min(x)
      k_cap <- min(k_cap, max(2L, as.integer(ceiling(rng_x / 0.05)) + 1L))
    }
    k <- min(k_cap, length(ux))
    if (k <= 1L) return(rep(1L, length(x)))
    probs <- if (isTRUE(include_tails)) {
      seq(0.02, 0.98, length.out = k)
    } else {
      seq(0.15, 0.85, length.out = k)
    }
    ctrs <- unique(as.numeric(stats::quantile(x, probs = probs, names = FALSE)))
    if (isTRUE(include_tails) && length(ctrs) < k) {
      ctrs <- unique(round(seq(min(x), max(x), length.out = k), 4L))
    }
    if (length(ctrs) <= 1L) return(rep(1L, length(x)))
    for (iter in 1:12) {
      d_mat <- abs(outer(x, ctrs, "-"))
      cl <- max.col(-d_mat, ties.method = "first")
      new_ctrs <- vapply(seq_along(ctrs), function(ci) {
        m_i <- which(cl == ci)
        if (length(m_i)) mean(x[m_i]) else ctrs[ci]
      }, numeric(1))
      if (max(abs(new_ctrs - ctrs)) < 1e-6) break
      ctrs <- new_ctrs
    }
    as.integer(factor(cl))
  }

  compute_gapfree_lsd <- function(vox_cols) {
    tri_idx <- if (n_mask >= 3L) which((rows[3:n_mask] - rows[1:(n_mask - 2L)]) == 2L) else integer(0)
    if (length(tri_idx) >= 8L) {
      d2_u <- y[rows[tri_idx + 2L], vox_cols, drop = FALSE] -
        2 * y[rows[tri_idx + 1L], vox_cols, drop = FALSE] +
        y[rows[tri_idx], vox_cols, drop = FALSE]
    } else {
      pair_idx <- if (n_mask >= 2L) which((rows[2:n_mask] - rows[1:(n_mask - 1L)]) == 1L) else integer(0)
      if (length(pair_idx) >= 6L) {
        d2_u <- y[rows[pair_idx + 1L], vox_cols, drop = FALSE] -
          y[rows[pair_idx], vox_cols, drop = FALSE]
      } else {
        d2_u <- y[rows[3:n_mask], vox_cols, drop = FALSE] -
          2 * y[rows[2:(n_mask - 1L)], vox_cols, drop = FALSE] +
          y[rows[1:(n_mask - 2L)], vox_cols, drop = FALSE]
      }
    }
    log(pmax(apply(d2_u, 2L, stats::mad, center = 0), 1e-12))
  }

  z_innov_all <- if (need_z_innov) {
    matrix(0.0, n_mask, n_vox)
  } else NULL
  df_res_accum <- numeric(0L)

  g_int_graph <- if (!is.null(spatial_graph) && !is.null(graph_groups) && length(graph_groups) == n_vox) {
    hrf_encode_group_ids(graph_groups, n_vox)
  } else {
    integer(0L)
  }

  for (g in unique(labels)) {
    index <- which(labels == g)
    usable <- index[colSums(is.finite(profiles[, index, drop = FALSE])) == length(rho_grid)]
    if (!length(usable)) next

    if (!is.null(spatial_graph) && length(usable) >= 6L && n_mask >= 36L) {
      ind_rhos_u <- vapply(usable, function(v) interp_rho(profiles[, v])$rho, numeric(1))
      df_est <- max(20.0, mean(vapply(md_list, `[[`, integer(1), "df_nominal")))
      r_med <- stats::median(ind_rhos_u)
      exp_sd <- sqrt(max(0.05, 1.0 - r_med * r_med) / df_est)
      raw_lsd_u <- compute_gapfree_lsd(usable)
      lsd_scale <- max(stats::sd(raw_lsd_u), 0.25)
      rho_scale <- max(0.12, 0.50 * stats::sd(ind_rhos_u))

      usable_flag <- rep(FALSE, n_vox)
      usable_flag[usable] <- TRUE
      ind_rhos_full <- rep(0.0, n_vox)
      ind_rhos_full[usable] <- ind_rhos_u
      raw_lsd_full <- rep(0.0, n_vox)
      raw_lsd_full[usable] <- raw_lsd_u

      # If labels partition discrete groups, keep smoothing strictly inside `index` unless single "run" group
      grp_vec_use <- if (length(unique(labels)) > 1L) {
        hrf_encode_group_ids(labels, n_vox)
      } else {
        g_int_graph
      }
      cw_use <- if (length(unique(labels)) > 1L) 0.0 else as.numeric(graph_cross_weight)

      sm_res <- csr_bilateral_noise_smooth_cpp(
        profiles = profiles,
        ind_rhos = ind_rhos_full,
        raw_lsd = raw_lsd_full,
        usable = usable_flag,
        row_ptr = spatial_graph$row_ptr,
        col_idx = spatial_graph$col_idx,
        dist = spatial_graph$dist,
        group_ids = grp_vec_use,
        sigma_sp = as.numeric(spatial_graph$sigma),
        cross_weight = cw_use,
        rho_scale = rho_scale,
        lsd_scale = lsd_scale
      )
      smooth_rhos <- vapply(usable, function(v) interp_rho(sm_res$smooth_profiles[, v])$rho, numeric(1))
      med_eff <- max(1.0, stats::median(sm_res$eff_n[usable], na.rm = TRUE))
      sp_98 <- as.numeric(diff(stats::quantile(smooth_rhos, probs = c(0.02, 0.98), names = FALSE)))
      sp_rng <- max(smooth_rhos) - min(smooth_rhos)
      thr_sp <- max(0.12, 2.6 * exp_sd / sqrt(max(1.0, 0.5 * med_eff)))

      if (is.finite(sp_98) && (sp_98 > thr_sp || sp_rng > max(0.18, 1.5 * thr_sp))) {
        cl <- cluster_1d_det(smooth_rhos, max_k = max_noise_strata, include_tails = TRUE)
        n_cl <- max(cl)
        cl_full <- integer(n_vox)
        cl_full[usable] <- cl
        w_cl_mat <- csr_cluster_member_weights_cpp(
          cluster_id = cl_full,
          n_clusters = n_cl,
          usable = usable_flag,
          row_ptr = spatial_graph$row_ptr,
          col_idx = spatial_graph$col_idx,
          edge_weights = sm_res$edge_weights
        )
        susan_pfx <- if (n_cl > 1L) "susan_" else ""
        for (cid in seq_len(n_cl)) {
          m_pos <- which(cl == cid)
          if (!length(m_pos)) next
          w_vec_c <- w_cl_mat[usable, cid]
          pos_w <- which(w_vec_c >= max(1e-6, 0.01 * max(w_vec_c)))
          fit_c <- if (length(pos_w)) {
            refine_ar_fit(usable[pos_w], w_vec_c[pos_w])
          } else {
            refine_ar_fit(usable[m_pos], rep(1.0, length(m_pos)))
          }
          v_target <- usable[m_pos]
          out$rho[v_target] <- fit_c$rho
          out$psi1[v_target] <- fit_c$psi[1L]
          out$psi2[v_target] <- fit_c$psi[2L]
          out$psi3[v_target] <- fit_c$psi[3L]
          out$psi4[v_target] <- fit_c$psi[4L]
          out$method[v_target] <- paste0(susan_pfx, if (fit_c$boundary) "pooled_reml_grid_boundary" else "pooled_reml")
          out$n_voxels[v_target] <- length(m_pos)
        }
        missing_v <- setdiff(index, usable)
        if (length(missing_v)) {
          maj_cid <- as.integer(names(which.max(table(cl))))
          ref_v <- usable[which(cl == maj_cid)[1L]]
          out$rho[missing_v] <- out$rho[ref_v]
          out$psi1[missing_v] <- out$psi1[ref_v]
          out$psi2[missing_v] <- out$psi2[ref_v]
          out$psi3[missing_v] <- out$psi3[ref_v]
          out$psi4[missing_v] <- out$psi4[ref_v]
          out$method[missing_v] <- out$method[ref_v]
          out$n_voxels[missing_v] <- out$n_voxels[ref_v]
        }
      } else {
        fit_g <- refine_ar_fit(usable, rep(1.0, length(usable)))
        out$rho[index] <- fit_g$rho
        out$psi1[index] <- fit_g$psi[1L]
        out$psi2[index] <- fit_g$psi[2L]
        out$psi3[index] <- fit_g$psi[3L]
        out$psi4[index] <- fit_g$psi[4L]
        out$method[index] <- if (fit_g$boundary) "pooled_reml_grid_boundary" else "pooled_reml"
        out$n_voxels[index] <- length(usable)
      }
    } else {
      # Pool likelihood PROFILES across concordant voxels, not BOLD signals.
      # Woolrich (2001) SUSAN principle: when a parcel mixes distinct noise regimes
      # (e.g., GM core vs CSF/pial-vein or WM partial-volume border voxels) whose
      # 10-90% profile-argmin spread exceeds homogeneous sampling noise, weight each
      # voxel's neighborhood by bilateral similarity + spatial affinity and quantize
      # into <= 3 strata so operator caching remains compact.
      use_susan_g <- FALSE
      if (length(usable) >= 10L && n_mask >= 36L) {
        ind_rhos <- vapply(usable, function(v) interp_rho(profiles[, v])$rho, numeric(1))
        df_est <- max(20.0, mean(vapply(md_list, `[[`, integer(1), "df_nominal")))
        r_med <- stats::median(ind_rhos)
        exp_sd <- sqrt(max(0.05, 1.0 - r_med * r_med) / df_est)
        q90_10 <- as.numeric(diff(stats::quantile(ind_rhos, probs = c(0.10, 0.90), names = FALSE)))
        if (is.finite(q90_10) && q90_10 > max(0.22, 3.8 * exp_sd)) {
          use_susan_g <- TRUE
        }
      }

      if (!use_susan_g) {
        fit_g <- refine_ar_fit(usable, rep(1.0, length(usable)))
        out$rho[index] <- fit_g$rho
        out$psi1[index] <- fit_g$psi[1L]
        out$psi2[index] <- fit_g$psi[2L]
        out$psi3[index] <- fit_g$psi[3L]
        out$psi4[index] <- fit_g$psi[4L]
        out$method[index] <- if (fit_g$boundary) "pooled_reml_grid_boundary" else "pooled_reml"
        out$n_voxels[index] <- length(usable)
      } else {
        raw_lsd <- compute_gapfree_lsd(usable)
        lsd_scale <- max(stats::sd(raw_lsd), 0.25)
        rho_scale <- max(0.12, 0.50 * stats::sd(ind_rhos))
        n_u <- length(usable)
        aff_mat <- matrix(0.0, n_u, n_u)
        smooth_rhos <- numeric(n_u)
        for (ii in seq_len(n_u)) {
          d_rho <- (ind_rhos - ind_rhos[ii]) / rho_scale
          d_lsd <- (raw_lsd - raw_lsd[ii]) / lsd_scale
          w_sim <- exp(-0.5 * (d_rho * d_rho + 0.5 * d_lsd * d_lsd))
          if (!is.null(spatial_w)) {
            w_sim <- w_sim * (0.25 + 0.75 * spatial_w[usable[ii], usable])
          }
          if (sum(w_sim) <= 1e-12) w_sim[ii] <- 1.0
          aff_mat[ii, ] <- w_sim / sum(w_sim)
          smooth_rhos[ii] <- interp_rho(as.numeric(profiles[, usable, drop = FALSE] %*% aff_mat[ii, ]))$rho
        }
        cl <- cluster_1d_det(smooth_rhos, max_k = max_noise_strata)
        susan_pfx <- if (length(unique(cl)) > 1L) "susan_" else ""
        for (cid in unique(cl)) {
          m_pos <- which(cl == cid)
          w_cl <- colMeans(aff_mat[m_pos, , drop = FALSE])
          fit_c <- refine_ar_fit(usable, w_cl)
          v_target <- usable[m_pos]
          out$rho[v_target] <- fit_c$rho
          out$psi1[v_target] <- fit_c$psi[1L]
          out$psi2[v_target] <- fit_c$psi[2L]
          out$psi3[v_target] <- fit_c$psi[3L]
          out$psi4[v_target] <- fit_c$psi[4L]
          out$method[v_target] <- paste0(susan_pfx, if (fit_c$boundary) "pooled_reml_grid_boundary" else "pooled_reml")
          out$n_voxels[v_target] <- length(m_pos)
        }
        missing_v <- setdiff(index, usable)
        if (length(missing_v)) {
          maj_cid <- as.integer(names(which.max(table(cl))))
          ref_v <- usable[which(cl == maj_cid)[1L]]
          out$rho[missing_v] <- out$rho[ref_v]
          out$psi1[missing_v] <- out$psi1[ref_v]
          out$psi2[missing_v] <- out$psi2[ref_v]
          out$psi3[missing_v] <- out$psi3[ref_v]
          out$psi4[missing_v] <- out$psi4[ref_v]
          out$method[missing_v] <- out$method[ref_v]
          out$n_voxels[missing_v] <- out$n_voxels[ref_v]
        }
      }
    }

    # Populate coherent spike-robust GLS residual noise_sd per (HRF, AR) stratum
    g_keys <- paste(k_sigs[index], out$rho[index], out$psi1[index], out$psi2[index],
                    out$psi3[index], out$psi4[index], sep = ":")
    for (gk in unique(g_keys)) {
      sub_v <- index[g_keys == gk]
      v1 <- sub_v[1L]
      ar_v <- c(out$rho[v1], out$psi1[v1], out$psi2[v1], out$psi3[v1], out$psi4[v1])
      df_nom <- md_list[[k_sigs[v1]]]$df_nominal
      df_res_accum <- c(df_res_accum, max(1.0, as.numeric(df_nom)))
      gls_g <- hrf_estimate_gls_noise_sd(
        y[rows, sub_v, drop = FALSE],
        des_list[[k_sigs[v1]]][rows, , drop = FALSE],
        rows,
        ar_v,
        return_innovations = !is.null(z_innov_all)
      )
      out$noise_sd[sub_v] <- if (is.list(gls_g)) gls_g$sd else gls_g
      if (!is.null(z_innov_all) && is.list(gls_g) && !is.null(gls_g$z_unit)) {
        z_innov_all[, sub_v] <- gls_g$z_unit
      }
    }
  }
  df_use <- if (length(df_res_accum)) mean(df_res_accum) else max(1.0, n_mask - ncol(z))
  attach_v_quad(out, z_innov_all, df_use)
}

# Estimate the effective number of independent candidate tests across the HRF
# library. For compact libraries (K <= 9, such as "standard"), exact Bonferroni
# multiplicity (K for task, K - 1 for non-fallback shape departure) is used.
# For extended or dense libraries (K > 9) where candidate HRF shapes are
# strongly positively correlated, standard Bonferroni is overly conservative;
# we compute the Galwey (2009) spectral effective number of tests
# M_eff = (sum_i sqrt(lambda_i))^2 / sum_i lambda_i from the eigenvalues of the
# absolute cosine-similarity matrix |H_norm' H_norm|, floored at 9.
hrf_effective_candidates <- function(kernels) {
  if (is.null(kernels) || !is.matrix(kernels) || ncol(kernels) <= 9L) {
    k <- if (is.matrix(kernels)) ncol(kernels) else 1L
    return(list(task = k, shape = max(1L, k - 1L)))
  }
  k <- ncol(kernels)
  norms <- sqrt(colSums(kernels^2))
  norms[norms <= 0] <- 1
  h_norm <- sweep(kernels, 2L, norms, "/")
  r_mat <- abs(crossprod(h_norm))
  eig <- eigen(r_mat, symmetric = TRUE, only.values = TRUE)$values
  eig <- pmax(eig, 0)
  # Galwey (2009) spectral effective number of independent tests for dense/extended libraries (K > 9).
  m_eff <- if (sum(eig) > 0) (sum(sqrt(eig))^2) / sum(eig) else k
  m_eff <- max(9, min(k, m_eff))
  list(task = m_eff, shape = max(1, m_eff - 1))
}

# Cauchy combination test (Liu & Xie, 2020) for pooling p-values across
# spatially correlated voxels within a parcel. Because the standard Cauchy tail
# is invariant to arbitrary dependency structures, this provides a valid
# parcel-level omnibus test without assuming spatial independence across voxels.
# Conservative upper-tail point masses (e.g., p = 1 from one-sided/clamped tests)
# are bounded at 0.90 so they do not generate -Inf tangents that mask true signals.
hrf_cauchy_p <- function(p, upper_bound = 0.90) {
  p <- p[is.finite(p)]
  if (!length(p)) return(NA_real_)
  p <- pmin(pmax(p, 1e-15), upper_bound)
  t_stat <- mean(tan((0.5 - p) * pi))
  0.5 - atan(t_stat) / pi
}

# Combine per-run calibrated GLS F-test p-values across training runs in a
# leave-target-run-out fold using Fisher's method (-2 * sum log p ~ chi^2_{2R}),
# apply library multiplicity correction (hrf_effective_candidates), and evaluate
# the three-gate calibrated decision rule:
#   1. compatible: no single training run rejects candidate k in favor of a
#      better-scoring alternative at alpha / (R * (K - 1)).
#   2. task_supported: Fisher-combined task F-test p <= alpha AND mean adjusted
#      partial R^2 >= min_r2 AND compatible across runs.
#   3. depart: Fisher-combined pairwise encompassing F-test vs fallback_kernel
#      has p <= alpha AND mean score strictly exceeds fallback_kernel (the
#      fallback kernel itself is exempt from the departure requirement).
hrf_fold_evidence <- function(results, train, counts, scores, fallback, alpha, min_r2,
                              kernels = NULL, subblock_mode = FALSE) {
  k <- nrow(scores); v <- ncol(scores)
  eff <- if (is.null(kernels)) list(task = k, shape = max(1L, k - 1L)) else hrf_effective_candidates(kernels)
  combine <- function(field, multiplicity) {
    logs <- matrix(0, k, v)
    cauchy_sum <- if (isTRUE(subblock_mode)) matrix(0, k, v) else NULL
    for (r in train) {
      p <- results[[r]][[field]]
      if (is.null(p)) next
      ok_p <- is.finite(p)
      p[!ok_p] <- 1
      logs <- logs + log(pmax(p, .Machine$double.xmin))
      if (isTRUE(subblock_mode)) {
        p_c <- pmin(pmax(p, 1e-15), 0.90)
        t_c <- tan((0.5 - p_c) * pi)
        t_c[!ok_p] <- 0
        cauchy_sum <- cauchy_sum + t_c
      }
    }
    raw <- matrix(stats::pchisq(as.vector(-2 * logs),
      df = rep(2 * pmax(counts, 1L), each = k), lower.tail = FALSE), k, v)
    if (isTRUE(subblock_mode)) {
      # Guard against residual within-run sub-block dependence (shared AR fit and
      # boundary carryover) using Brown-scaled Fisher + Cauchy tail bound.
      c_scale <- ifelse(counts >= 2L, 1.15, 1.0)
      raw_brown <- matrix(stats::pchisq(as.vector(sweep(-2 * logs, 2L, c_scale, "/")),
        df = rep(2 * pmax(counts, 1L) / c_scale, each = k), lower.tail = FALSE), k, v)
      raw_cauchy <- 0.5 - atan(sweep(cauchy_sum, 2L, pmax(counts, 1L), "/")) / pi
      raw <- pmax(raw, pmin(raw_brown, raw_cauchy))
    }
    raw[, counts == 0] <- NA_real_
    out <- pmin(raw * multiplicity, 1)
    out[, counts == 0] <- NA_real_
    list(adj = out, raw = raw)
  }
  task_res <- combine("task_p", eff$task)
  shape_res <- combine("fallback_p", eff$shape)
  task <- task_res$adj
  shape <- shape_res$adj
  compatible <- matrix(TRUE, k, v)
  # Simultaneous candidate compatibility within each voxel over all its usable
  # training runs. No averaging can hide one run's strong rejection of a shape.
  # Cap competitor multiplicity by 2 * eff$shape so dense/oversampled HRF grids
  # do not dilute incompatibility sensitivity through near-duplicate bases.
  m_incompat <- min(max(1L, k - 1L), max(eff$shape, as.integer(ceiling(2.0 * eff$shape))))
  cutoff <- alpha / (pmax(counts, 1L) * m_incompat)
  for (r in train) {
    p <- results[[r]]$incompatible_p
    if (is.null(p)) next
    p[!is.finite(p)] <- 1
    compatible <- compatible & sweep(p, 2L, cutoff, ">")
  }
  task_supported <- compatible & task <= alpha & scores >= min_r2
  task_supported[is.na(task_supported)] <- FALSE
  depart <- shape <= alpha & sweep(scores, 2L, scores[fallback, ], ">")
  depart[is.na(depart)] <- FALSE
  depart[fallback, ] <- TRUE
  list(task_p = task, shape_p = shape,
    task_p_raw = task_res$raw, shape_p_raw = shape_res$raw,
    eff_task = eff$task, eff_shape = eff$shape,
    compatible = compatible, task_supported = task_supported,
    depart = depart, supported = task_supported & depart)
}

#SNR-adaptive Empirical-Bayes canonical prior scaling:
# Keeps 100% of canonical_strength on null/weak voxels (z_ev <= 2.5) while
# decaying smoothly as O(1 / z_ev^2) when pooled task evidence is decisive.
hrf_adaptive_canonical_strength <- function(canonical_strength, z_ev) {
  if (!is.finite(canonical_strength) || canonical_strength <= 0 ||
      !is.finite(z_ev) || z_ev <= 2.5) {
    return(canonical_strength)
  }
  canonical_strength / (1.0 + 0.35 * ((z_ev - 2.5)^2))
}

# Pool calibrated evidence and compute Empirical-Bayes departure confidence omega
# across spatially correlated voxels within a parcel (or local neighborhood).
hrf_parcel_pooled_evidence <- function(gate, scores, eligible, fallback_kernel, alpha,
    min_r2, min_parcel_voxels, h_norm, kernel_cos, dist_fb, canonical_strength,
    rho_spatial = NULL) {
  scores_el <- scores[, eligible, drop = FALSE]
  K <- nrow(scores_el)
  V <- ncol(scores_el)
  pooled <- if (V >= 1L) rowMeans(scores_el, na.rm = TRUE) else rep(NA_real_, K)
  if (V < min_parcel_voxels || all(!is.finite(pooled)) || max(pooled, na.rm = TRUE) < min_r2) {
    return(list(
      admissible = if (is.null(gate)) rep(TRUE, K) else (seq_len(K) == fallback_kernel),
      comp_ok = is.finite(pooled) & (pooled >= min_r2),
      omega = if (is.null(gate) && any(is.finite(pooled)) && max(pooled, na.rm = TRUE) >= min_r2) 1.0 else 0.0,
      c_eff = canonical_strength,
      z_task_max = 0.0
    ))
  }
  if (is.null(gate)) {
    comp_ok <- is.finite(pooled) & (pooled >= min_r2)
    comp_ok[fallback_kernel] <- TRUE
    return(list(admissible = rep(TRUE, K), comp_ok = comp_ok, omega = 1.0,
      c_eff = canonical_strength, z_task_max = NA_real_))
  }

  # Estimate spatial NOISE correlation among parcel voxels from training-run
  # whitened task-free innovations (`rho_spatial`) when available; never correlate
  # raw `scores_el` across K directly when signal is present (shared HRF curve -> 0.95!).
  rho_bar <- if (!is.null(rho_spatial) && length(rho_spatial) == 1L && is.finite(rho_spatial)) {
    max(0.05, min(0.95, as.numeric(rho_spatial)))
  } else {
    scores_resid <- sweep(scores_el, 1L, pooled, "-")
    cor_mat <- tryCatch(stats::cor(scores_resid, use = "pairwise.complete.obs"),
      error = function(e) matrix(0.25, V, V))
    rb <- if (V > 1L) mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE) + 1.0 / max(1.0, V - 1.0) else 1.0
    if (!is.finite(rb)) rb <- 0.35
    max(0.10, min(0.90, rb))
  }
  denom_z <- sqrt(V + V * (V - 1) * rho_bar)

  eff_t <- if (!is.null(gate$eff_task)) gate$eff_task else K
  eff_s <- if (!is.null(gate$eff_shape)) gate$eff_shape else max(1L, K - 1L)
  p_task_raw <- if (!is.null(gate$task_p_raw)) gate$task_p_raw[, eligible, drop = FALSE] else
    (gate$task_p[, eligible, drop = FALSE] / eff_t)
  p_shape_raw <- if (!is.null(gate$shape_p_raw)) gate$shape_p_raw[, eligible, drop = FALSE] else
    (gate$shape_p[, eligible, drop = FALSE] / eff_s)
  p_t_clamped <- matrix(pmin(pmax(p_task_raw, 1e-15), 0.99), nrow = K, ncol = V)
  p_s_clamped <- matrix(pmin(pmax(p_shape_raw, 1e-15), 0.99), nrow = K, ncol = V)

  # Spatially-correlated Stouffer/Brown z-test combined with bounded Cauchy test
  z_task <- rowSums(stats::qnorm(1 - p_t_clamped), na.rm = TRUE) / denom_z
  z_shape <- rowSums(stats::qnorm(1 - p_s_clamped), na.rm = TRUE) / denom_z
  z_task_max <- if (any(is.finite(z_task))) max(z_task, na.rm = TRUE) else 0.0
  z_sh_nonfb <- if (K > 1L && any(is.finite(z_shape[-fallback_kernel]))) {
    max(0.0, max(z_shape[-fallback_kernel], na.rm = TRUE))
  } else 0.0
  v_breadth <- min(1.0, max(0.25, (as.numeric(V) - 2.0) / 3.0))
  z_deshrink <- v_breadth * min(z_task_max, 1.45 * z_sh_nonfb)
  c_eff <- hrf_adaptive_canonical_strength(canonical_strength, z_deshrink)

  p_task_stouffer <- pmin(1.0, stats::pnorm(z_task, lower.tail = FALSE) * eff_t)
  p_shape_stouffer <- pmin(1.0, stats::pnorm(z_shape, lower.tail = FALSE) * eff_s)
  p_task_cauchy <- pmin(1.0, apply(p_t_clamped, 1L, hrf_cauchy_p) * eff_t)
  p_shape_cauchy <- pmin(1.0, apply(p_s_clamped, 1L, hrf_cauchy_p) * eff_s)

  p_task_pool <- pmin(p_task_stouffer, p_task_cauchy, na.rm = TRUE)
  p_shape_pool <- pmin(p_shape_stouffer, p_shape_cauchy, na.rm = TRUE)
  p_shape_pool[fallback_kernel] <- 0.0

  comp_cnt <- rowSums(gate$compatible[, eligible, drop = FALSE])
  task_cnt <- rowSums(gate$task_supported[, eligible, drop = FALSE])
  beat_fb <- rowSums(sweep(scores_el, 2L, scores_el[fallback_kernel, ], ">="), na.rm = TRUE)

  adm_strict <- rowSums(gate$supported[, eligible, drop = FALSE]) >= min_parcel_voxels
  comp_ok <- (comp_cnt >= ceiling(0.6 * V)) &
    ((task_cnt >= ceiling(0.5 * min_parcel_voxels)) | (!is.na(p_task_pool) & p_task_pool <= alpha)) &
    is.finite(pooled) & (pooled >= min_r2)

  # Parcel-pooled admission: admit when pooled task AND shape tests pass candidate-corrected
  # alpha across concordant parcel voxels, even when single-voxel power is below Bonferroni.
  adm_pooled <- if (K > 1L) {
    (comp_cnt >= max(min_parcel_voxels, ceiling(0.6 * V))) &
      ((task_cnt >= ceiling(0.5 * min_parcel_voxels)) | (!is.na(p_task_pool) & p_task_pool <= alpha)) &
      (beat_fb >= max(min_parcel_voxels, ceiling(0.6 * V))) &
      (pooled > pooled[fallback_kernel]) &
      !is.na(p_task_pool) & (p_task_pool <= alpha) &
      !is.na(p_shape_pool) & (p_shape_pool <= alpha)
  } else {
    rep(FALSE, K)
  }
  adm_pooled[fallback_kernel] <- (comp_cnt[fallback_kernel] >= min_parcel_voxels) &
    !is.na(p_task_pool[fallback_kernel]) & (p_task_pool[fallback_kernel] <= alpha) &
    is.finite(pooled[fallback_kernel]) & (pooled[fallback_kernel] >= min_r2)
  admissible <- adm_strict | adm_pooled

  # Continuous Empirical-Bayes departure weight omega in [0, 1] for soft_average:
  cand_non_fb <- setdiff(which(comp_ok), fallback_kernel)
  if (!length(cand_non_fb) || K <= 1L) {
    omega <- 0.0
  } else {
    obj_non_fb <- pooled[cand_non_fb] - c_eff * dist_fb[cand_non_fb]
    best_c <- cand_non_fb[which.max(obj_non_fb)]
    tau_best <- h_norm[, best_c] - kernel_cos[fallback_kernel, best_c] * h_norm[, fallback_kernel]
    tau_best <- tau_best / max(sqrt(sum(tau_best^2)), 1e-8)
    opp_cands <- vapply(cand_non_fb, function(k) {
      tk <- h_norm[, k] - kernel_cos[fallback_kernel, k] * h_norm[, fallback_kernel]
      tk <- tk / max(sqrt(sum(tk^2)), 1e-8)
      sum(tau_best * tk)
    }, numeric(1))
    spread_all <- max(pooled[comp_ok], na.rm = TRUE) - min(pooled[comp_ok], na.rm = TRUE)
    max_d_fb <- max(dist_fb[comp_ok], 1e-6)
    expected_drop <- 0.35 * spread_all * (dist_fb[best_c] / max_d_fb)
    diffs <- (scores_el[best_c, ] - c_eff * dist_fb[best_c]) -
      (scores_el[fallback_kernel, ] - expected_drop)
    if (any(opp_cands < -0.5)) {
      opp_c <- cand_non_fb[which.min(opp_cands)]
      diffs_asym <- scores_el[best_c, ] - scores_el[opp_c, ]
      if (mean(diffs_asym, na.rm = TRUE) <= c_eff * dist_fb[best_c]) {
        diffs[] <- -1
      }
    }
    if (!any(is.finite(diffs)) || mean(diffs, na.rm = TRUE) <= 0) {
      omega <- 0.0
    } else {
      sd_d <- stats::sd(diffs, na.rm = TRUE)
      eff_v <- V / (1 + (V - 1) * rho_bar * 0.65)
      t_diff <- if (is.finite(sd_d) && sd_d > 1e-8) mean(diffs, na.rm = TRUE) / (sd_d / sqrt(eff_v)) else 0.0
      # Convert finite-parcel Student-t_{V-1} statistic to exact Gaussian Z equivalent
      # so small parcels (e.g. V = 3, df = 2) never suffer heavy-tail t_2 inflation
      z_t_diff <- if (V >= 2L && is.finite(t_diff) && t_diff > 0) {
        stats::qnorm(pmin(1.0 - 1e-12, stats::pt(t_diff, df = max(1.0, as.numeric(V - 1L)))))
      } else {
        0.0
      }
      cone_sh <- cand_non_fb[kernel_cos[best_c, cand_non_fb] >= 0.92 &
                             is.finite(pooled[cand_non_fb]) & (pooled[cand_non_fb] >= pooled[fallback_kernel])]
      if (!length(cone_sh)) cone_sh <- best_c
      z_sh_cone <- max(z_shape[cone_sh], na.rm = TRUE)
      strong_task_concord <- (V >= 5L) && !is.na(p_task_pool[best_c]) && (p_task_pool[best_c] <= alpha) &&
        (max(beat_fb[cone_sh], na.rm = TRUE) >= ceiling(0.50 * V))
      t_mult <- if (strong_task_concord) 0.95 else 0.80
      z_eff <- max(z_sh_cone, t_mult * z_t_diff, na.rm = TRUE)
      z_lo <- if (strong_task_concord && is.finite(z_deshrink) && z_deshrink > 3.0) {
        max(1.15, 1.35 - 0.10 * (z_deshrink - 3.0))
      } else {
        1.35
      }
      z_hi <- if (strong_task_concord && is.finite(z_deshrink) && z_deshrink > 3.0) {
        max(1.85, 2.50 - 0.30 * (z_deshrink - 3.0))
      } else {
        2.50
      }
      omega <- max(0.0, min(1.0, (z_eff - z_lo) / max(z_hi - z_lo, 0.25)))
    }
  }
  comp_ok[fallback_kernel] <- TRUE
  list(admissible = admissible, comp_ok = comp_ok, omega = omega,
    c_eff = c_eff, z_task_max = z_deshrink,
    p_task_pool = p_task_pool, p_shape_pool = p_shape_pool)
}

# Spread-normalized local Boltzmann evidence weights blended with Empirical-Bayes
# fallback prior weight (1 - omega). When calibrated `z_task_max` is strong (> 3.0),
# concentrates the local Boltzmann temperature and cone size so high-SNR composite
# HRFs converge cleanly onto the true sub-grid waveform without boundary pull.
hrf_soft_boltzmann_weights <- function(scores_vec, kernel_cos, dist_fb,
    fallback_kernel, comp_ok, omega, canonical_strength = 0.01,
    top_k = NULL, temp = 0.20, quantize_step = 0.02, z_task_max = NULL) {
  K <- length(scores_vec)
  w_fb <- rep(0.0, K)
  w_fb[fallback_kernel] <- 1.0
  if (K <= 1L || !is.finite(omega) || omega <= 1e-4) return(w_fb)
  comp_ok <- !is.na(comp_ok) & comp_ok & is.finite(scores_vec)
  if (!any(comp_ok)) return(w_fb)
  has_z <- !is.null(z_task_max) && length(z_task_max) == 1L && is.finite(z_task_max)
  c_eff <- if (has_z) hrf_adaptive_canonical_strength(canonical_strength, z_task_max) else canonical_strength
  temp_eff <- if (has_z && z_task_max > 3.5) {
    max(0.14, temp / (1.0 + 0.06 * (z_task_max - 3.5)))
  } else {
    temp
  }

  obj_adm <- scores_vec - c_eff * dist_fb
  obj_adm[!comp_ok] <- -Inf
  if (!any(is.finite(obj_adm))) return(w_fb)
  k_top <- which.max(obj_adm)

  # Retain immediate local sub-grid bracket neighbors (cos >= 0.90, R2 >= 0.90 * best)
  # around the compatible winner k_top so high-SNR sub-grid HRFs interpolate smoothly
  # even when pairwise encompassing tests reject individual grid endpoints alone.
  local_cone <- comp_ok | (is.finite(scores_vec) & (kernel_cos[k_top, ] >= 0.90) &
    (scores_vec >= max(0.01, 0.90 * scores_vec[k_top])))
  if (is.null(top_k)) {
    top_k <- max(4L, min(8L, as.integer(ceiling(0.30 * sum(local_cone)))))
  }

  obj <- scores_vec - c_eff * dist_fb
  obj[!local_cone] <- -Inf
  ord <- order(obj, decreasing = TRUE)
  ord <- ord[is.finite(obj[ord]) & (kernel_cos[k_top, ord] >= 0.85)]
  if (!length(ord)) return(w_fb)
  top_idx <- utils::head(ord, min(as.integer(top_k), length(ord)))

  sc_ref <- scores_vec[is.finite(scores_vec)]
  spread <- max(sc_ref, na.rm = TRUE) - stats::median(sc_ref, na.rm = TRUE)
  spread <- max(spread, 1e-4)
  delta <- (obj[top_idx] - max(obj[top_idx])) / (temp_eff * spread)
  w_sub <- exp(pmax(-20, delta))
  w_data <- rep(0.0, K)
  w_data[top_idx] <- w_sub / sum(w_sub)

  w_fin <- omega * w_data + (1.0 - omega) * w_fb
  if (quantize_step > 0) {
    w_q <- round(w_fin / quantize_step) * quantize_step
    if (sum(w_q) > 0) w_fin <- w_q / sum(w_q)
  }
  w_fin / sum(w_fin)
}

# Resolve optional Step-1 parametric trial modulators (`hrf_modulators`) into a
# run-aligned list of numeric event-by-modulator matrices. Accepts:
#   - character vector of column names in `events` (e.g., `c("rpe", "rt")`)
#   - run-named or positional list of numeric vectors / matrices / data.frames
#   - single numeric vector / matrix / data.frame for single-run inputs
hrf_resolve_modulators_input <- function(hrf_modulators, events, bold) {
  if (is.null(hrf_modulators)) return(NULL)
  n_runs <- length(bold)
  run_names <- names(bold)

  coerce_mod_mat <- function(x, m_ev, r_lbl) {
    if (is.null(x) || m_ev == 0L) return(matrix(0.0, nrow = m_ev, ncol = 0L))
    if (is.data.frame(x)) {
      num_ok <- vapply(x, is.numeric, logical(1L))
      if (!all(num_ok)) {
        stop(sprintf("'hrf_modulators' for run '%s' contains non-numeric columns.", r_lbl), call. = FALSE)
      }
      mat <- as.matrix(x)
    } else if (is.matrix(x)) {
      if (!is.numeric(x)) {
        stop(sprintf("'hrf_modulators' matrix for run '%s' must be numeric.", r_lbl), call. = FALSE)
      }
      mat <- x
    } else if (is.numeric(x) && is.null(dim(x))) {
      mat <- matrix(as.numeric(x), ncol = 1L)
    } else {
      stop(sprintf("Unsupported 'hrf_modulators' entry for run '%s'.", r_lbl), call. = FALSE)
    }
    if (nrow(mat) != m_ev) {
      stop(sprintf(
        "'hrf_modulators' for run '%s' has %d rows, but events has %d rows.",
        r_lbl, nrow(mat), m_ev
      ), call. = FALSE)
    }
    storage.mode(mat) <- "double"
    mat[!is.finite(mat)] <- 0.0
    mat
  }

  out <- stats::setNames(vector("list", n_runs), run_names)
  if (is.character(hrf_modulators)) {
    mod_cols <- unique(hrf_modulators[nzchar(hrf_modulators)])
    if (!length(mod_cols)) return(NULL)
    for (r in seq_len(n_runs)) {
      ev_r <- events[[r]]
      miss_c <- setdiff(mod_cols, names(ev_r))
      if (length(miss_c)) {
        stop(sprintf(
          "Column(s) %s specified in 'hrf_modulators' not found in events for run '%s'.",
          paste(miss_c, collapse = ", "), run_names[r]
        ), call. = FALSE)
      }
      out[[r]] <- coerce_mod_mat(ev_r[, mod_cols, drop = FALSE], nrow(ev_r), run_names[r])
    }
    return(out)
  }

  if (is.list(hrf_modulators) && !is.data.frame(hrf_modulators)) {
    if (!is.null(names(hrf_modulators)) && any(nzchar(names(hrf_modulators)))) {
      if (!all(run_names %in% names(hrf_modulators))) {
        stop("Named 'hrf_modulators' list must contain entries for all runs in 'bold'.", call. = FALSE)
      }
      src <- hrf_modulators[run_names]
    } else if (length(hrf_modulators) == n_runs) {
      src <- hrf_modulators
    } else {
      stop("'hrf_modulators' list length must match length(bold).", call. = FALSE)
    }
    for (r in seq_len(n_runs)) {
      out[[r]] <- coerce_mod_mat(src[[r]], nrow(events[[r]]), run_names[r])
    }
    return(out)
  }

  if (n_runs == 1L && (is.numeric(hrf_modulators) || is.data.frame(hrf_modulators) || is.matrix(hrf_modulators))) {
    out[[1L]] <- coerce_mod_mat(hrf_modulators, nrow(events[[1L]]), run_names[1L])
    return(out)
  }

  stop(
    "'hrf_modulators' must be NULL, event column name(s), or a per-run list of numeric vectors/matrices/data.frames.",
    call. = FALSE
  )
}

# Build a K-list of candidate HRF design matrices combining unpenalized condition
# intercepts with within-condition centered parametric modulator columns so zero-mean
# modulators (e.g., signed RPE with beta_0 = 0) inform HRF shape selection cleanly.
hrf_build_modulated_design <- function(n_obs, ev, mod_mat, tr, offset, kernels, mask) {
  m_ev <- nrow(ev)
  labels <- as.character(ev$condition)
  durations <- if (is.null(ev$duration)) rep(0.0, m_ev) else as.numeric(ev$duration)
  u_cond <- unique(labels)
  n_cond <- length(u_cond)
  cond_cols <- match(labels, u_cond)

  if (is.null(mod_mat) || !is.matrix(mod_mat) || ncol(mod_mat) == 0L || m_ev < 2L) {
    return(hrf_event_design_cpp(n_obs, tr, offset, ev$onset, durations, cond_cols, n_cond, kernels))
  }

  # Event-to-regressor weight matrix W_ev: [B_cond, M_mod_centered]
  w_cols <- list()
  for (ic in seq_len(n_cond)) {
    v_c <- as.numeric(cond_cols == ic)
    w_cols[[length(w_cols) + 1L]] <- v_c
  }
  n_mask <- sum(mask)
  per_cond_ok <- (n_cond * (1L + ncol(mod_mat))) <= max(4L, as.integer(floor(0.25 * n_mask)))
  for (jm in seq_len(ncol(mod_mat))) {
    mj <- as.numeric(mod_mat[, jm])
    mj[!is.finite(mj)] <- 0.0
    if (per_cond_ok) {
      for (ic in seq_len(n_cond)) {
        idx_c <- which(cond_cols == ic)
        if (length(idx_c) >= 3L) {
          mc <- mj[idx_c] - mean(mj[idx_c])
          sd_c <- stats::sd(mc)
          if (is.finite(sd_c) && sd_c > 1e-8) {
            col_m <- numeric(m_ev)
            col_m[idx_c] <- mc / sd_c
            w_cols[[length(w_cols) + 1L]] <- col_m
          }
        }
      }
    } else {
      mc <- mj - ave(mj, cond_cols, FUN = mean)
      sd_c <- stats::sd(mc)
      if (is.finite(sd_c) && sd_c > 1e-8) {
        w_cols[[length(w_cols) + 1L]] <- mc / sd_c
      }
    }
  }
  w_ev <- do.call(cbind, w_cols)
  if (ncol(w_ev) <= n_cond) {
    return(hrf_event_design_cpp(n_obs, tr, offset, ev$onset, durations, cond_cols, n_cond, kernels))
  }

  des_trial <- hrf_event_design_cpp(n_obs, tr, offset, ev$onset, durations, seq_len(m_ev), m_ev, kernels)
  # Ensure columns visible on `mask` across ALL K candidates are retained so candidate ranks stay comparable
  col_keep <- vapply(seq_len(ncol(w_ev)), function(jc) {
    wj <- w_ev[, jc, drop = FALSE]
    nrms <- vapply(des_trial, function(xt) sqrt(sum((xt[mask, , drop = FALSE] %*% wj)^2)), numeric(1L))
    all(is.finite(nrms)) && min(nrms) > 1e-5
  }, logical(1L))
  if (!any(col_keep)) {
    return(hrf_event_design_cpp(n_obs, tr, offset, ev$onset, durations, cond_cols, n_cond, kernels))
  }
  w_keep <- w_ev[, col_keep, drop = FALSE]
  lapply(des_trial, function(xt) xt %*% w_keep)
}

# SVD Fractional-Ridge Trial HRF Scoring (Rokem & Kay 2020 style matched-DF shrinkage):
#   - Retains the condition-mean subspace (and any supplied `mod_mat`) unpenalized (`w = 1`),
#     while decomposing orthogonal within-condition trial-to-trial deviations via SVD and
#     solving by bisection for candidate-specific ridge parameter alpha_k so EVERY candidate
#     HRF k in 1..K has the EXACT same effective degrees of freedom (`nu_1 = r_c + df_dev`).
#   - Prevents rapid-ITI collinearity blowup and derivative-spanning overfit of unregularized
#     `amplitude_model = "trial"`, while cleanly recovering HRF shapes for unmodeled zero-mean
#     parametric / endogenous trial variations (`beta_0 = 0, Var(beta_trial) > 0`).
#   - Computes independent orthogonal-complement residual variance (`rss_orth`), Satterthwaite
#     `task_p`, and Vuong quadratic-form difference `fallback_p` & `incompatible_p`.
hrf_score_run_fracridge <- function(y_blk, des_cond, des_trial, z_score, rows,
                                    ar_vec = 0.0, fracridge = 0.25,
                                    evidence = TRUE, fallback_index = 1L) {
  K <- length(des_trial)
  n_vox <- ncol(y_blk)
  m <- length(rows)
  na_mat <- matrix(NA_real_, K, n_vox)
  fail_out <- function(status_lbl) {
    out <- list(
      scores = na_mat,
      signal_status = rep(status_lbl, n_vox),
      task_rank = rep(NA_integer_, K),
      nuisance_rank = NA_integer_,
      n_observations = m,
      comparable = FALSE
    )
    if (isTRUE(evidence)) {
      out$task_p <- na_mat
      out$fallback_p <- na_mat
      out$incompatible_p <- na_mat
    }
    out
  }
  if (m < 12L || K < 1L || n_vox < 1L) return(fail_out("insufficient_observations"))

  st_y <- decon_signal_status(y_blk[rows, , drop = FALSE])
  ok_v <- which(st_y == "ok")

  z_obs <- z_score[rows, , drop = FALSE]
  wz <- hrf_apply_whitening(z_obs, rows, ar_vec, unit_marginal = FALSE)
  nrm_z <- sqrt(colSums(wz^2)); nrm_z[nrm_z == 0] <- 1.0
  sv_z <- hrf_safe_svd(sweep(wz, 2L, nrm_z, "/"), nv = 0L)
  rk_z <- if (length(sv_z$d) && sv_z$d[1L] > 0) sum(sv_z$d > 1e-8 * sv_z$d[1L]) else 0L
  qz <- if (rk_z > 0L) sv_z$u[, seq_len(rk_z), drop = FALSE] else matrix(0.0, m, 0L)
  df_nuis <- m - rk_z
  if (df_nuis < 16L || !length(ok_v)) {
    out <- fail_out(if (df_nuis < 16L) "insufficient_observations" else "ok")
    out$signal_status <- st_y
    out$nuisance_rank <- as.integer(rk_z)
    return(out)
  }

  wy <- hrf_apply_whitening(y_blk[rows, ok_v, drop = FALSE], rows, ar_vec, unit_marginal = FALSE)
  if (rk_z > 0L) wy <- wy - qz %*% crossprod(qz, wy)
  tss <- colSums(wy^2)
  nuis_only <- tss <= 1e-12
  if (all(nuis_only)) {
    st_y[ok_v[nuis_only]] <- "nuisance_only"
    out <- fail_out("nuisance_only")
    out$signal_status <- st_y
    out$nuisance_rank <- as.integer(rk_z)
    return(out)
  }

  # Step 1: Orthonormalize unpenalized condition basis Q_c(k) and SVD of orthogonal trial deviations
  qc_list <- ud_list <- sd_list <- s2_norm_list <- vector("list", K)
  rc_vec <- rd_vec <- integer(K)
  for (k in seq_len(K)) {
    xc_k <- hrf_apply_whitening(des_cond[[k]][rows, , drop = FALSE], rows, ar_vec, unit_marginal = FALSE)
    if (rk_z > 0L) xc_k <- xc_k - qz %*% crossprod(qz, xc_k)
    nc_k <- sqrt(colSums(xc_k^2)); pos_c <- nc_k > 1e-8
    if (any(pos_c)) {
      sv_c <- hrf_safe_svd(sweep(xc_k[, pos_c, drop = FALSE], 2L, nc_k[pos_c], "/"), nv = 0L)
      rc_k <- if (length(sv_c$d) && sv_c$d[1L] > 0) sum(sv_c$d > 1e-7 * sv_c$d[1L]) else 0L
      qc_k <- if (rc_k > 0L) sv_c$u[, seq_len(rc_k), drop = FALSE] else matrix(0.0, m, 0L)
    } else {
      rc_k <- 0L
      qc_k <- matrix(0.0, m, 0L)
    }
    rc_vec[k] <- rc_k
    qc_list[[k]] <- qc_k

    xt_k <- hrf_apply_whitening(des_trial[[k]][rows, , drop = FALSE], rows, ar_vec, unit_marginal = FALSE)
    if (rk_z > 0L) xt_k <- xt_k - qz %*% crossprod(qz, xt_k)
    if (rc_k > 0L) xt_k <- xt_k - qc_k %*% crossprod(qc_k, xt_k)
    col_ms_k <- if (ncol(xt_k) > 0L) max(mean(colSums(xt_k^2)), 1e-12) else 1.0
    sv_d <- hrf_safe_svd(xt_k, nv = 0L)
    rd_k <- if (length(sv_d$d) && sv_d$d[1L] > 0) sum(sv_d$d > 1e-5 * sv_d$d[1L]) else 0L
    rd_vec[k] <- rd_k
    ud_list[[k]] <- if (rd_k > 0L) sv_d$u[, seq_len(rd_k), drop = FALSE] else matrix(0.0, m, 0L)
    sd_list[[k]] <- if (rd_k > 0L) sv_d$d[seq_len(rd_k)] else numeric(0L)
    s2_norm_list[[k]] <- if (rd_k > 0L) (sv_d$d[seq_len(rd_k)]^2) / col_ms_k else numeric(0L)
  }

  if (length(unique(rc_vec)) != 1L || all(rc_vec + rd_vec == 0L)) {
    return(fail_out("noncomparable_candidates"))
  }
  rc <- rc_vec[1L]
  # Retain an independent orthogonal residual subspace of >= max(12, 0.42 * (df_nuis - rc)) DF
  max_keep <- max(0L, min(min(rd_vec), as.integer(floor(0.58 * (df_nuis - rc))), df_nuis - rc - 12L))
  if (max_keep == 0L && rc == 0L) {
    return(fail_out("insufficient_observations"))
  }

  frac_clamped <- min(0.85, max(0.05, as.numeric(fracridge)))
  df_dev_target <- if (max_keep > 0L) {
    min(as.numeric(max_keep) * 0.92, max(1.0, frac_clamped * stats::median(rd_vec)))
  } else {
    0.0
  }

  solve_ridge_weights <- function(s_vec, target_df) {
    if (!length(s_vec) || target_df <= 0) return(numeric(0L))
    s2 <- s_vec^2
    if (length(s2) <= target_df + 1e-6) return(rep(1.0, length(s2)))
    lo <- 0.0
    hi <- sum(s2) / max(target_df, 1e-4) * 4.0 + 1e-8
    while (sum(s2 / (s2 + hi)) > target_df) hi <- hi * 4.0
    for (it in seq_len(28L)) {
      mid <- 0.5 * (lo + hi)
      if (sum(s2 / (s2 + mid)) > target_df) lo <- mid else hi <- mid
    }
    s2 / (s2 + 0.5 * (lo + hi))
  }

  # Shared monotone singular-mode weights for overall task-presence Satterthwaite F-test
  sd_ref <- if (max_keep > 0L) {
    s2_mat <- vapply(seq_len(K), function(k) {
      s2k <- sd_list[[k]][seq_len(max_keep)]^2
      s2k / max(sum(s2k), 1e-16)
    }, numeric(max_keep))
    if (is.null(dim(s2_mat))) s2_mat <- matrix(s2_mat, nrow = max_keep)
    sqrt(pmax(0.0, rowMeans(s2_mat)))
  } else {
    numeric(0L)
  }
  w_shared <- solve_ridge_weights(sd_ref, df_dev_target)

  # Empirical-Bayes variance-ratio estimate gamma_hat per voxel for exact Gaussian
  # marginal (REML) fractional-ridge HRF discrimination across fast & slow Toeplitz spectra
  fb_idx <- min(max(1L, as.integer(fallback_index)), K)
  gamma_min <- max(0.12, 0.50 * (frac_clamped / max(0.10, 1.0 - frac_clamped)))
  qc_fb <- qc_list[[fb_idx]]
  ud_fb_keep <- if (max_keep > 0L) ud_list[[fb_idx]][, seq_len(max_keep), drop = FALSE] else matrix(0.0, m, 0L)
  ssr_c_fb <- if (rc > 0L) colSums(crossprod(qc_fb, wy)^2) else rep(0.0, length(ok_v))
  ssr_d_fb <- if (ncol(ud_fb_keep) > 0L) colSums(crossprod(ud_fb_keep, wy)^2) else rep(0.0, length(ok_v))
  sig2_eb <- pmax((tss - ssr_c_fb - ssr_d_fb) / max(1.0, as.numeric(df_nuis - rc - ncol(ud_fb_keep))), 1e-12)
  tr_s2_fb <- if (max_keep > 0L && length(s2_norm_list[[fb_idx]]) >= max_keep) {
    max(sum(s2_norm_list[[fb_idx]][seq_len(max_keep)]), 1e-8)
  } else 1.0
  g_hat_v <- if (max_keep > 0L) {
    pmax(gamma_min, pmin(35.0, pmax(0.0, ssr_d_fb - ncol(ud_fb_keep) * sig2_eb) / pmax(tr_s2_fb * sig2_eb, 1e-12)))
  } else {
    rep(gamma_min, length(ok_v))
  }
  active_g <- g_hat_v[g_hat_v > gamma_min * 1.15]
  gamma_pair <- if (length(active_g)) max(0.40, min(25.0, stats::median(active_g))) else max(0.40, gamma_min * 2.0)

  n_ok <- length(ok_v)
  ll_pair_mat <- scores_ok <- matrix(NA_real_, K, n_ok)
  task_p_ok <- matrix(1.0, K, n_ok)
  nu1_vec <- numeric(K)
  df_rem <- max(1.0, as.numeric(df_nuis - rc))

  for (k in seq_len(K)) {
    qc_k <- qc_list[[k]]
    ud_k_full <- ud_list[[k]]
    s2_k_full <- s2_norm_list[[k]]
    ssr_c_k <- if (rc > 0L) colSums(crossprod(qc_k, wy)^2) else rep(0.0, n_ok)
    rss_c_k <- pmax(1e-12, tss - ssr_c_k)
    proj_d_full <- if (ncol(ud_k_full) > 0L) crossprod(ud_k_full, wy)^2 else matrix(0.0, 0L, n_ok)

    # Exact Gaussian marginal / REML profile score per voxel
    if (length(s2_k_full) > 0L) {
      gs_mat <- outer(s2_k_full, g_hat_v)
      w_eb_mat <- gs_mat / (1.0 + gs_mat)
      qf_eb <- pmax(1e-12, rss_c_k - colSums(w_eb_mat * proj_d_full))
      ld_eb <- colSums(log1p(gs_mat))
    } else {
      qf_eb <- rss_c_k
      ld_eb <- rep(0.0, n_ok)
    }
    ll_gain_k <- -0.5 * (as.numeric(df_nuis) * log(qf_eb / pmax(tss, 1e-12)) + as.numeric(rc) + ld_eb)
    scores_ok[k, ] <- pmax(-0.25, 1.0 - exp(-2.0 * ll_gain_k / as.numeric(df_nuis)))

    # Shared-spectrum Satterthwaite task F-test + fixed-gamma_pair REML log-likelihood for exact pairwise LR tests
    ud_k <- if (max_keep > 0L) ud_k_full[, seq_len(max_keep), drop = FALSE] else matrix(0.0, m, 0L)
    d_k <- c(rep(1.0, rc), w_shared)
    nu1_k <- sum(d_k)
    nu2_k <- sum(d_k^2)
    nu1_vec[k] <- nu1_k

    if (isTRUE(evidence)) {
      proj_d_keep <- if (max_keep > 0L) proj_d_full[seq_len(max_keep), , drop = FALSE] else matrix(0.0, 0L, n_ok)
      ssr_full_k <- ssr_c_k + colSums(proj_d_keep)
      ssr1_k <- ssr_c_k + if (length(w_shared)) as.numeric(crossprod(w_shared, proj_d_keep)) else 0.0
      rss_orth_k <- pmax(0.0, tss - ssr_full_k)
      df2_k <- max(1.0, as.numeric(df_nuis - rc - ncol(ud_k)))
      df1_satt <- max(1.0, (nu1_k^2) / max(nu2_k, 1e-8))
      f_stat <- (ssr1_k / max(nu1_k, 1e-8)) / pmax(rss_orth_k / df2_k, 1e-16)
      task_p_ok[k, ] <- stats::pf(pmax(0.0, f_stat), df1 = df1_satt, df2 = df2_k, lower.tail = FALSE)

      w_p <- if (length(s2_k_full)) (gamma_pair * s2_k_full) / (1.0 + gamma_pair * s2_k_full) else numeric(0L)
      w_c_p <- 25.0 / 26.0
      qf_p <- pmax(1e-12, tss - w_c_p * ssr_c_k - if (length(w_p)) as.numeric(crossprod(w_p, proj_d_full)) else 0.0)
      ld_p <- rc * log(26.0) + sum(log1p(gamma_pair * s2_k_full))
      ll_pair_mat[k, ] <- -0.5 * (as.numeric(df_nuis) * log(qf_p) + ld_p)
    }
  }

  scores_ok[, nuis_only] <- NA_real_
  st_y[ok_v[nuis_only]] <- "nuisance_only"

  scores_all <- matrix(NA_real_, K, n_vox)
  scores_all[, ok_v] <- scores_ok
  task_rk_int <- as.integer(round(nu1_vec))

  out <- list(
    scores = scores_all,
    signal_status = st_y,
    task_rank = pmax(1L, task_rk_int),
    nuisance_rank = as.integer(rk_z),
    n_observations = m,
    comparable = TRUE
  )

  if (isTRUE(evidence)) {
    # Orthonormal basis for Qz_perp (m x df_nuis) so pairwise Gaussian marginal LR tests
    # evaluate exact scale-free two-chisquare F-ratios under H_{k0}: Y ~ N(0, sigma^2 Sigma_{k0})
    p_perp <- diag(m) - if (rk_z > 0L) tcrossprod(qz) else 0.0
    q_perp <- hrf_safe_svd(p_perp, nv = 0L)$u[, seq_len(df_nuis), drop = FALSE]
    u_perp_list <- eig_s_list <- vector("list", K)
    ld_pair_vec <- numeric(K)
    for (k in seq_len(K)) {
      uk <- cbind(qc_list[[k]], ud_list[[k]])
      u_perp_list[[k]] <- if (ncol(uk) > 0L) crossprod(q_perp, uk) else matrix(0.0, df_nuis, 0L)
      eig_s_list[[k]] <- c(rep(26.0, rc), 1.0 + gamma_pair * s2_norm_list[[k]])
      ld_pair_vec[k] <- sum(log(eig_s_list[[k]]))
    }

    eta_cache <- new.env(parent = emptyenv())
    get_pair_eta <- function(j_cand, k0_base) {
      key_p <- paste0(j_cand, "_", k0_base)
      if (exists(key_p, envir = eta_cache, inherits = FALSE)) {
        return(get(key_p, envir = eta_cache))
      }
      u0 <- u_perp_list[[k0_base]]; e0 <- eig_s_list[[k0_base]]
      uj <- u_perp_list[[j_cand]];  ej <- eig_s_list[[j_cand]]
      s0_half <- diag(df_nuis) + if (ncol(u0) > 0L) u0 %*% ((sqrt(e0) - 1.0) * t(u0)) else 0.0
      sj_inv  <- diag(df_nuis) - if (ncol(uj) > 0L) uj %*% (((ej - 1.0) / ej) * t(uj)) else 0.0
      mat_m <- s0_half %*% sj_inv %*% s0_half
      eta_v <- eigen(0.5 * (mat_m + t(mat_m)), symmetric = TRUE, only.values = TRUE)$values
      assign(key_p, eta_v, envir = eta_cache)
      eta_v
    }

    eval_pair_p <- function(j_cand, k0_base, v_sub = seq_len(n_ok)) {
      p_out <- rep(1.0, length(v_sub))
      if (j_cand == k0_base || !length(v_sub)) return(p_out)
      v_win <- v_sub[which(is.finite(scores_ok[j_cand, v_sub]) & is.finite(scores_ok[k0_base, v_sub]) &
                             (scores_ok[j_cand, v_sub] > scores_ok[k0_base, v_sub]) &
                             (ll_pair_mat[j_cand, v_sub] > ll_pair_mat[k0_base, v_sub]))]
      if (!length(v_win)) return(p_out)
      eta_v <- get_pair_eta(j_cand, k0_base)
      delta_jk <- (ld_pair_vec[k0_base] - ld_pair_vec[j_cand]) / as.numeric(df_nuis)
      t_vec <- ll_pair_mat[j_cand, v_win] - ll_pair_mat[k0_base, v_win]
      alpha_vec <- exp(delta_jk - 2.0 * t_vec / as.numeric(df_nuis))
      p_win <- vapply(alpha_vec, function(al) {
        pos <- al - eta_v[eta_v < al]
        neg <- eta_v[eta_v > al] - al
        if (!length(pos)) return(1.0)
        if (!length(neg)) return(0.0)
        mu_p <- sum(pos); v_p <- 2.0 * sum(pos * pos); df_p <- max(1.0, 2.0 * mu_p * mu_p / max(v_p, 1e-16))
        mu_n <- sum(neg); v_n <- 2.0 * sum(neg * neg); df_n <- max(1.0, 2.0 * mu_n * mu_n / max(v_n, 1e-16))
        stats::pf(mu_n / max(mu_p, 1e-16), df1 = df_p, df2 = df_n, lower.tail = FALSE)
      }, numeric(1L))
      p_out[match(v_win, v_sub)] <- p_win
      p_out
    }

    fb_p_ok <- matrix(1.0, K, n_ok)
    for (k in seq_len(K)) {
      if (k != fb_idx) fb_p_ok[k, ] <- eval_pair_p(k, fb_idx)
    }

    incomp_p_ok <- matrix(1.0, K, n_ok)
    if (K > 1L) {
      best_v <- max.col(t(ifelse(is.finite(scores_ok), scores_ok, -Inf)), ties.method = "first")
      for (u_best in unique(best_v)) {
        v_hit <- which(best_v == u_best)
        for (k in seq_len(K)) {
          if (k != u_best) {
            p_ub <- eval_pair_p(u_best, k, v_hit)
            incomp_p_ok[k, v_hit] <- pmin(incomp_p_ok[k, v_hit], p_ub)
          }
        }
      }
      for (k in seq_len(K)) {
        if (k != fb_idx) {
          v_fb_better <- which(scores_ok[fb_idx, ] > scores_ok[k, ])
          if (length(v_fb_better)) {
            incomp_p_ok[k, v_fb_better] <- pmin(incomp_p_ok[k, v_fb_better], eval_pair_p(fb_idx, k, v_fb_better))
          }
        }
      }
    }

    task_p_ok[, nuis_only] <- NA_real_
    fb_p_ok[, nuis_only] <- NA_real_
    incomp_p_ok[, nuis_only] <- NA_real_

    tp_all <- fb_all <- ic_all <- matrix(NA_real_, K, n_vox)
    tp_all[, ok_v] <- task_p_ok
    fb_all[, ok_v] <- fb_p_ok
    ic_all[, ok_v] <- incomp_p_ok
    out$task_p <- tp_all
    out$fallback_p <- fb_all
    out$incompatible_p <- ic_all
  }

  out
}

# Build stratum keys and per-voxel HRF kernel accessor for deconvolve_crossrun_hrf()
# and decon_epoch_regression(). In "hard" mode, preserves exact integer kernel_index
# strata; in "soft_average" mode, groups voxels that share identical composite weight
# vectors and returns their normalized composite HRF kernel.
hrf_stratum_kernels <- function(hrf_fit, r) {
  selected <- hrf_fit$selections[[r]]$kernel_index
  noise_r <- if (is.null(hrf_fit$noise)) NULL else hrf_fit$noise[[r]]
  rho <- if (is.null(noise_r)) rep(0, length(selected)) else noise_r$rho
  psi1 <- if (!is.null(noise_r) && "psi1" %in% names(noise_r)) noise_r$psi1 else rep(0, length(selected))
  psi2 <- if (!is.null(noise_r) && "psi2" %in% names(noise_r)) noise_r$psi2 else rep(0, length(selected))
  psi3 <- if (!is.null(noise_r) && "psi3" %in% names(noise_r)) noise_r$psi3 else rep(0, length(selected))
  psi4 <- if (!is.null(noise_r) && "psi4" %in% names(noise_r)) noise_r$psi4 else rep(0, length(selected))
  has_p34 <- abs(psi3) > 1e-10 | abs(psi4) > 1e-10
  has_p12 <- abs(psi1) > 1e-10 | abs(psi2) > 1e-10 | has_p34
  psi_sfx <- ifelse(
    has_p34,
    paste0(":", sprintf("%.6g", psi1), ":", sprintf("%.6g", psi2),
           ":", sprintf("%.6g", psi3), ":", sprintf("%.6g", psi4)),
    ifelse(
      has_p12,
      paste0(":", sprintf("%.6g", psi1), ":", sprintf("%.6g", psi2)),
      ""
    )
  )
  get_ar <- function(voxel_idx) {
    idx <- voxel_idx[1L]
    if (has_p34[idx]) {
      c(rho[idx], psi1[idx], psi2[idx], psi3[idx], psi4[idx])
    } else if (has_p12[idx]) {
      c(rho[idx], psi1[idx], psi2[idx])
    } else {
      rho[idx]
    }
  }
  is_soft <- identical(hrf_fit$settings$selection_mode, "soft_average") &&
    !is.null(hrf_fit$selected_kernels) && !is.null(hrf_fit$selected_kernels[[r]])
  if (!is_soft) {
    keys <- paste0(paste(selected, sprintf("%.17g", rho), sep = ":"), psi_sfx)
    get_kernel <- function(voxel_idx) {
      hrf_fit$kernels[, selected[voxel_idx], drop = FALSE]
    }
    return(list(
      keys = keys, rho = rho, psi1 = psi1, psi2 = psi2, psi3 = psi3, psi4 = psi4,
      selected = selected, is_soft = FALSE, get_kernel = get_kernel, get_ar = get_ar
    ))
  }
  w_sig <- if (!is.null(hrf_fit$hrf_strata) && !is.null(hrf_fit$hrf_strata[[r]]) &&
               length(hrf_fit$hrf_strata[[r]]) == length(selected)) {
    paste0("p", hrf_fit$hrf_strata[[r]])
  } else {
    w_mat <- if (!is.null(hrf_fit$kernel_weights) && !is.null(hrf_fit$kernel_weights[[r]])) {
      hrf_fit$kernel_weights[[r]]
    } else {
      hrf_fit$selected_kernels[[r]]
    }
    hrf_fast_col_keys(w_mat, digits = 4L)
  }
  keys <- paste0(paste(w_sig, sprintf("%.17g", rho), sep = ":"), psi_sfx)
  sel_k_mat <- hrf_fit$selected_kernels[[r]]
  get_kernel <- function(voxel_idx) {
    sel_k_mat[, voxel_idx, drop = FALSE]
  }
  list(
    keys = keys, rho = rho, psi1 = psi1, psi2 = psi2, psi3 = psi3, psi4 = psi4,
    selected = selected, is_soft = TRUE, get_kernel = get_kernel, get_ar = get_ar
  )
}

# Construct the exact gap-aware AR(1) or cascaded AR(1)xAR(q) Cholesky whitening
# matrix W on observed acquisition indices `rows`.
# - When `unit_marginal = FALSE` (default), W %*% Sigma_cor %*% t(W) = sig_u2 * I_m
#   (reducing identically to (1 - rho^2) * I_m under pure AR(1) when psi = 0).
# - When `unit_marginal = TRUE`, W %*% Sigma_cor %*% t(W) = I_m for exact GLS/REML
#   variance estimation.
# On contiguous grids (all(diff(rows) == 1)), builds the exact O(m * p) stationary
# AR(p) Cholesky factor with exact p x p Toeplitz initial block; across arbitrary
# scrub/censor gaps with higher-order psi != 0, builds the exact subset Toeplitz
# Cholesky factor so zero off-diagonal correlation leaks across censor boundaries.
hrf_whitening_matrix <- function(rows, rho, unit_marginal = FALSE) {
  m <- length(rows)
  w <- diag(m)
  if (!length(rho) || m <= 1L) return(w)
  has_psi <- length(rho) >= 2L && any(is.finite(rho[-1L]) & abs(rho[-1L]) > 1e-10)
  if (!has_psi) {
    rho1 <- rho[1L]
    if (!is.finite(rho1) || abs(rho1) <= 1e-12) return(w)
    rho1 <- min(max(rho1, -0.95), 0.95)
    s_u <- sqrt(max(1.0 - rho1 * rho1, 1e-8))
    w[1L, 1L] <- s_u
    for (i in seq.int(2L, m)) {
      a <- rho1^(rows[i] - rows[i - 1L])
      multiplier <- s_u / sqrt(max(1.0 - a * a, 1e-12))
      w[i, i] <- multiplier
      w[i, i - 1L] <- -a * multiplier
    }
    return(if (isTRUE(unit_marginal)) w / s_u else w)
  }

  is_contig <- all(diff(rows) == 1L)
  lag_span <- if (is_contig) max(8L, length(rho) + 2L) else max(8L, max(rows) - min(rows))
  info <- hrf_ar_cascade_info(rho, max_lag = lag_span)
  phi <- info$phi
  p <- length(phi)
  if (p == 0L) return(w)
  sig_u <- sqrt(info$sig_u2)

  if (is_contig && m > p) {
    w_unit <- matrix(0.0, m, m)
    Gamma_p <- matrix(info$acf[abs(outer(seq_len(p), seq_len(p), "-")) + 1L], p, p)
    L_p <- tryCatch(t(chol(Gamma_p)), error = function(e) NULL)
    if (!is.null(L_p)) {
      w_unit[seq_len(p), seq_len(p)] <- forwardsolve(L_p, diag(p))
      inv_su <- 1.0 / sig_u
      for (i in seq.int(p + 1L, m)) {
        w_unit[i, i] <- inv_su
        for (k in seq_len(p)) {
          w_unit[i, i - k] <- -phi[k] * inv_su
        }
      }
      return(if (isTRUE(unit_marginal)) w_unit else w_unit * sig_u)
    }
  }

  idx <- abs(outer(rows, rows, "-")) + 1L
  Sigma_sub <- matrix(info$acf[idx], m, m)
  R_chol <- tryCatch(chol(Sigma_sub), error = function(e) {
    diag(Sigma_sub) <- diag(Sigma_sub) + 1e-6
    chol(Sigma_sub)
  })
  w_unit <- t(backsolve(R_chol, diag(m)))
  if (isTRUE(unit_marginal)) w_unit else w_unit * sig_u
}

# Apply exact AR(1) or cascaded AR(p) Cholesky whitening to a matrix (or list of
# matrices) on observed acquisition indices `rows` via O(m * p) C++ ArWhitener,
# matching `hrf_whitening_matrix(rows, ar_vec, unit_marginal) %*% X` without
# allocating a dense m x m matrix in R.
hrf_apply_whitening <- function(target, rows, ar_vec, unit_marginal = FALSE) {
  rows <- as.integer(rows)
  m <- length(rows)
  as_mat_target <- function(obj) {
    if (is.list(obj)) {
      lapply(obj, function(el) if (is.null(el)) NULL else as.matrix(el))
    } else {
      as.matrix(obj)
    }
  }
  if (m <= 1L || !length(ar_vec)) return(as_mat_target(target))
  has_psi <- length(ar_vec) >= 2L && any(is.finite(ar_vec[-1L]) & abs(ar_vec[-1L]) > 1e-10)
  if (!has_psi) {
    rho1 <- ar_vec[1L]
    if (!is.finite(rho1) || abs(rho1) <= 1e-12) return(as_mat_target(target))
    rho1 <- min(max(rho1, -0.95), 0.95)
    return(hrf_whiten_apply_cpp(
      target = as_mat_target(target), rows = rows, rho = rho1,
      unit_marginal = isTRUE(unit_marginal)
    ))
  }
  is_contig <- all(diff(rows) == 1L)
  lag_span <- if (is_contig) max(8L, length(ar_vec) + 2L) else max(8L, max(rows) - min(rows))
  info <- hrf_ar_cascade_info(ar_vec, max_lag = lag_span)
  if (!length(info$phi)) return(as_mat_target(target))
  hrf_whiten_apply_cpp(
    target = as_mat_target(target), rows = rows, rho = 0.0,
    ar_phi = info$phi, ar_acf = info$acf, sig_u2 = info$sig_u2,
    unit_marginal = isTRUE(unit_marginal)
  )
}

hrf_score_run_ar <- function(y_sub, design, z, rows, batch_size, ar_vec,
    evidence = FALSE, fallback_index = 1L) {
  rho1 <- if (length(ar_vec)) ar_vec[1L] else 0.0
  has_psi <- length(ar_vec) >= 2L && any(is.finite(ar_vec[-1L]) & abs(ar_vec[-1L]) > 1e-10)
  if (!has_psi || length(rows) < 3L) {
    return(hrf_score_run_cpp(y_sub, design, z, rows, as.integer(batch_size),
      rho = rho1, evidence = evidence, fallback_index = as.integer(fallback_index)))
  }
  is_contig <- all(diff(rows) == 1L)
  lag_span <- if (is_contig) max(8L, length(ar_vec) + 2L) else max(8L, max(rows) - min(rows))
  info <- hrf_ar_cascade_info(ar_vec, max_lag = lag_span)
  if (!length(info$phi)) {
    return(hrf_score_run_cpp(y_sub, design, z, rows, as.integer(batch_size),
      rho = 0.0, evidence = evidence, fallback_index = as.integer(fallback_index)))
  }
  hrf_score_run_cpp(y_sub, design, z, rows, as.integer(batch_size),
    rho = 0.0, evidence = evidence, fallback_index = as.integer(fallback_index),
    ar_phi = info$phi, ar_acf = info$acf, sig_u2 = info$sig_u2)
}

# Estimate run-level temporal observation weights w in (0, 1] on `rows` by
# pooling standardized task+trial-projected residuals across usable voxels.
# Uses Winsorized-SVD Liu-Tang-Zhang (2009) 3-cumulant quadratic-form null
# calibration (exact under arbitrary spatial covariance / global signal modes)
# paired with Sparse Bayesian Learning (Tipping 2001) Woodbury coordinate updates
# for additive transient excess variance Sigma = Sigma_AR + diag((1 - w) / w),
# or Diedrichsen & Shadmehr (2005) innovation-space Empirical Bayes weights.
hrf_estimate_rwls_weights <- function(z_raw, r_obs = NULL, rows, ar_vec = 0.0,
                                      mode = c("auto", "woodbury", "innovation"),
                                      task_env = NULL) {
  mode <- match.arg(mode)
  rows <- as.integer(rows)
  m <- length(rows)
  if (!is.matrix(z_raw) || nrow(z_raw) != m) return(rep(1.0, m))
  ok_v <- which(colSums(!is.finite(z_raw)) == 0L &
    apply(z_raw, 2L, function(x) sum(abs(x)) > .Machine$double.eps))
  v_ok <- length(ok_v)
  min_v <- if (mode == "auto") 6L else 3L
  min_m <- if (mode == "auto") 24L else 16L
  if (m < min_m || v_ok < min_v) return(rep(1.0, m))

  Z_std <- z_raw[, ok_v, drop = FALSE]
  has_robs <- is.matrix(r_obs) && nrow(r_obs) == m && ncol(r_obs) == ncol(z_raw)
  R_obs <- if (has_robs) r_obs[, ok_v, drop = FALSE] else NULL
  if (!is.null(R_obs) && any(!is.finite(R_obs))) R_obs <- NULL

  # Guard against systematic task-locked HRF shape/latency mismatch or unmodeled
  # single-trial neural variance across trials masquerading as transient noise spikes
  # (Diedrichsen & Shadmehr 2005 Sec. 3.3). Comparing the 75th percentile of task vs
  # baseline frames has a 25% breakdown point in task windows (preserving full power
  # for genuine isolated motion glitches even when they collide with task trials)
  # while rescaling systematic multi-trial HRF mismatch.
  if (!is.null(task_env) && length(task_env) == m && all(is.finite(task_env))) {
    te <- pmax(as.numeric(task_env), 0.0)
    te_ref <- unname(stats::quantile(te, probs = 0.90, names = FALSE))
    if (is.finite(te_ref) && te_ref > 1e-12) {
      te <- pmin(te / te_ref, 2.5)
      hi <- which(te >= 0.20)
      lo <- which(te < 0.08)
      if (length(hi) < 8L || length(lo) < 8L) {
        med_te <- stats::median(te)
        hi <- which(te > med_te)
        lo <- which(te <= med_te)
      }
      if (length(hi) >= 8L && length(lo) >= 8L) {
        qin0 <- rowMeans(pmin(Z_std^2, 100.0))
        q_hi <- unname(stats::quantile(qin0[hi], probs = 0.75, names = FALSE))
        q_lo <- max(unname(stats::quantile(qin0[lo], probs = 0.75, names = FALSE)), 1e-6)
        if (is.finite(q_hi) && (q_hi / q_lo) > 1.30) {
          gam_in <- (q_hi / q_lo) - 1.10
          Z_std <- sweep(Z_std, 1L, sqrt(1.0 + gam_in * te), "/")
        }
        if (!is.null(R_obs)) {
          qob0 <- rowMeans(pmin(R_obs^2, 100.0))
          qb_hi <- unname(stats::quantile(qob0[hi], probs = 0.75, names = FALSE))
          qb_lo <- max(unname(stats::quantile(qob0[lo], probs = 0.75, names = FALSE)), 1e-6)
          if (is.finite(qb_hi) && (qb_hi / qb_lo) > 1.30) {
            gam_ob <- (qb_hi / qb_lo) - 1.10
            R_obs <- sweep(R_obs, 1L, sqrt(1.0 + gam_ob * te), "/")
          }
        }
      }
    }
  }

  cw <- 2.8
  bw <- (2.0 * stats::pnorm(cw) - 1.0) - 2.0 * cw * stats::dnorm(cw) +
    2.0 * (cw^2) * (1.0 - stats::pnorm(cw))
  sub_cols <- if (v_ok > 512L) unique(round(seq(1, v_ok, length.out = 512L))) else seq_len(v_ok)

  ltz_null <- function(M_mat) {
    if (is.null(M_mat)) return(NULL)
    M_sub <- M_mat[, sub_cols, drop = FALSE]
    Mw <- pmin(pmax(M_sub, -cw), cw) / sqrt(bw)
    sv <- tryCatch(hrf_safe_svd(Mw, nu = 0, nv = 0)$d, error = function(e) numeric(0))
    ev <- (sv^2) / (m * ncol(Mw))
    ev <- ev[is.finite(ev) & ev > 1e-12]
    if (!length(ev)) return(NULL)
    s1 <- sum(ev); s2 <- sum(ev^2); s3 <- sum(ev^3)
    c1 <- s3 / max(s2, 1e-12)
    nu_s <- max((s2^3) / max(s3^2, 1e-18), 0.5)
    c0 <- s1 - c1 * nu_s
    nu_eff <- min(max((s1^2) / max(s2, 1e-12), 1.0), as.numeric(v_ok))
    list(
      s1 = s1, c0 = c0, c1 = c1, nu_s = nu_s, nu_eff = nu_eff,
      qfun = function(p) c0 + c1 * stats::qchisq(p, df = nu_s)
    )
  }

  ltz_in <- ltz_null(Z_std)
  if (is.null(ltz_in) || !is.finite(ltz_in$s1) || ltz_in$s1 <= 0) return(rep(1.0, m))
  q_in <- rowMeans(pmin(Z_std^2, 100.0))

  if (mode == "auto") {
    gate_thr <- max(ltz_in$qfun(1.0 - 0.01 / m), 1.8 * ltz_in$s1)
    if (max(q_in) <= gate_thr) return(rep(1.0, m))
  } else {
    if (max(q_in) <= ltz_in$qfun(0.95)) return(rep(1.0, m))
  }

  compute_innov_weights <- function() {
    v_norm <- q_in / max(ltz_in$s1, 0.1)
    thr_98 <- max(ltz_in$qfun(0.98) / max(ltz_in$s1, 0.1), 1.15)
    nu_e <- ltz_in$nu_eff
    nu0_in <- 4.0
    v_eb <- (nu0_in + nu_e * v_norm) / (nu0_in + nu_e)
    gate <- pmin(pmax((v_norm - 0.85 * thr_98) / max(0.30 * thr_98, 0.15), 0.0), 1.0)
    v_innov <- 1.0 + gate * pmax(v_eb - 1.0, 0.0)
    pmin(pmax(1.0 / v_innov, 0.005), 1.0)
  }

  if (mode == "innovation") {
    return(compute_innov_weights())
  }

  q_ob <- if (!is.null(R_obs)) rowMeans(pmin(R_obs^2, 100.0)) else rep(0.0, m)
  ltz_ob <- ltz_null(R_obs)
  thr_in <- ltz_in$qfun(0.985)
  thr_ob <- if (!is.null(ltz_ob)) ltz_ob$qfun(0.985) else Inf
  cand <- which(q_in > thr_in | q_ob > thr_ob)
  if (!length(cand)) return(rep(1.0, m))

  score_cand <- pmax(
    q_ob[cand] / max(if (!is.null(ltz_ob)) ltz_ob$s1 else 1.0, 0.1),
    q_in[cand] / max(ltz_in$s1, 0.1)
  )
  cand <- cand[order(score_cand, decreasing = TRUE)]
  max_cand <- max(16L, as.integer(0.35 * m))
  if (length(cand) > max_cand) cand <- cand[seq_len(max_cand)]

  Z_norm <- Z_std / sqrt(max(ltz_in$s1, 0.1))
  k <- length(cand)
  E_c <- matrix(0.0, m, k)
  E_c[cbind(cand, seq_len(k))] <- 1.0
  F_c <- hrf_apply_whitening(E_c, rows, ar_vec, unit_marginal = TRUE)
  G_f <- crossprod(F_c)
  FZ  <- crossprod(F_c, Z_norm)

  nu_eff <- ltz_in$nu_eff
  d_vec <- rep(0.0, k)
  nu0 <- 3.0
  eb_shrink <- nu_eff / (nu0 + nu_eff)
  q_mult <- max(ltz_in$qfun(0.985) / max(ltz_in$s1, 0.1), 1.10)

  for (sweep_i in 1:3) {
    for (r in seq_len(k)) {
      act <- which(d_vec > 1e-6 & seq_len(k) != r)
      if (!length(act)) {
        q_mr <- max(G_f[r, r], 1e-8)
        p_mr <- mean(FZ[r, ]^2)
      } else {
        M_a <- G_f[act, act, drop = FALSE] + diag(1.0 / d_vec[act], length(act))
        g_ra <- G_f[act, r]
        w_ra <- tryCatch(solve(M_a, g_ra), error = function(e) rep(0.0, length(act)))
        q_mr <- max(G_f[r, r] - sum(g_ra * w_ra), 1e-8)
        proj_z <- FZ[r, ] - as.numeric(crossprod(w_ra, FZ[act, , drop = FALSE]))
        p_mr <- mean(proj_z^2)
      }
      if (is.finite(p_mr) && p_mr > q_mr * q_mult) {
        d_ml <- (p_mr - q_mr) / (q_mr * q_mr)
        d_vec[r] <- min(max(eb_shrink * d_ml, 0.0), 199.0)
      } else {
        d_vec[r] <- 0.0
      }
    }
  }

  delta2 <- rep(0.0, m)
  delta2[cand] <- d_vec
  w_wood <- pmin(pmax(1.0 / (1.0 + delta2), 0.005), 1.0)

  if (mode == "auto" && !is.null(ltz_ob) && k >= 1L) {
    # Distinguish isolated observation glitches (where W_unit amplifies spike energy
    # by g_imp = ||W_unit e_t||_2^2 = (1 + sum(phi^2)) / sigma_u^2 > 1) from sustained
    # recursive physiological/respiratory innovation bursts (where excess energy in
    # observation space is comparable to innovation space: ex_ratio << g_imp).
    g_imp <- mean(diag(G_f))
    if (is.finite(g_imp) && g_imp > 1.25) {
      ex_in <- sum(pmax(q_in / max(ltz_in$s1, 0.1) - 1.5, 0.0))
      ex_ob <- sum(pmax(q_ob / max(ltz_ob$s1, 0.1) - 1.5, 0.0))
      if (ex_ob > 1e-6 && (ex_in / ex_ob) < 0.72 * g_imp) {
        w_innov <- compute_innov_weights()
        if (!any(abs(w_innov - 1.0) > 1e-8)) return(rep(1.0, m))
        return(structure(w_innov, weight_space = "innovation"))
      }
    }
  }
  w_wood
}

hrf_inverse_operator <- function(n, kernel, rows, rho, nuisance, settings, obs_weights = NULL, lambda = NULL) {
  measured <- length(rows)
  lam_use <- if (is.null(lambda)) as.numeric(settings$lambda)[1L] else as.numeric(lambda)[1L]
  w_space <- if (!is.null(obs_weights) && !is.null(attr(obs_weights, "weight_space"))) {
    as.character(attr(obs_weights, "weight_space"))[1L]
  } else if (!is.null(settings$rwls_space)) {
    settings$rwls_space
  } else {
    "woodbury"
  }
  has_obs_w <- !is.null(obs_weights) && length(obs_weights) == measured &&
    any(abs(obs_weights - 1.0) > 1e-8)
  w_arg <- if (has_obs_w) pmin(pmax(as.numeric(obs_weights), 1e-8), 1.0) else NULL

  # Pre-projection mode: orthogonalize observed BOLD against nuisance/drift on `rows`
  # without stripping smooth low-frequency neural modes from the continuous forward matrix H.
  if (identical(settings$nuisance_mode, "preproject") && ncol(nuisance) > 0L) {
    settings_base <- settings
    settings_base$nuisance_mode <- "ignore"
    op_base <- hrf_inverse_operator(
      n = n, kernel = kernel, rows = rows, rho = rho,
      nuisance = matrix(numeric(), n, 0L), settings = settings_base,
      obs_weights = obs_weights, lambda = lambda
    )
    z_obs <- nuisance[rows, , drop = FALSE]
    sw <- if (has_obs_w) sqrt(pmax(w_arg, 1e-4)) else NULL
    z_w <- if (has_obs_w) sweep(z_obs, 1L, sw, "*") else z_obs
    nrm_z <- sqrt(colSums(z_w^2)); nrm_z[nrm_z == 0] <- 1.0
    sv_z <- hrf_safe_svd(sweep(z_w, 2L, nrm_z, "/"), nv = 0L)
    rk_z <- if (length(sv_z$d) && sv_z$d[1L] > 0) sum(sv_z$d > 1e-8 * sv_z$d[1L]) else 0L
    if (rk_z > 0L && measured - rk_z >= 8L) {
      qz <- sv_z$u[, seq_len(rk_z), drop = FALSE]
      if (has_obs_w) {
        op_qz <- op_base %*% sweep(qz, 1L, sw, "/")
        op_base <- op_base - op_qz %*% sweep(t(qz), 2L, sw, "*")
      } else {
        op_base <- op_base - (op_base %*% qz) %*% t(qz)
      }
    }
    return(op_base)
  }

  has_psi <- length(rho) >= 2L && any(is.finite(rho[-1L]) & abs(rho[-1L]) > 1e-10)
  if (settings$backend == "cpp") {
    if (!has_psi) {
      return(hrf_inverse_operator_cpp(as.integer(n), kernel,
        rows, rho[1L], nuisance, lam_use, settings$penalty, settings$ridge_floor,
        as.integer(settings$tail_padding),
        obs_weights = w_arg, weight_space = w_space))
    }
    all_rows <- c(rows, if (settings$tail_padding) n + seq_len(settings$tail_padding))
    is_contig <- all(diff(all_rows) == 1L)
    lag_span <- if (is_contig) max(8L, length(rho) + 2L) else max(8L, max(all_rows) - min(all_rows))
    info <- hrf_ar_cascade_info(rho, max_lag = lag_span)
    if (!length(info$phi)) {
      return(hrf_inverse_operator_cpp(as.integer(n), kernel,
        rows, 0.0, nuisance, lam_use, settings$penalty, settings$ridge_floor,
        as.integer(settings$tail_padding),
        obs_weights = w_arg, weight_space = w_space))
    }
    return(hrf_inverse_operator_cpp(as.integer(n), kernel,
      rows, 0.0, nuisance, lam_use, settings$penalty, settings$ridge_floor,
      as.integer(settings$tail_padding),
      ar_phi = info$phi, ar_acf = info$acf, sig_u2 = info$sig_u2,
      obs_weights = w_arg, weight_space = w_space))
  }
  # Exact masked, nuisance-adjusted GLS reference (used for backend = "r").
  all_rows <- c(rows, if (settings$tail_padding) n + seq_len(settings$tail_padding))
  h <- decon_convolution_matrix(n + settings$tail_padding, kernel)[all_rows, , drop = FALSE]
  if (!has_obs_w) {
    w <- hrf_whitening_matrix(all_rows, rho)
  } else if (identical(w_space, "innovation")) {
    w_std <- hrf_whitening_matrix(all_rows, rho)
    w_all <- c(w_arg, rep(1.0, settings$tail_padding))
    w <- sweep(w_std, 1L, sqrt(w_all), "*")
  } else {
    w_u <- hrf_whitening_matrix(all_rows, rho, unit_marginal = TRUE)
    w_s <- hrf_whitening_matrix(all_rows, rho, unit_marginal = FALSE)
    s0 <- if (nrow(w_u) > 0L && abs(w_u[1L, 1L]) > 1e-12) w_s[1L, 1L] / w_u[1L, 1L] else 1.0
    delta2 <- c((1.0 - w_arg) / w_arg, rep(0.0, settings$tail_padding))
    sigma_sub <- solve(crossprod(w_u))
    diag(sigma_sub) <- diag(sigma_sub) + delta2
    r_rwls <- chol(sigma_sub)
    w <- s0 * t(backsolve(r_rwls, diag(nrow(sigma_sub))))
  }
  z <- rbind(nuisance[rows, , drop = FALSE], matrix(0, settings$tail_padding, ncol(nuisance)))
  z <- w %*% z; h <- w %*% h
  rhs <- w[, seq_len(measured), drop = FALSE]
  if (ncol(z)) {
    norm <- sqrt(colSums(z^2)); norm[norm == 0] <- 1
    fit <- hrf_safe_svd(sweep(z, 2L, norm, "/"), nv = 0)
    rank <- sum(fit$d > 1e-8 * fit$d[1L])
    if (rank) {
      q <- fit$u[, seq_len(rank), drop = FALSE]
      h <- h - q %*% crossprod(q, h); rhs <- rhs - q %*% crossprod(q, rhs)
    }
  }
  lhs <- crossprod(h) + lam_use * decon_penalty_crossprod(ncol(h), settings$penalty)
  diag(lhs) <- diag(lhs) + settings$ridge_floor
  solve(lhs, crossprod(h, rhs))[length(kernel) - 1L + seq_len(n), , drop = FALSE]
}

# Weighted Pool-Adjacent-Violators (PAVA) for monotone non-increasing sequence y_1 >= y_2 >= ... >= y_G.
hrf_pava_nonincreasing <- function(y, w = NULL) {
  n <- length(y)
  if (n <= 1L) return(as.numeric(y))
  if (is.null(w)) w <- rep(1.0, n)
  w <- pmax(as.numeric(w), 1e-8)
  # Non-increasing y is equivalent to non-decreasing -y; apply standard PAVA on -y.
  z <- -as.numeric(y)
  blk_val <- numeric(n)
  blk_w <- numeric(n)
  blk_len <- integer(n)
  m <- 0L
  for (i in seq_len(n)) {
    m <- m + 1L
    blk_val[m] <- z[i]
    blk_w[m] <- w[i]
    blk_len[m] <- 1L
    while (m >= 2L && blk_val[m - 1L] > blk_val[m]) {
      w_tot <- blk_w[m - 1L] + blk_w[m]
      blk_val[m - 1L] <- (blk_w[m - 1L] * blk_val[m - 1L] + blk_w[m] * blk_val[m]) / w_tot
      blk_w[m - 1L] <- w_tot
      blk_len[m - 1L] <- blk_len[m - 1L] + blk_len[m]
      m <- m - 1L
    }
  }
  out <- numeric(n)
  pos <- 0L
  for (b in seq_len(m)) {
    len_b <- blk_len[b]
    out[(pos + 1L):(pos + len_b)] <- -blk_val[b]
    pos <- pos + len_b
  }
  out
}

# Gain-calibrated, task-weighted adaptive prewhitened GCV selection for the
# deconvolution penalty lambda in deconvolve_crossrun_hrf().
#
# Standard observation-space GCV minimizes prediction error ||W H (s_hat - s)||_2^2
# in convolved BOLD space, which severely over-smooths event-related MEDuSA
# deconvolution for three reasons:
#   1. Unfiltered null/low-SNR voxels have strictly decreasing E[GCV(lambda)],
#      pulling a global lambda toward the upper grid boundary.
#   2. Because the passband singular values of W H are small (~0.15-0.30),
#      shrinking s_hat -> 0 reduces observation-space noise faster than signal
#      loss whenever per-TR event power is smaller than per-TR noise variance.
#   3. In Stage 3 (decon_epoch_regression), physical scale is restored by
#      dividing by the diagonal resolution gain g_lvl(lambda) = b_0' h_0.
#      Consequently, the gain-restored variance of both a single-lag readout
#      V_lvl(lambda) = ||b_0||_2^2 / g_lvl(lambda)^2 and a 1-TR temporal contrast
#      V_diff(lambda) = ||b_0 - b_1||_2^2 / g_diff(lambda)^2 are U-shaped in
#      lambda: beyond a finite variance minimum lambda_var*(TR, rho, penalty),
#      increasing lambda increases BOTH smoothing bias AND gain-restored variance.
#   4. Discrete derivative penalties scale with TR^4 (diff2/sobolev) and AR(1)
#      whitening scales signal energy by (1 - rho)^2, so the optimal penalty
#      shifts systematically across acquisition TRs and noise autocorrelation.
#
# This helper pools task-supported voxels/parcels from hrf_fit$selections,
# evaluates the prewhitened GCV trace alongside the exact interior gain-restored
# level/contrast variances and relative modulation transfer function C_rel(lambda)
# on a compact local window, and selects the lambda that minimizes the combined
# gain-calibrated prewhitened GCV criterion gcv_adaptive.
hrf_tune_decon_lambda <- function(bold, hrf_fit, lambda_grid, scaling = "zscore",
    penalty = "sobolev", ridge_floor = 1e-8, observation_mask = "fit",
    backend = "cpp", gcv_rule = c("min", "1se"), nuisance_mode = "ignore", tail_padding = 0L,
    frame_weights = NULL, precomputed_nsd = NULL,
    lambda_mode = c("auto", "global", "stratified"), max_lambda_strata = 4L,
    voxel_map = NULL, geometry = NULL, drift_active = TRUE,
    noise_pcs_list = NULL) {
  gcv_rule <- match.arg(gcv_rule)
  lambda_mode <- match.arg(lambda_mode)
  max_lambda_strata <- as.integer(max_lambda_strata)
  lambda_grid <- sort(unique(as.numeric(lambda_grid)))
  n_runs <- length(bold)
  n_vox <- ncol(bold[[1L]])

  if (length(lambda_grid) == 1L) {
    lam1 <- lambda_grid[1L]
    return(list(
      lambda = lam1,
      tuning = data.frame(
        kernel = 1L, lambda = lam1, df = NA_real_, rss = NA_real_,
        gcv = NA_real_, gcv_se = NA_real_, var_level = NA_real_,
        var_contrast = NA_real_, mtf_rel = 1.0, gcv_adaptive = 1.0, snr_task = NA_real_
      ),
      lambda_by_voxel = rep(lam1, n_vox),
      lambda_strata = data.frame(
        tier = 1L, label = "global", n_voxels = n_vox,
        mean_r2 = NA_real_, snr_task = NA_real_, lambda = lam1,
        stringsAsFactors = FALSE
      ),
      lambda_mode_used = "global"
    ))
  }

  r2_mat <- matrix(0.0, nrow = n_vox, ncol = n_runs)
  r2_sel_mat <- matrix(0.0, nrow = n_vox, ncol = n_runs)
  tp_mat <- matrix(1.0, nrow = n_vox, ncol = n_runs)
  has_r2 <- FALSE
  if (!is.null(hrf_fit$selections)) {
    for (r in seq_len(min(n_runs, length(hrf_fit$selections)))) {
      sel_r <- hrf_fit$selections[[r]]
      if (is.data.frame(sel_r) && nrow(sel_r) == n_vox) {
        if ("best_adjusted_partial_r2" %in% names(sel_r)) {
          v_r2 <- as.numeric(sel_r$best_adjusted_partial_r2)
          v_r2[!is.finite(v_r2)] <- 0.0
          r2_mat[, r] <- pmax(0.0, v_r2)
          r2_sel_mat[, r] <- r2_mat[, r]
          has_r2 <- TRUE
        }
        if ("selected_adjusted_partial_r2" %in% names(sel_r)) {
          v_sr2 <- as.numeric(sel_r$selected_adjusted_partial_r2)
          v_sr2[!is.finite(v_sr2)] <- 0.0
          r2_sel_mat[, r] <- pmax(0.0, v_sr2)
          has_r2 <- TRUE
        }
        if ("selected_task_p" %in% names(sel_r)) {
          v_tp <- as.numeric(sel_r$selected_task_p)
          v_tp[!is.finite(v_tp)] <- 1.0
          tp_mat[, r] <- pmin(1.0, pmax(0.0, v_tp))
        }
      }
    }
  }
  r2_v <- rowMeans(r2_mat, na.rm = TRUE)
  r2_v[!is.finite(r2_v)] <- 0.0
  tp_v <- apply(tp_mat, 1L, min, na.rm = TRUE)
  tp_v[!is.finite(tp_v)] <- 1.0

  tr <- if (!is.null(hrf_fit$settings$tr) && is.finite(hrf_fit$settings$tr) && hrf_fit$settings$tr > 0) {
    as.numeric(hrf_fit$settings$tr)
  } else {
    1.0
  }

  # Pre-compute valid rows and ok signal columns per run once
  rows_by_run <- vector("list", n_runs)
  ok_by_run <- vector("list", n_runs)
  for (r_tune in seq_len(n_runs)) {
    y_t <- bold[[r_tune]]
    n_t <- nrow(y_t)
    mask_t <- if (observation_mask == "fit" && !is.null(hrf_fit$specification$valid_rows)) {
      hrf_fit$specification$valid_rows[[r_tune]]
    } else {
      rep(TRUE, n_t)
    }
    rows_t <- which(!is.na(mask_t) & mask_t)
    rows_by_run[[r_tune]] <- rows_t
    if (length(rows_t) >= 3L) {
      st_t <- decon_signal_status(y_t[rows_t, , drop = FALSE])
      ok_by_run[[r_tune]] <- which(st_t == "ok")
    } else {
      ok_by_run[[r_tune]] <- integer(0L)
    }
  }

  tune_subset_cols <- function(cand_cols = NULL, mode_type = "global_active", r2_tier_score = NULL) {
    run_tuning_list <- list()
    run_weights <- numeric(0)
    wiener_lambdas <- numeric(0)
    wiener_tunings <- list()

    for (r_tune in seq_len(n_runs)) {
      y_t <- bold[[r_tune]]
      n_t <- nrow(y_t)
      rows_t <- rows_by_run[[r_tune]]
      ok_cols <- ok_by_run[[r_tune]]
      if (length(rows_t) < 3L || !length(ok_cols)) next

      r2_r_vec <- if (has_r2 && ncol(r2_mat) >= r_tune) r2_mat[, r_tune] else r2_v
      tp_r_vec <- if (ncol(tp_mat) >= r_tune) tp_mat[, r_tune] else tp_v

      if (identical(mode_type, "global_active")) {
        r2_ok <- r2_r_vec[ok_cols]
        tp_ok <- tp_r_vec[ok_cols]
        pos_r2 <- r2_ok[r2_ok > 0]
        if (length(pos_r2)) {
          r2_cut <- max(0.005, as.numeric(stats::quantile(pos_r2, probs = 0.35, na.rm = TRUE)))
          active_cols <- ok_cols[r2_ok >= r2_cut | tp_ok <= 0.15]
        } else {
          active_cols <- integer(0)
        }
        if (length(active_cols) < min(8L, length(ok_cols))) {
          ord_ok <- order(r2_ok, -tp_ok, decreasing = TRUE)
          active_cols <- ok_cols[utils::head(ord_ok, min(64L, length(ok_cols)))]
        } else if (length(active_cols) > 128L) {
          ord_act <- order(r2_r_vec[active_cols], decreasing = TRUE)
          active_cols <- active_cols[utils::head(ord_act, 128L)]
        }
      } else {
        sub_ok <- intersect(as.integer(cand_cols), ok_cols)
        if (!length(sub_ok)) next
        r2_ok <- r2_r_vec[sub_ok]
        if (length(sub_ok) > 96L) {
          if (identical(mode_type, "tier_background")) {
            # Evenly stride across background voxels for spatially balanced representation
            stride_idx <- unique(as.integer(round(seq(1, length(sub_ok), length.out = 96L))))
            active_cols <- sub_ok[stride_idx]
          } else {
            sc_sub <- if (!is.null(r2_tier_score)) r2_tier_score[sub_ok] else r2_r_vec[sub_ok]
            ord_sub <- order(sc_sub, decreasing = TRUE)
            active_cols <- sub_ok[utils::head(ord_sub, 96L)]
          }
        } else {
          active_cols <- sub_ok
        }
      }

      w_act <- if (identical(mode_type, "tier_background")) {
        rep(1.0 / length(active_cols), length(active_cols))
      } else {
        wa <- pmax(r2_r_vec[active_cols], 1e-3)
        wa / sum(wa)
      }
      sel_t <- hrf_fit$selections[[r_tune]]$kernel_index[active_cols]
      k_tab <- tapply(w_act, sel_t, sum)
      k_mode <- as.integer(names(k_tab)[which.max(k_tab)])
      if (!length(k_mode) || is.na(k_mode)) k_mode <- 1L

      norm_mode_tune <- if (!is.null(hrf_fit$settings$kernel_normalize)) {
        hrf_fit$settings$kernel_normalize
      } else {
        "signed_sum"
      }
      k_tune <- if (identical(hrf_fit$settings$selection_mode, "soft_average") &&
          !is.null(hrf_fit$selected_kernels) && !is.null(hrf_fit$selected_kernels[[r_tune]])) {
        k_avg <- as.numeric(hrf_fit$selected_kernels[[r_tune]][, active_cols, drop = FALSE] %*% w_act)
        tryCatch(
          normalize_decon_kernel(k_avg, mode = norm_mode_tune),
          error = function(e) as.numeric(hrf_fit$kernels[, k_mode])
        )
      } else {
        as.numeric(hrf_fit$kernels[, k_mode])
      }
      attr(k_tune, "normalize") <- norm_mode_tune

      noise_rt <- if (is.null(hrf_fit$noise)) NULL else hrf_fit$noise[[r_tune]]
      rho_t <- if (is.null(noise_rt)) 0.0 else stats::median(noise_rt$rho[active_cols], na.rm = TRUE)
      if (!is.finite(rho_t)) rho_t <- 0.0
      psi1_t <- if (!is.null(noise_rt) && "psi1" %in% names(noise_rt)) stats::median(noise_rt$psi1[active_cols], na.rm = TRUE) else 0.0
      psi2_t <- if (!is.null(noise_rt) && "psi2" %in% names(noise_rt)) stats::median(noise_rt$psi2[active_cols], na.rm = TRUE) else 0.0
      psi3_t <- if (!is.null(noise_rt) && "psi3" %in% names(noise_rt)) stats::median(noise_rt$psi3[active_cols], na.rm = TRUE) else 0.0
      psi4_t <- if (!is.null(noise_rt) && "psi4" %in% names(noise_rt)) stats::median(noise_rt$psi4[active_cols], na.rm = TRUE) else 0.0
      if (!is.finite(psi1_t)) psi1_t <- 0.0
      if (!is.finite(psi2_t)) psi2_t <- 0.0
      if (!is.finite(psi3_t)) psi3_t <- 0.0
      if (!is.finite(psi4_t)) psi4_t <- 0.0
      ar_t <- if (abs(psi3_t) > 1e-10 || abs(psi4_t) > 1e-10) {
        c(rho_t, psi1_t, psi2_t, psi3_t, psi4_t)
      } else if (abs(psi1_t) > 1e-10 || abs(psi2_t) > 1e-10) {
        c(rho_t, psi1_t, psi2_t)
      } else {
        rho_t
      }

      w_rt <- if (!is.null(frame_weights) && length(frame_weights) >= r_tune &&
                  length(frame_weights[[r_tune]]) == length(rows_t)) {
        as.numeric(frame_weights[[r_tune]])
      } else {
        NULL
      }
      has_w_rt <- !is.null(w_rt) && any(is.finite(w_rt)) && any(w_rt < 1.0 - 1e-8)

      nsd_tune <- NULL
      if (identical(scaling, "noise_sd")) {
        if (!is.null(precomputed_nsd) && length(precomputed_nsd) >= r_tune &&
            length(precomputed_nsd[[r_tune]]) >= max(active_cols)) {
          cand_nsd <- as.numeric(precomputed_nsd[[r_tune]][active_cols])
          if (all(is.finite(cand_nsd) & cand_nsd > 0)) {
            nsd_tune <- cand_nsd
          }
        }
        if (is.null(nsd_tune) && !is.null(hrf_fit$specification$events[[r_tune]])) {
          mask_tune_bool <- seq_len(n_t) %in% rows_t
          nuis_raw_t <- if (!is.null(hrf_fit$specification$nuisance) && length(hrf_fit$specification$nuisance) >= r_tune) {
            hrf_fit$specification$nuisance[[r_tune]]
          } else NULL
          pcs_t <- if (!is.null(noise_pcs_list) && length(noise_pcs_list) >= r_tune) noise_pcs_list[[r_tune]] else NULL
          nuis_t <- hrf_combine_run_nuisance(nuis_raw_t, pcs_t, n_t)
          z_t_base <- cbind(intercept = 1, trend = seq(-1, 1, length.out = n_t), nuis_t)
          off_t <- if (length(hrf_fit$settings$time_offset) >= r_tune) hrf_fit$settings$time_offset[r_tune] else 0.0
          md_t <- hrf_build_noise_mean_design(n_t, hrf_fit$specification$events[[r_tune]], tr, off_t, k_tune, z_t_base, mask_tune_bool)
          if (md_t$df_nominal >= 16L) {
            nsd_tune <- hrf_estimate_gls_noise_sd(
              y_t[rows_t, active_cols, drop = FALSE],
              md_t$design[rows_t, , drop = FALSE],
              rows_t,
              ar_t
            )
          }
        }
      }

      sc_vox <- scale_decon_input(
        y_t[rows_t, active_cols, drop = FALSE],
        scaling = scaling,
        noise_rho = rho_t,
        noise_psi1 = psi1_t,
        noise_psi2 = psi2_t,
        noise_psi3 = psi3_t,
        noise_psi4 = psi4_t,
        noise_sd_override = nsd_tune,
        rows = rows_t,
        obs_weights = w_rt
      )$y
      if (nuisance_mode %in% c("joint", "preproject")) {
        hp_tune <- if (!is.null(hrf_fit$settings$highpass)) hrf_fit$settings$highpass else 128.0
        nuis_raw_j <- if (!is.null(hrf_fit$specification$nuisance) && length(hrf_fit$specification$nuisance) >= r_tune) {
          hrf_fit$specification$nuisance[[r_tune]]
        } else NULL
        pcs_j <- if (!is.null(noise_pcs_list) && length(noise_pcs_list) >= r_tune) noise_pcs_list[[r_tune]] else NULL
        nuis_j <- hrf_combine_run_nuisance(nuis_raw_j, pcs_j, n_t)
        z_j <- hrf_build_decon_nuisance_matrix(
          n_obs = n_t, tr = hrf_fit$settings$tr, rows = rows_t,
          nuisance_r = nuis_j, nuisance_mode = nuisance_mode, highpass_sec = hp_tune,
          drift_active = drift_active
        )[rows_t, , drop = FALSE]
        if (ncol(z_j) > 0L) {
          z_jw <- if (has_w_rt) sweep(z_j, 1L, sqrt(pmax(w_rt, 1e-4)), "*") else z_j
          nrm_j <- sqrt(colSums(z_jw^2)); nrm_j[nrm_j == 0] <- 1
          sv_j <- hrf_safe_svd(sweep(z_jw, 2L, nrm_j, "/"), nv = 0)
          rk_j <- if (length(sv_j$d) && sv_j$d[1L] > 0) sum(sv_j$d > 1e-8 * sv_j$d[1L]) else 0L
          if (rk_j > 0L && length(rows_t) - rk_j >= 8L) {
            qj <- sv_j$u[, seq_len(rk_j), drop = FALSE]
            if (has_w_rt) {
              sw_rt <- sqrt(pmax(w_rt, 1e-4))
              coef_j <- crossprod(qj, sweep(sc_vox, 1L, sw_rt, "*"))
              sc_vox <- sc_vox - sweep(qj %*% coef_j, 1L, sw_rt, "/")
            } else {
              sc_vox <- sc_vox - qj %*% crossprod(qj, sc_vox)
            }
          }
        }
      }
      if (has_w_rt) {
        sc_vox <- sweep(sc_vox, 1L, sqrt(pmin(pmax(w_rt, 1e-4), 1.0)), "*")
      }
      y_tune <- sc_vox
      groups_all <- hrf_fit$specification$spatial_groups
      if (!identical(mode_type, "tier_background") && !is.null(groups_all) && length(groups_all) == n_vox) {
        g_act <- as.character(groups_all[active_cols])
        u_groups <- unique(g_act[!is.na(g_act)])
        parcel_cols <- list()
        for (g in u_groups) {
          idx_g <- which(g_act == g)
          if (length(idx_g) >= 2L) {
            wg <- sqrt(pmax(r2_r_vec[active_cols[idx_g]], 1e-3))
            wg <- wg / sum(wg)
            sub_g <- sc_vox[, idx_g, drop = FALSE]
            # Align voxel signs deterministically via leading SVD direction so mixed-polarity
            # task voxels inside a parcel do not cancel out during GCV pooling.
            sgn_g <- tryCatch({
              vg <- hrf_safe_svd(sub_g, nu = 0, nv = 1)$v[, 1L]
              if (vg[which.max(wg)] < 0) vg <- -vg
              ifelse(vg >= 0, 1.0, -1.0)
            }, error = function(e) rep(1.0, length(idx_g)))
            p_series <- as.numeric(sub_g %*% (wg * sgn_g))
            if (all(is.finite(p_series)) && stats::sd(p_series) > 1e-8) {
              parcel_cols[[g]] <- p_series
            }
          }
        }
        if (length(parcel_cols)) {
          p_mat <- do.call(cbind, parcel_cols)
          p_scaled <- scale_decon_input(
            p_mat,
            scaling = "noise_sd",
            noise_rho = rho_t,
            noise_psi1 = psi1_t,
            noise_psi2 = psi2_t,
            noise_psi3 = psi3_t,
            noise_psi4 = psi4_t,
            rows = rows_t,
            obs_weights = w_rt
          )$y
          y_tune <- cbind(p_scaled, sc_vox[, seq_len(min(32L, ncol(sc_vox))), drop = FALSE])
        }
      }
      if (tail_padding > 0L) {
        y_tune <- rbind(y_tune, matrix(0.0, as.integer(tail_padding), ncol(y_tune)))
      }

      tune_res <- deconvolve_reglin(
        BOLDobs = y_tune,
        kernel = k_tune,
        lambda = NULL,
        lambda_grid = lambda_grid,
        penalty = penalty,
        tune_by = "global",
        normalize = FALSE,
        demean = FALSE,
        ridge_floor = ridge_floor,
        return_diagnostics = TRUE,
        prewhiten_gcv = rho_t != 0,
        prewhiten_rho = rho_t,
        backend = if (backend == "wiener") "wiener" else "cpp",
        gcv_rule = gcv_rule,
        scaling = "none"
      )
      tun <- tune_res$tuning
      if (backend == "wiener" || !is.data.frame(tun) || nrow(tun) != length(lambda_grid)) {
        lam_out <- as.numeric(tune_res$lambda[1L])
        if (!is.finite(lam_out)) lam_out <- lambda_grid[1L]
        wiener_lambdas <- c(wiener_lambdas, lam_out)
        wiener_tunings[[length(wiener_tunings) + 1L]] <- tun
        run_weights <- c(run_weights, as.numeric(length(rows_t)))
        next
      }

      # Evaluate exact interior gain-restored variance and modulation transfer on a
      # compact window using a fast 2-column normal-equation solve.
      l_k <- length(k_tune)
      n_win <- min(length(rows_t), max(64L, 2L * l_k + 16L))
      if (n_win < l_k + 8L) n_win <- max(length(rows_t), l_k + 8L)
      t0 <- max(3L, min(n_win - 2L, as.integer(floor(n_win * 0.45))))
      t1 <- min(n_win, t0 + 1L)

      h_full <- decon_convolution_matrix(n_win, k_tune)[seq_len(n_win), , drop = FALSE]
      wh <- hrf_apply_whitening(h_full, seq_len(n_win), ar_t, unit_marginal = FALSE)
      hth <- crossprod(wh)
      p_pen <- decon_penalty_crossprod(ncol(wh), penalty)
      idx_lat <- (l_k - 1L) + seq_len(n_win)
      idx0 <- idx_lat[t0]
      idx1 <- idx_lat[t1]
      h0 <- wh[, idx0]
      h1 <- wh[, idx1]
      h_diff <- h0 - h1

      e_rhs <- matrix(0.0, nrow = ncol(wh), ncol = 2L)
      e_rhs[idx0, 1L] <- 1.0
      e_rhs[idx0, 2L] <- 1.0
      e_rhs[idx1, 2L] <- -1.0

      n_grid <- length(lambda_grid)
      v_lvl <- v_diff <- c_mtf <- g_rel <- numeric(n_grid)
      rf <- max(ridge_floor, 1e-8)
      d_isi <- max(2L, as.integer(round(2.5 / tr)))
      far_mask <- abs(seq_len(n_win) - t0) >= d_isi
      for (i in seq_len(n_grid)) {
        lhs <- hth + lambda_grid[i] * p_pen
        diag(lhs) <- diag(lhs) + rf
        q_mat <- tryCatch(solve(lhs, e_rhs), error = function(e) MASS::ginv(lhs) %*% e_rhs)
        b0 <- as.numeric(wh %*% q_mat[, 1L])
        bdiff <- as.numeric(wh %*% q_mat[, 2L])
        r_row <- as.numeric(crossprod(q_mat[, 1L], hth[, idx_lat, drop = FALSE]))

        g_lvl <- max(sum(b0 * h0), 1e-4)
        g_diff <- max(sum(bdiff * h_diff), 1e-4)
        leak_isi <- if (any(far_mask)) {
          min(0.85, sqrt(sum(r_row[far_mask]^2)) / g_lvl)
        } else {
          0.0
        }
        g_eff <- g_lvl * max(0.15, 1.0 - 0.65 * (leak_isi^2))
        v_lvl[i] <- sum(b0^2) / (g_eff^2)
        v_diff[i] <- sum(bdiff^2) / (g_diff^2)
        g_rel[i] <- g_lvl

        conc_0 <- g_lvl / sqrt(max(sum(r_row^2), 1e-12))
        step_rel <- g_diff / (2.0 * g_lvl)
        c_mtf[i] <- 0.5 * conc_0 + 0.5 * step_rel
      }

      g_rel <- pmin(1.0, pmax(0.0, g_rel / max(g_rel[1L], 1e-6)))
      c_rel <- pmin(1.0, pmax(0.0, 0.65 * (c_mtf / max(c_mtf[1L], 1e-6)) + 0.35 * g_rel))
      p_exp <- min(1.0, tr^0.65)
      v_norm <- 0.50 * (v_lvl / min(v_lvl))^p_exp + 0.50 * (v_diff / min(v_diff))^p_exp

      gcv_vec <- as.numeric(tun$gcv)
      gcv_min <- min(gcv_vec, na.rm = TRUE)
      gcv_end <- gcv_vec[n_grid]
      snr_gcv <- if (is.finite(gcv_min) && gcv_min > 0 && is.finite(gcv_end)) {
        max(0.0, (gcv_end - gcv_min) / gcv_min)
      } else {
        0.0
      }

      if (identical(mode_type, "tier_background")) {
        snr_task <- min(0.20, 0.15 * snr_gcv)
        w_bias <- 0.04 + 0.20 * snr_task
        v_use <- 0.75 * (v_lvl / min(v_lvl))^p_exp + 0.25 * (gcv_vec / max(gcv_min, 1e-8))
      } else if (identical(mode_type, "tier_active")) {
        sc_ref <- if (!is.null(r2_tier_score)) r2_tier_score[active_cols] else r2_ok
        r2_rep <- if (has_r2 && length(sc_ref)) {
          as.numeric(stats::quantile(sc_ref, probs = 0.75, names = FALSE, na.rm = TRUE))
        } else {
          0.0
        }
        if (!is.finite(r2_rep)) r2_rep <- 0.0
        snr_r2 <- (8.0 / max(tr, 0.25)) * r2_rep / max(1.0 - r2_rep, 0.10)
        snr_task <- min(6.0, max(snr_r2, if (has_r2) 0.35 * snr_gcv else 2.0 * snr_gcv))
        w_bias <- 0.25 + 0.35 * snr_task
        v_use <- v_norm
      } else {
        r2_top <- if (has_r2 && length(r2_ok)) {
          mean(utils::head(sort(r2_ok, decreasing = TRUE), min(12L, length(r2_ok))), na.rm = TRUE)
        } else {
          0.0
        }
        if (!is.finite(r2_top)) r2_top <- 0.0
        snr_r2 <- (8.0 / max(tr, 0.25)) * r2_top / max(1.0 - r2_top, 0.10)
        snr_task <- min(6.0, max(snr_r2, if (has_r2) 0.35 * snr_gcv else 2.0 * snr_gcv))
        w_bias <- 0.25 + 0.35 * snr_task
        v_use <- v_norm
      }

      idx_gmin <- which.min(gcv_vec)
      gcv_excess <- rep(0.0, n_grid)
      if (length(idx_gmin) == 1L && is.finite(gcv_min) && gcv_min > 0) {
        s_hi <- seq_len(n_grid) >= idx_gmin
        gcv_excess[s_hi] <- pmax(0.0, (gcv_vec[s_hi] - gcv_min) / gcv_min)
      }

      gcv_adaptive <- v_use + w_bias * (1.0 - c_rel)^2 + 0.30 * (gcv_excess^2)

      tun$var_level <- v_lvl
      tun$var_contrast <- v_diff
      tun$mtf_rel <- c_rel
      tun$gcv_adaptive <- gcv_adaptive
      tun$snr_task <- rep(snr_task, n_grid)
      run_tuning_list[[length(run_tuning_list) + 1L]] <- tun
      run_weights <- c(run_weights, as.numeric(length(rows_t)))
    }

    if (length(wiener_lambdas) > 0L && !length(run_tuning_list)) {
      w_norm <- run_weights / sum(run_weights)
      ord_w <- order(w_norm, wiener_lambdas)
      lam_pool <- exp(sum(w_norm[ord_w] * log(pmax(wiener_lambdas[ord_w], 1e-12))))
      best_w <- which.min(abs(log(pmax(lambda_grid, 1e-12)) - log(max(lam_pool, 1e-12))))
      return(list(lambda = lambda_grid[best_w], tuning = wiener_tunings[[ord_w[1L]]]))
    }

    if (length(run_tuning_list) > 0L) {
      n_grid <- length(lambda_grid)
      w_norm <- run_weights / sum(run_weights)
      # Sort run contributions canonically so floating-point summation is strictly
      # permutation-invariant across arbitrary run orderings.
      sig_sort <- vapply(run_tuning_list, function(df_r) as.numeric(df_r$gcv_adaptive[1L]), numeric(1L))
      ord_r <- order(w_norm, sig_sort)
      tun_pool <- run_tuning_list[[ord_r[1L]]]
      avg_cols <- intersect(
        c("df", "rss", "gcv", "gcv_se", "var_level", "var_contrast", "mtf_rel", "gcv_adaptive", "snr_task"),
        names(tun_pool)
      )
      for (cn in avg_cols) {
        acc_c <- numeric(n_grid)
        for (k in ord_r) {
          acc_c <- acc_c + w_norm[k] * as.numeric(run_tuning_list[[k]][[cn]])
        }
        tun_pool[[cn]] <- acc_c
      }
      gcv_adaptive_pool <- tun_pool$gcv_adaptive
      best_idx <- which.min(gcv_adaptive_pool)
      if (!length(best_idx) || !is.finite(best_idx)) best_idx <- 1L
      if (gcv_rule == "1se" && best_idx < n_grid) {
        var_min_idx <- max(which.min(tun_pool$var_level), which.min(tun_pool$var_contrast))
        if (best_idx < var_min_idx) {
          best_idx <- min(var_min_idx, best_idx + 1L)
        }
      }
      return(list(lambda = lambda_grid[best_idx], tuning = tun_pool))
    }

    list(
      lambda = lambda_grid[1L],
      tuning = data.frame(
        kernel = 1L, lambda = lambda_grid, df = NA_real_, rss = NA_real_,
        gcv = NA_real_, gcv_se = NA_real_, var_level = NA_real_,
        var_contrast = NA_real_, mtf_rel = NA_real_, gcv_adaptive = NA_real_, snr_task = NA_real_
      )
    )
  }

  tuned_global <- tune_subset_cols(mode_type = "global_active")
  lam_ref <- tuned_global$lambda
  tun_ref <- tuned_global$tuning
  snr_ref <- if (is.data.frame(tun_ref) && "snr_task" %in% names(tun_ref) && nrow(tun_ref) >= 1L) {
    as.numeric(tun_ref$snr_task[1L])
  } else {
    NA_real_
  }

  global_fallback <- function() {
    list(
      lambda = lam_ref,
      tuning = tun_ref,
      lambda_by_voxel = rep(lam_ref, n_vox),
      lambda_strata = data.frame(
        tier = 1L, label = "global", n_voxels = n_vox,
        mean_r2 = mean(r2_v, na.rm = TRUE), snr_task = snr_ref, lambda = lam_ref,
        stringsAsFactors = FALSE
      ),
      lambda_mode_used = "global"
    )
  }

  if (identical(lambda_mode, "global") || max_lambda_strata <= 1L ||
      identical(backend, "wiener") || !has_r2 || n_vox < 12L) {
    return(global_fallback())
  }

  # Compute out-of-sample spatially stabilized task informativeness r2_stab
  r2_oos <- rowMeans(0.65 * r2_sel_mat + 0.35 * r2_mat, na.rm = TRUE)
  r2_oos[!is.finite(r2_oos)] <- 0.0
  tp_med <- apply(tp_mat, 1L, stats::median, na.rm = TRUE)
  tp_med[!is.finite(tp_med)] <- 1.0

  sg_lam <- NULL
  if (inherits(hrf_fit$specification$spatial_graph, "medusa_spatial_graph") &&
      hrf_fit$specification$spatial_graph$n_voxels == n_vox) {
    sg_lam <- hrf_fit$specification$spatial_graph
  } else if (!is.null(voxel_map) && is.data.frame(voxel_map) && nrow(voxel_map) == n_vox) {
    sg_lam <- tryCatch(
      build_spatial_graph(
        list(voxel_map = voxel_map, geometry = geometry),
        spatial_groups = hrf_fit$specification$spatial_groups,
        parcel_boundary = "soft",
        voxel_ids = hrf_fit$voxel_ids
      ),
      error = function(e) NULL
    )
  } else if (!is.null(hrf_fit$specification$spatial_coords)) {
    sg_lam <- tryCatch(
      build_spatial_graph(
        hrf_fit$specification$spatial_coords,
        spatial_fwhm = hrf_fit$settings$spatial_fwhm,
        spatial_groups = hrf_fit$specification$spatial_groups,
        parcel_boundary = "soft",
        voxel_ids = hrf_fit$voxel_ids
      ),
      error = function(e) NULL
    )
  }

  if (!is.null(sg_lam)) {
    pos_oos <- r2_oos[r2_oos > 0]
    rng_bw <- if (length(pos_oos) >= 4L) {
      max(0.025, 1.5 * stats::mad(pos_oos, constant = 1.4826))
    } else {
      0.035
    }
    g_ids_lam <- hrf_encode_group_ids(hrf_fit$specification$spatial_groups, n_vox)
    cw_lam <- if (!is.null(hrf_fit$settings$parcel_cross_weight)) {
      as.numeric(hrf_fit$settings$parcel_cross_weight)
    } else {
      0.35
    }
    r2_stab <- as.numeric(smooth_scalar_csr_cpp(
      sg_lam$row_ptr, sg_lam$col_idx, sg_lam$dist, sg_lam$sigma,
      r2_oos, rng_bw, g_ids_lam, cw_lam
    ))
    r2_stab[!is.finite(r2_stab)] <- r2_oos[!is.finite(r2_stab)]
  } else if (!is.null(hrf_fit$specification$spatial_groups) &&
             length(hrf_fit$specification$spatial_groups) == n_vox) {
    g_chr <- as.character(hrf_fit$specification$spatial_groups)
    g_means <- tapply(r2_oos, g_chr, mean, na.rm = TRUE)
    r2_stab <- 0.50 * r2_oos + 0.50 * as.numeric(g_means[g_chr])
  } else {
    r2_stab <- r2_oos
  }

  # Cross-run Fisher combined task p-value and BH q-value to distinguish genuine focal
  # task subpopulations from random multiple-comparisons fluctuations on pure-null maps
  p_comb_v <- stats::pchisq(
    -2.0 * rowSums(log(pmin(pmax(tp_mat, 1e-16), 1.0))),
    df = 2.0 * ncol(tp_mat),
    lower.tail = FALSE
  )
  p_comb_v[!is.finite(p_comb_v)] <- 1.0
  q_comb_v <- stats::p.adjust(p_comb_v, method = "BH")

  bg_mask <- (q_comb_v > 0.20) & (r2_stab < 0.025)
  bg_idx <- which(bg_mask)
  act_idx <- which(!bg_mask)
  n_sig_fdr <- sum((q_comb_v < 0.10) & (r2_stab >= 0.015))
  r2_span <- diff(range(r2_stab, na.rm = TRUE))

  if (identical(lambda_mode, "auto")) {
    # Only stratify when the map contains genuine FDR-supported task activity alongside
    # either silent background cortex or broad SNR heterogeneity across >= 24 voxels.
    if (n_sig_fdr < 4L || length(act_idx) < 4L ||
        (length(bg_idx) < 4L && (r2_span < 0.030 || n_vox < 24L))) {
      return(global_fallback())
    }
  }

  tier_cols <- list()
  tier_is_bg <- logical(0L)
  tier_labels <- character(0L)

  if (length(act_idx) < 4L) {
    # Forced "stratified" on a nearly uniform/low-SNR map: split at median if informative
    if (n_vox >= 8L && r2_span > 1e-4) {
      ord_all <- order(r2_stab, seq_len(n_vox))
      half_n <- as.integer(floor(n_vox / 2L))
      tier_cols <- list(ord_all[seq_len(half_n)], ord_all[(half_n + 1L):n_vox])
      tier_is_bg <- c(TRUE, FALSE)
      tier_labels <- c("background", "low_snr")
    } else {
      return(global_fallback())
    }
  } else {
    has_bg_tier <- length(bg_idx) >= 4L
    if (has_bg_tier) {
      tier_cols[[1L]] <- bg_idx
      tier_is_bg <- c(tier_is_bg, TRUE)
      tier_labels <- c(tier_labels, "background")
      k_act_budget <- max(1L, max_lambda_strata - 1L)
    } else {
      act_idx <- seq_len(n_vox)
      k_act_budget <- max_lambda_strata
    }

    min_act_sz <- max(8L, min(16L, as.integer(floor(length(act_idx) / max(1L, k_act_budget)))))
    k_act <- min(k_act_budget, max(1L, as.integer(floor(length(act_idx) / min_act_sz))))
    ord_a <- order(r2_stab[act_idx], act_idx)
    act_sorted <- act_idx[ord_a]
    bin_id <- as.integer(ceiling(seq_along(act_sorted) * k_act / length(act_sorted)))
    act_bins <- split(act_sorted, bin_id)

    # Merge adjacent active bins if their mean stabilized R^2 differs by < 0.012
    merged_bins <- list(act_bins[[1L]])
    if (length(act_bins) >= 2L) {
      for (bi in 2:length(act_bins)) {
        prev_b <- merged_bins[[length(merged_bins)]]
        cur_b <- act_bins[[bi]]
        if (abs(mean(r2_stab[cur_b]) - mean(r2_stab[prev_b])) < 0.012) {
          merged_bins[[length(merged_bins)]] <- c(prev_b, cur_b)
        } else {
          merged_bins[[length(merged_bins) + 1L]] <- cur_b
        }
      }
    }

    n_ab <- length(merged_bins)
    act_names <- if (n_ab == 1L) {
      "high_snr"
    } else if (n_ab == 2L) {
      c("low_snr", "high_snr")
    } else if (n_ab == 3L) {
      c("low_snr", "mid_snr", "high_snr")
    } else {
      paste0("snr_tier_", seq_len(n_ab))
    }
    for (bi in seq_len(n_ab)) {
      tier_cols[[length(tier_cols) + 1L]] <- merged_bins[[bi]]
      tier_is_bg <- c(tier_is_bg, FALSE)
      tier_labels <- c(tier_labels, act_names[bi])
    }
  }

  G <- length(tier_cols)
  if (G <= 1L) {
    return(global_fallback())
  }

  # Order tiers strictly by ascending mean r2_stab (background first, highest SNR last)
  tier_mean_r2 <- vapply(tier_cols, function(idx) mean(r2_stab[idx], na.rm = TRUE), numeric(1L))
  ord_t <- order(tier_is_bg, decreasing = TRUE)
  if (sum(!tier_is_bg) >= 2L) {
    act_pos <- which(!tier_is_bg)
    ord_act_t <- act_pos[order(tier_mean_r2[act_pos])]
    ord_t <- c(which(tier_is_bg), ord_act_t)
  }
  tier_cols <- tier_cols[ord_t]
  tier_is_bg <- tier_is_bg[ord_t]
  tier_labels <- tier_labels[ord_t]
  tier_mean_r2 <- tier_mean_r2[ord_t]
  tier_n <- vapply(tier_cols, length, integer(1L))

  tier_lam_raw <- numeric(G)
  tier_snr <- numeric(G)
  for (g in seq_len(G)) {
    mt_g <- if (tier_is_bg[g]) "tier_background" else "tier_active"
    res_g <- tune_subset_cols(cand_cols = tier_cols[[g]], mode_type = mt_g, r2_tier_score = r2_stab)
    tier_lam_raw[g] <- res_g$lambda
    tier_snr[g] <- if (is.data.frame(res_g$tuning) && "snr_task" %in% names(res_g$tuning) && nrow(res_g$tuning) >= 1L) {
      as.numeric(res_g$tuning$snr_task[1L])
    } else {
      NA_real_
    }
  }

  # Anchor background >= lam_ref >= highest-SNR focal tier, then enforce monotone
  # non-increasing log(lambda) across tiers via weighted PAVA and snap to lambda_grid.
  if (tier_is_bg[1L]) {
    tier_lam_raw[1L] <- max(tier_lam_raw[1L], lam_ref)
  }
  if (!tier_is_bg[G]) {
    tier_lam_raw[G] <- min(tier_lam_raw[G], lam_ref)
  }
  log_grid <- log(pmax(lambda_grid, 1e-12))
  log_iso <- hrf_pava_nonincreasing(log(pmax(tier_lam_raw, 1e-12)), w = sqrt(tier_n))
  tier_lam <- vapply(log_iso, function(lg) {
    lambda_grid[which.min(abs(log_grid - lg))]
  }, numeric(1L))
  # Guarantee exact discrete monotonicity after grid snapping
  if (G >= 2L) {
    for (g in 2:G) {
      if (tier_lam[g] > tier_lam[g - 1L]) {
        tier_lam[g] <- tier_lam[g - 1L]
      }
    }
  }

  lambda_by_voxel <- rep(lam_ref, n_vox)
  for (g in seq_len(G)) {
    lambda_by_voxel[tier_cols[[g]]] <- tier_lam[g]
  }

  lambda_strata <- data.frame(
    tier = seq_len(G),
    label = tier_labels,
    n_voxels = tier_n,
    mean_r2 = tier_mean_r2,
    snr_task = tier_snr,
    lambda = tier_lam,
    stringsAsFactors = FALSE
  )

  list(
    lambda = lam_ref,
    tuning = tun_ref,
    lambda_by_voxel = lambda_by_voxel,
    lambda_strata = lambda_strata,
    lambda_mode_used = if (length(unique(tier_lam)) > 1L) "stratified" else "global"
  )
}
