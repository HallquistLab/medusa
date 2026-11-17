#' Diagnose amplitude attenuation and event leakage for a frozen inverse
#'
#' @param deconvolution Result of deconvolve_crossrun_hrf().
#' @param events Optional run-named event tables (or a combined multi-run
#'   \code{data.frame} with a \code{run_col} column) with onset and optional
#'   condition and duration. Defaults to the HRF training event records. Every
#'   record is a separate unit-amplitude source, regardless of its condition label.
#' @param lags Readout times relative to each event onset, in seconds. Default 0.
#' @param weights Linear readout weights, default one per lag. For example,
#'   lags = -1:3 produces a five-sample window sum at TR = 1.
#' @param source_kernel_grid Optional normalized HRF candidates used to generate
#'   synthetic BOLD. Default uses each stratum's fitted HRF. Honors
#'   \code{attr(source_kernel_grid, "normalize")} when present, otherwise matching
#'   the fitted normalization mode. Supplying alternatives evaluates sensitivity
#'   to HRF mismatch; it does not estimate the true HRF.
#' @param return_matrices Return event-by-event recovered, ideal, and excess
#'   resolution matrices, default FALSE.
#' @param max_events Explicit limit per run on the quadratic diagnostic matrix,
#'   default 2000. Use a smaller prespecified subset or explicitly raise the limit.
#' @param run_col Column name in \code{events} identifying runs when \code{events}
#'   is supplied as a combined multi-run \code{data.frame}. Default \code{"run"}.
#' @return A list with per-event diagnostics, stratum summaries, voxel-to-stratum
#'   mappings, optional resolution matrices, and the conditional assumptions.
#' @details For source events B, convolution H, frozen inverse D, and readout S,
#'   recovered = S D H B and ideal = S B. Gain is their diagonal ratio. Excess
#'   leakage is recovered minus ideal, normalized by each receiving event's ideal
#'   diagonal. Subtracting ideal retains natural overlap of neural readout windows.
#'   Previous/next refer to chronological source records, including other event
#'   types. Zero ideal diagonals and readouts outside the acquisition grid are NA.
#'
#'   The calculation reuses one operator per fitted HRF/noise stratum and uses
#'   BLAS across source events, not a simulation or fit per voxel. It includes the
#'   actual masks, nuisance projection, padding, centering, and inverse penalty.
#'   Scaling divisors cancel in the gain ratios and are omitted. Their estimated
#'   values, noise parameters, selected HRFs, and lambda are held fixed: these are
#'   conditional linear diagnostics, not uncertainty intervals or a claim about
#'   physiological neural units. Source events of positive duration have unit
#'   height (fractional acquisition-bin occupancy), not unit integrated area.
#'   No automatic correction is applied; a scalar cannot undo a full matrix of
#'   attenuation and cross-talk. Wiener diagnostics are not supported here.
#' @export
diagnose_crossrun_recovery <- function(deconvolution, events = NULL, lags = 0,
    weights = NULL, source_kernel_grid = NULL, return_matrices = FALSE, max_events = 2000L,
    run_col = "run") {
  fit <- deconvolution$hrf_fit; settings <- deconvolution$settings
  if (!inherits(fit, "crossrun_hrf") || is.null(deconvolution$scaling) ||
      is.null(settings$nuisance_mode)) stop("Use a current deconvolve_crossrun_hrf result")
  if (settings$backend == "wiener") stop("Recovery diagnostics require the cpp or r GLS inverse")
  checkmate::assert_numeric(lags, min.len = 1L, any.missing = FALSE, finite = TRUE)
  if (is.null(weights)) weights <- rep(1, length(lags))
  checkmate::assert_numeric(weights, len = length(lags), any.missing = FALSE, finite = TRUE)
  checkmate::assert_flag(return_matrices)
  checkmate::assert_integerish(max_events, lower = 1L, len = 1L)
  if (is.null(events)) events <- fit$specification$events
  if (is.data.frame(events)) {
    events <- hrf_resolve_events_input(
      events,
      stats::setNames(vector("list", length(fit$run_ids)), fit$run_ids),
      run_col = run_col,
      allow_empty_runs = FALSE
    )
  }
  if (!is.list(events) || !identical(names(events), fit$run_ids))
    stop("events must be a list in fitted run order")
  sources <- if (is.null(source_kernel_grid)) {
    NULL
  } else {
    raw_src <- if (is.list(source_kernel_grid)) {
      lapply(source_kernel_grid, as.numeric)
    } else {
      src_mat <- as.matrix(source_kernel_grid)
      lapply(seq_len(ncol(src_mat)), function(i) as.numeric(src_mat[, i]))
    }
    source_norm_mode <- if (!is.null(attr(source_kernel_grid, "normalize"))) {
      attr(source_kernel_grid, "normalize")
    } else if (infer_kernel_grid_normalize(raw_src) == "positive_lobe") {
      "positive_lobe"
    } else if (!is.null(fit$settings$kernel_normalize)) {
      fit$settings$kernel_normalize
    } else {
      "signed_sum"
    }
    do.call(cbind, build_decon_kernel_grid(c(0, 1), source_kernel_grid, normalize = source_norm_mode))
  }
  if (!is.null(sources) && !is.null(attr(source_kernel_grid, "tr")) &&
      !isTRUE(all.equal(attr(source_kernel_grid, "tr"), fit$settings$tr)))
    stop("source_kernel_grid sampling TR does not match the fit")
  all_events <- summaries <- mappings <- matrices <- list()
  finite_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  for (r in seq_along(fit$run_ids)) {
    n <- fit$n_time[r]; ev <- events[[r]]; tr <- fit$settings$tr; offset <- fit$settings$time_offset[r]
    checkmate::assert_data_frame(ev, min.rows = 1L)
    checkmate::assert_numeric(ev$onset, any.missing = FALSE, finite = TRUE)
    if (nrow(ev) > max_events) stop("Diagnostic event matrix exceeds max_events")
    duration <- if (is.null(ev$duration)) rep(0, nrow(ev)) else ev$duration
    checkmate::assert_numeric(duration, lower = 0, len = nrow(ev), any.missing = FALSE, finite = TRUE)
    m <- nrow(ev); order <- order(ev$onset, seq_len(m)); position <- match(seq_len(m), order)
    # Identify chronologically preceding and succeeding events (across all condition types)
    # so we can measure immediate neighbor cross-talk in the resolution matrix.
    previous <- ifelse(position > 1L, order[pmax(1L, position - 1L)], NA_integer_)
    next_event <- ifelse(position < m, order[pmin(m, position + 1L)], NA_integer_)
    # Build linear readout matrix S (m_events x n_time) using two-bin fractional
    # bracket interpolation at each event onset + lag.
    readout <- matrix(0, m, n); observed <- rep(TRUE, m)
    for (j in seq_along(lags)) {
      at <- (ev$onset + lags[j] - offset) / tr
      lo <- floor(at); fraction <- at - lo
      ok <- at >= 0 & at <= n - 1L
      if (weights[j] == 0) next
      observed <- observed & ok
      index <- which(ok)
      readout[cbind(index, lo[index] + 1L)] <- readout[cbind(index, lo[index] + 1L)] + weights[j] * (1 - fraction[index])
      index <- which(ok & fraction > 0 & lo + 2L <= n)
      readout[cbind(index, lo[index] + 2L)] <- readout[cbind(index, lo[index] + 2L)] + weights[j] * fraction[index]
    }
    readout[!observed, ] <- NA_real_
    # Construct the unblurred neural design B (using a delta kernel [1, 0]) and
    # compute ideal = S * B so that natural overlap of neural readout windows is
    # separated from HRF + regularization blur (excess = S * D * H * B - S * B).
    neural <- hrf_event_design_cpp(n, tr, offset, ev$onset, duration,
      seq_len(m), m, matrix(c(1, 0), ncol = 1L))[[1L]]
    ideal <- readout %*% neural
    ideal_diagonal <- diag(ideal)
    mask <- if (settings$observation_mask == "fit") fit$specification$valid_rows[[r]] else rep(TRUE, n)
    rows <- which(!is.na(mask) & mask)
    if (length(rows) < 3L) next
    fw_r <- resolve_run_obs_weights(deconvolution, r, length(rows))
    hp_diag <- if (!is.null(settings$highpass)) settings$highpass else if (!is.null(fit$settings$highpass)) fit$settings$highpass else 128.0
    z <- hrf_build_decon_nuisance_matrix(
      n_obs = n, tr = tr, rows = rows,
      nuisance_r = fit$specification$nuisance[[r]],
      nuisance_mode = settings$nuisance_mode, highpass_sec = hp_diag,
      drift_active = !isFALSE(settings$drift_preproject_active)
    )
    strat <- hrf_stratum_kernels(fit, r)
    selected <- fit$selections[[r]]$kernel_index
    rho <- strat$rho
    valid <- deconvolution$scaling[[r]]$status == "ok"
    lam_vec_r <- if (!is.null(deconvolution$scaling[[r]]$lambda)) {
      as.numeric(deconvolution$scaling[[r]]$lambda)
    } else {
      rep(as.numeric(settings$lambda)[1L], length(fit$voxel_ids))
    }
    lam_vec_r[!is.finite(lam_vec_r) | lam_vec_r < 0] <- as.numeric(settings$lambda)[1L]
    # Group voxels by their (selected_hrf_kernel / weight_vector, ar1_rho, lambda) stratum so the
    # inverse operator D is constructed and multiplied only once per unique stratum.
    keys <- if (length(unique(lam_vec_r[valid])) > 1L) {
      paste(strat$keys, sprintf("%.12g", lam_vec_r), sep = "|lam:")
    } else {
      strat$keys
    }
    mappings[[r]] <- data.frame(run = fit$run_ids[r], voxel = fit$voxel_ids,
      stratum = keys, lambda = lam_vec_r, usable = valid)
    for (key in unique(keys[valid])) {
      voxels <- which(keys == key & valid); k <- selected[voxels[1L]]; ar <- rho[voxels[1L]]
      ar_op <- if (!is.null(strat$get_ar)) strat$get_ar(voxels[1L]) else ar
      hrf_k <- strat$get_kernel(voxels[1L])
      lam_k <- lam_vec_r[voxels[1L]]
      operator <- hrf_inverse_operator(n, as.numeric(hrf_k), rows, ar_op, z, settings, obs_weights = fw_r, lambda = lam_k)
      readout_operator <- readout %*% operator
      source <- if (is.null(sources)) hrf_k else sources
      designs <- hrf_event_design_cpp(n, tr, offset, ev$onset, duration, seq_len(m), m, source)
      for (s in seq_len(ncol(source))) {
        y <- designs[[s]][rows, , drop = FALSE]
        if (settings$scaling != "none") y <- center_cols_weighted(y, fw_r)
        recovered <- readout_operator %*% y
        excess <- recovered - ideal
        divisor <- ideal_diagonal
        divisor[!is.finite(divisor) | abs(divisor) < 1e-10] <- NA_real_
        normalized <- sweep(excess, 1L, divisor, "/")
        gain <- diag(recovered) / divisor
        off_diagonal <- normalized; diag(off_diagonal) <- NA_real_
        meta <- data.frame(run = fit$run_ids[r], stratum = key, fitted_kernel = colnames(fit$kernels)[k],
          source_kernel = if (is.null(sources)) "fitted" else paste0("source_", s), noise_rho = ar,
          lambda = lam_k, n_voxels = length(voxels))
        all_events[[length(all_events) + 1L]] <- cbind(meta[rep(1L, m), ],
          data.frame(event = seq_len(m), onset = ev$onset,
            condition = if (is.null(ev$condition)) NA_character_ else as.character(ev$condition),
            ideal_diagonal = ideal_diagonal, recovered_diagonal = diag(recovered), amplitude_gain = gain,
            previous_event = previous, previous_excess_leakage = normalized[cbind(seq_len(m), previous)],
            next_event = next_event, next_excess_leakage = normalized[cbind(seq_len(m), next_event)],
            mean_abs_excess_leakage = apply(abs(off_diagonal), 1L, finite_mean)))
        summaries[[length(summaries) + 1L]] <- cbind(meta,
          data.frame(n_events = m, n_valid_gains = sum(is.finite(gain)), mean_gain = finite_mean(gain),
            mean_abs_gain_error = finite_mean(abs(gain - 1)),
            mean_abs_excess_leakage = finite_mean(abs(off_diagonal))))
        if (return_matrices) matrices[[paste(fit$run_ids[r], key, s, sep = "/")]] <-
          list(recovered = recovered, ideal = ideal, excess = excess, normalized_excess = normalized)
      }
    }
  }
  list(events = if (length(all_events)) do.call(rbind, all_events) else data.frame(),
    summary = if (length(summaries)) do.call(rbind, summaries) else data.frame(),
    voxels = if (length(mappings)) do.call(rbind, mappings) else data.frame(),
    matrices = if (return_matrices) matrices else NULL,
    settings = list(lags = lags, weights = weights, inverse = settings,
      conditional_on = "HRF, covariance, lambda, masks, nuisance design, and scaling divisors",
      reference = paste0("unit-height event sources; ",
        if (!is.null(fit$settings$kernel_normalize)) fit$settings$kernel_normalize else "signed_sum",
        "-normalized HRFs; ratios within analysis units"),
      automatic_correction = FALSE))
}

# Build a trial-by-time epoch readout matrix S(tau) of dimension (m_events x n_time).
# For each event onset t_i and relative lag tau, S(tau) maps the continuous
# deconvolved time series (length n_time, sampled every TR seconds starting at
# time_offset) onto a single trial-level readout value at t_i + tau.
# - "linear" performs two-point linear bracket interpolation between adjacent TRs.
# - "boxcar", "triangle", and "gaussian" evaluate a 9-point symmetric quadrature
#   over [t_i + tau - 0.5*width, t_i + tau + 0.5*width], attenuating high-frequency
#   deconvolution noise while preserving local event timing.
# - When readout_width is a vector of length(onsets) (used by source_model = "duration"),
#   each trial i receives its own width_scale[i] = max(w_src, duration[i]).
build_epoch_readout_matrix <- function(onsets, lag, n_time, tr, offset,
    readout_kernel = "linear", readout_width = NULL) {
  m <- length(onsets)
  if (is.null(readout_width) || length(readout_width) <= 1L) {
    quad <- epoch_kernel_quadrature(readout_kernel, readout_width, tr)
    width_scale <- rep(1, m)
  } else {
    if (length(readout_width) != m) {
      stop("readout_width must be NULL, a scalar, or a vector of length(onsets)")
    }
    checkmate::assert_numeric(readout_width, lower = .Machine$double.eps, finite = TRUE, any.missing = FALSE)
    quad <- epoch_kernel_quadrature(readout_kernel, 1.0, tr)
    width_scale <- as.numeric(readout_width)
  }
  readout <- matrix(0, m, n_time)
  weight_sum <- rep(0, m)
  for (q in seq_along(quad$offsets)) {
    # Fractional grid index for quadrature node q of each trial
    at <- (onsets + lag + quad$offsets[q] * width_scale - offset) / tr
    lo <- floor(at)
    fraction <- at - lo
    ok <- is.finite(at) & at >= 0 & at <= n_time - 1L
    w_q <- quad$weights[q]
    if (!any(ok) || w_q == 0) next
    idx <- which(ok)
    readout[cbind(idx, lo[idx] + 1L)] <- readout[cbind(idx, lo[idx] + 1L)] +
      w_q * (1 - fraction[idx])
    idx2 <- which(ok & fraction > 0 & lo + 2L <= n_time)
    if (length(idx2)) {
      readout[cbind(idx2, lo[idx2] + 2L)] <- readout[cbind(idx2, lo[idx2] + 2L)] +
        w_q * fraction[idx2]
    }
    weight_sum[idx] <- weight_sum[idx] + w_q
  }
  # Normalize row weights for trials partially supported near run boundaries;
  # trials completely outside the run grid receive NA.
  valid <- weight_sum > 0
  if (any(valid)) {
    readout[valid, ] <- sweep(readout[valid, , drop = FALSE], 1L, weight_sum[valid], "/")
  }
  readout[!valid, ] <- NA_real_
  readout
}

# Construct the target-trial cross-talk de-mixing and gain-correction operator W(tau)
# for a given lag tau.
# Given readout matrix S_tau (m_ok x n_time), inverse deconvolution operator D (n_time x n_obs),
# uncentered forward trial design G_tau_unc (n_obs x m_ok), and column-centered
# design G_tau_cent = G_tau_unc - 1 * mu_G^T:
# 1. K_unc = S_tau * D * G_tau_unc is the uncentered trial-by-trial resolution matrix,
#    whose diagonal g_unc[i] is trial i's isolated amplitude gain and whose off-diagonals
#    measure adjacent-trial overlap leakage.
# 2. Row-normalizing R_unc = diag(g_safe)^{-1} K_unc places all trials on unit diagonal
#    even when trial durations vary (where g_unc[i] scales with duration d_i).
# 3. Penalizing only centered trial fluctuations via Q = I - (1/m)*1*1^T leaves a
#    constant baseline firing rate alpha_0 * 1 unattenuated (since Q * 1 = 0),
#    preventing duration-dependent shrinkage of the task baseline from confounding
#    covariates correlated with reaction time (RT).
# 4. When input BOLD was centered (scaling != "none"), G_tau_cent differs from
#    G_tau_unc by the rank-1 matrix -1 * mu_G^T. A Sherman-Morrison rank-1 update
#    adjusts W_unc into W_full so that constant baseline activity still maps to a
#    constant across trials of unequal duration.
solve_duration_demix_operator <- function(S_tau, operator, G_tau_unc, G_tau_cent,
    scaling_mode, demix_ridge, gain_floor, DG_tau_unc = NULL, DG_tau_cent = NULL) {
  m_ok <- nrow(S_tau)
  n_obs <- ncol(operator)
  K_norm <- S_tau %*% (if (!is.null(DG_tau_cent)) DG_tau_cent else (operator %*% G_tau_cent))
  K_unc <- S_tau %*% (if (!is.null(DG_tau_unc)) DG_tau_unc else (operator %*% G_tau_unc))
  g_unc <- diag(K_unc)
  g_safe <- sign(ifelse(g_unc == 0, 1, g_unc)) * pmax(abs(g_unc), gain_floor)
  R_unc <- sweep(K_unc, 1L, g_safe, "/")
  Q <- diag(m_ok) - matrix(1 / m_ok, m_ok, m_ok)
  lhs <- crossprod(R_unc) + demix_ridge * Q + 1e-8 * diag(m_ok)
  W_unc <- solve(lhs, sweep(t(R_unc), 2L, g_safe, "/"))
  if (scaling_mode != "none") {
    # Sherman-Morrison rank-1 adjustment for column-demeaned forward design G_tau_cent
    c0 <- as.numeric(S_tau %*% (operator %*% rep(1, n_obs)))
    mu_G <- if (nrow(G_tau_unc) >= 1L) {
      as.numeric(G_tau_unc[1L, ] - G_tau_cent[1L, ])
    } else {
      colMeans(G_tau_unc)
    }
    u0 <- as.numeric(W_unc %*% c0)
    denom <- 1 - sum(mu_G * u0)
    if (!is.finite(denom) || abs(denom) < 0.05) {
      denom <- sign(ifelse(!is.finite(denom) || denom == 0, 1, denom)) * 0.05
    }
    W_full <- W_unc + outer(u0, as.numeric(crossprod(mu_G, W_unc))) / denom
  } else {
    W_full <- W_unc
  }
  # Effective post-demixing gain in the centered trial subspace Q
  eff_gain <- if (m_ok > 1L) {
    sum(diag(Q %*% W_full %*% K_norm %*% Q)) / (m_ok - 1L)
  } else {
    1 / (1 + demix_ridge)
  }
  if (!is.finite(eff_gain) || abs(eff_gain) < gain_floor) {
    eff_gain <- 1 / (1 + demix_ridge)
  }
  eff_gain_vec <- if (m_ok > 1L) {
    g_d <- diag(Q %*% W_full %*% K_norm %*% Q) / diag(Q)
    g_d[!is.finite(g_d)] <- eff_gain
    sign(ifelse(g_d == 0, 1, g_d)) * pmax(abs(g_d), gain_floor)
  } else {
    rep(eff_gain, m_ok)
  }
  list(W = W_full, K_norm = K_norm, eff_gain = eff_gain, eff_gain_vec = eff_gain_vec)
}

resolve_run_obs_weights <- function(deconvolution, r, m_obs) {
  if (is.null(deconvolution) || is.null(deconvolution$frame_weights) ||
      length(deconvolution$frame_weights) < r) {
    return(NULL)
  }
  fw <- deconvolution$frame_weights[[r]]
  if (is.null(fw) || length(fw) != m_obs || !all(is.finite(fw))) {
    return(NULL)
  }
  ws_attr <- attr(fw, "weight_space")
  fw <- as.numeric(fw)
  if (!any(abs(fw - 1.0) > 1e-8)) return(NULL)
  if (!is.null(ws_attr)) attr(fw, "weight_space") <- ws_attr
  fw
}

center_cols_weighted <- function(mat, obs_weights = NULL) {
  if (nrow(mat) == 0L || ncol(mat) == 0L) return(mat)
  if (!is.null(obs_weights) && length(obs_weights) == nrow(mat) &&
      all(is.finite(obs_weights)) && any(abs(obs_weights - 1.0) > 1e-8)) {
    w_pos <- pmax(as.numeric(obs_weights), 1e-8)
    w_c <- w_pos / sum(w_pos)
    sweep(mat, 2L, as.numeric(crossprod(w_c, mat)), "-")
  } else {
    sweep(mat, 2L, colMeans(mat), "-")
  }
}

# Propagate observation-level AR(p) (+ optional rWLS transient) noise covariance
# through a linear trial readout/de-mixing operator B (M_trials x m_obs) in O(m_obs * p * M)
# operations without forming a dense m_obs x m_obs covariance matrix.
# Under Woodbury rWLS (Sigma = Sigma_sub + diag(delta^2)), U_0 = W_unit^{-T} B_eff^T
# is updated in O(m * k * M) via the thin SVD of A = W_unit E_S D_S^{1/2}:
#   U = U_0 + U_A diag(sqrt(1 + d_A^2) - 1) U_A^T U_0,
# satisfying U^T U = B_eff (Sigma_sub + diag(delta^2)) B_eff^T to machine precision.
compute_ar1_operator_factor <- function(B, rows, rho, center_obs = TRUE,
                                        obs_weights = NULL, weight_space = "woodbury") {
  if (!is.null(obs_weights) && !is.null(attr(obs_weights, "weight_space"))) {
    weight_space <- as.character(attr(obs_weights, "weight_space"))[1L]
  }
  if (length(rows) > 0L && ncol(B) > length(rows) && max(rows) <= ncol(B)) {
    B <- B[, rows, drop = FALSE]
  }
  m <- ncol(B)
  M <- nrow(B)
  if (m == 0L || M == 0L) return(matrix(0, m, M))
  has_obs_w <- !is.null(obs_weights) && length(obs_weights) == m &&
    all(is.finite(obs_weights)) && any(abs(obs_weights - 1.0) > 1e-8)
  w_use <- if (has_obs_w) pmin(pmax(as.numeric(obs_weights), 1e-12), 1.0) else NULL
  B_eff <- if (isTRUE(center_obs)) {
    if (has_obs_w) {
      w_pos <- pmax(as.numeric(obs_weights), 1e-8)
      w_c <- w_pos / sum(w_pos)
      B - tcrossprod(rowSums(B), w_c)
    } else {
      sweep(B, 1L, rowMeans(B), "-")
    }
  } else {
    B
  }
  Bt <- t(B_eff)
  rho1 <- if (length(rho)) rho[1L] else 0.0
  has_ar1 <- is.finite(rho1) && abs(rho1) >= 1e-10
  has_psi <- m >= 2L && length(rho) >= 2L && any(is.finite(rho[-1L]) & abs(rho[-1L]) >= 1e-10)
  rows_use <- if (length(rows) == m) as.integer(rows) else seq_len(m)
  U_0 <- if ((!has_ar1 && !has_psi) || m == 1L) {
    Bt
  } else if (has_psi) {
    is_contig <- all(diff(rows_use) == 1L)
    lag_span <- if (is_contig) max(8L, length(rho) + 2L) else max(8L, max(rows_use) - min(rows_use))
    info <- hrf_ar_cascade_info(rho, max_lag = lag_span)
    if (length(info$phi) == 0L) {
      Bt
    } else {
      hrf_ar_operator_factor_cpp(
        Bt, rows_use, rho = 0.0,
        ar_phi = info$phi, ar_acf = info$acf, sig_u2 = info$sig_u2
      )
    }
  } else {
    hrf_ar_operator_factor_cpp(Bt, rows_use, rho = as.numeric(rho1))
  }
  if (!has_obs_w) return(U_0)
  if (identical(weight_space, "innovation") || (!has_ar1 && !has_psi) || m == 1L) {
    return(sweep(U_0, 1L, 1.0 / sqrt(pmax(w_use, 1e-12)), "*"))
  }
  idx_s <- which(w_use < 1.0 - 1e-8)
  if (!length(idx_s)) return(U_0)
  delta_s <- sqrt(pmax(0.0, 1.0 / pmax(w_use[idx_s], 1e-12) - 1.0))
  E_s <- matrix(0.0, m, length(idx_s))
  E_s[cbind(idx_s, seq_along(idx_s))] <- delta_s
  A_mat <- hrf_apply_whitening(E_s, rows_use, rho, unit_marginal = TRUE)
  sv_a <- svd(A_mat, nv = 0)
  pos_a <- sv_a$d > 1e-12
  if (any(pos_a)) {
    ua <- sv_a$u[, pos_a, drop = FALSE]
    grow_a <- sqrt(1.0 + sv_a$d[pos_a]^2) - 1.0
    U_0 <- U_0 + ua %*% sweep(crossprod(ua, U_0), 1L, grow_a, "*")
  }
  U_0
}

compute_ar1_operator_cov <- function(B, rows, rho, center_obs = TRUE,
                                     obs_weights = NULL, weight_space = "woodbury") {
  m <- ncol(B)
  M <- nrow(B)
  if (m == 0L || M == 0L) return(matrix(0, M, M))
  U <- compute_ar1_operator_factor(
    B, rows, rho, center_obs = center_obs,
    obs_weights = obs_weights, weight_space = weight_space
  )
  crossprod(U)
}

# Construct linear operators that estimate non-target event amplitudes (L_other)
# and joint lag-0 target amplitudes (L_target_0) from deconvolved activity.
# In rapid multi-phase trials (e.g., Cue -> Decision -> Feedback), preceding or
# following non-target events leak into the target epoch through the HRF/inverse
# resolution matrix R_all_0 = S_all_0 * D * G_all_0.
# 1. First, condition-level mean amplitudes A_mean_op are estimated across all
#    conditions simultaneously so condition baselines are properly separated.
# 2. Next, residual trial-to-trial fluctuations Res_all_op are inverted using a
#    column-normalized, second-order Neumann-refined ridge inverse (W_2) to
#    capture trial-specific amplitude variation in non-target events (and, for
#    impulse models at lag 0, joint target-trial amplitudes).
build_other_demix_operator <- function(S_all_0, DG_all_0, ev, o_idx, t_idx, demix_ridge,
    source_model = "impulse", demix_target_trials = TRUE,
    other_event_model = c("auto", "trial", "condition"), T_whiten = NULL) {
  other_event_model <- match.arg(other_event_model)
  ok_all <- rowSums(!is.finite(S_all_0)) == 0L
  idx_ok <- which(ok_all)
  if (!length(idx_ok) || !length(o_idx)) return(NULL)
  m_all <- nrow(ev)
  cond_all <- if (is.null(ev$condition)) rep("event", m_all) else as.character(ev$condition)
  cond_all[is.na(cond_all) | !nzchar(cond_all)] <- "other_unknown"
  cond_levels <- unique(cond_all[idx_ok])
  M_cond_full <- vapply(cond_levels, function(cl) as.numeric(cond_all == cl), numeric(m_all))
  if (is.null(dim(M_cond_full))) M_cond_full <- matrix(M_cond_full, nrow = m_all, ncol = length(cond_levels))
  M_cond_ok <- M_cond_full[idx_ok, , drop = FALSE]
  S_ok <- S_all_0[idx_ok, , drop = FALSE]
  R_ok <- S_ok %*% DG_all_0[, idx_ok, drop = FALSE]
  R_cond <- R_ok %*% M_cond_ok
  rc_diag <- diag(crossprod(R_cond))
  rc_scale <- if (any(is.finite(rc_diag) & rc_diag > 0)) mean(rc_diag[is.finite(rc_diag) & rc_diag > 0]) else 1.0
  lhs_cond <- crossprod(R_cond) + max(1e-4 * rc_scale, 1e-8) * diag(ncol(R_cond))
  cond_mean_op <- solve(lhs_cond, crossprod(R_cond, S_ok))
  A_mean_op <- M_cond_full %*% cond_mean_op
  A_mean_ok <- M_cond_ok %*% cond_mean_op

  ok_o <- o_idx[o_idx %in% idx_ok]
  L_other_mean <- A_mean_op[o_idx, , drop = FALSE]
  if (!length(ok_o) || !isTRUE(demix_target_trials) || identical(other_event_model, "condition")) {
    return(list(
      L_other = L_other_mean,
      L_target_0 = NULL,
      mode = "condition"
    ))
  }

  if (!is.null(T_whiten) && is.matrix(T_whiten) && nrow(T_whiten) >= 8L && ncol(T_whiten) == nrow(DG_all_0)) {
    F_ok <- T_whiten %*% DG_all_0[, idx_ok, drop = FALSE]
    F_cond <- F_ok %*% M_cond_ok
    fc_diag <- colSums(F_cond^2)
    fc_scale <- if (any(is.finite(fc_diag) & fc_diag > 0)) mean(fc_diag[is.finite(fc_diag) & fc_diag > 0]) else 1.0
    cm_white <- solve(crossprod(F_cond) + max(1e-4 * fc_scale, 1e-8) * diag(ncol(F_cond)), t(F_cond))
    cm_act_op <- cm_white %*% T_whiten
    A_mean_ts_op <- M_cond_full %*% cm_act_op
    A_mean_ts_ok <- M_cond_ok %*% cm_act_op

    fn_ok <- pmax(sqrt(colSums(F_ok^2)), 1e-6)
    F_u <- sweep(F_ok, 2L, fn_ok, "/")
    Gu <- crossprod(F_u)
    m_ok_ev <- ncol(F_u)
    df_ratio_ev <- m_ok_ev / max(1.0, nrow(T_whiten))
    r_eff <- max(0.035, min(0.10, demix_ridge * max(0.90, df_ratio_ev)))
    H_1 <- solve(Gu + r_eff * diag(m_ok_ev))
    H_2 <- H_1 + 0.70 * H_1 %*% (diag(m_ok_ev) - Gu %*% H_1)
    W_u_white <- H_2 %*% t(F_u)
    W_u_act <- W_u_white %*% T_whiten - (W_u_white %*% F_cond) %*% cm_act_op
    g_eff <- pmax(diag(H_2 %*% Gu), 0.10)
    Delta_ok_op <- sweep(W_u_act, 1L, g_eff * fn_ok, "/")
    A_joint_op <- A_mean_ts_op
    A_joint_op[idx_ok, ] <- A_mean_ts_ok + Delta_ok_op
    proj_resid_info <- list(T_whiten = T_whiten, F_u = F_u, H_1 = H_1, eff_df = sum(diag(Gu %*% H_1)))
  } else {
    proj_resid_info <- NULL
    Res_ok_op <- S_ok - R_ok %*% A_mean_ok
    diag_ok <- pmax(diag(R_ok), 0.05)
    C_ok <- sweep(R_ok, 2L, diag_ok, "/")
    r_eff <- if (demix_ridge < 0.025) max(demix_ridge, 0.015) else min(max(demix_ridge, 0.025), 0.035)
    W_1 <- solve(crossprod(C_ok) + r_eff * diag(ncol(C_ok)), t(C_ok))
    W_2 <- W_1 + 0.80 * W_1 %*% (diag(nrow(C_ok)) - C_ok %*% W_1)
    g_eff <- pmax(diag(W_2 %*% C_ok), 0.10)
    Delta_ok_op <- sweep(W_2 %*% Res_ok_op, 1L, g_eff * diag_ok, "/")
    A_joint_op <- A_mean_op
    A_joint_op[idx_ok, ] <- A_mean_ok + Delta_ok_op
  }
  use_target_0 <- identical(source_model, "impulse") && length(t_idx) > 0L && all(t_idx %in% idx_ok)
  L_other_trial <- A_joint_op[o_idx, , drop = FALSE]
  L_target_trial <- if (use_target_0) A_joint_op[t_idx, , drop = FALSE] else NULL

  if (identical(other_event_model, "trial")) {
    return(list(
      L_other = L_other_trial,
      L_target_0 = L_target_trial,
      mode = "trial"
    ))
  }

  o_pos_ok <- match(ok_o, idx_ok)
  L_delta_o <- Delta_ok_op[o_pos_ok, , drop = FALSE]
  t_ok <- t_idx[t_idx %in% idx_ok]
  K_o_from_t <- if (length(t_ok)) L_delta_o %*% DG_all_0[, t_ok, drop = FALSE] else NULL
  L_t_pilot <- if (length(t_ok)) A_joint_op[t_ok, , drop = FALSE] else NULL
  list(
    L_other = L_other_mean,
    L_other_mean = L_other_mean,
    L_other_trial = L_other_trial,
    L_target_0 = NULL,
    L_target_trial = L_target_trial,
    L_delta_o = L_delta_o,
    ok_o = ok_o,
    cond_o = cond_all[ok_o],
    K_o_from_t = K_o_from_t,
    L_t_pilot = L_t_pilot,
    proj_resid_info = proj_resid_info,
    mode = "auto"
  )
}

# Compute a variance-preserving Zero-Phase Component Analysis (ZCA) whitening
# matrix Phi for a within-run trial covariance matrix Sigma_tot (M x M).
# Let D = diag(sqrt(diag(Sigma_tot))) and R_tot = D^{-1} Sigma_tot D^{-1} be the
# trial correlation matrix. Then Phi = D * R_tot^{-1/2} * D^{-1} (normalized to
# unit mean diagonal) removes off-diagonal trial-to-trial error correlations
# while preserving each trial's marginal error variance diag(Sigma_tot).
# This allows downstream mixed models (`mixed_by(..., weights = "precision_weight")`)
# to pair ZCA-decorrelated trial observations with diagonal precision weights.
compute_zca_whitening_matrix <- function(Sigma_tot, eig_floor = 0.08, offdiag_scale = 0.90) {
  M <- nrow(Sigma_tot)
  if (M <= 1L) return(diag(M))
  Sigma_sym <- 0.5 * (Sigma_tot + t(Sigma_tot))
  mean_mode_var <- sum(Sigma_sym) / M
  cent_mean_var <- max((sum(diag(Sigma_sym)) - mean_mode_var) / max(1L, M - 1L), .Machine$double.eps)
  # Stabilize the rank-1 run-mean null mode so run-centered covariance matrices
  # have exact centered reference scale along 1_M and are strictly invariant to c * (1_M 1_M' / M).
  if (is.finite(mean_mode_var) && mean_mode_var < cent_mean_var) {
    Sigma_sym <- Sigma_sym + (cent_mean_var - mean_mode_var) * matrix(1.0 / M, M, M)
  }
  d_tot <- sqrt(pmax(diag(Sigma_sym), .Machine$double.eps))
  R_tot <- sweep(sweep(Sigma_sym, 1L, d_tot, "/"), 2L, d_tot, "/")
  R_tot <- 0.5 * (R_tot + t(R_tot))
  if (offdiag_scale < 1.0) {
    R_tot <- (1.0 - offdiag_scale) * diag(M) + offdiag_scale * R_tot
  }
  eig <- eigen(R_tot, symmetric = TRUE)
  lam_pos <- pmax(eig$values, 0.0)
  inv_sqrt_vals <- 1.0 / ((lam_pos^2 + eig_floor^2)^0.25)
  R_inv_half <- eig$vectors %*% (inv_sqrt_vals * t(eig$vectors))
  Phi <- sweep(sweep(R_inv_half, 1L, d_tot, "*"), 2L, d_tot, "/")
  tr_mean <- mean(diag(Phi))
  if (is.finite(tr_mean) && abs(tr_mean) > 1e-6) {
    Phi <- Phi / tr_mean
  }
  Phi
}

# Resolve a sparse spatial graph for bilateral Empirical-Bayes variance shrinkage in Step 3.
resolve_epoch_spatial_graph <- function(deconvolution, fit, spatial_graph,
                                        aggregate_by_parcel, n_units, unit_ids,
                                        parcel_vox_list = NULL, parcel_weights_list = NULL) {
  if (!is.null(spatial_graph)) {
    if (inherits(spatial_graph, "medusa_spatial_graph")) {
      if (as.integer(spatial_graph$n_voxels) != as.integer(n_units)) {
        stop("spatial_graph node count (", spatial_graph$n_voxels,
             ") does not match expected unit count (", n_units, ").", call. = FALSE)
      }
      return(spatial_graph)
    }
    return(build_spatial_graph(
      spatial_graph,
      n_voxels = n_units,
      voxel_ids = unit_ids
    ))
  }
  if (isTRUE(aggregate_by_parcel)) {
    if (is.null(parcel_vox_list) || length(parcel_vox_list) != n_units || n_units < 2L) {
      return(NULL)
    }
    vox_coords <- NULL
    vm <- deconvolution$voxel_map
    n_vox <- length(fit$voxel_ids)
    if (!is.null(vm) && is.data.frame(vm) && nrow(vm) == n_vox) {
      if (all(c("x", "y", "z") %in% names(vm))) {
        vox_coords <- as.matrix(vm[, c("x", "y", "z"), drop = FALSE])
      } else if (all(c("i", "j", "k") %in% names(vm))) {
        sp <- hrf_infer_lattice_spacing(vm, deconvolution$geometry)
        ijk_m <- as.matrix(vm[, c("i", "j", "k"), drop = FALSE])
        storage.mode(ijk_m) <- "double"
        vox_coords <- sweep(ijk_m, 2L, as.numeric(sp), "*")
      }
    }
    if (is.null(vox_coords) && !is.null(fit$specification$spatial_coords) &&
        is.matrix(fit$specification$spatial_coords) &&
        nrow(fit$specification$spatial_coords) == n_vox) {
      vox_coords <- fit$specification$spatial_coords
    }
    if (is.null(vox_coords) || any(!is.finite(vox_coords))) {
      return(NULL)
    }
    storage.mode(vox_coords) <- "double"
    centroids <- matrix(
      NA_real_, nrow = n_units, ncol = ncol(vox_coords),
      dimnames = list(unit_ids, colnames(vox_coords))
    )
    for (g in seq_len(n_units)) {
      vg <- parcel_vox_list[[g]]
      wg <- if (!is.null(parcel_weights_list) && length(parcel_weights_list) >= g) {
        parcel_weights_list[[g]]
      } else {
        rep(1.0 / length(vg), length(vg))
      }
      wg <- wg / sum(wg)
      centroids[g, ] <- as.numeric(crossprod(wg, vox_coords[vg, , drop = FALSE]))
    }
    return(tryCatch(
      build_spatial_graph(
        centroids,
        metric = "euclidean",
        voxel_ids = unit_ids
      ),
      error = function(e) NULL
    ))
  }
  if (inherits(fit$specification$spatial_graph, "medusa_spatial_graph") &&
      fit$specification$spatial_graph$n_voxels == n_units) {
    return(fit$specification$spatial_graph)
  }
  if (!is.null(deconvolution$voxel_map) && is.data.frame(deconvolution$voxel_map) &&
      nrow(deconvolution$voxel_map) == n_units) {
    return(tryCatch(
      build_spatial_graph(
        list(voxel_map = deconvolution$voxel_map, geometry = deconvolution$geometry),
        spatial_groups = fit$specification$spatial_groups,
        parcel_boundary = "soft",
        voxel_ids = unit_ids
      ),
      error = function(e) NULL
    ))
  }
  if (!is.null(fit$specification$spatial_coords)) {
    return(tryCatch(
      build_spatial_graph(
        fit$specification$spatial_coords,
        spatial_fwhm = fit$settings$spatial_fwhm,
        spatial_groups = fit$specification$spatial_groups,
        parcel_boundary = "soft",
        voxel_ids = unit_ids
      ),
      error = function(e) NULL
    ))
  }
  NULL
}

# Non-spatial Global Empirical-Bayes shrinkage of untruncated dimensionless biological
# variance ratios theta_raw = tau2_raw / s2_obs across analysis units (with Isserlis
# residual-correlation adjustment and Efron-Morris limited-translation protection).
eb_shrink_tau2_global_r <- function(theta_raw, d_eff, U_resid = NULL, d_resid = 30.0) {
  V <- length(theta_raw)
  ok <- which(is.finite(theta_raw) & is.finite(d_eff) & d_eff > 0)
  if (length(ok) < 2L) {
    return(list(
      theta_eb = pmax(0.0, theta_raw),
      shrinkage_b = rep(0.0, V),
      kappa_minus = rep(1.0, V),
      c_self_prior = rep(0.0, V),
      weights = rep(0.0, V),
      ok = ok
    ))
  }
  d_res <- max(3.0, as.numeric(d_resid)[1L])
  w_raw <- rep(0.0, V)
  w_raw[ok] <- sqrt(pmax(2.0, d_eff[ok]))
  w_tot <- sum(w_raw)
  w_glob <- w_raw / w_tot
  th_all <- sum(w_glob[ok] * theta_raw[ok])

  # Leave-one-out prior mean and exact weights
  denom_loo <- pmax(w_tot - w_raw[ok], 1e-8)
  th0_ok <- (w_tot * th_all - w_raw[ok] * theta_raw[ok]) / denom_loo
  theta_prior <- rep(NA_real_, V)
  theta_prior[ok] <- th0_ok

  # Estimate global Isserlis squared residual correlation cleanly on up to 128 units
  kap_glob <- sum(w_glob[ok]^2)
  c_self_vec <- rep(0.0, V)
  if (!is.null(U_resid) && is.matrix(U_resid) && ncol(U_resid) == V && nrow(U_resid) >= 3L) {
    sub_idx <- if (length(ok) > 128L) {
      ok[unique(as.integer(round(seq(1, length(ok), length.out = 128L))))]
    } else {
      ok
    }
    U_sub <- U_resid[, sub_idx, drop = FALSE]
    nrm2_sub <- colSums(U_sub^2)
    good_s <- which(nrm2_sub > 0.25)
    if (length(good_s) >= 2L) {
      u_s <- sub_idx[good_s]
      ws <- w_raw[u_s] / sum(w_raw[u_s])
      G_s <- crossprod(U_resid[, u_s, drop = FALSE])
      R2_s <- matrix(
        pmin(pmax((d_res * (G_s * G_s) - 1.0) / (d_res - 1.0), 0.0), 1.0),
        length(u_s), length(u_s)
      )
      diag(R2_s) <- 1.0
      kap_glob <- min(1.0, max(sum(ws^2), as.numeric(crossprod(ws, R2_s %*% ws))))

      # Mean squared correlation of every valid unit with the reference pool
      C_all <- crossprod(U_resid[, ok, drop = FALSE], U_resid[, u_s, drop = FALSE])
      R2_all <- matrix(
        pmin(pmax((d_res * (C_all * C_all) - 1.0) / (d_res - 1.0), 0.0), 1.0),
        length(ok), length(u_s)
      )
      c_ok <- as.numeric(R2_all %*% ws)
      # Remove self-contribution when unit is in u_s
      in_sub <- match(ok, u_s)
      hit <- which(!is.na(in_sub))
      if (length(hit)) {
        ws_hit <- ws[in_sub[hit]]
        c_ok[hit] <- pmax(0.0, (c_ok[hit] - ws_hit) / pmax(1.0 - ws_hit, 1e-6))
      }
      c_self_vec[ok] <- pmin(1.0, pmax(0.0, c_ok))
    }
  }
  kappa_minus <- rep(kap_glob, V)

  de_ok <- pmax(2.0, d_eff[ok])
  th_lvl <- pmax(0.0, th0_ok)
  vt_ok <- (2.0 / de_ok) * (1.0 + th_lvl)^2
  vd_ok <- pmax(0.05 * vt_ok, vt_ok * pmax(0.0, 1.0 + kap_glob - 2.0 * c_self_vec[ok]))
  diff_ok <- theta_raw[ok] - th0_ok
  om2_glob <- max(0.0, stats::median((diff_ok^2) / 0.4549364 - vd_ok, na.rm = TRUE))
  om_eff <- pmax(om2_glob, (diff_ok^2) / 6.25 - vd_ok)
  b_ok <- pmin(0.90, pmax(0.0, (vt_ok * pmax(0.0, 1.0 - c_self_vec[ok])) / pmax(vd_ok + om_eff, 1e-12)))

  theta_eb <- rep(NA_real_, V)
  theta_eb[ok] <- pmax(0.0, (1.0 - b_ok) * theta_raw[ok] + b_ok * th0_ok)
  shrinkage_b <- rep(0.0, V)
  shrinkage_b[ok] <- b_ok

  list(
    theta_eb = theta_eb,
    theta_prior = theta_prior,
    shrinkage_b = shrinkage_b,
    kappa_minus = kappa_minus,
    c_self_prior = c_self_vec,
    weights = w_raw,
    ok = ok
  )
}

eb_moderate_scale_global_r <- function(scale_raw, d_eff, df_nom, eb_fit) {
  V <- nrow(scale_raw)
  S <- ncol(scale_raw)
  scale_mod <- matrix(pmax(scale_raw, 0.25), V, S)
  df_mod <- as.numeric(df_nom)
  ok <- eb_fit$ok
  if (length(ok) < 2L) return(list(scale_mod = scale_mod, df_mod = df_mod))

  w_raw <- eb_fit$weights
  kap <- eb_fit$kappa_minus[ok]
  c_u0 <- eb_fit$c_self_prior[ok]
  de_ok <- pmax(2.0, d_eff[ok])
  psi_ref <- rep(1.0, length(ok))

  for (col in seq_len(S)) {
    sc_col <- scale_raw[, col]
    ok_c <- ok[is.finite(sc_col[ok])]
    if (length(ok_c) < 2L) next
    wc <- w_raw[ok_c]
    wtot <- sum(wc)
    s_all <- sum((wc / wtot) * sc_col[ok_c])
    s0 <- (wtot * s_all - wc * sc_col[ok_c]) / pmax(wtot - wc, 1e-8)

    kap_c <- eb_fit$kappa_minus[ok_c]
    cu0_c <- eb_fit$c_self_prior[ok_c]
    de_c <- pmax(2.0, d_eff[ok_c])
    vs_c <- (2.0 / de_c) * (pmax(0.25, s0)^2)
    vd_c <- pmax(0.05 * vs_c, vs_c * pmax(0.0, 1.0 + kap_c - 2.0 * cu0_c))
    diff_c <- sc_col[ok_c] - s0
    om2_c <- max(0.0, stats::median((diff_c^2) / 0.4549364 - vd_c, na.rm = TRUE))
    om_eff <- pmax(om2_c, (diff_c^2) / 6.25 - vd_c)
    bc <- pmin(0.85, pmax(0.0, (vs_c * pmax(0.0, 1.0 - cu0_c)) / pmax(vd_c + om_eff, 1e-12)))
    scale_mod[ok_c, col] <- pmax(0.25, (1.0 - bc) * sc_col[ok_c] + bc * s0)
    if (col == S) {
      psi_c <- pmin(1.0, pmax(0.15, (1.0 - bc)^2 + 2.0 * bc * (1.0 - bc) * cu0_c + (bc^2) * kap_c))
      df_gain <- de_c / pmax(0.25, psi_c)
      df_mod[ok_c] <- pmax(df_nom[ok_c], pmin(2.5 * df_nom[ok_c], df_gain))
    }
  }
  list(scale_mod = scale_mod, df_mod = df_mod)
}

#' Operator-based epoch de-mixing and gain-corrected trial regression across lags
#'
#' Estimate trial-by-trial correlations and regression slopes across epoch lags
#' relative to a target event while removing preceding/neighboring event leakage,
#' correcting regularization-induced amplitude attenuation, and propagating
#' per-trial operator precision weights and within-run error covariance.
#'
#' @param deconvolution Result of \code{deconvolve_crossrun_hrf()}.
#' @param covariates Optional trial-level predictors for target trials, supplied
#'   as: (a) a character vector of numeric/logical column names present in
#'   \code{events} (e.g., \code{c("rpe", "abs_rpe")}), extracted automatically at
#'   the target event rows in guaranteed lockstep with \code{event_index}; (b) a
#'   one-sided formula evaluated on target event rows (e.g.,
#'   \code{~ rpe + abs_rpe}); (c) a numeric vector, matrix, or \code{data.frame}
#'   concatenated across target trials in fitted run order; or (d) a run-named
#'   list matching \code{deconvolution$hrf_fit$run_ids}.
#'   When \code{NULL} (default), performs covariate-agnostic epoch readout,
#'   de-mixing, biological variance estimation (\code{tau2_mode = "marginal"}),
#'   Empirical-Bayes precision weighting, and optional within-run ZCA whitening
#'   (\code{whiten_panel = TRUE}), returning the trial panel (\code{return_panel}
#'   defaults to \code{TRUE} when \code{covariates = NULL}) for downstream
#'   two-stage \code{\link{mixed_by}} behavioral modeling or inter-region
#'   connectivity (\code{R2 ~ R1_seed * x}).
#' @param events Optional combined multi-run \code{data.frame} (split by
#'   \code{run_col}) or run-named list of event \code{data.frame}s with
#'   \code{onset} and optional \code{condition}, \code{duration}, and behavioral
#'   covariate columns. Defaults to the HRF fit's event tables.
#' @param target_condition Optional character scalar naming the target event
#'   condition in \code{events[[r]]$condition} (e.g., \code{"feedback"}). When
#'   supplied, non-target events in \code{events[[r]]} are treated as interfering
#'   events whose leakage into the target epoch is removed when
#'   \code{demix_other_events = TRUE}.
#' @param lags Numeric vector of epoch readout lags (in seconds) relative to each
#'   target event onset.
#' @param readout_kernel Epoch readout kernel: \code{"linear"}, \code{"boxcar"},
#'   \code{"triangle"}, or \code{"gaussian"}. Default: \code{"triangle"}.
#' @param readout_width Optional kernel full width in seconds for non-linear
#'   readouts (defaults to \code{1.5 * tr}).
#' @param demix_other_events Logical; whether to estimate non-target event
#'   contributions at their own onsets and subtract their operator-predicted
#'   leakage into each target epoch lag. When \code{demix_ridge >= 0.03}
#'   (default \code{0.05}), uses a condition-mean non-target de-mixing operator
#'   that removes preceding cue/event bleed without inflating target-trial
#'   variance; when \code{demix_ridge < 0.03} (e.g. \code{0.02}), activates a
#'   2nd-order iterated Tikhonov joint \eqn{(T, O)} operator to de-mix strongly
#'   trial-varying non-target confounds. Default: \code{TRUE}.
#' @param demix_target_trials Logical; whether to use the stratum resolution
#'   matrix at each lag to de-mix adjacent target-trial overlap and restore
#'   regularization-attenuated amplitude gain via 2nd-order iterated Tikhonov
#'   inversion (\code{solve_duration_demix_operator()}). Default: \code{TRUE}.
#' @param demix_within_trial Logical; whether to apply an additional multi-lag
#'   point-spread de-mixing step across the requested \code{lags} (when at least
#'   3 lags are evaluated) to sharpen temporal localization across lags within
#'   each target trial. Recommended (\code{TRUE}) when estimating within-trial
#'   temporal waveforms \eqn{\beta(\tau)} across a multi-lag grid. Default:
#'   \code{FALSE}.
#' @param demix_ridge Positive ridge stabilization parameter for target and
#'   non-target cross-talk inversion. Values \code{>= 0.03} (default \code{0.05})
#'   use condition-mean non-target de-mixing; values \code{< 0.03} activate
#'   joint trial-varying non-target de-mixing. Default: \code{0.05}.
#' @param within_trial_ridge Positive ridge stabilization parameter for
#'   multi-lag within-trial point-spread inversion when
#'   \code{demix_within_trial = TRUE}. Default: \code{0.06}.
#' @param gain_floor Minimum absolute operator gain used when correcting
#'   amplitude slopes. Default: \code{0.05}.
#' @param restore_scale Logical; if \code{TRUE}, multiplies each run and voxel by
#'   its recorded \code{deconvolution$scaling[[r]]$scale} divisor so that
#'   \code{beta_raw} and \code{beta_corrected} (and panel readouts) are reported
#'   in unscaled physical input units rather than standardized
#'   (\code{zscore}/\code{noise_sd}) units. Strongly recommended (\code{TRUE})
#'   when pooling across multiple runs in \code{\link{mixed_by}}. Default:
#'   \code{FALSE}.
#' @param source_model Forward neural source model at each target lag:
#'   \code{"impulse"} (default unit impulse at \code{onset + lag}),
#'   \code{"window"} (local sustained boxcar of width \code{readout_width}
#'   centered at \code{onset + lag}), or \code{"duration"} (trial-varying
#'   sustained boxcar over \code{[onset + lag, onset + lag + d_i]} using each
#'   target trial's duration \code{d_i} from \code{trial_durations} or
#'   \code{events[[r]]$duration}, with duration-matched readout and unpenalized
#'   baseline-rate de-confounding).
#' @param trial_durations Optional single column name in \code{events} (e.g.,
#'   \code{"rt"}), numeric vector (concatenated across runs), or run-named list
#'   of non-negative target trial durations (in seconds), such as per-trial
#'   reaction times. Overrides \code{events[[r]]$duration} for target trials and
#'   defaults \code{source_model} to \code{"duration"} when \code{source_model}
#'   is not explicitly specified.
#' @param duration_scale Scaling convention for positive-duration target trials
#'   when \code{source_model = "duration"}: \code{"rate"} (default; trial
#'   amplitude represents neural firing rate / boxcar height per TR bin, so
#'   integrated neural activity scales with duration and constant baseline firing
#'   rates are de-confounded from duration variability) or \code{"integral"}
#'   (each trial boxcar is normalized to unit time-integrated neural mass across
#'   its duration).
#' @param weighting Trial regression weighting scheme across target trials and
#'   runs: \code{"precision"} (default; Empirical-Bayes / floor-stabilized
#'   inverse error-variance weights propagating per-trial operator standard
#'   errors and run noise scale), \code{"gls"} (null-mode-damped,
#'   Bartlett-adaptive within-run ZCA correlation whitening plus precision
#'   weighting accounting for off-diagonal trial-to-trial error correlations),
#'   or \code{"none"} (unweighted OLS). All three slope estimates
#'   (\code{beta_unweighted}, \code{beta_weighted}, \code{beta_gls}) and their
#'   standard errors are always included in \code{coefficients}.
#' @param tau2_mode Estimator for biological trial-to-trial variance
#'   \eqn{\hat{\tau}^2_{\mathrm{bio}}} used in Empirical-Bayes trial precision
#'   weights and within-run ZCA whitening: \code{"auto"} (default; uses
#'   \code{"marginal"} when \code{covariates = NULL} and an Empirical-Bayes
#'   blend of residual and marginal estimators when covariates are supplied),
#'   \code{"residual"} (method-of-moments estimator from OLS covariate regression
#'   residuals minus projected operator variance \eqn{\mathrm{tr}((\mathbf{I}-\mathbf{P}_X)\boldsymbol{\Omega})}),
#'   or \code{"marginal"} (covariate-agnostic estimator combining run-centered
#'   marginal variance, lag-1 trial difference variance, and robust MAD scale
#'   minus analytical operator variance expectations).
#' @param tau2_shrinkage Empirical-Bayes shrinkage across analysis units for
#'   biological trial variance \eqn{\hat{\tau}^2_{\mathrm{bio},v}} and sandwich
#'   residual scale factors: \code{"auto"} (default; enables bilateral
#'   \code{"graph_eb"} automatically for voxelwise/vertexwise map workflows with
#'   \eqn{\ge 6} units whenever spatial coordinates, \code{voxel_map}, or
#'   \code{spatial_graph} are available, and uses \code{"none"} for parcel-level
#'   runs unless \code{spatial_graph} is supplied explicitly),
#'   \code{"graph_eb"} (2-pass bilateral Empirical Bayes over a sparse spatial
#'   neighborhood graph with Isserlis spatial-correlation accounting and
#'   Efron-Morris limited translation; for \code{aggregate_by_parcel = TRUE},
#'   automatically constructs a precision-weighted 3D parcel-centroid Euclidean
#'   graph when \code{spatial_graph} is \code{NULL}), \code{"eb"} (non-spatial
#'   global Empirical Bayes across units with residual-correlation accounting),
#'   or \code{"none"} (independent single-unit method-of-moments estimates).
#' @param spatial_graph Optional \code{"medusa_spatial_graph"} object (from
#'   \code{\link{build_spatial_graph}}) or spatial input accepted by
#'   \code{\link{build_spatial_graph}} used when \code{tau2_shrinkage} resolves
#'   to \code{"graph_eb"}. When \code{NULL}, inferred automatically from
#'   \code{deconvolution$hrf_fit$specification$spatial_graph},
#'   \code{deconvolution$voxel_map}, or \code{spatial_coords} when present
#'   (using precision-weighted 3D parcel centroids when
#'   \code{aggregate_by_parcel = TRUE}).
#' @param aggregate_by_parcel Logical; if \code{TRUE}, aggregates deconvolved
#'   signals across voxels within each spatial parcel (using inverse-variance
#'   voxel precision weights and empirical spatial noise correlation propagation)
#'   prior to Stage-2 epoch de-mixing and Stage-3 regression/ZCA whitening.
#'   Replaces voxel-by-voxel Stage-3/4 fitting with parcel-level inference while
#'   calibrating \code{se_demixed}, \code{tau2_bio}, \code{precision_weight},
#'   and ZCA whitening directly on the high-SNR parcel series. Default:
#'   \code{FALSE}.
#' @param parcel_groups Optional character or factor vector of parcel/region
#'   labels for each voxel (either unnamed of length \code{n_voxels} or named by
#'   \code{voxel_ids}). Defaults to
#'   \code{deconvolution$hrf_fit$specification$spatial_groups} when
#'   \code{aggregate_by_parcel = TRUE}. Supplying a non-\code{NULL}
#'   \code{parcel_groups} automatically enables \code{aggregate_by_parcel = TRUE}.
#' @param whiten_panel Logical; if \code{TRUE} and \code{return_panel = TRUE},
#'   adds within-run null-mode-damped, Bartlett-adaptive ZCA-whitened outcome
#'   (\code{decon_demixed_whitened}) and covariate (\code{<var>_whitened})
#'   columns to \code{panel} so downstream
#'   \code{mixed_by(..., weights = "precision_weight")} models account for both
#'   heteroscedastic trial precision and off-diagonal within-run trial error
#'   correlations. Recommended (\code{TRUE}) for group mixed-effects inference.
#'   Default: \code{FALSE}.
#' @param return_epochs Logical; whether to include raw and de-mixed trial-by-voxel
#'   (or trial-by-parcel) epoch matrices per lag. Default: \code{FALSE}.
#' @param return_panel Logical; whether to return a tidy MEDUSA-compatible
#'   \code{data.table} (\code{panel}) with columns \code{run}, \code{trial},
#'   \code{event_index}, \code{onset}, \code{evt_time}, \code{voxel} (plus
#'   \code{id}/\code{session} when present on \code{events}, and \code{parcel}
#'   when \code{aggregate_by_parcel = TRUE}),
#'   \code{decon_raw}, \code{decon_demixed}, \code{decon_eb}, \code{se_raw},
#'   \code{se_demixed}, \code{precision_weight}, \code{precision_weight_cohort},
#'   \code{precision_weight_op}, \code{precision_weight_raw}, \code{error_ac1},
#'   \code{tau2_bio}, the covariate columns, and (under
#'   \code{shared_error_method = "poet_covreg"}) Block-PAF list-columns
#'   \code{shared_noise_vec}, \code{lag_noise_vec} (multi-lag fits), and
#'   \code{trial_cor_vec} (plus legacy \code{shared_noise_loading}) for
#'   downstream \code{\link{mixed_by}}, \code{\link{medusa_multilag_infer}},
#'   \code{\link{aggregate_decon_panel}}, or \code{\link{decon_connectivity_repair}}
#'   modeling. Defaults to \code{TRUE} when \code{covariates = NULL} and
#'   \code{return_epochs = FALSE}.
#' @param run_col Column name in \code{events} identifying runs when \code{events}
#'   is supplied as a combined multi-run \code{data.frame}. Default \code{"run"}.
#' @param variance_confounds Optional frame-level variance confound envelopes
#'   (e.g., framewise displacement, \code{std_dvars}/\code{dvars}, RVT) for
#'   LTZ-gated heteroskedastic covariance regression with exact Schur-complement
#'   restoration when confound columns are gated out. Accepts a
#'   \code{"medusa_confounds"} object from \code{\link{read_fmriprep_confounds}},
#'   a \code{"medusa_nifti_ts"} object with \code{$confounds} from
#'   \code{\link{read_fmriprep_subject}}, a per-run list of numeric
#'   matrices/\code{data.frame}s/vectors, a multi-run \code{data.frame} containing
#'   \code{run_col}, or (for single-run fits) a numeric matrix/\code{data.frame}/vector.
#' @param shared_error_method Shared cross-unit and dynamic error covariance
#'   estimator: \code{"poet_covreg"} (default; v3 Subspace-Debiased Adaptive-POET
#'   + Multi-Lag Toeplitz & General-\eqn{W_s} Covariance Regression),
#'   \code{"legacy"} (stationary innovation shrinkage), or \code{"none"}.
#' @param bold Optional list of per-run BOLD matrices (not required when
#'   \code{deconvolution} comes from current \code{\link{deconvolve_crossrun_hrf}},
#'   which caches observed/sketched prewhitened series automatically).
#' @param other_event_model Parameterization of non-target events during
#'   \code{demix_other_events = TRUE}: \code{"auto"} (default; tests for
#'   trial-level amplitude modulation on adjacent non-target events via a
#'   spatially correlated Satterthwaite \eqn{F}-screen and switches from
#'   condition-mean to per-trial projection only when supported), \code{"trial"}
#'   (one nuisance column per non-target trial), or \code{"condition"} (one
#'   nuisance column per non-target condition).
#' @return A list with \code{coefficients} (per voxel/parcel, lag, and predictor
#'   term), \code{summary} (averaged across usable units per lag and term),
#'   optional \code{epochs}, optional \code{panel}, \code{lag_resolution},
#'   \code{lag_error_cor}, \code{unit_error_cor}, \code{shared_error_model}
#'   (when \code{shared_error_method = "poet_covreg"}), optional
#'   \code{parcel_weights}, and \code{settings}.
#' @details
#'   \code{decon_epoch_regression()} bridges continuous regularized deconvolution
#'   and trial-level cognitive neuroscience models in five stages:
#'   \enumerate{
#'     \item \strong{Kernel Quadrature Readout:} At each requested lag \eqn{\tau},
#'       continuous deconvolved activity is integrated around each target event
#'       using \code{readout_kernel}. For \code{source_model = "duration"}, the
#'       readout window centers on the trial's midpoint
#'       \eqn{t_i + \tau + 0.5 d_i} and expands to \eqn{\max(w_{\mathrm{src}}, d_i)}.
#'     \item \strong{Non-Target Event De-Mixing (\code{demix_other_events}):}
#'       Using the frozen inverse operator \eqn{\mathbf{D}_r} and forward HRF
#'       design \eqn{\mathbf{G}_{r,0}}, the function subtracts resolution-predicted
#'       leakage from non-target events (e.g., cues or motor responses) using a
#'       condition-mean operator (\code{demix_ridge >= 0.03}) or a 2nd-order
#'       Tikhonov joint \eqn{(T, O)} operator (\code{demix_ridge < 0.03}).
#'     \item \strong{Target Cross-Talk & Duration De-Mixing (\code{demix_target_trials}):}
#'       Inverts the lag-specific target resolution matrix via a 2nd-order
#'       iterated Tikhonov operator (\code{solve_duration_demix_operator()}) with
#'       per-trial effective gain restoration to remove trial-to-trial overlap
#'       and divide out regularization-induced amplitude shrinkage. In
#'       \code{source_model = "duration"}, a centered-subspace penalty and
#'       rank-1 baseline adjustment keep constant baseline firing rates invariant
#'       to trial duration (preventing reaction-time confounding).
#'     \item \strong{Optional Multi-Lag Point-Spread De-Mixing (\code{demix_within_trial}):}
#'       When enabled across \eqn{\ge 3} lags, sharpens temporal specificity across
#'       lags using a 2nd-order Tikhonov lag-by-lag resolution inverse while
#'       preserving cross-trial operator covariance structure.
#'     \item \strong{Analytical Uncertainty, Adaptive-POET Covariance Regression & ZCA Whitening:}
#'       Propagates observation-level AR(1) noise covariance through the composite
#'       linear readout/de-mixing operator \eqn{\mathbf{B}_{r,\tau}} to obtain exact
#'       per-trial standard errors (\code{se_demixed}), adjacent-trial error
#'       autocorrelation (\code{error_ac1}), Empirical-Bayes precision weights
#'       (\code{precision_weight}, \code{precision_weight_cohort}), covariance-propagated
#'       sandwich SEs (\code{se_weighted}), and null-mode-damped, Bartlett-adaptive
#'       within-run ZCA-whitened variables for \code{\link{mixed_by}}. Under
#'       \code{shared_error_method = "poet_covreg"}, fits a subspace-debiased
#'       Adaptive-POET + Tukey-Hanning multi-lag Toeplitz and LTZ-gated dynamic
#'       covariance-regression model (automatically switching above 512 units to
#'       a voxel-safe dual-space factor representation that avoids \eqn{V \times V}
#'       allocations) and embeds per-trial shared-error factors onto \code{panel}
#'       for \code{\link{decon_connectivity_repair}} and \code{\link{aggregate_decon_panel}}.
#'   }
#' @section Recommended Best-Practice Configuration:
#'   Across the MEDuSA Monte Carlo validation suite (amplitude recovery,
#'   temporal alignment, inter-region connectivity, factorial power surfaces, and
#'   biophysical Windkessel-Balloon + \code{fmrisim} scanner realism), the
#'   recommended end-to-end configuration follows one shared Steps 0--3 backbone
#'   with estimand-specific downstream modes:
#'   \enumerate{
#'     \item \strong{Stage 0--1 (Ingestion & HRF Estimation):} Ingest BIDS /
#'       \code{fMRIPrep} data via \code{\link{read_fmriprep_subject}} (or
#'       \code{\link{read_fmriprep_confounds}}, which extracts frame-level
#'       \code{$variance_confounds} automatically). Use
#'       \code{hrf_library(tr, preset = "medusa")} and
#'       \code{\link{estimate_crossrun_hrf}} with \code{noise_model = "ar1"},
#'       \code{evidence = "calibrated"}, \code{selection_mode = "soft_average"},
#'       and parcel labels in \code{spatial_groups}.
#'     \item \strong{Stage 2 (Deconvolution):} Use \code{\link{deconvolve_crossrun_hrf}}
#'       with \code{lambda = NULL} (adaptive prewhitened GCV tuning),
#'       \code{scaling = "noise_sd"}, \code{penalty = "sobolev"}, and
#'       \code{rwls = "auto"} (automatically caching \code{bold_obs} up to 16,384
#'       voxels or a compact \code{whitened_sketch} above 16,384 voxels so
#'       \code{bold} need not be passed again in Stage 3).
#'     \item \strong{Stage 3 (Epoch De-Mixing, Parcel Aggregation & Shared-Error Engine):} Use
#'       \code{decon_epoch_regression(..., demix_other_events = TRUE,
#'       demix_target_trials = TRUE, restore_scale = TRUE, weighting = "gls",
#'       aggregate_by_parcel = TRUE, shared_error_method = "poet_covreg",
#'       variance_confounds = conf, whiten_panel = TRUE, return_panel = TRUE)}.
#'       If you first fit voxelwise panels (\code{aggregate_by_parcel = FALSE}),
#'       you can collapse them post-hoc to any parcellation with exact shared-error
#'       propagation via \code{\link{aggregate_decon_panel}}. Leave
#'       \code{demix_within_trial = FALSE} when feeding \code{\link{medusa_multilag_infer}}
#'       (enable \code{TRUE} only for single-subject lag sharpening), and set
#'       \code{source_model = "duration"} (\code{duration_scale = "rate"}) when
#'       modeling variable-duration (RT) epochs.
#'     \item \strong{Stage 4--6 (Group Mixed-Effects, Multi-Lag & Connectivity Inference):}
#'       For scalar behavioral models, pass \code{weights = "precision_weight_cohort"}
#'       (and matched \code{*_whitened} columns when desired) to \code{\link{mixed_by}}.
#'       For peri-event trajectories, use \code{\link{medusa_multilag_infer}} with
#'       \code{contrast_source = "unsmoothed"}. For seed-to-target connectivity and
#'       PPI, pass \code{subject_epochs} or a stacked/CSV/Parquet \code{panel}
#'       (serialized via \code{medusa_write_table()} to preserve pipe-delimited
#'       shared-error list-columns) to
#'       \code{decon_connectivity_repair(..., method = "mmoms", shared_error = "auto",
#'       use_whitened = FALSE, crosstrial_se = TRUE)}.
#'   }
#' @seealso \code{\link{deconvolve_crossrun_hrf}}, \code{\link{aggregate_decon_panel}},
#'   \code{\link{decon_connectivity_repair}}, \code{\link{medusa_multilag_infer}},
#'   \code{\link{mixed_by}}
#' @export
decon_epoch_regression <- function(
  deconvolution,
  covariates = NULL,
  events = NULL,
  target_condition = NULL,
  lags = seq(-2, 4, by = 1),
  readout_kernel = c("triangle", "linear", "boxcar", "gaussian"),
  readout_width = NULL,
  demix_other_events = TRUE,
  demix_target_trials = TRUE,
  demix_within_trial = FALSE,
  demix_ridge = 0.05,
  within_trial_ridge = 0.06,
  gain_floor = 0.05,
  restore_scale = FALSE,
  source_model = c("impulse", "window", "duration"),
  trial_durations = NULL,
  duration_scale = c("rate", "integral"),
  weighting = c("precision", "gls", "none"),
  tau2_mode = c("auto", "residual", "marginal"),
  tau2_shrinkage = c("auto", "graph_eb", "eb", "none"),
  spatial_graph = NULL,
  aggregate_by_parcel = FALSE,
  parcel_groups = NULL,
  whiten_panel = FALSE,
  return_epochs = FALSE,
  return_panel = FALSE,
  run_col = "run",
  variance_confounds = NULL,
  shared_error_method = c("poet_covreg", "legacy", "none"),
  bold = NULL,
  other_event_model = c("auto", "trial", "condition")
) {
  fit <- deconvolution$hrf_fit
  settings <- deconvolution$settings
  if (!inherits(fit, "crossrun_hrf") || is.null(deconvolution$scaling) ||
      is.null(settings$nuisance_mode)) {
    stop("Use a current deconvolve_crossrun_hrf result")
  }
  if (settings$backend == "wiener") {
    stop("decon_epoch_regression requires the cpp or r GLS inverse")
  }
  source_model_supplied <- !missing(source_model)
  within_ridge_supplied <- !missing(within_trial_ridge)
  demix_ridge_supplied <- !missing(demix_ridge)
  return_panel_supplied <- !missing(return_panel)
  other_model_supplied <- !missing(other_event_model)
  readout_kernel <- match.arg(readout_kernel)
  source_model <- match.arg(source_model)
  duration_scale <- match.arg(duration_scale)
  weighting <- match.arg(weighting)
  tau2_mode <- match.arg(tau2_mode)
  tau2_shrinkage <- match.arg(tau2_shrinkage)
  shared_error_method <- match.arg(shared_error_method)
  other_event_model <- match.arg(other_event_model)
  if (!is.null(trial_durations) && !source_model_supplied) {
    source_model <- "duration"
  }
  if (!is.null(parcel_groups)) {
    aggregate_by_parcel <- TRUE
  }
  checkmate::assert_numeric(lags, min.len = 1L, any.missing = FALSE, finite = TRUE)
  checkmate::assert_string(target_condition, null.ok = TRUE)
  checkmate::assert_flag(demix_other_events)
  checkmate::assert_flag(demix_target_trials)
  checkmate::assert_flag(demix_within_trial)
  checkmate::assert_number(demix_ridge, lower = 1e-6, finite = TRUE)
  if (!other_model_supplied && demix_ridge_supplied) {
    other_event_model <- if (demix_ridge < 0.03) "trial" else "condition"
  }
  checkmate::assert_number(within_trial_ridge, lower = 1e-6, finite = TRUE)
  checkmate::assert_number(gain_floor, lower = 1e-4, finite = TRUE)
  checkmate::assert_flag(restore_scale)
  checkmate::assert_flag(aggregate_by_parcel)
  checkmate::assert_flag(whiten_panel)
  checkmate::assert_flag(return_epochs)
  checkmate::assert_flag(return_panel)
  if (is.null(covariates) && !return_panel_supplied && !isTRUE(return_epochs)) {
    return_panel <- TRUE
  }

  if (is.null(events)) events <- fit$specification$events
  if (is.data.frame(events)) {
    events <- hrf_resolve_events_input(
      events,
      stats::setNames(vector("list", length(fit$run_ids)), fit$run_ids),
      run_col = run_col,
      allow_empty_runs = FALSE
    )
  }
  if (!is.list(events) || !identical(names(events), fit$run_ids)) {
    stop("events must be a list in fitted run order")
  }
  if (!is.null(variance_confounds)) {
    variance_confounds <- .resolve_variance_confounds_input(
      variance_confounds,
      run_ids = fit$run_ids,
      run_col = run_col
    )
  }

  n_runs <- length(fit$run_ids)
  target_idx_list <- vector("list", n_runs)
  other_idx_list <- vector("list", n_runs)
  for (r in seq_len(n_runs)) {
    ev <- events[[r]]
    checkmate::assert_data_frame(ev, min.rows = 1L)
    checkmate::assert_numeric(ev$onset, any.missing = FALSE, finite = TRUE)
    if (is.null(target_condition)) {
      target_idx_list[[r]] <- seq_len(nrow(ev))
      other_idx_list[[r]] <- integer(0L)
    } else {
      if (is.null(ev$condition)) {
        stop("target_condition supplied, but events[[", r, "]] has no 'condition' column")
      }
      cond <- as.character(ev$condition)
      target_idx_list[[r]] <- which(!is.na(cond) & cond == target_condition)
      other_idx_list[[r]] <- which(is.na(cond) | cond != target_condition)
      if (!length(target_idx_list[[r]])) {
        stop("No events match target_condition = '", target_condition, "' in run ", fit$run_ids[r])
      }
    }
  }

  target_counts <- vapply(target_idx_list, length, integer(1L))
  ends <- cumsum(target_counts)
  starts <- ends - target_counts + 1L
  reserved_panel_cols <- c(
    "run", "trial", "event_index", "onset", "evt_time", "voxel", "parcel",
    "decon_raw", "decon_demixed", "decon_eb", "se_raw", "se_demixed",
    "reliability_trial", "precision_weight", "precision_weight_cohort",
    "precision_weight_op", "precision_weight_raw", "error_ac1", "tau2_bio",
    "shared_noise_loading", "shared_noise_vec", "trial_cor_vec", "lag_noise_vec",
    "decon_demixed_whitened", "se_whitened", "decon_eb_whitened"
  )
  cov_list <- vector("list", n_runs)
  if (is.null(covariates)) {
    for (r in seq_len(n_runs)) {
      cov_list[[r]] <- matrix(numeric(0L), nrow = target_counts[r], ncol = 0L)
    }
    term_names <- character(0L)
    n_terms <- 0L
    eff_tau2_mode <- "marginal"
  } else if (is.character(covariates) && is.null(dim(covariates))) {
    checkmate::assert_character(covariates, min.len = 1L, any.missing = FALSE, min.chars = 1L, unique = TRUE)
    bad_res <- intersect(covariates, reserved_panel_cols)
    if (length(bad_res)) {
      stop("Covariate column name(s) collide with reserved panel output columns: ",
           paste(bad_res, collapse = ", "), call. = FALSE)
    }
    for (r in seq_len(n_runs)) {
      ev_r <- events[[r]]
      miss_c <- setdiff(covariates, names(ev_r))
      if (length(miss_c)) {
        stop(
          "Covariate column(s) not found in events[['", fit$run_ids[r], "']]: ",
          paste(miss_c, collapse = ", "),
          call. = FALSE
        )
      }
      sub_df <- ev_r[target_idx_list[[r]], covariates, drop = FALSE]
      num_cols <- vector("list", length(covariates))
      for (ci in seq_along(covariates)) {
        col_v <- sub_df[[ci]]
        if (!is.numeric(col_v) && !is.logical(col_v)) {
          stop(
            "Covariate column '", covariates[ci], "' in events[['", fit$run_ids[r],
            "']] must be numeric or logical (got ", class(col_v)[1L], ").",
            call. = FALSE
          )
        }
        num_cols[[ci]] <- as.double(col_v)
      }
      x_r <- do.call(cbind, num_cols)
      colnames(x_r) <- covariates
      cov_list[[r]] <- x_r
    }
    term_names <- covariates
    n_terms <- length(term_names)
    eff_tau2_mode <- tau2_mode
  } else if (inherits(covariates, "formula")) {
    if (length(covariates) != 2L) {
      stop("When covariates is a formula, it must be one-sided (e.g., ~ rpe + abs_rpe).", call. = FALSE)
    }
    f_vars <- all.vars(covariates)
    for (r in seq_len(n_runs)) {
      miss_v <- setdiff(f_vars, names(events[[r]]))
      if (length(miss_v)) {
        stop(
          "Formula variable(s) not found in events[['", fit$run_ids[r], "']]: ",
          paste(miss_v, collapse = ", "),
          call. = FALSE
        )
      }
    }
    target_ev_all <- do.call(
      rbind,
      lapply(seq_len(n_runs), function(r) events[[r]][target_idx_list[[r]], f_vars, drop = FALSE])
    )
    mf_all <- stats::model.frame(covariates, data = target_ev_all, na.action = stats::na.pass)
    mm_all <- stats::model.matrix(covariates, data = mf_all)
    int_col <- match("(Intercept)", colnames(mm_all), nomatch = 0L)
    if (int_col > 0L) mm_all <- mm_all[, -int_col, drop = FALSE]
    if (ncol(mm_all) == 0L) {
      stop("Formula 'covariates' produced zero predictor columns after removing intercept.", call. = FALSE)
    }
    bad_res <- intersect(colnames(mm_all), reserved_panel_cols)
    if (length(bad_res)) {
      stop("Formula term name(s) collide with reserved panel output columns: ",
           paste(bad_res, collapse = ", "), call. = FALSE)
    }
    for (r in seq_len(n_runs)) {
      cov_list[[r]] <- mm_all[starts[r]:ends[r], , drop = FALSE]
    }
    term_names <- colnames(mm_all)
    n_terms <- length(term_names)
    eff_tau2_mode <- tau2_mode
  } else if (is.list(covariates) && !is.data.frame(covariates)) {
    if (length(covariates) != n_runs) {
      stop("When covariates is a list, it must have one element per run")
    }
    if (!is.null(names(covariates)) && any(nzchar(names(covariates)))) {
      if (!setequal(names(covariates), fit$run_ids)) {
        stop("Named covariates list must match fitted run IDs: ",
             paste(fit$run_ids, collapse = ", "), call. = FALSE)
      }
      covariates <- covariates[fit$run_ids]
    }
    for (r in seq_len(n_runs)) {
      if (is.data.frame(covariates[[r]])) {
        bad_cols <- names(covariates[[r]])[
          !vapply(covariates[[r]], function(col) is.numeric(col) || is.logical(col), logical(1L))
        ]
        if (length(bad_cols)) {
          stop("Non-numeric column(s) in covariates[[", r, "]]: ",
               paste(bad_cols, collapse = ", "), call. = FALSE)
        }
      }
      x_r <- as.matrix(covariates[[r]])
      storage.mode(x_r) <- "double"
      if (nrow(x_r) != target_counts[r]) {
        stop("covariates[[", r, "]] row count (", nrow(x_r),
          ") does not match target event count (", target_counts[r], ")")
      }
      if (is.null(colnames(x_r))) colnames(x_r) <- paste0("x", seq_len(ncol(x_r)))
      cov_list[[r]] <- x_r
    }
    term_names <- colnames(cov_list[[1L]])
    n_terms <- length(term_names)
    eff_tau2_mode <- tau2_mode
  } else {
    if (is.data.frame(covariates)) {
      bad_cols <- names(covariates)[
        !vapply(covariates, function(col) is.numeric(col) || is.logical(col), logical(1L))
      ]
      if (length(bad_cols)) {
        stop(
          "Non-numeric column(s) in covariates data.frame (",
          paste(bad_cols, collapse = ", "),
          "); select numeric predictor columns or pass covariates = c(...) column names from events.",
          call. = FALSE
        )
      }
    }
    x_all <- as.matrix(covariates)
    storage.mode(x_all) <- "double"
    if (nrow(x_all) != sum(target_counts)) {
      stop("Total rows in covariates (", nrow(x_all),
        ") must equal total target events across runs (", sum(target_counts), ")")
    }
    if (is.null(colnames(x_all))) colnames(x_all) <- paste0("x", seq_len(ncol(x_all)))
    for (r in seq_len(n_runs)) {
      cov_list[[r]] <- x_all[starts[r]:ends[r], , drop = FALSE]
    }
    term_names <- colnames(cov_list[[1L]])
    n_terms <- length(term_names)
    eff_tau2_mode <- tau2_mode
  }
  for (r in seq_len(n_runs)) {
    rownames(cov_list[[r]]) <- NULL
  }

  n_voxels <- length(fit$voxel_ids)
  n_lags <- length(lags)
  tr <- fit$settings$tr

  # Resolve spatial parcel grouping and precompute cross-run inverse-variance voxel weights
  parcel_vox_list <- NULL
  parcel_weights_list <- NULL
  if (isTRUE(aggregate_by_parcel)) {
    if (!is.null(parcel_groups)) {
      if (!is.null(names(parcel_groups)) && all(fit$voxel_ids %in% names(parcel_groups))) {
        p_groups <- as.character(unname(parcel_groups[fit$voxel_ids]))
      } else {
        if (length(parcel_groups) != n_voxels) {
          stop("parcel_groups length (", length(parcel_groups),
            ") must equal number of voxels (", n_voxels, ")")
        }
        p_groups <- as.character(parcel_groups)
      }
    } else if (!is.null(fit$specification$spatial_groups)) {
      p_groups <- as.character(fit$specification$spatial_groups)
    } else {
      p_groups <- rep("parcel_1", n_voxels)
    }
    if (anyNA(p_groups) || any(!nzchar(p_groups))) {
      stop("parcel_groups must not contain NA or empty labels")
    }
    unit_ids <- unique(p_groups)
    n_units <- length(unit_ids)
    parcel_vox_list <- lapply(unit_ids, function( pid ) which(p_groups == pid))
    names(parcel_vox_list) <- unit_ids

    abs_s2_mat <- matrix(NA_real_, n_runs, n_voxels)
    rel_s2_mat <- matrix(NA_real_, n_runs, n_voxels)
    cv_r2_mat <- matrix(0.0, n_runs, n_voxels)
    for (r in seq_len(n_runs)) {
      sc_df_r <- deconvolution$scaling[[r]]
      ok_r <- sc_df_r$status == "ok"
      div_r0 <- sc_df_r$scale
      div_r0[!is.finite(div_r0) | div_r0 <= 0] <- 1
      nsd_r0 <- if (!is.null(sc_df_r$noise_sd)) sc_df_r$noise_sd else div_r0
      need_n0 <- ok_r & (!is.finite(nsd_r0) | nsd_r0 <= 0)
      if (any(need_n0)) nsd_r0[need_n0] <- div_r0[need_n0]
      rel_r0 <- nsd_r0 / div_r0
      rel_r0[!is.finite(rel_r0) | rel_r0 <= 0] <- 1
      abs_s2_mat[r, ok_r] <- nsd_r0[ok_r]^2
      rel_s2_mat[r, ok_r] <- rel_r0[ok_r]^2
      sel_r <- fit$selections[[r]]
      sc_raw <- if (!is.null(sel_r[["selected_adjusted_partial_r2"]])) {
        sel_r[["selected_adjusted_partial_r2"]]
      } else if (!is.null(sel_r[["best_adjusted_partial_r2"]])) {
        sel_r[["best_adjusted_partial_r2"]]
      } else {
        sel_r[["score"]]
      }
      if (!is.null(sc_raw) && is.numeric(sc_raw)) {
        sc_v <- as.numeric(sc_raw)
        sc_v[!is.finite(sc_v) | sc_v < 0] <- 0
        cv_r2_mat[r, ] <- sc_v
      }
    }
    parcel_weights_list <- vector("list", n_units)
    names(parcel_weights_list) <- unit_ids
    for (g in seq_len(n_units)) {
      vg <- parcel_vox_list[[g]]
      mean_abs_s2 <- colMeans(abs_s2_mat[, vg, drop = FALSE], na.rm = TRUE)
      mean_rel_s2 <- colMeans(rel_s2_mat[, vg, drop = FALSE], na.rm = TRUE)
      mean_cv_r2 <- colMeans(cv_r2_mat[, vg, drop = FALSE], na.rm = TRUE)
      mean_abs_s2[!is.finite(mean_abs_s2) | mean_abs_s2 <= 0] <- 1.0
      mean_rel_s2[!is.finite(mean_rel_s2) | mean_rel_s2 <= 0] <- 1.0
      med_abs <- stats::median(mean_abs_s2, na.rm = TRUE)
      med_rel <- stats::median(mean_rel_s2, na.rm = TRUE)
      s2_abs_fl <- pmax(mean_abs_s2, 0.10 * max(med_abs, 1e-8))
      s2_rel_fl <- pmax(mean_rel_s2, 0.10 * max(med_rel, 1e-8))
      med_r2 <- stats::median(mean_cv_r2, na.rm = TRUE)
      ev_mult <- if (any(mean_cv_r2 > 0)) {
        1.0 + mean_cv_r2 / (0.05 + max(med_r2, 1e-4))
      } else {
        rep(1.0, length(vg))
      }
      w_g <- if (isTRUE(restore_scale) || settings$scaling %in% c("demean", "none")) {
        ev_mult / s2_abs_fl
      } else {
        ev_mult / s2_rel_fl
      }
      w_g[!is.finite(w_g) | w_g <= 0] <- 1.0
      parcel_weights_list[[g]] <- w_g / sum(w_g)
    }
  } else {
    unit_ids <- fit$voxel_ids
    n_units <- n_voxels
  }

  dur_target_list <- vector("list", n_runs)
  if (!is.null(trial_durations)) {
    if (is.character(trial_durations) && is.null(dim(trial_durations))) {
      checkmate::assert_string(trial_durations, min.chars = 1L)
      for (r in seq_len(n_runs)) {
        if (!(trial_durations %in% names(events[[r]]))) {
          stop("trial_durations column '", trial_durations, "' not found in events[['",
               fit$run_ids[r], "']].", call. = FALSE)
        }
        d_r <- as.numeric(events[[r]][[trial_durations]][target_idx_list[[r]]])
        checkmate::assert_numeric(d_r, lower = 0, any.missing = FALSE, finite = TRUE)
        dur_target_list[[r]] <- d_r
      }
    } else if (is.list(trial_durations) && !is.data.frame(trial_durations)) {
      if (length(trial_durations) != n_runs) {
        stop("When trial_durations is a list, it must have one element per run")
      }
      if (!is.null(names(trial_durations)) && any(nzchar(names(trial_durations)))) {
        if (!setequal(names(trial_durations), fit$run_ids)) {
          stop("Named trial_durations list must match fitted run IDs: ",
               paste(fit$run_ids, collapse = ", "), call. = FALSE)
        }
        trial_durations <- trial_durations[fit$run_ids]
      }
      for (r in seq_len(n_runs)) {
        d_r <- as.numeric(trial_durations[[r]])
        if (length(d_r) != target_counts[r]) {
          stop("trial_durations[[", r, "]] length (", length(d_r),
            ") does not match target event count (", target_counts[r], ")")
        }
        checkmate::assert_numeric(d_r, lower = 0, any.missing = FALSE, finite = TRUE)
        dur_target_list[[r]] <- d_r
      }
    } else {
      d_all <- as.numeric(trial_durations)
      if (length(d_all) != sum(target_counts)) {
        stop("Total length of trial_durations (", length(d_all),
          ") must equal total target events across runs (", sum(target_counts), ")")
      }
      checkmate::assert_numeric(d_all, lower = 0, any.missing = FALSE, finite = TRUE)
      for (r in seq_len(n_runs)) {
        dur_target_list[[r]] <- d_all[starts[r]:ends[r]]
      }
    }
  }

  total_trials <- sum(target_counts)
  raw_epochs <- demixed_epochs <- se_raw_epochs <- se_demixed_epochs <- vector("list", n_lags)
  pw_epochs <- pw_cohort_epochs <- pw_op_epochs <- pw_raw_epochs <- error_ac1_epochs <- tau2_bio_epochs <- vector("list", n_lags)
  rel_epochs <- eb_y_epochs <- vector("list", n_lags)
  whitened_y_epochs <- se_whitened_epochs <- eb_whitened_y_epochs <- if (isTRUE(whiten_panel)) vector("list", n_lags) else NULL
  whitened_x_epochs <- if (isTRUE(whiten_panel) && n_terms > 0L) vector("list", n_lags) else NULL
  V_demix_blocks <- vector("list", n_lags)
  lag_names <- sprintf("lag_%g", lags)
  names(raw_epochs) <- names(demixed_epochs) <- names(se_raw_epochs) <- names(se_demixed_epochs) <- lag_names
  names(pw_epochs) <- names(pw_cohort_epochs) <- names(pw_op_epochs) <- names(pw_raw_epochs) <- names(error_ac1_epochs) <- names(tau2_bio_epochs) <- lag_names
  names(rel_epochs) <- names(eb_y_epochs) <- lag_names
  if (isTRUE(whiten_panel)) {
    names(whitened_y_epochs) <- names(se_whitened_epochs) <- names(eb_whitened_y_epochs) <- lag_names
  }
  names(V_demix_blocks) <- lag_names
  for (li in seq_len(n_lags)) {
    na_mat <- matrix(NA_real_, total_trials, n_units, dimnames = list(NULL, unit_ids))
    raw_epochs[[li]] <- na_mat
    demixed_epochs[[li]] <- na_mat
    se_raw_epochs[[li]] <- na_mat
    se_demixed_epochs[[li]] <- na_mat
    pw_epochs[[li]] <- na_mat
    pw_cohort_epochs[[li]] <- na_mat
    pw_op_epochs[[li]] <- na_mat
    pw_raw_epochs[[li]] <- na_mat
    error_ac1_epochs[[li]] <- na_mat
    tau2_bio_epochs[[li]] <- na_mat
    rel_epochs[[li]] <- na_mat
    eb_y_epochs[[li]] <- na_mat
    if (isTRUE(whiten_panel)) {
      whitened_y_epochs[[li]] <- na_mat
      se_whitened_epochs[[li]] <- na_mat
      eb_whitened_y_epochs[[li]] <- na_mat
      if (n_terms > 0L) {
        whitened_x_epochs[[li]] <- stats::setNames(
          replicate(n_terms, na_mat, simplify = FALSE),
          term_names
        )
      }
    }
    V_demix_blocks[[li]] <- stats::setNames(vector("list", n_runs), fit$run_ids)
    for (r in seq_len(n_runs)) V_demix_blocks[[li]][[r]] <- list()
  }
  gain_by_lag_voxel <- matrix(NA_real_, n_lags, n_units, dimnames = list(NULL, unit_ids))
  x_centered_list <- if (n_terms > 0L) {
    lapply(cov_list, function(x_r) sweep(x_r, 2L, colMeans(x_r, na.rm = TRUE), "-"))
  } else {
    cov_list
  }
  X_raw_pool <- do.call(rbind, cov_list)

  gain_accum <- matrix(0, n_lags, n_units)
  gain_count <- matrix(0L, n_lags, n_units)
  stratum_map <- vector("list", n_runs)
  sig_eff_mat <- matrix(NA_real_, n_runs, n_units, dimnames = list(fit$run_ids, unit_ids))
  center_obs <- !identical(settings$scaling, "none")
  lag_cov_accum <- stats::setNames(
    replicate(n_units, matrix(0, n_lags, n_lags), simplify = FALSE),
    unit_ids
  )
  lag_cov_trials <- integer(n_units)
  lag_res_accum <- stats::setNames(
    replicate(n_units, matrix(0, n_lags, n_lags), simplify = FALSE),
    unit_ids
  )
  lag_res_trials <- integer(n_units)
  track_unit_cor <- identical(shared_error_method, "legacy") && (n_units >= 2L && n_units <= 256L)
  unit_err_cor_accum <- if (track_unit_cor) {
    matrix(0.0, n_units, n_units, dimnames = list(unit_ids, unit_ids))
  } else NULL
  unit_err_df_tot <- 0.0
  li_ref_cor <- which.min(abs(lags))
  run_stratum_ops <- vector("list", n_runs)
  global_op_cache <- new.env(parent = emptyenv())
  global_tw_cache <- new.env(parent = emptyenv())
  other_model_eff <- if (identical(other_event_model, "auto") &&
                         (any(lags < -1e-6) || identical(source_model, "window") ||
                          (isTRUE(return_panel) && isFALSE(whiten_panel)))) {
    "condition"
  } else {
    other_event_model
  }
  need_crossrun_haus <- (n_runs >= 2L) && identical(other_model_eff, "auto") &&
    isTRUE(demix_other_events) && (n_terms > 0L) && any(lengths(other_idx_list) > 0L)
  gamma_haus_override <- NULL

  for (haus_pass in if (need_crossrun_haus) 1:2 else 1L) {
    if (haus_pass == 2L) {
      gain_accum[] <- 0
      gain_count[] <- 0L
      stratum_map <- vector("list", n_runs)
      sig_eff_mat[] <- NA_real_
      for (u_i in seq_len(n_units)) {
        lag_cov_accum[[u_i]][] <- 0
        lag_res_accum[[u_i]][] <- 0
      }
      lag_cov_trials[] <- 0L
      lag_res_trials[] <- 0L
      if (!is.null(unit_err_cor_accum)) unit_err_cor_accum[] <- 0.0
      unit_err_df_tot <- 0.0
      run_stratum_ops <- vector("list", n_runs)
      for (li_r in seq_len(n_lags)) {
        for (rr_i in seq_len(n_runs)) V_demix_blocks[[li_r]][[rr_i]] <- list()
      }
    }
    haus_W <- if (need_crossrun_haus && haus_pass == 1L) {
      replicate(n_units, matrix(0.0, n_terms, n_terms), simplify = FALSE)
    } else NULL
    haus_rhs <- if (need_crossrun_haus && haus_pass == 1L) {
      replicate(n_units, numeric(n_terms), simplify = FALSE)
    } else NULL
    haus_df2 <- if (need_crossrun_haus && haus_pass == 1L) numeric(n_units) else NULL
    haus_gam_p1 <- if (need_crossrun_haus && haus_pass == 1L) matrix(NA_real_, n_runs, n_units) else NULL

  for (r in seq_len(n_runs)) {
    n <- fit$n_time[r]
    ev <- events[[r]]
    offset <- fit$settings$time_offset[r]
    t_idx <- target_idx_list[[r]]
    o_idx <- other_idx_list[[r]]
    m_all <- nrow(ev)
    m_target <- length(t_idx)
    row_span <- starts[r]:ends[r]
    duration <- if (is.null(ev$duration)) rep(0, m_all) else as.numeric(ev$duration)
    checkmate::assert_numeric(duration, lower = 0, len = m_all, any.missing = FALSE, finite = TRUE)
    if (!is.null(trial_durations)) {
      duration[t_idx] <- dur_target_list[[r]]
    }
    dur_target <- duration[t_idx]
    w_src <- if (is.null(readout_width)) 1.5 * tr else readout_width[1L]

    mask <- if (settings$observation_mask == "fit") fit$specification$valid_rows[[r]] else rep(TRUE, n)
    rows <- which(!is.na(mask) & mask)
    if (length(rows) < 3L) next
    hp_ep <- if (!is.null(settings$highpass)) settings$highpass else if (!is.null(fit$settings$highpass)) fit$settings$highpass else 128.0
    z <- hrf_build_decon_nuisance_matrix(
      n_obs = n, tr = tr, rows = rows,
      nuisance_r = fit$specification$nuisance[[r]],
      nuisance_mode = settings$nuisance_mode, highpass_sec = hp_ep,
      drift_active = !isFALSE(settings$drift_preproject_active)
    )

    strat <- hrf_stratum_kernels(fit, r)
    rho <- strat$rho
    sc_df <- deconvolution$scaling[[r]]
    valid <- sc_df$status == "ok"
    act_r <- deconvolution$activity[[r]]
    lam_default_r <- as.numeric(settings$lambda)[1L]
    if (!is.finite(lam_default_r) || lam_default_r < 0) lam_default_r <- 0.1
    lam_vec_r <- if (!is.null(sc_df$lambda) && length(sc_df$lambda) == ncol(act_r)) {
      lv <- as.numeric(sc_df$lambda)
      lv[!is.finite(lv) | lv < 0] <- lam_default_r
      lv
    } else {
      rep(lam_default_r, ncol(act_r))
    }
    keys <- if (any(valid) && length(unique(lam_vec_r[valid])) > 1L) {
      paste(strat$keys, sprintf("%.12g", lam_vec_r), sep = "|lam:")
    } else {
      strat$keys
    }
    div_r <- sc_df$scale
    div_r[!is.finite(div_r) | div_r <= 0] <- 1
    noise_sd_r <- if (!is.null(sc_df$noise_sd)) sc_df$noise_sd else rep(NA_real_, ncol(act_r))
    need_nsd <- valid & (!is.finite(noise_sd_r) | noise_sd_r <= 0)
    if (any(need_nsd)) {
      noise_sd_r[need_nsd] <- div_r[need_nsd]
    }
    # Effective observation-noise SD in the working analysis units:
    # - If restore_scale = TRUE, deconvolved activity is multiplied back by div_r,
    #   so observation noise in restored units has SD = noise_sd_r.
    # - Otherwise, deconvolved activity is in standardized units (divided by div_r),
    #   so observation noise has SD = noise_sd_r / div_r (equals 1 under scaling = "noise_sd").
    sig_eff_r <- if (isTRUE(restore_scale)) noise_sd_r else (noise_sd_r / div_r)
    sig_eff_r[!is.finite(sig_eff_r) | sig_eff_r <= 0] <- 1

    if (isTRUE(restore_scale)) {
      act_r[, valid] <- sweep(act_r[, valid, drop = FALSE], 2L, div_r[valid], "*")
    }

    fw_r <- resolve_run_obs_weights(deconvolution, r, length(rows))
    w_space <- if (!is.null(settings$rwls_space)) settings$rwls_space else "woodbury"
    Z_vox_r <- if (!is.null(deconvolution$noise_innovations) &&
                   length(deconvolution$noise_innovations) >= r) {
      deconvolution$noise_innovations[[r]]
    } else NULL
    has_z_vox <- is.matrix(Z_vox_r) && nrow(Z_vox_r) == length(rows) && ncol(Z_vox_r) == n_voxels
    if (has_z_vox && !is.null(fw_r) && !isTRUE(attr(Z_vox_r, "rwls_whitened"))) {
      Z_vox_r <- sweep(Z_vox_r, 1L, sqrt(pmax(fw_r, 1e-8)), "*")
      zn_col <- sqrt(colSums(Z_vox_r^2))
      pos_zn <- is.finite(zn_col) & zn_col > 1e-12
      if (any(pos_zn)) {
        Z_vox_r[, pos_zn] <- sweep(Z_vox_r[, pos_zn, drop = FALSE], 2L, zn_col[pos_zn], "/")
      }
    }
    mix_by_key <- list()

    if (isTRUE(aggregate_by_parcel)) {
      act_parcel_r <- matrix(NA_real_, n, n_units, dimnames = list(NULL, unit_ids))
      valid_parcel <- logical(n_units)
      sig_eff_parcel_r <- rep(1.0, n_units)
      rho_parcel <- rep(0.0, n_units)
      lam_parcel <- rep(lam_default_r, n_units)
      keys_parcel <- character(n_units)
      kernels_parcel <- vector("list", n_units)
      ar_parcel <- vector("list", n_units)
      sketch_vox_r <- if (!is.null(deconvolution$whitened_sketch) &&
                          length(deconvolution$whitened_sketch) >= r) {
        deconvolution$whitened_sketch[[r]]
      } else NULL
      has_sketch_vox <- is.list(sketch_vox_r) && is.matrix(sketch_vox_r$U_sk) &&
        is.matrix(sketch_vox_r$B_sk) && nrow(sketch_vox_r$U_sk) == length(rows) &&
        ncol(sketch_vox_r$B_sk) == n_voxels

      z_unit_parcel_r <- if ((has_z_vox || has_sketch_vox) && track_unit_cor) {
        matrix(0.0, length(rows), n_units, dimnames = list(NULL, unit_ids))
      } else NULL

      # Precompute low-rank neural trial projection fallback when Z_vox_r & whitened_sketch are absent
      Q_task_fb <- NULL
      if (!has_z_vox && !has_sketch_vox && length(rows) >= 12L) {
        G_fb_act <- hrf_event_design_cpp(
          n, tr, offset, ev$onset, duration, seq_len(m_all), m_all, matrix(c(1.0, 0.0), 2L, 1L)
        )[[1L]][rows, , drop = FALSE]
        G_fb_all <- cbind(1.0, seq(-1.0, 1.0, length.out = length(rows)), G_fb_act)
        nrm_fb <- sqrt(colSums(G_fb_all^2)); nrm_fb[nrm_fb == 0] <- 1.0
        sv_fb <- svd(sweep(G_fb_all, 2L, nrm_fb, "/"), nv = 0)
        rk_max_fb <- max(0L, length(rows) - 10L)
        rk_fb <- if (length(sv_fb$d) && sv_fb$d[1L] > 0 && rk_max_fb > 0L) {
          min(rk_max_fb, sum(sv_fb$d > 0.05 * sv_fb$d[1L]))
        } else 0L
        if (rk_fb > 0L) Q_task_fb <- sv_fb$u[, seq_len(rk_fb), drop = FALSE]
      }

      for (g in seq_len(n_units)) {
        vox_g <- parcel_vox_list[[g]]
        ok_g <- valid[vox_g]
        if (!any(ok_g)) next
        v_ok <- vox_g[ok_g]
        valid_parcel[g] <- TRUE
        u_g <- parcel_weights_list[[g]][ok_g]
        u_g <- u_g / sum(u_g)
        act_parcel_r[, g] <- as.numeric(act_r[, v_ok, drop = FALSE] %*% u_g)

        seff_v <- sig_eff_r[v_ok]
        s2_indep <- sum((u_g^2) * (seff_v^2))
        s2_ceil <- (sum(u_g * seff_v))^2
        if (length(v_ok) == 1L || length(rows) < 4L) {
          sig_eff_parcel_r[g] <- sqrt(s2_indep)
          if (!is.null(z_unit_parcel_r)) {
            if (has_z_vox) {
              z_unit_parcel_r[, g] <- Z_vox_r[, v_ok[1L]]
            } else if (has_sketch_vox) {
              zg1 <- as.numeric(sketch_vox_r$U_sk %*% sketch_vox_r$B_sk[, v_ok[1L]])
              ng1 <- sqrt(sum(zg1^2))
              if (ng1 > 1e-12) z_unit_parcel_r[, g] <- zg1 / ng1
            }
          }
        } else if (has_z_vox) {
          zeta_g <- as.numeric(Z_vox_r[, v_ok, drop = FALSE] %*% (u_g * seff_v))
          z_norm2 <- sum(zeta_g^2)
          if (is.finite(z_norm2) && z_norm2 > 1e-20) {
            s2_eff_g <- min(s2_ceil, max(0.05 * s2_indep, z_norm2))
            if (!is.null(z_unit_parcel_r)) {
              z_unit_parcel_r[, g] <- zeta_g / sqrt(z_norm2)
            }
          } else {
            s2_eff_g <- s2_indep
          }
          sig_eff_parcel_r[g] <- sqrt(max(s2_eff_g, 1e-12))
        } else if (has_sketch_vox) {
          wt_sk <- u_g * (if (isTRUE(restore_scale)) 1.0 else (1.0 / div_r[v_ok]))
          b_sk_g <- as.numeric(sketch_vox_r$B_sk[, v_ok, drop = FALSE] %*% wt_sk)
          df_v_ok <- pmax(1.0, ifelse(is.finite(sketch_vox_r$df_vox[v_ok]), sketch_vox_r$df_vox[v_ok], length(rows) - 4.0))
          s2_cap_v <- colSums(sketch_vox_r$B_sk[, v_ok, drop = FALSE]^2) / df_v_ok
          s2_unexp <- pmax(0.0, ifelse(is.finite(sketch_vox_r$s2_raw[v_ok]), sketch_vox_r$s2_raw[v_ok], 0.0) - s2_cap_v)
          s2_sk_est <- sum(b_sk_g^2) / mean(df_v_ok) + sum((wt_sk^2) * s2_unexp)
          s2_eff_g <- min(s2_ceil, max(0.05 * s2_indep, s2_sk_est))
          if (!is.null(z_unit_parcel_r)) {
            zeta_g <- as.numeric(sketch_vox_r$U_sk %*% b_sk_g)
            ng_sk <- sqrt(sum(zeta_g^2))
            if (is.finite(ng_sk) && ng_sk > 1e-12) z_unit_parcel_r[, g] <- zeta_g / ng_sk
          }
          sig_eff_parcel_r[g] <- sqrt(max(s2_eff_g, 1e-12))
        } else {
          A_sub <- act_r[rows, v_ok, drop = FALSE]
          if (!is.null(Q_task_fb)) {
            A_sub <- A_sub - Q_task_fb %*% crossprod(Q_task_fb, A_sub)
          }
          A_sub <- sweep(A_sub, 2L, colMeans(A_sub), "-")
          v_sub <- sqrt(pmax(colSums(A_sub^2), 1e-24))
          U_sub <- sweep(A_sub, 2L, v_sub, "/")
          zeta_fb <- as.numeric(U_sub %*% (u_g * seff_v))
          s2_eff_g <- min(s2_ceil, max(0.05 * s2_indep, sum(zeta_fb^2)))
          sig_eff_parcel_r[g] <- sqrt(max(s2_eff_g, 1e-12))
        }

        u_keys <- unique(keys[v_ok])
        if (length(u_keys) == 1L) {
          keys_parcel[g] <- u_keys[1L]
          rho_parcel[g] <- rho[v_ok[1L]]
          lam_parcel[g] <- lam_vec_r[v_ok[1L]]
          kernels_parcel[[g]] <- strat$get_kernel(v_ok[1L])
          ar_parcel[[g]] <- if (!is.null(strat$get_ar)) strat$get_ar(v_ok[1L]) else rho[v_ok[1L]]
        } else {
          rho_g <- sum(u_g * rho[v_ok])
          lam_parcel[g] <- exp(sum(u_g * log(pmax(lam_vec_r[v_ok], 1e-12))))
          psi1_g <- if (!is.null(strat$psi1)) sum(u_g * strat$psi1[v_ok]) else 0.0
          psi2_g <- if (!is.null(strat$psi2)) sum(u_g * strat$psi2[v_ok]) else 0.0
          psi3_g <- if (!is.null(strat$psi3)) sum(u_g * strat$psi3[v_ok]) else 0.0
          psi4_g <- if (!is.null(strat$psi4)) sum(u_g * strat$psi4[v_ok]) else 0.0
          k_mat <- vapply(v_ok, function(vi) as.numeric(strat$get_kernel(vi)), numeric(nrow(fit$kernels)))
          k_g <- as.numeric(k_mat %*% u_g)
          norm_mode_g <- if (!is.null(fit$settings$kernel_normalize)) fit$settings$kernel_normalize else "signed_sum"
          k_g <- tryCatch(normalize_decon_kernel(k_g, mode = norm_mode_g), error = function(e) k_g)
          rho_parcel[g] <- rho_g
          kernels_parcel[[g]] <- matrix(k_g, ncol = 1L)
          pkey_g <- paste0("parcel_mix_", unit_ids[g], "_", r)
          keys_parcel[g] <- pkey_g
          has_p34_g <- abs(psi3_g) > 1e-10 || abs(psi4_g) > 1e-10
          has_p12_g <- abs(psi1_g) > 1e-10 || abs(psi2_g) > 1e-10 || has_p34_g
          ar_parcel[[g]] <- if (has_p34_g) {
            c(rho_g, psi1_g, psi2_g, psi3_g, psi4_g)
          } else if (has_p12_g) {
            c(rho_g, psi1_g, psi2_g)
          } else {
            rho_g
          }
          if (length(u_keys) <= 32L) {
            mix_by_key[[pkey_g]] <- lapply(u_keys, function(sk) {
              idx_sk <- which(keys[v_ok] == sk)
              v_rep <- v_ok[idx_sk[1L]]
              list(
                key = sk,
                weight = sum(u_g[idx_sk]),
                kernel = strat$get_kernel(v_rep),
                ar = if (!is.null(strat$get_ar)) strat$get_ar(v_rep) else rho[v_rep],
                lambda = lam_vec_r[v_rep]
              )
            })
          }
        }
      }
      act_r <- act_parcel_r
      valid <- valid_parcel
      sig_eff_r <- sig_eff_parcel_r
      rho <- rho_parcel
      keys <- keys_parcel
      get_unit_kernel <- function(idx) kernels_parcel[[idx]]
      get_unit_ar <- function(idx) ar_parcel[[idx]]
      get_unit_lambda <- function(idx) lam_parcel[idx]
      Z_unit_r <- z_unit_parcel_r
    } else {
      get_unit_kernel <- strat$get_kernel
      get_unit_ar <- if (!is.null(strat$get_ar)) strat$get_ar else function(idx) rho[idx]
      get_unit_lambda <- function(idx) lam_vec_r[idx]
      Z_unit_r <- if (has_z_vox && track_unit_cor) Z_vox_r else NULL
    }
    if (!is.null(gamma_haus_override) && any(valid)) {
      gam_q <- round(gamma_haus_override * 8.0) / 8.0
      if (length(unique(gam_q[valid])) > 1L) {
        keys <- paste(keys, sprintf("%.3f", gam_q), sep = "|gam:")
      }
    }
    stratum_map[[r]] <- keys
    sig_eff_mat[r, ] <- sig_eff_r

    # Build readout matrices:
    # - S_lags[[li]]: standard readout centered at onset + lag (used for raw readouts
    #   and for impulse/window models).
    # - S_lags_model[[li]]: for source_model = "duration", centers the readout at the
    #   trial boxcar midpoint (onset + lag + 0.5 * d_i) with trial-specific width
    #   max(w_src, d_i), aligning sustained processes of varying duration at lag = 0.
    dur_kernel <- if (readout_kernel == "linear" && source_model == "duration" && any(dur_target > 0)) {
      "boxcar"
    } else {
      readout_kernel
    }
    if (source_model == "duration" && any(duration > 0)) {
      S_all_0 <- build_epoch_readout_matrix(ev$onset + 0.5 * duration, 0, n, tr, offset,
        readout_kernel = dur_kernel, readout_width = pmax(w_src, duration))
    } else {
      S_all_0 <- build_epoch_readout_matrix(ev$onset, 0, n, tr, offset,
        readout_kernel = readout_kernel, readout_width = readout_width)
    }
    S_lags <- lapply(lags, function(lag_val) {
      build_epoch_readout_matrix(ev$onset[t_idx], lag_val, n, tr, offset,
        readout_kernel = readout_kernel, readout_width = readout_width)
    })
    S_lags_model <- if (source_model == "duration") {
      lapply(lags, function(lag_val) {
        build_epoch_readout_matrix(ev$onset[t_idx] + 0.5 * dur_target, lag_val, n, tr, offset,
          readout_kernel = dur_kernel, readout_width = pmax(w_src, dur_target))
      })
    } else {
      S_lags
    }

    for (li in seq_len(n_lags)) {
      raw_epochs[[li]][row_span, valid] <- S_lags[[li]] %*% act_r[, valid, drop = FALSE]
      demixed_epochs[[li]][row_span, valid] <- S_lags_model[[li]] %*% act_r[, valid, drop = FALSE]
    }

    get_cached_sub_op <- function(s_key, s_kern, s_ar, s_lam) {
      g_key <- paste0(r, "|", sub("\\|gam:.*$", "", s_key))
      if (exists(g_key, envir = global_op_cache, inherits = FALSE)) {
        get(g_key, envir = global_op_cache)
      } else {
        op_s <- hrf_inverse_operator(
          n, s_kern[, 1L], rows, s_ar, z, settings,
          obs_weights = fw_r, lambda = s_lam
        )
        assign(g_key, op_s, envir = global_op_cache)
        op_s
      }
    }
    U_ref_by_key <- if (!is.null(Z_unit_r)) list() else NULL
    u_skeys_r <- unique(keys[valid])
    stratum_specs_r <- vector("list", length(u_skeys_r))

    for (si_r in seq_along(u_skeys_r)) {
      key <- u_skeys_r[si_r]
      base_key <- sub("\\|gam:.*$", "", key)
      voxels <- which(keys == key & valid)
      ar <- get_unit_ar(voxels[1L])
      hrf_k <- get_unit_kernel(voxels[1L])
      lam_k <- get_unit_lambda(voxels[1L])
      mix_spec <- mix_by_key[[base_key]]
      if (!is.null(mix_spec)) {
        op_parts <- lapply(mix_spec, function(ms) get_cached_sub_op(ms$key, ms$kernel, ms$ar, ms$lambda))
        operator <- matrix(0.0, n, length(rows))
        for (mi in seq_along(mix_spec)) {
          operator <- operator + mix_spec[[mi]]$weight * op_parts[[mi]]
        }
      } else {
        op_parts <- NULL
        operator <- get_cached_sub_op(key, hrf_k, ar, lam_k)
      }

      eval_stratum_forward <- function(ons_v, dur_v, scale_int_idx = integer(0L),
                                       scale_int_div = numeric(0L), win_div = 1.0) {
        m_ev <- length(ons_v)
        build_one <- function(k_mat) {
          G_u <- hrf_event_design_cpp(n, tr, offset, ons_v, dur_v,
            seq_len(m_ev), m_ev, k_mat)[[1L]][rows, , drop = FALSE]
          if (length(scale_int_idx) > 0L) {
            G_u[, scale_int_idx] <- sweep(G_u[, scale_int_idx, drop = FALSE], 2L, scale_int_div, "/")
          }
          if (win_div != 1.0 && win_div > 0) {
            G_u <- G_u / win_div
          }
          G_c <- if (settings$scaling != "none") center_cols_weighted(G_u, fw_r) else G_u
          list(G_unc = G_u, G_cent = G_c)
        }
        base_g <- build_one(hrf_k)
        if (is.null(mix_spec)) {
          list(
            G_unc = base_g$G_unc,
            G_cent = base_g$G_cent,
            DG_unc = operator %*% base_g$G_unc,
            DG_cent = operator %*% base_g$G_cent
          )
        } else {
          DG_u <- matrix(0.0, n, m_ev)
          DG_c <- matrix(0.0, n, m_ev)
          for (mi in seq_along(mix_spec)) {
            gm <- build_one(mix_spec[[mi]]$kernel)
            wm <- mix_spec[[mi]]$weight
            DG_u <- DG_u + wm * (op_parts[[mi]] %*% gm$G_unc)
            DG_c <- DG_c + wm * (op_parts[[mi]] %*% gm$G_cent)
          }
          list(G_unc = base_g$G_unc, G_cent = base_g$G_cent, DG_unc = DG_u, DG_cent = DG_c)
        }
      }

      # Unclipped causal BOLD reference energy per target trial for boundary/censor gating
      if (source_model == "duration") {
        n_ref_g <- nrow(hrf_k) + max(4L, as.integer(ceiling(max(c(0.0, dur_target)) / tr)) + 4L)
        G_ref_t <- hrf_event_design_cpp(n_ref_g, tr, 0.0, rep(0.0, m_target), dur_target,
          seq_len(m_target), m_target, hrf_k)[[1L]]
        if (duration_scale == "integral" && any(dur_target > 0)) {
          dp_ref <- which(dur_target > 0)
          G_ref_t[, dp_ref] <- sweep(G_ref_t[, dp_ref, drop = FALSE], 2L, dur_target[dp_ref] / tr, "/")
        }
        g2_full_target <- pmax(colSums(G_ref_t^2), 1e-12)
      } else {
        dur_ref1 <- if (source_model == "window") w_src else 0.0
        n_ref_g <- nrow(hrf_k) + max(4L, as.integer(ceiling(dur_ref1 / tr)) + 4L)
        G_ref_1 <- hrf_event_design_cpp(n_ref_g, tr, 0.0, 0.0, dur_ref1, 1L, 1L, hrf_k)[[1L]]
        if (source_model == "window" && w_src > 0) G_ref_1 <- G_ref_1 / (w_src / tr)
        g2_full_target <- rep(max(sum(G_ref_1^2), 1e-12), m_target)
      }

      # Forward design for all events at their own onsets (for non-target leakage removal)
      dur_all_model <- duration
      if (source_model == "impulse") dur_all_model[t_idx] <- 0
      pos_t_all <- if (source_model == "duration" && duration_scale == "integral") t_idx[dur_target > 0] else integer(0L)
      div_t_all <- if (length(pos_t_all)) (dur_target[dur_target > 0] / tr) else numeric(0L)
      fw_all_0 <- eval_stratum_forward(ev$onset, dur_all_model, pos_t_all, div_t_all)
      G_all_0 <- fw_all_0$G_cent
      DG_all_0 <- fw_all_0$DG_cent

      T_whiten_cached <- NULL
      get_stratum_T_whiten <- function() {
        if (!is.null(T_whiten_cached)) return(T_whiten_cached)
        tw_key <- paste0(r, "|", base_key)
        if (exists(tw_key, envir = global_tw_cache, inherits = FALSE)) {
          T_whiten_cached <<- get(tw_key, envir = global_tw_cache)
          return(T_whiten_cached)
        }
        U_op <- compute_ar1_operator_factor(
          operator, rows, ar, center_obs = center_obs,
          obs_weights = fw_r, weight_space = w_space
        )
        sv_op <- hrf_safe_svd(U_op, nu = 0L, nv = ncol(U_op))
        k_sv <- min(length(sv_op$d), ncol(sv_op$v))
        if (k_sv < 1L) return(NULL)
        d_sv <- sv_op$d[seq_len(k_sv)]
        v_sv <- sv_op$v[, seq_len(k_sv), drop = FALSE]
        if (max(d_sv) <= 1e-12) return(NULL)
        keep_sv <- d_sv > 1e-5 * max(d_sv)
        if (!any(keep_sv)) return(NULL)
        T_whiten_cached <<- sweep(t(v_sv[, keep_sv, drop = FALSE]), 1L, d_sv[keep_sv], "/")
        assign(tw_key, T_whiten_cached, envir = global_tw_cache)
        T_whiten_cached
      }

      L_joint <- if (demix_other_events && length(o_idx) > 0L) {
        tw_oth <- if (!identical(other_model_eff, "condition")) get_stratum_T_whiten() else NULL
        build_other_demix_operator(S_all_0, DG_all_0, ev, o_idx, t_idx, demix_ridge,
          source_model = source_model, demix_target_trials = demix_target_trials,
          other_event_model = other_model_eff, T_whiten = tw_oth)
      } else {
        NULL
      }
      gamma_oth <- 0.0
      if (!is.null(L_joint) && identical(L_joint$mode, "auto") && !is.null(L_joint$L_delta_o) &&
          nrow(L_joint$L_delta_o) >= 4L) {
        if (!is.null(gamma_haus_override)) {
          gamma_oth <- mean(gamma_haus_override[voxels], na.rm = TRUE)
          if (!is.finite(gamma_oth)) gamma_oth <- 0.0
        } else {
          # Satterthwaite & cross-run pooled non-target confounding screen
          cond_o_f <- as.character(L_joint$cond_o)
          u_co <- unique(cond_o_f)
          m_o_ok <- nrow(L_joint$L_delta_o)
          Q_o <- diag(m_o_ok)
          for (uc in u_co) {
            ic <- which(cond_o_f == uc)
            if (length(ic) >= 1L) {
              Q_o[ic, ic] <- Q_o[ic, ic] - matrix(1.0 / length(ic), length(ic), length(ic))
            }
          }
          df_o_nom <- m_o_ok - length(u_co)
          if (df_o_nom >= 3L) {
            w_v <- 1.0 / pmax(sig_eff_r[voxels]^2, 1e-12)
            w_v <- w_v / sum(w_v)
            a_bar <- as.numeric(act_r[, voxels, drop = FALSE] %*% w_v)
            sig2_bar <- if (!is.null(Z_unit_r) && ncol(Z_unit_r) >= max(voxels)) {
              sum(as.numeric(Z_unit_r[, voxels, drop = FALSE] %*% (w_v * sig_eff_r[voxels]))^2)
            } else if (!isTRUE(aggregate_by_parcel) && has_z_vox && ncol(Z_vox_r) >= max(voxels)) {
              sum(as.numeric(Z_vox_r[, voxels, drop = FALSE] %*% (w_v * sig_eff_r[voxels]))^2)
            } else {
              sum((w_v * sig_eff_r[voxels])^2)
            }
            sig2_bar_raw <- max(sig2_bar, 1e-16)
            emp_mult_o <- 1.0
            if (!is.null(L_joint$proj_resid_info)) {
              pri <- L_joint$proj_resid_info
              y_w_bar <- as.numeric(pri$T_whiten %*% a_bar)
              res_w_bar <- y_w_bar - as.numeric(pri$F_u %*% (pri$H_1 %*% crossprod(pri$F_u, y_w_bar)))
              df_w_res <- max(8.0, as.numeric(length(y_w_bar)) - as.numeric(pri$eff_df))
              sig2_emp <- sum(res_w_bar^2) / df_w_res
              if (is.finite(sig2_emp) && sig2_emp > sig2_bar) {
                emp_mult_o <- sig2_emp / sig2_bar_raw
                sig2_bar <- sig2_emp
              }
            }
            sig2_bar <- max(sig2_bar, 1e-16)
            QL_do <- Q_o %*% L_joint$L_delta_o
            if (!is.null(L_joint$K_o_from_t) && !is.null(L_joint$L_t_pilot) &&
                ncol(L_joint$K_o_from_t) == nrow(L_joint$L_t_pilot)) {
              QL_do <- QL_do - (Q_o %*% L_joint$K_o_from_t) %*% L_joint$L_t_pilot
            }
            B_do <- QL_do %*% operator
            V_null_unit_o <- compute_ar1_operator_cov(
              B_do, rows, ar, center_obs = center_obs,
              obs_weights = fw_r, weight_space = w_space
            )
            V_null_o <- sig2_bar * V_null_unit_o
            if (n_lags >= 2L && m_target >= 3L && any(abs(lags) >= 1.5 * tr)) {
              ext_lags <- unique(c(min(lags), max(lags)))
              ext_lags <- ext_lags[abs(ext_lags) >= 1.5 * tr]
              for (lg_ex in ext_lags) {
                li_ex <- which.min(abs(lags - lg_ex))[1L]
                S_ex <- S_lags_model[[li_ex]]
                ok_ex <- rowSums(!is.finite(S_ex)) == 0L
                if (sum(ok_ex) >= 3L) {
                  dur_ex <- if (source_model == "duration") dur_target else rep(0.0, m_target)
                  fw_ex <- eval_stratum_forward(ev$onset[t_idx] + lg_ex, dur_ex)
                  S_ex_ok <- S_ex[ok_ex, , drop = FALSE]
                  DG_ex_ok <- fw_ex$DG_cent[, ok_ex, drop = FALSE]
                  g_ex <- max(0.20, mean(rowSums(S_ex_ok * t(DG_ex_ok)), na.rm = TRUE))
                  y_ex_pil <- as.numeric(S_ex_ok %*% a_bar)
                  B_ex_pil <- S_ex_ok %*% operator
                  v_err_ex <- mean(diag(sig2_bar * compute_ar1_operator_cov(
                    B_ex_pil, rows, ar, center_obs = center_obs,
                    obs_weights = fw_r, weight_space = w_space
                  )))
                  tau2_ex_up <- max(0.0, (stats::var(y_ex_pil) - v_err_ex) / (g_ex^2))
                  if (is.finite(tau2_ex_up) && tau2_ex_up > 0) {
                    QK_oex <- QL_do %*% DG_ex_ok
                    V_null_o <- V_null_o + tau2_ex_up * tcrossprod(QK_oex)
                  }
                }
              }
            }
            d_hat_o <- as.numeric(QL_do %*% a_bar)
            ss_obs_o <- sum(d_hat_o^2)
            ss_null_o <- sum(diag(V_null_o))
            df1_o <- max(1.0, (ss_null_o^2) / max(sum(V_null_o^2), 1e-16))
            df2_o <- max(16.0, as.numeric(length(rows) - nrow(ev)))
            f_rat_o <- ss_obs_o / max(ss_null_o, 1e-16)
            p_oth <- stats::pf(f_rat_o, df1 = df1_o, df2 = df2_o, lower.tail = FALSE)

            p_cov_o <- 1.0
            use_haus_1lag <- (n_lags == 1L) && identical(source_model, "impulse") &&
              !isTRUE(demix_within_trial) && (abs(lags[1L]) <= 1e-6) &&
              isTRUE(demix_target_trials) && !is.null(L_joint$L_target_trial)
            if (n_terms > 0L && is.matrix(cov_list[[r]]) && nrow(cov_list[[r]]) == m_target &&
                !is.null(L_joint$ok_o) && length(L_joint$ok_o) == m_o_ok) {
              S_chk_0 <- if (use_haus_1lag) S_lags_model[[1L]] else S_all_0[t_idx, , drop = FALSE]
              t_ok_m <- (rowSums(!is.finite(S_chk_0)) == 0L) &
                (rowSums(!is.finite(cov_list[[r]])) == 0L)
              if (use_haus_1lag) {
                fw_h1 <- eval_stratum_forward(ev$onset[t_idx], rep(0.0, m_target))
                obs_f_h1 <- colSums(fw_h1$G_unc^2) / g2_full_target
                t_ok_m <- t_ok_m & is.finite(obs_f_h1) & (obs_f_h1 >= 0.20)
              }
              m_t_ok <- sum(t_ok_m)
              if (m_t_ok >= n_terms + 3L) {
                Q_t_ok <- diag(m_t_ok) - matrix(1.0 / m_t_ok, m_t_ok, m_t_ok)
                X_r_c <- Q_t_ok %*% cov_list[[r]][t_ok_m, , drop = FALSE]
                keep_xc <- which(colSums(X_r_c^2) > 1e-12)
                if (length(keep_xc) > 0L) {
                  X_r_c <- X_r_c[, keep_xc, drop = FALSE]
                  q_xc <- ncol(X_r_c)
                  XtX_c <- crossprod(X_r_c)
                  reg_xtx <- 1e-10 * max(mean(diag(XtX_c)), 1e-12)
                  XtX_inv_c <- solve(XtX_c + reg_xtx * diag(q_xc))
                  df_lk <- max(1.0, as.numeric(m_t_ok - q_xc - 1L))

                  if (use_haus_1lag) {
                    R_oth_1 <- S_chk_0[t_ok_m, , drop = FALSE] %*% DG_all_0[, o_idx, drop = FALSE]
                    S_cd_1 <- S_chk_0[t_ok_m, , drop = FALSE] - R_oth_1 %*% L_joint$L_other_mean
                    dem_cd1 <- solve_duration_demix_operator(
                      S_cd_1, operator,
                      fw_h1$G_unc[, t_ok_m, drop = FALSE], fw_h1$G_cent[, t_ok_m, drop = FALSE],
                      settings$scaling, demix_ridge, gain_floor,
                      DG_tau_unc = fw_h1$DG_unc[, t_ok_m, drop = FALSE],
                      DG_tau_cent = fw_h1$DG_cent[, t_ok_m, drop = FALSE]
                    )
                    G_gain_1 <- sweep(Q_t_ok, 1L, dem_cd1$eff_gain_vec, "/") + (diag(m_t_ok) - Q_t_ok)
                    M_cd_1 <- G_gain_1 %*% (dem_cd1$W %*% S_cd_1)
                    M_tr_1 <- L_joint$L_target_trial[t_ok_m, , drop = FALSE]
                    Delta_M_1 <- Q_t_ok %*% (M_cd_1 - M_tr_1)
                    Delta_B_1 <- Delta_M_1 %*% operator
                    V_unit_1 <- compute_ar1_operator_cov(
                      Delta_B_1, rows, ar, center_obs = center_obs,
                      obs_weights = fw_r, weight_space = w_space
                    )
                    V_c_unit <- crossprod(X_r_c, V_unit_1 %*% X_r_c)
                    ss_lk_unit <- sum(diag(V_unit_1))
                    tr_xtx_vc_unit <- sum(XtX_inv_c * V_c_unit)
                    leak_t_v <- as.numeric(Delta_M_1 %*% a_bar)
                    V_c <- sig2_bar * V_c_unit
                    ss_lk_null <- sig2_bar * ss_lk_unit
                    L_vox_op <- Delta_M_1
                  } else {
                    M_leak_c <- Q_t_ok %*% (S_all_0[t_idx[t_ok_m], , drop = FALSE] %*% DG_all_0[, L_joint$ok_o, drop = FALSE])
                    MX_c <- crossprod(M_leak_c, X_r_c)
                    V_c_unit <- crossprod(MX_c, V_null_unit_o %*% MX_c)
                    ss_lk_unit <- sum(M_leak_c * (M_leak_c %*% V_null_unit_o))
                    tr_xtx_vc_unit <- sum(XtX_inv_c * V_c_unit)
                    leak_t_v <- as.numeric(M_leak_c %*% d_hat_o)
                    V_c <- crossprod(MX_c, V_null_o %*% MX_c)
                    ss_lk_null <- sum(M_leak_c * (M_leak_c %*% V_null_o))
                    L_vox_op <- M_leak_c %*% QL_do
                  }

                  num_c <- as.numeric(crossprod(X_r_c, leak_t_v))
                  b_lk <- as.numeric(XtX_inv_c %*% num_c)
                  res_lk <- leak_t_v - as.numeric(X_r_c %*% b_lk)
                  ss_res_null <- max(0.0, ss_lk_null - sum(XtX_inv_c * V_c))
                  ss_res_obs <- sum(res_lk^2)
                  s2_lk_ex <- max(0.0, (ss_res_obs - ss_res_null) / df_lk)
                  V_tot_c <- V_c + s2_lk_ex * XtX_c
                  w_obs_df <- ss_res_null / max(ss_res_obs, ss_res_null, 1e-16)
                  df2_cov <- max(4.0, w_obs_df * df2_o + (1.0 - w_obs_df) * df_lk)
                  reg_vc <- 1e-10 * max(mean(diag(V_tot_c)), 1e-12)
                  f_cov_o <- tryCatch(
                    as.numeric(crossprod(num_c, solve(V_tot_c + reg_vc * diag(q_xc), num_c))) / q_xc,
                    error = function(e) NA_real_
                  )
                  if (is.finite(f_cov_o) && f_cov_o >= 0) {
                    p_cov_o <- stats::pf(f_cov_o, df1 = q_xc, df2 = df2_cov, lower.tail = FALSE)
                  }

                  if (!is.null(haus_W)) {
                    if (length(voxels) == 1L) {
                      V_b_r <- XtX_inv_c %*% V_tot_c %*% XtX_inv_c
                      reg_vb <- 1e-10 * max(mean(diag(V_b_r)), 1e-12)
                      Prec_b_r <- tryCatch(solve(V_b_r + reg_vb * diag(q_xc)), error = function(e) NULL)
                      if (!is.null(Prec_b_r) && all(is.finite(Prec_b_r))) {
                        vi_h <- voxels[1L]
                        haus_W[[vi_h]][keep_xc, keep_xc] <- haus_W[[vi_h]][keep_xc, keep_xc] + Prec_b_r
                        haus_rhs[[vi_h]][keep_xc] <- haus_rhs[[vi_h]][keep_xc] + as.numeric(Prec_b_r %*% b_lk)
                        haus_df2[vi_h] <- haus_df2[vi_h] + df2_cov
                      }
                    } else {
                      L_vox <- L_vox_op %*% act_r[, voxels, drop = FALSE]
                      Num_vox <- crossprod(X_r_c, L_vox)
                      B_vox <- XtX_inv_c %*% Num_vox
                      Res_vox <- L_vox - X_r_c %*% B_vox
                      ss_res_vox <- colSums(Res_vox^2)
                      for (jj_v in seq_along(voxels)) {
                        vi_h <- voxels[jj_v]
                        s2_u <- max(sig_eff_r[vi_h]^2, 1e-16) * emp_mult_o
                        ss_rn_v <- max(0.0, s2_u * (ss_lk_unit - tr_xtx_vc_unit))
                        s2_ex_v <- max(0.0, (ss_res_vox[jj_v] - ss_rn_v) / df_lk)
                        V_tot_v <- s2_u * V_c_unit + s2_ex_v * XtX_c
                        V_b_v <- XtX_inv_c %*% V_tot_v %*% XtX_inv_c
                        reg_vb <- 1e-10 * max(mean(diag(V_b_v)), 1e-12)
                        Prec_b_v <- tryCatch(solve(V_b_v + reg_vb * diag(q_xc)), error = function(e) NULL)
                        if (!is.null(Prec_b_v) && all(is.finite(Prec_b_v))) {
                          w_df_v <- ss_rn_v / max(ss_res_vox[jj_v], ss_rn_v, 1e-16)
                          df2_v <- max(4.0, w_df_v * df2_o + (1.0 - w_df_v) * df_lk)
                          haus_W[[vi_h]][keep_xc, keep_xc] <- haus_W[[vi_h]][keep_xc, keep_xc] + Prec_b_v
                          haus_rhs[[vi_h]][keep_xc] <- haus_rhs[[vi_h]][keep_xc] + as.numeric(Prec_b_v %*% B_vox[, jj_v])
                          haus_df2[vi_h] <- haus_df2[vi_h] + df2_v
                        }
                      }
                    }
                  }
                }
              }
            }

            if (n_terms > 0L) {
              if (is.finite(p_cov_o) && p_cov_o < 0.045) {
                gamma_oth <- .poet_firm_gate(
                  -log10(max(p_cov_o, 1e-16)),
                  lo = -log10(0.045),
                  hi = -log10(0.008)
                )
              }
            } else if (is.finite(f_rat_o) && f_rat_o > 3.0 && is.finite(p_oth) && p_oth < 0.001) {
              sig2_unit_ref <- stats::median(sig_eff_r[voxels]^2) * emp_mult_o
              var_mult_u <- max(1.0, sig2_unit_ref / sig2_bar)
              gamma_oth <- max(0.0, min(1.0, (ss_obs_o - ss_null_o) / max(ss_obs_o + (var_mult_u - 1.0) * ss_null_o, 1e-16)))
            }
            if (!is.null(haus_gam_p1)) haus_gam_p1[r, voxels] <- gamma_oth
          }
        }
        if (gamma_oth > 0.0) {
          L_joint$L_other <- (1.0 - gamma_oth) * L_joint$L_other_mean + gamma_oth * L_joint$L_other_trial
        }
      }
      L_other <- if (!is.null(L_joint)) L_joint$L_other else NULL
      L_target_0 <- if (!is.null(L_joint)) L_joint$L_target_0 else NULL
      a_other_hat <- if (!is.null(L_other)) {
        L_other %*% act_r[, voxels, drop = FALSE]
      } else {
        NULL
      }

      # For duration source_model, estimate baseline firing rate alpha0_hat at lag = 0
      # so differential time-on-task baseline leakage is removed across all lags
      alpha0_hat <- NULL
      b_alpha0 <- NULL
      s_alpha0 <- NULL
      DG_target_0_sum <- NULL
      if (source_model == "duration") {
        S_0 <- build_epoch_readout_matrix(ev$onset[t_idx] + 0.5 * dur_target, 0, n, tr, offset,
          readout_kernel = dur_kernel, readout_width = pmax(w_src, dur_target))
        dp0_idx <- if (duration_scale == "integral" && any(dur_target > 0)) which(dur_target > 0) else integer(0L)
        dp0_div <- if (length(dp0_idx)) (dur_target[dp0_idx] / tr) else numeric(0L)
        fw_0 <- eval_stratum_forward(ev$onset[t_idx], dur_target, dp0_idx, dp0_div)
        G_0_unc <- fw_0$G_unc
        G_0_cent <- fw_0$G_cent
        obs_frac_0 <- colSums(G_0_unc^2) / g2_full_target
        ok_0 <- (rowSums(!is.finite(S_0)) == 0L) & is.finite(obs_frac_0) & (obs_frac_0 >= 0.20)
        if (any(ok_0)) {
          DG_target_0_sum <- rowSums(fw_0$DG_cent[, ok_0, drop = FALSE])
          Y_0 <- S_0[ok_0, , drop = FALSE] %*% act_r[, voxels, drop = FALSE]
          S_clean_0 <- S_0[ok_0, , drop = FALSE]
          if (!is.null(L_other)) {
            R_other_0 <- S_0[ok_0, , drop = FALSE] %*% DG_all_0[, o_idx, drop = FALSE]
            Y_0 <- Y_0 - R_other_0 %*% a_other_hat
            S_clean_0 <- S_clean_0 - R_other_0 %*% L_other
          }
          demix_0 <- solve_duration_demix_operator(
            S_clean_0, operator,
            G_0_unc[, ok_0, drop = FALSE], G_0_cent[, ok_0, drop = FALSE],
            settings$scaling, demix_ridge, gain_floor,
            DG_tau_unc = fw_0$DG_unc[, ok_0, drop = FALSE],
            DG_tau_cent = fw_0$DG_cent[, ok_0, drop = FALSE]
          )
          alpha0_hat <- colMeans(demix_0$W %*% Y_0)
          s_alpha0 <- colMeans(demix_0$W %*% S_clean_0)
          b_alpha0 <- as.numeric(s_alpha0 %*% operator)
        }
      }

      ok_tau_list <- vector("list", n_lags)
      B_demix_list <- vector("list", n_lags)
      M_act_list <- vector("list", n_lags)
      DG_tau_list <- vector("list", n_lags)

      for (li in seq_len(n_lags)) {
        # Raw readout operator error variance
        S_raw_l <- S_lags[[li]]
        ok_raw <- rowSums(!is.finite(S_raw_l)) == 0L
        if (any(ok_raw)) {
          m_raw <- sum(ok_raw)
          B_raw <- S_raw_l[ok_raw, , drop = FALSE] %*% operator
          V_raw <- compute_ar1_operator_cov(
            B_raw, rows, ar, center_obs = center_obs,
            obs_weights = fw_r, weight_space = w_space
          )
          if (m_raw > 1L) {
            Q_raw <- diag(m_raw) - matrix(1 / m_raw, m_raw, m_raw)
            V_raw_cent <- (Q_raw %*% V_raw %*% Q_raw) / (1 - 1 / m_raw)
          } else {
            V_raw_cent <- V_raw
          }
          d_raw_unit <- sqrt(pmax(0, diag(V_raw_cent)))
          se_raw_epochs[[li]][row_span[ok_raw], voxels] <- outer(d_raw_unit, sig_eff_r[voxels])
        }

        if (source_model == "duration") {
          shifted_onsets <- ev$onset[t_idx] + lags[li]
          dp_idx <- if (duration_scale == "integral" && any(dur_target > 0)) which(dur_target > 0) else integer(0L)
          dp_div <- if (length(dp_idx)) (dur_target[dp_idx] / tr) else numeric(0L)
          fw_tau <- eval_stratum_forward(shifted_onsets, dur_target, dp_idx, dp_div)
        } else {
          # Local target resolution matrix K_tau for unit sources at causal target onset + lag
          dur_tau <- if (source_model == "window") rep(w_src, m_target) else rep(0, m_target)
          shift_center <- if (source_model == "window") lags[li] - 0.5 * w_src else lags[li]
          shifted_onsets <- ev$onset[t_idx] + shift_center
          win_div_tau <- if (source_model == "window" && w_src > 0) (w_src / tr) else 1.0
          fw_tau <- eval_stratum_forward(shifted_onsets, dur_tau, win_div = win_div_tau)
        }
        G_tau_unc <- fw_tau$G_unc
        G_tau <- fw_tau$G_cent
        DG_tau_list[[li]] <- fw_tau$DG_cent

        S_tau <- S_lags_model[[li]]
        obs_frac_tau <- colSums(G_tau_unc^2) / g2_full_target
        ok_tau <- (rowSums(!is.finite(S_tau)) == 0L) & is.finite(obs_frac_tau) & (obs_frac_tau >= 0.20)
        ok_tau_list[[li]] <- ok_tau
        Y_clean <- demixed_epochs[[li]][row_span, voxels, drop = FALSE]
        if (any(!ok_tau)) {
          Y_clean[!ok_tau, ] <- NA_real_
          demixed_epochs[[li]][row_span[!ok_tau], voxels] <- NA_real_
        }
        if (!any(ok_tau)) next
        m_ok <- sum(ok_tau)
        S_clean_tau <- S_tau[ok_tau, , drop = FALSE]

        if (!is.null(L_other)) {
          R_other_tau <- S_tau[ok_tau, , drop = FALSE] %*% DG_all_0[, o_idx, drop = FALSE]
          Y_clean[ok_tau, ] <- Y_clean[ok_tau, , drop = FALSE] - R_other_tau %*% a_other_hat
          S_clean_tau <- S_clean_tau - R_other_tau %*% L_other
        }
        B_clean_tau <- S_clean_tau %*% operator

        demix_tau <- solve_duration_demix_operator(
          S_clean_tau, operator,
          G_tau_unc[, ok_tau, drop = FALSE], G_tau[, ok_tau, drop = FALSE],
          settings$scaling, demix_ridge, gain_floor,
          DG_tau_unc = fw_tau$DG_unc[, ok_tau, drop = FALSE],
          DG_tau_cent = fw_tau$DG_cent[, ok_tau, drop = FALSE]
        )
        K_norm <- demix_tau$K_norm
        g_vec <- diag(K_norm)
        mean_g <- mean(g_vec[is.finite(g_vec)])
        if (is.finite(mean_g)) {
          gain_accum[li, voxels] <- gain_accum[li, voxels] + mean_g
          gain_count[li, voxels] <- gain_count[li, voxels] + 1L
        }
        if (demix_target_trials) {
          if (identical(source_model, "impulse") && !isTRUE(demix_within_trial) &&
              abs(lags[li]) <= 1e-6 && !is.null(L_target_0)) {
            L_t0_ok <- L_target_0[ok_tau, , drop = FALSE]
            Y_clean[ok_tau, ] <- L_t0_ok %*% act_r[, voxels, drop = FALSE]
            M_act_tau <- L_t0_ok
            B_demix_tau <- L_t0_ok %*% operator
          } else if (n_lags == 1L && identical(source_model, "impulse") && !isTRUE(demix_within_trial) &&
                     abs(lags[li]) <= 1e-6 && !is.null(L_joint) &&
                     !is.null(L_joint$L_target_trial) && gamma_oth > 0.0) {
            R_oth_0b <- S_tau[ok_tau, , drop = FALSE] %*% DG_all_0[, o_idx, drop = FALSE]
            S_cd_0b <- S_tau[ok_tau, , drop = FALSE] - R_oth_0b %*% L_joint$L_other_mean
            dem_cd_0b <- solve_duration_demix_operator(
              S_cd_0b, operator,
              G_tau_unc[, ok_tau, drop = FALSE], G_tau[, ok_tau, drop = FALSE],
              settings$scaling, demix_ridge, gain_floor,
              DG_tau_unc = fw_tau$DG_unc[, ok_tau, drop = FALSE],
              DG_tau_cent = fw_tau$DG_cent[, ok_tau, drop = FALSE]
            )
            Q_tau <- diag(m_ok) - matrix(1 / m_ok, m_ok, m_ok)
            G_gain_cd <- sweep(Q_tau, 1L, dem_cd_0b$eff_gain_vec, "/") + (diag(m_ok) - Q_tau)
            M_cd_tau <- G_gain_cd %*% (dem_cd_0b$W %*% S_cd_0b)
            L_t0_ok <- L_joint$L_target_trial[ok_tau, , drop = FALSE]
            M_act_tau <- (1.0 - gamma_oth) * M_cd_tau + gamma_oth * L_t0_ok
            Y_clean[ok_tau, ] <- M_act_tau %*% act_r[, voxels, drop = FALSE]
            B_demix_tau <- M_act_tau %*% operator
          } else {
            Y_sub <- Y_clean[ok_tau, , drop = FALSE]
            B_sub_tau <- B_clean_tau
            S_sub_act <- S_clean_tau
            if (source_model == "duration" && !is.null(alpha0_hat) && !is.null(DG_target_0_sum) && abs(lags[li]) > 1e-8) {
              k_base_0 <- as.numeric(S_tau[ok_tau, , drop = FALSE] %*% DG_target_0_sum)
              k_base_tau <- rowSums(K_norm)
              # Retain only the physical temporal overlap fraction gamma_tau in [0, 1]
              # between shifted windows [onset + lag, onset + lag + d_i] and [onset, onset + d_i]
              # so pre-onset (lag < 0) and post-offset lags never receive artificial +alpha0_hat
              # baseline shifts, while W_tau %*% k_base_tau = 1_m keeps centered trial contrasts invariant.
              gamma_tau <- mean(pmax(0.0, pmin(1.0, 1.0 - abs(lags[li]) / pmax(dur_target[ok_tau], w_src, 1e-6))))
              delta_base <- k_base_0 - gamma_tau * k_base_tau
              Y_sub <- Y_sub - outer(delta_base, alpha0_hat)
              if (!is.null(b_alpha0)) {
                B_sub_tau <- B_sub_tau - tcrossprod(delta_base, b_alpha0)
              }
              if (!is.null(s_alpha0)) {
                S_sub_act <- S_sub_act - tcrossprod(delta_base, s_alpha0)
              }
            }
            A_hat <- demix_tau$W %*% Y_sub
            A_mean <- colMeans(A_hat)
            Y_clean[ok_tau, ] <- sweep(
              sweep(sweep(A_hat, 2L, A_mean, "-"), 1L, demix_tau$eff_gain_vec, "/"),
              2L, A_mean, "+"
            )
            Q_tau <- diag(m_ok) - matrix(1 / m_ok, m_ok, m_ok)
            G_gain_tau <- sweep(Q_tau, 1L, demix_tau$eff_gain_vec, "/") + (diag(m_ok) - Q_tau)
            M_act_tau <- G_gain_tau %*% (demix_tau$W %*% S_sub_act)
            B_demix_tau <- G_gain_tau %*% (demix_tau$W %*% B_sub_tau)
          }
        } else {
          M_act_tau <- S_clean_tau
          B_demix_tau <- B_clean_tau
        }
        demixed_epochs[[li]][row_span, voxels] <- Y_clean
        B_demix_list[[li]] <- B_demix_tau
        M_act_list[[li]] <- M_act_tau
      }

      # Cross-lag point-spread resolution matrix C_sym and optional within-trial de-mixing
      ok_common <- NULL
      m_common <- 0L
      if (n_lags >= 2L) {
        has_all_lags <- all(vapply(ok_tau_list, function(x) !is.null(x) && length(x) == m_target, logical(1L)))
        if (has_all_lags) {
          ok_common <- Reduce("&", ok_tau_list)
          m_common <- sum(ok_common)
          if (m_common >= 2L) {
            K_lag <- matrix(0, n_lags, n_lags)
            for (l1 in seq_len(n_lags)) {
              S_l1_sub <- S_lags_model[[l1]][ok_common, , drop = FALSE]
              for (l2 in seq_len(n_lags)) {
                DG_l2_sub <- DG_tau_list[[l2]][, ok_common, drop = FALSE]
                k_diag <- rowSums(S_l1_sub * t(DG_l2_sub))
                K_lag[l1, l2] <- mean(k_diag[is.finite(k_diag)])
              }
            }
            d_lag <- pmax(diag(K_lag), gain_floor)
            D_inv_sqrt <- diag(1 / sqrt(d_lag), n_lags)
            C_lag <- D_inv_sqrt %*% K_lag %*% D_inv_sqrt
            C_sym <- 0.5 * (C_lag + t(C_lag))
            for (vi in voxels) {
              lag_res_accum[[vi]] <- lag_res_accum[[vi]] + m_common * C_sym
              lag_res_trials[vi] <- lag_res_trials[vi] + m_common
            }

            dt_lag <- if (n_lags > 1L) min(diff(sort(unique(lags)))) else tr
            dt_ref <- max(tr, 1.0)
            if (isTRUE(demix_within_trial) && n_lags >= 2L && is.finite(dt_lag) && dt_lag >= dt_ref) {
              idx_common <- which(ok_common)
              Q_c <- diag(m_common) - matrix(1.0 / m_common, m_common, m_common)
              tw_wt <- get_stratum_T_whiten()
              act_sub_wt <- act_r[, voxels, drop = FALSE]

              # 1. Remove ITI-varying cross-trial/cross-lag baseline spillover on sparse causal multi-lag grids
              #    using the joint GLS lag-mean estimator (e.g., Trial(i+1, 0s) mean landing at Trial(i, +4s)).
              pos_lags_wt <- which(lags >= -1e-6)
              if (dt_lag >= 1.5 * dt_ref && length(pos_lags_wt) >= 2L &&
                  !is.null(tw_wt) && nrow(tw_wt) >= length(pos_lags_wt) + 4L && ncol(tw_wt) == n) {
                T_clean_wt <- if (!is.null(L_other)) {
                  tw_wt - (tw_wt %*% DG_all_0[, o_idx, drop = FALSE]) %*% L_other
                } else {
                  tw_wt
                }
                DG_lag_sums <- vapply(pos_lags_wt, function(li_p) {
                  rowSums(DG_tau_list[[li_p]][, idx_common, drop = FALSE])
                }, numeric(n))
                if (is.null(dim(DG_lag_sums))) DG_lag_sums <- matrix(DG_lag_sums, ncol = length(pos_lags_wt))
                F_lag_means <- T_clean_wt %*% DG_lag_sums
                fm_diag <- colSums(F_lag_means^2)
                fm_sc <- if (any(is.finite(fm_diag) & fm_diag > 0)) mean(fm_diag[is.finite(fm_diag) & fm_diag > 0]) else 1.0
                mu_lag_act_op <- solve(
                  crossprod(F_lag_means) + max(0.04 * fm_sc, 1e-6) * diag(length(pos_lags_wt)),
                  crossprod(F_lag_means, T_clean_wt)
                )
                mu_hat_wt <- mu_lag_act_op %*% act_sub_wt
                for (li in seq_len(n_lags)) {
                  ok_sub_li <- ok_common[ok_tau_list[[li]]]
                  M_li_old <- M_act_list[[li]][ok_sub_li, , drop = FALSE]
                  spill_li <- Q_c %*% (M_li_old %*% DG_lag_sums)
                  M_li_new <- M_li_old - spill_li %*% mu_lag_act_op
                  demixed_epochs[[li]][row_span[idx_common], voxels] <-
                    demixed_epochs[[li]][row_span[idx_common], voxels, drop = FALSE] - spill_li %*% mu_hat_wt
                  M_act_list[[li]][ok_sub_li, ] <- M_li_new
                  B_demix_list[[li]][ok_sub_li, ] <- M_li_new %*% operator
                }
              }

              # 2. Invert the exact post-trial-demixing n_lags x n_lags cross-lag resolution matrix
              #    R_post using first-power Neumann refinement (avoiding crossprod squaring of R_post).
              R_post <- matrix(0.0, n_lags, n_lags)
              for (l1 in seq_len(n_lags)) {
                ok_sub_l1 <- ok_common[ok_tau_list[[l1]]]
                QM_l1 <- Q_c %*% M_act_list[[l1]][ok_sub_l1, , drop = FALSE]
                for (l2 in seq_len(n_lags)) {
                  DGQ_l2 <- DG_tau_list[[l2]][, idx_common, drop = FALSE] %*% Q_c
                  R_post[l1, l2] <- sum(QM_l1 * t(DGQ_l2)) / max(1.0, m_common - 1.0)
                }
              }
              d_post <- pmax(abs(diag(R_post)), gain_floor) * ifelse(diag(R_post) < 0, -1.0, 1.0)
              C_asym <- sweep(R_post, 1L, d_post, "/")
              eig_sym_c <- eigen(0.5 * (C_asym + t(C_asym)), symmetric = TRUE, only.values = TRUE)$values
              if (all(is.finite(eig_sym_c)) && min(eig_sym_c) > 0.02) {
                r_1 <- if (!within_ridge_supplied) 0.065 else within_trial_ridge
                W_1 <- solve(C_asym + r_1 * diag(n_lags))
                W_2 <- W_1 + 0.65 * W_1 %*% (diag(n_lags) - C_asym %*% W_1)
              } else {
                eff_wt_ridge <- if (!within_ridge_supplied && dt_lag >= 1.5 * dt_ref) 0.015 else within_trial_ridge
                W_1 <- solve(crossprod(C_asym) + eff_wt_ridge * diag(n_lags), t(C_asym))
                W_2 <- W_1 + 0.80 * W_1 %*% (diag(n_lags) - C_asym %*% W_1)
              }
              W_raw <- sweep(W_2, 2L, d_post, "/")
              W_lag <- sweep(W_raw, 1L, pmax(diag(W_raw %*% R_post), gain_floor), "/")
              Y_common_old <- lapply(seq_len(n_lags), function(li) demixed_epochs[[li]][row_span[idx_common], voxels, drop = FALSE])
              B_common_old <- lapply(seq_len(n_lags), function(li) B_demix_list[[li]][ok_common[ok_tau_list[[li]]], , drop = FALSE])
              M_common_old <- lapply(seq_len(n_lags), function(li) M_act_list[[li]][ok_common[ok_tau_list[[li]]], , drop = FALSE])
              n_obs_b <- ncol(B_common_old[[1L]])
              for (l1 in seq_len(n_lags)) {
                Y_new_l1 <- matrix(0, m_common, length(voxels))
                B_new_l1 <- matrix(0, m_common, n_obs_b)
                M_new_l1 <- matrix(0, m_common, n)
                for (l2 in seq_len(n_lags)) {
                  w_12 <- W_lag[l1, l2]
                  if (abs(w_12) > 0) {
                    Y_new_l1 <- Y_new_l1 + w_12 * Y_common_old[[l2]]
                    B_new_l1 <- B_new_l1 + w_12 * B_common_old[[l2]]
                    M_new_l1 <- M_new_l1 + w_12 * M_common_old[[l2]]
                  }
                }
                demixed_epochs[[l1]][row_span[idx_common], voxels] <- Y_new_l1
                ok_sub_l1 <- ok_common[ok_tau_list[[l1]]]
                B_demix_list[[l1]][ok_sub_l1, ] <- B_new_l1
                M_act_list[[l1]][ok_sub_l1, ] <- M_new_l1
              }
            } else if (isTRUE(demix_within_trial) && n_lags >= 3L) {
              sub_tr_ratio <- if (is.finite(dt_lag) && dt_lag > 0 && dt_lag < dt_ref) {
                (dt_ref / max(dt_lag, 0.1 * dt_ref))^2
              } else {
                1.0
              }
              eff_wt_ridge <- if (!within_ridge_supplied && sub_tr_ratio > 1.0) {
                within_trial_ridge * sub_tr_ratio
              } else {
                within_trial_ridge
              }
              neumann_wt <- 0.6 / sub_tr_ratio
              W_1 <- solve(crossprod(C_sym) + eff_wt_ridge * diag(n_lags), t(C_sym))
              W_sym <- W_1 + neumann_wt * W_1 %*% (diag(n_lags) - C_sym %*% W_1)
              R_mat <- W_sym %*% C_sym
              diag_base <- pmax(diag(R_mat), gain_floor)
              if (is.finite(dt_lag) && dt_lag > 0 && dt_lag < dt_ref) {
                sigma_lobe <- 0.72 * sqrt(max(0.0, 1.0 - (dt_lag / dt_ref)^2))
                lag_dist <- abs(outer(lags, lags, "-"))
                L_lobe <- exp(-0.5 * (lag_dist / max(sigma_lobe, 1e-6))^2)
                lobe_ratios <- rowSums(R_mat * L_lobe) / diag_base
                mid_idx <- if (n_lags >= 3L) seq.int(2L, n_lags - 1L) else seq_len(n_lags)
                lobe_scale <- max(1.0, stats::median(lobe_ratios[mid_idx], na.rm = TRUE))
                diag_gain <- diag_base * lobe_scale
              } else {
                diag_gain <- diag_base
              }
              W_lag <- sweep(W_sym, 1L, diag_gain, "/")

              idx_common <- which(ok_common)
              Y_common_old <- lapply(seq_len(n_lags), function(li) {
                demixed_epochs[[li]][row_span[idx_common], voxels, drop = FALSE]
              })
              B_common_old <- lapply(seq_len(n_lags), function(li) {
                ok_sub_li <- ok_common[ok_tau_list[[li]]]
                B_demix_list[[li]][ok_sub_li, , drop = FALSE]
              })
              n_obs_b <- ncol(B_common_old[[1L]])
              for (l1 in seq_len(n_lags)) {
                Y_new_l1 <- matrix(0, m_common, length(voxels))
                B_new_l1 <- matrix(0, m_common, n_obs_b)
                for (l2 in seq_len(n_lags)) {
                  w_12 <- W_lag[l1, l2]
                  if (abs(w_12) > 0) {
                    Y_new_l1 <- Y_new_l1 + w_12 * Y_common_old[[l2]]
                    B_new_l1 <- B_new_l1 + w_12 * B_common_old[[l2]]
                  }
                }
                demixed_epochs[[l1]][row_span[idx_common], voxels] <- Y_new_l1
                ok_sub_l1 <- ok_common[ok_tau_list[[l1]]]
                B_demix_list[[l1]][ok_sub_l1, ] <- B_new_l1
              }
            }
          }
        }
      }

      stratum_specs_r[[si_r]] <- list(
        kernel = hrf_k,
        ar = ar,
        lambda = lam_k,
        mix = mix_spec,
        ok_tau_list = ok_tau_list,
        B_demix_list = B_demix_list
      )

      U_demix_list <- vector("list", n_lags)
      for (li in seq_len(n_lags)) {
        ok_tau <- ok_tau_list[[li]]
        if (is.null(ok_tau) || !any(ok_tau)) next
        m_ok <- sum(ok_tau)
        B_demix_tau <- B_demix_list[[li]]
        M_act_tau <- M_act_list[[li]]

        # Compute dimensionless trial innovation factor and error covariance for stratum (r, key, li)
        U_demix <- compute_ar1_operator_factor(
          B_demix_tau, rows, ar, center_obs = center_obs,
          obs_weights = fw_r, weight_space = w_space
        )
        U_demix_list[[li]] <- U_demix
        V_demix <- crossprod(U_demix)
        if (m_ok > 1L) {
          Q_tau <- diag(m_ok) - matrix(1 / m_ok, m_ok, m_ok)
          V_demix_cent <- (Q_tau %*% V_demix %*% Q_tau) / (1 - 1 / m_ok)
        } else {
          Q_tau <- matrix(1.0, 1L, 1L)
          V_demix_cent <- V_demix
        }
        if (m_ok <= 1L) {
          M_post <- matrix(1.0, m_ok, m_ok)
        } else if (!isTRUE(demix_target_trials) || isTRUE(demix_within_trial)) {
          # Under multi-lag within-trial trajectory de-mixing (demix_within_trial = TRUE),
          # W_lag combines readouts across multiple lags for a within-trial waveform rather
          # than a single-lag Dirac source, and keeps X localized to trial i so off-peak
          # boundary lags cannot couple trial i+1 non-target events into trial i's regressor.
          M_post <- Q_tau
        } else {
          M_post <- Q_tau %*% (M_act_tau %*% DG_tau_list[[li]][, ok_tau, drop = FALSE]) %*% Q_tau
        }
        if (!is.null(U_ref_by_key) && li == li_ref_cor && m_ok >= 2L) {
          U_c_ref <- sweep(U_demix, 1L, rowMeans(U_demix), "-")
          cn_ref <- sqrt(pmax(colSums(U_c_ref^2), 1e-24))
          U_ref_by_key[[key]] <- list(ok_tau = ok_tau, U_norm = sweep(U_c_ref, 2L, cn_ref, "/"))
        }
        d_demix_unit <- sqrt(pmax(0, diag(V_demix_cent)))
        se_demixed_epochs[[li]][row_span[ok_tau], voxels] <- outer(d_demix_unit, sig_eff_r[voxels])
        if (m_ok > 1L) {
          num_ac1 <- V_demix_cent[cbind(2:m_ok, 1:(m_ok - 1L))]
          den_ac1 <- pmax(d_demix_unit[2:m_ok] * d_demix_unit[1:(m_ok - 1L)], .Machine$double.eps)
          ac1_vec <- c(0, pmin(0.999, pmax(-0.999, num_ac1 / den_ac1)))
        } else {
          ac1_vec <- 0
        }
        error_ac1_epochs[[li]][row_span[ok_tau], voxels] <- matrix(ac1_vec, nrow = m_ok, ncol = length(voxels))
        V_demix_blocks[[li]][[r]][[key]] <- list(
          ok_tau = ok_tau,
          V_raw = V_demix,
          V_cent = V_demix_cent,
          M_post = M_post
        )
      }

      # Exact L x L cross-lag trial error covariance across lags in a single BLAS crossprod
      if (n_lags >= 2L && m_common >= 2L && !any(vapply(U_demix_list, is.null, logical(1L)))) {
        n_obs_u <- nrow(U_demix_list[[1L]])
        U_lag_mat <- matrix(0, n_obs_u * m_common, n_lags)
        for (li in seq_len(n_lags)) {
          ok_sub_li <- ok_common[ok_tau_list[[li]]]
          U_sub <- U_demix_list[[li]][, ok_sub_li, drop = FALSE]
          U_cent <- sweep(U_sub, 1L, rowMeans(U_sub), "-")
          U_lag_mat[, li] <- as.numeric(U_cent)
        }
        V_lag_unit <- crossprod(U_lag_mat) / max(1L, m_common - 1L)
        for (vi in voxels) {
          lag_cov_accum[[vi]] <- lag_cov_accum[[vi]] + m_common * (sig_eff_r[vi]^2) * V_lag_unit
          lag_cov_trials[vi] <- lag_cov_trials[vi] + m_common
        }
      }
    }

    # Accumulate cross-unit trial error correlation across runs from task-free innovations
    if (!is.null(unit_err_cor_accum) && !is.null(Z_unit_r)) {
      R_innov_r <- crossprod(Z_unit_r)
      R_innov_r <- 0.5 * (R_innov_r + t(R_innov_r))
      diag(R_innov_r) <- 1.0
      u_keys_r <- unique(keys[valid])
      if (length(u_keys_r) >= 2L && !is.null(U_ref_by_key) && length(U_ref_by_key) >= 2L) {
        O_unit_r <- matrix(1.0, n_units, n_units)
        for (ia in seq_along(u_keys_r)) {
          ka <- u_keys_r[ia]
          ua <- U_ref_by_key[[ka]]
          if (is.null(ua)) next
          va <- which(keys == ka & valid)
          for (ib in seq_len(ia - 1L)) {
            kb <- u_keys_r[ib]
            ub <- U_ref_by_key[[kb]]
            if (is.null(ub)) next
            vb <- which(keys == kb & valid)
            ok_ab <- ua$ok_tau & ub$ok_tau
            if (any(ok_ab)) {
              pa <- ok_ab[ua$ok_tau]
              pb <- ok_ab[ub$ok_tau]
              ov_ab <- mean(colSums(ua$U_norm[, pa, drop = FALSE] * ub$U_norm[, pb, drop = FALSE]))
              O_unit_r[va, vb] <- ov_ab
              O_unit_r[vb, va] <- ov_ab
            }
          }
        }
        R_innov_r <- R_innov_r * O_unit_r
        diag(R_innov_r) <- 1.0
      }
      w_df_r <- max(10.0, length(rows) - m_all - 2.0)
      unit_err_cor_accum <- unit_err_cor_accum + w_df_r * R_innov_r
      unit_err_df_tot <- unit_err_df_tot + w_df_r
    }

    if (length(u_skeys_r) > 0L) {
      s_idx_u <- match(keys, u_skeys_r)
      s_idx_u[!valid] <- NA_integer_
      run_stratum_ops[[r]] <- list(
        rows = rows,
        fw_r = fw_r,
        u_skeys = u_skeys_r,
        s_idx_by_unit = s_idx_u,
        stratum_specs = stratum_specs_r
      )
    }
  }

    if (need_crossrun_haus && haus_pass == 1L) {
      gamma_haus_override <- numeric(n_units)
      for (u_h in seq_len(n_units)) {
        W_uh <- haus_W[[u_h]]
        act_kh <- which(diag(W_uh) > 1e-12)
        if (length(act_kh) > 0L) {
          q_uh <- length(act_kh)
          W_sub <- W_uh[act_kh, act_kh, drop = FALSE]
          reg_wh <- 1e-10 * max(mean(diag(W_sub)), 1e-12)
          b_pool <- tryCatch(
            solve(W_sub + reg_wh * diag(q_uh), haus_rhs[[u_h]][act_kh]),
            error = function(e) NULL
          )
          if (!is.null(b_pool) && all(is.finite(b_pool))) {
            f_pool <- as.numeric(crossprod(b_pool, haus_rhs[[u_h]][act_kh])) / q_uh
            if (is.finite(f_pool) && f_pool >= 0) {
              p_pool <- stats::pf(f_pool, df1 = q_uh, df2 = max(8.0, haus_df2[u_h]), lower.tail = FALSE)
              if (is.finite(p_pool) && p_pool < 0.15) {
                gamma_haus_override[u_h] <- .poet_firm_gate(
                  -log10(max(p_pool, 1e-16)),
                  lo = -log10(0.15),
                  hi = -log10(0.035)
                )
              }
            }
          }
        }
      }
      diff_p1 <- sweep(haus_gam_p1, 2L, gamma_haus_override, "-")
      if (all(is.na(diff_p1) | abs(diff_p1) < 1e-6)) {
        break
      }
    }
  }

  unit_error_cor <- NULL
  shared_loading_vec <- rep(0.0, n_units)
  names(shared_loading_vec) <- unit_ids
  if (!is.null(unit_err_cor_accum) && unit_err_df_tot > 0) {
    R_u_raw <- unit_err_cor_accum / unit_err_df_tot
    R_u_raw <- 0.5 * (R_u_raw + t(R_u_raw))
    df_u <- max(20.0, unit_err_df_tot)
    R_u_eb <- matrix(0.0, n_units, n_units, dimnames = list(unit_ids, unit_ids))
    off_mask <- row(R_u_raw) != col(R_u_raw)
    r_off <- pmin(0.995, pmax(-0.995, R_u_raw[off_mask]))
    keep_off <- is.finite(r_off) & (abs(r_off) > 0.25 / sqrt(df_u))
    if (any(keep_off)) {
      r_k <- r_off[keep_off]
      var_rk <- ((1.0 - r_k^2)^2) / df_u
      r_sh <- r_k * (r_k^2 / pmax(r_k^2 + var_rk, 1e-12))
      r_off_new <- numeric(length(r_off))
      r_off_new[keep_off] <- pmin(0.98, pmax(-0.98, r_sh))
      R_u_eb[off_mask] <- r_off_new
    }
    diag(R_u_eb) <- 1.0
    eig_ue <- eigen(R_u_eb, symmetric = TRUE)
    R_u_psd <- eig_ue$vectors %*% (pmax(eig_ue$values, 1e-4) * t(eig_ue$vectors))
    du_psd <- sqrt(pmax(diag(R_u_psd), .Machine$double.eps))
    unit_error_cor <- R_u_psd / outer(du_psd, du_psd)
    unit_error_cor <- 0.5 * (unit_error_cor + t(unit_error_cor))
    diag(unit_error_cor) <- 1.0
    dimnames(unit_error_cor) <- list(unit_ids, unit_ids)

    max_off <- max(abs(unit_error_cor[off_mask]))
    if (is.finite(max_off) && max_off > 1e-6) {
      if (n_units == 2L) {
        r12 <- unit_error_cor[1L, 2L]
        if (r12 >= 0) {
          shared_loading_vec[] <- c(sqrt(r12), sqrt(r12))
        } else {
          shared_loading_vec[] <- c(sqrt(-r12), -sqrt(-r12))
        }
      } else {
        R_pa <- unit_error_cor
        diag(R_pa) <- apply(abs(unit_error_cor - diag(n_units)), 1L, max)
        for (it_pa in seq_len(4L)) {
          eg_pa <- eigen(0.5 * (R_pa + t(R_pa)), symmetric = TRUE)
          lv_pa <- eg_pa$vectors[, 1L] * sqrt(max(0.0, eg_pa$values[1L]))
          if (sum(lv_pa) < 0) lv_pa <- -lv_pa
          lv_pa <- pmin(0.995, pmax(-0.995, lv_pa))
          diag(R_pa) <- lv_pa^2
        }
        shared_loading_vec[] <- lv_pa
      }
    }
  }

  has_gain <- gain_count > 0L
  gain_by_lag_voxel[has_gain] <- gain_accum[has_gain] / gain_count[has_gain]

  # Center covariates within run for pooled regression
  X_pool <- do.call(rbind, x_centered_list)
  x_row_finite <- if (n_terms > 0L) (rowSums(!is.finite(X_pool)) == 0L) else rep(TRUE, total_trials)

  need_epoch_graph <- !identical(tau2_shrinkage, "none") &&
    !identical(tau2_shrinkage, "eb") &&
    n_units >= 2L &&
    (!isTRUE(aggregate_by_parcel) || !is.null(spatial_graph) || identical(tau2_shrinkage, "graph_eb"))
  epoch_graph <- if (need_epoch_graph) {
    resolve_epoch_spatial_graph(
      deconvolution, fit, spatial_graph, aggregate_by_parcel, n_units, unit_ids,
      parcel_vox_list = parcel_vox_list, parcel_weights_list = parcel_weights_list
    )
  } else {
    NULL
  }
  eff_tau2_shrinkage <- if (identical(tau2_shrinkage, "none") || n_units < 2L) {
    "none"
  } else if (identical(tau2_shrinkage, "auto")) {
    if ((!isTRUE(aggregate_by_parcel) || !is.null(spatial_graph)) &&
        !is.null(epoch_graph) && n_units >= 6L) {
      "graph_eb"
    } else {
      "none"
    }
  } else if (identical(tau2_shrinkage, "graph_eb")) {
    if (is.null(epoch_graph)) {
      stop(
        "tau2_shrinkage = 'graph_eb' requires spatial coordinates, voxel_map, or spatial_graph ",
        "(or use tau2_shrinkage = 'eb' for unstructured non-spatial units).",
        call. = FALSE
      )
    }
    "graph_eb"
  } else {
    "eb"
  }

  task_gate_vec <- numeric(n_units)
  eb_group_ids <- integer(0L)
  eb_cross_w <- 1.0
  if (identical(eff_tau2_shrinkage, "graph_eb")) {
    if (length(fit$selections) == n_runs) {
      r2_mat_eb <- matrix(0.0, n_runs, n_units)
      has_r2 <- FALSE
      for (r in seq_len(n_runs)) {
        sel_r <- fit$selections[[r]]
        sc_raw <- if (!is.null(sel_r[["selected_adjusted_partial_r2"]])) {
          sel_r[["selected_adjusted_partial_r2"]]
        } else if (!is.null(sel_r[["best_adjusted_partial_r2"]])) {
          sel_r[["best_adjusted_partial_r2"]]
        } else {
          sel_r[["score"]]
        }
        if (!is.null(sc_raw) && is.numeric(sc_raw)) {
          sv <- as.numeric(sc_raw)
          sv[!is.finite(sv) | sv < 0] <- 0.0
          if (!isTRUE(aggregate_by_parcel) && length(sv) == n_units) {
            r2_mat_eb[r, ] <- sv
            has_r2 <- TRUE
          } else if (isTRUE(aggregate_by_parcel) && length(sv) == n_voxels &&
                     !is.null(parcel_vox_list) && length(parcel_vox_list) == n_units) {
            for (g in seq_len(n_units)) {
              vg <- parcel_vox_list[[g]]
              wg <- parcel_weights_list[[g]]
              r2_mat_eb[r, g] <- sum(wg * sv[vg]) / max(sum(wg), 1e-12)
            }
            has_r2 <- TRUE
          }
        }
      }
      if (has_r2) {
        r2_mean_eb <- colMeans(r2_mat_eb)
        pos_r2 <- r2_mean_eb[r2_mean_eb > 1e-4]
        r2_half <- if (length(pos_r2)) max(0.02, stats::median(pos_r2)) else 0.05
        task_gate_vec <- r2_mean_eb / (r2_mean_eb + r2_half)
      } else {
        task_gate_vec <- rep(0.5, n_units)
      }
    } else {
      task_gate_vec <- rep(0.5, n_units)
    }
    if (!isTRUE(aggregate_by_parcel) && !is.null(fit$specification$spatial_groups) &&
        !identical(epoch_graph$parcel_boundary, "ignore")) {
      eb_group_ids <- hrf_encode_group_ids(fit$specification$spatial_groups, n_units)
      eb_cross_w <- if (identical(epoch_graph$parcel_boundary, "hard")) 0.0 else 0.35
    }
  }

  n_coef_alloc <- n_lags * n_units * max(1L, n_terms)
  cf_voxel <- character(n_coef_alloc)
  cf_lag <- rep(NA_real_, n_coef_alloc)
  cf_term <- rep(NA_character_, n_coef_alloc)
  cf_beta_raw <- rep(NA_real_, n_coef_alloc)
  cf_beta_corr <- rep(NA_real_, n_coef_alloc)
  cf_beta_unw <- rep(NA_real_, n_coef_alloc)
  cf_beta_pw <- rep(NA_real_, n_coef_alloc)
  cf_beta_gls <- rep(NA_real_, n_coef_alloc)
  cf_se_corr <- rep(NA_real_, n_coef_alloc)
  cf_se_unw <- rep(NA_real_, n_coef_alloc)
  cf_se_pw <- rep(NA_real_, n_coef_alloc)
  cf_se_gls <- rep(NA_real_, n_coef_alloc)
  cf_r_raw <- rep(NA_real_, n_coef_alloc)
  cf_r_dem <- rep(NA_real_, n_coef_alloc)
  cf_t_stat <- rep(NA_real_, n_coef_alloc)
  cf_p_val <- rep(NA_real_, n_coef_alloc)
  cf_t_unw <- rep(NA_real_, n_coef_alloc)
  cf_t_pw <- rep(NA_real_, n_coef_alloc)
  cf_t_gls <- rep(NA_real_, n_coef_alloc)
  cf_df <- rep(NA_real_, n_coef_alloc)
  cf_tau2 <- rep(NA_real_, n_coef_alloc)
  cf_gain <- rep(NA_real_, n_coef_alloc)
  cf_mean_se <- rep(NA_real_, n_coef_alloc)
  cf_mean_ac1 <- rep(NA_real_, n_coef_alloc)
  row_k <- 0L

  lag_key_names <- sprintf("lag_%g", lags)
  run_trial_operators <- if (isTRUE(return_panel)) {
    stats::setNames(
      lapply(unit_ids, function(uid) stats::setNames(vector("list", n_lags), lag_key_names)),
      unit_ids
    )
  } else NULL

  for (li in seq_len(n_lags)) {
    Y_raw_lag <- raw_epochs[[li]]
    Y_demix_lag <- demixed_epochs[[li]]
    # Center epoch readouts within each run to eliminate run baseline offsets
    for (r in seq_len(n_runs)) {
      rs <- starts[r]:ends[r]
      Y_raw_lag[rs, ] <- sweep(Y_raw_lag[rs, , drop = FALSE], 2L,
        colMeans(Y_raw_lag[rs, , drop = FALSE], na.rm = TRUE), "-")
      Y_demix_lag[rs, ] <- sweep(Y_demix_lag[rs, , drop = FALSE], 2L,
        colMeans(Y_demix_lag[rs, , drop = FALSE], na.rm = TRUE), "-")
    }

    # Stratum-level unit-variance covariance/design cache per run at this lag
    tpl_cache <- lapply(seq_len(n_runs), function(r) new.env(parent = emptyenv()))
    get_stratum_tpl <- function(r, key_rv, ok_r_local, rs_r, m_target_r) {
      if (is.null(key_rv) || !nzchar(key_rv)) return(NULL)
      blk_r <- V_demix_blocks[[li]][[r]][[key_rv]]
      if (is.null(blk_r)) return(NULL)
      m_r <- length(ok_r_local)
      ckey <- if (m_r == m_target_r) {
        key_rv
      } else {
        paste0(key_rv, "|", paste(ok_r_local, collapse = ","))
      }
      if (exists(ckey, envir = tpl_cache[[r]], inherits = FALSE)) {
        return(get(ckey, envir = tpl_cache[[r]], inherits = FALSE))
      }
      ok_tau_idx <- which(blk_r$ok_tau)
      pos_in_blk <- match(ok_r_local, ok_tau_idx)
      if (anyNA(pos_in_blk)) return(NULL)

      Q_mr <- if (m_r > 1L) (diag(m_r) - matrix(1.0 / m_r, m_r, m_r)) else matrix(0.0, 1L, 1L)
      if (!is.null(blk_r$V_raw)) {
        V_raw_u <- blk_r$V_raw[pos_in_blk, pos_in_blk, drop = FALSE]
        Omega_cent_u <- if (m_r > 1L) (Q_mr %*% V_raw_u %*% Q_mr) else matrix(0.0, 1L, 1L)
        Omega_u <- if (m_r > 1L) Omega_cent_u / (1.0 - 1.0 / m_r) else V_raw_u
      } else {
        Omega_u <- blk_r$V_cent[pos_in_blk, pos_in_blk, drop = FALSE]
        Omega_cent_u <- if (m_r > 1L) ((m_r - 1.0) / m_r) * (Q_mr %*% Omega_u %*% Q_mr) else matrix(0.0, 1L, 1L)
      }
      if (!is.null(blk_r$M_post)) {
        if (length(pos_in_blk) == nrow(blk_r$M_post)) {
          M_post_r <- blk_r$M_post
        } else if (m_r > 1L) {
          M_post_r <- Q_mr %*% blk_r$M_post[pos_in_blk, pos_in_blk, drop = FALSE] %*% Q_mr
        } else {
          M_post_r <- matrix(0.0, 1L, 1L)
        }
      } else {
        M_post_r <- Q_mr
      }
      M2_r <- tcrossprod(M_post_r)
      tr_Om_u <- sum(diag(Omega_cent_u))
      normF2_Om_u <- sum(Omega_cent_u * Omega_cent_u)
      tr_M2 <- sum(diag(M2_r))

      exp_diff_Om_u <- 0.0
      exp_diff_M2 <- 0.0
      if (m_r >= 2L) {
        diag_Om <- diag(Omega_cent_u)
        off1_Om <- Omega_cent_u[cbind(2:m_r, 1:(m_r - 1L))]
        exp_diff_Om_u <- sum(diag_Om[-1L] + diag_Om[-m_r] - 2.0 * off1_Om)
        diag_M2 <- diag(M2_r)
        off1_M2 <- M2_r[cbind(2:m_r, 1:(m_r - 1L))]
        exp_diff_M2 <- sum(diag_M2[-1L] + diag_M2[-m_r] - 2.0 * off1_M2)
      }

      Xr_tilde_sub <- NULL
      xOx_u <- NULL
      xM2x <- NULL
      if (n_terms > 0L) {
        Xr_raw_sub <- X_raw_pool[rs_r, , drop = FALSE]
        Xr_cent_sub <- if (m_r > 1L) {
          sweep(Xr_raw_sub, 2L, colMeans(Xr_raw_sub), "-")
        } else {
          matrix(0.0, 1L, n_terms)
        }
        Xr_tilde_sub <- M_post_r %*% Xr_cent_sub
        xOx_u <- crossprod(Xr_tilde_sub, Omega_cent_u %*% Xr_tilde_sub)
        xM2x <- crossprod(Xr_tilde_sub, M2_r %*% Xr_tilde_sub)
      }

      tpl <- list(
        ckey = ckey,
        m_r = m_r,
        rs_r = rs_r,
        ok_r_local = ok_r_local,
        Omega_u = Omega_u,
        Omega_cent_u = Omega_cent_u,
        M_post_r = M_post_r,
        M2_r = M2_r,
        tr_Om_u = tr_Om_u,
        normF2_Om_u = normF2_Om_u,
        tr_M2 = tr_M2,
        exp_diff_Om_u = exp_diff_Om_u,
        exp_diff_M2 = exp_diff_M2,
        Xr_tilde_sub = Xr_tilde_sub,
        xOx_u = xOx_u,
        xM2x = xM2x
      )
      assign(ckey, tpl, envir = tpl_cache[[r]])
      tpl
    }

    # Pass-A storage across units at lag li
    u_usable <- logical(n_units)
    u_n_ok <- integer(n_units)
    u_df_v <- rep(NA_real_, n_units)
    u_df_tau2 <- rep(NA_real_, n_units)
    u_d_eff_tau2 <- rep(NA_real_, n_units)
    u_s2_med <- rep(1.0, n_units)
    u_s2_obs_scale <- rep(1.0, n_units)
    u_theta_raw <- rep(NA_real_, n_units)
    u_tau2_unshrunk <- rep(NA_real_, n_units)
    u_mean_se <- rep(NA_real_, n_units)
    u_mean_ac1 <- rep(0.0, n_units)
    ckey_mat <- matrix("", n_runs, n_units)
    fb_blocks <- vector("list", n_units)

    u_beta_raw <- if (n_terms > 0L) matrix(NA_real_, n_terms, n_units) else NULL
    u_beta_unw <- if (n_terms > 0L) matrix(NA_real_, n_terms, n_units) else NULL
    u_ss_resid_unw <- if (n_terms > 0L) rep(NA_real_, n_units) else NULL
    u_xtx_inv_tilde <- if (n_terms > 0L) vector("list", n_units) else NULL
    U_resid_mat <- if (!identical(eff_tau2_shrinkage, "none")) {
      matrix(0.0, total_trials, n_units)
    } else {
      NULL
    }

    # Pass A: Per-unit raw method-of-moments tau2_raw, OLS slopes, and unit residuals
    for (vi in seq_len(n_units)) {
      y_raw_v <- Y_raw_lag[, vi]
      y_dem_v <- Y_demix_lag[, vi]
      s_raw_v <- se_raw_epochs[[li]][, vi]
      s_dem_v <- se_demixed_epochs[[li]][, vi]
      s2_raw_v <- s_raw_v * s_raw_v
      s2_dem_v <- s_dem_v * s_dem_v

      ok_s_dem <- is.finite(s2_dem_v) & (s2_dem_v > 0)
      if (any(ok_s_dem)) {
        floor_op <- 0.10 * stats::median(s2_dem_v[ok_s_dem])
        w_op_vec <- 1 / (floor_op + s2_dem_v[ok_s_dem])
        pw_op_epochs[[li]][ok_s_dem, vi] <- w_op_vec / mean(w_op_vec)
      }
      ok_s_raw <- is.finite(s2_raw_v) & (s2_raw_v > 0)
      if (any(ok_s_raw)) {
        floor_raw <- 0.10 * stats::median(s2_raw_v[ok_s_raw])
        w_raw_vec <- 1 / (floor_raw + s2_raw_v[ok_s_raw])
        pw_raw_epochs[[li]][ok_s_raw, vi] <- w_raw_vec / mean(w_raw_vec)
      }

      ok_v <- if (n_terms > 0L) {
        is.finite(y_raw_v) & is.finite(y_dem_v) & x_row_finite
      } else {
        is.finite(y_raw_v) & is.finite(y_dem_v)
      }
      n_ok <- sum(ok_v)
      u_n_ok[vi] <- n_ok
      mean_se_v <- if (any(ok_s_dem)) mean(s_dem_v[ok_s_dem]) else NA_real_
      ac1_v <- error_ac1_epochs[[li]][, vi]
      mean_ac1_v <- if (any(is.finite(ac1_v))) mean(ac1_v[is.finite(ac1_v) & ac1_v != 0]) else NA_real_
      if (!is.finite(mean_ac1_v)) mean_ac1_v <- 0
      u_mean_se[vi] <- mean_se_v
      u_mean_ac1[vi] <- mean_ac1_v

      if (n_ok <= n_terms + 2L) {
        if (any(ok_s_dem)) {
          pw_epochs[[li]][ok_s_dem, vi] <- pw_op_epochs[[li]][ok_s_dem, vi]
          fl_c <- 0.10 * stats::median(s2_dem_v[ok_s_dem])
          pw_cohort_epochs[[li]][ok_s_dem, vi] <- 1 / (fl_c + s2_dem_v[ok_s_dem])
        }
        next
      }

      u_usable[vi] <- TRUE
      X_tilde_pool <- if (n_terms > 0L) matrix(NA_real_, total_trials, n_terms) else NULL
      xOx_sum <- if (n_terms > 0L) matrix(0, n_terms, n_terms) else NULL
      xM2x_sum <- if (n_terms > 0L) matrix(0, n_terms, n_terms) else NULL
      tr_Omega_cent <- 0.0
      normF2_Omega_cent <- 0.0
      tr_M2_marg <- 0.0
      ss_marg <- 0.0
      df_marg <- 0L
      ss_diff <- 0.0
      exp_diff_Om <- 0.0
      exp_diff_M2 <- 0.0
      n_runs_active <- 0L
      u_marg_vec <- if (!is.null(U_resid_mat) && identical(eff_tau2_mode, "marginal")) {
        numeric(total_trials)
      } else {
        NULL
      }

      for (r in seq_len(n_runs)) {
        rs <- starts[r]:ends[r]
        m_target_r <- length(rs)
        ok_r_local <- which(ok_v[rs])
        m_r <- length(ok_r_local)
        if (!m_r) next
        n_runs_active <- n_runs_active + 1L
        rs_r <- rs[ok_r_local]
        key_rv <- if (!is.null(stratum_map[[r]])) stratum_map[[r]][vi] else NULL
        tpl <- get_stratum_tpl(r, key_rv, ok_r_local, rs_r, m_target_r)

        if (!is.null(tpl)) {
          ckey_mat[r, vi] <- tpl$ckey
          s2_eff_rv <- sig_eff_mat[r, vi]^2
          tr_Omega_cent <- tr_Omega_cent + s2_eff_rv * tpl$tr_Om_u
          normF2_Omega_cent <- normF2_Omega_cent + (s2_eff_rv * s2_eff_rv) * tpl$normF2_Om_u
          tr_M2_marg <- tr_M2_marg + tpl$tr_M2
          if (m_r >= 2L) {
            yr_dem_r <- demixed_epochs[[li]][rs_r, vi]
            yr_c_r <- yr_dem_r - mean(yr_dem_r)
            if (!is.null(u_marg_vec)) u_marg_vec[rs_r] <- yr_c_r
            ss_marg <- ss_marg + sum(yr_c_r^2)
            df_marg <- df_marg + (m_r - 1L)
            dy_r <- diff(yr_c_r)
            ss_diff <- ss_diff + sum(dy_r^2)
            exp_diff_Om <- exp_diff_Om + s2_eff_rv * tpl$exp_diff_Om_u
            exp_diff_M2 <- exp_diff_M2 + tpl$exp_diff_M2
          }
          if (n_terms > 0L) {
            X_tilde_pool[rs_r, ] <- tpl$Xr_tilde_sub
            xOx_sum <- xOx_sum + s2_eff_rv * tpl$xOx_u
            xM2x_sum <- xM2x_sum + tpl$xM2x
          }
        } else {
          ckey_mat[r, vi] <- "__FB__"
          Q_mr <- if (m_r > 1L) (diag(m_r) - matrix(1.0 / m_r, m_r, m_r)) else matrix(0.0, 1L, 1L)
          s2_sub <- s2_dem_v[rs_r]
          s2_sub[!is.finite(s2_sub)] <- 1
          Omega_r <- diag(s2_sub, m_r)
          Omega_cent_r <- if (m_r > 1L) (Q_mr %*% Omega_r %*% Q_mr) else matrix(0.0, 1L, 1L)
          M_post_r <- Q_mr
          M2_r <- tcrossprod(M_post_r)

          tr_Omega_cent <- tr_Omega_cent + sum(diag(Omega_cent_r))
          normF2_Omega_cent <- normF2_Omega_cent + sum(Omega_cent_r * Omega_cent_r)
          tr_M2_marg <- tr_M2_marg + sum(diag(M2_r))
          if (m_r >= 2L) {
            yr_dem_r <- demixed_epochs[[li]][rs_r, vi]
            yr_c_r <- yr_dem_r - mean(yr_dem_r)
            if (!is.null(u_marg_vec)) u_marg_vec[rs_r] <- yr_c_r
            ss_marg <- ss_marg + sum(yr_c_r^2)
            df_marg <- df_marg + (m_r - 1L)
            dy_r <- diff(yr_c_r)
            ss_diff <- ss_diff + sum(dy_r^2)
            diag_Om <- diag(Omega_cent_r)
            off1_Om <- Omega_cent_r[cbind(2:m_r, 1:(m_r - 1L))]
            exp_diff_Om <- exp_diff_Om + sum(diag_Om[-1L] + diag_Om[-m_r] - 2.0 * off1_Om)
            diag_M2 <- diag(M2_r)
            off1_M2 <- M2_r[cbind(2:m_r, 1:(m_r - 1L))]
            exp_diff_M2 <- exp_diff_M2 + sum(diag_M2[-1L] + diag_M2[-m_r] - 2.0 * off1_M2)
          }
          Xr_tilde_sub <- NULL
          if (n_terms > 0L) {
            Xr_raw_sub <- X_raw_pool[rs_r, , drop = FALSE]
            Xr_cent_sub <- if (m_r > 1L) {
              sweep(Xr_raw_sub, 2L, colMeans(Xr_raw_sub), "-")
            } else {
              matrix(0.0, 1L, n_terms)
            }
            Xr_tilde_sub <- M_post_r %*% Xr_cent_sub
            X_tilde_pool[rs_r, ] <- Xr_tilde_sub
            xOx_sum <- xOx_sum + crossprod(Xr_tilde_sub, Omega_cent_r %*% Xr_tilde_sub)
            xM2x_sum <- xM2x_sum + crossprod(Xr_tilde_sub, M2_r %*% Xr_tilde_sub)
          }
          if (is.null(fb_blocks[[vi]])) fb_blocks[[vi]] <- vector("list", n_runs)
          fb_blocks[[vi]][[r]] <- list(
            m_r = m_r,
            rs_r = rs_r,
            ok_r_local = ok_r_local,
            Omega_r = Omega_r,
            Omega_cent_r = Omega_cent_r,
            M_post_r = M_post_r,
            M2_r = M2_r,
            Xr_tilde_sub = Xr_tilde_sub
          )
        }
      }

      s2_ok <- s2_dem_v[ok_v]
      s2_med <- if (any(is.finite(s2_ok) & s2_ok > 0)) stats::median(s2_ok[is.finite(s2_ok) & s2_ok > 0]) else 1.0
      u_s2_med[vi] <- s2_med
      tau2_var_raw <- if (tr_M2_marg > 1e-6) (ss_marg - tr_Omega_cent) / tr_M2_marg else 0.0
      tau2_diff_raw <- if (exp_diff_M2 > 1e-6) (ss_diff - exp_diff_Om) / exp_diff_M2 else 0.0
      mad_y <- stats::mad(y_dem_v[ok_v], constant = 1.4826, na.rm = TRUE)
      s2_cent_med <- s2_med * max(0.5, df_marg / max(1, n_ok))
      m2_cent_scale <- max(0.5, tr_M2_marg / max(1, n_ok))
      tau2_mad_raw <- if (is.finite(mad_y)) (mad_y^2 - s2_cent_med) / m2_cent_scale else 0.0
      tau2_marginal_raw <- 0.45 * tau2_var_raw + 0.35 * tau2_diff_raw + 0.20 * tau2_mad_raw
      tau2_marginal <- max(0.0, tau2_marginal_raw)
      s2_obs_marg <- if (tr_Omega_cent > 1e-12 && tr_M2_marg > 1e-6) {
        tr_Omega_cent / tr_M2_marg
      } else {
        s2_med
      }

      if (n_terms > 0L) {
        X_ok <- X_pool[ok_v, , drop = FALSE]
        X_tilde_ok <- X_tilde_pool[ok_v, , drop = FALSE]
        xtx_raw <- crossprod(X_ok) + 1e-8 * diag(n_terms)
        xtx_inv_raw <- solve(xtx_raw)
        beta_raw_v <- as.numeric(xtx_inv_raw %*% crossprod(X_ok, y_raw_v[ok_v]))

        xtx_tilde <- crossprod(X_tilde_ok) + 1e-8 * diag(n_terms)
        xtx_inv_tilde <- solve(xtx_tilde)
        beta_unw_v <- as.numeric(xtx_inv_tilde %*% crossprod(X_tilde_ok, y_dem_v[ok_v]))

        resid_unw <- y_dem_v[ok_v] - as.numeric(X_tilde_ok %*% beta_unw_v)
        df_v <- max(1L, n_ok - n_terms - max(1L, n_runs_active))

        tr_proj_Omega <- max(0.0, tr_Omega_cent - sum(diag(xtx_inv_tilde %*% xOx_sum)))
        denom_resid <- max(1.0, tr_M2_marg - sum(diag(xtx_inv_tilde %*% xM2x_sum)))
        tau2_resid_raw <- (sum(resid_unw^2) - tr_proj_Omega) / denom_resid
        tau2_resid <- max(0.0, tau2_resid_raw)
        s2_obs_resid <- if (tr_proj_Omega > 1e-12) (tr_proj_Omega / denom_resid) else s2_obs_marg

        if (identical(eff_tau2_mode, "marginal")) {
          tau2_raw_untrunc <- tau2_marginal_raw
          tau2_bio_unshrunk <- tau2_marginal
          s2_obs_scale <- s2_obs_marg
        } else if (identical(eff_tau2_mode, "residual")) {
          tau2_raw_untrunc <- tau2_resid_raw
          tau2_bio_unshrunk <- tau2_resid
          s2_obs_scale <- s2_obs_resid
        } else {
          tau2_raw_untrunc <- 0.65 * tau2_resid_raw + 0.35 * tau2_marginal_raw
          tau2_bio_unshrunk <- max(0.0, tau2_raw_untrunc)
          s2_obs_scale <- 0.65 * s2_obs_resid + 0.35 * s2_obs_marg
        }

        u_beta_raw[, vi] <- beta_raw_v
        u_beta_unw[, vi] <- beta_unw_v
        u_ss_resid_unw[vi] <- sum(resid_unw^2)
        u_xtx_inv_tilde[[vi]] <- xtx_inv_tilde

        if (!is.null(U_resid_mat) && !identical(eff_tau2_mode, "marginal")) {
          u_res_full <- numeric(total_trials)
          u_res_full[ok_v] <- resid_unw
          for (r in seq_len(n_runs)) {
            rs <- starts[r]:ends[r]
            ok_r_idx <- rs[ok_v[rs]]
            if (length(ok_r_idx) >= 2L) {
              u_res_full[ok_r_idx] <- u_res_full[ok_r_idx] - mean(u_res_full[ok_r_idx])
            } else if (length(ok_r_idx) == 1L) {
              u_res_full[ok_r_idx] <- 0.0
            }
          }
          nrm_u <- sqrt(sum(u_res_full^2))
          if (is.finite(nrm_u) && nrm_u > 1e-12) {
            U_resid_mat[, vi] <- u_res_full / nrm_u
          }
        }
      } else {
        df_v <- max(1L, n_ok - max(1L, n_runs_active))
        tau2_raw_untrunc <- tau2_marginal_raw
        tau2_bio_unshrunk <- tau2_marginal
        s2_obs_scale <- s2_obs_marg
      }

      if (!is.null(u_marg_vec)) {
        nrm_m <- sqrt(sum(u_marg_vec^2))
        if (is.finite(nrm_m) && nrm_m > 1e-12) {
          U_resid_mat[, vi] <- u_marg_vec / nrm_m
        }
      }

      s2_obs_scale <- if (is.finite(s2_obs_scale) && s2_obs_scale > 1e-12) s2_obs_scale else max(s2_med, 1e-8)
      df_tau2_v <- if (identical(eff_tau2_mode, "marginal")) max(1L, df_marg) else df_v
      d_satt <- if (normF2_Omega_cent > 1e-24 && tr_Omega_cent > 1e-12) {
        (tr_Omega_cent * tr_Omega_cent) / normF2_Omega_cent
      } else {
        as.numeric(df_tau2_v)
      }
      u_df_v[vi] <- as.numeric(df_v)
      u_df_tau2[vi] <- as.numeric(df_tau2_v)
      u_d_eff_tau2[vi] <- max(2.0, min(as.numeric(df_tau2_v), d_satt))
      u_s2_obs_scale[vi] <- s2_obs_scale
      u_theta_raw[vi] <- tau2_raw_untrunc / s2_obs_scale
      u_tau2_unshrunk[vi] <- tau2_bio_unshrunk
    }

    # Empirical-Bayes shrinkage of untruncated biological variance ratio theta_raw = tau2_raw / s2_obs
    u_tau2_final <- u_tau2_unshrunk
    eb_fit_lag <- NULL
    n_usable_lag <- sum(u_usable)
    if (identical(eff_tau2_shrinkage, "graph_eb") && n_usable_lag >= 2L) {
      d_res_eb <- stats::median(u_df_tau2[u_usable], na.rm = TRUE)
      if (!is.finite(d_res_eb) || d_res_eb < 3.0) d_res_eb <- 3.0
      eb_fit_lag <- graph_eb_tau2_csr_cpp(
        row_ptr = epoch_graph$row_ptr,
        col_idx = epoch_graph$col_idx,
        dist = epoch_graph$dist,
        sigma_spatial = epoch_graph$sigma,
        theta_raw = u_theta_raw,
        d_eff = u_d_eff_tau2,
        task_gate = task_gate_vec,
        U_resid = U_resid_mat,
        d_resid = d_res_eb,
        group_ids = eb_group_ids,
        cross_parcel_weight = eb_cross_w
      )
      th_eb <- as.numeric(eb_fit_lag$theta_eb)
      ok_eb <- u_usable & is.finite(th_eb)
      u_tau2_final[ok_eb] <- pmax(0.0, th_eb[ok_eb] * u_s2_obs_scale[ok_eb])
    } else if (identical(eff_tau2_shrinkage, "eb") && n_usable_lag >= 2L) {
      d_res_eb <- stats::median(u_df_tau2[u_usable], na.rm = TRUE)
      if (!is.finite(d_res_eb) || d_res_eb < 3.0) d_res_eb <- 3.0
      eb_fit_lag <- eb_shrink_tau2_global_r(
        theta_raw = u_theta_raw,
        d_eff = u_d_eff_tau2,
        U_resid = U_resid_mat,
        d_resid = d_res_eb
      )
      th_eb <- as.numeric(eb_fit_lag$theta_eb)
      ok_eb <- u_usable & is.finite(th_eb)
      u_tau2_final[ok_eb] <- pmax(0.0, th_eb[ok_eb] * u_s2_obs_scale[ok_eb])
    }

    # Pass B: Precision weights, ZCA-whitened GLS, and raw sandwich scale factors
    u_beta_pw <- if (n_terms > 0L) matrix(NA_real_, n_terms, n_units) else NULL
    u_beta_gls <- if (n_terms > 0L) matrix(NA_real_, n_terms, n_units) else NULL
    u_v0_unw <- if (n_terms > 0L) matrix(NA_real_, n_terms, n_units) else NULL
    u_v0_pw <- if (n_terms > 0L) matrix(NA_real_, n_terms, n_units) else NULL
    u_v0_gls <- if (n_terms > 0L) matrix(NA_real_, n_terms, n_units) else NULL
    u_r_raw <- if (n_terms > 0L) matrix(NA_real_, n_terms, n_units) else NULL
    u_r_dem <- if (n_terms > 0L) matrix(NA_real_, n_terms, n_units) else NULL
    scale_raw_mat <- if (n_terms > 0L) matrix(NA_real_, n_units, 3L) else NULL
    u_d_eff_gls <- if (n_terms > 0L) rep(NA_real_, n_units) else NULL

    for (vi in which(u_usable)) {
      y_raw_v <- Y_raw_lag[, vi]
      y_dem_v <- Y_demix_lag[, vi]
      s_dem_v <- se_demixed_epochs[[li]][, vi]
      s2_dem_v <- s_dem_v * s_dem_v
      ok_s_dem <- is.finite(s2_dem_v) & (s2_dem_v > 0)
      ok_v <- if (n_terms > 0L) {
        is.finite(y_raw_v) & is.finite(y_dem_v) & x_row_finite
      } else {
        is.finite(y_raw_v) & is.finite(y_dem_v)
      }
      n_ok <- u_n_ok[vi]
      s2_med <- u_s2_med[vi]
      tau2_bio <- u_tau2_final[vi]
      tau2_bio_epochs[[li]][, vi] <- tau2_bio
      tau2_eff <- max(tau2_bio, 0.10 * s2_med)

      if (any(ok_s_dem)) {
        w_eb_all <- 1 / (tau2_eff + s2_dem_v[ok_s_dem])
        pw_cohort_epochs[[li]][ok_s_dem, vi] <- w_eb_all
        pw_epochs[[li]][ok_s_dem, vi] <- w_eb_all / mean(w_eb_all)
        rel_epochs[[li]][ok_s_dem, vi] <- tau2_eff / (tau2_eff + s2_dem_v[ok_s_dem])
      } else {
        pw_cohort_epochs[[li]][ok_v, vi] <- 1
        pw_epochs[[li]][ok_v, vi] <- 1
        rel_epochs[[li]][ok_v, vi] <- 1
      }
      w_v_all <- pw_epochs[[li]][, vi]
      w_v_all[!is.finite(w_v_all) | w_v_all <= 0] <- 1
      rel_v_all <- rel_epochs[[li]][, vi]
      rel_v_all[!is.finite(rel_v_all) | rel_v_all < 0] <- 1

      for (r in seq_len(n_runs)) {
        ck_r <- ckey_mat[r, vi]
        if (!nzchar(ck_r)) next
        rs_ok <- if (identical(ck_r, "__FB__")) {
          fb_blocks[[vi]][[r]]$rs_r
        } else {
          get(ck_r, envir = tpl_cache[[r]], inherits = FALSE)$rs_r
        }
        yr_raw <- demixed_epochs[[li]][rs_ok, vi]
        mu_r <- mean(yr_raw)
        eb_y_epochs[[li]][rs_ok, vi] <- mu_r + rel_v_all[rs_ok] * (yr_raw - mu_r)
      }

      if (n_terms == 0L) {
        if (isTRUE(whiten_panel) || !is.null(run_trial_operators)) {
          run_ops_ul <- stats::setNames(vector("list", n_runs), fit$run_ids)
          for (r in seq_len(n_runs)) {
            ck_r <- ckey_mat[r, vi]
            if (!nzchar(ck_r)) next
            if (identical(ck_r, "__FB__")) {
              fb_r <- fb_blocks[[vi]][[r]]
              rs_ok <- fb_r$rs_r
              ok_r_local <- fb_r$ok_r_local
              m_r <- fb_r$m_r
              Omega_r <- fb_r$Omega_r
              Omega_cent_r <- fb_r$Omega_cent_r
              M_post_r <- fb_r$M_post_r
              M2_r <- fb_r$M2_r
            } else {
              tpl <- get(ck_r, envir = tpl_cache[[r]], inherits = FALSE)
              s2_eff_rv <- sig_eff_mat[r, vi]^2
              rs_ok <- tpl$rs_r
              ok_r_local <- tpl$ok_r_local
              m_r <- tpl$m_r
              Omega_r <- s2_eff_rv * tpl$Omega_u
              Omega_cent_r <- s2_eff_rv * tpl$Omega_cent_u
              M_post_r <- tpl$M_post_r
              M2_r <- tpl$M2_r
            }
            yr_raw <- demixed_epochs[[li]][rs_ok, vi]
            yr_cent <- yr_raw - mean(yr_raw)
            Phi_r <- NULL
            if (isTRUE(whiten_panel)) {
              tau2_zca <- max(tau2_eff, 0.15 * s2_med)
              Sigma_tot_r <- Omega_cent_r + tau2_zca * M2_r
              Phi_r <- compute_zca_whitening_matrix(Sigma_tot_r, eig_floor = 0.08, offdiag_scale = 0.90)
              yr_zca <- as.numeric(Phi_r %*% yr_cent)
              cent_div_r <- if (m_r > 1L) (1.0 - 1.0 / m_r) else 1.0
              s2_zca_r <- if (m_r > 1L) {
                pmax(0, rowSums((Phi_r %*% Omega_cent_r) * Phi_r)) / cent_div_r
              } else {
                pmax(0, rowSums((Phi_r %*% Omega_r) * Phi_r))
              }
              se_zca_r <- sqrt(s2_zca_r)
              rel_zca_r <- tau2_eff / (tau2_eff + s2_zca_r)
              whitened_y_epochs[[li]][rs_ok, vi] <- mean(yr_raw) + yr_zca
              se_whitened_epochs[[li]][rs_ok, vi] <- se_zca_r
              eb_whitened_y_epochs[[li]][rs_ok, vi] <- mean(yr_raw) + rel_zca_r * yr_zca
            }
            if (!is.null(run_trial_operators)) {
              run_ops_ul[[r]] <- list(
                ok_r_local = ok_r_local,
                M_post = M_post_r,
                Phi_zca = Phi_r
              )
            }
          }
          if (!is.null(run_trial_operators)) {
            run_trial_operators[[vi]][[li]] <- run_ops_ul
          }
        }
        next
      }

      # Resolution-restored precision-weighted WLS and within-run ZCA-whitened GLS
      X_ok <- X_pool[ok_v, , drop = FALSE]
      X_tilde_ok <- matrix(0, n_ok, n_terms)
      X_w_ok <- matrix(0, n_ok, n_terms)
      y_w_ok <- numeric(n_ok)
      X_gls_ok <- matrix(0, n_ok, n_terms)
      y_gls_ok <- numeric(n_ok)
      w_ok <- numeric(n_ok)
      Sigma_tot_list <- vector("list", n_runs)
      run_ops_ul <- if (!is.null(run_trial_operators)) stats::setNames(vector("list", n_runs), fit$run_ids) else NULL
      pos_cursor <- 0L

      for (r in seq_len(n_runs)) {
        ck_r <- ckey_mat[r, vi]
        if (!nzchar(ck_r)) next
        if (identical(ck_r, "__FB__")) {
          fb_r <- fb_blocks[[vi]][[r]]
          rs_ok <- fb_r$rs_r
          ok_r_local <- fb_r$ok_r_local
          m_r <- fb_r$m_r
          Omega_r <- fb_r$Omega_r
          Omega_cent_r <- fb_r$Omega_cent_r
          M_post_r <- fb_r$M_post_r
          M2_r <- fb_r$M2_r
          Xr_tilde <- fb_r$Xr_tilde_sub
        } else {
          tpl <- get(ck_r, envir = tpl_cache[[r]], inherits = FALSE)
          s2_eff_rv <- sig_eff_mat[r, vi]^2
          rs_ok <- tpl$rs_r
          ok_r_local <- tpl$ok_r_local
          m_r <- tpl$m_r
          Omega_r <- s2_eff_rv * tpl$Omega_u
          Omega_cent_r <- s2_eff_rv * tpl$Omega_cent_u
          M_post_r <- tpl$M_post_r
          M2_r <- tpl$M2_r
          Xr_tilde <- tpl$Xr_tilde_sub
        }
        idx_span <- (pos_cursor + 1L):(pos_cursor + m_r)
        pos_cursor <- pos_cursor + m_r

        wr <- w_v_all[rs_ok]
        w_ok[idx_span] <- wr
        X_tilde_ok[idx_span, ] <- Xr_tilde

        yr_raw <- demixed_epochs[[li]][rs_ok, vi]
        yr_cent <- yr_raw - mean(yr_raw)
        Xr_raw <- X_raw_pool[rs_ok, , drop = FALSE]

        Q_wr <- if (m_r > 1L) {
          diag(m_r) - outer(rep(1.0, m_r), wr / sum(wr))
        } else {
          matrix(0.0, 1L, 1L)
        }
        y_w_ok[idx_span] <- as.numeric(Q_wr %*% yr_cent)
        X_w_ok[idx_span, ] <- Q_wr %*% Xr_tilde

        tau2_zca <- max(tau2_eff, 0.15 * s2_med)
        Sigma_tot_r <- Omega_cent_r + tau2_zca * M2_r
        Sigma_w_r <- Q_wr %*% Sigma_tot_r %*% t(Q_wr)

        Phi_r <- compute_zca_whitening_matrix(Sigma_tot_r, eig_floor = 0.08, offdiag_scale = 0.90)
        yr_zca <- as.numeric(Phi_r %*% yr_cent)
        Xr_zca <- Phi_r %*% Xr_tilde
        Sigma_zca_tot_r <- Phi_r %*% Sigma_tot_r %*% t(Phi_r)
        Sigma_gls_r <- Q_wr %*% Sigma_zca_tot_r %*% t(Q_wr)

        Sigma_tot_list[[r]] <- list(
          idx = idx_span,
          Sigma_unw = Sigma_tot_r,
          Sigma_pw = Sigma_w_r,
          Sigma_gls = Sigma_gls_r,
          wr = wr
        )
        if (!is.null(run_ops_ul)) {
          run_ops_ul[[r]] <- list(
            ok_r_local = ok_r_local,
            M_post = M_post_r,
            Phi_zca = if (isTRUE(whiten_panel)) Phi_r else NULL
          )
        }

        if (isTRUE(whiten_panel)) {
          cent_div_r <- if (m_r > 1L) (1.0 - 1.0 / m_r) else 1.0
          s2_zca_r <- if (m_r > 1L) {
            pmax(0, rowSums((Phi_r %*% Omega_cent_r) * Phi_r)) / cent_div_r
          } else {
            pmax(0, rowSums((Phi_r %*% Omega_r) * Phi_r))
          }
          se_zca_r <- sqrt(s2_zca_r)
          rel_zca_r <- tau2_eff / (tau2_eff + s2_zca_r)
          whitened_y_epochs[[li]][rs_ok, vi] <- mean(yr_raw) + yr_zca
          se_whitened_epochs[[li]][rs_ok, vi] <- se_zca_r
          eb_whitened_y_epochs[[li]][rs_ok, vi] <- mean(yr_raw) + rel_zca_r * yr_zca
          for (ti in seq_len(n_terms)) {
            whitened_x_epochs[[li]][[ti]][rs_ok, vi] <- mean(Xr_raw[, ti]) + Xr_zca[, ti]
          }
        }
        y_gls_ok[idx_span] <- as.numeric(Q_wr %*% yr_zca)
        X_gls_ok[idx_span, ] <- Q_wr %*% Xr_zca
      }
      if (!is.null(run_trial_operators)) {
        run_trial_operators[[vi]][[li]] <- run_ops_ul
      }

      xtx_inv_tilde <- u_xtx_inv_tilde[[vi]]
      meat_unw <- matrix(0, n_terms, n_terms)
      tr_Sigma_unw <- 0.0
      for (r in seq_len(n_runs)) {
        st_r <- Sigma_tot_list[[r]]
        if (is.null(st_r)) next
        idx_r <- st_r$idx
        X_u_r <- X_tilde_ok[idx_r, , drop = FALSE]
        Sigma_u_r <- st_r$Sigma_unw
        meat_u_r <- crossprod(X_u_r, Sigma_u_r %*% X_u_r)
        meat_unw <- meat_unw + meat_u_r
        tr_Sigma_unw <- tr_Sigma_unw + sum(diag(Sigma_u_r)) - sum(diag(xtx_inv_tilde %*% meat_u_r))
      }
      scale_unw_raw <- u_ss_resid_unw[vi] / max(tr_Sigma_unw, 1e-8)
      u_v0_unw[, vi] <- pmax(diag(xtx_inv_tilde %*% meat_unw %*% xtx_inv_tilde), 0.0)

      w_ok <- w_ok / mean(w_ok)
      WX_w <- X_w_ok * w_ok
      xwx <- crossprod(X_w_ok, WX_w) + 1e-8 * diag(n_terms)
      xwx_inv <- solve(xwx)
      beta_pw_v <- as.numeric(xwx_inv %*% crossprod(WX_w, y_w_ok))
      resid_pw <- y_w_ok - as.numeric(X_w_ok %*% beta_pw_v)

      meat_pw <- matrix(0, n_terms, n_terms)
      tr_W_Sigma <- 0.0
      for (r in seq_len(n_runs)) {
        st_r <- Sigma_tot_list[[r]]
        if (is.null(st_r)) next
        idx_r <- st_r$idx
        WX_r <- WX_w[idx_r, , drop = FALSE]
        Sigma_r <- st_r$Sigma_pw
        wr <- w_ok[idx_r]
        meat_r <- crossprod(WX_r, Sigma_r %*% WX_r)
        meat_pw <- meat_pw + meat_r
        tr_W_Sigma <- tr_W_Sigma + sum(wr * diag(Sigma_r)) - sum(diag(xwx_inv %*% meat_r))
      }
      scale_pw_raw <- sum(w_ok * (resid_pw^2)) / max(tr_W_Sigma, 1e-8)
      u_v0_pw[, vi] <- pmax(diag(xwx_inv %*% meat_pw %*% xwx_inv), 0.0)

      WX_gls <- X_gls_ok * w_ok
      xwx_gls <- crossprod(X_gls_ok, WX_gls) + 1e-8 * diag(n_terms)
      xwx_gls_inv <- solve(xwx_gls)
      beta_gls_v <- as.numeric(xwx_gls_inv %*% crossprod(WX_gls, y_gls_ok))
      resid_gls <- y_gls_ok - as.numeric(X_gls_ok %*% beta_gls_v)
      meat_gls <- matrix(0, n_terms, n_terms)
      tr_W_Sigma_gls <- 0.0
      normF2_W_Sigma_gls <- 0.0
      for (r in seq_len(n_runs)) {
        st_r <- Sigma_tot_list[[r]]
        if (is.null(st_r)) next
        idx_r <- st_r$idx
        WX_g_r <- WX_gls[idx_r, , drop = FALSE]
        Sigma_g_r <- st_r$Sigma_gls
        wr_g <- w_ok[idx_r]
        meat_g_r <- crossprod(WX_g_r, Sigma_g_r %*% WX_g_r)
        meat_gls <- meat_gls + meat_g_r
        tr_W_Sigma_gls <- tr_W_Sigma_gls + sum(wr_g * diag(Sigma_g_r)) - sum(diag(xwx_gls_inv %*% meat_g_r))
        WS_r <- sweep(Sigma_g_r, 1L, wr_g, "*")
        normF2_W_Sigma_gls <- normF2_W_Sigma_gls + sum(WS_r * t(WS_r))
      }
      scale_gls_raw <- sum(w_ok * (resid_gls^2)) / max(tr_W_Sigma_gls, 1e-8)
      u_v0_gls[, vi] <- pmax(diag(xwx_gls_inv %*% meat_gls %*% xwx_gls_inv), 0.0)

      u_beta_pw[, vi] <- beta_pw_v
      u_beta_gls[, vi] <- beta_gls_v
      scale_raw_mat[vi, 1L] <- scale_unw_raw
      scale_raw_mat[vi, 2L] <- scale_pw_raw
      scale_raw_mat[vi, 3L] <- scale_gls_raw
      d_gls_satt <- if (normF2_W_Sigma_gls > 1e-24 && tr_W_Sigma_gls > 1e-12) {
        (tr_W_Sigma_gls * tr_W_Sigma_gls) / normF2_W_Sigma_gls
      } else {
        u_df_v[vi]
      }
      u_d_eff_gls[vi] <- max(2.0, min(u_df_v[vi], d_gls_satt))

      y_raw_ok <- y_raw_v[ok_v]
      y_dem_ok <- y_dem_v[ok_v]
      sd_y_raw <- stats::sd(y_raw_ok)
      sd_y_dem <- stats::sd(y_dem_ok)
      for (ti in seq_len(n_terms)) {
        xk <- X_ok[, ti]
        xk_tilde <- X_tilde_ok[, ti]
        u_r_raw[ti, vi] <- if (stats::sd(xk) > 0 && sd_y_raw > 0) stats::cor(xk, y_raw_ok) else NA_real_
        u_r_dem[ti, vi] <- if (stats::sd(xk_tilde) > 0 && sd_y_dem > 0) stats::cor(xk_tilde, y_dem_ok) else NA_real_
      }
    }

    # Empirical-Bayes moderation of sandwich residual scale factors [scale_unw, scale_pw, scale_gls]
    scale_mod_mat <- if (n_terms > 0L) matrix(pmax(scale_raw_mat, 0.25), n_units, 3L) else NULL
    u_df_out <- u_df_v
    if (n_terms > 0L && !is.null(eb_fit_lag) && n_usable_lag >= 2L) {
      if (identical(eff_tau2_shrinkage, "graph_eb")) {
        mod_res <- moderate_scale_csr_cpp(
          alpha_row_ptr = eb_fit_lag$alpha_row_ptr,
          alpha_col_idx = eb_fit_lag$alpha_col_idx,
          alpha_w = eb_fit_lag$alpha_w,
          kappa_minus = eb_fit_lag$kappa_minus,
          c_self_prior = eb_fit_lag$c_self_prior,
          scale_raw = scale_raw_mat,
          d_eff = u_d_eff_gls,
          df_nom = u_df_v
        )
        scale_mod_mat <- mod_res$scale_mod
        u_df_out[u_usable] <- as.numeric(mod_res$df_mod)[u_usable]
      } else if (identical(eff_tau2_shrinkage, "eb")) {
        mod_res <- eb_moderate_scale_global_r(
          scale_raw = scale_raw_mat,
          d_eff = u_d_eff_gls,
          df_nom = u_df_v,
          eb_fit = eb_fit_lag
        )
        scale_mod_mat <- mod_res$scale_mod
        u_df_out[u_usable] <- as.numeric(mod_res$df_mod)[u_usable]
      }
    }

    # Populate pre-allocated coefficient arrays in exact unit/term order
    for (vi in seq_len(n_units)) {
      if (!u_usable[vi]) {
        n_ok <- u_n_ok[vi]
        if (n_terms > 0L) {
          for (ti in seq_len(n_terms)) {
            row_k <- row_k + 1L
            cf_voxel[row_k] <- unit_ids[vi]
            cf_lag[row_k] <- lags[li]
            cf_term[row_k] <- term_names[ti]
            cf_df[row_k] <- n_ok - n_terms
            cf_gain[row_k] <- gain_by_lag_voxel[li, vi]
            cf_mean_se[row_k] <- u_mean_se[vi]
            cf_mean_ac1[row_k] <- u_mean_ac1[vi]
          }
        } else {
          row_k <- row_k + 1L
          cf_voxel[row_k] <- unit_ids[vi]
          cf_lag[row_k] <- lags[li]
          cf_df[row_k] <- max(0L, n_ok - n_runs)
          cf_gain[row_k] <- gain_by_lag_voxel[li, vi]
          cf_mean_se[row_k] <- u_mean_se[vi]
          cf_mean_ac1[row_k] <- u_mean_ac1[vi]
        }
        next
      }

      if (n_terms == 0L) {
        row_k <- row_k + 1L
        cf_voxel[row_k] <- unit_ids[vi]
        cf_lag[row_k] <- lags[li]
        cf_df[row_k] <- u_df_v[vi]
        cf_tau2[row_k] <- u_tau2_final[vi]
        cf_gain[row_k] <- gain_by_lag_voxel[li, vi]
        cf_mean_se[row_k] <- u_mean_se[vi]
        cf_mean_ac1[row_k] <- u_mean_ac1[vi]
        next
      }

      beta_raw_v <- u_beta_raw[, vi]
      beta_unw_v <- u_beta_unw[, vi]
      beta_pw_v <- u_beta_pw[, vi]
      beta_gls_v <- u_beta_gls[, vi]
      se_unw_v <- sqrt(pmax(scale_mod_mat[vi, 1L] * u_v0_unw[, vi], .Machine$double.eps))
      se_pw_v <- sqrt(pmax(scale_mod_mat[vi, 2L] * u_v0_pw[, vi], .Machine$double.eps))
      se_gls_v <- sqrt(pmax(scale_mod_mat[vi, 3L] * u_v0_gls[, vi], .Machine$double.eps))
      t_unw_vec <- beta_unw_v / se_unw_v
      t_pw_vec <- beta_pw_v / se_pw_v
      t_gls_vec <- beta_gls_v / se_gls_v

      if (identical(weighting, "precision")) {
        beta_corr_v <- beta_pw_v
        se_corr_v <- se_pw_v
        t_corr_vec <- t_pw_vec
      } else if (identical(weighting, "gls")) {
        beta_corr_v <- beta_gls_v
        se_corr_v <- se_gls_v
        t_corr_vec <- t_gls_vec
      } else {
        beta_corr_v <- beta_unw_v
        se_corr_v <- se_unw_v
        t_corr_vec <- t_unw_vec
      }
      df_out_v <- u_df_out[vi]
      p_corr_vec <- 2 * stats::pt(abs(t_corr_vec), df = df_out_v, lower.tail = FALSE)

      for (ti in seq_len(n_terms)) {
        row_k <- row_k + 1L
        cf_voxel[row_k] <- unit_ids[vi]
        cf_lag[row_k] <- lags[li]
        cf_term[row_k] <- term_names[ti]
        cf_beta_raw[row_k] <- beta_raw_v[ti]
        cf_beta_corr[row_k] <- beta_corr_v[ti]
        cf_beta_unw[row_k] <- beta_unw_v[ti]
        cf_beta_pw[row_k] <- beta_pw_v[ti]
        cf_beta_gls[row_k] <- beta_gls_v[ti]
        cf_se_corr[row_k] <- se_corr_v[ti]
        cf_se_unw[row_k] <- se_unw_v[ti]
        cf_se_pw[row_k] <- se_pw_v[ti]
        cf_se_gls[row_k] <- se_gls_v[ti]
        cf_r_raw[row_k] <- u_r_raw[ti, vi]
        cf_r_dem[row_k] <- u_r_dem[ti, vi]
        cf_t_stat[row_k] <- t_corr_vec[ti]
        cf_p_val[row_k] <- p_corr_vec[ti]
        cf_t_unw[row_k] <- t_unw_vec[ti]
        cf_t_pw[row_k] <- t_pw_vec[ti]
        cf_t_gls[row_k] <- t_gls_vec[ti]
        cf_df[row_k] <- df_out_v
        cf_tau2[row_k] <- u_tau2_final[vi]
        cf_gain[row_k] <- gain_by_lag_voxel[li, vi]
        cf_mean_se[row_k] <- u_mean_se[vi]
        cf_mean_ac1[row_k] <- u_mean_ac1[vi]
      }
    }
  }

  keep_cf <- seq_len(row_k)
  coef_df <- data.frame(
    voxel = cf_voxel[keep_cf],
    lag = cf_lag[keep_cf],
    term = cf_term[keep_cf],
    beta_raw = cf_beta_raw[keep_cf],
    beta_corrected = cf_beta_corr[keep_cf],
    beta_unweighted = cf_beta_unw[keep_cf],
    beta_weighted = cf_beta_pw[keep_cf],
    beta_gls = cf_beta_gls[keep_cf],
    se_corrected = cf_se_corr[keep_cf],
    se_unweighted = cf_se_unw[keep_cf],
    se_weighted = cf_se_pw[keep_cf],
    se_gls = cf_se_gls[keep_cf],
    r_raw = cf_r_raw[keep_cf],
    r_demixed = cf_r_dem[keep_cf],
    t_stat = cf_t_stat[keep_cf],
    p_value = cf_p_val[keep_cf],
    t_unweighted = cf_t_unw[keep_cf],
    t_weighted = cf_t_pw[keep_cf],
    t_gls = cf_t_gls[keep_cf],
    df = cf_df[keep_cf],
    tau2_bio = cf_tau2[keep_cf],
    mean_gain = cf_gain[keep_cf],
    mean_se_trial = cf_mean_se[keep_cf],
    mean_error_ac1 = cf_mean_ac1[keep_cf],
    stringsAsFactors = FALSE
  )
  if (isTRUE(aggregate_by_parcel)) {
    coef_df$parcel <- coef_df$voxel
  }

  lag_names <- as.character(lags)
  lag_error_cov <- stats::setNames(vector("list", n_units), unit_ids)
  lag_error_cor <- stats::setNames(vector("list", n_units), unit_ids)
  lag_resolution <- stats::setNames(vector("list", n_units), unit_ids)
  for (vi in seq_len(n_units)) {
    if (n_lags >= 2L && lag_cov_trials[vi] > 0L) {
      V_u <- lag_cov_accum[[vi]] / lag_cov_trials[vi]
      V_u <- 0.5 * (V_u + t(V_u))
      d_u <- sqrt(pmax(diag(V_u), .Machine$double.eps))
      R_u <- V_u / outer(d_u, d_u)
      R_u[!is.finite(R_u)] <- 0
      diag(R_u) <- 1
    } else {
      V_u <- diag(1, n_lags)
      R_u <- diag(1, n_lags)
    }
    if (n_lags >= 2L && lag_res_trials[vi] > 0L) {
      C_u <- lag_res_accum[[vi]] / lag_res_trials[vi]
      C_u <- 0.5 * (C_u + t(C_u))
    } else {
      C_u <- diag(1, n_lags)
    }
    dimnames(V_u) <- dimnames(R_u) <- dimnames(C_u) <- list(lag_names, lag_names)
    lag_error_cov[[vi]] <- V_u
    lag_error_cor[[vi]] <- R_u
    lag_resolution[[vi]] <- C_u
  }

  if (nrow(coef_df) > 0L) {
    coef_df$beta_multilag <- coef_df$beta_corrected
    coef_df$se_multilag <- coef_df$se_corrected
    if (n_lags >= 3L && n_terms > 0L) {
      for (vi in seq_len(n_units)) {
        u_id <- unit_ids[vi]
        C_u <- lag_resolution[[vi]]
        R_u <- lag_error_cor[[vi]]
        for (ti in seq_len(n_terms)) {
          idx_vt <- which(coef_df$voxel == u_id & coef_df$term == term_names[ti])
          if (length(idx_vt) == n_lags) {
            ord_l <- match(lags, coef_df$lag[idx_vt])
            idx_vt <- idx_vt[ord_l]
            b_v <- coef_df$beta_corrected[idx_vt]
            se_v <- coef_df$se_corrected[idx_vt]
            if (all(is.finite(b_v)) && all(is.finite(se_v) & se_v > 0)) {
              Sigma_b <- outer(se_v, se_v) * R_u
              sm_res <- medusa_smooth_deblur_trajectory(
                beta_vec = b_v,
                Sigma_beta = Sigma_b,
                lags = lags,
                C_sym = C_u,
                tr = fit$settings$tr,
                demix_within_trial = isTRUE(demix_within_trial) || !isTRUE(demix_target_trials),
                within_trial_ridge = within_trial_ridge,
                trajectory_mode = "auto"
              )
              coef_df$beta_multilag[idx_vt] <- sm_res$beta
              coef_df$se_multilag[idx_vt] <- sm_res$se
            }
          }
        }
      }
    }
    attr(coef_df, "lag_error_cov") <- lag_error_cov
    attr(coef_df, "lag_error_cor") <- lag_error_cor
    attr(coef_df, "lag_resolution") <- lag_resolution
  }

  fmean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
  summary_rows <- list()
  for (li in seq_len(n_lags)) {
    if (n_terms > 0L) {
      for (ti in seq_len(n_terms)) {
        sub_df <- coef_df[coef_df$lag == lags[li] & coef_df$term == term_names[ti], , drop = FALSE]
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          lag = lags[li], term = term_names[ti],
          mean_beta_raw = fmean(sub_df$beta_raw),
          mean_beta_corrected = fmean(sub_df$beta_corrected),
          mean_beta_multilag = fmean(sub_df$beta_multilag),
          mean_beta_unweighted = fmean(sub_df$beta_unweighted),
          mean_beta_weighted = fmean(sub_df$beta_weighted),
          mean_beta_gls = fmean(sub_df$beta_gls),
          mean_r_raw = fmean(sub_df$r_raw),
          mean_r_demixed = fmean(sub_df$r_demixed),
          mean_tau2_bio = fmean(sub_df$tau2_bio),
          mean_gain = fmean(sub_df$mean_gain),
          mean_se_trial = fmean(sub_df$mean_se_trial),
          mean_error_ac1 = fmean(sub_df$mean_error_ac1),
          n_voxels = sum(is.finite(sub_df$beta_corrected)),
          stringsAsFactors = FALSE
        )
      }
    } else {
      sub_df <- coef_df[coef_df$lag == lags[li], , drop = FALSE]
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        lag = lags[li], term = NA_character_,
        mean_beta_raw = NA_real_,
        mean_beta_corrected = NA_real_,
        mean_beta_multilag = NA_real_,
        mean_beta_unweighted = NA_real_,
        mean_beta_weighted = NA_real_,
        mean_beta_gls = NA_real_,
        mean_r_raw = NA_real_,
        mean_r_demixed = NA_real_,
        mean_tau2_bio = fmean(sub_df$tau2_bio),
        mean_gain = fmean(sub_df$mean_gain),
        mean_se_trial = fmean(sub_df$mean_se_trial),
        mean_error_ac1 = fmean(sub_df$mean_error_ac1),
        n_voxels = sum(is.finite(sub_df$tau2_bio)),
        stringsAsFactors = FALSE
      )
    }
  }

  # v3 Subspace-Debiased Adaptive-POET + Multi-Lag Toeplitz & General-W_s Covariance-Regression
  # Shared-Error Engine (propagated after univariate Pass A / Graph-EB so marginal tau2_bio and
  # precision_weight invariance are preserved identically while Step-6 connectivity & panel SEs
  # gain exact cross-parcel, dynamic burst, multi-lag, and cross-trial error moments).
  shared_error_model <- NULL
  sem_fit <- NULL
  if (identical(shared_error_method, "poet_covreg") && n_units >= 1L) {
    sem_fit <- medusa_fit_shared_error_covreg(
      deconvolution = deconvolution,
      events = events,
      target_idx_list = target_idx_list,
      dur_target_list = if (!is.null(trial_durations)) dur_target_list else NULL,
      cov_list = cov_list,
      lags = lags,
      readout_width = readout_width,
      restore_scale = restore_scale,
      source_model = source_model,
      aggregate_by_parcel = aggregate_by_parcel,
      unit_ids = unit_ids,
      parcel_vox_list = parcel_vox_list,
      parcel_weights_list = parcel_weights_list,
      run_stratum_ops = run_stratum_ops,
      variance_confounds = variance_confounds,
      bold = bold,
      return_panel = return_panel,
      return_epochs = return_epochs
    )
    unit_error_cor <- sem_fit$unit_error_cor
    shared_loading_vec <- sem_fit$shared_loading_vec
    shared_error_model <- sem_fit$shared_error_model
    if (!is.null(sem_fit$run_panel_blocks)) {
      for (blk in sem_fit$run_panel_blocks) {
        rs_b <- starts[blk$orig_r]:ends[blk$orig_r]
        for (li in seq_len(n_lags)) {
          se_demixed_epochs[[li]][rs_b, ] <- blk$se[, , li]
        }
      }
    }
  }

  panel_dt <- NULL
  if (isTRUE(return_panel)) {
    run_vec <- rep(fit$run_ids, times = target_counts)
    trial_vec <- unlist(lapply(target_counts, seq_len), use.names = FALSE)
    evt_idx_vec <- unlist(target_idx_list, use.names = FALSE)
    onset_vec <- unlist(lapply(seq_len(n_runs), function(r) events[[r]]$onset[target_idx_list[[r]]]), use.names = FALSE)
    has_id_col <- all(vapply(events, function(ev) "id" %in% names(ev), logical(1L)))
    id_vec <- if (has_id_col) {
      unlist(lapply(seq_len(n_runs), function(r) as.character(events[[r]]$id[target_idx_list[[r]]])), use.names = FALSE)
    } else NULL
    has_ses_col <- all(vapply(events, function(ev) "session" %in% names(ev), logical(1L)))
    ses_vec <- if (has_ses_col) {
      unlist(lapply(seq_len(n_runs), function(r) as.character(events[[r]]$session[target_idx_list[[r]]])), use.names = FALSE)
    } else NULL
    cov_full <- if (n_terms > 0L) as.data.frame(do.call(rbind, cov_list)) else NULL
    shared_load_rep <- rep(as.numeric(shared_loading_vec), each = total_trials)
    has_v3_blocks <- !is.null(sem_fit) && !is.null(sem_fit$run_panel_blocks) && length(sem_fit$run_panel_blocks) >= 1L
    n_rows_per_lag <- total_trials * n_units

    panel_list <- vector("list", n_lags)
    for (li in seq_len(n_lags)) {
      base_meta <- data.table::data.table(
        run = rep(run_vec, times = n_units),
        trial = rep(trial_vec, times = n_units),
        event_index = rep(evt_idx_vec, times = n_units),
        onset = rep(onset_vec, times = n_units),
        evt_time = lags[li],
        voxel = rep(unit_ids, each = total_trials),
        decon_raw = as.numeric(raw_epochs[[li]]),
        decon_demixed = as.numeric(demixed_epochs[[li]]),
        decon_eb = as.numeric(eb_y_epochs[[li]]),
        se_raw = as.numeric(se_raw_epochs[[li]]),
        se_demixed = as.numeric(se_demixed_epochs[[li]]),
        reliability_trial = as.numeric(rel_epochs[[li]]),
        precision_weight = as.numeric(pw_epochs[[li]]),
        precision_weight_cohort = as.numeric(pw_cohort_epochs[[li]]),
        precision_weight_op = as.numeric(pw_op_epochs[[li]]),
        precision_weight_raw = as.numeric(pw_raw_epochs[[li]]),
        error_ac1 = as.numeric(error_ac1_epochs[[li]]),
        tau2_bio = as.numeric(tau2_bio_epochs[[li]]),
        shared_noise_loading = shared_load_rep
      )
      if (has_v3_blocks) {
        d_svec <- dim(sem_fit$run_panel_blocks[[1L]]$svec)[1L]
        d_tcor <- dim(sem_fit$run_panel_blocks[[1L]]$tcor)[1L]
        svec_arr_li <- array(0.0, dim = c(d_svec, total_trials, n_units))
        tcor_arr_li <- array(0.0, dim = c(d_tcor, total_trials, n_units))
        lvec_arr_li <- if (n_lags >= 2L && !is.null(sem_fit$run_panel_blocks[[1L]]$lvec)) {
          array(0.0, dim = c(dim(sem_fit$run_panel_blocks[[1L]]$lvec)[1L], total_trials, n_units))
        } else NULL
        for (blk in sem_fit$run_panel_blocks) {
          rs_b <- starts[blk$orig_r]:ends[blk$orig_r]
          svec_arr_li[, rs_b, ] <- blk$svec[, , , li]
          tcor_arr_li[, rs_b, ] <- blk$tcor[, , , li]
          if (!is.null(lvec_arr_li) && !is.null(blk$lvec)) {
            lvec_arr_li[, rs_b, ] <- blk$lvec[, , , li]
          }
        }
        dim(svec_arr_li) <- c(d_svec, n_rows_per_lag)
        dim(tcor_arr_li) <- c(d_tcor, n_rows_per_lag)
        base_meta[["shared_noise_vec"]] <- asplit(svec_arr_li, 2L)
        base_meta[["trial_cor_vec"]] <- asplit(tcor_arr_li, 2L)
        if (!is.null(lvec_arr_li)) {
          dim(lvec_arr_li) <- c(dim(lvec_arr_li)[1L], n_rows_per_lag)
          base_meta[["lag_noise_vec"]] <- asplit(lvec_arr_li, 2L)
        }
      }
      if (!is.null(id_vec) && (is.null(cov_full) || !("id" %in% colnames(cov_full)))) {
        base_meta[["id"]] <- rep(id_vec, times = n_units)
      }
      if (!is.null(ses_vec) && (is.null(cov_full) || !("session" %in% colnames(cov_full)))) {
        base_meta[["session"]] <- rep(ses_vec, times = n_units)
      }
      if (isTRUE(aggregate_by_parcel)) {
        base_meta[["parcel"]] <- base_meta[["voxel"]]
      }
      if (!is.null(cov_full)) {
        for (cn in colnames(cov_full)) {
          base_meta[[cn]] <- rep(cov_full[[cn]], times = n_units)
        }
      }
      if (isTRUE(whiten_panel)) {
        base_meta[["decon_demixed_whitened"]] <- as.numeric(whitened_y_epochs[[li]])
        base_meta[["se_whitened"]] <- as.numeric(se_whitened_epochs[[li]])
        base_meta[["decon_eb_whitened"]] <- as.numeric(eb_whitened_y_epochs[[li]])
        if (n_terms > 0L) {
          for (ti in seq_len(n_terms)) {
            cn_w <- paste0(term_names[ti], "_whitened")
            base_meta[[cn_w]] <- as.numeric(whitened_x_epochs[[li]][[ti]])
          }
        }
      }
      panel_list[[li]] <- base_meta
    }
    panel_dt <- data.table::rbindlist(panel_list)
    attr(panel_dt, "lag_error_cov") <- lag_error_cov
    attr(panel_dt, "lag_error_cor") <- lag_error_cor
    attr(panel_dt, "lag_resolution") <- lag_resolution
    attr(panel_dt, "unit_error_cor") <- unit_error_cor
    attr(panel_dt, "shared_error_model") <- shared_error_model
    attr(panel_dt, "run_trial_operators") <- run_trial_operators
    attr(panel_dt, "tr") <- fit$settings$tr
    attr(panel_dt, "raw_gain_by_lag_unit") <- gain_by_lag_voxel
    attr(panel_dt, "demix_target_trials") <- isTRUE(demix_target_trials)
  }

  list(
    coefficients = coef_df,
    summary = do.call(rbind, summary_rows),
    epochs = if (return_epochs) {
      list(
        raw = raw_epochs,
        demixed = demixed_epochs,
        eb = eb_y_epochs,
        se_raw = se_raw_epochs,
        se_demixed = se_demixed_epochs,
        reliability_trial = rel_epochs,
        precision_weight = pw_epochs,
        precision_weight_cohort = pw_cohort_epochs,
        precision_weight_op = pw_op_epochs,
        error_ac1 = error_ac1_epochs,
        tau2_bio = tau2_bio_epochs,
        whitened = whitened_y_epochs,
        se_whitened = se_whitened_epochs,
        eb_whitened = eb_whitened_y_epochs
      )
    } else NULL,
    panel = panel_dt,
    parcel_weights = parcel_weights_list,
    lag_error_cov = lag_error_cov,
    lag_error_cor = lag_error_cor,
    lag_resolution = lag_resolution,
    unit_error_cor = unit_error_cor,
    shared_error_model = shared_error_model,
    run_trial_operators = run_trial_operators,
    settings = list(
      target_condition = target_condition,
      lags = lags,
      readout_kernel = readout_kernel,
      readout_width = readout_width,
      demix_other_events = demix_other_events,
      demix_target_trials = demix_target_trials,
      demix_within_trial = demix_within_trial,
      demix_ridge = demix_ridge,
      within_trial_ridge = within_trial_ridge,
      gain_floor = gain_floor,
      restore_scale = restore_scale,
      source_model = source_model,
      duration_scale = duration_scale,
      weighting = weighting,
      tau2_mode = eff_tau2_mode,
      tau2_shrinkage = tau2_shrinkage,
      tau2_shrinkage_used = eff_tau2_shrinkage,
      aggregate_by_parcel = aggregate_by_parcel,
      whiten_panel = whiten_panel,
      shared_error_method = shared_error_method
    )
  )
}

# Smallest finite real generalized eigenvalue solving det(M - lambda * Psi) = 0
# over the error-prone subspace (projecting out error-free zero columns of Psi via Schur complement).
medusa_smallest_det_root <- function(M, Psi) {
  if (requireNamespace("geigen", quietly = TRUE)) {
    roots <- tryCatch(
      geigen::geigen(M, Psi, symmetric = FALSE, only.values = TRUE)$values,
      error = function(e) NULL
    )
    if (!is.null(roots)) {
      if (is.complex(roots)) {
        sc <- suppressWarnings(max(abs(roots), na.rm = TRUE))
        if (!is.finite(sc)) sc <- 1
        roots <- ifelse(abs(Im(roots)) <= sqrt(.Machine$double.eps) * max(1, sc), Re(roots), NA_real_)
      }
      roots <- roots[is.finite(roots)]
      if (length(roots)) return(min(roots))
    }
  }
  d_psi <- diag(Psi)
  E_idx <- which(is.finite(d_psi) & d_psi > 1e-14)
  if (!length(E_idx)) return(NA_real_)
  F_idx <- setdiff(seq_len(nrow(M)), E_idx)
  M_EF <- if (length(F_idx)) {
    M_FF <- M[F_idx, F_idx, drop = FALSE]
    M[E_idx, E_idx, drop = FALSE] -
      M[E_idx, F_idx, drop = FALSE] %*% solve(M_FF + 1e-12 * diag(length(F_idx)), M[F_idx, E_idx, drop = FALSE])
  } else {
    M[E_idx, E_idx, drop = FALSE]
  }
  Psi_EE <- 0.5 * (Psi[E_idx, E_idx, drop = FALSE] + t(Psi[E_idx, E_idx, drop = FALSE]))
  eig_P <- eigen(Psi_EE, symmetric = TRUE)
  pos <- eig_P$values > 1e-12 * max(1, max(eig_P$values))
  if (!any(pos)) return(NA_real_)
  U_pos <- eig_P$vectors[, pos, drop = FALSE]
  d_inv_half <- 1 / sqrt(eig_P$values[pos])
  H <- sweep(sweep(crossprod(U_pos, M_EF %*% U_pos), 1L, d_inv_half, "*"), 2L, d_inv_half, "*")
  H <- 0.5 * (H + t(H))
  vals <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  vals <- vals[is.finite(vals)]
  if (length(vals)) min(vals) else NA_real_
}

medusa_mmom_correction <- function(lambda, nobs, alpha) {
  if (is.finite(lambda) && lambda <= 1 + 1 / nobs) {
    lambda - (1 + alpha) / nobs
  } else {
    1 - alpha / nobs
  }
}

# Built-in fallback implementation of Vig & Hallquist (2026) M-MOM / M-MOMS
# matching repair::fit_mmom() when the external `repair` package is not installed.
medusa_fit_mmom_fallback <- function(
  data_list,
  preliminary_moment_control = list(correction = "step_up", alpha = NULL),
  final_moment_control = list(correction = "step_up", alpha = NULL)
) {
  y <- as.numeric(data_list$outcome)
  X_pred <- as.matrix(data_list$predictors)
  storage.mode(X_pred) <- "double"
  nobs <- length(y)
  q <- ncol(X_pred)
  p <- q + 1L
  if (is.null(colnames(X_pred))) colnames(X_pred) <- paste0("x", seq_len(q))
  design <- cbind(`(Intercept)` = 1.0, X_pred)
  pred_idx <- 2L:p
  err_x <- data_list$err_x
  err_y <- if (is.null(data_list$err_y)) rep(0.0, nobs) else as.numeric(data_list$err_y)
  err_xy <- if (is.null(data_list$err_xy)) matrix(0.0, nobs, q) else as.matrix(data_list$err_xy)

  alpha_min <- q + 2.0
  norm_ctrl <- function(ctrl) {
    corr <- if (is.null(ctrl$correction)) "minimum" else ctrl$correction
    alp <- if (corr == "custom" && !is.null(ctrl$alpha)) as.numeric(ctrl$alpha) else alpha_min
    list(correction = corr, alpha = alp)
  }
  c_pre <- norm_ctrl(preliminary_moment_control)
  c_fin <- norm_ctrl(final_moment_control)

  sum_err_x <- function(w) {
    out <- matrix(0.0, q, q)
    for (i in seq_len(nobs)) {
      out <- out + w[i] * err_x[, , i]
    }
    0.5 * (out + t(out))
  }
  build_err_moment <- function(ex_sum, ey_sum, exy_sum_design) {
    em <- matrix(0.0, p + 1L, p + 1L)
    em[1L, 1L] <- ey_sum
    em[1L, 2L:(p + 1L)] <- exy_sum_design
    em[2L:(p + 1L), 1L] <- exy_sum_design
    em[2L:(p + 1L), 2L:(p + 1L)] <- ex_sum
    em
  }

  ex_sum_1 <- sum_err_x(rep(1.0, nobs))
  err_cov_design_1 <- matrix(0.0, p, p)
  err_cov_design_1[pred_idx, pred_idx] <- ex_sum_1
  exy_sum_1 <- c(0.0, colSums(err_xy))
  err_mom_1 <- build_err_moment(err_cov_design_1, sum(err_y), exy_sum_1)
  obs_mom_1 <- crossprod(cbind(y, design))
  lam_1 <- medusa_smallest_det_root(obs_mom_1, err_mom_1)
  des_mom_1 <- crossprod(design)
  out_mom_1 <- as.numeric(crossprod(design, y))

  psd_latent_gram <- function(M_obs, E_cov, c_val) {
    # Subtract full unbiased error covariance (including off-diagonal shared predictor error),
    # project latent Gram matrix onto PSD cone if small-sample noise induces negative curvature,
    # and apply Fuller's (1 - c) modification strictly to the diagonal ridge so shared error
    # covariance between predictors never leaks across columns when c < 1.
    Sig_lat <- 0.5 * ((M_obs - E_cov) + t(M_obs - E_cov))
    if (p > 2L) {
      eg_s <- eigen(Sig_lat[pred_idx, pred_idx, drop = FALSE], symmetric = TRUE)
      if (min(eg_s$values) < 0.0) {
        Sig_lat[pred_idx, pred_idx] <- eg_s$vectors %*% (pmax(eg_s$values, 0.0) * t(eg_s$vectors))
      }
    }
    D_err <- diag(diag(E_cov), p)
    Sig_lat + (1.0 - c_val) * D_err
  }

  eval_fixed <- function(alpha_1, corr_1, alpha_2, corr_2) {
    c1 <- if (identical(corr_1, "none") || max(abs(err_mom_1)) <= .Machine$double.eps) {
      1.0
    } else {
      max(0.0, medusa_mmom_correction(lam_1, nobs, alpha_1))
    }
    S1_x <- psd_latent_gram(des_mom_1, err_cov_design_1, c1)
    s1_xy <- out_mom_1 - exy_sum_1
    coef0 <- tryCatch(as.numeric(qr.solve(S1_x, s1_xy)), error = function(e) NULL)
    if (is.null(coef0) || any(!is.finite(coef0))) {
      return(list(converged = FALSE, S1_x = S1_x, c1 = c1, lam_1 = lam_1))
    }
    names(coef0) <- colnames(design)
    g0_pred <- coef0[pred_idx]
    err_resid_var <- numeric(nobs)
    for (i in seq_len(nobs)) {
      err_resid_var[i] <- err_y[i] - 2.0 * sum(g0_pred * err_xy[i, ]) +
        as.numeric(crossprod(g0_pred, err_x[, , i] %*% g0_pred))
    }
    err_resid_var <- pmax(0.0, err_resid_var)
    resid1 <- y - as.numeric(design %*% coef0)
    sigma2_eps <- if (is.finite(lam_1) && lam_1 < 1.0) {
      0.0
    } else {
      max(0.0, sum(resid1^2) / max(1L, nobs - p) - mean(err_resid_var))
    }
    case_var <- sigma2_eps + err_resid_var
    if (any(!is.finite(case_var)) || any(case_var <= sqrt(.Machine$double.eps))) {
      return(list(converged = FALSE, S1_x = S1_x, c1 = c1, lam_1 = lam_1))
    }
    w_case <- 1.0 / case_var
    w_des <- design * sqrt(w_case)
    w_out <- y * sqrt(w_case)
    w_mom <- crossprod(cbind(w_out, w_des))
    ex_sum_2 <- sum_err_x(w_case)
    err_cov_design_2 <- matrix(0.0, p, p)
    err_cov_design_2[pred_idx, pred_idx] <- ex_sum_2
    exy_sum_2 <- c(0.0, colSums(err_xy * w_case))
    err_mom_2 <- build_err_moment(err_cov_design_2, sum(w_case * err_y), exy_sum_2)
    lam_2 <- medusa_smallest_det_root(w_mom, err_mom_2)
    c2 <- if (identical(corr_2, "none") || max(abs(err_mom_2)) <= .Machine$double.eps) {
      1.0
    } else {
      max(0.0, medusa_mmom_correction(lam_2, nobs, alpha_2))
    }
    S2_x <- psd_latent_gram(w_mom[2L:(p + 1L), 2L:(p + 1L), drop = FALSE], err_cov_design_2, c2)
    s2_xy <- w_mom[2L:(p + 1L), 1L] - exy_sum_2
    coef2 <- tryCatch(as.numeric(qr.solve(S2_x, s2_xy)), error = function(e) NULL)
    if (is.null(coef2) || any(!is.finite(coef2))) {
      return(list(converged = FALSE, S1_x = S1_x, S2_x = S2_x, c1 = c1, c2 = c2, lam_1 = lam_1, lam_2 = lam_2))
    }
    names(coef2) <- colnames(design)
    bread_inv <- tryCatch(qr.solve(S2_x / nobs, diag(p)), error = function(e) NULL)
    if (is.null(bread_inv) || any(!is.finite(bread_inv))) {
      return(list(converged = FALSE, S1_x = S1_x, S2_x = S2_x, c1 = c1, c2 = c2, lam_1 = lam_1, lam_2 = lam_2))
    }
    tilde_cp <- matrix(0.0, p, p)
    for (i in seq_len(nobs)) {
      om_i <- c(0.0, err_xy[i, ] - as.numeric(err_x[, , i] %*% g0_pred))
      tilde_cp <- tilde_cp + (w_case[i]^2) * tcrossprod(om_i)
    }
    meat <- crossprod(w_des) + tilde_cp
    vcov_mat <- 0.5 * ((bread_inv %*% meat %*% bread_inv) / (nobs^2) +
      t((bread_inv %*% meat %*% bread_inv) / (nobs^2)))
    dimnames(vcov_mat) <- list(colnames(design), colnames(design))
    if (any(!is.finite(diag(vcov_mat))) || any(diag(vcov_mat) < 0)) {
      return(list(converged = FALSE, S1_x = S1_x, S2_x = S2_x, c1 = c1, c2 = c2, lam_1 = lam_1, lam_2 = lam_2))
    }
    se_vec <- sqrt(diag(vcov_mat))
    z_vec <- coef2 / se_vec
    p_vec <- 2.0 * stats::pnorm(abs(z_vec), lower.tail = FALSE)
    coef_df <- data.frame(
      term = colnames(design),
      estimate = coef2,
      se = se_vec,
      zstat = z_vec,
      pval = p_vec,
      ci_low = coef2 - 1.959964 * se_vec,
      ci_high = coef2 + 1.959964 * se_vec,
      stringsAsFactors = FALSE
    )
    list(
      converged = TRUE,
      coefficients = coef_df,
      vcov = vcov_mat,
      preliminary_coef = coef0,
      corrected_design_moment = S2_x,
      analysis_weights = w_case,
      S1_x = S1_x,
      S2_x = S2_x,
      w_des_mom = crossprod(w_des),
      c1 = c1,
      c2 = c2,
      lam_1 = lam_1,
      lam_2 = lam_2,
      alpha_prelim = alpha_1,
      alpha_final = alpha_2,
      resid_var = sigma2_eps
    )
  }

  check_guard <- function(fit_cand, stage) {
    if (!isTRUE(fit_cand$converged)) return(FALSE)
    if (stage == "preliminary") {
      S_corr <- fit_cand$S1_x[pred_idx, pred_idx, drop = FALSE] / nobs
      S_obs <- des_mom_1[pred_idx, pred_idx, drop = FALSE] / nobs
    } else {
      S_corr <- fit_cand$S2_x[pred_idx, pred_idx, drop = FALSE]
      S_obs <- fit_cand$w_des_mom[pred_idx, pred_idx, drop = FALSE]
    }
    eig_c <- eigen(0.5 * (S_corr + t(S_corr)), symmetric = TRUE, only.values = TRUE)$values
    eig_o <- eigen(0.5 * (S_obs + t(S_obs)), symmetric = TRUE, only.values = TRUE)$values
    min_c <- min(eig_c)
    max_c <- max(eig_c)
    max_o <- max(eig_o)
    is.finite(min_c) && is.finite(max_c) && is.finite(max_o) &&
      (min_c >= sqrt(.Machine$double.eps)) &&
      (min_c >= 0.05 * max_o) &&
      ((max_c / max(min_c, .Machine$double.eps)) <= 1e5)
  }

  alpha_1 <- c_pre$alpha
  alpha_2 <- c_fin$alpha
  fit_cur <- eval_fixed(
    alpha_1, if (c_pre$correction == "step_up") "custom" else c_pre$correction,
    alpha_2, if (c_fin$correction == "step_up") "custom" else c_fin$correction
  )

  if (c_pre$correction == "step_up" && !check_guard(fit_cur, "preliminary")) {
    upper_1 <- if (is.finite(lam_1)) {
      if (lam_1 <= 1 + 1 / nobs) lam_1 * nobs - 1 else nobs
    } else alpha_min
    if (upper_1 > alpha_min) {
      grid_1 <- seq(alpha_min, upper_1, length.out = 9L)
      pass_hi <- NA_real_
      fail_lo <- alpha_min
      fit_hi <- NULL
      for (a_cand in grid_1[-1L]) {
        f_c <- eval_fixed(a_cand, "custom", alpha_2, if (c_fin$correction == "step_up") "custom" else c_fin$correction)
        if (check_guard(f_c, "preliminary")) {
          pass_hi <- a_cand
          fit_hi <- f_c
          break
        } else {
          fail_lo <- a_cand
        }
      }
      if (is.finite(pass_hi)) {
        lo <- fail_lo; hi <- pass_hi
        for (iter in seq_len(12L)) {
          if (hi - lo <= 0.005) break
          mid <- 0.5 * (lo + hi)
          f_m <- eval_fixed(mid, "custom", alpha_2, if (c_fin$correction == "step_up") "custom" else c_fin$correction)
          if (check_guard(f_m, "preliminary")) {
            hi <- mid; fit_hi <- f_m
          } else {
            lo <- mid
          }
        }
        alpha_1 <- hi
        fit_cur <- fit_hi
      }
    }
  }

  if (c_fin$correction == "step_up" && !check_guard(fit_cur, "final")) {
    lam_2_cur <- fit_cur$lam_2
    upper_2 <- if (!is.null(lam_2_cur) && is.finite(lam_2_cur)) {
      if (lam_2_cur <= 1 + 1 / nobs) lam_2_cur * nobs - 1 else nobs
    } else alpha_min
    if (upper_2 > alpha_min) {
      grid_2 <- seq(alpha_min, upper_2, length.out = 9L)
      pass_hi <- NA_real_
      fail_lo <- alpha_min
      fit_hi <- NULL
      for (a_cand in grid_2[-1L]) {
        f_c <- eval_fixed(alpha_1, if (c_pre$correction == "step_up") "custom" else c_pre$correction, a_cand, "custom")
        if (check_guard(f_c, "final")) {
          pass_hi <- a_cand
          fit_hi <- f_c
          break
        } else {
          fail_lo <- a_cand
        }
      }
      if (is.finite(pass_hi)) {
        lo <- fail_lo; hi <- pass_hi
        for (iter in seq_len(12L)) {
          if (hi - lo <= 0.005) break
          mid <- 0.5 * (lo + hi)
          f_m <- eval_fixed(alpha_1, if (c_pre$correction == "step_up") "custom" else c_pre$correction, mid, "custom")
          if (check_guard(f_m, "final")) {
            hi <- mid; fit_hi <- f_m
          } else {
            lo <- mid
          }
        }
        alpha_2 <- hi
        fit_cur <- fit_hi
      }
    }
  }

  fit_cur$diagnostics <- list(
    nobs = nobs,
    npredictors = q,
    lambda_prelim = fit_cur$lam_1,
    lambda_final = fit_cur$lam_2,
    alpha_prelim = alpha_1,
    alpha_final = alpha_2,
    correction_prelim = fit_cur$c1,
    correction_final = fit_cur$c2,
    resid_var = fit_cur$resid_var
  )
  fit_cur
}

#' Fit a REPAIR modified method-of-moments (M-MOM / M-MOMS) measurement-error model
#'
#' Estimates a linear regression model when predictors, the outcome, or both are
#' likelihood-only predictions measured with heteroscedastic sampling error
#' (Vig & Hallquist, 2026). Delegates directly to \code{repair::fit_mmom()} when
#' the \pkg{repair} package is installed and uses an identical built-in fallback
#' otherwise.
#'
#' @param outcome Numeric outcome vector of length \eqn{N}.
#' @param predictors Numeric predictor matrix (\eqn{N \times q}, excluding the
#'   intercept).
#' @param err_x Numeric 3D array (\eqn{q \times q \times N}) of per-observation
#'   predictor sampling-error covariance matrices \eqn{\hat{\boldsymbol{\Psi}}_{x,j}}.
#' @param err_y Non-negative numeric vector of length \eqn{N} giving per-observation
#'   outcome sampling-error variances \eqn{\hat{\psi}_{y,j}}.
#' @param err_xy Optional numeric matrix (\eqn{N \times q}) of per-observation
#'   outcome-predictor sampling-error covariances \eqn{\hat{\boldsymbol{\psi}}_{xy,j}^\top}.
#' @param regularization Regularization policy: \code{"step_up"} (default;
#'   adaptive M-MOMS step-up search), \code{"minimum"} (fixed minimum
#'   \eqn{\alpha = q + 2}), \code{"none"}, or \code{"custom"}.
#' @param alpha Optional custom regularization parameter when
#'   \code{regularization = "custom"}.
#' @param backend Computation backend passed to \code{repair::fit_mmom()} when
#'   available (\code{"auto"}, \code{"r"}, or \code{"rcpp"}).
#' @return A list with \code{converged}, \code{coefficients} (data.frame with
#'   \code{term}, \code{estimate}, \code{se}, \code{zstat}, \code{pval},
#'   \code{ci_low}, \code{ci_high}), \code{vcov}, \code{preliminary_coef},
#'   \code{corrected_design_moment}, \code{analysis_weights},
#'   \code{diagnostics}, and \code{repair_fit} (the raw \pkg{repair} object when
#'   \pkg{repair} is installed).
#' @export
medusa_fit_mmom <- function(
  outcome,
  predictors,
  err_x,
  err_y = NULL,
  err_xy = NULL,
  regularization = c("step_up", "minimum", "none", "custom"),
  alpha = NULL,
  backend = c("auto", "r", "rcpp")
) {
  regularization <- match.arg(regularization)
  backend <- match.arg(backend)
  predictors <- as.matrix(predictors)
  storage.mode(predictors) <- "double"
  nobs <- length(outcome)
  q <- ncol(predictors)
  if (is.null(colnames(predictors))) colnames(predictors) <- paste0("x", seq_len(q))
  if (is.null(err_y)) err_y <- rep(0.0, nobs)
  if (is.null(err_xy)) err_xy <- matrix(0.0, nobs, q)
  ctrl <- list(correction = regularization, alpha = alpha)
  data_list <- list(
    outcome = as.numeric(outcome),
    predictors = predictors,
    err_x = err_x,
    err_y = as.numeric(err_y),
    err_xy = as.matrix(err_xy)
  )

  has_offdiag_err_x <- FALSE
  if (q >= 2L) {
    ex_tot <- apply(abs(err_x), c(1L, 2L), sum)
    has_offdiag_err_x <- any(ex_tot[row(ex_tot) != col(ex_tot)] > 1e-12)
  }
  if (max(abs(err_xy)) <= 1e-12 && !has_offdiag_err_x && requireNamespace("repair", quietly = TRUE)) {
    rfit <- tryCatch(
      repair::fit_mmom(
        data_list = data_list,
        preliminary_moment_control = ctrl,
        final_moment_control = ctrl,
        backend = backend
      ),
      error = function(e) NULL
    )
    if (!is.null(rfit) && isTRUE(rfit@converged)) {
      diag_list <- rfit@diagnostics
      return(list(
        converged = TRUE,
        coefficients = as.data.frame(rfit@coefficients),
        vcov = as.matrix(rfit@vcov),
        preliminary_coef = rfit@preliminary_coef,
        corrected_design_moment = as.matrix(rfit@corrected_design_moment),
        analysis_weights = as.numeric(rfit@analysis_weights),
        diagnostics = list(
          nobs = nobs,
          npredictors = q,
          lambda_prelim = diag_list$lambda_prelim,
          lambda_final = diag_list$lambda_final,
          alpha_prelim = if (!is.null(diag_list$alpha_prelim)) diag_list$alpha_prelim else q + 2,
          alpha_final = if (!is.null(diag_list$alpha_final)) diag_list$alpha_final else q + 2
        ),
        repair_fit = rfit
      ))
    }
  }

  fb <- medusa_fit_mmom_fallback(
    data_list = data_list,
    preliminary_moment_control = ctrl,
    final_moment_control = ctrl
  )
  fb$repair_fit <- NULL
  fb
}

#' Measurement-error-corrected inter-region connectivity and PPI via REPAIR
#'
#' Estimates trial-by-trial seed-to-target neural connectivity, directed lagged
#' connectivity (controlling for baseline target activity at the seed lag), and
#' task-modulated psychophysiological interactions (PPI) from a
#' \code{\link{decon_epoch_regression}} trial panel using the REPAIR modified
#' method-of-moments (\code{M-MOMS} / \code{M-MOM}) framework of Vig & Hallquist
#' (2026). Propagates per-trial operator standard errors (\code{se_demixed}) into
#' exact heteroscedastic predictor and outcome error covariances,including the
#' rank-1 interaction error block for PPI (\code{seed * x}) and dual predictor
#' error variances for Granger-style auto-lag controls (\code{targ_t0}).
#'
#' @param panel_dt A \code{\link{decon_epoch_regression}} result list, a
#'   trial-level \code{data.frame}/\code{data.table} (including CSV, \code{.csv.gz},
#'   Parquet, or RDS panels where pipe-delimited list columns are automatically
#'   restored), a list of per-subject results/panels, or a file path (or character
#'   vector of file paths) readable by \code{medusa_read_table()}.
#' @param seed_region Character scalar naming the seed region/parcel/voxel.
#' @param target_region Character scalar naming the target region/parcel/voxel.
#' @param seed_lag Numeric lag (in seconds) for the seed readout (default \code{0}).
#' @param target_lag Numeric lag (in seconds) for the target readout (default \code{0}).
#' @param covariates Optional character vector of error-free trial covariate
#'   columns in \code{panel_dt}. Defaults to \code{"x"} when present in
#'   \code{panel_dt}, or \code{character(0)} otherwise.
#' @param test_ppi Logical; if \code{TRUE}, includes the interaction between
#'   the seed readout and the first covariate in \code{covariates}
#'   (\code{seed_c:x_c}) with exact per-trial interaction error covariance
#'   \eqn{\tilde{\psi}_{S,j} [1, x_{c,j}; x_{c,j}, x_{c,j}^2]}. Default: \code{FALSE}.
#' @param control_target_seed_lag Logical; if \code{TRUE} and
#'   \code{abs(target_lag - seed_lag) > 1e-6}, includes the target region's
#'   readout at \code{seed_lag} (\code{targ_t0_val}) as an error-prone control
#'   covariate with measurement-error variance \code{se_demixed^2}. Default:
#'   \code{FALSE}.
#' @param orthogonalize_seed_x Logical; if \code{TRUE} (default), residualizes
#'   the run-centered seed readout on the error-free covariates within each
#'   subject-run and scales the seed error variance by \eqn{1 - h_{x,r,j}}
#'   prior to forming PPI interactions, preventing quadratic task drive
#'   (\eqn{x_c^2}) from confounding physiological-by-psychological interactions.
#' @param method Estimation method: \code{"mmoms"} (default; adaptive step-up
#'   REPAIR \code{M-MOMS}), \code{"mmom"} (fixed minimum-regularization REPAIR
#'   \code{M-MOM}), \code{"eb_predictor"} (Skrondal & Laake 2001 role-specific
#'   Empirical-Bayes shrunk predictor on likelihood-only outcome),
#'   \code{"naive_gls"} (precision-weighted WLS without predictor error
#'   subtraction), or \code{"naive_unweighted"} (unweighted OLS).
#' @param use_whitened Logical; if \code{TRUE}, uses \code{decon_demixed_whitened}
#'   and \code{se_whitened} when present. Defaults to \code{FALSE} for
#'   \code{"mmoms"} and \code{"mmom"} (so elementwise PPI products are formed on
#'   unmixed trial states before error correction) and \code{TRUE} for
#'   \code{"naive_gls"}.
#' @param region_col Optional column name identifying regions/parcels/voxels in
#'   \code{panel_dt} (auto-detected from \code{"region"}, \code{"parcel"},
#'   \code{"voxel"}).
#' @param id_col Column name identifying subjects in \code{panel_dt} (default
#'   \code{"id"}).
#' @param shared_error Policy for cross-parcel shared observation-error covariance
#'   between seed and target readouts: \code{"auto"} (default; propagates the
#'   Stage-3 Adaptive-POET + dynamic covariance-regression
#'   \code{shared_error_model} when operator objects are supplied, or the
#'   embedded per-trial Block-PAF factor columns \code{shared_noise_vec},
#'   \code{lag_noise_vec}, and \code{trial_cor_vec}---falling back to stationary
#'   \code{unit_error_cor} / \code{shared_noise_loading} on legacy panels),
#'   \code{"zero"} (assumes zero cross-parcel noise correlation while still
#'   propagating same-parcel cross-lag error correlation when
#'   \code{control_target_seed_lag = TRUE}), or \code{"supplied"} (reads per-trial
#'   raw seed--target error covariance from \code{err_xy_col}).
#' @param err_xy_col Optional column name supplying per-trial raw seed--target error
#'   covariance when \code{shared_error = "supplied"}.
#' @param run_col Run identifier column name (default \code{"run"}).
#' @param trial_col Within-run trial identifier column name (default \code{"trial"}).
#' @param lag_col Within-trial lag column name (default \code{"evt_time"}).
#' @param y_col Optional readout column override.
#' @param se_col Optional trial standard-error column override.
#' @param eb_col Optional Empirical-Bayes trial readout column override.
#' @param modulator Optional alias for \code{covariates}.
#' @param return_augmented_panel Logical; if \code{TRUE}, returns the merged
#'   trial panel augmented with moment-preserving REPAIR predictor columns
#'   (\code{seed_repair}, \code{ppi_repair}, \code{repair_weight}) that reproduce
#'   subject-level \code{M-MOMS} slopes in weighted least squares / \code{lmer()}.
#' @param crosstrial_se Logical; if \code{TRUE} (default), computes cross-trial
#'   (\eqn{i \neq j}) hybrid operator + gated onset-kernel sandwich standard errors
#'   and Satterthwaite effective degrees of freedom (\code{df}) for single-subject
#'   (\eqn{N = 1}) and subject-level inference while keeping diagonal M-MOMS
#'   estimation weights.
#' @param control_seed_lags Optional numeric vector of other within-trial seed
#'   lags (in seconds, present in \code{panel_dt}) to partial out as
#'   measurement-error-prone control predictors alongside \code{seed_lag}.
#' @param multilag_control Logical; if \code{TRUE}, automatically partials out
#'   other seed lags present in a multi-lag panel (subject to a latent-variance
#'   identifiability guard so silent lags do not induce M-MOMS attenuation) and
#'   enables causal target auto-lag control (\code{control_target_seed_lag = TRUE})
#'   strictly when \code{target_lag > seed_lag}. Default: \code{FALSE}.
#' @return A list of class \code{"medusa_connectivity_repair"} containing
#'   \code{group} (group-level HKSJ random-effects meta-analytic inference table,
#'   or single-subject cross-trial Satterthwaite \eqn{t}-inference table when
#'   \eqn{N = 1}), \code{subjects} (per-subject estimates, SEs, Satterthwaite
#'   \code{df}, and REPAIR diagnostics), optional \code{panel_repair}, and
#'   \code{settings}.
#' @seealso \code{\link{decon_epoch_regression}}, \code{\link{aggregate_decon_panel}},
#'   \code{\link{medusa_fit_mmom}}, \code{\link{medusa_between_subject_repair}}
#' @export
decon_connectivity_repair <- function(
  panel_dt,
  seed_region,
  target_region,
  seed_lag = 0,
  target_lag = 0,
  covariates = NULL,
  test_ppi = FALSE,
  control_target_seed_lag = FALSE,
  orthogonalize_seed_x = TRUE,
  method = c("mmoms", "mmom", "eb_predictor", "naive_gls", "naive_unweighted"),
  shared_error = c("auto", "zero", "supplied"),
  err_xy_col = NULL,
  use_whitened = NULL,
  region_col = NULL,
  id_col = "id",
  run_col = "run",
  trial_col = "trial",
  lag_col = "evt_time",
  y_col = NULL,
  se_col = NULL,
  eb_col = NULL,
  modulator = NULL,
  return_augmented_panel = FALSE,
  crosstrial_se = TRUE,
  control_seed_lags = NULL,
  multilag_control = FALSE
) {
  ctrl_t0_missing <- missing(control_target_seed_lag)
  method <- match.arg(method)
  shared_error <- match.arg(shared_error)
  if (is.null(covariates) && !is.null(modulator)) {
    covariates <- modulator
  }

  # Unpack panel_dt (decon_epoch_regression result, list of per-subject results/panels, data.frame, or file path(s))
  subj_unit_cor_map <- list()
  subj_lag_cor_map <- list()
  subj_sem_map <- list()
  unit_cor_global <- NULL
  lag_cor_global <- NULL
  sem_global <- NULL
  raw_gain_global <- NULL
  demix_flag_global <- TRUE
  tr_ref <- 1.0

  if (is.character(panel_dt)) {
    if (length(panel_dt) == 1L) {
      panel_dt <- medusa_read_table(panel_dt)
    } else if (length(panel_dt) > 1L) {
      panel_dt <- lapply(panel_dt, medusa_read_table)
    } else {
      stop("Empty character vector supplied to `decon_connectivity_repair()`.")
    }
  }

  if (is.list(panel_dt) && !is.data.frame(panel_dt)) {
    if (!is.null(panel_dt$panel) && is.data.frame(panel_dt$panel)) {
      p0 <- panel_dt$panel
      unit_cor_global <- if (!is.null(panel_dt$unit_error_cor)) panel_dt$unit_error_cor else attr(p0, "unit_error_cor")
      lag_cor_global <- if (!is.null(panel_dt$lag_error_cor)) panel_dt$lag_error_cor else attr(p0, "lag_error_cor")
      sem_global <- if (!is.null(panel_dt$shared_error_model)) panel_dt$shared_error_model else attr(p0, "shared_error_model")
      raw_gain_global <- attr(p0, "raw_gain_by_lag_unit")
      if (!is.null(panel_dt$settings$demix_target_trials)) {
        demix_flag_global <- isTRUE(panel_dt$settings$demix_target_trials)
      } else if (!is.null(attr(p0, "demix_target_trials"))) {
        demix_flag_global <- isTRUE(attr(p0, "demix_target_trials"))
      }
      if (!is.null(attr(p0, "tr")) && is.finite(attr(p0, "tr"))) tr_ref <- as.numeric(attr(p0, "tr"))
      dt <- restore_panel_vec_cols(data.table::copy(data.table::as.data.table(p0)))
    } else if (length(panel_dt) >= 1L) {
      sub_names <- names(panel_dt)
      if (is.null(sub_names)) sub_names <- sprintf("sub%03d", seq_along(panel_dt))
      dt_parts <- vector("list", length(panel_dt))
      for (idx_p in seq_along(panel_dt)) {
        el <- panel_dt[[idx_p]]
        if (is.character(el) && length(el) == 1L) {
          el <- medusa_read_table(el)
        }
        if (is.list(el) && !is.data.frame(el) && !is.null(el$panel)) {
          p_i <- el$panel
          uc_i <- if (!is.null(el$unit_error_cor)) el$unit_error_cor else attr(p_i, "unit_error_cor")
          lc_i <- if (!is.null(el$lag_error_cor)) el$lag_error_cor else attr(p_i, "lag_error_cor")
          sem_i <- if (!is.null(el$shared_error_model)) el$shared_error_model else attr(p_i, "shared_error_model")
          if (!is.null(el$settings$demix_target_trials)) demix_flag_global <- isTRUE(el$settings$demix_target_trials)
        } else if (is.data.frame(el)) {
          p_i <- el
          uc_i <- attr(p_i, "unit_error_cor")
          lc_i <- attr(p_i, "lag_error_cor")
          sem_i <- attr(p_i, "shared_error_model")
          if (!is.null(attr(p_i, "demix_target_trials"))) demix_flag_global <- isTRUE(attr(p_i, "demix_target_trials"))
        } else {
          stop("Unsupported element in `panel_dt` list passed to `decon_connectivity_repair()`.")
        }
        if (!is.null(attr(p_i, "tr")) && is.finite(attr(p_i, "tr"))) tr_ref <- as.numeric(attr(p_i, "tr"))
        if (is.null(raw_gain_global) && !is.null(attr(p_i, "raw_gain_by_lag_unit"))) {
          raw_gain_global <- attr(p_i, "raw_gain_by_lag_unit")
        }
        dt_i <- restore_panel_vec_cols(data.table::copy(data.table::as.data.table(p_i)))
        sid_i <- if (id_col %in% names(dt_i)) as.character(dt_i[[id_col]][1L]) else sub_names[idx_p]
        if (!(id_col %in% names(dt_i))) dt_i[, (id_col) := sid_i]
        if (!is.null(uc_i)) subj_unit_cor_map[[sid_i]] <- uc_i
        if (!is.null(lc_i)) subj_lag_cor_map[[sid_i]] <- lc_i
        if (!is.null(sem_i)) subj_sem_map[[sid_i]] <- sem_i
        dt_parts[[idx_p]] <- dt_i
      }
      dt <- data.table::rbindlist(dt_parts, fill = TRUE)
    } else {
      stop("Empty list supplied to `decon_connectivity_repair()`.")
    }
  } else {
    unit_cor_global <- attr(panel_dt, "unit_error_cor")
    lag_cor_global <- attr(panel_dt, "lag_error_cor")
    sem_global <- attr(panel_dt, "shared_error_model")
    raw_gain_global <- attr(panel_dt, "raw_gain_by_lag_unit")
    if (!is.null(attr(panel_dt, "demix_target_trials"))) {
      demix_flag_global <- isTRUE(attr(panel_dt, "demix_target_trials"))
    }
    if (!is.null(attr(panel_dt, "tr")) && is.finite(attr(panel_dt, "tr"))) {
      tr_ref <- as.numeric(attr(panel_dt, "tr"))
    }
    dt <- restore_panel_vec_cols(data.table::copy(data.table::as.data.table(panel_dt)))
  }

  if (is.null(region_col)) {
    for (cand in c("region", "parcel", "voxel")) {
      if (cand %in% names(dt)) {
        region_col <- cand
        break
      }
    }
  }
  if (is.null(region_col) || !(region_col %in% names(dt))) {
    stop("Could not find region/parcel/voxel column in panel_dt")
  }
  if (!(id_col %in% names(dt))) dt[, (id_col) := "sub001"]
  if (is.data.frame(panel_dt) && length(unique(dt[[id_col]])) > 1L) {
    # Note: base `rbind()` on data.tables copies attributes from the first subject's table;
    # never apply a single subject's operator/correlation attribute across multiple subjects!
    sem_global <- NULL
    unit_cor_global <- NULL
    lag_cor_global <- NULL
  }
  if (!(run_col %in% names(dt))) dt[, (run_col) := "run1"]
  if (!(lag_col %in% names(dt))) {
    alt_lag <- intersect(c("evt_time", "lag"), names(dt))[1L]
    if (!is.na(alt_lag)) lag_col <- alt_lag else dt[, (lag_col) := 0.0]
  }
  if (!(trial_col %in% names(dt))) {
    dt[, (trial_col) := seq_len(.N), by = c(id_col, run_col, region_col, lag_col)]
  }
  if (identical(shared_error, "supplied")) {
    if (is.null(err_xy_col) || !(err_xy_col %in% names(dt))) {
      stop("`shared_error = 'supplied'` requires a valid `err_xy_col` present in `panel_dt`.")
    }
  }
  if (is.null(covariates)) {
    covariates <- if ("x" %in% names(dt)) "x" else character(0L)
  } else {
    covariates <- intersect(covariates, names(dt))
  }
  if (!is.null(modulator)) {
    mod_chr <- as.character(modulator)[1L]
    if (!(mod_chr %in% covariates) && (mod_chr %in% names(dt))) {
      covariates <- c(covariates, mod_chr)
    }
  }
  if (isTRUE(test_ppi) && !length(covariates)) {
    stop("test_ppi = TRUE requires at least one covariate column in panel_dt")
  }
  mod_idx <- 1L
  if (isTRUE(test_ppi) && !is.null(modulator)) {
    mod_chr <- as.character(modulator)[1L]
    mod_idx <- match(mod_chr, covariates)
    if (is.na(mod_idx)) {
      stop("modulator '", mod_chr, "' not found in covariates / panel_dt.")
    }
  }
  ppi_term_name <- if (isTRUE(test_ppi) && length(covariates) > 1L && !is.null(modulator) && as.character(modulator)[1L] != "x") {
    paste0("seed_c:", covariates[mod_idx], "_c")
  } else {
    "seed_c:x_c"
  }

  if (is.null(use_whitened)) {
    use_whitened <- identical(method, "naive_gls")
  }

  if (is.null(y_col) || !(y_col %in% names(dt))) {
    y_col <- if (isTRUE(use_whitened) && "decon_demixed_whitened" %in% names(dt)) {
      "decon_demixed_whitened"
    } else if ("decon_demixed" %in% names(dt)) {
      "decon_demixed"
    } else if ("beta_lss" %in% names(dt)) {
      "beta_lss"
    } else {
      "decon_raw"
    }
  }
  if (is.null(se_col) || !(se_col %in% names(dt))) {
    se_col <- if (isTRUE(use_whitened) && "se_whitened" %in% names(dt)) {
      "se_whitened"
    } else if ("se_demixed" %in% names(dt)) {
      "se_demixed"
    } else if (identical(y_col, "decon_raw") && "se_raw" %in% names(dt)) {
      "se_raw"
    } else {
      NULL
    }
  }
  if (is.null(eb_col) || !(eb_col %in% names(dt))) {
    eb_col <- if (isTRUE(use_whitened) && "decon_eb_whitened" %in% names(dt)) {
      "decon_eb_whitened"
    } else if ("decon_eb" %in% names(dt)) {
      "decon_eb"
    } else {
      NULL
    }
  }

  all_lags_panel <- sort(unique(as.numeric(dt[[lag_col]])))
  if (isTRUE(multilag_control)) {
    if (is.null(control_seed_lags) && length(all_lags_panel) > 1L) {
      control_seed_lags <- all_lags_panel[abs(all_lags_panel - seed_lag) > 1e-6]
    }
    if (ctrl_t0_missing) {
      control_target_seed_lag <- (target_lag - seed_lag > 1e-6)
    }
    if (target_lag <= seed_lag + 1e-6) {
      control_target_seed_lag <- FALSE
    }
  }
  csl_lags <- numeric(0L)
  if (!is.null(control_seed_lags) && length(control_seed_lags) > 0L) {
    for (cl_in in as.numeric(control_seed_lags)) {
      if (!is.finite(cl_in) || abs(cl_in - seed_lag) <= 1e-6) next
      hit_cl <- all_lags_panel[abs(all_lags_panel - cl_in) < 1e-6]
      if (length(hit_cl) > 0L && !any(abs(csl_lags - hit_cl[1L]) < 1e-6)) {
        csl_lags <- c(csl_lags, hit_cl[1L])
      }
    }
  }

  # Equalize gain ONLY when operating on un-demixed raw readouts (decon_demixed is already unit-gain restored)
  get_raw_gain_scale <- function(reg_name, lag_val) {
    if (!(identical(y_col, "decon_raw") || isFALSE(demix_flag_global))) return(1.0)
    if (is.null(raw_gain_global) || !is.matrix(raw_gain_global)) return(1.0)
    c_idx <- match(as.character(reg_name), colnames(raw_gain_global))
    r_idx <- which(abs(all_lags_panel - lag_val) < 1e-6)[1L]
    if (is.na(c_idx) || is.na(r_idx) || r_idx > nrow(raw_gain_global)) return(1.0)
    gv <- as.numeric(raw_gain_global[r_idx, c_idx])
    if (is.finite(gv) && gv > 0.05) gv else 1.0
  }

  seed_sub <- dt[as.character(dt[[region_col]]) == as.character(seed_region) &
    abs(as.numeric(dt[[lag_col]]) - seed_lag) < 1e-6]
  targ_sub <- dt[as.character(dt[[region_col]]) == as.character(target_region) &
    abs(as.numeric(dt[[lag_col]]) - target_lag) < 1e-6]
  if (!nrow(seed_sub)) stop("No rows match seed_region = '", seed_region, "' at seed_lag = ", seed_lag)
  if (!nrow(targ_sub)) stop("No rows match target_region = '", target_region, "' at target_lag = ", target_lag)

  g_seed_raw <- get_raw_gain_scale(seed_region, seed_lag)
  g_targ_raw <- get_raw_gain_scale(target_region, target_lag)
  g_t0_raw <- 1.0

  has_vec <- ("shared_noise_vec" %in% names(dt)) && !isTRUE(use_whitened)
  has_lvec <- ("lag_noise_vec" %in% names(dt)) && !isTRUE(use_whitened)
  has_tcor <- ("trial_cor_vec" %in% names(dt)) && !isTRUE(use_whitened)
  has_ons <- "onset" %in% names(dt)
  has_ac1 <- "error_ac1" %in% names(dt) && !isTRUE(use_whitened)

  seed_keep <- data.table::data.table(
    id = seed_sub[[id_col]],
    run = seed_sub[[run_col]],
    trial = seed_sub[[trial_col]],
    seed_val = as.numeric(seed_sub[[y_col]]) / g_seed_raw,
    se_seed = (if (!is.null(se_col)) as.numeric(seed_sub[[se_col]]) else rep(0.0, nrow(seed_sub))) / g_seed_raw,
    seed_eb = if (!is.null(eb_col)) as.numeric(seed_sub[[eb_col]]) / g_seed_raw else NA_real_,
    tau2_seed = if ("tau2_bio" %in% names(seed_sub)) as.numeric(seed_sub[["tau2_bio"]]) else NA_real_,
    load_seed = if ("shared_noise_loading" %in% names(seed_sub)) as.numeric(seed_sub[["shared_noise_loading"]]) else NA_real_,
    ac1_seed = if (has_ac1) as.numeric(seed_sub[["error_ac1"]]) else rep(0.0, nrow(seed_sub)),
    vec_seed = if (has_vec) lapply(seed_sub[["shared_noise_vec"]], function(v) as.numeric(v) / g_seed_raw) else replicate(nrow(seed_sub), numeric(0L), simplify = FALSE),
    lvec_seed = if (has_lvec) lapply(seed_sub[["lag_noise_vec"]], function(v) as.numeric(v) / g_seed_raw) else replicate(nrow(seed_sub), numeric(0L), simplify = FALSE),
    tcor_seed = if (has_tcor) seed_sub[["trial_cor_vec"]] else replicate(nrow(seed_sub), numeric(0L), simplify = FALSE)
  )
  for (cv in covariates) {
    cv_use <- if (isTRUE(use_whitened) && paste0(cv, "_whitened") %in% names(seed_sub)) {
      paste0(cv, "_whitened")
    } else {
      cv
    }
    seed_keep[[cv]] <- as.numeric(seed_sub[[cv_use]])
  }

  targ_keep <- data.table::data.table(
    id = targ_sub[[id_col]],
    run = targ_sub[[run_col]],
    trial = targ_sub[[trial_col]],
    onset_val = if (has_ons) as.numeric(targ_sub[["onset"]]) else as.numeric(targ_sub[[trial_col]]) * 4.5,
    targ_val = as.numeric(targ_sub[[y_col]]) / g_targ_raw,
    se_targ = (if (!is.null(se_col)) as.numeric(targ_sub[[se_col]]) else rep(0.0, nrow(targ_sub))) / g_targ_raw,
    w_targ = if ("precision_weight" %in% names(targ_sub)) as.numeric(targ_sub[["precision_weight"]]) else rep(1.0, nrow(targ_sub)),
    load_targ = if ("shared_noise_loading" %in% names(targ_sub)) as.numeric(targ_sub[["shared_noise_loading"]]) else NA_real_,
    ac1_targ = if (has_ac1) as.numeric(targ_sub[["error_ac1"]]) else rep(0.0, nrow(targ_sub)),
    vec_targ = if (has_vec) lapply(targ_sub[["shared_noise_vec"]], function(v) as.numeric(v) / g_targ_raw) else replicate(nrow(targ_sub), numeric(0L), simplify = FALSE),
    lvec_targ = if (has_lvec) lapply(targ_sub[["lag_noise_vec"]], function(v) as.numeric(v) / g_targ_raw) else replicate(nrow(targ_sub), numeric(0L), simplify = FALSE),
    tcor_targ = if (has_tcor) targ_sub[["trial_cor_vec"]] else replicate(nrow(targ_sub), numeric(0L), simplify = FALSE),
    err_xy_raw = if (identical(shared_error, "supplied")) {
      as.numeric(targ_sub[[err_xy_col]]) / (g_seed_raw * g_targ_raw)
    } else {
      rep(NA_real_, nrow(targ_sub))
    }
  )

  m_dt <- merge(targ_keep, seed_keep, by = c("id", "run", "trial"), sort = FALSE)
  has_t0 <- isTRUE(control_target_seed_lag) && abs(target_lag - seed_lag) > 1e-6
  if (has_t0) {
    t0_sub <- dt[as.character(dt[[region_col]]) == as.character(target_region) &
      abs(as.numeric(dt[[lag_col]]) - seed_lag) < 1e-6]
    if (nrow(t0_sub) > 0L) {
      g_t0_raw <- get_raw_gain_scale(target_region, seed_lag)
      t0_keep <- data.table::data.table(
        id = t0_sub[[id_col]],
        run = t0_sub[[run_col]],
        trial = t0_sub[[trial_col]],
        targ_t0_val = as.numeric(t0_sub[[y_col]]) / g_t0_raw,
        se_targ_t0 = (if (!is.null(se_col)) as.numeric(t0_sub[[se_col]]) else rep(0.0, nrow(t0_sub))) / g_t0_raw,
        targ_t0_eb = if (!is.null(eb_col)) as.numeric(t0_sub[[eb_col]]) / g_t0_raw else NA_real_,
        tau2_targ_t0 = if ("tau2_bio" %in% names(t0_sub)) as.numeric(t0_sub[["tau2_bio"]]) else NA_real_,
        ac1_t0 = if (has_ac1) as.numeric(t0_sub[["error_ac1"]]) else rep(0.0, nrow(t0_sub)),
        vec_t0 = if (has_vec) lapply(t0_sub[["shared_noise_vec"]], function(v) as.numeric(v) / g_t0_raw) else replicate(nrow(t0_sub), numeric(0L), simplify = FALSE),
        lvec_t0 = if (has_lvec) lapply(t0_sub[["lag_noise_vec"]], function(v) as.numeric(v) / g_t0_raw) else replicate(nrow(t0_sub), numeric(0L), simplify = FALSE),
        tcor_t0 = if (has_tcor) t0_sub[["trial_cor_vec"]] else replicate(nrow(t0_sub), numeric(0L), simplify = FALSE)
      )
      m_dt <- merge(m_dt, t0_keep, by = c("id", "run", "trial"), sort = FALSE)
    } else {
      has_t0 <- FALSE
    }
  }

  csl_lags_eff <- numeric(0L)
  g_csl_raw <- numeric(0L)
  if (length(csl_lags) > 0L) {
    for (cl_val in csl_lags) {
      csl_sub <- dt[as.character(dt[[region_col]]) == as.character(seed_region) &
        abs(as.numeric(dt[[lag_col]]) - cl_val) < 1e-6]
      if (nrow(csl_sub) > 0L) {
        k_csl <- length(csl_lags_eff) + 1L
        g_ck <- get_raw_gain_scale(seed_region, cl_val)
        csl_keep <- data.table::data.table(
          id = csl_sub[[id_col]],
          run = csl_sub[[run_col]],
          trial = csl_sub[[trial_col]]
        )
        csl_keep[[paste0("csl_val_", k_csl)]] <- as.numeric(csl_sub[[y_col]]) / g_ck
        csl_keep[[paste0("se_csl_", k_csl)]] <- (if (!is.null(se_col)) as.numeric(csl_sub[[se_col]]) else rep(0.0, nrow(csl_sub))) / g_ck
        csl_keep[[paste0("ac1_csl_", k_csl)]] <- if (has_ac1) as.numeric(csl_sub[["error_ac1"]]) else rep(0.0, nrow(csl_sub))
        csl_keep[[paste0("vec_csl_", k_csl)]] <- if (has_vec) lapply(csl_sub[["shared_noise_vec"]], function(v) as.numeric(v) / g_ck) else replicate(nrow(csl_sub), numeric(0L), simplify = FALSE)
        csl_keep[[paste0("lvec_csl_", k_csl)]] <- if (has_lvec) lapply(csl_sub[["lag_noise_vec"]], function(v) as.numeric(v) / g_ck) else replicate(nrow(csl_sub), numeric(0L), simplify = FALSE)
        csl_keep[[paste0("tcor_csl_", k_csl)]] <- if (has_tcor) csl_sub[["trial_cor_vec"]] else replicate(nrow(csl_sub), numeric(0L), simplify = FALSE)
        m_dt <- merge(m_dt, csl_keep, by = c("id", "run", "trial"), sort = FALSE)
        csl_lags_eff <- c(csl_lags_eff, cl_val)
        g_csl_raw <- c(g_csl_raw, g_ck)
      }
    }
  }
  n_csl <- length(csl_lags_eff)

  keep_finite <- is.finite(m_dt$seed_val) & is.finite(m_dt$targ_val)
  if (has_t0) keep_finite <- keep_finite & is.finite(m_dt$targ_t0_val)
  if (n_csl > 0L) {
    for (k_csl in seq_len(n_csl)) {
      keep_finite <- keep_finite & is.finite(m_dt[[paste0("csl_val_", k_csl)]])
    }
  }
  for (cv in covariates) keep_finite <- keep_finite & is.finite(m_dt[[cv]])
  m_dt <- m_dt[keep_finite]
  if (!nrow(m_dt)) {
    stop("No complete finite paired (seed, target) trials remain after filtering.")
  }
  data.table::setorderv(m_dt, c("id", "run", "trial"))

  m_dt[, targ_c := targ_val - mean(targ_val), by = c("id", "run")]
  m_dt[, seed_c_raw := seed_val - mean(seed_val), by = c("id", "run")]
  if (has_t0) m_dt[, targ_t0_c := targ_t0_val - mean(targ_t0_val), by = c("id", "run")]
  if (n_csl > 0L) {
    for (k_csl in seq_len(n_csl)) {
      vn_k <- paste0("csl_val_", k_csl)
      m_dt[, (paste0("csl_c_raw_", k_csl)) := get(vn_k) - mean(get(vn_k)), by = c("id", "run")]
    }
  }
  for (cv in covariates) {
    m_dt[, (paste0(cv, "_c")) := get(cv) - mean(get(cv)), by = c("id", "run")]
  }
  cov_c_names <- if (length(covariates)) paste0(covariates, "_c") else character(0L)

  # Helper to look up cross-lag operator error correlation between any two lags (lag_a, lag_b)
  ref_lag_cor_cache <- new.env(parent = emptyenv())
  get_cross_lag_cor <- function(sid_key, reg_key, lag_a = seed_lag, lag_b = target_lag) {
    if (abs(lag_b - lag_a) < 1e-6) return(1.0)
    lc_obj <- if (!is.null(subj_lag_cor_map[[sid_key]])) subj_lag_cor_map[[sid_key]] else lag_cor_global
    R_mat <- NULL
    if (is.matrix(lc_obj)) {
      R_mat <- lc_obj
    } else if (is.list(lc_obj) && !is.null(lc_obj[[as.character(reg_key)]])) {
      R_mat <- lc_obj[[as.character(reg_key)]]
    }
    if (is.matrix(R_mat) && nrow(R_mat) == ncol(R_mat)) {
      rn_num <- suppressWarnings(as.numeric(rownames(R_mat)))
      i1 <- if (length(rn_num) == nrow(R_mat) && all(is.finite(rn_num))) {
        which(abs(rn_num - lag_a) < 1e-6)[1L]
      } else if (nrow(R_mat) == length(all_lags_panel)) {
        which(abs(all_lags_panel - lag_a) < 1e-6)[1L]
      } else NA_integer_
      i2 <- if (length(rn_num) == nrow(R_mat) && all(is.finite(rn_num))) {
        which(abs(rn_num - lag_b) < 1e-6)[1L]
      } else if (nrow(R_mat) == length(all_lags_panel)) {
        which(abs(all_lags_panel - lag_b) < 1e-6)[1L]
      } else NA_integer_
      if (!is.na(i1) && !is.na(i2)) {
        rv <- as.numeric(R_mat[i1, i2])
        if (is.finite(rv)) return(min(0.98, max(-0.98, rv)))
      }
    }
    lags_pair <- sort(unique(c(lag_a, lag_b)))
    cache_key <- paste(format(lags_pair, digits = 8L), collapse = "_")
    if (!exists(cache_key, envir = ref_lag_cor_cache, inherits = FALSE)) {
      ref_ops_pair <- medusa_reference_lag_operators(lags = lags_pair, tr = tr_ref)
      assign(cache_key, as.numeric(ref_ops_pair$R_op[1L, 2L]), envir = ref_lag_cor_cache)
    }
    min(0.98, max(-0.98, get(cache_key, envir = ref_lag_cor_cache, inherits = FALSE)))
  }

  # Helper to look up zero-lag cross-parcel correlation fallback when neither sem nor shared_noise_vec is present
  get_cross_unit_cor <- function(sid_key, sub_df) {
    if (identical(shared_error, "zero")) return(0.0)
    if (identical(as.character(seed_region), as.character(target_region))) return(1.0)
    uc_obj <- if (!is.null(subj_unit_cor_map[[sid_key]])) subj_unit_cor_map[[sid_key]] else unit_cor_global
    if (is.matrix(uc_obj) &&
        (as.character(seed_region) %in% rownames(uc_obj)) &&
        (as.character(target_region) %in% colnames(uc_obj))) {
      rv <- as.numeric(uc_obj[as.character(seed_region), as.character(target_region)])
      if (is.finite(rv)) return(min(0.98, max(-0.98, rv)))
    }
    if ("load_seed" %in% names(sub_df) && "load_targ" %in% names(sub_df)) {
      prod_load <- sub_df$load_seed * sub_df$load_targ
      if (any(is.finite(prod_load))) {
        rv <- stats::median(prod_load[is.finite(prod_load)])
        if (is.finite(rv)) return(min(0.98, max(-0.98, rv)))
      }
    }
    0.0
  }

  build_stacked_trial_cor <- function(tcor_list, tr_nums, ons_vec, ac1_vec) {
    m_r <- length(tr_nums)
    R_tr <- diag(1.0, m_r)
    if (m_r <= 1L) return(R_tr)
    lens <- vapply(tcor_list, length, integer(1L))
    if (any(lens > 0L)) {
      for (a in seq_len(m_r - 1L)) {
        vc <- as.numeric(unlist(tcor_list[[a]], use.names = FALSE))
        if (!length(vc)) next
        for (b in seq.int(a + 1L, m_r)) {
          d_tr <- as.integer(round(tr_nums[b] - tr_nums[a]))
          if (d_tr >= 1L && d_tr <= length(vc) && is.finite(vc[d_tr])) {
            R_tr[a, b] <- R_tr[b, a] <- vc[d_tr]
          }
        }
      }
    } else if (any(is.finite(ac1_vec) & abs(ac1_vec) > 1e-6)) {
      r1_m <- min(0.85, max(-0.85, stats::median(ac1_vec[is.finite(ac1_vec)])))
      med_iti <- stats::median(diff(ons_vec))
      if (!is.finite(med_iti) || med_iti <= 0.5) med_iti <- 4.5
      for (a in seq_len(m_r - 1L)) {
        for (b in seq.int(a + 1L, m_r)) {
          d_tr <- as.integer(round(tr_nums[b] - tr_nums[a]))
          if (d_tr == 1L) {
            dt_s <- abs(ons_vec[b] - ons_vec[a])
            R_tr[a, b] <- R_tr[b, a] <- r1_m * exp(-max(0.0, dt_s - med_iti) / 4.5)
          }
        }
      }
    }
    .poet_project_cor_psd(R_tr, eig_floor = 1e-4)
  }

  harmonize_run_joint_psd <- function(J_raw, diag_target) {
    N_j <- nrow(J_raw)
    if (N_j == 0L) return(J_raw)
    J_sym <- 0.5 * (J_raw + t(J_raw))
    sd_raw <- sqrt(pmax(diag(J_sym), 1e-14))
    sd_tgt <- sqrt(pmax(diag_target, 0.0))
    R_j <- J_sym / outer(sd_raw, sd_raw)
    R_j[!is.finite(R_j)] <- 0.0
    R_j[row(R_j) != col(R_j)] <- pmin(0.98, pmax(-0.98, R_j[row(R_j) != col(R_j)]))
    diag(R_j) <- 1.0
    eg_j <- eigen(0.5 * (R_j + t(R_j)), symmetric = TRUE)
    R_psd <- eg_j$vectors %*% (pmax(eg_j$values, 0.0) * t(eg_j$vectors))
    d_r <- sqrt(pmax(diag(R_psd), 1e-12))
    R_psd <- R_psd / outer(d_r, d_r)
    diag(R_psd) <- 1.0
    J_out <- R_psd * outer(sd_tgt, sd_tgt)
    0.5 * (J_out + t(J_out))
  }

  subj_ids <- unique(m_dt$id)
  subj_rows <- vector("list", length(subj_ids))
  aug_panels <- if (isTRUE(return_augmented_panel)) vector("list", length(subj_ids)) else NULL

  cand_terms <- c(
    "seed_c",
    if (has_t0) "targ_t0_val" else character(0L),
    if (n_csl > 0L) sprintf("seed_lag_%s_c", format(csl_lags_eff, trim = TRUE)) else character(0L)
  )
  cand_regs <- c(
    as.character(seed_region),
    if (has_t0) as.character(target_region) else character(0L),
    rep(as.character(seed_region), n_csl)
  )
  cand_lags <- c(seed_lag, if (has_t0) seed_lag else numeric(0L), csl_lags_eff)
  cand_g_raw <- c(g_seed_raw, if (has_t0) g_t0_raw else numeric(0L), g_csl_raw)
  cand_use_Ps <- c(TRUE, if (has_t0) FALSE else logical(0L), rep(TRUE, n_csl))
  cand_screen <- c(FALSE, if (has_t0) isTRUE(multilag_control) else logical(0L), rep(TRUE, n_csl))
  cand_raw_cols <- c(
    "seed_c_raw",
    if (has_t0) "targ_t0_c" else character(0L),
    if (n_csl > 0L) paste0("csl_c_raw_", seq_len(n_csl)) else character(0L)
  )
  cand_se_cols <- c(
    "se_seed",
    if (has_t0) "se_targ_t0" else character(0L),
    if (n_csl > 0L) paste0("se_csl_", seq_len(n_csl)) else character(0L)
  )
  cand_ac1_cols <- c(
    "ac1_seed",
    if (has_t0) "ac1_t0" else character(0L),
    if (n_csl > 0L) paste0("ac1_csl_", seq_len(n_csl)) else character(0L)
  )
  cand_vec_cols <- c(
    "vec_seed",
    if (has_t0) "vec_t0" else character(0L),
    if (n_csl > 0L) paste0("vec_csl_", seq_len(n_csl)) else character(0L)
  )
  cand_lvec_cols <- c(
    "lvec_seed",
    if (has_t0) "lvec_t0" else character(0L),
    if (n_csl > 0L) paste0("lvec_csl_", seq_len(n_csl)) else character(0L)
  )
  cand_tcor_cols <- c(
    "tcor_seed",
    if (has_t0) "tcor_t0" else character(0L),
    if (n_csl > 0L) paste0("tcor_csl_", seq_len(n_csl)) else character(0L)
  )
  K_cand <- length(cand_terms)

  for (si in seq_along(subj_ids)) {
    sid <- subj_ids[si]
    sid_chr <- as.character(sid)
    sub_dt <- m_dt[m_dt$id == sid]
    n_i <- nrow(sub_dt)
    sem <- if (!is.null(subj_sem_map[[sid_chr]])) subj_sem_map[[sid_chr]] else sem_global
    can_use_sem <- !is.null(sem) && identical(y_col, "decon_demixed") && !isTRUE(use_whitened) &&
                   !identical(shared_error, "supplied")

    cand_mat <- matrix(0.0, nrow = n_i, ncol = K_cand, dimnames = list(NULL, cand_terms))
    psi_y <- numeric(n_i)
    psi_x_arr <- array(0.0, dim = c(K_cand, K_cand, n_i))
    psi_xy_mat <- matrix(0.0, nrow = n_i, ncol = K_cand)
    df_s_tot <- 0.0
    df_t0_tot <- 0.0
    u_runs_sub <- unique(sub_dt$run)
    run_idx_list <- vector("list", length(u_runs_sub))
    run_joint_list <- vector("list", length(u_runs_sub))

    r_ST_0_fb <- get_cross_unit_cor(sid_chr, sub_dt)

    for (ri_sub in seq_along(u_runs_sub)) {
      r_lbl <- u_runs_sub[ri_sub]
      idx_r <- which(sub_dt$run == r_lbl)
      run_idx_list[[ri_sub]] <- idx_r
      m_r <- length(idx_r)
      if (m_r <= 1L) {
        cand_mat[idx_r, ] <- 0.0
        next
      }

      M_0 <- matrix(1.0 / sqrt(m_r), nrow = m_r, ncol = 1L)
      Q_x <- NULL
      if (isTRUE(orthogonalize_seed_x) && length(cov_c_names) > 0L && m_r > length(cov_c_names) + 2L) {
        X_r <- as.matrix(sub_dt[idx_r, cov_c_names, with = FALSE])
        csd <- apply(X_r, 2L, stats::sd)
        X_act <- X_r[, is.finite(csd) & csd > 1e-10, drop = FALSE]
        if (ncol(X_act) > 0L && m_r > ncol(X_act) + 2L) {
          qr_x <- qr(X_act); rk_x <- qr_x$rank
          Q_x <- qr.Q(qr_x)[, seq_len(rk_x), drop = FALSE]
          M_s <- cbind(M_0, Q_x)
        } else {
          M_s <- M_0
        }
      } else {
        M_s <- M_0
      }
      for (k_c in seq_len(K_cand)) {
        v_raw_r <- sub_dt[[cand_raw_cols[k_c]]][idx_r]
        if (cand_use_Ps[k_c] && !is.null(Q_x)) {
          cand_mat[idx_r, k_c] <- v_raw_r - as.numeric(Q_x %*% crossprod(Q_x, v_raw_r))
        } else {
          cand_mat[idx_r, k_c] <- v_raw_r
        }
      }
      df_s_tot <- df_s_tot + max(0.0, m_r - ncol(M_s))
      df_t0_tot <- df_t0_tot + max(0.0, m_r - ncol(M_0))
      P_s <- diag(m_r) - tcrossprod(M_s)
      P_0 <- diag(m_r) - tcrossprod(M_0)

      used_sem <- FALSE
      if (can_use_sem && (method %in% c("mmoms", "mmom", "eb_predictor"))) {
        r_pos <- match(as.character(r_lbl), as.character(sem$run_ids))
        g_T <- match(as.character(target_region), sem$unit_ids)
        l_T <- which(abs(sem$lags - target_lag) < 1e-6)[1L]
        g_cand_sem <- match(cand_regs, sem$unit_ids)
        l_cand_sem <- vapply(cand_lags, function(lv) which(abs(sem$lags - lv) < 1e-6)[1L], integer(1L))
        if (!is.na(r_pos) && !is.na(g_T) && !is.na(l_T) &&
            !any(is.na(g_cand_sem)) && !any(is.na(l_cand_sem))) {
          rs_sem <- sem$runs[[r_pos]]
          tr_idx <- sub_dt$trial[idx_r]
          if (all(tr_idx >= 1L & tr_idx <= rs_sem$m_target)) {
            sT <- rs_sem$s_idx[g_T]
            s_cand_sem <- rs_sem$s_idx[g_cand_sem]
            if (!is.na(sT) && !any(is.na(s_cand_sem))) {
              n_blk <- 1L + K_cand
              U_blk <- vector("list", n_blk)
              V_blk <- vector("list", n_blk)
              g_blk <- c(g_T, g_cand_sem)
              U_blk[[1L]] <- (rs_sem$U_strata[[sT]][[l_T]][, tr_idx, drop = FALSE] %*% P_0) / g_targ_raw
              V_blk[[1L]] <- if (!is.null(rs_sem$V_obs_strata)) (rs_sem$V_obs_strata[[sT]][[l_T]][, tr_idx, drop = FALSE] %*% P_0) / g_targ_raw else NULL
              for (k_c in seq_len(K_cand)) {
                P_kc <- if (cand_use_Ps[k_c]) P_s else P_0
                sk <- s_cand_sem[k_c]; lk <- l_cand_sem[k_c]; gk_raw <- cand_g_raw[k_c]
                U_blk[[1L + k_c]] <- (rs_sem$U_strata[[sk]][[lk]][, tr_idx, drop = FALSE] %*% P_kc) / gk_raw
                V_blk[[1L + k_c]] <- if (!is.null(rs_sem$V_obs_strata)) (rs_sem$V_obs_strata[[sk]][[lk]][, tr_idx, drop = FALSE] %*% P_kc) / gk_raw else NULL
              }

              fl_mat <- matrix(0.0, nrow = m_r, ncol = n_blk)
              raw_arr <- array(0.0, dim = c(n_blk, n_blk, m_r))
              for (a_b in seq_len(n_blk)) {
                fl_mat[, a_b] <- 0.15 * pmax(0.0, medusa_eval_sem_bilinear(
                  U_blk[[a_b]], U_blk[[a_b]], V_blk[[a_b]], V_blk[[a_b]], rs_sem, sem, g_blk[a_b], g_blk[a_b], stat_floor = TRUE
                ))
                for (b_b in seq.int(a_b, n_blk)) {
                  v_ab <- if (identical(shared_error, "zero") && (g_blk[a_b] != g_blk[b_b])) {
                    rep(0.0, m_r)
                  } else {
                    medusa_eval_sem_bilinear(
                      U_blk[[a_b]], U_blk[[b_b]], V_blk[[a_b]], V_blk[[b_b]], rs_sem, sem, g_blk[a_b], g_blk[b_b]
                    )
                  }
                  raw_arr[a_b, b_b, ] <- v_ab
                  raw_arr[b_b, a_b, ] <- v_ab
                }
              }
              for (jj in seq_len(m_r)) {
                Cp <- .project_trial_err_psd(raw_arr[, , jj], var_floors = fl_mat[jj, ])
                psi_y[idx_r[jj]] <- Cp[1L, 1L]
                psi_xy_mat[idx_r[jj], ] <- Cp[1L, 2L:n_blk]
                psi_x_arr[, , idx_r[jj]] <- Cp[2L:n_blk, 2L:n_blk, drop = FALSE]
              }
              if (isTRUE(crosstrial_se) && (method %in% c("mmoms", "mmom"))) {
                J_r <- matrix(0.0, nrow = n_blk * m_r, ncol = n_blk * m_r)
                diag_tgt <- numeric(n_blk * m_r)
                diag_tgt[seq_len(m_r)] <- psi_y[idx_r]
                for (k_c in seq_len(K_cand)) {
                  diag_tgt[(k_c * m_r + 1L):((k_c + 1L) * m_r)] <- psi_x_arr[k_c, k_c, idx_r]
                }
                for (a_b in seq_len(n_blk)) {
                  ra <- ((a_b - 1L) * m_r + 1L):(a_b * m_r)
                  for (b_b in seq.int(a_b, n_blk)) {
                    rb <- ((b_b - 1L) * m_r + 1L):(b_b * m_r)
                    P_ab <- if (identical(shared_error, "zero") && (g_blk[a_b] != g_blk[b_b])) {
                      matrix(0.0, m_r, m_r)
                    } else {
                      medusa_eval_sem_bilinear_mat(
                        U_blk[[a_b]], U_blk[[b_b]], V_blk[[a_b]], V_blk[[b_b]], rs_sem, sem, g_blk[a_b], g_blk[b_b]
                      )
                    }
                    if (a_b != b_b) {
                      d_ab <- if (a_b == 1L) psi_xy_mat[idx_r, b_b - 1L] else psi_x_arr[a_b - 1L, b_b - 1L, idx_r]
                      diag(P_ab) <- d_ab
                      J_r[ra, rb] <- t(P_ab)
                      J_r[rb, ra] <- P_ab
                    } else {
                      J_r[ra, ra] <- P_ab
                    }
                  }
                }
                run_joint_list[[ri_sub]] <- harmonize_run_joint_psd(J_r, diag_tgt)
              }
              used_sem <- TRUE
            }
          }
        }
      }

      if (!used_sem && (method %in% c("mmoms", "mmom", "eb_predictor"))) {
        n_blk <- 1L + K_cand
        blk_regs <- c(as.character(target_region), cand_regs)
        blk_lags <- c(target_lag, cand_lags)
        blk_se_cols <- c("se_targ", cand_se_cols)
        blk_ac1_cols <- c("ac1_targ", cand_ac1_cols)
        blk_vec_cols <- c("vec_targ", cand_vec_cols)
        blk_lvec_cols <- c("lvec_targ", cand_lvec_cols)
        blk_tcor_cols <- c("tcor_targ", cand_tcor_cols)
        blk_use_Ps <- c(FALSE, cand_use_Ps)

        s2_blk <- matrix(0.0, nrow = m_r, ncol = n_blk)
        for (a_b in seq_len(n_blk)) {
          se_v <- sub_dt[[blk_se_cols[a_b]]][idx_r]
          s2_blk[, a_b] <- ifelse(is.finite(se_v), se_v^2, 0.0)
        }
        c_unproj_arr <- array(0.0, dim = c(n_blk, n_blk, m_r))
        raw_arr <- array(0.0, dim = c(n_blk, n_blk, m_r))

        for (a_b in seq_len(n_blk)) {
          M_a <- if (blk_use_Ps[a_b]) M_s else M_0
          c_unproj_arr[a_b, a_b, ] <- s2_blk[, a_b]
          raw_arr[a_b, a_b, ] <- pmax(0.0, .proj_cross_diag(s2_blk[, a_b], M_a, M_a))
          if (a_b < n_blk) {
            for (b_b in seq.int(a_b + 1L, n_blk)) {
              M_b <- if (blk_use_Ps[b_b]) M_s else M_0
              same_reg <- identical(blk_regs[a_b], blk_regs[b_b])
              same_lag <- abs(blk_lags[a_b] - blk_lags[b_b]) < 1e-6
              if (same_reg) {
                c_ab <- vapply(idx_r, function(ii) {
                  la <- as.numeric(unlist(sub_dt[[blk_lvec_cols[a_b]]][[ii]], use.names = FALSE))
                  lb <- as.numeric(unlist(sub_dt[[blk_lvec_cols[b_b]]][[ii]], use.names = FALSE))
                  if (length(la) > 0L && length(la) == length(lb)) return(sum(la * lb))
                  va <- as.numeric(unlist(sub_dt[[blk_vec_cols[a_b]]][[ii]], use.names = FALSE))
                  vb <- as.numeric(unlist(sub_dt[[blk_vec_cols[b_b]]][[ii]], use.names = FALSE))
                  if (length(va) > 0L && length(va) == length(vb)) sum(va * vb) else 0.0
                }, numeric(1L))
                if (all(abs(c_ab) < 1e-12)) {
                  r_lg <- get_cross_lag_cor(sid_chr, blk_regs[a_b], blk_lags[a_b], blk_lags[b_b])
                  c_ab <- r_lg * sqrt(pmax(0.0, s2_blk[, a_b] * s2_blk[, b_b]))
                }
              } else if (identical(shared_error, "zero")) {
                c_ab <- rep(0.0, m_r)
              } else if (a_b == 1L && b_b == 2L && identical(shared_error, "supplied")) {
                c_ab <- ifelse(is.finite(sub_dt$err_xy_raw[idx_r]), sub_dt$err_xy_raw[idx_r], 0.0)
              } else {
                has_ab_vec <- same_lag && has_vec &&
                  any(vapply(sub_dt[[blk_vec_cols[a_b]]][idx_r], length, integer(1L)) > 0L) &&
                  any(vapply(sub_dt[[blk_vec_cols[b_b]]][idx_r], length, integer(1L)) > 0L)
                if (has_ab_vec) {
                  c_ab <- vapply(idx_r, function(ii) {
                    va <- as.numeric(unlist(sub_dt[[blk_vec_cols[a_b]]][[ii]], use.names = FALSE))
                    vb <- as.numeric(unlist(sub_dt[[blk_vec_cols[b_b]]][[ii]], use.names = FALSE))
                    if (length(va) > 0L && length(va) == length(vb)) sum(va * vb) else 0.0
                  }, numeric(1L))
                } else {
                  r_lg_ab <- if (same_lag) {
                    1.0
                  } else {
                    0.5 * (
                      get_cross_lag_cor(sid_chr, blk_regs[a_b], blk_lags[a_b], blk_lags[b_b]) +
                      get_cross_lag_cor(sid_chr, blk_regs[b_b], blk_lags[a_b], blk_lags[b_b])
                    )
                  }
                  c_ab <- (r_ST_0_fb * r_lg_ab) * sqrt(pmax(0.0, s2_blk[, a_b] * s2_blk[, b_b]))
                }
              }
              c_unproj_arr[a_b, b_b, ] <- c_ab
              c_unproj_arr[b_b, a_b, ] <- c_ab
              p_ab <- .proj_cross_diag(c_ab, M_b, M_a)
              raw_arr[a_b, b_b, ] <- p_ab
              raw_arr[b_b, a_b, ] <- p_ab
            }
          }
        }
        for (jj in seq_len(m_r)) {
          Cp <- .project_trial_err_psd(raw_arr[, , jj])
          psi_y[idx_r[jj]] <- Cp[1L, 1L]
          psi_xy_mat[idx_r[jj], ] <- Cp[1L, 2L:n_blk]
          psi_x_arr[, , idx_r[jj]] <- Cp[2L:n_blk, 2L:n_blk, drop = FALSE]
        }
        if (isTRUE(crosstrial_se) && (method %in% c("mmoms", "mmom"))) {
          tr_n <- sub_dt$trial[idx_r]; ons_n <- sub_dt$onset_val[idx_r]
          R_tr_list <- vector("list", n_blk)
          sd_u_list <- vector("list", n_blk)
          P_proj_list <- vector("list", n_blk)
          for (a_b in seq_len(n_blk)) {
            R_tr_list[[a_b]] <- build_stacked_trial_cor(sub_dt[[blk_tcor_cols[a_b]]][idx_r], tr_n, ons_n, sub_dt[[blk_ac1_cols[a_b]]][idx_r])
            sd_u_list[[a_b]] <- sqrt(pmax(s2_blk[, a_b], 1e-12))
            P_proj_list[[a_b]] <- if (blk_use_Ps[a_b]) P_s else P_0
          }
          J_r <- matrix(0.0, nrow = n_blk * m_r, ncol = n_blk * m_r)
          diag_tgt <- numeric(n_blk * m_r)
          diag_tgt[seq_len(m_r)] <- psi_y[idx_r]
          for (k_c in seq_len(K_cand)) {
            diag_tgt[(k_c * m_r + 1L):((k_c + 1L) * m_r)] <- psi_x_arr[k_c, k_c, idx_r]
          }
          for (a_b in seq_len(n_blk)) {
            ra <- ((a_b - 1L) * m_r + 1L):(a_b * m_r)
            P_aa <- P_proj_list[[a_b]] %*% (outer(sd_u_list[[a_b]], sd_u_list[[a_b]]) * R_tr_list[[a_b]]) %*% P_proj_list[[a_b]]
            J_r[ra, ra] <- P_aa
            if (a_b < n_blk) {
              for (b_b in seq.int(a_b + 1L, n_blk)) {
                rb <- ((b_b - 1L) * m_r + 1L):(b_b * m_r)
                if (identical(shared_error, "zero") && !identical(blk_regs[a_b], blk_regs[b_b])) {
                  P_ba <- matrix(0.0, m_r, m_r)
                } else {
                  c_ba <- c_unproj_arr[b_b, a_b, ]
                  r_ba <- pmin(0.98, pmax(-0.98, c_ba / pmax(sd_u_list[[b_b]] * sd_u_list[[a_b]], 1e-12)))
                  P_ba_u <- outer(sd_u_list[[b_b]], sd_u_list[[a_b]]) *
                    (0.5 * outer(r_ba, rep(1.0, m_r)) + 0.5 * outer(rep(1.0, m_r), r_ba)) *
                    (0.5 * (R_tr_list[[b_b]] + R_tr_list[[a_b]]))
                  P_ba <- P_proj_list[[b_b]] %*% P_ba_u %*% P_proj_list[[a_b]]
                }
                d_ab <- if (a_b == 1L) psi_xy_mat[idx_r, b_b - 1L] else psi_x_arr[b_b - 1L, a_b - 1L, idx_r]
                diag(P_ba) <- d_ab
                J_r[rb, ra] <- P_ba
                J_r[ra, rb] <- t(P_ba)
              }
            }
          }
          run_joint_list[[ri_sub]] <- harmonize_run_joint_psd(J_r, diag_tgt)
        }
      }
    }

    seed_c_vec <- cand_mat[, 1L]
    psi_s <- psi_x_arr[1L, 1L, ]
    psi_sy <- psi_xy_mat[, 1L]

    # Screen nuisance control lags so pure-noise silent lags cannot induce M-MOMS step-up attenuation on seed_c
    keep_cand <- rep(TRUE, K_cand)
    obs_v1 <- sum(seed_c_vec^2)
    for (k_c in seq_len(K_cand)) {
      if (!cand_screen[k_c]) next
      obs_vk <- sum(cand_mat[, k_c]^2)
      err_vk <- sum(psi_x_arr[k_c, k_c, ])
      n_cur_pred <- sum(keep_cand[seq_len(k_c - 1L)]) + 1L + length(covariates) + as.integer(isTRUE(test_ppi))
      if (!is.finite(obs_vk) || !is.finite(err_vk) ||
          (n_i - n_cur_pred - 1L < 6L) ||
          ((obs_vk - err_vk) <= max(0.55 * err_vk, 0.10 * obs_v1))) {
        keep_cand[k_c] <- FALSE
      }
    }
    act_cand <- which(keep_cand)
    K_act <- length(act_cand)

    use_eb <- identical(method, "eb_predictor")
    cand_use_mat <- cand_mat[, act_cand, drop = FALSE]
    if (use_eb) {
      for (ka in seq_len(K_act)) {
        k_orig <- act_cand[ka]
        psi_k <- psi_x_arr[k_orig, k_orig, ]
        ok_k <- is.finite(psi_k) & (psi_k > 0)
        if (any(ok_k)) {
          df_k_eff <- max(1.0, if (cand_use_Ps[k_orig]) df_s_tot else df_t0_tot)
          fl_mult_k <- max(0.05, 2.0 / max(df_k_eff, 4.0))
          vk_raw <- cand_use_mat[, ka]
          t2_est_k <- max(fl_mult_k * stats::median(psi_k[ok_k]), (sum(vk_raw[ok_k]^2) - sum(psi_k[ok_k])) / df_k_eff)
          cand_use_mat[, ka] <- (t2_est_k / (t2_est_k + psi_k)) * vk_raw
        }
      }
    }

    pred_list <- list()
    for (ka in seq_len(K_act)) {
      pred_list[[cand_terms[act_cand[ka]]]] <- cand_use_mat[, ka]
    }
    s_vec <- cand_use_mat[, 1L]
    for (k in seq_along(covariates)) {
      pred_list[[paste0(covariates[k], "_c")]] <- sub_dt[[cov_c_names[k]]]
    }
    if (isTRUE(test_ppi)) {
      x_mod_c <- sub_dt[[cov_c_names[mod_idx]]]
      pred_list[[ppi_term_name]] <- s_vec * x_mod_c
    }
    X_i <- do.call(cbind, pred_list)
    q_i <- ncol(X_i)

    if (method %in% c("mmoms", "mmom")) {
      err_x_i <- array(0.0, dim = c(q_i, q_i, n_i))
      err_xy_i <- matrix(0.0, nrow = n_i, ncol = q_i)
      err_x_i[seq_len(K_act), seq_len(K_act), ] <- psi_x_arr[act_cand, act_cand, , drop = FALSE]
      err_xy_i[, seq_len(K_act)] <- psi_xy_mat[, act_cand, drop = FALSE]
      if (isTRUE(test_ppi)) {
        x_mod_c <- sub_dt[[cov_c_names[mod_idx]]]
        idx_ppi <- q_i
        for (ka in seq_len(K_act)) {
          err_x_i[ka, idx_ppi, ] <- x_mod_c * err_x_i[1L, ka, ]
          err_x_i[idx_ppi, ka, ] <- x_mod_c * err_x_i[1L, ka, ]
        }
        err_x_i[idx_ppi, idx_ppi, ] <- (x_mod_c^2) * psi_s
        err_xy_i[, idx_ppi] <- x_mod_c * psi_sy
      }

      reg_pol <- if (identical(method, "mmoms")) "step_up" else "minimum"
      mfit <- medusa_fit_mmom(
        outcome = sub_dt$targ_c,
        predictors = X_i,
        err_x = err_x_i,
        err_y = psi_y,
        err_xy = err_xy_i,
        regularization = reg_pol
      )
      r_eff_summary <- if (sum(psi_s) * sum(psi_y) > 1e-12) {
        min(0.98, max(-0.98, sum(psi_sy) / sqrt(sum(psi_s) * sum(psi_y))))
      } else 0.0

      if (isTRUE(mfit$converged)) {
        cf <- mfit$coefficients[mfit$coefficients$term != "(Intercept)", , drop = FALSE]
        se_out <- as.numeric(cf$se)
        stat_out <- as.numeric(cf$zstat)
        pval_out <- as.numeric(cf$pval)
        df_out <- rep(Inf, nrow(cf))

        # 5A: Cross-Trial (i != j) Hybrid Operator + Onset-Kernel Sandwich Variance & Satterthwaite df_eff
        if (isTRUE(crosstrial_se) && all(!vapply(run_joint_list, is.null, logical(1L)))) {
          ct_res <- tryCatch({
            p_full <- q_i + 1L
            des_full <- cbind(`(Intercept)` = 1.0, X_i)
            b_fin <- mfit$coefficients$estimate[match(colnames(des_full), mfit$coefficients$term)]
            b_pre <- as.numeric(mfit$preliminary_coef)
            g0_p <- b_pre[2:p_full]
            w_c <- as.numeric(mfit$analysis_weights)
            S2_raw <- as.matrix(mfit$corrected_design_moment)
            S2_sum <- S2_raw * (sum(w_c) / max(S2_raw[1L, 1L], 1e-12))
            B_inv <- qr.solve(S2_sum, diag(p_full))

            h_lev <- pmin(0.80, pmax(0.0, w_c * rowSums((des_full %*% B_inv) * des_full)))
            u_hat <- as.numeric(sub_dt$targ_c - des_full %*% b_fin)
            u_hc2 <- u_hat / sqrt(1.0 - h_lev)
            disp_2 <- mean(pmin(9.0, w_c * (u_hc2^2))) / 0.9707
            delta_s2_eps <- if (is.finite(disp_2) && disp_2 > 1.05) (disp_2 - 1.0) / max(mean(w_c), 1e-8) else 0.0

            n_runs_sub <- length(run_idx_list)
            Sig_u0_runs <- vector("list", n_runs_sub)
            C_om_runs <- vector("list", n_runs_sub)
            num_ex <- 0.0; den_var_ex <- 0.0; sum_k2 <- 0.0

            for (rk_i in seq_len(n_runs_sub)) {
              idx_r <- run_idx_list[[rk_i]]
              m_r <- length(idx_r)
              J_r_full <- run_joint_list[[rk_i]]
              if (K_act < K_cand) {
                keep_j_idx <- c(seq_len(m_r), unlist(lapply(act_cand, function(kc) (kc * m_r + 1L):((kc + 1L) * m_r))))
                J_r <- J_r_full[keep_j_idx, keep_j_idx, drop = FALSE]
              } else {
                J_r <- J_r_full
              }
              w_r <- w_c[idx_r]

              J_k_list <- vector("list", K_act)
              Jw_k_list <- vector("list", K_act)
              alp_k_list <- vector("list", K_act)
              for (ka in seq_len(K_act)) {
                J_ka <- matrix(0.0, m_r, p_full)
                J_ka[, 1L + ka] <- 1.0
                if (ka == 1L && isTRUE(test_ppi)) J_ka[, p_full] <- sub_dt[[cov_c_names[mod_idx]]][idx_r]
                J_k_list[[ka]] <- J_ka
                Jw_k_list[[ka]] <- J_ka * w_r
                alp_k_list[[ka]] <- as.numeric(J_ka[, 2:p_full, drop = FALSE] %*% g0_p)
              }
              L_nu <- do.call(cbind, c(list(diag(1.0, m_r)), lapply(alp_k_list, function(av) diag(-av, m_r))))
              Psi_nu <- L_nu %*% J_r %*% t(L_nu)
              Psi_nu <- 0.5 * (Psi_nu + t(Psi_nu))
              C_nu_list <- lapply(seq_len(K_act), function(ka) {
                J_r[(ka * m_r + 1L):((ka + 1L) * m_r), , drop = FALSE] %*% t(L_nu)
              })
              C_om <- matrix(0.0, p_full, p_full)
              for (ka in seq_len(K_act)) {
                for (kb in seq_len(K_act)) {
                  C_om <- C_om + crossprod(Jw_k_list[[ka]], (C_nu_list[[ka]] * t(C_nu_list[[kb]])) %*% Jw_k_list[[kb]])
                }
              }
              C_om_runs[[rk_i]] <- 0.5 * (C_om + t(C_om))

              s2_eps_r <- pmax(0.0, 1.0 / pmax(w_r, 1e-12) - diag(Psi_nu)) + delta_s2_eps
              Sig_u0 <- Psi_nu + diag(s2_eps_r, m_r)
              Sig_u0_runs[[rk_i]] <- Sig_u0

              if (m_r >= 2L) {
                ons_r <- sub_dt$onset_val[idx_r]
                u_r <- u_hc2[idx_r]
                sd_u0 <- sqrt(pmax(diag(Sig_u0), 1e-12))
                dt_m <- abs(outer(ons_r, ons_r, "-"))
                up_mask <- (row(dt_m) < col(dt_m)) & is.finite(dt_m) & (dt_m > 0.0) & (dt_m <= 24.0)
                if (any(up_mask)) {
                  k_v <- exp(-dt_m[up_mask] / 6.0)
                  R_u0 <- Sig_u0 / outer(sd_u0, sd_u0)
                  rex_v <- ((outer(u_r, u_r) - Sig_u0) / outer(sd_u0, sd_u0))[up_mask]
                  r0_v <- R_u0[up_mask]
                  num_ex <- num_ex + sum(k_v * rex_v)
                  den_var_ex <- den_var_ex + sum((k_v^2) * (1.0 + r0_v^2))
                  sum_k2 <- sum_k2 + sum(k_v^2)
                }
              }
            }

            z_eps <- if (den_var_ex > 1e-12) num_ex / sqrt(den_var_ex) else 0.0
            w_eps <- .poet_firm_gate(abs(z_eps), lo = 1.75, hi = 3.00)
            rho_eps <- if (sum_k2 > 1e-12) min(0.80, max(-0.60, num_ex / sum_k2)) else 0.0

            meat_tot <- matrix(0.0, p_full, p_full)
            meat_diag_ref <- crossprod(des_full * sqrt(w_c))
            tr_R2_sum <- 0.0

            for (rk_i in seq_len(n_runs_sub)) {
              idx_r <- run_idx_list[[rk_i]]
              m_r <- length(idx_r)
              w_r <- w_c[idx_r]
              D_r <- des_full[idx_r, , drop = FALSE]
              Sig_u <- Sig_u0_runs[[rk_i]]
              if (w_eps > 1e-6 && m_r >= 2L) {
                ons_r <- sub_dt$onset_val[idx_r]
                sd_u0 <- sqrt(pmax(diag(Sig_u), 1e-12))
                dt_m <- abs(outer(ons_r, ons_r, "-"))
                K_ex <- w_eps * rho_eps * exp(-dt_m / 6.0)
                diag(K_ex) <- 0.0; K_ex[!is.finite(K_ex)] <- 0.0
                Sig_try <- Sig_u + K_ex * outer(sd_u0, sd_u0)
                eg_u <- eigen(0.5 * (Sig_try + t(Sig_try)), symmetric = TRUE)
                Sig_u <- eg_u$vectors %*% (pmax(eg_u$values, 0.15 * min(sd_u0^2)) * t(eg_u$vectors))
              }
              sd_u <- sqrt(pmax(diag(Sig_u), 1e-12))
              R_u <- Sig_u / outer(sd_u, sd_u)
              tr_R2_sum <- tr_R2_sum + sum(R_u^2)
              G_r <- outer(w_r, w_r) * Sig_u
              meat_tot <- meat_tot + crossprod(D_r, G_r %*% D_r) + C_om_runs[[rk_i]]
            }

            meat_tot <- 0.5 * (meat_tot + t(meat_tot))
            eg_m <- eigen(meat_tot, symmetric = TRUE)
            floor_m <- 0.15 * min(pmax(eigen(meat_diag_ref, symmetric = TRUE, only.values = TRUE)$values, 1e-12))
            if (min(eg_m$values) < floor_m) {
              meat_tot <- eg_m$vectors %*% (pmax(eg_m$values, floor_m) * t(eg_m$vectors))
            }
            V_ct <- 0.5 * ((B_inv %*% meat_tot %*% B_inv) + t(B_inv %*% meat_tot %*% B_inv))
            dimnames(V_ct) <- list(colnames(des_full), colnames(des_full))
            n_eff_satt <- (as.numeric(n_i)^2) / max(as.numeric(n_i), tr_R2_sum)
            df_satt <- max(4.0, min(as.numeric(n_i - p_full), n_eff_satt - as.numeric(p_full)))
            list(V = V_ct, df = df_satt)
          }, error = function(e) NULL)

          if (!is.null(ct_res) && all(is.finite(diag(ct_res$V))) && all(diag(ct_res$V) > 0)) {
            idx_cf <- match(cf$term, colnames(ct_res$V))
            se_out <- sqrt(diag(ct_res$V)[idx_cf])
            df_out <- rep(ct_res$df, nrow(cf))
            stat_out <- as.numeric(cf$estimate) / se_out
            pval_out <- 2.0 * stats::pt(abs(stat_out), df = ct_res$df, lower.tail = FALSE)
          }
        }

        subj_rows[[si]] <- data.table::data.table(
          id = sid,
          term = cf$term,
          estimate = cf$estimate,
          std.error = se_out,
          statistic = stat_out,
          df = df_out,
          p.value = pval_out,
          n_trials = n_i,
          r_shared_error = r_eff_summary,
          lambda_prelim = mfit$diagnostics$lambda_prelim,
          lambda_final = mfit$diagnostics$lambda_final,
          alpha_prelim = mfit$diagnostics$alpha_prelim,
          alpha_final = mfit$diagnostics$alpha_final,
          converged = TRUE
        )

        if (isTRUE(return_augmented_panel)) {
          w_raw <- mfit$analysis_weights
          w_rep <- w_raw / mean(w_raw)
          S2_full <- mfit$corrected_design_moment
          S2_pred <- S2_full[-1L, -1L, drop = FALSE] -
            tcrossprod(S2_full[-1L, 1L], S2_full[1L, -1L]) / max(S2_full[1L, 1L], .Machine$double.eps)
          mu_X_w <- colSums(X_i * w_raw) / sum(w_raw)
          X_i_wcent <- sweep(X_i, 2L, mu_X_w, "-")
          X_w_obs <- crossprod(X_i_wcent, X_i_wcent * w_raw)
          b_mmom_vec <- as.numeric(cf$estimate)
          y_wcent <- sub_dt$targ_c - sum(w_raw * sub_dt$targ_c) / sum(w_raw)
          R_mat <- tryCatch(qr.solve(X_w_obs, S2_pred), error = function(e) diag(q_i))
          X_rep <- X_i_wcent %*% R_mat
          if (any(abs(err_xy_i) > 0)) {
            b_wls_vec <- tryCatch(
              as.numeric(qr.solve(X_w_obs, crossprod(X_i_wcent, w_raw * y_wcent))),
              error = function(e) NULL
            )
            b_norm2 <- sum(b_mmom_vec^2)
            if (!is.null(b_wls_vec) && all(is.finite(b_wls_vec)) && is.finite(b_norm2) && b_norm2 > 1e-12) {
              R_mat_adj <- R_mat + tcrossprod(b_wls_vec - as.numeric(R_mat %*% b_mmom_vec), b_mmom_vec) / b_norm2
              X_rep_try <- X_i_wcent %*% R_mat_adj
              b_check <- tryCatch(
                as.numeric(qr.solve(crossprod(X_rep_try, X_rep_try * w_raw), crossprod(X_rep_try, w_raw * y_wcent))),
                error = function(e) rep(NA_real_, q_i)
              )
              if (all(is.finite(b_check)) && max(abs(b_check - b_mmom_vec)) < 1e-6) {
                X_rep <- X_rep_try
              }
            }
          }
          colnames(X_rep) <- paste0(colnames(X_i), "_repair")
          sub_copy <- data.table::copy(sub_dt)
          data.table::set(sub_copy, j = "seed_c", value = seed_c_vec)
          for (col_j in seq_len(q_i)) {
            data.table::set(sub_copy, j = colnames(X_rep)[col_j], value = X_rep[, col_j])
          }
          data.table::set(sub_copy, j = "repair_weight", value = w_rep)
          aug_panels[[si]] <- sub_copy
        }
        next
      }
    }

    # Fallback / WLS / EB predictor / OLS path
    w_use <- if (identical(method, "naive_unweighted")) {
      rep(1.0, n_i)
    } else {
      w0 <- sub_dt$w_targ
      w0[!is.finite(w0) | w0 <= 0] <- 1.0
      w0 / mean(w0)
    }
    df_lm <- as.data.frame(X_i)
    df_lm$targ_c <- sub_dt$targ_c
    lm_fit <- stats::lm(targ_c ~ ., data = df_lm, weights = w_use)
    sm_cf <- summary(lm_fit)$coefficients
    keep_rows <- setdiff(rownames(sm_cf), "(Intercept)")
    clean_terms <- gsub("^`|`$", "", keep_rows)
    subj_rows[[si]] <- data.table::data.table(
      id = sid,
      term = clean_terms,
      estimate = as.numeric(sm_cf[keep_rows, "Estimate"]),
      std.error = as.numeric(sm_cf[keep_rows, "Std. Error"]),
      statistic = as.numeric(sm_cf[keep_rows, "t value"]),
      df = rep(as.numeric(lm_fit$df.residual), length(keep_rows)),
      p.value = as.numeric(sm_cf[keep_rows, "Pr(>|t|)"]),
      n_trials = n_i,
      r_shared_error = 0.0,
      lambda_prelim = NA_real_,
      lambda_final = NA_real_,
      alpha_prelim = NA_real_,
      alpha_final = NA_real_,
      converged = TRUE
    )
    if (isTRUE(return_augmented_panel)) {
      sub_copy <- data.table::copy(sub_dt)
      data.table::set(sub_copy, j = "seed_c", value = seed_c_vec)
      for (col_j in seq_len(q_i)) {
        data.table::set(sub_copy, j = paste0(colnames(X_i)[col_j], "_repair"), value = X_i[, col_j])
      }
      data.table::set(sub_copy, j = "repair_weight", value = w_use)
      aug_panels[[si]] <- sub_copy
    }
  }

  subj_dt <- data.table::rbindlist(subj_rows, fill = TRUE)

  # Group-level random-effects pooling with ancillary sample-size weights and unequal-weight HC2/HKSJ calibration
  terms_u <- unique(subj_dt$term)
  grp_rows <- vector("list", length(terms_u))
  for (ti in seq_along(terms_u)) {
    tm <- terms_u[ti]
    st <- subj_dt[term == tm & is.finite(estimate) & is.finite(std.error) & std.error > 0]
    S_k <- nrow(st)
    if (S_k == 0L) next
    if (S_k == 1L) {
      df_1 <- if ("df" %in% names(st) && is.finite(st$df[1L]) && st$df[1L] > 0) st$df[1L] else Inf
      t_crit_1 <- if (is.finite(df_1)) stats::qt(0.975, df = df_1) else 1.959964
      p_1 <- if (is.finite(df_1)) 2.0 * stats::pt(abs(st$statistic[1L]), df = df_1, lower.tail = FALSE) else st$p.value[1L]
      grp_rows[[ti]] <- data.table::data.table(
        term = tm,
        estimate = st$estimate[1L],
        std.error = st$std.error[1L],
        statistic = st$statistic[1L],
        df = df_1,
        p.value = p_1,
        conf.low = st$estimate[1L] - t_crit_1 * st$std.error[1L],
        conf.high = st$estimate[1L] + t_crit_1 * st$std.error[1L],
        tau2_between = 0.0,
        mean_unweighted = st$estimate[1L],
        median_subj = st$estimate[1L],
        n_subjects = 1L
      )
      next
    }
    theta_raw <- st$estimate
    theta_s <- theta_raw
    if (S_k >= 6L) {
      med_t <- stats::median(theta_raw)
      mad_t <- stats::mad(theta_raw, constant = 1.4826)
      if (is.finite(mad_t) && mad_t > 1e-8) {
        theta_s <- pmin(med_t + 3.0 * mad_t, pmax(med_t - 3.0 * mad_t, theta_raw))
      }
    }
    v_s <- st$std.error^2
    med_v <- stats::median(v_s)
    v_s_stab <- pmax(v_s, 0.15 * med_v)
    w0 <- 1.0 / v_s_stab
    mu0 <- sum(w0 * theta_s) / sum(w0)
    Q_stat <- sum(w0 * ((theta_s - mu0)^2))
    denom_dl <- sum(w0) - sum(w0^2) / sum(w0)
    tau2_dl <- if (denom_dl > 1e-12) max(0.0, (Q_stat - (S_k - 1.0)) / denom_dl) else 0.0

    w_anc <- if ("n_trials" %in% names(st) && all(is.finite(st$n_trials)) && all(st$n_trials > 0)) {
      as.numeric(st$n_trials)
    } else {
      rep(1.0, S_k)
    }
    w_anc <- w_anc / sum(w_anc)
    sum_w2 <- sum(w_anc^2)
    mu_re <- sum(w_anc * theta_s)
    denom_hksj <- pmin(1.0 - w_anc, pmax(1.0 - 2.0 * w_anc + sum_w2, 1.0 / max(2.0, S_k)))
    denom_hksj <- pmax(denom_hksj, 1.0 / max(2.0, 2.0 * S_k))
    var_hksj <- sum((w_anc^2) * ((theta_s - mu_re)^2) / denom_hksj)
    se_re <- sqrt(max(var_hksj, .Machine$double.eps))
    df_re <- max(1.0, (1.0 / max(sum_w2, .Machine$double.eps)) - 1.0)
    t_re <- mu_re / se_re
    p_re <- 2.0 * stats::pt(abs(t_re), df = df_re, lower.tail = FALSE)
    t_crit <- stats::qt(0.975, df = df_re)
    grp_rows[[ti]] <- data.table::data.table(
      term = tm,
      estimate = mu_re,
      std.error = se_re,
      statistic = t_re,
      df = df_re,
      p.value = p_re,
      conf.low = mu_re - t_crit * se_re,
      conf.high = mu_re + t_crit * se_re,
      tau2_between = tau2_dl,
      mean_unweighted = mean(theta_raw),
      median_subj = stats::median(theta_raw),
      n_subjects = S_k
    )
  }

  grp_dt <- data.table::rbindlist(grp_rows, fill = TRUE)
  if (nrow(grp_dt) > 0L) {
    grp_dt[, `:=`(
      seed_region = seed_region,
      target_region = target_region,
      seed_lag = seed_lag,
      target_lag = target_lag,
      method = method
    )]
  }

  structure(
    list(
      group = grp_dt,
      subjects = subj_dt,
      panel_repair = if (isTRUE(return_augmented_panel)) data.table::rbindlist(aug_panels, fill = TRUE) else NULL,
      settings = list(
        seed_region = seed_region,
        target_region = target_region,
        seed_lag = seed_lag,
        target_lag = target_lag,
        covariates = covariates,
        modulator = modulator,
        test_ppi = test_ppi,
        control_target_seed_lag = control_target_seed_lag,
        control_seed_lags = csl_lags_eff,
        multilag_control = multilag_control,
        orthogonalize_seed_x = orthogonalize_seed_x,
        crosstrial_se = crosstrial_se,
        method = method,
        shared_error = shared_error,
        use_whitened = use_whitened
      )
    ),
    class = "medusa_connectivity_repair"
  )
}

#' Between-subject REPAIR structural regression on MEDuSA mixed models
#'
#' Connects fitted \code{lme4::lmer()} (or \code{\link{mixed_by}} model objects)
#' from MEDuSA Stage 4 directly to \code{repair::repair()} (Vig & Hallquist,
#' 2026) when the \pkg{repair} package is installed, with an exact built-in
#' likelihood-only WLS + \code{\link{medusa_fit_mmom}} fallback when \pkg{repair}
#' is unavailable. Estimates measurement-error-corrected between-subject
#' regressions relating subject-specific neural random intercepts/slopes across
#' brain regions or to external behavioral/clinical variables.
#'
#' @param formula Structural regression formula using \code{re()} markers
#'   (e.g., \code{re("x_whitened", model = "target") ~ re("x_whitened", model = "seed") + age}).
#' @param models Named list of fitted \code{lme4::lmer} or \code{nlme::lme}
#'   models (or a single fitted model).
#' @param data Subject-level \code{data.frame} with one row per subject ID.
#' @param id Character scalar naming the subject ID column in \code{data}
#'   (default \code{"id"}).
#' @param regularization Regularization policy passed to \code{repair::repair()}
#'   or \code{\link{medusa_fit_mmom}}: \code{"step_up"} (default; adaptive
#'   M-MOMS), \code{"minimum"}, \code{"none"}, or \code{"custom"}.
#' @param cross_model_error Cross-model prediction-error covariance policy when
#'   \code{models} has multiple entries (defaults to \code{"independent"} when
#'   \code{length(models) > 1L}).
#' @param medusa_mlm Optional alias for \code{models} (e.g., a single fitted
#'   \code{lme4::lmer} model).
#' @param subject_data Optional alias for \code{data}.
#' @param id_col Optional alias for \code{id}.
#' @param ... Additional arguments passed to \code{repair::repair()} (or
#'   \code{alpha} for \code{\link{medusa_fit_mmom}}).
#' @return A list with \code{coefficients} (tidy \code{data.frame} with
#'   \code{term}, \code{estimate}, \code{std.error}, \code{statistic},
#'   \code{p.value}, \code{conf.low}, \code{conf.high}), \code{reliability},
#'   \code{predictions} (likelihood-only random-effect predictions and error
#'   variances), and \code{fit}.
#' @export
medusa_between_subject_repair <- function(
  formula,
  models = NULL,
  data = NULL,
  id = "id",
  regularization = c("step_up", "minimum", "none", "custom"),
  cross_model_error = NULL,
  medusa_mlm = NULL,
  subject_data = NULL,
  id_col = NULL,
  ...
) {
  if (is.null(models) && !is.null(medusa_mlm)) models <- medusa_mlm
  if (is.null(data) && !is.null(subject_data)) data <- subject_data
  if (!is.null(id_col)) id <- id_col
  if (is.null(models)) stop("`models` (or `medusa_mlm`) must be provided.")
  if (is.null(data)) stop("`data` (or `subject_data`) must be provided.")
  regularization <- match.arg(regularization)
  extra_dots <- list(...)

  if (requireNamespace("repair", quietly = TRUE)) {
    if (is.null(cross_model_error) && is.list(models) && length(models) > 1L) {
      cross_model_error <- "independent"
    }
    r_args <- list(
      formula = formula,
      models = models,
      data = as.data.frame(data),
      id = id,
      regularization = regularization
    )
    if (!is.null(cross_model_error)) r_args$cross_model_error <- cross_model_error
    if (length(extra_dots) > 0L) r_args <- c(r_args, extra_dots)
    r_mod <- do.call(repair::repair, r_args)
    return(list(
      coefficients = repair::repair_coefficients(r_mod),
      reliability = repair::repair_reliability(r_mod),
      predictions = repair::repair_predictions(r_mod),
      fit = r_mod
    ))
  }

  # Built-in fallback using exact lme4 subject-level likelihood-only WLS + medusa_fit_mmom()
  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Either `repair` or `lme4` is required for `medusa_between_subject_repair()`.")
  }
  df_sub <- as.data.frame(data)
  if (!(id %in% names(df_sub))) {
    stop("Subject ID column '", id, "' not found in `data`.")
  }
  subj_ids <- as.character(df_sub[[id]])
  is_single_mod <- inherits(models, c("merMod", "lmerMod")) || !is.list(models)
  has_named_models <- !is_single_mod && !is.null(names(models))
  mod_list <- if (is_single_mod) list(model = models) else models
  if (is.null(names(mod_list))) names(mod_list) <- paste0("model", seq_along(mod_list))

  re_cache <- list()
  extract_model_re <- function(mod_nm, grp_override = NULL) {
    key_c <- paste0(mod_nm, "||", if (is.null(grp_override)) "" else grp_override)
    if (!is.null(re_cache[[key_c]])) return(re_cache[[key_c]])
    if (!(mod_nm %in% names(mod_list))) {
      stop("Model '", mod_nm, "' referenced in `re()` was not found in `models`.")
    }
    mod <- mod_list[[mod_nm]]
    if (!inherits(mod, "merMod")) {
      stop("Built-in fallback for `medusa_between_subject_repair()` supports `lme4::lmer` (`merMod`) models; install `repair` for other model classes.")
    }
    flist <- lme4::getME(mod, "flist")
    cnms_list <- lme4::getME(mod, "cnms")
    zt_list <- lme4::getME(mod, "Ztlist")
    assign_vec <- attr(flist, "assign")
    grp_nm <- if (!is.null(grp_override) && (grp_override %in% names(flist))) {
      grp_override
    } else if (id %in% names(flist)) {
      id
    } else {
      names(flist)[1L]
    }
    grp_pos <- match(grp_nm, names(flist))
    bar_idx <- if (!is.null(assign_vec)) which(assign_vec == grp_pos) else which(names(cnms_list) == grp_nm)
    if (!length(bar_idx)) stop("Grouping factor '", grp_nm, "' not found in random effects of model '", mod_nm, "'.")

    zt_offsets <- cumsum(c(0L, lengths(cnms_list)))
    zt_cols <- list()
    re_terms <- character(0L)
    for (bi in bar_idx) {
      idx_zt <- seq.int(zt_offsets[bi] + 1L, zt_offsets[bi + 1L])
      zt_cols <- c(zt_cols, zt_list[idx_zt])
      re_terms <- c(re_terms, cnms_list[[bi]])
    }
    q_re <- length(re_terms)
    y_obs <- as.numeric(lme4::getME(mod, "y"))
    X_obs <- as.matrix(lme4::getME(mod, "X"))
    b_fix <- as.numeric(lme4::fixef(mod))
    r_obs <- y_obs - as.numeric(X_obs %*% b_fix)
    w_obs <- stats::weights(mod, type = "prior")
    if (is.null(w_obs) || length(w_obs) != length(y_obs)) w_obs <- rep(1.0, length(y_obs))
    w_obs <- pmax(as.numeric(w_obs), 1e-8)
    sig2 <- stats::sigma(mod)^2

    g_fac <- as.character(flist[[grp_nm]])
    g_levs <- levels(flist[[grp_nm]])
    n_s <- length(subj_ids)
    B_lik <- matrix(NA_real_, nrow = n_s, ncol = q_re, dimnames = list(subj_ids, re_terms))
    Psi_arr <- array(NA_real_, dim = c(q_re, q_re, n_s), dimnames = list(re_terms, re_terms, subj_ids))

    for (si in seq_len(n_s)) {
      sid <- subj_ids[si]
      lev_i <- match(sid, g_levs)
      idx_s <- which(g_fac == sid)
      if (is.na(lev_i) || length(idx_s) < q_re) next
      Z_s <- vapply(seq_len(q_re), function(j) as.numeric(zt_cols[[j]][lev_i, idx_s]), numeric(length(idx_s)))
      if (!is.matrix(Z_s)) Z_s <- matrix(Z_s, nrow = length(idx_s), ncol = q_re)
      w_s <- w_obs[idx_s]
      ZWZ <- crossprod(Z_s, Z_s * w_s)
      ZWZ_inv <- tryCatch(
        solve(ZWZ),
        error = function(e) tryCatch(qr.solve(ZWZ + 1e-8 * diag(q_re), diag(q_re)), error = function(e2) NULL)
      )
      if (is.null(ZWZ_inv) || any(!is.finite(ZWZ_inv))) next
      ZWZ_inv <- 0.5 * (ZWZ_inv + t(ZWZ_inv))
      B_lik[si, ] <- as.numeric(ZWZ_inv %*% crossprod(Z_s, w_s * r_obs[idx_s]))
      Psi_arr[, , si] <- sig2 * ZWZ_inv
    }

    vc_list <- lme4::VarCorr(mod)
    tau2_reml <- stats::setNames(rep(NA_real_, q_re), re_terms)
    for (vc_el in vc_list[bar_idx]) {
      sd_el <- attr(vc_el, "stddev")
      if (!is.null(sd_el)) {
        nm_hit <- intersect(names(sd_el), re_terms)
        tau2_reml[nm_hit] <- as.numeric(sd_el[nm_hit])^2
      }
    }

    res_obj <- list(B_lik = B_lik, Psi_arr = Psi_arr, re_terms = re_terms, tau2_reml = tau2_reml, grp_nm = grp_nm)
    re_cache[[key_c]] <<- res_obj
    res_obj
  }

  f_env <- environment(formula)
  if (is.null(f_env)) f_env <- parent.frame()
  parse_re_call <- function(expr) {
    mc <- match.call(function(term, model = NULL, group = NULL) NULL, expr)
    t_arg <- mc$term
    t_str <- if (is.character(t_arg)) {
      t_arg[1L]
    } else if (identical(deparse(t_arg), "1")) {
      "(Intercept)"
    } else {
      gsub("^\"|\"$", "", deparse(t_arg))
    }
    m_nm <- if (!is.null(mc$model)) as.character(eval(mc$model, f_env))[1L] else names(mod_list)[1L]
    g_nm <- if (!is.null(mc$group)) as.character(eval(mc$group, f_env))[1L] else NULL
    lbl <- if (!is.null(mc$model) || has_named_models) {
      sprintf("re(%s, model = \"%s\")", t_str, m_nm)
    } else {
      sprintf("re(%s)", t_str)
    }
    list(term = t_str, model = m_nm, group = g_nm, label = lbl)
  }

  if (length(formula) < 3L) stop("`formula` must be a two-sided formula (`lhs ~ rhs`).")
  lhs_expr <- formula[[2L]]
  rhs_expr <- formula[[3L]]

  lhs_re <- NULL
  if (is.call(lhs_expr) && identical(as.character(lhs_expr[[1L]]), "re")) {
    lhs_re <- parse_re_call(lhs_expr)
  }

  rhs_re_specs <- list()
  rewrite_rhs <- function(expr) {
    if (is.call(expr)) {
      if (identical(as.character(expr[[1L]]), "re")) {
        sp <- parse_re_call(expr)
        idx_sp <- length(rhs_re_specs) + 1L
        col_syn <- sprintf(".medusa_re_rhs_%d", idx_sp)
        sp$syn_col <- col_syn
        rhs_re_specs[[idx_sp]] <<- sp
        return(as.name(col_syn))
      }
      for (k in seq_along(expr)[-1L]) {
        expr[[k]] <- rewrite_rhs(expr[[k]])
      }
    }
    expr
  }
  rhs_clean <- rewrite_rhs(rhs_expr)

  df_work <- df_sub
  pred_out <- data.frame(id = subj_ids, stringsAsFactors = FALSE)
  names(pred_out)[1L] <- id
  rel_rows <- list()

  resolve_re_term_idx <- function(ext, t_str, m_nm) {
    j <- match(t_str, ext$re_terms)
    if (is.na(j) && t_str %in% c("1", "intercept", "Intercept")) {
      j <- match("(Intercept)", ext$re_terms)
    }
    if (is.na(j)) {
      stop("Random-effect term '", t_str, "' not found in model '", m_nm,
           "' (available: ", paste(ext$re_terms, collapse = ", "), ").")
    }
    j
  }

  if (!is.null(lhs_re)) {
    ext_y <- extract_model_re(lhs_re$model, lhs_re$group)
    jy <- resolve_re_term_idx(ext_y, lhs_re$term, lhs_re$model)
    y_all <- as.numeric(ext_y$B_lik[, jy])
    ey_all <- as.numeric(ext_y$Psi_arr[jy, jy, ])
    df_work$.medusa_re_lhs <- y_all
    pred_out[[paste0(lhs_re$model, "_", lhs_re$term, "_est")]] <- y_all
    pred_out[[paste0(lhs_re$model, "_", lhs_re$term, "_err_var")]] <- ey_all
    t2_y <- ext_y$tau2_reml[jy]
    if (!is.finite(t2_y)) t2_y <- max(0.0, stats::var(y_all, na.rm = TRUE) - mean(ey_all, na.rm = TRUE))
    rel_rows[[length(rel_rows) + 1L]] <- data.frame(
      role = "outcome", model = lhs_re$model, term = lhs_re$term,
      tau2 = t2_y, mean_err_var = mean(ey_all, na.rm = TRUE),
      reliability = t2_y / max(t2_y + mean(ey_all, na.rm = TRUE), 1e-12),
      stringsAsFactors = FALSE
    )
    f_work <- stats::as.formula(call("~", as.name(".medusa_re_lhs"), rhs_clean), env = f_env)
  } else {
    y_all <- as.numeric(eval(lhs_expr, df_work, f_env))
    ey_all <- rep(0.0, nrow(df_work))
    f_work <- stats::as.formula(call("~", lhs_expr, rhs_clean), env = f_env)
  }

  rhs_ext_list <- vector("list", length(rhs_re_specs))
  rhs_j_idx <- integer(length(rhs_re_specs))
  for (ki in seq_along(rhs_re_specs)) {
    sp <- rhs_re_specs[[ki]]
    ext_k <- extract_model_re(sp$model, sp$group)
    jk <- resolve_re_term_idx(ext_k, sp$term, sp$model)
    rhs_ext_list[[ki]] <- ext_k
    rhs_j_idx[ki] <- jk
    xk_all <- as.numeric(ext_k$B_lik[, jk])
    exk_all <- as.numeric(ext_k$Psi_arr[jk, jk, ])
    df_work[[sp$syn_col]] <- xk_all
    pred_out[[paste0(sp$model, "_", sp$term, "_est")]] <- xk_all
    pred_out[[paste0(sp$model, "_", sp$term, "_err_var")]] <- exk_all
    t2_k <- ext_k$tau2_reml[jk]
    if (!is.finite(t2_k)) t2_k <- max(0.0, stats::var(xk_all, na.rm = TRUE) - mean(exk_all, na.rm = TRUE))
    rel_rows[[length(rel_rows) + 1L]] <- data.frame(
      role = "predictor", model = sp$model, term = sp$term,
      tau2 = t2_k, mean_err_var = mean(exk_all, na.rm = TRUE),
      reliability = t2_k / max(t2_k + mean(exk_all, na.rm = TRUE), 1e-12),
      stringsAsFactors = FALSE
    )
  }

  mf <- stats::model.frame(f_work, data = df_work, na.action = stats::na.pass)
  mm_full <- stats::model.matrix(f_work, data = mf)
  keep_rows <- stats::complete.cases(mf) & (rowSums(!is.finite(mm_full)) == 0L) &
               is.finite(y_all) & is.finite(ey_all)
  for (ki in seq_along(rhs_re_specs)) {
    jk <- rhs_j_idx[ki]
    keep_rows <- keep_rows & is.finite(rhs_ext_list[[ki]]$Psi_arr[jk, jk, ])
  }
  idx_keep <- which(keep_rows)
  if (length(idx_keep) < 4L) {
    stop("Fewer than 4 complete finite subjects remain for `medusa_between_subject_repair()`.")
  }

  mm_keep <- mm_full[idx_keep, , drop = FALSE]
  int_col <- match("(Intercept)", colnames(mm_keep))
  X_pred <- if (!is.na(int_col)) mm_keep[, -int_col, drop = FALSE] else mm_keep
  q_pred <- ncol(X_pred)
  n_keep <- length(idx_keep)
  if (q_pred == 0L) stop("`formula` must include at least one predictor.")

  err_x <- array(0.0, dim = c(q_pred, q_pred, n_keep))
  err_xy <- matrix(0.0, nrow = n_keep, ncol = q_pred)
  syn_names <- vapply(rhs_re_specs, `[[`, character(1L), "syn_col")
  syn_labels <- vapply(rhs_re_specs, `[[`, character(1L), "label")

  pred_col_names <- colnames(X_pred)
  for (ki in seq_along(rhs_re_specs)) {
    c_ki <- match(syn_names[ki], pred_col_names)
    if (is.na(c_ki)) next
    ext_ki <- rhs_ext_list[[ki]]
    j_ki <- rhs_j_idx[ki]
    err_x[c_ki, c_ki, ] <- as.numeric(ext_ki$Psi_arr[j_ki, j_ki, idx_keep])
    if (ki < length(rhs_re_specs)) {
      for (kj in seq.int(ki + 1L, length(rhs_re_specs))) {
        c_kj <- match(syn_names[kj], pred_col_names)
        if (is.na(c_kj)) next
        if (identical(rhs_re_specs[[ki]]$model, rhs_re_specs[[kj]]$model) &&
            identical(ext_ki$grp_nm, rhs_ext_list[[kj]]$grp_nm)) {
          cov_ij <- as.numeric(ext_ki$Psi_arr[j_ki, rhs_j_idx[kj], idx_keep])
          err_x[c_ki, c_kj, ] <- cov_ij
          err_x[c_kj, c_ki, ] <- cov_ij
        }
      }
    }
    if (!is.null(lhs_re) && !identical(cross_model_error, "independent") &&
        identical(lhs_re$model, rhs_re_specs[[ki]]$model)) {
      ext_y <- extract_model_re(lhs_re$model, lhs_re$group)
      jy <- resolve_re_term_idx(ext_y, lhs_re$term, lhs_re$model)
      if (identical(ext_y$grp_nm, ext_ki$grp_nm)) {
        err_xy[, c_ki] <- as.numeric(ext_y$Psi_arr[jy, j_ki, idx_keep])
      }
    }
    pred_col_names[c_ki] <- syn_labels[ki]
  }
  colnames(X_pred) <- pred_col_names

  mfit <- medusa_fit_mmom(
    outcome = y_all[idx_keep],
    predictors = X_pred,
    err_x = err_x,
    err_y = ey_all[idx_keep],
    err_xy = err_xy,
    regularization = regularization,
    alpha = extra_dots$alpha
  )
  if (!isTRUE(mfit$converged)) {
    stop("`medusa_fit_mmom()` failed to converge in `medusa_between_subject_repair()`.")
  }

  cf_raw <- mfit$coefficients
  df_l2 <- max(1.0, as.numeric(n_keep - q_pred - 1L))
  t_crit <- stats::qt(0.975, df = df_l2)
  est_v <- as.numeric(cf_raw$estimate)
  se_v <- as.numeric(cf_raw$se)
  stat_v <- est_v / pmax(se_v, .Machine$double.eps)
  pval_v <- 2.0 * stats::pt(abs(stat_v), df = df_l2, lower.tail = FALSE)
  coef_out <- data.frame(
    term = as.character(cf_raw$term),
    estimate = est_v,
    std.error = se_v,
    statistic = stat_v,
    df = rep(df_l2, nrow(cf_raw)),
    p.value = pval_v,
    conf.low = est_v - t_crit * se_v,
    conf.high = est_v + t_crit * se_v,
    stringsAsFactors = FALSE
  )
  rel_df <- if (length(rel_rows) > 0L) do.call(rbind, rel_rows) else data.frame()

  list(
    coefficients = coef_out,
    reliability = rel_df,
    predictions = pred_out,
    fit = mfit
  )
}

# Build canonical reference cross-lag resolution matrix C_sym and operator error
# correlation R_op when only a lag grid (and not a decon_epoch_regression object)
# is supplied to Stage-4 multi-lag inference.
medusa_reference_lag_operators <- function(
  lags,
  tr = 1.0,
  readout_kernel = "triangle",
  readout_width = 1.5,
  lambda = 0.8,
  rho = 0.30,
  penalty = "sobolev",
  gain_floor = 0.15
) {
  lags <- as.numeric(lags)
  n_lags <- length(lags)
  lag_names <- as.character(lags)
  if (n_lags < 2L) {
    I1 <- diag(1, max(1L, n_lags))
    if (n_lags == 1L) dimnames(I1) <- list(lag_names, lag_names)
    return(list(C_sym = I1, R_op = I1, d_lag = rep(1, n_lags)))
  }
  tr_eff <- if (is.finite(tr) && tr > 0.05) tr else 1.0
  w_eff <- if (!is.null(readout_width) && is.finite(readout_width[1L]) && readout_width[1L] > 0) {
    readout_width[1L]
  } else {
    1.5 * tr_eff
  }
  hrf_k <- hrf_library(tr = tr_eff, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = 6.0)[, 1L, drop = FALSE]
  span_lags <- max(lags) - min(lags)
  n_ref <- max(96L, as.integer(ceiling((span_lags + 80.0) / tr_eff)))
  rows_ref <- seq_len(n_ref)
  z_ref <- cbind(1, seq(-1, 1, length.out = n_ref))
  settings_ref <- list(
    lambda = lambda,
    scaling = "noise_sd",
    penalty = penalty,
    backend = "cpp",
    tail_padding = 0L,
    ridge_floor = 1e-8,
    nuisance_mode = "joint",
    observation_mask = "all"
  )
  op_ref <- hrf_inverse_operator(n_ref, hrf_k[, 1L], rows_ref, rho, z_ref, settings_ref)
  t_start <- max(12.0, abs(min(lags)) + 6.0)
  t_end <- max(t_start + 32.0, (n_ref - 16L) * tr_eff - max(0, max(lags)))
  onsets_ref <- seq(t_start, t_end, by = max(6.3, 4.3 * tr_eff))
  if (length(onsets_ref) < 4L) {
    onsets_ref <- seq(t_start, t_start + 30.0 * tr_eff, length.out = 6L)
  }
  m_ref <- length(onsets_ref)
  max_run_time <- (n_ref - 1L) * tr_eff

  S_lags_ref <- lapply(lags, function(lv) {
    build_epoch_readout_matrix(
      onsets_ref, lv, n_ref, tr_eff, 0,
      readout_kernel = readout_kernel, readout_width = w_eff
    )
  })
  ok_list <- lapply(S_lags_ref, function(S_l) rowSums(!is.finite(S_l)) == 0L)
  ok_common <- Reduce("&", ok_list)
  if (sum(ok_common) < 2L) {
    dist_mat <- abs(outer(lags, lags, "-"))
    C_fb <- exp(-0.5 * (dist_mat / (1.35 * tr_eff))^2)
    R_fb <- exp(-0.5 * (dist_mat / (2.10 * tr_eff))^2)
    dimnames(C_fb) <- dimnames(R_fb) <- list(lag_names, lag_names)
    return(list(C_sym = C_fb, R_op = R_fb, d_lag = rep(0.28, n_lags)))
  }
  m_c <- sum(ok_common)
  DG_lags_ref <- lapply(lags, function(lv) {
    sh_onsets <- pmin(max_run_time, pmax(0, onsets_ref + lv))
    G_u <- hrf_event_design_cpp(n_ref, tr_eff, 0, sh_onsets, rep(0, m_ref), seq_len(m_ref), m_ref, hrf_k)[[1L]]
    G_c <- sweep(G_u, 2L, colMeans(G_u), "-")
    op_ref %*% G_c
  })
  K_ref <- matrix(0, n_lags, n_lags)
  U_lag_mat <- matrix(0, n_ref * m_c, n_lags)
  for (l1 in seq_len(n_lags)) {
    S_l1 <- S_lags_ref[[l1]][ok_common, , drop = FALSE]
    for (l2 in seq_len(n_lags)) {
      DG_l2 <- DG_lags_ref[[l2]][, ok_common, drop = FALSE]
      k_diag <- rowSums(S_l1 * t(DG_l2))
      K_ref[l1, l2] <- mean(k_diag[is.finite(k_diag)])
    }
    B_l1 <- S_l1 %*% op_ref
    U_l1 <- compute_ar1_operator_factor(B_l1, rows_ref, rho, center_obs = TRUE)
    U_cent <- sweep(U_l1, 1L, rowMeans(U_l1), "-")
    U_lag_mat[, l1] <- as.numeric(U_cent)
  }
  d_ref <- pmax(diag(K_ref), gain_floor)
  D_inv <- diag(1 / sqrt(d_ref), n_lags)
  C_lag <- D_inv %*% K_ref %*% D_inv
  C_sym <- 0.5 * (C_lag + t(C_lag))
  V_lag <- crossprod(U_lag_mat) / max(1L, m_c - 1L)
  d_v <- sqrt(pmax(diag(V_lag), .Machine$double.eps))
  R_op <- V_lag / outer(d_v, d_v)
  R_op[!is.finite(R_op)] <- 0
  diag(R_op) <- 1
  dimnames(C_sym) <- dimnames(R_op) <- list(lag_names, lag_names)
  list(C_sym = C_sym, R_op = R_op, d_lag = d_ref)
}

# Build lobe-normalized spectral Lavrentiev-Tikhonov cross-lag de-blurring operator W_deb
medusa_build_deblur_operator <- function(
  C_sym,
  lags,
  tr = 1.0,
  alpha1 = 0.025,
  alpha0_sq = 0.0012,
  sigma_lobe = 1.05
) {
  n_lags <- length(lags)
  if (n_lags < 3L || is.null(C_sym) || !all(dim(C_sym) == c(n_lags, n_lags))) {
    return(diag(1, n_lags))
  }
  tr_eff <- if (is.finite(tr) && tr > 0.05) tr else 1.0
  dt_lag <- if (n_lags > 1L) stats::median(diff(sort(unique(lags)))) else tr_eff
  density_scale <- max(1.0, tr_eff / max(dt_lag, 0.1 * tr_eff))
  C_sym <- 0.5 * (C_sym + t(C_sym))
  eig_C <- eigen(C_sym, symmetric = TRUE)
  mu_sc <- pmax(eig_C$values, 0) / density_scale
  f_spec <- (mu_sc / (mu_sc^2 + alpha1 * mu_sc + alpha0_sq)) / density_scale
  W_raw <- eig_C$vectors %*% (f_spec * t(eig_C$vectors))
  R_raw <- W_raw %*% C_sym
  lag_dist <- abs(outer(lags, lags, "-"))
  L_lobe <- exp(-0.5 * (lag_dist / max(sigma_lobe, 0.25))^2)
  g_lobe <- pmax(rowSums(R_raw * L_lobe), 0.20)
  sweep(W_raw, 1L, g_lobe, "/")
}

# Joint cross-lag spectral de-blurring + correlated-SURE 2nd-order boosted P-spline smoother
medusa_smooth_deblur_trajectory <- function(
  beta_vec,
  Sigma_beta,
  lags,
  C_sym = NULL,
  tr = 1.0,
  demix_within_trial = FALSE,
  within_trial_ridge = 0.15,
  recalibrate_wt_lobe = FALSE,
  trajectory_mode = c("auto", "deblur_pspline", "pspline", "deblur", "none"),
  alpha1 = 0.025,
  alpha0_sq = 0.0012,
  sigma_lobe = 1.05
) {
  trajectory_mode <- match.arg(trajectory_mode)
  L <- length(beta_vec)
  if (L < 3L || identical(trajectory_mode, "none")) {
    se_out <- sqrt(pmax(diag(Sigma_beta), .Machine$double.eps))
    return(list(
      beta = as.numeric(beta_vec),
      se = se_out,
      Sigma = Sigma_beta,
      beta_deblur = as.numeric(beta_vec),
      Sigma_deblur = Sigma_beta,
      W_pre = diag(1, L),
      S_smooth = diag(1, L),
      M_op = diag(1, L),
      eta_opt = 0,
      df_eff = L
    ))
  }
  do_deblur <- if (identical(trajectory_mode, "auto")) {
    !isTRUE(demix_within_trial) && !is.null(C_sym)
  } else {
    trajectory_mode %in% c("deblur_pspline", "deblur")
  }
  do_pspline <- trajectory_mode %in% c("auto", "deblur_pspline", "pspline")

  if (do_deblur && !is.null(C_sym)) {
    W_pre <- medusa_build_deblur_operator(
      C_sym = C_sym,
      lags = lags,
      tr = tr,
      alpha1 = alpha1,
      alpha0_sq = alpha0_sq,
      sigma_lobe = sigma_lobe
    )
  } else if (isTRUE(demix_within_trial) && isTRUE(recalibrate_wt_lobe) && !is.null(C_sym)) {
    W_1 <- solve(crossprod(C_sym) + within_trial_ridge * diag(L), t(C_sym))
    W_sym <- W_1 + 0.6 * W_1 %*% (diag(L) - C_sym %*% W_1)
    R_mat <- W_sym %*% C_sym
    diag_base <- pmax(diag(R_mat), 0.15)
    L_lobe <- exp(-0.5 * (outer(lags, lags, "-") / max(sigma_lobe, 0.25))^2)
    g_lobe_wt <- pmax(rowSums(R_mat * L_lobe), 0.20)
    W_pre <- diag(diag_base / g_lobe_wt, L)
  } else {
    W_pre <- diag(1, L)
  }

  b_pre <- as.numeric(W_pre %*% beta_vec)
  Sigma_pre <- W_pre %*% Sigma_beta %*% t(W_pre)
  Sigma_pre <- 0.5 * (Sigma_pre + t(Sigma_pre))

  if (!do_pspline || L < 4L) {
    se_out <- sqrt(pmax(diag(Sigma_pre), .Machine$double.eps))
    return(list(
      beta = b_pre,
      se = se_out,
      Sigma = Sigma_pre,
      beta_deblur = b_pre,
      Sigma_deblur = Sigma_pre,
      W_pre = W_pre,
      S_smooth = diag(1, L),
      M_op = W_pre,
      eta_opt = 0,
      df_eff = sum(diag(W_pre))
    ))
  }

  # Second-order difference penalty on potentially unequally spaced lags
  D2 <- matrix(0, L - 2L, L)
  for (i in seq_len(L - 2L)) {
    dt1 <- max(abs(lags[i + 1L] - lags[i]), 1e-4)
    dt2 <- max(abs(lags[i + 2L] - lags[i + 1L]), 1e-4)
    dt_bar <- 0.5 * (dt1 + dt2)
    D2[i, i:(i + 2L)] <- c(1 / (dt1 * dt_bar), -(1 / dt1 + 1 / dt2) / dt_bar, 1 / (dt2 * dt_bar))
  }
  P2 <- crossprod(D2)

  s2_diag <- pmax(diag(Sigma_pre), 1e-8)
  w_prec <- 1 / s2_diag
  w_prec <- w_prec / mean(w_prec)
  W_p <- diag(w_prec, L)

  eta_max <- if (do_deblur) 0.08 else 1.25
  eta_min <- if (do_deblur) 1e-4 else 2e-3
  eta_grid <- c(0, exp(seq(log(eta_min), log(eta_max), length.out = 28L)))
  IL <- diag(1, L)

  best_sure <- Inf
  best_S <- IL
  best_eta <- 0
  for (eta in eta_grid) {
    if (eta <= 0) {
      S_eta <- IL
    } else {
      S0 <- tryCatch(solve(W_p + eta * P2, W_p), error = function(e) IL)
      # Twicing / 2nd-order boosted P-spline reduces quadratic peak curvature bias from O(eta) to O(eta^2)
      S_eta <- S0 %*% (2 * IL - S0)
    }
    ImS <- IL - S_eta
    diff_b <- as.numeric(ImS %*% b_pre)
    rss_w <- sum(w_prec * (diff_b^2))
    # Exact generalized Stein's Unbiased Risk Estimate (SURE) under correlated cross-lag covariance Sigma_pre
    tr_1 <- sum(rowSums((W_p %*% ImS) * (ImS %*% Sigma_pre)))
    tr_2 <- sum(rowSums((W_p %*% S_eta) * (S_eta %*% Sigma_pre)))
    sure_val <- rss_w - tr_1 + tr_2
    if (is.finite(sure_val) && sure_val < best_sure) {
      best_sure <- sure_val
      best_S <- S_eta
      best_eta <- eta
    }
  }

  M_total <- best_S %*% W_pre
  b_post <- as.numeric(best_S %*% b_pre)
  Sigma_post <- best_S %*% Sigma_pre %*% t(best_S)
  Sigma_post <- 0.5 * (Sigma_post + t(Sigma_post))
  se_post <- sqrt(pmax(diag(Sigma_post), .Machine$double.eps))
  list(
    beta = b_post,
    se = se_post,
    Sigma = Sigma_post,
    beta_deblur = b_pre,
    Sigma_deblur = Sigma_pre,
    W_pre = W_pre,
    S_smooth = best_S,
    M_op = M_total,
    eta_opt = best_eta,
    df_eff = sum(diag(best_S))
  )
}

# Vectorized 1D Threshold-Free Cluster Enhancement (Smith & Nichols, 2009) across a matrix
# of lag t-statistics T_mat (L x B). Uses a fixed integration step dh so observed and
# permutation/Monte-Carlo null columns are evaluated on the exact same threshold grid.
medusa_tfce_1d_matrix <- function(T_mat, dt = 1.0, E = 0.5, H = 2.0, dh = 0.10, two_sided = TRUE) {
  if (is.null(dim(T_mat))) T_mat <- matrix(T_mat, ncol = 1L)
  L <- nrow(T_mat)
  B <- ncol(T_mat)
  tfce_side <- function(M_pos) {
    max_val <- max(M_pos, na.rm = TRUE)
    if (!is.finite(max_val) || max_val < dh) return(matrix(0, L, B))
    h_max <- min(max_val, 14.0)
    h_grid <- seq(dh, h_max, by = dh)
    acc <- matrix(0, L, B)
    for (h in h_grid) {
      A <- M_pos >= h
      if (!any(A)) break
      # Forward accumulation of contiguous supra-threshold run lengths along lags
      F_len <- matrix(0L, L, B)
      F_len[1L, ] <- as.integer(A[1L, ])
      if (L >= 2L) {
        for (l in 2L:L) {
          al <- as.integer(A[l, ])
          F_len[l, ] <- al * (F_len[l - 1L, ] + 1L)
        }
        # Backward propagation so every lag in a contiguous cluster receives the full cluster extent
        for (l in (L - 1L):1L) {
          al <- as.integer(A[l, ])
          F_len[l, ] <- al * pmax(F_len[l, ], F_len[l + 1L, ])
        }
      }
      acc <- acc + ((F_len * dt)^E) * (h^H) * dh
    }
    acc
  }
  if (isTRUE(two_sided)) {
    tfce_pos <- tfce_side(pmax(T_mat, 0))
    tfce_neg <- tfce_side(pmax(-T_mat, 0))
    pmax(tfce_pos, tfce_neg)
  } else {
    tfce_side(pmax(T_mat, 0))
  }
}

# Continuous sub-grid natural cubic spline latency, center-of-mass, and FWHM with
# parametric bootstrap SEs and confidence intervals from the joint L x L covariance
medusa_spline_latency <- function(
  beta_vec,
  Sigma_beta,
  lags,
  fine_dt = 0.01,
  n_boot = 400L,
  conf_level = 0.95,
  polarity = c("positive", "negative", "abs")
) {
  polarity <- match.arg(polarity)
  L <- length(lags)
  if (L < 2L || any(!is.finite(beta_vec))) {
    pk <- if (L >= 1L && is.finite(beta_vec[1L])) lags[1L] else NA_real_
    return(list(
      peak_lag = pk,
      peak_lag_grid = pk,
      peak_lag_se = NA_real_,
      peak_lag_conf.low = NA_real_,
      peak_lag_conf.high = NA_real_,
      peak_val = if (L >= 1L) beta_vec[1L] else NA_real_,
      peak_val_se = NA_real_,
      com_lag = pk,
      com_lag_se = NA_real_,
      fwhm = NA_real_,
      fwhm_se = NA_real_
    ))
  }
  tau_fine <- seq(min(lags), max(lags), by = fine_dt)
  G <- length(tau_fine)
  # Construct linear natural cubic spline interpolation basis M_fine (G x L) once
  IL <- diag(1, L)
  M_fine <- vapply(seq_len(L), function(j) {
    stats::spline(x = lags, y = IL[, j], xout = tau_fine, method = "natural")$y
  }, numeric(G))
  if (is.null(dim(M_fine))) M_fine <- matrix(M_fine, nrow = G, ncol = L)

  extract_metrics <- function(f_curve) {
    s_curve <- if (identical(polarity, "positive")) {
      f_curve
    } else if (identical(polarity, "negative")) {
      -f_curve
    } else {
      abs(f_curve)
    }
    idx_pk <- which.max(s_curve)
    t_pk <- tau_fine[idx_pk]
    v_pk <- f_curve[idx_pk]
    s_pk <- s_curve[idx_pk]
    pos_w <- pmax(s_curve, 0)^2
    t_com <- if (sum(pos_w) > 1e-12) sum(tau_fine * pos_w) / sum(pos_w) else t_pk
    if (is.finite(s_pk) && s_pk > 1e-8) {
      above_half <- s_curve >= (0.5 * s_pk)
      # Contiguous lobe containing idx_pk
      left_idx <- idx_pk
      while (left_idx > 1L && above_half[left_idx - 1L]) left_idx <- left_idx - 1L
      right_idx <- idx_pk
      while (right_idx < G && above_half[right_idx + 1L]) right_idx <- right_idx + 1L
      fwhm_val <- tau_fine[right_idx] - tau_fine[left_idx]
    } else {
      fwhm_val <- NA_real_
    }
    c(peak_lag = t_pk, peak_val = v_pk, com_lag = t_com, fwhm = fwhm_val)
  }

  f_hat <- as.numeric(M_fine %*% beta_vec)
  obs_m <- extract_metrics(f_hat)
  grid_s <- if (identical(polarity, "positive")) {
    beta_vec
  } else if (identical(polarity, "negative")) {
    -beta_vec
  } else {
    abs(beta_vec)
  }
  peak_lag_grid <- lags[which.max(grid_s)]

  if (n_boot <= 0L || any(!is.finite(Sigma_beta))) {
    return(list(
      peak_lag = unname(obs_m["peak_lag"]),
      peak_lag_grid = peak_lag_grid,
      peak_lag_se = NA_real_,
      peak_lag_conf.low = NA_real_,
      peak_lag_conf.high = NA_real_,
      peak_val = unname(obs_m["peak_val"]),
      peak_val_se = NA_real_,
      com_lag = unname(obs_m["com_lag"]),
      com_lag_se = NA_real_,
      fwhm = unname(obs_m["fwhm"]),
      fwhm_se = NA_real_
    ))
  }

  Sigma_sym <- 0.5 * (Sigma_beta + t(Sigma_beta))
  eig_S <- eigen(Sigma_sym, symmetric = TRUE)
  L_half <- eig_S$vectors %*% (sqrt(pmax(eig_S$values, 0)) * t(eig_S$vectors))
  Z_draw <- matrix(stats::rnorm(L * n_boot), nrow = L, ncol = n_boot)
  B_draw <- beta_vec + L_half %*% Z_draw
  F_draw <- M_fine %*% B_draw
  boot_mat <- apply(F_draw, 2L, extract_metrics)

  alpha_tail <- (1.0 - conf_level) / 2.0
  pk_draws <- boot_mat["peak_lag", ]
  val_draws <- boot_mat["peak_val", ]
  com_draws <- boot_mat["com_lag", ]
  fwhm_draws <- boot_mat["fwhm", ]
  ci_pk <- stats::quantile(pk_draws, probs = c(alpha_tail, 1.0 - alpha_tail), na.rm = TRUE, names = FALSE)

  list(
    peak_lag = unname(obs_m["peak_lag"]),
    peak_lag_grid = peak_lag_grid,
    peak_lag_se = stats::sd(pk_draws, na.rm = TRUE),
    peak_lag_conf.low = ci_pk[1L],
    peak_lag_conf.high = ci_pk[2L],
    peak_val = unname(obs_m["peak_val"]),
    peak_val_se = stats::sd(val_draws, na.rm = TRUE),
    com_lag = unname(obs_m["com_lag"]),
    com_lag_se = stats::sd(com_draws, na.rm = TRUE),
    fwhm = unname(obs_m["fwhm"]),
    fwhm_se = stats::sd(fwhm_draws, na.rm = TRUE)
  )
}

#' Joint Multi-Lag Trajectory Modeling and Operator-Informed Inference for MEDuSA
#'
#' Replaces independent lag-by-lag post-hoc testing in MEDuSA Stage 4 with a
#' unified multi-lag estimation and inference framework that exploits the exact
#' Stage-2/3 cross-lag point-spread resolution matrix \eqn{C_{\mathrm{sym}}} and
#' cross-lag operator error correlation \eqn{R_{\mathrm{op}}}:
#' \enumerate{
#'   \item \strong{Lobe-Normalized Spectral Lavrentiev--Tikhonov Cross-Lag De-Blurring
#'     + Correlated-SURE 2nd-Order Boosted P-Spline Smoothing}: sharpens multiphasic
#'     and ramped within-trial neural waveforms without boundary distortion or
#'     interior gain inflation.
#'   \item \strong{Multivariate Empirical Bayes (EB) Subject Trajectory Shrinkage}:
#'     borrows strength across all \eqn{L} lags simultaneously via the between-subject
#'     trajectory covariance \eqn{\hat{\Sigma}_b} and subject-specific cross-lag
#'     error covariance \eqn{\hat{\Sigma}_{e,s}}.
#'   \item \strong{Operator-Informed Simultaneous Cross-Lag Inference}: calibrates
#'     epoch-wide family-wise error rate (FWER) via exact correlated multivariate-\eqn{t}
#'     Max-\eqn{|T|} simultaneous confidence bands (\code{p_maxT}, \code{conf.low.sim},
#'     \code{conf.high.sim}) and 1D Threshold-Free Cluster Enhancement (\code{tfce_stat},
#'     \code{p_tfce}), avoiding the severe conservatism of Bonferroni/BY-FDR under
#'     strong positive cross-lag correlation.
#'   \item \strong{Continuous Sub-Grid Spline Latency & Shape Estimation}: estimates
#'     continuous peak latency (\code{peak_lag}), center-of-mass latency (\code{com_lag}),
#'     and full-width at half-maximum (\code{fwhm}) with parametric-bootstrap SEs and
#'     confidence intervals.
#' }
#'
#' @param data Either: (a) a \code{\link{decon_epoch_regression}} result list (or a
#'   list of per-subject \code{decon_epoch_regression} results), (b) a trial-level
#'   MEDuSA panel \code{data.frame}/\code{data.table} (single-subject or multi-subject),
#'   or (c) a lag-level coefficient table with columns \code{lag} (or \code{evt_time}),
#'   \code{estimate} (or \code{beta_corrected}), and \code{std.error} (or \code{se_corrected}).
#' @param covariates Optional character vector of covariate names when \code{data} is a
#'   trial-level panel. If \code{NULL}, inferred automatically from panel columns.
#' @param lag_resolution Optional named list (by region/parcel/voxel) of \eqn{L \times L}
#'   cross-lag resolution matrices \eqn{C_{\mathrm{sym}}}, or a single \eqn{L \times L} matrix.
#' @param lag_error_cor Optional named list (by region/parcel/voxel) of \eqn{L \times L}
#'   operator cross-lag error correlation matrices \eqn{R_{\mathrm{op}}}, or a single matrix.
#' @param trajectory_mode Trajectory estimation mode: \code{"auto"} (default; applies
#'   spectral cross-lag de-blurring when \code{!demix_within_trial} followed by
#'   correlated-SURE boosted P-spline smoothing), \code{"deblur_pspline"},
#'   \code{"pspline"}, \code{"deblur"}, or \code{"none"}.
#' @param use_whitened Logical; if \code{TRUE} (default) and whitened columns
#'   (\code{decon_demixed_whitened}, \code{*_whitened}) are present in a trial panel,
#'   uses them for within-subject lag regressions.
#' @param tr Acquisition TR in seconds (default \code{1.0}), used when constructing
#'   reference cross-lag operators if \code{lag_resolution} is not attached.
#' @param readout_kernel Readout kernel (\code{"triangle"}, \code{"boxcar"},
#'   \code{"gaussian"}, \code{"linear"}) for reference operator construction.
#' @param readout_width Readout kernel width in seconds (default \code{1.5}).
#' @param decon_lambda Reference deconvolution penalty \eqn{\lambda} (default \code{0.8}).
#' @param rho_ar1 Reference AR(1) parameter (default \code{0.30}).
#' @param demix_within_trial Logical indicating whether Stage-3 within-trial cross-lag
#'   de-mixing was already applied to the input panel/coefficients (default \code{FALSE}).
#' @param recalibrate_wt_lobe Logical; if \code{TRUE} and \code{demix_within_trial = TRUE},
#'   recalibrates interior vs. boundary lobe gain on \code{W_lag}.
#' @param contrasts Optional named list of lag contrast specifications. Each entry can be
#'   either a numeric weight vector of length \eqn{L} or a list \code{list(early = c(t1, t2),
#'   late = c(t3, t4))} specifying two lag windows to contrast (\code{late - early}).
#' @param contrast_source Trajectory source for linear window contrasts:
#'   \code{"unsmoothed"} (default; evaluates contrasts on the cross-lag de-blurred
#'   trajectory prior to P-spline smoothing so contrast point estimates and CIs are
#'   free of P-spline shrinkage bias), \code{"raw"} (raw lag-by-lag estimates), or
#'   \code{"multilag"} (P-spline smoothed trajectory).
#' @param n_perm Number of correlated Monte Carlo studentized draws for
#'   simultaneous Max-\eqn{|T|}, 1D TFCE, and omnibus calibration (default \code{4000L}).
#' @param n_boot_latency Number of parametric bootstrap draws for continuous spline
#'   latency SEs and confidence intervals (default \code{400L}).
#' @param fine_dt Sub-grid spline evaluation step in seconds (default \code{0.01}).
#' @param conf_level Simultaneous and pointwise confidence level (default \code{0.95}).
#' @param polarity Peak polarity for latency extraction (\code{"positive"},
#'   \code{"negative"}, or \code{"abs"}).
#' @param seed Optional random seed for reproducible Monte Carlo / bootstrap draws.
#' @return An object of class \code{"medusa_multilag"} containing:
#'   \describe{
#'     \item{\code{coefficients}}{\code{data.table} of group/unit lag estimates
#'       (\code{estimate}, \code{std.error}, \code{estimate_multilag},
#'       \code{std.error_multilag}), pointwise \eqn{t}/\eqn{p}-values, FDR-adjusted
#'       \eqn{p}-values (\code{p_fdr_bh}, \code{p_fdr_by}), simultaneous Max-\eqn{|T|}
#'       \eqn{p}-values and confidence bands (\code{p_maxT}, \code{conf.low.sim},
#'       \code{conf.high.sim}), 1D TFCE statistics and family-wise \eqn{p}-values
#'       (\code{tfce_stat}, \code{p_tfce}), and epoch omnibus tests
#'       (\code{omnibus_stat}, \code{omnibus_p}).}
#'     \item{\code{latency}}{\code{data.table} of continuous sub-grid spline peak
#'       latency (\code{peak_lag}, \code{peak_lag_se}, \code{peak_lag_conf.low},
#'       \code{peak_lag_conf.high}), discrete grid peak (\code{peak_lag_grid}),
#'       center-of-mass latency (\code{com_lag}, \code{com_lag_se}), and FWHM
#'       (\code{fwhm}, \code{fwhm_se}) per unit and term.}
#'     \item{\code{contrasts}}{\code{data.table} of covariance-propagated multi-lag
#'       window contrasts (\code{estimate}, \code{std.error}, \code{statistic},
#'       \code{df}, \code{p.value}, \code{conf.low}, \code{conf.high}).}
#'     \item{\code{subjects}}{\code{data.table} of subject-specific lag trajectories
#'       (\code{estimate_raw}, \code{estimate_smooth}, and Multivariate Empirical
#'       Bayes \code{estimate_eb}) when multi-subject data are supplied.}
#'     \item{\code{subject_latency}}{\code{data.table} of subject-specific continuous
#'       spline peak and center-of-mass latencies extracted from \code{estimate_eb}.}
#'     \item{\code{covariance}}{Nested list of \eqn{L \times L} raw and multi-lag
#'       trajectory covariance and correlation matrices per unit and term.}
#'   }
#' @export
medusa_multilag_infer <- function(
  data,
  covariates = NULL,
  lag_resolution = NULL,
  lag_error_cor = NULL,
  trajectory_mode = c("auto", "deblur_pspline", "pspline", "deblur", "none"),
  use_whitened = TRUE,
  tr = 1.0,
  readout_kernel = "triangle",
  readout_width = 1.5,
  decon_lambda = 0.8,
  rho_ar1 = 0.30,
  demix_within_trial = FALSE,
  recalibrate_wt_lobe = FALSE,
  contrasts = NULL,
  contrast_source = c("unsmoothed", "raw", "multilag"),
  n_perm = 4000L,
  n_boot_latency = 400L,
  fine_dt = 0.01,
  conf_level = 0.95,
  polarity = c("positive", "negative", "abs"),
  seed = NULL
) {
  trajectory_mode <- match.arg(trajectory_mode)
  contrast_source <- match.arg(contrast_source)
  polarity <- match.arg(polarity)
  if (!is.null(seed)) set.seed(as.integer(seed))

  subj_lag_ops <- list()
  subj_trial_ops <- list()
  global_trial_ops <- NULL

  avg_named_mat_list <- function(mat_lists) {
    mat_lists <- Filter(Negate(is.null), mat_lists)
    if (!length(mat_lists)) return(NULL)
    if (all(vapply(mat_lists, is.matrix, logical(1L)))) {
      d0 <- dim(mat_lists[[1L]])
      ok_m <- vapply(mat_lists, function(m) all(dim(m) == d0) && all(is.finite(m)), logical(1L))
      if (!any(ok_m)) return(mat_lists[[1L]])
      return(Reduce("+", mat_lists[ok_m]) / sum(ok_m))
    }
    all_keys <- unique(unlist(lapply(mat_lists, names)))
    if (!length(all_keys)) return(mat_lists[[1L]])
    out <- list()
    for (k in all_keys) {
      k_mats <- lapply(mat_lists, function(x) if (is.list(x)) x[[k]] else NULL)
      k_mats <- Filter(is.matrix, k_mats)
      if (!length(k_mats)) next
      d0 <- dim(k_mats[[1L]])
      ok_m <- vapply(k_mats, function(m) all(dim(m) == d0) && all(is.finite(m)), logical(1L))
      out[[k]] <- if (any(ok_m)) Reduce("+", k_mats[ok_m]) / sum(ok_m) else k_mats[[1L]]
    }
    out
  }

  # 1. Unpack input object (decon_epoch_regression result, mixed_by result, list of results/panels, panel, or coef table)
  if (is.list(data) && !is.data.frame(data)) {
    if (!is.null(data$coefficients) && is.data.frame(data$coefficients)) {
      # Single decon_epoch_regression output
      if (is.null(lag_resolution)) lag_resolution <- data$lag_resolution
      if (is.null(lag_error_cor)) lag_error_cor <- data$lag_error_cor
      global_trial_ops <- if (!is.null(data$run_trial_operators)) {
        data$run_trial_operators
      } else if (!is.null(data$panel)) {
        attr(data$panel, "run_trial_operators")
      } else NULL
      if (!is.null(data$settings$demix_within_trial)) {
        demix_within_trial <- isTRUE(data$settings$demix_within_trial)
      }
      if (!is.null(data$settings$readout_kernel)) readout_kernel <- data$settings$readout_kernel
      if (!is.null(data$settings$readout_width)) readout_width <- data$settings$readout_width
      if (!is.null(data$panel) && is.data.frame(data$panel) && nrow(data$panel) > 0L) {
        dt_in <- data.table::as.data.table(data$panel)
        input_kind <- "panel"
      } else {
        dt_in <- data.table::as.data.table(data$coefficients)
        input_kind <- "coef"
      }
    } else if ((!is.null(data$coef_df_reml) && is.data.frame(data$coef_df_reml)) ||
               (!is.null(data$coef_df_ml) && is.data.frame(data$coef_df_ml))) {
      # mixed_by() output object
      cf_mb <- if (!is.null(data$coef_df_reml) && is.data.frame(data$coef_df_reml)) data$coef_df_reml else data$coef_df_ml
      dt_in <- data.table::as.data.table(cf_mb)
      input_kind <- "coef"
    } else if (length(data) >= 1L && (is.data.frame(data[[1L]]) ||
               (is.list(data[[1L]]) && (!is.null(data[[1L]]$panel) || !is.null(data[[1L]]$coefficients))))) {
      # List of per-subject decon_epoch_regression outputs or per-subject panel data.tables
      sub_names <- names(data)
      if (is.null(sub_names)) sub_names <- sprintf("sub_%03d", seq_along(data))
      is_df_list <- all(vapply(data, is.data.frame, logical(1L)))
      if (is_df_list) {
        c_res_list <- lapply(data, attr, "lag_resolution")
        r_cor_list <- lapply(data, attr, "lag_error_cor")
        if (is.null(lag_resolution)) lag_resolution <- avg_named_mat_list(c_res_list)
        if (is.null(lag_error_cor)) lag_error_cor <- avg_named_mat_list(r_cor_list)
        dt_parts <- vector("list", length(data))
        for (i in seq_along(data)) {
          p_i <- data.table::copy(data.table::as.data.table(data[[i]]))
          sc_i <- intersect(c("id", "subject", "sub"), names(p_i))[1L]
          sid_i <- if (!is.na(sc_i)) as.character(p_i[[sc_i]][1L]) else sub_names[i]
          if (is.na(sc_i)) p_i[, id := sid_i]
          subj_lag_ops[[sid_i]] <- list(C = attr(data[[i]], "lag_resolution"), R = attr(data[[i]], "lag_error_cor"))
          subj_trial_ops[[sid_i]] <- attr(data[[i]], "run_trial_operators")
          dt_parts[[i]] <- p_i
        }
        dt_in <- data.table::rbindlist(dt_parts, fill = TRUE)
        input_kind <- if (any(c("decon_demixed", "decon_interp", "decon_raw", "decon_demixed_whitened") %in% names(dt_in))) "panel" else "coef"
      } else {
        c_res_list <- lapply(data, function(x) if (!is.null(x$lag_resolution)) x$lag_resolution else if (!is.null(x$panel)) attr(x$panel, "lag_resolution") else NULL)
        r_cor_list <- lapply(data, function(x) if (!is.null(x$lag_error_cor)) x$lag_error_cor else if (!is.null(x$panel)) attr(x$panel, "lag_error_cor") else NULL)
        if (is.null(lag_resolution)) lag_resolution <- avg_named_mat_list(c_res_list)
        if (is.null(lag_error_cor)) lag_error_cor <- avg_named_mat_list(r_cor_list)
        if (!is.null(data[[1L]]$settings$demix_within_trial)) {
          demix_within_trial <- isTRUE(data[[1L]]$settings$demix_within_trial)
        }
        if (!is.null(data[[1L]]$settings$readout_kernel)) readout_kernel <- data[[1L]]$settings$readout_kernel
        if (!is.null(data[[1L]]$settings$readout_width)) readout_width <- data[[1L]]$settings$readout_width
        has_panels <- all(vapply(data, function(x) !is.null(x$panel) && is.data.frame(x$panel), logical(1L)))
        if (has_panels) {
          dt_parts <- vector("list", length(data))
          for (i in seq_along(data)) {
            p_i <- data.table::copy(data.table::as.data.table(data[[i]]$panel))
            sc_i <- intersect(c("id", "subject", "sub"), names(p_i))[1L]
            sid_i <- if (!is.na(sc_i)) as.character(p_i[[sc_i]][1L]) else sub_names[i]
            if (is.na(sc_i)) p_i[, id := sid_i]
            subj_lag_ops[[sid_i]] <- list(C = c_res_list[[i]], R = r_cor_list[[i]])
            subj_trial_ops[[sid_i]] <- if (!is.null(data[[i]]$run_trial_operators)) {
              data[[i]]$run_trial_operators
            } else {
              attr(data[[i]]$panel, "run_trial_operators")
            }
            dt_parts[[i]] <- p_i
          }
          dt_in <- data.table::rbindlist(dt_parts, fill = TRUE)
          input_kind <- "panel"
        } else {
          dt_parts <- vector("list", length(data))
          for (i in seq_along(data)) {
            c_i <- data.table::copy(data.table::as.data.table(data[[i]]$coefficients))
            sc_i <- intersect(c("id", "subject", "sub"), names(c_i))[1L]
            sid_i <- if (!is.na(sc_i)) as.character(c_i[[sc_i]][1L]) else sub_names[i]
            if (is.na(sc_i)) c_i[, id := sid_i]
            subj_lag_ops[[sid_i]] <- list(C = c_res_list[[i]], R = r_cor_list[[i]])
            dt_parts[[i]] <- c_i
          }
          dt_in <- data.table::rbindlist(dt_parts, fill = TRUE)
          input_kind <- "coef"
        }
      }
    } else {
      stop("Unsupported list input to `medusa_multilag_infer()`.")
    }
  } else if (is.data.frame(data)) {
    if (is.null(lag_resolution)) lag_resolution <- attr(data, "lag_resolution")
    if (is.null(lag_error_cor)) lag_error_cor <- attr(data, "lag_error_cor")
    global_trial_ops <- attr(data, "run_trial_operators")
    dt_in <- data.table::as.data.table(data)
    if (any(c("decon_demixed", "decon_interp", "decon_raw", "decon_demixed_whitened") %in% names(dt_in))) {
      input_kind <- "panel"
    } else if (any(c("beta_corrected", "beta_gls", "estimate", "beta_weighted", "beta_unweighted", "beta_raw") %in% names(dt_in))) {
      input_kind <- "coef"
    } else {
      stop("Data frame supplied to `medusa_multilag_infer()` must be a MEDuSA trial panel or coefficient table.")
    }
  } else {
    stop("`data` must be a `decon_epoch_regression` result, list of results, `mixed_by` output, or `data.frame`/`data.table`.")
  }

  if (identical(input_kind, "coef") && ("effect" %in% names(dt_in))) {
    dt_in <- dt_in[is.na(effect) | as.character(effect) == "fixed"]
  }

  # Identify unit column (region / parcel / voxel)
  unit_col <- intersect(c("region", "parcel", "voxel", "unit"), names(dt_in))[1L]
  if (is.na(unit_col)) {
    dt_in <- data.table::copy(dt_in)
    dt_in[, unit := "unit_1"]
    unit_col <- "unit"
  }
  lag_col <- intersect(c("evt_time", "lag"), names(dt_in))[1L]
  if (is.na(lag_col)) stop("Could not find `evt_time` or `lag` column in `data`.")
  subj_col <- intersect(c("id", "subject", "sub"), names(dt_in))[1L]

  lags <- sort(unique(as.numeric(dt_in[[lag_col]])))
  L <- length(lags)
  if (L < 2L) stop("`medusa_multilag_infer()` requires at least 2 lags.")
  dt_lag <- stats::median(diff(lags))

  # Construct or retrieve reference operators C_sym and R_op (unit-level and optional subject-specific)
  ref_ops <- NULL
  extract_op_mat <- function(op_obj, u_id) {
    if (is.matrix(op_obj) && all(dim(op_obj) == c(L, L))) return(op_obj)
    if (is.list(op_obj) && !is.null(op_obj[[as.character(u_id)]]) &&
        all(dim(op_obj[[as.character(u_id)]]) == c(L, L))) {
      return(op_obj[[as.character(u_id)]])
    }
    NULL
  }
  get_unit_operators <- function(u_id, sid_opt = NULL) {
    C_u <- NULL
    R_u <- NULL
    if (!is.null(sid_opt) && !is.null(subj_lag_ops[[as.character(sid_opt)]])) {
      C_u <- extract_op_mat(subj_lag_ops[[as.character(sid_opt)]]$C, u_id)
      R_u <- extract_op_mat(subj_lag_ops[[as.character(sid_opt)]]$R, u_id)
    }
    if (is.null(C_u)) C_u <- extract_op_mat(lag_resolution, u_id)
    if (is.null(R_u)) R_u <- extract_op_mat(lag_error_cor, u_id)
    if (is.null(C_u) || is.null(R_u)) {
      if (is.null(ref_ops)) {
        ref_ops <<- medusa_reference_lag_operators(
          lags = lags, tr = tr, readout_kernel = readout_kernel,
          readout_width = readout_width, lambda = decon_lambda, rho = rho_ar1
        )
      }
      if (is.null(C_u)) C_u <- ref_ops$C_sym
      if (is.null(R_u)) R_u <- ref_ops$R_op
    }
    list(C_sym = C_u, R_op = R_u)
  }

  # 2. Identify covariates and columns
  if (identical(input_kind, "panel")) {
    has_w_y <- isTRUE(use_whitened) && ("decon_demixed_whitened" %in% names(dt_in)) &&
      any(is.finite(dt_in$decon_demixed_whitened))
    y_col_unw <- if ("decon_demixed" %in% names(dt_in) && any(is.finite(dt_in$decon_demixed))) {
      "decon_demixed"
    } else if ("decon_interp" %in% names(dt_in) && any(is.finite(dt_in$decon_interp))) {
      "decon_interp"
    } else {
      "decon_raw"
    }
    y_col <- if (has_w_y) "decon_demixed_whitened" else y_col_unw
    w_col <- if ("precision_weight" %in% names(dt_in)) "precision_weight" else NULL
    w_cohort_col <- if ("precision_weight_cohort" %in% names(dt_in)) "precision_weight_cohort" else w_col

    if (is.null(covariates)) {
      reserved_cols <- c(
        "id", "subject", "sub", "run", "trial", "event_index", "onset", "evt_time", "lag",
        "voxel", "parcel", "region", "unit", "decon_raw", "decon_demixed", "decon_eb",
        "decon_interp", "se_raw", "se_demixed", "reliability_trial", "precision_weight",
        "precision_weight_cohort", "precision_weight_op", "precision_weight_raw",
        "error_ac1", "tau2_bio", "shared_noise_loading",
        "decon_demixed_whitened", "se_whitened", "decon_eb_whitened"
      )
      cand_covs <- setdiff(names(dt_in), reserved_cols)
      cand_covs <- cand_covs[!grepl("_whitened$", cand_covs)]
      cand_covs <- cand_covs[vapply(cand_covs, function(cn) is.numeric(dt_in[[cn]]), logical(1L))]
      covariates <- if (length(cand_covs) > 0L) cand_covs else "intercept"
    } else {
      keep_int <- intersect(covariates, "intercept")
      keep_num <- intersect(setdiff(covariates, "intercept"), names(dt_in))
      covariates <- unique(c(keep_int, keep_num))
      if (!length(covariates)) stop("None of the requested `covariates` were found in `data`.")
    }
  } else {
    est_col <- intersect(c("beta_corrected", "beta_gls", "estimate", "beta_weighted", "beta_unweighted", "beta_raw"), names(dt_in))[1L]
    se_col <- intersect(c("se_corrected", "se_gls", "std.error", "se_weighted", "se_unweighted"), names(dt_in))[1L]
    if (is.na(est_col)) stop("Could not find coefficient estimate column in `data`.")
    if (!("term" %in% names(dt_in))) {
      dt_in <- data.table::copy(dt_in)
      dt_in[, term := "x"]
    }
    u_terms <- unique(as.character(dt_in$term))
    u_terms <- u_terms[!is.na(u_terms) & nzchar(u_terms)]
    if (is.null(covariates)) {
      non_int <- setdiff(u_terms, "(Intercept)")
      covariates <- if (length(non_int) > 0L) non_int else u_terms
    } else {
      covariates <- intersect(covariates, u_terms)
    }
  }

  units <- unique(as.character(dt_in[[unit_col]]))

  # Helper to look up run_trial_operators for (sid, u_id, li)
  lag_key_names <- sprintf("lag_%g", lags)
  get_run_trial_ops <- function(sid_chr, u_id_chr, li_idx) {
    ops_root <- if (!is.null(subj_trial_ops[[sid_chr]])) subj_trial_ops[[sid_chr]] else global_trial_ops
    if (is.null(ops_root) || !is.list(ops_root)) return(NULL)
    ops_u <- ops_root[[u_id_chr]]
    if (is.null(ops_u) && length(ops_root) == 1L) ops_u <- ops_root[[1L]]
    if (is.null(ops_u) || !is.list(ops_u)) return(NULL)
    if (!is.null(ops_u[[lag_key_names[li_idx]]])) return(ops_u[[lag_key_names[li_idx]]])
    if (length(ops_u) >= li_idx) return(ops_u[[li_idx]])
    NULL
  }

  coef_out_list <- list()
  lat_out_list <- list()
  con_out_list <- list()
  subj_out_list <- list()
  subj_lat_list <- list()
  cov_out_list <- list()

  for (u_id in units) {
    dt_u <- dt_in[as.character(dt_in[[unit_col]]) == u_id]
    ops_u <- get_unit_operators(u_id)
    C_u <- ops_u$C_sym
    R_op_u <- ops_u$R_op
    CC_u <- tcrossprod(C_u)
    d_cc_u <- sqrt(pmax(diag(CC_u), .Machine$double.eps))
    R_res_u <- CC_u / outer(d_cc_u, d_cc_u)
    R_res_u <- 0.5 * (R_res_u + t(R_res_u))
    diag(R_res_u) <- 1.0
    cov_out_list[[u_id]] <- list()

    # Precompute panel regressions jointly across covariates so correlated within-run predictors
    # are partialed in multiple regression and Level-2 between-subject covariates are estimated across subjects
    panel_fit_cache <- list()
    if (identical(input_kind, "panel")) {
      s_ids_u <- if (!is.na(subj_col)) unique(as.character(dt_u[[subj_col]])) else "sub_1"
      S_u <- length(s_ids_u)
      dt_u_l1 <- dt_u[abs(as.numeric(dt_u[[lag_col]]) - lags[1L]) < 1e-6]
      if (nrow(dt_u_l1) == 0L) dt_u_l1 <- dt_u

      non_int_covs <- setdiff(covariates, "intercept")
      cov_within <- character(0L)
      cov_between_subj <- character(0L)
      cov_between_run <- character(0L)

      for (cv in non_int_covs) {
        xv_all <- as.numeric(dt_u_l1[[cv]])
        tot_sd <- stats::sd(xv_all, na.rm = TRUE)
        if (!is.finite(tot_sd) || tot_sd <= 1e-12) next
        # Compute within-(subject, run) variation
        if ("run" %in% names(dt_u_l1)) {
          sr_grp <- if (!is.na(subj_col)) paste(dt_u_l1[[subj_col]], dt_u_l1$run, sep = "||") else as.character(dt_u_l1$run)
        } else {
          sr_grp <- if (!is.na(subj_col)) as.character(dt_u_l1[[subj_col]]) else rep("g1", nrow(dt_u_l1))
        }
        xv_w <- xv_all - stats::ave(xv_all, sr_grp, FUN = function(z) mean(z, na.rm = TRUE))
        w_sd <- stats::sd(xv_w, na.rm = TRUE)
        if (is.finite(w_sd) && w_sd > 1e-7 * max(1.0, tot_sd)) {
          cov_within <- c(cov_within, cv)
        } else if (!is.na(subj_col) && S_u >= 2L) {
          s_grp <- as.character(dt_u_l1[[subj_col]])
          xv_ws <- xv_all - stats::ave(xv_all, s_grp, FUN = function(z) mean(z, na.rm = TRUE))
          ws_sd <- stats::sd(xv_ws, na.rm = TRUE)
          if (!is.finite(ws_sd) || ws_sd <= 1e-7 * max(1.0, tot_sd)) {
            cov_between_subj <- c(cov_between_subj, cv)
          } else {
            cov_between_run <- c(cov_between_run, cv)
          }
        } else {
          cov_between_run <- c(cov_between_run, cv)
        }
      }

      # (A) Intercept and/or Level-2 between-subject mean extraction per (subject, lag)
      need_subj_means <- ("intercept" %in% covariates) || (length(cov_between_subj) > 0L)
      if (need_subj_means) {
        B_int <- matrix(NA_real_, S_u, L, dimnames = list(s_ids_u, as.character(lags)))
        SE_int <- matrix(NA_real_, S_u, L, dimnames = list(s_ids_u, as.character(lags)))
        df_int <- numeric(S_u)
        sc_int_accum <- matrix(0, L, L)
        sc_int_n <- 0L
        for (si in seq_len(S_u)) {
          sid_i <- s_ids_u[si]
          dt_su <- if (!is.na(subj_col)) dt_u[as.character(dt_u[[subj_col]]) == sid_i] else dt_u
          has_rt <- all(c("run", "trial") %in% names(dt_su))
          rt_keys <- if (has_rt) unique(paste(dt_su$run, dt_su$trial, sep = "_")) else NULL
          score_mat <- if (has_rt) matrix(NA_real_, length(rt_keys), L) else NULL
          df_l_vec <- numeric(L)
          for (li in seq_len(L)) {
            lag_val <- lags[li]
            dt_sul <- dt_su[abs(as.numeric(dt_su[[lag_col]]) - lag_val) < 1e-6]
            if (nrow(dt_sul) < 3L) next
            yv <- as.numeric(dt_sul[[y_col_unw]])
            wv <- if (!is.null(w_cohort_col)) as.numeric(dt_sul[[w_cohort_col]]) else rep(1.0, length(yv))
            wv[!is.finite(wv) | wv <= 0] <- 1.0
            ok_i <- is.finite(yv)
            if (sum(ok_i) < 3L) next
            w_ok <- wv[ok_i] / mean(wv[ok_i])
            y_ok <- yv[ok_i]
            b_hat <- sum(w_ok * y_ok) / sum(w_ok)
            res_i <- y_ok - b_hat
            df_i <- max(1L, sum(ok_i) - 1L)
            # Account for lag-1 within-run autocorrelation of trial residuals when estimating mean SE
            ac1_i <- if ("error_ac1" %in% names(dt_sul)) stats::median(as.numeric(dt_sul$error_ac1[ok_i]), na.rm = TRUE) else 0.0
            if (!is.finite(ac1_i)) ac1_i <- 0.0
            ac_infl <- max(0.5, min(3.0, (1.0 + ac1_i) / max(0.15, 1.0 - ac1_i)))
            se_hat <- sqrt(sum((w_ok^2) * (res_i^2)) / (df_i * (sum(w_ok)^2)) * sum(ok_i) * ac_infl)
            B_int[si, li] <- b_hat
            SE_int[si, li] <- se_hat
            df_l_vec[li] <- df_i
            if (!is.null(score_mat)) {
              k_ok <- paste(dt_sul$run[ok_i], dt_sul$trial[ok_i], sep = "_")
              score_mat[match(k_ok, rt_keys), li] <- w_ok * res_i
            }
          }
          df_int[si] <- if (any(df_l_vec > 0)) stats::median(df_l_vec[df_l_vec > 0]) else 10
          if (!is.null(score_mat)) {
            ok_rows <- rowSums(!is.finite(score_mat)) == 0L
            if (sum(ok_rows) >= 5L) {
              R_sc <- suppressWarnings(stats::cor(score_mat[ok_rows, , drop = FALSE]))
              if (all(is.finite(R_sc))) {
                sc_int_accum <- sc_int_accum + R_sc
                sc_int_n <- sc_int_n + 1L
              }
            }
          }
        }
        R_within_int <- if (sc_int_n >= 1L) 0.65 * (sc_int_accum / sc_int_n) + 0.35 * R_op_u else R_op_u
        diag(R_within_int) <- 1.0
        if ("intercept" %in% covariates) {
          panel_fit_cache[["intercept"]] <- list(
            kind = "within",
            B_mat = B_int,
            SE_mat = SE_int,
            df_sub = df_int,
            R_within = R_within_int
          )
        }

        # Fit Level-2 between-subject multiple regression jointly across cov_between_subj
        if (length(cov_between_subj) > 0L && S_u >= 3L) {
          p_b <- length(cov_between_subj)
          Z_subj <- matrix(NA_real_, S_u, p_b, dimnames = list(s_ids_u, cov_between_subj))
          for (si in seq_len(S_u)) {
            sid_i <- s_ids_u[si]
            dt_su1 <- dt_u_l1[as.character(dt_u_l1[[subj_col]]) == sid_i]
            if (nrow(dt_su1) == 0L) dt_su1 <- dt_u[as.character(dt_u[[subj_col]]) == sid_i]
            for (jb in seq_len(p_b)) {
              zv <- as.numeric(dt_su1[[cov_between_subj[jb]]])
              zv <- zv[is.finite(zv)]
              if (length(zv)) Z_subj[si, jb] <- mean(zv)
            }
          }
          B_bs <- matrix(NA_real_, p_b, L, dimnames = list(cov_between_subj, as.character(lags)))
          SE_bs <- matrix(NA_real_, p_b, L, dimnames = list(cov_between_subj, as.character(lags)))
          Res_bs <- matrix(NA_real_, S_u, L, dimnames = list(s_ids_u, as.character(lags)))
          S_eff_bs <- integer(L)
          err_share_bs <- rep(NA_real_, L)
          for (li in seq_len(L)) {
            yi <- B_int[, li]
            s2i <- SE_int[, li]^2
            s2i[!is.finite(s2i) | s2i <= 0] <- stats::var(yi, na.rm = TRUE)
            ok_s <- is.finite(yi) & (rowSums(!is.finite(Z_subj)) == 0L)
            n_s_ok <- sum(ok_s)
            if (n_s_ok <= p_b + 1L) next
            S_eff_bs[li] <- n_s_ok
            y_ok <- yi[ok_s]
            Z_ok <- sweep(Z_subj[ok_s, , drop = FALSE], 2L, colMeans(Z_subj[ok_s, , drop = FALSE]), "-")
            X_des <- cbind(1.0, Z_ok)
            # Preliminary OLS to estimate residual between-subject variance tau2_b
            b_ols <- tryCatch(qr.solve(X_des, y_ok), error = function(e) NULL)
            if (is.null(b_ols)) next
            e_ols <- y_ok - as.numeric(X_des %*% b_ols)
            df_l2 <- max(1L, n_s_ok - p_b - 1L)
            v_res <- sum(e_ols^2) / df_l2
            mean_s2_ok <- mean(s2i[ok_s], na.rm = TRUE)
            tau2_b <- max(0.0, v_res - mean_s2_ok)
            err_share_bs[li] <- mean_s2_ok / max(tau2_b + mean_s2_ok, 1e-12)
            w_s <- 1.0 / pmax(tau2_b + s2i[ok_s], 1e-8)
            w_s <- w_s / mean(w_s)
            WX <- X_des * w_s
            xwx <- crossprod(X_des, WX) + 1e-10 * diag(p_b + 1L)
            xwx_inv <- tryCatch(solve(xwx), error = function(e) NULL)
            if (is.null(xwx_inv)) next
            b_wls <- as.numeric(xwx_inv %*% crossprod(WX, y_ok))
            e_wls <- y_ok - as.numeric(X_des %*% b_wls)
            Res_bs[ok_s, li] <- e_wls
            # HC2 leverage-adjusted sandwich variance
            h_ii <- pmin(0.85, rowSums((X_des %*% xwx_inv) * WX))
            meat_hc2 <- crossprod(WX, (WX * (e_wls^2) / pmax(1.0 - h_ii, 0.15)))
            cov_hc2 <- xwx_inv %*% meat_hc2 %*% xwx_inv
            meat_mod <- crossprod(WX, WX * (tau2_b + s2i[ok_s]))
            cov_mod <- xwx_inv %*% meat_mod %*% xwx_inv
            cov_tot <- if (df_l2 >= 3L) cov_hc2 else (0.70 * cov_hc2 + 0.30 * cov_mod)
            B_bs[, li] <- b_wls[-1L]
            SE_bs[, li] <- sqrt(pmax(diag(cov_tot)[-1L], .Machine$double.eps))
          }
          theta_err_bs <- stats::median(err_share_bs[is.finite(err_share_bs)], na.rm = TRUE)
          if (!is.finite(theta_err_bs)) theta_err_bs <- 0.50
          theta_err_bs <- min(1.0, max(0.0, theta_err_bs))
          R_prior_bs <- theta_err_bs * R_within_int + (1.0 - theta_err_bs) * R_res_u
          R_emp_bs <- suppressWarnings(stats::cor(Res_bs, use = "pairwise.complete.obs"))
          if (is.null(dim(R_emp_bs)) || any(dim(R_emp_bs) != c(L, L))) R_emp_bs <- R_prior_bs
          bad_r <- !is.finite(R_emp_bs)
          if (any(bad_r)) R_emp_bs[bad_r] <- R_prior_bs[bad_r]
          S_med_bs <- max(3L, as.integer(stats::median(S_eff_bs[S_eff_bs > 0], na.rm = TRUE)))
          sh_bs <- min(0.75, max(0.20, (S_med_bs - 2.0) / (S_med_bs + 6.0)))
          R_grp_bs <- 0.5 * ((sh_bs * R_emp_bs + (1.0 - sh_bs) * R_prior_bs) +
                             t(sh_bs * R_emp_bs + (1.0 - sh_bs) * R_prior_bs))
          diag(R_grp_bs) <- 1.0
          eig_rbs <- eigen(R_grp_bs, symmetric = TRUE)
          R_grp_bs <- eig_rbs$vectors %*% (pmax(eig_rbs$values, 0.02) * t(eig_rbs$vectors))
          d_rbs <- sqrt(pmax(diag(R_grp_bs), .Machine$double.eps))
          R_grp_bs <- R_grp_bs / outer(d_rbs, d_rbs)
          diag(R_grp_bs) <- 1.0
          df_bs <- max(2, S_med_bs - p_b - 1)
          for (jb in seq_len(p_b)) {
            cv_name <- cov_between_subj[jb]
            panel_fit_cache[[cv_name]] <- list(
              kind = "between_subj",
              b_raw_grp = as.numeric(B_bs[jb, ]),
              se_raw_grp = as.numeric(SE_bs[jb, ]),
              R_grp = R_grp_bs,
              R_within = R_within_int,
              shrink_alpha = sh_bs,
              df_grp = df_bs,
              n_subjects = S_med_bs
            )
          }
        }
      }

      # (B) Joint within-run trial-level multiple regression across ALL cov_within
      if (length(cov_within) > 0L) {
        p_w <- length(cov_within)
        all_have_w <- has_w_y && all(paste0(cov_within, "_whitened") %in% names(dt_u))
        B_w_list <- stats::setNames(
          lapply(cov_within, function(z) matrix(NA_real_, S_u, L, dimnames = list(s_ids_u, as.character(lags)))),
          cov_within
        )
        SE_w_list <- stats::setNames(
          lapply(cov_within, function(z) matrix(NA_real_, S_u, L, dimnames = list(s_ids_u, as.character(lags)))),
          cov_within
        )
        df_w_sub <- numeric(S_u)
        sc_w_accum <- stats::setNames(lapply(cov_within, function(z) matrix(0, L, L)), cov_within)
        sc_w_n <- stats::setNames(integer(p_w), cov_within)

        for (si in seq_len(S_u)) {
          sid_cur <- s_ids_u[si]
          dt_su <- if (!is.na(subj_col)) dt_u[as.character(dt_u[[subj_col]]) == sid_cur] else dt_u
          has_rt <- all(c("run", "trial") %in% names(dt_su))
          rt_keys <- if (has_rt) unique(paste(dt_su$run, dt_su$trial, sep = "_")) else NULL
          score_mats <- if (has_rt) lapply(cov_within, function(z) matrix(NA_real_, length(rt_keys), L)) else NULL
          df_l_vec <- numeric(L)

          for (li in seq_len(L)) {
            lag_val <- lags[li]
            dt_sul <- dt_su[abs(as.numeric(dt_su[[lag_col]]) - lag_val) < 1e-6]
            if (nrow(dt_sul) < p_w + 3L) next
            run_ops_ul <- if (!all_have_w) get_run_trial_ops(sid_cur, u_id, li) else NULL
            has_phi_ops <- !is.null(run_ops_ul) && any(vapply(run_ops_ul, function(ro) !is.null(ro$Phi_zca), logical(1L)))
            use_w_here <- all_have_w || (has_w_y && has_phi_ops)
            y_col_here <- if (use_w_here) "decon_demixed_whitened" else y_col_unw

            yv <- as.numeric(dt_sul[[y_col_here]])
            wv <- if (!is.null(w_col)) as.numeric(dt_sul[[w_col]]) else rep(1.0, length(yv))
            wv[!is.finite(wv) | wv <= 0] <- 1.0

            X_mat <- matrix(NA_real_, nrow(dt_sul), p_w)
            for (jw in seq_len(p_w)) {
              cn_j <- if (all_have_w) paste0(cov_within[jw], "_whitened") else cov_within[jw]
              X_mat[, jw] <- as.numeric(dt_sul[[cn_j]])
            }

            ok_i <- is.finite(yv) & (rowSums(!is.finite(X_mat)) == 0L)
            if (sum(ok_i) < p_w + 3L) next
            w_ok <- wv[ok_i] / mean(wv[ok_i])
            y_ok <- yv[ok_i]
            X_ok <- X_mat[ok_i, , drop = FALSE]

            # Apply run-level post-demixing operator M_post (+ optional Phi_zca) when covariates were merged post-Stage-3,
            # followed by exact weighted within-run centering Q_wr on both y and X
            if ("run" %in% names(dt_sul)) {
              r_ok <- as.character(dt_sul$run[ok_i])
              tr_ok <- if ("trial" %in% names(dt_sul)) as.integer(dt_sul$trial[ok_i]) else NULL
              u_runs <- unique(r_ok)
              for (rk in u_runs) {
                idx_r <- which(r_ok == rk)
                m_r <- length(idx_r)
                if (m_r < 2L) {
                  y_ok[idx_r] <- 0.0
                  X_ok[idx_r, ] <- 0.0
                  next
                }
                if (!all_have_w && !is.null(run_ops_ul) && !is.null(run_ops_ul[[rk]])) {
                  ro <- run_ops_ul[[rk]]
                  ok_loc_idx <- if (is.logical(ro$ok_r_local)) which(ro$ok_r_local) else as.integer(ro$ok_r_local)
                  if (!is.null(tr_ok) && length(ok_loc_idx) == m_r && all(tr_ok[idx_r] == ok_loc_idx) &&
                      is.matrix(ro$M_post) && all(dim(ro$M_post) == c(m_r, m_r))) {
                    Xr_c0 <- sweep(X_ok[idx_r, , drop = FALSE], 2L, colMeans(X_ok[idx_r, , drop = FALSE]), "-")
                    Xr_til <- ro$M_post %*% Xr_c0
                    if (use_w_here && is.matrix(ro$Phi_zca) && all(dim(ro$Phi_zca) == c(m_r, m_r))) {
                      X_ok[idx_r, ] <- ro$Phi_zca %*% Xr_til
                    } else {
                      X_ok[idx_r, ] <- Xr_til
                    }
                  }
                }
                wr <- w_ok[idx_r]
                sw <- sum(wr)
                y_ok[idx_r] <- y_ok[idx_r] - sum(wr * y_ok[idx_r]) / sw
                mu_xr <- colSums(X_ok[idx_r, , drop = FALSE] * wr) / sw
                X_ok[idx_r, ] <- sweep(X_ok[idx_r, , drop = FALSE], 2L, mu_xr, "-")
              }
              n_runs_i <- length(u_runs)
            } else {
              sw <- sum(w_ok)
              y_ok <- y_ok - sum(w_ok * y_ok) / sw
              mu_x <- colSums(X_ok * w_ok) / sw
              X_ok <- sweep(X_ok, 2L, mu_x, "-")
              n_runs_i <- 1L
            }

            w_ok <- w_ok / mean(w_ok)
            WX <- X_ok * w_ok
            xwx <- crossprod(X_ok, WX)
            col_e <- diag(xwx)
            if (any(!is.finite(col_e)) || all(col_e <= 1e-10)) next
            reg_diag <- 1e-10 * max(mean(col_e[col_e > 0]), 1.0)
            xwx_inv <- tryCatch(solve(xwx + reg_diag * diag(p_w)), error = function(e) NULL)
            if (is.null(xwx_inv)) next

            b_hat <- as.numeric(xwx_inv %*% crossprod(WX, y_ok))
            res_i <- y_ok - as.numeric(X_ok %*% b_hat)
            df_i <- max(1L, sum(ok_i) - p_w - n_runs_i)
            s2_ols <- sum(w_ok * (res_i^2)) / df_i
            cov_model <- s2_ols * xwx_inv

            h_ii <- pmin(0.85, rowSums((X_ok %*% xwx_inv) * WX))
            u_sc <- WX * (res_i / sqrt(pmax(1.0 - h_ii, 0.15)))
            meat_hc <- crossprod(u_sc)
            # When fitting on unwhitened readouts, add lag-1 Bartlett autocorrelation cross-term within run
            if (!use_w_here && sum(ok_i) >= 6L) {
              if ("run" %in% names(dt_sul)) {
                r_ok <- as.character(dt_sul$run[ok_i])
                for (rk in unique(r_ok)) {
                  ir <- which(r_ok == rk)
                  if (length(ir) >= 3L) {
                    U_r <- u_sc[ir, , drop = FALSE]
                    G1_r <- crossprod(U_r[-1L, , drop = FALSE], U_r[-length(ir), , drop = FALSE])
                    meat_hc <- meat_hc + 0.5 * (G1_r + t(G1_r))
                  }
                }
              } else {
                G1 <- crossprod(u_sc[-1L, , drop = FALSE], u_sc[-nrow(u_sc), , drop = FALSE])
                meat_hc <- meat_hc + 0.5 * (G1 + t(G1))
              }
            }
            cov_hc <- xwx_inv %*% meat_hc %*% xwx_inv
            cov_comb <- 0.5 * cov_model + 0.5 * cov_hc
            se_hat <- sqrt(pmax(diag(cov_comb), .Machine$double.eps))

            for (jw in seq_len(p_w)) {
              if (col_e[jw] > 1e-10) {
                B_w_list[[jw]][si, li] <- b_hat[jw]
                SE_w_list[[jw]][si, li] <- se_hat[jw]
              }
            }
            df_l_vec[li] <- df_i
            if (!is.null(score_mats)) {
              k_ok <- paste(dt_sul$run[ok_i], dt_sul$trial[ok_i], sep = "_")
              m_idx <- match(k_ok, rt_keys)
              infl_mat <- (WX * res_i) %*% xwx_inv
              for (jw in seq_len(p_w)) {
                score_mats[[jw]][m_idx, li] <- infl_mat[, jw]
              }
            }
          }
          df_w_sub[si] <- if (any(df_l_vec > 0)) stats::median(df_l_vec[df_l_vec > 0]) else 12
          if (!is.null(score_mats)) {
            for (jw in seq_len(p_w)) {
              sm_j <- score_mats[[jw]]
              ok_rows <- rowSums(!is.finite(sm_j)) == 0L
              if (sum(ok_rows) >= 5L) {
                R_sc <- suppressWarnings(stats::cor(sm_j[ok_rows, , drop = FALSE]))
                if (all(is.finite(R_sc))) {
                  sc_w_accum[[jw]] <- sc_w_accum[[jw]] + R_sc
                  sc_w_n[jw] <- sc_w_n[jw] + 1L
                }
              }
            }
          }
        }
        for (jw in seq_len(p_w)) {
          cv_name <- cov_within[jw]
          R_w_j <- if (sc_w_n[jw] >= 1L) 0.65 * (sc_w_accum[[jw]] / sc_w_n[jw]) + 0.35 * R_op_u else R_op_u
          diag(R_w_j) <- 1.0
          panel_fit_cache[[cv_name]] <- list(
            kind = "within",
            B_mat = B_w_list[[jw]],
            SE_mat = SE_w_list[[jw]],
            df_sub = df_w_sub,
            R_within = R_w_j
          )
        }
      }

      # (C) Between-run covariates (constant within run, varying across runs within subject)
      if (length(cov_between_run) > 0L && ("run" %in% names(dt_u))) {
        for (cv_br in cov_between_run) {
          B_br <- matrix(NA_real_, S_u, L, dimnames = list(s_ids_u, as.character(lags)))
          SE_br <- matrix(NA_real_, S_u, L, dimnames = list(s_ids_u, as.character(lags)))
          df_br <- rep(4, S_u)
          for (si in seq_len(S_u)) {
            sid_i <- s_ids_u[si]
            dt_su <- if (!is.na(subj_col)) dt_u[as.character(dt_u[[subj_col]]) == sid_i] else dt_u
            for (li in seq_len(L)) {
              lag_val <- lags[li]
              dt_sul <- dt_su[abs(as.numeric(dt_su[[lag_col]]) - lag_val) < 1e-6]
              if (nrow(dt_sul) < 4L) next
              rm_dt <- dt_sul[, .(
                ym = mean(as.numeric(.SD[[y_col_unw]]), na.rm = TRUE),
                xm = mean(as.numeric(.SD[[cv_br]]), na.rm = TRUE)
              ), by = run, .SDcols = c(y_col_unw, cv_br)]
              rm_dt <- rm_dt[is.finite(ym) & is.finite(xm)]
              if (nrow(rm_dt) >= 3L && stats::sd(rm_dt$xm) > 1e-10) {
                fit_r <- stats::lm(ym ~ xm, data = rm_dt)
                cf_r <- summary(fit_r)$coefficients
                if ("xm" %in% rownames(cf_r)) {
                  B_br[si, li] <- cf_r["xm", "Estimate"]
                  SE_br[si, li] <- cf_r["xm", "Std. Error"]
                  df_br[si] <- max(1L, nrow(rm_dt) - 2L)
                }
              }
            }
          }
          panel_fit_cache[[cv_br]] <- list(
            kind = "within",
            B_mat = B_br,
            SE_mat = SE_br,
            df_sub = df_br,
            R_within = R_op_u
          )
        }
      }
    }

    for (tm in covariates) {
      is_between_subj_term <- FALSE
      if (identical(input_kind, "panel")) {
        pf <- panel_fit_cache[[tm]]
        if (is.null(pf)) next
        if (identical(pf$kind, "between_subj")) {
          is_between_subj_term <- TRUE
          b_raw_grp <- pf$b_raw_grp
          se_raw_grp <- pf$se_raw_grp
          if (any(!is.finite(b_raw_grp))) {
            bad_b <- !is.finite(b_raw_grp)
            b_raw_grp[bad_b] <- 0.0
            med_se_b <- stats::median(se_raw_grp[is.finite(se_raw_grp) & se_raw_grp > 0], na.rm = TRUE)
            se_raw_grp[bad_b | !is.finite(se_raw_grp) | se_raw_grp <= 0] <- if (is.finite(med_se_b)) 2.0 * med_se_b else 0.25
          }
          R_within <- pf$R_within
          R_grp <- pf$R_grp
          shrink_alpha <- pf$shrink_alpha
          Sigma_raw_grp <- outer(se_raw_grp, se_raw_grp) * R_grp
          df_grp <- pf$df_grp
          S_n <- pf$n_subjects
          B_mat_orig <- NULL
        } else {
          B_mat <- pf$B_mat
          SE_mat <- pf$SE_mat
          df_sub <- pf$df_sub
          R_within <- pf$R_within
        }
      } else {
        # Coefficient table input
        dt_ut <- dt_u[as.character(dt_u$term) == tm]
        if (nrow(dt_ut) == 0L) next
        s_ids <- if (!is.na(subj_col)) unique(as.character(dt_ut[[subj_col]])) else "sub_1"
        S_n <- length(s_ids)
        B_mat <- matrix(NA_real_, S_n, L, dimnames = list(s_ids, as.character(lags)))
        SE_mat <- matrix(NA_real_, S_n, L, dimnames = list(s_ids, as.character(lags)))
        df_sub <- rep(24, S_n)
        for (si in seq_len(S_n)) {
          sid_i <- s_ids[si]
          dt_s <- if (!is.na(subj_col)) dt_ut[as.character(dt_ut[[subj_col]]) == sid_i] else dt_ut
          ord <- match(lags, as.numeric(dt_s[[lag_col]]))
          B_mat[si, ] <- as.numeric(dt_s[[est_col]])[ord]
          if (!is.na(se_col)) {
            SE_mat[si, ] <- as.numeric(dt_s[[se_col]])[ord]
          } else {
            SE_mat[si, ] <- 0.10
          }
          if ("df" %in% names(dt_s)) {
            df_vals <- as.numeric(dt_s$df)[ord]
            if (any(is.finite(df_vals) & df_vals > 0)) {
              df_sub[si] <- stats::median(df_vals[is.finite(df_vals) & df_vals > 0])
            }
          }
        }
        R_within <- R_op_u
      }

      if (!is_between_subj_term) {
        # Retain all subjects with at least 1 observed lag (available-case pooling instead of listwise deletion)
        ok_subj <- rowSums(is.finite(B_mat)) >= 1L
        if (!any(ok_subj)) next
        B_mat <- B_mat[ok_subj, , drop = FALSE]
        SE_mat <- SE_mat[ok_subj, , drop = FALSE]
        df_sub <- df_sub[ok_subj]
        s_ids_ok <- rownames(B_mat)
        S_n <- nrow(B_mat)
        B_mat_orig <- B_mat
        SE_mat_orig <- SE_mat

        for (li in seq_len(L)) {
          ok_b_li <- is.finite(B_mat[, li])
          bad_se <- ok_b_li & (!is.finite(SE_mat[, li]) | SE_mat[, li] <= 0)
          if (any(bad_se)) {
            sd_li <- if (sum(ok_b_li) >= 2L) stats::sd(B_mat[ok_b_li, li]) else NA_real_
            SE_mat[bad_se, li] <- if (is.finite(sd_li) && sd_li > 0) sd_li else 0.10
          }
        }

        # 3. Group / Unit Mean Trajectory & Full L x L Cross-Lag Covariance (available-case per lag)
        if (S_n >= 2L) {
          b_raw_grp <- rep(NA_real_, L)
          se_raw_grp <- rep(NA_real_, L)
          S_per_lag <- integer(L)
          err_share_grp <- rep(NA_real_, L)
          for (li in seq_len(L)) {
            ok_li <- is.finite(B_mat[, li]) & is.finite(SE_mat[, li]) & (SE_mat[, li] > 0)
            n_li <- sum(ok_li)
            S_per_lag[li] <- n_li
            if (n_li >= 2L) {
              yi <- B_mat[ok_li, li]
              s2i <- SE_mat[ok_li, li]^2
              v_tot <- stats::var(yi)
              mean_s2_li <- mean(s2i)
              tau2_b <- max(0, v_tot - mean_s2_li)
              err_share_grp[li] <- mean_s2_li / max(tau2_b + mean_s2_li, 1e-12)
              wi <- 1 / pmax(tau2_b + s2i, 1e-8)
              wi <- wi / sum(wi)
              b_raw_grp[li] <- sum(wi * yi)
              res_y <- yi - b_raw_grp[li]
              sum_w2_li <- sum(wi^2)
              denom_li <- pmin(1.0 - wi, pmax(1.0 - 2.0 * wi + sum_w2_li, 1.0 / n_li))
              v_emp <- sum((wi^2) * (res_y^2) / pmax(denom_li, 1.0 / (2.0 * n_li)))
              v_mod <- sum((wi^2) * (tau2_b + s2i))
              v_tot_grp <- if (n_li >= 3L) v_emp else (0.50 * v_emp + 0.50 * v_mod)
              se_raw_grp[li] <- sqrt(max(v_tot_grp, .Machine$double.eps))
            } else if (n_li == 1L) {
              b_raw_grp[li] <- B_mat[ok_li, li]
              se_raw_grp[li] <- SE_mat[ok_li, li]
            }
          }
          if (!any(is.finite(b_raw_grp))) next
          if (any(!is.finite(b_raw_grp))) {
            bad_l <- !is.finite(b_raw_grp)
            b_raw_grp[bad_l] <- 0.0
            med_se_g <- stats::median(se_raw_grp[is.finite(se_raw_grp) & se_raw_grp > 0], na.rm = TRUE)
            se_raw_grp[bad_l] <- if (is.finite(med_se_g)) 2.5 * med_se_g else 0.25
          }
          theta_err_grp <- stats::median(err_share_grp[is.finite(err_share_grp)], na.rm = TRUE)
          if (!is.finite(theta_err_grp)) theta_err_grp <- 0.50
          theta_err_grp <- min(1.0, max(0.0, theta_err_grp))
          R_prior_grp <- theta_err_grp * R_within + (1.0 - theta_err_grp) * R_res_u
          B_cent <- sweep(B_mat, 2L, b_raw_grp, "-")
          R_emp_grp <- suppressWarnings(stats::cor(B_cent, use = "pairwise.complete.obs"))
          if (is.null(dim(R_emp_grp)) || any(dim(R_emp_grp) != c(L, L))) R_emp_grp <- R_prior_grp
          bad_rg <- !is.finite(R_emp_grp)
          if (any(bad_rg)) R_emp_grp[bad_rg] <- R_prior_grp[bad_rg]
          S_eff_med <- max(2L, as.integer(stats::median(S_per_lag[S_per_lag > 0], na.rm = TRUE)))
          shrink_alpha <- min(0.75, max(0.20, (S_eff_med - 2.0) / (S_eff_med + 6.0)))
          R_grp <- shrink_alpha * R_emp_grp + (1.0 - shrink_alpha) * R_prior_grp
          R_grp <- 0.5 * (R_grp + t(R_grp))
          diag(R_grp) <- 1.0
          eig_R <- eigen(R_grp, symmetric = TRUE)
          R_grp <- eig_R$vectors %*% (pmax(eig_R$values, 0.02) * t(eig_R$vectors))
          d_rg <- sqrt(pmax(diag(R_grp), .Machine$double.eps))
          R_grp <- R_grp / outer(d_rg, d_rg)
          diag(R_grp) <- 1.0
          Sigma_raw_grp <- outer(se_raw_grp, se_raw_grp) * R_grp
          df_grp <- max(2, S_eff_med - 1)
        } else {
          b_raw_grp <- as.numeric(B_mat[1L, ])
          se_raw_grp <- as.numeric(SE_mat[1L, ])
          if (!any(is.finite(b_raw_grp))) next
          if (any(!is.finite(b_raw_grp))) {
            bad_l <- !is.finite(b_raw_grp)
            b_raw_grp[bad_l] <- 0.0
            med_se_1 <- stats::median(se_raw_grp[is.finite(se_raw_grp) & se_raw_grp > 0], na.rm = TRUE)
            se_raw_grp[bad_l | !is.finite(se_raw_grp) | se_raw_grp <= 0] <- if (is.finite(med_se_1)) 2.5 * med_se_1 else 0.25
          }
          shrink_alpha <- 0.0
          R_grp <- R_within
          Sigma_raw_grp <- outer(se_raw_grp, se_raw_grp) * R_grp
          df_grp <- if (is.finite(df_sub[1L]) && df_sub[1L] >= 2) df_sub[1L] else 20
        }
      }

      # 4. Joint Cross-Lag Spectral De-Blurring + Correlated-SURE Boosted P-Spline
      sm_grp <- medusa_smooth_deblur_trajectory(
        beta_vec = b_raw_grp,
        Sigma_beta = Sigma_raw_grp,
        lags = lags,
        C_sym = C_u,
        tr = tr,
        demix_within_trial = demix_within_trial,
        recalibrate_wt_lobe = recalibrate_wt_lobe,
        trajectory_mode = trajectory_mode
      )
      b_ml_grp <- sm_grp$beta
      se_ml_grp <- sm_grp$se
      Sigma_ml_grp <- sm_grp$Sigma
      b_deblur_grp <- sm_grp$beta_deblur
      Sigma_deblur_grp <- sm_grp$Sigma_deblur

      # 5. Subject-Level Smooth and Multivariate Empirical Bayes (EB) Trajectories
      if (!is_between_subj_term && S_n >= 1L) {
        B_smooth_subj <- matrix(0, S_n, L)
        Sigma_e_subj <- vector("list", S_n)
        for (si in seq_len(S_n)) {
          ops_si <- get_unit_operators(u_id, sid_opt = s_ids_ok[si])
          b_si_in <- B_mat[si, ]
          se_si_in <- SE_mat[si, ]
          miss_si <- !is.finite(b_si_in) | !is.finite(se_si_in) | (se_si_in <= 0)
          if (any(miss_si)) {
            b_si_in[miss_si] <- b_raw_grp[miss_si]
            max_se_si <- max(se_si_in[!miss_si], na.rm = TRUE)
            se_si_in[miss_si] <- if (is.finite(max_se_si) && max_se_si > 0) 3.0 * max_se_si else 0.50
          }
          R_si_use <- if (identical(ops_si$R_op, R_op_u)) R_within else (0.5 * R_within + 0.5 * ops_si$R_op)
          diag(R_si_use) <- 1.0
          Sigma_raw_si <- outer(se_si_in, se_si_in) * R_si_use
          sm_si <- medusa_smooth_deblur_trajectory(
            beta_vec = b_si_in,
            Sigma_beta = Sigma_raw_si,
            lags = lags,
            C_sym = ops_si$C_sym,
            tr = tr,
            demix_within_trial = demix_within_trial,
            recalibrate_wt_lobe = recalibrate_wt_lobe,
            trajectory_mode = trajectory_mode
          )
          B_smooth_subj[si, ] <- sm_si$beta
          Sigma_e_subj[[si]] <- sm_si$Sigma
        }

        if (S_n >= 2L) {
          Sigma_tot_sub <- stats::cov(B_smooth_subj)
          Sigma_e_bar <- Reduce("+", Sigma_e_subj) / S_n
          Diff_b <- 0.5 * ((Sigma_tot_sub - Sigma_e_bar) + t(Sigma_tot_sub - Sigma_e_bar))
          eig_b <- eigen(Diff_b, symmetric = TRUE)
          floor_b <- 0.05 * stats::median(diag(Sigma_e_bar))
          Sigma_b_hat <- eig_b$vectors %*% (pmax(eig_b$values, floor_b) * t(eig_b$vectors))

          B_eb_subj <- matrix(0, S_n, L)
          for (si in seq_len(S_n)) {
            S_tot_i <- Sigma_b_hat + Sigma_e_subj[[si]] + 1e-6 * diag(L)
            K_eb_i <- Sigma_b_hat %*% solve(S_tot_i)
            B_eb_subj[si, ] <- b_ml_grp + as.numeric(K_eb_i %*% (B_smooth_subj[si, ] - b_ml_grp))
          }
        } else {
          B_eb_subj <- B_smooth_subj
        }

        for (si in seq_len(S_n)) {
          df_s_rows <- data.table::data.table(
            id = s_ids_ok[si],
            unit = u_id,
            term = tm,
            lag = lags,
            estimate_raw = B_mat_orig[si, ],
            std.error_raw = SE_mat_orig[si, ],
            estimate_smooth = B_smooth_subj[si, ],
            estimate_eb = B_eb_subj[si, ]
          )
          df_s_rows[[unit_col]] <- u_id
          subj_out_list[[length(subj_out_list) + 1L]] <- df_s_rows

          lat_si <- medusa_spline_latency(
            beta_vec = B_eb_subj[si, ],
            Sigma_beta = Sigma_e_subj[[si]],
            lags = lags,
            fine_dt = fine_dt,
            n_boot = 0L,
            polarity = polarity
          )
          dt_lat_si <- data.table::data.table(
            id = s_ids_ok[si],
            unit = u_id,
            term = tm,
            peak_lag = lat_si$peak_lag,
            peak_lag_grid = lat_si$peak_lag_grid,
            peak_val = lat_si$peak_val,
            com_lag = lat_si$com_lag,
            fwhm = lat_si$fwhm
          )
          dt_lat_si[[unit_col]] <- u_id
          subj_lat_list[[length(subj_lat_list) + 1L]] <- dt_lat_si
        }
      }

      # 6. Operator-Informed Simultaneous Cross-Lag Inference (Studentized Max-|T|, 1D TFCE, & Pivotal Omnibus)
      t_obs <- b_raw_grp / pmax(se_raw_grp, .Machine$double.eps)
      p_point <- 2.0 * stats::pt(abs(t_obs), df = df_grp, lower.tail = FALSE)
      p_bh <- stats::p.adjust(p_point, method = "BH")
      p_by <- stats::p.adjust(p_point, method = "BY")

      t_ml_obs <- b_ml_grp / pmax(se_ml_grp, .Machine$double.eps)
      p_ml_point <- 2.0 * stats::pt(abs(t_ml_obs), df = df_grp, lower.tail = FALSE)

      eig_Rg <- eigen(R_grp, symmetric = TRUE)
      L_Rg <- eig_Rg$vectors %*% (sqrt(pmax(eig_Rg$values, 0)) * t(eig_Rg$vectors))

      # Studentized finite-sample null simulation: re-estimates per-lag SEs and shrunk correlation R_grp^(b)
      # under H0 so Max-|T|, TFCE, and Hotelling/Mahalanobis omnibus_p account for Wishart covariance uncertainty
      m_df_sim <- max(3L, min(80L, as.integer(round(df_grp)) + 1L))
      Z_rep <- matrix(stats::rnorm(m_df_sim * n_perm * L), nrow = m_df_sim * n_perm, ncol = L) %*% t(L_Rg)
      dim(Z_rep) <- c(m_df_sim, n_perm, L)
      # Vectorized sample means and sample variances across m_df_sim rows for all (n_perm x L)
      Z_perm_L <- matrix(aperm(Z_rep, c(1L, 3L, 2L)), nrow = m_df_sim, ncol = L * n_perm)
      mu_null_vec <- colMeans(Z_perm_L)
      Z_cent_L <- sweep(Z_perm_L, 2L, mu_null_vec, "-")
      se_null_vec <- sqrt(pmax(colSums(Z_cent_L^2) / (m_df_sim * (m_df_sim - 1.0)), .Machine$double.eps))
      T_null <- matrix(mu_null_vec / se_null_vec, nrow = L, ncol = n_perm)

      max_t_null <- apply(abs(T_null), 2L, max)
      crit_maxT <- as.numeric(stats::quantile(max_t_null, probs = conf_level, na.rm = TRUE, names = FALSE))
      p_maxT <- vapply(abs(t_obs), function(tv) (1.0 + sum(max_t_null >= tv)) / (n_perm + 1.0), numeric(1L))

      # 1D Threshold-Free Cluster Enhancement (TFCE) calibrated under correlated studentized T_null
      T_all <- cbind(t_obs, T_null)
      tfce_all <- medusa_tfce_1d_matrix(T_all, dt = dt_lag, E = 0.5, H = 2.0, dh = 0.10, two_sided = TRUE)
      tfce_obs <- tfce_all[, 1L]
      max_tfce_null <- apply(tfce_all[, -1L, drop = FALSE], 2L, max)
      p_tfce <- vapply(tfce_obs, function(sv) (1.0 + sum(max_tfce_null >= sv)) / (n_perm + 1.0), numeric(1L))

      # Pivotal regularized studentized Mahalanobis / Hotelling epoch omnibus test against H0: beta(tau) == 0
      compute_omni_stat <- function(t_vec, R_mat) {
        eg <- eigen(R_mat, symmetric = TRUE)
        l_fl <- 0.05 * max(eg$values, 1e-8)
        R_inv <- eg$vectors %*% ((1.0 / pmax(eg$values, l_fl)) * t(eg$vectors))
        sum(t_vec * as.numeric(R_inv %*% t_vec)) / L
      }
      q_obs <- compute_omni_stat(t_obs, R_grp)
      B_omni <- min(as.integer(n_perm), 2000L)
      q_null <- numeric(B_omni)
      if (shrink_alpha > 0 && m_df_sim <= 60L) {
        for (bi in seq_len(B_omni)) {
          idx_bi <- ((bi - 1L) * L + 1L):(bi * L)
          Zc_b <- Z_cent_L[, idx_bi, drop = FALSE]
          sd_b <- sqrt(pmax(colSums(Zc_b^2), .Machine$double.eps))
          Zn_b <- sweep(Zc_b, 2L, sd_b, "/")
          R_emp_b <- crossprod(Zn_b)
          R_b <- shrink_alpha * R_emp_b + (1.0 - shrink_alpha) * R_grp
          R_b <- 0.5 * (R_b + t(R_b))
          diag(R_b) <- 1.0
          eg_rb <- eigen(R_b, symmetric = TRUE)
          if (min(eg_rb$values) < 0.02) {
            R_b <- eg_rb$vectors %*% (pmax(eg_rb$values, 0.02) * t(eg_rb$vectors))
            d_b <- sqrt(pmax(diag(R_b), .Machine$double.eps))
            R_b <- R_b / outer(d_b, d_b)
            diag(R_b) <- 1.0
          }
          q_null[bi] <- compute_omni_stat(T_null[, bi], R_b)
        }
      } else {
        eg_r0 <- eigen(R_grp, symmetric = TRUE)
        l_fl0 <- 0.05 * max(eg_r0$values, 1e-8)
        R_inv0 <- eg_r0$vectors %*% ((1.0 / pmax(eg_r0$values, l_fl0)) * t(eg_r0$vectors))
        T_sub <- T_null[, seq_len(B_omni), drop = FALSE]
        q_null <- colSums(T_sub * (R_inv0 %*% T_sub)) / L
      }
      omnibus_p <- (1.0 + sum(q_null >= q_obs)) / (B_omni + 1.0)

      # Marra & Wood (2012) smoothing-bias-aware SE for simultaneous confidence bands around estimate_multilag
      IL_mat <- diag(1, L)
      ImS_grp <- IL_mat - sm_grp$S_smooth
      diff_smooth_grp <- as.numeric(ImS_grp %*% b_deblur_grp)
      var_diff_null <- pmax(0, diag(ImS_grp %*% Sigma_deblur_grp %*% t(ImS_grp)))
      bias2_mw <- pmax(0, (diff_smooth_grp^2) - 0.5 * var_diff_null) +
        0.25 * pmax(0, diag(Sigma_deblur_grp) - diag(Sigma_ml_grp))
      se_ml_band <- sqrt(pmax(se_ml_grp^2 + bias2_mw, .Machine$double.eps))

      t_crit_pt <- stats::qt(1.0 - (1.0 - conf_level) / 2.0, df = df_grp)
      dt_coef <- data.table::data.table(
        unit = u_id,
        term = tm,
        lag = lags,
        evt_time = lags,
        estimate = b_raw_grp,
        std.error = se_raw_grp,
        estimate_multilag = b_ml_grp,
        std.error_multilag = se_ml_grp,
        statistic = t_obs,
        statistic_multilag = t_ml_obs,
        df = df_grp,
        p.value = p_point,
        p.value_multilag = p_ml_point,
        p_fdr_bh = p_bh,
        p_fdr_by = p_by,
        p_maxT = p_maxT,
        tfce_stat = tfce_obs,
        p_tfce = p_tfce,
        crit_maxT = crit_maxT,
        conf.low = b_raw_grp - t_crit_pt * se_raw_grp,
        conf.high = b_raw_grp + t_crit_pt * se_raw_grp,
        conf.low.sim = b_raw_grp - crit_maxT * se_raw_grp,
        conf.high.sim = b_raw_grp + crit_maxT * se_raw_grp,
        conf.low.multilag.sim = b_ml_grp - crit_maxT * se_ml_band,
        conf.high.multilag.sim = b_ml_grp + crit_maxT * se_ml_band,
        omnibus_stat = q_obs,
        omnibus_p = omnibus_p,
        n_subjects = S_n
      )
      dt_coef[[unit_col]] <- u_id
      coef_out_list[[length(coef_out_list) + 1L]] <- dt_coef

      # 7. Continuous Sub-Grid Spline Peak Latency, Center-of-Mass, & FWHM
      lat_grp <- medusa_spline_latency(
        beta_vec = b_ml_grp,
        Sigma_beta = Sigma_ml_grp,
        lags = lags,
        fine_dt = fine_dt,
        n_boot = n_boot_latency,
        conf_level = conf_level,
        polarity = polarity
      )
      dt_lat <- data.table::data.table(
        unit = u_id,
        term = tm,
        peak_lag = lat_grp$peak_lag,
        peak_lag_grid = lat_grp$peak_lag_grid,
        peak_lag_se = lat_grp$peak_lag_se,
        peak_lag_conf.low = lat_grp$peak_lag_conf.low,
        peak_lag_conf.high = lat_grp$peak_lag_conf.high,
        peak_val = lat_grp$peak_val,
        peak_val_se = lat_grp$peak_val_se,
        com_lag = lat_grp$com_lag,
        com_lag_se = lat_grp$com_lag_se,
        fwhm = lat_grp$fwhm,
        fwhm_se = lat_grp$fwhm_se,
        omnibus_p = omnibus_p,
        min_p_maxT = min(p_maxT, na.rm = TRUE),
        min_p_tfce = min(p_tfce, na.rm = TRUE)
      )
      dt_lat[[unit_col]] <- u_id
      lat_out_list[[length(lat_out_list) + 1L]] <- dt_lat

      # 8. Covariance-Propagated Multi-Lag Contrasts (evaluated on unsmoothed deblurred trajectory by default)
      b_con_src <- if (identical(contrast_source, "unsmoothed")) {
        b_deblur_grp
      } else if (identical(contrast_source, "raw")) {
        b_raw_grp
      } else {
        b_ml_grp
      }
      Sig_con_src <- if (identical(contrast_source, "unsmoothed")) {
        Sigma_deblur_grp
      } else if (identical(contrast_source, "raw")) {
        Sigma_raw_grp
      } else {
        Sigma_ml_grp
      }

      con_specs <- contrasts
      if (is.null(con_specs) && L >= 4L) {
        mid_tau <- stats::median(lags)
        early_idx <- which(lags < mid_tau)
        late_idx <- which(lags > mid_tau)
        if (length(early_idx) >= 1L && length(late_idx) >= 1L) {
          c_default <- numeric(L)
          c_default[early_idx] <- -1.0 / length(early_idx)
          c_default[late_idx] <- 1.0 / length(late_idx)
          con_specs <- list(late_minus_early = c_default)
        }
      }
      if (!is.null(con_specs) && length(con_specs) > 0L) {
        c_names <- names(con_specs)
        if (is.null(c_names)) c_names <- paste0("contrast_", seq_along(con_specs))
        for (ci in seq_along(con_specs)) {
          cs <- con_specs[[ci]]
          if (is.numeric(cs) && length(cs) == L) {
            c_vec <- as.numeric(cs)
          } else if (is.list(cs) && !is.null(cs$early) && !is.null(cs$late)) {
            c_vec <- numeric(L)
            e_idx <- which(lags >= min(cs$early) - 1e-6 & lags <= max(cs$early) + 1e-6)
            l_idx <- which(lags >= min(cs$late) - 1e-6 & lags <= max(cs$late) + 1e-6)
            if (length(e_idx)) c_vec[e_idx] <- -1.0 / length(e_idx)
            if (length(l_idx)) c_vec[l_idx] <- 1.0 / length(l_idx)
          } else {
            next
          }
          est_c <- sum(c_vec * b_con_src)
          var_c <- as.numeric(crossprod(c_vec, Sig_con_src %*% c_vec))
          se_c <- sqrt(max(var_c, .Machine$double.eps))
          t_c <- est_c / se_c
          p_c <- 2.0 * stats::pt(abs(t_c), df = df_grp, lower.tail = FALSE)
          dt_c <- data.table::data.table(
            unit = u_id,
            term = tm,
            contrast = c_names[ci],
            estimate = est_c,
            std.error = se_c,
            statistic = t_c,
            df = df_grp,
            p.value = p_c,
            conf.low = est_c - t_crit_pt * se_c,
            conf.high = est_c + t_crit_pt * se_c
          )
          dt_c[[unit_col]] <- u_id
          con_out_list[[length(con_out_list) + 1L]] <- dt_c
        }
      }

      cov_out_list[[u_id]][[tm]] <- list(
        Sigma_raw = Sigma_raw_grp,
        Sigma_deblur = Sigma_deblur_grp,
        Sigma_multilag = Sigma_ml_grp,
        R_group = R_grp,
        C_sym = C_u,
        M_op = sm_grp$M_op,
        eta_opt = sm_grp$eta_opt
      )
    }
  }

  structure(
    list(
      coefficients = data.table::rbindlist(coef_out_list, fill = TRUE),
      latency = data.table::rbindlist(lat_out_list, fill = TRUE),
      contrasts = if (length(con_out_list)) data.table::rbindlist(con_out_list, fill = TRUE) else data.table::data.table(),
      subjects = if (length(subj_out_list)) data.table::rbindlist(subj_out_list, fill = TRUE) else NULL,
      subject_latency = if (length(subj_lat_list)) data.table::rbindlist(subj_lat_list, fill = TRUE) else NULL,
      covariance = cov_out_list,
      settings = list(
        lags = lags,
        trajectory_mode = trajectory_mode,
        contrast_source = contrast_source,
        demix_within_trial = demix_within_trial,
        n_perm = n_perm,
        n_boot_latency = n_boot_latency,
        fine_dt = fine_dt,
        conf_level = conf_level,
        polarity = polarity
      )
    ),
    class = "medusa_multilag"
  )
}


