# Subspace-Debiased Adaptive-POET + Multi-Lag Toeplitz & General-W_s
# LTZ-Gated Continuous Covariance Regression (v3 Architecture)
# and Post-Hoc Voxel-to-Parcel Panel Aggregation

.poet_firm_threshold <- function(x, t1, t2) {
  ax <- abs(x)
  out <- numeric(length(x))
  if (t2 <= t1) {
    keep <- is.finite(ax) & (ax > t1)
    out[keep] <- x[keep]
    return(out)
  }
  mid <- is.finite(ax) & (ax > t1) & (ax < t2)
  hi <- is.finite(ax) & (ax >= t2)
  if (any(mid)) {
    out[mid] <- sign(x[mid]) * t2 * ((ax[mid] - t1) / (t2 - t1))
  }
  if (any(hi)) {
    out[hi] <- x[hi]
  }
  out
}

.poet_firm_gate <- function(w_stat, lo = 2.2, hi = 6.5) {
  if (!is.finite(w_stat) || w_stat <= lo) return(0.0)
  if (w_stat >= hi) return(1.0)
  (w_stat - lo) / max(hi - lo, 1e-8)
}

.poet_project_cor_psd <- function(R_in, eig_floor = 1e-3, max_off = 0.985) {
  P <- nrow(R_in)
  R_sym <- 0.5 * (R_in + t(R_in))
  R_sym[!is.finite(R_sym)] <- 0.0
  off <- row(R_sym) != col(R_sym)
  R_sym[off] <- pmin(max_off, pmax(-max_off, R_sym[off]))
  diag(R_sym) <- 1.0
  dimnames(R_sym) <- dimnames(R_in)
  if (P <= 1L) return(R_sym)
  eg <- eigen(R_sym, symmetric = TRUE)
  if (min(eg$values) >= eig_floor) {
    return(R_sym)
  }
  R_psd <- eg$vectors %*% (pmax(eg$values, eig_floor) * t(eg$vectors))
  d <- sqrt(pmax(diag(R_psd), .Machine$double.eps))
  R_out <- R_psd / outer(d, d)
  R_out <- 0.5 * (R_out + t(R_out))
  R_out[off] <- pmin(max_off, pmax(-max_off, R_out[off]))
  diag(R_out) <- 1.0
  dimnames(R_out) <- dimnames(R_in)
  R_out
}

poet_shrink_spikes <- function(vals, vecs, P, df_eff, k_max = 12L, poet_mode = "adaptive_logp") {
  gamma_aspect <- P / df_eff
  tw_margin <- if (identical(poet_mode, "adaptive_logp")) {
    max(0.12, 2.2 * ((1.0 + sqrt(gamma_aspect))^(4.0 / 3.0)) * (max(gamma_aspect, 1e-4)^(1.0 / 6.0)) / (df_eff^(2.0 / 3.0)))
  } else {
    0.12
  }
  mp_edge <- (1.0 + sqrt(gamma_aspect))^2 + tw_margin
  k_cap <- min(as.integer(k_max), length(vals), as.integer(floor(P / 3L)))
  min_ipr <- if (identical(poet_mode, "adaptive_logp")) max(2.8, min(6.0, 0.12 * P)) else 2.25

  keep_f <- integer(0L)
  if (k_cap >= 1L) {
    for (k in seq_len(k_cap)) {
      if (vals[k] <= max(1.55, mp_edge)) break
      vk2 <- vecs[, k]^2
      ipr_n <- 1.0 / max(sum(vk2^2), 1e-12)
      top2_share <- if (P >= 6L) sum(utils::head(sort(vk2, decreasing = TRUE), 2L)) else 0.0
      top2_cap <- min(0.48, max(0.20, 3.5 / sqrt(P)))
      is_pervasive <- if (identical(poet_mode, "adaptive_logp") && P >= 8L) {
        (ipr_n > min_ipr) && (top2_share < top2_cap)
      } else {
        (ipr_n > 2.25)
      }
      if (is_pervasive) {
        keep_f <- c(keep_f, k)
      }
    }
  }

  K_sel <- length(keep_f)
  if (K_sel > 0L) {
    sum_kept <- sum(vals[keep_f])
    rem_cnt <- max(1L, P - K_sel)
    bulk_mean <- min(1.0, max(0.05, (P - sum_kept) / rem_cnt))
    if (length(vals) == P) {
      rem_idx <- setdiff(seq_len(P), keep_f)
      if (length(rem_idx)) bulk_mean <- min(1.0, max(0.05, mean(vals[rem_idx])))
    }
    if (identical(poet_mode, "adaptive_logp") && gamma_aspect > 0.02) {
      shrunk_vals <- numeric(K_sel)
      sp_edge <- 1.0 + sqrt(gamma_aspect)
      for (ki in seq_len(K_sel)) {
        ell_k <- vals[keep_f[ki]] / bulk_mean
        if (ell_k > sp_edge + 1e-6) {
          b_term <- ell_k + 1.0 - gamma_aspect
          disc <- max(0.0, b_term^2 - 4.0 * ell_k)
          tau_k <- 0.5 * (b_term + sqrt(disc)) - 1.0
          cos2_k <- if (tau_k > sqrt(gamma_aspect)) {
            (1.0 - gamma_aspect / (tau_k^2)) / (1.0 + gamma_aspect / tau_k)
          } else 0.0
          shrunk_vals[ki] <- max(0.0, bulk_mean * tau_k * cos2_k)
        } else {
          shrunk_vals[ki] <- max(0.0, vals[keep_f[ki]] - bulk_mean)
        }
      }
    } else {
      shrunk_vals <- pmax(0.0, vals[keep_f] - bulk_mean)
    }
    L_mat <- sweep(vecs[, keep_f, drop = FALSE], 2L, sqrt(shrunk_vals), "*")
    comm <- rowSums(L_mat^2)
    over <- comm > 0.96
    if (any(over)) {
      L_mat[over, ] <- sweep(L_mat[over, , drop = FALSE], 1L, sqrt(0.96 / comm[over]), "*")
    }
  } else {
    L_mat <- matrix(0.0, P, 0L)
  }
  list(K = K_sel, L = L_mat, keep_f = keep_f)
}

poet_regularize_cor <- function(
  R_raw,
  df_u,
  k_max = 12L,
  poet_mode = c("adaptive_logp", "fixed_v1")
) {
  poet_mode <- match.arg(poet_mode)
  P <- nrow(R_raw)
  dn <- dimnames(R_raw)
  R_sym <- 0.5 * (R_raw + t(R_raw))
  R_sym[!is.finite(R_sym)] <- 0.0
  off <- row(R_sym) != col(R_sym)
  R_sym[off] <- pmin(0.985, pmax(-0.985, R_sym[off]))
  diag(R_sym) <- 1.0
  df_eff <- max(20.0, as.numeric(df_u)[1L])
  df_z <- max(10.0, df_eff - 3.0)

  if (P <= 3L) {
    R_reg <- R_sym
    if (identical(poet_mode, "fixed_v1")) {
      R_reg[off] <- .poet_firm_threshold(R_sym[off], 1.65 / sqrt(df_eff), 3.00 / sqrt(df_eff))
    } else {
      z_off <- sqrt(df_z) * atanh(pmin(0.999, pmax(-0.999, R_sym[off])))
      z_thr <- .poet_firm_threshold(z_off, 1.68, 3.00)
      R_reg[off] <- tanh(z_thr / sqrt(df_z))
    }
    R_psd <- .poet_project_cor_psd(R_reg)
    return(list(R = R_psd, K = 0L, L = matrix(0.0, P, 0L), S_local = R_psd - diag(P)))
  }

  eg <- eigen(R_sym, symmetric = TRUE)
  vals <- pmax(eg$values, 0.0)
  vecs <- eg$vectors
  sp_fit <- poet_shrink_spikes(vals, vecs, P = P, df_eff = df_eff, k_max = k_max, poet_mode = poet_mode)
  K_sel <- sp_fit$K
  L_mat <- sp_fit$L
  LLt <- if (K_sel > 0L) tcrossprod(L_mat) else matrix(0.0, P, P)

  E_res <- R_sym - LLt
  d_res <- sqrt(pmax(0.04, 1.0 - diag(LLt)))
  C_res <- E_res / outer(d_res, d_res)
  C_res[!is.finite(C_res)] <- 0.0
  C_res[off] <- pmin(0.999, pmax(-0.999, C_res[off]))

  C_thr <- matrix(0.0, P, P)
  if (identical(poet_mode, "fixed_v1")) {
    t1 <- 1.80 / sqrt(df_eff)
    t2 <- 3.20 / sqrt(df_eff)
    C_thr[off] <- .poet_firm_threshold(C_res[off], t1, t2)
  } else {
    z1 <- max(1.85, sqrt(max(0.0, 3.0 * log(P) - 0.80)))
    z2 <- z1 + 1.35
    z_off <- sqrt(df_z) * atanh(C_res[off])
    z_thr <- .poet_firm_threshold(z_off, z1, z2)
    C_thr[off] <- tanh(z_thr / sqrt(df_z))
  }

  S_loc <- C_thr * outer(d_res, d_res)
  S_loc <- 0.5 * (S_loc + t(S_loc))
  diag(S_loc) <- 0.0

  R_comb <- LLt + S_loc
  diag(R_comb) <- 1.0
  dimnames(R_comb) <- dn
  R_psd <- .poet_project_cor_psd(R_comb)
  list(R = R_psd, K = K_sel, L = L_mat, S_local = S_loc)
}

.ltz_from_cor_quad <- function(z_vec, R_spatial, p_lo = 0.25, p_hi = 0.01) {
  ok <- is.finite(z_vec)
  z <- z_vec[ok]
  m_p <- length(z)
  if (m_p == 0L) return(0.0)
  R_sub <- R_spatial[ok, ok, drop = FALSE]
  Sigma_z <- R_sub * R_sub
  diag(Sigma_z) <- 1.0
  ev <- eigen(0.5 * (Sigma_z + t(Sigma_z)), symmetric = TRUE, only.values = TRUE)$values
  ev <- ev[is.finite(ev) & ev > 1e-10]
  if (!length(ev)) return(0.0)
  s1 <- sum(ev); s2 <- sum(ev^2); s3 <- sum(ev^3)
  c1 <- s3 / max(s2, 1e-12)
  nu_s <- max(0.5, (s2^3) / max(s3^2, 1e-18))
  c0 <- s1 - c1 * nu_s
  q_lo <- c0 + c1 * stats::qchisq(1.0 - p_lo, df = nu_s)
  q_hi <- c0 + c1 * stats::qchisq(1.0 - p_hi, df = nu_s)
  Q_val <- sum(z^2)
  .poet_firm_gate(Q_val, lo = q_lo, hi = max(q_lo + 0.5, q_hi))
}

.ltz_from_factor_quad <- function(z_vec, L_poet, p_lo = 0.25, p_hi = 0.01) {
  ok <- is.finite(z_vec)
  z <- z_vec[ok]
  V <- length(z)
  if (V == 0L) return(0.0)
  L_sub <- if (!is.null(L_poet) && nrow(L_poet) == length(z_vec) && ncol(L_poet) > 0L) {
    L_poet[ok, , drop = FALSE]
  } else {
    matrix(0.0, V, 0L)
  }
  K <- ncol(L_sub)
  s1 <- as.numeric(V)
  if (K == 0L) {
    s2 <- as.numeric(V)
    s3 <- as.numeric(V)
  } else {
    comm <- rowSums(L_sub^2)
    d_diag <- 1.0 - comm^2
    # Exact trace(Sigma_z^2) and trace(Sigma_z^3) in O(V K^2) without forming V x V!
    # Sigma_z = (L L') o (L L') + diag(1 - comm^2), where A o A = B B' with B_{v, (a,b)} = L_{v,a} L_{v,b}.
    B_kr <- matrix(0.0, V, K * K)
    col_i <- 0L
    for (a in seq_len(K)) {
      for (b in seq_len(K)) {
        col_i <- col_i + 1L
        B_kr[, col_i] <- L_sub[, a] * L_sub[, b]
      }
    }
    G_b <- crossprod(B_kr)
    tr_A2 <- sum(G_b^2)
    s2 <- tr_A2 + 2.0 * sum(d_diag * (comm^2)) + sum(d_diag^2)
    G_b_d <- crossprod(B_kr, d_diag * B_kr)
    tr_A3 <- sum(G_b * (G_b %*% G_b))
    tr_A2_D <- sum(G_b * G_b_d)
    s3 <- tr_A3 + 3.0 * tr_A2_D + 3.0 * sum((d_diag^2) * (comm^2)) + sum(d_diag^3)
  }
  c1 <- s3 / max(s2, 1e-12)
  nu_s <- max(0.5, (s2^3) / max(s3^2, 1e-18))
  c0 <- s1 - c1 * nu_s
  q_lo <- c0 + c1 * stats::qchisq(1.0 - p_lo, df = nu_s)
  q_hi <- c0 + c1 * stats::qchisq(1.0 - p_hi, df = nu_s)
  .poet_firm_gate(sum(z^2), lo = q_lo, hi = max(q_lo + 0.5, q_hi))
}

.project_trial_err_psd <- function(C_mat, var_floors = NULL, max_cor = 0.98) {
  k <- nrow(C_mat)
  C_sym <- 0.5 * (C_mat + t(C_mat))
  C_sym[!is.finite(C_sym)] <- 0.0
  d_raw <- diag(C_sym)
  d_targ <- if (!is.null(var_floors)) pmax(var_floors, d_raw, 0.0) else pmax(d_raw, 0.0)
  pos <- d_targ > 1e-16
  if (!any(pos)) return(matrix(0.0, k, k))
  C_out <- matrix(0.0, k, k)
  if (sum(pos) == 1L) {
    C_out[pos, pos] <- d_targ[pos]
    return(C_out)
  }
  idx <- which(pos)
  sd_v <- sqrt(d_targ[idx])
  R_sub <- C_sym[idx, idx, drop = FALSE] / outer(sd_v, sd_v)
  off <- row(R_sub) != col(R_sub)
  R_sub[off] <- pmin(max_cor, pmax(-max_cor, R_sub[off]))
  diag(R_sub) <- 1.0
  eg <- eigen(0.5 * (R_sub + t(R_sub)), symmetric = TRUE)
  R_psd <- eg$vectors %*% (pmax(eg$values, 1e-4) * t(eg$vectors))
  d_r <- sqrt(pmax(diag(R_psd), .Machine$double.eps))
  R_psd <- R_psd / outer(d_r, d_r)
  R_psd[off] <- pmin(max_cor, pmax(-max_cor, R_psd[off]))
  diag(R_psd) <- 1.0
  C_out[idx, idx] <- R_psd * outer(sd_v, sd_v)
  0.5 * (C_out + t(C_out))
}

.proj_cross_diag <- function(c_vec, M_L, M_R = M_L) {
  c_use <- ifelse(is.finite(c_vec), c_vec, 0.0)
  if (is.null(M_L) || ncol(M_L) == 0L || is.null(M_R) || ncol(M_R) == 0L) {
    return(c_use)
  }
  h_L <- rowSums(M_L^2)
  if (identical(M_L, M_R)) {
    cross_term <- rowSums((M_L %*% crossprod(M_L, c_use * M_L)) * M_L)
    return(c_use * (1.0 - 2.0 * h_L) + cross_term)
  }
  h_R <- rowSums(M_R^2)
  p_L <- ncol(M_L); p_R <- ncol(M_R)
  cross_term <- numeric(length(c_use))
  for (a in seq_len(p_L)) {
    la <- M_L[, a]
    for (b in seq_len(p_R)) {
      rb <- M_R[, b]
      cross_term <- cross_term + sum(la * rb * c_use) * (la * rb)
    }
  }
  c_use * (1.0 - h_L - h_R) + cross_term
}

.extract_psd_factors <- function(C_in, k_rank, nC = NULL, n_lags = 1L, probe_cache = NULL) {
  q <- nrow(C_in)
  C_sym <- 0.5 * (C_in + t(C_in))
  C_sym[!is.finite(C_sym)] <- 0.0
  d_ci <- diag(C_sym)
  pos_c <- which(d_ci > 1e-14)
  F_out <- matrix(0.0, k_rank, q)
  np <- length(pos_c)
  if (np == 0L) return(F_out)
  if (np == 1L) {
    F_out[1L, pos_c] <- sqrt(d_ci[pos_c])
    return(F_out)
  }

  sd_c <- sqrt(d_ci[pos_c])
  R_ci <- C_sym[pos_c, pos_c, drop = FALSE] / outer(sd_c, sd_c)
  off_c <- row(R_ci) != col(R_ci)
  R_ci[off_c] <- pmin(0.98, pmax(-0.98, R_ci[off_c]))
  diag(R_ci) <- 1.0
  R_sym_c <- 0.5 * (R_ci + t(R_ci))

  if (np <= k_rank) {
    eg_c <- eigen(R_sym_c, symmetric = TRUE)
    if (min(eg_c$values) < -1e-8) {
      R_p <- eg_c$vectors %*% (pmax(eg_c$values, 0.0) * t(eg_c$vectors))
      dr_p <- sqrt(pmax(diag(R_p), .Machine$double.eps))
      R_p <- R_p / outer(dr_p, dr_p)
      eg_c <- eigen(0.5 * (R_p + t(R_p)), symmetric = TRUE)
    }
    k_use <- min(k_rank, sum(eg_c$values > 1e-12))
    if (k_use >= 1L) {
      F_cor <- sweep(eg_c$vectors[, seq_len(k_use), drop = FALSE], 2L, sqrt(eg_c$values[seq_len(k_use)]), "*")
      F_cov <- sweep(F_cor, 1L, sd_c, "*")
      F_out[seq_len(k_use), pos_c] <- t(F_cov)
    }
    return(F_out)
  }

  R_p <- R_sym_c
  R_inv <- tryCatch(
    chol2inv(chol(R_sym_c + 1e-3 * diag(np))),
    error = function(e) {
      eg_c <- eigen(R_sym_c, symmetric = TRUE)
      R_p <<- eg_c$vectors %*% (pmax(eg_c$values, 1e-4) * t(eg_c$vectors))
      dr_p <- sqrt(pmax(diag(R_p), .Machine$double.eps))
      R_p <<- R_p / outer(dr_p, dr_p)
      R_p[off_c] <<- pmin(0.98, pmax(-0.98, R_p[off_c]))
      diag(R_p) <<- 1.0
      eg_p <- eigen(0.5 * (R_p + t(R_p)), symmetric = TRUE)
      eg_p$vectors %*% ((1.0 / pmax(eg_p$values, 1e-4)) * t(eg_p$vectors))
    }
  )
  unit_of_pos <- if (!is.null(nC) && nC >= 2L && q == nC * n_lags) {
    ((pos_c - 1L) %% nC) + 1L
  } else {
    seq_len(np)
  }
  u_blocks <- split(seq_len(np), unit_of_pos)

  R_work <- R_p
  for (idx_b in u_blocks) {
    if (length(idx_b) == 1L) {
      cond_b <- 1.0 / max(R_inv[idx_b, idx_b], 1.0)
      R_work[idx_b, idx_b] <- max(0.0, 1.0 - 0.98 * cond_b)
    } else {
      Prec_bb <- 0.5 * (R_inv[idx_b, idx_b, drop = FALSE] + t(R_inv[idx_b, idx_b, drop = FALSE]))
      Cond_bb <- solve(Prec_bb + 1e-6 * diag(length(idx_b)))
      R_work[idx_b, idx_b] <- R_p[idx_b, idx_b, drop = FALSE] - 0.98 * Cond_bb
    }
  }

  F_cor <- matrix(0.0, np, 0L)
  k_use <- 0L
  for (it in 1:3) {
    if (np <= 512L) {
      eg_s <- eigen(0.5 * (R_work + t(R_work)), symmetric = TRUE)
      k_use <- min(k_rank, sum(eg_s$values > 1e-6))
      if (k_use < 1L) break
      F_cor <- sweep(eg_s$vectors[, seq_len(k_use), drop = FALSE], 2L, sqrt(eg_s$values[seq_len(k_use)]), "*")
    } else {
      R_shift <- 0.5 * (R_work + t(R_work)) + 1.0 * diag(np)
      r_sub <- min(np, k_rank + 12L)
      Omega <- if (!is.null(probe_cache) && nrow(probe_cache) >= np && ncol(probe_cache) >= r_sub) {
        probe_cache[seq_len(np), seq_len(r_sub), drop = FALSE]
      } else {
        qr.Q(qr(outer(seq_len(np), seq_len(r_sub), function(i, j) sin((i * 1.6180339887 + j) * j))))
      }
      Q_it <- qr.Q(qr(R_shift %*% Omega))
      Q_it <- qr.Q(qr(R_shift %*% Q_it))
      Q_it <- qr.Q(qr(R_shift %*% Q_it))
      B_sm <- crossprod(Q_it, R_shift %*% Q_it)
      eg_s <- eigen(0.5 * (B_sm + t(B_sm)), symmetric = TRUE)
      vals_un <- eg_s$values - 1.0
      k_use <- min(k_rank, sum(vals_un > 1e-6))
      if (k_use < 1L) break
      V_big <- Q_it %*% eg_s$vectors[, seq_len(k_use), drop = FALSE]
      F_cor <- sweep(V_big, 2L, sqrt(vals_un[seq_len(k_use)]), "*")
    }
    rn <- sqrt(rowSums(F_cor^2))
    cap <- rn > 0.995
    if (any(cap)) F_cor[cap, ] <- F_cor[cap, , drop = FALSE] * (0.995 / rn[cap])
    if (it < 3L) {
      for (idx_b in u_blocks) {
        R_work[idx_b, idx_b] <- tcrossprod(F_cor[idx_b, , drop = FALSE])
      }
    }
  }

  if (k_use >= 1L) {
    F_cov <- sweep(F_cor, 1L, sd_c, "*")
    F_out[seq_len(k_use), pos_c] <- t(F_cov)
  }
  F_out
}

compute_shared_loading_from_cor <- function(R_cor, unit_ids) {
  P <- length(unit_ids)
  out <- stats::setNames(numeric(P), unit_ids)
  if (is.null(R_cor) || !is.matrix(R_cor) || nrow(R_cor) != P || P < 2L) return(out)
  off_mask <- row(R_cor) != col(R_cor)
  max_off <- max(abs(R_cor[off_mask]))
  if (!is.finite(max_off) || max_off <= 1e-6) return(out)
  if (P == 2L) {
    r12 <- R_cor[1L, 2L]
    out[] <- if (r12 >= 0) c(sqrt(r12), sqrt(r12)) else c(sqrt(-r12), -sqrt(-r12))
    return(out)
  }
  R_pa <- R_cor
  diag(R_pa) <- apply(abs(R_cor - diag(P)), 1L, max)
  for (it_pa in seq_len(4L)) {
    eg_pa <- eigen(0.5 * (R_pa + t(R_pa)), symmetric = TRUE)
    lv_pa <- eg_pa$vectors[, 1L] * sqrt(max(0.0, eg_pa$values[1L]))
    if (sum(lv_pa) < 0) lv_pa <- -lv_pa
    lv_pa <- pmin(0.995, pmax(-0.995, lv_pa))
    diag(R_pa) <- lv_pa^2
  }
  out[] <- lv_pa
  out
}

# Resolve (stat_vec, dyn_vec, rob_pair) for pair (g, h) from exact or lowrank_voxel shared_error_model
.sem_pair_params <- function(sem, g, h) {
  if (identical(sem$mode, "lowrank_voxel")) {
    if (g == h) {
      return(list(
        gam_s = as.numeric(sem$Gamma_stat_diag[, g]),
        b_d = as.numeric(sem$B_dyn_diag[, g])
      ))
    }
    L_p <- sem$L_poet
    if (is.null(L_p) || ncol(L_p) == 0L) {
      return(list(
        gam_s = numeric(sem$L_max + 1L),
        b_d = numeric(sem$M_dyn)
      ))
    }
    sd_g <- sqrt(max(sem$Gamma_stat_diag[1L, g], 1e-10))
    sd_h <- sqrt(max(sem$Gamma_stat_diag[1L, h], 1e-10))
    lg <- L_p[g, ] * sd_g
    lh <- L_p[h, ] * sd_h
    L1 <- sem$L_max + 1L
    gam_s <- numeric(L1)
    for (idx_l in seq_len(L1)) {
      gam_s[idx_l] <- as.numeric(crossprod(lg, matrix(sem$Gamma_stat_fac[idx_l, , ], ncol(L_p), ncol(L_p)) %*% lh))
    }
    b_d <- numeric(sem$M_dyn)
    for (j in seq_len(sem$M_dyn)) {
      b_d[j] <- as.numeric(crossprod(lg, matrix(sem$B_dyn_fac[j, , ], ncol(L_p), ncol(L_p)) %*% lh))
    }
    return(list(gam_s = gam_s, b_d = b_d))
  }
  list(
    gam_s = as.numeric(sem$Gamma_stat_gated[, g, h]),
    b_d = as.numeric(sem$B_dyn_gated[, g, h])
  )
}

medusa_eval_sem_bilinear <- function(U_L, U_R, V_L, V_R, rs_sem, sem, g, h, stat_floor = FALSE) {
  sc_gh <- if (!is.null(rs_sem$run_scale)) (rs_sem$run_scale[g] * rs_sem$run_scale[h]) else 1.0
  if (stat_floor) {
    sig_fl <- if (!is.null(sem$sig0_rob)) sem$sig0_rob else sem$sig0_sd
    return(sc_gh * (sig_fl[g] * sig_fl[h]) * colSums(U_L * U_R))
  }
  if (is.null(sem$Gamma_stat_gated) && is.null(sem$Gamma_stat_diag) && !is.null(sem$Theta_gated)) {
    m <- rs_sem$m
    K_r <- sem$K_reg; B_t <- sem$B_tot
    th_vec <- sem$Theta_gated[, g, h]
    f_mat <- cbind(1.0, rs_sem$e_task, rs_sem$E_cov_env)
    v0 <- as.numeric(f_mat %*% th_vec[1:K_r])
    f_mat1 <- 0.5 * (f_mat[2:m, , drop = FALSE] + f_mat[1:(m - 1L), , drop = FALSE])
    v1 <- as.numeric(f_mat1 %*% th_vec[(K_r + 1L):B_t])
    return(sc_gh * (colSums(U_L * (v0 * U_R)) +
                      colSums(U_L[2:m, , drop = FALSE] * (v1 * U_R[1:(m - 1L), , drop = FALSE])) +
                      colSums(U_L[1:(m - 1L), , drop = FALSE] * (v1 * U_R[2:m, , drop = FALSE]))))
  }

  pars <- .sem_pair_params(sem, g, h)
  gam_s <- pars$gam_s
  b_d <- pars$b_d
  if (!is.null(sem$kappa_eff_v1)) {
    m <- rs_sem$m
    f_dyn <- cbind(rs_sem$e_task, rs_sem$E_cov_env)
    v0 <- gam_s[1L] + as.numeric(f_dyn %*% b_d)
    f_d1 <- 0.5 * (f_dyn[2:m, , drop = FALSE] + f_dyn[1:(m - 1L), , drop = FALSE])
    v1 <- sem$kappa_eff_v1[g, h] * as.numeric(f_d1 %*% b_d)
    return(sc_gh * (colSums(U_L * (v0 * U_R)) +
                      colSums(U_L[2:m, , drop = FALSE] * (v1 * U_R[1:(m - 1L), , drop = FALSE])) +
                      colSums(U_L[1:(m - 1L), , drop = FALSE] * (v1 * U_R[2:m, , drop = FALSE]))))
  }

  ans <- gam_s[1L] * colSums(U_L * U_R)
  L_m <- sem$L_max
  for (d in seq_len(L_m)) {
    gd <- gam_s[d + 1L]
    if (abs(gd) > 1e-14) {
      pd <- rs_sem$lag_pairs[[d]]
      if (nrow(pd) > 0L) {
        ans <- ans + gd * colSums(U_L[pd[, 1L], , drop = FALSE] * U_R[pd[, 2L], , drop = FALSE])
      }
    }
  }
  if (any(abs(b_d) > 1e-14)) {
    f_dyn <- cbind(rs_sem$e_task, rs_sem$E_cov_env)
    gm_vec <- rep(sem$gamma_obs, length(b_d))
    if (!is.null(sem$j_rwls) && length(sem$j_rwls) == 1L) gm_vec[sem$j_rwls] <- 1.0
    v_inn <- as.numeric(f_dyn %*% ((1.0 - gm_vec) * b_d))
    v_obs <- as.numeric(f_dyn %*% (gm_vec * b_d))
    if (any(abs(v_inn) > 1e-14)) ans <- ans + colSums(U_L * (v_inn * U_R))
    if (any(abs(v_obs) > 1e-14) && !is.null(V_L) && !is.null(V_R)) ans <- ans + colSums(V_L * (v_obs * V_R))
  }
  sc_gh * ans
}

medusa_eval_sem_bilinear_mat <- function(U_L, U_R, V_L, V_R, rs_sem, sem, g, h) {
  sc_gh <- if (!is.null(rs_sem$run_scale)) (rs_sem$run_scale[g] * rs_sem$run_scale[h]) else 1.0
  if (is.null(sem$Gamma_stat_gated) && is.null(sem$Gamma_stat_diag) && !is.null(sem$Theta_gated)) {
    m <- rs_sem$m
    K_r <- sem$K_reg; B_t <- sem$B_tot
    th_vec <- sem$Theta_gated[, g, h]
    f_mat <- cbind(1.0, rs_sem$e_task, rs_sem$E_cov_env)
    v0 <- as.numeric(f_mat %*% th_vec[1:K_r])
    f_mat1 <- 0.5 * (f_mat[2:m, , drop = FALSE] + f_mat[1:(m - 1L), , drop = FALSE])
    v1 <- as.numeric(f_mat1 %*% th_vec[(K_r + 1L):B_t])
    return(sc_gh * (crossprod(U_L, v0 * U_R) +
                      crossprod(U_L[2:m, , drop = FALSE], v1 * U_R[1:(m - 1L), , drop = FALSE]) +
                      crossprod(U_L[1:(m - 1L), , drop = FALSE], v1 * U_R[2:m, , drop = FALSE])))
  }

  pars <- .sem_pair_params(sem, g, h)
  gam_s <- pars$gam_s
  b_d <- pars$b_d
  if (!is.null(sem$kappa_eff_v1)) {
    m <- rs_sem$m
    f_dyn <- cbind(rs_sem$e_task, rs_sem$E_cov_env)
    v0 <- gam_s[1L] + as.numeric(f_dyn %*% b_d)
    f_d1 <- 0.5 * (f_dyn[2:m, , drop = FALSE] + f_dyn[1:(m - 1L), , drop = FALSE])
    v1 <- sem$kappa_eff_v1[g, h] * as.numeric(f_d1 %*% b_d)
    return(sc_gh * (crossprod(U_L, v0 * U_R) +
                      crossprod(U_L[2:m, , drop = FALSE], v1 * U_R[1:(m - 1L), , drop = FALSE]) +
                      crossprod(U_L[1:(m - 1L), , drop = FALSE], v1 * U_R[2:m, , drop = FALSE])))
  }

  ans <- gam_s[1L] * crossprod(U_L, U_R)
  L_m <- sem$L_max
  for (d in seq_len(L_m)) {
    gd <- gam_s[d + 1L]
    if (abs(gd) > 1e-14) {
      pd <- rs_sem$lag_pairs[[d]]
      if (nrow(pd) > 0L) {
        ans <- ans + gd * crossprod(U_L[pd[, 1L], , drop = FALSE], U_R[pd[, 2L], , drop = FALSE])
      }
    }
  }
  if (any(abs(b_d) > 1e-14)) {
    f_dyn <- cbind(rs_sem$e_task, rs_sem$E_cov_env)
    gm_vec <- rep(sem$gamma_obs, length(b_d))
    if (!is.null(sem$j_rwls) && length(sem$j_rwls) == 1L) gm_vec[sem$j_rwls] <- 1.0
    v_inn <- as.numeric(f_dyn %*% ((1.0 - gm_vec) * b_d))
    v_obs <- as.numeric(f_dyn %*% (gm_vec * b_d))
    if (any(abs(v_inn) > 1e-14)) ans <- ans + crossprod(U_L, v_inn * U_R)
    if (any(abs(v_obs) > 1e-14) && !is.null(V_L) && !is.null(V_R)) ans <- ans + crossprod(V_L, v_obs * V_R)
  }
  sc_gh * ans
}

# Resolve variance_confounds input (NULL, medusa_confounds, medusa_nifti_ts, data.frame, or run list)
.resolve_variance_confounds_input <- function(
  variance_confounds,
  run_ids,
  run_col = "run"
) {
  if (is.null(variance_confounds)) return(NULL)
  run_ids <- as.character(run_ids)
  n_runs <- length(run_ids)

  coerce_vc_mat <- function(el, r_lbl) {
    if (is.null(el)) return(NULL)
    if (is.data.frame(el)) {
      if (!is.null(run_col) && nzchar(run_col) && (run_col %in% names(el))) {
        el <- el[, setdiff(names(el), run_col), drop = FALSE]
      }
      if (ncol(el) == 0L) {
        return(matrix(numeric(0L), nrow = nrow(el), ncol = 0L))
      }
      num_ok <- vapply(el, function(col) is.numeric(col) || is.logical(col), logical(1L))
      if (any(!num_ok)) {
        stop(
          "Non-numeric column(s) in variance_confounds[['", r_lbl, "']]: ",
          paste(names(el)[!num_ok], collapse = ", "),
          call. = FALSE
        )
      }
      m_out <- as.matrix(el)
      storage.mode(m_out) <- "double"
      return(m_out)
    }
    if (is.numeric(el) || is.logical(el)) {
      m_out <- as.matrix(el)
      storage.mode(m_out) <- "double"
      return(m_out)
    }
    stop(
      "variance_confounds[['", r_lbl, "']] must be a numeric vector, matrix, or data.frame.",
      call. = FALSE
    )
  }

  if (inherits(variance_confounds, "medusa_nifti_ts")) {
    if (!is.null(variance_confounds$confounds)) {
      variance_confounds <- variance_confounds$confounds
    } else if (!is.null(variance_confounds$variance_confounds)) {
      variance_confounds <- variance_confounds$variance_confounds
    } else {
      stop(
        "'variance_confounds' is a 'medusa_nifti_ts' object without '$confounds'; ",
        "pass a 'medusa_confounds' object from read_fmriprep_confounds() or read_fmriprep_subject().",
        call. = FALSE
      )
    }
  }

  if (inherits(variance_confounds, "medusa_confounds")) {
    vc_src <- variance_confounds$variance_confounds
    if (is.null(vc_src) || !is.list(vc_src)) {
      # Fallback reconstruction for legacy serialized medusa_confounds objects
      fd_src <- variance_confounds$framewise_displacement
      raw_nuis <- attr(variance_confounds, "raw_nuisance", exact = TRUE)
      nuis_src <- if (is.list(raw_nuis)) raw_nuis else variance_confounds$nuisance
      r_names_src <- names(variance_confounds$valid_rows)
      if (is.null(r_names_src)) r_names_src <- names(fd_src)
      n_src <- length(variance_confounds$valid_rows)
      vc_src <- stats::setNames(vector("list", n_src), r_names_src)
      for (ii in seq_len(n_src)) {
        n_v <- length(variance_confounds$valid_rows[[ii]])
        cols_ii <- list()
        if (is.list(fd_src) && length(fd_src) >= ii && !is.null(fd_src[[ii]])) {
          fd_v <- as.numeric(fd_src[[ii]])
          if (length(fd_v) == n_v && any(is.finite(fd_v))) {
            if (!is.finite(fd_v[1L])) fd_v[1L] <- 0.0
            bad_f <- !is.finite(fd_v)
            if (any(bad_f)) fd_v[bad_f] <- stats::median(fd_v[!bad_f])
            fd_v <- pmax(0.0, fd_v)
            if (n_v >= 2L && is.finite(stats::sd(fd_v)) && stats::sd(fd_v) > 1e-10) {
              cols_ii[["framewise_displacement"]] <- fd_v
            }
          }
        }
        if (is.list(nuis_src) && length(nuis_src) >= ii && is.matrix(nuis_src[[ii]])) {
          zm <- nuis_src[[ii]]
          for (cn in intersect(c("std_dvars", "dvars", "rvt", "respiratory_volume_per_time"), colnames(zm))) {
            cv <- as.numeric(zm[, cn])
            if (any(is.finite(cv))) {
              cv[!is.finite(cv)] <- stats::median(cv[is.finite(cv)])
              if (cn %in% c("std_dvars", "dvars")) cv <- pmax(0.0, cv)
              if (n_v >= 2L && is.finite(stats::sd(cv)) && stats::sd(cv) > 1e-10) {
                cols_ii[[cn]] <- cv
              }
            }
          }
        }
        vc_src[[ii]] <- if (length(cols_ii) > 0L) {
          do.call(cbind, cols_ii)
        } else {
          matrix(numeric(0L), nrow = n_v, ncol = 0L)
        }
      }
    }
    src_names <- names(vc_src)
    if (!is.null(src_names) && all(run_ids %in% src_names)) {
      vc_aligned <- vc_src[run_ids]
    } else if (length(vc_src) == n_runs && (is.null(src_names) || !any(nzchar(src_names)) || n_runs == 1L)) {
      vc_aligned <- vc_src
    } else {
      stop(
        sprintf(
          "Run mismatch in 'variance_confounds': medusa_confounds has run(s) [%s], expected [%s].",
          paste(if (is.null(src_names)) seq_along(vc_src) else src_names, collapse = ", "),
          paste(run_ids, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    return(stats::setNames(
      lapply(seq_len(n_runs), function(i) coerce_vc_mat(vc_aligned[[i]], run_ids[i])),
      run_ids
    ))
  }

  if (is.data.frame(variance_confounds)) {
    if (!is.null(run_col) && nzchar(run_col) && (run_col %in% names(variance_confounds))) {
      r_chr <- as.character(variance_confounds[[run_col]])
      if (all(run_ids %in% unique(r_chr))) {
        return(stats::setNames(
          lapply(run_ids, function(rn) {
            sub_df <- variance_confounds[!is.na(r_chr) & r_chr == rn, , drop = FALSE]
            coerce_vc_mat(sub_df, rn)
          }),
          run_ids
        ))
      }
    }
    if (n_runs == 1L) {
      return(stats::setNames(list(coerce_vc_mat(variance_confounds, run_ids[1L])), run_ids))
    }
    stop(
      "Multi-run 'variance_confounds' data.frame must contain a '", run_col,
      "' column matching fitted run_ids (or pass a 'medusa_confounds' object or per-run list).",
      call. = FALSE
    )
  }

  if (is.list(variance_confounds)) {
    v_names <- names(variance_confounds)
    if (!is.null(v_names) && any(nzchar(v_names))) {
      if (!all(run_ids %in% v_names)) {
        stop(
          "Named 'variance_confounds' list is missing run(s): ",
          paste(setdiff(run_ids, v_names), collapse = ", "),
          call. = FALSE
        )
      }
      vc_sub <- variance_confounds[run_ids]
    } else if (length(variance_confounds) == n_runs) {
      vc_sub <- variance_confounds
    } else {
      stop(
        sprintf(
          "'variance_confounds' list length (%d) must match number of runs (%d).",
          length(variance_confounds), n_runs
        ),
        call. = FALSE
      )
    }
    return(stats::setNames(
      lapply(seq_len(n_runs), function(i) coerce_vc_mat(vc_sub[[i]], run_ids[i])),
      run_ids
    ))
  }

  if (is.numeric(variance_confounds) || is.logical(variance_confounds)) {
    if (n_runs == 1L) {
      return(stats::setNames(list(coerce_vc_mat(variance_confounds, run_ids[1L])), run_ids))
    }
    stop(
      "Multi-run 'variance_confounds' must be a 'medusa_confounds' object or a per-run list of matrices/vectors.",
      call. = FALSE
    )
  }

  stop(
    "Unsupported 'variance_confounds' input; supply NULL, a 'medusa_confounds' object, ",
    "a per-run list, or (for single-run fits) a numeric matrix/data.frame/vector.",
    call. = FALSE
  )
}

# Core v3 POET + Multi-Lag Toeplitz & General-W_s Covariance-Regression Engine
medusa_fit_shared_error_covreg <- function(
  deconvolution,
  events,
  target_idx_list,
  dur_target_list,
  cov_list = NULL,
  lags = 0,
  readout_width = NULL,
  restore_scale = FALSE,
  source_model = "impulse",
  aggregate_by_parcel = FALSE,
  unit_ids,
  parcel_vox_list = NULL,
  parcel_weights_list = NULL,
  run_stratum_ops,
  variance_confounds = NULL,
  bold = NULL,
  return_panel = TRUE,
  return_epochs = FALSE,
  max_exact_units = 512L,
  embed_rank = 24L,
  k_max_poet = 12L
) {
  fit <- deconvolution$hrf_fit
  settings <- deconvolution$settings
  tr <- fit$settings$tr
  n_runs <- length(fit$run_ids)
  lags <- as.numeric(lags)
  n_lags <- length(lags)
  P <- length(unit_ids)
  n_vox_tot <- length(fit$voxel_ids)
  L_max <- min(12L, max(8L, as.integer(round(8.0 / tr))))
  L1 <- L_max + 1L
  center_obs <- !identical(settings$scaling, "none")

  diff1_cols <- function(a) {
    if (ncol(a) == 0L || nrow(a) < 2L) return(a)
    rbind(a[2:nrow(a), , drop = FALSE], 0.0) - rbind(0.0, a[1:(nrow(a) - 1L), , drop = FALSE])
  }
  std_env <- function(v) {
    s <- stats::sd(v)
    if (!is.finite(s) || s <= 1e-10) return(rep(0.0, length(v)))
    as.numeric((v - mean(v)) / s)
  }

  q_cov <- if (!is.null(cov_list) && length(cov_list) >= 1L && is.matrix(cov_list[[1L]])) ncol(cov_list[[1L]]) else 0L
  cov_names <- if (q_cov > 0L) colnames(cov_list[[1L]]) else character(0L)
  exact_mode <- (P <= as.integer(max_exact_units))
  vc_resolved <- .resolve_variance_confounds_input(variance_confounds, fit$run_ids)

  run_store <- vector("list", n_runs)
  any_saturated <- FALSE
  df_u_tot <- 0.0

  for (r in seq_len(n_runs)) {
    r_op <- run_stratum_ops[[r]]
    if (is.null(r_op)) next
    n <- fit$n_time[r]
    rows <- r_op$rows
    m <- length(rows)
    ev <- events[[r]]
    offset <- if (length(fit$settings$time_offset) >= r) fit$settings$time_offset[r] else 0.0
    t_idx <- target_idx_list[[r]]
    m_all <- nrow(ev)
    m_target <- length(t_idx)
    duration <- if (is.null(ev$duration)) rep(0.0, m_all) else as.numeric(ev$duration)
    if (!is.null(dur_target_list) && length(dur_target_list) >= r && !is.null(dur_target_list[[r]])) {
      duration[t_idx] <- as.numeric(dur_target_list[[r]])
    }
    dur_target <- duration[t_idx]
    w_src <- if (is.null(readout_width)) 1.5 * tr else readout_width[1L]

    fw_r <- r_op$fw_r
    nuis_r <- if (!is.null(fit$specification$nuisance) &&
                  length(fit$specification$nuisance) >= r &&
                  !is.null(fit$specification$nuisance[[r]])) {
      fit$specification$nuisance[[r]]
    } else {
      matrix(numeric(), n, 0L)
    }
    z_r_base <- cbind(intercept = 1, trend = seq(-1, 1, length.out = n), nuis_r)
    mask_bool <- seq_len(n) %in% rows

    sc_df <- deconvolution$scaling[[r]]
    valid_vox <- sc_df$status == "ok"
    div_r <- sc_df$scale
    div_r[!is.finite(div_r) | div_r <= 0] <- 1.0
    scale_mult_vox <- if (isTRUE(restore_scale)) rep(1.0, n_vox_tot) else (1.0 / div_r)

    # Resolve observed BOLD or Step-2 prewhitened residual sketch without requiring raw `bold`
    y_obs_r <- NULL
    if (!is.null(bold) && length(bold) >= r && is.matrix(bold[[r]])) {
      y_obs_r <- if (nrow(bold[[r]]) == n) bold[[r]][rows, , drop = FALSE] else bold[[r]]
    } else if (!is.null(deconvolution$bold_obs) && length(deconvolution$bold_obs) >= r &&
               is.matrix(deconvolution$bold_obs[[r]]) && nrow(deconvolution$bold_obs[[r]]) == m) {
      y_obs_r <- deconvolution$bold_obs[[r]]
    }
    sketch_r <- if (is.null(y_obs_r) && !is.null(deconvolution$whitened_sketch) &&
                    length(deconvolution$whitened_sketch) >= r) {
      deconvolution$whitened_sketch[[r]]
    } else NULL

    # Build unit-level observed matrix Y_unit when y_obs_r is available and exact_mode (or moderate V)
    Y_unit <- NULL
    if (!is.null(y_obs_r)) {
      if (isTRUE(aggregate_by_parcel)) {
        Y_unit <- matrix(0.0, m, P, dimnames = list(NULL, unit_ids))
        for (g in seq_len(P)) {
          vg <- parcel_vox_list[[g]]
          ok_g <- valid_vox[vg] & (colSums(!is.finite(y_obs_r[, vg, drop = FALSE])) == 0L)
          if (!any(ok_g)) next
          v_ok <- vg[ok_g]
          ug <- parcel_weights_list[[g]][ok_g]
          ug <- ug / sum(ug)
          y_sub <- y_obs_r[, v_ok, drop = FALSE]
          if (!isTRUE(restore_scale)) y_sub <- sweep(y_sub, 2L, div_r[v_ok], "/")
          Y_unit[, g] <- as.numeric(y_sub %*% ug)
        }
      } else if (exact_mode) {
        Y_unit <- y_obs_r
        bad_v <- !valid_vox | (colSums(!is.finite(Y_unit)) > 0L)
        if (any(bad_v)) Y_unit[, bad_v] <- 0.0
        if (!isTRUE(restore_scale)) Y_unit <- sweep(Y_unit, 2L, div_r, "/")
        colnames(Y_unit) <- unit_ids
      }
    }

    u_skeys <- r_op$u_skeys
    S_cnt <- length(u_skeys)
    s_idx_by_unit <- r_op$s_idx_by_unit
    stratum_specs <- r_op$stratum_specs
    rho_unit <- vapply(seq_len(P), function(g) {
      si <- s_idx_by_unit[g]
      if (is.na(si) || si < 1L) return(0.0)
      as.numeric(stratum_specs[[si]]$ar[1L])
    }, numeric(1L))

    X_cov_c <- if (q_cov > 0L && nrow(cov_list[[r]]) == m_target) {
      xc <- scale(cov_list[[r]], center = TRUE, scale = FALSE)
      xc[!is.finite(xc)] <- 0.0
      matrix(as.numeric(xc), nrow = m_target, ncol = q_cov, dimnames = list(NULL, cov_names))
    } else {
      matrix(0.0, nrow = m_target, ncol = 0L)
    }

    M_strata <- vector("list", S_cnt)
    G_inn_strata <- vector("list", S_cnt)
    H_strata <- vector("list", S_cnt)
    U_strata <- vector("list", S_cnt)
    V_obs_strata <- vector("list", S_cnt)
    ok_strata <- vector("list", S_cnt)
    sat_strata <- logical(S_cnt)
    df_strata <- numeric(S_cnt)
    c_strata <- numeric(S_cnt)
    kpeak_strata <- numeric(S_cnt)
    h_task_accum <- numeric(m)
    x_env_accum <- matrix(0.0, m, ncol(X_cov_c))

    E_unit <- if (exact_mode) matrix(0.0, m, P, dimnames = list(NULL, unit_ids)) else NULL
    E_rob_unit <- if (exact_mode) matrix(0.0, m, P, dimnames = list(NULL, unit_ids)) else NULL
    parcel_idio_s2 <- numeric(P)

    # Low-rank voxel accumulators when P > max_exact_units
    s2_raw_u <- if (!exact_mode) numeric(P) else NULL
    s2_rob_u <- if (!exact_mode) numeric(P) else NULL
    s_blk_u <- if (!exact_mode) numeric(P) else NULL
    df_unit_v <- if (!exact_mode) rep(max(1.0, m - 4.0), P) else NULL
    E_vox_full <- if (!exact_mode && !is.null(y_obs_r)) matrix(0.0, m, P) else NULL
    E_rob_vox_full <- if (!exact_mode && !is.null(y_obs_r)) matrix(0.0, m, P) else NULL

    fw_eff <- if (!is.null(fw_r) && length(fw_r) == m) {
      fw_num <- as.numeric(fw_r)
      fw_out <- ifelse(fw_num < 0.03, fw_num, 1.0)
      if (!is.null(attr(fw_r, "weight_space"))) attr(fw_out, "weight_space") <- attr(fw_r, "weight_space")
      fw_out
    } else NULL
    has_rwls_w <- !is.null(fw_eff) && any(fw_eff < 1.0 - 1e-6)
    rwls_w_sp <- if (has_rwls_w && identical(attr(fw_eff, "weight_space"), "innovation")) "innovation" else "woodbury"
    e_rwls_act <- NULL
    if (has_rwls_w && rwls_w_sp == "woodbury") {
      w_cl <- pmin(pmax(as.numeric(fw_eff), 1e-4), 1.0)
      act_rw <- which(w_cl < 1.0 - 1e-6)
      if (length(act_rw) > 0L) {
        e_rwls_act <- matrix(0.0, m, length(act_rw))
        e_rwls_act[cbind(act_rw, seq_along(act_rw))] <- sqrt(pmax(1.0 / w_cl[act_rw] - 1.0, 0.0))
      }
    }
    sup_r <- if (!is.null(fw_r) && length(fw_r) == m) which(as.numeric(fw_r) < 0.85) else integer(0L)
    w_sup2_r <- if (length(sup_r)) pmin(1.0, pmax(0.0, as.numeric(fw_r[sup_r])^2)) else numeric(0L)

    for (si in seq_len(S_cnt)) {
      sp <- stratum_specs[[si]]
      cols_s <- which(s_idx_by_unit == si)
      w_share <- length(cols_s) / max(1L, P)
      k_vec <- as.numeric(sp$kernel[, 1L])
      ar_op <- sp$ar
      kpeak_strata[si] <- which.max(k_vec)[1L] * tr

      k_pos <- matrix(pmax(k_vec, 0.0), ncol = 1L)
      if (max(k_pos) > 0) k_pos <- k_pos / max(k_pos)
      G_fb_pos <- hrf_event_design_cpp(
        n, tr, offset, ev$onset[t_idx], dur_target, seq_len(m_target), m_target, k_pos
      )[[1L]][rows, , drop = FALSE]
      h_task_accum <- h_task_accum + w_share * rowSums(G_fb_pos)
      if (ncol(X_cov_c) > 0L) {
        x_env_accum <- x_env_accum + w_share * (G_fb_pos %*% X_cov_c)
      }

      # 4A: Duration-aware & df_res-guarded sustained / multi-lag neural protection in M_s
      ev_eff <- ev
      dur_proj <- duration
      if (source_model == "window" && all(dur_proj[t_idx] <= 0.0)) dur_proj[t_idx] <- w_src
      ev_eff$duration <- dur_proj

      md_s <- hrf_build_noise_mean_design(n, ev_eff, tr, offset, k_vec, z_r_base, mask_bool)
      x_gls <- md_s$design[rows, , drop = FALSE]
      is_extra_only <- FALSE
      x_tr_in <- if (!isTRUE(md_s$used_trial_design) && !is.null(md_s$trial_design) && ncol(md_s$trial_design) > 0L) {
        td_full <- md_s$trial_design
        cbind(td_full, diff1_cols(td_full))[rows, , drop = FALSE]
      } else if (isTRUE(md_s$used_trial_design) && m_target >= 1L) {
        extra_parts <- list()
        if (max(lags) >= 1.5 * tr || min(lags) <= -1.5 * tr) {
          lag_shift_ms <- if (max(lags) >= 1.5 * tr) max(lags) else min(lags)
          G_late <- hrf_event_design_cpp(
            n, tr, offset, ev$onset[t_idx] + lag_shift_ms, dur_proj[t_idx],
            seq_len(m_target), m_target, matrix(k_vec, ncol = 1L)
          )[[1L]][rows, , drop = FALSE]
          extra_parts[[length(extra_parts) + 1L]] <- G_late
        }
        if (any(dur_proj[t_idx] > 0.0) && !is.null(md_s$trial_design) && ncol(md_s$trial_design) >= max(t_idx)) {
          td_targ <- md_s$trial_design[, t_idx, drop = FALSE]
          extra_parts[[length(extra_parts) + 1L]] <- diff1_cols(td_targ)[rows, , drop = FALSE]
        }
        if (length(extra_parts) > 0L) {
          is_extra_only <- TRUE
          do.call(cbind, extra_parts)
        } else NULL
      } else NULL

      wy_in <- if (!is.null(Y_unit)) {
        Y_unit[, cols_s, drop = FALSE]
      } else if (!exact_mode && !is.null(y_obs_r) && length(cols_s) > 0L) {
        ys <- y_obs_r[, cols_s, drop = FALSE]
        bad_c <- !valid_vox[cols_s] | (colSums(!is.finite(ys)) > 0L)
        if (any(bad_c)) ys[, bad_c] <- 0.0
        if (!isTRUE(restore_scale)) ys <- sweep(ys, 2L, div_r[cols_s], "/")
        ys
      } else NULL

      wb <- hrf_apply_whitening(
        list(wx = x_gls, wy = wy_in, wt = x_tr_in, we = e_rwls_act),
        rows = rows, ar_vec = ar_op, unit_marginal = TRUE
      )
      apply_rwls_whiten <- function(mat_in) mat_in
      if (has_rwls_w) {
        if (rwls_w_sp == "innovation") {
          sw_rw <- sqrt(pmax(as.numeric(fw_eff), 1e-8))
          apply_rwls_whiten <- function(mat_in) if (is.null(mat_in)) NULL else sweep(mat_in, 1L, sw_rw, "*")
        } else if (!is.null(wb$we) && ncol(wb$we) > 0L) {
          sv_we <- hrf_safe_svd(wb$we, nv = 0L)
          pos_we <- if (length(sv_we$d)) sv_we$d > 1e-8 else logical(0L)
          if (any(pos_we)) {
            u_we <- sv_we$u[, pos_we, drop = FALSE]
            fac_we <- 1.0 / sqrt(1.0 + sv_we$d[pos_we]^2) - 1.0
            apply_rwls_whiten <- function(mat_in) {
              if (is.null(mat_in)) return(NULL)
              mat_in + u_we %*% sweep(crossprod(u_we, mat_in), 1L, fac_we, "*")
            }
          }
        }
        wb$wx <- apply_rwls_whiten(wb$wx)
        if (!is.null(wb$wy)) wb$wy <- apply_rwls_whiten(wb$wy)
        if (!is.null(wb$wt)) wb$wt <- apply_rwls_whiten(wb$wt)
      }

      nrm <- sqrt(colSums(wb$wx^2)); nrm[nrm == 0] <- 1.0
      sv <- hrf_safe_svd(sweep(wb$wx, 2L, nrm, "/"), nv = 0L)
      rk <- if (length(sv$d) && sv$d[1L] > 0) sum(sv$d > 1e-8 * sv$d[1L]) else 0L
      df_res <- m - rk
      Q_tot <- if (rk > 0L) sv$u[, seq_len(rk), drop = FALSE] else matrix(0.0, m, 0L)
      is_sat <- !isTRUE(md_s$used_trial_design)
      df_floor_ms <- if (is_extra_only) max(120L, as.integer(floor(0.60 * m))) else 24L
      if (!is.null(wb$wt) && df_res > df_floor_ms) {
        wt <- if (rk > 0L) wb$wt - Q_tot %*% crossprod(Q_tot, wb$wt) else wb$wt
        nt <- sqrt(colSums(wt^2)); pos_t <- nt > 1e-8
        if (any(pos_t)) {
          sv_t <- hrf_safe_svd(sweep(wt[, pos_t, drop = FALSE], 2L, nt[pos_t], "/"), nv = 0L)
          rk_max <- if (is_extra_only) min(m_target, max(0L, df_res - df_floor_ms)) else max(0L, df_res - 24L)
          sv_tol_t <- if (is_extra_only) 0.12 else 0.06
          n_want <- if (length(sv_t$d) && sv_t$d[1L] > 0) sum(sv_t$d > sv_tol_t * sv_t$d[1L]) else 0L
          rk_t <- min(rk_max, n_want)
          if (rk_t > 0L) Q_tot <- cbind(Q_tot, sv_t$u[, seq_len(rk_t), drop = FALSE])
          if (!is_extra_only) is_sat <- (n_want > rk_max)
        }
      }
      sat_strata[si] <- is_sat
      if (is_sat) any_saturated <- TRUE

      M_s <- diag(m) - tcrossprod(Q_tot)
      df_s <- max(1.0, sum(diag(M_s)))
      df_strata[si] <- df_s
      M_strata[[si]] <- M_s
      G_inn_strata[[si]] <- if (has_rwls_w) t(apply_rwls_whiten(M_s)) else M_s

      W_s <- hrf_whitening_matrix(rows, ar_op, unit_marginal = TRUE)
      c_s <- sqrt(max(mean(rowSums(W_s^2)), 1e-12))
      W_tilde_0 <- W_s / c_s
      W_tilde <- if (has_rwls_w) apply_rwls_whiten(W_s) / c_s else W_tilde_0
      c_strata[si] <- c_s
      H_strata[[si]] <- M_s %*% W_tilde

      if (length(cols_s) > 0L) {
        if (!is.null(wb$wy)) {
          lev_s <- pmin(0.95, pmax(0.0, 1.0 - diag(M_s)))
          ok_lev <- (1.0 - lev_s) > 0.10
          E_s <- M_s %*% wb$wy
          E_rob_s <- E_s
          if (sum(ok_lev) >= 16L) {
            inv_sd <- 1.0 / sqrt(1.0 - lev_s[ok_lev])
            for (jj in seq_along(cols_s)) {
              zl <- E_s[ok_lev, jj] * inv_sd
              s0 <- stats::mad(zl, center = 0.0)
              if (is.finite(s0) && s0 > .Machine$double.eps) {
                E_rob_s[ok_lev, jj] <- (sign(zl) * pmin(abs(zl), 4.0 * s0) / sqrt(0.99944)) / inv_sd
                E_s[ok_lev, jj] <- (sign(zl) * pmin(abs(zl), 12.0 * s0)) / inv_sd
              }
            }
          }
          if (exact_mode) {
            E_unit[, cols_s] <- E_s
            E_rob_unit[, cols_s] <- E_rob_s
          } else {
            E_vox_full[, cols_s] <- E_s
            E_rob_vox_full[, cols_s] <- E_rob_s
            s2_raw_u[cols_s] <- colSums(E_s^2) / df_s
            s2_rob_u[cols_s] <- colSums(E_rob_s^2) / df_s
            df_unit_v[cols_s] <- df_s
            if (length(sup_r) > 0L) {
              ex_s <- matrix(pmax(0.0, E_s[sup_r, , drop = FALSE]^2 - E_rob_s[sup_r, , drop = FALSE]^2),
                             nrow = length(sup_r), ncol = length(cols_s))
              s_blk_u[cols_s] <- colSums((1.0 - w_sup2_r) * ex_s)
            }
          }
        } else if (!is.null(sketch_r)) {
          # Reconstruct from Step-2 whole-brain prewhitened residual sketch without raw `bold`
          if (isTRUE(aggregate_by_parcel)) {
            for (jj in seq_along(cols_s)) {
              g <- cols_s[jj]
              vg <- parcel_vox_list[[g]]
              ok_g <- valid_vox[vg] & is.finite(sketch_r$s2_raw[vg])
              if (!any(ok_g)) next
              v_ok <- vg[ok_g]
              ug <- parcel_weights_list[[g]][ok_g]
              ug <- ug / sum(ug)
              wt_v <- ug * scale_mult_vox[v_ok]
              b_g <- as.numeric(sketch_r$B_sk[, v_ok, drop = FALSE] %*% wt_v)
              br_g <- as.numeric(sketch_r$B_rob_sk[, v_ok, drop = FALSE] %*% wt_v)
              E_unit[, g] <- as.numeric(M_s %*% (sketch_r$U_sk %*% b_g))
              E_rob_unit[, g] <- as.numeric(M_s %*% (sketch_r$U_sk %*% br_g))
              s2_cap_v <- colSums(sketch_r$B_sk[, v_ok, drop = FALSE]^2) / pmax(sketch_r$df_vox[v_ok], 1.0)
              s2_unexp <- pmax(0.0, sketch_r$s2_raw[v_ok] - s2_cap_v)
              parcel_idio_s2[g] <- parcel_idio_s2[g] + sum((wt_v^2) * s2_unexp)
            }
          } else if (exact_mode) {
            wt_c <- scale_mult_vox[cols_s]
            E_unit[, cols_s] <- M_s %*% (sketch_r$U_sk %*% sweep(sketch_r$B_sk[, cols_s, drop = FALSE], 2L, wt_c, "*"))
            E_rob_unit[, cols_s] <- M_s %*% (sketch_r$U_sk %*% sweep(sketch_r$B_rob_sk[, cols_s, drop = FALSE], 2L, wt_c, "*"))
          } else {
            wt_c <- scale_mult_vox[cols_s]
            s2_raw_u[cols_s] <- ifelse(is.finite(sketch_r$s2_raw[cols_s]), sketch_r$s2_raw[cols_s] * (wt_c^2), 0.0)
            s2_rob_u[cols_s] <- ifelse(is.finite(sketch_r$s2_rob[cols_s]), sketch_r$s2_rob[cols_s] * (wt_c^2), 0.0)
            s_blk_u[cols_s] <- sketch_r$s_blocked[cols_s] * (wt_c^2)
            df_unit_v[cols_s] <- ifelse(is.finite(sketch_r$df_vox[cols_s]), sketch_r$df_vox[cols_s], df_s)
          }
        } else if (!is.null(deconvolution$noise_innovations) && length(deconvolution$noise_innovations) >= r &&
                   is.matrix(deconvolution$noise_innovations[[r]])) {
          # Legacy fallback when crossrun_decon came from an older serialized object without bold_obs
          Z_in <- deconvolution$noise_innovations[[r]]
          nsd_v <- ifelse(is.finite(sc_df$noise_sd) & sc_df$noise_sd > 0, sc_df$noise_sd, div_r) * scale_mult_vox
          if (isTRUE(aggregate_by_parcel)) {
            for (jj in seq_along(cols_s)) {
              g <- cols_s[jj]
              vg <- parcel_vox_list[[g]]
              ok_g <- valid_vox[vg]
              if (!any(ok_g)) next
              v_ok <- vg[ok_g]
              ug <- parcel_weights_list[[g]][ok_g]; ug <- ug / sum(ug)
              e_g <- as.numeric(Z_in[, v_ok, drop = FALSE] %*% (ug * nsd_v[v_ok] * sqrt(df_s)))
              E_unit[, g] <- e_g
              E_rob_unit[, g] <- e_g
            }
          } else if (exact_mode) {
            E_unit[, cols_s] <- sweep(Z_in[, cols_s, drop = FALSE], 2L, nsd_v[cols_s] * sqrt(df_s), "*")
            E_rob_unit[, cols_s] <- E_unit[, cols_s]
          }
        }
      }

      if (isTRUE(return_panel) || isTRUE(return_epochs)) {
        U_l <- vector("list", n_lags)
        V_l <- vector("list", n_lags)
        ok_tau_l <- sp$ok_tau_list
        B_demix_l <- sp$B_demix_list
        for (li in seq_len(n_lags)) {
          ok_t <- ok_tau_l[[li]]
          U_full <- matrix(0.0, m, m_target)
          V_full <- matrix(0.0, m, m_target)
          if (!is.null(ok_t) && any(ok_t) && !is.null(B_demix_l[[li]])) {
            U_sub_0 <- compute_ar1_operator_factor(
              B_demix_l[[li]], rows, ar_op, center_obs = center_obs, obs_weights = NULL
            )
            U_sub <- if (has_rwls_w) {
              compute_ar1_operator_factor(
                B_demix_l[[li]], rows, ar_op, center_obs = center_obs, obs_weights = fw_eff, weight_space = rwls_w_sp
              )
            } else {
              U_sub_0
            }
            U_full[, ok_t] <- U_sub
            if (is.null(sp$mix)) {
              V_full[, ok_t] <- crossprod(W_tilde_0, U_sub_0)
            } else {
              V_mix <- matrix(0.0, m, ncol(U_sub_0))
              for (mi in seq_along(sp$mix)) {
                W_m <- hrf_whitening_matrix(rows, sp$mix[[mi]]$ar, unit_marginal = TRUE)
                V_mix <- V_mix + sp$mix[[mi]]$weight * crossprod(W_m / c_s, U_sub_0)
              }
              V_full[, ok_t] <- V_mix
            }
          }
          U_l[[li]] <- U_full
          V_l[[li]] <- V_full
        }
        U_strata[[si]] <- U_l
        V_obs_strata[[si]] <- V_l
        ok_strata[[si]] <- ok_tau_l
      }
    }

    e_task <- std_env(h_task_accum)
    E_cov_env <- if (ncol(x_env_accum) > 0L) {
      matrix(vapply(seq_len(ncol(x_env_accum)), function(j) std_env(x_env_accum[, j]), numeric(m)),
             nrow = m, ncol = ncol(x_env_accum))
    } else {
      matrix(0.0, m, 0L)
    }

    vc_raw_mat <- matrix(0.0, m, 0L)
    if (!is.null(vc_resolved) && length(vc_resolved) >= r && !is.null(vc_resolved[[r]])) {
      V_in <- as.matrix(vc_resolved[[r]])
      storage.mode(V_in) <- "double"
      if (nrow(V_in) == n) V_in <- V_in[rows, , drop = FALSE]
      if (ncol(V_in) > 0L && nrow(V_in) != m) {
        stop(
          sprintf(
            "variance_confounds[['%s']] has %d rows, expected %d (full run volumes) or %d (usable uncensored volumes).",
            fit$run_ids[r], nrow(vc_resolved[[r]]), n, m
          ),
          call. = FALSE
        )
      }
      if (nrow(V_in) == m && ncol(V_in) > 0L) {
        for (jc in seq_len(ncol(V_in))) {
          bad_c <- !is.finite(V_in[, jc])
          if (any(bad_c)) {
            mc <- mean(V_in[!bad_c, jc])
            V_in[bad_c, jc] <- if (is.finite(mc)) mc else 0.0
          }
        }
        vc_raw_mat <- V_in
      }
    }

    d_mat_r <- abs(outer(rows, rows, "-"))
    adj_idx1 <- if (m >= 2L) which(abs(rows[2:m] - rows[1:(m - 1L)]) == 1L) else integer(0L)
    lag_pairs_r <- lapply(seq_len(L_max), function(d) which(d_mat_r == d, arr.ind = TRUE))
    fw_vec_r <- if (!is.null(fw_r) && length(fw_r) == m) as.numeric(fw_r) else rep(1.0, m)

    # For lowrank_voxel when y_obs_r was used (512 < P <= 16384), compress run residuals to rank <= 64 sketch
    vox_sketch_r <- sketch_r
    if (!exact_mode && is.null(vox_sketch_r) && !is.null(E_vox_full)) {
      K_sk_r <- min(m, 64L)
      nrm_r_v <- sqrt(pmax(colSums(E_rob_vox_full^2), 1e-16))
      Z_r_tmp <- sweep(E_rob_vox_full, 2L, nrm_r_v, "/")
      G_r_tmp <- tcrossprod(Z_r_tmp)
      eg_r_tmp <- eigen(0.5 * (G_r_tmp + t(G_r_tmp)), symmetric = TRUE)
      U_sk_r <- eg_r_tmp$vectors[, seq_len(K_sk_r), drop = FALSE]
      vox_sketch_r <- list(
        U_sk = U_sk_r,
        B_sk = crossprod(U_sk_r, E_vox_full),
        B_rob_sk = crossprod(U_sk_r, E_rob_vox_full)
      )
    } else if (!exact_mode && !is.null(vox_sketch_r) && !isTRUE(restore_scale)) {
      vox_sketch_r <- list(
        U_sk = vox_sketch_r$U_sk,
        B_sk = sweep(vox_sketch_r$B_sk, 2L, scale_mult_vox, "*"),
        B_rob_sk = sweep(vox_sketch_r$B_rob_sk, 2L, scale_mult_vox, "*")
      )
    }

    df_u_tot <- df_u_tot + mean(df_strata)
    run_store[[r]] <- list(
      m = m, rows = rows, m_target = m_target,
      E = E_unit, E_rob = E_rob_unit, E_vox_full = E_vox_full, vox_sketch = vox_sketch_r,
      s2_raw_u = s2_raw_u, s2_rob_u = s2_rob_u, s_blk_u = s_blk_u, df_unit_v = df_unit_v,
      parcel_idio_s2 = parcel_idio_s2,
      s_idx = s_idx_by_unit, S_cnt = S_cnt,
      M_strata = M_strata, G_inn_strata = G_inn_strata, H_strata = H_strata,
      U_strata = U_strata, V_obs_strata = V_obs_strata, ok_strata = ok_strata,
      e_task = e_task, E_cov_env = E_cov_env,
      fw_vec = fw_vec_r, vc_raw_mat = vc_raw_mat, saturated = any(sat_strata),
      rho_unit = rho_unit, df_strata = df_strata, c_strata = c_strata, kpeak_strata = kpeak_strata,
      adj_idx1 = adj_idx1, lag_pairs = lag_pairs_r
    )
  }

  keep_runs <- which(!vapply(run_store, is.null, logical(1L)))
  if (!length(keep_runs)) {
    return(list(unit_error_cor = NULL, shared_loading_vec = stats::setNames(numeric(P), unit_ids), shared_error_model = NULL))
  }
  run_store <- run_store[keep_runs]
  run_ids_active <- fit$run_ids[keep_runs]
  n_runs_act <- length(run_store)

  # 3A: Append optional run-level variance_confounds with pooled SD and quadratic quiet-anchor for non-negative profiles
  j_rwls <- integer(0L)
  vc_col_idx <- integer(0L)
  vc_pos_idx <- integer(0L)
  k_vc_in <- min(vapply(run_store, function(rs) ncol(rs$vc_raw_mat), integer(1L)))
  if (is.finite(k_vc_in) && k_vc_in >= 1L) {
    for (jc in seq_len(k_vc_in)) {
      raw_jc_all <- unlist(lapply(run_store, function(rs) rs$vc_raw_mat[, jc]), use.names = FALSE)
      is_nonneg_jc <- is.finite(min(raw_jc_all)) && (min(raw_jc_all) >= -1e-8)
      if (is_nonneg_jc) {
        z2_all <- raw_jc_all^2
        q0_jc <- as.numeric(stats::quantile(z2_all, probs = 0.15, names = FALSE))
        vc_prof <- lapply(run_store, function(rs) pmax(0.0, (rs$vc_raw_mat[, jc]^2) - q0_jc))
        sd_vc_pool <- stats::sd(unlist(vc_prof, use.names = FALSE))
      } else {
        vc_prof <- lapply(run_store, function(rs) rs$vc_raw_mat[, jc] - mean(rs$vc_raw_mat[, jc]))
        sd_vc_pool <- stats::sd(unlist(vc_prof, use.names = FALSE))
      }
      if (is.finite(sd_vc_pool) && sd_vc_pool > 1e-10) {
        f_exist <- do.call(rbind, lapply(run_store, function(rs) cbind(rs$e_task, rs$E_cov_env)))
        vc_all <- unlist(vc_prof, use.names = FALSE) / sd_vc_pool
        max_cor_vc <- if (ncol(f_exist) > 0L) max(abs(stats::cor(vc_all, f_exist)), na.rm = TRUE) else 0.0
        if (!is.finite(max_cor_vc) || max_cor_vc < 0.992) {
          for (r in seq_len(n_runs_act)) {
            run_store[[r]]$E_cov_env <- cbind(run_store[[r]]$E_cov_env, vc_prof[[r]] / sd_vc_pool)
          }
          j_new_vc <- 1L + ncol(run_store[[1L]]$E_cov_env)
          vc_col_idx <- c(vc_col_idx, j_new_vc)
          if (is_nonneg_jc) vc_pos_idx <- c(vc_pos_idx, j_new_vc)
        }
      }
    }
  }

  # Empirical-Bayes run-scale heterogeneity standardization across multi-run datasets:
  # When runs have homogeneous noise variance, firm gating yields run_scale == 1.0 identically;
  # when runs differ substantially in scanner/physiological noise SD, standardizing each run
  # prevents noisy runs from dominating POET/covreg and propagates exact per-run SE scales.
  for (r in seq_len(n_runs_act)) run_store[[r]]$run_scale <- rep(1.0, P)
  if (n_runs_act >= 2L && P >= 1L) {
    df_run_m <- matrix(1.0, nrow = n_runs_act, ncol = P)
    vrob_run_m <- matrix(0.0, nrow = n_runs_act, ncol = P)
    for (r in seq_len(n_runs_act)) {
      rs <- run_store[[r]]
      if (exact_mode) {
        s_sf <- ifelse(is.na(rs$s_idx) | rs$s_idx < 1L, 1L, rs$s_idx)
        df_rg <- pmax(1.0, rs$df_strata[s_sf])
        df_run_m[r, ] <- df_rg
        vrob_run_m[r, ] <- colSums(rs$E_rob^2) / df_rg + rs$parcel_idio_s2
      } else {
        df_rg <- pmax(1.0, rs$df_unit_v)
        df_run_m[r, ] <- df_rg
        vrob_run_m[r, ] <- rs$s2_rob_u
      }
    }
    df_col_sum <- pmax(colSums(df_run_m), 1.0)
    v_bar_g <- colSums(df_run_m * vrob_run_m) / df_col_sum
    ok_vg <- is.finite(v_bar_g) & (v_bar_g > 1e-12)
    if (any(ok_vg)) {
      l_mat <- matrix(0.0, nrow = n_runs_act, ncol = P)
      l_mat[, ok_vg] <- 0.5 * log(
        pmax(vrob_run_m[, ok_vg, drop = FALSE], 1e-14) /
          outer(rep(1.0, n_runs_act), pmax(v_bar_g[ok_vg], 1e-14))
      )
      se_l_mat <- sqrt(pmax((1.0 - 1.0 / n_runs_act) / (2.0 * pmax(df_run_m, 4.0)), 1e-8))
      z_run_mat <- abs(l_mat) / se_l_mat
      sh_run <- sqrt(2.0 * log(max(1.0, as.numeric(n_runs_act) * as.numeric(P))))
      w_het_glob <- .poet_firm_gate(max(z_run_mat[, ok_vg, drop = FALSE], na.rm = TRUE) - sh_run, lo = 1.25, hi = 2.60)
      if (w_het_glob > 1e-8) {
        w_eb_m <- w_het_glob * pmax(0.0, 1.0 - (1.5 * (se_l_mat^2)) / pmax(l_mat^2, 1e-12))
        sc_raw <- exp(w_eb_m * l_mat)
        norm_g <- sqrt(pmax(colSums(df_run_m * (sc_raw^2)) / df_col_sum, 1e-8))
        run_sc_m <- sweep(sc_raw, 2L, norm_g, "/")
        run_sc_m[, !ok_vg] <- 1.0
        for (r in seq_len(n_runs_act)) {
          sc_r <- pmin(20.0, pmax(0.05, run_sc_m[r, ]))
          inv_r <- 1.0 / sc_r
          run_store[[r]]$run_scale <- sc_r
          if (exact_mode) {
            run_store[[r]]$E <- sweep(run_store[[r]]$E, 2L, inv_r, "*")
            run_store[[r]]$E_rob <- sweep(run_store[[r]]$E_rob, 2L, inv_r, "*")
            run_store[[r]]$parcel_idio_s2 <- run_store[[r]]$parcel_idio_s2 * (inv_r^2)
          } else {
            run_store[[r]]$s2_raw_u <- run_store[[r]]$s2_raw_u * (inv_r^2)
            run_store[[r]]$s2_rob_u <- run_store[[r]]$s2_rob_u * (inv_r^2)
            run_store[[r]]$s_blk_u <- run_store[[r]]$s_blk_u * (inv_r^2)
            if (!is.null(run_store[[r]]$E_vox_full)) {
              run_store[[r]]$E_vox_full <- sweep(run_store[[r]]$E_vox_full, 2L, inv_r, "*")
            }
            if (!is.null(run_store[[r]]$vox_sketch)) {
              run_store[[r]]$vox_sketch$B_sk <- sweep(run_store[[r]]$vox_sketch$B_sk, 2L, inv_r, "*")
              run_store[[r]]$vox_sketch$B_rob_sk <- sweep(run_store[[r]]$vox_sketch$B_rob_sk, 2L, inv_r, "*")
            }
          }
        }
      }
    }
  }

  # =========================================================================
  # PATH 1: EXACT v3 POET + COVARIANCE REGRESSION ENGINE (P <= max_exact_units)
  # =========================================================================
  if (exact_mode) {
    S_num_1lag <- matrix(0.0, P, P, dimnames = list(unit_ids, unit_ids))
    S_rob_mat <- matrix(0.0, P, P, dimnames = list(unit_ids, unit_ids))
    S_blocked_diag <- numeric(P)
    S_den_1lag <- matrix(0.0, P, P)
    idio_s2_accum <- numeric(P)
    any_rwls_active <- FALSE

    for (r in seq_len(n_runs_act)) {
      rs <- run_store[[r]]
      S_num_1lag <- S_num_1lag + crossprod(rs$E)
      S_rob_mat <- S_rob_mat + crossprod(rs$E_rob)
      sup_idx <- which(rs$fw_vec < 0.85)
      if (length(sup_idx) > 0L) {
        any_rwls_active <- TRUE
        w_sup2 <- pmin(1.0, pmax(0.0, rs$fw_vec[sup_idx]^2))
        ex_sup <- matrix(pmax(0.0, rs$E[sup_idx, , drop = FALSE]^2 - rs$E_rob[sup_idx, , drop = FALSE]^2),
                         nrow = length(sup_idx), ncol = P)
        S_blocked_diag <- S_blocked_diag + colSums((1.0 - w_sup2) * ex_sup)
      }
      M_flat <- vapply(rs$M_strata, as.numeric, numeric(rs$m * rs$m))
      if (!is.matrix(M_flat)) M_flat <- matrix(M_flat, ncol = rs$S_cnt)
      tr_MM <- crossprod(M_flat)
      s_safe <- ifelse(is.na(rs$s_idx) | rs$s_idx < 1L, 1L, rs$s_idx)
      S_den_1lag <- S_den_1lag + tr_MM[s_safe, s_safe, drop = FALSE]
      if (any(rs$parcel_idio_s2 > 0)) {
        idio_s2_accum <- idio_s2_accum + rs$parcel_idio_s2 * diag(tr_MM)[s_safe]
      }
    }
    if (any(idio_s2_accum > 0)) {
      diag(S_num_1lag) <- diag(S_num_1lag) + idio_s2_accum
      diag(S_rob_mat) <- diag(S_rob_mat) + idio_s2_accum
    }

    Sig0_1lag_raw <- S_num_1lag / pmax(S_den_1lag, 1e-8)
    Sig0_1lag_raw <- 0.5 * (Sig0_1lag_raw + t(Sig0_1lag_raw))
    sig0_raw <- sqrt(pmax(diag(Sig0_1lag_raw), 1e-10))
    sig0_rob <- pmin(sig0_raw, sqrt(pmax(diag(S_rob_mat) / pmax(diag(S_den_1lag), 1e-8), 1e-10)))
    sig0_sd <- if (any_rwls_active) {
      sqrt(pmax(sig0_rob^2, (diag(S_num_1lag) - S_blocked_diag) / pmax(diag(S_den_1lag), 1e-8)))
    } else {
      sig0_raw
    }
    if (any_rwls_active) {
      S_rob_cov <- 0.5 * ((S_rob_mat / pmax(S_den_1lag, 1e-8)) + t(S_rob_mat / pmax(S_den_1lag, 1e-8)))
      sd_rob_u <- sqrt(pmax(diag(S_rob_cov), 1e-10))
      R0_debiased <- S_rob_cov / outer(sd_rob_u, sd_rob_u)
    } else {
      R0_debiased <- Sig0_1lag_raw / outer(sig0_raw, sig0_raw)
    }
    poet_fit <- poet_regularize_cor(R0_debiased, df_u = df_u_tot, k_max = k_max_poet, poet_mode = "adaptive_logp")
    R_poet <- poet_fit$R
    shared_loading_vec <- compute_shared_loading_from_cor(R_poet, unit_ids)

    # Fast exit when neither panel nor epochs are requested (univariate summary call)
    if (!isTRUE(return_panel) && !isTRUE(return_epochs)) {
      return(list(
        unit_error_cor = R_poet,
        shared_loading_vec = shared_loading_vec,
        shared_error_model = NULL
      ))
    }

    r_supp_tol <- 1e-3
    R_off_abs <- abs(R_poet - diag(diag(R_poet)))
    r_max_by_unit <- if (P >= 2L) apply(R_off_abs, 1L, max) else numeric(P)
    ltz_units <- which(r_max_by_unit > r_supp_tol)

    visited_u <- logical(P)
    comp_list <- list()
    for (g in seq_len(P)) {
      if (visited_u[g]) next
      if (r_max_by_unit[g] <= r_supp_tol) {
        visited_u[g] <- TRUE
        next
      }
      q_bfs <- g
      visited_u[g] <- TRUE
      idx_head <- 1L
      while (idx_head <= length(q_bfs)) {
        cur <- q_bfs[idx_head]; idx_head <- idx_head + 1L
        nbrs <- which(!visited_u & abs(R_poet[cur, ]) > r_supp_tol)
        if (length(nbrs)) {
          visited_u[nbrs] <- TRUE
          q_bfs <- c(q_bfs, nbrs)
        }
      }
      if (length(q_bfs) >= 2L) comp_list[[length(comp_list) + 1L]] <- sort(q_bfs)
    }

    M_dyn <- 1L + ncol(run_store[[1L]]$E_cov_env)
    K_reg <- 1L + M_dyn
    B_tot <- L1 + 2L * M_dyn

    s_mat <- matrix(
      vapply(run_store, function(rs) ifelse(is.na(rs$s_idx) | rs$s_idx < 1L, 1L, rs$s_idx), integer(P)),
      nrow = P, ncol = n_runs_act
    )
    unit_sig <- apply(s_mat, 1L, paste, collapse = ":")
    u_sigs <- unique(unit_sig)
    G_sig <- length(u_sigs)
    sig_idx <- match(unit_sig, u_sigs)
    rep_unit_for_sig <- match(u_sigs, unit_sig)

    active_off <- if (P >= 2L) which(row(R_poet) < col(R_poet) & abs(R_poet) > r_supp_tol, arr.ind = TRUE) else matrix(0L, 0L, 2L)
    n_act_off <- nrow(active_off)

    X_diag_cache <- lapply(seq_len(n_runs_act), function(r) vector("list", G_sig))
    XtX_diag_sig <- vector("list", G_sig)
    for (a in seq_len(G_sig)) {
      ua <- rep_unit_for_sig[a]
      XtX_aa <- matrix(0.0, B_tot, B_tot)
      for (r in seq_len(n_runs_act)) {
        rs <- run_store[[r]]; sa <- s_mat[ua, r]
        f_dyn <- cbind(rs$e_task, rs$E_cov_env)
        if (rs$saturated) f_dyn[] <- 0.0
        X_ra <- poet_build_pair_X_cpp(
          rs$M_strata[[sa]], rs$M_strata[[sa]],
          rs$G_inn_strata[[sa]], rs$G_inn_strata[[sa]],
          rs$H_strata[[sa]], rs$H_strata[[sa]],
          rs$adj_idx1, rs$lag_pairs, f_dyn, L_max
        )
        X_diag_cache[[r]][[a]] <- X_ra
        XtX_aa <- XtX_aa + crossprod(X_ra)
      }
      XtX_diag_sig[[a]] <- XtX_aa
    }

    E_runs_list <- lapply(run_store, function(rs) rs$E)
    adj_runs_list <- lapply(run_store, function(rs) as.integer(rs$adj_idx1))
    lp_runs_list <- lapply(run_store, function(rs) rs$lag_pairs)
    Xty_diag <- poet_extract_Xty_diag_cpp(
      E_runs_list, X_diag_cache, as.integer(sig_idx),
      adj_runs_list, lp_runs_list, L_max, B_tot
    )

    gamma_obs <- 0.80
    eval_units_gam <- if (length(ltz_units) >= 1L) ltz_units else seq_len(P)
    if (!any_saturated && length(eval_units_gam) >= 1L) {
      keep_gam_dyn <- seq_len(min(M_dyn, 1L + q_cov))
      n_gam_dyn <- length(keep_gam_dyn)
      grid_gam <- seq(0.0, 1.0, by = 0.05)
      prior_wt <- 0.0025 * sum(sig0_sd[eval_units_gam]^4)
      obj_gam <- vapply(grid_gam, function(gm) {
        A_gm <- matrix(0.0, B_tot, L1 + n_gam_dyn)
        A_gm[1:L1, 1:L1] <- diag(L1)
        for (ii in seq_len(n_gam_dyn)) {
          j <- keep_gam_dyn[ii]
          A_gm[L1 + j, L1 + ii] <- 1.0 - gm
          A_gm[L1 + M_dyn + j, L1 + ii] <- gm
        }
        ssr_val <- 0.0
        for (g in eval_units_gam) {
          XtX_gg <- XtX_diag_sig[[sig_idx[g]]]
          XtX_r <- crossprod(A_gm, XtX_gg %*% A_gm)
          act_c <- which(diag(XtX_r) > 1e-8)
          if (length(act_c)) {
            rhs_r <- as.numeric(crossprod(A_gm, Xty_diag[, g]))[act_c]
            b_r <- solve(XtX_r[act_c, act_c, drop = FALSE] + 1e-6 * diag(length(act_c)), rhs_r)
            ssr_val <- ssr_val + sum(rhs_r * b_r)
          }
        }
        ssr_val - prior_wt * ((gm - 0.80)^2)
      }, numeric(1L))
      if (any(is.finite(obj_gam))) {
        gamma_obs <- grid_gam[which.max(ifelse(is.finite(obj_gam), obj_gam, -Inf))[1L]]
      }
    }

    A_tie_full <- matrix(0.0, B_tot, L1 + M_dyn)
    A_tie_full[1:L1, 1:L1] <- diag(L1)
    for (j in seq_len(M_dyn)) {
      gm_j <- if (length(j_rwls) == 1L && j == j_rwls) 1.0 else gamma_obs
      A_tie_full[L1 + j, L1 + j] <- 1.0 - gm_j
      A_tie_full[L1 + M_dyn + j, L1 + j] <- gm_j
    }
    A_tie_1lag <- matrix(0.0, B_tot, 1L + M_dyn)
    A_tie_1lag[1L, 1L] <- 1.0
    for (j in seq_len(M_dyn)) {
      gm_j <- if (length(j_rwls) == 1L && j == j_rwls) 1.0 else gamma_obs
      A_tie_1lag[L1 + j, 1L + j] <- 1.0 - gm_j
      A_tie_1lag[L1 + M_dyn + j, 1L + j] <- gm_j
    }

    Q_port <- numeric(P)
    z_tail_g <- numeric(P)
    w_bart_lf <- numeric(L_max)
    if (L_max >= 2L) {
      idx_tl <- 2:min(6L, L_max)
      w_bart_lf[idx_tl] <- 1.0 - (seq_along(idx_tl) - 1.0) / length(idx_tl)
    } else {
      w_bart_lf[1L] <- 1.0
    }
    sig0_mlag_diag <- numeric(P)
    for (g in seq_len(P)) {
      XtX_gg <- XtX_diag_sig[[sig_idx[g]]]
      XtX_tf <- crossprod(A_tie_full, XtX_gg %*% A_tie_full)
      act_f <- which(diag(XtX_tf) > 1e-8)
      b_tf <- numeric(L1 + M_dyn)
      if (length(act_f)) {
        rhs_tf <- as.numeric(crossprod(A_tie_full, Xty_diag[, g]))[act_f]
        b_tf[act_f] <- solve(XtX_tf[act_f, act_f, drop = FALSE] + 1e-6 * diag(length(act_f)), rhs_tf)
      }
      XtX_s <- XtX_gg[1:L1, 1:L1, drop = FALSE]
      gam_u <- solve(XtX_s + 1e-6 * diag(L1), Xty_diag[1:L1, g])
      sig0_mlag_diag[g] <- sqrt(max(gam_u[1L], 0.15 * sig0_sd[g]^2))

      s2_cg <- max(b_tf[1L], 0.15 * sig0_sd[g]^2)
      acf_g_raw <- pmin(0.95, pmax(-0.95, b_tf[2:L1] / s2_cg))
      V_inv_sub <- XtX_s[2:L1, 2:L1, drop = FALSE] - XtX_s[2:L1, 1L, drop = FALSE] %*% solve(XtX_s[1L, 1L, drop = FALSE], XtX_s[1L, 2:L1, drop = FALSE])
      gam_sub <- s2_cg * acf_g_raw
      Q_port[g] <- as.numeric(crossprod(gam_sub, V_inv_sub %*% gam_sub)) / (2.0 * (s2_cg^2))
      V_sub_cov_w <- as.numeric(crossprod(w_bart_lf, solve(V_inv_sub + 1e-6 * diag(L_max), w_bart_lf)))
      z_tail_g[g] <- sum(w_bart_lf * acf_g_raw) / sqrt(max(2.0 * V_sub_cov_w, 1e-12))
    }

    r2_bar <- if (P >= 2L) mean((R_poet^2)[row(R_poet) != col(R_poet)]) else 0.0
    p_eff_port <- P / (1.0 + (P - 1.0) * max(0.0, r2_bar))
    z_omni_glob <- (mean(Q_port) - L_max) / max(sqrt(2.0 * L_max / max(1.0, p_eff_port)), 1e-6)
    z_signed_glob <- sqrt(max(1.0, p_eff_port)) * mean(z_tail_g)
    z_comp_max <- -Inf
    if (length(comp_list) >= 1L) {
      shift_comp_port <- sqrt(max(0.0, 2.0 * log(max(1.0, as.numeric(length(comp_list))))))
      for (Cu_p in comp_list) {
        nCp <- length(Cu_p)
        R_Cp <- R_poet[Cu_p, Cu_p, drop = FALSE]
        r2_Cp <- if (nCp >= 2L) mean((R_Cp^2)[row(R_Cp) != col(R_Cp)]) else 0.0
        p_eff_Cp <- nCp / (1.0 + (nCp - 1.0) * max(0.0, r2_Cp))
        z_Cp_omni <- (mean(Q_port[Cu_p]) - L_max) / max(sqrt(2.0 * L_max / max(1.0, p_eff_Cp)), 1e-6) - shift_comp_port
        z_Cp_signed <- abs(sqrt(max(1.0, p_eff_Cp)) * mean(z_tail_g[Cu_p])) - shift_comp_port
        z_comp_max <- max(z_comp_max, z_Cp_omni - 0.35, z_Cp_signed - 0.45)
      }
    }
    z_port <- max(z_omni_glob, abs(z_signed_glob) - 0.15, z_comp_max)
    w_res_global <- .poet_firm_gate(z_port, lo = 1.55, hi = 2.75)

    use_mlag <- (w_res_global > 1e-6)
    A_tie_use <- if (use_mlag) A_tie_full else A_tie_1lag
    L_act <- if (use_mlag) L1 else 1L

    off_sig_cache <- new.env(parent = emptyenv())
    get_exact_sig_pair <- function(a, b) {
      if (a == b) return(list(XtX = XtX_diag_sig[[a]], X_runs = lapply(seq_len(n_runs_act), function(r) X_diag_cache[[r]][[a]])))
      ia <- min(a, b); ib <- max(a, b)
      key <- paste0(ia, ":", ib)
      if (exists(key, envir = off_sig_cache, inherits = FALSE)) return(get(key, envir = off_sig_cache))
      ua <- rep_unit_for_sig[ia]; ub <- rep_unit_for_sig[ib]
      XtX_ab <- matrix(0.0, B_tot, B_tot)
      X_runs <- vector("list", n_runs_act)
      for (r in seq_len(n_runs_act)) {
        rs <- run_store[[r]]; sa <- s_mat[ua, r]; sb <- s_mat[ub, r]
        f_dyn <- cbind(rs$e_task, rs$E_cov_env)
        if (rs$saturated) f_dyn[] <- 0.0
        Xr <- poet_build_pair_X_cpp(
          rs$M_strata[[sa]], rs$M_strata[[sb]],
          rs$G_inn_strata[[sa]], rs$G_inn_strata[[sb]],
          rs$H_strata[[sa]], rs$H_strata[[sb]],
          rs$adj_idx1, rs$lag_pairs, f_dyn, L_max
        )
        X_runs[[r]] <- Xr
        XtX_ab <- XtX_ab + crossprod(Xr)
      }
      res <- list(XtX = XtX_ab, X_runs = X_runs)
      assign(key, res, envir = off_sig_cache)
      res
    }

    u_off_keys <- if (n_act_off > 0L) {
      sa_v <- pmin(sig_idx[active_off[, 1L]], sig_idx[active_off[, 2L]])
      sb_v <- pmax(sig_idx[active_off[, 1L]], sig_idx[active_off[, 2L]])
      unique(paste0(sa_v[sa_v != sb_v], ":", sb_v[sa_v != sb_v]))
    } else character(0L)
    dense_rep_mode <- (length(u_off_keys) > 240L)

    if (dense_rep_mode) {
      feat_u <- cbind(
        scale(vapply(seq_len(P), function(g) mean(vapply(seq_len(n_runs_act), function(r) run_store[[r]]$df_strata[s_mat[g, r]], numeric(1L))), numeric(1L))),
        scale(vapply(seq_len(P), function(g) mean(vapply(seq_len(n_runs_act), function(r) run_store[[r]]$c_strata[s_mat[g, r]], numeric(1L))), numeric(1L))),
        scale(vapply(seq_len(P), function(g) mean(vapply(seq_len(n_runs_act), function(r) run_store[[r]]$kpeak_strata[s_mat[g, r]], numeric(1L))), numeric(1L)))
      )
      feat_u[!is.finite(feat_u)] <- 0.0
      G_rep <- min(8L, G_sig)
      ord_f <- order(feat_u[, 1L] + 0.5 * feat_u[, 2L] + 0.25 * feat_u[, 3L])
      rep_grp <- integer(P)
      rep_grp[ord_f] <- as.integer(ceiling(seq_len(P) * G_rep / P))
      rep_u_vec <- vapply(seq_len(G_rep), function(k) {
        idx_k <- which(rep_grp == k)
        idx_k[ceiling(length(idx_k) / 2.0)]
      }, integer(1L))
      X_rep_runs <- lapply(seq_len(n_runs_act), function(r) vector("list", G_rep * G_rep))
      XtX_rep <- vector("list", G_rep * G_rep)
      for (ka in seq_len(G_rep)) {
        ua <- rep_u_vec[ka]
        for (kb in seq.int(ka, G_rep)) {
          ub <- rep_u_vec[kb]
          XtX_ab <- matrix(0.0, B_tot, B_tot)
          for (r in seq_len(n_runs_act)) {
            rs <- run_store[[r]]; sa <- s_mat[ua, r]; sb <- s_mat[ub, r]
            f_dyn <- cbind(rs$e_task, rs$E_cov_env)
            if (rs$saturated) f_dyn[] <- 0.0
            Xr <- poet_build_pair_X_cpp(
              rs$M_strata[[sa]], rs$M_strata[[sb]],
              rs$G_inn_strata[[sa]], rs$G_inn_strata[[sb]],
              rs$H_strata[[sa]], rs$H_strata[[sb]],
              rs$adj_idx1, rs$lag_pairs, f_dyn, L_max
            )
            X_rep_runs[[r]][[(ka - 1L) * G_rep + kb]] <- Xr
            XtX_ab <- XtX_ab + crossprod(Xr)
          }
          XtX_rep[[(ka - 1L) * G_rep + kb]] <- XtX_ab
        }
      }
      Xty_bulk <- array(0.0, dim = c(B_tot, P, P))
      for (ka in seq_len(G_rep)) {
        ca <- which(rep_grp == ka)
        for (kb in seq.int(ka, G_rep)) {
          cb <- which(rep_grp == kb)
          for (r in seq_len(n_runs_act)) {
            rs <- run_store[[r]]; m <- rs$m; idx1 <- rs$adj_idx1
            Xr <- X_rep_runs[[r]][[(ka - 1L) * G_rep + kb]]
            EA <- rs$E[, ca, drop = FALSE]; EB <- rs$E[, cb, drop = FALSE]
            for (ell in seq_len(B_tot)) {
              w0 <- Xr[1:m, ell]
              M_ell <- crossprod(EA, w0 * EB)
              if (length(idx1)) {
                w1 <- Xr[m + idx1, ell] / sqrt(2.0)
                M_ell <- M_ell + crossprod(EA[idx1 + 1L, , drop = FALSE], w1 * EB[idx1, , drop = FALSE]) +
                  crossprod(EA[idx1, , drop = FALSE], w1 * EB[idx1 + 1L, , drop = FALSE])
              }
              if (L_max >= 2L) {
                for (d in 2:L_max) {
                  wd <- Xr[2L * m - 1L + d - 1L, ell]
                  if (abs(wd) > 1e-14) {
                    pd <- rs$lag_pairs[[d]]; N_d <- 0.5 * nrow(pd)
                    if (N_d > 0) {
                      M_ell <- M_ell + (wd / sqrt(2.0 * N_d)) * crossprod(EA[pd[, 1L], , drop = FALSE], EB[pd[, 2L], , drop = FALSE])
                    }
                  }
                }
              }
              Xty_bulk[ell, ca, cb] <- Xty_bulk[ell, ca, cb] + M_ell
              if (ka != kb) Xty_bulk[ell, cb, ca] <- Xty_bulk[ell, cb, ca] + t(M_ell)
            }
          }
        }
      }
      for (g in seq_len(P)) Xty_bulk[, g, g] <- Xty_diag[, g]
    }

    # Stage D: C++-Accelerated Component-Coherent Gating + Schur-Complement VC Restoration + PSD Cone
    B_dyn_gated <- array(0.0, dim = c(M_dyn, P, P))
    w_ltz_global <- numeric(M_dyn)
    all_f_frames <- do.call(rbind, lapply(run_store, function(rs) cbind(rs$e_task, rs$E_cov_env)))

    sig0_use <- if (use_mlag) {
      sqrt((1.0 - w_res_global) * (sig0_sd^2) + w_res_global * (sig0_mlag_diag^2))
    } else {
      sig0_sd
    }
    Gamma_cond_diag <- matrix(0.0, L_act, P)
    for (g in seq_len(P)) {
      Gamma_cond_diag[1L, g] <- sig0_use[g]^2
      if (use_mlag) {
        XtX_s <- XtX_diag_sig[[sig_idx[g]]][1:L1, 1:L1, drop = FALSE]
        gam_u <- solve(XtX_s + 1e-6 * diag(L1), Xty_diag[1:L1, g])
        Gamma_cond_diag[, g] <- gam_u
        Gamma_cond_diag[1L, g] <- sig0_use[g]^2
      }
    }
    Delta0_1lag_all <- matrix(0.0, P, P)

    for (C_idx in seq_along(comp_list)) {
      C_u <- comp_list[[C_idx]]
      nC <- length(C_u)
      R_C <- R_poet[C_u, C_u, drop = FALSE]
      supp_C <- abs(R_C) > r_supp_tol
      Sig0_C_base <- R_C * outer(sig0_use[C_u], sig0_use[C_u])

      off_pairs_C <- which(row(R_C) < col(R_C) & supp_C, arr.ind = TRUE)
      hac_mask_C <- matrix(FALSE, nC, nC)
      diag(hac_mask_C) <- TRUE
      if (nrow(off_pairs_C) <= 120L) {
        hac_mask_C[supp_C] <- TRUE
      } else {
        r_vals_off <- abs(R_C[off_pairs_C])
        keep_top <- utils::head(order(r_vals_off, decreasing = TRUE), 120L)
        for (kk in keep_top) {
          ia_k <- off_pairs_C[kk, 1L]; ib_k <- off_pairs_C[kk, 2L]
          hac_mask_C[ia_k, ib_k] <- hac_mask_C[ib_k, ia_k] <- TRUE
        }
      }

      M_core_c <- min(M_dyn, 1L + q_cov)
      n_vc_c <- length(vc_col_idx)
      B_raw_C <- array(0.0, dim = c(M_dyn, nC, nC))
      Z_C <- array(0.0, dim = c(M_dyn, nC, nC))
      Z_core_al_C <- array(0.0, dim = c(M_core_c, nC, nC))
      R_cv_C <- if (n_vc_c >= 1L) array(0.0, dim = c(M_core_c, n_vc_c, nC, nC)) else NULL
      R_sd_C <- array(0.0, dim = c(L_act, M_dyn, nC, nC))
      gam_unc_C <- array(0.0, dim = c(L_act, nC, nC))

      # Build unique design slots for all diagonal & active off-diagonal pairs in C_u
      slot_map <- new.env(parent = emptyenv())
      slot_XtX_list <- list()
      slot_X_runs_list <- list()
      get_slot_id <- function(key_str, sp_info) {
        if (exists(key_str, envir = slot_map, inherits = FALSE)) return(get(key_str, envir = slot_map))
        sid <- length(slot_XtX_list) + 1L
        slot_XtX_list[[sid]] <<- sp_info$XtX
        slot_X_runs_list[[sid]] <<- sp_info$X_runs
        assign(key_str, sid, envir = slot_map)
        sid
      }

      n_off_act <- nrow(off_pairs_C)
      n_batch <- nC + n_off_act
      pair_ia <- integer(n_batch)
      pair_ib <- integer(n_batch)
      pair_g <- integer(n_batch)
      pair_h <- integer(n_batch)
      pair_slot <- integer(n_batch)
      pair_hac <- logical(n_batch)
      pair_has_rhs <- logical(n_batch)
      rhs_over <- matrix(0.0, B_tot, n_batch)

      for (ia in seq_len(nC)) {
        g <- C_u[ia]
        a_sig <- sig_idx[g]
        sp_d <- list(XtX = XtX_diag_sig[[a_sig]], X_runs = lapply(seq_len(n_runs_act), function(r) X_diag_cache[[r]][[a_sig]]))
        pair_ia[ia] <- ia; pair_ib[ia] <- ia
        pair_g[ia] <- g; pair_h[ia] <- g
        pair_slot[ia] <- get_slot_id(paste0("d:", a_sig), sp_d)
        pair_hac[ia] <- TRUE
        pair_has_rhs[ia] <- TRUE
        rhs_over[, ia] <- Xty_diag[, g]
      }

      if (n_off_act > 0L) {
        for (kk in seq_len(n_off_act)) {
          idx_p <- nC + kk
          ia <- off_pairs_C[kk, 1L]; ib <- off_pairs_C[kk, 2L]
          g <- C_u[ia]; h <- C_u[ib]
          pair_ia[idx_p] <- ia; pair_ib[idx_p] <- ib
          pair_g[idx_p] <- g; pair_h[idx_p] <- h
          pair_hac[idx_p] <- hac_mask_C[ia, ib]
          if (!dense_rep_mode || abs(poet_fit$S_local[g, h]) > 1e-5) {
            sa_i <- min(sig_idx[g], sig_idx[h]); sb_i <- max(sig_idx[g], sig_idx[h])
            sp_o <- get_exact_sig_pair(sa_i, sb_i)
            pair_slot[idx_p] <- get_slot_id(paste0("e:", sa_i, ":", sb_i), sp_o)
            pair_has_rhs[idx_p] <- FALSE
          } else {
            ka <- min(rep_grp[g], rep_grp[h]); kb <- max(rep_grp[g], rep_grp[h])
            id_rep <- (ka - 1L) * G_rep + kb
            sp_r <- list(XtX = XtX_rep[[id_rep]], X_runs = lapply(seq_len(n_runs_act), function(r) X_rep_runs[[r]][[id_rep]]))
            pair_slot[idx_p] <- get_slot_id(paste0("r:", id_rep), sp_r)
            pair_has_rhs[idx_p] <- TRUE
            rhs_over[, idx_p] <- Xty_bulk[, g, h]
          }
        }
      }

      batch_fit <- poet_fit_pair_batch_cpp(
        E_runs_list, adj_runs_list, lp_runs_list,
        slot_XtX_list, slot_X_runs_list, A_tie_use,
        L_max, L_act, M_dyn, M_core_c, as.integer(vc_col_idx),
        sig0_sd, R_poet,
        pair_g, pair_h, pair_slot, pair_hac, pair_has_rhs, rhs_over
      )

      vif_hac_diag <- rep(1.25, M_dyn)
      diag_ratios <- matrix(NA_real_, nC, M_dyn)
      for (ia in seq_len(nC)) {
        B_raw_C[, ia, ia] <- batch_fit$b_dyn[, ia]
        R_sd_C[, , ia, ia] <- matrix(batch_fit$R_sd_flat[, ia], nrow = L_act, ncol = M_dyn)
        if (n_vc_c >= 1L) R_cv_C[, , ia, ia] <- matrix(batch_fit$R_cv_flat[, ia], nrow = M_core_c, ncol = n_vc_c)
        gam_unc_C[, ia, ia] <- batch_fit$gam_unc[, ia]
        for (j in seq_len(M_dyn)) {
          vj <- batch_fit$v_dyn[j, ia]
          vg0 <- batch_fit$v_gauss[j, ia]
          Z_C[j, ia, ia] <- if (is.finite(vj) && vj > 1e-16) batch_fit$b_dyn[j, ia] / sqrt(vj) else 0.0
          if (is.finite(vj) && is.finite(vg0) && vg0 > 1e-16) diag_ratios[ia, j] <- vj / vg0
        }
        if (n_vc_c >= 1L) {
          for (jc0 in seq_len(M_core_c)) {
            va0 <- batch_fit$v_core_alone[jc0, ia]
            Z_core_al_C[jc0, ia, ia] <- if (is.finite(va0) && va0 > 1e-16) batch_fit$b_core_alone[jc0, ia] / sqrt(va0) else Z_C[jc0, ia, ia]
          }
        }
      }
      for (j in seq_len(M_dyn)) {
        ok_r <- is.finite(diag_ratios[, j]) & diag_ratios[, j] > 0
        if (any(ok_r)) vif_hac_diag[j] <- max(1.0, stats::median(diag_ratios[ok_r, j]))
      }

      if (n_off_act > 0L) {
        for (kk in seq_len(n_off_act)) {
          idx_p <- nC + kk
          ia <- pair_ia[idx_p]; ib <- pair_ib[idx_p]
          use_h <- pair_hac[idx_p]
          B_raw_C[, ia, ib] <- B_raw_C[, ib, ia] <- batch_fit$b_dyn[, idx_p]
          rsd_m <- matrix(batch_fit$R_sd_flat[, idx_p], nrow = L_act, ncol = M_dyn)
          R_sd_C[, , ia, ib] <- rsd_m; R_sd_C[, , ib, ia] <- rsd_m
          if (n_vc_c >= 1L) {
            rcv_m <- matrix(batch_fit$R_cv_flat[, idx_p], nrow = M_core_c, ncol = n_vc_c)
            R_cv_C[, , ia, ib] <- rcv_m; R_cv_C[, , ib, ia] <- rcv_m
          }
          gam_unc_C[, ia, ib] <- gam_unc_C[, ib, ia] <- batch_fit$gam_unc[, idx_p]
          for (j in seq_len(M_dyn)) {
            vj <- batch_fit$v_dyn[j, idx_p] * (if (use_h) 1.0 else vif_hac_diag[j])
            zk <- if (is.finite(vj) && vj > 1e-16) batch_fit$b_dyn[j, idx_p] / sqrt(vj) else 0.0
            Z_C[j, ia, ib] <- Z_C[j, ib, ia] <- zk
          }
          if (n_vc_c >= 1L) {
            for (jc0 in seq_len(M_core_c)) {
              va0 <- batch_fit$v_core_alone[jc0, idx_p] * (if (use_h) 1.0 else vif_hac_diag[jc0])
              za0 <- if (is.finite(va0) && va0 > 1e-16) batch_fit$b_core_alone[jc0, idx_p] / sqrt(va0) else Z_C[jc0, ia, ib]
              Z_core_al_C[jc0, ia, ib] <- Z_core_al_C[jc0, ib, ia] <- za0
            }
          }
        }
      }

      off_cand_mask <- (row(R_C) < col(R_C)) & (abs(R_C) >= 0.10)
      n_cand_off <- sum(off_cand_mask)
      r_eff_C <- max(2.0, (sum(diag(R_C))^2) / max(sum(R_C^2), 1e-8))
      m_eff_off <- max(1.0, min(as.numeric(max(1L, n_cand_off)), 0.5 * r_eff_C * (r_eff_C - 1.0)))
      shift_mult <- 2.0 * log(m_eff_off)

      B_gated_C <- vector("list", M_dyn)
      if (!any_saturated) {
        w_vc_max_C <- 0.0
        if (n_vc_c >= 1L) {
          k_vc_mult <- max(1.0, as.numeric(n_vc_c))
          shift_vc <- 2.0 * log(k_vc_mult)
          for (j in vc_col_idx) {
            Bj_sym <- 0.5 * (B_raw_C[j, , ] + t(B_raw_C[j, , ]))
            is_pos_j <- j %in% vc_pos_idx
            sh_j <- shift_mult + shift_vc
            if (is_pos_j) {
              egj <- eigen(Bj_sym, symmetric = TRUE)
              Bj_psd <- egj$vectors %*% (pmax(egj$values, 0.0) * t(egj$vectors))
              Bj_psd[!supp_C] <- 0.0
              zj_diag_pos <- pmax(0.0, diag(Z_C[j, , ]))
              w_ltz_j <- .ltz_from_cor_quad(zj_diag_pos, R_C, p_lo = 0.12 / k_vc_mult, p_hi = 0.004 / k_vc_mult)
              zj_off_pos <- if (n_cand_off > 0L) max(0.0, max((sign(R_C) * Z_C[j, , ])[off_cand_mask])) else 0.0
              w_j <- max(w_ltz_j, .poet_firm_gate(zj_off_pos^2, lo = 2.2 + sh_j, hi = 5.8 + sh_j))
              w_ltz_global[j] <- max(w_ltz_global[j], w_j)
              w_vc_max_C <- max(w_vc_max_C, w_j)
              B_gated_C[[j]] <- w_j * Bj_psd
            } else {
              Bj_sym[!supp_C] <- 0.0
              zj_diag <- diag(Z_C[j, , ])
              w_ltz_j <- .ltz_from_cor_quad(zj_diag, R_C, p_lo = 0.08 / k_vc_mult, p_hi = 0.002 / k_vc_mult)
              z2_off_max <- if (n_cand_off > 0L) max(0.0, max((Z_C[j, , ]^2)[off_cand_mask])) else 0.0
              w_j <- max(w_ltz_j, .poet_firm_gate(z2_off_max, lo = 2.6 + sh_j, hi = 6.2 + sh_j))
              w_ltz_global[j] <- max(w_ltz_global[j], w_j)
              w_vc_max_C <- max(w_vc_max_C, w_j)
              B_gated_C[[j]] <- w_j * Bj_sym
            }
          }
          if (!is.null(R_cv_C)) {
            for (ia in seq_len(nC)) {
              for (ib in seq_len(nC)) {
                if (!supp_C[ia, ib]) next
                b_vc_raw <- B_raw_C[vc_col_idx, ia, ib]
                b_vc_gat <- vapply(vc_col_idx, function(jv) B_gated_C[[jv]][ia, ib], numeric(1L))
                R_cv_ab <- matrix(R_cv_C[, , ia, ib], nrow = M_core_c, ncol = n_vc_c)
                B_raw_C[1:M_core_c, ia, ib] <- B_raw_C[1:M_core_c, ia, ib] + as.numeric(R_cv_ab %*% (b_vc_raw - b_vc_gat))
                Z_C[1:M_core_c, ia, ib] <- (1.0 - w_vc_max_C) * Z_core_al_C[1:M_core_c, ia, ib] + w_vc_max_C * Z_C[1:M_core_c, ia, ib]
              }
            }
          }
        }

        B1_sym <- 0.5 * (B_raw_C[1L, , ] + t(B_raw_C[1L, , ]))
        eg1 <- eigen(B1_sym, symmetric = TRUE)
        B1_psd <- eg1$vectors %*% (pmax(eg1$values, 0.0) * t(eg1$vectors))
        B1_psd[!supp_C] <- 0.0
        z1_diag_pos <- pmax(0.0, diag(Z_C[1L, , ]))
        w_ltz_1 <- .ltz_from_cor_quad(z1_diag_pos, R_C, p_lo = 0.25, p_hi = 0.01)
        z1_off_pos <- if (n_cand_off > 0L) max(0.0, max((sign(R_C) * Z_C[1L, , ])[off_cand_mask])) else 0.0
        w_1 <- max(w_ltz_1, .poet_firm_gate(z1_off_pos^2, lo = 1.6 + shift_mult, hi = 5.0 + shift_mult))
        w_ltz_global[1L] <- max(w_ltz_global[1L], w_1)
        B_gated_C[[1L]] <- w_1 * B1_psd

        core_extra <- setdiff(seq_len(M_dyn), c(1L, vc_col_idx))
        if (length(core_extra) >= 1L) {
          for (j in core_extra) {
            Bj_sym <- 0.5 * (B_raw_C[j, , ] + t(B_raw_C[j, , ]))
            Bj_sym[!supp_C] <- 0.0
            zj_diag <- diag(Z_C[j, , ])
            w_ltz_j <- .ltz_from_cor_quad(zj_diag, R_C, p_lo = 0.20, p_hi = 0.005)
            z2_off_max <- if (n_cand_off > 0L) max(0.0, max((Z_C[j, , ]^2)[off_cand_mask])) else 0.0
            w_j <- max(w_ltz_j, .poet_firm_gate(z2_off_max, lo = 1.8 + shift_mult, hi = 5.5 + shift_mult))
            w_ltz_global[j] <- max(w_ltz_global[j], w_j)
            B_gated_C[[j]] <- w_j * Bj_sym
          }
        }
      } else {
        for (j in seq_len(M_dyn)) B_gated_C[[j]] <- matrix(0.0, nC, nC)
      }

      j_zero_anc <- union(j_rwls, vc_pos_idx)
      Delta0_za <- matrix(0.0, nC, nC)
      Delta0_oth <- matrix(0.0, nC, nC)
      for (j in seq_len(M_dyn)) {
        if (j %in% j_zero_anc) {
          Delta0_za <- Delta0_za - matrix(R_sd_C[1L, j, , ], nC, nC) * B_gated_C[[j]]
        } else {
          Delta0_oth <- Delta0_oth - matrix(R_sd_C[1L, j, , ], nC, nC) * B_gated_C[[j]]
        }
      }
      Delta0_za <- 0.5 * (Delta0_za + t(Delta0_za)); Delta0_za[!supp_C] <- 0.0
      Delta0_oth <- 0.5 * (Delta0_oth + t(Delta0_oth)); Delta0_oth[!supp_C] <- 0.0
      Delta0_C <- Delta0_za + Delta0_oth

      if (nrow(all_f_frames) > 0L) {
        Sig0_rob_C <- R_C * outer(sig0_rob[C_u], sig0_rob[C_u])
        eg0 <- eigen(0.5 * (Sig0_rob_C + t(Sig0_rob_C)), symmetric = TRUE)
        S0_isqrt <- eg0$vectors %*% ((1.0 / sqrt(pmax(eg0$values, 1e-10))) * t(eg0$vectors))
        S0_sqrt <- eg0$vectors %*% (sqrt(pmax(eg0$values, 1e-10)) * t(eg0$vectors))
        Base_s <- S0_isqrt %*% Sig0_C_base %*% S0_isqrt
        if (length(j_zero_anc) >= 1L && any(abs(Delta0_za) > 1e-12)) {
          M_za <- 0.5 * ((Base_s + S0_isqrt %*% Delta0_za %*% S0_isqrt) + t(Base_s + S0_isqrt %*% Delta0_za %*% S0_isqrt))
          eg_za <- eigen(M_za, symmetric = TRUE)
          min_base_ev <- min(0.92, min(eigen(0.5 * (Base_s + t(Base_s)), symmetric = TRUE, only.values = TRUE)$values))
          if (min(eg_za$values) < min_base_ev) {
            M_za_cl <- eg_za$vectors %*% (pmax(eg_za$values, min_base_ev) * t(eg_za$vectors))
            Delta0_za <- 0.5 * (((S0_sqrt %*% M_za_cl %*% S0_sqrt) - Sig0_C_base) + t((S0_sqrt %*% M_za_cl %*% S0_sqrt) - Sig0_C_base))
            Delta0_za[!supp_C] <- 0.0
          }
        }
        Delta0_C <- Delta0_za + Delta0_oth
        H_head <- 0.5 * (Base_s + t(Base_s)) - 0.35 * diag(nC)
        h_min_head <- max(0.10, min(eigen(H_head, symmetric = TRUE, only.values = TRUE)$values))
        D0_s <- S0_isqrt %*% Delta0_C %*% S0_isqrt
        Dk_s <- lapply(seq_len(M_dyn), function(j) S0_isqrt %*% B_gated_C[[j]] %*% S0_isqrt)
        if (nC == 2L) {
          a_h <- rep(H_head[1L, 1L] + D0_s[1L, 1L], nrow(all_f_frames))
          d_h <- rep(H_head[2L, 2L] + D0_s[2L, 2L], nrow(all_f_frames))
          b_h <- rep(H_head[1L, 2L] + D0_s[1L, 2L], nrow(all_f_frames))
          a_v <- rep(D0_s[1L, 1L], nrow(all_f_frames))
          d_v <- rep(D0_s[2L, 2L], nrow(all_f_frames))
          b_v <- rep(D0_s[1L, 2L], nrow(all_f_frames))
          for (j in seq_len(M_dyn)) {
            fj <- all_f_frames[, j]; Mj <- Dk_s[[j]]
            a_h <- a_h + fj * Mj[1L, 1L]; d_h <- d_h + fj * Mj[2L, 2L]; b_h <- b_h + fj * Mj[1L, 2L]
            a_v <- a_v + fj * Mj[1L, 1L]; d_v <- d_v + fj * Mj[2L, 2L]; b_v <- b_v + fj * Mj[1L, 2L]
          }
          lam_tot_min <- min(0.5 * ((a_h + d_h) - sqrt((a_h - d_h)^2 + 4.0 * (b_h^2))))
          lam_dyn_min <- min(0.5 * ((a_v + d_v) - sqrt((a_v - d_v)^2 + 4.0 * (b_v^2))))
        } else {
          chk_idx <- unique(as.integer(unlist(lapply(seq_len(M_dyn), function(j) {
            ord_j <- order(all_f_frames[, j])
            c(utils::head(ord_j, 4L), utils::tail(ord_j, 4L))
          }))))
          lam_tot_min <- 0.0; lam_dyn_min <- 0.0
          for (ti_f in chk_idx) {
            Mt <- D0_s
            for (j in seq_len(M_dyn)) Mt <- Mt + all_f_frames[ti_f, j] * Dk_s[[j]]
            ev_v <- eigen(0.5 * (Mt + t(Mt)), symmetric = TRUE, only.values = TRUE)$values
            ev_h <- eigen(0.5 * ((H_head + Mt) + t(H_head + Mt)), symmetric = TRUE, only.values = TRUE)$values
            lam_dyn_min <- min(lam_dyn_min, min(ev_v))
            lam_tot_min <- min(lam_tot_min, min(ev_h))
          }
        }
        if (is.finite(lam_tot_min) && lam_tot_min < 0.0 && is.finite(lam_dyn_min) && lam_dyn_min < -1e-8) {
          sc_cone <- min(1.0, h_min_head / abs(lam_dyn_min))
          Delta0_C <- sc_cone * Delta0_C
          for (j in seq_len(M_dyn)) B_gated_C[[j]] <- sc_cone * B_gated_C[[j]]
        }
      }

      for (ia in seq_len(nC)) {
        g <- C_u[ia]
        b_gg <- vapply(seq_len(M_dyn), function(j) B_gated_C[[j]][ia, ia], numeric(1L))
        rsd_gg <- matrix(R_sd_C[, , ia, ia], nrow = L_act, ncol = M_dyn)
        Gamma_cond_diag[, g] <- gam_unc_C[, ia, ia] - as.numeric(rsd_gg %*% b_gg)
        for (ib in seq_len(nC)) {
          h <- C_u[ib]
          if (!supp_C[ia, ib]) next
          Delta0_1lag_all[g, h] <- Delta0_C[ia, ib]
          for (j in seq_len(M_dyn)) B_dyn_gated[j, g, h] <- B_gated_C[[j]][ia, ib]
        }
      }
    }

    Gamma_stat_gated <- array(0.0, dim = c(L1, P, P))
    Sig0_1lag_cond <- R_poet * outer(sig0_sd, sig0_sd) + Delta0_1lag_all
    diag(Sig0_1lag_cond) <- pmax(0.15 * (sig0_rob^2), diag(Sig0_1lag_cond))

    if (!use_mlag) {
      Gamma_stat_gated[1L, , ] <- Sig0_1lag_cond
    } else {
      N_freq <- 256L
      omega_grid <- 2.0 * pi * (0:(N_freq - 1L)) / N_freq
      cos_grid <- outer(omega_grid, 0:L_max, function(w, d) cos(w * d))
      w_hann <- 0.5 * (1.0 + cos(pi * seq_len(L_max) / (L_max + 1.0)))
      wt_lag <- c(1.0, 2.0 * w_hann)

      sig0_cond_mlag <- sqrt(pmax(Gamma_cond_diag[1L, ], 0.15 * (sig0_rob^2)))
      acf_c_mat <- matrix(0.0, L_max, P)
      Q_port_c <- numeric(P); Q_lf_c <- numeric(P)
      for (g in seq_len(P)) {
        s2_g <- sig0_cond_mlag[g]^2
        acf_c_mat[, g] <- pmin(0.95, pmax(-0.95, Gamma_cond_diag[2:L1, g] / s2_g))
        XtX_s <- XtX_diag_sig[[sig_idx[g]]][1:L1, 1:L1, drop = FALSE]
        V_inv_sub <- XtX_s[2:L1, 2:L1, drop = FALSE] - XtX_s[2:L1, 1L, drop = FALSE] %*% solve(XtX_s[1L, 1L, drop = FALSE], XtX_s[1L, 2:L1, drop = FALSE])
        gam_cg <- s2_g * acf_c_mat[, g]
        Q_port_c[g] <- as.numeric(crossprod(gam_cg, V_inv_sub %*% gam_cg)) / (2.0 * (s2_g^2))
        V_cov_w <- as.numeric(crossprod(w_bart_lf, solve(V_inv_sub + 1e-6 * diag(L_max), w_bart_lf)))
        Q_lf_c[g] <- (sum(w_bart_lf * gam_cg)^2) / max(2.0 * (s2_g^2) * V_cov_w, 1e-12)
      }
      acf_c_mean <- rowMeans(acf_c_mat)
      sqrt_spec_unit <- matrix(1.0, N_freq, P)
      Gamma_mlag <- array(0.0, dim = c(L1, P, P))

      for (g in seq_len(P)) {
        q_snr_g <- max(Q_port_c[g] / L_max, Q_lf_c[g], 1e-6)
        w_dev_g <- max(0.0, 1.0 - 1.35 / q_snr_g)
        acf_g <- acf_c_mean + w_dev_g * (acf_c_mat[, g] - acf_c_mean)
        sp_g <- pmax(as.numeric(cos_grid %*% (wt_lag * c(1.0, acf_g))), 0.008)
        sp_g <- sp_g / mean(sp_g)
        sqrt_spec_unit[, g] <- sqrt(sp_g)
        r_g_psd <- as.numeric(crossprod(cos_grid, sp_g)) / N_freq
        Gamma_mlag[, g, g] <- (sig0_cond_mlag[g]^2) * r_g_psd
      }

      if (nrow(active_off) > 0L) {
        for (kk in seq_len(nrow(active_off))) {
          g <- active_off[kk, 1L]; h <- active_off[kk, 2L]
          cross_sp <- sqrt_spec_unit[, g] * sqrt_spec_unit[, h]
          r_gh_psd <- as.numeric(crossprod(cos_grid, cross_sp)) / N_freq
          gam_gh <- R_poet[g, h] * sig0_cond_mlag[g] * sig0_cond_mlag[h] * r_gh_psd
          Gamma_mlag[, g, h] <- gam_gh
          Gamma_mlag[, h, g] <- gam_gh
        }
      }

      Gamma_stat_gated <- w_res_global * Gamma_mlag
      Gamma_stat_gated[1L, , ] <- Gamma_stat_gated[1L, , ] + (1.0 - w_res_global) * Sig0_1lag_cond
    }

    shared_error_model <- list(
      version = 3L,
      mode = "exact",
      unit_ids = unit_ids,
      lags = lags,
      run_ids = run_ids_active,
      K_reg = K_reg,
      M_dyn = M_dyn,
      L_max = L_max,
      j_rwls = j_rwls,
      vc_col_idx = vc_col_idx,
      gamma_obs = gamma_obs,
      w_res_global = w_res_global,
      w_ltz_global = w_ltz_global,
      R_poet = R_poet,
      poet_K = poet_fit$K,
      sig0_sd = sig0_sd,
      sig0_rob = sig0_rob,
      Gamma_stat_gated = Gamma_stat_gated,
      B_dyn_gated = B_dyn_gated,
      runs = lapply(seq_len(n_runs_act), function(r) {
        rs <- run_store[[r]]
        list(
          m = rs$m, rows = rs$rows, m_target = rs$m_target, s_idx = s_mat[, r],
          U_strata = rs$U_strata, V_obs_strata = rs$V_obs_strata, ok_strata = rs$ok_strata,
          e_task = rs$e_task, E_cov_env = rs$E_cov_env, lag_pairs = rs$lag_pairs,
          run_scale = rs$run_scale
        )
      })
    )

    # Stage E: Exact trial variances, length-4 cross-trial correlations, lag_noise_vec, and Block-PAF shared_noise_vec
    n_comp <- length(comp_list)
    if (n_comp == 0L) {
      D_emb_tot <- 1L
      comp_slices <- list()
    } else {
      d_want <- vapply(comp_list, function(Cu) min(as.integer(embed_rank), max(2L, length(Cu) * n_lags)), integer(1L))
      D_cap <- max(as.integer(embed_rank), min(64L, sum(d_want)))
      d_alloc <- if (sum(d_want) <= D_cap) d_want else pmax(2L, as.integer(floor(d_want * (D_cap / sum(d_want)))))
      D_emb_tot <- sum(d_alloc)
      ends_c <- cumsum(d_alloc)
      starts_c <- c(1L, utils::head(ends_c, -1L) + 1L)
      comp_slices <- Map(seq.int, starts_c, ends_c)
    }
    D_lag <- min(6L, max(1L, n_lags))
    K_tcor <- 4L

    run_panel_blocks <- vector("list", n_runs_act)
    for (r in seq_len(n_runs_act)) {
      rs <- run_store[[r]]
      rs_sem <- shared_error_model$runs[[r]]
      m_t <- rs$m_target

      center_op_list <- function(op_strata_list) {
        lapply(seq_len(rs$S_cnt), function(si) {
          lapply(seq_len(n_lags), function(li) {
            ok_t <- rs$ok_strata[[si]][[li]]
            U_out <- matrix(0.0, rs$m, m_t)
            if (!is.null(ok_t) && sum(ok_t) > 1L) {
              m_ok <- sum(ok_t)
              Q_t <- diag(m_ok) - matrix(1.0 / m_ok, m_ok, m_ok)
              U_out[, ok_t] <- (op_strata_list[[si]][[li]][, ok_t, drop = FALSE] %*% Q_t) / sqrt(1.0 - 1.0 / m_ok)
            } else if (!is.null(ok_t) && sum(ok_t) == 1L) {
              U_out[, ok_t] <- op_strata_list[[si]][[li]][, ok_t, drop = FALSE]
            }
            U_out
          })
        })
      }

      U_cent_s <- center_op_list(rs$U_strata)
      V_cent_s <- center_op_list(rs$V_obs_strata)

      se_arr_r <- array(NA_real_, dim = c(m_t, P, n_lags))
      s2_trial_arr <- array(NA_real_, dim = c(m_t, P, n_lags))
      tcor_arr_r <- if (isTRUE(return_panel)) array(0.0, dim = c(K_tcor, m_t, P, n_lags)) else NULL
      lvec_arr_r <- if (isTRUE(return_panel) && n_lags >= 2L) array(0.0, dim = c(D_lag, m_t, P, n_lags)) else NULL
      svec_arr_r <- if (isTRUE(return_panel)) array(0.0, dim = c(D_emb_tot, m_t, P, n_lags)) else NULL

      for (g in seq_len(P)) {
        si <- rs_sem$s_idx[g]
        for (li in seq_len(n_lags)) {
          ok_t <- rs$ok_strata[[si]][[li]]
          if (is.null(ok_t) || !any(ok_t)) next
          Uc <- U_cent_s[[si]][[li]][, ok_t, drop = FALSE]
          Vc <- V_cent_s[[si]][[li]][, ok_t, drop = FALSE]
          v_base <- pmax(1e-10, medusa_eval_sem_bilinear(Uc, Uc, Vc, Vc, rs_sem, shared_error_model, g, g, stat_floor = TRUE))
          v_dyn <- pmax(0.15 * v_base, medusa_eval_sem_bilinear(Uc, Uc, Vc, Vc, rs_sem, shared_error_model, g, g, stat_floor = FALSE))
          s2_full <- rep(NA_real_, m_t)
          s2_full[ok_t] <- v_dyn
          s2_trial_arr[, g, li] <- s2_full
          se_arr_r[, g, li] <- sqrt(s2_full)

          if (isTRUE(return_panel) && m_t >= 2L) {
            Uu_full <- rs$U_strata[[si]][[li]]
            Vu_full <- rs$V_obs_strata[[si]][[li]]
            vu_base <- pmax(1e-12, medusa_eval_sem_bilinear(Uu_full, Uu_full, Vu_full, Vu_full, rs_sem, shared_error_model, g, g, stat_floor = TRUE))
            vu_dyn <- pmax(0.15 * vu_base, medusa_eval_sem_bilinear(Uu_full, Uu_full, Vu_full, Vu_full, rs_sem, shared_error_model, g, g, stat_floor = FALSE))
            for (kd in seq_len(min(K_tcor, m_t - 1L))) {
              i_a <- seq_len(m_t - kd); i_b <- i_a + kd
              ok_ab <- ok_t[i_a] & ok_t[i_b]
              if (any(ok_ab)) {
                c_kd <- medusa_eval_sem_bilinear(
                  Uu_full[, i_a[ok_ab], drop = FALSE], Uu_full[, i_b[ok_ab], drop = FALSE],
                  Vu_full[, i_a[ok_ab], drop = FALSE], Vu_full[, i_b[ok_ab], drop = FALSE],
                  rs_sem, shared_error_model, g, g, stat_floor = FALSE
                )
                den_kd <- sqrt(pmax(1e-20, vu_dyn[i_a[ok_ab]] * vu_dyn[i_b[ok_ab]]))
                tcor_arr_r[kd, i_a[ok_ab], g, li] <- pmin(0.95, pmax(-0.95, c_kd / den_kd))
              }
            }
          }
        }

        if (isTRUE(return_panel) && n_lags >= 2L) {
          Cg_tr <- array(0.0, dim = c(n_lags, n_lags, m_t))
          for (la in seq_len(n_lags)) {
            Cg_tr[la, la, ] <- ifelse(is.finite(s2_trial_arr[, g, la]), s2_trial_arr[, g, la], 0.0)
            if (la < n_lags) {
              Ua <- U_cent_s[[si]][[la]]; Va <- V_cent_s[[si]][[la]]
              for (lb in seq.int(la + 1L, n_lags)) {
                Ub <- U_cent_s[[si]][[lb]]; Vb <- V_cent_s[[si]][[lb]]
                c_ab <- medusa_eval_sem_bilinear(Ua, Ub, Va, Vb, rs_sem, shared_error_model, g, g, stat_floor = FALSE)
                Cg_tr[la, lb, ] <- c_ab; Cg_tr[lb, la, ] <- c_ab
              }
            }
          }
          for (ti in seq_len(m_t)) {
            Ci <- .project_trial_err_psd(Cg_tr[, , ti])
            eg_l <- eigen(Ci, symmetric = TRUE)
            k_l <- min(D_lag, sum(eg_l$values > 1e-12))
            if (k_l >= 1L) {
              Fl <- sweep(eg_l$vectors[, seq_len(k_l), drop = FALSE], 2L, sqrt(eg_l$values[seq_len(k_l)]), "*")
              lvec_arr_r[seq_len(k_l), ti, g, ] <- t(Fl)
            }
          }
        }
      }

      if (isTRUE(return_panel) && n_comp >= 1L) {
        f_dyn_r <- cbind(rs$e_task, rs$E_cov_env)
        ones_lag <- matrix(1.0, n_lags, n_lags)

        for (cm in seq_len(n_comp)) {
          C_u <- comp_list[[cm]]
          nC <- length(C_u)
          q_m <- nC * n_lags
          sl_idx <- comp_slices[[cm]]
          d_m <- length(sl_idx)

          s_vec_C <- rs_sem$s_idx[C_u]
          u_s_C <- unique(s_vec_C)
          n_sC <- length(u_s_C)
          q_sC <- n_sC * n_lags
          map_al_to_sl <- integer(q_m)
          for (la in seq_len(n_lags)) {
            map_al_to_sl[(la - 1L) * nC + seq_len(nC)] <- (la - 1L) * n_sC + match(s_vec_C, u_s_C)
          }

          probe_m <- if (q_m > 512L) {
            r_sub <- min(q_m, d_m + 12L)
            qr.Q(qr(outer(seq_len(q_m), seq_len(r_sub), function(i, j) sin((i * 1.6180339887 + j) * j))))
          } else NULL
          G_stat_C <- lapply(seq_len(L1), function(idx_l) {
            kronecker(ones_lag, matrix(Gamma_stat_gated[idx_l, C_u, C_u], nC, nC))
          })
          B_dyn_C <- lapply(seq_len(M_dyn), function(j) {
            kronecker(ones_lag, matrix(B_dyn_gated[j, C_u, C_u], nC, nC))
          })

          for (ti in seq_len(m_t)) {
            U_ti <- matrix(0.0, rs$m, q_sC)
            V_ti <- matrix(0.0, rs$m, q_sC)
            for (la in seq_len(n_lags)) {
              for (isi in seq_len(n_sC)) {
                col_k <- (la - 1L) * n_sC + isi
                si_orig <- u_s_C[isi]
                U_ti[, col_k] <- U_cent_s[[si_orig]][[la]][, ti]
                V_ti[, col_k] <- V_cent_s[[si_orig]][[la]][, ti]
              }
            }

            T0_s <- crossprod(U_ti)
            Ci <- G_stat_C[[1L]] * T0_s[map_al_to_sl, map_al_to_sl, drop = FALSE]

            for (d in seq_len(L_max)) {
              Gd_mat <- G_stat_C[[d + 1L]]
              if (any(abs(Gd_mat) > 1e-14)) {
                pd <- rs$lag_pairs[[d]]
                if (nrow(pd) > 0L) {
                  Td_s <- crossprod(U_ti[pd[, 1L], , drop = FALSE], U_ti[pd[, 2L], , drop = FALSE])
                  Ci <- Ci + Gd_mat * Td_s[map_al_to_sl, map_al_to_sl, drop = FALSE]
                }
              }
            }
            for (j in seq_len(M_dyn)) {
              Bj_mat <- B_dyn_C[[j]]
              if (any(abs(Bj_mat) > 1e-14)) {
                fj <- f_dyn_r[, j]
                gm_j <- if (length(j_rwls) == 1L && j == j_rwls) 1.0 else gamma_obs
                D_comb <- matrix(0.0, q_sC, q_sC)
                if (gm_j < 1.0) D_comb <- D_comb + (1.0 - gm_j) * crossprod(U_ti, fj * U_ti)
                if (gm_j > 0.0) D_comb <- D_comb + gm_j * crossprod(V_ti, fj * V_ti)
                Ci <- Ci + Bj_mat * D_comb[map_al_to_sl, map_al_to_sl, drop = FALSE]
              }
            }

            if (!is.null(rs_sem$run_scale)) {
              sc_Cu <- rep(rs_sem$run_scale[C_u], times = n_lags)
              Ci <- Ci * tcrossprod(sc_Cu)
            }
            for (la in seq_len(n_lags)) {
              idx_diag_la <- (la - 1L) * nC + seq_len(nC)
              s2_d <- as.numeric(s2_trial_arr[ti, C_u, la])
              diag(Ci)[idx_diag_la] <- ifelse(is.finite(s2_d), s2_d, 0.0)
            }

            F_ti <- .extract_psd_factors(Ci, d_m, nC = nC, n_lags = n_lags, probe_cache = probe_m)
            for (la in seq_len(n_lags)) {
              cols_la <- (la - 1L) * nC + seq_len(nC)
              svec_arr_r[sl_idx, ti, C_u, la] <- F_ti[, cols_la, drop = FALSE]
            }
          }
        }
      }

      run_panel_blocks[[r]] <- list(
        orig_r = keep_runs[r],
        se = se_arr_r,
        tcor = tcor_arr_r,
        lvec = lvec_arr_r,
        svec = svec_arr_r
      )
    }

    return(list(
      unit_error_cor = R_poet,
      shared_loading_vec = shared_loading_vec,
      shared_error_model = shared_error_model,
      run_panel_blocks = run_panel_blocks
    ))
  }

  # =========================================================================
  # PATH 2: VOXEL-SAFE DUAL-SPACE LOW-RANK ENGINE (P > max_exact_units; zero V x V!)
  # =========================================================================
  df_v_tot <- numeric(P)
  s_num_v <- numeric(P)
  s_rob_v <- numeric(P)
  s_blk_tot <- numeric(P)
  any_rwls_vox <- FALSE

  for (r in seq_len(n_runs_act)) {
    rs <- run_store[[r]]
    df_r <- pmax(1.0, rs$df_unit_v)
    df_v_tot <- df_v_tot + df_r
    s_num_v <- s_num_v + rs$s2_raw_u * df_r
    s_rob_v <- s_rob_v + rs$s2_rob_u * df_r
    if (any(rs$s_blk_u > 0)) {
      any_rwls_vox <- TRUE
      s_blk_tot <- s_blk_tot + rs$s_blk_u
    }
  }
  df_v_tot <- pmax(df_v_tot, 1e-8)
  sig0_raw <- sqrt(pmax(s_num_v / df_v_tot, 1e-10))
  sig0_rob <- pmin(sig0_raw, sqrt(pmax(s_rob_v / df_v_tot, 1e-10)))
  sig0_sd <- if (any_rwls_vox) {
    sqrt(pmax(sig0_rob^2, (s_num_v - s_blk_tot) / df_v_tot))
  } else {
    sig0_raw
  }

  # Dual-space SVD over stacked sketch rows (sum_r K_sk_r x P, <= 64*n_runs x P; NEVER P x P!)
  B_rob_stack <- do.call(rbind, lapply(run_store, function(rs) {
    if (!is.null(rs$vox_sketch)) rs$vox_sketch$B_rob_sk else matrix(0.0, 0L, P)
  }))
  col_nrm_rob <- sqrt(pmax(s_rob_v, colSums(B_rob_stack^2), 1e-16))
  Z_rob_low <- if (nrow(B_rob_stack) > 0L) sweep(B_rob_stack, 2L, col_nrm_rob, "/") else matrix(0.0, 0L, P)

  L_poet <- matrix(0.0, P, 0L, dimnames = list(unit_ids, NULL))
  K_poet <- 0L
  if (nrow(Z_rob_low) >= 2L) {
    G_dual <- tcrossprod(Z_rob_low)
    eg_d <- eigen(0.5 * (G_dual + t(G_dual)), symmetric = TRUE)
    pos_d <- which(eg_d$values > 1e-8)
    if (length(pos_d) >= 1L) {
      vals_d <- eg_d$values[pos_d]
      u_d <- eg_d$vectors[, pos_d, drop = FALSE]
      v_d <- crossprod(Z_rob_low, sweep(u_d, 2L, 1.0 / sqrt(vals_d), "*"))
      sp_v <- poet_shrink_spikes(vals_d, v_d, P = P, df_eff = max(20.0, df_u_tot), k_max = k_max_poet, poet_mode = "adaptive_logp")
      K_poet <- sp_v$K
      L_poet <- sp_v$L
      rownames(L_poet) <- unit_ids
    }
  }
  shared_loading_vec <- if (K_poet >= 1L) {
    lv1 <- L_poet[, 1L]
    if (sum(lv1) < 0) lv1 <- -lv1
    stats::setNames(pmin(0.995, pmax(-0.995, lv1)), unit_ids)
  } else {
    stats::setNames(numeric(P), unit_ids)
  }

  if (!isTRUE(return_panel) && !isTRUE(return_epochs)) {
    return(list(unit_error_cor = NULL, shared_loading_vec = shared_loading_vec, shared_error_model = NULL))
  }

  M_dyn <- 1L + ncol(run_store[[1L]]$E_cov_env)
  K_reg <- 1L + M_dyn
  B_tot <- L1 + 2L * M_dyn

  s_mat <- matrix(
    vapply(run_store, function(rs) ifelse(is.na(rs$s_idx) | rs$s_idx < 1L, 1L, rs$s_idx), integer(P)),
    nrow = P, ncol = n_runs_act
  )
  unit_sig <- apply(s_mat, 1L, paste, collapse = ":")
  u_sigs <- unique(unit_sig)
  G_sig <- length(u_sigs)
  sig_idx <- match(unit_sig, u_sigs)
  rep_unit_for_sig <- match(u_sigs, unit_sig)

  X_diag_cache <- lapply(seq_len(n_runs_act), function(r) vector("list", G_sig))
  XtX_diag_sig <- vector("list", G_sig)
  for (a in seq_len(G_sig)) {
    ua <- rep_unit_for_sig[a]
    XtX_aa <- matrix(0.0, B_tot, B_tot)
    for (r in seq_len(n_runs_act)) {
      rs <- run_store[[r]]; sa <- s_mat[ua, r]
      f_dyn <- cbind(rs$e_task, rs$E_cov_env)
      if (rs$saturated) f_dyn[] <- 0.0
      X_ra <- poet_build_pair_X_cpp(
        rs$M_strata[[sa]], rs$M_strata[[sa]],
        rs$G_inn_strata[[sa]], rs$G_inn_strata[[sa]],
        rs$H_strata[[sa]], rs$H_strata[[sa]],
        rs$adj_idx1, rs$lag_pairs, f_dyn, L_max
      )
      X_diag_cache[[r]][[a]] <- X_ra
      XtX_aa <- XtX_aa + crossprod(X_ra)
    }
    XtX_diag_sig[[a]] <- XtX_aa
  }

  # Compute Xty_diag (B_tot x P): exact from E_vox_full when available, or fast O(P K_sk^2) from sketch!
  has_full_E_vox <- all(vapply(run_store, function(rs) !is.null(rs$E_vox_full), logical(1L)))
  if (has_full_E_vox) {
    Xty_diag <- poet_extract_Xty_diag_cpp(
      lapply(run_store, function(rs) rs$E_vox_full), X_diag_cache, as.integer(sig_idx),
      lapply(run_store, function(rs) as.integer(rs$adj_idx1)),
      lapply(run_store, function(rs) rs$lag_pairs), L_max, B_tot
    )
  } else {
    Xty_diag <- matrix(0.0, B_tot, P)
    for (r in seq_len(n_runs_act)) {
      rs <- run_store[[r]]; m <- rs$m; idx1 <- rs$adj_idx1
      sk_r <- rs$vox_sketch
      if (is.null(sk_r)) next
      U_sk <- sk_r$U_sk; B_sk <- sk_r$B_sk
      # Unexplained idiosyncratic lag-0 variance adjustment per voxel
      df_r <- pmax(1.0, rs$df_unit_v)
      s2_unexp_r <- pmax(0.0, rs$s2_raw_u - colSums(B_sk^2) / df_r)
      for (a in seq_len(G_sig)) {
        cols_a <- which(sig_idx == a)
        if (!length(cols_a)) next
        Xr <- X_diag_cache[[r]][[a]]
        B_sub <- B_sk[, cols_a, drop = FALSE]
        for (ell in seq_len(B_tot)) {
          w0 <- Xr[1:m, ell]
          K_ell <- crossprod(U_sk, w0 * U_sk)
          if (length(idx1)) {
            w1 <- Xr[m + idx1, ell] / sqrt(2.0)
            K1 <- crossprod(U_sk[idx1 + 1L, , drop = FALSE], w1 * U_sk[idx1, , drop = FALSE])
            K_ell <- K_ell + K1 + t(K1)
          }
          if (L_max >= 2L) {
            for (d in 2:L_max) {
              wd <- Xr[2L * m - 1L + d - 1L, ell]
              if (abs(wd) > 1e-14) {
                pd <- rs$lag_pairs[[d]]; N_d <- 0.5 * nrow(pd)
                if (N_d > 0) {
                  K_ell <- K_ell + (wd / sqrt(2.0 * N_d)) * crossprod(U_sk[pd[, 1L], , drop = FALSE], U_sk[pd[, 2L], , drop = FALSE])
                }
              }
            }
          }
          Xty_diag[ell, cols_a] <- Xty_diag[ell, cols_a] + colSums(B_sub * (K_ell %*% B_sub)) + sum(w0) * s2_unexp_r[cols_a]
        }
      }
    }
  }
  # Free full m x P voxel residual matrices if allocated
  for (r in seq_len(n_runs_act)) run_store[[r]]$E_vox_full <- NULL

  # Estimate gamma_obs and Portmanteau w_res_global across voxels
  gamma_obs <- 0.80
  if (!any_saturated) {
    keep_gam_dyn <- seq_len(min(M_dyn, 1L + q_cov))
    n_gam_dyn <- length(keep_gam_dyn)
    grid_gam <- seq(0.0, 1.0, by = 0.05)
    prior_wt <- 0.0025 * sum(sig0_sd^4)
    obj_gam <- vapply(grid_gam, function(gm) {
      A_gm <- matrix(0.0, B_tot, L1 + n_gam_dyn)
      A_gm[1:L1, 1:L1] <- diag(L1)
      for (ii in seq_len(n_gam_dyn)) {
        j <- keep_gam_dyn[ii]
        A_gm[L1 + j, L1 + ii] <- 1.0 - gm
        A_gm[L1 + M_dyn + j, L1 + ii] <- gm
      }
      ssr_val <- 0.0
      for (a in seq_len(G_sig)) {
        cols_a <- which(sig_idx == a)
        if (!length(cols_a)) next
        XtX_r <- crossprod(A_gm, XtX_diag_sig[[a]] %*% A_gm)
        act_c <- which(diag(XtX_r) > 1e-8)
        if (length(act_c)) {
          rhs_a <- crossprod(A_gm, Xty_diag[, cols_a, drop = FALSE])[act_c, , drop = FALSE]
          b_a <- solve(XtX_r[act_c, act_c, drop = FALSE] + 1e-6 * diag(length(act_c)), rhs_a)
          ssr_val <- ssr_val + sum(rhs_a * b_a)
        }
      }
      ssr_val - prior_wt * ((gm - 0.80)^2)
    }, numeric(1L))
    if (any(is.finite(obj_gam))) {
      gamma_obs <- grid_gam[which.max(ifelse(is.finite(obj_gam), obj_gam, -Inf))[1L]]
    }
  }

  A_tie_full <- matrix(0.0, B_tot, L1 + M_dyn)
  A_tie_full[1:L1, 1:L1] <- diag(L1)
  for (j in seq_len(M_dyn)) {
    A_tie_full[L1 + j, L1 + j] <- 1.0 - gamma_obs
    A_tie_full[L1 + M_dyn + j, L1 + j] <- gamma_obs
  }
  A_tie_1lag <- matrix(0.0, B_tot, 1L + M_dyn)
  A_tie_1lag[1L, 1L] <- 1.0
  for (j in seq_len(M_dyn)) {
    A_tie_1lag[L1 + j, 1L + j] <- 1.0 - gamma_obs
    A_tie_1lag[L1 + M_dyn + j, 1L + j] <- gamma_obs
  }

  w_bart_lf <- numeric(L_max)
  if (L_max >= 2L) {
    idx_tl <- 2:min(6L, L_max)
    w_bart_lf[idx_tl] <- 1.0 - (seq_along(idx_tl) - 1.0) / length(idx_tl)
  } else {
    w_bart_lf[1L] <- 1.0
  }

  Q_port <- numeric(P); z_tail_g <- numeric(P)
  sig0_mlag_diag <- numeric(P)
  gam_u_all <- matrix(0.0, L1, P)
  for (a in seq_len(G_sig)) {
    cols_a <- which(sig_idx == a)
    if (!length(cols_a)) next
    XtX_gg <- XtX_diag_sig[[a]]
    XtX_tf <- crossprod(A_tie_full, XtX_gg %*% A_tie_full)
    act_f <- which(diag(XtX_tf) > 1e-8)
    b_tf <- matrix(0.0, L1 + M_dyn, length(cols_a))
    if (length(act_f)) {
      rhs_tf <- crossprod(A_tie_full, Xty_diag[, cols_a, drop = FALSE])[act_f, , drop = FALSE]
      b_tf[act_f, ] <- solve(XtX_tf[act_f, act_f, drop = FALSE] + 1e-6 * diag(length(act_f)), rhs_tf)
    }
    XtX_s <- XtX_gg[1:L1, 1:L1, drop = FALSE]
    gam_u <- solve(XtX_s + 1e-6 * diag(L1), Xty_diag[1:L1, cols_a, drop = FALSE])
    gam_u_all[, cols_a] <- gam_u
    sig0_mlag_diag[cols_a] <- sqrt(pmax(gam_u[1L, ], 0.15 * (sig0_sd[cols_a]^2)))

    s2_cg <- pmax(sig0_sd[cols_a]^2, 0.25 * (sig0_rob[cols_a]^2), 1e-10)
    acf_raw <- matrix(
      pmin(0.95, pmax(-0.95, sweep(b_tf[2:L1, , drop = FALSE], 2L, s2_cg, "/"))),
      nrow = L_max, ncol = length(cols_a)
    )
    V_inv_sub <- XtX_s[2:L1, 2:L1, drop = FALSE] - XtX_s[2:L1, 1L, drop = FALSE] %*% solve(XtX_s[1L, 1L, drop = FALSE], XtX_s[1L, 2:L1, drop = FALSE])
    V_sub_cov_w <- as.numeric(crossprod(w_bart_lf, solve(V_inv_sub + 1e-6 * diag(L_max), w_bart_lf)))
    Q_port[cols_a] <- 0.5 * colSums(acf_raw * (V_inv_sub %*% acf_raw))
    z_tail_g[cols_a] <- as.numeric(crossprod(w_bart_lf, acf_raw)) / sqrt(max(2.0 * V_sub_cov_w, 1e-12))
  }

  r2_bar_v <- if (K_poet >= 1L) {
    G_L <- crossprod(L_poet)
    max(0.0, (sum(G_L^2) - sum(rowSums(L_poet^2)^2)) / max(1.0, as.numeric(P) * (P - 1.0)))
  } else 0.0
  p_eff_port <- P / (1.0 + (P - 1.0) * r2_bar_v)
  z_omni_glob <- (mean(Q_port) - L_max) / max(sqrt(2.0 * L_max / max(1.0, p_eff_port)), 1e-6)
  z_signed_glob <- sqrt(max(1.0, p_eff_port)) * mean(z_tail_g)
  w_res_global <- .poet_firm_gate(max(z_omni_glob, abs(z_signed_glob) - 0.15), lo = 1.55, hi = 2.75)

  use_mlag <- (w_res_global > 1e-6)
  A_tie_use <- if (use_mlag) A_tie_full else A_tie_1lag
  L_act <- if (use_mlag) L1 else 1L
  sig0_use <- if (use_mlag) sqrt((1.0 - w_res_global) * (sig0_sd^2) + w_res_global * (sig0_mlag_diag^2)) else sig0_sd

  # Vectorized diagonal dynamic slope fit + LTZ gating across all P voxels
  b_dyn_raw_v <- matrix(0.0, M_dyn, P)
  z_dyn_v <- matrix(0.0, M_dyn, P)
  r_sd1_v <- matrix(0.0, M_dyn, P)
  R_sd_sig <- vector("list", G_sig)
  M_core_c <- min(M_dyn, 1L + q_cov)
  n_vc_c <- length(vc_col_idx)
  b_core_al_v <- matrix(0.0, M_core_c, P)
  z_core_al_v <- matrix(0.0, M_core_c, P)
  R_cv_sig <- vector("list", G_sig)

  for (a in seq_len(G_sig)) {
    cols_a <- which(sig_idx == a)
    if (!length(cols_a)) next
    XtX_t <- crossprod(A_tie_use, XtX_diag_sig[[a]] %*% A_tie_use)
    n_p <- L_act + M_dyn
    act_k <- which(diag(XtX_t) > 1e-8)
    if (!length(act_k)) next
    XtX_t_inv <- matrix(0.0, n_p, n_p)
    XtX_t_inv[act_k, act_k] <- solve(XtX_t[act_k, act_k, drop = FALSE] + 1e-8 * diag(length(act_k)))
    rhs_t <- crossprod(A_tie_use, Xty_diag[, cols_a, drop = FALSE])
    b_t <- XtX_t_inv %*% rhs_t
    b_dyn_raw_v[, cols_a] <- b_t[(L_act + 1L):n_p, , drop = FALSE]
    R_sd_a <- solve(XtX_t[1:L_act, 1:L_act, drop = FALSE] + 1e-8 * diag(L_act), XtX_t[1:L_act, (L_act + 1L):n_p, drop = FALSE])
    R_sd_sig[[a]] <- R_sd_a
    r_sd1_v[, cols_a] <- R_sd_a[1L, ]

    v_sc_cols <- 2.0 * (sig0_sd[cols_a]^4) * 1.25
    for (j in seq_len(M_dyn)) {
      vj <- v_sc_cols * XtX_t_inv[L_act + j, L_act + j]
      z_dyn_v[j, cols_a] <- ifelse(vj > 1e-16, b_dyn_raw_v[j, cols_a] / sqrt(vj), 0.0)
    }
    if (n_vc_c >= 1L) {
      idx_sc <- seq_len(L_act + M_core_c)
      idx_vc <- L_act + vc_col_idx
      inv_sc <- solve(XtX_t[idx_sc, idx_sc, drop = FALSE] + 1e-8 * diag(length(idx_sc)))
      R_cv_sig[[a]] <- inv_sc[(L_act + 1L):(L_act + M_core_c), , drop = FALSE] %*% XtX_t[idx_sc, idx_vc, drop = FALSE]
      b_sc <- inv_sc %*% rhs_t[idx_sc, , drop = FALSE]
      b_core_al_v[, cols_a] <- b_sc[(L_act + 1L):(L_act + M_core_c), , drop = FALSE]
      for (jc0 in seq_len(M_core_c)) {
        vj0 <- v_sc_cols * inv_sc[L_act + jc0, L_act + jc0]
        z_core_al_v[jc0, cols_a] <- ifelse(vj0 > 1e-16, b_core_al_v[jc0, cols_a] / sqrt(vj0), 0.0)
      }
    }
  }

  w_ltz_global <- numeric(M_dyn)
  B_dyn_diag <- matrix(0.0, M_dyn, P)
  if (!any_saturated) {
    w_vc_max <- 0.0
    if (n_vc_c >= 1L) {
      k_vc_mult <- max(1.0, as.numeric(n_vc_c))
      for (j in vc_col_idx) {
        is_pos_j <- j %in% vc_pos_idx
        w_j <- if (is_pos_j) {
          .ltz_from_factor_quad(pmax(0.0, z_dyn_v[j, ]), L_poet, p_lo = 0.12 / k_vc_mult, p_hi = 0.004 / k_vc_mult)
        } else {
          .ltz_from_factor_quad(z_dyn_v[j, ], L_poet, p_lo = 0.08 / k_vc_mult, p_hi = 0.002 / k_vc_mult)
        }
        w_ltz_global[j] <- w_j
        w_vc_max <- max(w_vc_max, w_j)
        B_dyn_diag[j, ] <- w_j * (if (is_pos_j) pmax(0.0, b_dyn_raw_v[j, ]) else b_dyn_raw_v[j, ])
      }
      for (a in seq_len(G_sig)) {
        cols_a <- which(sig_idx == a)
        if (!length(cols_a) || is.null(R_cv_sig[[a]])) next
        diff_vc <- b_dyn_raw_v[vc_col_idx, cols_a, drop = FALSE] - B_dyn_diag[vc_col_idx, cols_a, drop = FALSE]
        b_dyn_raw_v[1:M_core_c, cols_a] <- b_dyn_raw_v[1:M_core_c, cols_a, drop = FALSE] + R_cv_sig[[a]] %*% diff_vc
        z_dyn_v[1:M_core_c, cols_a] <- (1.0 - w_vc_max) * z_core_al_v[, cols_a, drop = FALSE] + w_vc_max * z_dyn_v[1:M_core_c, cols_a, drop = FALSE]
      }
    }
    w_1 <- .ltz_from_factor_quad(pmax(0.0, z_dyn_v[1L, ]), L_poet, p_lo = 0.25, p_hi = 0.01)
    w_ltz_global[1L] <- w_1
    B_dyn_diag[1L, ] <- w_1 * pmax(0.0, b_dyn_raw_v[1L, ])
    core_extra <- setdiff(seq_len(M_dyn), c(1L, vc_col_idx))
    for (j in core_extra) {
      w_j <- .ltz_from_factor_quad(z_dyn_v[j, ], L_poet, p_lo = 0.20, p_hi = 0.005)
      w_ltz_global[j] <- w_j
      B_dyn_diag[j, ] <- w_j * b_dyn_raw_v[j, ]
    }
  }

  # Clamp per-voxel dynamic trajectory so quiet baseline never drops below 0.92 * sig0_rob^2 or < 0.35 * sig0_use^2
  all_f_frames <- do.call(rbind, lapply(run_store, function(rs) cbind(rs$e_task, rs$E_cov_env)))
  delta0_v <- -colSums(r_sd1_v * B_dyn_diag)
  if (length(vc_pos_idx) >= 1L) {
    d_za <- -colSums(r_sd1_v[vc_pos_idx, , drop = FALSE] * B_dyn_diag[vc_pos_idx, , drop = FALSE])
    d_oth <- delta0_v - d_za
    d_za <- pmax(d_za, 0.92 * (sig0_rob^2) - (sig0_use^2))
    delta0_v <- d_za + d_oth
  }
  if (nrow(all_f_frames) > 0L && any(abs(B_dyn_diag) > 1e-12)) {
    v_traj <- outer(rep(1.0, nrow(all_f_frames)), delta0_v) + all_f_frames %*% B_dyn_diag
    v_min <- apply(v_traj, 2L, min)
    head_v <- pmax(0.10 * (sig0_rob^2), (sig0_use^2) - 0.35 * (sig0_rob^2))
    over_v <- which(v_min < -head_v & v_min < -1e-10)
    if (length(over_v)) {
      sc_v <- pmin(1.0, head_v[over_v] / abs(v_min[over_v]))
      delta0_v[over_v] <- delta0_v[over_v] * sc_v
      B_dyn_diag[, over_v] <- sweep(B_dyn_diag[, over_v, drop = FALSE], 2L, sc_v, "*")
    }
  }

  Gamma_stat_diag <- matrix(0.0, L1, P)
  s2_cond_0 <- pmax(0.15 * (sig0_rob^2), (sig0_use^2) + delta0_v)
  if (!use_mlag) {
    Gamma_stat_diag[1L, ] <- s2_cond_0
  } else {
    gam_cond_all <- gam_u_all
    for (a in seq_len(G_sig)) {
      cols_a <- which(sig_idx == a)
      if (!length(cols_a) || is.null(R_sd_sig[[a]])) next
      gam_cond_all[, cols_a] <- gam_u_all[, cols_a, drop = FALSE] - R_sd_sig[[a]] %*% B_dyn_diag[, cols_a, drop = FALSE]
    }
    gam_norm <- sweep(gam_cond_all, 2L, pmax(sig0_sd^2, 1e-10), "/")
    acf_pool <- pmin(0.95, pmax(-0.95, rowMeans(gam_norm[2:L1, , drop = FALSE]) / max(mean(gam_norm[1L, ]), 0.25)))
    N_freq <- 256L
    omega_grid <- 2.0 * pi * (0:(N_freq - 1L)) / N_freq
    cos_grid <- outer(omega_grid, 0:L_max, function(w, d) cos(w * d))
    w_hann <- 0.5 * (1.0 + cos(pi * seq_len(L_max) / (L_max + 1.0)))
    wt_lag <- c(1.0, 2.0 * w_hann)
    sp_pool <- pmax(as.numeric(cos_grid %*% (wt_lag * c(1.0, acf_pool))), 0.008)
    sp_pool <- sp_pool / mean(sp_pool)
    r_pool_psd <- as.numeric(crossprod(cos_grid, sp_pool)) / N_freq
    for (idx_l in seq_len(L1)) {
      Gamma_stat_diag[idx_l, ] <- w_res_global * s2_cond_0 * r_pool_psd[idx_l]
    }
    Gamma_stat_diag[1L, ] <- Gamma_stat_diag[1L, ] + (1.0 - w_res_global) * s2_cond_0
  }
  sig0_stat_v <- sqrt(pmax(Gamma_stat_diag[1L, ], 1e-10))

  # K_poet x K_poet factor-space cross-voxel dynamic & stationary covariance
  Gamma_stat_fac <- array(0.0, dim = c(L1, max(1L, K_poet), max(1L, K_poet)))
  B_dyn_fac <- array(0.0, dim = c(M_dyn, max(1L, K_poet), max(1L, K_poet)))
  if (K_poet >= 1L) {
    Gamma_stat_fac[1L, , ] <- diag(1.0, K_poet)
    if (use_mlag) {
      for (idx_l in 2:L1) Gamma_stat_fac[idx_l, , ] <- r_pool_psd[idx_l] * diag(1.0, K_poet)
    }
    # Project voxel dynamic slopes into factor space aligned with L_sc = L_poet * sig0_stat_v
    L_sc <- sweep(L_poet, 1L, sig0_stat_v, "*")
    LtL_inv <- solve(crossprod(L_sc) + 1e-8 * diag(K_poet))
    for (j in seq_len(M_dyn)) {
      if (w_ltz_global[j] > 1e-6) {
        # Estimate factor dynamic excess from sketch residuals
        B_j_acc <- matrix(0.0, K_poet, K_poet)
        den_j_acc <- 0.0
        for (r in seq_len(n_runs_act)) {
          rs <- run_store[[r]]
          sk_r <- rs$vox_sketch
          if (is.null(sk_r)) next
          fj <- cbind(rs$e_task, rs$E_cov_env)[, j]
          if (all(fj == 0.0)) next
          # Factor time series F_r = U_sk (B_sk L_sc LtL_inv): m x K_poet
          F_r <- sk_r$U_sk %*% (sk_r$B_sk %*% (L_sc %*% LtL_inv))
          B_j_acc <- B_j_acc + crossprod(F_r, fj * F_r)
          den_j_acc <- den_j_acc + sum(fj^2)
        }
        if (den_j_acc > 1e-8) {
          Bj_s <- 0.5 * (B_j_acc + t(B_j_acc)) / den_j_acc
          if (j == 1L || (j %in% vc_pos_idx)) {
            eg_f <- eigen(Bj_s, symmetric = TRUE)
            Bj_s <- eg_f$vectors %*% (pmax(eg_f$values, 0.0) * t(eg_f$vectors))
          }
          B_dyn_fac[j, , ] <- w_ltz_global[j] * Bj_s
        }
      }
    }
  }

  shared_error_model <- list(
    version = 3L,
    mode = "lowrank_voxel",
    unit_ids = unit_ids,
    lags = lags,
    run_ids = run_ids_active,
    K_reg = K_reg,
    M_dyn = M_dyn,
    L_max = L_max,
    j_rwls = j_rwls,
    vc_col_idx = vc_col_idx,
    gamma_obs = gamma_obs,
    w_res_global = w_res_global,
    w_ltz_global = w_ltz_global,
    L_poet = L_poet,
    poet_K = K_poet,
    df_v_tot = df_v_tot,
    sig0_sd = sig0_sd,
    sig0_rob = sig0_rob,
    Gamma_stat_diag = Gamma_stat_diag,
    Gamma_stat_fac = Gamma_stat_fac,
    B_dyn_diag = B_dyn_diag,
    B_dyn_fac = B_dyn_fac,
    runs = lapply(seq_len(n_runs_act), function(r) {
      rs <- run_store[[r]]
      list(
        m = rs$m, rows = rs$rows, m_target = rs$m_target, s_idx = s_mat[, r],
        U_strata = rs$U_strata, V_obs_strata = rs$V_obs_strata, ok_strata = rs$ok_strata,
        e_task = rs$e_task, E_cov_env = rs$E_cov_env, lag_pairs = rs$lag_pairs,
        vox_sketch = rs$vox_sketch, run_scale = rs$run_scale
      )
    })
  )

  # Vectorized Stage E over strata S_cnt << P (fast even at P = 50,000!)
  D_emb_tot <- max(1L, min(as.integer(embed_rank), max(1L, K_poet * n_lags)))
  D_lag <- min(6L, max(1L, n_lags))
  K_tcor <- 4L
  L_sc_global <- if (K_poet >= 1L) sweep(L_poet, 1L, sig0_stat_v, "*") else NULL

  run_panel_blocks <- vector("list", n_runs_act)
  for (r in seq_len(n_runs_act)) {
    rs <- run_store[[r]]
    rs_sem <- shared_error_model$runs[[r]]
    m_t <- rs$m_target
    f_dyn_r <- cbind(rs$e_task, rs$E_cov_env)

    center_op_list <- function(op_strata_list) {
      lapply(seq_len(rs$S_cnt), function(si) {
        lapply(seq_len(n_lags), function(li) {
          ok_t <- rs$ok_strata[[si]][[li]]
          U_out <- matrix(0.0, rs$m, m_t)
          if (!is.null(ok_t) && sum(ok_t) > 1L) {
            m_ok <- sum(ok_t)
            Q_t <- diag(m_ok) - matrix(1.0 / m_ok, m_ok, m_ok)
            U_out[, ok_t] <- (op_strata_list[[si]][[li]][, ok_t, drop = FALSE] %*% Q_t) / sqrt(1.0 - 1.0 / m_ok)
          } else if (!is.null(ok_t) && sum(ok_t) == 1L) {
            U_out[, ok_t] <- op_strata_list[[si]][[li]][, ok_t, drop = FALSE]
          }
          U_out
        })
      })
    }

    U_cent_s <- center_op_list(rs$U_strata)
    V_cent_s <- center_op_list(rs$V_obs_strata)

    se_arr_r <- array(NA_real_, dim = c(m_t, P, n_lags))
    s2_trial_arr <- array(NA_real_, dim = c(m_t, P, n_lags))
    tcor_arr_r <- if (isTRUE(return_panel)) array(0.0, dim = c(K_tcor, m_t, P, n_lags)) else NULL
    lvec_arr_r <- if (isTRUE(return_panel) && n_lags >= 2L) array(0.0, dim = c(D_lag, m_t, P, n_lags)) else NULL
    svec_arr_r <- if (isTRUE(return_panel)) array(0.0, dim = c(D_emb_tot, m_t, P, n_lags)) else NULL

    # Vectorized over S_cnt strata (avoids O(P) inner operator crossprods!)
    for (si in seq_len(rs$S_cnt)) {
      cols_si <- which(rs_sem$s_idx == si)
      if (!length(cols_si)) next
      for (li in seq_len(n_lags)) {
        ok_t <- rs$ok_strata[[si]][[li]]
        if (is.null(ok_t) || !any(ok_t)) next
        Uc <- U_cent_s[[si]][[li]][, ok_t, drop = FALSE]
        Vc <- V_cent_s[[si]][[li]][, ok_t, drop = FALSE]
        m_ok <- sum(ok_t)

        # Precompute stratum temporal basis profiles (m_ok x (L1 + M_dyn))
        T_stat <- matrix(0.0, m_ok, L1)
        T_stat[, 1L] <- colSums(Uc * Uc)
        for (d in seq_len(L_max)) {
          pd <- rs$lag_pairs[[d]]
          if (nrow(pd) > 0L) T_stat[, d + 1L] <- colSums(Uc[pd[, 1L], , drop = FALSE] * Uc[pd[, 2L], , drop = FALSE])
        }
        T_dyn <- matrix(0.0, m_ok, M_dyn)
        for (j in seq_len(M_dyn)) {
          fj <- f_dyn_r[, j]
          gm_j <- if (length(j_rwls) == 1L && j == j_rwls) 1.0 else gamma_obs
          T_dyn[, j] <- (1.0 - gm_j) * colSums(Uc * (fj * Uc)) + gm_j * colSums(Vc * (fj * Vc))
        }

        rsc2_si <- if (!is.null(rs_sem$run_scale)) (rs_sem$run_scale[cols_si]^2) else rep(1.0, length(cols_si))
        v_base_m <- outer(T_stat[, 1L], (sig0_rob[cols_si]^2) * rsc2_si)
        v_dyn_m <- sweep(
          T_stat %*% Gamma_stat_diag[, cols_si, drop = FALSE] + T_dyn %*% B_dyn_diag[, cols_si, drop = FALSE],
          2L, rsc2_si, "*"
        )
        v_fin_m <- matrix(pmax(1e-10, 0.15 * v_base_m, v_dyn_m), nrow = m_ok, ncol = length(cols_si))
        s2_trial_arr[ok_t, cols_si, li] <- v_fin_m
        se_arr_r[ok_t, cols_si, li] <- sqrt(v_fin_m)

        if (isTRUE(return_panel) && m_t >= 2L) {
          Uu <- rs$U_strata[[si]][[li]]
          Vu <- rs$V_obs_strata[[si]][[li]]
          Tu_s <- matrix(0.0, m_t, L1)
          Tu_s[, 1L] <- colSums(Uu * Uu)
          for (d in seq_len(L_max)) {
            pd <- rs$lag_pairs[[d]]
            if (nrow(pd) > 0L) Tu_s[, d + 1L] <- colSums(Uu[pd[, 1L], , drop = FALSE] * Uu[pd[, 2L], , drop = FALSE])
          }
          Tu_d <- matrix(0.0, m_t, M_dyn)
          for (j in seq_len(M_dyn)) {
            fj <- f_dyn_r[, j]
            gm_j <- if (length(j_rwls) == 1L && j == j_rwls) 1.0 else gamma_obs
            Tu_d[, j] <- (1.0 - gm_j) * colSums(Uu * (fj * Uu)) + gm_j * colSums(Vu * (fj * Vu))
          }
          vu_fin <- matrix(
            pmax(1e-12, 0.15 * outer(Tu_s[, 1L], sig0_rob[cols_si]^2),
                 Tu_s %*% Gamma_stat_diag[, cols_si, drop = FALSE] + Tu_d %*% B_dyn_diag[, cols_si, drop = FALSE]),
            nrow = m_t, ncol = length(cols_si)
          )
          for (kd in seq_len(min(K_tcor, m_t - 1L))) {
            i_a <- seq_len(m_t - kd); i_b <- i_a + kd
            ok_ab <- ok_t[i_a] & ok_t[i_b]
            if (any(ok_ab)) {
              ia_k <- i_a[ok_ab]; ib_k <- i_b[ok_ab]
              Cu_s <- matrix(0.0, length(ia_k), L1)
              Cu_s[, 1L] <- colSums(Uu[, ia_k, drop = FALSE] * Uu[, ib_k, drop = FALSE])
              for (d in seq_len(L_max)) {
                pd <- rs$lag_pairs[[d]]
                if (nrow(pd) > 0L) Cu_s[, d + 1L] <- colSums(Uu[pd[, 1L], ia_k, drop = FALSE] * Uu[pd[, 2L], ib_k, drop = FALSE])
              }
              Cu_d <- matrix(0.0, length(ia_k), M_dyn)
              for (j in seq_len(M_dyn)) {
                fj <- f_dyn_r[, j]
                gm_j <- if (length(j_rwls) == 1L && j == j_rwls) 1.0 else gamma_obs
                Cu_d[, j] <- (1.0 - gm_j) * colSums(Uu[, ia_k, drop = FALSE] * (fj * Uu[, ib_k, drop = FALSE])) +
                  gm_j * colSums(Vu[, ia_k, drop = FALSE] * (fj * Vu[, ib_k, drop = FALSE]))
              }
              c_mat <- Cu_s %*% Gamma_stat_diag[, cols_si, drop = FALSE] + Cu_d %*% B_dyn_diag[, cols_si, drop = FALSE]
              den_m <- matrix(sqrt(pmax(1e-20, vu_fin[ia_k, , drop = FALSE] * vu_fin[ib_k, , drop = FALSE])), nrow = length(ia_k), ncol = length(cols_si))
              tcor_arr_r[kd, ia_k, cols_si, li] <- matrix(pmin(0.95, pmax(-0.95, c_mat / den_m)), nrow = length(ia_k), ncol = length(cols_si))
            }
          }
        }
      }

      if (isTRUE(return_panel) && n_lags >= 2L) {
        # Stratum unit-variance lag correlation factors scaled per voxel by sqrt(s2_trial_arr)
        R_lag_si <- array(0.0, dim = c(n_lags, n_lags, m_t))
        for (la in seq_len(n_lags)) {
          R_lag_si[la, la, ] <- 1.0
          if (la < n_lags) {
            Ua <- U_cent_s[[si]][[la]]
            na_v <- sqrt(pmax(colSums(Ua * Ua), 1e-16))
            for (lb in seq.int(la + 1L, n_lags)) {
              Ub <- U_cent_s[[si]][[lb]]
              nb_v <- sqrt(pmax(colSums(Ub * Ub), 1e-16))
              rab <- pmin(0.98, pmax(-0.98, colSums(Ua * Ub) / (na_v * nb_v)))
              R_lag_si[la, lb, ] <- rab; R_lag_si[lb, la, ] <- rab
            }
          }
        }
        for (ti in seq_len(m_t)) {
          Ri <- .project_trial_err_psd(R_lag_si[, , ti])
          eg_l <- eigen(Ri, symmetric = TRUE)
          k_l <- min(D_lag, sum(eg_l$values > 1e-12))
          if (k_l >= 1L) {
            Fl_u <- t(sweep(eg_l$vectors[, seq_len(k_l), drop = FALSE], 2L, sqrt(eg_l$values[seq_len(k_l)]), "*"))
            for (la in seq_len(n_lags)) {
              sd_v_ti <- sqrt(pmax(0.0, ifelse(is.finite(s2_trial_arr[ti, cols_si, la]), s2_trial_arr[ti, cols_si, la], 0.0)))
              lvec_arr_r[seq_len(k_l), ti, cols_si, la] <- tcrossprod(Fl_u[, la], sd_v_ti)
            }
          }
        }
      }
    }

    # Factor-space O(P K) shared_noise_vec per trial (never touches P x P!)
    if (isTRUE(return_panel) && K_poet >= 1L) {
      # Map each factor k = 1..K_poet to its dominant stratum (or stratum-weighted readout)
      W_str_fac <- matrix(0.0, rs$S_cnt, K_poet)
      for (si in seq_len(rs$S_cnt)) {
        cols_si <- which(rs_sem$s_idx == si)
        if (length(cols_si)) W_str_fac[si, ] <- colSums(L_poet[cols_si, , drop = FALSE]^2)
      }
      W_str_fac <- sweep(W_str_fac, 2L, pmax(colSums(W_str_fac), 1e-12), "/")
      ones_lag <- matrix(1.0, n_lags, n_lags)
      G_fac_L <- lapply(seq_len(L1), function(idx_l) kronecker(ones_lag, matrix(Gamma_stat_fac[idx_l, , ], K_poet, K_poet)))
      B_fac_L <- lapply(seq_len(M_dyn), function(j) kronecker(ones_lag, matrix(B_dyn_fac[j, , ], K_poet, K_poet)))
      q_f <- K_poet * n_lags

      for (ti in seq_len(m_t)) {
        U_f <- matrix(0.0, rs$m, q_f)
        V_f <- matrix(0.0, rs$m, q_f)
        for (la in seq_len(n_lags)) {
          U_s_la <- vapply(seq_len(rs$S_cnt), function(si) U_cent_s[[si]][[la]][, ti], numeric(rs$m))
          V_s_la <- vapply(seq_len(rs$S_cnt), function(si) V_cent_s[[si]][[la]][, ti], numeric(rs$m))
          if (!is.matrix(U_s_la)) U_s_la <- matrix(U_s_la, nrow = rs$m, ncol = rs$S_cnt)
          if (!is.matrix(V_s_la)) V_s_la <- matrix(V_s_la, nrow = rs$m, ncol = rs$S_cnt)
          cols_la <- (la - 1L) * K_poet + seq_len(K_poet)
          U_f[, cols_la] <- U_s_la %*% W_str_fac
          V_f[, cols_la] <- V_s_la %*% W_str_fac
        }
        Cf <- G_fac_L[[1L]] * crossprod(U_f)
        for (d in seq_len(L_max)) {
          Gd_f <- G_fac_L[[d + 1L]]
          if (any(abs(Gd_f) > 1e-14)) {
            pd <- rs$lag_pairs[[d]]
            if (nrow(pd) > 0L) Cf <- Cf + Gd_f * crossprod(U_f[pd[, 1L], , drop = FALSE], U_f[pd[, 2L], , drop = FALSE])
          }
        }
        for (j in seq_len(M_dyn)) {
          Bj_f <- B_fac_L[[j]]
          if (any(abs(Bj_f) > 1e-14)) {
            fj <- f_dyn_r[, j]
            gm_j <- if (length(j_rwls) == 1L && j == j_rwls) 1.0 else gamma_obs
            Cf <- Cf + Bj_f * ((1.0 - gm_j) * crossprod(U_f, fj * U_f) + gm_j * crossprod(V_f, fj * V_f))
          }
        }
        eg_cf <- eigen(0.5 * (Cf + t(Cf)), symmetric = TRUE)
        k_cf <- min(D_emb_tot, sum(eg_cf$values > 1e-12))
        if (k_cf >= 1L) {
          # F_fac: k_cf x (K_poet * n_lags)
          F_fac <- t(sweep(eg_cf$vectors[, seq_len(k_cf), drop = FALSE], 2L, sqrt(eg_cf$values[seq_len(k_cf)]), "*"))
          L_sc_r <- if (!is.null(rs_sem$run_scale)) sweep(L_sc_global, 1L, rs_sem$run_scale, "*") else L_sc_global
          for (la in seq_len(n_lags)) {
            cols_la <- (la - 1L) * K_poet + seq_len(K_poet)
            emb_v <- F_fac[, cols_la, drop = FALSE] %*% t(L_sc_r)
            # Cap embedding norm at 0.99 * se_arr_r so shared variance never exceeds total trial error variance
            nrm_emb <- sqrt(colSums(emb_v^2))
            se_cap <- 0.99 * ifelse(is.finite(se_arr_r[ti, , la]), se_arr_r[ti, , la], 0.0)
            over_e <- nrm_emb > se_cap & nrm_emb > 1e-14
            if (any(over_e)) {
              emb_v[, over_e] <- sweep(emb_v[, over_e, drop = FALSE], 2L, se_cap[over_e] / nrm_emb[over_e], "*")
            }
            svec_arr_r[seq_len(k_cf), ti, , la] <- emb_v
          }
        }
      }
    }

    run_panel_blocks[[r]] <- list(
      orig_r = keep_runs[r],
      se = se_arr_r,
      tcor = tcor_arr_r,
      lvec = lvec_arr_r,
      svec = svec_arr_r
    )
  }

  list(
    unit_error_cor = NULL,
    shared_loading_vec = shared_loading_vec,
    shared_error_model = shared_error_model,
    run_panel_blocks = run_panel_blocks
  )
}

# Collapse an exact or lowrank_voxel shared_error_model from voxels to parcels
.collapse_sem_to_parcels <- function(sem, parcel_ids, vox_idx_list, w_vox_list) {
  if (is.null(sem)) return(NULL)
  P_par <- length(parcel_ids)
  V <- length(sem$unit_ids)
  n_lags <- length(sem$lags)
  L_max <- sem$L_max
  L1 <- L_max + 1L
  M_dyn <- sem$M_dyn

  # Build V x P_par weight matrix W_vp
  W_vp <- matrix(0.0, V, P_par, dimnames = list(sem$unit_ids, parcel_ids))
  for (g in seq_len(P_par)) {
    vg <- vox_idx_list[[g]]
    if (length(vg)) W_vp[vg, g] <- w_vox_list[[g]]
  }

  if (identical(sem$mode, "exact") && !is.null(sem$Gamma_stat_gated)) {
    Sig0_rob_v <- sem$R_poet * outer(sem$sig0_rob, sem$sig0_rob)
    Sig0_sd_v <- sem$R_poet * outer(sem$sig0_sd, sem$sig0_sd)
    C_rob_p <- crossprod(W_vp, Sig0_rob_v %*% W_vp)
    C_sd_p <- crossprod(W_vp, Sig0_sd_v %*% W_vp)
    sig0_rob_p <- sqrt(pmax(diag(C_rob_p), 1e-10))
    sig0_sd_p <- sqrt(pmax(diag(C_sd_p), 1e-10))
    R_poet_p <- .poet_project_cor_psd(C_sd_p / outer(sig0_sd_p, sig0_sd_p))

    Gamma_p <- array(0.0, dim = c(L1, P_par, P_par))
    for (idx_l in seq_len(L1)) {
      Gamma_p[idx_l, , ] <- crossprod(W_vp, matrix(sem$Gamma_stat_gated[idx_l, , ], V, V) %*% W_vp)
    }
    B_dyn_p <- array(0.0, dim = c(M_dyn, P_par, P_par))
    for (j in seq_len(M_dyn)) {
      B_dyn_p[j, , ] <- crossprod(W_vp, matrix(sem$B_dyn_gated[j, , ], V, V) %*% W_vp)
    }
  } else {
    # Collapse lowrank_voxel sem (using run vox_sketch for exact within-parcel local spatial correlation when present!)
    has_sk <- all(vapply(sem$runs, function(rs) !is.null(rs$vox_sketch), logical(1L)))
    if (has_sk) {
      B_rob_all <- do.call(rbind, lapply(sem$runs, function(rs) rs$vox_sketch$B_rob_sk))
      df_vec <- if (!is.null(sem$df_v_tot) && length(sem$df_v_tot) == V) {
        pmax(1.0, sem$df_v_tot)
      } else {
        rep(max(1.0, sum(vapply(sem$runs, function(rs) max(1.0, rs$m - 30.0), numeric(1L)))), V)
      }
      B_rob_vn <- sweep(B_rob_all, 2L, 1.0 / sqrt(df_vec), "*")
      C_rob_sk <- crossprod(B_rob_vn %*% W_vp)
      idio_v <- pmax(0.0, (sem$sig0_rob^2) - colSums(B_rob_vn^2))
      diag(C_rob_sk) <- diag(C_rob_sk) + as.numeric(crossprod(W_vp^2, idio_v))
      sig0_rob_p <- sqrt(pmax(diag(C_rob_sk), 1e-10))
      R0_p <- C_rob_sk / outer(sig0_rob_p, sig0_rob_p)
      R_poet_p <- poet_regularize_cor(R0_p, df_u = mean(df_vec))$R
      ratio_sd <- sem$sig0_sd / pmax(sem$sig0_rob, 1e-10)
      sig0_sd_p <- sig0_rob_p * sqrt(pmax(1.0, as.numeric(crossprod(sweep(W_vp, 2L, pmax(colSums(W_vp), 1e-12), "/"), ratio_sd^2))))
    } else {
      K_p <- if (!is.null(sem$L_poet)) ncol(sem$L_poet) else 0L
      L_sc_sd <- if (K_p > 0L) sweep(sem$L_poet, 1L, sem$sig0_sd, "*") else matrix(0.0, V, 0L)
      L_sc_rob <- if (K_p > 0L) sweep(sem$L_poet, 1L, sem$sig0_rob, "*") else matrix(0.0, V, 0L)
      idio_sd <- pmax(0.0, (sem$sig0_sd^2) - rowSums(L_sc_sd^2))
      idio_rob <- pmax(0.0, (sem$sig0_rob^2) - rowSums(L_sc_rob^2))
      C_sd_p <- tcrossprod(crossprod(W_vp, L_sc_sd)) + diag(as.numeric(crossprod(W_vp^2, idio_sd)), P_par)
      C_rob_p <- tcrossprod(crossprod(W_vp, L_sc_rob)) + diag(as.numeric(crossprod(W_vp^2, idio_rob)), P_par)
      sig0_sd_p <- sqrt(pmax(diag(C_sd_p), 1e-10))
      sig0_rob_p <- sqrt(pmax(diag(C_rob_p), 1e-10))
      R_poet_p <- .poet_project_cor_psd(C_sd_p / outer(sig0_sd_p, sig0_sd_p))
    }

    K_p <- if (!is.null(sem$L_poet)) ncol(sem$L_poet) else 0L
    sig0_stat_v <- sqrt(pmax(sem$Gamma_stat_diag[1L, ], 1e-10))
    ratio_cond <- sig0_stat_v / pmax(sem$sig0_rob, 1e-10)
    sig0_cond_p <- sig0_rob_p * sqrt(pmax(0.15, as.numeric(crossprod(sweep(W_vp, 2L, pmax(colSums(W_vp), 1e-12), "/"), ratio_cond^2))))
    L_sc <- if (K_p > 0L) sweep(sem$L_poet, 1L, sig0_stat_v, "*") else matrix(0.0, V, 0L)
    WL <- crossprod(W_vp, L_sc)
    Gamma_p <- array(0.0, dim = c(L1, P_par, P_par))
    Gamma_p[1L, , ] <- R_poet_p * outer(sig0_cond_p, sig0_cond_p)
    for (idx_l in seq_len(L1)) {
      if (idx_l >= 2L && K_p > 0L) {
        Gamma_p[idx_l, , ] <- WL %*% matrix(sem$Gamma_stat_fac[idx_l, , ], K_p, K_p) %*% t(WL)
      }
      diag_gl <- as.numeric(crossprod(W_vp^2, sem$Gamma_stat_diag[idx_l, ] - if (K_p > 0L) rowSums((L_sc %*% matrix(sem$Gamma_stat_fac[idx_l, , ], K_p, K_p)) * L_sc) else 0.0))
      if (idx_l >= 2L) diag(Gamma_p[idx_l, , ]) <- diag(Gamma_p[idx_l, , ]) + diag_gl
    }
    B_dyn_p <- array(0.0, dim = c(M_dyn, P_par, P_par))
    for (j in seq_len(M_dyn)) {
      Bj_f_part <- if (K_p > 0L) WL %*% matrix(sem$B_dyn_fac[j, , ], K_p, K_p) %*% t(WL) else matrix(0.0, P_par, P_par)
      diag_bj <- as.numeric(crossprod(W_vp^2, sem$B_dyn_diag[j, ] - if (K_p > 0L) rowSums((L_sc %*% matrix(sem$B_dyn_fac[j, , ], K_p, K_p)) * L_sc) else 0.0))
      diag(Bj_f_part) <- diag(Bj_f_part) + diag_bj
      B_dyn_p[j, , ] <- 0.5 * (Bj_f_part + t(Bj_f_part))
    }
  }
  dimnames(R_poet_p) <- list(parcel_ids, parcel_ids)

  # Collapse run U_strata and V_obs_strata per parcel
  runs_p <- lapply(sem$runs, function(rs) {
    U_p <- vector("list", P_par)
    V_p <- vector("list", P_par)
    ok_p <- vector("list", P_par)
    for (g in seq_len(P_par)) {
      vg <- vox_idx_list[[g]]
      wg <- w_vox_list[[g]]
      s_vg <- rs$s_idx[vg]
      u_sg <- unique(s_vg)
      U_gl <- vector("list", n_lags)
      V_gl <- vector("list", n_lags)
      ok_gl <- vector("list", n_lags)
      for (li in seq_len(n_lags)) {
        if (length(u_sg) == 1L) {
          U_gl[[li]] <- rs$U_strata[[u_sg[1L]]][[li]]
          V_gl[[li]] <- rs$V_obs_strata[[u_sg[1L]]][[li]]
          ok_gl[[li]] <- rs$ok_strata[[u_sg[1L]]][[li]]
        } else {
          U_acc <- matrix(0.0, rs$m, rs$m_target)
          V_acc <- matrix(0.0, rs$m, rs$m_target)
          ok_acc <- rep(TRUE, rs$m_target)
          for (su in u_sg) {
            w_su <- sum(wg[s_vg == su])
            U_acc <- U_acc + w_su * rs$U_strata[[su]][[li]]
            V_acc <- V_acc + w_su * rs$V_obs_strata[[su]][[li]]
            ok_acc <- ok_acc & rs$ok_strata[[su]][[li]]
          }
          U_gl[[li]] <- U_acc
          V_gl[[li]] <- V_acc
          ok_gl[[li]] <- ok_acc
        }
      }
      U_p[[g]] <- U_gl; V_p[[g]] <- V_gl; ok_p[[g]] <- ok_gl
    }
    rsc_p <- if (!is.null(rs$run_scale) && length(rs$run_scale) == V) {
      W_nrm <- sweep(W_vp, 2L, pmax(colSums(W_vp), 1e-12), "/")
      sqrt(pmax(1e-8, as.numeric(crossprod(W_nrm, rs$run_scale^2))))
    } else rep(1.0, P_par)
    list(
      m = rs$m, rows = rs$rows, m_target = rs$m_target, s_idx = seq_len(P_par),
      U_strata = U_p, V_obs_strata = V_p, ok_strata = ok_p,
      e_task = rs$e_task, E_cov_env = rs$E_cov_env, lag_pairs = rs$lag_pairs,
      run_scale = rsc_p
    )
  })

  list(
    version = 3L,
    mode = "exact",
    unit_ids = parcel_ids,
    lags = sem$lags,
    run_ids = sem$run_ids,
    K_reg = sem$K_reg,
    M_dyn = sem$M_dyn,
    L_max = L_max,
    j_rwls = sem$j_rwls,
    vc_col_idx = sem$vc_col_idx,
    gamma_obs = sem$gamma_obs,
    w_res_global = sem$w_res_global,
    w_ltz_global = sem$w_ltz_global,
    R_poet = R_poet_p,
    poet_K = sem$poet_K,
    sig0_sd = sig0_sd_p,
    sig0_rob = sig0_rob_p,
    Gamma_stat_gated = Gamma_p,
    B_dyn_gated = B_dyn_p,
    runs = runs_p
  )
}

#' Aggregate a Voxel-Level MEDuSA Trial Panel to Parcels / ROIs
#'
#' Collapses a voxel-level \code{\link{decon_epoch_regression}} output (or stacked
#' trial-by-voxel \code{data.table}) into a parcel-level panel while propagating
#' shared cross-voxel and cross-parcel measurement-error vectors
#' (\code{shared_noise_vec}), cross-trial correlations (\code{trial_cor_vec}),
#' cross-lag embeddings (\code{lag_noise_vec}), exact parcel \code{se_demixed},
#' Empirical-Bayes \code{tau2_bio} / \code{precision_weight}, and any attached
#' \code{shared_error_model} for downstream \code{\link{decon_connectivity_repair}}
#' or \code{\link{mixed_by}} analyses.
#'
#' @param x A voxel-level \code{\link{decon_epoch_regression}} result list, a
#'   trial-level \code{data.frame}/\code{data.table} (including CSV/Parquet panels
#'   read via \code{data.table::fread} or \code{medusa_read_table}), a character
#'   file path (or vector of paths), or a list of per-subject results.
#' @param parcel_groups Either a named character/factor vector mapping voxel IDs
#'   to parcel labels, an unnamed vector matching \code{voxel_ids} order, or a
#'   \code{data.frame} containing voxel and parcel mapping columns.
#' @param parcel_weights Weighting rule within each parcel: \code{"auto"}
#'   (inverse squared \code{se_demixed} with 10\% median floor), \code{"equal"},
#'   or \code{"supplied"}.
#' @param weights Optional named numeric vector of voxel weights when
#'   \code{parcel_weights = "supplied"}.
#' @param voxel_col Name of the voxel column in the input panel (default \code{"voxel"}).
#' @param parcel_col Name of the output parcel column (default \code{"parcel"}).
#' @param id_col,run_col,trial_col,lag_col Subject, run, trial, and lag column names.
#'
#' @return If \code{x} is a \code{decon_epoch_regression} result list, an updated
#'   parcel-level result list with \code{panel}, \code{parcel_weights},
#'   \code{unit_error_cor}, and \code{shared_error_model}. If \code{x} is a
#'   \code{data.frame}/\code{data.table}, a parcel-aggregated \code{data.table}.
#' @export
aggregate_decon_panel <- function(
  x,
  parcel_groups,
  parcel_weights = c("auto", "equal", "supplied"),
  weights = NULL,
  voxel_col = "voxel",
  parcel_col = "parcel",
  id_col = "id",
  run_col = "run",
  trial_col = "trial",
  lag_col = "evt_time"
) {
  parcel_weights <- match.arg(parcel_weights)

  # Handle file paths (single path or vector of per-subject paths)
  if (is.character(x)) {
    if (length(x) == 1L) {
      x <- medusa_read_table(x)
    } else if (length(x) > 1L) {
      return(lapply(x, aggregate_decon_panel,
        parcel_groups = parcel_groups, parcel_weights = parcel_weights, weights = weights,
        voxel_col = voxel_col, parcel_col = parcel_col, id_col = id_col,
        run_col = run_col, trial_col = trial_col, lag_col = lag_col
      ))
    }
  }

  # Handle list of per-subject results/panels recursively
  if (is.list(x) && !is.data.frame(x) && is.null(x$panel)) {
    return(lapply(x, aggregate_decon_panel,
      parcel_groups = parcel_groups, parcel_weights = parcel_weights, weights = weights,
      voxel_col = voxel_col, parcel_col = parcel_col, id_col = id_col,
      run_col = run_col, trial_col = trial_col, lag_col = lag_col
    ))
  }

  is_fit_obj <- is.list(x) && !is.data.frame(x) && !is.null(x$panel)
  p_in <- if (is_fit_obj) x$panel else x
  if (!is.data.frame(p_in)) {
    stop("`x` must be a `decon_epoch_regression(..., return_panel = TRUE)` result, a panel data.frame, or a file path.")
  }
  sem_in <- if (is_fit_obj && !is.null(x$shared_error_model)) {
    x$shared_error_model
  } else {
    attr(p_in, "shared_error_model")
  }
  dt <- restore_panel_vec_cols(data.table::copy(data.table::as.data.table(p_in)))
  if (!(voxel_col %in% names(dt))) {
    stop("Column '", voxel_col, "' not found in panel.")
  }

  vox_u <- if (!is.null(sem_in)) as.character(sem_in$unit_ids) else unique(as.character(dt[[voxel_col]]))

  # Resolve voxel -> parcel map
  supplied_w_map <- NULL
  if (is.data.frame(parcel_groups)) {
    v_c <- intersect(c(voxel_col, "voxel", "unit"), names(parcel_groups))[1L]
    p_c <- intersect(c(parcel_col, "parcel", "region", "group"), names(parcel_groups))[1L]
    if (is.na(v_c) || is.na(p_c)) {
      stop("`parcel_groups` data.frame must contain voxel and parcel/region columns.")
    }
    pg_map <- stats::setNames(as.character(parcel_groups[[p_c]]), as.character(parcel_groups[[v_c]]))
    if ("weight" %in% names(parcel_groups)) {
      supplied_w_map <- stats::setNames(as.numeric(parcel_groups[["weight"]]), as.character(parcel_groups[[v_c]]))
    }
  } else {
    pg_chr <- as.character(parcel_groups)
    if (!is.null(names(parcel_groups)) && any(nzchar(names(parcel_groups)))) {
      pg_map <- stats::setNames(pg_chr, names(parcel_groups))
    } else if (length(pg_chr) == length(vox_u)) {
      pg_map <- stats::setNames(pg_chr, vox_u)
    } else {
      stop("Unnamed `parcel_groups` must have length equal to number of voxels (", length(vox_u), ").")
    }
  }
  if (!is.null(weights)) {
    if (!is.null(names(weights))) {
      supplied_w_map <- stats::setNames(as.numeric(weights), names(weights))
    } else if (length(weights) == length(vox_u)) {
      supplied_w_map <- stats::setNames(as.numeric(weights), vox_u)
    }
  }

  dt[, .parcel_tmp := unname(pg_map[as.character(get(voxel_col))])]
  dt <- dt[!is.na(.parcel_tmp) & nzchar(.parcel_tmp)]
  if (!nrow(dt)) stop("No panel rows matched `parcel_groups`.")
  parcel_ids <- unique(dt$.parcel_tmp)
  P_par <- length(parcel_ids)

  has_id <- id_col %in% names(dt)
  if (!has_id) dt[, (id_col) := ".single_subj"]

  # Compute normalized within-parcel voxel weights per subject
  if (identical(parcel_weights, "equal") || !("se_demixed" %in% names(dt))) {
    w_dt <- unique(dt[, c(id_col, voxel_col, ".parcel_tmp"), with = FALSE])
    w_dt[, .w_raw := 1.0]
  } else if (identical(parcel_weights, "supplied")) {
    if (is.null(supplied_w_map)) stop("`parcel_weights = 'supplied'` requires `weights` or a `weight` column in `parcel_groups`.")
    w_dt <- unique(dt[, c(id_col, voxel_col, ".parcel_tmp"), with = FALSE])
    w_dt[, .w_raw := pmax(0.0, unname(supplied_w_map[as.character(get(voxel_col))]))]
    w_dt[!is.finite(.w_raw) | .w_raw <= 0, .w_raw := 1.0]
  } else {
    w_dt <- dt[, list(
      .s2_bar = mean(se_demixed^2, na.rm = TRUE)
    ), by = c(id_col, voxel_col, ".parcel_tmp")]
    w_dt[!is.finite(.s2_bar) | .s2_bar <= 0, .s2_bar := 1.0]
    w_dt[, .w_raw := {
      med_s2 <- stats::median(.s2_bar, na.rm = TRUE)
      1.0 / pmax(.s2_bar, 0.10 * max(med_s2, 1e-8))
    }, by = c(id_col, ".parcel_tmp")]
  }
  w_dt[, .w_norm := .w_raw / sum(.w_raw), by = c(id_col, ".parcel_tmp")]
  dt <- merge(dt, w_dt[, c(id_col, voxel_col, ".w_norm"), with = FALSE], by = c(id_col, voxel_col), sort = FALSE)

  # Aggregate numeric & list-vector columns per (id, run, trial, lag, parcel)
  grp_keys <- c(id_col, run_col, trial_col, lag_col, ".parcel_tmp")
  sig_cols <- intersect(
    c("decon_raw", "decon_demixed", "decon_eb", "decon_demixed_whitened", "decon_eb_whitened"),
    names(dt)
  )
  has_svec <- "shared_noise_vec" %in% names(dt)
  has_lvec <- "lag_noise_vec" %in% names(dt)
  has_tcor <- "trial_cor_vec" %in% names(dt)
  has_se <- "se_demixed" %in% names(dt)
  has_se_raw <- "se_raw" %in% names(dt)
  has_ac1 <- "error_ac1" %in% names(dt)

  carry_cols <- setdiff(
    names(dt),
    c(grp_keys, voxel_col, parcel_col, ".w_norm", sig_cols,
      "se_raw", "se_demixed", "se_whitened", "reliability_trial",
      "precision_weight", "precision_weight_cohort", "precision_weight_op",
      "precision_weight_raw", "error_ac1", "tau2_bio", "shared_noise_loading",
      "shared_noise_vec", "lag_noise_vec", "trial_cor_vec")
  )

  agg_vec_weighted <- function(vlist, w) {
    lens <- vapply(vlist, length, integer(1L))
    d <- max(c(0L, lens))
    if (d == 0L) return(numeric(0L))
    acc <- numeric(d)
    for (i in seq_along(vlist)) {
      vi <- as.numeric(vlist[[i]])
      if (length(vi) == d && all(is.finite(vi))) acc <- acc + w[i] * vi
    }
    acc
  }

  out_dt <- dt[, {
    wv <- .w_norm / sum(.w_norm)
    res <- list()
    for (sc in sig_cols) {
      vals <- as.numeric(get(sc))
      ok <- is.finite(vals)
      res[[sc]] <- if (any(ok)) sum(wv[ok] * vals[ok]) / sum(wv[ok]) else NA_real_
    }
    u_bar <- if (has_svec) agg_vec_weighted(shared_noise_vec, wv) else numeric(0L)
    if (has_se) {
      s2_v <- as.numeric(se_demixed)^2
      u2_v <- if (has_svec && length(u_bar) > 0L) {
        vapply(shared_noise_vec, function(uv) sum(as.numeric(uv)^2), numeric(1L))
      } else numeric(length(s2_v))
      idio_v <- pmax(0.0, ifelse(is.finite(s2_v), s2_v, 0.0) - u2_v)
      s2_par <- sum(u_bar^2) + sum((wv^2) * idio_v)
      res[["se_demixed"]] <- sqrt(max(s2_par, 1e-12))
    }
    if (has_se_raw) {
      sr_v <- as.numeric(se_raw)^2
      res[["se_raw"]] <- sqrt(max(sum((wv^2) * ifelse(is.finite(sr_v), sr_v, 0.0)), 1e-12))
    }
    if (has_ac1) {
      ac_v <- as.numeric(error_ac1)
      ok_a <- is.finite(ac_v)
      res[["error_ac1"]] <- if (any(ok_a)) sum(wv[ok_a] * ac_v[ok_a]) / sum(wv[ok_a]) else 0.0
    }
    if (has_svec) res[["shared_noise_vec"]] <- list(u_bar)
    if (has_tcor) res[["trial_cor_vec"]] <- list(agg_vec_weighted(trial_cor_vec, wv))
    if (has_lvec) {
      l_bar <- agg_vec_weighted(lag_noise_vec, wv)
      if (length(l_bar) > 0L && has_se && sum(l_bar^2) > 1e-16) {
        l_bar <- l_bar * (res[["se_demixed"]] / sqrt(sum(l_bar^2)))
      }
      res[["lag_noise_vec"]] <- list(l_bar)
    }
    for (cc in carry_cols) {
      res[[cc]] <- get(cc)[1L]
    }
    res
  }, by = grp_keys]

  data.table::setnames(out_dt, ".parcel_tmp", parcel_col)
  if (parcel_col != "voxel" && !("voxel" %in% names(out_dt))) {
    out_dt[, voxel := get(parcel_col)]
  }
  data.table::setorderv(out_dt, c(id_col, lag_col, parcel_col, run_col, trial_col))

  # If shared_error_model is attached (single-subject fit or panel), collapse it and refresh exact parcel SEs/vectors!
  sem_par <- NULL
  R_poet_par <- NULL
  p_weights_out <- NULL
  if (!is.null(sem_in) && length(unique(out_dt[[id_col]])) == 1L) {
    vox_idx_l <- stats::setNames(vector("list", P_par), parcel_ids)
    w_vox_l <- stats::setNames(vector("list", P_par), parcel_ids)
    for (g in seq_len(P_par)) {
      sub_w <- w_dt[.parcel_tmp == parcel_ids[g]]
      m_v <- match(as.character(sub_w[[voxel_col]]), sem_in$unit_ids)
      ok_m <- !is.na(m_v)
      vox_idx_l[[g]] <- m_v[ok_m]
      wv_g <- sub_w$.w_norm[ok_m]
      w_vox_l[[g]] <- wv_g / sum(wv_g)
    }
    p_weights_out <- w_vox_l
    sem_par <- .collapse_sem_to_parcels(sem_in, parcel_ids, vox_idx_l, w_vox_l)
    R_poet_par <- sem_par$R_poet
    sh_load_par <- compute_shared_loading_from_cor(R_poet_par, parcel_ids)
    out_dt[, shared_noise_loading := unname(sh_load_par[as.character(get(parcel_col))])]

    # Refresh parcel se_demixed, trial_cor_vec, and shared_noise_vec from exact collapsed sem_par!
    n_lags_p <- length(sem_par$lags)
    K_tcor <- 4L
    D_emb_p <- min(24L, max(2L, P_par * n_lags_p))
    key_p <- paste(as.character(out_dt[[run_col]]), as.character(out_dt[[parcel_col]]),
                   sprintf("%.6f", as.numeric(out_dt[[lag_col]])), sep = "\r")
    row_lk <- split(seq_len(nrow(out_dt)), key_p)
    se_p_col <- as.numeric(out_dt$se_demixed)
    tcor_p_col <- if ("trial_cor_vec" %in% names(out_dt)) out_dt$trial_cor_vec else replicate(nrow(out_dt), numeric(K_tcor), simplify = FALSE)
    svec_p_col <- replicate(nrow(out_dt), numeric(D_emb_p), simplify = FALSE)

    for (ri in seq_along(sem_par$runs)) {
      rs_p <- sem_par$runs[[ri]]
      r_lbl <- as.character(sem_par$run_ids[ri])
      m_t <- rs_p$m_target
      U_cs <- lapply(seq_len(P_par), function(g) {
        lapply(seq_len(n_lags_p), function(li) {
          ok_t <- rs_p$ok_strata[[g]][[li]]
          Uo <- matrix(0.0, rs_p$m, m_t)
          if (!is.null(ok_t) && sum(ok_t) > 1L) {
            mok <- sum(ok_t); Qt <- diag(mok) - matrix(1.0 / mok, mok, mok)
            Uo[, ok_t] <- (rs_p$U_strata[[g]][[li]][, ok_t, drop = FALSE] %*% Qt) / sqrt(1.0 - 1.0 / mok)
          } else if (!is.null(ok_t) && sum(ok_t) == 1L) Uo[, ok_t] <- rs_p$U_strata[[g]][[li]][, ok_t, drop = FALSE]
          Uo
        })
      })
      V_cs <- lapply(seq_len(P_par), function(g) {
        lapply(seq_len(n_lags_p), function(li) {
          ok_t <- rs_p$ok_strata[[g]][[li]]
          Vo <- matrix(0.0, rs_p$m, m_t)
          if (!is.null(ok_t) && sum(ok_t) > 1L) {
            mok <- sum(ok_t); Qt <- diag(mok) - matrix(1.0 / mok, mok, mok)
            Vo[, ok_t] <- (rs_p$V_obs_strata[[g]][[li]][, ok_t, drop = FALSE] %*% Qt) / sqrt(1.0 - 1.0 / mok)
          } else if (!is.null(ok_t) && sum(ok_t) == 1L) Vo[, ok_t] <- rs_p$V_obs_strata[[g]][[li]][, ok_t, drop = FALSE]
          Vo
        })
      })
      s2_p_arr <- array(NA_real_, dim = c(m_t, P_par, n_lags_p))
      for (g in seq_len(P_par)) {
        for (li in seq_len(n_lags_p)) {
          ok_t <- rs_p$ok_strata[[g]][[li]]
          if (is.null(ok_t) || !any(ok_t)) next
          Uc <- U_cs[[g]][[li]][, ok_t, drop = FALSE]; Vc <- V_cs[[g]][[li]][, ok_t, drop = FALSE]
          vb <- pmax(1e-10, medusa_eval_sem_bilinear(Uc, Uc, Vc, Vc, rs_p, sem_par, g, g, stat_floor = TRUE))
          vd <- pmax(0.15 * vb, medusa_eval_sem_bilinear(Uc, Uc, Vc, Vc, rs_p, sem_par, g, g, stat_floor = FALSE))
          s2f <- rep(NA_real_, m_t); s2f[ok_t] <- vd
          s2_p_arr[, g, li] <- s2f
          tc_m <- matrix(0.0, K_tcor, m_t)
          if (m_t >= 2L) {
            Uu <- rs_p$U_strata[[g]][[li]]; Vu <- rs_p$V_obs_strata[[g]][[li]]
            vub <- pmax(1e-12, medusa_eval_sem_bilinear(Uu, Uu, Vu, Vu, rs_p, sem_par, g, g, stat_floor = TRUE))
            vud <- pmax(0.15 * vub, medusa_eval_sem_bilinear(Uu, Uu, Vu, Vu, rs_p, sem_par, g, g, stat_floor = FALSE))
            for (kd in seq_len(min(K_tcor, m_t - 1L))) {
              ia <- seq_len(m_t - kd); ib <- ia + kd; okab <- ok_t[ia] & ok_t[ib]
              if (any(okab)) {
                ck <- medusa_eval_sem_bilinear(Uu[, ia[okab], drop = FALSE], Uu[, ib[okab], drop = FALSE],
                                               Vu[, ia[okab], drop = FALSE], Vu[, ib[okab], drop = FALSE], rs_p, sem_par, g, g)
                tc_m[kd, ia[okab]] <- pmin(0.95, pmax(-0.95, ck / sqrt(pmax(1e-20, vud[ia[okab]] * vud[ib[okab]]))))
              }
            }
          }
          lk <- paste(r_lbl, parcel_ids[g], sprintf("%.6f", sem_par$lags[li]), sep = "\r")
          hr <- row_lk[[lk]]
          if (!is.null(hr) && length(hr) == m_t) {
            otr <- match(out_dt[[trial_col]][hr], seq_len(m_t))
            se_p_col[hr] <- sqrt(s2f[otr])
            tcor_p_col[hr] <- lapply(otr, function(ti) as.numeric(tc_m[, ti]))
          }
        }
      }
      if (P_par >= 2L) {
        q_tot <- P_par * n_lags_p
        ones_l <- matrix(1.0, n_lags_p, n_lags_p)
        G_s <- lapply(seq_len(sem_par$L_max + 1L), function(il) kronecker(ones_l, matrix(sem_par$Gamma_stat_gated[il, , ], P_par, P_par)))
        B_d <- lapply(seq_len(sem_par$M_dyn), function(j) kronecker(ones_l, matrix(sem_par$B_dyn_gated[j, , ], P_par, P_par)))
        f_dr <- cbind(rs_p$e_task, rs_p$E_cov_env)
        emb_p <- array(0.0, dim = c(D_emb_p, q_tot, m_t))
        for (ti in seq_len(m_t)) {
          Ut <- matrix(0.0, rs_p$m, q_tot); Vt <- matrix(0.0, rs_p$m, q_tot)
          for (la in seq_len(n_lags_p)) {
            for (g in seq_len(P_par)) {
              ck <- (la - 1L) * P_par + g
              Ut[, ck] <- U_cs[[g]][[la]][, ti]; Vt[, ck] <- V_cs[[g]][[la]][, ti]
            }
          }
          Ci <- G_s[[1L]] * crossprod(Ut)
          for (d in seq_len(sem_par$L_max)) {
            if (any(abs(G_s[[d + 1L]]) > 1e-14)) {
              pd <- rs_p$lag_pairs[[d]]
              if (nrow(pd) > 0L) Ci <- Ci + G_s[[d + 1L]] * crossprod(Ut[pd[, 1L], , drop = FALSE], Ut[pd[, 2L], , drop = FALSE])
            }
          }
          for (j in seq_len(sem_par$M_dyn)) {
            if (any(abs(B_d[[j]]) > 1e-14)) {
              fj <- f_dr[, j]
              gm_j <- if (!is.null(sem_par$j_rwls) && length(sem_par$j_rwls) == 1L && j == sem_par$j_rwls) 1.0 else sem_par$gamma_obs
              Ci <- Ci + B_d[[j]] * ((1.0 - gm_j) * crossprod(Ut, fj * Ut) + gm_j * crossprod(Vt, fj * Vt))
            }
          }
          if (!is.null(rs_p$run_scale)) {
            Ci <- Ci * tcrossprod(rep(rs_p$run_scale, times = n_lags_p))
          }
          for (la in seq_len(n_lags_p)) {
            id_d <- (la - 1L) * P_par + seq_len(P_par)
            diag(Ci)[id_d] <- ifelse(is.finite(s2_p_arr[ti, , la]), s2_p_arr[ti, , la], 0.0)
          }
          emb_p[, , ti] <- .extract_psd_factors(Ci, D_emb_p, nC = P_par, n_lags = n_lags_p)
        }
        for (g in seq_len(P_par)) {
          for (la in seq_len(n_lags_p)) {
            ra <- (la - 1L) * P_par + g
            lk <- paste(r_lbl, parcel_ids[g], sprintf("%.6f", sem_par$lags[la]), sep = "\r")
            hr <- row_lk[[lk]]
            if (!is.null(hr) && length(hr) == m_t) {
              otr <- match(out_dt[[trial_col]][hr], seq_len(m_t))
              svec_p_col[hr] <- lapply(otr, function(ti) as.numeric(emb_p[, ra, ti]))
            }
          }
        }
      }
    }
    data.table::set(out_dt, j = "se_demixed", value = se_p_col)
    data.table::set(out_dt, j = "trial_cor_vec", value = tcor_p_col)
    data.table::set(out_dt, j = "shared_noise_vec", value = svec_p_col)
  } else if ("shared_noise_vec" %in% names(out_dt)) {
    # Populate approximate scalar shared_noise_loading on bare stacked data.table
    out_dt[, shared_noise_loading := {
      u2 <- vapply(shared_noise_vec, function(uv) sum(as.numeric(uv)^2), numeric(1L))
      s2 <- pmax(se_demixed^2, 1e-12)
      sqrt(pmin(0.99, pmax(0.0, mean(u2 / s2, na.rm = TRUE))))
    }, by = c(id_col, parcel_col)]
  }

  # Refresh parcel-level tau2_bio, reliability_trial, decon_eb, and precision_weight
  if (all(c("decon_demixed", "se_demixed") %in% names(out_dt))) {
    out_dt[, c("tau2_bio", "reliability_trial", "decon_eb", "precision_weight") := {
      y_c <- decon_demixed - ave(decon_demixed, get(run_col), FUN = function(z) mean(z, na.rm = TRUE))
      s2 <- se_demixed^2
      ok <- is.finite(y_c) & is.finite(s2) & (s2 > 0)
      if (any(ok)) {
        t2 <- max(0.0, mean(y_c[ok]^2) - mean(s2[ok]))
        fl <- 0.10 * stats::median(s2[ok])
        rel <- t2 / pmax(t2 + s2, 1e-12)
        mu_r <- ave(decon_demixed, get(run_col), FUN = function(z) mean(z, na.rm = TRUE))
        eb_v <- mu_r + rel * y_c
        pw <- 1.0 / pmax(t2 + fl + s2, 1e-12)
        pw <- pw / mean(pw[ok])
        list(rep(t2, .N), rel, eb_v, pw)
      } else {
        list(rep(0.0, .N), rep(0.0, .N), decon_demixed, rep(1.0, .N))
      }
    }, by = c(id_col, parcel_col, lag_col)]
  }

  if (!has_id) out_dt[, (id_col) := NULL]

  for (at_nm in c("tr", "demix_target_trials")) {
    if (!is.null(attr(p_in, at_nm))) attr(out_dt, at_nm) <- attr(p_in, at_nm)
  }
  if (!is.null(R_poet_par)) attr(out_dt, "unit_error_cor") <- R_poet_par
  if (!is.null(sem_par)) attr(out_dt, "shared_error_model") <- sem_par

  if (is_fit_obj) {
    res_obj <- x
    res_obj$panel <- out_dt
    res_obj$parcel_weights <- p_weights_out
    res_obj$unit_error_cor <- R_poet_par
    res_obj$shared_error_model <- sem_par
    if (!is.null(res_obj$settings)) res_obj$settings$aggregate_by_parcel <- TRUE
    return(res_obj)
  }
  out_dt
}
