#include "medusa.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace {

inline arma::vec extract_obs_vec_impl(
    const arma::mat& C_mat,
    const arma::uvec& adj_idx0,
    const std::vector<arma::umat>& lag_pairs0,
    const arma::uword L_max
) {
  const arma::uword m = C_mat.n_rows;
  if (m == 0) return arma::vec();
  const arma::uword n_hi = (L_max >= 2) ? (L_max - 1) : 0;
  const arma::uword n_obs = 2 * m - 1 + n_hi;
  arma::vec out(n_obs, arma::fill::zeros);

  for (arma::uword i = 0; i < m; ++i) {
    out(i) = C_mat(i, i);
  }

  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
  for (arma::uword k = 0; k < adj_idx0.n_elem; ++k) {
    const arma::uword i0 = adj_idx0(k);
    if (i0 + 1 < m) {
      out(m + i0) = (C_mat(i0 + 1, i0) + C_mat(i0, i0 + 1)) * inv_sqrt2;
    }
  }

  if (L_max >= 2) {
    for (arma::uword d = 2; d <= L_max; ++d) {
      const arma::umat& pd = lag_pairs0[d - 1];
      const double N_d = 0.5 * static_cast<double>(pd.n_rows);
      if (N_d > 0.0) {
        double acc = 0.0;
        for (arma::uword r = 0; r < pd.n_rows; ++r) {
          const arma::uword i1 = pd(r, 0);
          const arma::uword i2 = pd(r, 1);
          if (i1 < m && i2 < m) {
            acc += C_mat(i1, i2);
          }
        }
        out(2 * m - 1 + (d - 2)) = acc / std::sqrt(2.0 * N_d);
      }
    }
  }

  return out;
}

inline arma::vec extract_y_vec_impl(
    const double* eg,
    const double* eh,
    const arma::uword m,
    const arma::uvec& adj_idx0,
    const std::vector<arma::umat>& lag_pairs0,
    const arma::uword L_max
) {
  if (m == 0) return arma::vec();
  const arma::uword n_hi = (L_max >= 2) ? (L_max - 1) : 0;
  const arma::uword n_obs = 2 * m - 1 + n_hi;
  arma::vec out(n_obs, arma::fill::zeros);

  for (arma::uword i = 0; i < m; ++i) {
    out(i) = eg[i] * eh[i];
  }

  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
  for (arma::uword k = 0; k < adj_idx0.n_elem; ++k) {
    const arma::uword i0 = adj_idx0(k);
    if (i0 + 1 < m) {
      out(m + i0) = (eg[i0 + 1] * eh[i0] + eg[i0] * eh[i0 + 1]) * inv_sqrt2;
    }
  }

  if (L_max >= 2) {
    for (arma::uword d = 2; d <= L_max; ++d) {
      const arma::umat& pd = lag_pairs0[d - 1];
      const double N_d = 0.5 * static_cast<double>(pd.n_rows);
      if (N_d > 0.0) {
        double acc = 0.0;
        for (arma::uword r = 0; r < pd.n_rows; ++r) {
          const arma::uword i1 = pd(r, 0);
          const arma::uword i2 = pd(r, 1);
          if (i1 < m && i2 < m) {
            acc += eg[i1] * eh[i2];
          }
        }
        out(2 * m - 1 + (d - 2)) = acc / std::sqrt(2.0 * N_d);
      }
    }
  }

  return out;
}

inline arma::uvec to_zero_based_uvec(const Rcpp::IntegerVector& iv) {
  arma::uvec out(iv.size());
  for (R_xlen_t i = 0; i < iv.size(); ++i) {
    out(static_cast<arma::uword>(i)) = static_cast<arma::uword>(std::max(1, iv[i]) - 1);
  }
  return out;
}

inline std::vector<arma::umat> to_zero_based_lag_pairs(const Rcpp::List& lag_pairs_r, arma::uword L_max) {
  std::vector<arma::umat> out(L_max);
  for (arma::uword d = 0; d < L_max; ++d) {
    if (static_cast<R_xlen_t>(d) < lag_pairs_r.size() && !Rf_isNull(lag_pairs_r[d])) {
      Rcpp::IntegerMatrix pm = Rcpp::as<Rcpp::IntegerMatrix>(lag_pairs_r[d]);
      arma::umat um(static_cast<arma::uword>(pm.nrow()), 2);
      for (int r = 0; r < pm.nrow(); ++r) {
        um(static_cast<arma::uword>(r), 0) = static_cast<arma::uword>(std::max(1, pm(r, 0)) - 1);
        um(static_cast<arma::uword>(r), 1) = static_cast<arma::uword>(std::max(1, pm(r, 1)) - 1);
      }
      out[d] = std::move(um);
    } else {
      out[d] = arma::umat(0, 2);
    }
  }
  return out;
}

inline arma::mat damped_sym_inv(const arma::mat& A, double ridge) {
  arma::mat As = 0.5 * (A + A.t());
  As.diag() += ridge;
  arma::mat invA;
  if (arma::inv_sympd(invA, As)) {
    return invA;
  }
  return arma::solve(As, arma::eye<arma::mat>(As.n_rows, As.n_cols), arma::solve_opts::likely_sympd);
}

} // namespace

// [[Rcpp::export]]
arma::mat poet_build_pair_X_cpp(
    const arma::mat& M_A,
    const arma::mat& M_B,
    const arma::mat& G_A,
    const arma::mat& G_B,
    const arma::mat& H_A,
    const arma::mat& H_B,
    const Rcpp::IntegerVector& adj_idx1,
    const Rcpp::List& lag_pairs,
    const arma::mat& f_dyn,
    int L_max_in
) {
  const arma::uword m = M_A.n_rows;
  const arma::uword L_max = static_cast<arma::uword>(std::max(1, L_max_in));
  const arma::uword L1 = L_max + 1;
  const arma::uword M_dyn = f_dyn.n_cols;
  const arma::uword B_tot = L1 + 2 * M_dyn;
  if (m == 0) return arma::mat(0, B_tot);
  const arma::uword n_hi = (L_max >= 2) ? (L_max - 1) : 0;
  const arma::uword n_obs = 2 * m - 1 + n_hi;

  const arma::uvec adj_idx0 = to_zero_based_uvec(adj_idx1);
  const std::vector<arma::umat> lag_pairs0 = to_zero_based_lag_pairs(lag_pairs, L_max);

  arma::mat X_out(n_obs, B_tot, arma::fill::zeros);

  const arma::mat C0 = M_A * M_B;
  X_out.col(0) = extract_obs_vec_impl(C0, adj_idx0, lag_pairs0, L_max);

  for (arma::uword d = 1; d <= L_max; ++d) {
    const arma::umat& pd = lag_pairs0[d - 1];
    arma::mat Jd_MB(m, M_B.n_cols, arma::fill::zeros);
    for (arma::uword r = 0; r < pd.n_rows; ++r) {
      const arma::uword i1 = pd(r, 0);
      const arma::uword i2 = pd(r, 1);
      if (i1 != i2 && i1 < m && i2 < M_B.n_rows) {
        Jd_MB.row(i1) += M_B.row(i2);
      }
    }
    const arma::mat Cd = M_A * Jd_MB;
    X_out.col(d) = extract_obs_vec_impl(Cd, adj_idx0, lag_pairs0, L_max);
  }

  if (M_dyn > 0) {
    const arma::mat G_Bt = G_B.t();
    const arma::mat H_Bt = H_B.t();
    for (arma::uword j = 0; j < M_dyn; ++j) {
      const arma::vec fk = f_dyn.col(j);
      if (arma::all(fk == 0.0)) continue;

      arma::mat GB_scaled = G_Bt;
      GB_scaled.each_col() %= fk;
      const arma::mat C_inn = G_A * GB_scaled;
      X_out.col(L1 + j) = extract_obs_vec_impl(C_inn, adj_idx0, lag_pairs0, L_max);

      arma::mat HB_scaled = H_Bt;
      HB_scaled.each_col() %= fk;
      const arma::mat C_obs = H_A * HB_scaled;
      X_out.col(L1 + M_dyn + j) = extract_obs_vec_impl(C_obs, adj_idx0, lag_pairs0, L_max);
    }
  }

  return X_out;
}

// [[Rcpp::export]]
arma::mat poet_extract_Xty_diag_cpp(
    const Rcpp::List& E_runs,
    const Rcpp::List& X_diag_runs,
    const Rcpp::IntegerVector& sig_idx,
    const Rcpp::List& adj_idx1_runs,
    const Rcpp::List& lag_pairs_runs,
    int L_max_in,
    int B_tot_in
) {
  const arma::uword n_runs = static_cast<arma::uword>(E_runs.size());
  const arma::uword P = static_cast<arma::uword>(sig_idx.size());
  const arma::uword L_max = static_cast<arma::uword>(std::max(1, L_max_in));
  const arma::uword B_tot = static_cast<arma::uword>(B_tot_in);

  arma::mat Xty_diag(B_tot, P, arma::fill::zeros);
  if (n_runs == 0 || P == 0) return Xty_diag;

  for (arma::uword r = 0; r < n_runs; ++r) {
    const arma::mat E_r = Rcpp::as<arma::mat>(E_runs[r]);
    const Rcpp::List X_list_r = Rcpp::as<Rcpp::List>(X_diag_runs[r]);
    const arma::uvec adj_idx0 = to_zero_based_uvec(Rcpp::as<Rcpp::IntegerVector>(adj_idx1_runs[r]));
    const std::vector<arma::umat> lag_pairs0 = to_zero_based_lag_pairs(Rcpp::as<Rcpp::List>(lag_pairs_runs[r]), L_max);
    const arma::uword m = E_r.n_rows;

    const arma::uword G_sig = static_cast<arma::uword>(X_list_r.size());
    std::vector<arma::mat> X_t_cache(G_sig);
    for (arma::uword a = 0; a < G_sig; ++a) {
      X_t_cache[a] = Rcpp::as<arma::mat>(X_list_r[a]).t();
    }

    for (arma::uword g = 0; g < P; ++g) {
      const arma::uword ag = static_cast<arma::uword>(std::max(1, sig_idx[g]) - 1);
      if (ag >= G_sig || g >= E_r.n_cols || m == 0) continue;
      const double* eg = E_r.colptr(g);
      const arma::vec yr = extract_y_vec_impl(eg, eg, m, adj_idx0, lag_pairs0, L_max);
      Xty_diag.col(g) += X_t_cache[ag] * yr;
    }
  }

  return Xty_diag;
}

// [[Rcpp::export]]
Rcpp::List poet_fit_pair_batch_cpp(
    const Rcpp::List& E_runs,
    const Rcpp::List& adj_idx1_runs,
    const Rcpp::List& lag_pairs_runs,
    const Rcpp::List& slot_XtX_list,
    const Rcpp::List& slot_X_runs_list,
    const arma::mat& A_tie_use,
    int L_max_in,
    int L_act_in,
    int M_dyn_in,
    int M_core_in,
    const Rcpp::IntegerVector& vc_col_idx,
    const arma::vec& sig0_sd,
    const arma::mat& R_poet,
    const Rcpp::IntegerVector& pair_g,
    const Rcpp::IntegerVector& pair_h,
    const Rcpp::IntegerVector& pair_slot,
    const Rcpp::LogicalVector& pair_do_hac,
    const Rcpp::LogicalVector& pair_has_rhs,
    const arma::mat& rhs_override
) {
  const arma::uword n_runs = static_cast<arma::uword>(E_runs.size());
  const arma::uword L_max = static_cast<arma::uword>(std::max(1, L_max_in));
  const arma::uword L_act = static_cast<arma::uword>(L_act_in);
  const arma::uword M_dyn = static_cast<arma::uword>(M_dyn_in);
  const arma::uword M_core = static_cast<arma::uword>(M_core_in);
  const arma::uword n_p = L_act + M_dyn;
  const arma::uword B_tot = A_tie_use.n_rows;
  const arma::uword n_vc = static_cast<arma::uword>(vc_col_idx.size());
  const arma::uword n_slots = static_cast<arma::uword>(slot_XtX_list.size());
  const arma::uword n_pairs = static_cast<arma::uword>(pair_g.size());

  std::vector<arma::mat> E_vec(n_runs);
  std::vector<arma::uvec> adj_vec(n_runs);
  std::vector<std::vector<arma::umat>> lp_vec(n_runs);
  for (arma::uword r = 0; r < n_runs; ++r) {
    E_vec[r] = Rcpp::as<arma::mat>(E_runs[r]);
    adj_vec[r] = to_zero_based_uvec(Rcpp::as<Rcpp::IntegerVector>(adj_idx1_runs[r]));
    lp_vec[r] = to_zero_based_lag_pairs(Rcpp::as<Rcpp::List>(lag_pairs_runs[r]), L_max);
  }

  arma::uvec idx_sc(L_act + M_core);
  for (arma::uword i = 0; i < L_act + M_core; ++i) idx_sc(i) = i;
  arma::uvec idx_vc(n_vc);
  for (arma::uword j = 0; j < n_vc; ++j) {
    idx_vc(j) = L_act + static_cast<arma::uword>(std::max(1, vc_col_idx[j]) - 1);
  }

  struct SlotCache {
    bool has_act = false;
    arma::uword n_act = 0;
    arma::mat XtX_t_inv;
    arma::mat R_sd;
    arma::mat XtX_stat_inv;
    arma::mat inv_sc;
    arma::mat R_cv;
    std::vector<arma::mat> X_runs;
    std::vector<arma::mat> X_runs_t;
    std::vector<arma::mat> X_rt_runs;
  };

  std::vector<SlotCache> slots(n_slots);
  for (arma::uword s = 0; s < n_slots; ++s) {
    SlotCache& sc = slots[s];
    const arma::mat XtX = Rcpp::as<arma::mat>(slot_XtX_list[s]);
    const arma::mat XtX_t = A_tie_use.t() * XtX * A_tie_use;
    sc.XtX_stat_inv = damped_sym_inv(XtX.submat(0, 0, L_act - 1, L_act - 1), 1e-8);

    const arma::uvec act_k = arma::find(XtX_t.diag() > 1e-8);
    sc.n_act = act_k.n_elem;
    sc.has_act = (sc.n_act > 0);
    sc.XtX_t_inv.zeros(n_p, n_p);
    sc.R_sd.zeros(L_act, M_dyn);
    sc.R_cv.zeros(M_core, n_vc);

    if (sc.has_act) {
      sc.XtX_t_inv.submat(act_k, act_k) = damped_sym_inv(XtX_t.submat(act_k, act_k), 1e-8);
      if (M_dyn > 0) {
        const arma::mat XtX_ss_inv = damped_sym_inv(XtX_t.submat(0, 0, L_act - 1, L_act - 1), 1e-8);
        sc.R_sd = XtX_ss_inv * XtX_t.submat(0, L_act, L_act - 1, n_p - 1);
      }

      if (n_vc >= 1 && M_core > 0) {
        sc.inv_sc = damped_sym_inv(XtX_t.submat(idx_sc, idx_sc), 1e-8);
        const arma::mat R_sc_vc = sc.inv_sc * XtX_t.submat(idx_sc, idx_vc);
        sc.R_cv = R_sc_vc.rows(L_act, L_act + M_core - 1);
      }
    }

    const Rcpp::List X_r_list = Rcpp::as<Rcpp::List>(slot_X_runs_list[s]);
    sc.X_runs.resize(n_runs);
    sc.X_runs_t.resize(n_runs);
    sc.X_rt_runs.resize(n_runs);
    for (arma::uword r = 0; r < n_runs; ++r) {
      sc.X_runs[r] = Rcpp::as<arma::mat>(X_r_list[r]);
      sc.X_runs_t[r] = sc.X_runs[r].t();
      sc.X_rt_runs[r] = sc.X_runs[r] * A_tie_use;
    }
  }

  arma::mat b_dyn_mat(M_dyn, n_pairs, arma::fill::zeros);
  arma::mat v_dyn_mat(M_dyn, n_pairs, arma::fill::zeros);
  arma::mat v_gauss_mat(M_dyn, n_pairs, arma::fill::zeros);
  arma::mat gam_unc_mat(L_act, n_pairs, arma::fill::zeros);
  arma::mat R_sd_flat(L_act * M_dyn, n_pairs, arma::fill::zeros);
  arma::mat R_cv_flat(M_core * n_vc, n_pairs, arma::fill::zeros);
  arma::mat b_core_al_mat(M_core, n_pairs, arma::fill::zeros);
  arma::mat v_core_al_mat(M_core, n_pairs, arma::fill::zeros);

  const arma::mat A_tie_t = A_tie_use.t();

  for (arma::uword p = 0; p < n_pairs; ++p) {
    const arma::uword g = static_cast<arma::uword>(std::max(1, pair_g[p]) - 1);
    const arma::uword h = static_cast<arma::uword>(std::max(1, pair_h[p]) - 1);
    const arma::uword s = static_cast<arma::uword>(std::max(1, pair_slot[p]) - 1);
    if (s >= n_slots) continue;
    const bool do_hac = (pair_do_hac[p] == TRUE);
    const bool has_rhs = (pair_has_rhs[p] == TRUE);
    const SlotCache& sc = slots[s];

    arma::vec rhs_full(B_tot, arma::fill::zeros);
    if (has_rhs) {
      rhs_full = rhs_override.col(p);
    } else {
      for (arma::uword r = 0; r < n_runs; ++r) {
        if (g >= E_vec[r].n_cols || h >= E_vec[r].n_cols) continue;
        const arma::uword m = E_vec[r].n_rows;
        const arma::vec yr = extract_y_vec_impl(
          E_vec[r].colptr(g), E_vec[r].colptr(h), m, adj_vec[r], lp_vec[r], L_max
        );
        rhs_full += sc.X_runs_t[r] * yr;
      }
    }

    gam_unc_mat.col(p) = sc.XtX_stat_inv * rhs_full.head(L_act);
    if (!sc.has_act || M_dyn == 0) continue;

    const arma::vec rhs_t = A_tie_t * rhs_full;
    const arma::vec b_t = sc.XtX_t_inv * rhs_t;
    b_dyn_mat.col(p) = b_t.subvec(L_act, n_p - 1);
    R_sd_flat.col(p) = arma::vectorise(sc.R_sd);

    if (M_core > 0) {
      arma::vec b_core_al = b_t.subvec(L_act, L_act + M_core - 1);
      if (n_vc >= 1) {
        R_cv_flat.col(p) = arma::vectorise(sc.R_cv);
        const arma::vec b_sc = sc.inv_sc * rhs_t.head(L_act + M_core);
        b_core_al = b_sc.subvec(L_act, L_act + M_core - 1);
      }
      b_core_al_mat.col(p) = b_core_al;
    }

    const double r_gh = R_poet(g, h);
    const double v_scale = (1.0 + r_gh * r_gh) * (sig0_sd(g) * sig0_sd(g)) * (sig0_sd(h) * sig0_sd(h));
    for (arma::uword j = 0; j < M_dyn; ++j) {
      v_gauss_mat(j, p) = v_scale * sc.XtX_t_inv(L_act + j, L_act + j);
    }

    if (!do_hac || sc.n_act <= L_act) {
      v_dyn_mat.col(p) = v_gauss_mat.col(p);
      if (n_vc >= 1) {
        for (arma::uword jc = 0; jc < M_core; ++jc) {
          v_core_al_mat(jc, p) = v_scale * sc.inv_sc(L_act + jc, L_act + jc);
        }
      } else {
        v_core_al_mat.col(p) = v_gauss_mat.col(p).head(M_core);
      }
      continue;
    }

    const arma::vec th_full = A_tie_use * b_t;
    arma::mat meat_t(n_p, n_p, arma::fill::zeros);
    for (arma::uword r = 0; r < n_runs; ++r) {
      const arma::uword m = E_vec[r].n_rows;
      const arma::mat& X_r = sc.X_runs[r];
      const arma::mat& X_rt = sc.X_rt_runs[r];
      const arma::vec y_r = extract_y_vec_impl(
        E_vec[r].colptr(g), E_vec[r].colptr(h), m, adj_vec[r], lp_vec[r], L_max
      );
      const arma::vec u_r = y_r - X_r * th_full;

      arma::mat S_t = X_rt.rows(0, m - 1);
      S_t.each_col() %= u_r.head(m);
      if (m >= 2) {
        arma::mat S_l1 = X_rt.rows(m, 2 * m - 2);
        S_l1.each_col() %= u_r.subvec(m, 2 * m - 2);
        S_t.rows(1, m - 1) += S_l1;
      }
      meat_t += S_t.t() * S_t;

      if (X_rt.n_rows >= 2 * m) {
        for (arma::uword rr = 2 * m - 1; rr < X_rt.n_rows; ++rr) {
          const arma::rowvec xr_u = X_rt.row(rr) * u_r(rr);
          meat_t += xr_u.t() * xr_u;
        }
      }

      for (arma::uword lh = 1; lh <= 3; ++lh) {
        if (m > lh + 2) {
          const double wb_h = 1.0 - static_cast<double>(lh) / 4.0;
          const arma::mat Gm = S_t.rows(lh, m - 1).t() * S_t.rows(0, m - lh - 1);
          meat_t += wb_h * (Gm + Gm.t());
        }
      }
    }

    const arma::mat V_t = sc.XtX_t_inv * meat_t * sc.XtX_t_inv;
    for (arma::uword j = 0; j < M_dyn; ++j) {
      v_dyn_mat(j, p) = V_t(L_act + j, L_act + j);
    }
    if (n_vc >= 1) {
      const arma::mat V_sc = sc.inv_sc * meat_t.submat(idx_sc, idx_sc) * sc.inv_sc;
      for (arma::uword jc = 0; jc < M_core; ++jc) {
        v_core_al_mat(jc, p) = V_sc(L_act + jc, L_act + jc);
      }
    } else {
      v_core_al_mat.col(p) = v_dyn_mat.col(p).head(M_core);
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("b_dyn") = b_dyn_mat,
    Rcpp::Named("v_dyn") = v_dyn_mat,
    Rcpp::Named("v_gauss") = v_gauss_mat,
    Rcpp::Named("gam_unc") = gam_unc_mat,
    Rcpp::Named("R_sd_flat") = R_sd_flat,
    Rcpp::Named("R_cv_flat") = R_cv_flat,
    Rcpp::Named("b_core_alone") = b_core_al_mat,
    Rcpp::Named("v_core_alone") = v_core_al_mat
  );
}
