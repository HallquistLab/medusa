#include "medusa.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

namespace {

struct CsrMat {
  arma::uword n_rows = 0;
  arma::uword n_cols = 0;
  std::vector<int> row_ptr;
  std::vector<int> col_idx;
  std::vector<double> val;
};

inline CsrMat parse_csr_components(
    const Rcpp::IntegerVector& row_ptr_r,
    const Rcpp::IntegerVector& col_idx_r,
    const Rcpp::NumericVector& val_r,
    arma::uword n_rows,
    arma::uword n_cols
) {
  CsrMat out;
  out.n_rows = n_rows;
  out.n_cols = n_cols;
  out.row_ptr.assign(row_ptr_r.begin(), row_ptr_r.end());
  out.col_idx.assign(col_idx_r.begin(), col_idx_r.end());
  out.val.assign(val_r.begin(), val_r.end());
  if (out.row_ptr.size() < n_rows + 1) {
    out.row_ptr.resize(n_rows + 1, static_cast<int>(out.val.size()));
  }
  return out;
}

inline arma::mat csr_mult_dense(const CsrMat& W, const arma::mat& X) {
  arma::mat out(W.n_rows, X.n_cols, arma::fill::zeros);
  if (W.n_rows == 0 || X.n_cols == 0 || W.val.empty()) return out;
  for (arma::uword c = 0; c < X.n_cols; ++c) {
    const double* xc = X.colptr(c);
    double* yc = out.colptr(c);
    for (arma::uword i = 0; i < W.n_rows; ++i) {
      double acc = 0.0;
      const int p_end = W.row_ptr[i + 1];
      for (int p = W.row_ptr[i]; p < p_end; ++p) {
        acc += W.val[p] * xc[W.col_idx[p]];
      }
      yc[i] = acc;
    }
  }
  return out;
}

inline arma::mat csr_trans_mult_dense(const CsrMat& W, const arma::mat& X) {
  arma::mat out(W.n_cols, X.n_cols, arma::fill::zeros);
  if (W.n_rows == 0 || W.n_cols == 0 || X.n_cols == 0 || W.val.empty()) return out;
  const arma::uword r_lim = std::min(W.n_rows, X.n_rows);
  for (arma::uword c = 0; c < X.n_cols; ++c) {
    const double* xc = X.colptr(c);
    double* yc = out.colptr(c);
    for (arma::uword i = 0; i < r_lim; ++i) {
      const double xi = xc[i];
      if (xi == 0.0) continue;
      const int p_end = W.row_ptr[i + 1];
      for (int p = W.row_ptr[i]; p < p_end; ++p) {
        yc[W.col_idx[p]] += W.val[p] * xi;
      }
    }
  }
  return out;
}

inline void accum_csr_f_csr_entry(
    const CsrMat& W_A,
    const CsrMat& W_B,
    arma::uword i1,
    arma::uword i2,
    double weight,
    arma::uword obs_row,
    arma::uword col_obs_base,
    const std::vector<arma::uword>& act_j,
    const arma::mat& f_dyn_obs,
    arma::mat& X_out
) {
  if (i1 >= W_A.n_rows || i2 >= W_B.n_rows) return;
  int pa = W_A.row_ptr[i1];
  const int pa_end = W_A.row_ptr[i1 + 1];
  int pb = W_B.row_ptr[i2];
  const int pb_end = W_B.row_ptr[i2 + 1];
  if (pa >= pa_end || pb >= pb_end) return;
  if (W_A.col_idx[pa_end - 1] < W_B.col_idx[pb] ||
      W_B.col_idx[pb_end - 1] < W_A.col_idx[pa]) {
    return;
  }
  const arma::uword nj = act_j.size();
  while (pa < pa_end && pb < pb_end) {
    const int ca = W_A.col_idx[pa];
    const int cb = W_B.col_idx[pb];
    if (ca == cb) {
      const double prod = weight * W_A.val[pa] * W_B.val[pb];
      const arma::uword t = static_cast<arma::uword>(ca);
      for (arma::uword jj = 0; jj < nj; ++jj) {
        const arma::uword j = act_j[jj];
        X_out(obs_row, col_obs_base + j) += prod * f_dyn_obs(t, j);
      }
      ++pa;
      ++pb;
    } else if (ca < cb) {
      ++pa;
    } else {
      ++pb;
    }
  }
}

struct StratumFactors {
  arma::uword m = 0;
  arma::mat Q;
  arma::mat UG;
  arma::mat VG;
  arma::mat UH;
  arma::mat VH;
  CsrMat W_b;
  bool has_woodbury = false;
};

inline StratumFactors parse_stratum_factors(const Rcpp::List& fac_r) {
  StratumFactors sf;
  sf.m = static_cast<arma::uword>(Rcpp::as<int>(fac_r["m"]));
  Rcpp::NumericMatrix q_nm = Rcpp::as<Rcpp::NumericMatrix>(fac_r["Q"]);
  sf.Q = arma::mat(q_nm.begin(), static_cast<arma::uword>(q_nm.nrow()),
                   static_cast<arma::uword>(q_nm.ncol()), false, true);

  sf.W_b = parse_csr_components(
    Rcpp::as<Rcpp::IntegerVector>(fac_r["W_row_ptr"]),
    Rcpp::as<Rcpp::IntegerVector>(fac_r["W_col_idx"]),
    Rcpp::as<Rcpp::NumericVector>(fac_r["W_val"]),
    sf.m, sf.m
  );

  arma::mat VH_Q;
  if (fac_r.containsElementNamed("VH_Q") && !Rf_isNull(fac_r["VH_Q"])) {
    Rcpp::NumericMatrix vh_nm = Rcpp::as<Rcpp::NumericMatrix>(fac_r["VH_Q"]);
    VH_Q = arma::mat(vh_nm.begin(), static_cast<arma::uword>(vh_nm.nrow()),
                     static_cast<arma::uword>(vh_nm.ncol()), false, true);
  } else if (sf.Q.n_cols > 0) {
    VH_Q = csr_trans_mult_dense(sf.W_b, sf.Q);
  } else {
    VH_Q.zeros(sf.m, 0);
  }

  sf.has_woodbury = fac_r.containsElementNamed("has_woodbury") &&
                    Rcpp::as<bool>(fac_r["has_woodbury"]);
  if (sf.has_woodbury) {
    const arma::mat u_we = Rcpp::as<arma::mat>(fac_r["u_we"]);
    const arma::mat MU_we = Rcpp::as<arma::mat>(fac_r["MU_we"]);
    const arma::vec fac_we = Rcpp::as<arma::vec>(fac_r["fac_we"]);
    const arma::mat P_we = Rcpp::as<arma::mat>(fac_r["P_we"]);
    const arma::uword r_q = sf.Q.n_cols;
    const arma::uword k_w = u_we.n_cols;

    arma::mat neg_MU_scaled = -MU_we;
    neg_MU_scaled.each_row() %= fac_we.t();

    sf.UG.set_size(sf.m, r_q + k_w);
    sf.VG.set_size(sf.m, r_q + k_w);
    sf.UH.set_size(sf.m, r_q + k_w);
    sf.VH.set_size(sf.m, r_q + k_w);
    if (r_q > 0) {
      sf.UG.cols(0, r_q - 1) = sf.Q;
      sf.VG.cols(0, r_q - 1) = sf.Q;
      sf.UH.cols(0, r_q - 1) = sf.Q;
      sf.VH.cols(0, r_q - 1) = VH_Q;
    }
    if (k_w > 0) {
      sf.UG.cols(r_q, r_q + k_w - 1) = neg_MU_scaled;
      sf.VG.cols(r_q, r_q + k_w - 1) = u_we;
      sf.UH.cols(r_q, r_q + k_w - 1) = -MU_we;
      sf.VH.cols(r_q, r_q + k_w - 1) = P_we;
    }
  } else {
    sf.VH = VH_Q;
  }
  return sf;
}

inline std::vector<arma::mat> build_Jt_U_list(
    const arma::mat& U,
    const std::vector<arma::umat>& lag_pairs0,
    arma::uword m,
    arma::uword L_max,
    bool skip_diag
) {
  std::vector<arma::mat> out(L_max);
  for (arma::uword d = 1; d <= L_max; ++d) {
    const arma::umat& pd = lag_pairs0[d - 1];
    arma::mat Jt_U(m, U.n_cols, arma::fill::zeros);
    if (U.n_cols > 0 && pd.n_rows > 0) {
      for (arma::uword r = 0; r < pd.n_rows; ++r) {
        const arma::uword i1 = pd(r, 0);
        const arma::uword i2 = pd(r, 1);
        if (i1 < m && i2 < m && (!skip_diag || i1 != i2)) {
          Jt_U.row(i2) += U.row(i1);
        }
      }
    }
    out[d - 1] = std::move(Jt_U);
  }
  return out;
}

inline std::vector<arma::mat> build_J_V_list(
    const arma::mat& V,
    const std::vector<arma::umat>& lag_pairs0,
    arma::uword m,
    arma::uword L_max,
    bool skip_diag
) {
  std::vector<arma::mat> out(L_max);
  for (arma::uword d = 1; d <= L_max; ++d) {
    const arma::umat& pd = lag_pairs0[d - 1];
    arma::mat J_V(m, V.n_cols, arma::fill::zeros);
    if (V.n_cols > 0 && pd.n_rows > 0) {
      for (arma::uword r = 0; r < pd.n_rows; ++r) {
        const arma::uword i1 = pd(r, 0);
        const arma::uword i2 = pd(r, 1);
        if (i1 < m && i2 < m && (!skip_diag || i1 != i2)) {
          J_V.row(i1) += V.row(i2);
        }
      }
    }
    out[d - 1] = std::move(J_V);
  }
  return out;
}

// Subtracts extract_obs_vec(X * Y^T) into col_ptr using precomputed J_{d'}^T X (when mode_JtX = true)
// or J_{d'} Y (when mode_JtX = false) for d' = 2..L_max.
inline void subtract_uv_obs(
    double* col_ptr,
    const arma::mat& X,
    const arma::mat& Y,
    const std::vector<arma::mat>& J_cache,
    bool cache_is_JtX,
    const arma::uvec& adj_idx0,
    const std::vector<arma::umat>& lag_pairs0,
    arma::uword m,
    arma::uword L_max
) {
  const arma::uword k = X.n_cols;
  if (k == 0 || Y.n_cols != k || m == 0) return;

  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
  const arma::uword n_adj = adj_idx0.n_elem;

  for (arma::uword c = 0; c < k; ++c) {
    const double* xc = X.colptr(c);
    const double* yc = Y.colptr(c);
    for (arma::uword i = 0; i < m; ++i) {
      col_ptr[i] -= xc[i] * yc[i];
    }
    for (arma::uword idx = 0; idx < n_adj; ++idx) {
      const arma::uword i0 = adj_idx0(idx);
      if (i0 + 1 < m) {
        col_ptr[m + i0] -= (xc[i0 + 1] * yc[i0] + xc[i0] * yc[i0 + 1]) * inv_sqrt2;
      }
    }
  }

  if (L_max >= 2) {
    for (arma::uword dp = 2; dp <= L_max; ++dp) {
      const arma::umat& pd = lag_pairs0[dp - 1];
      const double N_d = 0.5 * static_cast<double>(pd.n_rows);
      if (N_d > 0.0) {
        const double acc = cache_is_JtX
          ? arma::dot(J_cache[dp - 1], Y)
          : arma::dot(X, J_cache[dp - 1]);
        col_ptr[2 * m - 1 + (dp - 2)] -= acc / std::sqrt(2.0 * N_d);
      }
    }
  }
}

// Selective extraction of extract_obs_vec(A * B_right) from A^T (m x m) and B_right (m x m)
// without forming the m x m product matrix C = A * B_right.
inline arma::vec extract_prod_obs_selective(
    const arma::mat& A_t,
    const arma::mat& B_right,
    const arma::uvec& adj_idx0,
    const std::vector<arma::umat>& lag_pairs0,
    arma::uword L_max
) {
  const arma::uword m = A_t.n_cols;
  if (m == 0) return arma::vec();
  const arma::uword n_hi = (L_max >= 2) ? (L_max - 1) : 0;
  const arma::uword n_obs = 2 * m - 1 + n_hi;
  arma::vec out(n_obs, arma::fill::zeros);

  auto dot_col = [m](const double* p1, const double* p2) -> double {
    double s = 0.0;
    for (arma::uword k = 0; k < m; ++k) s += p1[k] * p2[k];
    return s;
  };

  for (arma::uword i = 0; i < m; ++i) {
    out(i) = dot_col(A_t.colptr(i), B_right.colptr(i));
  }

  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
  for (arma::uword k = 0; k < adj_idx0.n_elem; ++k) {
    const arma::uword i0 = adj_idx0(k);
    if (i0 + 1 < m) {
      const double c_10 = dot_col(A_t.colptr(i0 + 1), B_right.colptr(i0));
      const double c_01 = dot_col(A_t.colptr(i0), B_right.colptr(i0 + 1));
      out(m + i0) = (c_10 + c_01) * inv_sqrt2;
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
            acc += dot_col(A_t.colptr(i1), B_right.colptr(i2));
          }
        }
        out(2 * m - 1 + (d - 2)) = acc / std::sqrt(2.0 * N_d);
      }
    }
  }
  return out;
}

// Selective extraction of extract_obs_vec(A * diag(fk) * B^T) from A^T (m x m), B^T (m x m), and fk (m).
inline arma::vec extract_weighted_obs_selective(
    const arma::mat& A_t,
    const arma::mat& B_t,
    const arma::vec& fk,
    const arma::uvec& adj_idx0,
    const std::vector<arma::umat>& lag_pairs0,
    arma::uword L_max
) {
  const arma::uword m = A_t.n_cols;
  if (m == 0) return arma::vec();
  const arma::uword n_hi = (L_max >= 2) ? (L_max - 1) : 0;
  const arma::uword n_obs = 2 * m - 1 + n_hi;
  arma::vec out(n_obs, arma::fill::zeros);
  const double* fp = fk.memptr();

  auto wdot_col = [m, fp](const double* p1, const double* p2) -> double {
    double s = 0.0;
    for (arma::uword k = 0; k < m; ++k) s += p1[k] * fp[k] * p2[k];
    return s;
  };

  for (arma::uword i = 0; i < m; ++i) {
    out(i) = wdot_col(A_t.colptr(i), B_t.colptr(i));
  }

  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
  for (arma::uword k = 0; k < adj_idx0.n_elem; ++k) {
    const arma::uword i0 = adj_idx0(k);
    if (i0 + 1 < m) {
      const double c_10 = wdot_col(A_t.colptr(i0 + 1), B_t.colptr(i0));
      const double c_01 = wdot_col(A_t.colptr(i0), B_t.colptr(i0 + 1));
      out(m + i0) = (c_10 + c_01) * inv_sqrt2;
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
            acc += wdot_col(A_t.colptr(i1), B_t.colptr(i2));
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

inline std::vector<uint64_t> encode_sorted_pair_keys(
    const arma::umat& pd,
    arma::uword m,
    bool skip_diag
) {
  std::vector<uint64_t> keys;
  keys.reserve(pd.n_rows);
  for (arma::uword r = 0; r < pd.n_rows; ++r) {
    const arma::uword i1 = pd(r, 0);
    const arma::uword i2 = pd(r, 1);
    if (i1 < m && i2 < m && (!skip_diag || i1 != i2)) {
      keys.push_back((static_cast<uint64_t>(i1) << 32) | static_cast<uint64_t>(i2));
    }
  }
  std::sort(keys.begin(), keys.end());
  return keys;
}

inline double count_sorted_key_overlap(
    const std::vector<uint64_t>& a,
    const std::vector<uint64_t>& b
) {
  std::size_t ia = 0, ib = 0;
  double cnt = 0.0;
  while (ia < a.size() && ib < b.size()) {
    if (a[ia] == b[ib]) {
      const uint64_t kv = a[ia];
      std::size_t ca = 0, cb = 0;
      while (ia < a.size() && a[ia] == kv) { ++ca; ++ia; }
      while (ib < b.size() && b[ib] == kv) { ++cb; ++ib; }
      cnt += static_cast<double>(ca) * static_cast<double>(cb);
    } else if (a[ia] < b[ib]) {
      ++ia;
    } else {
      ++ib;
    }
  }
  return cnt;
}

} // namespace

// [[Rcpp::export]]
arma::mat poet_csr_transpose_mult_cpp(
    const Rcpp::IntegerVector& row_ptr,
    const Rcpp::IntegerVector& col_idx,
    const Rcpp::NumericVector& val,
    int n_rows,
    int n_cols,
    const arma::mat& X
) {
  const CsrMat W = parse_csr_components(
    row_ptr, col_idx, val,
    static_cast<arma::uword>(std::max(0, n_rows)),
    static_cast<arma::uword>(std::max(0, n_cols))
  );
  return csr_trans_mult_dense(W, X);
}

// [[Rcpp::export]]
arma::mat poet_tr_MM_factored_cpp(
    const Rcpp::List& Q_list,
    int m_in
) {
  const arma::uword S = static_cast<arma::uword>(Q_list.size());
  const arma::uword m = static_cast<arma::uword>(std::max(0, m_in));
  arma::mat tr_MM(S, S, arma::fill::zeros);
  if (S == 0 || m == 0) return tr_MM;

  std::vector<arma::mat> Q_vec(S);
  arma::uword max_r = 0;
  for (arma::uword s = 0; s < S; ++s) {
    Rcpp::NumericMatrix qm = Rcpp::as<Rcpp::NumericMatrix>(Q_list[s]);
    Q_vec[s] = arma::mat(qm.begin(), static_cast<arma::uword>(qm.nrow()),
                         static_cast<arma::uword>(qm.ncol()), false, true);
    if (Q_vec[s].n_cols > max_r) max_r = Q_vec[s].n_cols;
  }

  // Fast O(S^2 m r^2) direct trace path when S is small or rank r is small
  if (S <= 8 || (max_r * max_r * 2 <= m)) {
    arma::vec q_norm2(S, arma::fill::zeros);
    for (arma::uword s = 0; s < S; ++s) {
      q_norm2(s) = arma::accu(arma::square(Q_vec[s]));
    }
    for (arma::uword a = 0; a < S; ++a) {
      for (arma::uword b = a; b < S; ++b) {
        double cross_f2 = 0.0;
        if (Q_vec[a].n_cols > 0 && Q_vec[b].n_cols > 0) {
          const arma::mat QtQ = Q_vec[a].t() * Q_vec[b];
          cross_f2 = arma::accu(arma::square(QtQ));
        }
        const double val = static_cast<double>(m) - q_norm2(a) - q_norm2(b) + cross_f2;
        tr_MM(a, b) = val;
        tr_MM(b, a) = val;
      }
    }
    return tr_MM;
  }

  // Blocked upper-triangular BLAS dsyrk path: O(0.5 * m^2 * S^2) with <= ~112 MiB buffer
  const arma::uword row_blk = 32;
  const double sqrt2 = std::sqrt(2.0);
  const arma::uword max_pairs = row_blk * m;
  arma::mat B_chunk(max_pairs, S);

  for (arma::uword i0 = 0; i0 < m; i0 += row_blk) {
    const arma::uword i1 = std::min(m, i0 + row_blk);
    arma::uword n_pairs_blk = 0;
    for (arma::uword r = i0; r < i1; ++r) {
      n_pairs_blk += (m - r);
    }
    arma::subview<double> B_sub = B_chunk.rows(0, n_pairs_blk - 1);

    for (arma::uword s = 0; s < S; ++s) {
      double* col_p = B_sub.colptr(s);
      const arma::mat& Qs = Q_vec[s];
      arma::mat M_rows(i1 - i0, m, arma::fill::zeros);
      if (Qs.n_cols > 0) {
        M_rows = -(Qs.rows(i0, i1 - 1) * Qs.t());
      }
      for (arma::uword r = i0; r < i1; ++r) {
        M_rows(r - i0, r) += 1.0;
      }
      arma::uword cursor = 0;
      for (arma::uword r = i0; r < i1; ++r) {
        const arma::uword rr = r - i0;
        col_p[cursor++] = M_rows(rr, r);
        for (arma::uword c = r + 1; c < m; ++c) {
          col_p[cursor++] = sqrt2 * M_rows(rr, c);
        }
      }
    }
    tr_MM += B_sub.t() * B_sub;
  }

  return 0.5 * (tr_MM + tr_MM.t());
}

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

  const arma::mat M_At = M_A.t();
  X_out.col(0) = extract_prod_obs_selective(M_At, M_B, adj_idx0, lag_pairs0, L_max);

  arma::mat Jd_MB(m, M_B.n_cols);
  for (arma::uword d = 1; d <= L_max; ++d) {
    const arma::umat& pd = lag_pairs0[d - 1];
    Jd_MB.zeros();
    for (arma::uword r = 0; r < pd.n_rows; ++r) {
      const arma::uword i1 = pd(r, 0);
      const arma::uword i2 = pd(r, 1);
      if (i1 != i2 && i1 < m && i2 < M_B.n_rows) {
        Jd_MB.row(i1) += M_B.row(i2);
      }
    }
    X_out.col(d) = extract_prod_obs_selective(M_At, Jd_MB, adj_idx0, lag_pairs0, L_max);
  }

  if (M_dyn > 0) {
    const arma::mat G_At = G_A.t();
    const arma::mat G_Bt = G_B.t();
    const arma::mat H_At = H_A.t();
    const arma::mat H_Bt = H_B.t();
    for (arma::uword j = 0; j < M_dyn; ++j) {
      const arma::vec fk = f_dyn.col(j);
      if (arma::all(fk == 0.0)) continue;
      X_out.col(L1 + j) = extract_weighted_obs_selective(G_At, G_Bt, fk, adj_idx0, lag_pairs0, L_max);
      X_out.col(L1 + M_dyn + j) = extract_weighted_obs_selective(H_At, H_Bt, fk, adj_idx0, lag_pairs0, L_max);
    }
  }

  return X_out;
}

// [[Rcpp::export]]
arma::mat poet_build_pair_X_factored_cpp(
    const Rcpp::List& fac_A_r,
    const Rcpp::List& fac_B_r,
    const Rcpp::IntegerVector& adj_idx1,
    const Rcpp::List& lag_pairs,
    const arma::mat& f_dyn_inn,
    const arma::mat& f_dyn_obs,
    int L_max_in
) {
  const StratumFactors A = parse_stratum_factors(fac_A_r);
  const StratumFactors B = parse_stratum_factors(fac_B_r);
  const arma::uword m = A.m;
  const arma::uword L_max = static_cast<arma::uword>(std::max(1, L_max_in));
  const arma::uword L1 = L_max + 1;
  const arma::uword M_dyn = f_dyn_obs.n_cols;
  const arma::uword B_tot = L1 + 2 * M_dyn;
  if (m == 0) return arma::mat(0, B_tot);
  const arma::uword n_hi = (L_max >= 2) ? (L_max - 1) : 0;
  const arma::uword n_obs = 2 * m - 1 + n_hi;

  const arma::uvec adj_idx0 = to_zero_based_uvec(adj_idx1);
  const std::vector<arma::umat> lag_pairs0 = to_zero_based_lag_pairs(lag_pairs, L_max);
  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);

  arma::mat X_out(n_obs, B_tot, arma::fill::zeros);

  // J_{d'}^T Q_A and J_{d'} Q_B: Note stationary J_d skips diagonal (i1 != i2), while
  // extract_obs_vec sums all valid (i1, i2) in pd_{d'}. Build un-skipped caches for extract_obs_vec,
  // and skip_diag = true vectors W_d, Z_d for stationary operator J_d (identical whenever lag_pairs has i1 != i2).
  const std::vector<arma::mat> Jt_QA_obs = build_Jt_U_list(A.Q, lag_pairs0, m, L_max, false);
  const std::vector<arma::mat> J_QB_obs = build_J_V_list(B.Q, lag_pairs0, m, L_max, false);

  // Column 0: C0 = M_A * M_B = I_m - Q_A Q_A^T - (Q_B - Q_A (Q_A^T Q_B)) Q_B^T
  {
    double* c0 = X_out.colptr(0);
    for (arma::uword i = 0; i < m; ++i) c0[i] = 1.0;
    if (L_max >= 2) {
      for (arma::uword dp = 2; dp <= L_max; ++dp) {
        const arma::umat& pd = lag_pairs0[dp - 1];
        const double N_d = 0.5 * static_cast<double>(pd.n_rows);
        if (N_d > 0.0) {
          double diag_cnt = 0.0;
          for (arma::uword r = 0; r < pd.n_rows; ++r) {
            if (pd(r, 0) == pd(r, 1) && pd(r, 0) < m) diag_cnt += 1.0;
          }
          if (diag_cnt > 0.0) c0[2 * m - 1 + (dp - 2)] += diag_cnt / std::sqrt(2.0 * N_d);
        }
      }
    }
    subtract_uv_obs(c0, A.Q, A.Q, Jt_QA_obs, true, adj_idx0, lag_pairs0, m, L_max);
    arma::mat X2_0 = B.Q;
    if (A.Q.n_cols > 0 && B.Q.n_cols > 0) {
      X2_0 -= A.Q * (A.Q.t() * B.Q);
    }
    subtract_uv_obs(c0, X2_0, B.Q, J_QB_obs, false, adj_idx0, lag_pairs0, m, L_max);
  }

  // Columns 1..L_max: Cd = M_A * (J_d * M_B) = J_d - Q_A W_d^T - (Z_d - Q_A (Q_A^T Z_d)) Q_B^T
  std::vector<std::vector<uint64_t>> obs_keys(L_max);
  if (L_max >= 2) {
    for (arma::uword dp = 2; dp <= L_max; ++dp) {
      obs_keys[dp - 1] = encode_sorted_pair_keys(lag_pairs0[dp - 1], m, false);
    }
  }
  arma::vec jd_adj(m);
  for (arma::uword d = 1; d <= L_max; ++d) {
    const arma::umat& pd = lag_pairs0[d - 1];
    double* cd = X_out.colptr(d);
    jd_adj.zeros();
    arma::mat W_d(m, A.Q.n_cols, arma::fill::zeros);
    arma::mat Z_d(m, B.Q.n_cols, arma::fill::zeros);

    for (arma::uword r = 0; r < pd.n_rows; ++r) {
      const arma::uword i1 = pd(r, 0);
      const arma::uword i2 = pd(r, 1);
      if (i1 != i2 && i1 < m && i2 < m) {
        if (i1 == i2 + 1) jd_adj(i2) += 1.0;
        if (i2 == i1 + 1) jd_adj(i1) += 1.0;
        if (A.Q.n_cols > 0) W_d.row(i2) += A.Q.row(i1);
        if (B.Q.n_cols > 0) Z_d.row(i1) += B.Q.row(i2);
      }
    }

    for (arma::uword k = 0; k < adj_idx0.n_elem; ++k) {
      const arma::uword i0 = adj_idx0(k);
      if (i0 + 1 < m && jd_adj(i0) != 0.0) {
        cd[m + i0] += jd_adj(i0) * inv_sqrt2;
      }
    }
    if (L_max >= 2 && pd.n_rows > 0) {
      const std::vector<uint64_t> jd_keys = encode_sorted_pair_keys(pd, m, true);
      for (arma::uword dp = 2; dp <= L_max; ++dp) {
        const double N_dp = 0.5 * static_cast<double>(lag_pairs0[dp - 1].n_rows);
        if (N_dp > 0.0) {
          const double ov = count_sorted_key_overlap(jd_keys, obs_keys[dp - 1]);
          if (ov > 0.0) cd[2 * m - 1 + (dp - 2)] += ov / std::sqrt(2.0 * N_dp);
        }
      }
    }

    subtract_uv_obs(cd, A.Q, W_d, Jt_QA_obs, true, adj_idx0, lag_pairs0, m, L_max);
    if (A.Q.n_cols > 0 && Z_d.n_cols > 0) {
      Z_d -= A.Q * (A.Q.t() * Z_d);
    }
    subtract_uv_obs(cd, Z_d, B.Q, J_QB_obs, false, adj_idx0, lag_pairs0, m, L_max);
  }

  // Dynamic columns: innovation (L1..L1+M_dyn-1) and observation (L1+M_dyn..B_tot-1)
  if (M_dyn > 0) {
    std::vector<arma::uword> act_j_inn;
    std::vector<arma::uword> act_j_obs;
    act_j_inn.reserve(M_dyn);
    act_j_obs.reserve(M_dyn);
    for (arma::uword j = 0; j < M_dyn; ++j) {
      if (j < f_dyn_inn.n_cols && !arma::all(f_dyn_inn.col(j) == 0.0)) act_j_inn.push_back(j);
      if (!arma::all(f_dyn_obs.col(j) == 0.0)) act_j_obs.push_back(j);
    }

    // 1. Innovation: C_inn = G_A * diag(f_inn) * G_B^T
    if (!act_j_inn.empty()) {
      const arma::mat& UGA = A.has_woodbury ? A.UG : A.Q;
      const arma::mat& VGA = A.has_woodbury ? A.VG : A.Q;
      const arma::mat& UGB = B.has_woodbury ? B.UG : B.Q;
      const arma::mat& VGB = B.has_woodbury ? B.VG : B.Q;

      const std::vector<arma::mat> Jt_UGA_local = A.has_woodbury
        ? build_Jt_U_list(UGA, lag_pairs0, m, L_max, false)
        : std::vector<arma::mat>();
      const std::vector<arma::mat> J_UGB_local = B.has_woodbury
        ? build_J_V_list(UGB, lag_pairs0, m, L_max, false)
        : std::vector<arma::mat>();
      const std::vector<arma::mat>& Jt_UGA = A.has_woodbury ? Jt_UGA_local : Jt_QA_obs;
      const std::vector<arma::mat>& J_UGB = B.has_woodbury ? J_UGB_local : J_QB_obs;

      for (arma::uword jj = 0; jj < act_j_inn.size(); ++jj) {
        const arma::uword j = act_j_inn[jj];
        const arma::vec fk = f_dyn_inn.col(j);
        double* ci = X_out.colptr(L1 + j);
        for (arma::uword i = 0; i < m; ++i) ci[i] = fk(i);
        if (L_max >= 2) {
          for (arma::uword dp = 2; dp <= L_max; ++dp) {
            const arma::umat& pd = lag_pairs0[dp - 1];
            const double N_d = 0.5 * static_cast<double>(pd.n_rows);
            if (N_d > 0.0) {
              double acc_d = 0.0;
              for (arma::uword r = 0; r < pd.n_rows; ++r) {
                if (pd(r, 0) == pd(r, 1) && pd(r, 0) < m) acc_d += fk(pd(r, 0));
              }
              if (acc_d != 0.0) ci[2 * m - 1 + (dp - 2)] += acc_d / std::sqrt(2.0 * N_d);
            }
          }
        }

        arma::mat FA = VGA;
        if (FA.n_cols > 0) FA.each_col() %= fk;
        subtract_uv_obs(ci, UGA, FA, Jt_UGA, true, adj_idx0, lag_pairs0, m, L_max);

        arma::mat FB = VGB;
        if (FB.n_cols > 0) {
          FB.each_col() %= fk;
          if (UGA.n_cols > 0) {
            FB -= UGA * (VGA.t() * FB);
          }
        }
        subtract_uv_obs(ci, FB, UGB, J_UGB, false, adj_idx0, lag_pairs0, m, L_max);
      }
    }

    // 2. Observation: C_obs = H_A * diag(f_obs) * H_B^T, where H_s = W_b_s - U_H_s * V_H_s^T
    if (!act_j_obs.empty()) {
      const arma::uword col_obs_base = L1 + M_dyn;

      // Sparse base term W_A * diag(f_obs) * W_B^T evaluated across all active j simultaneously
      for (arma::uword i = 0; i < m; ++i) {
        accum_csr_f_csr_entry(A.W_b, B.W_b, i, i, 1.0, i, col_obs_base, act_j_obs, f_dyn_obs, X_out);
      }
      for (arma::uword k = 0; k < adj_idx0.n_elem; ++k) {
        const arma::uword i0 = adj_idx0(k);
        if (i0 + 1 < m) {
          accum_csr_f_csr_entry(A.W_b, B.W_b, i0 + 1, i0, inv_sqrt2, m + i0, col_obs_base, act_j_obs, f_dyn_obs, X_out);
          accum_csr_f_csr_entry(A.W_b, B.W_b, i0, i0 + 1, inv_sqrt2, m + i0, col_obs_base, act_j_obs, f_dyn_obs, X_out);
        }
      }
      if (L_max >= 2) {
        for (arma::uword dp = 2; dp <= L_max; ++dp) {
          const arma::umat& pd = lag_pairs0[dp - 1];
          const double N_d = 0.5 * static_cast<double>(pd.n_rows);
          if (N_d > 0.0) {
            const double wt_d = 1.0 / std::sqrt(2.0 * N_d);
            const arma::uword obs_r = 2 * m - 1 + (dp - 2);
            for (arma::uword r = 0; r < pd.n_rows; ++r) {
              const arma::uword i1 = pd(r, 0);
              const arma::uword i2 = pd(r, 1);
              if (i1 < m && i2 < m) {
                accum_csr_f_csr_entry(A.W_b, B.W_b, i1, i2, wt_d, obs_r, col_obs_base, act_j_obs, f_dyn_obs, X_out);
              }
            }
          }
        }
      }

      const arma::mat& UHA = A.has_woodbury ? A.UH : A.Q;
      const arma::mat& VHA = A.VH;
      const arma::mat& UHB = B.has_woodbury ? B.UH : B.Q;
      const arma::mat& VHB = B.VH;

      const std::vector<arma::mat> Jt_UHA_local = A.has_woodbury
        ? build_Jt_U_list(UHA, lag_pairs0, m, L_max, false)
        : std::vector<arma::mat>();
      const std::vector<arma::mat> J_UHB_local = B.has_woodbury
        ? build_J_V_list(UHB, lag_pairs0, m, L_max, false)
        : std::vector<arma::mat>();
      const std::vector<arma::mat>& Jt_UHA = A.has_woodbury ? Jt_UHA_local : Jt_QA_obs;
      const std::vector<arma::mat>& J_UHB = B.has_woodbury ? J_UHB_local : J_QB_obs;

      for (arma::uword jj = 0; jj < act_j_obs.size(); ++jj) {
        const arma::uword j = act_j_obs[jj];
        const arma::vec fk = f_dyn_obs.col(j);
        double* co = X_out.colptr(col_obs_base + j);

        if (UHA.n_cols > 0) {
          arma::mat EA = VHA;
          EA.each_col() %= fk;
          const arma::mat WB_EA = csr_mult_dense(B.W_b, EA);
          subtract_uv_obs(co, UHA, WB_EA, Jt_UHA, true, adj_idx0, lag_pairs0, m, L_max);
        }

        if (UHB.n_cols > 0) {
          arma::mat EB = VHB;
          EB.each_col() %= fk;
          arma::mat X2_obs = csr_mult_dense(A.W_b, EB);
          if (UHA.n_cols > 0) {
            X2_obs -= UHA * (VHA.t() * EB);
          }
          subtract_uv_obs(co, X2_obs, UHB, J_UHB, false, adj_idx0, lag_pairs0, m, L_max);
        }
      }
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
    Rcpp::NumericMatrix E_nm = Rcpp::as<Rcpp::NumericMatrix>(E_runs[r]);
    const arma::mat E_r(E_nm.begin(), static_cast<arma::uword>(E_nm.nrow()),
                        static_cast<arma::uword>(E_nm.ncol()), false, true);
    const Rcpp::List X_list_r = Rcpp::as<Rcpp::List>(X_diag_runs[r]);
    const arma::uvec adj_idx0 = to_zero_based_uvec(Rcpp::as<Rcpp::IntegerVector>(adj_idx1_runs[r]));
    const std::vector<arma::umat> lag_pairs0 = to_zero_based_lag_pairs(Rcpp::as<Rcpp::List>(lag_pairs_runs[r]), L_max);
    const arma::uword m = E_r.n_rows;

    const arma::uword G_sig = static_cast<arma::uword>(X_list_r.size());
    std::vector<arma::mat> X_views(G_sig);
    for (arma::uword a = 0; a < G_sig; ++a) {
      Rcpp::NumericMatrix X_nm = Rcpp::as<Rcpp::NumericMatrix>(X_list_r[a]);
      X_views[a] = arma::mat(X_nm.begin(), static_cast<arma::uword>(X_nm.nrow()),
                             static_cast<arma::uword>(X_nm.ncol()), false, true);
    }

    for (arma::uword g = 0; g < P; ++g) {
      const arma::uword ag = static_cast<arma::uword>(std::max(1, sig_idx[g]) - 1);
      if (ag >= G_sig || g >= E_r.n_cols || m == 0) continue;
      const double* eg = E_r.colptr(g);
      const arma::vec yr = extract_y_vec_impl(eg, eg, m, adj_idx0, lag_pairs0, L_max);
      Xty_diag.col(g) += X_views[ag].t() * yr;
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
    Rcpp::NumericMatrix E_nm = Rcpp::as<Rcpp::NumericMatrix>(E_runs[r]);
    E_vec[r] = arma::mat(E_nm.begin(), static_cast<arma::uword>(E_nm.nrow()),
                         static_cast<arma::uword>(E_nm.ncol()), false, true);
    adj_vec[r] = to_zero_based_uvec(Rcpp::as<Rcpp::IntegerVector>(adj_idx1_runs[r]));
    lp_vec[r] = to_zero_based_lag_pairs(Rcpp::as<Rcpp::List>(lag_pairs_runs[r]), L_max);
  }

  arma::uvec idx_sc(L_act + M_core);
  for (arma::uword i = 0; i < L_act + M_core; ++i) idx_sc(i) = i;
  arma::uvec idx_vc(n_vc);
  for (arma::uword j = 0; j < n_vc; ++j) {
    idx_vc(j) = L_act + static_cast<arma::uword>(std::max(1, vc_col_idx[j]) - 1);
  }

  // Determine per-slot usage so X_runs views are only bound when needed and X_rt_runs
  // is only materialized for slots referenced by at least one HAC pair.
  std::vector<bool> slot_needs_X(n_slots, false);
  std::vector<bool> slot_needs_hac(n_slots, false);
  for (arma::uword p = 0; p < n_pairs; ++p) {
    const arma::uword s = static_cast<arma::uword>(std::max(1, pair_slot[p]) - 1);
    if (s >= n_slots) continue;
    const bool do_hac = (pair_do_hac[p] == TRUE);
    const bool has_rhs = (pair_has_rhs[p] == TRUE);
    if (!has_rhs || do_hac) slot_needs_X[s] = true;
    if (do_hac) slot_needs_hac[s] = true;
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

    if (slot_needs_X[s] && static_cast<R_xlen_t>(s) < slot_X_runs_list.size() && !Rf_isNull(slot_X_runs_list[s])) {
      const Rcpp::List X_r_list = Rcpp::as<Rcpp::List>(slot_X_runs_list[s]);
      sc.X_runs.resize(n_runs);
      const bool build_rt = slot_needs_hac[s] && sc.has_act && (sc.n_act > L_act) && (M_dyn > 0);
      if (build_rt) sc.X_rt_runs.resize(n_runs);
      for (arma::uword r = 0; r < n_runs; ++r) {
        Rcpp::NumericMatrix X_nm = Rcpp::as<Rcpp::NumericMatrix>(X_r_list[r]);
        sc.X_runs[r] = arma::mat(X_nm.begin(), static_cast<arma::uword>(X_nm.nrow()),
                                 static_cast<arma::uword>(X_nm.ncol()), false, true);
        if (build_rt) {
          sc.X_rt_runs[r] = sc.X_runs[r] * A_tie_use;
        }
      }
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
        rhs_full += sc.X_runs[r].t() * yr;
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

// [[Rcpp::export]]
arma::mat poet_build_trial_Ci_cpp(
    const arma::mat& U_ti,
    const arma::mat& V_ti,
    const Rcpp::IntegerVector& map_s_idx1,
    const Rcpp::List& G_stat_sub,
    const Rcpp::List& B_dyn_sub,
    const Rcpp::List& lag_pairs,
    const arma::mat& f_dyn_r,
    const arma::vec& gm_vec,
    const arma::vec& run_scale_Cu,
    const arma::mat& s2_diag_mat,
    int n_lags_in,
    int n_sC_in,
    int L_max_in
) {
  const arma::uword nC = static_cast<arma::uword>(map_s_idx1.size());
  const arma::uword n_lags = static_cast<arma::uword>(std::max(1, n_lags_in));
  const arma::uword n_sC = static_cast<arma::uword>(std::max(1, n_sC_in));
  const arma::uword L_max = static_cast<arma::uword>(std::max(1, L_max_in));
  const arma::uword q_m = nC * n_lags;
  const arma::uword q_sC = n_sC * n_lags;
  const arma::uword m = U_ti.n_rows;

  arma::mat Ci(q_m, q_m, arma::fill::zeros);
  if (nC == 0 || m == 0) return Ci;

  arma::uvec map_s0 = to_zero_based_uvec(map_s_idx1);
  const std::vector<arma::umat> lag_pairs0 = to_zero_based_lag_pairs(lag_pairs, L_max);

  auto add_tiled_term = [&](const arma::mat& M_coeff, const arma::mat& T_s) {
    for (arma::uword lb = 0; lb < n_lags; ++lb) {
      for (arma::uword ib = 0; ib < nC; ++ib) {
        const arma::uword col_ci = lb * nC + ib;
        const arma::uword col_ts = lb * n_sC + map_s0(ib);
        const double* mc_col = M_coeff.colptr(ib);
        for (arma::uword la = 0; la < n_lags; ++la) {
          double* ci_blk = Ci.colptr(col_ci) + la * nC;
          const double* ts_blk = T_s.colptr(col_ts) + la * n_sC;
          for (arma::uword ia = 0; ia < nC; ++ia) {
            ci_blk[ia] += mc_col[ia] * ts_blk[map_s0(ia)];
          }
        }
      }
    }
  };

  arma::mat T_buf = U_ti.t() * U_ti;
  if (G_stat_sub.size() >= 1) {
    Rcpp::NumericMatrix g0_nm = Rcpp::as<Rcpp::NumericMatrix>(G_stat_sub[0]);
    const arma::mat G0(g0_nm.begin(), nC, nC, false, true);
    add_tiled_term(G0, T_buf);
  }

  for (arma::uword d = 1; d <= L_max; ++d) {
    if (static_cast<R_xlen_t>(d) >= G_stat_sub.size()) break;
    Rcpp::NumericMatrix gd_nm = Rcpp::as<Rcpp::NumericMatrix>(G_stat_sub[d]);
    const arma::mat Gd(gd_nm.begin(), nC, nC, false, true);
    if (!arma::any(arma::abs(arma::vectorise(Gd)) > 1e-14)) continue;
    const arma::umat& pd = lag_pairs0[d - 1];
    if (pd.n_rows > 0) {
      T_buf = U_ti.rows(pd.col(0)).t() * U_ti.rows(pd.col(1));
      add_tiled_term(Gd, T_buf);
    }
  }

  const arma::uword M_dyn = static_cast<arma::uword>(B_dyn_sub.size());
  for (arma::uword j = 0; j < M_dyn; ++j) {
    if (j >= f_dyn_r.n_cols || j >= gm_vec.n_elem) break;
    Rcpp::NumericMatrix bj_nm = Rcpp::as<Rcpp::NumericMatrix>(B_dyn_sub[j]);
    const arma::mat Bj(bj_nm.begin(), nC, nC, false, true);
    if (!arma::any(arma::abs(arma::vectorise(Bj)) > 1e-14)) continue;

    const arma::vec fj = f_dyn_r.col(j);
    const double gm_j = gm_vec(j);
    T_buf.zeros(q_sC, q_sC);
    if (gm_j < 1.0) {
      arma::mat U_sc = U_ti;
      U_sc.each_col() %= fj;
      T_buf += (1.0 - gm_j) * (U_ti.t() * U_sc);
    }
    if (gm_j > 0.0) {
      arma::mat V_sc = V_ti;
      V_sc.each_col() %= fj;
      T_buf += gm_j * (V_ti.t() * V_sc);
    }
    add_tiled_term(Bj, T_buf);
  }

  if (run_scale_Cu.n_elem == nC) {
    for (arma::uword lb = 0; lb < n_lags; ++lb) {
      for (arma::uword ib = 0; ib < nC; ++ib) {
        const double sc_b = run_scale_Cu(ib);
        double* ci_col = Ci.colptr(lb * nC + ib);
        for (arma::uword la = 0; la < n_lags; ++la) {
          double* ci_blk = ci_col + la * nC;
          for (arma::uword ia = 0; ia < nC; ++ia) {
            ci_blk[ia] *= run_scale_Cu(ia) * sc_b;
          }
        }
      }
    }
  }

  if (s2_diag_mat.n_rows == nC && s2_diag_mat.n_cols == n_lags) {
    for (arma::uword la = 0; la < n_lags; ++la) {
      for (arma::uword ia = 0; ia < nC; ++ia) {
        const arma::uword idx_d = la * nC + ia;
        Ci(idx_d, idx_d) = s2_diag_mat(ia, la);
      }
    }
  }

  return Ci;
}

namespace {

inline void csr_trans_mult_vec(const CsrMat& W, const double* x, double* y) {
  std::memset(y, 0, static_cast<std::size_t>(W.n_cols) * sizeof(double));
  if (W.n_rows == 0 || W.n_cols == 0 || W.val.empty()) return;
  for (arma::uword i = 0; i < W.n_rows; ++i) {
    const double xi = x[i];
    if (xi == 0.0) continue;
    const int p_end = W.row_ptr[i + 1];
    for (int p = W.row_ptr[i]; p < p_end; ++p) {
      y[W.col_idx[p]] += W.val[p] * xi;
    }
  }
}

inline void symmetrize_inplace(arma::mat& M) {
  const arma::uword n = M.n_rows;
  for (arma::uword j = 0; j < n; ++j) {
    double* cj = M.colptr(j);
    if (!std::isfinite(cj[j])) cj[j] = 0.0;
    for (arma::uword i = j + 1; i < n; ++i) {
      double a = cj[i];
      double b = M(j, i);
      if (!std::isfinite(a)) a = 0.0;
      if (!std::isfinite(b)) b = 0.0;
      const double avg = 0.5 * (a + b);
      cj[i] = avg;
      M(j, i) = avg;
    }
  }
}

inline arma::mat project_trial_err_psd_impl(const arma::mat& C_mat, double max_cor = 0.98) {
  const arma::uword k = C_mat.n_rows;
  arma::mat C_sym = C_mat;
  symmetrize_inplace(C_sym);
  arma::vec d_targ = arma::clamp(C_sym.diag(), 0.0, arma::datum::inf);
  arma::uvec idx = arma::find(d_targ > 1e-16);
  arma::mat C_out(k, k, arma::fill::zeros);
  if (idx.n_elem == 0) return C_out;
  if (idx.n_elem == 1) {
    C_out(idx(0), idx(0)) = d_targ(idx(0));
    return C_out;
  }

  const arma::uword np = idx.n_elem;
  const arma::vec sd_v = arma::sqrt(d_targ.elem(idx));
  arma::mat R_sub = C_sym.submat(idx, idx);
  for (arma::uword j = 0; j < np; ++j) {
    R_sub(j, j) = 1.0;
    for (arma::uword i = j + 1; i < np; ++i) {
      const double v = std::min(max_cor, std::max(-max_cor, R_sub(i, j) / (sd_v(i) * sd_v(j))));
      R_sub(i, j) = v;
      R_sub(j, i) = v;
    }
  }

  arma::vec eg_val;
  arma::mat eg_vec;
  arma::eig_sym(eg_val, eg_vec, R_sub);
  arma::vec ev_cl = arma::clamp(eg_val, 1e-4, arma::datum::inf);
  arma::mat V_sc = eg_vec;
  V_sc.each_row() %= ev_cl.t();
  arma::mat R_psd = V_sc * eg_vec.t();
  const double eps = std::numeric_limits<double>::epsilon();
  arma::vec d_r = arma::sqrt(arma::clamp(R_psd.diag(), eps, arma::datum::inf));
  for (arma::uword j = 0; j < np; ++j) {
    R_psd(j, j) = sd_v(j) * sd_v(j);
    for (arma::uword i = j + 1; i < np; ++i) {
      const double rc = std::min(max_cor, std::max(-max_cor, R_psd(i, j) / (d_r(i) * d_r(j))));
      const double cv = rc * sd_v(i) * sd_v(j);
      R_psd(i, j) = cv;
      R_psd(j, i) = cv;
    }
  }
  C_out.submat(idx, idx) = R_psd;
  return C_out;
}

// In-place Block-PAF factor extraction matching .extract_psd_factors() with minimal memory footprint.
// Note: mutates C_work in-place when np == q to avoid allocating extra q x q matrices.
inline arma::mat extract_psd_factors_inplace(
    arma::mat& C_work,
    arma::uword k_rank,
    arma::uword nC,
    arma::uword n_lags,
    const arma::mat* probe_ptr
) {
  const arma::uword q = C_work.n_rows;
  arma::mat F_out(k_rank, q, arma::fill::zeros);
  if (q == 0 || k_rank == 0) return F_out;

  symmetrize_inplace(C_work);
  const arma::vec d_ci = C_work.diag();
  const arma::uvec pos_c = arma::find(d_ci > 1e-14);
  const arma::uword np = pos_c.n_elem;
  if (np == 0) return F_out;
  if (np == 1) {
    F_out(0, pos_c(0)) = std::sqrt(d_ci(pos_c(0)));
    return F_out;
  }

  const arma::vec sd_c = arma::sqrt(d_ci.elem(pos_c));
  arma::mat R_sub;
  const bool use_c_inplace = (np == q);
  if (!use_c_inplace) {
    R_sub = C_work.submat(pos_c, pos_c);
  }
  arma::mat& R_work = use_c_inplace ? C_work : R_sub;

  for (arma::uword j = 0; j < np; ++j) {
    double* cj = R_work.colptr(j);
    cj[j] = 1.0;
    const double sj = sd_c(j);
    for (arma::uword i = j + 1; i < np; ++i) {
      const double rv = std::min(0.98, std::max(-0.98, cj[i] / (sd_c(i) * sj)));
      cj[i] = rv;
      R_work(j, i) = rv;
    }
  }

  const double eps = std::numeric_limits<double>::epsilon();

  if (np <= k_rank) {
    arma::vec eg_val;
    arma::mat eg_vec;
    arma::eig_sym(eg_val, eg_vec, R_work);
    if (eg_val.n_elem > 0 && eg_val(0) < -1e-8) {
      arma::vec ev_p = arma::clamp(eg_val, 0.0, arma::datum::inf);
      arma::mat V_sc = eg_vec;
      V_sc.each_row() %= ev_p.t();
      arma::mat R_p = V_sc * eg_vec.t();
      arma::vec dr_p = arma::sqrt(arma::clamp(R_p.diag(), eps, arma::datum::inf));
      for (arma::uword j = 0; j < np; ++j) {
        R_p(j, j) = 1.0;
        for (arma::uword i = j + 1; i < np; ++i) {
          const double v = R_p(i, j) / (dr_p(i) * dr_p(j));
          R_p(i, j) = v;
          R_p(j, i) = v;
        }
      }
      arma::eig_sym(eg_val, eg_vec, R_p);
    }
    arma::uword k_use = 0;
    for (arma::uword j = 0; j < eg_val.n_elem; ++j) {
      if (eg_val(eg_val.n_elem - 1 - j) > 1e-12 && k_use < k_rank) ++k_use;
    }
    for (arma::uword j = 0; j < k_use; ++j) {
      const arma::uword src_col = eg_val.n_elem - 1 - j;
      const double sc_ev = std::sqrt(eg_val(src_col));
      const double* v_ptr = eg_vec.colptr(src_col);
      for (arma::uword r = 0; r < np; ++r) {
        F_out(j, pos_c(r)) = v_ptr[r] * sc_ev * sd_c(r);
      }
    }
    return F_out;
  }

  arma::mat R_inv = R_work;
  R_inv.diag() += 1e-3;
  if (!arma::inv_sympd(R_inv, R_inv)) {
    arma::vec eg_val;
    arma::mat eg_vec;
    arma::eig_sym(eg_val, eg_vec, R_work);
    arma::vec ev_cl = arma::clamp(eg_val, 1e-4, arma::datum::inf);
    arma::mat V_sc = eg_vec;
    V_sc.each_row() %= ev_cl.t();
    R_work = V_sc * eg_vec.t();
    arma::vec dr_p = arma::sqrt(arma::clamp(R_work.diag(), eps, arma::datum::inf));
    for (arma::uword j = 0; j < np; ++j) {
      double* cj = R_work.colptr(j);
      cj[j] = 1.0;
      for (arma::uword i = j + 1; i < np; ++i) {
        const double rv = std::min(0.98, std::max(-0.98, cj[i] / (dr_p(i) * dr_p(j))));
        cj[i] = rv;
        R_work(j, i) = rv;
      }
    }
    arma::mat R_try = R_work;
    R_try.diag() += 1e-4;
    if (!arma::inv_sympd(R_inv, R_try)) {
      arma::eig_sym(eg_val, eg_vec, R_work);
      arma::vec inv_ev = 1.0 / arma::clamp(eg_val, 1e-4, arma::datum::inf);
      V_sc = eg_vec;
      V_sc.each_row() %= inv_ev.t();
      R_inv = V_sc * eg_vec.t();
    }
  }

  std::vector<arma::uvec> u_blocks;
  if (nC >= 2 && q == nC * n_lags) {
    std::vector<std::vector<arma::uword>> buckets(nC);
    for (arma::uword r = 0; r < np; ++r) {
      buckets[pos_c(r) % nC].push_back(r);
    }
    u_blocks.reserve(nC);
    for (arma::uword u = 0; u < nC; ++u) {
      if (!buckets[u].empty()) {
        u_blocks.emplace_back(buckets[u]);
      }
    }
  } else {
    u_blocks.resize(np);
    for (arma::uword r = 0; r < np; ++r) {
      u_blocks[r] = arma::uvec({r});
    }
  }

  for (std::size_t b = 0; b < u_blocks.size(); ++b) {
    const arma::uvec& idx_b = u_blocks[b];
    if (idx_b.n_elem == 1) {
      const arma::uword ib = idx_b(0);
      const double cond_b = 1.0 / std::max(R_inv(ib, ib), 1.0);
      R_work(ib, ib) = std::max(0.0, 1.0 - 0.98 * cond_b);
    } else {
      const arma::mat R_bb = R_inv.submat(idx_b, idx_b);
      const arma::mat Prec_bb = 0.5 * (R_bb + R_bb.t());
      const arma::mat Cond_bb = damped_sym_inv(Prec_bb, 1e-6);
      R_work.submat(idx_b, idx_b) -= 0.98 * Cond_bb;
    }
  }
  R_inv.reset();

  arma::mat Omega;
  if (np > 512) {
    const arma::uword r_sub = std::min(np, k_rank + static_cast<arma::uword>(12));
    if (probe_ptr != nullptr && probe_ptr->n_rows >= np && probe_ptr->n_cols >= r_sub) {
      Omega = probe_ptr->submat(0, 0, np - 1, r_sub - 1);
    } else {
      arma::mat raw_om(np, r_sub);
      for (arma::uword j = 0; j < r_sub; ++j) {
        const double jd = static_cast<double>(j + 1);
        double* col_o = raw_om.colptr(j);
        for (arma::uword i = 0; i < np; ++i) {
          const double id = static_cast<double>(i + 1);
          col_o[i] = std::sin((id * 1.6180339887 + jd) * jd);
        }
      }
      arma::mat R_qr;
      arma::qr_econ(Omega, R_qr, raw_om);
    }
  }

  arma::mat F_cor;
  arma::uword k_use = 0;
  for (int it = 0; it < 3; ++it) {
    symmetrize_inplace(R_work);
    if (np <= 512) {
      arma::vec eg_val;
      arma::mat eg_vec;
      arma::eig_sym(eg_val, eg_vec, R_work);
      k_use = 0;
      for (arma::uword j = 0; j < eg_val.n_elem; ++j) {
        if (eg_val(eg_val.n_elem - 1 - j) > 1e-6 && k_use < k_rank) ++k_use;
      }
      if (k_use < 1) break;
      F_cor.set_size(np, k_use);
      for (arma::uword j = 0; j < k_use; ++j) {
        const arma::uword src_col = eg_val.n_elem - 1 - j;
        F_cor.col(j) = eg_vec.col(src_col) * std::sqrt(eg_val(src_col));
      }
    } else {
      arma::mat Q_it, R_qr;
      arma::mat Y_it = R_work * Omega + Omega;
      arma::qr_econ(Q_it, R_qr, Y_it);
      Y_it = R_work * Q_it + Q_it;
      arma::qr_econ(Q_it, R_qr, Y_it);
      Y_it = R_work * Q_it + Q_it;
      arma::qr_econ(Q_it, R_qr, Y_it);
      Y_it = R_work * Q_it + Q_it;
      arma::mat B_sm = Q_it.t() * Y_it;
      B_sm = 0.5 * (B_sm + B_sm.t());
      arma::vec eg_val;
      arma::mat eg_vec;
      arma::eig_sym(eg_val, eg_vec, B_sm);
      k_use = 0;
      for (arma::uword j = 0; j < eg_val.n_elem; ++j) {
        const double vu = eg_val(eg_val.n_elem - 1 - j) - 1.0;
        if (vu > 1e-6 && k_use < k_rank) ++k_use;
      }
      if (k_use < 1) break;
      arma::mat evecs_top(eg_vec.n_rows, k_use);
      arma::rowvec scales_top(k_use);
      for (arma::uword j = 0; j < k_use; ++j) {
        const arma::uword src_col = eg_val.n_elem - 1 - j;
        evecs_top.col(j) = eg_vec.col(src_col);
        scales_top(j) = std::sqrt(eg_val(src_col) - 1.0);
      }
      F_cor = Q_it * evecs_top;
      F_cor.each_row() %= scales_top;
    }

    for (arma::uword r = 0; r < np; ++r) {
      double rn2 = 0.0;
      for (arma::uword j = 0; j < k_use; ++j) {
        const double v = F_cor(r, j);
        rn2 += v * v;
      }
      const double rn = std::sqrt(rn2);
      if (rn > 0.995) {
        const double mult = 0.995 / rn;
        for (arma::uword j = 0; j < k_use; ++j) F_cor(r, j) *= mult;
      }
    }

    if (it < 2) {
      for (std::size_t b = 0; b < u_blocks.size(); ++b) {
        const arma::uvec& idx_b = u_blocks[b];
        const arma::mat F_b = F_cor.rows(idx_b);
        R_work.submat(idx_b, idx_b) = F_b * F_b.t();
      }
    }
  }

  if (k_use >= 1) {
    for (arma::uword j = 0; j < k_use; ++j) {
      const double* fc_col = F_cor.colptr(j);
      for (arma::uword r = 0; r < np; ++r) {
        F_out(j, pos_c(r)) = fc_col[r] * sd_c(r);
      }
    }
  }

  return F_out;
}

struct LagSpan {
  arma::uword s0;
  arma::uword s1;
  arma::uword len;
};

inline std::vector<std::vector<LagSpan>> build_lag_spans(
    const std::vector<arma::umat>& lag_pairs0,
    arma::uword m
) {
  std::vector<std::vector<LagSpan>> out(lag_pairs0.size());
  for (std::size_t d = 0; d < lag_pairs0.size(); ++d) {
    const arma::umat& pd = lag_pairs0[d];
    const arma::uword n_p = pd.n_rows;
    if (n_p == 0) continue;
    const arma::uword* p0 = pd.colptr(0);
    const arma::uword* p1 = pd.colptr(1);
    std::vector<LagSpan>& sp = out[d];
    bool active = false;
    LagSpan cur{0, 0, 0};
    for (arma::uword r = 0; r < n_p; ++r) {
      const arma::uword i0 = p0[r];
      const arma::uword i1 = p1[r];
      if (i0 >= m || i1 >= m) continue;
      if (active && i0 == cur.s0 + cur.len && i1 == cur.s1 + cur.len) {
        ++cur.len;
      } else {
        if (active) sp.push_back(cur);
        cur = LagSpan{i0, i1, 1};
        active = true;
      }
    }
    if (active) sp.push_back(cur);
  }
  return out;
}

inline void eval_sem_basis_profiles(
    const arma::mat& UL,
    const arma::mat& UR,
    const arma::mat& VL,
    const arma::mat& VR,
    const std::vector<std::vector<LagSpan>>& lag_spans,
    const arma::mat& f_u,
    const arma::mat& f_v,
    bool any_fu,
    bool any_fv,
    bool need_mlag,
    arma::uword L_max,
    arma::mat& W_prod,
    arma::mat& T_stat,
    arma::mat& T_dyn
) {
  const arma::uword m = UL.n_rows;
  const arma::uword k_c = UL.n_cols;
  const arma::uword L1 = L_max + 1;
  const arma::uword M_dyn = f_u.n_cols;
  T_stat.zeros(k_c, L1);
  T_dyn.zeros(k_c, M_dyn);
  if (m == 0 || k_c == 0) return;

  W_prod.set_size(m, k_c);
  for (arma::uword c = 0; c < k_c; ++c) {
    const double* ul = UL.colptr(c);
    const double* ur = UR.colptr(c);
    double* wp = W_prod.colptr(c);
    double s0 = 0.0;
    for (arma::uword t = 0; t < m; ++t) {
      const double p = ul[t] * ur[t];
      wp[t] = p;
      s0 += p;
    }
    T_stat(c, 0) = s0;
  }

  if (need_mlag) {
    for (arma::uword d = 1; d <= L_max && (d - 1) < lag_spans.size(); ++d) {
      const std::vector<LagSpan>& spans_d = lag_spans[d - 1];
      if (spans_d.empty()) continue;
      for (arma::uword c = 0; c < k_c; ++c) {
        const double* ul = UL.colptr(c);
        const double* ur = UR.colptr(c);
        double acc = 0.0;
        for (std::size_t s = 0; s < spans_d.size(); ++s) {
          const double* a = ul + spans_d[s].s0;
          const double* b = ur + spans_d[s].s1;
          const arma::uword len = spans_d[s].len;
          for (arma::uword k = 0; k < len; ++k) acc += a[k] * b[k];
        }
        T_stat(c, d) = acc;
      }
    }
  }

  if (M_dyn > 0) {
    if (any_fu) {
      T_dyn = W_prod.t() * f_u;
    }
    const bool has_V = (VL.n_rows == m && VR.n_rows == m && VL.n_cols == k_c && VR.n_cols == k_c);
    if (any_fv && has_V) {
      for (arma::uword c = 0; c < k_c; ++c) {
        const double* vl = VL.colptr(c);
        const double* vr = VR.colptr(c);
        double* wp = W_prod.colptr(c);
        for (arma::uword t = 0; t < m; ++t) wp[t] = vl[t] * vr[t];
      }
      T_dyn += W_prod.t() * f_v;
    }
  }
}

inline arma::mat center_trial_op_cols(const arma::mat& U_in, const arma::uvec& ok_idx) {
  const arma::uword m = U_in.n_rows;
  const arma::uword m_t = U_in.n_cols;
  arma::mat U_out(m, m_t, arma::fill::zeros);
  const arma::uword m_ok = ok_idx.n_elem;
  if (m == 0 || m_ok == 0) return U_out;
  if (m_ok == 1) {
    U_out.col(ok_idx(0)) = U_in.col(ok_idx(0));
    return U_out;
  }
  arma::vec mu(m, arma::fill::zeros);
  for (arma::uword k = 0; k < m_ok; ++k) {
    mu += U_in.col(ok_idx(k));
  }
  mu /= static_cast<double>(m_ok);
  const double inv_sc = 1.0 / std::sqrt(1.0 - 1.0 / static_cast<double>(m_ok));
  for (arma::uword k = 0; k < m_ok; ++k) {
    const arma::uword c = ok_idx(k);
    U_out.col(c) = (U_in.col(c) - mu) * inv_sc;
  }
  return U_out;
}

} // namespace

// [[Rcpp::export]]
arma::cube poet_build_Xty_bulk_cpp(
    const Rcpp::List& E_runs,
    const Rcpp::List& X_rep_runs,
    const Rcpp::IntegerVector& rep_grp,
    const arma::mat& Xty_diag,
    const Rcpp::List& adj_idx1_runs,
    const Rcpp::List& lag_pairs_runs,
    int G_rep_in,
    int L_max_in,
    int B_tot_in
) {
  const arma::uword n_runs = static_cast<arma::uword>(E_runs.size());
  const arma::uword P = static_cast<arma::uword>(rep_grp.size());
  const arma::uword G_rep = static_cast<arma::uword>(std::max(1, G_rep_in));
  const arma::uword L_max = static_cast<arma::uword>(std::max(1, L_max_in));
  const arma::uword B_tot = static_cast<arma::uword>(B_tot_in);

  // Store with slice = P, col = P, row = B_tot so layout matches R array(dim = c(B_tot, P, P))!
  arma::cube Xty_bulk(B_tot, P, P, arma::fill::zeros);
  if (n_runs == 0 || P == 0 || B_tot == 0) return Xty_bulk;

  std::vector<arma::uvec> grp_cols(G_rep);
  for (arma::uword k = 0; k < G_rep; ++k) {
    std::vector<arma::uword> tmp;
    for (arma::uword g = 0; g < P; ++g) {
      if (static_cast<arma::uword>(std::max(1, rep_grp[g]) - 1) == k) tmp.push_back(g);
    }
    grp_cols[k] = arma::uvec(tmp);
  }

  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);

  for (arma::uword r = 0; r < n_runs; ++r) {
    Rcpp::NumericMatrix E_nm = Rcpp::as<Rcpp::NumericMatrix>(E_runs[r]);
    const arma::mat E_r(E_nm.begin(), static_cast<arma::uword>(E_nm.nrow()),
                        static_cast<arma::uword>(E_nm.ncol()), false, true);
    const arma::uword m = E_r.n_rows;
    if (m == 0) continue;

    const arma::uvec adj_idx0 = to_zero_based_uvec(Rcpp::as<Rcpp::IntegerVector>(adj_idx1_runs[r]));
    const std::vector<arma::umat> lag_pairs0 = to_zero_based_lag_pairs(Rcpp::as<Rcpp::List>(lag_pairs_runs[r]), L_max);
    const Rcpp::List X_list_r = Rcpp::as<Rcpp::List>(X_rep_runs[r]);

    // Precompute normalized lag-d residual cross-product matrices C_d (P x P) once per run!
    std::vector<arma::mat> C_lag(L_max + 1);
    if (L_max >= 2) {
      for (arma::uword d = 2; d <= L_max; ++d) {
        const arma::umat& pd = lag_pairs0[d - 1];
        const double N_d = 0.5 * static_cast<double>(pd.n_rows);
        if (N_d > 0.0) {
          C_lag[d] = (1.0 / std::sqrt(2.0 * N_d)) * (E_r.rows(pd.col(0)).t() * E_r.rows(pd.col(1)));
        }
      }
    }

    for (arma::uword ka = 0; ka < G_rep; ++ka) {
      const arma::uvec& ca = grp_cols[ka];
      if (ca.n_elem == 0) continue;
      const arma::mat EA = E_r.cols(ca);

      for (arma::uword kb = ka; kb < G_rep; ++kb) {
        const arma::uvec& cb = grp_cols[kb];
        if (cb.n_elem == 0) continue;
        const arma::mat EB = E_r.cols(cb);
        const arma::uword slot_idx = ka * G_rep + kb;
        Rcpp::NumericMatrix X_nm = Rcpp::as<Rcpp::NumericMatrix>(X_list_r[slot_idx]);
        const arma::mat Xr(X_nm.begin(), static_cast<arma::uword>(X_nm.nrow()),
                           static_cast<arma::uword>(X_nm.ncol()), false, true);

        arma::mat Z_ell(m, cb.n_elem);
        for (arma::uword ell = 0; ell < B_tot; ++ell) {
          const double* x_col = Xr.colptr(ell);
          Z_ell = EB;
          for (arma::uword c = 0; c < cb.n_elem; ++c) {
            double* zc = Z_ell.colptr(c);
            for (arma::uword t = 0; t < m; ++t) zc[t] *= x_col[t];
          }
          for (arma::uword k = 0; k < adj_idx0.n_elem; ++k) {
            const arma::uword i0 = adj_idx0(k);
            if (i0 + 1 < m) {
              const double w1 = x_col[m + i0] * inv_sqrt2;
              if (w1 != 0.0) {
                for (arma::uword c = 0; c < cb.n_elem; ++c) {
                  Z_ell(i0 + 1, c) += w1 * EB(i0, c);
                  Z_ell(i0, c) += w1 * EB(i0 + 1, c);
                }
              }
            }
          }

          arma::mat M_ell = EA.t() * Z_ell;
          if (L_max >= 2) {
            for (arma::uword d = 2; d <= L_max; ++d) {
              const double wd = x_col[2 * m - 1 + (d - 2)];
              if (std::abs(wd) > 1e-14 && C_lag[d].n_rows == P) {
                M_ell += wd * C_lag[d].submat(ca, cb);
              }
            }
          }

          for (arma::uword jb = 0; jb < cb.n_elem; ++jb) {
            const arma::uword h = cb(jb);
            for (arma::uword ia = 0; ia < ca.n_elem; ++ia) {
              const arma::uword g = ca(ia);
              const double val = M_ell(ia, jb);
              Xty_bulk(ell, g, h) += val;
              if (ka != kb) Xty_bulk(ell, h, g) += val;
            }
          }
        }
      }
    }
  }

  for (arma::uword g = 0; g < P; ++g) {
    if (g < Xty_diag.n_cols) {
      for (arma::uword ell = 0; ell < B_tot; ++ell) {
        Xty_bulk(ell, g, g) = Xty_diag(ell, g);
      }
    }
  }

  return Xty_bulk;
}

// [[Rcpp::export]]
arma::mat poet_extract_psd_factors_cpp(
    const arma::mat& C_in,
    int k_rank_in,
    int nC_in = 0,
    int n_lags_in = 1,
    Rcpp::Nullable<Rcpp::NumericMatrix> probe_cache_r = R_NilValue
) {
  arma::mat C_work = C_in;
  arma::mat probe_mat;
  const arma::mat* probe_ptr = nullptr;
  if (probe_cache_r.isNotNull()) {
    Rcpp::NumericMatrix pm(probe_cache_r);
    probe_mat = arma::mat(pm.begin(), static_cast<arma::uword>(pm.nrow()),
                          static_cast<arma::uword>(pm.ncol()), false, true);
    probe_ptr = &probe_mat;
  }
  return extract_psd_factors_inplace(
    C_work,
    static_cast<arma::uword>(std::max(0, k_rank_in)),
    static_cast<arma::uword>(std::max(0, nC_in)),
    static_cast<arma::uword>(std::max(1, n_lags_in)),
    probe_ptr
  );
}

// [[Rcpp::export]]
Rcpp::List poet_stage_e_stratum_moments_cpp(
    const Rcpp::List& U_si,
    const Rcpp::List& V_si,
    const Rcpp::List& ok_si,
    const Rcpp::List& lag_pairs,
    const arma::mat& f_dyn_r,
    const arma::vec& gm_vec,
    const arma::mat& Gamma_stat_cols,
    const arma::mat& B_dyn_cols,
    const arma::vec& sig0_rob_cols,
    const arma::vec& rsc2_cols,
    int L_max_in,
    int K_tcor_in,
    int D_lag_in,
    bool return_panel
) {
  const arma::uword n_lags = static_cast<arma::uword>(U_si.size());
  const arma::uword n_u = Gamma_stat_cols.n_cols;
  const arma::uword L_max = static_cast<arma::uword>(std::max(1, L_max_in));
  const arma::uword K_tcor = static_cast<arma::uword>(std::max(1, K_tcor_in));
  const arma::uword D_lag = static_cast<arma::uword>(std::max(1, D_lag_in));

  arma::uword m_t = 0;
  if (n_lags > 0) {
    Rcpp::NumericMatrix u0 = Rcpp::as<Rcpp::NumericMatrix>(U_si[0]);
    m_t = static_cast<arma::uword>(u0.ncol());
  }

  Rcpp::NumericVector s2_r(static_cast<R_xlen_t>(m_t * n_u * n_lags), R_NaReal);
  s2_r.attr("dim") = Rcpp::IntegerVector::create(m_t, n_u, n_lags);
  Rcpp::NumericVector se_r(static_cast<R_xlen_t>(m_t * n_u * n_lags), R_NaReal);
  se_r.attr("dim") = Rcpp::IntegerVector::create(m_t, n_u, n_lags);
  double* s2_ptr = REAL(s2_r);
  double* se_ptr = REAL(se_r);

  Rcpp::NumericVector tcor_r;
  double* tcor_ptr = nullptr;
  if (return_panel) {
    tcor_r = Rcpp::NumericVector(static_cast<R_xlen_t>(K_tcor * m_t * n_u * n_lags), 0.0);
    tcor_r.attr("dim") = Rcpp::IntegerVector::create(K_tcor, m_t, n_u, n_lags);
    tcor_ptr = REAL(tcor_r);
  }

  Rcpp::NumericVector lvec_r;
  double* lvec_ptr = nullptr;
  if (return_panel && n_lags >= 2) {
    lvec_r = Rcpp::NumericVector(static_cast<R_xlen_t>(D_lag * m_t * n_u * n_lags), 0.0);
    lvec_r.attr("dim") = Rcpp::IntegerVector::create(D_lag, m_t, n_u, n_lags);
    lvec_ptr = REAL(lvec_r);
  }

  if (n_lags == 0 || n_u == 0 || m_t == 0) {
    return Rcpp::List::create(
      Rcpp::Named("s2_trial") = s2_r,
      Rcpp::Named("se") = se_r,
      Rcpp::Named("tcor") = return_panel ? SEXP(tcor_r) : R_NilValue,
      Rcpp::Named("lvec") = (return_panel && n_lags >= 2) ? SEXP(lvec_r) : R_NilValue
    );
  }

  Rcpp::NumericMatrix u0_nm = Rcpp::as<Rcpp::NumericMatrix>(U_si[0]);
  const arma::uword m_rows = static_cast<arma::uword>(u0_nm.nrow());
  const std::vector<arma::umat> lag_pairs0 = to_zero_based_lag_pairs(lag_pairs, L_max);
  const std::vector<std::vector<LagSpan>> lag_spans = build_lag_spans(lag_pairs0, m_rows);
  const bool need_mlag = (Gamma_stat_cols.n_rows > 1) &&
    arma::any(arma::abs(arma::vectorise(Gamma_stat_cols.rows(1, Gamma_stat_cols.n_rows - 1))) > 1e-14);

  const arma::uword M_dyn = f_dyn_r.n_cols;
  arma::mat f_u(m_rows, M_dyn, arma::fill::zeros);
  arma::mat f_v(m_rows, M_dyn, arma::fill::zeros);
  bool any_fu = false;
  bool any_fv = false;
  const bool need_dyn = (B_dyn_cols.n_elem > 0) &&
    arma::any(arma::abs(arma::vectorise(B_dyn_cols)) > 1e-14);
  if (need_dyn && m_rows > 0) {
    for (arma::uword j = 0; j < M_dyn && j < gm_vec.n_elem; ++j) {
      const double gm_j = gm_vec(j);
      const double wu = 1.0 - gm_j;
      if (std::abs(wu) > 1e-14) {
        f_u.col(j) = wu * f_dyn_r.col(j);
        any_fu = true;
      }
      if (std::abs(gm_j) > 1e-14) {
        f_v.col(j) = gm_j * f_dyn_r.col(j);
        any_fv = true;
      }
    }
  }

  std::vector<arma::mat> U_cent_l(n_lags);
  std::vector<arma::mat> V_cent_l(n_lags);
  arma::mat W_prod, T_stat, T_dyn;

  for (arma::uword li = 0; li < n_lags; ++li) {
    Rcpp::NumericMatrix u_nm = Rcpp::as<Rcpp::NumericMatrix>(U_si[li]);
    Rcpp::NumericMatrix v_nm = Rcpp::as<Rcpp::NumericMatrix>(V_si[li]);
    const arma::mat Uu(u_nm.begin(), static_cast<arma::uword>(u_nm.nrow()),
                       static_cast<arma::uword>(u_nm.ncol()), false, true);
    const arma::mat Vu(v_nm.begin(), static_cast<arma::uword>(v_nm.nrow()),
                       static_cast<arma::uword>(v_nm.ncol()), false, true);

    std::vector<arma::uword> ok_vec;
    if (!Rf_isNull(ok_si[li])) {
      Rcpp::LogicalVector okl = Rcpp::as<Rcpp::LogicalVector>(ok_si[li]);
      for (int c = 0; c < okl.size(); ++c) {
        if (okl[c] == TRUE) ok_vec.push_back(static_cast<arma::uword>(c));
      }
    }
    const arma::uvec ok_idx(ok_vec);
    U_cent_l[li] = center_trial_op_cols(Uu, ok_idx);
    V_cent_l[li] = center_trial_op_cols(Vu, ok_idx);

    if (ok_idx.n_elem == 0) continue;

    const arma::mat Uc = U_cent_l[li].cols(ok_idx);
    const arma::mat Vc = V_cent_l[li].cols(ok_idx);
    eval_sem_basis_profiles(Uc, Uc, Vc, Vc, lag_spans, f_u, f_v, any_fu, any_fv, need_mlag, L_max, W_prod, T_stat, T_dyn);
    const arma::mat raw_m = T_stat * Gamma_stat_cols + T_dyn * B_dyn_cols;

    for (arma::uword u = 0; u < n_u; ++u) {
      const double rsc2 = rsc2_cols(u);
      const double s2_rob = sig0_rob_cols(u) * sig0_rob_cols(u) * rsc2;
      for (arma::uword k = 0; k < ok_idx.n_elem; ++k) {
        const arma::uword ti = ok_idx(k);
        const double v_base = std::max(1e-10, s2_rob * T_stat(k, 0));
        const double v_dyn = std::max(0.15 * v_base, rsc2 * raw_m(k, u));
        const R_xlen_t idx_s = static_cast<R_xlen_t>(ti) + m_t * (u + n_u * li);
        s2_ptr[idx_s] = v_dyn;
        se_ptr[idx_s] = std::sqrt(v_dyn);
      }
    }

    if (return_panel && m_t >= 2) {
      eval_sem_basis_profiles(Uu, Uu, Vu, Vu, lag_spans, f_u, f_v, any_fu, any_fv, need_mlag, L_max, W_prod, T_stat, T_dyn);
      const arma::mat raw_u = T_stat * Gamma_stat_cols + T_dyn * B_dyn_cols;
      arma::mat vu_dyn(m_t, n_u);
      for (arma::uword u = 0; u < n_u; ++u) {
        const double rsc2 = rsc2_cols(u);
        const double s2_rob = sig0_rob_cols(u) * sig0_rob_cols(u) * rsc2;
        for (arma::uword ti = 0; ti < m_t; ++ti) {
          const double vub = std::max(1e-12, s2_rob * T_stat(ti, 0));
          vu_dyn(ti, u) = std::max(0.15 * vub, rsc2 * raw_u(ti, u));
        }
      }

      std::vector<bool> ok_mask(m_t, false);
      for (arma::uword k = 0; k < ok_idx.n_elem; ++k) ok_mask[ok_idx(k)] = true;

      const arma::uword kd_max = std::min(K_tcor, m_t - 1);
      for (arma::uword kd = 1; kd <= kd_max; ++kd) {
        std::vector<arma::uword> ia_v, ib_v;
        for (arma::uword i = 0; i + kd < m_t; ++i) {
          if (ok_mask[i] && ok_mask[i + kd]) {
            ia_v.push_back(i);
            ib_v.push_back(i + kd);
          }
        }
        if (ia_v.empty()) continue;
        const arma::uvec ia_u(ia_v), ib_u(ib_v);
        eval_sem_basis_profiles(
          Uu.cols(ia_u), Uu.cols(ib_u), Vu.cols(ia_u), Vu.cols(ib_u),
          lag_spans, f_u, f_v, any_fu, any_fv, need_mlag, L_max, W_prod, T_stat, T_dyn
        );
        const arma::mat c_mat = T_stat * Gamma_stat_cols + T_dyn * B_dyn_cols;
        for (arma::uword u = 0; u < n_u; ++u) {
          const double rsc2 = rsc2_cols(u);
          for (arma::uword k = 0; k < ia_u.n_elem; ++k) {
            const arma::uword ia = ia_u(k);
            const arma::uword ib = ib_u(k);
            const double c_kd = rsc2 * c_mat(k, u);
            const double den = std::sqrt(std::max(1e-20, vu_dyn(ia, u) * vu_dyn(ib, u)));
            const double tc = std::min(0.95, std::max(-0.95, c_kd / den));
            const R_xlen_t idx_tc = static_cast<R_xlen_t>(kd - 1) + K_tcor * (ia + m_t * (u + n_u * li));
            tcor_ptr[idx_tc] = tc;
          }
        }
      }
    }
  }

  if (return_panel && n_lags >= 2) {
    // Cross-lag covariances for each pair (la, lb), then per-(u, ti) PSD projection + eigendecomposition
    std::vector<arma::cube> Cg_u(n_u, arma::cube(n_lags, n_lags, m_t, arma::fill::zeros));
    for (arma::uword la = 0; la < n_lags; ++la) {
      for (arma::uword u = 0; u < n_u; ++u) {
        for (arma::uword ti = 0; ti < m_t; ++ti) {
          const double v = s2_ptr[ti + m_t * (u + n_u * la)];
          Cg_u[u](la, la, ti) = std::isfinite(v) ? v : 0.0;
        }
      }
      for (arma::uword lb = la + 1; lb < n_lags; ++lb) {
        eval_sem_basis_profiles(
          U_cent_l[la], U_cent_l[lb], V_cent_l[la], V_cent_l[lb],
          lag_spans, f_u, f_v, any_fu, any_fv, need_mlag, L_max, W_prod, T_stat, T_dyn
        );
        const arma::mat cab_m = T_stat * Gamma_stat_cols + T_dyn * B_dyn_cols;
        for (arma::uword u = 0; u < n_u; ++u) {
          const double rsc2 = rsc2_cols(u);
          for (arma::uword ti = 0; ti < m_t; ++ti) {
            const double cab = rsc2 * cab_m(ti, u);
            Cg_u[u](la, lb, ti) = cab;
            Cg_u[u](lb, la, ti) = cab;
          }
        }
      }
    }

    arma::vec eg_val;
    arma::mat eg_vec;
    for (arma::uword u = 0; u < n_u; ++u) {
      for (arma::uword ti = 0; ti < m_t; ++ti) {
        const arma::mat Ci = project_trial_err_psd_impl(Cg_u[u].slice(ti), 0.98);
        arma::eig_sym(eg_val, eg_vec, Ci);
        arma::uword k_l = 0;
        for (arma::uword j = 0; j < eg_val.n_elem; ++j) {
          if (eg_val(eg_val.n_elem - 1 - j) > 1e-12 && k_l < D_lag) ++k_l;
        }
        for (arma::uword j = 0; j < k_l; ++j) {
          const arma::uword src_col = eg_val.n_elem - 1 - j;
          const double sc_ev = std::sqrt(eg_val(src_col));
          const double* v_ptr = eg_vec.colptr(src_col);
          for (arma::uword la = 0; la < n_lags; ++la) {
            const R_xlen_t idx_lv = static_cast<R_xlen_t>(j) + D_lag * (ti + m_t * (u + n_u * la));
            lvec_ptr[idx_lv] = v_ptr[la] * sc_ev;
          }
        }
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("s2_trial") = s2_r,
    Rcpp::Named("se") = se_r,
    Rcpp::Named("tcor") = return_panel ? SEXP(tcor_r) : R_NilValue,
    Rcpp::Named("lvec") = (return_panel && n_lags >= 2) ? SEXP(lvec_r) : R_NilValue
  );
}

// [[Rcpp::export]]
Rcpp::NumericVector poet_stage_e_component_svec_cpp(
    const Rcpp::List& U_sub_strata,
    const Rcpp::List& ok_sub_strata,
    const Rcpp::List& vsp_sub_strata,
    const Rcpp::IntegerVector& map_s_idx1,
    const Rcpp::List& G_stat_sub,
    const Rcpp::List& B_dyn_sub,
    const Rcpp::List& lag_pairs,
    const arma::mat& f_dyn_r,
    const arma::vec& gm_vec,
    const arma::vec& run_scale_Cu,
    const Rcpp::NumericVector& s2_comp_arr,
    Rcpp::Nullable<Rcpp::NumericMatrix> probe_cache_r,
    int d_m_in,
    int m_t_in,
    int n_lags_in,
    int L_max_in
) {
  const arma::uword n_sC = static_cast<arma::uword>(U_sub_strata.size());
  const arma::uword nC = static_cast<arma::uword>(map_s_idx1.size());
  const arma::uword d_m = static_cast<arma::uword>(std::max(1, d_m_in));
  const arma::uword m_t = static_cast<arma::uword>(std::max(1, m_t_in));
  const arma::uword n_lags = static_cast<arma::uword>(std::max(1, n_lags_in));
  const arma::uword L_max = static_cast<arma::uword>(std::max(1, L_max_in));
  const arma::uword q_m = nC * n_lags;
  const arma::uword q_sC = n_sC * n_lags;

  Rcpp::NumericVector svec_out(static_cast<R_xlen_t>(d_m * m_t * nC * n_lags), 0.0);
  svec_out.attr("dim") = Rcpp::IntegerVector::create(d_m, m_t, nC, n_lags);
  if (nC == 0 || n_sC == 0) return svec_out;

  Rcpp::List u0_list = Rcpp::as<Rcpp::List>(U_sub_strata[0]);
  Rcpp::NumericMatrix u00 = Rcpp::as<Rcpp::NumericMatrix>(u0_list[0]);
  const arma::uword m = static_cast<arma::uword>(u00.nrow());
  if (m == 0) return svec_out;

  struct VspInfo {
    bool has_rwls = false;
    int mode = 0; // 0 = none, 1 = scale, 2 = woodbury
    arma::vec sqrt_w;
    arma::mat ua;
    arma::vec inv_fac;
    bool has_udiff = false;
    arma::vec u_diff;
    std::vector<arma::vec> b_row_sums;
    CsrMat W_obs;
  };

  std::vector<VspInfo> vsp(n_sC);
  for (arma::uword isi = 0; isi < n_sC; ++isi) {
    Rcpp::List vl = Rcpp::as<Rcpp::List>(vsp_sub_strata[isi]);
    VspInfo& info = vsp[isi];
    Rcpp::List w_csr = Rcpp::as<Rcpp::List>(vl["W_obs_csr"]);
    info.W_obs = parse_csr_components(
      Rcpp::as<Rcpp::IntegerVector>(w_csr["row_ptr"]),
      Rcpp::as<Rcpp::IntegerVector>(w_csr["col_idx"]),
      Rcpp::as<Rcpp::NumericVector>(w_csr["val"]),
      m, m
    );
    info.has_rwls = vl.containsElementNamed("has_rwls") && Rcpp::as<bool>(vl["has_rwls"]);
    if (info.has_rwls) {
      const std::string rmode = Rcpp::as<std::string>(vl["rwls_mode"]);
      if (rmode == "scale" && !Rf_isNull(vl["sqrt_w_use"])) {
        info.mode = 1;
        info.sqrt_w = Rcpp::as<arma::vec>(vl["sqrt_w_use"]);
      } else if (rmode == "woodbury" && !Rf_isNull(vl["ua"])) {
        info.mode = 2;
        info.ua = Rcpp::as<arma::mat>(vl["ua"]);
        info.inv_fac = Rcpp::as<arma::vec>(vl["inv_fac_a"]);
      }
      if (!Rf_isNull(vl["u_diff"]) && !Rf_isNull(vl["b_row_sums"])) {
        info.has_udiff = true;
        info.u_diff = Rcpp::as<arma::vec>(vl["u_diff"]);
        Rcpp::List brs = Rcpp::as<Rcpp::List>(vl["b_row_sums"]);
        info.b_row_sums.resize(n_lags);
        for (arma::uword la = 0; la < n_lags; ++la) {
          if (static_cast<R_xlen_t>(la) < brs.size() && !Rf_isNull(brs[la])) {
            info.b_row_sums[la] = Rcpp::as<arma::vec>(brs[la]);
          }
        }
      }
    }
  }

  // Bind zero-copy views to U_sub_strata and precompute per-(la, isi) centering means & scales
  std::vector<arma::mat> U_views(q_sC);
  std::vector<std::vector<bool>> ok_masks(q_sC, std::vector<bool>(m_t, false));
  arma::uvec m_ok_vec(q_sC, arma::fill::zeros);
  arma::vec inv_sc_vec(q_sC, arma::fill::ones);
  arma::mat mu_U(m, q_sC, arma::fill::zeros);
  arma::vec mu_s(q_sC, arma::fill::zeros);

  for (arma::uword isi = 0; isi < n_sC; ++isi) {
    Rcpp::List u_l = Rcpp::as<Rcpp::List>(U_sub_strata[isi]);
    Rcpp::List ok_l = Rcpp::as<Rcpp::List>(ok_sub_strata[isi]);
    for (arma::uword la = 0; la < n_lags; ++la) {
      const arma::uword col_k = la * n_sC + isi;
      Rcpp::NumericMatrix unm = Rcpp::as<Rcpp::NumericMatrix>(u_l[la]);
      U_views[col_k] = arma::mat(unm.begin(), m, m_t, false, true);
      arma::uword cnt = 0;
      if (!Rf_isNull(ok_l[la])) {
        Rcpp::LogicalVector okl = Rcpp::as<Rcpp::LogicalVector>(ok_l[la]);
        for (arma::uword ti = 0; ti < m_t && static_cast<R_xlen_t>(ti) < okl.size(); ++ti) {
          if (okl[ti] == TRUE) {
            ok_masks[col_k][ti] = true;
            ++cnt;
            mu_U.col(col_k) += U_views[col_k].col(ti);
            if (vsp[isi].has_udiff && vsp[isi].b_row_sums[la].n_elem > ti) {
              mu_s(col_k) += vsp[isi].b_row_sums[la](ti);
            }
          }
        }
      }
      m_ok_vec(col_k) = cnt;
      if (cnt > 1) {
        mu_U.col(col_k) /= static_cast<double>(cnt);
        mu_s(col_k) /= static_cast<double>(cnt);
        inv_sc_vec(col_k) = 1.0 / std::sqrt(1.0 - 1.0 / static_cast<double>(cnt));
      }
    }
  }

  arma::mat probe_mat;
  const arma::mat* probe_ptr = nullptr;
  if (probe_cache_r.isNotNull()) {
    Rcpp::NumericMatrix pm(probe_cache_r);
    probe_mat = arma::mat(pm.begin(), static_cast<arma::uword>(pm.nrow()),
                          static_cast<arma::uword>(pm.ncol()), false, true);
    probe_ptr = &probe_mat;
  }

  const arma::uvec map_s0 = to_zero_based_uvec(map_s_idx1);
  bool identity_map = (n_sC == nC);
  if (identity_map) {
    for (arma::uword i = 0; i < nC; ++i) {
      if (map_s0(i) != i) {
        identity_map = false;
        break;
      }
    }
  }
  const std::vector<arma::umat> lag_pairs0 = to_zero_based_lag_pairs(lag_pairs, L_max);
  const std::vector<std::vector<LagSpan>> lag_spans = build_lag_spans(lag_pairs0, m);
  const double* s2_in_ptr = REAL(s2_comp_arr);
  double* svec_ptr = REAL(svec_out);

  // Pre-parse and pre-scale (nC x nC) coefficient matrices once outside the trial loop!
  const bool has_rsc = (run_scale_Cu.n_elem == nC);
  auto scale_coeff_mat = [&](const arma::mat& M_in) {
    if (!has_rsc) return M_in;
    arma::mat M_sc = M_in;
    M_sc.each_col() %= run_scale_Cu;
    M_sc.each_row() %= run_scale_Cu.t();
    return M_sc;
  };
  const arma::uword n_gstat = static_cast<arma::uword>(G_stat_sub.size());
  std::vector<arma::mat> G_mats(n_gstat);
  std::vector<bool> G_active(n_gstat, false);
  for (arma::uword d = 0; d < n_gstat; ++d) {
    Rcpp::NumericMatrix g_nm = Rcpp::as<Rcpp::NumericMatrix>(G_stat_sub[d]);
    const arma::mat Gm(g_nm.begin(), nC, nC, false, true);
    if (d == 0 || arma::any(arma::abs(arma::vectorise(Gm)) > 1e-14)) {
      G_mats[d] = scale_coeff_mat(Gm);
      G_active[d] = true;
    }
  }
  const arma::uword M_dyn = static_cast<arma::uword>(B_dyn_sub.size());
  std::vector<arma::mat> B_mats(M_dyn);
  std::vector<bool> B_active(M_dyn, false);
  for (arma::uword j = 0; j < M_dyn; ++j) {
    if (j >= f_dyn_r.n_cols || j >= gm_vec.n_elem) break;
    Rcpp::NumericMatrix b_nm = Rcpp::as<Rcpp::NumericMatrix>(B_dyn_sub[j]);
    const arma::mat Bm(b_nm.begin(), nC, nC, false, true);
    if (arma::any(arma::abs(arma::vectorise(Bm)) > 1e-14)) {
      B_mats[j] = scale_coeff_mat(Bm);
      B_active[j] = true;
    }
  }

  // Single reusable trial buffers across all m_t trials (zero per-trial heap allocations!)
  arma::mat U_ti(m, q_sC, arma::fill::zeros);
  arma::mat V_ti(m, q_sC, arma::fill::zeros);
  arma::mat W_sc(m, q_sC, arma::fill::zeros);
  arma::mat U_p0, U_p1;
  arma::vec u0_buf(m);
  arma::mat Ci(q_m, q_m, arma::fill::zeros);
  arma::mat T_buf(q_sC, q_sC, arma::fill::zeros);

  auto add_tiled_term = [&](const arma::mat& M_coeff, const arma::mat& T_s) {
    if (identity_map) {
      for (arma::uword lb = 0; lb < n_lags; ++lb) {
        for (arma::uword ib = 0; ib < nC; ++ib) {
          const arma::uword col_ci = lb * nC + ib;
          const double* mc_col = M_coeff.colptr(ib);
          for (arma::uword la = 0; la < n_lags; ++la) {
            double* ci_blk = Ci.colptr(col_ci) + la * nC;
            const double* ts_blk = T_s.colptr(col_ci) + la * nC;
            for (arma::uword ia = 0; ia < nC; ++ia) {
              ci_blk[ia] += mc_col[ia] * ts_blk[ia];
            }
          }
        }
      }
    } else {
      for (arma::uword lb = 0; lb < n_lags; ++lb) {
        for (arma::uword ib = 0; ib < nC; ++ib) {
          const arma::uword col_ci = lb * nC + ib;
          const arma::uword col_ts = lb * n_sC + map_s0(ib);
          const double* mc_col = M_coeff.colptr(ib);
          for (arma::uword la = 0; la < n_lags; ++la) {
            double* ci_blk = Ci.colptr(col_ci) + la * nC;
            const double* ts_blk = T_s.colptr(col_ts) + la * n_sC;
            for (arma::uword ia = 0; ia < nC; ++ia) {
              ci_blk[ia] += mc_col[ia] * ts_blk[map_s0(ia)];
            }
          }
        }
      }
    }
  };

  auto scale_cols_by_vec = [&](const arma::mat& X_in, const double* fv) {
    for (arma::uword c = 0; c < q_sC; ++c) {
      const double* xp = X_in.colptr(c);
      double* wp = W_sc.colptr(c);
      for (arma::uword r = 0; r < m; ++r) wp[r] = xp[r] * fv[r];
    }
  };

  for (arma::uword ti = 0; ti < m_t; ++ti) {
    for (arma::uword la = 0; la < n_lags; ++la) {
      for (arma::uword isi = 0; isi < n_sC; ++isi) {
        const arma::uword col_k = la * n_sC + isi;
        double* u_dst = U_ti.colptr(col_k);
        double* v_dst = V_ti.colptr(col_k);
        if (!ok_masks[col_k][ti] || m_ok_vec(col_k) == 0) {
          std::memset(u_dst, 0, static_cast<std::size_t>(m) * sizeof(double));
          std::memset(v_dst, 0, static_cast<std::size_t>(m) * sizeof(double));
          continue;
        }
        const double* u_src = U_views[col_k].colptr(ti);
        if (m_ok_vec(col_k) == 1) {
          std::memcpy(u_dst, u_src, static_cast<std::size_t>(m) * sizeof(double));
        } else {
          const double* mu_p = mu_U.colptr(col_k);
          const double sc = inv_sc_vec(col_k);
          for (arma::uword r = 0; r < m; ++r) u_dst[r] = (u_src[r] - mu_p[r]) * sc;
        }

        const VspInfo& info = vsp[isi];
        if (!info.has_rwls) {
          csr_trans_mult_vec(info.W_obs, u_dst, v_dst);
        } else {
          double* u0p = u0_buf.memptr();
          if (info.mode == 1) {
            const double* sw = info.sqrt_w.memptr();
            for (arma::uword r = 0; r < m; ++r) u0p[r] = u_dst[r] * sw[r];
          } else if (info.mode == 2 && info.ua.n_cols > 0) {
            std::memcpy(u0p, u_dst, static_cast<std::size_t>(m) * sizeof(double));
            for (arma::uword k = 0; k < info.ua.n_cols; ++k) {
              const double* uak = info.ua.colptr(k);
              double dot_v = 0.0;
              for (arma::uword r = 0; r < m; ++r) dot_v += uak[r] * u_dst[r];
              const double mult = dot_v * info.inv_fac(k);
              for (arma::uword r = 0; r < m; ++r) u0p[r] += mult * uak[r];
            }
          } else {
            std::memcpy(u0p, u_dst, static_cast<std::size_t>(m) * sizeof(double));
          }
          if (info.has_udiff && info.b_row_sums[la].n_elem > ti) {
            const double s_raw = info.b_row_sums[la](ti);
            const double s_c = (m_ok_vec(col_k) > 1) ? ((s_raw - mu_s(col_k)) * inv_sc_vec(col_k)) : s_raw;
            const double* udp = info.u_diff.memptr();
            for (arma::uword r = 0; r < m; ++r) u0p[r] -= s_c * udp[r];
          }
          csr_trans_mult_vec(info.W_obs, u0p, v_dst);
        }
      }
    }

    Ci.zeros();
    T_buf = U_ti.t() * U_ti;
    if (n_gstat >= 1 && G_active[0]) {
      add_tiled_term(G_mats[0], T_buf);
    }
    for (arma::uword d = 1; d <= L_max && d < n_gstat; ++d) {
      if (!G_active[d] || (d - 1) >= lag_spans.size()) continue;
      const std::vector<LagSpan>& spans_d = lag_spans[d - 1];
      if (spans_d.empty()) continue;
      if (spans_d.size() <= 16) {
        T_buf.zeros(q_sC, q_sC);
        for (std::size_t s = 0; s < spans_d.size(); ++s) {
          const LagSpan& sp = spans_d[s];
          T_buf += U_ti.rows(sp.s0, sp.s0 + sp.len - 1).t() * U_ti.rows(sp.s1, sp.s1 + sp.len - 1);
        }
      } else {
        const arma::umat& pd = lag_pairs0[d - 1];
        const arma::uword n_p = pd.n_rows;
        U_p0.set_size(n_p, q_sC);
        U_p1.set_size(n_p, q_sC);
        const arma::uword* p0 = pd.colptr(0);
        const arma::uword* p1 = pd.colptr(1);
        for (arma::uword c = 0; c < q_sC; ++c) {
          const double* uc = U_ti.colptr(c);
          double* d0 = U_p0.colptr(c);
          double* d1 = U_p1.colptr(c);
          for (arma::uword r = 0; r < n_p; ++r) {
            d0[r] = (p0[r] < m) ? uc[p0[r]] : 0.0;
            d1[r] = (p1[r] < m) ? uc[p1[r]] : 0.0;
          }
        }
        T_buf = U_p0.t() * U_p1;
      }
      add_tiled_term(G_mats[d], T_buf);
    }
    for (arma::uword j = 0; j < M_dyn; ++j) {
      if (!B_active[j]) continue;
      const double* fj = f_dyn_r.colptr(j);
      const double gm_j = gm_vec(j);
      T_buf.zeros(q_sC, q_sC);
      if (gm_j < 1.0) {
        scale_cols_by_vec(U_ti, fj);
        T_buf += (1.0 - gm_j) * (U_ti.t() * W_sc);
      }
      if (gm_j > 0.0) {
        scale_cols_by_vec(V_ti, fj);
        T_buf += gm_j * (V_ti.t() * W_sc);
      }
      add_tiled_term(B_mats[j], T_buf);
    }
    for (arma::uword la = 0; la < n_lags; ++la) {
      for (arma::uword ia = 0; ia < nC; ++ia) {
        const double v = s2_in_ptr[ti + m_t * (ia + nC * la)];
        Ci(la * nC + ia, la * nC + ia) = std::isfinite(v) ? v : 0.0;
      }
    }

    const arma::mat F_ti = extract_psd_factors_inplace(Ci, d_m, nC, n_lags, probe_ptr);
    for (arma::uword la = 0; la < n_lags; ++la) {
      for (arma::uword ia = 0; ia < nC; ++ia) {
        const double* src_f = F_ti.colptr(la * nC + ia);
        double* dst_f = svec_ptr + static_cast<R_xlen_t>(d_m) * (ti + m_t * (ia + nC * la));
        std::memcpy(dst_f, src_f, static_cast<std::size_t>(d_m) * sizeof(double));
      }
    }
  }

  return svec_out;
}

// [[Rcpp::export]]
Rcpp::List poet_cols_to_list_cpp(const Rcpp::NumericMatrix& M) {
  const R_xlen_t nr = M.nrow();
  const R_xlen_t nc = M.ncol();
  Rcpp::List out(nc);
  const double* src = M.begin();
  const std::size_t nbytes = static_cast<std::size_t>(nr) * sizeof(double);
  for (R_xlen_t j = 0; j < nc; ++j) {
    Rcpp::NumericVector v = Rcpp::no_init(nr);
    if (nr > 0) std::memcpy(REAL(v), src + j * nr, nbytes);
    out[j] = v;
  }
  return out;
}
