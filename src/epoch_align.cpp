#include "medusa.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

enum SupportCode : unsigned char {
  SUP_OBSERVED = 0,
  SUP_INTERPOLATED = 1,
  SUP_CENSORED = 2,
  SUP_MISSING_VALUE = 3,
  SUP_MISSING_EPOCH = 4,
  SUP_GAP = 5,
  SUP_OUTSIDE = 6
};

enum LagBracketMode : unsigned char {
  MODE_MISSING_EPOCH = 0,
  MODE_OUTSIDE = 1,
  MODE_EXACT = 2,
  MODE_BRACKET = 3
};

struct LagSpec {
  LagBracketMode mode;
  bool is_gap;
  int j_lo;       // 0-based index within trial window [0 .. n_win - 1]
  int j_hi;       // 0-based index within trial window [0 .. n_win - 1]
  int r_lo;       // 0-based index within run clock [0 .. T - 1]
  int r_hi;       // 0-based index within run clock [0 .. T - 1]
  int src_left;   // 1-based source_id (r_lo + 1), or NA_INTEGER
  int src_right;  // 1-based source_id (r_hi + 1), or NA_INTEGER
  double weight;  // interpolation weight in [0, 1]
};

inline bool same_clock_entry(double a, double b) {
  const bool fa = std::isfinite(a);
  const bool fb = std::isfinite(b);
  if (fa != fb) return false;
  return !fa || (a == b);
}

inline bool is_cell_censored(const int* cens_ptr, R_xlen_t cens_len, bool cens_is_matrix,
                             int r, int v, int n_time) {
  if (cens_len == 0) return false;
  if (!cens_is_matrix) {
    const int val = cens_ptr[r];
    return (val == NA_LOGICAL) || (val != 0);
  }
  const R_xlen_t idx = static_cast<R_xlen_t>(r) + static_cast<R_xlen_t>(v) * static_cast<R_xlen_t>(n_time);
  const int val = cens_ptr[idx];
  return (val == NA_LOGICAL) || (val != 0);
}

inline void compute_mean_median_sd(std::vector<double>& vals, double& out_mean,
                                   double& out_median, double& out_sd) {
  const size_t n = vals.size();
  if (n == 0) {
    out_mean = NA_REAL;
    out_median = NA_REAL;
    out_sd = NA_REAL;
    return;
  }
  double sum = 0.0;
  for (size_t i = 0; i < n; ++i) {
    sum += vals[i];
  }
  const double mean_val = sum / static_cast<double>(n);
  out_mean = mean_val;

  if (n == 1) {
    out_median = vals[0];
    out_sd = NA_REAL;
    return;
  }

  double ss = 0.0;
  for (size_t i = 0; i < n; ++i) {
    const double d = vals[i] - mean_val;
    ss += d * d;
  }
  out_sd = std::sqrt(ss / static_cast<double>(n - 1));

  const size_t mid = n / 2;
  std::nth_element(vals.begin(), vals.begin() + mid, vals.end());
  if ((n & 1U) == 1U) {
    out_median = vals[mid];
  } else {
    const double hi = vals[mid];
    const double lo = *std::max_element(vals.begin(), vals.begin() + mid);
    out_median = 0.5 * (lo + hi);
  }
}

inline std::string join_sorted_ids(const std::vector<int>& ids) {
  if (ids.empty()) return std::string();
  std::string out = std::to_string(ids[0]);
  for (size_t i = 1; i < ids.size(); ++i) {
    out.push_back(';');
    out += std::to_string(ids[i]);
  }
  return out;
}

} // namespace

// [[Rcpp::export]]
bool check_shared_clock_cpp(NumericVector time_vec, int n_time, int n_units) {
  if (n_time < 0 || n_units < 1) return false;
  const R_xlen_t expected = static_cast<R_xlen_t>(n_time) * static_cast<R_xlen_t>(n_units);
  if (time_vec.size() != expected) return false;
  if (n_units == 1 || n_time == 0) return true;

  const double* ptr = REAL(time_vec);
  for (int u = 1; u < n_units; ++u) {
    const double* u_ptr = ptr + static_cast<R_xlen_t>(u) * static_cast<R_xlen_t>(n_time);
    for (int t = 0; t < n_time; ++t) {
      if (!same_clock_entry(ptr[t], u_ptr[t])) {
        return false;
      }
    }
  }
  return true;
}

// [[Rcpp::export]]
List fast_pivot_decon_cpp(NumericVector time_vec, NumericVector vnum_vec,
                          IntegerVector group_idx, NumericVector value_vec,
                          Nullable<LogicalVector> valid_vec) {
  const R_xlen_t n = time_vec.size();
  if (vnum_vec.size() != n || group_idx.size() != n || value_vec.size() != n) {
    stop("Mismatched column lengths when pivoting deconvolution table");
  }
  const bool has_valid = valid_vec.isNotNull();
  LogicalVector valid_in;
  if (has_valid) {
    valid_in = LogicalVector(valid_vec);
    if (valid_in.size() != n) {
      stop("valid column length does not match deconvolution table");
    }
  }

  const double* t_ptr = REAL(time_vec);
  const double* v_ptr = REAL(vnum_vec);
  const int* g_ptr = INTEGER(group_idx);
  const double* y_ptr = REAL(value_vec);
  const int* val_ptr = has_valid ? LOGICAL(valid_in) : nullptr;

  // Count finite timestamps first; handle all-NA time tables cleanly
  bool any_finite_time = false;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (std::isfinite(t_ptr[i])) {
      any_finite_time = true;
      break;
    }
  }
  if (!any_finite_time) {
    // Collect unique (vnum, group_idx) units
    std::vector<std::pair<double, int>> units;
    for (R_xlen_t i = 0; i < n; ++i) {
      units.emplace_back(v_ptr[i], g_ptr[i]);
    }
    std::sort(units.begin(), units.end());
    units.erase(std::unique(units.begin(), units.end()), units.end());
    const int n_u = std::max(1, static_cast<int>(units.size()));
    IntegerVector u_grp(n_u);
    for (int u = 0; u < static_cast<int>(units.size()); ++u) {
      u_grp[u] = units[u].second;
    }
    NumericMatrix empty_y(0, n_u);
    LogicalVector empty_cens(0);
    return List::create(
      Named("rectangular") = true,
      Named("y") = empty_y,
      Named("clock") = NumericVector(0),
      Named("censored") = empty_cens,
      Named("cens_is_matrix") = false,
      Named("unit_group") = u_grp
    );
  }

  // Fast-path A: Time-major V x T layout emitted by mat2df(t(deconv_tv))
  if (n > 0 && std::isfinite(t_ptr[0])) {
    R_xlen_t v_cand = 1;
    const double t0 = t_ptr[0];
    while (v_cand < n && same_clock_entry(t_ptr[v_cand], t0)) {
      ++v_cand;
    }
    if (v_cand >= 1 && (n % v_cand) == 0) {
      const int n_u = static_cast<int>(v_cand);
      const int n_t_raw = static_cast<int>(n / v_cand);
      // Verify unique (vnum, group) in first block and exact repetition across all T blocks
      std::vector<std::pair<double, int>> first_units(n_u);
      for (int u = 0; u < n_u; ++u) {
        first_units[u] = std::make_pair(v_ptr[u], g_ptr[u]);
      }
      std::vector<std::pair<double, int>> sorted_u = first_units;
      std::sort(sorted_u.begin(), sorted_u.end());
      bool units_unique = (std::adjacent_find(sorted_u.begin(), sorted_u.end()) == sorted_u.end());

      bool layout_ok = units_unique;
      std::vector<int> fin_t_rows;
      fin_t_rows.reserve(n_t_raw);
      for (int t = 0; t < n_t_raw && layout_ok; ++t) {
        const R_xlen_t base = static_cast<R_xlen_t>(t) * v_cand;
        const double tt = t_ptr[base];
        if (std::isfinite(tt)) fin_t_rows.push_back(t);
        for (int u = 0; u < n_u; ++u) {
          if (!same_clock_entry(t_ptr[base + u], tt) ||
              v_ptr[base + u] != v_ptr[u] ||
              g_ptr[base + u] != g_ptr[u]) {
            layout_ok = false;
            break;
          }
        }
      }

      if (layout_ok) {
        const int n_t = static_cast<int>(fin_t_rows.size());
        std::vector<std::pair<double, int>> clock_ord(n_t);
        for (int k = 0; k < n_t; ++k) {
          const int t_orig = fin_t_rows[k];
          clock_ord[k] = std::make_pair(t_ptr[static_cast<R_xlen_t>(t_orig) * v_cand], t_orig);
        }
        std::sort(clock_ord.begin(), clock_ord.end(),
                  [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                    return a.first < b.first;
                  });
        for (int k = 1; k < n_t; ++k) {
          if (clock_ord[k].first == clock_ord[k - 1].first) {
            stop("Duplicate source times within an epoch unit");
          }
        }

        NumericVector clock_out(n_t);
        for (int k = 0; k < n_t; ++k) clock_out[k] = clock_ord[k].first;

        NumericMatrix y_out = no_init_matrix(n_t, n_u);
        double* y_out_ptr = REAL(y_out);

        // Cache-blocked transpose from (V x T_raw) into (T_fin x V)
        constexpr int TILE = 32;
        for (int k0 = 0; k0 < n_t; k0 += TILE) {
          const int k_lim = std::min(n_t, k0 + TILE);
          for (int u0 = 0; u0 < n_u; u0 += TILE) {
            const int u_lim = std::min(n_u, u0 + TILE);
            for (int k = k0; k < k_lim; ++k) {
              const R_xlen_t src_base = static_cast<R_xlen_t>(clock_ord[k].second) * v_cand;
              for (int u = u0; u < u_lim; ++u) {
                y_out_ptr[static_cast<R_xlen_t>(k) + static_cast<R_xlen_t>(u) * n_t] = y_ptr[src_base + u];
              }
            }
          }
        }

        LogicalVector cens_out;
        bool cens_is_mat = false;
        if (has_valid) {
          // Check if valid is uniform within each timepoint
          bool uniform_t = true;
          for (int k = 0; k < n_t && uniform_t; ++k) {
            const R_xlen_t src_base = static_cast<R_xlen_t>(clock_ord[k].second) * v_cand;
            const bool c0 = (val_ptr[src_base] == NA_LOGICAL) || (val_ptr[src_base] == 0);
            for (int u = 1; u < n_u; ++u) {
              const bool cu = (val_ptr[src_base + u] == NA_LOGICAL) || (val_ptr[src_base + u] == 0);
              if (cu != c0) {
                uniform_t = false;
                break;
              }
            }
          }
          if (uniform_t) {
            cens_out = LogicalVector(n_t);
            for (int k = 0; k < n_t; ++k) {
              const R_xlen_t src_base = static_cast<R_xlen_t>(clock_ord[k].second) * v_cand;
              cens_out[k] = (val_ptr[src_base] == NA_LOGICAL) || (val_ptr[src_base] == 0);
            }
          } else {
            cens_is_mat = true;
            cens_out = LogicalVector(static_cast<R_xlen_t>(n_t) * n_u);
            int* co_ptr = LOGICAL(cens_out);
            for (int k = 0; k < n_t; ++k) {
              const R_xlen_t src_base = static_cast<R_xlen_t>(clock_ord[k].second) * v_cand;
              for (int u = 0; u < n_u; ++u) {
                co_ptr[static_cast<R_xlen_t>(k) + static_cast<R_xlen_t>(u) * n_t] =
                    (val_ptr[src_base + u] == NA_LOGICAL) || (val_ptr[src_base + u] == 0);
              }
            }
          }
        } else {
          cens_out = LogicalVector(0);
        }

        IntegerVector u_grp(n_u);
        for (int u = 0; u < n_u; ++u) u_grp[u] = g_ptr[u];

        return List::create(
          Named("rectangular") = true,
          Named("y") = y_out,
          Named("clock") = clock_out,
          Named("censored") = cens_out,
          Named("cens_is_matrix") = cens_is_mat,
          Named("unit_group") = u_grp
        );
      }
    }
  }

  // General rectangular grid check for arbitrary row order
  std::vector<double> clock_vals;
  std::vector<std::pair<double, int>> units;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (std::isfinite(t_ptr[i])) {
      clock_vals.push_back(t_ptr[i]);
    }
    units.emplace_back(v_ptr[i], g_ptr[i]);
  }
  std::sort(clock_vals.begin(), clock_vals.end());
  clock_vals.erase(std::unique(clock_vals.begin(), clock_vals.end()), clock_vals.end());

  std::sort(units.begin(), units.end());
  units.erase(std::unique(units.begin(), units.end()), units.end());

  const int n_t = static_cast<int>(clock_vals.size());
  const int n_u = static_cast<int>(units.size());
  const R_xlen_t total_cells = static_cast<R_xlen_t>(n_t) * static_cast<R_xlen_t>(n_u);
  if (n_u == 0 || total_cells > n) {
    // Missing unit-by-time combinations: ragged per-unit clock
    return List::create(Named("rectangular") = false);
  }

  NumericMatrix y_out(n_t, n_u);
  double* y_out_ptr = REAL(y_out);
  std::vector<unsigned char> seen(total_cells, 0);
  LogicalVector cens_mat;
  int* cens_mat_ptr = nullptr;
  if (has_valid) {
    cens_mat = LogicalVector(total_cells);
    cens_mat_ptr = LOGICAL(cens_mat);
  }

  R_xlen_t finite_seen = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    const double tv = t_ptr[i];
    if (!std::isfinite(tv)) continue;
    const int ti = static_cast<int>(std::lower_bound(clock_vals.begin(), clock_vals.end(), tv) - clock_vals.begin());
    const std::pair<double, int> key(v_ptr[i], g_ptr[i]);
    const int ui = static_cast<int>(std::lower_bound(units.begin(), units.end(), key) - units.begin());
    const R_xlen_t cell = static_cast<R_xlen_t>(ti) + static_cast<R_xlen_t>(ui) * static_cast<R_xlen_t>(n_t);
    if (seen[cell] != 0) {
      stop("Duplicate source times within an epoch unit");
    }
    seen[cell] = 1;
    ++finite_seen;
    y_out_ptr[cell] = y_ptr[i];
    if (has_valid) {
      cens_mat_ptr[cell] = (val_ptr[i] == NA_LOGICAL) || (val_ptr[i] == 0);
    }
  }

  if (finite_seen != total_cells) {
    return List::create(Named("rectangular") = false);
  }

  NumericVector clock_out(clock_vals.begin(), clock_vals.end());
  IntegerVector u_grp(n_u);
  for (int u = 0; u < n_u; ++u) u_grp[u] = units[u].second;

  return List::create(
    Named("rectangular") = true,
    Named("y") = y_out,
    Named("clock") = clock_out,
    Named("censored") = has_valid ? cens_mat : LogicalVector(0),
    Named("cens_is_matrix") = has_valid,
    Named("unit_group") = u_grp
  );
}

// [[Rcpp::export]]
List align_interpolate_matrix_cpp(
    List signals,
    NumericVector clock,
    LogicalVector censored,
    bool cens_is_matrix,
    IntegerVector unit_group,
    int n_groups,
    NumericVector onsets,
    NumericVector evt_before,
    NumericVector evt_after,
    IntegerVector overlap_n,
    std::string overlap_policy,
    double time_before,
    double time_after,
    double pad_before,
    double pad_after,
    NumericVector tout,
    double max_gap,
    std::string readout_kernel,
    NumericVector quad_offsets,
    NumericVector quad_weights,
    bool trial_grouped,
    bool spatial_first
) {
  const int n_sig = signals.size();
  const int n_time = clock.size();
  const int n_units = unit_group.size();
  const int n_trials = onsets.size();
  const int n_lags = tout.size();

  if (n_sig < 1) stop("At least one signal matrix is required");
  if (n_units < 1) stop("At least one spatial unit is required");
  if (n_groups < 1) stop("At least one spatial group is required");
  if (n_trials < 1) stop("No trial events supplied for alignment");
  if (n_lags < 1) stop("Output time grid must be nonempty");

  const R_xlen_t expected_cells = static_cast<R_xlen_t>(n_time) * static_cast<R_xlen_t>(n_units);
  std::vector<const double*> sig_ptrs(n_sig);
  for (int s = 0; s < n_sig; ++s) {
    SEXP s_sexp = signals[s];
    if (TYPEOF(s_sexp) != REALSXP || Rf_xlength(s_sexp) != expected_cells) {
      stop("Signal matrix dimensions do not match length(clock) * n_units");
    }
    sig_ptrs[s] = REAL(s_sexp);
  }

  const R_xlen_t cens_len = censored.size();
  if (cens_len != 0) {
    if (!cens_is_matrix && cens_len != n_time) {
      stop("Frame censor vector length must match length(clock)");
    }
    if (cens_is_matrix && cens_len != expected_cells) {
      stop("Matrix censor dimensions must match length(clock) * n_units");
    }
  }
  const int* cens_ptr = (cens_len > 0) ? LOGICAL(censored) : nullptr;
  const double* clock_ptr = REAL(clock);
  const double* tout_ptr = REAL(tout);
  const double tol = 1e-8 * std::max(1.0, max_gap);
  const bool use_quad = (readout_kernel != "linear") && (quad_offsets.size() > 0);
  const int n_quad = quad_offsets.size();
  const double* qo_ptr = use_quad ? REAL(quad_offsets) : nullptr;
  const double* qw_ptr = use_quad ? REAL(quad_weights) : nullptr;

  // Partition units by spatial group g in 0 .. n_groups - 1
  std::vector<std::vector<int>> group_units(n_groups);
  for (int u = 0; u < n_units; ++u) {
    const int g = unit_group[u];
    if (g < 0 || g >= n_groups) stop("Invalid unit_group index");
    group_units[g].push_back(u);
  }

  // Step 1: Precompute trial window bounds [k_start[t], k_end[t]) and initial trial statuses
  const double rel_min = time_before + pad_before;
  const double rel_max = time_after + pad_after;
  const bool exclude_overlap = (overlap_policy == "exclude");

  std::vector<int> k_start(n_trials, 0);
  std::vector<int> k_end(n_trials, 0);
  std::vector<std::string> trial_init_status(n_trials, "ok");

  for (int t = 0; t < n_trials; ++t) {
    const double onset = onsets[t];
    if (!std::isfinite(onset)) {
      trial_init_status[t] = "missing_onset";
      continue;
    }
    if (exclude_overlap && overlap_n[t] != NA_INTEGER && overlap_n[t] > 0) {
      trial_init_status[t] = "overlap_excluded";
      continue;
    }

    const double eb = evt_before[t];
    const double ea = evt_after[t];
    if (n_time > 0) {
      const double approx_lo = std::max(onset + rel_min, eb);
      const double approx_hi = std::min(onset + rel_max, ea);
      int lo = static_cast<int>(std::lower_bound(clock_ptr, clock_ptr + n_time, approx_lo) - clock_ptr);
      lo = std::max(0, lo - 2);
      while (lo < n_time) {
        const double ct = clock_ptr[lo];
        const double rt = ct - onset;
        if (rt >= rel_min && ct > eb) break;
        ++lo;
      }

      int hi = static_cast<int>(std::upper_bound(clock_ptr, clock_ptr + n_time, approx_hi) - clock_ptr);
      hi = std::min(n_time, hi + 2);
      while (hi > lo) {
        const double ct = clock_ptr[hi - 1];
        const double rt = ct - onset;
        if (rt <= rel_max && ct < ea) break;
        --hi;
      }
      if (lo < hi) {
        k_start[t] = lo;
        k_end[t] = hi;
      }
    }

    if (k_start[t] >= k_end[t]) {
      if (n_time == 0) {
        trial_init_status[t] = "no_data";
      } else if ((onset + rel_max < clock_ptr[0]) || (onset + rel_min > clock_ptr[n_time - 1])) {
        trial_init_status[t] = "outside_run";
      } else {
        trial_init_status[t] = "censored";
      }
    }
  }

  // Step 2: Precompute shared lag bracket geometry and quadrature brackets per trial
  std::vector<std::vector<double>> trial_rel_t(n_trials);
  std::vector<std::vector<LagSpec>> trial_lags(n_trials, std::vector<LagSpec>(n_lags));
  std::vector<std::vector<int>> trial_quad_left;
  if (use_quad) {
    trial_quad_left.resize(n_trials);
  }

  for (int t = 0; t < n_trials; ++t) {
    const int ks = k_start[t];
    const int ke = k_end[t];
    const int n_win = ke - ks;
    if (n_win <= 0) {
      for (int l = 0; l < n_lags; ++l) {
        LagSpec& ls = trial_lags[t][l];
        ls.mode = MODE_MISSING_EPOCH;
        ls.is_gap = false;
        ls.j_lo = ls.j_hi = ls.r_lo = ls.r_hi = 0;
        ls.src_left = NA_INTEGER;
        ls.src_right = NA_INTEGER;
        ls.weight = 0.0;
      }
      continue;
    }

    const double onset = onsets[t];
    std::vector<double>& rel_t = trial_rel_t[t];
    rel_t.resize(n_win);
    for (int j = 0; j < n_win; ++j) {
      rel_t[j] = clock_ptr[ks + j] - onset;
    }

    for (int l = 0; l < n_lags; ++l) {
      const double x = tout_ptr[l];
      const int left = static_cast<int>(std::upper_bound(rel_t.begin(), rel_t.end(), x) - rel_t.begin());
      const bool near_left = (left > 0) && (std::abs(x - rel_t[left - 1]) <= tol);
      const bool near_right = (left < n_win) && (std::abs(x - rel_t[left]) <= tol);
      LagSpec& ls = trial_lags[t][l];

      if (near_left || near_right) {
        const int j0 = near_left ? (left - 1) : left;
        const int r0 = ks + j0;
        ls.mode = MODE_EXACT;
        ls.is_gap = false;
        ls.j_lo = j0;
        ls.j_hi = j0;
        ls.r_lo = r0;
        ls.r_hi = r0;
        ls.src_left = r0 + 1;
        ls.src_right = NA_INTEGER;
        ls.weight = 0.0;
      } else if (left > 0 && left < n_win) {
        const int j0 = left - 1;
        const int j1 = left;
        const int r0 = ks + j0;
        const int r1 = ks + j1;
        ls.mode = MODE_BRACKET;
        ls.is_gap = (rel_t[j1] - rel_t[j0]) > (max_gap + tol);
        ls.j_lo = j0;
        ls.j_hi = j1;
        ls.r_lo = r0;
        ls.r_hi = r1;
        ls.src_left = r0 + 1;
        ls.src_right = r1 + 1;
        ls.weight = (x - rel_t[j0]) / (rel_t[j1] - rel_t[j0]);
      } else {
        ls.mode = MODE_OUTSIDE;
        ls.is_gap = false;
        ls.j_lo = ls.j_hi = ls.r_lo = ls.r_hi = 0;
        ls.src_left = NA_INTEGER;
        ls.src_right = NA_INTEGER;
        ls.weight = 0.0;
      }
    }

    if (use_quad && n_win > 1) {
      std::vector<int>& ql = trial_quad_left[t];
      ql.resize(static_cast<size_t>(n_lags) * static_cast<size_t>(n_quad));
      for (int l = 0; l < n_lags; ++l) {
        for (int m = 0; m < n_quad; ++m) {
          const double qt = tout_ptr[l] + qo_ptr[m];
          ql[l * n_quad + m] = static_cast<int>(
              std::upper_bound(rel_t.begin(), rel_t.end(), qt) - rel_t.begin());
        }
      }
    }
  }

  // Output row count
  const R_xlen_t n_rows = trial_grouped
      ? (static_cast<R_xlen_t>(n_trials) * static_cast<R_xlen_t>(n_groups) * static_cast<R_xlen_t>(n_lags))
      : (static_cast<R_xlen_t>(n_groups) * static_cast<R_xlen_t>(n_lags));

  auto row_index_tg = [&](int t, int g, int l) -> R_xlen_t {
    const R_xlen_t block = spatial_first
        ? (static_cast<R_xlen_t>(g) * static_cast<R_xlen_t>(n_trials) + static_cast<R_xlen_t>(t))
        : (static_cast<R_xlen_t>(t) * static_cast<R_xlen_t>(n_groups) + static_cast<R_xlen_t>(g));
    return block * static_cast<R_xlen_t>(n_lags) + static_cast<R_xlen_t>(l);
  };

  // Helper to evaluate one unit's estimate and support at (t, l)
  auto eval_unit_lag = [&](const double* y_ptr, int u, int t, int l, double& est_out, SupportCode& sup_out) {
    const LagSpec& ls = trial_lags[t][l];
    if (ls.mode == MODE_MISSING_EPOCH) {
      est_out = NA_REAL;
      sup_out = SUP_MISSING_EPOCH;
      return;
    }
    if (ls.mode == MODE_OUTSIDE) {
      est_out = NA_REAL;
      sup_out = SUP_OUTSIDE;
      return;
    }

    const double* col_u = y_ptr + static_cast<R_xlen_t>(u) * static_cast<R_xlen_t>(n_time);
    int center_j = ls.j_lo;

    if (ls.mode == MODE_EXACT) {
      if (is_cell_censored(cens_ptr, cens_len, cens_is_matrix, ls.r_lo, u, n_time)) {
        est_out = NA_REAL;
        sup_out = SUP_CENSORED;
        return;
      }
      const double v0 = col_u[ls.r_lo];
      if (!std::isfinite(v0)) {
        est_out = NA_REAL;
        sup_out = SUP_MISSING_VALUE;
        return;
      }
      est_out = v0;
      sup_out = SUP_OBSERVED;
    } else {
      const bool c_lo = is_cell_censored(cens_ptr, cens_len, cens_is_matrix, ls.r_lo, u, n_time);
      const bool c_hi = is_cell_censored(cens_ptr, cens_len, cens_is_matrix, ls.r_hi, u, n_time);
      if (c_lo || c_hi) {
        est_out = NA_REAL;
        sup_out = SUP_CENSORED;
        return;
      }
      const double v_lo = col_u[ls.r_lo];
      const double v_hi = col_u[ls.r_hi];
      if (!std::isfinite(v_lo) || !std::isfinite(v_hi)) {
        est_out = NA_REAL;
        sup_out = SUP_MISSING_VALUE;
        return;
      }
      if (ls.is_gap) {
        est_out = NA_REAL;
        sup_out = SUP_GAP;
        return;
      }
      est_out = (1.0 - ls.weight) * v_lo + ls.weight * v_hi;
      sup_out = SUP_INTERPOLATED;
    }

    // Optional local window-integrated readout quadrature within contiguous uncensored segment
    if (use_quad) {
      const int ks = k_start[t];
      const int n_win = k_end[t] - ks;
      if (n_win > 1) {
        const std::vector<double>& rel_t = trial_rel_t[t];
        auto is_usable = [&](int j) -> bool {
          const int r = ks + j;
          return !is_cell_censored(cens_ptr, cens_len, cens_is_matrix, r, u, n_time) &&
                 std::isfinite(col_u[r]);
        };
        int seg_lo = center_j;
        while (seg_lo > 0 && (rel_t[seg_lo] - rel_t[seg_lo - 1] <= max_gap + tol) && is_usable(seg_lo - 1)) {
          --seg_lo;
        }
        int seg_hi = center_j;
        while (seg_hi + 1 < n_win && (rel_t[seg_hi + 1] - rel_t[seg_hi] <= max_gap + tol) && is_usable(seg_hi + 1)) {
          ++seg_hi;
        }
        if (seg_lo == seg_hi) {
          est_out = col_u[ks + seg_lo];
        } else {
          const double t_min = rel_t[seg_lo];
          const double t_max = rel_t[seg_hi];
          const std::vector<int>& ql = trial_quad_left[t];
          double vw_sum = 0.0;
          double w_sum = 0.0;
          for (int m = 0; m < n_quad; ++m) {
            const double qt = tout_ptr[l] + qo_ptr[m];
            if (qt >= t_min - tol && qt <= t_max + tol) {
              double qv;
              if (qt <= t_min) {
                qv = col_u[ks + seg_lo];
              } else if (qt >= t_max) {
                qv = col_u[ks + seg_hi];
              } else {
                const int qb = ql[l * n_quad + m] - 1;
                const double alpha = (qt - rel_t[qb]) / (rel_t[qb + 1] - rel_t[qb]);
                const double y0 = col_u[ks + qb];
                qv = y0 + (col_u[ks + qb + 1] - y0) * alpha;
              }
              const double wm = qw_ptr[m];
              vw_sum += wm * qv;
              w_sum += wm;
            }
          }
          if (w_sum > 0.0) {
            est_out = vw_sum / w_sum;
          }
        }
      }
    }
  };

  // Allocate per-signal output columns
  List sig_results(n_sig);
  // Track whether at least one signal has finite mean at each row (for epoch_status refinement)
  std::vector<unsigned char> any_sig_usable;
  if (trial_grouped) {
    any_sig_usable.assign(n_rows, 0);
  }

  std::vector<double> stat_buf;
  std::vector<int> epoch_count(std::max(1, n_time + 1), 0);
  std::vector<int> last_epoch(std::max(1, n_time + 1), -1);

  for (int s = 0; s < n_sig; ++s) {
    const double* y_ptr = sig_ptrs[s];

    NumericVector col_mean = no_init_vector(n_rows);
    NumericVector col_median = no_init_vector(n_rows);
    NumericVector col_sd = no_init_vector(n_rows);
    IntegerVector col_n = no_init_vector(n_rows);
    IntegerVector col_n_total = no_init_vector(n_rows);
    IntegerVector col_n_censored = no_init_vector(n_rows);
    IntegerVector col_n_missing = no_init_vector(n_rows);
    IntegerVector col_n_gap = no_init_vector(n_rows);
    IntegerVector col_n_outside = no_init_vector(n_rows);
    CharacterVector col_source_ids(n_rows);
    IntegerVector col_source_max_epochs(n_rows);

    double* p_mean = REAL(col_mean);
    double* p_median = REAL(col_median);
    double* p_sd = REAL(col_sd);
    int* p_n = INTEGER(col_n);
    int* p_n_total = INTEGER(col_n_total);
    int* p_n_cens = INTEGER(col_n_censored);
    int* p_n_miss = INTEGER(col_n_missing);
    int* p_n_gap = INTEGER(col_n_gap);
    int* p_n_out = INTEGER(col_n_outside);
    int* p_max_ep = INTEGER(col_source_max_epochs);

    if (trial_grouped) {
      for (int g = 0; g < n_groups; ++g) {
        const std::vector<int>& u_vec = group_units[g];
        const int v_g = static_cast<int>(u_vec.size());
        stat_buf.reserve(v_g);
        std::fill(epoch_count.begin(), epoch_count.end(), 0);
        std::fill(last_epoch.begin(), last_epoch.end(), -1);

        for (int t = 0; t < n_trials; ++t) {
          for (int l = 0; l < n_lags; ++l) {
            const R_xlen_t r_idx = row_index_tg(t, g, l);
            stat_buf.clear();
            int n_cens = 0;
            int n_miss = 0;
            int n_gap_c = 0;
            int n_out_c = 0;

            for (int idx_u = 0; idx_u < v_g; ++idx_u) {
              double est = NA_REAL;
              SupportCode sup = SUP_MISSING_EPOCH;
              eval_unit_lag(y_ptr, u_vec[idx_u], t, l, est, sup);
              if (std::isfinite(est)) {
                stat_buf.push_back(est);
              } else {
                switch (sup) {
                  case SUP_CENSORED: ++n_cens; break;
                  case SUP_MISSING_VALUE:
                  case SUP_MISSING_EPOCH: ++n_miss; break;
                  case SUP_GAP: ++n_gap_c; break;
                  case SUP_OUTSIDE: ++n_out_c; break;
                  default: break;
                }
              }
            }

            const int n_ok = static_cast<int>(stat_buf.size());
            double m_val, med_val, sd_val;
            compute_mean_median_sd(stat_buf, m_val, med_val, sd_val);

            p_mean[r_idx] = m_val;
            p_median[r_idx] = med_val;
            p_sd[r_idx] = sd_val;
            p_n[r_idx] = n_ok;
            p_n_total[r_idx] = v_g;
            p_n_cens[r_idx] = n_cens;
            p_n_miss[r_idx] = n_miss;
            p_n_gap[r_idx] = n_gap_c;
            p_n_out[r_idx] = n_out_c;

            if (n_ok > 0) {
              any_sig_usable[r_idx] = 1;
              const LagSpec& ls = trial_lags[t][l];
              const int sl = ls.src_left;
              const int sr = ls.src_right;
              if (sl != NA_INTEGER) {
                if (last_epoch[sl] != t) {
                  last_epoch[sl] = t;
                  ++epoch_count[sl];
                }
              }
              if (sr != NA_INTEGER) {
                if (last_epoch[sr] != t) {
                  last_epoch[sr] = t;
                  ++epoch_count[sr];
                }
                std::string s_id = std::to_string(sl);
                s_id.push_back(';');
                s_id += std::to_string(sr);
                col_source_ids[r_idx] = s_id;
              } else if (sl != NA_INTEGER) {
                col_source_ids[r_idx] = std::to_string(sl);
              }
            }
          }
        }

        // Populate exact distinct-epoch source reuse count for group g
        for (int t = 0; t < n_trials; ++t) {
          for (int l = 0; l < n_lags; ++l) {
            const R_xlen_t r_idx = row_index_tg(t, g, l);
            if (p_n[r_idx] > 0) {
              const LagSpec& ls = trial_lags[t][l];
              int mx = (ls.src_left != NA_INTEGER) ? epoch_count[ls.src_left] : 0;
              if (ls.src_right != NA_INTEGER && epoch_count[ls.src_right] > mx) {
                mx = epoch_count[ls.src_right];
              }
              p_max_ep[r_idx] = mx;
            } else {
              p_max_ep[r_idx] = 0;
            }
          }
        }
      }
    } else {
      // !trial_grouped: aggregate across all (trial, unit) pairs within each spatial group g
      for (int g = 0; g < n_groups; ++g) {
        const std::vector<int>& u_vec = group_units[g];
        const int v_g = static_cast<int>(u_vec.size());
        const int n_tot_g = n_trials * v_g;
        std::fill(epoch_count.begin(), epoch_count.end(), 0);
        std::fill(last_epoch.begin(), last_epoch.end(), -1);

        std::vector<std::vector<double>> lag_vals(n_lags);
        for (int l = 0; l < n_lags; ++l) lag_vals[l].reserve(n_tot_g);
        std::vector<int> lag_cens(n_lags, 0);
        std::vector<int> lag_miss(n_lags, 0);
        std::vector<int> lag_gap(n_lags, 0);
        std::vector<int> lag_out(n_lags, 0);
        std::vector<std::vector<int>> lag_sources(n_lags);

        for (int t = 0; t < n_trials; ++t) {
          for (int l = 0; l < n_lags; ++l) {
            int ok_in_trial = 0;
            for (int idx_u = 0; idx_u < v_g; ++idx_u) {
              double est = NA_REAL;
              SupportCode sup = SUP_MISSING_EPOCH;
              eval_unit_lag(y_ptr, u_vec[idx_u], t, l, est, sup);
              if (std::isfinite(est)) {
                lag_vals[l].push_back(est);
                ++ok_in_trial;
              } else {
                switch (sup) {
                  case SUP_CENSORED: ++lag_cens[l]; break;
                  case SUP_MISSING_VALUE:
                  case SUP_MISSING_EPOCH: ++lag_miss[l]; break;
                  case SUP_GAP: ++lag_gap[l]; break;
                  case SUP_OUTSIDE: ++lag_out[l]; break;
                  default: break;
                }
              }
            }
            if (ok_in_trial > 0) {
              const LagSpec& ls = trial_lags[t][l];
              if (ls.src_left != NA_INTEGER) {
                lag_sources[l].push_back(ls.src_left);
                if (last_epoch[ls.src_left] != t) {
                  last_epoch[ls.src_left] = t;
                  ++epoch_count[ls.src_left];
                }
              }
              if (ls.src_right != NA_INTEGER) {
                lag_sources[l].push_back(ls.src_right);
                if (last_epoch[ls.src_right] != t) {
                  last_epoch[ls.src_right] = t;
                  ++epoch_count[ls.src_right];
                }
              }
            }
          }
        }

        for (int l = 0; l < n_lags; ++l) {
          const R_xlen_t r_idx = static_cast<R_xlen_t>(g) * static_cast<R_xlen_t>(n_lags) + static_cast<R_xlen_t>(l);
          const int n_ok = static_cast<int>(lag_vals[l].size());
          double m_val, med_val, sd_val;
          compute_mean_median_sd(lag_vals[l], m_val, med_val, sd_val);

          p_mean[r_idx] = m_val;
          p_median[r_idx] = med_val;
          p_sd[r_idx] = sd_val;
          p_n[r_idx] = n_ok;
          p_n_total[r_idx] = n_tot_g;
          p_n_cens[r_idx] = lag_cens[l];
          p_n_miss[r_idx] = lag_miss[l];
          p_n_gap[r_idx] = lag_gap[l];
          p_n_out[r_idx] = lag_out[l];

          std::vector<int>& srcs = lag_sources[l];
          if (!srcs.empty()) {
            std::sort(srcs.begin(), srcs.end());
            srcs.erase(std::unique(srcs.begin(), srcs.end()), srcs.end());
            col_source_ids[r_idx] = join_sorted_ids(srcs);
            int mx = 0;
            for (int sid : srcs) {
              if (epoch_count[sid] > mx) mx = epoch_count[sid];
            }
            p_max_ep[r_idx] = mx;
          } else {
            p_max_ep[r_idx] = 0;
          }
        }
      }
    }

    sig_results[s] = List::create(
      Named("mean") = col_mean,
      Named("median") = col_median,
      Named("sd") = col_sd,
      Named("n") = col_n,
      Named("n_total") = col_n_total,
      Named("n_censored") = col_n_censored,
      Named("n_missing") = col_n_missing,
      Named("n_gap") = col_n_gap,
      Named("n_outside") = col_n_outside,
      Named("source_ids") = col_source_ids,
      Named("source_max_epochs") = col_source_max_epochs
    );
  }

  CharacterVector out_status;
  if (trial_grouped) {
    out_status = CharacterVector(n_rows);
    for (int g = 0; g < n_groups; ++g) {
      for (int t = 0; t < n_trials; ++t) {
        std::string st = trial_init_status[t];
        if (st == "ok") {
          int usable_lags = 0;
          for (int l = 0; l < n_lags; ++l) {
            if (any_sig_usable[row_index_tg(t, g, l)] != 0) ++usable_lags;
          }
          if (usable_lags == 0) {
            st = "insufficient_data";
          } else if (usable_lags < n_lags) {
            st = "partial";
          }
        }
        for (int l = 0; l < n_lags; ++l) {
          out_status[row_index_tg(t, g, l)] = st;
        }
      }
    }
  }

  return List::create(
    Named("signals") = sig_results,
    Named("epoch_status") = out_status
  );
}
