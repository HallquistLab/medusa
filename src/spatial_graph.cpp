#include "medusa.h"
#include <Rmath.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

inline double cauchy_combine_two(double p1, double p2, double w1 = 0.80) {
  const double c1 = std::min(std::max(p1, 1e-15), 0.90);
  const double c2 = std::min(std::max(p2, 1e-15), 0.90);
  const double t_mean = w1 * std::tan((0.5 - c1) * M_PI) + (1.0 - w1) * std::tan((0.5 - c2) * M_PI);
  return std::min(1.0, std::max(1e-16, 0.5 - std::atan(t_mean) / M_PI));
}

inline uint64_t pack_grid_key(int ix, int iy, int iz) {
  const uint64_t ux = static_cast<uint64_t>(static_cast<uint32_t>(ix)) & 0x1FFFFFULL;
  const uint64_t uy = static_cast<uint64_t>(static_cast<uint32_t>(iy)) & 0x1FFFFFULL;
  const uint64_t uz = static_cast<uint64_t>(static_cast<uint32_t>(iz)) & 0x1FFFFFULL;
  return (ux << 42) | (uy << 21) | uz;
}

struct HeapNode {
  int64_t dist_q;
  int node;
  double dist;
  bool operator>(const HeapNode& other) const {
    if (dist_q != other.dist_q) return dist_q > other.dist_q;
    return node > other.node;
  }
};

} // namespace

// Fast deterministic median nearest-neighbor Euclidean distance without O(V^2) memory.
// [[Rcpp::export]]
double estimate_median_nn_dist_cpp(const arma::mat& coords) {
  const arma::uword V = coords.n_rows;
  const arma::uword D = coords.n_cols;
  if (V <= 1 || D == 0) return 1.0;

  const arma::uword M = std::min<arma::uword>(V, 1024);
  std::vector<double> nn_dists;
  nn_dists.reserve(M);

  for (arma::uword q = 0; q < M; ++q) {
    const arma::uword i = (M == V) ? q : ((q * (V - 1)) / (M - 1));
    double best_d2 = std::numeric_limits<double>::infinity();
    for (arma::uword j = 0; j < V; ++j) {
      if (i == j) continue;
      double d2 = 0.0;
      for (arma::uword d = 0; d < D; ++d) {
        const double diff = coords(i, d) - coords(j, d);
        d2 += diff * diff;
      }
      if (d2 > 0.0 && d2 < best_d2) {
        best_d2 = d2;
      }
    }
    if (std::isfinite(best_d2)) {
      nn_dists.push_back(std::sqrt(best_d2));
    }
  }

  if (nn_dists.empty()) return 1.0;
  std::sort(nn_dists.begin(), nn_dists.end());
  const size_t n = nn_dists.size();
  if (n % 2 == 1) {
    return std::max(1e-6, nn_dists[n / 2]);
  }
  return std::max(1e-6, 0.5 * (nn_dists[n / 2 - 1] + nn_dists[n / 2]));
}

// Build a 26-connected in-mask geodesic CSR graph on 3D integer lattice coordinates (i, j, k).
// Out-of-mask gaps (e.g., sulcal CSF or background) act as natural barriers because paths
// can only step through in-mask lattice nodes.
// [[Rcpp::export]]
Rcpp::List build_lattice_geodesic_graph_cpp(
    const Rcpp::IntegerMatrix& ijk,
    const Rcpp::NumericVector& spacing,
    double max_radius,
    int max_neighbors,
    const Rcpp::IntegerVector& group_ids,
    bool hard_boundary) {
  const int V = ijk.nrow();
  if (ijk.ncol() < 3) Rcpp::stop("ijk must have at least 3 columns (i, j, k)");
  if (spacing.size() < 3) Rcpp::stop("spacing must have at least 3 positive elements");
  if (V == 0) {
    return Rcpp::List::create(
      Rcpp::Named("row_ptr") = Rcpp::IntegerVector::create(1),
      Rcpp::Named("col_idx") = Rcpp::IntegerVector(0),
      Rcpp::Named("dist") = Rcpp::NumericVector(0)
    );
  }
  const double sx = spacing[0], sy = spacing[1], sz = spacing[2];
  if (sx <= 0.0 || sy <= 0.0 || sz <= 0.0 || !std::isfinite(sx) || !std::isfinite(sy) || !std::isfinite(sz)) {
    Rcpp::stop("spacing values must be positive and finite");
  }
  const bool check_group = hard_boundary && (group_ids.size() == V);
  const int k_cap = std::max(0, max_neighbors);
  const double r_cap = std::max(0.0, max_radius);

  int min_i = ijk(0, 0), max_i = ijk(0, 0);
  int min_j = ijk(0, 1), max_j = ijk(0, 1);
  int min_k = ijk(0, 2), max_k = ijk(0, 2);
  for (int u = 1; u < V; ++u) {
    const int iu = ijk(u, 0), ju = ijk(u, 1), ku = ijk(u, 2);
    if (iu < min_i) min_i = iu;
    if (iu > max_i) max_i = iu;
    if (ju < min_j) min_j = ju;
    if (ju > max_j) max_j = ju;
    if (ku < min_k) min_k = ku;
    if (ku > max_k) max_k = ku;
  }

  const int64_t nx = static_cast<int64_t>(max_i) - static_cast<int64_t>(min_i) + 1;
  const int64_t ny = static_cast<int64_t>(max_j) - static_cast<int64_t>(min_j) + 1;
  const int64_t nz = static_cast<int64_t>(max_k) - static_cast<int64_t>(min_k) + 1;
  const bool use_flat = (nx > 0 && ny > 0 && nz > 0 &&
                         nx <= 4096 && ny <= 4096 && nz <= 4096 &&
                         (nx * ny) <= 30000000LL &&
                         (nx * ny * nz) <= 30000000LL);

  std::vector<int> flat_grid;
  std::unordered_map<uint64_t, int> hash_grid;
  if (use_flat) {
    flat_grid.assign(static_cast<size_t>(nx * ny * nz), -1);
    for (int u = 0; u < V; ++u) {
      const int64_t idx = (ijk(u, 0) - min_i) +
                          nx * ((ijk(u, 1) - min_j) + ny * (ijk(u, 2) - min_k));
      flat_grid[static_cast<size_t>(idx)] = u;
    }
  } else {
    hash_grid.reserve(static_cast<size_t>(V) * 2);
    for (int u = 0; u < V; ++u) {
      hash_grid[pack_grid_key(ijk(u, 0) - min_i, ijk(u, 1) - min_j, ijk(u, 2) - min_k)] = u;
    }
  }

  int off_i[26], off_j[26], off_k[26];
  double step_len[26];
  int n_off = 0;
  for (int dk = -1; dk <= 1; ++dk) {
    for (int dj = -1; dj <= 1; ++dj) {
      for (int di = -1; di <= 1; ++di) {
        if (di == 0 && dj == 0 && dk == 0) continue;
        off_i[n_off] = di;
        off_j[n_off] = dj;
        off_k[n_off] = dk;
        step_len[n_off] = std::sqrt((di * sx) * (di * sx) +
                                    (dj * sy) * (dj * sy) +
                                    (dk * sz) * (dk * sz));
        ++n_off;
      }
    }
  }

  auto lookup_voxel = [&](int ii, int jj, int kk) -> int {
    if (ii < min_i || ii > max_i || jj < min_j || jj > max_j || kk < min_k || kk > max_k) {
      return -1;
    }
    if (use_flat) {
      const int64_t idx = (ii - min_i) + nx * ((jj - min_j) + ny * (kk - min_k));
      return flat_grid[static_cast<size_t>(idx)];
    }
    const auto it = hash_grid.find(pack_grid_key(ii - min_i, jj - min_j, kk - min_k));
    return (it == hash_grid.end()) ? -1 : it->second;
  };

  const double INF = std::numeric_limits<double>::infinity();
  std::vector<double> dist_buf(static_cast<size_t>(V), INF);
  std::vector<int> touched;
  touched.reserve(256);

  std::vector<int> row_ptr(static_cast<size_t>(V) + 1, 1);
  std::vector<int> col_idx;
  std::vector<double> dist_out;
  col_idx.reserve(static_cast<size_t>(V) * static_cast<size_t>(std::min(k_cap + 1, 27)));
  dist_out.reserve(static_cast<size_t>(V) * static_cast<size_t>(std::min(k_cap + 1, 27)));

  std::vector<std::pair<int, double>> settled;
  settled.reserve(static_cast<size_t>(k_cap + 2));

  for (int u = 0; u < V; ++u) {
    if ((u & 4095) == 0) Rcpp::checkUserInterrupt();
    settled.clear();
    std::priority_queue<HeapNode, std::vector<HeapNode>, std::greater<HeapNode>> pq;
    dist_buf[static_cast<size_t>(u)] = 0.0;
    touched.push_back(u);
    pq.push(HeapNode{0LL, u, 0.0});
    const int gu = check_group ? group_ids[u] : 0;

    while (!pq.empty()) {
      const HeapNode top = pq.top();
      pq.pop();
      const int curr = top.node;
      const double d = top.dist;
      if (d > dist_buf[static_cast<size_t>(curr)] + 1e-12) continue;

      settled.emplace_back(curr, d);
      if (static_cast<int>(settled.size()) >= k_cap + 1) break;
      if (r_cap <= 0.0) break;

      const int ci = ijk(curr, 0), cj = ijk(curr, 1), ck = ijk(curr, 2);
      for (int m = 0; m < 26; ++m) {
        const double nd = d + step_len[m];
        if (nd > r_cap + 1e-9) continue;
        const int nb = lookup_voxel(ci + off_i[m], cj + off_j[m], ck + off_k[m]);
        if (nb < 0) continue;
        if (check_group && group_ids[nb] != gu) continue;
        if (nd < dist_buf[static_cast<size_t>(nb)] - 1e-12) {
          if (!std::isfinite(dist_buf[static_cast<size_t>(nb)])) {
            touched.push_back(nb);
          }
          dist_buf[static_cast<size_t>(nb)] = nd;
          const int64_t dq = static_cast<int64_t>(std::llround(nd * 1e8));
          pq.push(HeapNode{dq, nb, nd});
        }
      }
    }

    for (int idx : touched) {
      dist_buf[static_cast<size_t>(idx)] = INF;
    }
    touched.clear();

    // Ensure self is strictly at position 0 (guaranteed since dist == 0, but explicit)
    if (settled.empty() || settled[0].first != u) {
      settled.insert(settled.begin(), std::make_pair(u, 0.0));
      if (static_cast<int>(settled.size()) > k_cap + 1) {
        settled.resize(static_cast<size_t>(k_cap + 1));
      }
    }
    for (const auto& pr : settled) {
      col_idx.push_back(pr.first + 1); // 1-based index for R
      dist_out.push_back(pr.second);
    }
    row_ptr[static_cast<size_t>(u) + 1] = static_cast<int>(col_idx.size()) + 1;
  }

  return Rcpp::List::create(
    Rcpp::Named("row_ptr") = Rcpp::wrap(row_ptr),
    Rcpp::Named("col_idx") = Rcpp::wrap(col_idx),
    Rcpp::Named("dist") = Rcpp::wrap(dist_out)
  );
}

// Build a sparse radius + k-NN Euclidean CSR graph on arbitrary coordinates using
// a 3D spatial-hash cell linked-list (with exact fallback for D > 3 or tiny V).
// [[Rcpp::export]]
Rcpp::List build_coords_knn_graph_cpp(
    const arma::mat& coords,
    double max_radius,
    int max_neighbors,
    const Rcpp::IntegerVector& group_ids,
    bool hard_boundary) {
  const int V = static_cast<int>(coords.n_rows);
  const int D = static_cast<int>(coords.n_cols);
  if (V == 0) {
    return Rcpp::List::create(
      Rcpp::Named("row_ptr") = Rcpp::IntegerVector::create(1),
      Rcpp::Named("col_idx") = Rcpp::IntegerVector(0),
      Rcpp::Named("dist") = Rcpp::NumericVector(0)
    );
  }
  const bool check_group = hard_boundary && (group_ids.size() == V);
  const int k_cap = std::max(0, max_neighbors);
  const double r_cap = std::max(0.0, max_radius);
  const double r2_cap = r_cap * r_cap;

  std::vector<int> row_ptr(static_cast<size_t>(V) + 1, 1);
  std::vector<int> col_idx;
  std::vector<double> dist_out;
  col_idx.reserve(static_cast<size_t>(V) * static_cast<size_t>(std::min(k_cap + 1, 27)));
  dist_out.reserve(static_cast<size_t>(V) * static_cast<size_t>(std::min(k_cap + 1, 27)));

  const bool use_hash = (V > 512 && D >= 1 && D <= 3 && r_cap > 1e-9);

  std::vector<double> cx, cy, cz;
  double min_x = 0.0, min_y = 0.0, min_z = 0.0;
  double inv_cell = 1.0;
  std::unordered_map<uint64_t, int> cell_head;
  std::vector<int> next_in_cell;

  if (use_hash) {
    cx.assign(static_cast<size_t>(V), 0.0);
    cy.assign(static_cast<size_t>(V), 0.0);
    cz.assign(static_cast<size_t>(V), 0.0);
    for (int u = 0; u < V; ++u) {
      cx[static_cast<size_t>(u)] = coords(u, 0);
      if (D >= 2) cy[static_cast<size_t>(u)] = coords(u, 1);
      if (D >= 3) cz[static_cast<size_t>(u)] = coords(u, 2);
    }
    min_x = *std::min_element(cx.begin(), cx.end());
    min_y = *std::min_element(cy.begin(), cy.end());
    min_z = *std::min_element(cz.begin(), cz.end());
    inv_cell = 1.0 / r_cap;
    cell_head.reserve(static_cast<size_t>(V));
    next_in_cell.assign(static_cast<size_t>(V), -1);
    for (int u = 0; u < V; ++u) {
      const int ix = static_cast<int>(std::floor((cx[static_cast<size_t>(u)] - min_x) * inv_cell));
      const int iy = static_cast<int>(std::floor((cy[static_cast<size_t>(u)] - min_y) * inv_cell));
      const int iz = static_cast<int>(std::floor((cz[static_cast<size_t>(u)] - min_z) * inv_cell));
      const uint64_t key = pack_grid_key(ix, iy, iz);
      auto it = cell_head.find(key);
      if (it == cell_head.end()) {
        cell_head[key] = u;
      } else {
        next_in_cell[static_cast<size_t>(u)] = it->second;
        it->second = u;
      }
    }
  }

  std::vector<HeapNode> cand;
  cand.reserve(128);

  for (int u = 0; u < V; ++u) {
    if ((u & 4095) == 0) Rcpp::checkUserInterrupt();
    cand.clear();
    const int gu = check_group ? group_ids[u] : 0;

    if (k_cap > 0 && r_cap > 0.0) {
      if (use_hash) {
        const double ux = cx[static_cast<size_t>(u)];
        const double uy = cy[static_cast<size_t>(u)];
        const double uz = cz[static_cast<size_t>(u)];
        const int ix = static_cast<int>(std::floor((ux - min_x) * inv_cell));
        const int iy = static_cast<int>(std::floor((uy - min_y) * inv_cell));
        const int iz = static_cast<int>(std::floor((uz - min_z) * inv_cell));
        const int dz_lo = (D >= 3) ? -1 : 0;
        const int dz_hi = (D >= 3) ? 1 : 0;
        const int dy_lo = (D >= 2) ? -1 : 0;
        const int dy_hi = (D >= 2) ? 1 : 0;
        for (int dz = dz_lo; dz <= dz_hi; ++dz) {
          if (iz + dz < 0) continue;
          for (int dy = dy_lo; dy <= dy_hi; ++dy) {
            if (iy + dy < 0) continue;
            for (int dx = -1; dx <= 1; ++dx) {
              if (ix + dx < 0) continue;
              const auto it = cell_head.find(pack_grid_key(ix + dx, iy + dy, iz + dz));
              if (it == cell_head.end()) continue;
              for (int v = it->second; v >= 0; v = next_in_cell[static_cast<size_t>(v)]) {
                if (v == u) continue;
                if (check_group && group_ids[v] != gu) continue;
                const double ddx = ux - cx[static_cast<size_t>(v)];
                const double ddy = uy - cy[static_cast<size_t>(v)];
                const double ddz = uz - cz[static_cast<size_t>(v)];
                const double d2 = ddx * ddx + ddy * ddy + ddz * ddz;
                if (d2 <= r2_cap + 1e-12) {
                  const double d = std::sqrt(std::max(0.0, d2));
                  cand.push_back(HeapNode{static_cast<int64_t>(std::llround(d * 1e8)), v, d});
                }
              }
            }
          }
        }
      } else {
        for (int v = 0; v < V; ++v) {
          if (v == u) continue;
          if (check_group && group_ids[v] != gu) continue;
          double d2 = 0.0;
          for (int d = 0; d < D; ++d) {
            const double diff = coords(u, d) - coords(v, d);
            d2 += diff * diff;
          }
          if (d2 <= r2_cap + 1e-12) {
            const double dist = std::sqrt(std::max(0.0, d2));
            cand.push_back(HeapNode{static_cast<int64_t>(std::llround(dist * 1e8)), v, dist});
          }
        }
      }
    }

    if (static_cast<int>(cand.size()) > k_cap) {
      std::partial_sort(cand.begin(), cand.begin() + k_cap, cand.end(),
                        [](const HeapNode& a, const HeapNode& b) { return b > a; });
      cand.resize(static_cast<size_t>(k_cap));
    } else {
      std::sort(cand.begin(), cand.end(),
                [](const HeapNode& a, const HeapNode& b) { return b > a; });
    }

    col_idx.push_back(u + 1);
    dist_out.push_back(0.0);
    for (const auto& c : cand) {
      col_idx.push_back(c.node + 1);
      dist_out.push_back(c.dist);
    }
    row_ptr[static_cast<size_t>(u) + 1] = static_cast<int>(col_idx.size()) + 1;
  }

  return Rcpp::List::create(
    Rcpp::Named("row_ptr") = Rcpp::wrap(row_ptr),
    Rcpp::Named("col_idx") = Rcpp::wrap(col_idx),
    Rcpp::Named("dist") = Rcpp::wrap(dist_out)
  );
}

// Sparse bilateral smoothing of REML AR(1) log-likelihood profiles over CSR neighborhoods.
// [[Rcpp::export]]
Rcpp::List csr_bilateral_noise_smooth_cpp(
    const arma::mat& profiles,
    const Rcpp::NumericVector& ind_rhos,
    const Rcpp::NumericVector& raw_lsd,
    const Rcpp::LogicalVector& usable,
    const Rcpp::IntegerVector& row_ptr,
    const Rcpp::IntegerVector& col_idx,
    const Rcpp::NumericVector& dist,
    const Rcpp::IntegerVector& group_ids,
    double sigma_sp,
    double cross_weight,
    double rho_scale,
    double lsd_scale) {
  const int G = static_cast<int>(profiles.n_rows);
  const int V = static_cast<int>(profiles.n_cols);
  const int nnz = col_idx.size();
  const bool has_groups = (group_ids.size() == V);
  const double inv_s2 = 1.0 / std::max(sigma_sp * sigma_sp, 1e-12);
  const double inv_rho = 1.0 / std::max(rho_scale, 1e-6);
  const double inv_lsd = 1.0 / std::max(lsd_scale, 1e-6);

  arma::mat smooth_prof(G, V, arma::fill::value(NA_REAL));
  Rcpp::NumericVector edge_weights(nnz, 0.0);
  Rcpp::NumericVector eff_n(V, 0.0);

  std::vector<double> w_tmp;
  w_tmp.reserve(64);

  for (int u = 0; u < V; ++u) {
    if (!usable[u]) continue;
    const int p_start = row_ptr[u] - 1;
    const int p_end = row_ptr[u + 1] - 1;
    const int deg = p_end - p_start;
    if (deg <= 0) continue;

    w_tmp.assign(static_cast<size_t>(deg), 0.0);
    double w_sum = 0.0;
    const double ru = ind_rhos[u];
    const double lu = raw_lsd[u];
    const int gu = has_groups ? group_ids[u] : 0;

    for (int off = 0; off < deg; ++off) {
      const int p = p_start + off;
      const int v = col_idx[p] - 1;
      if (v < 0 || v >= V || !usable[v]) continue;
      double b_mult = 1.0;
      if (has_groups && gu > 0 && group_ids[v] != gu) {
        b_mult = cross_weight;
        if (b_mult <= 0.0) continue;
      }
      const double dp = dist[p];
      const double w_sp = std::exp(-0.5 * dp * dp * inv_s2) * b_mult;
      const double dr = (ind_rhos[v] - ru) * inv_rho;
      const double dl = (raw_lsd[v] - lu) * inv_lsd;
      const double w_sim = std::exp(-0.5 * (dr * dr + 0.5 * dl * dl));
      const double w = w_sp * w_sim;
      w_tmp[static_cast<size_t>(off)] = w;
      w_sum += w;
    }

    if (w_sum <= 1e-15) {
      smooth_prof.col(u) = profiles.col(u);
      edge_weights[p_start] = 1.0;
      eff_n[u] = 1.0;
      continue;
    }

    const double inv_sum = 1.0 / w_sum;
    double w_sq_sum = 0.0;
    arma::vec acc(G, arma::fill::zeros);
    for (int off = 0; off < deg; ++off) {
      const double wn = w_tmp[static_cast<size_t>(off)] * inv_sum;
      if (wn <= 0.0) continue;
      const int p = p_start + off;
      const int v = col_idx[p] - 1;
      edge_weights[p] = wn;
      w_sq_sum += wn * wn;
      acc += wn * profiles.col(v);
    }
    smooth_prof.col(u) = acc;
    eff_n[u] = (w_sq_sum > 0.0) ? (1.0 / w_sq_sum) : 1.0;
  }

  return Rcpp::List::create(
    Rcpp::Named("smooth_profiles") = smooth_prof,
    Rcpp::Named("edge_weights") = edge_weights,
    Rcpp::Named("eff_n") = eff_n
  );
}

// Aggregate CSR normalized bilateral weights per noise stratum without O(V^2) affinity matrices.
// [[Rcpp::export]]
arma::mat csr_cluster_member_weights_cpp(
    const Rcpp::IntegerVector& cluster_id,
    int n_clusters,
    const Rcpp::LogicalVector& usable,
    const Rcpp::IntegerVector& row_ptr,
    const Rcpp::IntegerVector& col_idx,
    const Rcpp::NumericVector& edge_weights) {
  const int V = usable.size();
  arma::mat out(V, std::max(1, n_clusters), arma::fill::zeros);
  std::vector<int> counts(static_cast<size_t>(std::max(1, n_clusters)), 0);

  for (int u = 0; u < V; ++u) {
    if (!usable[u]) continue;
    const int cid = cluster_id[u] - 1;
    if (cid < 0 || cid >= n_clusters) continue;
    counts[static_cast<size_t>(cid)] += 1;
    const int p_start = row_ptr[u] - 1;
    const int p_end = row_ptr[u + 1] - 1;
    for (int p = p_start; p < p_end; ++p) {
      const double w = edge_weights[p];
      if (w <= 0.0) continue;
      const int v = col_idx[p] - 1;
      if (v >= 0 && v < V && usable[v]) {
        out(v, cid) += w;
      }
    }
  }
  for (int c = 0; c < n_clusters; ++c) {
    if (counts[static_cast<size_t>(c)] > 0) {
      out.col(c) /= static_cast<double>(counts[static_cast<size_t>(c)]);
    }
  }
  return out;
}

// Debiased quadratic spatial innovation correlation v_quad(u) = sum_{v,v'} w_{uv} w_{uv'} R_{sp,vv'}^2
// computed per run in O(m * k^2 * V) via tiny k x k cross-products of unit-norm whitened residuals.
// [[Rcpp::export]]
Rcpp::NumericVector csr_quad_cor_cpp(
    const arma::mat& Z_unit,
    double d_res,
    const Rcpp::IntegerVector& row_ptr,
    const Rcpp::IntegerVector& col_idx,
    const Rcpp::NumericVector& dist,
    const Rcpp::IntegerVector& group_ids,
    const Rcpp::NumericVector& rho_vec,
    const Rcpp::NumericVector& lsd_vec,
    const Rcpp::LogicalVector& usable_col,
    double sigma_sp,
    double cross_weight,
    double rho_scale,
    double lsd_scale) {
  const int m = static_cast<int>(Z_unit.n_rows);
  const int V = static_cast<int>(Z_unit.n_cols);
  Rcpp::NumericVector out(V, 1.0);
  if (m < 4 || d_res <= 2.0 || V == 0) return out;

  const bool has_groups = (group_ids.size() == V);
  const bool has_noise = (rho_vec.size() == V && lsd_vec.size() == V);
  const double inv_s2 = 1.0 / std::max(sigma_sp * sigma_sp, 1e-12);
  const double inv_rho = 1.0 / std::max(rho_scale, 1e-6);
  const double inv_lsd = 1.0 / std::max(lsd_scale, 1e-6);

  std::vector<int> nbrs;
  std::vector<double> w_loc;
  nbrs.reserve(64);
  w_loc.reserve(64);

  for (int u = 0; u < V; ++u) {
    if (!usable_col[u]) {
      out[u] = 1.0;
      continue;
    }
    const int p_start = row_ptr[u] - 1;
    const int p_end = row_ptr[u + 1] - 1;
    nbrs.clear();
    w_loc.clear();
    const int gu = has_groups ? group_ids[u] : 0;
    const double ru = has_noise ? rho_vec[u] : 0.0;
    const double lu = has_noise ? lsd_vec[u] : 0.0;
    const bool ok_nu = has_noise && std::isfinite(ru) && std::isfinite(lu);
    double w_sum = 0.0;

    for (int p = p_start; p < p_end; ++p) {
      const int v = col_idx[p] - 1;
      if (v < 0 || v >= V || !usable_col[v]) continue;
      double b_mult = 1.0;
      if (has_groups && gu > 0 && group_ids[v] != gu) {
        b_mult = cross_weight;
        if (b_mult <= 0.0) continue;
      }
      const double dp = dist[p];
      double w = std::exp(-0.5 * dp * dp * inv_s2) * b_mult;
      if (ok_nu && std::isfinite(rho_vec[v]) && std::isfinite(lsd_vec[v])) {
        const double dr = (rho_vec[v] - ru) * inv_rho;
        const double dl = (lsd_vec[v] - lu) * inv_lsd;
        w *= std::exp(-0.5 * (dr * dr + 0.5 * dl * dl));
      }
      if (w > 1e-12) {
        nbrs.push_back(v);
        w_loc.push_back(w);
        w_sum += w;
      }
    }

    const int ku = static_cast<int>(nbrs.size());
    if (ku <= 1 || w_sum <= 1e-12) {
      out[u] = 1.0;
      continue;
    }

    double sw2 = 0.0;
    arma::mat Z_u(m, ku);
    for (int a = 0; a < ku; ++a) {
      w_loc[static_cast<size_t>(a)] /= w_sum;
      const double wa = w_loc[static_cast<size_t>(a)];
      sw2 += wa * wa;
      Z_u.col(a) = std::sqrt(wa) * Z_unit.col(nbrs[static_cast<size_t>(a)]);
    }
    const arma::mat G_u = Z_u.t() * Z_u;
    const double raw_f2 = arma::accu(G_u % G_u);
    double diag_f2 = 0.0;
    for (int a = 0; a < ku; ++a) {
      const double gaa = G_u(a, a);
      diag_f2 += gaa * gaa;
    }
    const double off_f2 = std::max(0.0, raw_f2 - diag_f2);
    const double v_est = sw2 + (d_res / (d_res - 1.0)) * off_f2 - (1.0 - sw2) / (d_res - 1.0);
    out[u] = std::min(1.0, std::max(sw2, v_est));
  }

  return out;
}

// Sparse bilateral leave-target-run-out HRF evidence gating and Boltzmann soft/hard pooling.
// [[Rcpp::export]]
Rcpp::List csr_hrf_graph_pool_cpp(
    const arma::mat& scores,
    const Rcpp::List& train_scores_list,
    bool has_gate,
    const arma::mat& task_p_raw,
    const arma::mat& shape_p_raw,
    const Rcpp::LogicalMatrix& compatible,
    const Rcpp::LogicalMatrix& task_supported,
    const Rcpp::LogicalMatrix& supported,
    double eff_task,
    double eff_shape,
    const Rcpp::IntegerVector& counts,
    int min_training_runs,
    const Rcpp::IntegerVector& row_ptr,
    const Rcpp::IntegerVector& col_idx,
    const Rcpp::NumericVector& dist,
    const Rcpp::IntegerVector& group_ids,
    int boundary_mode,
    double cross_weight,
    double sigma_sp,
    const Rcpp::NumericVector& rho_train,
    const Rcpp::NumericVector& lsd_train,
    const Rcpp::NumericVector& v_quad_fold,
    const arma::mat& kernel_cos,
    const Rcpp::NumericVector& dist_fb,
    int fallback_kernel,
    double evidence_alpha,
    double min_r2,
    double canonical_strength,
    bool soft_average) {
  const int K = static_cast<int>(scores.n_rows);
  const int V = static_cast<int>(scores.n_cols);
  const int fb = fallback_kernel - 1;
  const int S = train_scores_list.size();
  const bool has_groups = (group_ids.size() == V);
  const bool has_noise = (rho_train.size() == V && lsd_train.size() == V);
  const double inv_s2 = 1.0 / std::max(sigma_sp * sigma_sp, 1e-12);

  // Noise scale for bilateral downweighting across CSF/vein vs GM boundaries
  double sd_lsd = 0.50;
  if (has_noise && V > 2) {
    double s1 = 0.0, s2 = 0.0;
    int cnt_l = 0;
    for (int v = 0; v < V; ++v) {
      if (std::isfinite(lsd_train[v])) {
        s1 += lsd_train[v];
        s2 += lsd_train[v] * lsd_train[v];
        ++cnt_l;
      }
    }
    if (cnt_l > 2) {
      const double mu = s1 / cnt_l;
      const double vr = std::max(0.0, (s2 - cnt_l * mu * mu) / (cnt_l - 1));
      sd_lsd = std::max(0.50, std::sqrt(vr));
    }
  }
  const double inv_rho = 1.0 / 0.25;
  const double inv_lsd = 1.0 / sd_lsd;

  // Precompute centered, unit-norm per-run candidate score profiles d_tilde^(s) (K x V)
  std::vector<arma::mat> d_runs(static_cast<size_t>(S));
  std::vector<std::vector<uint8_t>> d_valid(static_cast<size_t>(S), std::vector<uint8_t>(static_cast<size_t>(V), 0));
  for (int s = 0; s < S; ++s) {
    arma::mat ms = Rcpp::as<arma::mat>(train_scores_list[s]);
    arma::mat ds(K, V, arma::fill::zeros);
    if (K >= 2) {
      for (int v = 0; v < V; ++v) {
        double sum_k = 0.0;
        bool ok = true;
        for (int k = 0; k < K; ++k) {
          const double val = ms(k, v);
          if (!std::isfinite(val)) {
            ok = false;
            break;
          }
          sum_k += val;
        }
        if (!ok) continue;
        const double mean_k = sum_k / static_cast<double>(K);
        double nrm2 = 0.0;
        for (int k = 0; k < K; ++k) {
          const double c = ms(k, v) - mean_k;
          ds(k, v) = c;
          nrm2 += c * c;
        }
        if (nrm2 > 1e-16) {
          ds.col(v) *= (1.0 / std::sqrt(nrm2));
          d_valid[static_cast<size_t>(s)][static_cast<size_t>(v)] = 1;
        }
      }
    }
    d_runs[static_cast<size_t>(s)] = std::move(ds);
  }

  // Precompute voxelwise task informativeness eta_v in [0, 1] and cross-run
  // HRF profile self-reliability r_self_sm[v] from training runs only.
  std::vector<double> eta(static_cast<size_t>(V), 0.0);
  std::vector<double> r_self(static_cast<size_t>(V), 0.0);
  std::vector<uint8_t> r_has(static_cast<size_t>(V), 0);
  std::vector<uint8_t> elig_mask(static_cast<size_t>(V), 0);
  arma::mat z_task_mat, z_shape_mat, t_task_mat, t_shape_mat;
  if (has_gate) {
    z_task_mat.set_size(K, V);
    z_shape_mat.set_size(K, V);
    t_task_mat.set_size(K, V);
    t_shape_mat.set_size(K, V);
  }

  const double r2_floor = std::max(min_r2, 1e-4);
  arma::vec sum_self_vec(K);
  for (int v = 0; v < V; ++v) {
    if (counts[v] < min_training_runs || !std::isfinite(scores(0, v))) continue;
    elig_mask[static_cast<size_t>(v)] = 1;
    double r_max = -std::numeric_limits<double>::infinity();
    double p_min = 1.0;
    for (int k = 0; k < K; ++k) {
      const double sc = scores(k, v);
      if (std::isfinite(sc) && sc > r_max) r_max = sc;
      if (has_gate) {
        const double pt_in = task_p_raw(k, v);
        const double ps_in = shape_p_raw(k, v);
        const double pt = std::isfinite(pt_in) ? std::min(std::max(pt_in, 1e-15), 0.99) : 0.99;
        const double ps = std::isfinite(ps_in) ? std::min(std::max(ps_in, 1e-15), 0.99) : 0.99;
        if (pt < p_min) p_min = pt;
        z_task_mat(k, v) = R::qnorm5(1.0 - pt, 0.0, 1.0, 1, 0);
        z_shape_mat(k, v) = R::qnorm5(1.0 - ps, 0.0, 1.0, 1, 0);
        const double pt_c = std::min(pt, 0.90);
        const double ps_c = std::min(ps, 0.90);
        t_task_mat(k, v) = std::tan((0.5 - pt_c) * M_PI);
        t_shape_mat(k, v) = std::tan((0.5 - ps_c) * M_PI);
      }
    }
    const double eta_r2 = std::min(1.0, std::max(0.0, (r_max - 0.5 * r2_floor) / (2.0 * r2_floor)));
    if (has_gate) {
      const double neg_log_p = -std::log10(std::max(p_min, 1e-12));
      const double eta_p = std::min(1.0, std::max(0.0, (neg_log_p - 1.30) / 2.20));
      eta[static_cast<size_t>(v)] = eta_r2 * eta_p;
    } else {
      eta[static_cast<size_t>(v)] = eta_r2;
    }

    if (K >= 2 && S >= 2) {
      int n_sv = 0;
      sum_self_vec.zeros();
      for (int s = 0; s < S; ++s) {
        if (d_valid[static_cast<size_t>(s)][static_cast<size_t>(v)]) {
          sum_self_vec += d_runs[static_cast<size_t>(s)].col(v);
          ++n_sv;
        }
      }
      if (n_sv >= 2) {
        const double norm2 = arma::dot(sum_self_vec, sum_self_vec);
        r_self[static_cast<size_t>(v)] = (norm2 - static_cast<double>(n_sv)) /
          static_cast<double>(n_sv * (n_sv - 1));
        r_has[static_cast<size_t>(v)] = 1;
      }
    }
  }

  // Spatially stabilize cross-run self-reliability across baseline geometry/noise neighbors
  // (positive on both sides of any functional boundary, ~0 under pure noise)
  std::vector<double> r_self_sm(static_cast<size_t>(V), 0.0);
  std::vector<double> gamma_sh(static_cast<size_t>(V), 0.0);
  for (int u = 0; u < V; ++u) {
    if (!elig_mask[static_cast<size_t>(u)]) continue;
    if (K >= 2 && S >= 2) {
      const int p_s = row_ptr[u] - 1;
      const int p_e = row_ptr[u + 1] - 1;
      const double ru = has_noise ? rho_train[u] : 0.0;
      const double lu = has_noise ? lsd_train[u] : 0.0;
      const bool ok_nu = has_noise && std::isfinite(ru) && std::isfinite(lu);
      double acc_r = 0.0, acc_w = 0.0;
      for (int p = p_s; p < p_e; ++p) {
        const int v = col_idx[p] - 1;
        if (v < 0 || v >= V || !r_has[static_cast<size_t>(v)]) continue;
        const double dp = dist[p];
        double w = std::exp(-0.5 * dp * dp * inv_s2);
        if (ok_nu && std::isfinite(rho_train[v]) && std::isfinite(lsd_train[v])) {
          const double dr = (rho_train[v] - ru) * inv_rho;
          const double dl = (lsd_train[v] - lu) * inv_lsd;
          w *= std::exp(-0.5 * (dr * dr + 0.5 * dl * dl));
        }
        acc_r += w * r_self[static_cast<size_t>(v)];
        acc_w += w;
      }
      const double r_nb = (acc_w > 1e-12) ? (acc_r / acc_w) : r_self[static_cast<size_t>(u)];
      const double r_comb = r_has[static_cast<size_t>(u)]
        ? (0.35 * r_self[static_cast<size_t>(u)] + 0.65 * r_nb)
        : r_nb;
      r_self_sm[static_cast<size_t>(u)] = std::max(0.0, std::min(0.95, r_comb));
      gamma_sh[static_cast<size_t>(u)] = eta[static_cast<size_t>(u)] *
        std::min(1.0, std::max(0.0, (r_self_sm[static_cast<size_t>(u)] - 0.04) / 0.24));
    } else if (K >= 2 && S == 1) {
      r_self_sm[static_cast<size_t>(u)] = 0.55;
      gamma_sh[static_cast<size_t>(u)] = 0.70 * eta[static_cast<size_t>(u)];
    }
  }

  // Empirical-Bayes departure ramp thresholds calibrated to library multiplicity eff_shape
  const double z_eb_0 = 1.40 + 0.15 * std::log(std::max(1.0, eff_shape / 4.0));
  const double z_eb_1 = z_eb_0 + 1.10;
  const double z_eb_span = z_eb_1 - z_eb_0;

  arma::mat w_fold(K, V, arma::fill::zeros);
  w_fold.row(fb).ones();
  Rcpp::IntegerVector chosen(V, fb + 1);
  Rcpp::CharacterVector status(V, "fallback_prespecified");
  Rcpp::NumericVector eff_n_out(V, 0.0);
  Rcpp::NumericVector omega_out(V, 0.0);

  std::vector<int> nbrs;
  std::vector<double> w_base, w_bilat;
  nbrs.reserve(64);
  w_base.reserve(64);
  w_bilat.reserve(64);

  arma::vec pooled_bilat(K);
  std::vector<double> p_task_pool(static_cast<size_t>(K), 1.0);
  std::vector<double> p_shape_pool(static_cast<size_t>(K), 1.0);
  std::vector<double> p_shape_raw_comb(static_cast<size_t>(K), 1.0);
  std::vector<double> w_compat(static_cast<size_t>(K), 0.0);
  std::vector<double> w_task_sup(static_cast<size_t>(K), 0.0);
  std::vector<double> w_beat_fb(static_cast<size_t>(K), 0.0);
  std::vector<uint8_t> admissible(static_cast<size_t>(K), 0);
  std::vector<uint8_t> comp_ok(static_cast<size_t>(K), 0);

  for (int u = 0; u < V; ++u) {
    if ((u & 4095) == 0) Rcpp::checkUserInterrupt();
    if (!elig_mask[static_cast<size_t>(u)]) continue;

    nbrs.clear();
    w_base.clear();
    w_bilat.clear();
    const int p_start = row_ptr[u] - 1;
    const int p_end = row_ptr[u + 1] - 1;
    const int gu = has_groups ? group_ids[u] : 0;
    const double ru = has_noise ? rho_train[u] : 0.0;
    const double lu = has_noise ? lsd_train[u] : 0.0;
    const bool ok_nu = has_noise && std::isfinite(ru) && std::isfinite(lu);
    const double eta_u = eta[static_cast<size_t>(u)];
    const double gamma_u = gamma_sh[static_cast<size_t>(u)];
    const double r_uu = r_self_sm[static_cast<size_t>(u)];

    double sum_wb = 0.0, sum_wbl = 0.0;

    for (int p = p_start; p < p_end; ++p) {
      const int v = col_idx[p] - 1;
      if (v < 0 || v >= V || !elig_mask[static_cast<size_t>(v)]) continue;

      double b_base = 1.0;
      const bool cross_grp = (has_groups && gu > 0 && group_ids[v] != gu);
      if (cross_grp) {
        if (boundary_mode == 2) continue; // hard parcel boundary
        if (boundary_mode == 1) b_base = cross_weight;
      }
      if (b_base <= 0.0) continue;

      const double dp = dist[p];
      const double w_geo = std::exp(-0.5 * dp * dp * inv_s2);
      double w_ns = 1.0;
      if (ok_nu && std::isfinite(rho_train[v]) && std::isfinite(lsd_train[v])) {
        const double dr = (rho_train[v] - ru) * inv_rho;
        const double dl = (lsd_train[v] - lu) * inv_lsd;
        w_ns = std::exp(-0.5 * (dr * dr + 0.5 * dl * dl));
      }

      const double wb = w_geo * w_ns * b_base;
      if (wb <= 1e-12) continue;

      double wbl = wb;
      if (v != u) {
        double c_uv = 0.0;
        int n_sh = 0;
        if (K >= 2 && S >= 1) {
          double sum_diag = 0.0;
          arma::vec sum_u(K, arma::fill::zeros);
          arma::vec sum_v(K, arma::fill::zeros);
          for (int s = 0; s < S; ++s) {
            if (d_valid[static_cast<size_t>(s)][static_cast<size_t>(u)] &&
                d_valid[static_cast<size_t>(s)][static_cast<size_t>(v)]) {
              const arma::mat& ds = d_runs[static_cast<size_t>(s)];
              const double dot_s = arma::dot(ds.col(u), ds.col(v));
              sum_diag += dot_s;
              sum_u += ds.col(u);
              sum_v += ds.col(v);
              ++n_sh;
            }
          }
          if (n_sh >= 2) {
            const double sum_all = arma::dot(sum_u, sum_v);
            c_uv = (sum_all - sum_diag) / static_cast<double>(n_sh * (n_sh - 1));
          } else if (n_sh == 1) {
            c_uv = sum_diag;
          }
        } else if (K == 1) {
          c_uv = 1.0;
        }
        const double r_vv = r_self_sm[static_cast<size_t>(v)];
        const double denom = (K >= 2 && n_sh >= 1)
          ? std::max(0.22, std::sqrt(std::max(0.06, r_uu) * std::max(0.06, r_vv)))
          : 1.0;
        const double c_tilde = std::max(-1.0, std::min(1.0, c_uv / denom));
        const double c_pos = std::max(0.0, c_tilde);
        const double c_neg = std::max(0.0, -c_tilde);
        const double eta_v = eta[static_cast<size_t>(v)];
        const double gamma_v = gamma_sh[static_cast<size_t>(v)];
        const double gamma_pair = std::sqrt(gamma_u * gamma_v);
        const double s_act = (eta_u > eta_v)
          ? ((1.0 - eta_u) + eta_u * ((eta_v + 0.05) / (eta_u + 0.05)))
          : 1.0;
        const double s_conc = ((1.0 - gamma_pair) + gamma_pair * c_pos) *
                              std::exp(-3.5 * gamma_pair * c_neg * c_neg);
        double b_eff = 1.0;
        if (cross_grp && boundary_mode == 1) {
          b_eff = cross_weight + (1.0 - cross_weight) * gamma_pair * c_pos;
          b_eff *= std::exp(-2.5 * gamma_pair * c_neg * c_neg);
        }
        wbl = w_geo * w_ns * b_eff * s_act * s_conc;
      }

      nbrs.push_back(v);
      w_base.push_back(wb);
      w_bilat.push_back(wbl);
      sum_wb += wb;
      sum_wbl += wbl;
    }

    const int ku = static_cast<int>(nbrs.size());
    if (ku == 0 || sum_wb <= 1e-15 || sum_wbl <= 1e-15) continue;

    double swb2 = 0.0, swbl2 = 0.0;
    for (int a = 0; a < ku; ++a) {
      w_base[static_cast<size_t>(a)] /= sum_wb;
      w_bilat[static_cast<size_t>(a)] /= sum_wbl;
      swb2 += w_base[static_cast<size_t>(a)] * w_base[static_cast<size_t>(a)];
      swbl2 += w_bilat[static_cast<size_t>(a)] * w_bilat[static_cast<size_t>(a)];
    }
    eff_n_out[u] = (swbl2 > 0.0) ? (1.0 / swbl2) : 1.0;

    pooled_bilat.zeros();
    for (int a = 0; a < ku; ++a) {
      pooled_bilat += w_bilat[static_cast<size_t>(a)] * scores.col(nbrs[static_cast<size_t>(a)]);
    }

    if (!pooled_bilat.is_finite() || pooled_bilat.max() < min_r2) {
      continue;
    }

    if (!has_gate) {
      for (int k = 0; k < K; ++k) {
        admissible[static_cast<size_t>(k)] = (pooled_bilat(k) >= min_r2) ? 1 : 0;
        comp_ok[static_cast<size_t>(k)] = (k == fb || pooled_bilat(k) >= min_r2) ? 1 : 0;
      }
    } else {
      const double vq = (v_quad_fold.size() == V && std::isfinite(v_quad_fold[u]))
        ? std::min(1.0, std::max(swb2, v_quad_fold[u]))
        : std::min(1.0, std::max(swb2, 0.50));
      const double sd_st = std::sqrt(std::max(1e-8, vq));
      const double sd_st_bl = std::sqrt(std::max(1e-8, std::min(1.0, std::max(swbl2, vq))));

      for (int k = 0; k < K; ++k) {
        double zt_acc = 0.0, zs_acc = 0.0, tt_acc = 0.0, ts_acc = 0.0;
        double wc_acc = 0.0, wts_acc = 0.0, wpl_acc = 0.0, wbf_acc = 0.0;
        for (int a = 0; a < ku; ++a) {
          const int v = nbrs[static_cast<size_t>(a)];
          const double wb = w_base[static_cast<size_t>(a)];
          const double wbl = w_bilat[static_cast<size_t>(a)];
          zt_acc += wb * z_task_mat(k, v);
          zs_acc += wbl * z_shape_mat(k, v);
          tt_acc += wb * t_task_mat(k, v);
          ts_acc += wbl * t_shape_mat(k, v);
          if (compatible(k, v)) wc_acc += wb;
          if (task_supported(k, v)) wts_acc += wb;
          if ((scores(k, v) >= 0.5 * min_r2) ||
              (std::isfinite(task_p_raw(k, v)) && task_p_raw(k, v) <= 0.20)) {
            wpl_acc += wb;
          }
          if (scores(k, v) > scores(fb, v)) wbf_acc += wbl;
        }
        const double pt_st = R::pnorm5(zt_acc / sd_st, 0.0, 1.0, 0, 0);
        const double ps_st = R::pnorm5(zs_acc / sd_st_bl, 0.0, 1.0, 0, 0);
        const double pt_ca = std::min(1.0, std::max(1e-16, 0.5 - std::atan(tt_acc) / M_PI));
        const double ps_ca = std::min(1.0, std::max(1e-16, 0.5 - std::atan(ts_acc) / M_PI));
        const double pt_raw = cauchy_combine_two(pt_st, pt_ca, 0.80);
        const double ps_raw = (k == fb) ? 0.0 : cauchy_combine_two(ps_st, ps_ca, 0.80);

        p_task_pool[static_cast<size_t>(k)] = std::min(1.0, pt_raw * eff_task);
        p_shape_raw_comb[static_cast<size_t>(k)] = ps_raw;
        p_shape_pool[static_cast<size_t>(k)] = (k == fb) ? 0.0 : std::min(1.0, ps_raw * eff_shape);
        w_compat[static_cast<size_t>(k)] = wc_acc;
        w_task_sup[static_cast<size_t>(k)] = wts_acc;
        w_beat_fb[static_cast<size_t>(k)] = wbf_acc;

        const bool comp_u = static_cast<bool>(compatible(k, u));
        const bool self_plausible = (scores(k, u) >= 0.5 * min_r2) ||
                                    (std::isfinite(task_p_raw(k, u)) && task_p_raw(k, u) <= 0.20);
        const bool halo_ok = self_plausible || ((wts_acc >= 0.50 || wpl_acc >= 0.65) && wbf_acc >= 0.55);
        const bool adm_self = static_cast<bool>(supported(k, u)) && (scores(k, u) >= min_r2);

        if (k == fb) {
          admissible[static_cast<size_t>(k)] = (comp_u && pooled_bilat(fb) >= min_r2 &&
            (p_task_pool[static_cast<size_t>(fb)] <= evidence_alpha || task_supported(fb, u))) ? 1 : 0;
          comp_ok[static_cast<size_t>(k)] = 1;
        } else {
          const bool adm_gr = comp_u && (wc_acc >= 0.60) && halo_ok &&
                              (pooled_bilat(k) >= min_r2) &&
                              (pooled_bilat(k) > pooled_bilat(fb)) &&
                              (wbf_acc >= 0.35) &&
                              (p_task_pool[static_cast<size_t>(k)] <= evidence_alpha) &&
                              (p_shape_pool[static_cast<size_t>(k)] <= evidence_alpha);
          admissible[static_cast<size_t>(k)] = (adm_self || adm_gr) ? 1 : 0;
          comp_ok[static_cast<size_t>(k)] = (adm_self ||
            (comp_u && (wc_acc >= 0.60) && halo_ok &&
             (pooled_bilat(k) >= min_r2) &&
             (p_task_pool[static_cast<size_t>(k)] <= evidence_alpha || wts_acc >= 0.50))) ? 1 : 0;
        }
      }
    }

    // Hard candidate selection
    int best_hard = fb;
    double best_hard_obj = -std::numeric_limits<double>::infinity();
    for (int k = 0; k < K; ++k) {
      if (!admissible[static_cast<size_t>(k)]) continue;
      const double obj_k = pooled_bilat(k) - canonical_strength * dist_fb[k];
      if (obj_k > best_hard_obj) {
        best_hard_obj = obj_k;
        best_hard = k;
      }
    }

    if (!soft_average || K <= 1) {
      chosen[u] = best_hard + 1;
      w_fold(fb, u) = 0.0;
      w_fold(best_hard, u) = 1.0;
      if (best_hard != fb || admissible[static_cast<size_t>(fb)]) {
        status[u] = (eff_n_out[u] > 1.05) ? "pooled_graph" : "voxel";
      }
      omega_out[u] = (best_hard != fb) ? 1.0 : 0.0;
      continue;
    }

    // Check if center voxel has no compatible shapes at all
    if (has_gate) {
      bool any_comp_u = false;
      for (int k = 0; k < K; ++k) {
        if (compatible(k, u)) {
          any_comp_u = true;
          break;
        }
      }
      if (!any_comp_u) continue;
    }

    // Empirical-Bayes departure weight omega_u for soft_average
    double omega_u = 0.0;
    if (!has_gate) {
      omega_u = 1.0;
    } else {
      int best_c = -1;
      double best_c_obj = -std::numeric_limits<double>::infinity();
      for (int k = 0; k < K; ++k) {
        if (k == fb || !comp_ok[static_cast<size_t>(k)]) continue;
        const double obj_k = pooled_bilat(k) - canonical_strength * dist_fb[k];
        if (obj_k > best_c_obj) {
          best_c_obj = obj_k;
          best_c = k;
        }
      }
      if (best_c >= 0 && best_c_obj > pooled_bilat(fb) && w_beat_fb[static_cast<size_t>(best_c)] >= 0.35) {
        if (supported(best_c, u) && scores(best_c, u) >= min_r2) {
          omega_u = 1.0;
        } else {
          const double ps_raw = std::min(0.999999, std::max(1e-15, p_shape_raw_comb[static_cast<size_t>(best_c)]));
          const double z_s = R::qnorm5(1.0 - ps_raw, 0.0, 1.0, 1, 0);

          // Bilateral neighborhood score-margin t-evidence (handles high-rank trial
          // designs where per-run encompassing F-tests have diluted single-voxel power)
          double z_tdiff = 0.0;
          if (ku >= 2) {
            const double c_fb_best = kernel_cos(fb, best_c);
            const double norm_best = std::sqrt(std::max(0.0, 1.0 - c_fb_best * c_fb_best));
            int opp_c = -1;
            double min_cos_tau = 1.0;
            if (norm_best > 1e-6) {
              for (int k = 0; k < K; ++k) {
                if (k == fb || !std::isfinite(pooled_bilat(k))) continue;
                const double c_fb_k = kernel_cos(fb, k);
                const double norm_k = std::sqrt(std::max(0.0, 1.0 - c_fb_k * c_fb_k));
                if (norm_k <= 1e-6) continue;
                const double cos_tau = (kernel_cos(best_c, k) - c_fb_best * c_fb_k) / (norm_best * norm_k);
                if (cos_tau < min_cos_tau) {
                  min_cos_tau = cos_tau;
                  opp_c = k;
                }
              }
            }
            const bool sym_blocked = (opp_c >= 0 && min_cos_tau < -0.5 &&
              (pooled_bilat(best_c) - pooled_bilat(opp_c)) <= 2.0 * canonical_strength * dist_fb[best_c]);
            if (!sym_blocked) {
              double mu_d = 0.0, w_d_sum = 0.0;
              for (int a = 0; a < ku; ++a) {
                const int v = nbrs[static_cast<size_t>(a)];
                const double sc_b = scores(best_c, v);
                const double sc_f = scores(fb, v);
                if (std::isfinite(sc_b) && std::isfinite(sc_f)) {
                  const double wbl = w_bilat[static_cast<size_t>(a)];
                  const double d_a = (sc_b - canonical_strength * dist_fb[best_c]) - sc_f;
                  mu_d += wbl * d_a;
                  w_d_sum += wbl;
                }
              }
              if (w_d_sum > 0.5 && (mu_d / w_d_sum) > 0.0) {
                mu_d /= w_d_sum;
                double s2_d = 0.0;
                for (int a = 0; a < ku; ++a) {
                  const int v = nbrs[static_cast<size_t>(a)];
                  const double sc_b = scores(best_c, v);
                  const double sc_f = scores(fb, v);
                  if (std::isfinite(sc_b) && std::isfinite(sc_f)) {
                    const double wbl = w_bilat[static_cast<size_t>(a)] / w_d_sum;
                    const double d_a = (sc_b - canonical_strength * dist_fb[best_c]) - sc_f;
                    const double dev = d_a - mu_d;
                    s2_d += wbl * dev * dev;
                  }
                }
                const double vq_u = (v_quad_fold.size() == V && std::isfinite(v_quad_fold[u]))
                  ? std::min(0.90, std::max(swb2, v_quad_fold[u]))
                  : std::min(0.90, std::max(swb2, 0.50));
                const double vq_score = std::min(0.88, std::max(swbl2, 0.70 * vq_u));
                const double var_mean = s2_d * (vq_score / std::max(0.12, 1.0 - vq_score));
                if (var_mean > 1e-16) {
                  const double t_diff = mu_d / std::sqrt(var_mean);
                  const double df_eff = std::max(2.5, (1.0 / vq_score) - 1.0);
                  const double p_td = std::min(0.999999, std::max(1e-15, R::pt(t_diff, df_eff, 0, 0)));
                  z_tdiff = R::qnorm5(1.0 - p_td, 0.0, 1.0, 1, 0);
                }
              }
            }
          }

          const double z_eff = std::max(z_s, z_tdiff);
          omega_u = std::min(1.0, std::max(0.0, (z_eff - z_eb_0) / z_eb_span));
          if (admissible[static_cast<size_t>(best_c)]) {
            omega_u = std::max(0.85, omega_u);
          }
        }
      }
    }
    omega_out[u] = omega_u;

    if (omega_u <= 1e-4) {
      chosen[u] = fb + 1;
      if (admissible[static_cast<size_t>(fb)]) {
        status[u] = (eff_n_out[u] > 1.05) ? "pooled_graph" : "voxel";
      }
      continue;
    }

    // Boltzmann composite weights within cosine cone around top candidate
    int k_top = fb;
    double top_obj = -std::numeric_limits<double>::infinity();
    int n_comp_ok = 0;
    std::vector<double> sc_ok_vals;
    sc_ok_vals.reserve(static_cast<size_t>(K));
    for (int k = 0; k < K; ++k) {
      if (!comp_ok[static_cast<size_t>(k)] || !std::isfinite(pooled_bilat(k))) continue;
      ++n_comp_ok;
      sc_ok_vals.push_back(pooled_bilat(k));
      const double obj_k = pooled_bilat(k) - canonical_strength * dist_fb[k];
      if (obj_k > top_obj) {
        top_obj = obj_k;
        k_top = k;
      }
    }

    if (n_comp_ok == 0 || !std::isfinite(top_obj)) continue;

    std::vector<std::pair<double, int>> cone_cands;
    cone_cands.reserve(static_cast<size_t>(n_comp_ok));
    for (int k = 0; k < K; ++k) {
      if (!comp_ok[static_cast<size_t>(k)] || !std::isfinite(pooled_bilat(k))) continue;
      if (kernel_cos(k_top, k) >= 0.85) {
        cone_cands.emplace_back(pooled_bilat(k) - canonical_strength * dist_fb[k], k);
      }
    }
    if (cone_cands.empty()) continue;

    std::sort(cone_cands.begin(), cone_cands.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                if (a.first != b.first) return a.first > b.first;
                return a.second < b.second;
              });
    const int top_k_cap = std::max(4, std::min(10, static_cast<int>(std::ceil(0.25 * n_comp_ok))));
    if (static_cast<int>(cone_cands.size()) > top_k_cap) {
      cone_cands.resize(static_cast<size_t>(top_k_cap));
    }

    std::sort(sc_ok_vals.begin(), sc_ok_vals.end());
    const size_t n_sc = sc_ok_vals.size();
    const double med_sc = (n_sc % 2 == 1)
      ? sc_ok_vals[n_sc / 2]
      : 0.5 * (sc_ok_vals[n_sc / 2 - 1] + sc_ok_vals[n_sc / 2]);
    const double spread = std::max(1e-4, sc_ok_vals.back() - med_sc);
    const double max_cone_obj = cone_cands[0].first;

    arma::vec w_data(K, arma::fill::zeros);
    double w_sub_sum = 0.0;
    for (const auto& pr : cone_cands) {
      const double delta = std::max(-20.0, (pr.first - max_cone_obj) / (0.20 * spread));
      const double ws = std::exp(delta);
      w_data(pr.second) = ws;
      w_sub_sum += ws;
    }
    if (w_sub_sum <= 0.0) continue;
    w_data /= w_sub_sum;

    arma::vec w_fin = omega_u * w_data;
    w_fin(fb) += (1.0 - omega_u);
    const double fin_sum = arma::accu(w_fin);
    if (fin_sum > 0.0) w_fin /= fin_sum;

    w_fold.col(u) = w_fin;
    int max_k = 0;
    double max_w = w_fin(0);
    for (int k = 1; k < K; ++k) {
      if (w_fin(k) > max_w) {
        max_w = w_fin(k);
        max_k = k;
      }
    }
    chosen[u] = max_k + 1;
    const double non_fb_mass = 1.0 - w_fin(fb);
    if (non_fb_mass > 0.05 || admissible[static_cast<size_t>(fb)]) {
      status[u] = (eff_n_out[u] > 1.05) ? "pooled_graph" : "voxel";
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("w_fold") = w_fold,
    Rcpp::Named("chosen") = chosen,
    Rcpp::Named("status") = status,
    Rcpp::Named("eff_n") = eff_n_out,
    Rcpp::Named("omega") = omega_out
  );
}
