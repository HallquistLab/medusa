#' Build a sparse spatial neighborhood graph for bilateral HRF and noise pooling
#'
#' Constructs a compressed-sparse-row (CSR) spatial neighborhood graph
#' (\code{"medusa_spatial_graph"}) for \code{\link{estimate_crossrun_hrf}} when
#' \code{pooling = "graph"}.
#'
#' @param x Spatial input: a \code{"medusa_nifti_ts"} object, a list containing
#'   \code{voxel_map} (and optional \code{geometry}), a \code{voxel_map}
#'   \code{data.frame} with lattice columns \code{i, j, k} (or world coordinates
#'   \code{x, y, z}), a numeric \eqn{V \times D} coordinate matrix, an edge
#'   \code{data.frame} with columns \code{from}, \code{to} (and optional positive
#'   \code{dist} together with \code{n_voxels}), a square numeric/sparse distance
#'   matrix (nonzero off-diagonal entries giving neighbor distances), or an
#'   existing \code{"medusa_spatial_graph"}.
#' @param spatial_fwhm Optional positive numeric Gaussian spatial kernel FWHM (in
#'   mm or coordinate units). Defaults to \code{2.5 * median(spacing)} for lattice
#'   inputs and \code{2.5} times the median nearest-neighbor distance for
#'   coordinate matrices.
#' @param max_neighbors Maximum number of non-self spatial neighbors retained per
#'   node (default \code{32L}). Every node always retains itself at index 1 of its
#'   CSR row with distance \code{0}.
#' @param max_radius Optional positive search radius cutoff. Defaults to
#'   \code{max(2.75 * sigma, 1.80 * step_scale)} where
#'   \code{sigma = spatial_fwhm / sqrt(8 * log(2))}.
#' @param spatial_groups Optional character or integer vector of parcel/group
#'   labels of length \eqn{V}. Only restricts graph construction when
#'   \code{parcel_boundary = "hard"}; under \code{"soft"} or \code{"ignore"},
#'   cross-border neighbors remain in the graph so bilateral pooling can modulate
#'   cross-border weights cleanly.
#' @param parcel_boundary Boundary handling when \code{spatial_groups} is present:
#'   \code{"soft"} (default), \code{"hard"} (omit all cross-parcel edges and
#'   geodesic paths), or \code{"ignore"}.
#' @param metric Distance metric for volumetric/coordinate inputs: \code{"auto"}
#'   (uses 26-connected in-mask lattice geodesic distance whenever unique integer
#'   lattice indices \code{i, j, k} are available, so out-of-mask CSF/sulcal gaps
#'   act as natural barriers, and falls back to 3D Euclidean radius/\eqn{k}-NN
#'   otherwise), \code{"geodesic"}, or \code{"euclidean"}.
#' @param n_voxels Optional total node count \eqn{V} when \code{x} is an edge
#'   \code{data.frame} that may contain isolated nodes.
#' @param voxel_ids Optional character vector of node IDs of length \eqn{V}.
#'
#' @return An S3 object of class \code{"medusa_spatial_graph"} containing 1-based
#'   CSR arrays \code{row_ptr} (length \eqn{V + 1}), \code{col_idx}, and
#'   \code{dist}, along with \code{n_voxels}, \code{voxel_ids},
#'   \code{spatial_fwhm}, \code{sigma}, \code{max_radius}, \code{max_neighbors},
#'   \code{metric}, and \code{parcel_boundary}.
#' @export
build_spatial_graph <- function(
    x,
    spatial_fwhm = NULL,
    max_neighbors = 32L,
    max_radius = NULL,
    spatial_groups = NULL,
    parcel_boundary = c("soft", "hard", "ignore"),
    metric = c("auto", "geodesic", "euclidean"),
    n_voxels = NULL,
    voxel_ids = NULL) {
  parcel_boundary <- match.arg(parcel_boundary)
  metric <- match.arg(metric)
  checkmate::assert_integerish(max_neighbors, lower = 1L, len = 1L, any.missing = FALSE)
  max_neighbors <- as.integer(max_neighbors)
  if (!is.null(spatial_fwhm)) {
    checkmate::assert_number(spatial_fwhm, lower = .Machine$double.eps, finite = TRUE)
  }
  if (!is.null(max_radius)) {
    checkmate::assert_number(max_radius, lower = .Machine$double.eps, finite = TRUE)
  }

  fwhm_to_sigma <- function(fwhm) max(as.numeric(fwhm) / sqrt(8.0 * log(2.0)), 1e-6)

  # 1. Existing medusa_spatial_graph
  if (inherits(x, "medusa_spatial_graph")) {
    V <- as.integer(x$n_voxels)
    if (!is.null(n_voxels) && as.integer(n_voxels) != V) {
      stop("spatial_graph node count (", V, ") does not match expected voxel count (", n_voxels, ").", call. = FALSE)
    }
    out <- x
    if (!is.null(spatial_fwhm)) {
      out$spatial_fwhm <- as.numeric(spatial_fwhm)
      out$sigma <- fwhm_to_sigma(spatial_fwhm)
    }
    if (!is.null(voxel_ids) && length(voxel_ids) == V) {
      out$voxel_ids <- as.character(voxel_ids)
    }
    if (parcel_boundary == "hard" && !is.null(spatial_groups) && !identical(x$parcel_boundary, "hard")) {
      g_ids <- hrf_encode_group_ids(spatial_groups, V)
      out <- hrf_filter_csr_hard_groups(out, g_ids)
    }
    out$parcel_boundary <- parcel_boundary
    return(out)
  }

  # Unpack medusa_nifti_ts, attributed matrix, or list(voxel_map, geometry)
  vm <- NULL
  geom <- NULL
  if (inherits(x, "medusa_nifti_ts")) {
    unp <- hrf_unpack_bold_input(x)
    vm <- unp$voxel_map
    geom <- unp$geometry
    if (is.null(voxel_ids) && length(unp$bold) >= 1L && !is.null(colnames(unp$bold[[1L]]))) {
      voxel_ids <- colnames(unp$bold[[1L]])
    }
    if (is.null(spatial_groups) && !is.null(unp$parcel_ids)) {
      spatial_groups <- unp$parcel_ids
    }
  } else if (is.matrix(x) && !is.null(attr(x, "voxel_map", exact = TRUE))) {
    vm <- attr(x, "voxel_map", exact = TRUE)
    geom <- attr(x, "geometry", exact = TRUE)
  } else if (is.list(x) && !is.data.frame(x) && !is.null(x$voxel_map)) {
    vm <- x$voxel_map
    geom <- x$geometry
  } else if (is.data.frame(x) && all(c("i", "j", "k") %in% names(x)) && !all(c("from", "to") %in% names(x))) {
    vm <- x
  }

  # 2. Volumetric lattice / voxel_map input
  if (!is.null(vm) && is.data.frame(vm)) {
    V <- nrow(vm)
    if (is.null(voxel_ids) && "voxel" %in% names(vm)) {
      voxel_ids <- as.character(vm$voxel)
    }
    g_ids <- hrf_encode_group_ids(spatial_groups, V)
    hard_bnd <- (parcel_boundary == "hard") && (length(g_ids) == V)

    has_ijk <- all(c("i", "j", "k") %in% names(vm)) &&
      all(is.finite(as.numeric(vm$i))) &&
      all(is.finite(as.numeric(vm$j))) &&
      all(is.finite(as.numeric(vm$k)))
    ijk_unique <- FALSE
    if (has_ijk) {
      ijk_mat <- cbind(as.integer(round(vm$i)), as.integer(round(vm$j)), as.integer(round(vm$k)))
      ijk_keys <- paste(ijk_mat[, 1L], ijk_mat[, 2L], ijk_mat[, 3L], sep = ":")
      ijk_unique <- !anyDuplicated(ijk_keys)
    }

    if (metric %in% c("auto", "geodesic") && has_ijk && ijk_unique) {
      spacing <- hrf_infer_lattice_spacing(vm, geom)
      step_med <- stats::median(spacing)
      eff_fwhm <- if (is.null(spatial_fwhm)) 2.5 * step_med else as.numeric(spatial_fwhm)
      sigma <- fwhm_to_sigma(eff_fwhm)
      eff_rad <- if (is.null(max_radius)) max(2.75 * sigma, 1.80 * max(spacing)) else as.numeric(max_radius)
      csr <- build_lattice_geodesic_graph_cpp(
        ijk = ijk_mat,
        spacing = as.numeric(spacing),
        max_radius = eff_rad,
        max_neighbors = max_neighbors,
        group_ids = g_ids,
        hard_boundary = hard_bnd
      )
      return(hrf_new_spatial_graph(
        csr = csr,
        n_voxels = V,
        voxel_ids = voxel_ids,
        spatial_fwhm = eff_fwhm,
        sigma = sigma,
        max_radius = eff_rad,
        max_neighbors = max_neighbors,
        metric = "lattice_geodesic",
        parcel_boundary = parcel_boundary
      ))
    }

    if (metric == "geodesic") {
      stop("metric = 'geodesic' requires unique finite integer lattice columns 'i', 'j', 'k'.", call. = FALSE)
    }
    if (all(c("x", "y", "z") %in% names(vm))) {
      x <- as.matrix(vm[, c("x", "y", "z"), drop = FALSE])
    } else if (has_ijk) {
      x <- matrix(as.numeric(ijk_mat), ncol = 3L)
    } else {
      stop("voxel_map must contain either ('i', 'j', 'k') or ('x', 'y', 'z') columns.", call. = FALSE)
    }
  }

  # 3. Edge list data.frame (from, to, optional dist)
  if (is.data.frame(x) && all(c("from", "to") %in% names(x))) {
    return(hrf_build_custom_edge_graph(
      edges = x,
      n_voxels = n_voxels,
      voxel_ids = voxel_ids,
      spatial_fwhm = spatial_fwhm,
      max_neighbors = max_neighbors,
      max_radius = max_radius,
      spatial_groups = spatial_groups,
      parcel_boundary = parcel_boundary
    ))
  }

  # 4. Square sparse Matrix (e.g., Matrix::dgCMatrix)
  if (inherits(x, "Matrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package 'Matrix' is required to convert sparse Matrix spatial_graph inputs.", call. = FALSE)
    }
    if (nrow(x) != ncol(x)) {
      stop("Sparse Matrix spatial_graph must be square (V x V).", call. = FALSE)
    }
    sm <- Matrix::summary(x)
    ed <- data.frame(from = as.integer(sm$i), to = as.integer(sm$j), dist = as.numeric(sm$x))
    ed <- ed[is.finite(ed$dist) & ed$dist > 0 & ed$from != ed$to, , drop = FALSE]
    if (is.null(voxel_ids) && !is.null(rownames(x))) voxel_ids <- rownames(x)
    return(hrf_build_custom_edge_graph(
      edges = ed,
      n_voxels = nrow(x),
      voxel_ids = voxel_ids,
      spatial_fwhm = spatial_fwhm,
      max_neighbors = max_neighbors,
      max_radius = max_radius,
      spatial_groups = spatial_groups,
      parcel_boundary = parcel_boundary
    ))
  }

  # 5. Coordinate matrix (or data.frame of coordinates)
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.atomic(x) && is.null(dim(x)) && is.numeric(x)) {
    x <- matrix(as.numeric(x), ncol = 1L, dimnames = list(names(x), NULL))
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(
      "Unsupported input 'x' for build_spatial_graph(): supply a 'medusa_nifti_ts', ",
      "voxel_map, numeric coordinate matrix, or edge data.frame.",
      call. = FALSE
    )
  }
  if (metric == "geodesic") {
    stop("metric = 'geodesic' requires a lattice voxel_map with columns 'i', 'j', 'k'.", call. = FALSE)
  }
  if (any(!is.finite(x))) stop("Spatial coordinates must be finite.", call. = FALSE)

  V <- nrow(x)
  if (is.null(voxel_ids) && !is.null(rownames(x))) voxel_ids <- rownames(x)
  g_ids <- hrf_encode_group_ids(spatial_groups, V)
  hard_bnd <- (parcel_boundary == "hard") && (length(g_ids) == V)

  med_nn <- estimate_median_nn_dist_cpp(x)
  eff_fwhm <- if (is.null(spatial_fwhm)) 2.5 * med_nn else as.numeric(spatial_fwhm)
  sigma <- fwhm_to_sigma(eff_fwhm)
  eff_rad <- if (is.null(max_radius)) max(2.75 * sigma, 1.80 * med_nn) else as.numeric(max_radius)

  csr <- build_coords_knn_graph_cpp(
    coords = x,
    max_radius = eff_rad,
    max_neighbors = max_neighbors,
    group_ids = g_ids,
    hard_boundary = hard_bnd
  )
  hrf_new_spatial_graph(
    csr = csr,
    n_voxels = V,
    voxel_ids = voxel_ids,
    spatial_fwhm = eff_fwhm,
    sigma = sigma,
    max_radius = eff_rad,
    max_neighbors = max_neighbors,
    metric = "euclidean",
    parcel_boundary = parcel_boundary
  )
}

#' @export
print.medusa_spatial_graph <- function(x, ...) {
  deg <- diff(x$row_ptr) - 1L
  med_deg <- if (length(deg)) stats::median(deg) else 0
  cat(sprintf(
    "<medusa_spatial_graph> %d node%s | metric: %s | FWHM: %.3g (sigma = %.3g) | median neighbors: %.1f (max %d) | boundary: %s\n",
    x$n_voxels,
    if (x$n_voxels == 1L) "" else "s",
    x$metric,
    x$spatial_fwhm,
    x$sigma,
    med_deg,
    x$max_neighbors,
    x$parcel_boundary
  ))
  invisible(x)
}

hrf_new_spatial_graph <- function(csr, n_voxels, voxel_ids, spatial_fwhm, sigma,
    max_radius, max_neighbors, metric, parcel_boundary) {
  structure(
    list(
      n_voxels = as.integer(n_voxels),
      voxel_ids = if (is.null(voxel_ids)) NULL else as.character(voxel_ids),
      row_ptr = as.integer(csr$row_ptr),
      col_idx = as.integer(csr$col_idx),
      dist = as.numeric(csr$dist),
      spatial_fwhm = as.numeric(spatial_fwhm),
      sigma = as.numeric(sigma),
      max_radius = as.numeric(max_radius),
      max_neighbors = as.integer(max_neighbors),
      metric = as.character(metric),
      parcel_boundary = as.character(parcel_boundary)
    ),
    class = "medusa_spatial_graph"
  )
}

hrf_encode_group_ids <- function(spatial_groups, n_voxels) {
  if (is.null(spatial_groups)) return(integer(0L))
  if (length(spatial_groups) != n_voxels) {
    stop("spatial_groups length (", length(spatial_groups), ") must match voxel count (", n_voxels, ").", call. = FALSE)
  }
  g_chr <- as.character(spatial_groups)
  if (anyNA(g_chr) || any(!nzchar(g_chr))) {
    stop("spatial_groups must contain nonmissing, nonempty labels.", call. = FALSE)
  }
  as.integer(match(g_chr, unique(g_chr)))
}

hrf_infer_lattice_spacing <- function(vm, geom = NULL) {
  if (!is.null(geom) && !is.null(geom$spacing) && length(geom$spacing) >= 3L) {
    sp <- abs(as.numeric(geom$spacing[1:3]))
    if (all(is.finite(sp)) && all(sp > 0)) return(sp)
  }
  if (!is.null(geom) && !is.null(geom$affine) && is.matrix(geom$affine) &&
      nrow(geom$affine) >= 3L && ncol(geom$affine) >= 3L) {
    sp <- sqrt(colSums(geom$affine[1:3, 1:3, drop = FALSE]^2))
    if (all(is.finite(sp)) && all(sp > 0)) return(as.numeric(sp))
  }
  if (all(c("x", "y", "z") %in% names(vm)) && nrow(vm) >= 2L) {
    sp <- c(1.0, 1.0, 1.0)
    cols_ijk <- c("i", "j", "k")
    xyz <- as.matrix(vm[, c("x", "y", "z"), drop = FALSE])
    for (ax in 1:3) {
      d_ijk <- abs(diff(vm[[cols_ijk[ax]]]))
      other1 <- abs(diff(vm[[cols_ijk[(ax %% 3L) + 1L]]]))
      other2 <- abs(diff(vm[[cols_ijk[((ax + 1L) %% 3L) + 1L]]]))
      hit <- which(d_ijk == 1L & other1 == 0L & other2 == 0L)
      if (length(hit)) {
        step_d <- sqrt(rowSums((xyz[hit + 1L, , drop = FALSE] - xyz[hit, , drop = FALSE])^2))
        step_d <- step_d[is.finite(step_d) & step_d > 0]
        if (length(step_d)) sp[ax] <- stats::median(step_d)
      }
    }
    return(sp)
  }
  c(1.0, 1.0, 1.0)
}

hrf_filter_csr_hard_groups <- function(graph, group_ids) {
  V <- graph$n_voxels
  if (length(group_ids) != V) return(graph)
  rp_in <- graph$row_ptr
  ci_in <- graph$col_idx
  ds_in <- graph$dist
  rp_out <- integer(V + 1L)
  rp_out[1L] <- 1L
  ci_chunks <- vector("list", V)
  ds_chunks <- vector("list", V)
  tot <- 0L
  for (u in seq_len(V)) {
    p_s <- rp_in[u]
    p_e <- rp_in[u + 1L] - 1L
    if (p_e < p_s) {
      ci_chunks[[u]] <- u
      ds_chunks[[u]] <- 0.0
      tot <- tot + 1L
    } else {
      nbrs <- ci_in[p_s:p_e]
      dists <- ds_in[p_s:p_e]
      keep <- (nbrs == u) | (group_ids[nbrs] == group_ids[u])
      ci_u <- nbrs[keep]
      ds_u <- dists[keep]
      if (!length(ci_u) || ci_u[1L] != u) {
        ci_u <- c(u, setdiff(ci_u, u))
        ds_u <- c(0.0, ds_u[nbrs[keep] != u])
      }
      ci_chunks[[u]] <- ci_u
      ds_chunks[[u]] <- ds_u
      tot <- tot + length(ci_u)
    }
    rp_out[u + 1L] <- tot + 1L
  }
  graph$row_ptr <- rp_out
  graph$col_idx <- unlist(ci_chunks, use.names = FALSE)
  graph$dist <- unlist(ds_chunks, use.names = FALSE)
  graph
}

hrf_build_custom_edge_graph <- function(edges, n_voxels, voxel_ids, spatial_fwhm,
    max_neighbors, max_radius, spatial_groups, parcel_boundary) {
  from <- edges$from
  to <- edges$to
  if (is.character(from) || is.character(to)) {
    if (is.null(voxel_ids)) {
      voxel_ids <- unique(c(as.character(from), as.character(to)))
    }
    from <- match(as.character(from), voxel_ids)
    to <- match(as.character(to), voxel_ids)
    if (anyNA(from) || anyNA(to)) {
      stop("Edge list contains node IDs not present in voxel_ids.", call. = FALSE)
    }
  } else {
    from <- as.integer(from)
    to <- as.integer(to)
  }
  V <- if (!is.null(n_voxels)) {
    as.integer(n_voxels)
  } else if (!is.null(voxel_ids)) {
    length(voxel_ids)
  } else if (length(from) || length(to)) {
    max(c(from, to), na.rm = TRUE)
  } else {
    0L
  }
  if (any(from < 1L | from > V | to < 1L | to > V)) {
    stop("Edge indices in 'from'/'to' must lie within 1..n_voxels.", call. = FALSE)
  }
  d_vec <- if ("dist" %in% names(edges)) as.numeric(edges$dist) else rep(1.0, nrow(edges))
  if (any(!is.finite(d_vec) | d_vec < 0)) {
    stop("Edge 'dist' values must be finite and nonnegative.", call. = FALSE)
  }
  keep_e <- (from != to) & (d_vec > 0)
  from <- from[keep_e]
  to <- to[keep_e]
  d_vec <- d_vec[keep_e]

  # Symmetrize undirected edges and keep minimum distance on duplicate pairs
  all_from <- c(from, to)
  all_to <- c(to, from)
  all_d <- c(d_vec, d_vec)

  g_ids <- hrf_encode_group_ids(spatial_groups, V)
  if (parcel_boundary == "hard" && length(g_ids) == V && length(all_from)) {
    same_g <- g_ids[all_from] == g_ids[all_to]
    all_from <- all_from[same_g]
    all_to <- all_to[same_g]
    all_d <- all_d[same_g]
  }

  pos_d <- all_d[all_d > 0]
  med_d <- if (length(pos_d)) stats::median(pos_d) else 1.0
  eff_fwhm <- if (is.null(spatial_fwhm)) 2.5 * med_d else as.numeric(spatial_fwhm)
  sigma <- max(eff_fwhm / sqrt(8.0 * log(2.0)), 1e-6)
  eff_rad <- if (is.null(max_radius)) max(2.75 * sigma, 1.80 * med_d) else as.numeric(max_radius)

  within_r <- all_d <= (eff_rad + 1e-9)
  all_from <- all_from[within_r]
  all_to <- all_to[within_r]
  all_d <- all_d[within_r]

  split_idx <- if (length(all_from)) split(seq_along(all_from), all_from) else list()
  rp <- integer(V + 1L)
  rp[1L] <- 1L
  ci_list <- vector("list", V)
  ds_list <- vector("list", V)
  tot <- 0L
  for (u in seq_len(V)) {
    e_u <- split_idx[[as.character(u)]]
    if (!length(e_u)) {
      ci_list[[u]] <- u
      ds_list[[u]] <- 0.0
      tot <- tot + 1L
    } else {
      nb <- all_to[e_u]
      dd <- all_d[e_u]
      ord <- order(round(dd, 8L), nb)
      nb <- nb[ord]
      dd <- dd[ord]
      uniq <- !duplicated(nb)
      nb <- nb[uniq]
      dd <- dd[uniq]
      if (length(nb) > max_neighbors) {
        nb <- nb[seq_len(max_neighbors)]
        dd <- dd[seq_len(max_neighbors)]
      }
      ci_list[[u]] <- c(u, nb)
      ds_list[[u]] <- c(0.0, dd)
      tot <- tot + 1L + length(nb)
    }
    rp[u + 1L] <- tot + 1L
  }

  hrf_new_spatial_graph(
    csr = list(
      row_ptr = rp,
      col_idx = unlist(ci_list, use.names = FALSE),
      dist = unlist(ds_list, use.names = FALSE)
    ),
    n_voxels = V,
    voxel_ids = voxel_ids,
    spatial_fwhm = eff_fwhm,
    sigma = sigma,
    max_radius = eff_rad,
    max_neighbors = max_neighbors,
    metric = "custom",
    parcel_boundary = parcel_boundary
  )
}

# Fast vectorized column signature without per-column R apply() overhead.
hrf_fast_col_keys <- function(mat, digits = 4L) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  if (nc == 0L) return(character(0L))
  if (nr == 0L) return(rep("", nc))
  rmat <- round(mat, digits)
  row_lists <- lapply(seq_len(nr), function(i) rmat[i, ])
  do.call(paste, c(row_lists, list(sep = ",")))
}

# Deterministic HRF waveform prototype quantizer for scalable operator caching in
# whole-brain soft_average mode. Preserves exact fallback waveforms as prototype 1,
# keeps exact unique waveforms whenever total unique shapes <= max_strata, and
# otherwise runs deterministic maximin seeding + count-weighted Lloyd updates.
hrf_quantize_prototypes <- function(k_mat, fallback_vec, kernel_normalize = "signed_sum",
    max_strata = 48L, rel_tol = 0.0075) {
  n_lags <- nrow(k_mat)
  n_vox <- ncol(k_mat)
  max_strata <- max(2L, as.integer(max_strata))
  if (n_vox == 0L) {
    return(list(kernels = k_mat, stratum_id = integer(0L), n_strata = 0L))
  }

  fb_diff <- colSums(abs(sweep(k_mat, 1L, fallback_vec, "-")))
  is_fb <- is.finite(fb_diff) & (fb_diff <= 1e-7)
  stratum_id <- integer(n_vox)
  stratum_id[is_fb] <- 1L
  k_out <- k_mat
  if (any(is_fb)) {
    k_out[, is_fb] <- fallback_vec
  }

  rem_idx <- which(!is_fb)
  if (!length(rem_idx)) {
    return(list(kernels = k_out, stratum_id = stratum_id, n_strata = 1L))
  }

  sigs <- hrf_fast_col_keys(k_mat[, rem_idx, drop = FALSE], digits = 4L)
  u_sigs <- unique(sigs)
  grp_map <- match(sigs, u_sigs)
  n_u <- length(u_sigs)
  first_pos <- rem_idx[match(seq_len(n_u), grp_map)]
  U <- k_mat[, first_pos, drop = FALSE]

  n_non_fb_cap <- max_strata - 1L
  if (n_u <= n_non_fb_cap) {
    offset <- if (any(is_fb)) 1L else 0L
    stratum_id[rem_idx] <- grp_map + offset
    for (g in seq_len(n_u)) {
      k_out[, rem_idx[grp_map == g]] <- U[, g]
    }
    if (!any(is_fb)) {
      return(list(kernels = k_out, stratum_id = stratum_id, n_strata = n_u))
    }
    return(list(kernels = k_out, stratum_id = stratum_id, n_strata = n_u + 1L))
  }

  # Deterministic count-weighted maximin + Lloyd prototype quantization
  cnts <- as.numeric(tabulate(grp_map, nbins = n_u))
  u_norms <- pmax(sqrt(colSums(U^2)), 1e-8)

  # Seed 1: highest-frequency non-fallback shape (ties broken by first index)
  seed_1 <- which.max(cnts)
  proto_idx <- seed_1
  d2_min <- colSums((U - U[, seed_1])^2)

  while (length(proto_idx) < n_non_fb_cap) {
    rel_err <- sqrt(pmax(d2_min, 0.0)) / u_norms
    if (max(rel_err) <= rel_tol) break
    nxt <- which.max(rel_err)
    proto_idx <- c(proto_idx, nxt)
    d2_new <- colSums((U - U[, nxt])^2)
    d2_min <- pmin(d2_min, d2_new)
  }

  C <- U[, proto_idx, drop = FALSE]
  K_p <- ncol(C)
  assign_u <- rep(1L, n_u)

  for (iter in seq_len(10L)) {
    # Exact Euclidean distance matrix between U (L x n_u) and C (L x K_p)
    u2 <- colSums(U^2)
    c2 <- colSums(C^2)
    d2_mat <- outer(u2, c2, "+") - 2.0 * crossprod(U, C)
    new_assign <- max.col(-d2_mat, ties.method = "first")
    if (iter > 1L && identical(new_assign, assign_u)) break
    assign_u <- new_assign
    for (p in seq_len(K_p)) {
      members <- which(assign_u == p)
      if (length(members)) {
        w_m <- cnts[members] / sum(cnts[members])
        c_raw <- as.numeric(U[, members, drop = FALSE] %*% w_m)
        C[, p] <- normalize_decon_kernel(c_raw, mode = kernel_normalize)
      }
    }
  }

  offset <- if (any(is_fb)) 1L else 0L
  vox_proto <- assign_u[grp_map]
  stratum_id[rem_idx] <- vox_proto + offset
  k_out[, rem_idx] <- C[, vox_proto, drop = FALSE]

  list(
    kernels = k_out,
    stratum_id = as.integer(stratum_id),
    n_strata = as.integer(K_p + offset)
  )
}
