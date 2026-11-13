# Fast C++-backed NIfTI <-> time-series matrix ingestion and export

resolve_nifti_spatial_geometry <- function(image, label = "image") {
  if (is.character(image)) {
    checkmate::assert_string(image, min.chars = 1L)
    checkmate::assert_file_exists(image)
    ref_file <- normalizePath(image, winslash = "/", mustWork = TRUE)
  } else {
    ref_file <- NA_character_
  }
  header <- RNifti::niftiHeader(image)
  ndim <- as.integer(header$dim[1L])
  if (is.na(ndim) || ndim < 2L || ndim > 7L) {
    stop("Invalid NIfTI dimension count in ", label, ".", call. = FALSE)
  }
  if (ndim > 4L && any(as.integer(header$dim[6L:(ndim + 1L)]) > 1L)) {
    stop(label, " must be a scalar 2D, 3D, or 4D NIfTI image (higher dimensions must be 1).", call. = FALSE)
  }
  raw_dim3 <- as.integer(header$dim[2:4])
  if (ndim == 2L) raw_dim3[3L] <- 1L
  if (any(is.na(raw_dim3) | raw_dim3 < 1L)) {
    stop("Invalid spatial dimensions in ", label, ".", call. = FALSE)
  }
  volumes <- if (ndim >= 4L) max(1L, as.integer(header$dim[5L])) else 1L

  raw_tr <- if (ndim >= 4L) as.numeric(header$pixdim[5L]) else NA_real_
  t_units <- bitwAnd(as.integer(header$xyzt_units), 56L)
  if (is.finite(raw_tr) && raw_tr > 0) {
    if (t_units == 16L) raw_tr <- raw_tr * 1e-3
    if (t_units == 24L) raw_tr <- raw_tr * 1e-6
  } else {
    raw_tr <- NA_real_
  }

  spacing <- as.numeric(header$pixdim[2:4])
  if (ndim == 2L && (!is.finite(spacing[3L]) || spacing[3L] <= 0)) spacing[3L] <- 1.0

  xform_mat <- RNifti::xform(header, useQuaternionFirst = FALSE)
  xform_code <- attr(xform_mat, "code")

  list(
    dim = raw_dim3,
    n_spatial = prod(as.numeric(raw_dim3)),
    volumes = volumes,
    header_tr = raw_tr,
    toffset = as.numeric(header$toffset),
    spacing = spacing,
    units = bitwAnd(as.integer(header$xyzt_units), 7L),
    affine = matrix(as.numeric(xform_mat), 4L, 4L),
    xform_code = if (is.null(xform_code)) NA_integer_ else as.integer(xform_code),
    reference_file = ref_file,
    nifti_header = header
  )
}

resolve_reference_from_geometry <- function(geom) {
  if (!is.list(geom)) return(NULL)
  ref_file <- geom$reference_file
  if (is.character(ref_file) && length(ref_file) == 1L && !is.na(ref_file) && file.exists(ref_file)) {
    return(ref_file)
  }
  if (!is.null(geom$nifti_header)) {
    return(geom$nifti_header)
  }
  if (is.character(ref_file) && length(ref_file) == 1L && !is.na(ref_file)) {
    return(ref_file)
  }
  NULL
}

validate_nifti_spatial_match <- function(ref, candidate, ref_label, candidate_label, tolerance = 1e-4) {
  close_num <- function(x, y) {
    length(x) == length(y) && all(is.finite(x)) && all(is.finite(y)) &&
      all(abs(x - y) <= tolerance)
  }
  compatible <- identical(ref$dim, candidate$dim) &&
    identical(ref$units, candidate$units) &&
    close_num(ref$spacing, candidate$spacing) &&
    close_num(ref$affine, candidate$affine)
  if (!isTRUE(compatible)) {
    stop(
      "Incompatible spatial NIfTI geometry (dimensions, voxel spacing, spatial units, or affine transform) ",
      "between ", ref_label, " and ", candidate_label, ". Resample/register onto a common grid first.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

read_3d_roi_values <- function(roi, label) {
  geom <- resolve_nifti_spatial_geometry(roi, label = label)
  if (geom$volumes != 1L) {
    stop(label, " must be a single-volume 3D NIfTI image (found ", geom$volumes, " volumes).", call. = FALSE)
  }
  vals <- if (is.character(roi)) {
    as.numeric(RNifti::readNifti(roi))
  } else {
    as.numeric(roi)
  }
  if (length(vals) != geom$n_spatial) {
    stop("Unexpected element count when reading ", label, ".", call. = FALSE)
  }
  list(geometry = geom, values = vals)
}

medusa_parallel_lapply <- function(X, FUN, ..., ncpus = 1L) {
  n_items <- length(X)
  if (ncpus <= 1L || n_items <= 1L) {
    return(lapply(X, FUN, ...))
  }
  n_workers <- min(as.integer(ncpus), n_items)
  can_fork <- identical(.Platform$OS.type, "unix") && !nzchar(Sys.getenv("RSTUDIO"))
  if (can_fork) {
    res <- parallel::mclapply(X, FUN, ..., mc.cores = n_workers, mc.preschedule = TRUE)
    errs <- vapply(res, inherits, logical(1L), "try-error")
    if (any(errs)) {
      first_err <- res[[which(errs)[1L]]]
      cond <- attr(first_err, "condition")
      stop(if (!is.null(cond)) conditionMessage(cond) else as.character(first_err), call. = FALSE)
    }
    return(res)
  }
  cl <- parallel::makeCluster(n_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  lib_paths <- .libPaths()
  parallel::clusterCall(cl, function(lp) .libPaths(lp), lib_paths)
  parallel::parLapply(cl, X, FUN, ...)
}

#' Read 3D/4D NIfTI runs directly into MEDuSA time-by-voxel matrices
#'
#' @description
#' Validates NIfTI spatial headers (`dim`, `pixdim`, units, and `sform`/`qform`
#' affine) across runs, optional mask, and optional parcellation atlas before
#' extracting selected voxels in C++ directly into `time x voxels` (`T x V`)
#' numeric matrices without allocating full 4D double-precision R arrays.
#'
#' The returned `"medusa_nifti_ts"` object supplies `$bold`, `$voxel_map`,
#' `$parcel_ids`, and `$tr` ready to pass directly to [estimate_crossrun_hrf()],
#' [deconvolve_crossrun_hrf()], [deconvolve_reglin()], and [write_nifti_ts()].
#'
#' @param niftis Character vector of paths to 3D or 4D scalar NIfTI files
#'   (`.nii` or `.nii.gz`). Named vectors preserve run names in `$bold`.
#' @param mask Optional path to a single-volume 3D NIfTI mask (or an in-memory
#'   `RNifti` image/array) retaining voxels strictly greater than `mask_threshold`.
#' @param atlas Optional path to a single-volume 3D NIfTI parcellation/atlas (or
#'   in-memory `RNifti` image/array) retaining voxels strictly greater than
#'   `atlas_threshold` and populating `$voxel_map$atlas_value` and `$parcel_ids`.
#' @param voxel_indices Optional numeric/integer vector of unique 1-based linear
#'   spatial indices to extract in exact order. When supplied, `mask` and `atlas`
#'   do not filter indices, though their values at `voxel_indices` are still
#'   recorded in `$voxel_map` if provided.
#' @param atlas_subset Optional numeric/integer vector of `atlas` parcel values
#'   to retain when `voxel_indices` is `NULL` and `atlas` is supplied.
#' @param mask_threshold Numeric threshold (`>` comparison) applied to `mask`
#'   when `voxel_indices` is `NULL`. Default `0`.
#' @param atlas_threshold Numeric threshold (`>` comparison) applied to `atlas`
#'   when `voxel_indices` is `NULL`. Default `0`.
#' @param tr Optional positive scalar repetition time in seconds. If `NULL`,
#'   inferred from 4D NIfTI headers (`pixdim[5]`, converting `msec`/`usec` units
#'   to seconds) when positive and consistent across runs.
#' @param run_names Optional character vector of unique run names matching
#'   `length(niftis)`. Overrides `names(niftis)`. Defaults to `"run1"`, `"run2"`,
#'   etc. when `niftis` is unnamed.
#' @param voxel_ids Optional character vector of unique voxel identifiers of
#'   length equal to the selected voxel count.
#' @param id_style Naming convention when `voxel_ids` is `NULL`: `"ijk"`
#'   (default, 1-based `"i_j_k"` grid coordinates) or `"vnum"` (`"v1"`, `"v2"`, ...).
#' @param orientation Matrix orientation: `"time_by_voxel"` (default, `T x V`
#'   matching all MEDuSA solvers) or `"voxel_by_time"` (`V x T`).
#' @param drop Logical; if `TRUE` and `length(niftis) == 1L`, returns the single
#'   extracted matrix directly with `voxel_map`, `parcel_ids`, `tr`, and
#'   `geometry` attached as attributes. Default `FALSE`.
#' @param ncpus Integer number of parallel worker processes used when extracting
#'   multiple NIfTI runs (`length(niftis) > 1L`). Default `1L`.
#'
#' @return By default (`drop = FALSE`), an S3 list of class `"medusa_nifti_ts"`
#'   with elements:
#'   \itemize{
#'     \item `bold`: Named list of numeric matrices (`T x V` by default) with
#'       voxel IDs on the voxel dimension.
#'     \item `voxel_map`: `data.frame` with columns `voxel`, `vnum`,
#'       `linear_index`, `i`, `j`, `k`, `x`, `y`, `z`, and optional `mask_value`
#'       and `atlas_value`.
#'     \item `parcel_ids`: Named vector mapping `voxel` IDs to `atlas_value` when
#'       `atlas` is supplied (`NULL` otherwise), directly usable as
#'       `spatial_groups` or `parcel_groups`.
#'     \item `tr`: Scalar repetition time in seconds (`NA_real_` if unavailable).
#'     \item `geometry`: Spatial geometry metadata list (`dim`, `spacing`,
#'       `units`, `affine`, `reference_file`).
#'     \item `orientation`: Selected matrix orientation string.
#'     \item `run_info`: `data.frame` summarizing `run`, `file`, `n_time`, and
#'       `header_tr`.
#'   }
#'
#' @seealso [write_nifti_ts()], [align_medusa_ts()], [estimate_crossrun_hrf()],
#'   [deconvolve_crossrun_hrf()], [decon_epoch_regression()]
#'
#' @examples
#' img <- array(rnorm(4 * 4 * 2 * 10), dim = c(4, 4, 2, 10))
#' mask <- array(0L, dim = c(4, 4, 2))
#' mask[1:2, 1:2, 1] <- 1L
#' f_img <- tempfile(fileext = ".nii.gz")
#' f_mask <- tempfile(fileext = ".nii.gz")
#' RNifti::writeNifti(img, f_img)
#' RNifti::writeNifti(mask, f_mask)
#'
#' ts <- read_nifti_ts(c(run1 = f_img), mask = f_mask, tr = 1.0)
#' dim(ts$bold$run1)
#' head(ts$voxel_map)
#' unlink(c(f_img, f_mask))
#'
#' @export
read_nifti_ts <- function(
  niftis,
  mask = NULL,
  atlas = NULL,
  voxel_indices = NULL,
  atlas_subset = NULL,
  mask_threshold = 0,
  atlas_threshold = 0,
  tr = NULL,
  run_names = NULL,
  voxel_ids = NULL,
  id_style = c("ijk", "vnum"),
  orientation = c("time_by_voxel", "voxel_by_time"),
  drop = FALSE,
  ncpus = 1L
) {
  checkmate::assert_character(niftis, min.len = 1L, any.missing = FALSE, min.chars = 1L)
  checkmate::assert_file_exists(niftis)
  checkmate::assert_number(mask_threshold, finite = TRUE)
  checkmate::assert_number(atlas_threshold, finite = TRUE)
  checkmate::assert_numeric(atlas_subset, any.missing = FALSE, min.len = 1L, null.ok = TRUE)
  if (!is.null(tr)) {
    checkmate::assert_number(tr, lower = 0.01, finite = TRUE)
  }
  id_style <- match.arg(id_style)
  orientation <- match.arg(orientation)
  time_first <- identical(orientation, "time_by_voxel")
  checkmate::assert_flag(drop)
  checkmate::assert_integerish(ncpus, lower = 1L, len = 1L, any.missing = FALSE)
  ncpus <- as.integer(ncpus)

  if (!is.null(run_names)) {
    checkmate::assert_character(run_names, len = length(niftis), any.missing = FALSE, min.chars = 1L, unique = TRUE)
  } else if (!is.null(names(niftis)) && all(nzchar(names(niftis))) && !anyDuplicated(names(niftis))) {
    run_names <- names(niftis)
  } else if (!is.null(names(niftis)) && any(nzchar(names(niftis)))) {
    stop("When niftis is named, all run names must be nonempty and unique.", call. = FALSE)
  } else {
    run_names <- paste0("run", seq_along(niftis))
  }

  norm_paths <- normalizePath(unname(niftis), winslash = "/", mustWork = TRUE)
  run_geoms <- lapply(seq_along(norm_paths), function(i) {
    resolve_nifti_spatial_geometry(norm_paths[i], label = paste0("niftis[", i, "] (", basename(norm_paths[i]), ")"))
  })
  ref_geom <- run_geoms[[1L]]
  ref_label <- basename(norm_paths[1L])

  if (length(run_geoms) > 1L) {
    for (i in seq.int(2L, length(run_geoms))) {
      validate_nifti_spatial_match(ref_geom, run_geoms[[i]], ref_label, basename(norm_paths[i]))
    }
  }

  mask_info <- NULL
  if (!is.null(mask)) {
    mask_label <- if (is.character(mask) && length(mask) == 1L) basename(mask) else "mask"
    mask_info <- read_3d_roi_values(mask, label = mask_label)
    validate_nifti_spatial_match(ref_geom, mask_info$geometry, ref_label, mask_label)
  }

  atlas_info <- NULL
  if (!is.null(atlas)) {
    atlas_label <- if (is.character(atlas) && length(atlas) == 1L) basename(atlas) else "atlas"
    atlas_info <- read_3d_roi_values(atlas, label = atlas_label)
    validate_nifti_spatial_match(ref_geom, atlas_info$geometry, ref_label, atlas_label)
  }

  n_spatial <- as.integer(ref_geom$n_spatial)
  if (!is.null(voxel_indices)) {
    checkmate::assert_numeric(voxel_indices, min.len = 1L, any.missing = FALSE, finite = TRUE)
    if (any(voxel_indices != floor(voxel_indices) | voxel_indices < 1 | voxel_indices > n_spatial)) {
      stop("voxel_indices must be finite, whole, one-based indices within 1..", n_spatial, ".", call. = FALSE)
    }
    voxel_indices <- as.integer(voxel_indices)
    if (anyDuplicated(voxel_indices)) {
      stop("voxel_indices must not contain duplicate indices.", call. = FALSE)
    }
  } else if (is.null(mask_info) && is.null(atlas_info)) {
    voxel_indices <- seq_len(n_spatial)
  } else {
    keep <- rep(TRUE, n_spatial)
    if (!is.null(mask_info)) {
      keep <- keep & is.finite(mask_info$values) & (mask_info$values > mask_threshold)
    }
    if (!is.null(atlas_info)) {
      keep <- keep & is.finite(atlas_info$values) & (atlas_info$values > atlas_threshold)
      if (!is.null(atlas_subset)) {
        keep <- keep & (atlas_info$values %in% atlas_subset)
      }
    }
    voxel_indices <- which(keep)
    if (length(voxel_indices) == 0L) {
      stop("No voxels remain after mask/atlas selection.", call. = FALSE)
    }
  }

  ijk <- arrayInd(voxel_indices, ref_geom$dim)
  if (ncol(ijk) == 2L) ijk <- cbind(ijk, 1L)
  xyz <- cbind(ijk - 1L, 1) %*% t(ref_geom$affine)

  if (!is.null(voxel_ids)) {
    checkmate::assert_character(
      voxel_ids, len = length(voxel_indices), any.missing = FALSE, min.chars = 1L, unique = TRUE
    )
  } else if (id_style == "ijk") {
    voxel_ids <- sprintf("%d_%d_%d", ijk[, 1L], ijk[, 2L], ijk[, 3L])
  } else {
    voxel_ids <- sprintf("v%d", seq_along(voxel_indices))
  }

  voxel_map <- data.frame(
    voxel = voxel_ids,
    vnum = seq_along(voxel_indices),
    linear_index = as.integer(voxel_indices),
    i = as.integer(ijk[, 1L]),
    j = as.integer(ijk[, 2L]),
    k = as.integer(ijk[, 3L]),
    x = as.numeric(xyz[, 1L]),
    y = as.numeric(xyz[, 2L]),
    z = as.numeric(xyz[, 3L]),
    stringsAsFactors = FALSE
  )
  if (!is.null(mask_info)) {
    voxel_map$mask_value <- mask_info$values[voxel_indices]
  }
  if (!is.null(atlas_info)) {
    voxel_map$atlas_value <- atlas_info$values[voxel_indices]
  }
  parcel_ids <- if (!is.null(atlas_info)) stats::setNames(voxel_map$atlas_value, voxel_map$voxel) else NULL

  header_trs <- vapply(run_geoms, "[[", numeric(1L), "header_tr")
  n_times <- vapply(run_geoms, "[[", integer(1L), "volumes")
  if (!is.null(tr)) {
    tr_resolved <- as.numeric(tr)
  } else {
    valid_trs <- header_trs[n_times > 1L & is.finite(header_trs) & header_trs > 0]
    if (length(valid_trs) > 0L) {
      if (max(abs(valid_trs - valid_trs[1L])) > 1e-4) {
        stop(
          "4D NIfTI headers report inconsistent repetition times across runs (",
          paste(unique(round(valid_trs, 5L)), collapse = ", "),
          " s); verify inputs or supply 'tr' explicitly.",
          call. = FALSE
        )
      }
      tr_resolved <- valid_trs[1L]
    } else {
      tr_resolved <- NA_real_
    }
  }

  bold_list <- medusa_parallel_lapply(
    norm_paths,
    function(p, idx, tf) {
      medusa_ns <- asNamespace("medusa")
      ext_fn <- get("extract_nifti_voxel_ts", envir = medusa_ns)
      ext_fn(p, idx, tf)
    },
    idx = voxel_indices,
    tf = time_first,
    ncpus = ncpus
  )

  for (i in seq_along(bold_list)) {
    if (time_first) {
      colnames(bold_list[[i]]) <- voxel_ids
    } else {
      rownames(bold_list[[i]]) <- voxel_ids
    }
  }
  names(bold_list) <- run_names

  run_info <- data.frame(
    run = run_names,
    file = norm_paths,
    n_time = as.integer(vapply(bold_list, if (time_first) nrow else ncol, integer(1L))),
    header_tr = header_trs,
    stringsAsFactors = FALSE
  )

  if (isTRUE(drop) && length(bold_list) == 1L) {
    mat <- bold_list[[1L]]
    attr(mat, "voxel_map") <- voxel_map
    attr(mat, "parcel_ids") <- parcel_ids
    attr(mat, "tr") <- tr_resolved
    attr(mat, "geometry") <- ref_geom
    attr(mat, "orientation") <- orientation
    return(mat)
  }

  structure(
    list(
      bold = bold_list,
      voxel_map = voxel_map,
      parcel_ids = parcel_ids,
      tr = tr_resolved,
      geometry = ref_geom,
      orientation = orientation,
      run_info = run_info
    ),
    class = "medusa_nifti_ts"
  )
}

#' @export
print.medusa_nifti_ts <- function(x, ...) {
  n_runs <- length(x$bold)
  n_vox <- nrow(x$voxel_map)
  tr_txt <- if (is.finite(x$tr)) sprintf("%.4g s", x$tr) else "unspecified"
  cat(sprintf(
    "<medusa_nifti_ts> %d run%s | %d voxel%s | spatial dim: %s | TR: %s | orientation: %s\n",
    n_runs, if (n_runs == 1L) "" else "s",
    n_vox, if (n_vox == 1L) "" else "s",
    paste(x$geometry$dim, collapse = "x"),
    tr_txt,
    x$orientation
  ))
  if (!is.null(x$parcel_ids)) {
    cat(sprintf("  Atlas parcels: %d unique non-zero values ($parcel_ids available)\n",
                length(unique(x$parcel_ids))))
  }
  if (!is.null(x$run_info) && nrow(x$run_info) > 0L) {
    show_n <- min(6L, nrow(x$run_info))
    for (i in seq_len(show_n)) {
      cat(sprintf("  - %s: %d volume%s (%s)\n",
                  x$run_info$run[i],
                  x$run_info$n_time[i],
                  if (x$run_info$n_time[i] == 1L) "" else "s",
                  basename(x$run_info$file[i])))
    }
    if (nrow(x$run_info) > show_n) {
      cat(sprintf("  ... and %d more run(s)\n", nrow(x$run_info) - show_n))
    }
  }
  invisible(x)
}

resolve_write_voxel_indices <- function(mat, time_first, voxel_indices, voxel_map, n_spatial) {
  n_vox <- if (time_first) ncol(mat) else nrow(mat)
  vox_names <- if (time_first) colnames(mat) else rownames(mat)

  if (!is.null(voxel_indices)) {
    checkmate::assert_numeric(voxel_indices, len = n_vox, any.missing = FALSE, finite = TRUE)
    return(as.numeric(voxel_indices))
  }

  if (!is.null(voxel_map)) {
    checkmate::assert_data_frame(voxel_map, min.rows = 1L)
    if (!"linear_index" %in% names(voxel_map)) {
      stop("voxel_map must contain a 'linear_index' column.", call. = FALSE)
    }
    if (!is.null(vox_names) && "voxel" %in% names(voxel_map)) {
      pos <- match(vox_names, as.character(voxel_map$voxel))
      if (anyNA(pos)) {
        missing_ids <- utils::head(vox_names[is.na(pos)], 5L)
        stop(
          "Matrix contains voxel IDs absent from voxel_map$voxel (e.g., ",
          paste(missing_ids, collapse = ", "), ").",
          call. = FALSE
        )
      }
      return(as.numeric(voxel_map$linear_index[pos]))
    }
    if (nrow(voxel_map) != n_vox) {
      stop(
        sprintf(
          "Matrix has %d voxels without matching names, while voxel_map has %d rows.",
          n_vox, nrow(voxel_map)
        ),
        call. = FALSE
      )
    }
    return(as.numeric(voxel_map$linear_index))
  }

  if (n_vox == n_spatial) {
    return(as.numeric(seq_len(n_spatial)))
  }

  stop(
    "Supply 'voxel_indices' or 'voxel_map' when writing a masked/subsetted matrix (",
    n_vox, " voxels vs ", n_spatial, " spatial voxels in reference).",
    call. = FALSE
  )
}

#' Write MEDuSA time-by-voxel matrices or voxel maps directly to 3D/4D NIfTI files
#'
#' @description
#' Streams a `time x voxels` (`T x V`) matrix, a `V`-length static map vector,
#' a list of run matrices (such as `deconvolve_crossrun_hrf()$activity`), or a
#' `"medusa_nifti_ts"` object directly into 3D or 4D NIfTI files (`.nii` or
#' `.nii.gz`) in C++, preserving the spatial geometry (`dim`, `pixdim`,
#' `sform`/`qform` affine) of `reference` without allocating full 4D R arrays.
#'
#' @param x Numeric matrix (`T x V` by default), numeric vector of length `V`
#'   (written as a single 3D volume), list of numeric matrices/vectors, a
#'   `"crossrun_decon"` result from [deconvolve_crossrun_hrf()], or a
#'   `"medusa_nifti_ts"` object.
#' @param file Output `.nii` or `.nii.gz` path, or a character vector of paths
#'   matching the number of runs/matrices in `x`.
#' @param reference Path to a reference 3D/4D NIfTI file, an in-memory `RNifti`
#'   image/header, or a `"medusa_nifti_ts"`, `"crossrun_decon"`, or
#'   `"crossrun_hrf"` object (which also supplies defaults for `voxel_map` and
#'   `tr`). Optional when `x` itself carries spatial geometry (`"medusa_nifti_ts"`,
#'   `"crossrun_decon"` fit from a `"medusa_nifti_ts"` input, or a matrix from
#'   `read_nifti_ts(..., drop = TRUE)`).
#' @param voxel_indices Optional numeric/integer vector of unique 1-based spatial
#'   linear indices.
#' @param voxel_map Optional `data.frame` containing `linear_index` (and optional
#'   `voxel` column IDs, matched automatically against `colnames(x)` so column
#'   subsets or reorderings map accurately to brain space).
#' @param tr Optional positive repetition time in seconds stored in `pixdim[5]`
#'   for 4D outputs. Defaults to `reference`'s TR when available, else `1.0`.
#' @param time_offset Scalar acquisition time offset stored in header `toffset`.
#'   Default `0`.
#' @param fill Numeric scalar written to non-indexed background voxels (e.g., `0`
#'   or `NA_real_`). Default `0`.
#' @param datatype On-disk NIfTI scalar type: `"float"` (32-bit IEEE float,
#'   default) or `"double"` (64-bit IEEE float).
#' @param compression Integer gzip compression level in `0..9` when writing
#'   `.nii.gz` files (`2L` by default for fast streaming deflate with compact
#'   output). Ignored for uncompressed `.nii`.
#' @param orientation Matrix orientation of `x`: `"time_by_voxel"` (default,
#'   `T x V`) or `"voxel_by_time"` (`V x T`).
#' @param ncpus Integer number of parallel workers when writing multiple files.
#'   Default `1L`.
#'
#' @return Invisibly returns a character vector of written NIfTI file paths.
#'
#' @seealso [read_nifti_ts()], [backproject_medusa()], [deconvolve_crossrun_hrf()]
#'
#' @examples
#' ref_arr <- array(0L, dim = c(4, 4, 2))
#' ref_arr[1:2, 1:2, 1] <- 1L
#' f_ref <- tempfile(fileext = ".nii.gz")
#' f_out <- tempfile(fileext = ".nii.gz")
#' RNifti::writeNifti(ref_arr, f_ref)
#'
#' idx <- which(ref_arr > 0)
#' mat <- matrix(seq_len(5 * length(idx)), nrow = 5, ncol = length(idx))
#' write_nifti_ts(mat, f_out, reference = f_ref, voxel_indices = idx, tr = 1.2)
#'
#' rt <- read_nifti_ts(f_out, mask = f_ref)
#' unlink(c(f_ref, f_out))
#'
#' @export
write_nifti_ts <- function(
  x,
  file,
  reference = NULL,
  voxel_indices = NULL,
  voxel_map = NULL,
  tr = NULL,
  time_offset = 0,
  fill = 0,
  datatype = c("float", "double"),
  compression = 2L,
  orientation = c("time_by_voxel", "voxel_by_time"),
  ncpus = 1L
) {
  datatype <- match.arg(datatype)
  orientation_supplied <- !missing(orientation)
  time_offset_supplied <- !missing(time_offset)
  orientation <- match.arg(orientation)
  checkmate::assert_integerish(compression, lower = 0L, upper = 9L, len = 1L, any.missing = FALSE)
  compression <- as.integer(compression)
  checkmate::assert_number(fill, na.ok = TRUE)
  checkmate::assert_integerish(ncpus, lower = 1L, len = 1L, any.missing = FALSE)
  ncpus <- as.integer(ncpus)

  if (inherits(x, "medusa_nifti_ts")) {
    if (is.null(voxel_map) && is.null(voxel_indices)) voxel_map <- x$voxel_map
    if (is.null(reference)) reference <- resolve_reference_from_geometry(x$geometry)
    if (is.null(tr) && is.finite(x$tr)) tr <- x$tr
    if (!orientation_supplied) orientation <- x$orientation
    x <- x$bold
  } else if (inherits(x, "crossrun_decon") ||
             (is.list(x) && !is.data.frame(x) && all(c("activity", "hrf_fit") %in% names(x)))) {
    vm_dec <- if (!is.null(x$voxel_map)) x$voxel_map else x$hrf_fit$voxel_map
    gm_dec <- if (!is.null(x$geometry)) x$geometry else x$hrf_fit$geometry
    if (is.null(voxel_map) && is.null(voxel_indices) && !is.null(vm_dec)) voxel_map <- vm_dec
    if (is.null(reference) && !is.null(gm_dec)) reference <- resolve_reference_from_geometry(gm_dec)
    if (is.null(tr) && !is.null(x$hrf_fit$settings$tr) && is.finite(x$hrf_fit$settings$tr)) {
      tr <- x$hrf_fit$settings$tr
    }
    if (!time_offset_supplied && !is.null(x$hrf_fit$settings$time_offset)) {
      time_offset <- as.numeric(x$hrf_fit$settings$time_offset)
    }
    x <- x$activity
  } else if (is.matrix(x)) {
    vm_attr <- attr(x, "voxel_map", exact = TRUE)
    if (is.null(voxel_map) && is.null(voxel_indices) && !is.null(vm_attr)) {
      voxel_map <- vm_attr
    }
    geom_attr <- attr(x, "geometry", exact = TRUE)
    if (is.null(reference) && !is.null(geom_attr)) {
      reference <- resolve_reference_from_geometry(geom_attr)
    }
    tr_attr <- attr(x, "tr", exact = TRUE)
    if (is.null(tr) && !is.null(tr_attr) && is.finite(tr_attr)) {
      tr <- tr_attr
    }
    off_attr <- attr(x, "time_offset", exact = TRUE)
    if (!time_offset_supplied && !is.null(off_attr) && all(is.finite(off_attr))) {
      time_offset <- as.numeric(off_attr)
    }
    orient_attr <- attr(x, "orientation", exact = TRUE)
    if (!orientation_supplied && !is.null(orient_attr)) {
      orientation <- orient_attr
    }
  } else if (is.list(x) && !is.data.frame(x) && length(x) > 0L && is.matrix(x[[1L]])) {
    vm_attr <- attr(x[[1L]], "voxel_map", exact = TRUE)
    if (is.null(voxel_map) && is.null(voxel_indices) && !is.null(vm_attr)) {
      voxel_map <- vm_attr
    }
    geom_attr <- attr(x[[1L]], "geometry", exact = TRUE)
    if (is.null(reference) && !is.null(geom_attr)) {
      reference <- resolve_reference_from_geometry(geom_attr)
    }
    tr_attr <- attr(x[[1L]], "tr", exact = TRUE)
    if (is.null(tr) && !is.null(tr_attr) && is.finite(tr_attr)) {
      tr <- tr_attr
    }
    if (!time_offset_supplied &&
        all(vapply(x, function(m) is.matrix(m) && !is.null(attr(m, "time_offset", exact = TRUE)), logical(1L)))) {
      time_offset <- vapply(x, function(m) as.numeric(attr(m, "time_offset", exact = TRUE))[1L], numeric(1L))
    }
  }

  if (inherits(reference, "medusa_nifti_ts")) {
    if (is.null(voxel_map) && is.null(voxel_indices)) voxel_map <- reference$voxel_map
    if (is.null(tr) && is.finite(reference$tr)) tr <- reference$tr
    reference <- resolve_reference_from_geometry(reference$geometry)
  } else if (inherits(reference, "crossrun_decon") ||
             (is.list(reference) && !is.data.frame(reference) && all(c("activity", "hrf_fit") %in% names(reference)))) {
    vm_ref <- if (!is.null(reference$voxel_map)) reference$voxel_map else reference$hrf_fit$voxel_map
    gm_ref <- if (!is.null(reference$geometry)) reference$geometry else reference$hrf_fit$geometry
    if (is.null(voxel_map) && is.null(voxel_indices) && !is.null(vm_ref)) voxel_map <- vm_ref
    if (is.null(tr) && !is.null(reference$hrf_fit$settings$tr) && is.finite(reference$hrf_fit$settings$tr)) {
      tr <- reference$hrf_fit$settings$tr
    }
    reference <- resolve_reference_from_geometry(gm_ref)
  } else if (inherits(reference, "crossrun_hrf")) {
    if (is.null(voxel_map) && is.null(voxel_indices) && !is.null(reference$voxel_map)) {
      voxel_map <- reference$voxel_map
    }
    if (is.null(tr) && !is.null(reference$settings$tr) && is.finite(reference$settings$tr)) {
      tr <- reference$settings$tr
    }
    reference <- resolve_reference_from_geometry(reference$geometry)
  }

  if (is.null(reference) || (is.character(reference) && length(reference) == 1L && is.na(reference))) {
    stop("A valid 'reference' NIfTI file path or NIfTI object is required.", call. = FALSE)
  }

  ref_geom <- resolve_nifti_spatial_geometry(reference, label = "reference")
  n_spatial <- as.integer(ref_geom$n_spatial)
  time_first <- identical(orientation, "time_by_voxel")

  mat_list <- if (is.list(x) && !is.data.frame(x)) {
    x
  } else {
    list(x)
  }
  if (length(mat_list) == 0L) {
    stop("x must contain at least one matrix or vector to write.", call. = FALSE)
  }

  checkmate::assert_character(file, len = length(mat_list), any.missing = FALSE, min.chars = 1L)

  if (!is.null(tr)) {
    checkmate::assert_number(tr, lower = 0.01, finite = TRUE)
    tr_val <- as.numeric(tr)
  } else if (is.finite(ref_geom$header_tr) && ref_geom$header_tr > 0) {
    tr_val <- ref_geom$header_tr
  } else {
    tr_val <- 1.0
  }

  checkmate::assert_numeric(time_offset, any.missing = FALSE, finite = TRUE)
  if (length(time_offset) == 1L) {
    offsets <- rep(as.numeric(time_offset), length(mat_list))
  } else if (length(time_offset) == length(mat_list)) {
    offsets <- as.numeric(time_offset)
  } else {
    stop("time_offset must be length 1 or match the number of output files.", call. = FALSE)
  }

  prepared <- vector("list", length(mat_list))
  for (i in seq_along(mat_list)) {
    item <- mat_list[[i]]
    if (is.atomic(item) && is.numeric(item) && is.null(dim(item))) {
      nm <- names(item)
      if (time_first) {
        item <- matrix(as.numeric(item), nrow = 1L)
        if (!is.null(nm)) colnames(item) <- nm
      } else {
        item <- matrix(as.numeric(item), ncol = 1L)
        if (!is.null(nm)) rownames(item) <- nm
      }
    } else {
      checkmate::assert_matrix(item, mode = "numeric", min.rows = 1L, min.cols = 1L)
      if (storage.mode(item) != "double") storage.mode(item) <- "double"
    }

    idx_i <- resolve_write_voxel_indices(item, time_first, voxel_indices, voxel_map, n_spatial)
    out_dir <- dirname(file[i])
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    prepared[[i]] <- list(
      file = file[i],
      mat = item,
      indices = idx_i,
      offset = offsets[i]
    )
  }

  ref_arg <- if (is.character(reference)) ref_geom$reference_file else RNifti::niftiHeader(reference)

  written_list <- medusa_parallel_lapply(
    prepared,
    function(task, ref_in, tr_in, fill_in, dt_in, comp_in, tf_in) {
      medusa_ns <- asNamespace("medusa")
      write_fn <- get("write_nifti_indexed_cpp", envir = medusa_ns)
      write_fn(ref_in, task$file, task$mat, task$indices, tr_in, task$offset, fill_in, dt_in, comp_in, tf_in)
    },
    ref_in = ref_arg,
    tr_in = tr_val,
    fill_in = as.numeric(fill),
    dt_in = datatype,
    comp_in = compression,
    tf_in = time_first,
    ncpus = ncpus
  )
  written <- unlist(written_list, use.names = FALSE)

  invisible(written)
}
