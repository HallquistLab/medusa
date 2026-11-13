#' Backproject parcel-level MEDuSA coefficients onto a 4D NIfTI brain mask
#'
#' @description
#' Maps parcel-by-time model statistics (for example from [mixed_by()]) back
#' into 4D NIfTI volumes (3D space across event-locked timepoints) using an
#' integer-valued 3D parcellation mask.
#'
#' @param coef_df A `data.frame` containing parcel, time, optional effect, and
#'   statistic columns to write.
#' @param brain_mask Path to a 3D integer-valued NIfTI parcellation mask.
#' @param plot_cols Character vector of numeric column names in `coef_df` to
#'   backproject.
#' @param parcel_col Column name in `coef_df` identifying integer mask parcel
#'   values. Default `"roi"`.
#' @param time_col Column name in `coef_df` identifying event-locked timepoints.
#'   Default `"time"`.
#' @param effect_col Column name in `coef_df` identifying model terms/effects,
#'   or `NULL` if `coef_df` contains a single effect. Default `"term"`.
#' @param output_dir Output directory where `.nii.gz` files will be written.
#'   Default `getwd()`.
#' @param ncpus Integer number of parallel worker processes. Default `4L`.
#'
#' @return Invisibly returns a character vector of written `.nii.gz` file paths.
#' @importFrom checkmate assert_data_frame assert_string assert_file_exists assert_character assert_subset assert_integerish
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom iterators iter
#' @importFrom foreach %dopar% foreach registerDoSEQ
#' @importFrom stats median
#' @export
backproject_medusa <- function(
  coef_df,
  brain_mask,
  plot_cols = NULL,
  parcel_col = "roi",
  time_col = "time",
  effect_col = "term",
  output_dir = getwd(),
  ncpus = 4L
) {
  checkmate::assert_data_frame(coef_df)
  checkmate::assert_string(brain_mask)
  checkmate::assert_file_exists(brain_mask)
  checkmate::assert_character(plot_cols, min.len = 1L, any.missing = FALSE)
  checkmate::assert_subset(plot_cols, names(coef_df))
  checkmate::assert_string(parcel_col)
  checkmate::assert_string(time_col)
  checkmate::assert_string(effect_col, null.ok = TRUE)
  checkmate::assert_subset(c(parcel_col, time_col), names(coef_df))
  if (!is.null(effect_col)) {
    checkmate::assert_subset(effect_col, names(coef_df))
  }
  checkmate::assert_string(output_dir)
  checkmate::assert_integerish(ncpus, lower = 1L, len = 1L)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  mask_hdr <- RNifti::niftiHeader(brain_mask)
  n_dim <- as.integer(mask_hdr$dim[1L])
  if (n_dim != 3L && !(n_dim == 4L && mask_hdr$dim[5L] <= 1L)) {
    stop("brain_mask does not appear to be a 3-D mask image.", call. = FALSE)
  }
  mask_img <- RNifti::readNifti(brain_mask)
  dim3 <- as.integer(mask_hdr$dim[2:4])
  mask_arr <- array(as.numeric(mask_img), dim = dim3)
  n_vox3 <- prod(dim3)

  mvals <- sort(unique(coef_df[[parcel_col]]))
  imgvals <- sort(unique(as.vector(mask_arr)))
  imgvals <- imgvals[is.finite(imgvals) & imgvals != 0]
  checkmate::assert_subset(mvals, imgvals)

  voxel_lookup <- stats::setNames(
    lapply(mvals, function(xx) which(mask_arr == xx)),
    as.character(mvals)
  )

  effects <- if (!is.null(effect_col)) sort(unique(as.character(coef_df[[effect_col]]))) else NA_character_
  to_loop <- expand.grid(pc = plot_cols, ee = effects, stringsAsFactors = FALSE)

  use_parallel <- ncpus > 1L && nrow(to_loop) > 1L
  if (use_parallel) {
    cl <- parallel::makeCluster(min(as.integer(ncpus), nrow(to_loop)))
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    doParallel::registerDoParallel(cl)
  } else {
    foreach::registerDoSEQ()
  }

  written <- foreach::foreach(
    ii = iterators::iter(to_loop, by = "row"),
    .packages = "RNifti"
  ) %dopar% {
    this_df <- if (!is.null(effect_col)) {
      coef_df[as.character(coef_df[[effect_col]]) == ii$ee, , drop = FALSE]
    } else {
      coef_df
    }

    timevals <- sort(unique(this_df[[time_col]]))
    n_t <- length(timevals)
    tr <- if (n_t > 1L) as.numeric(stats::median(diff(timevals))) else 1.0

    out_arr <- array(0, dim = c(dim3, n_t))
    p_vec <- this_df[[parcel_col]]
    t_vec <- this_df[[time_col]]
    v_vec <- as.numeric(this_df[[ii$pc]])

    for (mm in mvals) {
      vox_idx <- voxel_lookup[[as.character(mm)]]
      for (tt in seq_len(n_t)) {
        hit <- which(p_vec == mm & t_vec == timevals[tt])
        if (length(hit) != 1L) {
          eff_msg <- if (!is.null(effect_col)) paste0(", ", effect_col, "=", ii$ee) else ""
          stop(
            sprintf(
              "Expected exactly 1 row in coef_df for %s=%s, %s=%s%s, found %d.",
              parcel_col, as.character(mm), time_col, as.character(timevals[tt]), eff_msg, length(hit)
            ),
            call. = FALSE
          )
        }
        out_arr[vox_idx + (tt - 1L) * n_vox3] <- v_vec[hit]
      }
    }

    out_nii <- RNifti::asNifti(out_arr, reference = mask_img)
    pixdim_new <- mask_hdr$pixdim
    pixdim_new[5L] <- tr
    finite_vals <- out_arr[is.finite(out_arr)]
    cal_min <- if (length(finite_vals) > 0L) min(finite_vals) else 0
    cal_max <- if (length(finite_vals) > 0L) max(finite_vals) else 0
    out_nii <- RNifti::updateNifti(
      out_nii,
      template = list(pixdim = pixdim_new, cal_min = cal_min, cal_max = cal_max)
    )

    stem <- if (!is.null(effect_col)) {
      make.names(paste(ii$pc, ii$ee, sep = "_"))
    } else {
      make.names(ii$pc)
    }
    out_path <- file.path(output_dir, paste0(stem, ".nii.gz"))
    RNifti::writeNifti(out_nii, out_path, datatype = "float")
    out_path
  }

  invisible(unlist(written, use.names = FALSE))
}
