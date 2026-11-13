#' Backproject parcel-level MEDuSA coefficients onto a 4D NIfTI brain mask
#'
#' @description
#' Maps parcel-by-time model statistics (for example from [mixed_by()]) back
#' into 3D/4D NIfTI volumes (3D space across event-locked timepoints) using an
#' integer-valued 3D parcellation mask and a streaming C++ NIfTI writer.
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
#' @param compression Integer gzip compression level (`0`--`9`) when writing
#'   `.nii.gz` files. Default `2L` for fast deflate compression.
#' @param fill Numeric scalar written to non-parcel background voxels (for
#'   example `0` or `NA_real_`). Default `0`.
#'
#' @return Invisibly returns a character vector of written `.nii.gz` file paths.
#' @importFrom checkmate assert_data_frame assert_string assert_file_exists assert_character assert_subset assert_integerish assert_number
#' @importFrom parallel makeCluster stopCluster clusterCall parLapply mclapply
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
  ncpus = 4L,
  compression = 2L,
  fill = 0
) {
  checkmate::assert_data_frame(coef_df, min.rows = 1L)
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
  checkmate::assert_integerish(ncpus, lower = 1L, len = 1L, any.missing = FALSE)
  checkmate::assert_integerish(compression, lower = 0L, upper = 9L, len = 1L, any.missing = FALSE)
  checkmate::assert_number(fill, na.ok = TRUE)
  ncpus <- as.integer(ncpus)
  compression <- as.integer(compression)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  brain_mask_path <- normalizePath(brain_mask, winslash = "/", mustWork = TRUE)
  mask_hdr <- RNifti::niftiHeader(brain_mask_path)
  n_dim <- as.integer(mask_hdr$dim[1L])
  if (n_dim != 3L && !(n_dim == 4L && mask_hdr$dim[5L] <= 1L)) {
    stop("brain_mask does not appear to be a 3-D mask image.", call. = FALSE)
  }
  dim3 <- as.integer(mask_hdr$dim[2:4])
  n_vox3 <- prod(dim3)
  mask_vec <- as.numeric(RNifti::readNifti(brain_mask_path))
  if (length(mask_vec) != n_vox3) {
    stop("Unexpected voxel count when reading brain_mask.", call. = FALSE)
  }

  mvals <- sort(unique(coef_df[[parcel_col]]))
  if (length(mvals) == 0L) {
    stop("No non-missing parcel values found in coef_df[['", parcel_col, "']].", call. = FALSE)
  }
  imgvals <- sort(unique(mask_vec))
  imgvals <- imgvals[is.finite(imgvals) & imgvals != 0]
  checkmate::assert_subset(mvals, imgvals)

  parcel_slot <- as.integer(match(as.numeric(mask_vec), as.numeric(mvals), nomatch = 0L))
  n_m <- length(mvals)

  effects <- if (!is.null(effect_col)) sort(unique(as.character(coef_df[[effect_col]]))) else NA_character_
  tasks <- vector("list", length(effects) * length(plot_cols))
  task_idx <- 0L

  for (ee in effects) {
    this_df <- if (!is.null(effect_col)) {
      coef_df[!is.na(coef_df[[effect_col]]) & as.character(coef_df[[effect_col]]) == ee, , drop = FALSE]
    } else {
      coef_df
    }

    timevals <- sort(unique(this_df[[time_col]]))
    n_t <- length(timevals)
    if (n_t == 0L) {
      stop("No non-missing timepoints found in coef_df[['", time_col, "']].", call. = FALSE)
    }
    tr <- if (n_t > 1L) as.numeric(stats::median(diff(timevals))) else 1.0

    p_idx <- match(as.numeric(this_df[[parcel_col]]), as.numeric(mvals))
    t_idx <- match(this_df[[time_col]], timevals)
    valid_row <- !is.na(p_idx) & !is.na(t_idx)
    cell_idx <- p_idx[valid_row] + (t_idx[valid_row] - 1L) * n_m
    counts <- tabulate(cell_idx, nbins = n_m * n_t)
    if (any(counts != 1L) || !all(valid_row)) {
      bad_cells <- which(counts != 1L)
      if (length(bad_cells) > 0L) {
        first_bad <- bad_cells[1L]
        mm_i <- ((first_bad - 1L) %% n_m) + 1L
        tt_i <- ((first_bad - 1L) %/% n_m) + 1L
        eff_msg <- if (!is.null(effect_col)) paste0(", ", effect_col, "=", ee) else ""
        stop(
          sprintf(
            "Expected exactly 1 row in coef_df for %s=%s, %s=%s%s, found %d.",
            parcel_col, as.character(mvals[mm_i]),
            time_col, as.character(timevals[tt_i]),
            eff_msg, counts[first_bad]
          ),
          call. = FALSE
        )
      }
      stop("coef_df contains missing values in parcel or time columns.", call. = FALSE)
    }

    for (pc in plot_cols) {
      p_mat <- matrix(0, nrow = n_m, ncol = n_t)
      p_mat[cell_idx] <- as.numeric(this_df[[pc]])
      stem <- if (!is.null(effect_col)) {
        make.names(paste(pc, ee, sep = "_"))
      } else {
        make.names(pc)
      }
      out_path <- file.path(output_dir, paste0(stem, ".nii.gz"))
      task_idx <- task_idx + 1L
      tasks[[task_idx]] <- list(out_path = out_path, p_mat = p_mat, tr = tr)
    }
  }

  written_list <- medusa_parallel_lapply(
    tasks,
    function(task, mask_ref, slots, fill_val, comp_val) {
      medusa_ns <- asNamespace("medusa")
      write_fn <- get("write_nifti_parcel_cpp", envir = medusa_ns)
      write_fn(mask_ref, task$out_path, slots, task$p_mat, task$tr, fill_val, "float", comp_val)
      task$out_path
    },
    mask_ref = brain_mask_path,
    slots = parcel_slot,
    fill_val = as.numeric(fill),
    comp_val = compression,
    ncpus = ncpus
  )
  written <- unlist(written_list, use.names = FALSE)

  invisible(written)
}
