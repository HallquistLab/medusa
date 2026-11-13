#' Function to perform voxelwise deconvolution on an fMRI dataset using the fMRI model arguments object
#'
#' @param niftis A vector of processed fMRI timeseries images (4D files) to be deconvolved
#' @param add_metadata A data.frame with one row per value of \code{niftis}, saved
#'   as atlas-level metadata. The \code{id}, \code{run}, and optional \code{session}
#'   columns are also saved in per-file sidecars for event alignment.
#' @param out_dir Base output directory for deconvolved time series files. Default is \code{getwd()}.
#' @param out_file_expression Expression evaluated to resolve the filename for the deconvolved csv files. Default is
#'   to prefix \code{"decon_"}, replace slashes in \code{niftis} with periods, and add
#'   the atlas image name. Note that the suffix \code{_deconvolved.csv.gz} or \code{_original.csv.gz} will be added
#'   to the expression, so don't pass this piece.
#' @param log_file Name (and path) of log file for any deconvolution errors or messages
#' @param TR the repetition time of the sequence in seconds. Required
#' @param time_offset The number of seconds that will be subtracted or added to the time field. Default: 0. Useful if some number of volumes
#'   have been dropped from the NIfTI data prior to deconvolution.
#' @param atlas_files optional atlas file vector specifying voxels used in deconvolution. If omitted, perform whole-brain deconvolution
#' @param mask an optional character string specifying a mask that should be used to constrain bounds of deconvolution.
#' @param nprocs The number of processors to use simultaneously for deconvolution
#' @param save_original_ts Whether to save the voxelwise BOLD data prior to deconvolution (for comparison/diagnosis). Default: TRUE
#' @param algorithm Which deconvolution algorithm to use for deconvolving voxelwise time series. Default: "bush2011". Alternative is
#'   "bush2015", which implements a resampling approach as well, or "reglin"/"regularized_linear",
#'   which estimates a continuous latent activity time series using regularized linear deconvolution
#'   with optional HRF tuning.
#' @param decon_settings A list of settings passed to the deconvolution algorithm.
#' @param afni_dir Retained for compatibility; ignored. NIfTI extraction is internal
#'   and does not require AFNI.
#'
#' @details The Bush 2011 algorithm uses the package's internal RcppArmadillo
#'   implementation and does not require an external executable.
#'   Atlas, mask, and input images must share a spatial grid. Constant signals
#'   and signals containing nonfinite observations are excluded before fitting;
#'   their output rows are retained with missing estimates. The \code{decon_status}
#'   column identifies \code{"ok"}, \code{"constant"}, or \code{"nonfinite"} signals.
#'   Each deconvolution file has a \code{.metadata.csv} sidecar containing its
#'   source NIfTI path and any \code{id}, \code{run}, and \code{session} columns
#'   in \code{add_metadata}. These identify files during event alignment, including
#'   custom filenames. If absent, alignment falls back to the source path or
#'   requires an explicit \code{file_manifest}.
#'
#'   If you want to use subject metadata to name the output file, use \code{this_subj} in your \code{out_file_expression}, which will give you access to a one-row
#'   data.frame containing the metadata for the current subject in the loop.
#'
#'   Set \code{decon_settings$scaling} explicitly to \code{"zscore"},
#'   \code{"noise_sd"}, \code{"demean"}, \code{"none"}, or \code{"psc"}. The default is z-scoring.
#'   An explicit mode overrides \code{normalize}/\code{demean}; otherwise
#'   \code{normalize = FALSE} selects demeaning unless \code{demean = FALSE} too.
#'   The voxelwise Bush paths require z-scoring because their learning-rate and
#'   bounded-encoding/resampling behavior is sensitive to scale. Their output remains
#'   a bounded encoding. Use the linear algorithm for estimates in alternative units.
#'   For linear methods, z-scoring is recommended for standardized trial responses;
#'   demeaning retains input units, and percent signal change requires an appropriate
#'   positive baseline rather than residualized/zero-mean data. Raw scale does not by
#'   itself degrade the conditioning of the fixed linear inverse, but it changes
#'   signal weights during global tuning. See \code{\link{deconvolve_reglin}}.
#'
#'   Scaling uses the observed run once, before any zero tail is appended
#'   (\code{decon_settings$tail_padding}, defaulting to \code{0L} for
#'   \code{reglin}/\code{regularized_linear} and \code{length(decon_settings$kernel)}
#'   for \code{bush2011}/\code{bush2015});
#'   solver normalization is disabled, as is linear-solver demeaning of padded data. Each output has
#'   a \code{.scaling.csv} sidecar containing voxel number, mode, subtracted center,
#'   divisor, and validity status. Its \code{.metadata.csv} records algorithm, scaling,
#'   output units, and the scaling filename. Original-series files contain
#'   \code{BOLD_input}, the transformed observations supplied before padding, plus
#'   the legacy \code{BOLD_z} alias only for z-scored inputs. Existing output files
#'   are skipped; regenerate them to adopt the new scaling behavior and sidecars.
#' @examples
#'   \dontrun{
#'
#'     #name outputs according to subject metadata
#'     xx <- voxelwise_deconvolution(
#'       niftis="/proj/mnhallqlab/user/michael/test_nifti.nii.gz",
#'       out_dir="/proj/mnhallqlab/user/michael/decon_outputs",
#'       out_file_expression = expression(
#'         paste0(
#'           this_subj$subid,
#'           "_run",
#'           this_subj$run_num,
#'           "_",
#'           atlas_img_name
#'         )
#'       )
#'     )
#'   }
#'
#' @return Nothing (invisible NULL).
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster clusterCall
#' @importFrom checkmate assert_file_exists
#' @importFrom foreach %do% %dopar% foreach registerDoSEQ
#' @importFrom dplyr mutate across select left_join
#' @importFrom readr write_csv
#' @importFrom stats setNames
#' @export
voxelwise_deconvolution <- function(
  niftis, add_metadata=NULL, out_dir=getwd(), out_file_expression=NULL,
  log_file=file.path(out_dir, "deconvolve_errors"),
  TR = NULL, time_offset=0, atlas_files=NULL, mask=NULL, nprocs=20, save_original_ts=TRUE, algorithm="bush2011",
  decon_settings=list(nev_lr = .01, #neural events learning rate (default in algorithm)
    epsilon = .005, #convergence criterion (default)
    beta = 60, #best from Bush 2015 update
    kernel = spm_hrf(TR)$hrf, #canonical SPM difference of gammas
    n_resample = 25), afni_dir=NULL) { #for Bush 2015 only

  sapply(niftis, checkmate::assert_file_exists)
  checkmate::assert_data_frame(add_metadata, nrows=length(niftis), null.ok=TRUE)
  sapply(atlas_files, checkmate::assert_file_exists)
  checkmate::assert_numeric(TR, lower=0.01)
  algorithm <- match.arg(
    algorithm,
    c("bush2011", "bush2015", "reglin", "regularized_linear")
  )

  checkmate::assert_list(decon_settings)
  decon_settings <- utils::modifyList(list(nev_lr = 0.01, epsilon = 0.005,
    beta = 60, n_resample = 25), decon_settings)

  # Explicit scaling overrides legacy flags. The default remains z-scoring.
  # Bush's bounded encoding/fixed learning rate and resampling implementation
  # assume standardized inputs; alternative units are supported by reglin.
  input_scaling <- resolve_decon_scaling(decon_settings[["scaling"]],
    normalize = if (is.null(decon_settings[["normalize"]])) TRUE else decon_settings[["normalize"]],
    demean = if (is.null(decon_settings[["demean"]])) TRUE else decon_settings[["demean"]])
  if (algorithm %in% c("bush2011", "bush2015") && input_scaling != "zscore") {
    stop("Voxelwise Bush/Cisler deconvolution requires scaling = 'zscore'; ",
      "use reglin for alternative input units.")
  }

  if (is.null(decon_settings$kernel)) {
    stop("At present no alternative way to set hrf_pad if $kernel not passed in decon_settings")
  }
  if (is.null(decon_settings$tail_padding)) {
    hrf_pad <- if (algorithm %in% c("reglin", "regularized_linear")) 0L else length(decon_settings$kernel)
  } else {
    checkmate::assert_int(decon_settings$tail_padding, lower = 0L)
    hrf_pad <- as.integer(decon_settings$tail_padding)
  }

  # Configuration errors must be raised before workers can turn fit failures
  # into missing estimates in otherwise successful output files.
  if (algorithm %in% c("reglin", "regularized_linear") &&
      identical(decon_settings[["backend"]], "wiener") &&
      identical(decon_settings[["wiener_shrinkage"]], "psd")) {
    candidate_kernels <- build_decon_kernel_grid(decon_settings$kernel,
      decon_settings$kernel_grid,
      if (is.null(decon_settings$hrf_lags)) 0L else decon_settings$hrf_lags,
      normalize = decon_settings$kernel_normalize)
    n_lambda <- if (!is.null(decon_settings[["lambda"]])) 1L else
      if (is.null(decon_settings$lambda_grid)) 20L else length(unique(decon_settings$lambda_grid))
    validate_wiener_tuning("wiener", "psd", n_lambda, length(candidate_kernels))
  }

  if (is.null(atlas_files)) {
    message("Performing whole-brain voxelwise deconvolution")
    if (is.null(mask)) {
      stop("Cannot conduct voxelwise deconvolution without a relevant mask")
    } else {
      checkmate::assert_file_exists(mask)
    }

    atlas_files <- mask
  }
  
  zero_thresh <- 1e-4 # for binarizing/indexing
  atlas_imgs <- lapply(atlas_files, RNifti::readNifti)
  mask_file <- mask
  mask_geom <- NULL
  if (!is.null(mask)) {
    checkmate::assert_file_exists(mask)
    mask_geom <- roi_image_geometry(mask_file)
    mask <- RNifti::readNifti(mask)
    # ensure that mask is binary
    mask <- is.finite(mask) & mask >= zero_thresh
  }

  # Pre-validate all atlases once before reading or deconvolving NIfTI runs so
  # multi-atlas workflows decompress each 4D NIfTI only once per run.
  atlas_specs <- vector("list", length(atlas_files))
  for (ai in seq_along(atlas_files)) {
    cat("Working on atlas: ", atlas_files[ai], "\n")
    aimg <- atlas_imgs[[ai]]

    atlas_geometry <- roi_image_geometry(atlas_files[ai])
    if (atlas_geometry$volumes != 1L) stop("Deconvolution atlases must be single-volume images")
    if (!is.null(mask_file)) {
      validate_roi_geometry(atlas_geometry, mask_geom, atlas_files[ai], mask_file)
    }
    if (!is.null(mask)) {
      aimg <- aimg * mask
    }

    # One shared index vector defines both metadata and extraction order.
    # NIfTI and R arrays both store i first, then j, then k within each volume.
    voxel_indices <- which(is.finite(aimg) & aimg > zero_thresh)
    if (length(voxel_indices) == 0L) stop("No voxels remain in atlas after masking: ", atlas_files[ai])
    a_indices <- arrayInd(voxel_indices, dim(aimg)[seq_len(min(3L, length(dim(aimg))))])
    if (ncol(a_indices) == 2L) a_indices <- cbind(a_indices, 1L)
    xyz <- cbind(a_indices - 1L, 1) %*% t(atlas_geometry$affine)
    a_coordinates <- cbind(a_indices, xyz[, 1:3, drop = FALSE])

    a_coordinates <- as.data.frame(a_coordinates) %>%
      setNames(c("i", "j", "k", "x", "y", "z")) %>%
      dplyr::mutate(
        vnum = seq_along(voxel_indices),
        atlas_value = as.numeric(aimg[voxel_indices]),
        atlas_name = basename(atlas_files[ai])
      ) %>%
      mutate(across(c("x", "y", "z"), ~round(.x, 2))) %>%
      dplyr::select("vnum", "atlas_value", everything())

    atlas_img_name <- basename(sub("\\.nii(\\.gz)?$", "", atlas_files[ai], perl = TRUE))
    dir.create(file.path(out_dir, atlas_img_name, "deconvolved"), showWarnings = FALSE, recursive = TRUE)
    if (isTRUE(save_original_ts)) {
      dir.create(file.path(out_dir, atlas_img_name, "original"), showWarnings = FALSE, recursive = TRUE)
    }

    atlas_specs[[ai]] <- list(
      atlas_file = atlas_files[ai],
      atlas_img_name = atlas_img_name,
      atlas_geometry = atlas_geometry,
      voxel_indices = voxel_indices,
      a_coordinates = a_coordinates
    )
  }

  out_file_expr <- if (is.null(out_file_expression)) {
    expression(paste0("decon_", gsub("[:/\\]+", ".", niftis[si]), "_", atlas_img_name))
  } else {
    out_file_expression
  }

  if (!is.null(add_metadata)) {
    add_metadata$.nifti <- niftis
  }

  # Setup cluster only when multiple NIfTIs are processed in parallel; keep
  # single-worker execution in-process so testthat mocked bindings work cleanly.
  use_parallel <- nprocs > 1 && length(niftis) > 1L
  if (use_parallel) {
    n_workers <- min(as.integer(nprocs), length(niftis))
    can_fork <- identical(.Platform$OS.type, "unix") && !nzchar(Sys.getenv("RSTUDIO"))
    if (can_fork) {
      registerDoParallel(cores = n_workers)
      on.exit(try(registerDoSEQ(), silent = TRUE), add = TRUE)
    } else {
      cl <- makeCluster(n_workers)
      on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
      lib_paths <- .libPaths()
      parallel::clusterCall(cl, function(lp) .libPaths(lp), lib_paths)
      registerDoParallel(cl)
    }
  } else {
    registerDoSEQ()
  }

  dopar_pkgs <- c(
    "dplyr", "readr", "data.table",
    if ("medusa" %in% loadedNamespaces()) "medusa",
    "foreach", "iterators"
  )
  dopar_exports <- c(
    "mat2df", "extract_nifti_voxel_ts", "roi_image_geometry",
    "validate_roi_geometry", "decon_signal_status", "scale_decon_input",
    "resolve_decon_scaling", "estimate_decon_noise_sd",
    "deconvolve_reglin", "deconvolve_nlreg", "deconvolve_nlreg_resample"
  )

  ff <- foreach(
    si = seq_along(niftis),
    .packages = dopar_pkgs,
    .export = dopar_exports
  ) %dopar% {
    if (!is.null(add_metadata)) {
      this_subj <- add_metadata %>% dplyr::slice(si)
    }

    n_atlases <- length(atlas_specs)
    out_stems <- character(n_atlases)
    out_names <- character(n_atlases)
    needed <- logical(n_atlases)
    for (ai in seq_len(n_atlases)) {
      atlas_img_name <- atlas_specs[[ai]]$atlas_img_name
      out_stems[ai] <- eval(out_file_expr)
      out_names[ai] <- file.path(
        out_dir, atlas_img_name, "deconvolved",
        paste0(out_stems[ai], "_deconvolved.csv.gz")
      )
      if (file.exists(out_names[ai])) {
        message("Deconvolved file already exists: ", out_names[ai])
      } else {
        needed[ai] <- TRUE
      }
    }

    if (!any(needed)) {
      return(NULL)
    }

    cat("  Deconvolving subject: ", niftis[si], "\n")
    input_geometry <- roi_image_geometry(niftis[si])
    input_geometry$volumes <- 1L
    needed_idx <- which(needed)
    for (ai in needed_idx) {
      validate_roi_geometry(
        atlas_specs[[ai]]$atlas_geometry,
        input_geometry,
        atlas_specs[[ai]]$atlas_file,
        niftis[si]
      )
    }

    # Extract the union of needed voxel indices in a single pass directly into
    # time x voxels (T x V_union), avoiding repeated .nii.gz decompression and
    # eliminating intermediate V x T <-> T x V matrix transposes.
    if (length(needed_idx) == 1L) {
      union_indices <- atlas_specs[[needed_idx]]$voxel_indices
    } else {
      union_indices <- sort(unique(unlist(lapply(atlas_specs[needed_idx], "[[", "voxel_indices"), use.names = FALSE)))
    }
    union_tv <- extract_nifti_voxel_ts(niftis[si], union_indices, TRUE)
    n_time <- nrow(union_tv)

    for (ai in needed_idx) {
      spec <- atlas_specs[[ai]]
      atlas_img_name <- spec$atlas_img_name
      out_name <- out_names[ai]

      if (identical(spec$voxel_indices, union_indices)) {
        obs_tv <- union_tv
      } else {
        obs_tv <- union_tv[, match(spec$voxel_indices, union_indices), drop = FALSE]
      }
      n_vox <- ncol(obs_tv)

      # Screen before scaling in native time x voxels orientation: a constant
      # signal would otherwise become NaN and contaminate global tuning.
      decon_status <- decon_signal_status(obs_tv)
      valid_voxels <- decon_status == "ok"
      scaled_tv <- matrix(NA_real_, nrow = n_time, ncol = n_vox)
      centers <- divisors <- rep(NA_real_, n_vox)
      if (any(valid_voxels)) {
        scaled <- scale_decon_input(
          obs_tv[, valid_voxels, drop = FALSE],
          scaling = input_scaling,
          noise_rho = decon_settings$prewhiten_rho,
          noise_psi1 = decon_settings$noise_psi1,
          noise_psi2 = decon_settings$noise_psi2,
          noise_psi3 = decon_settings$noise_psi3,
          noise_psi4 = decon_settings$noise_psi4
        )
        scaled_tv[, valid_voxels] <- scaled$y
        centers[valid_voxels] <- scaled$center
        divisors[valid_voxels] <- scaled$scale
      }
      if (any(!valid_voxels)) {
        warning(sum(!valid_voxels), " constant or nonfinite voxel signals excluded from deconvolution: ", niftis[si])
      }

      if (!any(valid_voxels)) {
        deconv_tv <- matrix(NA_real_, nrow = n_time, ncol = n_vox)
      } else {
        alg_input <- scaled_tv[, valid_voxels, drop = FALSE]
        if (hrf_pad > 0L) {
          alg_input <- rbind(alg_input, matrix(0, nrow = hrf_pad, ncol = ncol(alg_input)))
        }
        alg_input <- as.matrix(alg_input)

        if (algorithm == "bush2015") {
          n_valid <- ncol(alg_input)
          n_padded <- nrow(alg_input)
          deconv_valid <- matrix(NA_real_, nrow = n_padded, ncol = n_valid)
          for (vj in seq_len(n_valid)) {
            vox_ts <- as.vector(alg_input[, vj])
            reg <- tryCatch(
              deconvolve_nlreg_resample(
                vox_ts,
                kernel = decon_settings$kernel,
                nev_lr = decon_settings$nev_lr,
                epsilon = decon_settings$epsilon,
                n_resample = decon_settings$n_resample,
                beta = decon_settings$beta,
                normalize = FALSE
              ),
              error = function(e) {
                cat("Problem deconvolving: ", niftis[si], as.character(e), "\n", file = log_file, append = TRUE)
                rep(NA_real_, n_padded)
              }
            )
            if (is.list(reg)) reg <- reg$NEVmean
            deconv_valid[, vj] <- as.numeric(reg)
          }
        } else if (algorithm %in% c("reglin", "regularized_linear")) {
          deconv_valid <- tryCatch(
            deconvolve_reglin(
              BOLDobs = alg_input,
              kernel = decon_settings$kernel,
              lambda = decon_settings[["lambda"]],
              lambda_grid = if (is.null(decon_settings$lambda_grid)) 10^seq(-4, 2, length.out = 20) else decon_settings$lambda_grid,
              penalty = if (is.null(decon_settings$penalty)) "diff2" else decon_settings$penalty,
              kernel_grid = decon_settings$kernel_grid,
              hrf_lags = if (is.null(decon_settings$hrf_lags)) 0L else decon_settings$hrf_lags,
              tune_by = if (is.null(decon_settings$tune_by)) "global" else decon_settings$tune_by,
              normalize = FALSE,
              demean = FALSE,
              prewhiten_gcv = if (is.null(decon_settings$prewhiten_gcv)) FALSE else decon_settings$prewhiten_gcv,
              prewhiten_rho = decon_settings$prewhiten_rho,
              backend = if (is.null(decon_settings$backend)) "cpp" else decon_settings$backend,
              wiener_shrinkage = if (is.null(decon_settings$wiener_shrinkage)) "scalar" else decon_settings$wiener_shrinkage,
              wiener_psd_smooth = if (is.null(decon_settings$wiener_psd_smooth)) 7L else decon_settings$wiener_psd_smooth,
              wiener_psd_floor = if (is.null(decon_settings$wiener_psd_floor)) 1e-8 else decon_settings$wiener_psd_floor,
              gcv_rule = if (is.null(decon_settings$gcv_rule)) "min" else decon_settings$gcv_rule,
              trim_kernel = TRUE,
              kernel_normalize = if (is.null(decon_settings$kernel_normalize)) "auto" else decon_settings$kernel_normalize
            ),
            error = function(e) {
              cat("Problem deconvolving: ", niftis[si], as.character(e), "\n", file = log_file, append = TRUE)
              matrix(NA_real_, nrow = nrow(alg_input), ncol = ncol(alg_input))
            }
          )
        } else if (algorithm == "bush2011") {
          deconv_valid <- tryCatch(
            deconvolve_nlreg(
              BOLDobs = alg_input,
              kernel = decon_settings$kernel,
              nev_lr = decon_settings$nev_lr,
              epsilon = decon_settings$epsilon,
              beta = decon_settings$beta,
              normalize = FALSE
            ),
            error = function(e) {
              cat("Problem deconvolving: ", niftis[si], as.character(e), "\n", file = log_file, append = TRUE)
              matrix(NA_real_, nrow = nrow(alg_input), ncol = ncol(alg_input))
            }
          )
        }

        deconv_valid <- deconv_valid[seq_len(n_time), , drop = FALSE]
        deconv_tv <- matrix(NA_real_, nrow = n_time, ncol = n_vox)
        deconv_tv[, valid_voxels] <- deconv_valid
      }

      # Melt to legacy (vnum, volume) long layout via a single final transpose
      deconv_melt <- mat2df(unname(t(deconv_tv)), value_name = "decon", varnames = c("vnum", "volume"))
      deconv_melt$time <- (deconv_melt$volume - 1) * TR + time_offset

      deconv_df <- deconv_melt %>%
        dplyr::mutate(vnum = as.numeric(.data$vnum)) %>%
        left_join(spec$a_coordinates, by = "vnum") %>%
        dplyr::select(-"i", -"j", -"k")

      deconv_df$decon_status <- decon_status[deconv_df$vnum]
      readr::write_csv(deconv_df, file = out_name)

      scaling_file <- paste0(out_name, ".scaling.csv")
      readr::write_csv(
        data.frame(
          vnum = seq_along(centers),
          scaling = input_scaling,
          center = centers,
          scale = divisors,
          decon_status = decon_status
        ),
        scaling_file
      )
      identity_metadata <- data.frame(
        .nifti = niftis[si],
        algorithm = algorithm,
        scaling = input_scaling,
        scaling_file = basename(scaling_file),
        decon_units = if (algorithm %in% c("bush2011", "bush2015")) {
          "bounded_encoding"
        } else {
          switch(input_scaling, zscore = "input_sd", noise_sd = "noise_sd", psc = "percent_signal_change", "input_units")
        }
      )
      if (!is.null(add_metadata)) {
        for (column in intersect(c("id", "run", "session"), names(this_subj))) {
          identity_metadata[[column]] <- this_subj[[column]]
        }
      }
      readr::write_csv(identity_metadata, file = paste0(out_name, ".metadata.csv"))

      if (isTRUE(save_original_ts)) {
        to_deconvolve_melt <- mat2df(unname(t(scaled_tv)), value_name = "BOLD_input", varnames = c("vnum", "time"))
        to_deconvolve_melt$time <- (to_deconvolve_melt$time - 1) * TR + time_offset

        orig_df <- to_deconvolve_melt %>%
          dplyr::mutate(vnum = as.numeric(.data$vnum)) %>%
          left_join(spec$a_coordinates, by = "vnum") %>%
          dplyr::select(-"i", -"j", -"k")

        if (input_scaling == "zscore") orig_df$BOLD_z <- orig_df$BOLD_input

        orig_out_name <- file.path(out_dir, atlas_img_name, "original", paste0(out_stems[ai], "_original.csv.gz"))
        readr::write_csv(orig_df, file = orig_out_name)
      }
    }
    NULL
  }

  # Write atlas-level metadata once per atlas, deduplicating across re-runs
  if (!is.null(add_metadata)) {
    for (ai in seq_along(atlas_specs)) {
      atlas_img_name <- atlas_specs[[ai]]$atlas_img_name
      out_name <- file.path(out_dir, atlas_img_name, paste0(atlas_img_name, "_metadata.csv"))
      if (file.exists(out_name)) {
        existing_meta <- tryCatch(
          as.data.frame(readr::read_csv(out_name, show_col_types = FALSE)),
          error = function(e) NULL
        )
        if (!is.null(existing_meta) && nrow(existing_meta) > 0L) {
          key_cols <- intersect(c(".nifti", "id", "run", "session"), intersect(names(existing_meta), names(add_metadata)))
          if (!length(key_cols)) key_cols <- intersect(names(existing_meta), names(add_metadata))
          combined_meta <- dplyr::bind_rows(
            dplyr::mutate(existing_meta, dplyr::across(dplyr::everything(), as.character)),
            dplyr::mutate(as.data.frame(add_metadata), dplyr::across(dplyr::everything(), as.character))
          )
          if (length(key_cols)) {
            combined_meta <- combined_meta[!duplicated(combined_meta[, key_cols, drop = FALSE], fromLast = TRUE), , drop = FALSE]
          } else {
            combined_meta <- unique(combined_meta)
          }
          readr::write_csv(combined_meta, file = out_name, append = FALSE, col_names = TRUE)
        } else {
          readr::write_csv(add_metadata, file = out_name, append = FALSE, col_names = TRUE)
        }
      } else {
        readr::write_csv(add_metadata, file = out_name, append = FALSE, col_names = TRUE)
      }
    }
  }

  return(invisible(NULL))
}
