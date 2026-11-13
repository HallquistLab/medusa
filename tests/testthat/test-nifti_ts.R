test_that("read_nifti_ts and write_nifti_ts round-trip 4D runs, masks, atlases, and headers", {
  tmp_dir <- tempfile("nifti_ts_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  dim3 <- c(5L, 4L, 3L)
  n_spatial <- prod(dim3)
  n_time <- 18L

  affine <- matrix(
    c(
      -2.0, 0.0, 0.0, 40.0,
       0.0, 2.0, 0.0, -50.0,
       0.0, 0.0, 2.5, -20.0,
       0.0, 0.0, 0.0, 1.0
    ),
    nrow = 4L, byrow = TRUE
  )
  attr(affine, "code") <- 4L

  make_nii <- function(arr, path, tr = NULL, xyzt_units = NULL, datatype = "float") {
    nii <- RNifti::asNifti(arr)
    RNifti::pixdim(nii) <- if (length(dim(arr)) == 4L) c(2, 2, 2.5, if (is.null(tr)) 1.25 else tr) else c(2, 2, 2.5)
    RNifti::sform(nii) <- affine
    if (!is.null(xyzt_units)) {
      nii <- RNifti::updateNifti(nii, template = list(xyzt_units = as.integer(xyzt_units)))
    }
    RNifti::writeNifti(nii, path, datatype = datatype)
    path
  }

  set.seed(42L)
  arr1 <- array(rnorm(n_spatial * n_time, mean = 100, sd = 5), dim = c(dim3, n_time))
  arr2 <- array(rnorm(n_spatial * n_time, mean = 100, sd = 5), dim = c(dim3, n_time))

  mask_arr <- array(0L, dim = dim3)
  mask_arr[2:4, 1:3, 1:2] <- 1L

  atlas_arr <- array(0L, dim = dim3)
  atlas_arr[2:3, 1:2, 1] <- 10L
  atlas_arr[4, 2:3, 2] <- 20L
  atlas_arr[1, 1, 1] <- 30L # outside mask

  f_run1 <- make_nii(arr1, file.path(tmp_dir, "run1.nii.gz"), tr = 1.25)
  f_run2 <- make_nii(arr2, file.path(tmp_dir, "run2.nii"), tr = 1.25)
  f_mask <- make_nii(mask_arr, file.path(tmp_dir, "mask.nii.gz"), datatype = "int16")
  f_atlas <- make_nii(atlas_arr, file.path(tmp_dir, "atlas.nii.gz"), datatype = "int16")

  ts <- read_nifti_ts(
    niftis = c(r1 = f_run1, r2 = f_run2),
    mask = f_mask,
    atlas = f_atlas
  )

  expect_s3_class(ts, "medusa_nifti_ts")
  expect_equal(names(ts$bold), c("r1", "r2"))
  expect_equal(ts$tr, 1.25, tolerance = 1e-6)
  expect_equal(nrow(ts$voxel_map), 6L) # 4 in parcel 10 + 2 in parcel 20 (parcel 30 masked out)
  expect_equal(sort(unique(ts$parcel_ids)), c(10, 20))
  expect_equal(names(ts$parcel_ids), colnames(ts$bold$r1))
  expect_output(print(ts), "<medusa_nifti_ts> 2 runs \\| 6 voxels")

  # Verify world coordinates match affine * (ijk - 1)
  ijk0 <- cbind(ts$voxel_map$i - 1, ts$voxel_map$j - 1, ts$voxel_map$k - 1, 1)
  xyz_exp <- ijk0 %*% t(affine)
  expect_equal(as.matrix(ts$voxel_map[, c("x", "y", "z")]), xyz_exp[, 1:3], ignore_attr = TRUE, tolerance = 1e-6)

  # atlas_subset filters to parcel 20 only
  ts_sub <- read_nifti_ts(f_run1, mask = f_mask, atlas = f_atlas, atlas_subset = 20L, id_style = "vnum")
  expect_equal(ncol(ts_sub$bold$run1), 2L)
  expect_equal(colnames(ts_sub$bold$run1), c("v1", "v2"))
  expect_true(all(ts_sub$parcel_ids == 20))

  # Write both runs back out (.nii from .nii.gz ref, and .nii.gz with double precision)
  out1 <- file.path(tmp_dir, "out_r1.nii")
  out2 <- file.path(tmp_dir, "out_r2.nii.gz")
  written <- write_nifti_ts(ts, c(out1, out2), datatype = "double", compression = 2L, fill = -5)
  expect_equal(length(written), 2L)
  expect_true(all(file.exists(c(out1, out2))))

  rt <- read_nifti_ts(c(r1 = out1, r2 = out2), mask = f_mask, atlas = f_atlas)
  expect_equal(rt$bold$r1, ts$bold$r1, tolerance = 1e-10)
  expect_equal(rt$bold$r2, ts$bold$r2, tolerance = 1e-10)
  expect_equal(rt$tr, 1.25, tolerance = 1e-6)
  expect_equal(rt$geometry$affine, ts$geometry$affine, tolerance = 1e-6)

  # Non-indexed background voxels received fill = -5
  raw_out1 <- as.numeric(RNifti::readNifti(out1))
  bg_idx <- setdiff(seq_len(n_spatial), ts$voxel_map$linear_index)
  expect_true(all(raw_out1[bg_idx] == -5))
})

test_that("write_nifti_ts handles drop=TRUE matrices, column subsets/permutations, 3D maps, NA fill, and parallel I/O", {
  tmp_dir <- tempfile("nifti_ts_modes_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  arr <- array(seq_len(4L * 3L * 2L * 10L) * 0.25, dim = c(4L, 3L, 2L, 10L))
  mask <- array(0L, dim = c(4L, 3L, 2L))
  mask[c(2L, 5L, 9L, 15L, 22L)] <- 1L

  f_run1 <- file.path(tmp_dir, "r1.nii.gz")
  f_run2 <- file.path(tmp_dir, "r2.nii.gz")
  f_mask <- file.path(tmp_dir, "mask.nii.gz")
  nii_r1 <- RNifti::asNifti(arr); RNifti::pixdim(nii_r1) <- c(2, 2, 2, 1.5); RNifti::writeNifti(nii_r1, f_run1)
  nii_r2 <- RNifti::asNifti(arr * 2); RNifti::pixdim(nii_r2) <- c(2, 2, 2, 1.5); RNifti::writeNifti(nii_r2, f_run2)
  nii_m <- RNifti::asNifti(mask); RNifti::pixdim(nii_m) <- c(2, 2, 2); RNifti::writeNifti(nii_m, f_mask)

  # 1. drop = TRUE with voxel_by_time orientation round-trips with zero boilerplate
  mat_vt <- read_nifti_ts(f_run1, mask = f_mask, orientation = "voxel_by_time", drop = TRUE)
  expect_equal(dim(mat_vt), c(5L, 10L))
  f_vt_out <- file.path(tmp_dir, "vt_out.nii.gz")
  write_nifti_ts(mat_vt, f_vt_out, datatype = "double")
  mat_vt_rt <- read_nifti_ts(f_vt_out, mask = f_mask, orientation = "voxel_by_time", drop = TRUE)
  expect_equal(mat_vt_rt, mat_vt, ignore_attr = TRUE, tolerance = 1e-10)

  # 2. Permuted column subset matched automatically via voxel_map$voxel + fill = NA_real_
  ts <- read_nifti_ts(c(r1 = f_run1, r2 = f_run2), mask = f_mask, ncpus = 2L)
  ts_seq <- read_nifti_ts(c(r1 = f_run1, r2 = f_run2), mask = f_mask, ncpus = 1L)
  expect_equal(ts$bold, ts_seq$bold)

  sub_cols <- rev(colnames(ts$bold$r1)[c(1L, 3L, 5L)])
  sub_mat <- ts$bold$r1[, sub_cols, drop = FALSE]
  f_sub <- file.path(tmp_dir, "sub_na.nii.gz")
  write_nifti_ts(sub_mat, f_sub, reference = ts, fill = NA_real_)

  img_sub <- RNifti::readNifti(f_sub)
  sub_lin <- ts$voxel_map$linear_index[match(sub_cols, ts$voxel_map$voxel)]
  raw_sub <- matrix(as.numeric(img_sub), nrow = 24L)
  expect_equal(t(raw_sub[sub_lin, , drop = FALSE]), unname(sub_mat), tolerance = 1e-6)
  other_lin <- setdiff(seq_len(24L), sub_lin)
  expect_true(all(is.na(raw_sub[other_lin, ])))

  # 3. 1D static named vector writes a 3D volume (ndim == 3)
  stat_vec <- setNames(c(3.5, -2.25, 8.0), sub_cols)
  f_stat <- file.path(tmp_dir, "stat3d.nii.gz")
  write_nifti_ts(stat_vec, f_stat, reference = ts, fill = 0)
  hdr_stat <- RNifti::niftiHeader(f_stat)
  expect_equal(as.integer(hdr_stat$dim[1:4]), c(3L, 4L, 3L, 2L))
  vec_img <- as.numeric(RNifti::readNifti(f_stat))
  expect_equal(vec_img[sub_lin], unname(stat_vec), tolerance = 1e-6)
  expect_true(all(vec_img[other_lin] == 0))

  # 4. Parallel multi-file write (ncpus = 2L)
  par_files <- file.path(tmp_dir, c("par1.nii.gz", "par2.nii.gz"))
  write_nifti_ts(ts$bold, par_files, reference = ts, datatype = "double", ncpus = 2L)
  ts_par_rt <- read_nifti_ts(c(r1 = par_files[1L], r2 = par_files[2L]), mask = f_mask)
  expect_equal(ts_par_rt$bold, ts$bold, tolerance = 1e-10)
})

test_that("read_nifti_ts converts msec TR headers and rejects spatial or TR mismatches", {
  tmp_dir <- tempfile("nifti_ts_checks_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # xyzt_units = 2 (mm) + 16 (msec) = 18; pixdim[5] = 1500 msec -> 1.5 s
  f_ms <- file.path(tmp_dir, "msec.nii.gz")
  nii_ms <- RNifti::asNifti(array(rnorm(8 * 5), c(2L, 2L, 2L, 5L)))
  RNifti::pixdim(nii_ms) <- c(2, 2, 2, 1500)
  nii_ms <- RNifti::updateNifti(nii_ms, template = list(xyzt_units = 18L))
  RNifti::writeNifti(nii_ms, f_ms)
  ts_ms <- read_nifti_ts(f_ms)
  expect_equal(ts_ms$tr, 1.5, tolerance = 1e-6)

  # Inconsistent TR across runs errors unless tr is supplied explicitly
  f_tr2 <- file.path(tmp_dir, "tr2.nii.gz")
  nii_tr2 <- RNifti::asNifti(array(rnorm(8 * 5), c(2L, 2L, 2L, 5L)))
  RNifti::pixdim(nii_tr2) <- c(2, 2, 2, 2.0)
  nii_tr2 <- RNifti::updateNifti(nii_tr2, template = list(xyzt_units = 10L)) # mm (2) + sec (8)
  RNifti::writeNifti(nii_tr2, f_tr2)
  expect_error(read_nifti_ts(c(f_ms, f_tr2)), "inconsistent repetition times")
  expect_equal(read_nifti_ts(c(f_ms, f_tr2), tr = 1.5)$tr, 1.5)

  # Mismatched affine with identical dimensions fails before extraction
  f_shift <- file.path(tmp_dir, "shifted.nii.gz")
  nii_shift <- RNifti::readNifti(f_tr2)
  aff <- RNifti::xform(nii_shift)
  aff[1L, 4L] <- aff[1L, 4L] + 10
  attr(aff, "code") <- 1L
  RNifti::sform(nii_shift) <- aff
  RNifti::writeNifti(nii_shift, f_shift)
  expect_error(read_nifti_ts(c(f_tr2, f_shift), tr = 2.0), "Incompatible spatial NIfTI geometry")
})

test_that("read_nifti_ts integrates directly with estimate_crossrun_hrf and deconvolve_crossrun_hrf", {
  set.seed(20261113L)
  tmp_dir <- tempfile("nifti_ts_pipeline_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  tr <- 1.0
  n_time <- 80L
  dim3 <- c(3L, 2L, 2L)
  lib <- hrf_library(tr, peak_delays = c(5.0, 6.0), dispersions = 1.0, undershoot_ratios = 6)

  atlas <- array(0L, dim = dim3)
  atlas[1:2, 1, 1] <- 1L
  atlas[1:2, 2, 1] <- 2L
  f_atlas <- file.path(tmp_dir, "atlas.nii.gz")
  nii_at <- RNifti::asNifti(atlas)
  RNifti::pixdim(nii_at) <- c(2, 2, 2)
  RNifti::writeNifti(nii_at, f_atlas)

  f_runs <- character(2L)
  events <- list()
  for (r in 1:2) {
    onsets <- seq(8, 68, by = 10)
    des <- hrf_event_design_cpp(n_time, tr, 0, onsets, rep(0, length(onsets)),
      seq_along(onsets), length(onsets), lib[, 2L, drop = FALSE])[[1L]]
    sig <- as.numeric(des %*% rep(1.5, length(onsets)))
    arr <- array(0, dim = c(dim3, n_time))
    for (i in which(atlas > 0L)) {
      arr[i + prod(dim3) * (0:(n_time - 1L))] <- 100 + sig + rnorm(n_time, sd = 0.2)
    }
    f_runs[r] <- file.path(tmp_dir, sprintf("run%d.nii.gz", r))
    nii_r <- RNifti::asNifti(arr)
    RNifti::pixdim(nii_r) <- c(2, 2, 2, tr)
    RNifti::writeNifti(nii_r, f_runs[r])
    events[[paste0("run", r)]] <- data.frame(onset = onsets, duration = 0, condition = "stim")
  }
  names(f_runs) <- c("run1", "run2")

  ts <- read_nifti_ts(f_runs, atlas = f_atlas)
  hrf_fit <- estimate_crossrun_hrf(
    ts$bold, events, ts$tr, lib,
    spatial_groups = ts$parcel_ids, pooling = "parcel", min_parcel_voxels = 2L,
    fallback_kernel = 2L, noise_model = "ols"
  )
  dec <- deconvolve_crossrun_hrf(ts$bold, hrf_fit, lambda = 0.5, scaling = "noise_sd")

  out_decon <- file.path(tmp_dir, c("decon_run1.nii.gz", "decon_run2.nii.gz"))
  write_nifti_ts(dec$activity, out_decon, reference = ts, datatype = "double")
  dec_rt <- read_nifti_ts(setNames(out_decon, c("run1", "run2")), atlas = f_atlas)
  expect_equal(dec_rt$bold$run1, dec$activity$run1, ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(dec_rt$bold$run2, dec$activity$run2, ignore_attr = TRUE, tolerance = 1e-10)
})
