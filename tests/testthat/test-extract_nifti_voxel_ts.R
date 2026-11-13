test_that("native extraction preserves requested voxel order and all volumes", {
  image <- array(seq_len(2L * 3L * 2L * 4L) / 7, c(2L, 3L, 2L, 4L))
  indices <- c(12L, 2L, 7L, 1L)
  for (extension in c(".nii", ".nii.gz")) {
    path <- tempfile(fileext = extension)
    RNifti::writeNifti(image, path)
    expected <- matrix(as.numeric(RNifti::readNifti(path)), nrow = 12L)[indices, , drop = FALSE]
    expect_equal(medusa:::extract_nifti_voxel_ts(path, indices), expected)
    expect_equal(dim(medusa:::extract_nifti_voxel_ts(path, 7L)), c(1L, 4L))
    unlink(path)
  }
})

test_that("native extraction applies on-disk integer scaling and retains missing data", {
  path <- tempfile(fileext = ".nii")
  on.exit(unlink(path))
  image <- array(seq_len(48L), c(2L, 2L, 2L, 6L))
  RNifti::writeNifti(image, path, datatype = "int16")
  # Set scaling on the stored integer data without rescaling the payload.
  # NIfTI-1 stores float32 scl_slope and scl_inter at byte offsets 112 and 116.
  con <- file(path, "r+b")
  seek(con, 112L, rw = "write")
  writeBin(c(2, -10), con, size = 4L, endian = .Platform$endian)
  close(con)
  expect_equal(RNifti::niftiHeader(path)$scl_slope, 2)
  expected <- matrix(2 * as.numeric(image) - 10, nrow = 8L)[c(8L, 2L), , drop = FALSE]
  expect_equal(medusa:::extract_nifti_voxel_ts(path, c(8L, 2L)), expected)

  image <- array(as.numeric(image), dim(image))
  image[c(2L, 8L, 16L)] <- c(NA_real_, NaN, Inf)
  RNifti::writeNifti(image, path)
  expected <- matrix(as.numeric(RNifti::readNifti(path)), nrow = 8L)[c(8L, 2L), , drop = FALSE]
  expect_equal(medusa:::extract_nifti_voxel_ts(path, c(8L, 2L)), expected)
})

test_that("native extraction validates indices and handles a single volume", {
  path <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(path))
  RNifti::writeNifti(array(seq_len(8L), c(2L, 2L, 2L)), path)
  expect_equal(medusa:::extract_nifti_voxel_ts(path, c(8L, 1L)), matrix(c(8, 1), ncol = 1L))
  expect_equal(medusa:::extract_nifti_voxel_ts(path, c(8L, 1L), TRUE), matrix(c(8, 1), nrow = 1L))
  for (bad in list(0, 9, -1, NA_real_, Inf, 1.5)) {
    expect_error(medusa:::extract_nifti_voxel_ts(path, bad), "one-based indices")
  }
})

test_that("time_first cache-blocked extraction matches RNifti across datatypes, scaling, and >64x64 tiles", {
  set.seed(20261113L)
  # 6 x 5 x 6 = 180 spatial voxels (> 2x64 tiles); 135 timepoints (> 2x64 tiles)
  dims <- c(6L, 5L, 6L, 135L)
  n_spatial <- prod(dims[1:3])
  n_time <- dims[4L]
  # Permuted subset of 145 voxels crossing tile boundaries at 64 and 128
  vox_idx <- sample(seq_len(n_spatial), size = 145L, replace = FALSE)

  for (dt in c("uint8", "int16", "int32", "float", "double")) {
    path <- tempfile(fileext = if (dt %in% c("int16", "float")) ".nii.gz" else ".nii")
    raw_vals <- if (dt == "uint8") {
      sample(0:250, n_spatial * n_time, replace = TRUE)
    } else if (dt %in% c("int16", "int32")) {
      sample(-500:500, n_spatial * n_time, replace = TRUE)
    } else {
      rnorm(n_spatial * n_time, mean = 100, sd = 15)
    }
    if (dt == "int32") {
      raw_vals[c(5L, 200L, 4000L)] <- NA_integer_
    }
    arr <- array(raw_vals, dim = dims)
    RNifti::writeNifti(arr, path, datatype = dt)

    if (dt == "int32" && grepl("\\.nii$", path)) {
      # Also test scaled int32 with NA_integer_ preserving NA_real_
      con <- file(path, "r+b")
      seek(con, 112L, rw = "write")
      writeBin(c(1.5, 10.25), con, size = 4L, endian = .Platform$endian)
      close(con)
    }

    ref_vt <- matrix(as.numeric(RNifti::readNifti(path)), nrow = n_spatial)[vox_idx, , drop = FALSE]
    ext_vt <- medusa:::extract_nifti_voxel_ts(path, vox_idx, FALSE)
    ext_tv <- medusa:::extract_nifti_voxel_ts(path, vox_idx, TRUE)

    expect_equal(dim(ext_vt), c(length(vox_idx), n_time))
    expect_equal(dim(ext_tv), c(n_time, length(vox_idx)))
    expect_equal(ext_vt, ref_vt, tolerance = 1e-7)
    expect_equal(ext_tv, t(ref_vt), tolerance = 1e-7)
    unlink(path)
  }
})

test_that("low-level C++ NIfTI writers validate arguments cleanly", {
  ref_path <- tempfile(fileext = ".nii.gz")
  out_path <- tempfile(fileext = ".nii.gz")
  on.exit(unlink(c(ref_path, out_path)), add = TRUE)
  RNifti::writeNifti(array(1L, c(2L, 2L, 2L)), ref_path)

  mat_tv <- matrix(as.numeric(1:12), nrow = 3L, ncol = 4L)
  expect_error(
    medusa:::write_nifti_indexed_cpp(ref_path, out_path, mat_tv, c(1, 2, 2, 4)),
    "duplicate"
  )
  expect_error(
    medusa:::write_nifti_indexed_cpp(ref_path, out_path, mat_tv, c(1, 2, 3, 9)),
    "one-based indices"
  )
  expect_error(
    medusa:::write_nifti_indexed_cpp(ref_path, out_path, mat_tv, c(1, 2, 3)),
    "must match the voxel dimension"
  )
  expect_error(
    medusa:::write_nifti_indexed_cpp(ref_path, out_path, mat_tv, c(1, 2, 3, 4), compression = 12L),
    "compression"
  )
  expect_error(
    medusa:::write_nifti_indexed_cpp(ref_path, out_path, mat_tv, c(1, 2, 3, 4), datatype = "int16"),
    "Unsupported output NIfTI datatype"
  )

  p_mat <- matrix(c(1, 2, 3, 4), nrow = 2L, ncol = 2L)
  expect_error(
    medusa:::write_nifti_parcel_cpp(ref_path, out_path, rep(1L, 7L), p_mat),
    "parcel_slot length"
  )
  expect_error(
    medusa:::write_nifti_parcel_cpp(ref_path, out_path, c(0L, 1L, 3L, rep(0L, 5L)), p_mat),
    "parcel_slot values"
  )
})
