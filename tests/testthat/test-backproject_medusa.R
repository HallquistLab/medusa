test_that("backproject_medusa writes 4D RNifti volumes with preserved geometry and TR", {
  out_dir <- tempfile("backproj_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  mask_arr <- array(0L, dim = c(3L, 2L, 2L))
  mask_arr[1:2, 1, 1] <- 1L
  mask_arr[3, 2, 2] <- 2L
  mask_path <- file.path(out_dir, "mask.nii.gz")
  affine <- diag(c(2, 2, 2, 1))
  affine[1:3, 4] <- c(-10, -20, 30)
  attr(affine, "code") <- 1L
  mask_nii <- RNifti::asNifti(mask_arr)
  RNifti::sform(mask_nii) <- affine
  RNifti::writeNifti(mask_nii, mask_path)

  coef_df <- expand.grid(
    roi = c(1L, 2L),
    time = c(-1.5, 0.0, 1.5),
    term = c("Intercept", "rpe"),
    stringsAsFactors = FALSE
  )
  coef_df$estimate <- seq_len(nrow(coef_df)) * 0.5
  coef_df$statistic <- -seq_len(nrow(coef_df))

  written <- backproject_medusa(
    coef_df = coef_df,
    brain_mask = mask_path,
    plot_cols = c("estimate", "statistic"),
    parcel_col = "roi",
    time_col = "time",
    effect_col = "term",
    output_dir = out_dir,
    ncpus = 1L
  )

  expect_length(written, 4L)
  expect_true(all(file.exists(written)))

  est_rpe_file <- file.path(out_dir, "estimate_rpe.nii.gz")
  expect_true(est_rpe_file %in% written)
  hdr <- RNifti::niftiHeader(est_rpe_file)
  expect_equal(as.integer(hdr$dim[1:5]), c(4L, 3L, 2L, 2L, 3L))
  expect_equal(hdr$pixdim[5], 1.5, tolerance = 1e-6)
  expect_equal(matrix(as.numeric(RNifti::xform(hdr)), 4L, 4L), matrix(as.numeric(affine), 4L, 4L), tolerance = 1e-6)

  img <- RNifti::readNifti(est_rpe_file)
  rpe_sub <- coef_df[coef_df$term == "rpe", ]
  expect_equal(img[1, 1, 1, ], rpe_sub$estimate[rpe_sub$roi == 1L][order(rpe_sub$time[rpe_sub$roi == 1L])], tolerance = 1e-6)
  expect_equal(img[2, 1, 1, ], img[1, 1, 1, ], tolerance = 1e-6)
  expect_equal(img[3, 2, 2, ], rpe_sub$estimate[rpe_sub$roi == 2L][order(rpe_sub$time[rpe_sub$roi == 2L])], tolerance = 1e-6)
  expect_equal(img[3, 1, 1, ], c(0, 0, 0))

  # Duplicate parcel-time-effect rows fail cleanly instead of entering browser()
  dup_df <- rbind(coef_df, coef_df[1L, , drop = FALSE])
  expect_error(
    backproject_medusa(dup_df, mask_path, plot_cols = "estimate", output_dir = out_dir, ncpus = 1L),
    "Expected exactly 1 row"
  )

  # Missing parcel-time-effect cell also fails with exact message
  miss_df <- coef_df[-2L, , drop = FALSE]
  expect_error(
    backproject_medusa(miss_df, mask_path, plot_cols = "estimate", output_dir = out_dir, ncpus = 1L),
    "Expected exactly 1 row.*found 0"
  )

  # Parallel execution (ncpus = 2L), fill = NA_real_, compression = 1L, effect_col = NULL, and single-timepoint 3D output
  out_par <- file.path(out_dir, "par")
  written_par <- backproject_medusa(
    coef_df = coef_df,
    brain_mask = mask_path,
    plot_cols = c("estimate", "statistic"),
    output_dir = out_par,
    ncpus = 2L,
    compression = 1L,
    fill = NA_real_
  )
  expect_length(written_par, 4L)
  img_par <- RNifti::readNifti(file.path(out_par, "estimate_rpe.nii.gz"))
  expect_equal(img_par[1, 1, 1, ], img[1, 1, 1, ], tolerance = 1e-6)
  expect_true(all(is.na(img_par[3, 1, 1, ])))

  single_t_df <- data.frame(roi = c(1L, 2L), time = c(0, 0), z = c(2.5, -3.25))
  out_3d <- file.path(out_dir, "single_t")
  w_3d <- backproject_medusa(
    coef_df = single_t_df,
    brain_mask = mask_path,
    plot_cols = "z",
    effect_col = NULL,
    output_dir = out_3d,
    ncpus = 1L
  )
  hdr_3d <- RNifti::niftiHeader(w_3d)
  expect_equal(as.integer(hdr_3d$dim[1:4]), c(3L, 3L, 2L, 2L))
  img_3d <- RNifti::readNifti(w_3d)
  expect_equal(img_3d[1, 1, 1], 2.5, tolerance = 1e-6)
  expect_equal(img_3d[3, 2, 2], -3.25, tolerance = 1e-6)
})

test_that("backproject_medusa coerces character/factor numeric parcel IDs and auto-detects parcel/evt_time/lag columns", {
  out_dir <- withr::local_tempdir()
  mask_arr <- array(0L, dim = c(3L, 2L, 2L))
  mask_arr[1:2, 1, 1] <- 10L
  mask_arr[3, 2, 2] <- 20L
  mask_path <- file.path(out_dir, "atlas.nii.gz")
  RNifti::writeNifti(RNifti::asNifti(mask_arr, pixdim = c(1, 2, 2, 2)), mask_path)

  # 1. mixed_by()-style output with character 'parcel' and 'evt_time' (no manual as.integer or parcel_col needed)
  mb_df <- data.frame(
    parcel = c("10", "20", "10", "20"),
    evt_time = c(0.0, 0.0, 1.5, 1.5),
    term = "rpe",
    statistic = c(3.2, -1.8, 4.1, -0.5),
    stringsAsFactors = FALSE
  )
  out_mb <- file.path(out_dir, "mb")
  f_mb <- backproject_medusa(mb_df, brain_mask = mask_path, plot_cols = "statistic", output_dir = out_mb, ncpus = 1L)
  expect_length(f_mb, 1L)
  img_mb <- RNifti::readNifti(f_mb)
  expect_equal(as.numeric(img_mb[1, 1, 1, ]), c(3.2, 4.1), tolerance = 1e-6)
  expect_equal(as.numeric(img_mb[3, 2, 2, ]), c(-1.8, -0.5), tolerance = 1e-6)
  # Input data.frame parcel column remains untouched as character
  expect_type(mb_df$parcel, "character")

  # 2. decon_epoch_regression()-style output with factor 'parcel', 'lag', and all-NA 'term' (marginal tau2_bio)
  ep_df <- data.frame(
    parcel = factor(c("10", "20")),
    lag = c(0.0, 0.0),
    term = c(NA_character_, NA_character_),
    tau2_bio = c(0.75, 1.40)
  )
  out_ep <- file.path(out_dir, "ep")
  f_ep <- backproject_medusa(ep_df, brain_mask = mask_path, plot_cols = "tau2_bio", output_dir = out_ep, ncpus = 1L)
  expect_equal(basename(f_ep), "tau2_bio.nii.gz")
  img_ep <- RNifti::readNifti(f_ep)
  expect_equal(img_ep[1, 1, 1], 0.75, tolerance = 1e-6)
  expect_equal(img_ep[3, 2, 2], 1.40, tolerance = 1e-6)

  # 3. Non-numeric string parcel labels fail with a clear diagnostic message
  bad_df <- data.frame(parcel = c("L_V1", "R_V1"), evt_time = c(0, 0), term = "rpe", statistic = c(1.0, 2.0))
  expect_error(
    backproject_medusa(bad_df, brain_mask = mask_path, plot_cols = "statistic", output_dir = out_dir, ncpus = 1L),
    "non-numeric label.*L_V1"
  )
})
