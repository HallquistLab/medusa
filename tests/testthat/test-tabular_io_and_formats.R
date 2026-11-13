make_io_decon_fixture <- function() {
  root <- tempfile("io_decon_")
  dir.create(root)
  atlas <- array(0L, c(2L, 2L, 2L))
  atlas[c(1L, 2L, 3L, 4L, 5L, 7L)] <- c(1L, 1L, 2L, 2L, 1L, 2L)
  mask <- array(0, dim(atlas))
  mask[c(1L, 2L, 3L, 4L, 5L, 7L)] <- 1

  n_t <- 24L
  bold <- array(0, c(dim(atlas), n_t))
  t_idx <- 0:(n_t - 1L)
  for (i in c(1L, 2L, 3L, 4L)) {
    bold[i + 8L * t_idx] <- 100 + 2 * sin(t_idx * (i + 1) / 5) + cos(t_idx / (i + 2))
  }
  bold[5L + 8L * t_idx] <- 100 # constant voxel in parcel 1
  bold[7L + 8L * t_idx] <- NA_real_ # nonfinite voxel in parcel 2

  atlas_path <- file.path(root, "schaefer_mini.nii.gz")
  mask_path <- file.path(root, "mask.nii.gz")
  bold_path <- file.path(root, "sub-01_run-01_bold.nii.gz")
  RNifti::writeNifti(atlas, atlas_path)
  RNifti::writeNifti(mask, mask_path)
  RNifti::writeNifti(bold, bold_path)

  trials <- data.frame(
    id = "01",
    run = 1L,
    trial = 1:3,
    onset = c(4.0, 11.0, 18.0),
    stringsAsFactors = FALSE
  )
  alignment <- list(
    evt_col = "onset",
    time_before = -1,
    time_after = 2,
    pad_before = -1,
    pad_after = 1,
    output_resolution = 1
  )
  list(
    root = root,
    atlas = atlas_path,
    mask = mask_path,
    bold = bold_path,
    trials = trials,
    alignment = alignment
  )
}

test_that("voxelwise_deconvolution csv.gz, rds, and nifti outputs are numerically equivalent and align identically", {
  f <- make_io_decon_fixture()
  on.exit(unlink(f$root, recursive = TRUE), add = TRUE)

  dset_dbl <- list(kernel = c(0.1, 0.6, 0.3), lambda = 0.05, nifti_datatype = "double")
  dset_flt <- list(kernel = c(0.1, 0.6, 0.3), lambda = 0.05)
  meta <- data.frame(id = "01", run = 1L, session = "A", stringsAsFactors = FALSE)

  dir_csv <- file.path(f$root, "out_csv")
  dir_rds <- file.path(f$root, "out_rds")
  dir_nii_dbl <- file.path(f$root, "out_nii_dbl")
  dir_nii_flt <- file.path(f$root, "out_nii_flt")

  expect_warning(
    voxelwise_deconvolution(
      f$bold, add_metadata = meta, out_dir = dir_csv,
      out_file_expression = expression("sub01_run1"),
      TR = 1.0, time_offset = 0.5, atlas_files = f$atlas, mask = f$mask,
      nprocs = 1L, save_original_ts = TRUE, algorithm = "reglin",
      decon_settings = dset_dbl, output_format = "csv.gz"
    ),
    "2 constant or nonfinite"
  )
  expect_warning(
    voxelwise_deconvolution(
      f$bold, add_metadata = meta, out_dir = dir_rds,
      out_file_expression = expression("sub01_run1"),
      TR = 1.0, time_offset = 0.5, atlas_files = f$atlas, mask = f$mask,
      nprocs = 1L, save_original_ts = TRUE, algorithm = "reglin",
      decon_settings = dset_dbl, output_format = "rds"
    ),
    "2 constant or nonfinite"
  )
  expect_warning(
    voxelwise_deconvolution(
      f$bold, add_metadata = meta, out_dir = dir_nii_dbl,
      out_file_expression = expression("sub01_run1"),
      TR = 1.0, time_offset = 0.5, atlas_files = f$atlas, mask = f$mask,
      nprocs = 1L, save_original_ts = TRUE, algorithm = "reglin",
      decon_settings = dset_dbl, output_format = "nifti"
    ),
    "2 constant or nonfinite"
  )
  expect_warning(
    voxelwise_deconvolution(
      f$bold, add_metadata = meta, out_dir = dir_nii_flt,
      out_file_expression = expression("sub01_run1"),
      TR = 1.0, time_offset = 0.5, atlas_files = f$atlas, mask = f$mask,
      nprocs = 1L, save_original_ts = FALSE, algorithm = "reglin",
      decon_settings = dset_flt, output_format = "nifti"
    ),
    "2 constant or nonfinite"
  )

  csv_file <- file.path(dir_csv, "schaefer_mini", "deconvolved", "sub01_run1_deconvolved.csv.gz")
  rds_file <- file.path(dir_rds, "schaefer_mini", "deconvolved", "sub01_run1_deconvolved.rds")
  nii_dbl_file <- file.path(dir_nii_dbl, "schaefer_mini", "deconvolved", "sub01_run1_deconvolved.nii.gz")
  nii_flt_file <- file.path(dir_nii_flt, "schaefer_mini", "deconvolved", "sub01_run1_deconvolved.nii.gz")

  expect_true(all(file.exists(c(csv_file, rds_file, nii_dbl_file, paste0(nii_dbl_file, ".coords.rds")))))

  dt_csv <- data.table::fread(csv_file)
  expect_identical(
    names(dt_csv),
    c("vnum", "volume", "decon", "time", "atlas_value", "x", "y", "z", "atlas_name", "decon_status")
  )
  orig_csv <- data.table::fread(file.path(dir_csv, "schaefer_mini", "original", "sub01_run1_original.csv.gz"))
  expect_identical(
    names(orig_csv),
    c("vnum", "time", "BOLD_input", "atlas_value", "x", "y", "z", "atlas_name", "BOLD_z")
  )

  rds_obj <- readRDS(rds_file)
  expect_equal(as.vector(t(rds_obj$decon)), dt_csv$decon, tolerance = 1e-12)
  expect_equal(rds_obj$time, seq(0.5, by = 1.0, length.out = 24L))
  expect_equal(rds_obj$decon_status, c("ok", "ok", "ok", "ok", "constant", "nonfinite"))

  orig_rds <- readRDS(file.path(dir_rds, "schaefer_mini", "original", "sub01_run1_original.rds"))
  expect_equal(as.vector(t(orig_rds$BOLD_input)), orig_csv$BOLD_input, tolerance = 1e-12)
  expect_equal(as.vector(t(orig_rds$BOLD_z)), orig_csv$BOLD_z, tolerance = 1e-12)

  trials_sess <- cbind(f$trials, session = "A", stringsAsFactors = FALSE)
  al_csv <- evt_align_decon_files(csv_file, trials_sess, f$alignment, tr = 1.0, ncpus = 1L, time_audit = TRUE)
  al_rds <- evt_align_decon_files(rds_file, trials_sess, f$alignment, tr = 1.0, ncpus = 1L, time_audit = TRUE)
  al_nii_dbl <- evt_align_decon_files(nii_dbl_file, trials_sess, f$alignment, tr = 1.0, ncpus = 1L, time_audit = TRUE)
  al_nii_flt <- evt_align_decon_files(nii_flt_file, trials_sess, f$alignment, tr = 1.0, ncpus = 1L, time_audit = TRUE)

  cmp_cols <- setdiff(names(al_csv), "alignment_file")
  expect_equal(as.data.frame(al_rds[, ..cmp_cols]), as.data.frame(al_csv[, ..cmp_cols]), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(as.data.frame(al_nii_dbl[, ..cmp_cols]), as.data.frame(al_csv[, ..cmp_cols]), tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(as.data.frame(al_nii_flt[, ..cmp_cols]), as.data.frame(al_csv[, ..cmp_cols]), tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(attr(al_rds, "source_clock")$time, attr(al_csv, "source_clock")$time, tolerance = 1e-12)
  expect_equal(attr(al_nii_dbl, "source_clock")$time, attr(al_csv, "source_clock")$time, tolerance = 1e-12)

  # Verify run_decon_alignment discovers rds and nifti formats and writes aligned outputs
  run_decon_alignment(
    f$atlas, dir_rds, trials_sess, list(stim = f$alignment),
    tr = 1.0, ncpus = 1L, scheduler = "local", out_format = "rds"
  )
  rds_aligned_path <- file.path(dir_rds, "schaefer_mini", "schaefer_mini_stim_decon_aligned.rds")
  expect_true(file.exists(rds_aligned_path))
  expect_true(file.exists(paste0(rds_aligned_path, ".files.csv")))
  expect_true(file.exists(paste0(rds_aligned_path, ".sources.csv")))

  run_decon_alignment(
    f$atlas, dir_nii_dbl, trials_sess, list(stim = f$alignment),
    tr = 1.0, ncpus = 1L, scheduler = "local", out_format = "csv.gz"
  )
  nii_aligned_path <- file.path(dir_nii_dbl, "schaefer_mini", "schaefer_mini_stim_decon_aligned.csv.gz")
  expect_true(file.exists(nii_aligned_path))
  rds_aligned_dt <- readRDS(rds_aligned_path)
  cmp_cols_no_audit <- setdiff(names(rds_aligned_dt), "alignment_file")
  expect_equal(
    as.data.frame(data.table::fread(nii_aligned_path, colClasses = list(character = "id"))[, ..cmp_cols_no_audit]),
    as.data.frame(rds_aligned_dt[, ..cmp_cols_no_audit]),
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
})

test_that("discover_decon_files excludes sidecars and deduplicates multi-format stems", {
  adir <- tempfile("discover_")
  dir.create(adir)
  on.exit(unlink(adir, recursive = TRUE), add = TRUE)

  files_to_touch <- c(
    "sub1_run1_deconvolved.csv.gz",
    "sub1_run1_deconvolved.csv.gz.metadata.csv",
    "sub1_run1_deconvolved.csv.gz.scaling.csv",
    "sub1_run1_deconvolved.rds",
    "sub1_run1_deconvolved.rds.metadata.csv",
    "sub1_run2_deconvolved.nii.gz",
    "sub1_run2_deconvolved.nii.gz.coords.rds",
    "sub1_run2_deconvolved.nii.gz.metadata.csv"
  )
  file.create(file.path(adir, files_to_touch))

  expect_message(
    found_auto <- medusa:::discover_decon_files(adir, input_format = "auto"),
    "Multiple deconvolution file formats"
  )
  expect_setequal(basename(found_auto), c("sub1_run1_deconvolved.rds", "sub1_run2_deconvolved.nii.gz"))
  expect_equal(basename(medusa:::discover_decon_files(adir, input_format = "csv.gz")), "sub1_run1_deconvolved.csv.gz")
  expect_equal(basename(medusa:::discover_decon_files(adir, input_format = "nifti")), "sub1_run2_deconvolved.nii.gz")
})

test_that("Parquet I/O, column projection, and split_filter pushdown work via Arrow adapter", {
  root <- tempfile("pq_mock_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  # Fail-fast message when arrow is absent
  local_mocked_bindings(arrow_is_available = function() FALSE, .package = "medusa")
  expect_error(
    medusa:::medusa_write_table(data.frame(a = 1), file.path(root, "x.parquet")),
    "Install it with install.packages\\('arrow'\\)"
  )
  expect_error(
    mixed_by(
      data = {
        p <- file.path(root, "dummy.parquet")
        file.create(p)
        p
      },
      outcomes = "y",
      rhs_model_formulae = list(m1 = ~ x + (1 | id))
    ),
    "Package 'arrow' must be installed"
  )

  # Mock the thin Arrow primitives with RDS-backed columnar storage to record exact
  # projected columns and pushed-down filters even when arrow is not installed locally.
  pq_calls <- list()
  local_mocked_bindings(
    arrow_is_available = function() TRUE,
    arrow_write_parquet_impl = function(df, path) {
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(as.data.frame(df), path)
      invisible(path)
    },
    arrow_dataset_names = function(path) {
      if (dir.exists(path)) {
        pfiles <- list.files(path, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
        names(readRDS(pfiles[1L]))
      } else {
        names(readRDS(path))
      }
    },
    arrow_read_parquet_impl = function(path, col_select = NULL, split_filter = NULL, as_data_table = TRUE) {
      pq_calls[[length(pq_calls) + 1L]] <<- list(
        path = path,
        col_select = col_select,
        split_filter = split_filter
      )
      raw <- if (dir.exists(path)) {
        pfiles <- list.files(path, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
        do.call(rbind, lapply(pfiles, readRDS))
      } else {
        readRDS(path)
      }
      avail <- names(raw)
      cols_use <- if (!is.null(col_select)) intersect(unique(col_select), avail) else avail
      raw <- medusa:::apply_split_filter_dt(raw, split_filter, cols = avail)
      raw <- raw[, cols_use, drop = FALSE]
      if (isTRUE(as_data_table)) data.table::as.data.table(raw) else as.data.frame(raw)
    },
    .package = "medusa"
  )

  # Build synthetic stacked MEDuSA panel with many unused wide audit columns
  set.seed(42)
  panel <- expand.grid(
    id = sprintf("sub%02d", 1:12),
    trial = 1:8,
    atlas_value = 1:3,
    evt_time = c(-1, 0, 1),
    stringsAsFactors = FALSE
  )
  behav <- unique(panel[, c("id", "trial")])
  behav$pe <- rnorm(nrow(behav))
  behav$cond <- factor(rep(c("low", "high"), length.out = nrow(behav)))

  subj_re <- stats::setNames(rnorm(12, sd = 0.6), sprintf("sub%02d", 1:12))
  panel <- merge(panel, behav, by = c("id", "trial"))
  panel$decon_mean <- 0.4 * panel$pe + subj_re[panel$id] + rnorm(nrow(panel), sd = 0.5)
  panel$precision_weight <- runif(nrow(panel), 0.5, 2.0)
  for (k in 1:12) {
    panel[[paste0("unused_audit_col_", k)]] <- paste0("extra_string_payload_", k)
  }
  panel_neural <- panel[, setdiff(names(panel), c("pe", "cond"))]

  pq_panel_path <- file.path(root, "aligned_panel.parquet")
  csvgz_panel_path <- file.path(root, "aligned_panel.csv.gz")
  pq_behav_path <- file.path(root, "behav.parquet")
  csv_behav_path <- file.path(root, "behav.csv.gz")

  medusa:::medusa_write_table(panel_neural, pq_panel_path)
  medusa:::medusa_write_table(panel_neural, csvgz_panel_path)
  medusa:::medusa_write_table(behav, pq_behav_path)
  medusa:::medusa_write_table(behav, csv_behav_path)

  # Partitioned Parquet directory support
  pq_dir <- file.path(root, "partitioned_panel")
  for (av in 1:3) {
    part_dir <- file.path(pq_dir, paste0("atlas_value=", av))
    dir.create(part_dir, recursive = TRUE)
    saveRDS(panel_neural[panel_neural$atlas_value == av, ], file.path(part_dir, "part-0.parquet"))
  }

  s_filt <- list(atlas_value = c(1L, 3L), evt_time = c(0, 1))
  ref_sub <- panel[panel$atlas_value %in% c(1L, 3L) & panel$evt_time %in% c(0, 1), ]

  fit_mem <- mixed_by(
    data = ref_sub,
    outcomes = "decon_mean",
    rhs_model_formulae = list(m1 = ~ pe + cond + (1 | id)),
    split_on = c("atlas_value", "evt_time"),
    weights = "precision_weight",
    calculate = c("parameter_estimates_reml", "fit_statistics"),
    emmeans_spec = list(cond_emm = list(outcome = "decon_mean", model_name = "m1", specs = ~ cond))
  )

  pq_calls <- list()
  fit_pq <- mixed_by(
    data = pq_panel_path,
    outcomes = "decon_mean",
    rhs_model_formulae = list(m1 = ~ pe + cond + (1 | id)),
    split_on = c("atlas_value", "evt_time"),
    external_df = pq_behav_path,
    external_merge_by = c("id", "trial"),
    weights = "precision_weight",
    split_filter = s_filt,
    calculate = c("parameter_estimates_reml", "fit_statistics"),
    emmeans_spec = list(cond_emm = list(outcome = "decon_mean", model_name = "m1", specs = ~ cond))
  )

  # Verify column projection excluded all 12 unused audit columns and pushed down split_filter
  panel_scan_call <- pq_calls[[2L]]
  expect_false(any(grepl("^unused_audit_col_", panel_scan_call$col_select)))
  expect_equal(panel_scan_call$split_filter, s_filt)

  cmp_reml <- intersect(names(fit_mem$coef_df_reml), names(fit_pq$coef_df_reml))
  cmp_reml <- setdiff(cmp_reml, "padj_BY_term")
  expect_equal(
    as.data.frame(fit_pq$coef_df_reml[, ..cmp_reml]),
    as.data.frame(fit_mem$coef_df_reml[, ..cmp_reml]),
    tolerance = 1e-8,
    ignore_attr = TRUE
  )
  expect_equal(
    fit_pq$emmeans_list$cond_emm$emmean,
    fit_mem$emmeans_list$cond_emm$emmean,
    tolerance = 1e-8
  )

  # Partitioned Parquet directory + CSV external_df must also match identically
  fit_pq_dir <- mixed_by(
    data = pq_dir,
    outcomes = "decon_mean",
    rhs_model_formulae = list(m1 = ~ pe + cond + (1 | id)),
    split_on = c("atlas_value", "evt_time"),
    external_df = csv_behav_path,
    external_merge_by = c("id", "trial"),
    weights = "precision_weight",
    split_filter = s_filt,
    calculate = c("parameter_estimates_reml", "fit_statistics")
  )
  expect_equal(
    as.data.frame(fit_pq_dir$coef_df_reml[, ..cmp_reml]),
    as.data.frame(fit_mem$coef_df_reml[, ..cmp_reml]),
    tolerance = 1e-7,
    ignore_attr = TRUE
  )

  # Delimited .csv.gz column projection + split_filter must also match identically
  fit_csvgz <- mixed_by(
    data = csvgz_panel_path,
    outcomes = "decon_mean",
    rhs_model_formulae = list(m1 = ~ pe + cond + (1 | id)),
    split_on = c("atlas_value", "evt_time"),
    external_df = csv_behav_path,
    external_merge_by = c("id", "trial"),
    weights = "precision_weight",
    split_filter = s_filt,
    calculate = c("parameter_estimates_reml", "fit_statistics")
  )
  expect_equal(
    as.data.frame(fit_csvgz$coef_df_reml[, ..cmp_reml]),
    as.data.frame(fit_mem$coef_df_reml[, ..cmp_reml]),
    tolerance = 1e-7,
    ignore_attr = TRUE
  )

  # Typo in split_filter column name fails loudly
  expect_error(
    mixed_by(
      data = ref_sub,
      outcomes = "decon_mean",
      rhs_model_formulae = list(m1 = ~ pe + (1 | id)),
      split_on = c("atlas_value", "evt_time"),
      split_filter = list(nonexistent_col = 1L)
    ),
    "nonexistent_col"
  )
})

test_that("real arrow Parquet roundtrip and Hive dataset projection work when arrow is installed", {
  skip_if_not_installed("arrow")
  f <- make_io_decon_fixture()
  on.exit(unlink(f$root, recursive = TRUE), add = TRUE)

  expect_warning(
    voxelwise_deconvolution(
      f$bold, out_dir = f$root, out_file_expression = expression("sub01_run1"),
      TR = 1.0, atlas_files = f$atlas, mask = f$mask, nprocs = 1L,
      save_original_ts = FALSE, algorithm = "reglin",
      decon_settings = list(kernel = c(0.1, 0.6, 0.3), lambda = 0.05),
      output_format = "parquet"
    ),
    "2 constant or nonfinite"
  )

  pq_decon <- file.path(f$root, "schaefer_mini", "deconvolved", "sub01_run1_deconvolved.parquet")
  expect_true(file.exists(pq_decon))

  pq_aligned <- file.path(f$root, "aligned.parquet")
  res <- evt_align_decon_files(
    pq_decon, f$trials, f$alignment,
    tr = 1.0, out_file = pq_aligned, ncpus = 1L,
    file_manifest = data.frame(file = pq_decon, id = "01", run = 1L)
  )
  expect_true(file.exists(pq_aligned))
  rt <- medusa:::medusa_read_table(pq_aligned, col_select = c("id", "trial", "atlas_value", "evt_time", "decon_mean"))
  expect_equal(ncol(rt), 5L)
  expect_equal(rt$decon_mean, res$decon_mean, tolerance = 1e-12)
})
