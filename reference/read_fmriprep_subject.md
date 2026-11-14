# Discover and ingest an fMRIPrep or postprocessed participant into a unified MEDuSA time-series bundle

Locates a participant's preprocessed or `BrainGnomes`-postprocessed BOLD
runs (`*_desc-preproc_bold.nii.gz` or `*_desc-<desc>_bold.nii.gz`),
brain or template masks (`*_desc-brain_mask.nii.gz`,
`*_templatemask.nii.gz`), confound tables
(`*_desc-confounds_timeseries.tsv`, `*_desc-<desc>_confounds.tsv`),
optional `*_censor.1D` masks, and optional BIDS event tables
(`*_events.tsv`), returning a `"medusa_nifti_ts"` object populated with
C++-extracted `T x V` `$bold` matrices, `$voxel_map`, `$parcel_ids`,
`$tr`, `$nuisance`, `$valid_rows`, `$confounds`, and `$events` in locked
run order, ready to pass directly into
[`estimate_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md)
and
[`deconvolve_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/deconvolve_crossrun_hrf.md).

## Usage

``` r
read_fmriprep_subject(
  fmriprep_dir,
  sub,
  bids_dir = NULL,
  ses = NULL,
  task = NULL,
  runs = NULL,
  space = NULL,
  cohort = NULL,
  res = NULL,
  desc = "preproc",
  mask = "auto",
  mask_mode = c("intersection", "union", "first"),
  atlas = NULL,
  atlas_subset = NULL,
  mask_threshold = 0,
  atlas_threshold = 0,
  tr = NULL,
  events = "auto",
  condition_col = "trial_type",
  condition_map = NULL,
  default_condition = NULL,
  exclude_conditions = NULL,
  na_duration = 0,
  strategy = c("24p", "6p", "12p", "24p_csf_wm", "27p", "36p", "compcor", "custom",
    "none"),
  columns = NULL,
  column_pattern = NULL,
  n_compcor = NULL,
  compcor_type = c("anatomical", "split", "temporal"),
  include_cosines = "auto",
  censor = "auto",
  fd_threshold = 0.5,
  fd_col = "auto",
  dvars_threshold = NULL,
  dvars_col = "std_dvars",
  include_non_steady_state = TRUE,
  include_motion_outliers = FALSE,
  drop_initial_volumes = 0L,
  censor_prev = 0L,
  censor_next = 0L,
  min_segment_length = 1L,
  id_style = c("ijk", "vnum"),
  ncpus = 1L
)
```

## Arguments

- fmriprep_dir:

  Path to the fMRIPrep or `BrainGnomes` postprocessing derivatives root
  directory (or directly to `sub-<id>`).

- sub:

  Participant identifier (with or without `"sub-"` prefix, e.g., `"01"`
  or `"sub-01"`).

- bids_dir:

  Optional path to the raw BIDS root directory containing `*_events.tsv`
  files (when `events = "auto"`).

- ses:

  Optional session filter (with or without `"ses-"` prefix).

- task:

  Optional BIDS task label (with or without `"task-"` prefix). Required
  when multiple tasks exist for the participant.

- runs:

  Optional vector of run identifiers to include (e.g., `1:3` or
  `c("01", "02")`).

- space:

  Optional template/output space filter (e.g., `"MNI152NLin2009cAsym"`).
  Auto-selected when only one space is present; required when multiple
  spaces exist.

- cohort:

  Optional TemplateFlow `cohort-*` entity filter (e.g., `"1"` or `"2"`
  for `MNIPediatricAsym`). Auto-selected when unique; required when
  multiple cohorts exist.

- res:

  Optional BIDS resolution entity (`res-*`). Auto-selected when unique.

- desc:

  BOLD derivative `desc-*` entity (default `"preproc"` for fMRIPrep; use
  e.g. `"postproc"` for `BrainGnomes` postprocessed outputs, or `NULL`
  to match un-described BOLD files).

- mask:

  Either `"auto"` (default; discovers per-run `*_desc-brain_mask.nii.gz`
  or `*_templatemask.nii.gz` files and combines them via `mask_mode`,
  falling back to anatomical mask, provenance spatial mask, or `atlas`),
  an explicit 3D NIfTI path / `RNifti` image, or `NULL`.

- mask_mode:

  Rule for combining per-run brain masks when `mask = "auto"`:
  `"intersection"` (default; retains voxels in-mask across every run),
  `"union"`, or `"first"`.

- atlas:

  Optional 3D parcellation NIfTI path or `RNifti` image on the BOLD
  spatial grid.

- atlas_subset:

  Optional integer/numeric vector of parcel labels to retain.

- mask_threshold:

  Numeric threshold (`>`) applied to `mask`. Default `0`.

- atlas_threshold:

  Numeric threshold (`>`) applied to `atlas`. Default `0`.

- tr:

  Optional repetition time in seconds (inferred from NIfTI headers or
  `*_bold.json` `RepetitionTime` when `NULL`).

- events:

  Either `"auto"` (default; discovers matching `*_events.tsv` in
  `bids_dir` / `fmriprep_dir`), explicit event file paths /
  `data.frame`s passed to
  [`read_bids_events`](https://hallquistlab.github.io/medusa/reference/read_bids_events.md),
  or `NULL` to skip event ingestion (e.g., resting-state).

- condition_col, condition_map, default_condition, exclude_conditions,
  na_duration:

  Event normalization options forwarded to
  [`read_bids_events`](https://hallquistlab.github.io/medusa/reference/read_bids_events.md).

- strategy, columns, column_pattern, n_compcor, compcor_type,
  include_cosines:

  Nuisance selection options forwarded to
  [`read_fmriprep_confounds`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_confounds.md).
  When `BrainGnomes` provenance indicates `confound_regression` was
  already applied to the BOLD derivative and `strategy` is not
  explicitly specified, `strategy` defaults automatically to `"none"` to
  avoid double-regressing unfiltered confounds.

- censor:

  Either `"auto"` (default; automatically uses matching `*_censor.1D`
  sidecars when present for postprocessed outputs), explicit censor file
  paths / vectors forwarded to
  [`read_fmriprep_confounds`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_confounds.md),
  or `NULL` to ignore `*_censor.1D` sidecars.

- fd_threshold, fd_col, dvars_threshold, dvars_col,
  include_non_steady_state, include_motion_outliers,
  drop_initial_volumes, censor_prev, censor_next, min_segment_length:

  Clock-preserving censoring options forwarded to
  [`read_fmriprep_confounds`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_confounds.md).

- id_style:

  Voxel ID convention forwarded to
  [`read_nifti_ts`](https://hallquistlab.github.io/medusa/reference/read_nifti_ts.md)
  (`"ijk"` or `"vnum"`).

- ncpus:

  Number of parallel workers forwarded to
  [`read_nifti_ts`](https://hallquistlab.github.io/medusa/reference/read_nifti_ts.md).

## Value

An S3 object of class `"medusa_nifti_ts"` enriched with `$nuisance`,
`$valid_rows`, `$confounds` (`"medusa_confounds"`), `$events` (when
discovered/supplied), and `$bids` provenance metadata.

## See also

[`read_nifti_ts`](https://hallquistlab.github.io/medusa/reference/read_nifti_ts.md),
[`read_fmriprep_confounds`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_confounds.md),
[`read_bids_events`](https://hallquistlab.github.io/medusa/reference/read_bids_events.md),
[`estimate_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md),
[`deconvolve_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/deconvolve_crossrun_hrf.md)

## Examples

``` r
root <- tempfile("medusa_bids_demo_")
fprep_func <- file.path(root, "derivatives", "fmriprep", "sub-01", "func")
bids_func <- file.path(root, "bids", "sub-01", "func")
dir.create(fprep_func, recursive = TRUE)
dir.create(bids_func, recursive = TRUE)

mask_img <- RNifti::asNifti(array(1L, dim = c(2L, 2L, 2L)))
RNifti::pixdim(mask_img) <- c(2, 2, 2)
for (r in 1:2) {
  stem <- sprintf("sub-01_task-mid_run-0%d", r)
  bold_arr <- array(rnorm(2L * 2L * 2L * 10L), dim = c(2L, 2L, 2L, 10L))
  bold_img <- RNifti::asNifti(bold_arr)
  RNifti::pixdim(bold_img) <- c(2, 2, 2, 1.0)
  RNifti::writeNifti(
    bold_img,
    file.path(fprep_func, paste0(stem, "_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"))
  )
  RNifti::writeNifti(
    mask_img,
    file.path(fprep_func, paste0(stem, "_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz"))
  )
  conf_df <- data.frame(
    trans_x = rnorm(10, sd = 0.02), trans_y = rnorm(10, sd = 0.02),
    trans_z = rnorm(10, sd = 0.02), rot_x = rnorm(10, sd = 0.002),
    rot_y = rnorm(10, sd = 0.002), rot_z = rnorm(10, sd = 0.002),
    framewise_displacement = c("n/a", rep("0.08", 9))
  )
  utils::write.table(
    conf_df,
    file.path(fprep_func, paste0(stem, "_desc-confounds_timeseries.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
  ev_df <- data.frame(onset = c(2, 6), duration = c(0.5, 0.5), trial_type = c("cue", "feedback"))
  utils::write.table(
    ev_df,
    file.path(bids_func, paste0(stem, "_events.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

sub01 <- read_fmriprep_subject(
  fmriprep_dir = file.path(root, "derivatives", "fmriprep"),
  bids_dir = file.path(root, "bids"),
  sub = "01",
  task = "mid",
  strategy = "6p"
)
sub01
#> <medusa_nifti_ts> 2 runs | 8 voxels | spatial dim: 2x2x2 | TR: 1 s | orientation: time_by_voxel
#>   Confounds attached: 0/20 volume(s) censored | nuisance: 6 col(s)
#>   BIDS events attached: 4 event row(s) ($events available)
#>   - run1: 10 volumes (sub-01_task-mid_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz)
#>   - run2: 10 volumes (sub-01_task-mid_run-02_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz)
unlink(root, recursive = TRUE)
```
