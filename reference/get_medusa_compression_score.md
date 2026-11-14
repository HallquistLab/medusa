# Compute event-locked multivariate PCA compression scores across spatial units

Compute event-locked multivariate PCA compression scores across spatial
units

## Usage

``` r
get_medusa_compression_score(
  fmri_obj,
  event = NULL,
  time_before = -3,
  time_after = 3,
  collide_before = NULL,
  collide_after = NULL,
  group_by = NULL,
  events = NULL,
  tr = NULL,
  time = NULL,
  time_offset = NULL,
  parcel_ids = NULL,
  unit_data = NULL,
  value_name = "value",
  trial_col = "trial",
  pexp_target = 0.9,
  scale_columns = TRUE,
  use_legacy = FALSE
)
```

## Arguments

- fmri_obj:

  an fmri_ts object, numeric `T x V` matrix, or single-run
  `"medusa_nifti_ts"` object

- event:

  the event to which the fmri time series should be aligned (column in
  `fmri_obj$event_data` or `events`)

- time_before:

  How many seconds before the `event` do we want data. Default: -3

- time_after:

  How many seconds after the `event` do we want data. Default: 3

- collide_before:

  An optional vector of column names in `trial_df` that set a boundary
  on the earliest time point used in interpolation. This effectively
  truncates the time series for interpolation to a smaller window than
  specified by `time_before`.

- collide_after:

  An optional vector of column names in `trial_df` that set a boundary
  on the latest time point used in interpolation. This effectively
  truncates the time series for interpolation to a smaller window than
  specified by `time_after`.

- group_by:

  return compression scores for each combination of group_by variables.

- events:

  Optional trial-level `data.frame` for matrix inputs.

- tr:

  Optional repetition time in seconds for matrix inputs.

- time:

  Optional numeric vector of acquisition timestamps for matrix inputs.

- time_offset:

  Optional scalar start-time offset in seconds for matrix inputs.

- parcel_ids:

  Optional vector assigning matrix columns to parcels/ROIs.

- unit_data:

  Optional `data.frame` of unit metadata for matrix inputs.

- value_name:

  Signal name prefix for single-signal matrix inputs. Default `"value"`.

- trial_col:

  Trial identifier column name for matrix inputs. Default `"trial"`.

- pexp_target:

  Target cumulative variance proportion(s) passed to
  [`compress_mts_pca()`](https://hallquistlab.github.io/medusa/reference/compress_mts_pca.md).
  Default `0.9`.

- scale_columns:

  Logical; standardize columns prior to SVD. Default `TRUE`.

- use_legacy:

  Logical; force legacy `fmri_ts` path when `TRUE`.

## Author

Michael Hallquist
