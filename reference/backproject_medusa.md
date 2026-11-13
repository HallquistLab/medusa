# Backproject parcel-level MEDuSA coefficients onto a 4D NIfTI brain mask

Maps parcel-by-time model statistics (for example from \[mixed_by()\])
back into 4D NIfTI volumes (3D space across event-locked timepoints)
using an integer-valued 3D parcellation mask.

## Usage

``` r
backproject_medusa(
  coef_df,
  brain_mask,
  plot_cols = NULL,
  parcel_col = "roi",
  time_col = "time",
  effect_col = "term",
  output_dir = getwd(),
  ncpus = 4L
)
```

## Arguments

- coef_df:

  A \`data.frame\` containing parcel, time, optional effect, and
  statistic columns to write.

- brain_mask:

  Path to a 3D integer-valued NIfTI parcellation mask.

- plot_cols:

  Character vector of numeric column names in \`coef_df\` to
  backproject.

- parcel_col:

  Column name in \`coef_df\` identifying integer mask parcel values.
  Default \`"roi"\`.

- time_col:

  Column name in \`coef_df\` identifying event-locked timepoints.
  Default \`"time"\`.

- effect_col:

  Column name in \`coef_df\` identifying model terms/effects, or
  \`NULL\` if \`coef_df\` contains a single effect. Default \`"term"\`.

- output_dir:

  Output directory where \`.nii.gz\` files will be written. Default
  \`getwd()\`.

- ncpus:

  Integer number of parallel worker processes. Default \`4L\`.

## Value

Invisibly returns a character vector of written \`.nii.gz\` file paths.
