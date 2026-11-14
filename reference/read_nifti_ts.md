# Read 3D/4D NIfTI runs directly into MEDuSA time-by-voxel matrices

Validates NIfTI spatial headers (\`dim\`, \`pixdim\`, units, and
\`sform\`/\`qform\` affine) across runs, optional mask, and optional
parcellation atlas before extracting selected voxels in C++ directly
into \`time x voxels\` (\`T x V\`) numeric matrices without allocating
full 4D double-precision R arrays.

The returned \`"medusa_nifti_ts"\` object supplies \`\$bold\`,
\`\$voxel_map\`, \`\$parcel_ids\`, and \`\$tr\` ready to pass directly
to \[estimate_crossrun_hrf()\], \[deconvolve_crossrun_hrf()\],
\[deconvolve_reglin()\], and \[write_nifti_ts()\].

## Usage

``` r
read_nifti_ts(
  niftis,
  mask = NULL,
  atlas = NULL,
  voxel_indices = NULL,
  atlas_subset = NULL,
  mask_threshold = 0,
  atlas_threshold = 0,
  tr = NULL,
  run_names = NULL,
  voxel_ids = NULL,
  id_style = c("ijk", "vnum"),
  orientation = c("time_by_voxel", "voxel_by_time"),
  drop = FALSE,
  ncpus = 1L
)
```

## Arguments

- niftis:

  Character vector of paths to 3D or 4D scalar NIfTI files (\`.nii\` or
  \`.nii.gz\`). Named vectors preserve run names in \`\$bold\`.

- mask:

  Optional path to a single-volume 3D NIfTI mask (or an in-memory
  \`RNifti\` image/array) retaining voxels strictly greater than
  \`mask_threshold\`.

- atlas:

  Optional path to a single-volume 3D NIfTI parcellation/atlas (or
  in-memory \`RNifti\` image/array) retaining voxels strictly greater
  than \`atlas_threshold\` and populating \`\$voxel_map\$atlas_value\`
  and \`\$parcel_ids\`.

- voxel_indices:

  Optional numeric/integer vector of unique 1-based linear spatial
  indices to extract in exact order. When supplied, \`mask\` and
  \`atlas\` do not filter indices, though their values at
  \`voxel_indices\` are still recorded in \`\$voxel_map\` if provided.

- atlas_subset:

  Optional numeric/integer vector of \`atlas\` parcel values to retain
  when \`voxel_indices\` is \`NULL\` and \`atlas\` is supplied.

- mask_threshold:

  Numeric threshold (\`\>\` comparison) applied to \`mask\` when
  \`voxel_indices\` is \`NULL\`. Default \`0\`.

- atlas_threshold:

  Numeric threshold (\`\>\` comparison) applied to \`atlas\` when
  \`voxel_indices\` is \`NULL\`. Default \`0\`.

- tr:

  Optional positive scalar repetition time in seconds. If \`NULL\`,
  inferred from 4D NIfTI headers (\`pixdim\[5\]\`, converting
  \`msec\`/\`usec\` units to seconds) when positive and consistent
  across runs.

- run_names:

  Optional character vector of unique run names matching
  \`length(niftis)\`. Overrides \`names(niftis)\`. Defaults to
  \`"run1"\`, \`"run2"\`, etc. when \`niftis\` is unnamed.

- voxel_ids:

  Optional character vector of unique voxel identifiers of length equal
  to the selected voxel count.

- id_style:

  Naming convention when \`voxel_ids\` is \`NULL\`: \`"ijk"\` (default,
  1-based \`"i_j_k"\` grid coordinates) or \`"vnum"\` (\`"v1"\`,
  \`"v2"\`, ...).

- orientation:

  Matrix orientation: \`"time_by_voxel"\` (default, \`T x V\` matching
  all MEDuSA solvers) or \`"voxel_by_time"\` (\`V x T\`).

- drop:

  Logical; if \`TRUE\` and \`length(niftis) == 1L\`, returns the single
  extracted matrix directly with \`voxel_map\`, \`parcel_ids\`, \`tr\`,
  and \`geometry\` attached as attributes. Default \`FALSE\`.

- ncpus:

  Integer number of parallel worker processes used when extracting
  multiple NIfTI runs (\`length(niftis) \> 1L\`). Default \`1L\`.

## Value

By default (\`drop = FALSE\`), an S3 list of class \`"medusa_nifti_ts"\`
with elements:

- \`bold\`: Named list of numeric matrices (\`T x V\` by default) with
  voxel IDs on the voxel dimension.

- \`voxel_map\`: \`data.frame\` with columns \`voxel\`, \`vnum\`,
  \`linear_index\`, \`i\`, \`j\`, \`k\`, \`x\`, \`y\`, \`z\`, and
  optional \`mask_value\` and \`atlas_value\`.

- \`parcel_ids\`: Named vector mapping \`voxel\` IDs to \`atlas_value\`
  when \`atlas\` is supplied (\`NULL\` otherwise), directly usable as
  \`spatial_groups\` or \`parcel_groups\`.

- \`tr\`: Scalar repetition time in seconds (\`NA_real\_\` if
  unavailable).

- \`geometry\`: Spatial geometry metadata list (\`dim\`, \`spacing\`,
  \`units\`, \`affine\`, \`reference_file\`).

- \`orientation\`: Selected matrix orientation string.

- \`run_info\`: \`data.frame\` summarizing \`run\`, \`file\`,
  \`n_time\`, and \`header_tr\`.

## See also

\[write_nifti_ts()\], \[align_medusa_ts()\],
\[estimate_crossrun_hrf()\], \[deconvolve_crossrun_hrf()\],
\[decon_epoch_regression()\]

## Examples

``` r
img <- array(rnorm(4 * 4 * 2 * 10), dim = c(4, 4, 2, 10))
mask <- array(0L, dim = c(4, 4, 2))
mask[1:2, 1:2, 1] <- 1L
f_img <- tempfile(fileext = ".nii.gz")
f_mask <- tempfile(fileext = ".nii.gz")
RNifti::writeNifti(img, f_img)
RNifti::writeNifti(mask, f_mask)

ts <- read_nifti_ts(c(run1 = f_img), mask = f_mask, tr = 1.0)
dim(ts$bold$run1)
#> [1] 10  4
head(ts$voxel_map)
#>   voxel vnum linear_index i j k x y z mask_value
#> 1 1_1_1    1            1 1 1 1 0 0 0          1
#> 2 2_1_1    2            2 2 1 1 1 0 0          1
#> 3 1_2_1    3            5 1 2 1 0 1 0          1
#> 4 2_2_1    4            6 2 2 1 1 1 0          1
unlink(c(f_img, f_mask))
```
