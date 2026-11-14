# Extract selected voxel time series in a specified order

Extract selected voxel time series in a specified order

## Usage

``` r
extract_nifti_voxel_ts(infile, voxel_indices, time_first = FALSE)
```

## Arguments

- infile:

  Path to a scalar NIfTI image, optionally gzipped.

- voxel_indices:

  One-based linear indices within the spatial volume.

- time_first:

  Logical; if `TRUE`, returns a `volumes x voxels` matrix directly via
  cache-blocked extraction without an intermediate `voxels x volumes`
  allocation. Defaults to `FALSE` for backward compatibility.

## Value

A numeric matrix (`voxels x volumes` when `time_first = FALSE`, or
`volumes x voxels` when `time_first = TRUE`) in the supplied voxel
order.

## Details

Uses the RNifti C++ API with single-dispatch typed pointer access,
including NIfTI slope/intercept scaling, without creating a full
double-precision R array or per-element virtual dispatch. Spatial
alignment with the atlas/mask must be checked by the caller.
