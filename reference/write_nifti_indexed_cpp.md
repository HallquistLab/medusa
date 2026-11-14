# Write indexed voxel time series or static maps directly to a 3D/4D NIfTI file

Write indexed voxel time series or static maps directly to a 3D/4D NIfTI
file

## Usage

``` r
write_nifti_indexed_cpp(
  reference,
  out_file,
  ts_mat,
  voxel_indices,
  tr = 1,
  time_offset = 0,
  fill_value = 0,
  datatype = "float",
  compression = 2L,
  time_first = TRUE
)
```

## Arguments

- reference:

  File path to a reference NIfTI image or an in-memory NIfTI object.

- out_file:

  Output file path (`.nii` or `.nii.gz`).

- ts_mat:

  Numeric matrix (`time x voxels` if `time_first = TRUE`, else
  `voxels x time`).

- voxel_indices:

  Unique one-based spatial linear indices matching the voxel dimension
  of `ts_mat`.

- tr:

  Repetition time in seconds stored in `pixdim[5]` for 4D outputs.

- time_offset:

  Acquisition start offset stored in `toffset`.

- fill_value:

  Value written to unindexed background voxels.

- datatype:

  Output NIfTI datatype (`"float"` or `"double"`).

- compression:

  Gzip compression level (`0`–`9`) for `.nii.gz`.

- time_first:

  Logical; whether `ts_mat` is `time x voxels`.
