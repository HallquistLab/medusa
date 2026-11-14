# Write parcel-level time series or statistics directly to a 3D/4D NIfTI file

Write parcel-level time series or statistics directly to a 3D/4D NIfTI
file

## Usage

``` r
write_nifti_parcel_cpp(
  reference,
  out_file,
  parcel_slot,
  parcel_values,
  tr = 1,
  fill_value = 0,
  datatype = "float",
  compression = 2L
)
```

## Arguments

- reference:

  File path to a reference 3D NIfTI mask or an in-memory NIfTI object.

- out_file:

  Output file path (`.nii` or `.nii.gz`).

- parcel_slot:

  Integer vector of length `n_spatial` with `0` for background voxels
  and `1..M` indexing rows of `parcel_values`.

- parcel_values:

  Numeric matrix of dimensions `M x T` (`parcels x time`).

- tr:

  Repetition time in seconds stored in `pixdim[5]` for 4D outputs.

- fill_value:

  Value written to background voxels (`parcel_slot == 0`).

- datatype:

  Output NIfTI datatype (`"float"` or `"double"`).

- compression:

  Gzip compression level (`0`–`9`) for `.nii.gz`.
