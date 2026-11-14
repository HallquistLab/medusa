# Align a set of deconvolved time series files to an event of interest function.

Align a set of deconvolved time series files to an event of interest
function.

## Usage

``` r
evt_align_decon_files(
  d_files,
  trial_df,
  alignment = list(),
  aggregate_by = "atlas_value",
  tr = NULL,
  atlas_cuts = NULL,
  atlas_subset = NULL,
  out_file = NULL,
  ncpus = 8,
  time_audit = FALSE,
  file_manifest = NULL,
  fail_on_error = TRUE
)
```

## Arguments

- d_files:

  A vector of deconvolved `.csv`/`.csv.gz`, `.rds`, `.nii`/`.nii.gz`, or
  `.parquet` files (created by `voxelwise_deconvolution`).

- trial_df:

  The trial-level data.frame containing id and run for each subject
  represented in `d_files`.

- alignment:

  A list containing alignment details passed to
  [`align_medusa_ts`](https://hallquistlab.github.io/medusa/reference/align_medusa_ts.md)
  /
  [`get_medusa_interpolated_ts`](https://hallquistlab.github.io/medusa/reference/get_medusa_interpolated_ts.md),
  including `evt_col`, `time_before`, `time_after`, `pad_before`,
  `pad_after`, `output_resolution`, `collide_before`, `collide_after`,
  `overlap_policy`, `valid_col`, `neighbor_events`, `max_gap`,
  `readout_kernel`, and `readout_width`.

- aggregate_by:

  The column name in the individual deconvolved files used for averaging
  repeated units (e.g., voxels) into single event-aligned time series.
  Most commonly, this is "atlas_value", which will lead to averaging of
  voxels within each parcel in the mask.

- tr:

  The repetition time of the sequence in seconds.

- atlas_cuts:

  For a continuous-valued atlas (e.g., containing a gradient of
  interest), a vector of cut points for binning values. These cuts are
  applied to the `aggregate_by` column, commonly "atlas_value"

- atlas_subset:

  An optional numeric vector containing values of `aggregate_by`

- out_file:

  Optional output file path for the event-aligned panel (`.csv`,
  `.csv.gz`, `.parquet`, or `.rds`). If NULL, nothing is written to disk
  (the event-aligned data are always returned).

- ncpus:

  Number of local worker processes to use.

- time_audit:

  If TRUE, additional columns will be added to the output showing how
  alignment is calculated vis-a-vis event timing

- file_manifest:

  Optional data.frame with `file`, `id`, `run`, and optional `session`
  columns. Overrides per-file metadata and legacy filename parsing;
  ambiguous identities are errors.

- fail_on_error:

  Stop if any file fails (default TRUE). FALSE permits partial output
  with a warning. A status row is retained for every file.

## Value

Aligned data with a `file_status` data.frame attribute. If `out_file` is
supplied, the manifest is saved to `paste0(out_file, ".files.csv")`,
including on failure. With `fail_on_error = TRUE`, a
`decon_alignment_error` condition includes the manifest and aligned
output is not written if any file fails. Successful results also carry a
`source_clock` table, saved to `paste0(out_file, ".sources.csv")` when
output is requested. The `alignment_file` column scopes acquisition-time
IDs to their input file.

## Details

This is intended to be used internally by `run_decon_alignment`, which
accepts a set of mask/atlas files and alignments, then processes these
in parallel. Long `.csv`/`.csv.gz` and `.parquet` tables sharing a run
clock are pivoted directly into `time x voxels` matrices in C++ without
intermediate `fmri_ts` allocation (automatically falling back to
`fmri_ts` for ragged per-unit clocks); compact `.rds` matrix bundles
(`list(decon = Y, coords = ..., time = ...)`) and 4D `.nii`/`.nii.gz`
deconvolved images (paired with `.coords.rds` sidecars) are also
supported directly.
