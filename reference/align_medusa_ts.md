# Event-lock and interpolate MEDuSA time-series matrices onto an epoch grid

Matrix-native counterpart to \[event_lock_ts()\] and
\[get_medusa_interpolated_ts()\] that aligns and interpolates directly
on \`time x voxels\` (\`T x V\`) numeric matrices without expanding
intermediate long tables or constructing \`fmri_ts\` R6 objects.

Preserves every audit, censoring, gap-fencing, collision/overlap, and
kernel quadrature guarantee of \[get_medusa_interpolated_ts()\], while
running in C++ across single-run matrices, \`"medusa_nifti_ts"\` bundles
from \[read_nifti_ts()\], deconvolution outputs from
\[deconvolve_crossrun_hrf()\], multi-run matrix lists, multivariate
signal lists, and legacy \`fmri_ts\` objects.

## Usage

``` r
align_medusa_ts(
  ts,
  events = NULL,
  event = NULL,
  time_before = -3,
  time_after = 3,
  collide_before = NULL,
  collide_after = NULL,
  pad_before = -1.5,
  pad_after = 1.5,
  output_resolution = NULL,
  tr = NULL,
  time = NULL,
  time_offset = NULL,
  parcel_ids = NULL,
  unit_data = NULL,
  group_by = "trial",
  value_name = "value",
  trial_col = "trial",
  run_col = "run",
  valid = NULL,
  overlap_policy = NULL,
  neighbor_events = NULL,
  max_gap = NULL,
  readout_kernel = "linear",
  readout_width = NULL,
  time_audit = FALSE,
  list_mode = c("auto", "runs", "signals"),
  orientation = NULL,
  ncpus = 1L
)
```

## Arguments

- ts:

  Time-series input: a \`T x V\` numeric \`matrix\`, a
  \`"medusa_nifti_ts"\` object from \[read_nifti_ts()\], a deconvolution
  result list from \[deconvolve_crossrun_hrf()\], a named \`list\` of
  \`T x V\` matrices (multi-run or multivariate signals; see
  \`list_mode\`), or a legacy \`fmri_ts\` object.

- events:

  Trial-level \`data.frame\` (or named \`list\` of per-run
  \`data.frame\`s for multi-run inputs). Optional when \`ts\` is an
  \`fmri_ts\` object or a \`deconvolve_crossrun_hrf()\` fit whose stored
  events include \`trial_col\` and \`event\`.

- event:

  Name of the target event onset column in \`events\`.

- time_before:

  Earliest epoch lag relative to \`event\` in seconds (\`\<= 0\`).
  Default \`-3.0\`.

- time_after:

  Latest epoch lag relative to \`event\` in seconds (\`\>= 0\`). Default
  \`3.0\`.

- collide_before:

  Optional character vector of column names in \`events\` truncating
  earlier observations (\`time \> evt_before\`) when \`overlap_policy =
  "truncate"\`.

- collide_after:

  Optional character vector of column names in \`events\` truncating
  later observations (\`time \< evt_after\`) when \`overlap_policy =
  "truncate"\`.

- pad_before:

  Extra window padding (in seconds, \`\<= 0\`) retained before
  \`time_before\` to supply interpolation brackets at the left boundary.
  Default \`-1.5\`.

- pad_after:

  Extra window padding (in seconds, \`\>= 0\`) retained after
  \`time_after\` to supply interpolation brackets at the right boundary.
  Default \`1.5\`.

- output_resolution:

  Epoch output step in seconds. Defaults to \`tr\`.

- tr:

  Optional positive scalar repetition time in seconds. Inferred from
  \`ts\` (\`ts\$tr\` or \`attr(ts, "tr")\`) when available.

- time:

  Optional numeric vector of length \`T\` (or per-run list) giving
  acquisition timestamps in seconds. Defaults to \`(0:(T - 1)) \* tr +
  time_offset\`.

- time_offset:

  Optional scalar (or per-run numeric vector) start-time offset in
  seconds when \`time\` is not supplied. Defaults to \`attr(mat,
  "time_offset")\` or \`0\`.

- parcel_ids:

  Optional vector of length \`V\` assigning each spatial column to a
  parcel/ROI (populates both \`parcel\` and \`atlas_value\` for
  \`group_by\`).

- unit_data:

  Optional \`V\`-row \`data.frame\` of spatial unit metadata (e.g.,
  \`ts\$voxel_map\`) whose columns may be referenced in \`group_by\`.

- group_by:

  Character vector of grouping variables, or \`NULL\`. For
  \`align_medusa_ts()\`, \`missing(group_by)\` defaults to
  \`c(trial_col, "parcel")\` when \`parcel_ids\` (or \`ts\$parcel_ids\`)
  is present and \`trial_col\` otherwise; passing explicit \`group_by =
  NULL\` aggregates across trials and units.

- value_name:

  Signal name prefix used for single-signal matrix inputs (e.g.,
  \`"value"\` produces \`value_mean\`, \`value_sd\`, ...). Default
  \`"value"\`.

- trial_col:

  Name of the unique within-run trial identifier column in \`events\`.
  Default \`"trial"\`.

- run_col:

  Name of the run identifier column used when \`events\` is a combined
  multi-run \`data.frame\` and attached to multi-run outputs. Default
  \`"run"\`.

- valid:

  Optional logical validity mask (\`TRUE\` = usable, \`FALSE\`/\`NA\` =
  censored): a length-\`T\` vector, a \`T x V\` logical matrix, a
  per-run list, or \`"observation_mask"\` to use \`attr(mat,
  "observation_mask")\`. For \`fmri_ts\` inputs, a single column name
  string is also accepted (matching \`valid_col\`).

- overlap_policy:

  One of \`"allow"\`, \`"truncate"\`, or \`"exclude"\`. \`NULL\` uses
  \`"truncate"\` when collision columns are supplied and \`"allow"\`
  otherwise.

- neighbor_events:

  Event columns used for \`align\_\*\` nearest-event annotations.
  Defaults to \`event\`.

- max_gap:

  Maximum allowed interpolation bracket duration in seconds (defaults to
  \`tr\`).

- readout_kernel:

  Readout kernel at each output lag: \`"linear"\` (default),
  \`"boxcar"\`, \`"triangle"\`, or \`"gaussian"\`.

- readout_width:

  Full window width in seconds for non-linear \`readout_kernel\` modes
  (defaults to \`1.5 \* max_gap\`).

- time_audit:

  Logical; if \`TRUE\`, appends \`evt_onset\`, \`collide_before\`,
  \`collide_after\`, and absolute \`time\` when \`trial_col\` is in
  \`group_by\`.

- list_mode:

  Interpretation when \`ts\` is a plain \`list\` of matrices: \`"auto"\`
  (default), \`"runs"\` (multi-run list), or \`"signals"\` (multivariate
  signals for one run).

- orientation:

  Optional matrix orientation override (\`"time_by_voxel"\` or
  \`"voxel_by_time"\`).

- ncpus:

  Integer number of parallel workers used across runs when aligning
  multi-run inputs. Default \`1L\`.

## Value

A \`data.table\` of aligned, interpolated epoch summaries and
\`align\_\*\` / \`epoch_status\` diagnostics with a \`source_clock\`
\`data.frame\` attribute.

## Details

\`align_medusa_ts()\` performs descriptive event-locked interpolation
and overlap/censoring auditing without assuming a forward generative
design. When \`ts\` is a \[deconvolve_crossrun_hrf()\] result and the
goal is to remove non-target event bleed, invert target-to-target
cross-talk, model variable durations, or propagate analytical per-trial
standard errors (\`se_demixed\`), use \[decon_epoch_regression()\]
instead.

## See also

\[read_nifti_ts()\], \[deconvolve_crossrun_hrf()\],
\[decon_epoch_regression()\], \[get_medusa_interpolated_ts()\],
\[evt_align_decon_files()\], \[get_medusa_compression_score()\]

## Examples

``` r
set.seed(1)
Y <- matrix(rnorm(60 * 8), nrow = 60, ncol = 8)
parcels <- rep(c("ROI_A", "ROI_B"), each = 4)
events <- data.frame(trial = 1:4, onset = c(8.5, 20.2, 33.7, 46.1))
aligned <- align_medusa_ts(
  Y, events = events, event = "onset", tr = 1.0,
  parcel_ids = parcels, time_before = -2, time_after = 4, output_resolution = 1.0
)
head(aligned)
#>    trial parcel evt_time  value_mean value_median   value_sd value_n
#>    <int> <char>    <num>       <num>        <num>      <num>   <int>
#> 1:     1  ROI_A       -2  0.18420709  0.146826684 0.36228740       4
#> 2:     1  ROI_A       -1  0.12716573  0.148702852 0.70929427       4
#> 3:     1  ROI_A        0  0.02925075 -0.183884447 0.83269516       4
#> 4:     1  ROI_A        1  0.31089884  0.235570737 0.82705471       4
#> 5:     1  ROI_A        2  0.17042016 -0.002382271 0.54287742       4
#> 6:     1  ROI_A        3 -0.08971901 -0.082654355 0.06232842       4
#>    value_n_total value_n_censored value_n_missing value_n_gap value_n_outside
#>            <int>            <int>           <int>       <int>           <int>
#> 1:             4                0               0           0               0
#> 2:             4                0               0           0               0
#> 3:             4                0               0           0               0
#> 4:             4                0               0           0               0
#> 5:             4                0               0           0               0
#> 6:             4                0               0           0               0
#>    value_source_ids value_source_max_epochs epoch_status align_onset
#>              <char>                   <int>       <char>       <num>
#> 1:              7;8                       1           ok         8.5
#> 2:              8;9                       1           ok         8.5
#> 3:             9;10                       1           ok         8.5
#> 4:            10;11                       1           ok         8.5
#> 5:            11;12                       1           ok         8.5
#> 6:            12;13                       1           ok         8.5
#>    align_previous_event align_next_event align_simultaneous_n
#>                   <num>            <num>                <int>
#> 1:                   NA             20.2                    0
#> 2:                   NA             20.2                    0
#> 3:                   NA             20.2                    0
#> 4:                   NA             20.2                    0
#> 5:                   NA             20.2                    0
#> 6:                   NA             20.2                    0
#>    align_window_overlap_n align_run_start_distance align_run_end_distance
#>                     <int>                    <num>                  <num>
#> 1:                      0                      8.5                   50.5
#> 2:                      0                      8.5                   50.5
#> 3:                      0                      8.5                   50.5
#> 4:                      0                      8.5                   50.5
#> 5:                      0                      8.5                   50.5
#> 6:                      0                      8.5                   50.5
#>    align_overlap_policy
#>                  <char>
#> 1:                allow
#> 2:                allow
#> 3:                allow
#> 4:                allow
#> 5:                allow
#> 6:                allow
```
