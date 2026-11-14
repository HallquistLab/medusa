# Read and normalize BIDS event files for MEDuSA HRF, deconvolution, and alignment workflows

Reads one or more BIDS `*_events.tsv` files (or normalizes in-memory
BIDS `data.frame`s) into MEDuSA's event table contract: finite numeric
`onset` and non-negative `duration` (in seconds), non-empty `condition`
labels (mapped by default from BIDS `trial_type`), run-matched `run`
labels in `bold` run order, and unique within-run `trial` indices.
Standard BIDS missing tokens (`"n/a"`) are parsed automatically without
coercing numeric behavioral columns to `character`.

## Usage

``` r
read_bids_events(
  events,
  bold = NULL,
  run_names = NULL,
  condition_col = "trial_type",
  condition_map = NULL,
  default_condition = NULL,
  exclude_conditions = NULL,
  na_duration = 0,
  drop_missing_onset = FALSE,
  drop_missing_condition = FALSE,
  time_shift = 0,
  trial_col = "trial",
  run_col = "run",
  sort_by_onset = FALSE,
  na_strings = c("n/a", "N/A", "NA", ".", "NaN", ""),
  output = c("data.frame", "list")
)
```

## Arguments

- events:

  Character vector of BIDS event file paths (`.tsv`, `.tsv.gz`, `.csv`,
  `.csv.gz`), a single `data.frame`, or a list of per-run `data.frame`s.

- bold:

  Optional `"medusa_nifti_ts"` object, named list of run BOLD matrices,
  single BOLD matrix, or character vector of target run names used to
  validate BIDS `sub-`/`ses-`/`task-` entities and order event runs
  identically to BOLD runs.

- run_names:

  Optional explicit character vector of unique run names when `bold` is
  `NULL`.

- condition_col:

  Column name in the BIDS event table supplying condition labels when
  `condition` is not already present (or when different from
  `"condition"`). Default `"trial_type"`.

- condition_map:

  Optional named character vector recoding raw condition values (e.g.,
  `c(win_fb = "feedback", loss_fb = "feedback")`). Unmatched condition
  labels are preserved unchanged.

- default_condition:

  Optional scalar string assigned to `condition` only when neither
  `"condition"` nor `condition_col` exists in a single-condition BIDS
  event table (e.g., `"event"`).

- exclude_conditions:

  Optional character vector of raw condition labels (prior to
  `condition_map`) to exclude from the returned table.

- na_duration:

  Replacement duration (in seconds, `>= 0`) used when `duration` is
  missing or `"n/a"` (default `0` for impulse events). Set `NULL` to
  error on missing durations.

- drop_missing_onset:

  Logical; if `TRUE`, rows with missing or non-finite `onset` are
  dropped instead of raising an error. Default `FALSE`.

- drop_missing_condition:

  Logical; if `TRUE`, rows with missing or empty `condition` labels are
  dropped instead of raising an error. Default `FALSE`.

- time_shift:

  Numeric scalar or run-ordered numeric vector added to `onset` (in
  seconds; default `0`).

- trial_col:

  Name of the within-run unique trial index column populated when not
  already present in `events`. Default `"trial"`.

- run_col:

  Name of the run identifier column populated on returned event tables.
  Default `"run"`.

- sort_by_onset:

  Logical; if `TRUE`, sorts events chronologically by `onset` within
  each run while recording original 1-based file row numbers in
  `bids_event_row`. Default `FALSE` preserves exact file row order for
  1-to-1 behavioral joins.

- na_strings:

  Character vector of missing-value tokens passed to
  [`data.table::fread`](https://rdrr.io/pkg/data.table/man/fread.html).
  Default `c("n/a", "N/A", "NA", ".", "NaN", "")`.

- output:

  Return container format: `"data.frame"` (default; single combined
  multi-run `data.frame` ready to pass directly to
  [`estimate_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md),
  [`decon_epoch_regression`](https://hallquistlab.github.io/medusa/reference/decon_epoch_regression.md),
  and
  [`align_medusa_ts`](https://hallquistlab.github.io/medusa/reference/align_medusa_ts.md))
  or `"list"` (run-named list of per-run `data.frame`s).

## Value

Either a combined multi-run `data.frame` (`output = "data.frame"`) or a
run-named list of `data.frame`s (`output = "list"`).

## See also

[`read_fmriprep_confounds`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_confounds.md),
[`read_fmriprep_subject`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_subject.md),
[`read_nifti_ts`](https://hallquistlab.github.io/medusa/reference/read_nifti_ts.md),
[`estimate_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md),
[`decon_epoch_regression`](https://hallquistlab.github.io/medusa/reference/decon_epoch_regression.md)

## Examples

``` r
f1 <- tempfile("sub-01_task-mid_run-01_", fileext = "_events.tsv")
f2 <- tempfile("sub-01_task-mid_run-02_", fileext = "_events.tsv")
writeLines(c(
  "onset\tduration\ttrial_type\tresponse_time\trpe",
  "2.0\t0.5\tcue\tn/a\tn/a",
  "6.0\t0.8\tfeedback\t0.52\t0.75"
), f1)
writeLines(c(
  "onset\tduration\ttrial_type\tresponse_time\trpe",
  "3.0\t0.5\tcue\tn/a\tn/a",
  "7.5\tn/a\tfeedback\t0.61\t-0.40"
), f2)

ev <- read_bids_events(c(f1, f2))
unlink(c(f1, f2))
```
