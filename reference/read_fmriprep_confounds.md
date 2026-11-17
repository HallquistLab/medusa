# Read fMRIPrep or postprocessed confound tables into clock-preserving MEDuSA nuisance and censoring inputs

Converts fMRIPrep `*_desc-confounds_timeseries.tsv` or `BrainGnomes`
`*_desc-<desc>_confounds.tsv` files (or in-memory confound
`data.frame`s) into run-matched `valid_rows` logical masks and finite
`nuisance` design matrices ready for
[`estimate_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md)
and
[`deconvolve_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/deconvolve_crossrun_hrf.md).
Enforces MEDuSA's clock-preservation invariant: censored volumes are
marked `FALSE` in `valid_rows` rather than deleted, initial `"n/a"`
derivative entries on volume 1 are zero-filled, and censored or missing
nuisance entries are safely imputed so `nuisance` is strictly finite on
all rows.

## Usage

``` r
read_fmriprep_confounds(
  confounds,
  bold = NULL,
  run_names = NULL,
  strategy = c("24p", "6p", "12p", "csf_wm", "csf_wm_8p", "24p_csf_wm", "27p", "36p",
    "compcor", "custom", "none"),
  columns = NULL,
  column_pattern = NULL,
  n_compcor = NULL,
  compcor_type = c("anatomical", "split", "temporal"),
  include_cosines = "auto",
  censor = NULL,
  fd_threshold = 0.5,
  fd_col = "auto",
  head_radius = 50,
  rot_units = c("rad", "deg"),
  tr = NULL,
  motion_bandstop_bpm = NULL,
  motion_lowpass_hz = NULL,
  motion_filter_order = NULL,
  motion_filter_pad = c("constant", "odd", "even", "zero"),
  motion_filter_target = c("fd_only", "fd_and_nuisance"),
  dvars_threshold = NULL,
  dvars_col = "std_dvars",
  include_non_steady_state = TRUE,
  include_motion_outliers = FALSE,
  drop_initial_volumes = 0L,
  censor_prev = 0L,
  censor_next = 0L,
  min_segment_length = 1L,
  impute = c("mean", "zero"),
  impute_censored = TRUE,
  center = TRUE,
  scale = FALSE,
  drop_constant = TRUE,
  aroma = NULL,
  aroma_mode = c("nonaggressive", "aggressive"),
  aroma_fit_valid = FALSE,
  na_strings = c("n/a", "N/A", "NA", ".", "NaN", "")
)
```

## Arguments

- confounds:

  Character vector of confound TSV paths
  (`*_desc-confounds_timeseries.tsv`, `*_desc-<desc>_confounds.tsv`,
  `.tsv.gz`, `.csv`), a single `data.frame`, or a list of per-run
  `data.frame`s.

- bold:

  Optional `"medusa_nifti_ts"` object, named list of run BOLD matrices,
  or single BOLD matrix used to match run ordering, validate BIDS
  `sub-`/`ses-`/`task-`/`space-`/`cohort-` entities, check `BrainGnomes`
  provenance sidecars against volume deletion, and verify exact volume
  counts (`nrow(confounds) == nrow(bold)`).

- run_names:

  Optional explicit character vector of unique run names when `bold` is
  `NULL`.

- strategy:

  Predefined nuisance regressor set:

  - `"24p"` (default): 6 realignment parameters (`trans_x`, `trans_y`,
    `trans_z`, `rot_x`, `rot_y`, `rot_z`), their temporal derivatives,
    and quadratic terms (Friston 24-parameter model).

  - `"12p"`: 6 motion parameters + 6 temporal derivatives.

  - `"6p"`: 6 rigid-body motion parameters.

  - `"csf_wm"`: mean `csf` and `white_matter` physiological signals (2p
    model; recommended alongside non-aggressive ICA-AROMA per Pruim et
    al., 2015 and Parkes et al., 2018).

  - `"csf_wm_8p"`: `csf` and `white_matter` with their temporal
    derivatives and quadratic terms (8p physiological model).

  - `"24p_csf_wm"`: `"24p"` + `csf` and `white_matter` with their
    derivatives and quadratic terms (32p physiological model without
    GSR).

  - `"27p"`: `"24p"` + mean `csf`, `white_matter`, and `global_signal`
    without physiological derivatives or quadratics.

  - `"36p"`: `"24p_csf_wm"` + `global_signal` with its derivative and
    quadratic terms (Satterthwaite 36p model).

  - `"compcor"`: `"12p"` motion + top `n_compcor` (default `5`) CompCor
    components + discrete cosine high-pass basis.

  - `"custom"`: only columns requested via `columns`, `column_pattern`,
    and/or `n_compcor`.

  - `"none"`: empty 0-column nuisance matrices (censoring only).

- columns:

  Optional character vector of explicit confound column names, shorthand
  sets (`"6p"`, `"12p"`, `"24p"`, `"27p"`, `"csf_wm"`, `"csf_wm_8p"`,
  `"36p"`), or zero-padded index range specifications such as
  `"a_comp_cor_<0-4>"` or `"cosine<0,2-3>"` (in addition to `strategy`,
  or exact set when `strategy = "custom"`).

- column_pattern:

  Optional regular expression selecting additional columns present in
  each confound table (e.g., `"^aroma_motion_"`).

- n_compcor:

  Optional integer count of leading CompCor components to add (defaults
  to `5L` when `strategy = "compcor"` and `0L` otherwise).

- compcor_type:

  CompCor family when `n_compcor > 0`: `"anatomical"` (`a_comp_cor_*`,
  falling back to split `c_comp_cor_*` + `w_comp_cor_*`), `"split"`
  (`c_comp_cor_*` + `w_comp_cor_*`), or `"temporal"` (`t_comp_cor_*`).

- include_cosines:

  Whether to include fMRIPrep discrete cosine high-pass columns
  (`cosine00`, `cosine01`, ...): `"auto"` (default; included
  automatically whenever CompCor columns are selected, per fMRIPrep
  recommendations), `TRUE`, or `FALSE`.

- censor:

  Optional external censoring mask(s): character vector of `*_censor.1D`
  (`1` = keep, `0` = censor) or `*_scrub.tsv` (one-hot spike indicator)
  paths, a logical or `0/1` vector (single run), or a list of per-run
  logical/`0-1` vectors. When `censor` is supplied and `fd_threshold` is
  not explicitly set, `fd_threshold` defaults to `NULL` so the external
  censor mask is not overridden by a second FD cutoff.

- fd_threshold:

  Framewise displacement censoring cutoff in mm (default `0.5`; volumes
  with `FD > fd_threshold` are marked `FALSE` in `valid_rows`). Set
  `NULL` or `Inf` to disable FD censoring.

- fd_col:

  Column name supplying framewise displacement, or `"auto"` (default),
  which prefers `"framewise_displacement_filtered"` when present (e.g.,
  respiratory-filtered `BrainGnomes` postprocessing outputs), then
  `"framewise_displacement"`, and falls back to computing Power (2012)
  FD from the 6 rigid-body motion parameters when missing or all `NA`.

- head_radius:

  Head radius in mm (default `50`) used when computing Power (2012)
  framewise displacement from the 6 rigid-body motion columns.

- rot_units:

  Rotation parameter units (`"rad"` default for fMRIPrep, or `"deg"`)
  used when computing Power (2012) FD from the 6 rigid-body motion
  columns.

- tr:

  Optional repetition time in seconds (scalar or per-run numeric
  vector), required only when `motion_bandstop_bpm` or
  `motion_lowpass_hz` is specified and TR cannot be inferred
  automatically from `bold` or a sibling `*_bold.json` /
  `*_bold.nii(.gz)` file.

- motion_bandstop_bpm:

  Optional length-2 numeric vector `c(min_bpm, max_bpm)` specifying a
  respiratory band-stop (notch) filter in breaths per minute (e.g.,
  `c(12, 18)` for typical adult breathing at 0.2–0.3 Hz, or
  `c(18.58, 25.73)` per Fair et al., 2020) applied zero-phase via
  forward-backward biquad filtering to the 6 rigid-body motion
  parameters before recomputing Power (2012) FD. Suppresses fast-TR
  multiband respiratory B0 pseudomotion that would otherwise over-censor
  clean volumes. Mutually exclusive with `motion_lowpass_hz`. Default
  `NULL`.

- motion_lowpass_hz:

  Optional numeric scalar specifying a Butterworth low-pass cutoff
  frequency in Hz (e.g., `0.2` Hz; Gratton et al., 2020) applied
  zero-phase via forward-backward filtering to the 6 rigid-body motion
  parameters before recomputing Power (2012) FD. Mutually exclusive with
  `motion_bandstop_bpm`. Default `NULL`.

- motion_filter_order:

  Optional integer effective zero-phase filter order (default `NULL`,
  which uses `4L`). Must be a positive multiple of `4` for
  `motion_bandstop_bpm` (each forward-backward 2nd-order biquad pass
  contributes order 4) or an even integer `>= 2` for `motion_lowpass_hz`
  (a forward-backward `N`-pole Butterworth low-pass has effective order
  `2 * N`).

- motion_filter_pad:

  Boundary reflection padding mode for zero-phase forward-backward
  filtering: `"constant"` (default), `"odd"`, `"even"`, or `"zero"`.

- motion_filter_target:

  Whether respiratory motion filtering updates only Power (2012) FD and
  censoring (`"fd_only"`, default; preserves unfiltered motion nuisance
  columns so respiratory B0 variance remaining in unfiltered BOLD is
  still modeled) or also replaces the 6 rigid-body motion columns and
  their Friston-24 derivatives/quadratics in `nuisance`
  (`"fd_and_nuisance"`).

- dvars_threshold:

  Optional DVARS censoring cutoff (e.g., `1.5` on `std_dvars`). Default
  `NULL` (disabled).

- dvars_col:

  Column name supplying DVARS when `dvars_threshold` is set. Default
  `"std_dvars"`.

- include_non_steady_state:

  Logical; if `TRUE` (default), marks any volume flagged in
  `non_steady_state_outlier*` columns as `FALSE` in `valid_rows`.

- include_motion_outliers:

  Logical; if `TRUE`, also marks volumes flagged in fMRIPrep's
  precomputed `motion_outlier*` columns as `FALSE` in `valid_rows`.
  Default `FALSE` when controlled via `fd_threshold`/`dvars_threshold`.

- drop_initial_volumes:

  Integer count of initial dummy volumes to mark `FALSE` in `valid_rows`
  without deleting rows (default `0L`).

- censor_prev:

  Integer count of preceding volumes to scrub before each
  FD/DVARS/motion-outlier spike (default `0L`).

- censor_next:

  Integer count of subsequent volumes to scrub after each
  FD/DVARS/motion-outlier spike (default `0L`).

- min_segment_length:

  Minimum contiguous uncensored segment length in volumes; shorter
  isolated islands between motion spikes are also marked `FALSE`
  (default `1L`).

- impute:

  Imputation rule for missing entries (beyond volume-1 derivative zeros)
  and censored rows: `"mean"` (usable-row column mean, default) or
  `"zero"`.

- impute_censored:

  Logical; if `TRUE` (default), replaces censored rows (`!valid_rows`)
  in `nuisance` with the imputed baseline value so motion spikes have
  zero leverage even if inspected or used unmasked.

- center:

  Logical; if `TRUE` (default), centers each nuisance column on
  `valid_rows`.

- scale:

  Logical; if `TRUE`, scales each nuisance column to unit SD on
  `valid_rows`. Default `FALSE`.

- drop_constant:

  Logical; if `TRUE` (default), removes columns with near-zero variance
  on `valid_rows` (such as all-zero outlier columns).

- aroma:

  Optional `"medusa_aroma"` object from
  [`read_fmriprep_aroma`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_aroma.md).
  When supplied, applies the matching ICA-AROMA projection
  (`aroma_mode`) to selected nuisance columns prior to final centering
  and censoring sanitization so subsequent joint nuisance regression in
  [`estimate_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md)
  and
  [`deconvolve_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/deconvolve_crossrun_hrf.md)
  cannot re-introduce AROMA motion variance into AROMA-denoised BOLD
  (Lindquist et al., 2019). Framewise displacement and `valid_rows`
  remain strictly derived from physical (or respiratory-filtered)
  head-motion parameters.

- aroma_mode:

  AROMA regression mode when `aroma` is supplied: `"nonaggressive"`
  (default; Pruim et al., 2015 partial regression) or `"aggressive"`.

- aroma_fit_valid:

  Logical; if `TRUE`, estimates AROMA component weights on uncensored
  `valid_rows` only rather than full MELODIC run timecourses (`FALSE`,
  default).

- na_strings:

  Character vector of missing-value tokens passed to
  [`data.table::fread`](https://rdrr.io/pkg/data.table/man/fread.html).
  Default `c("n/a", "N/A", "NA", ".", "NaN", "")`.

## Value

An S3 list of class `"medusa_confounds"` containing:

- `nuisance`: Run-named list of finite numeric `T x P` matrices.

- `valid_rows`: Run-named list of length-`T` logical observation masks
  (`TRUE` = usable, `FALSE` = censored).

- `summary`: Per-run QC `data.frame` reporting volume counts, censored
  counts/percentages, outlier breakdowns, mean/max FD, whether
  respiratory motion filtering was applied, and retained nuisance column
  counts.

- `framewise_displacement`: Run-named list of numeric FD vectors used
  for QC and censoring (respiratory-filtered when `motion_bandstop_bpm`
  or `motion_lowpass_hz` is active).

- `framewise_displacement_raw`: Run-named list of unfiltered FD vectors
  when `motion_bandstop_bpm` or `motion_lowpass_hz` is requested (`NULL`
  otherwise).

- `variance_confounds`: Run-named list of numeric frame-level
  variance-envelope matrices (`framewise_displacement`,
  `std_dvars`/`dvars`, and optional `rvt`) for
  [`decon_epoch_regression`](https://hallquistlab.github.io/medusa/reference/decon_epoch_regression.md).

The `"medusa_confounds"` object can be passed directly as
`nuisance = conf` to
[`estimate_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md)
(automatically supplying both `nuisance` and `valid_rows`), as
`variance_confounds = conf` to
[`decon_epoch_regression`](https://hallquistlab.github.io/medusa/reference/decon_epoch_regression.md),
or as `valid = conf` to
[`align_medusa_ts`](https://hallquistlab.github.io/medusa/reference/align_medusa_ts.md).

## See also

[`read_bids_events`](https://hallquistlab.github.io/medusa/reference/read_bids_events.md),
[`read_fmriprep_subject`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_subject.md),
[`read_fmriprep_aroma`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_aroma.md),
[`apply_aroma`](https://hallquistlab.github.io/medusa/reference/apply_aroma.md),
[`read_nifti_ts`](https://hallquistlab.github.io/medusa/reference/read_nifti_ts.md),
[`estimate_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md)

## Examples

``` r
f_conf <- tempfile("sub-01_task-mid_run-01_", fileext = "_desc-confounds_timeseries.tsv")
conf_df <- data.frame(
  trans_x = c(0.01, 0.02, 0.40, 0.03, 0.02),
  trans_y = c(0.00, -0.01, 0.25, 0.01, 0.00),
  trans_z = c(0.02, 0.03, 0.80, 0.04, 0.03),
  rot_x = c(0.001, 0.001, 0.01, 0.002, 0.001),
  rot_y = c(0.000, 0.001, 0.005, 0.001, 0.001),
  rot_z = c(-0.001, 0.000, 0.004, 0.000, -0.001),
  framewise_displacement = c("n/a", "0.08", "0.92", "0.12", "0.05"),
  non_steady_state_outlier00 = c(1L, 0L, 0L, 0L, 0L),
  stringsAsFactors = FALSE
)
utils::write.table(conf_df, f_conf, sep = "\t", row.names = FALSE, quote = FALSE)

conf <- read_fmriprep_confounds(f_conf, strategy = "12p", fd_threshold = 0.5)
conf$valid_rows$run1
#> [1] FALSE  TRUE FALSE  TRUE  TRUE
conf$summary
#>    run
#> 1 run1
#>                                                                                file
#> 1 /tmp/RtmpIDu8ce/sub-01_task-mid_run-01_19d61fdae2cd_desc-confounds_timeseries.tsv
#>   n_volumes n_valid n_censored pct_censored n_fd_outliers n_dvars_outliers
#> 1         5       3          2           40             1                0
#>   n_non_steady_state n_external_censor mean_fd max_fd motion_filtered
#> 1                  1                 0  0.2925   0.92           FALSE
#>   aroma_cleaned n_nuisance_cols
#> 1         FALSE              11
unlink(f_conf)
```
