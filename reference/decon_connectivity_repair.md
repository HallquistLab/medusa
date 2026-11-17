# Measurement-error-corrected inter-region connectivity and PPI via REPAIR

Estimates trial-by-trial seed-to-target neural connectivity, directed
lagged connectivity (controlling for baseline target activity at the
seed lag), and task-modulated psychophysiological interactions (PPI)
from a
[`decon_epoch_regression`](https://hallquistlab.github.io/medusa/reference/decon_epoch_regression.md)
trial panel using the REPAIR modified method-of-moments (`M-MOMS` /
`M-MOM`) framework of Vig & Hallquist (2026). Propagates per-trial
operator standard errors (`se_demixed`) into exact heteroscedastic
predictor and outcome error covariances,including the rank-1 interaction
error block for PPI (`seed * x`) and dual predictor error variances for
Granger-style auto-lag controls (`targ_t0`).

## Usage

``` r
decon_connectivity_repair(
  panel_dt,
  seed_region,
  target_region,
  seed_lag = 0,
  target_lag = 0,
  covariates = NULL,
  test_ppi = FALSE,
  control_target_seed_lag = FALSE,
  orthogonalize_seed_x = TRUE,
  method = c("mmoms", "mmom", "eb_predictor", "naive_gls", "naive_unweighted"),
  shared_error = c("auto", "zero", "supplied"),
  err_xy_col = NULL,
  use_whitened = NULL,
  region_col = NULL,
  id_col = "id",
  run_col = "run",
  trial_col = "trial",
  lag_col = "evt_time",
  y_col = NULL,
  se_col = NULL,
  eb_col = NULL,
  modulator = NULL,
  return_augmented_panel = FALSE,
  crosstrial_se = TRUE,
  control_seed_lags = NULL,
  multilag_control = FALSE,
  match_resolution = TRUE
)
```

## Arguments

- panel_dt:

  A
  [`decon_epoch_regression`](https://hallquistlab.github.io/medusa/reference/decon_epoch_regression.md)
  result list, a trial-level `data.frame`/`data.table` (including CSV,
  `.csv.gz`, Parquet, or RDS panels where pipe-delimited list columns
  are automatically restored), a list of per-subject results/panels, or
  a file path (or character vector of file paths) readable by
  `medusa_read_table()`.

- seed_region:

  Character scalar naming the seed region/parcel/voxel.

- target_region:

  Character scalar naming the target region/parcel/voxel.

- seed_lag:

  Numeric lag (in seconds) for the seed readout (default `0`).

- target_lag:

  Numeric lag (in seconds) for the target readout (default `0`).

- covariates:

  Optional character vector of error-free trial covariate columns in
  `panel_dt`. Defaults to `"x"` when present in `panel_dt`, or
  `character(0)` otherwise.

- test_ppi:

  Logical; if `TRUE`, includes the interaction between the seed readout
  and the first covariate in `covariates` (`seed_c:x_c`) with exact
  per-trial interaction error covariance \\\tilde{\psi}\_{S,j} \[1,
  x\_{c,j}; x\_{c,j}, x\_{c,j}^2\]\\. Default: `FALSE`.

- control_target_seed_lag:

  Logical; if `TRUE` and `abs(target_lag - seed_lag) > 1e-6`, includes
  the target region's readout at `seed_lag` (`targ_t0_val`) as an
  error-prone control covariate with measurement-error variance
  `se_demixed^2`. Default: `FALSE`.

- orthogonalize_seed_x:

  Logical; if `TRUE` (default), residualizes the run-centered seed
  readout on the error-free covariates within each subject-run and
  scales the seed error variance by \\1 - h\_{x,r,j}\\ prior to forming
  PPI interactions, preventing quadratic task drive (\\x_c^2\\) from
  confounding physiological-by-psychological interactions.

- method:

  Estimation method: `"mmoms"` (default; adaptive step-up REPAIR
  `M-MOMS`), `"mmom"` (fixed minimum-regularization REPAIR `M-MOM`),
  `"eb_predictor"` (Skrondal & Laake 2001 role-specific Empirical-Bayes
  shrunk predictor on likelihood-only outcome), `"naive_gls"`
  (precision-weighted WLS without predictor error subtraction), or
  `"naive_unweighted"` (unweighted OLS).

- shared_error:

  Policy for cross-parcel shared observation-error covariance between
  seed and target readouts: `"auto"` (default; propagates the Stage-3
  Adaptive-POET + dynamic covariance-regression `shared_error_model`
  when operator objects are supplied, or the embedded per-trial
  Block-PAF factor columns `shared_noise_vec`, `lag_noise_vec`, and
  `trial_cor_vec`—falling back to stationary `unit_error_cor` /
  `shared_noise_loading` on legacy panels), `"zero"` (assumes zero
  cross-parcel noise correlation while still propagating same-parcel
  cross-lag error correlation when `control_target_seed_lag = TRUE`), or
  `"supplied"` (reads per-trial raw seed–target error covariance from
  `err_xy_col`).

- err_xy_col:

  Optional column name supplying per-trial raw seed–target error
  covariance when `shared_error = "supplied"`.

- use_whitened:

  Logical; if `TRUE`, uses `decon_demixed_whitened` and `se_whitened`
  when present. Defaults to `FALSE` for `"mmoms"` and `"mmom"` (so
  elementwise PPI products are formed on unmixed trial states before
  error correction) and `TRUE` for `"naive_gls"`.

- region_col:

  Optional column name identifying regions/parcels/voxels in `panel_dt`
  (auto-detected from `"region"`, `"parcel"`, `"voxel"`).

- id_col:

  Column name identifying subjects in `panel_dt` (default `"id"`).

- run_col:

  Run identifier column name (default `"run"`).

- trial_col:

  Within-run trial identifier column name (default `"trial"`).

- lag_col:

  Within-trial lag column name (default `"evt_time"`).

- y_col:

  Optional readout column override.

- se_col:

  Optional trial standard-error column override.

- eb_col:

  Optional Empirical-Bayes trial readout column override.

- modulator:

  Optional alias for `covariates`.

- return_augmented_panel:

  Logical; if `TRUE`, returns the merged trial panel augmented with
  moment-preserving REPAIR predictor columns (`seed_repair`,
  `ppi_repair`, `repair_weight`) that reproduce subject-level `M-MOMS`
  slopes in weighted least squares / `lmer()`.

- crosstrial_se:

  Logical; if `TRUE` (default), computes cross-trial (\\i \neq j\\)
  hybrid operator + gated onset-kernel sandwich standard errors and
  Satterthwaite effective degrees of freedom (`df`) for single-subject
  (\\N = 1\\) and subject-level inference while keeping diagonal M-MOMS
  estimation weights.

- control_seed_lags:

  Optional numeric vector of other within-trial seed lags (in seconds,
  present in `panel_dt`) to partial out as measurement-error-prone
  control predictors alongside `seed_lag`.

- multilag_control:

  Logical; if `TRUE`, automatically partials out other seed lags present
  in a multi-lag panel (subject to a latent-variance identifiability
  guard so silent lags do not induce M-MOMS attenuation) and enables
  causal target auto-lag control (`control_target_seed_lag = TRUE`)
  strictly when `target_lag > seed_lag`. Default: `FALSE`.

- match_resolution:

  Logical; if `TRUE` (default), uses Stage-3 post-demixing trial
  resolution operators (`run_trial_operators` or embedded
  `trial_res_vec` bands on stacked/CSV panels) to match differential
  seed–target trial smoothing (\\M_S \neq M_T\\) and post-demixing PPI
  product transfer (\\(M_S \eta_S) \odot x\\ vs. \\M_T (\eta_S \odot
  x)\\) across EIV estimators (`"mmoms"`, `"mmom"`, `"eb_predictor"`).

## Value

A list of class `"medusa_connectivity_repair"` containing `group`
(group-level HKSJ random-effects meta-analytic inference table, or
single-subject cross-trial Satterthwaite \\t\\-inference table when \\N
= 1\\), `subjects` (per-subject estimates, SEs, Satterthwaite `df`, and
REPAIR diagnostics), optional `panel_repair`, and `settings`.

## See also

[`decon_epoch_regression`](https://hallquistlab.github.io/medusa/reference/decon_epoch_regression.md),
[`aggregate_decon_panel`](https://hallquistlab.github.io/medusa/reference/aggregate_decon_panel.md),
[`medusa_fit_mmom`](https://hallquistlab.github.io/medusa/reference/medusa_fit_mmom.md),
[`medusa_between_subject_repair`](https://hallquistlab.github.io/medusa/reference/medusa_between_subject_repair.md)
