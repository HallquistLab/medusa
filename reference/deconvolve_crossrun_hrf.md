# Deconvolve runs using frozen HRFs learned from other runs

Deconvolve runs using frozen HRFs learned from other runs

## Usage

``` r
deconvolve_crossrun_hrf(
  bold,
  hrf_fit,
  lambda = NULL,
  scaling = "zscore",
  penalty = c("diff2", "diff1", "ridge", "matern1", "sobolev"),
  backend = c("cpp", "r", "wiener"),
  batch_size = 2048L,
  tail_padding = 0L,
  ridge_floor = 1e-08,
  operator_cache_mb = 128,
  nuisance_mode = NULL,
  observation_mask = NULL,
  lambda_grid = 10^seq(-2, 1, length.out = 16),
  gcv_rule = c("min", "1se"),
  lambda_mode = c("auto", "global", "stratified"),
  max_lambda_strata = 4L,
  rwls = "auto",
  highpass = NULL
)
```

## Arguments

- bold:

  Either a `"medusa_nifti_ts"` object from
  [`read_nifti_ts`](https://hallquistlab.github.io/medusa/reference/read_nifti_ts.md)
  or named time-by-voxel run matrices with IDs/order/dimensions matching
  the HRF fit. These should be the same preprocessed BOLD used for
  estimation.

- hrf_fit:

  Object returned by estimate_crossrun_hrf().

- lambda:

  Nonnegative regularization strength (either a scalar or a per-voxel
  numeric vector of length `ncol(bold[[1]])`). If `NULL` (the default
  for modern `ar1`/`calibrated` fits), selected automatically from
  `lambda_grid` by gain-calibrated, task-weighted adaptive prewhitened
  generalized cross-validation using `gcv_rule` and `lambda_mode`.

- scaling:

  Input scaling: zscore (default), noise_sd, demean, none, or psc.

- penalty:

  Quadratic penalty: diff2, diff1, ridge, matern1, or sobolev.

- backend:

  Fixed linear inverse backend: cpp, r, or scalar wiener.

- batch_size:

  Maximum voxel columns per scaling/multiplication block.

- tail_padding:

  Number of appended zero BOLD observations, default zero. Use
  nrow(hrf_fit\$kernels) to reproduce the current voxelwise wrapper
  assumption.

- ridge_floor:

  Numerical diagonal ridge for cpp/r.

- operator_cache_mb:

  Maximum inverse-operator cache size in MiB, default 128. Zero disables
  reuse across runs. A single operator may exceed the cache budget while
  it is being applied, but is not retained in the cache.

- nuisance_mode:

  Nuisance handling for intercept, linear trend, DCT highpass drift
  (`highpass`), and supplied `nuisance` columns: `"preproject"` (default
  on modern AR/calibrated fits with `cpp` or `r` backends;
  observation-weighted orthogonal residualization before applying the
  continuous inverse operator), `"joint"` (simultaneous unpenalized
  nuisance normal equations), or `"ignore"` (default on legacy OLS or
  `wiener` fits).

- observation_mask:

  Use fit's valid_rows in the inverse (fit), or all observations (all).
  Defaults to fit for AR/calibrated, all for OLS/legacy.

- lambda_grid:

  Candidate regularization values evaluated when `lambda = NULL`.

- gcv_rule:

  Selection rule (`"min"` or `"1se"`) used when `lambda = NULL`.

- lambda_mode:

  Spatial regularization tuning mode when `lambda = NULL`: `"auto"`
  (default; stratifies `lambda` across out-of-sample
  task-informativeness tiers when the map exhibits SNR heterogeneity
  across active and background cortex, and falls back to `"global"` on
  homogeneous maps), `"global"` (single scalar `lambda` for all voxels),
  or `"stratified"` (monotone SNR-tier `lambda` up to
  `max_lambda_strata` strata).

- max_lambda_strata:

  Maximum number of SNR regularization strata when `lambda = NULL` and
  `lambda_mode` is `"auto"` or `"stratified"` (default `4L`).

- rwls:

  Robust weighted least-squares frame weighting mode for transient
  temporal artifacts (e.g., sub-censoring-threshold motion glitches or
  RF spikes): `"auto"` (default; enables `"woodbury"` on modern
  AR/calibrated fits and falls back to `"none"` on legacy OLS or
  `backend = "wiener"`), `"none"` (or `FALSE`; unweighted GLS),
  `"woodbury"` (or `TRUE`; additive-transient Sparse Bayesian Learning
  Woodbury covariance \\\Sigma\_{\mathrm{AR}} +
  \mathrm{diag}(\delta^2)\\), or `"innovation"` (Diedrichsen-style
  innovation-space downweighting).

- highpass:

  Optional DCT highpass cutoff period in seconds inherited from
  `hrf_fit$settings$highpass` when `NULL` (`128` s in modern mode; `Inf`
  / disabled when `Inf` or `<= 0`).

## Value

An S3 list of class `"crossrun_decon"` with `activity` matrices, per-run
`scaling`/status tables (including per-voxel `lambda`), the `hrf_fit`,
frozen inverse `settings` (including scalar reference `settings$lambda`,
the named per-voxel penalty vector `settings$lambda_by_voxel` ready for
1-line 3D NIfTI export via
[`write_nifti_ts`](https://hallquistlab.github.io/medusa/reference/write_nifti_ts.md),
and SNR tier summary table `settings$lambda_strata` when stratified),
standardized `noise_innovations`, cached uncensored `bold_obs` (for
\\\le 16{,}384\\ voxels) or compact \\m_r \times 64\\ `whitened_sketch`
factors (above \\16{,}384\\ voxels, enabling bold-free Stage-3
`"poet_covreg"` shared-error estimation without retaining full 4D BOLD
arrays), per-run observation `frame_weights`, and optional NIfTI
`voxel_map` and `geometry` metadata (enabling direct
[`write_nifti_ts`](https://hallquistlab.github.io/medusa/reference/write_nifti_ts.md)
export without an external reference file). Activity matrices retain
voxel IDs and carry `tr`, `time_offset`, and `observation_mask`
attributes for subsequent event alignment.

## Details

One inverse operator is constructed per selected HRF/noise/lambda
stratum and run length, cached for consecutive runs with identical
masks/nuisance designs, and applied using block matrix products for
cpp/r. Scalar Wiener uses direct block FFT fits. HRFs are frozen; no
voxelwise factorization or R fitting loop is required. AR/calibrated
fits use a masked GLS inverse with the frozen run-specific noise
covariance. A gap of d acquisitions uses rho^d and the corresponding
accumulated innovation variance. Censored rows contribute no
observations; latent activity is still estimated on the complete
original grid. Its observation_mask attribute identifies modeled versus
observed samples. Nonfinite/constant signals on usable observations
remain NA. Nonfinite censored observations do not invalidate a run.
Joint nuisance fitting avoids projecting only the BOLD and inadvertently
changing the task model. Wiener requires OLS, nuisance_mode = "ignore",
observation_mask = "all". Scaling is applied once before optional
padding, and its centers/divisors are returned. Different run z-scores
do not preserve absolute between-run units.
