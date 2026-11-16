# Apply non-aggressive (or aggressive) ICA-AROMA denoising to MEDuSA time-series or confound objects

Removes ICA-AROMA motion/artifact components in memory using vectorized
QR least-squares projection across all voxels simultaneously, without
writing intermediate 4D NIfTIs to disk or breaking MEDuSA's continuous
acquisition clock.

## Usage

``` r
apply_aroma(
  x,
  aroma = NULL,
  mixing = NULL,
  classification = NULL,
  mode = c("nonaggressive", "aggressive"),
  clean_nuisance = TRUE,
  fit_valid_rows = FALSE,
  valid_rows = NULL,
  preserve_mean = TRUE,
  force = FALSE
)
```

## Arguments

- x:

  Target object to denoise: a `"medusa_nifti_ts"` object (updates
  `$bold` and, when `clean_nuisance = TRUE`, attached
  `$nuisance`/`$confounds`), a `"medusa_confounds"` object (updates
  `$nuisance`), a run-named list of `T x V` numeric BOLD matrices, or a
  single `T x V` numeric matrix.

- aroma:

  Optional `"medusa_aroma"` object from
  [`read_fmriprep_aroma`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_aroma.md).
  Required unless `mixing` is supplied directly.

- mixing:

  Optional MELODIC mixing file path(s) or `T x K` matrix / list of
  matrices passed to
  [`read_fmriprep_aroma`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_aroma.md)
  when `aroma` is `NULL`.

- classification:

  Optional component classification path(s), metrics `data.frame`,
  logical mask, or 1-based integer vector of rejected ICs passed to
  [`read_fmriprep_aroma`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_aroma.md)
  when `aroma` is `NULL`.

- mode:

  Denoising projection mode: `"nonaggressive"` (default; Pruim et al.,
  2015 partial regression protecting accepted signal ICs) or
  `"aggressive"` (full orthogonal projection of rejected ICs only).

- clean_nuisance:

  Logical; if `TRUE` (default) and `x` is a `"medusa_nifti_ts"` object
  with `$confounds` or `$nuisance` attached, applies the identical AROMA
  projection to nuisance columns so downstream joint regression cannot
  re-introduce removed motion variance.

- fit_valid_rows:

  Logical; if `TRUE`, fits AROMA component weights using only uncensored
  `valid_rows` while projecting across all rows. Default `FALSE` fits
  across full continuous MELODIC run series (matching FSL `fsl_regfilt`,
  Pruim et al., 2015, and `BrainGnomes`).

- valid_rows:

  Optional logical vector or run-named list of logical masks used when
  `fit_valid_rows = TRUE` (defaults to `x$valid_rows` when present on
  `x`).

- preserve_mean:

  Logical; if `TRUE` (default), preserves each voxel's temporal mean
  across fitted timepoints. Ignored for `"medusa_confounds"` inputs
  (which follow `x$settings$center`).

- force:

  Logical; if `TRUE`, re-applies AROMA even when metadata or
  `BrainGnomes` provenance indicates AROMA was already executed. Default
  `FALSE`.

## Value

Object of the same class and structure as `x` with denoised signal
matrices and AROMA diagnostic metadata attached (`x$aroma` on
`"medusa_nifti_ts"`, `x$settings$aroma` on `"medusa_confounds"`, or
`attr(x, "aroma")` on plain matrices / lists).

## Details

Under `mode = "nonaggressive"` (the default; Pruim et al., 2015;
matching FSL `fsl_regfilt` and `BrainGnomes`), all accepted and rejected
MELODIC components are fit jointly with an intercept (\\Y =
\mathbf{1}\mu^\top + M\_{\mathrm{acc}} B\_{\mathrm{acc}} +
M\_{\mathrm{rej}} B\_{\mathrm{rej}} + E\\), and *only* the partial
contribution of rejected components (\\\tilde{M}\_{\mathrm{rej}}
\hat{B}\_{\mathrm{rej}}\\) is subtracted. Shared variance between
task-evoked neural components (\\M\_{\mathrm{acc}}\\) and
task-correlated motion components (\\M\_{\mathrm{rej}}\\) is therefore
retained rather than attenuated, and temporal degrees of freedom in Step
1
([`estimate_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md))
and Step 2
([`deconvolve_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/deconvolve_crossrun_hrf.md))
are not consumed by dozens of ICA columns. When `clean_nuisance = TRUE`,
attached nuisance regressors (`$nuisance` / `$confounds`) receive the
identical AROMA projection so downstream joint nuisance regression
cannot re-introduce AROMA motion variance (Lindquist et al., 2019).

## See also

[`read_fmriprep_aroma`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_aroma.md),
[`read_fmriprep_subject`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_subject.md),
[`read_fmriprep_confounds`](https://hallquistlab.github.io/medusa/reference/read_fmriprep_confounds.md),
[`estimate_crossrun_hrf`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md)

## Examples

``` r
set.seed(1)
s_task <- sin(seq(0, 6 * pi, length.out = 80))
m_motion <- 0.5 * s_task + rnorm(80, sd = 0.6)
mix <- cbind(IC_01 = s_task, IC_02 = m_motion, IC_03 = rnorm(80))
Y <- 500 + cbind(v1 = 2.0 * s_task + 1.5 * m_motion, v2 = 1.0 * s_task - 0.8 * m_motion) +
  matrix(rnorm(160, sd = 0.05), 80, 2)

ar <- read_fmriprep_aroma(mix, classification = 2L)
Y_clean <- apply_aroma(Y, ar, mode = "nonaggressive")
attr(Y_clean, "aroma")$summary
#>    run          mode n_volumes n_fit_volumes n_components n_accepted n_rejected
#> 1 run1 nonaggressive        80            80            3          2          1
#>   applied mean_pct_var_removed median_pct_var_removed
#> 1    TRUE             7.811599               7.811599
```
