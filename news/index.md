# Changelog

## medusa 0.1.0

### Fast C++ NIfTI matrix I/O and matrix-native event alignment

- Added
  [`read_nifti_ts()`](https://hallquistlab.github.io/medusa/reference/read_nifti_ts.md)
  and
  [`write_nifti_ts()`](https://hallquistlab.github.io/medusa/reference/write_nifti_ts.md)
  for streaming C++ conversion between 3D/4D `.nii`/`.nii.gz` images and
  `time x voxels` (`T x V`) numeric matrices (`"medusa_nifti_ts"`
  bundles carrying `voxel_map`, `parcel_ids`, `tr`, and validated
  `sform`/`qform` `geometry` without allocating full 4D R arrays).
  Updated
  [`voxelwise_deconvolution()`](https://hallquistlab.github.io/medusa/reference/voxelwise_deconvolution.md)
  and
  [`backproject_medusa()`](https://hallquistlab.github.io/medusa/reference/backproject_medusa.md)
  to use the streaming C++ NIfTI path.
- Added
  [`align_medusa_ts()`](https://hallquistlab.github.io/medusa/reference/align_medusa_ts.md)
  (`src/epoch_align.cpp`) for C++ matrix-native event-locked alignment
  and multi-kernel readout quadrature (`"linear"`, `"boxcar"`,
  `"triangle"`, `"gaussian"`) across `T x V` matrices,
  `"medusa_nifti_ts"` bundles,
  [`deconvolve_crossrun_hrf()`](https://hallquistlab.github.io/medusa/reference/deconvolve_crossrun_hrf.md)
  outputs, multi-run lists, and multivariate signal lists while
  preserving all `overlap_policy`, censoring, `max_gap`, `epoch_status`,
  `align_*`, `time_audit`, and `source_clock` invariants (`~80x–130x`
  faster and `~10x–135x` lower peak R heap memory than legacy long-table
  workflows).
- Accelerated
  [`get_medusa_interpolated_ts()`](https://hallquistlab.github.io/medusa/reference/get_medusa_interpolated_ts.md)
  and
  [`get_medusa_compression_score()`](https://hallquistlab.github.io/medusa/reference/get_medusa_compression_score.md)
  on shared-clock `fmri_ts` objects via a zero-copy C++ matrix view
  (with automatic fallback for ragged per-unit clocks), added
  `fmri_ts$to_matrix()`, and upgraded
  [`evt_align_decon_files()`](https://hallquistlab.github.io/medusa/reference/evt_align_decon_files.md)
  with `fast_pivot_decon_cpp()`, `.rds` matrix-bundle support, and fork
  parallelism on Unix CLI/HPC.
- Modernized file and tabular serialization across the pipeline:
  [`voxelwise_deconvolution()`](https://hallquistlab.github.io/medusa/reference/voxelwise_deconvolution.md)
  now supports `output_format = c("csv.gz", "rds", "nifti", "parquet")`
  (eliminating intermediate
  [`mat2df()`](https://hallquistlab.github.io/medusa/reference/mat2df.md) +
  [`dplyr::left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)
  expansions and writing delimited outputs via multi-threaded
  [`data.table::fwrite()`](https://rdrr.io/pkg/data.table/man/fwrite.html));
  [`run_decon_alignment()`](https://hallquistlab.github.io/medusa/reference/run_decon_alignment.md)
  and
  [`evt_align_decon_files()`](https://hallquistlab.github.io/medusa/reference/evt_align_decon_files.md)
  discover and align all four deconvolution formats, skip unused
  coordinate/status columns when reading long tables, and write aligned
  panels to `.csv.gz`, `.parquet`, or `.rds`; and
  [`mixed_by()`](https://hallquistlab.github.io/medusa/reference/mixed_by.md)
  now reads `.parquet` files and Hive-partitioned Parquet directories
  (via optional `arrow`) as well as `.rds` and delimited `external_df`
  tables, with automatic column projection and optional `split_filter`
  pushdown.
- Unified in-memory NIfTI, cross-run HRF, deconvolution, and
  epoch-regression ergonomics:
  [`estimate_crossrun_hrf()`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md)
  and
  [`deconvolve_crossrun_hrf()`](https://hallquistlab.github.io/medusa/reference/deconvolve_crossrun_hrf.md)
  now accept `"medusa_nifti_ts"` objects directly (auto-populating `tr`,
  `spatial_groups = ts$parcel_ids`,
  `kernel_grid = hrf_library(tr, preset = "medusa")`, and
  `spatial_coords` when `spatial_fwhm` is requested), split combined
  single-participant multi-run `events` `data.frame`s automatically
  across BIDS/MEDuSA run labels, and propagate `voxel_map` and
  `geometry` onto `"crossrun_hrf"` and `"crossrun_decon"` so
  `write_nifti_ts(decon, ...)` and
  `backproject_medusa(..., brain_mask = decon)` need no extra reference
  file;
  [`decon_epoch_regression()`](https://hallquistlab.github.io/medusa/reference/decon_epoch_regression.md)
  now accepts column-name `covariates` (e.g., `c("rpe", "abs_rpe")`),
  one-sided formulas (`~ rpe + abs_rpe`), column-name `trial_durations`
  (`"rt"`), and combined multi-run `events` tables.

### Standalone package extraction

- Extracted all MEDuSA deconvolution, cross-run HRF estimation,
  event-locked interpolation, operator-based epoch regression,
  connectivity repair, NIfTI backprojection
  ([`backproject_medusa()`](https://hallquistlab.github.io/medusa/reference/backproject_medusa.md)
  migrated to `RNifti`), and split-wise multilevel modeling
  ([`mixed_by()`](https://hallquistlab.github.io/medusa/reference/mixed_by.md),
  [`jlmer()`](https://hallquistlab.github.io/medusa/reference/jlmer.md),
  [`meta_mixed_by()`](https://hallquistlab.github.io/medusa/reference/meta_mixed_by.md))
  tools from `fmri.pipeline` into the standalone `medusa` package.

### MEDuSA learning sequence, precision-weighted epoch regression, and workflow documentation

- Add conceptual (`vignettes/medusa-concepts.Rmd`) and executable
  workflow (`vignettes/medusa-workflows.Rmd`) vignettes covering
  cross-run HRF learning, noise scaling, operator-based epoch de-mixing,
  duration models, per-trial precision weights, within-run ZCA
  whitening, NIfTI inputs, and the handoff to multilevel analysis
  (`mixed_by(..., weights = "precision_weight")`).
- Revise the technical companion
  (`medusa_deconvolution_review_prs_42_53.qmd`) through PR
  [\#55](https://github.com/HallquistLab/medusa/issues/55)+, align
  descriptions with the current implementation, distinguish legacy
  benchmark tables from live experiments, and add thorough inline
  documentation across all MEDuSA R functions.
- [`decon_epoch_regression()`](https://hallquistlab.github.io/medusa/reference/decon_epoch_regression.md)
  now propagates exact run-specific AR(1) observation noise covariance
  through the readout and unmixing operators to compute per-trial
  standard errors (`se_raw`, `se_demixed`), lag-1 trial error
  autocorrelation (`error_ac1`), Empirical-Bayes stabilized precision
  weights (`precision_weight`), ZCA-whitened GLS slopes (`beta_gls`,
  `se_gls`), within-trial multi-lag de-mixing
  (`demix_within_trial = TRUE`), and optional within-run ZCA-whitened
  panel columns (`whiten_panel = TRUE`).
- [`mixed_by()`](https://hallquistlab.github.io/medusa/reference/mixed_by.md)
  now accepts `weights = "precision_weight"` (normalizing positive
  precision weights within each split and forwarding to `lme4`,
  `rstanarm`, and `jlmer`).
- [`hrf_library()`](https://hallquistlab.github.io/medusa/reference/hrf_library.md)
  adds `preset = c("custom", "standard", "extended", "dense")`;
  [`estimate_crossrun_hrf()`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md)
  adds Galwey spectral effective candidate multiplicity
  (`hrf_effective_candidates()`) and Cauchy combination parcel pooling
  (`hrf_cauchy_p()`); and
  [`deconvolve_crossrun_hrf()`](https://hallquistlab.github.io/medusa/reference/deconvolve_crossrun_hrf.md)
  supports automatic prewhitened GCV tuning when `lambda = NULL`.

### Banded linear deconvolution solver and exact Takahashi trace

- `deconvolve_reglin_cpp()` now exploits the strict banded structure of
  $`H_k^\top H_k + \lambda P + \delta I`$ (half-bandwidth
  $`b = \max(\text{kernel\_len} - 1 + \mathbb{I}(\text{prewhiten}), \; b_P)`$),
  replacing dense $`O(n_{\text{latent}}^3)`$ normal-equation solves and
  dense matrix-RHS trace evaluations with an in-place
  $`O(n_{\text{latent}} b^2)`$ lower-banded Cholesky factorization,
  contiguous multi-signal forward/backward substitution, and the exact
  $`O(n_{\text{latent}} b^2)`$ Takahashi banded-inverse recursion for
  effective degrees of freedom
  $`\operatorname{tr}(A(\lambda)^{-1} H_k^\top H_k)`$.
- `hrf_inverse_operator_cpp()`, `hrf_score_run_cpp()`, and
  `hrf_noise_profile_cpp()` now eliminate dense $`m \times m`$ AR(1)
  whitening matrix construction and dense $`O(m^2 L)`$ / $`O(L m m_0)`$
  matrix multiplications in favor of in-place column-wise bidiagonal
  whitening, tridiagonal $`W^\top W`$ stencil accumulation, L1-tiled
  banded Cholesky solves (for `nuisance_mode = "ignore"`), and
  single-pass batch normalization across the REML noise grid.

### Coefficient recovery validation

- Test fixed window-sum and calibrated-template amplitude readouts, and
  regression through the frozen inverse’s transformed event design.
  Include active/confounded nulls, off-library HRFs, matched references
  and a direct forward GLS benchmark. See
  `COEFFICIENT_RECOVERY_VALIDATION.md`. These are validation candidates;
  package algorithms and defaults are unchanged.

### Parametric recovery validation

- Add a separate simulation study with genuine neural trial variability,
  task-active RPE nulls, correlated task confounds, and paired true-HRF,
  true-covariance, noiseless and matched-grid references. Audit
  coefficient bias, coverage, false positives and power across runs and
  independent datasets at eight fixed lambdas. Native inverse operators
  and trial regression factorizations are reused. See
  `PARAMETRIC_RECOVERY_VALIDATION.md`.

### Noise-aware HRF evidence and event recovery diagnostics

- Cross-run HRF fitting now defaults to pooled grid-REML AR(1)
  covariance, conditional GLS task/shape evidence with library
  corrections, training-run compatibility, and a fixed penalty toward
  the prespecified shape. Unsupported adaptation has an explicit
  parcel/prespecified fallback. Tests are conditional on the noise/mean
  model and are not whole-brain significance or posterior confidence.
  `noise_model = "ols", evidence = "legacy"` retains the prior fit.
- The native inverse shares operators across HRF/noise strata and uses
  original acquisition times across censor gaps. AR/calibrated fits use
  observation masks during scaling and fitting by default. Nuisance
  adjustment remains explicit; optional joint projection applies to both
  convolution and observations. Unidentifiable covariance is reported
  rather than silently assumed white.
- [`diagnose_crossrun_recovery()`](https://hallquistlab.github.io/medusa/reference/diagnose_crossrun_recovery.md)
  reports conditional event amplitude gain and excess leakage relative
  to the ideal readout, accounting for natural neural window overlap,
  fractional onsets, masks, and the frozen inverse. It does not
  automatically rescale amplitudes. See `NOISE_AWARE_HRF_RECOVERY.md`
  for review, independent reference checks, repeated validation, and
  throughput benchmarks.

### Cross-run HRF estimation and spatial regularization

- [`hrf_library()`](https://hallquistlab.github.io/medusa/reference/hrf_library.md)
  constructs a small unit-sum HRF library.
  [`estimate_crossrun_hrf()`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md)
  selects HRFs using other runs only, with batched Rcpp/BLAS training
  fits, nuisance regressors, observation masks, explicit support gates,
  and a recorded parcel/prespecified fallback hierarchy. Training can
  use condition means or free event amplitudes for signed trial
  responses that would otherwise average away, without supplying RPE
  values.
- Explicit parcel pooling and optional voxel shape penalties toward
  training-derived parcel HRFs provide spatial regularization. These
  constrain HRF selection without spatially smoothing recovered neural
  activity. Strength and effect-size thresholds are user-controlled. The
  original heuristic mode does not supply significance tests.
- [`deconvolve_crossrun_hrf()`](https://hallquistlab.github.io/medusa/reference/deconvolve_crossrun_hrf.md)
  applies frozen HRFs with required fixed lambda, grouped inverse
  operators, bounded caching, block multiplication, and recorded
  per-voxel scaling. Existing NIfTI/file deconvolution defaults are
  unchanged. See `CROSSRUN_HRF_GUIDE.md` and the reproducible
  review/benchmark scripts.
- `HRF_ESTIMATION_VALIDATION.md` reports 480 independent simulated
  datasets comparing learned, canonical, and true HRFs at three fixed
  inverse penalties. Learning helps with stable shape mismatch but can
  hurt with low SNR, zero-mean signed responses, run changes, or
  neural-duration mismatch. The original support gate accepted all null
  free-event fits under colored noise; support labels are not confidence
  measures. Estimator settings were not tuned to these results.
  Reproducible scripts, compact tables, and figures are included.

### Auditable event overlap and impulse recovery

- Alignment now interpolates within each unit’s valid acquisition
  segments before aggregation. Explicit invalid values/masks block
  interpolation, and brackets wider than the acquisition TR are missing
  by default. Off-grid median/SD now summarize interpolated units; means
  with changing contributors also change.
- `overlap_policy` supports `allow`, legacy `truncate`, and whole-window
  `exclude`. NULL retains collision-column truncation when explicitly
  configured, otherwise allows overlap. `valid_col`, `neighbor_events`,
  and `max_gap` pass through file and batch alignment settings. No
  policy claims physiological isolation.
- Outputs report per-lag contributor/reason counts, acquisition-time IDs
  and epoch reuse, neighboring events, simultaneous events, window
  overlap, and run edges. Source clocks are returned as attributes; file
  alignment adds `alignment_file` and writes a `.sources.csv` clock
  sidecar beside requested outputs.
