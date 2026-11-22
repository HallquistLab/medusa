# medusa 0.1.2

## Native heap-safety & LAPACK divide-and-conquer (`dgesdd` / `dsyevd`) hardening

* Eliminated native C heap corruption and `SIGABRT`/`SIGSEGV` crashes on Linux
  HPC systems running unpatched OpenBLAS / Reference-LAPACK (`< 3.12.1`,
  LAPACK Issues `#672` / `#1138`), where divide-and-conquer SVD (`dgesdd` /
  `DLASD2` $\to$ `DLASD3`) and symmetric eigendecomposition (`dsyevd` / `DLAED2`)
  can write out-of-bounds on the C heap when multi-trial Toeplitz designs or
  AR(1)-whitened sparse Woodbury outlier-spike matrices have tightly clustered
  singular/eigenvalues:
  - Added C++ helper `medusa_svd_cpp()` (`src/crossrun_hrf.cpp`) using LAPACK
    QR-iteration (`dgesvd`, `"std"`) with a polar-reorthogonalized symmetric QR
    (`dsyev`) fallback, and routed `hrf_safe_svd()` (`R/hrf_noise_evidence.R`)
    directly through `medusa_svd_cpp()` so `base::svd()` / `La.svd()` (`dgesdd`)
    is never invoked.
  - Updated `hrf_safe_svd_econ()` and `hrf_safe_svd_full()`
    (`src/crossrun_hrf.cpp`) to use `"std"` (`dgesvd` / `dsyev`) exclusively and
    routed all remaining Step-3 Woodbury whitening, fallback task-projection,
    response-matching, and time-series compression SVDs
    (`R/diagnose_crossrun_recovery.R`, `R/poet_covreg.R`,
    `R/event_lock_timeseries.R`) through `hrf_safe_svd()`.
  - Added `poet_safe_eig_sym()` (`src/poet_covreg.cpp`) using `"std"` (`dsyev`)
    with input symmetrization/finiteness sanitization and a guaranteed $n \times n$
    diagonal fallback across all 8 symmetric eigendecompositions in the POET
    shared-error engine, preventing unchecked `0 x 0` matrix indexing under
    `ARMA_NO_DEBUG` if a LAPACK decomposition ever fails.
  - Added `arma::solve_opts::no_approx` (with explicit `"std"` spectral
    pseudoinverse fallbacks) in `damped_sym_inv()` (`src/poet_covreg.cpp`),
    `solve_sympd_fallback()` (`src/deconvolve_reglin.cpp`), and Markov boundary
    triangular solves (`src/crossrun_hrf.cpp`) so near-singular linear systems do
    not implicitly trigger Armadillo's `svd_dc_econ` (`dgesdd`) fallback.
* Hardened C++ SEXP lifetime and CSR bounds safety:
  - Retained `Rcpp::NumericMatrix` storage handles inside `CompactStratumFac`
    (`src/poet_covreg.cpp`) with upfront row/column dimension checks so raw
    factor pointers remain GC-protected across any coerced R inputs.
  - Added upfront CSR pointer, neighborhood distance, kernel similarity, and
    evidence-gate dimension checks in `csr_bilateral_noise_smooth_cpp()`,
    `csr_cluster_member_weights_cpp()`, `csr_quad_cor_cpp()`, and
    `csr_hrf_graph_pool_cpp()` (`src/spatial_graph.cpp`).

# medusa 0.1.0

## Subspace-debiased Adaptive-POET + covariance-regression shared-error engine & connectivity inference

* Added the `v3` Subspace-Debiased Adaptive-POET + Multi-Lag Toeplitz &
  General-$W_s$ Covariance-Regression shared-error engine
  (`R/poet_covreg.R`, `src/poet_covreg.cpp`, `shared_error_method = "poet_covreg"`
  in `decon_epoch_regression()`):
  - **Bold-free Step 2 $\to$ Step 3 handoff:** `deconvolve_crossrun_hrf()` now
    caches uncensored `bold_obs` up to 16,384 voxels and switches automatically
    above 16,384 voxels to a compact `m x 64` `whitened_sketch`, while
    `medusa_fit_shared_error_covreg()` switches automatically above 512 units to
    a voxel-safe dual-space factor mode (`mode = "lowrank_voxel"`) that never
    allocates $V \times V$ matrices.
  - **Subspace-debiased stationary & dynamic error moments:** Projects out
    convolved trial/duration/first-difference designs ($M_s$) with exact
    $\operatorname{tr}(M_{s_g} M_{s_h})$ degrees-of-freedom debiasing, Tukey-Hanning
    tapered multi-lag Toeplitz PSD projection, run-scale Empirical-Bayes
    harmonization across heterogeneous-noise runs, and LTZ-gated dynamic task,
    covariate, Step-2 rWLS, and external `variance_confounds` variance
    envelopes with exact Schur-complement restoration when confound columns are
    gated out.
  - **Operator & bare-panel parity (`decon_connectivity_repair()`):** Propagates
    full-rank Block-PAF shared-error embeddings (`shared_noise_vec`), cross-lag
    embeddings (`lag_noise_vec`), and cross-trial error correlations
    (`trial_cor_vec`) onto `panel` (including automatic CSV/`.csv.gz`
    pipe-delimited list-column round-trip restoration across
    `medusa_write_table()`, `medusa_read_table()`, `data.table::fread()`,
    `decon_connectivity_repair()`, and `aggregate_decon_panel()`), and adds
    single-subject ($N = 1$) & group hybrid operator + gated onset-kernel
    cross-trial sandwich SEs and Satterthwaite degrees of freedom
    (`crosstrial_se = TRUE`).
  - **Post-hoc parcel aggregation & `medusa_confounds` integration:** Added
    exported `aggregate_decon_panel()` to collapse voxel-level panels or
    `decon_epoch_regression()` fits to parcels with exact shared-error
    propagation, and updated `read_fmriprep_confounds()` and
    `decon_epoch_regression(..., variance_confounds = ...)` to extract and
    accept `"medusa_confounds"` (and `"medusa_nifti_ts"`) objects directly.

## Fast C++ NIfTI matrix I/O and matrix-native event alignment

* Added `read_nifti_ts()` and `write_nifti_ts()` for streaming C++ conversion
  between 3D/4D `.nii`/`.nii.gz` images and `time x voxels` (`T x V`) numeric
  matrices (`"medusa_nifti_ts"` bundles carrying `voxel_map`, `parcel_ids`,
  `tr`, and validated `sform`/`qform` `geometry` without allocating full 4D R
  arrays). Updated `voxelwise_deconvolution()` and `backproject_medusa()` to
  use the streaming C++ NIfTI path.
* Added `align_medusa_ts()` (`src/epoch_align.cpp`) for C++ matrix-native
  event-locked alignment and multi-kernel readout quadrature (`"linear"`,
  `"boxcar"`, `"triangle"`, `"gaussian"`) across `T x V` matrices,
  `"medusa_nifti_ts"` bundles, `deconvolve_crossrun_hrf()` outputs, multi-run
  lists, and multivariate signal lists while preserving all `overlap_policy`,
  censoring, `max_gap`, `epoch_status`, `align_*`, `time_audit`, and
  `source_clock` invariants (`~80x–130x` faster and `~10x–135x` lower peak R
  heap memory than legacy long-table workflows).
* Accelerated `get_medusa_interpolated_ts()` and `get_medusa_compression_score()`
  on shared-clock `fmri_ts` objects via a zero-copy C++ matrix view (with
  automatic fallback for ragged per-unit clocks), added `fmri_ts$to_matrix()`,
  and upgraded `evt_align_decon_files()` with `fast_pivot_decon_cpp()`, `.rds`
  matrix-bundle support, and fork parallelism on Unix CLI/HPC.
* Modernized file and tabular serialization across the pipeline:
  `voxelwise_deconvolution()` now supports `output_format = c("csv.gz", "rds", "nifti", "parquet")`
  (eliminating intermediate `mat2df()` + `dplyr::left_join()` expansions and
  writing delimited outputs via multi-threaded `data.table::fwrite()`);
  `run_decon_alignment()` and `evt_align_decon_files()` discover and align all
  four deconvolution formats, skip unused coordinate/status columns when reading
  long tables, and write aligned panels to `.csv.gz`, `.parquet`, or `.rds`; and
  `mixed_by()` now reads `.parquet` files and Hive-partitioned Parquet
  directories (via optional `arrow`) as well as `.rds` and delimited `external_df`
  tables, with automatic column projection and optional `split_filter` pushdown.
* Unified in-memory NIfTI, cross-run HRF, deconvolution, and epoch-regression
  ergonomics: `estimate_crossrun_hrf()` and `deconvolve_crossrun_hrf()` now
  accept `"medusa_nifti_ts"` objects directly (auto-populating `tr`,
  `spatial_groups = ts$parcel_ids`, `kernel_grid = hrf_library(tr, preset = "medusa")`,
  and `spatial_coords` when `spatial_fwhm` is requested), split combined
  single-participant multi-run `events` `data.frame`s automatically across BIDS/MEDuSA
  run labels, and propagate `voxel_map` and `geometry` onto `"crossrun_hrf"` and
  `"crossrun_decon"` so `write_nifti_ts(decon, ...)` and
  `backproject_medusa(..., brain_mask = decon)` need no extra reference file;
  `decon_epoch_regression()` now accepts column-name `covariates` (e.g.,
  `c("rpe", "abs_rpe")`), one-sided formulas (`~ rpe + abs_rpe`), column-name
  `trial_durations` (`"rt"`), and combined multi-run `events` tables.
* Added zero-dependency BIDS, fMRIPrep, and `BrainGnomes` ingestion helpers:
  `read_bids_events()` normalizes `*_events.tsv` files (`trial_type` to
  `condition`, BIDS `"n/a"`/`"."`/`"NaN"` missing tokens without coercing
  numeric behavioral columns, `cohort-`/`hemi-` entity validation, and run
  ordering matched to `bold`); `read_fmriprep_confounds()` converts
  `*_desc-confounds_timeseries.tsv` and `*_desc-<desc>_confounds.tsv` tables
  into clock-preserving `valid_rows` masks (`fd_threshold`, `dvars_threshold`,
  optional zero-phase respiratory notch `motion_bandstop_bpm` or Butterworth
  low-pass `motion_lowpass_hz` filtering of the 6 rigid-body motion parameters
  prior to Power 2012 FD calculation, `framewise_displacement_filtered` / Power
  2012 FD fallback, external `*_censor.1D` / `*_scrub.tsv` sidecars,
  `non_steady_state_outlier*`, spike dilation, `min_segment_length`) and
  strictly finite `nuisance` matrices
  (`"24p"`, `"12p"`, `"6p"`, `"csf_wm"`, `"csf_wm_8p"`, `"24p_csf_wm"`,
  `"27p"`, `"36p"`, `"compcor"`, `"custom"` with `<0-k>` range/shorthand
  expansion, `"none"`) without deleting censored TRs; `read_fmriprep_aroma()`
  and `apply_aroma()` ingest `fMRIPost-AROMA` (`*_desc-melodic_mixing.tsv` +
  `*_desc-aroma_metrics.tsv`) and legacy fMRIPrep (`*_desc-MELODIC_mixing.tsv`
  + `*_AROMAnoiseICs.csv`) components to perform in-memory non-aggressive
  Pruim et al. (2015) partial regression on BOLD and matched nuisance
  regressors (preventing Lindquist et al., 2019 artifact re-introduction while
  preserving task-correlated neural amplitudes and temporal degrees of freedom);
  and `read_fmriprep_subject()` discovers raw fMRIPrep (`desc = "preproc"`,
  optional `aroma = "nonaggressive"` / `"auto"` and `aroma_dir`) or
  `BrainGnomes` postprocessed (`desc = "postproc"`, `*_templatemask.nii.gz`,
  `*_censor.1D`, provenance checks blocking physically scrubbed BOLD,
  duplicate AROMA regression, and double confound regression) BOLD runs, masks,
  confounds, and BIDS events into an enriched `"medusa_nifti_ts"` bundle
  accepted directly by `estimate_crossrun_hrf()`, `deconvolve_crossrun_hrf()`,
  and `align_medusa_ts()`.
* Added sparse bilateral HRF and noise graph pooling (`build_spatial_graph()`
  and `estimate_crossrun_hrf(..., pooling = "graph")`, while preserving
  `pooling = "voxel"` as the default for atlas-free NIfTI inputs): constructs
  `O(V)` CSR neighborhoods via 26-connected in-mask geodesic lattice Dijkstra
  (`metric = "geodesic"`, blocking cross-sulcal CSF/out-of-mask bleed), 3D
  spatial-hash radius/`k`-NN (`metric = "euclidean"`), or custom surface/voxel
  edge graphs; performs bilateral REML noise profile smoothing with adaptive
  stratum quantization (`max_noise_strata = 10L`) and empirical innovation
  spatial correlation forms (`v_quad`); pools out-of-sample training-run HRF
  evidence using Spearman-disattenuated cross-run profile concordance, signed
  anti-correlation boundary stopping, soft/hard/ignored atlas boundary priors
  (`parcel_boundary = c("soft", "hard", "ignore")`, `parcel_cross_weight = 0.35`),
  and prototype waveform quantization (`max_hrf_strata = 48L`) for whole-brain
  scalability; and hardens multi-trial GLS SVD projections (`hrf_safe_svd()`).
* Added SNR-stratified deconvolution regularization and bilateral graph
  Empirical-Bayes variance stabilization in Steps 2–3 (`deconvolve_crossrun_hrf()`
  and `decon_epoch_regression()`): `deconvolve_crossrun_hrf()` gains
  `lambda_mode = c("auto", "global", "stratified")`, `max_lambda_strata = 4L`,
  and per-voxel vector `lambda` support, stratifying prewhitened GCV
  regularization across spatially stabilized out-of-sample task-$R^2$ tiers
  (gated by cross-run Fisher + BH FDR significance with weighted PAVA
  monotonicity) and propagating heterogeneous `lambda` through
  `diagnose_crossrun_recovery()` and `decon_epoch_regression()`; and
  `decon_epoch_regression()` gains `tau2_shrinkage = c("auto", "graph_eb", "eb", "none")`
  and `spatial_graph = NULL` (`src/spatial_graph.cpp`), shrinking scale-free
  untruncated trial biological variance ratios prior to zero truncation with
  Isserlis spatial error-correlation accounting, bilateral boundary stopping,
  Efron-Morris hotspot protection, and spatial Limma/Smyth moderation of
  sandwich residual scales (`[scale_unw, scale_pw, scale_gls]`) without ever
  smoothing neural timecourses or trial regression slopes.
* Polished spatial/parcel pooling ergonomics and volumetric backprojection:
  `build_spatial_graph()` now accepts `stats::dist`, dense `Matrix`, and sparse
  symmetric `Matrix` distance objects directly alongside explicit
  `x_type = c("auto", "coords", "distance")` disambiguation;
  `deconvolve_crossrun_hrf()` exposes named per-voxel `settings$lambda_by_voxel`;
  `decon_epoch_regression(..., aggregate_by_parcel = TRUE, tau2_shrinkage = "graph_eb")`
  auto-constructs a precision-weighted 3D parcel-centroid spatial graph and
  propagates parcel-aggregated task-$R^2$ gating;
  `backproject_medusa()` auto-detects `parcel`, `evt_time`/`lag`, and all-`NA`
  `term` defaults while parsing character/factor numeric parcel IDs safely;
  and `book/medusa_spatial_parcel_pooling/` adds a full cortical-ribbon
  `ggbrain` walkthrough.

## Standalone package extraction

* Extracted all MEDuSA deconvolution, cross-run HRF estimation, event-locked
  interpolation, operator-based epoch regression, connectivity repair, NIfTI
  backprojection (`backproject_medusa()` migrated to `RNifti`), and split-wise
  multilevel modeling (`mixed_by()`, `jlmer()`, `meta_mixed_by()`) tools from
  `fmri.pipeline` into the standalone `medusa` package.

## MEDuSA learning sequence, precision-weighted epoch regression, and workflow documentation

* Add conceptual (`vignettes/medusa-concepts.Rmd`) and executable workflow
  (`vignettes/medusa-workflows.Rmd`) vignettes covering cross-run HRF
  learning, noise scaling, operator-based epoch de-mixing, duration models,
  per-trial precision weights, within-run ZCA whitening, NIfTI inputs, and the
  handoff to multilevel analysis (`mixed_by(..., weights = "precision_weight")`).
* Revise the technical companion (`medusa_deconvolution_review_prs_42_53.qmd`)
  through PR #55+, align descriptions with the current implementation,
  distinguish legacy benchmark tables from live experiments, and add thorough
  inline documentation across all MEDuSA R functions.
* `decon_epoch_regression()` now propagates exact run-specific AR(1) observation
  noise covariance through the readout and unmixing operators to compute
  per-trial standard errors (`se_raw`, `se_demixed`), lag-1 trial error
  autocorrelation (`error_ac1`), Empirical-Bayes stabilized precision weights
  (`precision_weight`), ZCA-whitened GLS slopes (`beta_gls`, `se_gls`),
  within-trial multi-lag de-mixing (`demix_within_trial = TRUE`), and optional
  within-run ZCA-whitened panel columns (`whiten_panel = TRUE`).
* `mixed_by()` now accepts `weights = "precision_weight"` (normalizing positive
  precision weights within each split and forwarding to `lme4`, `rstanarm`, and
  `jlmer`).
* `hrf_library()` adds `preset = c("custom", "standard", "extended", "dense")`;
  `estimate_crossrun_hrf()` adds Galwey spectral effective candidate multiplicity
  (`hrf_effective_candidates()`) and Cauchy combination parcel pooling
  (`hrf_cauchy_p()`); and `deconvolve_crossrun_hrf()` supports automatic
  prewhitened GCV tuning when `lambda = NULL`.

## Banded linear deconvolution solver and exact Takahashi trace

* `deconvolve_reglin_cpp()` now exploits the strict banded structure of
  $H_k^\top H_k + \lambda P + \delta I$ (half-bandwidth $b = \max(\text{kernel\_len} - 1 + \mathbb{I}(\text{prewhiten}), \; b_P)$),
  replacing dense $O(n_{\text{latent}}^3)$ normal-equation solves and dense
  matrix-RHS trace evaluations with an in-place $O(n_{\text{latent}} b^2)$
  lower-banded Cholesky factorization, contiguous multi-signal forward/backward
  substitution, and the exact $O(n_{\text{latent}} b^2)$ Takahashi banded-inverse
  recursion for effective degrees of freedom $\operatorname{tr}(A(\lambda)^{-1} H_k^\top H_k)$.
* `hrf_inverse_operator_cpp()`, `hrf_score_run_cpp()`, and `hrf_noise_profile_cpp()`
  now eliminate dense $m \times m$ AR(1) whitening matrix construction and
  dense $O(m^2 L)$ / $O(L m m_0)$ matrix multiplications in favor of
  in-place column-wise bidiagonal whitening, tridiagonal $W^\top W$ stencil
  accumulation, L1-tiled banded Cholesky solves (for `nuisance_mode = "ignore"`),
  and single-pass batch normalization across the REML noise grid.

## Coefficient recovery validation

* Test fixed window-sum and calibrated-template amplitude readouts, and
  regression through the frozen inverse's transformed event design. Include
  active/confounded nulls, off-library HRFs, matched references and a direct
  forward GLS benchmark. See `COEFFICIENT_RECOVERY_VALIDATION.md`. These are
  validation candidates; package algorithms and defaults are unchanged.

## Parametric recovery validation

* Add a separate simulation study with genuine neural trial variability,
  task-active RPE nulls, correlated task confounds, and paired true-HRF,
  true-covariance, noiseless and matched-grid references. Audit coefficient bias,
  coverage, false positives and power across runs and independent datasets at
  eight fixed lambdas. Native inverse operators and trial regression
  factorizations are reused. See `PARAMETRIC_RECOVERY_VALIDATION.md`.

## Noise-aware HRF evidence and event recovery diagnostics

* Cross-run HRF fitting now defaults to pooled grid-REML AR(1) covariance,
  conditional GLS task/shape evidence with library corrections, training-run
  compatibility, and a fixed penalty toward the prespecified shape. Unsupported
  adaptation has an explicit parcel/prespecified fallback. Tests are conditional
  on the noise/mean model and are not whole-brain significance or posterior
  confidence. `noise_model = "ols", evidence = "legacy"` retains the prior fit.
* The native inverse shares operators across HRF/noise strata and uses original
  acquisition times across censor gaps. AR/calibrated fits use observation masks
  during scaling and fitting by default. Nuisance adjustment remains explicit;
  optional joint projection applies to both convolution and observations.
  Unidentifiable covariance is reported rather than silently assumed white.
* `diagnose_crossrun_recovery()` reports conditional event amplitude gain and
  excess leakage relative to the ideal readout, accounting for natural neural
  window overlap, fractional onsets, masks, and the frozen inverse. It does not
  automatically rescale amplitudes. See `NOISE_AWARE_HRF_RECOVERY.md` for review,
  independent reference checks, repeated validation, and throughput benchmarks.

## Cross-run HRF estimation and spatial regularization

* `hrf_library()` constructs a small unit-sum HRF library.
  `estimate_crossrun_hrf()` selects HRFs using other runs only, with batched
  Rcpp/BLAS training fits, nuisance regressors, observation masks, explicit
  support gates, and a recorded parcel/prespecified fallback hierarchy.
  Training can use condition means or free event amplitudes for signed trial
  responses that would otherwise average away, without supplying RPE values.
* Explicit parcel pooling and optional voxel shape penalties toward
  training-derived parcel HRFs provide spatial regularization. These constrain
  HRF selection without spatially smoothing recovered neural activity. Strength
  and effect-size thresholds are user-controlled. The original heuristic mode
  does not supply significance tests.
* `deconvolve_crossrun_hrf()` applies frozen HRFs with required fixed lambda,
  grouped inverse operators, bounded caching, block multiplication, and recorded
  per-voxel scaling. Existing NIfTI/file deconvolution defaults are unchanged.
  See `CROSSRUN_HRF_GUIDE.md` and the reproducible review/benchmark scripts.
* `HRF_ESTIMATION_VALIDATION.md` reports 480 independent simulated datasets
  comparing learned, canonical, and true HRFs at three fixed inverse penalties.
  Learning helps with stable shape mismatch but can hurt with low SNR,
  zero-mean signed responses, run changes, or neural-duration mismatch. The
  original support gate accepted all null free-event fits under colored noise;
  support labels are not confidence measures. Estimator settings were not tuned
  to these results. Reproducible scripts, compact tables, and figures are included.

## Auditable event overlap and impulse recovery

* Alignment now interpolates within each unit's valid acquisition segments before
  aggregation. Explicit invalid values/masks block interpolation, and brackets
  wider than the acquisition TR are missing by default. Off-grid median/SD now
  summarize interpolated units; means with changing contributors also change.
* `overlap_policy` supports `allow`, legacy `truncate`, and whole-window `exclude`.
  NULL retains collision-column truncation when explicitly configured, otherwise
  allows overlap. `valid_col`, `neighbor_events`, and `max_gap` pass through file
  and batch alignment settings. No policy claims physiological isolation.
* Outputs report per-lag contributor/reason counts, acquisition-time IDs and epoch
  reuse, neighboring events, simultaneous events, window overlap, and run edges.
  Source clocks are returned as attributes; file alignment adds `alignment_file`
  and writes a `.sources.csv` clock sidecar beside requested outputs.
