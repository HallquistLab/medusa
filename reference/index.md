# Package index

## All functions

- [`backproject_medusa()`](https://hallquistlab.github.io/medusa/reference/backproject_medusa.md)
  : Backproject parcel-level MEDuSA coefficients onto a 4D NIfTI brain
  mask

- [`compress_mts_pca()`](https://hallquistlab.github.io/medusa/reference/compress_mts_pca.md)
  : Core function for computing multivariate time series compression
  scores by principal components

- [`decon_connectivity_repair()`](https://hallquistlab.github.io/medusa/reference/decon_connectivity_repair.md)
  : Measurement-error-corrected inter-region connectivity and PPI via
  REPAIR

- [`decon_epoch_regression()`](https://hallquistlab.github.io/medusa/reference/decon_epoch_regression.md)
  : Operator-based epoch de-mixing and gain-corrected trial regression
  across lags

- [`deconvolution_impulse_diagnostics()`](https://hallquistlab.github.io/medusa/reference/deconvolution_impulse_diagnostics.md)
  : Diagnose temporal spread of a fixed linear deconvolution operator

- [`deconvolve_crossrun_hrf()`](https://hallquistlab.github.io/medusa/reference/deconvolve_crossrun_hrf.md)
  : Deconvolve runs using frozen HRFs learned from other runs

- [`deconvolve_nlreg`](https://hallquistlab.github.io/medusa/reference/deconvolve_nlreg.md)
  : C++ port of Bush and Cisler 2013, Magnetic Resonance Imaging Adapted
  from the original provided by Keith Bush as well as C++ code from
  Jiang Bian

- [`deconvolve_nlreg_resample()`](https://hallquistlab.github.io/medusa/reference/deconvolve_nlreg_resample.md)
  : This function deconvolves the BOLD signal using Bush 2011 method,
  augmented by the resampling approach of Bush 2015.

- [`deconvolve_reglin()`](https://hallquistlab.github.io/medusa/reference/deconvolve_reglin.md)
  : Regularized linear deconvolution of BOLD time series

- [`diagnose_crossrun_recovery()`](https://hallquistlab.github.io/medusa/reference/diagnose_crossrun_recovery.md)
  : Diagnose amplitude attenuation and event leakage for a frozen
  inverse

- [`estimate_crossrun_hrf()`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md)
  : Estimate HRF shapes from other runs, with optional spatial
  regularization

- [`event_lock_ts()`](https://hallquistlab.github.io/medusa/reference/event_lock_ts.md)
  : function to get interpolated event locked data

- [`evt_align_decon_files()`](https://hallquistlab.github.io/medusa/reference/evt_align_decon_files.md)
  : Align a set of deconvolved time series files to an event of interest
  function.

- [`fmri_ts`](https://hallquistlab.github.io/medusa/reference/fmri_ts.md)
  : R6 class representing a multivariate time series object for fMRI
  analysis

- [`generate_feature`](https://hallquistlab.github.io/medusa/reference/generate_feature.md)
  : This function creates K shifts of a neural events vector according
  to the kernel length, K.

- [`generate_feature_armadillo`](https://hallquistlab.github.io/medusa/reference/generate_feature_armadillo.md)
  : This function creates K shifts of a neural events vector according
  to the kernel length, K.

- [`get_medusa_compression_score()`](https://hallquistlab.github.io/medusa/reference/get_medusa_compression_score.md)
  : front-end function for taking a list of windowed time series by mask
  value, interpolating them onto a time grid, and (optionally) averaging
  across voxels/units within a value to derive the mean interpolated
  time series

- [`get_medusa_interpolated_ts()`](https://hallquistlab.github.io/medusa/reference/get_medusa_interpolated_ts.md)
  : front-end function for taking a list of windowed time series by mask
  value, interpolating them onto a time grid, and (optionally) averaging
  across voxels/units within a value to derive the mean interpolated
  time series

- [`hrf_library()`](https://hallquistlab.github.io/medusa/reference/hrf_library.md)
  : Build a small, constrained library of causal HRFs

- [`jlmer()`](https://hallquistlab.github.io/medusa/reference/jlmer.md)
  : Fit a linear mixed model with Julia MixedModels

- [`medusa_between_subject_repair()`](https://hallquistlab.github.io/medusa/reference/medusa_between_subject_repair.md)
  : Between-subject REPAIR structural regression on MEDuSA mixed models

- [`medusa_fit_mmom()`](https://hallquistlab.github.io/medusa/reference/medusa_fit_mmom.md)
  : Fit a REPAIR modified method-of-moments (M-MOM / M-MOMS)
  measurement-error model

- [`medusa_multilag_infer()`](https://hallquistlab.github.io/medusa/reference/medusa_multilag_infer.md)
  : Joint Multi-Lag Trajectory Modeling and Operator-Informed Inference
  for MEDuSA

- [`meta_mixed_by()`](https://hallquistlab.github.io/medusa/reference/meta_mixed_by.md)
  : Function to run Bayesian random-effects meta-regression on
  coefficients from mixed_by

- [`mixed_by()`](https://hallquistlab.github.io/medusa/reference/mixed_by.md)
  : Mixed by runs a set of mixed-effects models for each combination of
  a set of factors. Its primary use is to run the same model on
  different splits of the data.

- [`run_decon_alignment()`](https://hallquistlab.github.io/medusa/reference/run_decon_alignment.md)
  :

  Align a set of deconvolved time series files from
  `voxelwise_deconvolution` to a set of events in a task-related design

- [`spm_hrf()`](https://hallquistlab.github.io/medusa/reference/spm_hrf.md)
  : Canonical SPM difference-of-gammas hemodynamic response function

- [`voxelwise_deconvolution()`](https://hallquistlab.github.io/medusa/reference/voxelwise_deconvolution.md)
  : Function to perform voxelwise deconvolution on an fMRI dataset using
  the fMRI model arguments object
