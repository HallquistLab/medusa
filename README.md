# medusa

`medusa` (**M**ultilevel **E**vent-Related **D**econvolved **S**ignal **A**nalysis) is an R package for hemodynamic response estimation, regularized BOLD deconvolution, event-locked alignment, operator-based epoch de-mixing, observation-noise propagation, connectivity disattenuation, and split-wise multilevel modeling (`mixed_by`, `jlmer`, `meta_mixed_by`) of fMRI time series.

`medusa` was extracted from [`fmri.pipeline`](https://github.com/HallquistLab/fmri.pipeline) so that event-locked deconvolution and multilevel trajectory workflows can be installed and maintained independently of voxelwise FSL/AFNI/SPM GLM cluster pipelines.

## Installation

```r
install.packages("remotes")
remotes::install_github("HallquistLab/medusa")
```

## Core Workflow

1. **NIfTI matrix I/O, cross-run HRF & AR($p$) noise estimation:** `read_nifti_ts()`, `write_nifti_ts()`, `hrf_library()`, `estimate_crossrun_hrf()`
2. **Deconvolution:** `deconvolve_crossrun_hrf()`, `deconvolve_reglin()`, `deconvolve_nlreg()`, `deconvolve_nlreg_resample()`, `voxelwise_deconvolution()`
3. **Event-locked alignment & epoch de-mixing:** `align_medusa_ts()`, `get_medusa_interpolated_ts()`, `evt_align_decon_files()`, `run_decon_alignment()`, `decon_epoch_regression()`, `fmri_ts`, `event_lock_ts()`
4. **Diagnostics & connectivity repair:** `diagnose_crossrun_recovery()`, `deconvolution_impulse_diagnostics()`, `decon_connectivity_repair()`, `medusa_between_subject_repair()`
5. **Multilevel trajectory & joint multi-lag inference:** `mixed_by()`, `jlmer()`, `meta_mixed_by()`, `medusa_fit_mmom()`, `medusa_multilag_infer()`, `backproject_medusa()`

## Documentation

MEDuSA documentation includes a comprehensive Quarto book alongside package vignettes:

1. **[MEDuSA Quarto Book](https://hallquistlab.github.io/medusa/book/)** (`book/`): End-to-end methodological guide covering conceptual foundations, a reproducible worked example, dedicated chapters for all six pipeline steps (goals, input/output architecture, recommended settings, event-informed continuous deconvolution vs. blind deconvolution, epoch de-mixing, parcel panels, group multilevel & multi-lag trajectory inference, and disattenuated connectivity), measurement-error theory, and calibration diagnostics. Render locally from `book/` with `quarto render`.
2. **[Conceptual foundations](vignettes/medusa-concepts.Rmd)**: Package vignette introducing convolution, regularization, HRF learning, event mixing, units, duration, and precision-weight propagation.
3. **[Worked workflows](vignettes/medusa-workflows.Rmd)**: Executable simulation and recipes for NIfTI inputs, file alignment, and downstream precision-weighted / ZCA-whitened multilevel models (`mixed_by`).
4. **[Technical companion](medusa_deconvolution_review_prs_42_53.qmd)**: Historical technical derivations and early validation notes.

See [rendering instructions](vignettes/medusa-docs.md) for producing the HTML vignettes and Quarto book from a checkout.
