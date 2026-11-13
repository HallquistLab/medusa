# medusa

`medusa` (**M**ultilevel **E**vent-Related **D**econvolved **S**ignal **A**nalysis) is an R package for hemodynamic response estimation, regularized BOLD deconvolution, event-locked alignment, operator-based epoch de-mixing, observation-noise propagation, connectivity disattenuation, and split-wise multilevel modeling (`mixed_by`, `jlmer`, `meta_mixed_by`) of fMRI time series.

`medusa` was extracted from [`fmri.pipeline`](https://github.com/HallquistLab/fmri.pipeline) so that event-locked deconvolution and multilevel trajectory workflows can be installed and maintained independently of voxelwise FSL/AFNI/SPM GLM cluster pipelines.

## Installation

```r
install.packages("remotes")
remotes::install_github("HallquistLab/medusa")
```

## Core Workflow

1. **Cross-run HRF & AR($p$) noise estimation:** `hrf_library()`, `estimate_crossrun_hrf()`
2. **Deconvolution:** `deconvolve_crossrun_hrf()`, `deconvolve_reglin()`, `deconvolve_nlreg()`, `deconvolve_nlreg_resample()`, `voxelwise_deconvolution()`
3. **Event-locked alignment & epoch de-mixing:** `fmri_ts`, `event_lock_ts()`, `get_medusa_interpolated_ts()`, `evt_align_decon_files()`, `run_decon_alignment()`, `decon_epoch_regression()`
4. **Diagnostics & connectivity repair:** `diagnose_crossrun_recovery()`, `deconvolution_impulse_diagnostics()`, `decon_connectivity_repair()`, `medusa_between_subject_repair()`
5. **Multilevel trajectory & joint multi-lag inference:** `mixed_by()`, `jlmer()`, `meta_mixed_by()`, `medusa_fit_mmom()`, `medusa_multilag_infer()`, `backproject_medusa()`

## Documentation

MEDuSA documentation follows a three-part learning sequence:

1. [Conceptual foundations](vignettes/medusa-concepts.Rmd): convolution, regularization, HRF learning, event mixing, units, duration, and precision-weight propagation, with an annotated reading guide.
2. [Worked workflows](vignettes/medusa-workflows.Rmd): executable simulation and recipes for NIfTI inputs, file alignment, and downstream precision-weighted / ZCA-whitened multilevel models (`mixed_by`).
3. [Technical companion](medusa_deconvolution_review_prs_42_53.qmd): derivations, simulation demonstrations, implementation details, and validation limits.

See [rendering instructions](vignettes/medusa-docs.md) for producing all three HTML documents from a checkout.
