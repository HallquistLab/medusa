# medusa

`medusa` (**M**ultilevel **E**vent-Related **D**econvolved **S**ignal
**A**nalysis) is an R package for hemodynamic response estimation,
regularized BOLD deconvolution, event-locked alignment, operator-based
epoch de-mixing, observation-noise propagation, connectivity
disattenuation, and split-wise multilevel modeling (`mixed_by`, `jlmer`,
`meta_mixed_by`) of fMRI time series.

`medusa` was extracted from
[`fmri.pipeline`](https://github.com/HallquistLab/fmri.pipeline) so that
event-locked deconvolution and multilevel trajectory workflows can be
installed and maintained independently of voxelwise FSL/AFNI/SPM GLM
cluster pipelines.

## Installation

``` r

install.packages("remotes")
remotes::install_github("HallquistLab/medusa")
```

## Core Workflow

1.  **NIfTI matrix I/O, cross-run HRF & AR($`p`$) noise estimation:**
    [`read_nifti_ts()`](https://hallquistlab.github.io/medusa/reference/read_nifti_ts.md),
    [`write_nifti_ts()`](https://hallquistlab.github.io/medusa/reference/write_nifti_ts.md),
    [`hrf_library()`](https://hallquistlab.github.io/medusa/reference/hrf_library.md),
    [`estimate_crossrun_hrf()`](https://hallquistlab.github.io/medusa/reference/estimate_crossrun_hrf.md)
2.  **Deconvolution:**
    [`deconvolve_crossrun_hrf()`](https://hallquistlab.github.io/medusa/reference/deconvolve_crossrun_hrf.md),
    [`deconvolve_reglin()`](https://hallquistlab.github.io/medusa/reference/deconvolve_reglin.md),
    [`deconvolve_nlreg()`](https://hallquistlab.github.io/medusa/reference/deconvolve_nlreg.md),
    [`deconvolve_nlreg_resample()`](https://hallquistlab.github.io/medusa/reference/deconvolve_nlreg_resample.md),
    [`voxelwise_deconvolution()`](https://hallquistlab.github.io/medusa/reference/voxelwise_deconvolution.md)
3.  **Event-locked alignment & epoch de-mixing:**
    [`align_medusa_ts()`](https://hallquistlab.github.io/medusa/reference/align_medusa_ts.md),
    [`get_medusa_interpolated_ts()`](https://hallquistlab.github.io/medusa/reference/get_medusa_interpolated_ts.md),
    [`evt_align_decon_files()`](https://hallquistlab.github.io/medusa/reference/evt_align_decon_files.md),
    [`run_decon_alignment()`](https://hallquistlab.github.io/medusa/reference/run_decon_alignment.md),
    [`decon_epoch_regression()`](https://hallquistlab.github.io/medusa/reference/decon_epoch_regression.md),
    `fmri_ts`,
    [`event_lock_ts()`](https://hallquistlab.github.io/medusa/reference/event_lock_ts.md)
4.  **Diagnostics & connectivity repair:**
    [`diagnose_crossrun_recovery()`](https://hallquistlab.github.io/medusa/reference/diagnose_crossrun_recovery.md),
    [`deconvolution_impulse_diagnostics()`](https://hallquistlab.github.io/medusa/reference/deconvolution_impulse_diagnostics.md),
    [`decon_connectivity_repair()`](https://hallquistlab.github.io/medusa/reference/decon_connectivity_repair.md),
    [`medusa_between_subject_repair()`](https://hallquistlab.github.io/medusa/reference/medusa_between_subject_repair.md)
5.  **Multilevel trajectory & joint multi-lag inference:**
    [`mixed_by()`](https://hallquistlab.github.io/medusa/reference/mixed_by.md),
    [`jlmer()`](https://hallquistlab.github.io/medusa/reference/jlmer.md),
    [`meta_mixed_by()`](https://hallquistlab.github.io/medusa/reference/meta_mixed_by.md),
    [`medusa_fit_mmom()`](https://hallquistlab.github.io/medusa/reference/medusa_fit_mmom.md),
    [`medusa_multilag_infer()`](https://hallquistlab.github.io/medusa/reference/medusa_multilag_infer.md),
    [`backproject_medusa()`](https://hallquistlab.github.io/medusa/reference/backproject_medusa.md)

## Documentation

MEDuSA documentation includes a comprehensive Quarto book alongside
package vignettes:

1.  **[MEDuSA Quarto
    Book](https://hallquistlab.github.io/medusa/book/)** (`book/`):
    End-to-end methodological guide covering conceptual foundations, a
    reproducible worked example, dedicated chapters for all six pipeline
    steps (goals, input/output architecture, recommended settings,
    event-informed continuous deconvolution vs. blind deconvolution,
    epoch de-mixing, parcel panels, group multilevel & multi-lag
    trajectory inference, and disattenuated connectivity),
    measurement-error theory, and calibration diagnostics. Render
    locally from `book/` with `quarto render`.
2.  **[Conceptual
    foundations](https://hallquistlab.github.io/medusa/vignettes/medusa-concepts.Rmd)**:
    Package vignette introducing convolution, regularization, HRF
    learning, event mixing, units, duration, and precision-weight
    propagation.
3.  **[Worked
    workflows](https://hallquistlab.github.io/medusa/vignettes/medusa-workflows.Rmd)**:
    Executable simulation and recipes for NIfTI inputs, file alignment,
    and downstream precision-weighted / ZCA-whitened multilevel models
    (`mixed_by`).
4.  **[Technical
    companion](https://hallquistlab.github.io/medusa/medusa_deconvolution_review_prs_42_53.qmd)**:
    Historical technical derivations and early validation notes.

See [rendering
instructions](https://hallquistlab.github.io/medusa/vignettes/medusa-docs.md)
for producing the HTML vignettes and Quarto book from a checkout.
