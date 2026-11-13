# MEDuSA Quarto Book

This directory contains the Quarto source, standalone simulation and figure scripts, and rendered HTML output (`_book/index.html`) for **MEDuSA: Multi-level Event-related Deconvolved Signal Analysis**.

## Rendering the Book

From `book/`, run:

```bash
quarto render
```

The rendered HTML book is written to `_book/index.html` (configured in [`_quarto.yml`](_quarto.yml)). Rendering uses the saved `.rds`, `.csv`, and `.png` outputs generated from the installed `medusa` package and completes in seconds.

## Structure

- **Introduction & Conceptual Foundations:** [`index.qmd`](index.qmd), [`foundations.qmd`](foundations.qmd), [`worked-example.qmd`](worked-example.qmd) (`worked/`)
- **Step 1 · Cross-run HRF and AR($p$) noise:** [`medusa_step1_hrf_noise/medusa_step1_hrf_noise.qmd`](medusa_step1_hrf_noise/medusa_step1_hrf_noise.qmd) and [`woolrich_review/woolrich_noise_assessment.qmd`](medusa_step1_hrf_noise/woolrich_review/woolrich_noise_assessment.qmd)
- **Step 2 · Event-informed continuous deconvolution:** [`medusa_step2_deconvolution/medusa_step2_deconvolution.qmd`](medusa_step2_deconvolution/medusa_step2_deconvolution.qmd)
- **Step 3 · Peri-event epoch readout and demixing:** [`medusa_step3_epochs/medusa_step3_epochs.qmd`](medusa_step3_epochs/medusa_step3_epochs.qmd)
- **Step 4 · Parcel panels and cohort precision weights:** [`medusa_step4_parcel_panels/medusa_step4_parcel_panels.qmd`](medusa_step4_parcel_panels/medusa_step4_parcel_panels.qmd)
- **Step 5 · Group mixed models and multi-lag trajectory inference:** [`medusa_step5_group_models/medusa_step5_group_models.qmd`](medusa_step5_group_models/medusa_step5_group_models.qmd) and [`medusa_pr62_calibration/medusa_pr62_calibration.qmd`](medusa_pr62_calibration/medusa_pr62_calibration.qmd)
- **Step 6 · Trial-level connectivity, M-MOMS, and PPI:** [`medusa_step6_connectivity/medusa_step6_connectivity.qmd`](medusa_step6_connectivity/medusa_step6_connectivity.qmd), [`medusa_standard_errors_mmom.qmd`](medusa_standard_errors_mmom.qmd), and [`medusa_connectivity_inference/medusa_connectivity_inference.qmd`](medusa_connectivity_inference/medusa_connectivity_inference.qmd)
- **Validation & Reproduction:** [`evidence.qmd`](evidence.qmd), [`reference.qmd`](reference.qmd), and [`library/medusa_validation_suite.html`](library/medusa_validation_suite.html)

See [`reference.qmd`](reference.qmd) for one-line `Rscript` commands to regenerate every chapter's simulations and figures from `library(medusa)`.
