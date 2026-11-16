# SNR-Stratified Deconvolution & Bilateral Graph-EB Epoch Validation (`lambda_mode = "auto"`, `tau2_shrinkage = "graph_eb"`)

Monte Carlo validation outputs from `inst/validation/run_track2_decon_eb_recovery.R`.

## Files

- `exp1_snr_stratified_lambda_summary.csv` & `exp1_snr_stratified_lambda_by_subject.csv`: Experiment 1 (`N = 20` subjects, `16 x 16` heterogeneous cortical slice with 10% focal biphasic core, 15% moderate periphery, and 75% silent spatially smooth `AR(2)` background) comparing `lambda_mode = "auto"` (FDR-gated SNR-stratified PAVA) against `lambda_mode = "global"` and naive unweighted GCV.
- `exp2_graph_eb_sulcal_summary.csv`: Experiment 2 (`N = 24` subjects, folded 3D cortical sulcus with CSF gap + isolated Efron-Morris hotspots under variable-ISI heteroscedasticity) comparing `tau2_shrinkage` modes (`"none"`, `"eb"`, `"graph_eb"`) paired with Step 1 `"voxel"` and `"graph"` pooling.
- `exp3_smooth_null_calibration.csv`: Experiment 3 (`N = 32` subjects, `6 x 6 x 6` spatially smooth 3D `AR(3)` pure-null volumes at `FWHM = 6 mm`) verifying subject-level and group-level Type-I error calibration and suppression of variance-underestimation `|t|` spikes.
- `fig1_track2a_snr_stratified_lambda.png` & `fig2_track2b_tau2_spatial_profile.png`: Diagnostic summary panels for Experiments 1 and 2.

## Reproducing

```sh
Rscript inst/validation/run_track2_decon_eb_recovery.R
```
