# Sparse Bilateral HRF & Noise Graph Validation (`pooling = "graph"`)

Monte Carlo and volumetric scalability validation outputs from `inst/validation/run_sparse_graph_hrf_recovery.R`.

## Files

- `scenario1_summary.csv` & `scenario1_subject_metrics.csv`: `N = 24` participants with person-specific `+/- 4-6 mm` shifted functional/neurovascular parcel borders (Kong et al., 2021), a `2 mm` out-of-mask sulcal CSF gap, and a high-`rho`/high-variance pial draining-vein border rim across 8 pooling strategies (`voxel`, `parcel`, Euclidean vs. 26-connected geodesic CSR graphs, soft parcel boundaries, and person-specific oracle parcellation).
- `scenario2_null_summary.csv` & `scenario2_null_subject_metrics.csv`: Spatially smooth 3D `AR(3)` pure-null calibration (`FWHM = 6 mm`, `N = 24`) for both `"condition"` and `"trial"` amplitude HRF screening plus Step 3 epoch-regression Type-I error control.
- `scenario3_scalability_summary.csv`: Whole-brain 3D volumetric runtime and memory scaling (`V = 2,000` to `50,000` cortical voxels) for geodesic CSR graph construction, bilateral HRF/noise pooling, `<= 48` prototype quantization, and Step 2 GLS deconvolution.
- `fig1_scenario1_accuracy.png` & `fig2_scenario3_scalability.png`: Diagnostic summary figures for Scenarios 1 and 3.
- `sparse_graph_hrf_benchmark.rds`: Bundled tabular outputs from all three scenarios.

## Reproducing

```sh
Rscript inst/validation/run_sparse_graph_hrf_recovery.R
```
