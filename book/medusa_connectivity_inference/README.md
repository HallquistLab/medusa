# Latent Connectivity Beyond Per-Trial Disattenuation

Companion simulation and derivation directory for [`medusa_connectivity_inference.qmd`](medusa_connectivity_inference.qmd) in Part IV (**Methodological deep dives**) of the MEDuSA Quarto book.

While `medusa::decon_connectivity_repair()` ([Step 6](../medusa_step6_connectivity/medusa_step6_connectivity.qmd)) implements per-trial M-MOMS predictor disattenuation, cross-parcel `shared_error = "auto"` subtraction, and HKSJ participant random-effects $t$-intervals, the scripts in this directory evaluate full cross-trial corrected-score estimating equations, finite-basis likelihood-only readouts ($B_g G_g = I$), and joint Gaussian resolution-retaining models across 400-replicate paired stress benchmarks ($K_S \neq K_T$, serial/shared Gaussian and $t_5$ errors, and small-$S$ cohorts).

## Reproduce

```sh
Rscript likelihood_only_bridge.R .
Rscript joint_gaussian_demo.R .
MEDUSA_REPS=400 Rscript connectivity_prototype.R .
Rscript prepare_figures.R .
quarto render medusa_connectivity_inference.qmd
```

- `measurement_error_models.R`: research helper functions for full cross-trial corrected-score estimation.
- `connectivity_prototype.R`: fixed-seed 400-replicate paired simulation driver (`simulation_summary.csv`, `simulation_draws.csv`).
- `likelihood_only_bridge.R` & `joint_gaussian_demo.R`: numerical audits (`likelihood_bridge_checks.csv`, `joint_gaussian_example.csv`).
- `prepare_figures.R`: generates `bias.png` and `coverage.png`.
