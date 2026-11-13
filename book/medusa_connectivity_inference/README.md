# MEDuSA connectivity/PPI inference exploration

Start with `medusa_connectivity_inference.html`. This is a separate research exploration, not a production MEDuSA API or a validated extension of REPAIR. Earlier primers and reports are preserved.

The proposed route is a finite-basis likelihood-only measurement model, full seed/target error covariance, and corrected participant-cluster scores. A fitted joint Gaussian latent-source/PPI example demonstrates an alternative retaining resolution. The simulations assume known measurement covariance; no full fMRI coverage guarantee is claimed.

## Reproduce

Base R runs the prototypes. Figure/report generation additionally needs ggplot2, scales, ragg, knitr, and Quarto.

```sh
Rscript likelihood_only_bridge.R .
Rscript joint_gaussian_demo.R .
MEDUSA_REPS=400 Rscript connectivity_prototype.R .
Rscript prepare_figures.R .
quarto render medusa_connectivity_inference.qmd
```

`measurement_error_models.R` contains reusable research functions without running the simulations. `connectivity_prototype.R` is the complete fixed-seed simulation driver, retaining its own exact helper definitions for reproducibility. The main experiment has 400 datasets in each of six scenarios, 24 trials per participant, and 60 participants except for the 12-participant scenario. Different methods on a dataset are paired. The 29,600 saved coefficient evaluations are not 29,600 independent datasets.

- `simulation_draws.csv`: all coefficient evaluations, including usability and coverage.
- `simulation_summary.csv`: all 74 scenario/method/effect summaries.
- `figure_data.csv`: plotted values and Monte Carlo intervals.
- `simulation_settings.rds`: settings and R session.
- `prototype_checks.txt`, `likelihood_bridge_checks.csv`, `joint_gaussian_checks.txt`: independent algebra/numerical checks.
- `source_manifest.json`: inspected MEDuSA revision and source hashes.

The group prototype targets common connectivity/PPI coefficients and uses normal critical values with a simple G/(G-1) score-covariance adjustment. It is neither a random-slope population-mean estimator nor a calibrated small-cluster procedure. Unknown error covariance, misspecified source means, informative censoring, and multiplicity remain separate validation requirements.
