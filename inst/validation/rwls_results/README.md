# Adaptive rWLS Temporal Frame Weighting Validation (`rwls = "auto"`)

Monte Carlo validation outputs from `inst/validation/run_rwls_decon_recovery.R` (`100` replications x `7` fMRI noise/artifact/signal regimes x `4` `rwls` modes = `2,800` end-to-end Step 1-4 fits), corresponding to **Study 14 (Pillar 14)** in `book/library/medusa_validation_suite.html`.

## Files

- `pillar14_rwls_regime_summary.csv`: Parcel and voxel trial-slope RMSE, trigger rate, Woodbury vs. Innovation Parseval routing share, single-trial correlation, Stage-3 analytical SE calibration (`SE / Emp SD`), 95% CI coverage, and null Type-I error across `R1`-`R7`.
- `pillar14_rwls_connectivity_summary.csv`: Stage-4 REPAIR (`M-MOMS`, `decon_connectivity_repair(..., shared_error = "auto")`) true coupling (`P1 -> P2`) and uncoupled null (`P1 -> P3`) bias, RMSE, and false-positive rates under synchronous cross-parcel sub-scrubbing glitches (`R7`).
- `pillar14_rwls_illustrative_traces.csv`: Single-run neural waveforms and learned frame weights `w_t` comparing sub-scrubbing motion glitches (`R2`) vs. `4.5x` heavy-tailed neural bursts + `1.1 s` HRF latency mismatch (`R5`).
- `fig14_pillar14_rwls_decon.png`: Multi-panel summary figure for Pillar 14.

## Reproducing

```sh
Rscript inst/validation/run_rwls_decon_recovery.R 100
```
