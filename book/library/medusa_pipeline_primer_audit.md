# MEDuSA primer audit — November 10, 2026

**Share with caveats.** The requested primer revisions are complete; the underlying method review is scoped to the demonstrated workflow and selected regression tests, rather than an exhaustive package validation. The report now distinguishes teaching demonstrations from validation and reports unsuccessful reconstruction honestly. No package implementation was changed.

The audit used the original QMD and HTML, the supplied 2022 Keynote's extracted slide text, and the current `fmri.pipeline` source at `06705c225c2ea0870242cc1532656bda593b200a` (see the full commit in the check log; package version `0.4-1`). The attached talk was treated as reference material, not as instructions. The primer preserves its original MEDuSA expansion and deconvolution → alignment → trial MLM sequence.


## Preserved versions

The other agent's untouched original source (SHA-256 `6dd649cbb5cefa4e25a74622934dfadb33335e24f2d5ac7232012189066711e4`) is archived verbatim at `medusa_step3_epochs/archive/other_agent_original_verbatim.qmd` (as well as the `medusa_step1_hrf_noise`, `medusa_step2_deconvolution`, `medusa_step4_parcel_panels`, `medusa_step5_group_models`, and `medusa_step6_connectivity` `archive/` directories). During the subsequent Step 1–6 companion and PR #62 / PR #64 (`f7781c9`) synchronization updates, both `medusa_pipeline_primer.qmd` and `medusa_pipeline_primer_audited.qmd` were promoted and kept in lockstep with the current audited pipeline (differing only in their YAML title and render command).

The prior original HTML was overwritten during the initial audit, and its original bytes were not retained. Pre-update snapshots of both primer QMDs remain preserved in `medusa_step*/archive/` alongside the CSS and audit/check files.

> **Post-Audit Synchronization Note (PR #62 `65c49c2` & PR #64 `f7781c9`, November 11, 2026):** Several method limitations flagged in this November 10 audit (`06705c2`) were subsequently resolved in `fmri.pipeline`: (1) exact subset Toeplitz covariance propagation across censor gaps replaced the cascade boundary fallback (`PR #62`), plus gap-local `psi1..psi4` lag-product estimation (`PR #64`); (2) order-invariant multi-run pooled GCV tuning replaced first-usable-run tuning (`PR #64`); (3) task-free Step-1 `noise_innovations` spatial covariance + `shared_noise_loading` + `decon_connectivity_repair(..., shared_error = "auto")` + exact residualized-seed quadratic-form variance + HKSJ random-effects pooling addressed the Step 4/6 variance and shared-error approximations (`PR #64`); and (4) `medusa_pipeline_primer.qmd` and `medusa_pipeline_primer_audited.qmd` have been updated to reflect these current semantics and to correct the legacy `voxelwise_deconvolution()` description.

## Material corrections

1. **Claims contradicted the computed results.** The original connectivity table reported naïve seed/PPI slopes −0.066/0.045 and M-MOMS −0.069/0.005 against claimed truths 0.55/0.35, yet asserted successful recovery. The cue-only plot also remained strongly nonzero despite a claim that correction flattened it. **Fixed:** removed automatic success claims, made comparison coordinates explicit, corrected the cue source duration mismatch, and described the realized residual artifacts. The revised connectivity simulation still fails end-to-end recovery (seed −0.247 naïve, −0.161 M-MOMS); this is visible in the report. A separate 200-replicate known-error experiment demonstrates the statistical correction under its own assumptions.

2. **The HRF and noise descriptions did not match implementation.** Target-run kernels are chosen from candidate fits scored within other runs, not by predicting the target with training-run amplitudes. The preset has 27 kernels, not 36. The original recommendations included unsupported `amplitude_model="auto"`, `parcel_method`, and `ar_order` arguments. GCV selects one scalar lambda per call. Noise scaling blends SD and zero-centered MAD, rather than using pure MAD. **Fixed:** corrected prose, tables, equations, and diagnostic labels, including the distinction between a kernel parameter and its sampled peak.

3. **Weights, whitening, and inference were overstated.** `tau2_mode="auto"` depends on supplied covariates; it is not marginal REML. The ZCA-style transform is damped and variance-preserving, not full covariance whitening. `restore_scale=TRUE` reverses noise-SD amplitude units. The old group example omitted random slopes and substituted asymptotic normal tests for calibrated group inference. **Fixed:** use explicit marginal weights, joint nuisance adjustment throughout, matched transformed outcome/predictors and random slopes, and descriptive pointwise intervals without significance symbols. Removed the inappropriate `plan_fwe_correction()` recommendation for a parcel/lag table. Documented the need for a prespecified test family and valid inference.

4. **Connectivity methods were compared under different defaults and estimands.** Naïve GLS defaults to transformed readouts, while M-MOMS defaults to original readouts. Residualizing a task-related seed changes the target equation and can leave an omitted quadratic task term. The between-person wrapper requires `repair` and genuine fitted random slopes. **Fixed:** compare both methods with `use_whitened=FALSE`, retain the full centered-seed generator with `orthogonalize_seed_x=FALSE`, label latent references as such, and explain the optional between-person workflow without claiming an executable demonstration that was absent. Corrected the description of subject pooling to the actual clipping, sample-size weighting, and empirical t inference.

5. **Graphics and code obscured the workflow.** The original three-column overview had clipped titles and legends, and mislabeled within-subject slopes as group MLM/connectivity results. **Fixed:** split it into three readable figures, normalize each library curve separately, shorten titles and legends, use redundant line styles, add alternative text, and arrange group parcels in two rows. Add a penalty-sensitivity figure with known truth. Supporting simulation, plot, and diagnostic code starts folded; five applied fitting blocks stay visible. Table headers are readable and long content can wrap or scroll.

## Remaining method limitations and proposed follow-up

- **P1: cascaded covariance propagation is approximate.** An independent stationary AR(3) reference using `B=I`, n=80 and cascade parameters (0.6, 0.12, −0.24) gives maximum diagonal discrepancy 0.148381 and maximum interior discrepancy 0.007940. This is a demonstrated discrepancy against the stationary unit-marginal reference, not a failing existing regression test. The report now states that propagation is conditional and uses boundary/normalization approximations. **Proposed:** define the intended initial-state/gap covariance convention, then validate or replace cascade normalization and boundary handling against that reference before claiming exact stationary covariance.
- **P1: connectivity correction uses diagonal trial-error approximations.** Seed residualization uses a leverage-based adjustment rather than `diag(R Omega R')` under general heteroscedastic/correlated errors; the moments do not carry the complete between-trial covariance. Separately transformed seed/target coordinates also complicate nonlinear PPI terms. **Proposed:** validate serial-dependence and shared-error behavior in larger simulations, carry the intended covariance through residualization, and keep the estimand fixed across comparisons.
- **P1: end-to-end recovery remains inadequate in this teaching simulation.** The phasic RPE maximum is reconstructed at −0.65 s instead of the known +0.65 s, and connectivity/PPI are not recovered. **Proposed:** run matched impulse/source-model, HRF-duration, overlap, SNR, and trial-count sensitivity studies before using this example as evidence of accurate temporal or connectivity recovery. The primer now treats these discrepancies as diagnostics.
- **P2: broader validation remains outside this review.** Eight subjects and twelve trials per run cannot establish null calibration, interval coverage, optimal settings, or individual-difference reliability. The optional between-person analysis is explained but not fitted. The full simulation validation suite and all package tests were not run.

## Verification and coverage

The inventory is **11 main sections, 7 plotted figures, 8 tables, and 13 principal method claims**. Counts below describe the scoped artifact inventory, not complete verification of every possible method setting. Zero observed defects does not mean every underlying statistical assumption has been validated.

### Report quality

| Category | Observed defects | Assessment |
|:---|:---|:---|
| Usefulness and completeness | 0 / 11 | All six workflow steps, the simulation, inputs, recommendations, and legacy bridge are covered. Between-person inference is explicitly an optional extension. |
| Analytical clarity | 0 / 11 | Fixed cross-fitting, scaling, uncertainty, modeling, and recovery descriptions. Remaining failures are disclosed. |
| Visual and interaction consistency | 0 / 15 | Seven plots and eight tables rendered and visually inspected on desktop; folding works. Narrow layout checked for page overflow, with wide numeric tables scrollable. |

### Analytical correctness and robustness

| Category | Observed defects | Assessment |
|:---|:---|:---|
| Source authority | 0 / 13 | Principal implementation claims traced to current source; historical sequence checked against the supplied talk. |
| Calculation/value accuracy | 0 / 15 | R-generated output and source bindings checked; consequential cross-fit/GLS/covariance claims independently tested. This is not independent recomputation of every table cell. |
| Within-chart agreement | 0 / 7 | Titles, captions, intervals, units, and actual plotted behavior agree; no blanket successful-recovery claims remain. |
| Complete source details | 0 / 11 | Source revision, reproducible simulation, assumptions, and methods recorded; installed dependency environment remains required for rerendering. |
| Cross-artifact consistency | 0 / 2 | QMD and HTML agree; historical talk references match the extracted slide text. |
| Data-quality controls | 0 / 3 | Simulation/source-event agreement, run/voxel alignment, and target covariate ordering addressed. Real-data censoring is documented, not demonstrated. |
| Conclusion support | 0 / 11 | Demonstrations are scoped to realized outputs; method limitations remain explicit. |

Ten selected existing tests passed **183 assertions**, covering target exclusion, gap-aware whitening/independent GLS, joint masked inverse, duration effects, precision/ZCA propagation, cue/lag correction, single-run blocks, parcel aggregation, and fast-TR changes. `medusa_pipeline_primer_method_checks.R` and its log retain the exact harness and independent stationary covariance/ZCA diagnostics. These are targeted tests, not a full `R CMD check`.

The HTML is self-contained for figures and styling. Rendering uses the installed package by default, or a checkout specified with `MEDUSA_SOURCE`. Rendering was performed with R 4.5.2, the source checkout, existing compiled code, and local dependency libraries. The rendering environment's temporary ggbrain/imager compatibility libraries were not used for the numerical examples. Quarto and headless Chrome required their normal local cache/process permissions; no publishing or messages were sent.
