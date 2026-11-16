# Possible Future Feature: Multi-Target & Connectome Batch Mode in `decon_connectivity_repair()`

> **Status:** Proposed Future Feature (Design Specification — Not Yet Implemented)  
> **Target Package / Module:** `medusa` (`R/diagnose_crossrun_recovery.R`, `R/poet_covreg.R`, `src/poet_covreg.cpp`)  
> **Date:** November 2026

---

## 1. Motivation & Bottleneck Analysis

`decon_connectivity_repair()` currently estimates measurement-error-corrected inter-region connectivity, directed lagged connectivity (`control_target_seed_lag = TRUE`), and task-modulated PPI (`test_ppi = TRUE`) for **one `(seed_region, target_region, seed_lag, target_lag)` pair per call**.

When users scale from hypothesis-driven ROI pairs to **seed-to-parcellation maps** ($1 \times 400$ parcels) or **whole-brain directed parcellated connectomes** ($200 \times 199 = 39{,}800$ directed edges; note REPAIR `M-MOMS` is directionally asymmetric because predictor and outcome error covariances enter differently), calling `decon_connectivity_repair()` inside an external R `lapply()` loop repeats three expensive operations on every edge $e = 1, \dots, E$:

1. **Repeated Panel Ingestion, Scanning & Merging ($O(E \cdot N_{\text{rows}})$):**
   - Each call copies `panel_dt`, checks/restores serialized list-columns (`shared_noise_vec`, `trial_cor_vec`, `lag_noise_vec`), performs two or three full-table scans (`dt[region == ... & lag == ...]`), and runs `merge()` on `(id, run, trial)`.
2. **Redundant Per-Node & Per-Run Projections ($O(P^2)$ instead of $O(P)$):**
   - Run-level covariate QR orthogonalization ($\mathbf{M}_0, \mathbf{M}_s, \mathbf{P}_0 = \mathbf{I} - \mathbf{M}_0 \mathbf{M}_0^\top, \mathbf{P}_s = \mathbf{I} - \mathbf{M}_s \mathbf{M}_s^\top$) and trial-onset distance kernels $\exp(-|\Delta t|/6)$ depend **only on `(subject, run)`**.
   - Node-level quantities — `seed_c`, `targ_c`, operator projections ($\mathbf{U}_S \mathbf{P}_s$, $\mathbf{V}_S \mathbf{P}_s$, $\mathbf{U}_T \mathbf{P}_0$, $\mathbf{V}_T \mathbf{P}_0$), marginal error variances ($\psi_S$, $\psi_Y$, $\psi_{T_0}$), and within-run trial correlation blocks ($\mathbf{P}_{SS}$, $\mathbf{P}_{YY}$, `build_stacked_trial_cor()`) — depend **only on `(subject, run, region, lag)`**, not on the partner region.
   - Across a $200 \times 200$ connectome, per-edge evaluation performs $\sim 200\times$ redundant operator and matrix setup.
3. **Inner R Loops Over Trials, Subjects & Edges:**
   - `.project_trial_err_psd()` currently calls R's `eigen()` on $2\times 2$ / $3\times 3$ matrices per trial, and both `medusa_fit_mmom()` and the `5A` cross-trial hybrid operator + onset-kernel sandwich SEs run in R per `(subject, edge)`.

---

## 2. Proposed Backward-Compatible API

Scalar `(seed_region, target_region, seed_lag, target_lag)` inputs continue to work identically to the current interface, while vector, `"all"`, or explicit `edges` inputs activate batch mode:

```r
decon_connectivity_repair(
  panel_dt,
  seed_region = NULL,          # scalar, character vector, or "all"
  target_region = NULL,        # scalar, character vector, or "all" (NULL => "all" when seed supplied)
  seed_lag = 0,                # scalar or numeric vector
  target_lag = 0,              # scalar or numeric vector
  edges = NULL,                # optional data.frame(seed_region, target_region[, seed_lag, target_lag])
  include_self = c("lagged", "FALSE", "TRUE"), # omit S == T when seed_lag == target_lag by default
  p_adjust = c("BH", "BY", "holm", "bonferroni", "none"),
  return_subjects = c("auto", "TRUE", "FALSE"), # TRUE when n_edges <= 2000; FALSE above unless requested
  return_matrices = FALSE,     # optional P x P array/list of group estimates, t-stats, & q-values
  n_threads = 1L,              # OpenMP threads for batched C++ edge solver
  ...
)
```

### Supported Topologies
- **Seed-to-Parcel Map ($1 \to P$):** `seed_region = "AMY", target_region = "all"`
- **Target Inflow Map ($P \to 1$):** `seed_region = "all", target_region = "vmPFC"`
- **Full Directed Connectome ($P \times (P - 1)$):** `seed_region = "all", target_region = "all"`
- **Lag-Grid / Directed Granger Panel:** vector `seed_lag` / `target_lag` with `control_target_seed_lag = TRUE`
- **Custom Hypothesis Subgraph:** `edges = data.frame(seed_region = ..., target_region = ...)`

---

## 3. Proposed Two-Stage Architecture

### Stage A — Single-Pass Panel Indexing & Unary Node Precomputation ($O(P)$)
1. **Canonical `(id, run, trial)` Grid Alignment:**
   - Unpack `panel_dt` once, restore list-columns once, and index all requested `(region, lag)` slices onto a shared per-`(id, run)` trial index so zero per-edge `data.table` filtering or `merge()` calls occur.
2. **Per-`(subject, run)` Design Cache:**
   - Precompute $\mathbf{M}_0$, $\mathbf{M}_s$, $\mathbf{P}_0$, $\mathbf{P}_s$, centered covariates $\mathbf{X}_{c,r}$, and inter-trial onset differences $|\Delta t_{ab}|$ once per `(subject, run)`.
3. **Per-`(subject, run, region, lag)` Unary Node Cache:**
   - Precompute **as-seed** quantities (`seed_c`, $\psi_{S,\text{raw}}$, $\psi_{S,\text{floor}}$, $\mathbf{U}_S \mathbf{P}_s$, $\mathbf{V}_S \mathbf{P}_s$, $\mathbf{P}_{SS}$) and **as-target** quantities (`targ_c`, $\psi_{Y,\text{raw}}$, $\psi_{Y,\text{floor}}$, $\mathbf{U}_T \mathbf{P}_0$, $\mathbf{V}_T \mathbf{P}_0$, $\mathbf{P}_{YY}$, plus baseline lag-$t_0$ blocks when `control_target_seed_lag = TRUE`).
   - **Operator-mode stratum deduplication:** Cache $\mathbf{U}_{\text{strata}}[[s]][[\ell]] \mathbf{P}_s$ and $\mathbf{U}_{\text{strata}}[[s]][[\ell]] \mathbf{P}_0$ by `(stratum_id, lag)` so parcels sharing an HRF/AR(1) stratum reuse the exact same projected operator matrices.
   - **Stacked-mode dense embedding packing:** Pack `shared_noise_vec` and `lag_noise_vec` once per node into contiguous $m_r \times K_{\text{emb}}$ matrices so cross-parcel trial error covariances reduce to fast matrix operations rather than per-trial R `unlist()` loops.

### Stage B — Batched C++/Armadillo Edge Solver (`src/poet_covreg.cpp`)
Implement a compiled C++ batch kernel (with optional OpenMP `n_threads` parallelism over edges):
1. **Closed-form $2\times 2$ and fast $3\times 3$ trial error PSD projection** (eliminating per-trial R `eigen()` allocations).
2. **Compiled `M-MOMS` (`step_up` / `minimum`) solver** matching `medusa_fit_mmom()` to machine precision.
3. **Compiled `5A` Cross-Trial Sandwich SE & Satterthwaite `df`**, assembling only the cross-block $\mathbf{P}_{SY}$ per edge and combining it with the cached unary blocks $\mathbf{P}_{SS}$ and $\mathbf{P}_{YY}$.
4. **Vectorized Group Random-Effects + HKSJ Pooling** across subjects, returning `group` with `p.value_adj` alongside `p.value`.

---

## 4. Numerical Parity & Memory Guardrails

- **Exact Single-Edge Parity:** Every edge in batch mode is guaranteed to match calling single-edge `decon_connectivity_repair()` to machine precision (`~1e-12`).
- **Memory Guardrails:**
  - `return_subjects = "auto"` retains the per-subject `subjects` table when `n_edges <= 2000` and omits it by default on large connectomes (e.g., $40{,}000$ edges $\times$ $50$ subjects) unless `return_subjects = TRUE`.
  - `return_augmented_panel = TRUE` is disallowed (or restricted with an explicit warning/limit) when `n_edges > 1` to avoid $O(E \cdot N_{\text{trials}})$ memory blow-up.

---

## 5. Open Design Decisions to Resolve Before Implementation

1. **Multiple-Testing Family (`p_adjust`):**
   - Default to adjusting `p.value_adj` across all tested edges **within each `(term, seed_lag, target_lag)` slice** (recommended for connectome/seed maps so main connectivity `seed_c` and PPI `seed_c:x_c` are controlled within their respective edge families), with an option to adjust jointly across lags when estimating a multi-lag surface.
2. **Matrix Extractor / S3 Interface:**
   - Decide whether to include both `return_matrices = TRUE` and a lightweight S3 helper `connectome_matrix(x, term = "seed_c", value = c("estimate", "statistic", "p.value", "p.value_adj"))` returning a `[seed_region, target_region]` matrix for plotting and network/NBS downstream tools.
3. **Optional Empirical-Bayes Variance Moderation Across Edges:**
   - Decide whether group-level inference should remain strictly edge-wise HKSJ by default (preserving exact single-edge invariance) with an optional `moderate_eb = TRUE` mode for small-$N$ connectome studies ($N < 20$).
