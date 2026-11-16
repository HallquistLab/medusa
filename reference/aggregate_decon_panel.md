# Aggregate a Voxel-Level MEDuSA Trial Panel to Parcels / ROIs

Collapses a voxel-level
[`decon_epoch_regression`](https://hallquistlab.github.io/medusa/reference/decon_epoch_regression.md)
output (or stacked trial-by-voxel `data.table`) into a parcel-level
panel while propagating shared cross-voxel and cross-parcel
measurement-error vectors (`shared_noise_vec`), cross-trial correlations
(`trial_cor_vec`), cross-lag embeddings (`lag_noise_vec`), exact parcel
`se_demixed`, Empirical-Bayes `tau2_bio` / `precision_weight`, and any
attached `shared_error_model` for downstream
[`decon_connectivity_repair`](https://hallquistlab.github.io/medusa/reference/decon_connectivity_repair.md)
or
[`mixed_by`](https://hallquistlab.github.io/medusa/reference/mixed_by.md)
analyses.

## Usage

``` r
aggregate_decon_panel(
  x,
  parcel_groups,
  parcel_weights = c("auto", "equal", "supplied"),
  weights = NULL,
  voxel_col = "voxel",
  parcel_col = "parcel",
  id_col = "id",
  run_col = "run",
  trial_col = "trial",
  lag_col = "evt_time"
)
```

## Arguments

- x:

  A voxel-level
  [`decon_epoch_regression`](https://hallquistlab.github.io/medusa/reference/decon_epoch_regression.md)
  result list, a trial-level `data.frame`/`data.table` (including
  CSV/Parquet panels read via
  [`data.table::fread`](https://rdrr.io/pkg/data.table/man/fread.html)
  or `medusa_read_table`), a character file path (or vector of paths),
  or a list of per-subject results.

- parcel_groups:

  Either a named character/factor vector mapping voxel IDs to parcel
  labels, an unnamed vector matching `voxel_ids` order, or a
  `data.frame` containing voxel and parcel mapping columns.

- parcel_weights:

  Weighting rule within each parcel: `"auto"` (inverse squared
  `se_demixed` with 10% median floor), `"equal"`, or `"supplied"`.

- weights:

  Optional named numeric vector of voxel weights when
  `parcel_weights = "supplied"`.

- voxel_col:

  Name of the voxel column in the input panel (default `"voxel"`).

- parcel_col:

  Name of the output parcel column (default `"parcel"`).

- id_col, run_col, trial_col, lag_col:

  Subject, run, trial, and lag column names.

## Value

If `x` is a `decon_epoch_regression` result list, an updated
parcel-level result list with `panel`, `parcel_weights`,
`unit_error_cor`, and `shared_error_model`. If `x` is a
`data.frame`/`data.table`, a parcel-aggregated `data.table`.
