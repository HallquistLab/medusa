# Between-subject REPAIR structural regression on MEDuSA mixed models

Connects fitted [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)
(or
[`mixed_by`](https://hallquistlab.github.io/medusa/reference/mixed_by.md)
model objects) from MEDuSA Stage 4 directly to
[`repair::repair()`](https://hallquistlab.github.io/repair/reference/repair.html)
(Vig & Hallquist, 2026) when the repair package is installed, with an
exact built-in likelihood-only WLS +
[`medusa_fit_mmom`](https://hallquistlab.github.io/medusa/reference/medusa_fit_mmom.md)
fallback when repair is unavailable. Estimates
measurement-error-corrected between-subject regressions relating
subject-specific neural random intercepts/slopes across brain regions or
to external behavioral/clinical variables.

## Usage

``` r
medusa_between_subject_repair(
  formula,
  models = NULL,
  data = NULL,
  id = "id",
  regularization = c("step_up", "minimum", "none", "custom"),
  cross_model_error = NULL,
  medusa_mlm = NULL,
  subject_data = NULL,
  id_col = NULL,
  ...
)
```

## Arguments

- formula:

  Structural regression formula using `re()` markers (e.g.,
  `re("x_whitened", model = "target") ~ re("x_whitened", model = "seed") + age`).

- models:

  Named list of fitted
  [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) or
  [`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html) models (or a
  single fitted model).

- data:

  Subject-level `data.frame` with one row per subject ID.

- id:

  Character scalar naming the subject ID column in `data` (default
  `"id"`).

- regularization:

  Regularization policy passed to
  [`repair::repair()`](https://hallquistlab.github.io/repair/reference/repair.html)
  or
  [`medusa_fit_mmom`](https://hallquistlab.github.io/medusa/reference/medusa_fit_mmom.md):
  `"step_up"` (default; adaptive M-MOMS), `"minimum"`, `"none"`, or
  `"custom"`.

- cross_model_error:

  Cross-model prediction-error covariance policy when `models` has
  multiple entries (defaults to `"independent"` when
  `length(models) > 1L`).

- medusa_mlm:

  Optional alias for `models` (e.g., a single fitted
  [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) model).

- subject_data:

  Optional alias for `data`.

- id_col:

  Optional alias for `id`.

- ...:

  Additional arguments passed to
  [`repair::repair()`](https://hallquistlab.github.io/repair/reference/repair.html)
  (or `alpha` for
  [`medusa_fit_mmom`](https://hallquistlab.github.io/medusa/reference/medusa_fit_mmom.md)).

## Value

A list with `coefficients` (tidy `data.frame` with `term`, `estimate`,
`std.error`, `statistic`, `p.value`, `conf.low`, `conf.high`),
`reliability`, `predictions` (likelihood-only random-effect predictions
and error variances), and `fit`.
