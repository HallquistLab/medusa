# Controlled trial-panel simulations for Step 5 using live library(medusa).
library(medusa)
library(data.table)
set.seed(5112026)

# Within- and between-participant effects can have opposite signs.
center <- function(x) x - mean(x)
m <- 18L
ids <- factor(sprintf("P%02d", seq_len(24)))
wb <- CJ(id = ids, trial = seq_len(m))
wb[, x_mean := rep(seq(-1.6, 1.6, length.out = length(ids)), each = m)]
wb[, x_within := center(rnorm(.N, sd = .52)), by = id]
wb[, x := x_mean + x_within]
wb[, activity := -.72 * x_mean + rep(rnorm(length(ids), sd = .13), each = m) +
     .55 * x_within + rnorm(.N, sd = .27)]
pooled <- lm(activity ~ x, wb)
separated <- lme4::lmer(activity ~ x_within + x_mean + (1 | id), wb)
within_slopes <- wb[, .(intercept = coef(lm(activity ~ x))[1L],
                        slope = coef(lm(activity ~ x))[2L],
                        x_min = min(x), x_max = max(x)), by = id]
wb_coefs <- c(pooled = unname(coef(pooled)["x"]), lme4::fixef(separated))

# A balanced population-null design. Conditional normal sampling of the slope
# estimates and residual sums of squares gives the same OLS and participant tests
# as explicitly generating every row with x = -1,+1.
set.seed(5112027)
ns <- 12L
nrep <- 3000L
sd_slope <- .45
sd_error <- .70
null <- rbindlist(lapply(c(12L, 48L, 192L), function(n_trial) {
  slopes <- matrix(rnorm(nrep * ns, sd = sqrt(sd_slope^2 + sd_error^2 / n_trial)), nrep)
  estimate <- rowMeans(slopes)
  ss_slopes <- rowSums((slopes - estimate)^2)
  rss_local <- sd_error^2 * rchisq(nrep, df = ns * (n_trial - 1L) - 1L)
  se_people <- sqrt(ss_slopes / ((ns - 1L) * ns))
  se_pooled <- sqrt((rss_local + n_trial * ss_slopes) / (ns * n_trial - 2L) / (ns * n_trial))
  data.table(n_trial, method = rep(c("Trial-row OLS", "Participant slopes"), each = nrep),
             estimate = rep(estimate, 2), se = c(se_pooled, se_people),
             df = rep(c(ns * n_trial - 2L, ns - 1L), each = nrep))
}))
null[, reject := abs(estimate / se) > qt(.975, df)]
null_summary <- null[, .(rejection = mean(reject), mean_se = mean(se),
                         empirical_sd = sd(estimate), replications = .N), by = .(n_trial, method)]
null_summary[, mcse := sqrt(rejection * (1 - rejection) / replications)]

# A synthetic parcel/lag panel with two behavioral terms and one null parcel.
# Trial residuals are independent within each readout after the random effects.
# Repeated readouts at neighboring lags have correlated errors.
set.seed(5112028)
lags <- -2:5
parcels <- c("Phasic", "Sustained", "Null")
base <- CJ(id = ids, run = factor(c("run1", "run2")), trial = seq_len(24))
base[, event_index := 2L * trial]
base[, onset := 10 + 8 * trial]
base[, value_mean := rep(rnorm(length(ids), sd = .6), each = 48L)]
base[, value_within := center(rnorm(.N)), by = id]
base[, entropy_within := center(rnorm(.N)), by = id]
latent_sd <- .16
base[, se_demixed := ifelse(run == "run1", .11, .37)]
base[, precision_weight_cohort := 1 / (latent_sd^2 + se_demixed^2)]
base[, precision_weight := {
  w <- 1 / (latent_sd^2 + se_demixed^2)
  w / mean(w)
}, by = id]
lag_factor <- chol(toeplitz(.68^(0:(length(lags) - 1L))))
draw_by_lag <- function(n, sd) {
  matrix(rnorm(n * length(lags), sd = sd), n) %*% lag_factor
}
known <- panel <- list()
for (parcel in parcels) {
  beta_value <- if (parcel == "Phasic") .27 * exp(-.5 * ((lags - 1) / .85)^2) else rep(0, length(lags))
  beta_entropy <- if (parcel == "Sustained") .23 * plogis((lags - .65) / .35) *
    plogis((4.6 - lags) / .55) else rep(0, length(lags))
  known[[parcel]] <- data.table(parcel, evt_time = rep(lags, 2),
    term = rep(c("value_within", "entropy_within"), each = length(lags)),
    truth = c(beta_value, beta_entropy))
  alpha_id <- draw_by_lag(length(ids), .27)
  slope_value <- draw_by_lag(length(ids), .09)
  slope_entropy <- draw_by_lag(length(ids), .085)
  alpha_run <- draw_by_lag(2L * length(ids), .15)
  E <- draw_by_lag(nrow(base), 1) * sqrt(latent_sd^2 + base$se_demixed^2)
  panel[[parcel]] <- rbindlist(lapply(seq_along(lags), function(k) {
    z <- copy(base)
    s <- as.integer(z$id)
    j <- 2L * (s - 1L) + as.integer(z$run)
    z[, `:=`(parcel = parcel, evt_time = lags[k],
      decon_demixed = .35 - .2 * value_mean + alpha_id[s, k] + alpha_run[j, k] +
        (beta_value[k] + slope_value[s, k]) * value_within +
        (beta_entropy[k] + slope_entropy[s, k]) * entropy_within + E[, k])]
    z
  }))
}
panel <- rbindlist(panel)
known <- rbindlist(known)
stopifnot(!anyDuplicated(panel[, .(id, run, event_index, evt_time, parcel)]),
          all(panel$precision_weight > 0))
model <- decon_demixed ~ value_within + value_mean + entropy_within +
  (1 + value_within + entropy_within || id) + (1 | id:run)
batch <- mixed_by(
  panel, model_formulae = list(trial_effects = model),
  split_on = c("parcel", "evt_time"), weights = "precision_weight",
  calculate = "parameter_estimates_reml", padjust_by = NULL,
  tidy_args = list(effects = "fixed", conf.int = FALSE),
  return_models = TRUE
)
fixed <- merge(batch$coef_df_reml, known, by = c("parcel", "evt_time", "term"))
stopifnot(nrow(fixed) == nrow(known), all(is.finite(fixed$p.value)),
          max(abs(fixed$p.value - 2 * pt(-abs(fixed$statistic), fixed$df))) < 1e-10)
fixed[, `:=`(lo = estimate - qt(.975, df) * std.error,
             hi = estimate + qt(.975, df) * std.error, p_by = p.adjust(p.value, method = "BY"))]
models <- lapply(batch$models$model, function(entry) entry$reml)
cat("Returned model entries:", paste(names(batch$models$model[[1L]]), collapse = ", "), "\n")
stopifnot(all(vapply(models, inherits, logical(1), "merMod")))
diagnostics <- batch$models[, .(parcel, evt_time)]
diagnostics[, `:=`(
  singular = vapply(models, lme4::isSingular, logical(1)),
  n_obs = vapply(models, function(x) as.integer(nobs(x)), integer(1)),
  n_subjects = vapply(models, function(x) nlevels(model.frame(x)$id), integer(1)),
  n_issues = vapply(models, function(x) length(x@optinfo$conv$lme4$messages), integer(1)))]
focal <- panel[parcel == "Phasic" & evt_time == 1]
direct <- lmerTest::lmer(model, data = focal, weights = precision_weight, REML = TRUE)
check_fit <- models[[which(diagnostics$parcel == "Phasic" & diagnostics$evt_time == 1)]]
stopifnot(max(abs(lme4::fixef(direct) - lme4::fixef(check_fit))) < 1e-8,
          max(abs(weights(check_fit) - focal$precision_weight)) < 1e-10)
unweighted <- lmerTest::lmer(model, data = focal, REML = TRUE)
weight_fit <- rbindlist(lapply(list(Weighted = direct, Unweighted = unweighted), function(fit) {
  data.table(run = focal$run, sd_model = sigma(fit) / sqrt(weights(fit)),
             sd_true = sqrt(latent_sd^2 + focal$se_demixed^2))
}), idcol = "model")
weight_fit <- unique(weight_fit)
sensitivity <- rbindlist(lapply(list(Weighted = direct, Unweighted = unweighted), function(fit) {
  x <- coef(summary(fit))["value_within", ]
  data.table(estimate = x[1L], se = x[2L], df = x[3L], p = x[5L])
}), idcol = "model")
focal_residuals <- focal[, .(id, run, trial)]
focal_residuals[, scaled_residual := sqrt(weights(direct)) * residuals(direct) / sigma(direct)]
setorder(focal_residuals, id, run, trial)
residual_acf <- focal_residuals[, .(ac1 = cor(head(scaled_residual, -1L),
                                            tail(scaled_residual, -1L))), by = .(id, run)]
pointwise <- fixed[, sum(p.value <= .05)]
family <- fixed[, sum(p_by <= .05)]
cat("Within/between comparison\n")
print(wb_coefs)
cat("Null-simulation rejection and Monte Carlo SE\n")
print(null_summary)
cat("Parcel/lag fits:", nrow(diagnostics), "Singular:", sum(diagnostics$singular),
    "Fixed tests:", nrow(fixed), "Pointwise:", pointwise, "BY:", family, "\n")
print(weight_fit)
print(sensitivity)
cat("Matched direct and mixed_by fixed effects and observation weights.\n")

fwrite(panel, "synthetic_panel.csv")
fwrite(fixed, "group_coefficients.csv")
fwrite(null_summary, "null_calibration.csv")
fwrite(diagnostics, "model_diagnostics.csv")
saveRDS(list(wb = wb, wb_coefs = wb_coefs, within_slopes = within_slopes,
             fixed = fixed, null = null_summary, diagnostics = diagnostics,
             weight_fit = weight_fit, sensitivity = sensitivity,
             focal_residuals = focal_residuals, residual_acf = residual_acf),
        "step5_examples.rds")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
