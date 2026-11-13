## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 8,
                      fig.height = 4.5, fig.align = "center")
library(medusa)
library(ggplot2)
needed <- c("estimate_crossrun_hrf", "deconvolve_crossrun_hrf", "decon_epoch_regression",
            "mixed_by", "medusa_multilag_infer", "decon_connectivity_repair")
stopifnot(all(vapply(needed, exists, logical(1),
                    envir = asNamespace("medusa"), inherits = FALSE)))


## ----simulate-----------------------------------------------------------------
set.seed(5401)
tr <- 1
n_time <- 220L
n_trials <- 20L
voxel_ids <- paste0("v", 1:6)
run_ids <- paste0("run", 1:3)
true_beta <- seq(0.5, 1.0, length.out = length(voxel_ids))
candidates <- hrf_library(tr, peak_delays = c(5, 6, 7),
                          dispersions = 1, undershoot_ratios = 6)
# Canonical (6, 1, 6) is placed first when present in hrf_library().
true_kernel_index <- c(1L, 1L, 1L, 2L, 2L, 2L)

source_column <- function(onset, duration, n, tr) {
  left <- (0:(n - 1L)) * tr
  if (duration > 0) {
    return(pmax(0, pmin(left + tr, onset + duration) - pmax(left, onset)) / tr)
  }
  pmax(0, 1 - abs(left / tr - onset / tr))
}
forward_bold <- function(source, hrf) {
  as.numeric(stats::convolve(source, rev(hrf), type = "open"))[seq_along(source)]
}
ar_noise <- function(n, rho, sigma) {
  as.numeric(stats::arima.sim(list(ar = rho), n = n,
                             sd = sigma * sqrt(1 - rho^2)))
}

simulate_participant <- function(beta = true_beta) {
  bold <- events <- covariates <- decision_covariates <- rt_by_run <- list()
  for (r in run_ids) {
    cue <- seq(6, by = 9.5, length.out = n_trials)
    decision <- cue + runif(n_trials, 1.0, 1.6)
    value <- as.numeric(scale(rnorm(n_trials)))
    rt <- pmin(2.5, pmax(0.5, 1.4 - 0.25 * value + rnorm(n_trials, sd = 0.2)))
    feedback <- decision + rt + runif(n_trials, 0.8, 1.5)
    rpe <- as.numeric(scale(rnorm(n_trials)))
    # Keep the event table sorted; event_index will refer to these exact rows.
    ev <- rbind(
      data.frame(onset = cue, duration = 0, condition = "cue", trial_id = 1:n_trials),
      data.frame(onset = decision, duration = rt, condition = "decision", trial_id = 1:n_trials),
      data.frame(onset = feedback, duration = 0, condition = "feedback", trial_id = 1:n_trials)
    )
    ev <- ev[order(ev$onset), ]
    rownames(ev) <- NULL
    E <- vapply(seq_len(nrow(ev)), function(i) {
      source_column(ev$onset[i], ev$duration[i], n_time, tr)
    }, numeric(n_time))
    Y <- vapply(seq_along(voxel_ids), function(v) {
      amplitude <- ifelse(ev$condition == "cue", 0.8,
                          ifelse(ev$condition == "decision",
                                 1.0 + 0.35 * value[ev$trial_id],
                                 0.6 + beta[v] * rpe[ev$trial_id]))
      source <- as.numeric(E %*% amplitude)
      100 + forward_bold(source, candidates[, true_kernel_index[v]]) + ar_noise(n_time, 0.35, 0.04)
    }, numeric(n_time))
    colnames(Y) <- voxel_ids
    bold[[r]] <- Y
    events[[r]] <- ev
    # Rows must match the selected event records, not all events in the run.
    fi <- ev$trial_id[ev$condition == "feedback"]
    di <- ev$trial_id[ev$condition == "decision"]
    covariates[[r]] <- data.frame(rpe = rpe[fi])
    decision_covariates[[r]] <- data.frame(value = value[di])
    rt_by_run[[r]] <- rt[di]
  }
  list(bold = bold, events = events, covariates = covariates,
       decision_covariates = decision_covariates, rt_by_run = rt_by_run)
}
sim <- simulate_participant()
bold <- sim$bold
events <- sim$events
covariates <- sim$covariates
decision_covariates <- sim$decision_covariates
rt_by_run <- sim$rt_by_run
stopifnot(identical(names(bold), names(events)),
          identical(names(bold), names(covariates)),
          all(vapply(bold, function(y) identical(colnames(y), voxel_ids), logical(1))))
head(events$run1, 6)


## ----observed-plot, fig.cap="Observed BOLD from one synthetic voxel. Event overlap makes the individual cue, decision, and feedback contributions difficult to read directly."----
ggplot(data.frame(time = (0:(n_time - 1L)) * tr, bold = bold$run1[, 1L]),
       aes(time, bold)) +
  geom_line() +
  geom_vline(xintercept = events$run1$onset[events$run1$condition == "feedback"],
             color = "#0072B2", alpha = 0.25) +
  labs(x = "Run time (s)", y = "BOLD input units") + theme_minimal(base_size = 12)


## ----estimate-hrf, class.source="fold-show"-----------------------------------
hrf_fit <- estimate_crossrun_hrf(
  bold = bold, events = events, tr = tr, kernel_grid = candidates,
  amplitude_model = "trial", noise_model = "ar1", evidence = "calibrated",
  pooling = "voxel", spatial_strength = 0, return_scores = TRUE
)
hrf_fit$training_runs$run1
head(hrf_fit$selections$run1)


## ----parcel-recipe, eval=FALSE------------------------------------------------
# medusa_candidates <- hrf_library(tr, preset = "medusa")
# parcel_fit <- estimate_crossrun_hrf(
#   bold, events, tr, medusa_candidates,
#   spatial_groups = parcel_ids, pooling = "parcel",
#   amplitude_model = "trial", noise_model = "ar1", evidence = "calibrated",
#   selection_mode = "soft_average"
# )


## ----deconvolve, class.source="fold-show"-------------------------------------
decon <- deconvolve_crossrun_hrf(
  bold, hrf_fit,
  lambda = 0.1, scaling = "noise_sd", penalty = "sobolev", backend = "cpp",
  nuisance_mode = "ignore", observation_mask = "fit", tail_padding = 0L
)
head(decon$scaling$run1)
dim(decon$activity$run1)


## ----operator-diagnostics-----------------------------------------------------
recovery <- diagnose_crossrun_recovery(decon, events = events, lags = 0)
knitr::kable(recovery$summary[, c("run", "fitted_kernel", "n_voxels",
                                  "mean_gain", "mean_abs_excess_leakage")], digits = 3)


## ----feedback-epochs, class.source="fold-show"--------------------------------
feedback_fit <- decon_epoch_regression(
  deconvolution = decon, covariates = covariates, events = events,
  target_condition = "feedback", lags = seq(-2, 4, by = 1),
  source_model = "impulse", readout_kernel = "triangle", readout_width = 1.5 * tr,
  demix_other_events = TRUE, demix_target_trials = TRUE,
  demix_within_trial = TRUE, weighting = "precision", whiten_panel = TRUE,
  restore_scale = TRUE, return_panel = TRUE
)
knitr::kable(feedback_fit$summary, digits = 3)


## ----recovery-comparison, fig.cap="Raw and corrected feedback RPE slopes. The dashed line marks the simulated onset effect averaged across voxels. Nonzero lag analyses shift the assumed source as well as the readout; they are not independent temporal measurements."----
s <- feedback_fit$summary
plot_df <- rbind(
  data.frame(lag = s$lag, beta = s$mean_beta_raw, estimate = "Raw readout"),
  data.frame(lag = s$lag, beta = s$mean_beta_corrected, estimate = "Corrected")
)
ggplot(plot_df, aes(lag, beta, color = estimate)) +
  geom_hline(yintercept = mean(true_beta), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted") + geom_line() + geom_point() +
  labs(x = "Lag relative to feedback (s)", y = "RPE slope (input units)", color = NULL) +
  theme_minimal(base_size = 12)

at_zero <- subset(feedback_fit$coefficients, lag == 0)
knitr::kable(data.frame(
  voxel = at_zero$voxel,
  truth = true_beta[match(at_zero$voxel, voxel_ids)],
  raw = at_zero$beta_raw,
  corrected = at_zero$beta_corrected,
  gls = at_zero$beta_gls,
  se_weighted = at_zero$se_weighted,
  se_gls = at_zero$se_gls
), digits = 3)


## ----canonical-comparison-----------------------------------------------------
canonical_fit <- estimate_crossrun_hrf(
  bold, events, tr, candidates[, 1L, drop = FALSE],
  amplitude_model = "trial", noise_model = "ar1", evidence = "calibrated"
)
canonical_decon <- deconvolve_crossrun_hrf(
  bold, canonical_fit, lambda = 0.1, scaling = "noise_sd", penalty = "sobolev",
  backend = "cpp", nuisance_mode = "ignore", observation_mask = "fit", tail_padding = 0L
)
canonical_feedback <- decon_epoch_regression(
  canonical_decon, covariates, events = events, target_condition = "feedback",
  lags = 0, source_model = "impulse", readout_kernel = "triangle",
  readout_width = 1.5 * tr, restore_scale = TRUE
)
comparison <- data.frame(
  voxel = voxel_ids,
  true_hrf_column = true_kernel_index,
  true_beta = true_beta,
  fixed_hrf_beta = canonical_feedback$coefficients$beta_corrected[
    match(voxel_ids, canonical_feedback$coefficients$voxel)],
  learned_hrf_beta = at_zero$beta_corrected[match(voxel_ids, at_zero$voxel)]
)
knitr::kable(comparison, digits = 3)


## ----duration-epochs----------------------------------------------------------
decision_fit <- decon_epoch_regression(
  deconvolution = decon, covariates = decision_covariates, events = events,
  target_condition = "decision", lags = c(0, 1),
  source_model = "duration", trial_durations = rt_by_run, duration_scale = "rate",
  readout_kernel = "triangle", readout_width = 1.5 * tr,
  restore_scale = TRUE, return_panel = TRUE
)
knitr::kable(decision_fit$summary, digits = 3)


## ----panel-identity-----------------------------------------------------------
panel <- as.data.frame(feedback_fit$panel)
panel$id <- "sim01"
# event_index indexes the original events table; trial is a target-row counter.
panel$trial_id <- vapply(seq_len(nrow(panel)), function(j) {
  events[[panel$run[j]]]$trial_id[panel$event_index[j]]
}, integer(1))
head(panel[, c("id", "run", "trial", "event_index", "trial_id", "evt_time",
               "voxel", "decon_raw", "decon_demixed", "se_demixed",
               "precision_weight", "precision_weight_cohort", "error_ac1",
               "decon_demixed_whitened", "rpe", "rpe_whitened")])
stopifnot(!anyDuplicated(panel[c("id", "run", "trial_id", "evt_time", "voxel")]),
          all(is.finite(panel$decon_demixed)),
          all(is.finite(panel$precision_weight) & panel$precision_weight > 0),
          all(is.finite(panel$precision_weight_cohort) & panel$precision_weight_cohort > 0))


## ----cohort-------------------------------------------------------------------
set.seed(5412)
parcel_ids <- setNames(rep(c("Parcel_A", "Parcel_B"), each = 3L), voxel_ids)
lags_group <- -2:4
participant_ids <- sprintf("sim%02d", seq_len(12L))
fit_one <- function(id) {
  subject <- simulate_participant(true_beta + rnorm(1L, sd = .12))
  fit <- estimate_crossrun_hrf(
    subject$bold, subject$events, tr, candidates, amplitude_model = "trial",
    noise_model = "ar1", evidence = "calibrated"
  )
  dec <- deconvolve_crossrun_hrf(
    subject$bold, fit, lambda = .1, scaling = "noise_sd", penalty = "sobolev",
    nuisance_mode = "ignore", observation_mask = "fit"
  )
  ep <- decon_epoch_regression(
    dec, subject$covariates, target_condition = "feedback", lags = lags_group,
    parcel_groups = parcel_ids, restore_scale = TRUE,
    demix_within_trial = FALSE, tau2_mode = "marginal", return_panel = TRUE
  )
  panel_i <- data.table::copy(ep$panel)
  panel_i[, id := id]
  for (run_i in names(subject$events)) {
    observed <- panel_i[run == run_i]
    records <- subject$events[[run_i]][observed$event_index, ]
    stopifnot(all(records$condition == "feedback"),
              isTRUE(all.equal(observed$onset, records$onset)))
  }
  stopifnot(!anyDuplicated(panel_i[, .(id, run, event_index, evt_time, parcel)]),
            all(panel_i$precision_weight > 0), all(is.finite(panel_i$se_demixed)))
  ep$panel <- panel_i
  ep
}
subject_epochs <- setNames(lapply(participant_ids, fit_one), participant_ids)
group_panel <- data.table::rbindlist(lapply(subject_epochs, "[[", "panel"))
group_panel[, rpe_within := rpe - mean(rpe), by = .(id, run, parcel, evt_time)]


## ----group-fit, class.source="fold-show"--------------------------------------
planned <- group_panel[parcel == "Parcel_A" & evt_time == 0]
stopifnot(data.table::uniqueN(planned$id) == length(participant_ids),
          nrow(planned) == length(participant_ids) * length(run_ids) * n_trials)
group_result <- mixed_by(
  planned, model_formulae = list(
    feedback = decon_demixed ~ rpe_within + (1 + rpe_within || id) + (1 | id:run)
  ),
  weights = "precision_weight_cohort", calculate = "parameter_estimates_reml",
  padjust_by = NULL, return_models = TRUE
)
group_fit <- group_result$models$model[[1L]]$reml
knitr::kable(subset(group_result$coef_df_reml, term == "rpe_within")[
  , c("estimate", "std.error", "df", "p.value"), with = FALSE], digits = 3)
lme4::isSingular(group_fit)


## ----group-trajectory, class.source="fold-show"-------------------------------
trajectory <- medusa_multilag_infer(
  subject_epochs, covariates = "rpe", use_whitened = FALSE,
  trajectory_mode = "none", contrast_source = "unsmoothed",
  tr = tr, n_perm = 1999L, n_boot_latency = 0L, seed = 5420L
)
knitr::kable(
  subset(trajectory$coefficients, parcel == "Parcel_A")[
    , c("lag", "estimate", "conf.low.sim", "conf.high.sim", "p_maxT"), with = FALSE],
  digits = 3
)


## ----save-recipe, eval=FALSE--------------------------------------------------
# saveRDS(list(hrf_fit = hrf_fit, deconvolution = decon,
#              feedback = feedback_fit, decision = decision_fit,
#              events = events, covariates = covariates,
#              session_info = sessionInfo()), "subject_medusa.rds")


## ----nifti-matrix-recipe, eval=FALSE------------------------------------------
# nifti_paths <- c(run1 = "/study/sub01/run1_preprocessed.nii.gz",
#                  run2 = "/study/sub01/run2_preprocessed.nii.gz",
#                  run3 = "/study/sub01/run3_preprocessed.nii.gz")
# mask_path <- "/study/sub01/analysis_mask.nii.gz"
# mask_img <- RNifti::readNifti(mask_path)
# mask_array <- as.array(mask_img)
# voxel_indices <- which(is.finite(mask_array) & mask_array > 0)
# ijk <- arrayInd(voxel_indices, dim(mask_array))
# voxel_ids <- apply(ijk, 1L, function(x) paste(x, collapse = "_"))
# 
# # Prefer a valid sform, otherwise qform. Resolve inconsistent headers first.
# mask_affine <- RNifti::xform(mask_img, useQuaternionFirst = FALSE)
# stopifnot(attr(mask_affine, "code") > 0)
# bold <- lapply(nifti_paths, function(path) {
#   img <- RNifti::readNifti(path)
#   affine <- RNifti::xform(img, useQuaternionFirst = FALSE)
#   stopifnot(identical(as.integer(dim(img)[1:3]), as.integer(dim(mask_img))),
#             attr(affine, "code") > 0,
#             isTRUE(all.equal(as.numeric(RNifti::pixdim(img)[1:3]),
#                              as.numeric(RNifti::pixdim(mask_img)[1:3]))),
#             isTRUE(all.equal(unname(affine), unname(mask_affine),
#                              check.attributes = FALSE)))
#   # Public RNifti/R array route: reshape as VOXEL x TIME, select, transpose.
#   # This materializes the image; use a modest mask/run size for this recipe.
#   volume <- matrix(as.numeric(img), nrow = prod(dim(img)[1:3]))
#   y <- t(volume[voxel_indices, , drop = FALSE])
#   colnames(y) <- voxel_ids
#   y
# })
# voxel_map <- data.frame(voxel = voxel_ids, linear_index = voxel_indices,
#                         i = ijk[, 1], j = ijk[, 2], k = ijk[, 3])


## ----file-workflow-recipe, eval=FALSE-----------------------------------------
# voxelwise_deconvolution(
#   niftis = unname(nifti_paths),
#   add_metadata = data.frame(id = "sub01", run = names(nifti_paths)),
#   out_dir = "/study/medusa_linear", TR = tr, atlas_files = atlas_path,
#   mask = mask_path, nprocs = 1L, algorithm = "reglin",
#   decon_settings = list(
#     kernel = hrf_library(tr, peak_delays = 6, dispersions = 1,
#                          undershoot_ratios = 6)[, 1L],
#     scaling = "noise_sd", penalty = "sobolev", lambda = 0.1,
#     backend = "cpp", prewhiten_gcv = TRUE, tail_padding = 0L
#   )
# )
# # Locate only the desired atlas's deconvolved files after the call above.
# d_files <- list.files("/study/medusa_linear/atlas/deconvolved",
#                       pattern = "_deconvolved\\.csv\\.gz$", full.names = TRUE)
# aligned <- evt_align_decon_files(
#   d_files = d_files, trial_df = trial_df,
#   alignment = list(evt_col = "feedback_onset", time_before = -2, time_after = 4,
#                    output_resolution = 1, overlap_policy = "allow",
#                    readout_kernel = "triangle", readout_width = 1.5 * tr),
#   aggregate_by = "atlas_value", tr = tr, ncpus = 1L, fail_on_error = TRUE
# )
# attr(aligned, "file_status")

