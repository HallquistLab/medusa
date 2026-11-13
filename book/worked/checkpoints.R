# Reproduce from worked/ using the installed medusa package (or MEDUSA_SOURCE).
source_dir <- Sys.getenv("MEDUSA_SOURCE")
if (nzchar(source_dir)) {
  compile <- Sys.getenv("MEDUSA_COMPILE", "true") != "false"
  pkgload::load_all(source_dir, compile = compile, quiet = TRUE)
}
library(medusa)
library(data.table)
purl_path <- tempfile(fileext = ".R")
invisible(knitr::purl("workflow-source.Rmd", output = purl_path, quiet = TRUE))
script <- readLines(purl_path)
unlink(purl_path)
capture <- new.env(parent = emptyenv())
beta_line <- match("  subject <- simulate_participant(true_beta + rnorm(1L, sd = .12))", script)
return_line <- match("  ep$panel <- panel_i", script)
stopifnot(!anyNA(c(beta_line, return_line)))
script[beta_line] <- "  beta_s <- true_beta + rnorm(1L, sd = .12); subject <- simulate_participant(beta_s)"
script[return_line] <- paste(script[return_line],
  "capture[[id]] <- list(subject = subject, fit = fit, dec = dec, beta = beta_s)", sep = "; ")
eval(parse(text = script))
parcel <- "Parcel_A"
sid <- "sim01"
run_focus <- "run1"
trial_focus <- 8L
rec <- capture[[sid]]
ep <- subject_epochs[[sid]]
y <- rec$subject$bold[[run_focus]]
onset <- ep$panel[run == run_focus & trial == trial_focus & parcel == "Parcel_A" & evt_time == 0, onset]
stopifnot(length(onset) == 1L)
parcel_w <- ep$parcel_weights[[parcel]]
idx <- which(parcel_ids == parcel)
time <- (seq_len(nrow(y)) - 1L) * tr - onset
activity <- rec$dec$activity[[run_focus]][, idx, drop = FALSE]
sc <- rec$dec$scaling[[run_focus]]
stage12 <- data.frame(time, bold = as.numeric((y[, idx] - matrix(colMeans(y[, idx]), nrow(y), length(idx), byrow=TRUE)) %*% parcel_w),
  recovered = as.numeric(sweep(activity, 2L, sc$scale[idx], "*") %*% parcel_w))
staged <- ep$panel[run == run_focus & trial == trial_focus & parcel == "Parcel_A",
  .(lag=evt_time, decon_raw, decon_demixed, se_demixed, precision_weight_cohort)]
event_row <- ep$panel[run == run_focus & trial == trial_focus & parcel == "Parcel_A" & evt_time == 0, event_index]
hrf_truth <- as.numeric(candidates[, true_kernel_index[idx], drop=FALSE] %*% parcel_w)
rpos <- match(run_focus, names(rec$subject$bold))
stratum <- medusa:::hrf_stratum_kernels(rec$fit, rpos)
hrf_fitted <- as.numeric(do.call(cbind,lapply(idx,stratum$get_kernel)) %*% parcel_w)
hrf <- data.frame(time = seq_along(hrf_truth)*tr-tr, generating = hrf_truth, fitted = hrf_fitted)

# Project only the generating feedback modulation through the frozen inverse.
# Preserve masks, scaling divisors, innovations and parcel weights.
make_epoch <- function(dec, rec) decon_epoch_regression(
  dec, rec$subject$covariates, target_condition = "feedback", lags = lags_group,
  parcel_groups = parcel_ids, restore_scale = TRUE, demix_within_trial = FALSE,
  tau2_mode = "marginal", return_panel = TRUE
)
project_signal <- function(rec, orig) {
  dec <- rec$dec
  for (r in seq_along(dec$activity)) {
    rn <- names(dec$activity)[r]
    n <- nrow(dec$activity[[r]])
    mask <- attr(dec$activity[[r]], "observation_mask")
    rows <- which(mask)
    fw_r <- medusa:::resolve_run_obs_weights(dec, r, length(rows))
    e <- rec$subject$events[[rn]]
    ti <- which(e$condition == "feedback")
    gen <- vapply(e$onset[ti], source_column, numeric(n), duration=0, n=n, tr=tr)
    source <- as.numeric(gen %*% rec$subject$covariates[[rn]]$rpe)
    stratum <- medusa:::hrf_stratum_kernels(rec$fit,r)
    no_nuisance <- matrix(numeric(),n,0L)
    for (v in seq_along(voxel_ids)) {
      inv <- medusa:::hrf_inverse_operator(n,as.numeric(stratum$get_kernel(v)),
        rows,stratum$get_ar(v),no_nuisance,dec$settings,obs_weights=fw_r)
      x <- matrix(rec$subject$bold[[rn]][rows,v], ncol=1L)
      scale <- dec$scaling[[rn]]$scale[v]
      reconstructed <- as.numeric(inv %*% (medusa:::center_cols_weighted(x, fw_r)/scale))
      stopifnot(max(abs(reconstructed-dec$activity[[r]][,v]))<1e-8)
      signal <- matrix(forward_bold(source*rec$beta[v],candidates[,true_kernel_index[v]])[rows], ncol=1L)
      dec$activity[[r]][,v] <- as.numeric(inv %*% (medusa:::center_cols_weighted(signal, fw_r)/scale))
    }
  }
  projected <- make_epoch(dec,rec)$panel
  key <- c("run","event_index","evt_time","parcel")
  stopifnot(identical(as.data.frame(projected[, ..key]),as.data.frame(orig$panel[, ..key])))
  orig$panel <- copy(orig$panel)
  orig$panel[,decon_demixed := projected$decon_demixed]
  orig
}
inference <- function(epochs, mode = "none", contrast = "unsmoothed") medusa_multilag_infer(
  epochs, covariates="rpe",use_whitened=FALSE,trajectory_mode=mode,
  contrast_source=contrast,contrasts=list(late_minus_early=list(early=c(0,1),late=c(2,4))),
  tr=tr,n_perm=399L,n_boot_latency=0L,seed=5420L
)
oracle <- lapply(participant_ids,function(id) project_signal(capture[[id]],subject_epochs[[id]]))
names(oracle) <- participant_ids
projected <- inference(oracle)$coefficients[parcel=="Parcel_A" & term=="rpe"]
truth <- mean(vapply(participant_ids,function(id) {
  sum(capture[[id]]$beta[idx] * subject_epochs[[id]]$parcel_weights[[parcel]])
},numeric(1)))

sources <- c("raw","unsmoothed","multilag")
compare <- lapply(sources,function(src) inference(subject_epochs,"deblur_pspline",src))
contrasts <- rbindlist(lapply(seq_along(sources),function(i) {
  cc <- compare[[i]]$contrasts[parcel=="Parcel_A" & term=="rpe"]
  stopifnot(nrow(cc)==1L)
  cc[, .(source=sources[i], estimate, std.error, conf.low, conf.high)]
}))
stopifnot(length(unique(round(contrasts$estimate,5)))==3L)
cc <- compare[[1]]$coefficients[parcel=="Parcel_A" & term=="rpe"]
unsmoothed_curve <- inference(subject_epochs,"deblur","unsmoothed")$coefficients[
  parcel=="Parcel_A" & term=="rpe",estimate_multilag]
curve_versions <- cc[, .(lag, raw=estimate, unsmoothed=unsmoothed_curve, multilag=estimate_multilag)]
stopifnot(abs(mean(unsmoothed_curve[cc$lag %in% 2:4]) -
              mean(unsmoothed_curve[cc$lag %in% 0:1]) -
              contrasts[source=="unsmoothed",estimate]) < 1e-8)

# Remove the final acquisition support in one run, keeping the acquisition grid.
last <- max(ep$panel[run==run_focus,onset])
trials_at_edge <- ep$panel[run==run_focus & onset==last & parcel=="Parcel_A", trial]
edge_trial <- unique(trials_at_edge)
fit_gap <- rec$fit
edge <- NULL
for (extra in 1:12) {
  fit_gap$specification$valid_rows[[run_focus]] <-
    (seq_len(nrow(y))-1)*tr <= last+extra
  dec_gap <- deconvolve_crossrun_hrf(rec$subject$bold,fit_gap,lambda=.1,scaling="noise_sd",
    penalty="sobolev",nuisance_mode="ignore",observation_mask="fit")
  p_gap <- make_epoch(dec_gap,rec)$panel
  edge <- p_gap[run==run_focus & trial==edge_trial & parcel=="Parcel_A", .(lag=evt_time,estimate=decon_demixed)]
  if (is.finite(edge[lag==0,estimate]) && any(!is.finite(edge$estimate))) break
}
stopifnot(is.finite(edge[lag==0,estimate]), anyNA(edge$estimate))
partial <- subject_epochs
partial[[sid]]$panel <- copy(partial[[sid]]$panel)
partial[[sid]]$panel[evt_time==max(lags_group), decon_demixed:=NA_real_]
partial_fit <- inference(partial)
support <- partial_fit$subjects[parcel=="Parcel_A" & term=="rpe",
  .(n_subjects=sum(is.finite(estimate_raw) & is.finite(std.error_raw))),by=lag]
support[, group_total := partial_fit$coefficients[parcel=="Parcel_A" & term=="rpe",n_subjects]]
setorder(support,lag)
stopifnot(support[lag==max(lags_group),n_subjects]==length(participant_ids)-1,
          all(support[lag<max(lags_group),n_subjects]==length(participant_ids)))

# Keep scalar loading columns in the round-trip CSV. A minimal CSV drops them.
csv_path <- tempfile(fileext=".csv")
fwrite(group_panel,csv_path)
csv <- fread(csv_path)
unlink(csv_path)
stopifnot(nrow(csv)==nrow(group_panel))
minimal <- copy(csv)
minimal[, shared_noise_loading := NULL]
views <- list("Participant fits"=subject_epochs,"Stacked panel"=group_panel,
              "Full CSV"=csv,"CSV without loadings"=minimal)
fits <- lapply(views,function(x) decon_connectivity_repair(
  x,seed_region="Parcel_A",target_region="Parcel_B",covariates="rpe",
  test_ppi=TRUE,use_whitened=FALSE,orthogonalize_seed_x=FALSE,shared_error="auto",
  method="mmoms",region_col="parcel"))
handoff <- rbindlist(lapply(seq_along(views),function(i) {
  x <- views[[i]]
  f <- fits[[i]]
  g <- f$group[term=="seed_c:x_c"]
  stopifnot(nrow(g)==1L)
  ss <- f$subjects[term=="seed_c:x_c" & converged==TRUE]
  list(input=names(views)[i],
       covariance=if (is.list(x) && !is.data.frame(x)) "Participant matrices"
         else if(!is.null(attr(x,"unit_error_cor"))) "Single panel attribute"
         else if("shared_noise_loading" %in% names(x)) "Loading products" else "Zero fallback",
       median_r=median(ss$r_shared_error), estimate=g$estimate,
       lower=g$conf.low,upper=g$conf.high,n_subjects=g$n_subjects)
}))
stopifnot(handoff[input=="CSV without loadings",median_r]==0)
current_curve <- trajectory$coefficients[parcel=="Parcel_A" & term=="rpe"]
current_focal <- group_result$coef_df_reml[term=="rpe_within"]
stopifnot(nrow(current_focal)==1L,nrow(current_curve)==length(lags_group),nrow(planned)==720L)
out <- list(group_curve=current_curve,focal=current_focal,
            n_trials=nrow(planned),subjects=length(participant_ids),singular=lme4::isSingular(group_fit),
            hrf=hrf,stage12=stage12,trial=staged,event_row=event_row,onset=onset,
            hrf_training=rec$fit$training_runs[[run_focus]], parcel_weights=parcel_w,
            edge=edge,edge_trial=edge_trial,edge_cut=extra,support=support,truth=truth,
            projected=projected[,.(lag,estimate)],contrasts=contrasts,
            curve_versions=curve_versions,handoff=handoff)
saveRDS(out,"checkpoints.rds")
print(list(edge=edge,support=support,contrasts=contrasts, handoff=handoff))
