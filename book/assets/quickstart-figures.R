# Generate publication-grade visual schematics and live signal walkthrough figures
# for book/quickstart.qmd
# Outputs:
#   - book/assets/quickstart-decision-tree.png
#   - book/assets/quickstart-step-mechanics.png
#   - book/assets/quickstart_metrics.rds

suppressPackageStartupMessages({
  library(medusa)
  library(ggplot2)
  library(grid)
  library(data.table)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
out_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
} else if (dir.exists("/Users/hallquist/Data_Analysis/r_packages/medusa/book/assets")) {
  "/Users/hallquist/Data_Analysis/r_packages/medusa/book/assets"
} else if (dir.exists("book/assets")) {
  "book/assets"
} else if (dir.exists("assets")) {
  "assets"
} else {
  "."
}

# Book palette
ink    <- "#173145"
blue   <- "#245D80"
teal   <- "#176E66"
orange <- "#A74D2B"
slate  <- "#526575"

# ==============================================================================
# FIGURE 1: Empirical Setup & Single-Decon Branching Architecture
# ==============================================================================
make_decision_tree_plot <- function() {
  stream_cards <- data.frame(
    xmin = 1.2, xmax = 29.2,
    ymin = c(43.5, 23.0, 2.5),
    ymax = c(58.5, 38.0, 17.5),
    hdr_fill = c("#245D80", "#176E66", "#5B6B79"),
    box_fill = c("#F3F7FA", "#EFF7F6", "#F6F8F9"),
    edge_col = c("#9BB8CD", "#96C8C3", "#C5D0D8"),
    title = c(
      "Scenario 1 · Raw fMRIPrep (+ AROMA)",
      "Scenario 2 · Postprocessed + Smoothed",
      "Scenario 3 · Custom 4D NIfTIs"
    ),
    body = c(
      "read_fmriprep_subject(desc = \"preproc\")\n• No AROMA: strategy = \"24p\"\n• With AROMA: aroma = \"nonaggressive\",\n  aroma_clean_nuisance = TRUE, strategy = \"csf_wm\"",
      "read_fmriprep_subject(desc = \"postproc\")\n• Already denoised: aroma = \"none\", strategy = \"none\"\n• Censor spikes via valid_rows (no volume deletion)\n• Smoothed BOLD: use spatial_fwhm = 4 in Step 1",
      "read_nifti_ts(niftis, mask, atlas, tr)\n• Attach read_bids_events() & confounds\n• Intersect run masks; carry voxel_map & parcel_ids"
    ),
    stringsAsFactors = FALSE
  )

  spine_bg <- data.frame(
    xmin = 33.8, xmax = 61.8, ymin = 6.5, ymax = 55.5
  )

  spine_cards <- data.frame(
    xmin = 35.4, xmax = 60.2,
    ymin = c(31.5, 9.0),
    ymax = c(50.5, 28.0),
    hdr_fill = c("#173145", "#173145"),
    box_fill = c("#FFFFFF", "#FFFFFF"),
    edge_col = c("#8FA6B8", "#8FA6B8"),
    title = c(
      "Step 1 · Cross-Run HRF & AR(1) Noise",
      "Step 2 · Continuous GLS Deconvolution"
    ),
    body = c(
      "estimate_crossrun_hrf(sub_obj, ...)\n• pooling = \"graph\", parcel_boundary = \"soft\"\n• amplitude_model = \"auto\" (catches zero-mean RPE)\n• noise_pcs = \"auto\", highpass = 128",
      "deconvolve_crossrun_hrf(sub_obj, hrf_fit, ...)\n• scaling = \"noise_sd\", penalty = \"sobolev\"\n• lambda_mode = \"auto\" (monotone SNR tiers)\n• rwls = \"auto\", nuisance_mode = \"preproject\""
    ),
    stringsAsFactors = FALSE
  )

  branch_cards <- data.frame(
    xmin = 66.2, xmax = 98.8,
    ymin = c(31.5, 2.5),
    ymax = c(58.5, 27.5),
    hdr_fill = c("#176E66", "#A74D2B"),
    box_fill = c("#EFF7F6", "#FBF5F0"),
    edge_col = c("#8EC3BD", "#DFB89E"),
    title = c(
      "Branch A · Parcelwise Panel (Steps 4–6)",
      "Branch B · Voxelwise Graph-EB (Steps 3 & 5)"
    ),
    body = c(
      "decon_epoch_regression(..., parcel_groups = sub$parcel_ids,\n  covariates = NULL, restore_scale = TRUE, whiten_panel = TRUE)\n→ Compact covariate-agnostic panel (.parquet / .rds)\n• Step 5 MLM: mixed_by(..., weights = \"precision_weight_cohort\")\n• Step 5 Curves: medusa_multilag_infer(trajectory_mode = \"auto\")\n• Step 6 PPI: decon_connectivity_repair() + backproject_medusa()",
      "decon_epoch_regression(..., aggregate_by_parcel = FALSE,\n  tau2_shrinkage = \"graph_eb\", restore_scale = TRUE)\n• Whole-brain maps (150k+ vox): return_panel = FALSE →\n  write_nifti_ts(t_gls, reference = decon) & 1/se_gls^2 pooling\n• Focal ROI single-voxel MLM: subset atlas_subset →\n  return_panel = TRUE, whiten_panel = TRUE → mixed_by()"
    ),
    stringsAsFactors = FALSE
  )

  all_cards <- rbind(
    transform(stream_cards, hdr_h = 3.9),
    transform(spine_cards, hdr_h = 3.9),
    transform(branch_cards, hdr_h = 4.1)
  )

  arrows_df <- data.frame(
    x    = c(29.4, 29.4, 29.4, 47.8, 62.0, 62.0),
    y    = c(51.0, 30.5, 10.0, 31.3, 39.0, 19.5),
    xend = c(33.5, 33.5, 33.5, 47.8, 65.9, 65.9),
    yend = c(42.0, 30.5, 19.0, 28.3, 45.0, 15.0)
  )

  ggplot() +
    annotate("text", x = 15.2, y = 61.5, label = "1. Choose Ingestion Recipe (Step 0)",
             fontface = "bold", size = 3.7, color = ink) +
    annotate("text", x = 47.8, y = 61.5, label = "2. Run Once per Subject (Steps 1–2)",
             fontface = "bold", size = 3.7, color = ink) +
    annotate("text", x = 82.5, y = 61.5, label = "3. Branch from Single decon Object (Steps 3–6)",
             fontface = "bold", size = 3.7, color = ink) +

    geom_rect(data = spine_bg,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#EEF3F7", color = "#B8C7D4", linewidth = 0.6, linetype = "dashed") +
    annotate("text", x = 47.8, y = 53.2,
             label = "SHARED CONTINUOUS FIT  (Single decon S3 Bundle)",
             fontface = "bold", size = 2.85, color = blue) +

    geom_rect(data = all_cards,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = all_cards$box_fill, color = all_cards$edge_col, linewidth = 0.55) +
    geom_rect(data = all_cards,
              aes(xmin = xmin, xmax = xmax, ymin = ymax - hdr_h, ymax = ymax),
              fill = all_cards$hdr_fill, color = NA) +
    geom_text(data = all_cards,
              aes(x = xmin + 0.9, y = ymax - 0.5 * hdr_h, label = title),
              hjust = 0, vjust = 0.5, color = "white", fontface = "bold", size = 3.0) +
    geom_text(data = all_cards,
              aes(x = xmin + 0.9, y = ymax - hdr_h - 1.1, label = body),
              hjust = 0, vjust = 1, color = ink, size = 2.62, lineheight = 1.22,
              family = "sans") +

    geom_segment(data = arrows_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 color = slate, linewidth = 0.75,
                 arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +

    coord_cartesian(xlim = c(0.5, 99.5), ylim = c(1.0, 63.2), expand = FALSE) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(6, 6, 6, 6)
    )
}

p_tree <- make_decision_tree_plot()
ggsave(
  filename = file.path(out_dir, "quickstart-decision-tree.png"),
  plot = p_tree, width = 12.4, height = 6.8, dpi = 260, bg = "white"
)

# ==============================================================================
# FIGURE 2: Live 4-Panel Signal Walkthrough Across Steps 0 -> 6
# ==============================================================================
set.seed(20261118L)
tr <- 1.0
n_time <- 195L
n_runs <- 3L
n_trials <- 20L

lib <- hrf_library(tr = tr, preset = "medusa")
pk_times <- apply(lib, 2L, function(k) (which.max(k) - 1L) * tr)
idx_fast <- which.min(abs(pk_times - 3.0))[1L]
idx_can  <- which.min(abs(pk_times - 6.0))[1L]
h_true   <- lib[, idx_fast, drop = FALSE]
h_can    <- lib[, idx_can, drop = FALSE]
t_hrf    <- seq(0, (nrow(lib) - 1L) * tr, by = tr)

# 6 parcels x 6 voxels = 36 voxels (P_Target, P_Seed, plus 4 background parcels)
# >= 5 voxels/parcel satisfies min_parcel_voxels = 5L and gives POET 6 parcels!
parcel_names <- c("P_Target", "P_Seed", "BG1", "BG2", "BG3", "BG4")
voxel_groups <- stats::setNames(rep(parcel_names, each = 6L), paste0("v", 1:36))
n_vox <- length(voxel_groups)

sim_subject_data <- function(n_r = 3L, inject_artifacts_r1 = FALSE) {
  bold_l <- ev_l <- val_l <- cov_l <- list()
  glitch_frame <- NA_integer_
  scrub_frames <- integer(0L)

  for (r in seq_len(n_r)) {
    rn <- paste0("run", r)
    itis <- stats::runif(n_trials - 1L, 5.2, 7.8)
    # Leave an explicit 13s quiet inter-trial gap after trial 10 in Run 1 so an isolated
    # non-task scanner glitch is clearly separated from neural bursts
    if (inject_artifacts_r1 && r == 1L) itis[10L] <- 14.0
    cue_on <- cumsum(c(8.0, itis))
    isi <- stats::runif(n_trials, 1.5, 3.8)
    fb_on <- cue_on + isi

    # Orthogonal cue value V_t and feedback RPE delta_t (E[delta_t | V_t] = 0)
    ev_val <- as.numeric(scale(stats::rnorm(n_trials)))
    rpe_raw <- stats::rnorm(n_trials)
    rpe <- as.numeric(scale(resid(stats::lm(rpe_raw ~ ev_val))))

    cue_amp <- 1.25 + 0.50 * ev_val
    seed_fb <- 0.25 + stats::rnorm(n_trials, sd = 0.65)
    targ_fb <- 0.35 + 0.70 * rpe + 0.50 * (seed_fb - mean(seed_fb)) + stats::rnorm(n_trials, sd = 0.18)

    ev_df <- rbind(
      data.frame(onset = cue_on, duration = 0.0, condition = "cue"),
      data.frame(onset = fb_on,  duration = 0.0, condition = "feedback")
    )
    ev_df <- ev_df[order(ev_df$onset), , drop = FALSE]
    rownames(ev_df) <- NULL

    G_cue <- medusa:::hrf_event_design_cpp(n_time, tr, 0, cue_on, rep(0.0, n_trials), seq_len(n_trials), n_trials, h_true)[[1L]]
    G_fb  <- medusa:::hrf_event_design_cpp(n_time, tr, 0, fb_on,  rep(0.0, n_trials), seq_len(n_trials), n_trials, h_true)[[1L]]
    sig_t <- as.numeric(G_cue %*% cue_amp + G_fb %*% targ_fb)
    sig_s <- as.numeric(G_cue %*% (0.90 + 0.25 * ev_val) + G_fb %*% seed_fb)

    t_norm <- seq(-1, 1, length.out = n_time)
    drift <- 0.35 * (1.5 * t_norm^2 - 0.5)
    shared_eps <- as.numeric(stats::arima.sim(model = list(ar = 0.45), n = n_time, sd = 0.24))

    Y_r <- matrix(0.0, n_time, n_vox, dimnames = list(NULL, names(voxel_groups)))
    for (v in seq_len(n_vox)) {
      eps_v <- as.numeric(stats::arima.sim(model = list(ar = 0.38), n = n_time, sd = 0.16))
      p_nm <- voxel_groups[v]
      sig_v <- if (p_nm == "P_Target") sig_t else if (p_nm == "P_Seed") sig_s else 0.0
      Y_r[, v] <- 100.0 + sig_v + drift + 0.85 * shared_eps + eps_v
    }

    vm <- rep(TRUE, n_time)
    if (inject_artifacts_r1 && r == 1L) {
      # Non-task transient scanner spike in the quiet gap after trial 10
      glitch_frame <- as.integer(round(fb_on[10L] + 7.5)) + 1L
      Y_r[glitch_frame, 1:12] <- Y_r[glitch_frame, 1:12] + 2.35
      # 2-TR high-motion burst partially shadowing trial 15's HRF window (inflating SE without full loss)
      scrub_start <- min(n_time - 6L, as.integer(round(fb_on[15L] + 4.0)) + 1L)
      scrub_frames <- scrub_start:(scrub_start + 1L)
      Y_r[scrub_frames, ] <- Y_r[scrub_frames, ] + 3.40
      vm[scrub_frames] <- FALSE
    }

    bold_l[[rn]] <- Y_r
    ev_l[[rn]]   <- ev_df
    val_l[[rn]]  <- vm
    cov_l[[rn]]  <- data.frame(rpe = rpe)
  }
  list(bold = bold_l, events = ev_l, valid_rows = val_l, covariates = cov_l,
       glitch_frame = glitch_frame, scrub_frames = scrub_frames)
}

ex1 <- sim_subject_data(n_r = 3L, inject_artifacts_r1 = TRUE)

hrf_fit <- estimate_crossrun_hrf(
  bold = ex1$bold, events = ex1$events, tr = tr, kernel_grid = lib,
  valid_rows = ex1$valid_rows, spatial_groups = voxel_groups, pooling = "parcel",
  amplitude_model = "auto", selection_mode = "soft_average",
  noise_model = "ar1", evidence = "calibrated", noise_pcs = "auto", highpass = 128
)

decon_rwls <- deconvolve_crossrun_hrf(
  bold = ex1$bold, hrf_fit = hrf_fit, scaling = "noise_sd",
  penalty = "sobolev", lambda_mode = "auto", rwls = "auto",
  nuisance_mode = "preproject", noise_pcs = "auto",
  observation_mask = "fit", highpass = 128
)
decon_norwls <- deconvolve_crossrun_hrf(
  bold = ex1$bold, hrf_fit = hrf_fit,
  lambda = decon_rwls$settings$lambda_by_voxel, scaling = "noise_sd",
  penalty = "sobolev", rwls = "none",
  nuisance_mode = "preproject", noise_pcs = "auto",
  observation_mask = "fit", highpass = 128
)

ep_par <- decon_epoch_regression(
  deconvolution = decon_rwls, covariates = ex1$covariates, events = ex1$events,
  target_condition = "feedback", lags = -2:4, parcel_groups = voxel_groups,
  other_event_model = "auto", demix_ridge = "auto", restore_scale = TRUE,
  whiten_panel = TRUE, return_panel = TRUE
)

# --- Panel A: HRF curves ---
h_est <- rowMeans(hrf_fit$selected_kernels$run1[, 1:6, drop = FALSE])
norm_peak <- function(x) x / max(abs(x))
df_pA <- rbind(
  data.frame(sec = t_hrf, val = norm_peak(h_true[, 1L]), series = "True regional HRF (3.0 s peak)"),
  data.frame(sec = t_hrf, val = norm_peak(h_est),        series = "Step 1 cross-fitted HRF (runs 2–3)"),
  data.frame(sec = t_hrf, val = norm_peak(h_can[, 1L]),  series = "Assumed canonical GLM HRF (6.0 s)")
)
df_pA <- subset(df_pA, sec <= 22)
df_pA$series <- factor(df_pA$series, levels = unique(df_pA$series))

pA <- ggplot(df_pA, aes(sec, val, color = series, linetype = series)) +
  geom_hline(yintercept = 0, color = "#CBD5E1", linewidth = 0.4) +
  geom_line(linewidth = 1.05) +
  scale_color_manual(values = c(orange, teal, "#7F8EA0")) +
  scale_linetype_manual(values = c("dashed", "solid", "dotted")) +
  labs(
    title = "A. Step 1 · Out-of-Sample HRF Identification",
    subtitle = "Learns fast 3.0 s regional timing from runs 2–3 without circularity",
    x = "Seconds post-impulse", y = "Peak-scaled HRF", color = NULL, linetype = NULL
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", color = ink, size = 11),
    plot.subtitle = element_text(color = slate, size = 9),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8.1),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE))

# --- Panel B: Window around glitch & scrub burst in Run 1 ---
gl_s <- ex1$glitch_frame - 1L
sc_s <- ex1$scrub_frames - 1L
win_lo <- max(1L, ex1$glitch_frame - 18L)
win_hi <- min(n_time, max(ex1$scrub_frames) + 12L)
win_t  <- win_lo:win_hi

y_obs_c <- rowMeans(ex1$bold$run1[, 1:6]) - mean(ex1$bold$run1[ex1$valid_rows$run1, 1:6])
a_rwls  <- rowMeans(decon_rwls$activity$run1[, 1:6])
a_unw   <- rowMeans(decon_norwls$activity$run1[, 1:6])

# Scale both decon traces by a shared reference SD so RWLS suppression is on identical units
sd_ref_a <- stats::sd(a_rwls[ex1$valid_rows$run1])
sd_ref_y <- stats::sd(y_obs_c[ex1$valid_rows$run1])

df_pB <- rbind(
  data.frame(sec = win_t - 1, val = y_obs_c[win_t] / sd_ref_y, series = "Observed BOLD (delayed + glitch + scrub)"),
  data.frame(sec = win_t - 1, val = a_unw[win_t] / sd_ref_a,   series = "Unweighted decon (rings at non-task glitch)"),
  data.frame(sec = win_t - 1, val = a_rwls[win_t] / sd_ref_a,  series = "Step 2 GLS + rwls = \"auto\" (glitch suppressed)")
)
df_pB$series <- factor(df_pB$series, levels = unique(df_pB$series))
fb_in_win <- ex1$events$run1$onset[ex1$events$run1$condition == "feedback" &
                                     ex1$events$run1$onset >= (win_lo - 1) & ex1$events$run1$onset <= (win_hi - 1)]

pB <- ggplot(df_pB, aes(sec, val, color = series, linetype = series)) +
  annotate("rect", xmin = min(sc_s) - 0.5, xmax = max(sc_s) + 0.5,
           ymin = -Inf, ymax = Inf, fill = "#FDE8E8", alpha = 0.85) +
  geom_vline(xintercept = fb_in_win, color = "#D0DCE5", linetype = "dashed", linewidth = 0.45) +
  geom_line(linewidth = 0.85) +
  annotate("text", x = gl_s, y = 3.55, label = "Non-task glitch\n(downweighted)",
           size = 2.65, color = orange, fontface = "bold", lineheight = 0.95) +
  annotate("text", x = mean(sc_s), y = -2.85, label = "valid_rows\nscrubbed",
           size = 2.65, color = "#9B2C2C", fontface = "bold", lineheight = 0.95) +
  scale_color_manual(values = c("#8A99A8", orange, teal)) +
  scale_linetype_manual(values = c("solid", "dotted", "solid")) +
  coord_cartesian(ylim = c(-3.4, 4.1)) +
  labs(
    title = "B. Step 0 & Step 2 · Intact-Clock Censoring & Continuous GLS",
    subtitle = "Re-aligns BOLD to onsets (dashed), suppresses glitches, & bridges scrubbed TRs",
    x = "Run 1 acquisition time (s)", y = "Standardized amplitude", color = NULL, linetype = NULL
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", color = ink, size = 11),
    plot.subtitle = element_text(color = slate, size = 9),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8.0),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE))

# --- Panel C: De-mixed trial readouts & scrub-shadowed SEs ---
pan_r1 <- as.data.table(ep_par$panel)[parcel == "P_Target" & run == "run1" & abs(evt_time) < 1e-6]
fb_r1 <- ex1$events$run1[ex1$events$run1$condition == "feedback", ]
pan_r1[, fb_onset := fb_r1$onset[trial]]
pan_r1 <- pan_r1[is.finite(decon_demixed) & is.finite(se_demixed)]
med_se <- stats::median(pan_r1$se_demixed, na.rm = TRUE)
pan_r1[, near_scrub := se_demixed > 1.06 * med_se]
pan_r1[, status := ifelse(near_scrub,
                          "HRF overlaps scrub gap (inflated SE, lower weight)",
                          "Clean BOLD support (tight SE, full weight)")]

pC <- ggplot(pan_r1, aes(x = fb_onset, y = decon_demixed, color = status)) +
  annotate("rect", xmin = min(sc_s) - 0.5, xmax = max(sc_s) + 0.5,
           ymin = -Inf, ymax = Inf, fill = "#FDE8E8", alpha = 0.85) +
  geom_hline(yintercept = 0, color = "#CBD5E1", linewidth = 0.4) +
  geom_errorbar(aes(ymin = decon_demixed - se_demixed, ymax = decon_demixed + se_demixed),
                width = 2.0, linewidth = 0.65) +
  geom_point(aes(size = precision_weight), alpha = 0.92) +
  scale_color_manual(values = c(
    "Clean BOLD support (tight SE, full weight)" = blue,
    "HRF overlaps scrub gap (inflated SE, lower weight)" = orange
  )) +
  scale_size_continuous(range = c(1.7, 3.6), guide = "none") +
  labs(
    title = "C. Steps 3–4 · De-Mixed Trial Readouts & Propagated SEs",
    subtitle = "Subtracts cue bleed (other_event_model = \"auto\") & weights trials by exact support",
    x = "Feedback trial onset in Run 1 (s)", y = "De-mixed readout (± 1 SE)", color = NULL
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", color = ink, size = 11),
    plot.subtitle = element_text(color = slate, size = 9),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8.0),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# --- Panel D: 14-subject cohort for Step 5 RPE lag profile & Step 6 M-MOMS connectivity ---
n_sub_coh <- 14L
sub_ep_list <- vector("list", n_sub_coh)
names(sub_ep_list) <- sprintf("sub%02d", seq_len(n_sub_coh))
naive_conn <- numeric(n_sub_coh)

for (s in seq_len(n_sub_coh)) {
  sim_s <- sim_subject_data(n_r = 3L, inject_artifacts_r1 = FALSE)
  hf_s <- estimate_crossrun_hrf(
    sim_s$bold, sim_s$events, tr, lib,
    spatial_groups = voxel_groups, pooling = "parcel",
    amplitude_model = "auto", noise_pcs = "auto", highpass = 128
  )
  dc_s <- deconvolve_crossrun_hrf(
    sim_s$bold, hf_s, scaling = "noise_sd", penalty = "sobolev",
    lambda_mode = "auto", nuisance_mode = "preproject",
    noise_pcs = "auto", highpass = 128
  )
  ep_s <- decon_epoch_regression(
    dc_s, covariates = sim_s$covariates, events = sim_s$events,
    target_condition = "feedback", lags = -2:4,
    parcel_groups = voxel_groups, other_event_model = "auto",
    demix_ridge = "auto", restore_scale = TRUE, return_panel = TRUE
  )
  ep_s$panel$id <- names(sub_ep_list)[s]
  sub_ep_list[[s]] <- ep_s

  p_s <- as.data.table(ep_s$panel)[abs(evt_time) < 1e-6]
  ps_seed <- p_s[parcel == "P_Seed", .(run, trial, seed = decon_demixed, rpe)]
  ps_targ <- p_s[parcel == "P_Target", .(run, trial, targ = decon_demixed)]
  m_st <- merge(ps_seed, ps_targ, by = c("run", "trial"))
  naive_conn[s] <- stats::coef(stats::lm(targ ~ seed + rpe + factor(run), data = m_st))[["seed"]]
}

ml_fit <- medusa_multilag_infer(
  sub_ep_list, covariates = "rpe", trajectory_mode = "deblur",
  use_whitened = FALSE, tr = tr, n_perm = 999L, n_boot_latency = 0L, seed = 42L
)

conn_fit <- decon_connectivity_repair(
  panel_dt = sub_ep_list,
  seed_region = "P_Seed",
  target_region = "P_Target",
  seed_lag = 0,
  target_lag = 0,
  covariates = "rpe",
  test_ppi = FALSE,
  use_whitened = FALSE,
  region_col = "parcel"
)

cf_curve <- as.data.table(ml_fit$coefficients)[parcel == "P_Target" & term == "rpe"]
y_col_curve <- if ("estimate_unsmoothed" %in% names(cf_curve)) "estimate_unsmoothed" else "estimate"
cf_curve[, y_plot := get(y_col_curve)]

pD <- ggplot(cf_curve, aes(x = lag, y = y_plot)) +
  geom_ribbon(aes(ymin = conf.low.sim, ymax = conf.high.sim), fill = "#CFE3EA", alpha = 0.75) +
  geom_hline(yintercept = 0, color = "#CBD5E1", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "#94A3B8", linetype = "dotted") +
  geom_line(aes(color = "Deblurred neural-time RPE slope (95% sim. band)"), linewidth = 1.05) +
  geom_point(aes(color = "Deblurred neural-time RPE slope (95% sim. band)"), size = 2.2) +
  geom_point(data = data.frame(lag = 0, y_plot = 0.70),
             aes(color = "True instantaneous RPE slope (0.70 at 0 s)"), shape = 18, size = 3.8) +
  scale_color_manual(values = c(
    "Deblurred neural-time RPE slope (95% sim. band)" = blue,
    "True instantaneous RPE slope (0.70 at 0 s)" = orange
  )) +
  scale_x_continuous(breaks = -2:4) +
  labs(
    title = "D. Step 5 · Deblurred Cohort Trajectory & Simultaneous Band",
    subtitle = "medusa_multilag_infer() sharpens cross-lag mixing onto true neural event time",
    x = "Lag from feedback onset (s)", y = "Group RPE slope", color = NULL
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", color = ink, size = 11),
    plot.subtitle = element_text(color = slate, size = 9),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8.0),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

png(file.path(out_dir, "quickstart-step-mechanics.png"), width = 3150, height = 2280, res = 260, bg = "white")
grid.newpage()
pushViewport(viewport(layout = grid.layout(2L, 2L, widths = unit(c(1, 1), "null"), heights = unit(c(1, 1), "null"))))
print(pA, vp = viewport(layout.pos.row = 1L, layout.pos.col = 1L))
print(pB, vp = viewport(layout.pos.row = 1L, layout.pos.col = 2L))
print(pC, vp = viewport(layout.pos.row = 2L, layout.pos.col = 1L))
print(pD, vp = viewport(layout.pos.row = 2L, layout.pos.col = 2L))
dev.off()

grp_dt <- as.data.table(conn_fit$group)
conn_seed_row <- grp_dt[grepl("seed", term, ignore.case = TRUE)][1L]
metrics_out <- list(
  hrf_true_peak = unname(pk_times[idx_fast]),
  hrf_est_peak  = (which.max(h_est) - 1L) * tr,
  lag0_rpe_est  = cf_curve[abs(lag) < 1e-6]$y_plot[1L],
  conn_true     = 0.50,
  conn_naive    = mean(naive_conn),
  conn_mmom     = conn_seed_row$estimate
)
saveRDS(metrics_out, file.path(out_dir, "quickstart_metrics.rds"))
cat("Generated quickstart figures and metrics successfully.\n")
print(metrics_out)
