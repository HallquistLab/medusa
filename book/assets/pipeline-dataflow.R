# Generate MEDuSA visual data-flow, operator, and uncertainty schematic
# Output: assets/pipeline-dataflow.png

library(ggplot2)
library(grid)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
out_dir <- if (length(script_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
} else if (dir.exists("assets")) {
  "assets"
} else {
  "."
}

# Swimlane background bands (Columns 1-4: participant-level stages)
bands <- data.frame(
  xmin  = c(0.8, 0.8, 0.8),
  xmax  = c(73.2, 73.2, 73.2),
  ymin  = c(40.2, 20.6, 1.0),
  ymax  = c(57.8, 38.2, 18.6),
  fill  = c("#F1F5F9", "#EDF7F6", "#FBF5EE"),
  edge  = c("#CBD5E1", "#B8DCD8", "#E9D5BE"),
  label = c(
    "DATA & S3 OBJECTS  (Representation & Matrix Dimensions)",
    "LINEAR OPERATORS & RESOLUTION  (Forward & Inverse Transformations)",
    "NOISE & PROPAGATED COVARIANCE  (Serial, Spatial, & Trial Uncertainty)"
  ),
  lbl_col = c("#1E3A5F", "#175954", "#8C431E"),
  stringsAsFactors = FALSE
)

# Stage column headers
col_x <- c(10.1, 28.0, 45.9, 63.8)
stage_hdrs <- data.frame(
  x = col_x,
  y = 59.8,
  title = c(
    "Step 0 · BIDS Ingestion",
    "Step 1 · Cross-Run HRF & Noise",
    "Step 2 · Continuous Decon",
    "Steps 3–4 · Epochs & Parcels"
  ),
  stringsAsFactors = FALSE
)

# 12 cards (4 columns x 3 lanes)
w_half <- 7.95
cards <- data.frame(
  x = rep(col_x, each = 3L),
  y = rep(c(47.8, 28.2, 8.6), times = 4L),
  lane = rep(c("data", "op", "noise"), times = 4L),
  head = c(
    # Step 0
    "read_fmriprep_subject()",
    "Intact Clock & Row Mask",
    "Confound & Censor Bundle",
    # Step 1
    "estimate_crossrun_hrf()",
    "Out-of-Sample HRF Filter",
    "Own-Run 2-Pass REML AR(p)",
    # Step 2
    "deconvolve_crossrun_hrf()",
    "Smoothness-Penalized GLS",
    "Event-Projected Noise Scale",
    # Steps 3-4
    "decon_epoch_regression()",
    "Readout & (T, O) De-Mixing",
    "Analytical Trial Covariance"
  ),
  body = c(
    # Step 0
    "S3 \"medusa_nifti_ts\"\nbold[[r]]: T_r × V matrices\nevents, parcel_ids, geometry",
    "Scanner clock t = 0..(T_r − 1)TR\nValid-frame selector M_{r,o}\n(m_r unmasked frames I_r)",
    "Nuisance Z_r (24p / CompCor /\nnon-aggressive ICA-AROMA)\nvalid_rows (FD & dummy TRs)",
    # Step 1
    "S3 \"crossrun_hrf\" (hrf_fit)\nLeave-one-run-out folds s ≠ r\nParcel or cortical-graph pool",
    "Selected kernels h^(−r) (L_hrf × V)\nCausal Toeplitz matrix H_r\n(T_r × T_r per stratum)",
    "Stationary AR(p) model Σ_cor\nExact gapped whitener W_r\n(m_r × m_r on retained I_r)",
    # Step 2
    "S3 \"crossrun_decon\" (decon)\nactivity[[r]]: â_r (T_r × V)\nMonotone SNR-tier penalty λ_v",
    "Linear inverse D_r (T_r × m_r)\nâ_r = D_r ỹ_{r,o}\nResolution R_r = D_r M_{r,o} H_r",
    "Outlier-trimmed GLS scale σ_obs\nTask-free noise_innovations\n(m_r × V) + glitch weights w_t",
    # Steps 3-4
    "epochs$panel (data.table)\ntrial × lag × parcel rows\ndecon_demixed: ẑ = θ_T + e",
    "Local readout S_τ (m × T_r)\nComposite trial map B_{r,τ}\nCross-lag matrix lag_resolution",
    "Ω_{r,τ} → se_demixed (s_{ri})\nprecision_weight_cohort (w_{ri})\nunit_error_cor / shared_loading"
  ),
  stringsAsFactors = FALSE
)

cards$xmin <- cards$x - w_half
cards$xmax <- cards$x + w_half
cards$ymin <- cards$y - 5.9
cards$ymax <- cards$y + 5.9
cards$border <- ifelse(cards$lane == "data", "#2B5C8F",
                ifelse(cards$lane == "op", "#1F7A73", "#B5622B"))

# Horizontal arrows between columns 1->2, 2->3, 3->4 across all 3 lanes
h_arrows <- expand.grid(
  idx = 1:3,
  y   = c(47.8, 28.2, 8.6)
)
h_arrows$x0 <- col_x[h_arrows$idx] + w_half + 0.16
h_arrows$x1 <- col_x[h_arrows$idx + 1L] - w_half - 0.16
h_arrows$col <- ifelse(h_arrows$y > 40, "#2B5C8F",
                ifelse(h_arrows$y > 20, "#1F7A73", "#B5622B"))

# Vertical within-stage arrows (placed on right half of columns 3-4, clear of left band titles)
v_arrows <- data.frame(
  x  = c(49.5, 49.5, 67.4, 67.4),
  y0 = c(14.6, 34.2, 14.6, 34.2),
  y1 = c(22.1, 41.7, 22.1, 41.7)
)

# Clean straight/routed connectors across the x = 72.0..78.0 gutter
gutter_segs <- data.frame(
  x0  = c(71.9, 71.9, 71.9, 71.9, 71.9),
  y0  = c(49.2, 28.8, 10.6, 46.0,  7.0),
  x1  = c(77.9, 77.9, 77.9, 77.9, 77.9),
  y1  = c(48.0, 41.5, 35.5, 22.0, 11.5),
  col = c("#1B6CA8", "#1F7A73", "#1B6CA8", "#A74D2B", "#A74D2B"),
  lty = c("solid", "solid", "dashed", "solid", "solid"),
  stringsAsFactors = FALSE
)

# Right-hand cohort structural stage boxes (Step 5 vs Step 6)
cohort_boxes <- data.frame(
  xmin = c(78.1, 78.1),
  xmax = c(99.3, 99.3),
  ymin = c(31.2, 1.0),
  ymax = c(57.8, 27.8),
  fill = c("#EEF5FB", "#FDF4EE"),
  edge = c("#1B6CA8", "#A74D2B"),
  title = c(
    "Step 5 · Neural as Outcome",
    "Step 6 · Neural as Predictor"
  ),
  subtitle = c(
    "ẑ_{sri} ~ Behavior + (1 + Behavior || id)",
    "ẑ_{T,i} ~ ẑ_{S,i} + x_i + ẑ_{S,i} : x_i"
  ),
  body = c(
    paste(
      "• mixed_by(): Split-wise REML mixed models",
      "  weighted by precision_weight_cohort (w_{ri});",
      "  supports disk-backed Parquet streaming",
      "• medusa_multilag_infer(): Spectral deblur",
      "  via lag_resolution + P-spline + Wishart",
      "  simultaneous Max-T bands, TFCE, & contrasts",
      "• backproject_medusa() / write_nifti_ts():",
      "  3D/4D NIfTI group statistic volumes",
      sep = "\n"
    ),
    paste(
      "• decon_connectivity_repair() &",
      "  medusa_fit_mmom() (regularized M-MOMS):",
      "• Subtracts trial error blocks Ψ_i(s_{S,i}²)",
      "  from X'X (removes predictor attenuation)",
      "• Subtracts shared noise c_i(r̂_{ST}) via",
      "  unit_error_cor from X'Y (shared_error='auto')",
      "• Pools S participant slopes via finite-sample",
      "  HKSJ random-effects t-intervals",
      sep = "\n"
    )
  ),
  stringsAsFactors = FALSE
)

p <- ggplot() +
  # Swimlane backgrounds
  geom_rect(
    data = bands,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = bands$fill, color = bands$edge, linewidth = 0.55
  ) +
  geom_text(
    data = bands,
    aes(x = xmin + 1.0, y = ymax - 1.9, label = label),
    color = bands$lbl_col, hjust = 0, vjust = 0.5,
    fontface = "bold", size = 3.20
  ) +
  # Stage column headers
  geom_label(
    data = stage_hdrs,
    aes(x = x, y = y, label = title),
    fill = "#173145", color = "white", fontface = "bold",
    size = 3.40, label.padding = unit(0.32, "lines"), linewidth = 0
  ) +
  # Cohort structural header
  annotate(
    "label", x = 88.7, y = 59.8,
    label = "Cohort Structural Inference (S Participants)",
    fill = "#236B70", color = "white", fontface = "bold",
    size = 3.40, label.padding = unit(0.32, "lines"), linewidth = 0
  ) +
  # 12 participant stage cards
  geom_rect(
    data = cards,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "white", color = cards$border, linewidth = 0.7
  ) +
  geom_rect(
    data = cards,
    aes(xmin = xmin, xmax = xmax, ymin = ymax - 2.6, ymax = ymax),
    fill = cards$border, color = NA
  ) +
  geom_text(
    data = cards,
    aes(x = x, y = ymax - 1.3, label = head),
    color = "white", fontface = "bold", size = 3.02
  ) +
  geom_text(
    data = cards,
    aes(x = xmin + 0.62, y = ymin + 4.45, label = body),
    color = "#1E293B", hjust = 0, vjust = 0.5, size = 2.82, lineheight = 1.13
  ) +
  # Horizontal flow arrows
  geom_segment(
    data = h_arrows,
    aes(x = x0, xend = x1, y = y, yend = y),
    color = h_arrows$col, linewidth = 0.9,
    arrow = arrow(length = unit(0.18, "cm"), type = "closed")
  ) +
  # Vertical within-stage coupling arrows
  geom_segment(
    data = v_arrows,
    aes(x = x, xend = x, y = y0, yend = y1),
    color = "#64748B", linewidth = 0.65, linetype = "dashed",
    arrow = arrow(length = unit(0.15, "cm"), type = "closed")
  ) +
  # Gutter connectors into Step 5 & Step 6
  geom_segment(
    data = gutter_segs,
    aes(x = x0, xend = x1, y = y0, yend = y1),
    color = gutter_segs$col, linetype = gutter_segs$lty, linewidth = 0.88,
    arrow = arrow(length = unit(0.19, "cm"), type = "closed")
  ) +
  # Cohort Step 5 & Step 6 boxes
  geom_rect(
    data = cohort_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = cohort_boxes$fill, color = cohort_boxes$edge, linewidth = 0.95
  ) +
  geom_rect(
    data = cohort_boxes,
    aes(xmin = xmin, xmax = xmax, ymin = ymax - 3.3, ymax = ymax),
    fill = cohort_boxes$edge, color = NA
  ) +
  geom_text(
    data = cohort_boxes,
    aes(x = (xmin + xmax) / 2, y = ymax - 1.65, label = title),
    color = "white", fontface = "bold", size = 3.32
  ) +
  geom_text(
    data = cohort_boxes,
    aes(x = (xmin + xmax) / 2, y = ymax - 5.1, label = subtitle),
    color = cohort_boxes$edge, fontface = "bold", size = 2.95
  ) +
  geom_text(
    data = cohort_boxes,
    aes(x = xmin + 0.72, y = ymin + 10.0, label = body),
    color = "#1E293B", hjust = 0, vjust = 0.5, size = 2.80, lineheight = 1.15
  ) +
  coord_cartesian(xlim = c(0.4, 99.6), ylim = c(0.4, 61.5), expand = FALSE) +
  labs(
    title = "MEDuSA Data-Structure, Linear-Operator, and Uncertainty Architecture",
    subtitle = "Participant stages (Steps 0–4) preserve explicit forward/inverse operators and observation-noise covariances before branching by the structural role of neural activity (Steps 5–6)"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, color = "#173145", margin = margin(b = 3)),
    plot.subtitle = element_text(size = 10.3, color = "#495E6E", margin = margin(b = 8)),
    plot.margin = margin(10, 10, 8, 10),
    plot.background = element_rect(fill = "white", color = NA)
  )

out_file <- file.path(out_dir, "pipeline-dataflow.png")
ggsave(out_file, p, width = 13.6, height = 7.6, dpi = 200, bg = "white")
cat("Saved schematic to", out_file, "\n")
