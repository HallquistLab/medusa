# Generate figures for Step 0: Task Design & BIDS/fMRIPrep Ingestion

library(ggplot2)
library(data.table)
library(patchwork)

d <- readRDS("step0_examples.rds")
db <- fread("design_benchmarks.csv")
dir.create("figures", showWarnings = FALSE)

theme_step0 <- theme_minimal(base_size = 12.5) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 13),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    plot.margin = margin(10, 15, 10, 10)
  )

# ------------------------------------------------------------------------------
# Figure 1: Experimental Design Feasibility (Cue-Feedback Jitter & Run Count)
# ------------------------------------------------------------------------------
jit_dt <- copy(db[grepl("^Cue-Feedback", dimension)])
jit_dt[, setting := factor(setting, levels = setting)]
jit_long <- melt(
  jit_dt,
  id.vars = "setting",
  measure.vars = c("target_rmse", "cue_bleed_abs"),
  variable.name = "metric",
  value.name = "value"
)
jit_long[, metric := factor(
  metric,
  levels = c("target_rmse", "cue_bleed_abs"),
  labels = c("Target Feedback Slope RMSE", "Cue-to-Feedback Bleed (abs. slope)")
)]

p1a <- ggplot(jit_long, aes(x = value, y = factor(setting, levels = rev(levels(setting))), color = metric, shape = metric)) +
  geom_segment(aes(x = 0, xend = value, yend = factor(setting, levels = rev(levels(setting)))),
               position = position_dodge(width = 0.45), linewidth = 0.85) +
  geom_point(position = position_dodge(width = 0.45), size = 3.0) +
  scale_color_manual(values = c("#1B6CA8", "#B35806")) +
  labs(
    title = "A. Phase Jitter Decouples Overlapping Cue & Feedback Epochs",
    x = "Error / Spillover (BOLD slope units; lower is better)",
    y = NULL, color = NULL, shape = NULL
  ) +
  theme_step0

run_dt <- copy(db[grepl("^Run Structure", dimension)])
run_dt[, setting := factor(setting, levels = setting)]
run_long <- melt(
  run_dt,
  id.vars = "setting",
  measure.vars = c("hrf_mae_sec", "target_rmse"),
  variable.name = "metric",
  value.name = "value"
)
run_long[, metric := factor(
  metric,
  levels = c("hrf_mae_sec", "target_rmse"),
  labels = c("Step 1 HRF Peak Latency MAE (s)", "Step 3 Feedback Slope RMSE")
)]

p1b <- ggplot(run_long, aes(x = value, y = factor(setting, levels = rev(levels(setting))), color = metric)) +
  geom_segment(aes(x = 0, xend = value, yend = factor(setting, levels = rev(levels(setting)))),
               linewidth = 0.9, show.legend = FALSE) +
  geom_point(size = 3.0, show.legend = FALSE) +
  facet_wrap(~metric, scales = "free_x") +
  scale_color_manual(values = c("#2A7F62", "#1B6CA8")) +
  labs(
    title = "B. Splitting Total Trials Across >= 3 Runs Stabilizes Out-of-Sample HRF Cross-Fitting",
    x = "Estimation Error (lower is better)",
    y = NULL
  ) +
  theme_step0

ggsave("figures/01-design-jitter-benchmarks.png", p1a / p1b, width = 10.8, height = 7.2, dpi = 175, bg = "white")

# ------------------------------------------------------------------------------
# Figure 2: Clock Alignment — Intact valid_rows vs Deleted Dummy Scans Bug
# ------------------------------------------------------------------------------
hc <- copy(d$hrf_curves)
hc[, curve := factor(curve, levels = c(
  "Generating HRF",
  "Intact clock (`valid_rows`)",
  "Deleted 4 dummy TRs (no `time_offset`)"
))]

p2a <- ggplot(hc[time <= 22], aes(x = time, y = value, color = curve, linetype = curve)) +
  geom_hline(yintercept = 0, color = "grey75") +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = c(
    "Generating HRF" = "#2A7F62",
    "Intact clock (`valid_rows`)" = "#1B6CA8",
    "Deleted 4 dummy TRs (no `time_offset`)" = "#C0392B"
  )) +
  scale_linetype_manual(values = c("dashed", "solid", "solid")) +
  guides(color = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3)) +
  labs(
    title = "A. Step 1 Cross-Fitted HRF Waveform",
    x = "Seconds after event onset",
    y = "Normalized HRF",
    color = NULL, linetype = NULL
  ) +
  theme_step0

cp <- copy(d$clock_profiles)
cp[, mode := factor(mode, levels = c(
  "Intact clock (`valid_rows`)",
  "Truncated + `time_offset = 4`",
  "Deleted 4 dummy TRs (no `time_offset`)"
))]

p2b <- ggplot(cp, aes(x = lag, y = beta, color = mode, linetype = mode)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  geom_hline(yintercept = 0.65, linetype = "dashed", color = "#2A7F62", alpha = 0.7) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c(
    "Intact clock (`valid_rows`)" = "#1B6CA8",
    "Truncated + `time_offset = 4`" = "#2A7F62",
    "Deleted 4 dummy TRs (no `time_offset`)" = "#C0392B"
  )) +
  scale_linetype_manual(values = c("solid", "dashed", "solid")) +
  guides(color = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3)) +
  scale_x_continuous(breaks = -5:4) +
  labs(
    title = "B. Step 3 Feedback RPE Slope Across Lags",
    x = "Peri-event lag (s)",
    y = "Estimated RPE slope (true = 0.65 at lag 0)",
    color = NULL, linetype = NULL
  ) +
  theme_step0

ggsave("figures/02-clock-alignment-pitfalls.png", (p2a | p2b) + plot_layout(widths = c(1, 1.15)),
       width = 11.2, height = 5.5, dpi = 175, bg = "white")

# ------------------------------------------------------------------------------
# Figure 3: Missed/Timeout Trial Handling & Safe event_index Joins
# ------------------------------------------------------------------------------
mt <- copy(d$missed_trace)
# Zoom around first missed trial in run-01 (~50 to 82 s)
t_miss1 <- d$run1_missed_time[1L]
mt_sub <- mt[time >= (t_miss1 - 12) & time <= (t_miss1 + 18)]
mt_long <- melt(mt_sub, id.vars = "time", variable.name = "model", value.name = "activity")
mt_long[, model := factor(
  model,
  levels = c("intact", "dropped"),
  labels = c("Keep 'missed_feedback' in events (Recommended)", "Drop missed trials from events.tsv")
)]

val_in_win <- d$run1_valid_time[d$run1_valid_time >= min(mt_sub$time) & d$run1_valid_time <= max(mt_sub$time)]

p3a <- ggplot(mt_long, aes(x = time, y = activity, color = model)) +
  geom_vline(xintercept = t_miss1, color = "#C0392B", linetype = "dashed", linewidth = 0.85) +
  geom_vline(xintercept = val_in_win, color = "grey65", linetype = "dotted") +
  geom_line(linewidth = 0.95) +
  annotate("text", x = t_miss1 + 0.4, y = max(mt_long$activity) * 0.94,
           label = "Timeout / missed trial", color = "#C0392B", hjust = 0, size = 3.3, fontface = "bold") +
  scale_color_manual(values = c("#1B6CA8", "#B35806")) +
  guides(color = guide_legend(nrow = 2)) +
  labs(
    title = "A. Parcel Activity Around a Missed Trial",
    x = "Run acquisition time (s; dotted lines = valid feedback trials)",
    y = "Restored activity",
    color = NULL
  ) +
  theme_step0

ms_tbl <- copy(d$missed_tbl)
ms_tbl[, short_lbl := factor(
  c(
    "Keep 'missed_feedback' +\nevent_index join (Recommended)",
    "Delete missed trials\nfrom events.tsv",
    "Join RL table on target\ncounter `trial == trial_number`"
  ),
  levels = rev(c(
    "Keep 'missed_feedback' +\nevent_index join (Recommended)",
    "Delete missed trials\nfrom events.tsv",
    "Join RL table on target\ncounter `trial == trial_number`"
  ))
)]

p3b <- ggplot(ms_tbl, aes(x = lag0_beta, y = short_lbl, color = short_lbl)) +
  geom_vline(xintercept = 0.65, linetype = "dashed", color = "#2A7F62", linewidth = 0.85) +
  geom_vline(xintercept = 0.00, linetype = "dotted", color = "grey60") +
  geom_errorbar(aes(xmin = lag0_beta - 1.96 * lag0_se, xmax = lag0_beta + 1.96 * lag0_se),
                orientation = "y", width = 0.2, linewidth = 0.95) +
  geom_point(size = 3.1) +
  scale_color_manual(values = c("#C0392B", "#B35806", "#1B6CA8"), guide = "none") +
  labs(
    title = "B. Recovered Lag-0 RPE Slope (True = 0.65)",
    x = "Estimated RPE Slope (95% CI)",
    y = NULL
  ) +
  theme_step0

ggsave("figures/03-missed-trials-handling.png", (p3a | p3b) + plot_layout(widths = c(1.2, 1)),
       width = 11.2, height = 4.9, dpi = 175, bg = "white")

# ------------------------------------------------------------------------------
# Figure 4: Non-Aggressive ICA-AROMA & Matched Nuisance Projection
# ------------------------------------------------------------------------------
ap <- copy(d$aroma_profiles)
aroma_lvls <- c(
  "Non-aggressive AROMA + matched csf_wm (Recommended)",
  "Non-aggressive AROMA + uncleaned csf_wm",
  "Aggressive AROMA (aroma = 'aggressive')",
  "No AROMA (csf_wm only)"
)
ap[, regime := factor(regime, levels = aroma_lvls)]
aroma_cols <- c(
  "Non-aggressive AROMA + matched csf_wm (Recommended)" = "#1B6CA8",
  "Non-aggressive AROMA + uncleaned csf_wm"             = "#2A7F62",
  "Aggressive AROMA (aroma = 'aggressive')"             = "#B35806",
  "No AROMA (csf_wm only)"                              = "#C0392B"
)

p4a <- ggplot(ap, aes(x = lag, y = beta, color = regime, linetype = regime)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey65") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey65") +
  geom_hline(yintercept = 0.65, linetype = "dashed", color = "#2A7F62", alpha = 0.65) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.2) +
  scale_color_manual(values = aroma_cols) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "solid")) +
  scale_x_continuous(breaks = -2:4) +
  guides(color = guide_legend(nrow = 4), linetype = guide_legend(nrow = 4)) +
  labs(
    title = "A. Step 3 RPE Trajectory (Task-Correlated Motion)",
    x = "Peri-event lag relative to feedback (s)",
    y = "Estimated RPE slope (true = 0.65 at lag 0)",
    color = NULL, linetype = NULL
  ) +
  theme_step0 +
  theme(legend.text = element_text(size = 10.5), legend.key.height = unit(0.42, "cm"))

at <- copy(d$aroma_tbl)
short_aroma <- c(
  "No AROMA\n(csf_wm only)",
  "Aggressive AROMA\n(aroma = 'aggressive')",
  "Non-aggr. AROMA +\nuncleaned csf_wm",
  "Non-aggr. AROMA +\nmatched csf_wm (Rec.)"
)
at[, short_lbl := factor(short_aroma, levels = short_aroma)]
at[, col_hex := c("#C0392B", "#B35806", "#2A7F62", "#1B6CA8")]

p4b1 <- ggplot(at, aes(x = est_noise_sd, y = short_lbl)) +
  geom_vline(xintercept = 0.030, linetype = "dashed", color = "#2A7F62", linewidth = 0.85) +
  geom_segment(aes(x = 0, xend = est_noise_sd, yend = short_lbl), color = at$col_hex, linewidth = 0.95) +
  geom_point(color = at$col_hex, size = 3.0) +
  geom_text(
    aes(label = sprintf("R2_mot = %.2f", resid_motion_r2)),
    hjust = -0.14, size = 3.1, color = "grey25"
  ) +
  scale_x_continuous(limits = c(0, max(at$est_noise_sd) * 1.38)) +
  labs(
    title = "B. Step 2 Noise SD (True = 0.030) & Post-Nuisance Motion R2",
    x = "Estimated observation noise_sd",
    y = NULL
  ) +
  theme_step0

p4b2 <- ggplot(at, aes(x = lag0_beta, y = short_lbl)) +
  geom_vline(xintercept = 0.65, linetype = "dashed", color = "#2A7F62", linewidth = 0.85) +
  geom_errorbar(
    aes(xmin = lag0_beta - 1.96 * lag0_se, xmax = lag0_beta + 1.96 * lag0_se),
    orientation = "y", width = 0.22, linewidth = 0.95, color = at$col_hex
  ) +
  geom_point(color = at$col_hex, size = 3.0) +
  labs(
    title = "C. Step 3 Lag-0 RPE Slope & 95% CI (True = 0.65)",
    x = "Lag-0 RPE slope (95% CI)",
    y = NULL
  ) +
  theme_step0

p4_full <- (p4a | (p4b1 / p4b2)) + plot_layout(widths = c(1.15, 1.0))
ggsave("figures/04-nonaggressive-aroma.png", p4_full, width = 12.0, height = 6.5, dpi = 175, bg = "white")

cat("Saved Step 0 figures.\n")
