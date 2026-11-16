# Generate side-by-side "Healthy vs. Pathological" diagnostic figures for the QC Gallery

library(ggplot2)
library(data.table)
library(patchwork)

d <- readRDS("qc_gallery_examples.rds")
dir.create("figures", showWarnings = FALSE)

theme_qc <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 12.4, color = "#173145"),
    plot.subtitle = element_text(size = 10.0, color = "#495E6E"),
    strip.text = element_text(face = "bold", size = 10.6),
    legend.position = "bottom",
    plot.margin = margin(10, 14, 10, 10)
  )

# ==============================================================================
# Figure 1: Steps 0 & 1 — Clock/Units Misalignment & AR(1) rho = 0.90 Ceiling
# ==============================================================================
cc <- copy(d$clock_curves)
cc[, state := factor(state, levels = c(
  "Healthy: Intact seconds clock (valid_rows)",
  "Pathology: Dropped 3 dummy TRs (time_offset = 0)",
  "Pathology: Onsets in TR indices instead of seconds"
))]

p1a <- ggplot(cc, aes(x = lag, y = beta, color = state, linetype = state)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  geom_hline(yintercept = 0.60, linetype = "dashed", color = "#1F7A73", alpha = 0.7) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c("#1B6CA8", "#C0392B", "#B35806")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed")) +
  scale_x_continuous(breaks = seq(-6, 6, by = 2)) +
  guides(color = guide_legend(nrow = 3), linetype = guide_legend(nrow = 3)) +
  labs(
    title = "A. Step 0/3 Event-Clock & Units Signature",
    subtitle = "True lag-0 feedback slope = 0.60 at TR = 2.0 s",
    x = "Peri-event lag (s)",
    y = "Estimated RPE slope",
    color = NULL, linetype = NULL
  ) +
  theme_qc

rd <- copy(d$rho_dt)
rd[, state := factor(state, levels = c(
  "Healthy: Clean / detrended BOLD",
  "Pathology: Unmodeled shim step & mid-band drift"
))]

p1b <- ggplot(rd, aes(x = rho, fill = state)) +
  geom_vline(xintercept = 0.90, color = "#C0392B", linetype = "dashed", linewidth = 0.85) +
  geom_histogram(binwidth = 0.035, position = "identity", alpha = 0.78, color = "white") +
  annotate("text", x = 0.875, y = 26, label = "0.90 stationarity\nceiling pile-up",
           color = "#C0392B", hjust = 1, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = c("#1B6CA8", "#C0392B")) +
  scale_x_continuous(breaks = seq(0.2, 0.9, by = 0.1)) +
  coord_cartesian(xlim = c(0.20, 0.95)) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(
    title = "B. Step 1 AR(1) Autocorrelation (rho) Distribution",
    subtitle = "36 parcels (true rho = 0.30–0.62); unmodeled drift pins rho to 0.90",
    x = "Estimated hrf_fit$noise$run1$rho",
    y = "Parcel count",
    fill = NULL
  ) +
  theme_qc

ggsave("figures/01-qc-clock-and-rho.png", (p1a | p1b) + plot_layout(widths = c(1.08, 1)),
       width = 11.4, height = 5.4, dpi = 180, bg = "white")

# ==============================================================================
# Figure 2: Step 2 — Oversmoothing Gain Collapse & Uncensored Spike Ringing
# ==============================================================================
s2t <- copy(d$step2_traces)
s2t[, mode := factor(mode, levels = c(
  "Healthy: valid_rows censor + auto lambda/rwls",
  "Pathology A: Oversmoothed (lambda = 8.0; gain = 0.04)",
  "Pathology B: Uncensored motion spike at 91 s (rwls = 'none')"
))]

p2a <- ggplot(s2t, aes(x = time, y = activity, color = mode)) +
  geom_hline(yintercept = 0, color = "grey75") +
  geom_vline(xintercept = d$on2_run1, linetype = "dotted", color = "#1F7A73", linewidth = 0.75) +
  geom_vline(xintercept = d$spike_time, linetype = "dashed", color = "#C0392B", linewidth = 0.8) +
  geom_line(linewidth = 0.95) +
  annotate("text", x = d$spike_time + 0.6, y = max(s2t$activity) * 0.93,
           label = "Uncensored spike (91 s)\n-> side-lobe ringing", color = "#C0392B",
           hjust = 0, size = 2.85, fontface = "bold") +
  annotate("text", x = d$on2_run1[1] - 0.8, y = max(s2t$activity) * 0.80,
           label = "True event\n(76 s)", color = "#1F7A73", hjust = 1, size = 3.0) +
  scale_color_manual(values = c("#1B6CA8", "#B35806", "#C0392B")) +
  guides(color = guide_legend(nrow = 3)) +
  labs(
    title = "A. Continuous Deconvolved Activity Trace (67–117 s)",
    subtitle = "Dotted green = true events (76, 108 s); dashed red = motion spike (91 s)",
    x = "Run acquisition time (s)",
    y = "Restored activity â(t)",
    color = NULL
  ) +
  theme_qc

g2 <- copy(d$step2_gain_tbl)
g2[, short := factor(
  c("Healthy\n(auto lambda + valid_rows)", "Oversmoothed\n(lambda = 8.0)", "Uncensored spike\n(rwls = 'none')"),
  levels = rev(c("Healthy\n(auto lambda + valid_rows)", "Oversmoothed\n(lambda = 8.0)", "Uncensored spike\n(rwls = 'none')"))
)]
g2_long <- melt(g2, id.vars = "short", measure.vars = c("mean_gain", "gain_err"))
g2_long[, variable := factor(variable,
  levels = c("mean_gain", "gain_err"),
  labels = c("Conditional Self-Gain (higher = sharper)", "Mean Absolute Gain Attenuation (1 - gain)")
)]

p2b <- ggplot(g2_long, aes(x = value, y = short, color = variable)) +
  geom_segment(aes(x = 0, xend = value, yend = short), linewidth = 1.0, show.legend = FALSE) +
  geom_point(size = 3.1, show.legend = FALSE) +
  facet_wrap(~variable, ncol = 1, scales = "free_x") +
  scale_color_manual(values = c("#1B6CA8", "#C0392B")) +
  labs(
    title = "B. diagnose_crossrun_recovery()",
    subtitle = "Oversmoothing collapses self-gain < 0.05",
    x = "Operator recovery diagnostic",
    y = NULL
  ) +
  theme_qc

ggsave("figures/02-qc-decon-gain-and-spikes.png", (p2a | p2b) + plot_layout(widths = c(1.32, 1)),
       width = 11.4, height = 5.4, dpi = 180, bg = "white")

# ==============================================================================
# Figure 3: Steps 3 & 4 — Censoring Shadow on se_demixed & decon_eb Role Misuse
# ==============================================================================
st <- copy(d$se_trace_dt)
st[, status := factor(status, levels = c(
  "Healthy: Uncensored run (uptick only at run-end Trial 16)",
  "Scrubbed: Partial gap spikes Trial 7 SE; full gap gates Trial 11 to NA"
))]
na_trials <- st[status %like% "Scrubbed" & is.na(se_demixed), trial]
t7_se <- st[status %like% "Scrubbed" & trial == 7, se_demixed]

p3a <- ggplot(st[!is.na(se_demixed)], aes(x = trial, y = se_demixed, color = status, shape = status)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  geom_vline(xintercept = na_trials, linetype = "dashed", color = "#C0392B", alpha = 0.7) +
  annotate("text", x = 7.3, y = t7_se,
           label = "Trial 7: partial scrub\nspikes se_demixed", color = "#C0392B",
           hjust = 0, vjust = 0.7, size = 3.0, fontface = "bold") +
  annotate("text", x = 11.25, y = min(st$se_demixed, na.rm = TRUE) + 0.006,
           label = "Trial 11: low BOLD\nsupport -> gated NA", color = "#C0392B",
           hjust = 0, size = 3.0) +
  scale_color_manual(values = c("#1B6CA8", "#C0392B")) +
  scale_x_continuous(breaks = seq(1, 16, by = 2)) +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  labs(
    title = "A. Step 3 Trial-by-Trial se_demixed Audit (Run 1)",
    subtitle = "Scrubbing inside a trial's causal HRF window spikes SE or gates to NA",
    x = "Target feedback trial index",
    y = "Propagated trial SE (se_demixed)",
    color = NULL, shape = NULL
  ) +
  theme_qc

ebc <- copy(d$eb_comp_dt)
ebc[, lbl := factor(
  c("Healthy:\ndecon_demixed +\nprecision_weight_cohort", "Pathology:\nShrunken decon_eb\nentered as outcome"),
  levels = rev(c("Healthy:\ndecon_demixed +\nprecision_weight_cohort", "Pathology:\nShrunken decon_eb\nentered as outcome"))
)]

p3b <- ggplot(ebc, aes(x = estimate, y = lbl, color = lbl)) +
  geom_vline(xintercept = 0.50, linetype = "dashed", color = "#1F7A73", linewidth = 0.9) +
  geom_vline(xintercept = 0.00, linetype = "dotted", color = "grey60") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0.18, linewidth = 1.0) +
  geom_point(size = 3.3) +
  annotate("text", x = 0.50, y = 2.38, label = "True slope = 0.50", color = "#1F7A73",
           hjust = 0.5, size = 3.1, fontface = "bold") +
  scale_color_manual(values = c("#C0392B", "#1B6CA8"), guide = "none") +
  coord_cartesian(xlim = c(-0.02, 0.66), ylim = c(0.6, 2.5)) +
  labs(
    title = "B. Step 4/5 Panel Column-Role Misuse",
    subtitle = "16-subject mixed_by(): decon_eb shrinks group slopes toward zero",
    x = "Estimated population RPE slope (95% CI)",
    y = NULL
  ) +
  theme_qc

ggsave("figures/03-qc-censor-shadow-and-eb.png", (p3a | p3b) + plot_layout(widths = c(1.22, 1)),
       width = 11.4, height = 5.2, dpi = 180, bg = "white")

# ==============================================================================
# Figure 4: Steps 5 & 6 — Centering/Random-Slope Pitfalls & Connectivity Shared Noise
# ==============================================================================
s5p <- copy(d$step5_pitfalls)
s5p[, spec := factor(spec, levels = rev(spec))]

p4a <- ggplot(s5p, aes(x = estimate, y = spec, color = spec)) +
  geom_vline(xintercept = 0.45, linetype = "dashed", color = "#1F7A73", linewidth = 0.85) +
  geom_vline(xintercept = 0.00, linetype = "dotted", color = "grey60") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0.20, linewidth = 0.95) +
  geom_point(size = 3.1) +
  scale_color_manual(values = c("#B35806", "#C0392B", "#1B6CA8"), guide = "none") +
  labs(
    title = "A. Step 5 Multilevel Specification Pitfalls",
    subtitle = "True within slope = +0.45 (between = -0.35, S = 36)",
    x = "Estimated within-person slope (95% CI)",
    y = NULL
  ) +
  theme_qc

ps <- copy(d$ppi_subj_dt)
ps[, method := factor(method, levels = rev(c(
  "Healthy: M-MOMS +\nshared_error = 'auto'",
  "Pathology A: Stripped metadata /\nshared_error = 'zero' (inflated)",
  "Pathology B: Naive OLS\n(predictor-noise attenuated)"
)))]

p4b <- ggplot(ps, aes(x = estimate, y = method, color = method, fill = method)) +
  geom_vline(xintercept = 0.35, linetype = "dashed", color = "#1F7A73", linewidth = 0.85) +
  geom_vline(xintercept = 0.00, linetype = "dotted", color = "grey60") +
  geom_boxplot(width = 0.34, alpha = 0.22, outlier.shape = NA) +
  geom_jitter(height = 0.12, size = 2.2, alpha = 0.85) +
  scale_color_manual(values = c("#B35806", "#C0392B", "#1B6CA8"), guide = "none") +
  scale_fill_manual(values = c("#B35806", "#C0392B", "#1B6CA8"), guide = "none") +
  labs(
    title = "B. Step 6 Seed-Target Connectivity (ppi$subjects)",
    subtitle = "True coupling = 0.35 (shared scanner noise r = 0.55, S = 16)",
    x = "Participant connectivity slope",
    y = NULL
  ) +
  theme_qc

ggsave("figures/04-qc-residual-ac1-and-mmoms.png", (p4a | p4b) + plot_layout(widths = c(1.05, 1.05)),
       width = 11.6, height = 5.1, dpi = 180, bg = "white")

cat("Saved all 4 QC gallery figures.\n")
