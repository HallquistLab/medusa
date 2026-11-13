library(ggplot2)
a <- read.csv("diagnostic_acfs.csv")
a <- subset(a, run == "run1" & method %in% c("Canonical_HRF", "Current_MEDuSA", "Tukey15_selected"))
t_df <- subset(a, method == "Current_MEDuSA")
t_df$method <- "Generating covariance"
t_df$estimated <- t_df$true
a <- rbind(a, t_df)
a$method <- factor(
  a$method,
  levels = c("Generating covariance", "Canonical_HRF", "Current_MEDuSA", "Tukey15_selected"),
  labels = c(
    "Generating covariance",
    "Fixed canonical HRF noise fit",
    "Current MEDuSA (2-pass REML refit)",
    "Pooled Tukey taper (M = 15)"
  )
)
theme_set(
  theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 16),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(12, 20, 12, 12)
    )
)
p <- ggplot(a, aes(lag * 0.65, estimated, color = method, linetype = method)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_line(linewidth = 0.95) +
  facet_wrap(~parcel, nrow = 1) +
  scale_color_manual(values = c("#222222", "#B35806", "#2166AC", "#2A7F62")) +
  scale_linetype_manual(values = c("solid", "longdash", "dotted", "dotdash")) +
  labs(
    x = "Lag (seconds)",
    y = "Implied noise autocorrelation",
    color = NULL,
    linetype = NULL,
    title = "Refitting AR(p) noise under the learned HRF eliminates residual task-mismatch bias"
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE))
ggsave("noise_comparison.png", ggplot_gtable(ggplot_build(p)), width = 11, height = 4.8, dpi = 170, bg = "white")
