library(ggplot2)
library(patchwork)

d <- read.csv("poet_covreg_step6_summary.csv")
dyn <- subset(d, grepl("Dynamic Shared-Noise", regime) & term == "seed_c")

ink    <- "#172937"
blue   <- "#245D80"
teal   <- "#167C74"
orange <- "#B75126"
amber  <- "#C88A2E"
gray   <- "#637282"

ord <- c(
  "Zero Shared Error",
  "Legacy Stationary (auto)",
  "POET + CovReg (CSV Panel)",
  "POET + CovReg (Operator)"
)
dyn$method <- factor(dyn$method, levels = ord)
dyn$ci_lo  <- dyn$mean - 1.96 * dyn$mean_se / sqrt(dyn$n)
dyn$ci_hi  <- dyn$mean + 1.96 * dyn$mean_se / sqrt(dyn$n)

cols <- c(
  "Zero Shared Error"        = orange,
  "Legacy Stationary (auto)" = amber,
  "POET + CovReg (CSV Panel)" = blue,
  "POET + CovReg (Operator)"  = teal
)

theme_s6 <- theme_minimal(base_size = 11.5) +
  theme(
    plot.title       = element_text(face = "bold", size = 12.0, color = ink),
    plot.subtitle    = element_text(size = 9.2, color = "#405362", margin = margin(b = 6)),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    plot.margin      = margin(8, 12, 8, 8)
  )

p_bias <- ggplot(dyn, aes(x = mean, y = method, color = method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#5C6E7C", linewidth = 0.7) +
  geom_segment(aes(x = 0, xend = mean, yend = method), linewidth = 1.1, alpha = 0.75) +
  geom_point(size = 3.5) +
  geom_text(aes(label = sprintf("%+.3f", mean)),
            hjust = ifelse(dyn$mean >= 0, -0.30, 1.30),
            size = 3.3, fontface = "bold", color = ink) +
  scale_color_manual(values = cols) +
  scale_x_continuous(limits = c(-0.06, 0.215), breaks = seq(-0.05, 0.20, by = 0.05)) +
  labs(
    title    = "A. Null connectivity bias (truth = 0)",
    subtitle = "Dynamic RVT-modulated shared noise (S = 16 cohorts)",
    x        = "Mean estimated seed-target slope",
    y        = NULL
  ) +
  theme_s6

p_fpr <- ggplot(dyn, aes(x = 100 * reject, y = method, fill = method)) +
  geom_vline(xintercept = 5, linetype = "dashed", color = ink, linewidth = 0.75) +
  geom_col(width = 0.56, alpha = 0.90) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * reject)),
            hjust = -0.18, size = 3.3, fontface = "bold", color = ink) +
  annotate("text", x = 7.2, y = 4.38, label = "5% nominal", hjust = 0,
           size = 2.9, color = ink, fontface = "italic") +
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(0, 90), breaks = c(0, 5, 25, 50, 75)) +
  labs(
    title    = "B. Null false-positive rate (p <= .05)",
    subtitle = "Operator & CSV panel modes match identically (4.2%)",
    x        = "False-positive rate (%)",
    y        = NULL
  ) +
  theme_s6 +
  theme(axis.text.y = element_blank())

fig <- p_bias + p_fpr + plot_layout(widths = c(1.15, 1.0))

ggsave(
  "poet_covreg_step6.png",
  fig,
  width = 9.8, height = 3.85, units = "in", dpi = 180, bg = "white"
)
cat("Saved poet_covreg_step6.png\n")
