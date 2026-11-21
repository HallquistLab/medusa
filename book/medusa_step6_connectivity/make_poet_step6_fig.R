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
# Each cohort contributes one rejection indicator for each estimator.
dyn$rejections <- round(dyn$reject * dyn$n)
intervals <- t(vapply(seq_len(nrow(dyn)), function(i) {
  stats::binom.test(dyn$rejections[i], dyn$n[i])$conf.int
}, numeric(2)))
dyn$lo <- 100 * intervals[, 1L]
dyn$hi <- 100 * intervals[, 2L]

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
    subtitle = "24 simulated cohorts of 16 participants each",
    x        = "Mean estimated seed-target slope",
    y        = NULL
  ) +
  theme_s6

p_fpr <- ggplot(dyn, aes(y = method, color = method)) +
  geom_vline(xintercept = 5, linetype = "dashed", color = ink) +
  geom_segment(aes(x = lo, xend = hi, yend = method)) +
  geom_point(aes(x = 100 * reject), size = 3) +
  geom_text(aes(x = hi, label = sprintf("%d/%d", rejections, n)),
            hjust = -0.15, size = 3.2, color = ink) +
  scale_color_manual(values = cols) +
  scale_x_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
  labs(
    title = "B. Cohorts rejecting the null",
    subtitle = "Exact 95% binomial intervals. Dashed line marks 5%.",
    x = "False-positive rate (%)", y = NULL
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
