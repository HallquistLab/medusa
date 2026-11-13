library(ggplot2)
library(data.table)
d <- readRDS("step5_examples.rds")
blue <- "#245D80"; teal <- "#167C74"; orange <- "#B75126"
theme_set(theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 17, color = "#172937"),
        plot.subtitle = element_text(size = 12, color = "#405362", margin = margin(b = 10)),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom", plot.margin = margin(12, 17, 12, 12)))
save_fig <- function(name, plot, height) {
  ggsave(file.path("figures", paste0(name, ".png")), plot,
         width = 8.4, height = height, units = "in", dpi = 180, bg = "white")
}

show_ids <- levels(d$wb$id)[c(1, 5, 10, 14, 19, 24)]
wb <- d$wb[id %in% show_ids]
ws <- d$within_slopes[id %in% show_ids]
means <- d$wb[, .(x = mean(x), activity = mean(activity)), by = id]
p <- ggplot(wb, aes(x, activity)) +
  geom_point(color = blue, alpha = .26, size = 1.7) +
  geom_segment(data = ws, aes(x = x_min, xend = x_max,
    y = intercept + slope * x_min, yend = intercept + slope * x_max),
    color = blue, linewidth = 1.05, inherit.aes = FALSE) +
  geom_smooth(data = means, method = "lm", formula = y ~ x, se = FALSE,
              color = orange, linetype = "dashed", linewidth = 1.05) +
  geom_point(data = means, shape = 23, fill = orange, color = "white", size = 3.15) +
  labs(title = "Trial effects and participant differences can disagree",
       subtitle = "Blue lines show six participants. Orange diamonds show all participant means.",
       x = "Computational predictor", y = "Trial activity")
save_fig("01-within-between", p, 4.4)

n <- copy(d$null)
n[, method := factor(method, levels = c("Participant slopes", "Trial-row OLS"))]
n[, trials := factor(n_trial)]
lo <- pmax(0, n$rejection - 1.96 * n$mcse)
hi <- pmin(1, n$rejection + 1.96 * n$mcse)
p <- ggplot(n, aes(trials, rejection, group = method, color = method, shape = method)) +
  geom_hline(yintercept = .05, color = "#798A94", linetype = "dotted", linewidth = .65) +
  geom_line(linewidth = .9) + geom_point(size = 3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .12, linewidth = .7) +
  scale_color_manual(values = c(blue, orange)) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1)) +
  labs(title = "More trials do not replace more participants",
       subtitle = "True population slope is zero. The dotted line marks the 5% test level.",
       x = "Trials per participant", y = "Null rejection rate", color = NULL, shape = NULL)
save_fig("02-null", p, 4.0)

w <- unique(d$weight_fit)
true <- unique(w[, .(run, sd_model = sd_true)])
true[, model := "Generating SD"]
w <- rbind(w[, .(run, sd_model, model)], true)
w[, model := factor(model, c("Generating SD", "Weighted", "Unweighted"))]
p <- ggplot(w, aes(run, sd_model, color = model, shape = model)) +
  geom_point(size = 3.75, position = position_dodge(width = .45)) +
  scale_color_manual(values = c("#263948", teal, orange)) +
  scale_shape_manual(values = c(4, 16, 17)) +
  scale_x_discrete(labels = c(run1 = "Lower measurement noise", run2 = "Higher measurement noise")) +
  scale_y_continuous(limits = c(0, .48), breaks = seq(0, .4, .1)) +
  labs(title = "Weights specify the model’s residual spread",
       subtitle = "One parcel and lag with known heteroskedasticity",
       x = NULL, y = "Conditional residual SD", color = NULL, shape = NULL)
save_fig("03-weights", p, 3.7)

f <- copy(d$fixed)
f[, term := factor(term, c("value_within", "entropy_within"), c("Value", "Entropy"))]
f[, parcel := factor(parcel, c("Phasic", "Sustained", "Null"))]
p <- ggplot(f, aes(evt_time, estimate, color = term, fill = term)) +
  geom_hline(yintercept = 0, color = "#9BAAB3", linewidth = .35) +
  geom_vline(xintercept = 0, color = "#9BAAB3", linetype = 3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .16, color = NA) +
  geom_line(aes(y = truth), linetype = "dashed", linewidth = .65) +
  geom_line(linewidth = .8) + geom_point(size = 1.8) +
  facet_wrap(~parcel, ncol = 1) +
  scale_color_manual(values = c(blue, orange)) +
  scale_fill_manual(values = c(blue, orange)) +
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 5)) +
  labs(title = "A group slope has a profile across lags",
       subtitle = "Dashed lines are generating slopes. Ribbons are pointwise intervals.",
       x = "Time from feedback (s)", y = "Population slope", color = NULL, fill = NULL)
save_fig("04-profiles", p, 6.55)

f[, status := factor(fifelse(p_by <= .05, "BY ≤ .05", fifelse(p.value <= .05, "Pointwise only", "Neither")),
                     levels = c("BY ≤ .05", "Pointwise only", "Neither"))]
f[, parcel := factor(parcel, c("Null", "Sustained", "Phasic"))]
p <- ggplot(f, aes(evt_time, parcel, fill = status)) +
  geom_tile(color = "white", linewidth = 1.1, width = .94, height = .86) +
  facet_wrap(~term, ncol = 1) +
  scale_fill_manual(values = c("#176E66", "#D29466", "#E5ECF0"), drop = FALSE) +
  scale_x_continuous(breaks = -2:5) +
  labs(title = "The hypothesis family changes which cells pass",
       subtitle = "One family of 48 tests across parcels, lags, and both terms",
       x = "Time from feedback (s)", y = NULL, fill = NULL)
save_fig("05-family", p, 4.6)

a <- copy(d$residual_acf)
a[, label := factor(id, levels = rev(levels(id)))]
p <- ggplot(a, aes(ac1, label, color = run, shape = run)) +
  geom_vline(xintercept = 0, color = "#82939D", linetype = "dotted") +
  geom_point(position = position_dodge(width = .55), size = 2.5, alpha = .85) +
  scale_color_manual(values = c(blue, orange), labels = c("Run 1", "Run 2")) +
  scale_shape_discrete(labels = c("Run 1", "Run 2")) +
  scale_y_discrete(breaks = sprintf("P%02d", c(1, 6, 12, 18, 24))) +
  coord_cartesian(xlim = c(-.65, .65)) +
  labs(title = "Check residual dependence within each run",
       subtitle = "Phasic parcel, +1 s. Adjacent trials after the weighted model.",
       x = "Lag-one correlation of scaled conditional residuals",
       y = "Participant", color = NULL, shape = NULL)
stopifnot(all(abs(a$ac1) < .65))
save_fig("06-residuals", p, 5.0)
