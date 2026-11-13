library(ggplot2)
library(data.table)
d <- readRDS("pr62_examples.rds")
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

cf <- d$result$coefficients
f <- cf[parcel == "Phasic" & term == "value_within"]
fine <- as.data.table(spline(f$lag, f$estimate_multilag, method = "natural", n = 701))
p <- ggplot(f, aes(lag, estimate)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            data = data.table(xmin = c(0, 3), xmax = c(1, 4)),
            inherit.aes = FALSE, fill = "#DAE8E1", alpha = .6) +
  geom_hline(yintercept = 0, color = "#8999A3", linetype = 3) +
  geom_ribbon(aes(ymin = conf.low.sim, ymax = conf.high.sim,
                  fill = "95% Max-|T| band"), alpha = .28) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .08, color = blue) +
  geom_point(aes(color = "Raw lag estimates"), size = 2.5) +
  geom_line(data = fine, aes(x, y, color = "Smoothed spline"), linewidth = 1.05, inherit.aes = FALSE) +
  scale_fill_manual(values = c("95% Max-|T| band" = "#8AAFC7")) +
  scale_color_manual(values = c("Raw lag estimates" = blue, "Smoothed spline" = orange)) +
  labs(title = "Raw lag tests and smoothed timing",
       subtitle = "Phasic value effect. Green shading marks the two contrast windows.",
       x = "Time from feedback (s)", y = "Population slope", color = NULL, fill = NULL)
save_fig("07-pr62-trajectory", p, 4.6)

pick <- function(p, tm, label) {
  z <- cf[parcel == p & term == tm]
  data.table(trajectory = label, test = c("Max-|T|", "TFCE", "Omnibus"),
             p = c(min(z$p_maxT), min(z$p_tfce), z$omnibus_p[1L]))
}
ps <- rbindlist(list(pick("Phasic", "value_within", "Phasic · value"),
                    pick("Sustained", "entropy_within", "Sustained · entropy"),
                    pick("Null", "entropy_within", "Null · entropy")))
ps[, trajectory := factor(trajectory, c("Null · entropy", "Sustained · entropy", "Phasic · value"))]
p <- ggplot(ps, aes(p, trajectory, color = test, shape = test)) +
  geom_vline(xintercept = .05, color = "#8999A3", linetype = 3) +
  geom_point(size = 3.7, position = position_dodge(width = .56)) +
  scale_x_log10(breaks = c(.00025, .001, .01, .05, .5, 1),
                labels = c("0.00025", "0.001", "0.01", "0.05", "0.5", "1")) +
  scale_color_manual(values = c(blue, orange, teal)) +
  labs(title = "Epoch tests detect different patterns",
       subtitle = "Minimum adjusted lag p-values for Max-|T| and TFCE, one p-value for omnibus",
       x = "Within-trajectory p-value (log scale)", y = NULL, color = NULL, shape = NULL)
save_fig("08-pr62-tests", p, 3.85)

b <- d$blur
z <- rbind(data.table(lag = b$lag, effect = b$estimate, curve = "Blurred input"),
           data.table(lag = b$lag, effect = b$estimate_multilag, curve = "After deblurring"),
           data.table(lag = b$lag, effect = d$truth, curve = "Known target"))
z[, marker_size := fifelse(curve == "Known target", 4.1, 2.6)]
p <- ggplot(z, aes(lag, effect, color = curve, linetype = curve)) +
  geom_hline(yintercept = 0, color = "#8999A3", linetype = 3) +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = curve, size = marker_size)) +
  scale_size_identity() +
  scale_shape_manual(values = c("After deblurring" = 16, "Blurred input" = 17, "Known target" = 1)) +
  scale_color_manual(values = c("After deblurring" = blue, "Blurred input" = orange, "Known target" = "#495761")) +
  scale_linetype_manual(values = c("After deblurring" = "solid", "Blurred input" = "longdash", "Known target" = "dotted")) +
  labs(title = "Deblurring changes the curve being summarized",
       subtitle = "Controlled coefficient-table example with a known resolution matrix",
       x = "Time from feedback (s)", y = "Effect amplitude", color = NULL, linetype = NULL, shape = NULL)
save_fig("09-pr62-deblur", p, 3.85)
