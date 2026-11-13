library(ggplot2)
d <- readRDS("step4_examples.rds")
blue <- "#245D80"; teal <- "#167C74"; orange <- "#B75126"; gray <- "#637282"
theme_set(theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 17, color = "#172937"),
        plot.subtitle = element_text(size = 12, color = "#405362", margin = margin(b = 10)),
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom", plot.margin = margin(12, 17, 12, 12)))
save_fig <- function(name, plot, height) {
  ggsave(file.path("figures", paste0(name, ".png")), plot,
         width = 8.4, height = height, units = "in", dpi = 180, bg = "white")
}

w <- d$weights
w$run <- factor(w$run, c("run1", "run2"), c("All four voxels", "Voxel 1 unavailable"))
p <- ggplot(w, aes(voxel, realized, fill = voxel)) +
  geom_col(width = .64) +
  geom_text(aes(label = scales::percent(realized, accuracy = 1)),
            vjust = -.48, size = 4.1, color = "#233B4D") +
  facet_wrap(~run) +
  scale_fill_manual(values = c(teal, blue, "#6793AC", "#A2B2BF"), guide = "none") +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, .62),
                     breaks = c(0, .2, .4, .6)) +
  labs(title = "When a voxel drops out, the parcel changes",
       subtitle = "Actual weights in the simulated four-voxel parcel",
       x = NULL, y = "Contribution to weighted average")
save_fig("01-weights", p, 3.9)

m <- d$metrics
m$parcel <- factor(m$parcel, c("Negative", "Independent", "Shared", "Nearly common"))
m$run <- factor(m$run, c("run1", "run2"), c("Run 1", "Run 2"))
method_names <- c(oracle = "True parcel noise", medusa = "MEDuSA (task-free innovations)", naive_d2 = "Naive activity Δ² proxy")
ml <- do.call(rbind, lapply(names(method_names), function(k) {
  data.frame(run = m$run, parcel = m$parcel, method = unname(method_names[k]), sd = m[[k]])
}))
ml$method <- factor(ml$method, unname(method_names))
p <- ggplot(m, aes(parcel)) +
  geom_linerange(aes(ymin = independence, ymax = ceiling), linewidth = 4,
                 color = "#CBD6DC", lineend = "round") +
  geom_point(data = ml, aes(y = sd, shape = method, color = method),
             position = position_dodge(width = .48), size = 3.3) +
  facet_wrap(~run, ncol = 1) +
  scale_color_manual(values = c("#1F3444", teal, orange)) +
  scale_shape_manual(values = c(16, 15, 17)) +
  scale_y_continuous(breaks = c(.005, .01, .015, .02), limits = c(.004, .020),
                     labels = scales::label_number(accuracy = .001)) +
  labs(title = "Task-free noise innovations isolate true spatial noise correlation",
       subtitle = "Gray spans mark independence to common-noise bounds; naive Δ² inflates under shared task signal",
       x = "True spatial noise regime", y = "Effective BOLD noise SD", color = NULL, shape = NULL)
save_fig("02-noise", p, 5.35)

a <- subset(d$mixed, lag >= -10 & lag <= 14)
am <- rbind(data.frame(lag = a$lag, activity = a$actual, method = "MEDuSA multi-stratum parcel operator"),
            data.frame(lag = a$lag, activity = a$representative, method = "Single averaged-HRF approximation"))
p <- ggplot(am, aes(lag, activity, color = method, linetype = method)) +
  geom_hline(yintercept = 0, color = "#A5B5BF", linewidth = .3) +
  geom_vline(xintercept = 0, color = "#A5B5BF", linetype = 3) +
  geom_line(linewidth = 1.15) + geom_point(size = 1.45) +
  scale_color_manual(values = c(blue, orange)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(title = "Multi-stratum parcel operators avoid single-HRF approximation error",
       subtitle = "Same unit impulse in two voxels with HRF peaks at 4 s and 8 s",
       x = "Time from impulse (s)", y = "Reconstructed activity", color = NULL, linetype = NULL)
save_fig("03-mixed-hrf", p, 4.1)

q <- subset(d$panel, parcel == "Nearly common")
r <- subset(d$panel_resid, parcel == "Nearly common")
cols <- c("se_demixed", "precision_weight")
q2 <- rbind(data.frame(q[cols], variance = "Marginal"), data.frame(r[cols], variance = "Residual"))
p <- ggplot(q2, aes(se_demixed, precision_weight, color = variance)) +
  geom_hline(yintercept = 1, color = "#98ABB5", linetype = 3) +
  geom_point(alpha = .74, size = 2.1) +
  scale_color_manual(values = c(blue, orange)) +
  labs(title = "The working trial variance changes the weights",
       subtitle = "Same measurements and SEs, two estimates of latent trial variance",
       x = "Trial measurement SE", y = "Relative precision weight", color = NULL)
save_fig("04-precision", p, 4.0)

q <- subset(q, run == "run1")
mu <- mean(q$decon_demixed)
p <- ggplot(q, aes(decon_demixed, decon_eb)) +
  geom_abline(color = "#98ABB5", linetype = 3) +
  geom_hline(yintercept = mu, color = "#98ABB5", linewidth = .35) +
  geom_point(aes(color = reliability_trial), size = 3) +
  scale_color_gradient(low = orange, high = blue, breaks = scales::breaks_pretty(n = 3),
                       guide = guide_colourbar(barwidth = grid::unit(6, "cm"))) +
  coord_equal() +
  labs(title = "Shrinkage pulls estimates toward the run mean",
       subtitle = "Nearly common noise parcel, Run 1",
       x = "Demixed trial amplitude", y = "Shrunk trial amplitude", color = "Reliability")
save_fig("05-shrinkage", p, 4.6)

cc <- d$transform_check
long <- function(x, label) {
  diag(x) <- NA_real_
  data.frame(trial_i = as.vector(row(x)), trial_j = as.vector(col(x)),
             cor = as.vector(x), transform = label)
}
cor_data <- rbind(long(cc$cor_before, "Before"), long(cc$cor_after, "After"))
cor_data$transform <- factor(cor_data$transform, c("Before", "After"))
mx <- max(abs(cor_data$cor), na.rm = TRUE)
p <- ggplot(cor_data, aes(trial_j, trial_i, fill = cor)) +
  geom_raster() + facet_wrap(~transform) + coord_equal() +
  scale_fill_gradient2(low = blue, mid = "white", high = orange, midpoint = 0,
    limits = c(-mx, mx), breaks = c(-.02, 0, .02), na.value = "#DEE6EB",
    guide = guide_colourbar(barwidth = grid::unit(6, "cm"))) +
  scale_x_continuous(breaks = c(1, 9, 18, 27)) +
  scale_y_continuous(breaks = c(1, 9, 18, 27)) +
  labs(title = "Regularized ZCA changes the working correlation",
       subtitle = "Nearly common noise parcel, Run 1. Diagonals omitted.",
       x = "Trial index", y = "Trial index", fill = "Correlation")
save_fig("06-zca", p, 4.95)
