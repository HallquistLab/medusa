library(ggplot2)
d <- readRDS("step2_examples.rds")
list2env(d, environment())
dir.create("figures", showWarnings = FALSE)

C <- c("#2166AC", "#B35806", "#2A7F62", "#80539B", "#333333")
theme_set(
  theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 16),
      strip.text = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      plot.margin = margin(12, 18, 12, 12)
    )
)
savep <- function(name, p, w = 11, h = 5) {
  ggsave(paste0("figures/", name, ".png"), ggplot_gtable(ggplot_build(p)), width = w, height = h, dpi = 170, bg = "white")
}

# 1. High-frequency forward ambiguity.
s0 <- truth$run1[, 1]
fast <- 0.45 * (-1)^(1:n)
s1 <- s0 + fast
amb <- rbind(
  data.frame(time = rep(tt, 2), value = c(s0, s1), series = rep(c("Source A", "Source A + fast fluctuation"), each = n), panel = "Candidate latent activity"),
  data.frame(time = rep(tt, 2), value = c(H %*% s0, H %*% s1), series = rep(c("Source A", "Source A + fast fluctuation"), each = n), panel = "Convolved BOLD")
)
savep(
  "01-ambiguity",
  ggplot(subset(amb, time >= 32 & time <= 60), aes(time, value, color = series, linetype = series)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") +
    scale_color_manual(values = C[1:2]) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(x = "Acquisition time (s)", y = "Arbitrary model units", color = NULL, linetype = NULL, title = "BOLD convolution suppresses rapid neural fluctuations"),
  h = 5.6
)

# 2. Regularization sweep.
labs_l <- c(weak = "Weak: lambda = 0.0001", moderate = "Moderate: lambda = 0.03", strong = "Strong: lambda = 3")
sw <- do.call(rbind, lapply(names(fixed), function(k) {
  rbind(
    data.frame(time = tt, value = truth$run1[, 1], series = "Generating source", panel = unname(labs_l[k])),
    data.frame(time = tt, value = fixed[[k]]$activity$run1[, 1], series = "Estimated activity", panel = unname(labs_l[k]))
  )
}))
sw$panel <- factor(sw$panel, levels = labs_l)
savep(
  "02-lambda",
  ggplot(subset(sw, time >= 28 & time <= 103), aes(time, value, color = series, linetype = series)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~panel, ncol = 1) +
    scale_color_manual(values = c("Estimated activity" = C[1], "Generating source" = C[5])) +
    scale_linetype_manual(values = c("Estimated activity" = "solid", "Generating source" = "dashed")) +
    labs(x = "Acquisition time (s)", y = "Activity in input model units", color = NULL, linetype = NULL, title = "Regularization trades noise amplification against smoothing and shrinkage"),
  h = 7.4
)

# 3. Distinct penalty geometries.
pp <- do.call(rbind, lapply(c("ridge", "diff1", "diff2", "matern1", "sobolev"), function(k) {
  r <- penalties[[k]]
  panel <- sprintf("%s: lambda = %.3g", k, r$settings$lambda)
  rbind(
    data.frame(time = tt, value = truth$run1[, 1], series = "Generating source", panel = panel),
    data.frame(time = tt, value = r$activity$run1[, 1], series = "Estimated activity", panel = panel)
  )
}))
pp$panel <- factor(pp$panel, levels = unique(pp$panel))
savep(
  "03-penalties",
  ggplot(subset(pp, time >= 30 & time <= 104), aes(time, value, color = series, linetype = series)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~panel, ncol = 1) +
    scale_color_manual(values = c("Estimated activity" = C[1], "Generating source" = C[5])) +
    scale_linetype_manual(values = c("Estimated activity" = "solid", "Generating source" = "dashed")) +
    labs(x = "Acquisition time (s)", y = "Activity in input model units", color = NULL, linetype = NULL, title = "Different quadratic penalties encode distinct smoothness priors"),
  h = 10
)

# 4. Raw GCV vs adaptive task-weighted criterion.
t_df <- auto$settings$lambda_tuning
cv <- rbind(
  data.frame(lambda = t_df$lambda, value = t_df$gcv / min(t_df$gcv), series = "Raw GCV / its minimum", panel = "Selection criteria (unitless)"),
  data.frame(lambda = t_df$lambda, value = t_df$gcv_adaptive / min(t_df$gcv_adaptive), series = "Adaptive score / its minimum", panel = "Selection criteria (unitless)"),
  data.frame(lambda = t_df$lambda, value = t_df$var_level / min(t_df$var_level), series = "Level variance / its minimum", panel = "Gain-restored variance proxies"),
  data.frame(lambda = t_df$lambda, value = t_df$var_contrast / min(t_df$var_contrast), series = "Contrast variance / its minimum", panel = "Gain-restored variance proxies")
)
cv$panel <- factor(cv$panel, levels = c("Selection criteria (unitless)", "Gain-restored variance proxies"))
savep(
  "04-tuning",
  ggplot(cv, aes(lambda, value, color = series, linetype = series)) +
    geom_vline(xintercept = auto$settings$lambda, linetype = 3, color = "grey40") +
    geom_line(linewidth = 0.95) +
    geom_point(size = 1.8) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") +
    scale_x_log10() +
    scale_color_manual(values = C[1:4]) +
    scale_linetype_manual(values = c("solid", "dashed", "dotdash", "longdash")) +
    labs(x = "Lambda (log scale); vertical line is selected lambda", y = "Relative criterion / proxy", color = NULL, linetype = NULL, title = "Task-weighted GCV balances fit, gain-restored variance, and modulation transfer") +
    guides(color = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)),
  h = 6.7
)

# 5. Conditional impulse responses and gain-restoration tradeoff.
irs <- do.call(rbind, lapply(c("weak", "moderate", "strong"), function(k) {
  x <- impulse_curves[[k]]
  g <- x[91]
  rbind(
    data.frame(lag = tt - 90, value = x, mode = "Raw recovery", method = unname(labs_l[k])),
    data.frame(lag = tt - 90, value = x / g, mode = "Divided by its diagonal gain", method = unname(labs_l[k]))
  )
}))
irs$method <- factor(irs$method, levels = labs_l)
irs$mode <- factor(irs$mode, levels = c("Raw recovery", "Divided by its diagonal gain"))
savep(
  "05-resolution",
  ggplot(subset(irs, abs(lag) <= 12), aes(lag, value, color = method, linetype = method)) +
    geom_hline(yintercept = 0, color = "grey65") +
    geom_vline(xintercept = 0, color = "grey75") +
    geom_line(linewidth = 0.9) +
    facet_wrap(~mode, ncol = 1, scales = "free_y") +
    scale_color_manual(values = C[1:3]) +
    scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
    labs(x = "Time relative to injected impulse (s)", y = "Recovered response / injected amplitude", color = NULL, linetype = NULL, title = "Scalar gain division restores peak height but leaves temporal spread intact"),
  h = 5.7
)

# 6. Event-resolution matrix.
mat <- diag_fit$matrices[[1]]
ne <- nrow(mat$recovered)
off <- mat$normalized_excess
diag(off) <- NA_real_
leak <- data.frame(receiver = rep(1:ne, ne), source = rep(1:ne, each = ne), value = as.vector(off))
mx <- max(abs(off), na.rm = TRUE)
savep(
  "06-event-leakage",
  ggplot(leak, aes(source, receiver, fill = value)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradient2(low = "#B35806", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-mx, mx), breaks = c(-0.1, 0, 0.1), labels = c("-0.10", "0", "0.10"), na.value = "grey85", guide = guide_colourbar(barwidth = grid::unit(5, "cm"), barheight = grid::unit(0.35, "cm"))) +
    scale_x_continuous(breaks = 1:15) +
    scale_y_reverse(breaks = 1:15) +
    coord_equal() +
    labs(x = "Injected event (column)", y = "Receiving event readout (row)", fill = "Excess response /\nideal self-response", title = "Raw point readouts retain directed cross-talk that Step 3 de-mixes"),
  w = 9, h = 6
)

# 7. Nuisance and censoring choices.
nu <- rbind(
  data.frame(time = tt, value = truth$run1[, 1], series = "Generating source", panel = "Nuisance handling"),
  data.frame(time = tt, value = auto$activity$run1[, 1], series = "Joint nuisance", panel = "Nuisance handling"),
  data.frame(time = tt, value = ignore$activity$run1[, 1], series = "Ignore supplied nuisance", panel = "Nuisance handling"),
  data.frame(time = tt, value = truth$run1[, 1], series = "Generating source", panel = "Artifact and observation mask"),
  data.frame(time = tt, value = masked$activity$run1[, 1], series = "Exclude artifact rows", panel = "Artifact and observation mask"),
  data.frame(time = tt, value = allrows$activity$run1[, 1], series = "Use artifact rows (unweighted)", panel = "Artifact and observation mask")
)
nu$panel <- factor(nu$panel, levels = c("Nuisance handling", "Artifact and observation mask"))
shade <- data.frame(panel = "Artifact and observation mask", xmin = 70, xmax = 80, ymin = -Inf, ymax = Inf)
shade$panel <- factor(shade$panel, levels = levels(nu$panel))
savep(
  "07-mask-nuisance",
  ggplot(subset(nu, time >= 50 & time <= 100), aes(time, value, color = series, linetype = series)) +
    geom_rect(data = shade, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey85", inherit.aes = FALSE) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c("Exclude artifact rows" = C[1], "Generating source" = C[5], "Ignore supplied nuisance" = C[2], "Joint nuisance" = C[1], "Use artifact rows (unweighted)" = C[2])) +
    scale_linetype_manual(values = c("Exclude artifact rows" = "solid", "Generating source" = "dashed", "Ignore supplied nuisance" = "dotdash", "Joint nuisance" = "solid", "Use artifact rows (unweighted)" = "dotdash")) +
    labs(x = "Original acquisition time (s)", y = "Activity in input model units", color = NULL, linetype = NULL, title = "Nuisance projection and observation masks operate on the continuous acquisition grid") +
    guides(color = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)),
  h = 6.0
)

# 8. Scaling divisors.
sc <- rbind(
  data.frame(voxel = ids, sd = scales$zscore$scaling$run1$scale, kind = "Whole-BOLD SD (zscore)"),
  data.frame(voxel = ids, sd = scales$noise_sd$scaling$run1$scale, kind = "Event-informed GLS noise SD"),
  data.frame(voxel = ids, sd = sigma, kind = "Generating observation-noise SD")
)
savep(
  "08-scaling",
  ggplot(sc, aes(voxel, sd, color = kind, shape = kind)) +
    geom_point(size = 3, position = position_dodge(width = 0.3)) +
    scale_color_manual(values = c(C[1], C[5], C[2])) +
    scale_shape_manual(values = c(17, 4, 16)) +
    labs(x = "Voxel", y = "BOLD divisor in input signal units", color = NULL, shape = NULL, title = "Event-informed GLS scaling isolates observation noise from evoked task variance"),
  h = 4.5
)

# 9. Boundary conditions and HRF mismatch.
sens <- rbind(
  data.frame(time = tt, value = truth$run1[, 1], series = "Generating source", panel = "Run-end assumptions"),
  data.frame(time = tt, value = auto$activity$run1[, 1], series = "No appended BOLD observations", panel = "Run-end assumptions"),
  data.frame(time = tt, value = padded$activity$run1[, 1], series = "Append 33 zero BOLD observations", panel = "Run-end assumptions"),
  data.frame(time = tt, value = truth$run1[, 1], series = "Generating source", panel = "Generating HRF delay 6 s"),
  data.frame(time = tt, value = auto$activity$run1[, 1], series = "Correct assumed HRF", panel = "Generating HRF delay 6 s"),
  data.frame(time = tt, value = badhrf$activity$run1[, 1], series = "Assumed HRF delay 8 s", panel = "Generating HRF delay 6 s")
)
sens <- subset(sens, (panel == "Run-end assumptions" & time >= 145) | (panel == "Generating HRF delay 6 s" & time >= 30 & time <= 55))
sens$panel <- factor(sens$panel, levels = c("Run-end assumptions", "Generating HRF delay 6 s"))
savep(
  "09-sensitivity",
  ggplot(sens, aes(time, value, color = series, linetype = series)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~panel, ncol = 1, scales = "free_x") +
    scale_color_manual(values = c("Append 33 zero BOLD observations" = C[2], "Assumed HRF delay 8 s" = C[2], "Correct assumed HRF" = C[1], "Generating source" = C[5], "No appended BOLD observations" = C[1])) +
    scale_linetype_manual(values = c("Append 33 zero BOLD observations" = "dotdash", "Assumed HRF delay 8 s" = "dotdash", "Correct assumed HRF" = "solid", "Generating source" = "dashed", "No appended BOLD observations" = "solid")) +
    labs(x = "Original acquisition time (s)", y = "Activity in input model units", color = NULL, linetype = NULL, title = "Boundary padding and HRF delay mismatch alter inferred neural timing") +
    guides(color = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)),
  h = 6.0
)

# 10. Event-informed continuous GLS vs blind deconvolution.
bve <- blind_vs_event
traj <- rbind(
  data.frame(time = bve$time, value = bve$truth, series = "Generating neural source", panel = "A. Heavy neural burst (65 s): event-informed GLS vs. blind methods"),
  data.frame(time = bve$time, value = bve$event_auto, series = "Event-informed GLS (rwls = 'auto')", panel = "A. Heavy neural burst (65 s): event-informed GLS vs. blind methods"),
  data.frame(time = bve$time, value = bve$blind_lin, series = "Blind linear GCV (zscore)", panel = "A. Heavy neural burst (65 s): event-informed GLS vs. blind methods"),
  data.frame(time = bve$time, value = bve$blind_bush, series = "Blind nonlinear sigmoid (Bush 2011)", panel = "A. Heavy neural burst (65 s): event-informed GLS vs. blind methods")
)
glitch_df <- rbind(
  data.frame(time = bve$time, value = bve$truth, series = "Generating neural source", panel = "B. Non-task scanner glitch (124 s): adaptive rwls = 'auto' vs. unweighted"),
  data.frame(time = bve$time, value = bve$event_auto, series = "Event-informed GLS (rwls = 'auto')", panel = "B. Non-task scanner glitch (124 s): adaptive rwls = 'auto' vs. unweighted"),
  data.frame(time = bve$time, value = bve$event_norwls, series = "Event-informed GLS (rwls = 'none')", panel = "B. Non-task scanner glitch (124 s): adaptive rwls = 'auto' vs. unweighted"),
  data.frame(time = bve$time, value = bve$frame_weights, series = "Learned frame weight w_t", panel = "C. Task-envelope-guarded frame weights w_t (protects 65 s burst, downweights 124 s glitch)")
)
p10_df <- rbind(
  subset(traj, time >= 55 & time <= 145),
  subset(glitch_df, time >= 55 & time <= 145)
)
p10_df$panel <- factor(p10_df$panel, levels = c(
  "A. Heavy neural burst (65 s): event-informed GLS vs. blind methods",
  "B. Non-task scanner glitch (124 s): adaptive rwls = 'auto' vs. unweighted",
  "C. Task-envelope-guarded frame weights w_t (protects 65 s burst, downweights 124 s glitch)"
))
savep(
  "10-event-informed-vs-blind",
  ggplot(p10_df, aes(time, value, color = series, linetype = series)) +
    geom_vline(xintercept = c(bve$burst_time, bve$glitch_time), linetype = 3, color = "grey55") +
    geom_line(linewidth = 0.9) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c(
      "Generating neural source" = C[5],
      "Event-informed GLS (rwls = 'auto')" = C[1],
      "Event-informed GLS (rwls = 'none')" = C[2],
      "Blind linear GCV (zscore)" = C[3],
      "Blind nonlinear sigmoid (Bush 2011)" = C[4],
      "Learned frame weight w_t" = "#167C74"
    )) +
    scale_linetype_manual(values = c(
      "Generating neural source" = "dashed",
      "Event-informed GLS (rwls = 'auto')" = "solid",
      "Event-informed GLS (rwls = 'none')" = "dotdash",
      "Blind linear GCV (zscore)" = "longdash",
      "Blind nonlinear sigmoid (Bush 2011)" = "twodash",
      "Learned frame weight w_t" = "solid"
    )) +
    labs(
      x = "Acquisition time (s)",
      y = "Recovered activity / frame weight",
      color = NULL,
      linetype = NULL,
      title = "Event-informed continuous GLS preserves heavy neural peaks while suppressing non-task glitches"
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)),
  w = 11, h = 8.4
)
