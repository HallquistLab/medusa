library(medusa)
library(ggplot2)
d <- readRDS("step1_examples.rds")
list2env(d, environment())
idx <- function(d_val, s_val) which(pars$peak_delay == d_val & pars$dispersion == s_val)
dir.create("figures", showWarnings = FALSE)

cols <- c("#2166AC", "#B35806", "#6A6A6A", "#2A7F62")
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
  g <- ggplot_gtable(ggplot_build(p))
  ggsave(paste0("figures/", name, ".png"), g, width = w, height = h, dpi = 170, bg = "white")
}
longcurve <- function(m, labels, tr_val = tr) {
  do.call(rbind, lapply(seq_len(ncol(m)), function(i) {
    data.frame(time = (0:(nrow(m) - 1L)) * tr_val, value = m[, i], curve = labels[i])
  }))
}

# 1. 3D physiological HRF manifold dimensions (positive-lobe normalized).
K_del <- hrf_library(tr, peak_delays = c(3.0, 6.0, 8.0), dispersions = 1.0, undershoot_ratios = 6, normalize = "positive_lobe")
p_del <- attr(K_del, "parameters")
K_disp <- hrf_library(tr, peak_delays = 6.0, dispersions = c(0.6, 1.0, 1.4), undershoot_ratios = 6, normalize = "positive_lobe")
p_disp <- attr(K_disp, "parameters")
K_und <- hrf_library(tr, peak_delays = 6.0, dispersions = 1.0, undershoot_ratios = c(3, 6, Inf), normalize = "positive_lobe")
p_und <- attr(K_und, "parameters")
lib <- rbind(
  transform(longcurve(K_del[, match(c(3.0, 6.0, 8.0), p_del$peak_delay)], c("3.0 s (fast)", "6.0 s (canonical)", "8.0 s (slow)")), panel = "Peak delay (disp = 1, ratio = 6)"),
  transform(longcurve(K_disp[, match(c(0.6, 1.0, 1.4), p_disp$dispersion)], c("0.6 (narrow)", "1.0 (standard)", "1.4 (broad)")), panel = "Dispersion (delay = 6 s, ratio = 6)"),
  transform(longcurve(K_und[, c(which(p_und$undershoot_ratio == 3), which(p_und$undershoot_ratio == 6), which(is.infinite(p_und$undershoot_ratio)))], c("Ratio 3 (strong dip)", "Ratio 6 (canonical)", "Inf (no undershoot)")), panel = "Undershoot (delay = 6 s, disp = 1)")
)
lib$panel <- factor(lib$panel, levels = unique(lib$panel))
savep(
  "01-library",
  ggplot(lib, aes(time, value, color = curve)) +
    geom_hline(yintercept = 0, color = "grey75") +
    geom_line(linewidth = 0.95) +
    facet_wrap(~panel, nrow = 1) +
    labs(
      x = "Time after neural impulse (s)",
      y = "Positive-lobe-normalized HRF",
      color = NULL,
      title = "The physiological HRF manifold spans delay, width, and undershoot depth"
    ) +
    scale_color_manual(values = rep(cols[c(1, 2, 4)], 3)) +
    guides(color = guide_legend(nrow = 3)),
  w = 11.5, h = 5.1
)

# 2. Forward model.
fr <- rbind(
  data.frame(time = events$run1$onset, value = amplitudes$run1, curve = "Neural amplitudes", panel = "Neural events (arbitrary amplitude)"),
  data.frame(time = rep(tt, 2), value = c(bold$run1[, 1], bold$run1[, 9]), curve = rep(c("Fast HRF region", "Slow HRF region"), each = n), panel = "Observed BOLD (arbitrary signal units)")
)
savep(
  "02-forward",
  ggplot(subset(fr, time <= 90), aes(time, value, color = curve)) +
    geom_line(data = function(x) subset(x, panel != "Neural events (arbitrary amplitude)"), linewidth = 0.7) +
    geom_segment(data = function(x) subset(x, panel == "Neural events (arbitrary amplitude)"), aes(xend = time, yend = 0), linewidth = 0.7) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") +
    scale_color_manual(values = c(cols[1], cols[3], cols[2])) +
    labs(x = "Acquisition time (s)", y = NULL, color = NULL, title = "The same neural events produce different BOLD timing across vascular regimes"),
  h = 5.4
)

# 3. Cross-fitting diagram.
f <- expand.grid(source = 1:3, target = 1:3)
f$role <- ifelse(f$source == f$target, "Own-run noise only", "HRF evidence")
f$label <- ifelse(f$source == f$target, "Noise refit here\nHRF evidence excluded", "Candidate scores\ncontribute to HRF")
savep(
  "03-folds",
  ggplot(f, aes(factor(source), factor(target, levels = 3:1), fill = role)) +
    geom_tile(color = "white", linewidth = 4) +
    geom_text(aes(label = label), size = 4.5, lineheight = 1.1) +
    scale_fill_manual(values = c("HRF evidence" = "#CEE2F0", "Own-run noise only" = "#F4DDC5")) +
    labs(x = "Run providing data", y = "Target run receiving HRF", fill = NULL, title = "Leave the target run out of its HRF decision") +
    theme(panel.grid = element_blank(), legend.position = "none"),
  w = 10, h = 5.5
)

# 4. Fold evidence; target 1 uses runs 2 and 3.
s <- soft$scores$run1
heat <- do.call(rbind, lapply(names(true_idx), function(g) {
  data.frame(pars, score = rowMeans(s[, groups == g, drop = FALSE]), region = g)
}))
truth <- data.frame(region = names(true_idx)[1:2], pars[true_idx[1:2], c("peak_delay", "dispersion")])
savep(
  "04-evidence",
  ggplot(heat, aes(peak_delay, factor(dispersion), fill = score)) +
    geom_tile(color = "white") +
    geom_point(data = truth, aes(peak_delay, factor(dispersion)), inherit.aes = FALSE, shape = 4, size = 4, stroke = 1.4) +
    facet_wrap(~region, nrow = 1) +
    scale_fill_gradient2(low = "#B35806", mid = "white", high = "#2166AC", midpoint = 0, limits = range(heat$score)) +
    scale_x_continuous(breaks = c(4, 5, 6, 7, 8)) +
    labs(x = "Library delay parameter (s)", y = "Dispersion", fill = "Mean adjusted\npartial R²", title = "Out-of-sample training runs recover regional HRF shape when task signal exists"),
  w = 11, h = 4.5
)

# 5. Signed example.
sc <- cbind(K[, idx(7.5, 1.1)], signed_condition$selected_kernels$run1[, 1], signed_trial$selected_kernels$run1[, 1])
sc <- sweep(sc, 2, apply(sc, 2, max), "/")
savep(
  "05-signed",
  ggplot(longcurve(sc, c("Generating HRF", "Condition amplitudes", "Free trial amplitudes")), aes(time, value, color = curve, linetype = curve)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c(cols[2], cols[1], cols[3])) +
    scale_linetype_manual(values = c("dashed", "solid", "dotdash")) +
    labs(x = "Time after neural impulse (s)", y = "HRF / its positive peak", color = NULL, linetype = NULL, title = "Zero-mean signed trial modulations cancel under a condition-mean model"),
  h = 4.5
)

# 6. Pooling stabilizes shape across heterogeneous-noise voxels.
peak <- function(m) apply(m, 2, which.max) * tr - tr
pk <- rbind(
  data.frame(sd = sd_vox, peak = peak(voxel$selected_kernels$run1), method = "Voxel hard selection", region = groups),
  data.frame(sd = sd_vox, peak = peak(hard$selected_kernels$run1), method = "Parcel hard selection", region = groups)
)
pk <- subset(pk, region != "Null")
pk$truth <- ifelse(pk$region == "Fast", peak(K[, true_idx["Fast"], drop = FALSE]), peak(K[, true_idx["Slow"], drop = FALSE]))
savep(
  "06-pooling",
  ggplot(pk, aes(sd, peak, color = method, shape = method)) +
    geom_hline(data = unique(pk[, c("region", "truth")]), aes(yintercept = truth), linetype = 2, color = "grey40") +
    geom_point(size = 3, position = position_dodge(width = 0.0018)) +
    facet_wrap(~region) +
    scale_color_manual(values = cols[1:2]) +
    labs(x = "Generating marginal noise SD", y = "Sampled HRF peak time (s)", color = NULL, shape = NULL, title = "Parcel evidence pooling stabilizes HRF shape without averaging raw BOLD"),
  h = 4.5
)

# 7. Actual composites and hard selections.
curves <- do.call(rbind, lapply(names(true_idx)[1:2], function(g) {
  v <- which(groups == g)[1]
  transform(longcurve(cbind(K[, true_idx[g]], hard$selected_kernels$run1[, v], soft$selected_kernels$run1[, v]), c("Generating HRF", "Hard library node", "Soft composite")), region = g)
}))
savep(
  "07-composites",
  ggplot(curves, aes(time, value, color = curve, linetype = curve)) +
    geom_line(linewidth = 1) +
    facet_wrap(~region) +
    scale_color_manual(values = c(cols[3], cols[2], cols[1])) +
    scale_linetype_manual(values = c("solid", "dashed", "dotdash")) +
    labs(x = "Time after neural impulse (s)", y = "Normalized HRF kernel", color = NULL, linetype = NULL, title = "Soft evidence averaging blends nearby library candidates smoothly"),
  h = 4.5
)

# 8. ACF using known simulation errors: post-selection AR(p) refit vs misspecified fixed HRF.
ar_cols <- intersect(c("rho", "psi1", "psi2", "psi3", "psi4"), names(soft$noise$run1))
aout <- do.call(rbind, lapply(c("Fast", "Null"), function(g) {
  v <- which(groups == g)
  ah_soft <- unlist(soft$noise$run1[v[1], ar_cols], use.names = FALSE)
  ah_can <- unlist(canonical_only$noise$run1[v[1], ar_cols], use.names = FALSE)
  e <- errors$run1[, v]
  rbind(
    data.frame(lag = 0:16, acf = rowMeans(sapply(1:8, function(i) as.numeric(acf(e[, i], lag.max = 16, plot = FALSE)$acf))), series = "Before whitening", region = g),
    data.frame(lag = 0:16, acf = rowMeans(sapply(1:8, function(i) as.numeric(acf(medusa:::hrf_whitening_matrix(1:n, ah_can) %*% e[, i], lag.max = 16, plot = FALSE)$acf))), series = "Fixed canonical HRF noise fit", region = g),
    data.frame(lag = 0:16, acf = rowMeans(sapply(1:8, function(i) as.numeric(acf(medusa:::hrf_whitening_matrix(1:n, ah_soft) %*% e[, i], lag.max = 16, plot = FALSE)$acf))), series = "Post-selection AR(p) refit (MEDuSA)", region = g)
  )
}))
aout$series <- factor(aout$series, levels = c("Before whitening", "Fixed canonical HRF noise fit", "Post-selection AR(p) refit (MEDuSA)"))
savep(
  "08-noise",
  ggplot(subset(aout, lag > 0), aes(lag * tr, acf, color = series, shape = series)) +
    geom_hline(yintercept = 0, color = "grey45") +
    geom_line(linewidth = 0.85) +
    geom_point(size = 2.2) +
    facet_wrap(~region) +
    scale_color_manual(values = c(cols[3], cols[2], cols[1])) +
    labs(x = "Lag (s)", y = "Mean error autocorrelation", color = NULL, shape = NULL, title = "Post-selection AR(p) refitting prevents HRF mismatch from distorting the whitener"),
  h = 4.7
)

# 9. Invariance of HRF weights after changing target BOLD.
wdat <- do.call(rbind, lapply(c("Original", "Perturbed target"), function(z) {
  fit <- if (z == "Original") soft else perturbed
  data.frame(candidate = 1:ncol(K), weight = fit$kernel_weights$run1[, 1], version = z)
}))
savep(
  "09-target-check",
  ggplot(wdat, aes(candidate, weight, color = version, shape = version)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 3) +
    scale_color_manual(values = cols[1:2]) +
    scale_shape_manual(values = c(1, 4)) +
    labs(x = "HRF library candidate index", y = "Weight for target run 1, voxel 1", color = NULL, shape = NULL, title = "Perturbing target-run BOLD leaves its cross-fitted HRF weights strictly invariant"),
  h = 4.2
)
write.csv(aout, "whitening_acf.csv", row.names = FALSE)
