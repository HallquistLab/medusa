library(medusa)
library(ggplot2)
library(patchwork)

set.seed(20261121L)

ink    <- "#172937"
blue   <- "#245D80"
teal   <- "#167C74"
orange <- "#B75126"
gray   <- "#637282"

theme_step1 <- theme_minimal(base_size = 11.5) +
  theme(
    plot.title       = element_text(face = "bold", size = 11.8, color = ink),
    plot.subtitle    = element_text(size = 9.0, color = "#405362", margin = margin(b = 6)),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    plot.margin      = margin(8, 12, 8, 10)
  )

# -------------------------------------------------------------------------
# Panel A: HRF Low-Pass Filter |H(f)| vs. Sobolev Inverse Gain |D(f)|
# -------------------------------------------------------------------------
tr <- 0.8
cand <- hrf_library(tr = tr, preset = "medusa")
h_can <- as.numeric(cand[, 1L])
N_fft <- 1024L
h_pad <- c(h_can, rep(0, N_fft - length(h_can)))
H_fft <- Mod(fft(h_pad))
freqs <- (0:(N_fft / 2L)) / (N_fft * tr)
H_mag <- H_fft[seq_along(freqs)]

s2 <- 4 * sin(pi * freqs * tr)^2
P_sob <- 0.05 + 0.45 * s2 + 0.50 * s2^2

gain_lam001 <- H_mag / (H_mag^2 + 0.01 * P_sob)
gain_lam010 <- H_mag / (H_mag^2 + 0.10 * P_sob)

keep_f <- freqs >= 0.005 & freqs <= 0.40
freq_df <- rbind(
  data.frame(freq = freqs[keep_f], gain = H_mag[keep_f],       curve = "Forward HRF |H(f)|"),
  data.frame(freq = freqs[keep_f], gain = gain_lam001[keep_f], curve = "Inverse gain |D(f)| (lambda = 0.01)"),
  data.frame(freq = freqs[keep_f], gain = gain_lam010[keep_f], curve = "Inverse gain |D(f)| (lambda = 0.10)")
)
freq_df$curve <- factor(freq_df$curve, levels = unique(freq_df$curve))

ymax_a <- max(gain_lam001[keep_f])
p_freq <- ggplot(freq_df, aes(x = freq, y = gain, color = curve, linetype = curve)) +
  annotate("rect", xmin = 0.02, xmax = 0.10, ymin = -Inf, ymax = Inf, fill = teal, alpha = 0.09) +
  annotate("rect", xmin = 0.14, xmax = 0.30, ymin = -Inf, ymax = Inf, fill = orange, alpha = 0.10) +
  geom_hline(yintercept = 1, color = "#95A5B0", linetype = "dotted") +
  geom_line(linewidth = 1.05) +
  annotate("text", x = 0.06, y = ymax_a * 1.13,
           label = "Evoked\ntask band", size = 2.8, color = teal, fontface = "bold", lineheight = 1.05) +
  annotate("text", x = 0.22, y = ymax_a * 1.13,
           label = "Respiratory & aliased\ncardiac (0.14-0.30 Hz)", size = 2.8, color = orange, fontface = "bold", lineheight = 1.05) +
  scale_color_manual(values = c(gray, orange, blue)) +
  scale_linetype_manual(values = c("dashed", "solid", "solid")) +
  scale_x_continuous(breaks = seq(0, 0.40, by = 0.10)) +
  scale_y_continuous(limits = c(0, ymax_a * 1.24), expand = expansion(mult = c(0.02, 0.02))) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    title    = "A. Smoothing limits inverse-filter gain",
    subtitle = "Frequency response for the simulated HRF on a regular grid",
    x        = "Frequency (Hz; TR = 0.8 s)",
    y        = "Amplitude gain",
    color    = NULL, linetype = NULL
  ) +
  theme_step1 +
  theme(legend.text = element_text(size = 8.8))

# -------------------------------------------------------------------------
# Panel B: Live medusa simulation comparing noise_pcs = "auto" vs. 0L
# -------------------------------------------------------------------------
n_t <- 260L
tt  <- (0:(n_t - 1L)) * tr
n_parcels <- 14L   # 2 active + 12 background (>= 8 bg parcels, >= 24 bg voxels)
vox_per_p <- 3L
V <- n_parcels * vox_per_p
parcel_ids <- setNames(rep(sprintf("P%02d", seq_len(n_parcels)), each = vox_per_p),
                       sprintf("v%02d", seq_len(V)))

# Two spatially orthogonal pervasive physiological modes across cortex
load_p1 <- seq(0.055, 0.080, length.out = V)
load_p2 <- 0.065 * cos(seq(0.1, 1.8 * pi, length.out = V))

make_run_data <- function(run_seed, with_physio = TRUE) {
  set.seed(run_seed)
  onsets <- sort(seq(16, 188, length.out = 13L) + runif(13L, -1.5, 1.5))
  amps   <- 1.1 + 0.35 * rnorm(length(onsets))
  ev <- data.frame(onset = onsets, duration = 0, condition = "event", stringsAsFactors = FALSE)

  a_true <- numeric(n_t)
  for (j in seq_along(onsets)) {
    w_imp <- pmax(0, 1 - abs(tt / tr - onsets[j] / tr))
    a_true <- a_true + amps[j] * w_imp
  }
  bold_sig <- as.numeric(stats::convolve(a_true, rev(h_can), type = "open"))[seq_len(n_t)]

  # Narrowband respiratory (0.17 Hz) + aliased cardiac (0.24 Hz) waves inside inverse-gain band
  p1 <- as.numeric(scale(sin(2 * pi * 0.17 * tt + runif(1, 0, 2 * pi)) + 0.25 * cos(2 * pi * 0.18 * tt)))
  p2 <- as.numeric(scale(sin(2 * pi * 0.24 * tt + runif(1, 0, 2 * pi)) + 0.25 * sin(2 * pi * 0.25 * tt)))

  Y <- matrix(0, n_t, V)
  colnames(Y) <- names(parcel_ids)
  for (v in seq_len(V)) {
    p_num <- ceiling(v / vox_per_p)
    eps <- as.numeric(stats::arima.sim(list(ar = 0.40), n = n_t, sd = 0.022 * sqrt(1 - 0.40^2)))
    phys <- if (with_physio) load_p1[v] * p1 + load_p2[v] * p2 else 0
    sig  <- if (p_num <= 2L) 0.70 * bold_sig else 0
    Y[, v] <- 100 + sig + phys + eps
  }
  list(Y = Y, ev = ev, a_true = 0.70 * a_true)
}

r1_p <- make_run_data(101L, with_physio = TRUE)
r2_p <- make_run_data(102L, with_physio = TRUE)
r3_p <- make_run_data(103L, with_physio = TRUE)
bold_p <- list(run1 = r1_p$Y, run2 = r2_p$Y, run3 = r3_p$Y)
ev_p   <- list(run1 = r1_p$ev, run2 = r2_p$ev, run3 = r3_p$ev)

r1_c <- make_run_data(101L, with_physio = FALSE)
r2_c <- make_run_data(102L, with_physio = FALSE)
r3_c <- make_run_data(103L, with_physio = FALSE)
bold_c <- list(run1 = r1_c$Y, run2 = r2_c$Y, run3 = r3_c$Y)
ev_c   <- list(run1 = r1_c$ev, run2 = r2_c$ev, run3 = r3_c$ev)

hf_auto <- estimate_crossrun_hrf(bold_p, ev_p, tr = tr, kernel_grid = cand,
                                 spatial_groups = parcel_ids, pooling = "parcel",
                                 selection_mode = "soft_average", noise_pcs = "auto")
dc_auto <- deconvolve_crossrun_hrf(bold_p, hf_auto, lambda = 0.015, scaling = "noise_sd",
                                   penalty = "sobolev", nuisance_mode = "preproject",
                                   noise_pcs = "auto", rwls = "none")

hf_zero <- estimate_crossrun_hrf(bold_p, ev_p, tr = tr, kernel_grid = cand,
                                 spatial_groups = parcel_ids, pooling = "parcel",
                                 selection_mode = "soft_average", noise_pcs = 0L)
dc_zero <- deconvolve_crossrun_hrf(bold_p, hf_zero, lambda = 0.015, scaling = "noise_sd",
                                   penalty = "sobolev", nuisance_mode = "preproject",
                                   noise_pcs = 0L, rwls = "none")

hf_clean <- estimate_crossrun_hrf(bold_c, ev_c, tr = tr, kernel_grid = cand,
                                  spatial_groups = parcel_ids, pooling = "parcel",
                                  selection_mode = "soft_average", noise_pcs = "auto")
dc_clean <- deconvolve_crossrun_hrf(bold_c, hf_clean, lambda = 0.015, scaling = "noise_sd",
                                    penalty = "sobolev", nuisance_mode = "preproject",
                                    noise_pcs = "auto", rwls = "none")

get_k_used <- function(dc_obj) {
  if (!is.null(dc_obj$settings$noise_pcs_used)) {
    u <- dc_obj$settings$noise_pcs_used
    if (is.list(u)) return(as.integer(u[[1L]][1L]))
    return(as.integer(u[1L]))
  }
  0L
}
k_phys  <- get_k_used(dc_auto)
k_clean <- get_k_used(dc_clean)

sc_auto <- dc_auto$scaling$run1$scale[1:3]
sc_zero <- dc_zero$scaling$run1$scale[1:3]
a_auto <- rowMeans(sweep(dc_auto$activity$run1[, 1:3, drop = FALSE], 2L, sc_auto, "*"))
a_zero <- rowMeans(sweep(dc_zero$activity$run1[, 1:3, drop = FALSE], 2L, sc_zero, "*"))
a_true <- r1_p$a_true

cor_zero <- stats::cor(a_zero, a_true)
cor_auto <- stats::cor(a_auto, a_true)
far_mask <- vapply(tt, function(t0) min(abs(t0 - r1_p$ev$onset)) > 3.0, logical(1))
sd_bg_zero <- stats::sd(a_zero[far_mask])
sd_bg_auto <- stats::sd(a_auto[far_mask])

win <- tt >= 40 & tt <= 106
on_win <- r1_p$ev$onset[r1_p$ev$onset >= 40 & r1_p$ev$onset <= 106]
trace_df <- rbind(
  data.frame(time = tt[win], val = a_zero[win], method = sprintf("noise_pcs = 0 (r = %.2f, off-event SD = %.3f)", cor_zero, sd_bg_zero)),
  data.frame(time = tt[win], val = a_auto[win], method = sprintf("noise_pcs = 'auto' (K* = %d, r = %.2f, off-event SD = %.3f)", k_phys, cor_auto, sd_bg_auto))
)
trace_df$method <- factor(trace_df$method, levels = unique(trace_df$method))

ymax_b <- max(trace_df$val)
p_trace <- ggplot(trace_df, aes(x = time, y = val, color = method)) +
  geom_vline(xintercept = on_win, linetype = "dotted", color = "#7B8C98", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "#B8C4CC", linewidth = 0.35) +
  geom_line(linewidth = 0.95) +
  annotate("label", x = 73, y = ymax_b * 1.13,
           label = sprintf("AR(1) control: %d PCs selected", k_clean),
           size = 2.85, color = ink, fill = "#F4F7F9", linewidth = 0.25) +
  scale_color_manual(values = c(orange, teal)) +
  scale_y_continuous(limits = c(min(trace_df$val) * 1.08, ymax_b * 1.24)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    title    = "B. Live medusa deconvolution (Parcel 1, Run 1)",
    subtitle = "Compare the two traces between event onsets",
    x        = "Run time (s; dotted lines = task events)",
    y        = "Restored activity (BOLD units)",
    color    = NULL
  ) +
  theme_step1 +
  theme(legend.text = element_text(size = 8.8))

fig <- p_freq | p_trace

ggsave(
  "figures/12-noise-pool-pcs.png",
  fig,
  width = 10.8, height = 4.75, units = "in", dpi = 180, bg = "white"
)

saveRDS(
  list(
    k_phys = k_phys, k_clean = k_clean,
    cor_zero = cor_zero, cor_auto = cor_auto,
    sd_bg_zero = sd_bg_zero, sd_bg_auto = sd_bg_auto,
    max_gain_001 = max(gain_lam001[freqs >= 0.14 & freqs <= 0.30])
  ),
  "noise_pcs_demo_metrics.rds"
)

cat(sprintf(
  "Saved figures/12-noise-pool-pcs.png | K_phys=%s, K_clean=%s | cor: %.3f -> %.3f | bg SD: %.3f -> %.3f | max gain=%.2f\n",
  as.character(k_phys), as.character(k_clean), cor_zero, cor_auto, sd_bg_zero, sd_bg_auto,
  max(gain_lam001[freqs >= 0.14 & freqs <= 0.30])
))
