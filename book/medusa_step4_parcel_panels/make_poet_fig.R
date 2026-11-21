library(medusa)
library(ggplot2)
library(patchwork)

set.seed(20261120L)

# Book palette
ink    <- "#172937"
blue   <- "#245D80"
teal   <- "#167C74"
orange <- "#B75126"
gray   <- "#637282"

theme_poet <- theme_minimal(base_size = 11.5) +
  theme(
    plot.title       = element_text(face = "bold", size = 12.2, color = ink),
    plot.subtitle    = element_text(size = 9.6, color = "#405362", margin = margin(b = 6)),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 9.5, color = ink),
    legend.position  = "bottom",
    plot.margin      = margin(8, 12, 8, 10)
  )

# -------------------------------------------------------------------------
# Part 1: Stationary Adaptive-POET Decomposition (P = 30 parcels, df = 180)
# -------------------------------------------------------------------------
P <- 30L
df_eff <- 180L

# 2 pervasive global physiological factors (global respiration + vascular drift)
f1 <- seq(0.42, 0.58, length.out = P)
f2 <- 0.34 * sin(seq(0.3, 2.2 * pi, length.out = P))
L_true <- cbind(f1, f2)
LL_true <- tcrossprod(L_true)

# 3 sparse local vascular neighbor clusters (adjacent parcels sharing draining veins)
S_true <- matrix(0, P, P)
local_pairs <- list(c(4L, 5L), c(14L, 15L), c(24L, 25L), c(24L, 26L), c(25L, 26L))
for (pr in local_pairs) {
  S_true[pr[1], pr[2]] <- 0.42
  S_true[pr[2], pr[1]] <- 0.42
}
R_true <- LL_true + S_true
diag(R_true) <- 1.0
eg_t <- eigen(0.5 * (R_true + t(R_true)), symmetric = TRUE)
R_true <- eg_t$vectors %*% (pmax(eg_t$values, 0.05) * t(eg_t$vectors))
d_t <- sqrt(diag(R_true))
R_true <- R_true / outer(d_t, d_t)
diag(R_true) <- 1.0

# Simulate df_eff Gaussian noise innovation draws with correlation R_true
chol_R <- chol(R_true)
E_sim <- matrix(rnorm(df_eff * P), nrow = df_eff, ncol = P) %*% chol_R
R_raw <- stats::cor(E_sim)

# Fit MEDuSA's Adaptive-POET estimator directly
poet_fit <- medusa:::poet_regularize_cor(R_raw, df_u = df_eff, k_max = 12L, poet_mode = "adaptive_logp")
R_poet   <- poet_fit$R
LL_hat   <- if (poet_fit$K > 0L) tcrossprod(poet_fit$L) else matrix(0, P, P)
S_hat    <- poet_fit$S_local

# Marchenko-Pastur + Tracy-Widom threshold used inside poet_shrink_spikes()
gamma_asp <- P / df_eff
tw_margin <- max(0.12, 2.2 * ((1 + sqrt(gamma_asp))^(4 / 3)) * (max(gamma_asp, 1e-4)^(1 / 6)) / (df_eff^(2 / 3)))
mp_edge   <- (1 + sqrt(gamma_asp))^2 + tw_margin
eig_raw   <- eigen(0.5 * (R_raw + t(R_raw)), symmetric = TRUE, only.values = TRUE)$values

off_mask  <- row(R_true) != col(R_true)
frob_raw  <- sqrt(mean((R_raw[off_mask]  - R_true[off_mask])^2))
frob_poet <- sqrt(mean((R_poet[off_mask] - R_true[off_mask])^2))
local_idx <- do.call(rbind, local_pairs)
mean_local_true <- mean(R_true[local_idx] - LL_true[local_idx])
mean_local_poet <- mean(S_hat[local_idx])

# Panel A: Eigenvalue spectrum vs. MP + TW bulk edge
eig_df <- data.frame(
  idx   = seq_len(P),
  val   = eig_raw,
  class = ifelse(seq_len(P) <= poet_fit$K,
                 sprintf("Pervasive physiological spike (K = %d)", poet_fit$K),
                 "Random-matrix bulk & local pairs")
)
eig_df$class <- factor(eig_df$class, levels = c(
  sprintf("Pervasive physiological spike (K = %d)", poet_fit$K),
  "Random-matrix bulk & local pairs"
))

p_scree <- ggplot(eig_df, aes(x = idx, y = val, color = class)) +
  geom_hline(yintercept = mp_edge, linetype = "dashed", color = orange, linewidth = 0.8) +
  geom_segment(aes(xend = idx, y = 0, yend = val), linewidth = 0.75, alpha = 0.75) +
  geom_point(size = 2.3) +
  annotate("label", x = 18.5, y = mp_edge + 1.15,
           label = sprintf("MP + TW edge = %.2f", mp_edge),
           size = 3.0, color = orange, fill = "white", linewidth = 0.25) +
  scale_color_manual(values = c(orange, blue)) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20, 25, 30)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    title    = "A. Random-matrix spike detection",
    subtitle = sprintf("P = %d parcels, df = %d; off-diag RMSE: %.3f (raw) -> %.3f (POET)", P, df_eff, frob_raw, frob_poet),
    x        = "Eigenvalue rank (k)",
    y        = expression("Sample eigenvalue " * lambda[k]),
    color    = NULL
  ) +
  theme_poet +
  theme(legend.text = element_text(size = 9.2))

# -------------------------------------------------------------------------
# Part 2: Dynamic Covariance Regression + LTZ Spatial Chi-Square Gating
# -------------------------------------------------------------------------
set.seed(20261125L)
n_tr <- 40L
env_fn <- function(x) 0.10 + 1.45 * exp(-0.5 * ((x - 12) / 2.2)^2) +
                             1.65 * exp(-0.5 * ((x - 30) / 2.4)^2)


load_g   <- seq(0.65, 0.85, length.out = P)
m_frames <- 300L
t_fr     <- seq(1, n_tr, length.out = m_frames)
w_fr_raw  <- env_fn(t_fr)
w_fr_rvt  <- as.numeric(scale(w_fr_raw))
w_fr_null <- as.numeric(scale(rnorm(m_frames)))
env_amp   <- pmax(0.05, 0.20 + 0.55 * (w_fr_raw - min(w_fr_raw)) / diff(range(w_fr_raw)))

f_dyn <- rnorm(m_frames) * sqrt(env_amp)
E_dyn <- outer(f_dyn, load_g) + matrix(rnorm(m_frames * P, sd = 0.60), m_frames, P)

E_sq   <- scale(E_dyn^2, center = TRUE, scale = TRUE)
z_rvt  <- as.numeric(crossprod(E_sq, w_fr_rvt))  / sqrt(m_frames)
z_null <- as.numeric(crossprod(E_sq, w_fr_null)) / sqrt(m_frames)
R_dyn_est <- medusa:::poet_regularize_cor(stats::cor(E_dyn), df_u = m_frames)$R

gate_rvt  <- medusa:::.ltz_from_cor_quad(z_rvt,  R_dyn_est, p_lo = 0.25, p_hi = 0.01)
gate_null <- medusa:::.ltz_from_cor_quad(z_null, R_dyn_est, p_lo = 0.25, p_hi = 0.01)

seed_i <- 1L; targ_i <- 2L
c_true <- env_amp * load_g[seed_i] * load_g[targ_i]
xy_fr  <- E_dyn[, seed_i] * E_dyn[, targ_i]
fit_cov <- lm(xy_fr ~ w_fr_rvt + w_fr_null)
b0     <- coef(fit_cov)[1]
b_rvt  <- coef(fit_cov)[2]
b_null <- coef(fit_cov)[3]

c_stat   <- rep(mean(xy_fr), m_frames)
c_covreg <- b0 + gate_rvt * b_rvt * w_fr_rvt + gate_null * b_null * w_fr_null

dyn_df <- rbind(
  data.frame(frame = seq_len(m_frames), cov = c_true,   method = "Generating covariance"),
  data.frame(frame = seq_len(m_frames), cov = c_covreg, method = "Gated pair regression"),
  data.frame(frame = seq_len(m_frames), cov = c_stat,   method = "Stationary-only average")
)
dyn_df$method <- factor(dyn_df$method, levels = c(
  "Generating covariance", "Gated pair regression", "Stationary-only average"
))

ymax_dyn <- max(c_true, c_covreg)
p_dyn <- ggplot(dyn_df, aes(x = frame, y = cov, color = method, linetype = method)) +
  geom_line(linewidth = 1.05) +
  annotate("label", x = 155, y = ymax_dyn * 0.96,
           label = sprintf("LTZ spatial gate:\nRVT burst = %.2f\nNull decoy = %.2f", gate_rvt, gate_null),
           size = 2.9, color = ink, fill = "#F4F7F9", linewidth = 0.25, lineheight = 1.12) +
  scale_color_manual(values = c(ink, teal, gray)) +
  scale_linetype_manual(values = c("dashed", "solid", "dotted")) +
  scale_y_continuous(limits = c(0, ymax_dyn * 1.20), expand = expansion(mult = c(0.02, 0.04))) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    title    = "B. When shared noise varies over a run",
    subtitle = "One parcel pair with a measured variance envelope",
    x        = "Acquisition frame",
    y        = "Framewise noise covariance",
    color    = NULL, linetype = NULL
  ) +
  theme_poet +
  theme(legend.text = element_text(size = 9.2))

# -------------------------------------------------------------------------
# Part 3: POET Matrix Decomposition Heatmaps (5 compact panels across Row 2)
# -------------------------------------------------------------------------
mat_to_df <- function(M, label, hide_diag = TRUE) {
  M2 <- M
  if (hide_diag) diag(M2) <- NA_real_
  data.frame(
    row = as.vector(row(M2)),
    col = as.vector(col(M2)),
    val = as.vector(M2),
    panel = label
  )
}

heat_df <- rbind(
  mat_to_df(R_true, "1. Ground truth\nR_true"),
  mat_to_df(R_raw,  "2. Sample corr\nR_hat (noisy)"),
  mat_to_df(LL_hat, sprintf("3. Global factors\nK = %d", poet_fit$K)),
  mat_to_df(S_hat,  "4. Sparse local\nS_loc (thresholded)"),
  mat_to_df(R_poet, "5. Adaptive-POET\nProj(L L' + S_loc)")
)
heat_df$panel <- factor(heat_df$panel, levels = unique(heat_df$panel))

p_heat <- ggplot(heat_df, aes(x = col, y = row, fill = val)) +
  geom_raster() +
  facet_wrap(~panel, nrow = 1) +
  scale_y_reverse(breaks = c(1, 15, 30)) +
  scale_x_continuous(breaks = c(1, 15, 30)) +
  coord_equal() +
  scale_fill_gradient2(
    low = blue, mid = "white", high = orange, midpoint = 0,
    limits = c(-0.75, 0.75), oob = scales::squish,
    na.value = "#E6EDF1",
    breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
    guide = guide_colourbar(barwidth = grid::unit(7.5, "cm"), barheight = grid::unit(0.35, "cm"))
  ) +
  labs(
    title    = "C. Step-by-step Adaptive-POET decomposition of cross-parcel noise correlation (diagonals omitted)",
    subtitle = "Each square represents one parcel pair. All five panels use the same color scale.",
    x        = "Parcel index",
    y        = "Parcel index",
    fill     = "Off-diagonal correlation"
  ) +
  theme_poet +
  theme(
    strip.text = element_text(size = 9.0, face = "bold", lineheight = 1.08),
    panel.grid = element_blank(),
    axis.text  = element_text(size = 8.5)
  )

fig <- (p_scree | p_dyn) / p_heat +
  plot_layout(heights = c(1.08, 1.0))

ggsave(
  "figures/07-poet-covreg.png",
  fig,
  width = 10.6, height = 7.6, units = "in", dpi = 180, bg = "white"
)

saveRDS(
  list(
    P = P, df_eff = df_eff, K_sel = poet_fit$K, mp_edge = mp_edge,
    frob_raw = frob_raw, frob_poet = frob_poet,
    gate_rvt = gate_rvt, gate_null = gate_null,
    mean_local_true = mean_local_true, mean_local_poet = mean_local_poet
  ),
  "poet_demo_metrics.rds"
)

cat(sprintf(
  "Saved figures/07-poet-covreg.png | K=%d | MP+TW=%.2f | RMSE raw=%.3f -> POET=%.3f | LTZ gates: RVT=%.2f, Null=%.2f\n",
  poet_fit$K, mp_edge, frob_raw, frob_poet, gate_rvt, gate_null
))
