# Generate ggbrain MNI152 brain maps & quantitative companion visuals for
# book/medusa_spatial_parcel_pooling/medusa_spatial_parcel_pooling.qmd.

if (!nzchar(Sys.getenv("KMP_USE_SHM"))) {
  Sys.setenv(KMP_USE_SHM = "0", OMP_NUM_THREADS = "1")
}
.args <- commandArgs(trailingOnly = FALSE)
.farg <- grep("^--file=", .args, value = TRUE)
if (length(.farg) == 1L) setwd(dirname(normalizePath(sub("^--file=", "", .farg))))

library(medusa)
library(ggbrain)
library(RNifti)
library(ggplot2)
library(patchwork)
library(data.table)

dir.create("figures", showWarnings = FALSE)
d <- readRDS("spatial_pooling_examples.rds")
vinfo <- copy(d$vinfo)
vmap  <- d$voxel_map
geom  <- d$geometry

# Book palette & clean theme
col_blue   <- "#2166AC"
col_navy   <- "#245D80"
col_teal   <- "#176E66"
col_orange <- "#B35806"
col_grey   <- "#5C6B73"

theme_book <- theme_minimal(base_size = 12.5) +
  theme(
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold", size = 13.5, color = "#172937"),
    plot.subtitle     = element_text(size = 11, color = "#405362", margin = margin(b = 8)),
    strip.text        = element_text(face = "bold", size = 11.5, color = "#172937"),
    legend.position   = "bottom",
    plot.margin       = margin(8, 12, 8, 10)
  )

# Locate bundled MNI152 2009c 2mm template & Schaefer-200 atlas
ext_dir    <- system.file("extdata", package = "ggbrain")
t1_path    <- file.path(ext_dir, "mni_template_2009c_2mm.nii.gz")
atlas_path <- file.path(ext_dir, "Schaefer_200_7networks_2009c_2mm.nii.gz")
t1_nii     <- readNifti(t1_path)
atlas_nii  <- readNifti(atlas_path)

# Temporary directory for 3D NIfTI stat maps passed to ggbrain
nii_dir <- tempfile("medusa_ggbrain_maps_", tmpdir = "/tmp")
dir.create(nii_dir, recursive = TRUE)
on.exit(unlink(nii_dir, recursive = TRUE), add = TRUE)

# Write slab-masked Schaefer atlas so outlines are drawn cleanly on z = 8 mm
atlas_slab_arr <- array(0L, dim = dim(atlas_nii))
atlas_slab_arr[, , d$slab_k] <- atlas_nii[, , d$slab_k]
atlas_slab_file <- file.path(nii_dir, "atlas_slab.nii.gz")
writeNifti(updateNifti(atlas_slab_arr, template = atlas_nii), atlas_slab_file)

# Reconstruct a lightweight MEDuSA reference object carrying voxel_map + geometry
# so write_nifti_ts() and backproject_medusa() write 3D NIfTIs directly.
ref_obj <- structure(
  list(voxel_map = vmap, geometry = geom, tr = d$tr),
  class = "medusa_nifti_ts"
)

write_vox_map <- function(vals, stem) {
  fpath <- file.path(nii_dir, paste0(stem, ".nii.gz"))
  mat <- matrix(as.numeric(vals), nrow = 1L, dimnames = list(NULL, vmap$voxel))
  write_nifti_ts(mat, fpath, reference = ref_obj, tr = d$tr, fill = NA_real_)
  fpath
}

# Helper to build one axial z = 8 mm ggbrain ggplot panel
brain_slice_panel <- function(overlay_file, title_txt, fill_scale, limits = NULL,
                              breaks = NULL,
                              show_legend = TRUE, bisided = FALSE,
                              sig_outline_file = NULL, sig_color = "#111111") {
  imgs <- c(underlay = t1_path, atlas = atlas_slab_file, stat = overlay_file)
  if (!is.null(sig_outline_file)) imgs <- c(imgs, sig = sig_outline_file)

  def_str <- if (bisided) "stat[abs(stat) > 1e-6]" else "stat[!is.na(stat)]"
  brk_fn  <- if (is.null(breaks)) {
    scales::breaks_extended(5)
  } else if (is.function(breaks)) {
    breaks
  } else {
    force(breaks)
    function(x) breaks
  }

  g <- ggbrain(bg_color = "white", text_color = "#172937", base_size = 12) +
    images(imgs) +
    slices(
      "z = 8",
      title = title_txt,
      bg_color = "white",
      text_color = "#172937",
      theme_custom = theme(
        plot.title = element_text(face = "bold", size = 12, hjust = 0.5, color = "#172937"),
        legend.position = "bottom",
        legend.key.width = grid::unit(0.95, "cm"),
        legend.key.height = grid::unit(0.30, "cm"),
        plot.margin = margin(4, 4, 4, 4)
      )
    ) +
    geom_brain(
      definition = "underlay[underlay > 8]",
      fill_scale = scale_fill_gradient(low = "#30363D", high = "#E6ECF0"),
      show_legend = FALSE
    ) +
    geom_brain(
      definition  = def_str,
      fill_scale  = fill_scale,
      limits      = limits,
      breaks      = brk_fn,
      show_legend = show_legend
    ) +
    geom_outline(
      definition  = "atlas[atlas > 0]",
      mapping     = aes(group = value),
      outline     = "#1F2933",
      size        = 1,
      alpha       = 0.35,
      blur_edge   = 0.4,
      show_legend = FALSE
    )

  if (!is.null(sig_outline_file)) {
    g <- g + geom_outline(
      definition  = "sig[!is.na(sig) & sig > 0]",
      mapping     = aes(group = value),
      outline     = sig_color,
      size        = 1,
      dil_ero     = 1L,
      alpha       = 0.98,
      show_legend = FALSE
    )
  }

  render(g)[[1L]]
}

save_fig <- function(stem, p, w = 12, h = 5.6) {
  ggsave(
    filename = file.path("figures", paste0(stem, ".png")),
    plot     = p,
    width    = w,
    height   = h,
    dpi      = 185,
    bg       = "white"
  )
}

# ==============================================================================
# FIGURE 1: Geodesic Ribbon Graph vs 3D Euclidean Shortcuts Across Sulcal CSF
# ==============================================================================
f_pk_true <- write_vox_map(d$hrf_peaks$truth,      "pk_truth")
f_pk_euc  <- write_vox_map(d$hrf_peaks$graph_euc,  "pk_euc")
f_pk_geo  <- write_vox_map(d$hrf_peaks$graph_soft, "pk_geo")

pk_scale <- function(nm = "HRF peak (s)") {
  scale_fill_viridis_c(name = nm, option = "plasma", oob = scales::squish)
}
pk_lim <- c(2.4, 7.2)
pk_brk <- c(3, 4, 5, 6, 7)

p1_a <- brain_slice_panel(f_pk_true, "Generating HRF peak", pk_scale(), pk_lim, breaks = pk_brk, show_legend = FALSE)
p1_b <- brain_slice_panel(f_pk_euc,  "3D Euclidean graph (bleeds across CSF)", pk_scale(), pk_lim, breaks = pk_brk, show_legend = TRUE)
p1_c <- brain_slice_panel(f_pk_geo,  "Geodesic ribbon graph (stops at CSF)",   pk_scale(), pk_lim, breaks = pk_brk, show_legend = FALSE)

# Zoomed cortical lattice inset on posterior interhemispheric fissure (Left vs Right Calcarine, z = 8 mm)
zoom_xlim <- c(-24, 24)
zoom_ylim <- c(-94, -52)

# Extract background T1 slice at k = 44 for anatomical CSF cleft context
ij_grid <- expand.grid(i = seq_len(dim(t1_nii)[1L]), j = seq_len(dim(t1_nii)[2L]))
w_z8 <- voxelToWorld(cbind(ij_grid$i, ij_grid$j, 44L), t1_nii)
t1_z8_dt <- data.table(x = w_z8[, 1], y = w_z8[, 2], t1 = as.numeric(t1_nii[, , 44L]))
t1_zoom  <- t1_z8_dt[x >= zoom_xlim[1] & x <= zoom_xlim[2] & y >= zoom_ylim[1] & y <= zoom_ylim[2]]

v_z8_zoom <- vinfo[k == 44L & x >= zoom_xlim[1] & x <= zoom_xlim[2] & y >= zoom_ylim[1] & y <= zoom_ylim[2]]
v_z8_zoom[, bank_type := fifelse(
  parcel == 10L,  "Active L Calcarine (2.8 s HRF)",
  fifelse(parcel == 109L, "Silent R Calcarine (5.1 s HRF)", "Adjacent ribbon cortex")
)]

sc_zoom <- d$shortcuts_z8[
  x1 >= zoom_xlim[1] & x1 <= zoom_xlim[2] & y1 >= zoom_ylim[1] & y1 <= zoom_ylim[2] &
  x2 >= zoom_xlim[1] & x2 <= zoom_xlim[2] & y2 >= zoom_ylim[1] & y2 <= zoom_ylim[2]
]

p1_d <- ggplot() +
  geom_raster(data = t1_zoom, aes(x = x, y = y, fill = t1), show.legend = FALSE, alpha = 0.42) +
  scale_fill_gradient(low = "#1E252B", high = "#EBEEF0") +
  geom_point(
    data = v_z8_zoom,
    aes(x = x, y = y, color = bank_type),
    size = 1.85, alpha = 0.92
  ) +
  geom_segment(
    data = sc_zoom,
    aes(x = x1, y = y1, xend = x2, yend = y2),
    color = col_orange, linewidth = 0.55, alpha = 0.78
  ) +
  geom_segment(
    data = d$ego_geo_calc,
    aes(x = x1, y = y1, xend = x2, yend = y2),
    color = col_teal, linewidth = 0.85, alpha = 0.95
  ) +
  geom_point(
    data = vinfo[vnum == d$seed_calc],
    aes(x = x, y = y),
    shape = 21, size = 3.6, stroke = 1.1, fill = "#FFD166", color = "#172937"
  ) +
  coord_fixed(xlim = zoom_xlim, ylim = zoom_ylim, expand = FALSE) +
  scale_color_manual(
    values = c(
      "Active L Calcarine (2.8 s HRF)" = col_blue,
      "Silent R Calcarine (5.1 s HRF)" = "#5A6B7C",
      "Adjacent ribbon cortex"         = "#9AA8B2"
    )
  ) +
  labs(
    title    = "Posterior midline zoom (z = 8 mm)",
    subtitle = "Orange: Euclidean CSF shortcuts\nTeal: Geodesic seed edges",
    x = "MNI x (mm)", y = "MNI y (mm)", color = NULL
  ) +
  guides(color = guide_legend(ncol = 1)) +
  theme_book +
  theme(legend.text = element_text(size = 9.5))

fig1 <- (p1_a | p1_b | p1_c | p1_d) +
  plot_layout(widths = c(1, 1, 1, 1.12)) +
  plot_annotation(
    title    = "In-mask geodesic lattice Dijkstra prevents cross-sulcal and cross-fissure smoothing",
    subtitle = "MNI152 z = 8 mm cortical ribbon with Schaefer-200 parcel boundaries (FWHM = 6 mm)",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15.5, color = "#172937"),
      plot.subtitle = element_text(size = 12, color = "#405362")
    )
  )
save_fig("01-geodesic-vs-euclidean-topology", fig1, w = 14.2, h = 5.4)

# ==============================================================================
# FIGURE 2: Step 1 HRF Peak Recovery Across Voxel, Parcel, & Bilateral Graph Pooling
# ==============================================================================
f_pk_vox <- write_vox_map(d$hrf_peaks$voxel_only,  "pk_vox")
f_pk_par <- write_vox_map(d$hrf_peaks$parcel_pool, "pk_par")

p2_a <- brain_slice_panel(f_pk_true, "1. Generating HRF peak",         pk_scale(), pk_lim, breaks = pk_brk, show_legend = FALSE)
p2_b <- brain_slice_panel(f_pk_vox,  "2. Single-voxel ('voxel')",      pk_scale(), pk_lim, breaks = pk_brk, show_legend = FALSE)
p2_c <- brain_slice_panel(f_pk_par,  "3. Group parcel ('parcel')",     pk_scale(), pk_lim, breaks = pk_brk, show_legend = TRUE)
p2_d <- brain_slice_panel(f_pk_geo,  "4. Bilateral graph ('soft')",    pk_scale(), pk_lim, breaks = pk_brk, show_legend = FALSE)

# Spatial transect along Left Posterior Temporal ribbon (Parcels 79 & 78, y in [-66, -24])
temp_dt <- d$hrf_peaks[parcel %in% c(78L, 79L)]
# Nominal Schaefer border y-coordinate between 79 and 78
y_atlas_border <- 0.5 * (max(temp_dt[parcel == 79L, y]) + min(temp_dt[parcel == 78L, y]))

transect_long <- melt(
  temp_dt,
  id.vars       = c("voxel", "parcel", "y", "truth"),
  measure.vars  = c("voxel_only", "parcel_pool", "graph_hard", "graph_soft"),
  variable.name = "method",
  value.name    = "peak_sec"
)
transect_long[, method := factor(
  method,
  levels = c("voxel_only", "parcel_pool", "graph_hard", "graph_soft"),
  labels = c("Single-voxel", "Group parcel", "Graph (hard parcel wall)", "Graph (soft parcel + bilateral)")
)]

transect_means <- transect_long[, .(mean_peak = mean(peak_sec)), by = .(method, y)]
truth_curve    <- unique(temp_dt[, .(y, truth)])[order(y)]

p2_transect <- ggplot(transect_long, aes(x = y, y = peak_sec)) +
  geom_vline(xintercept = y_atlas_border, linetype = "dotted", color = "#6C7A86", linewidth = 0.8) +
  geom_vline(xintercept = -43, linetype = "dashed", color = "#172937", linewidth = 0.85) +
  annotate("text", x = y_atlas_border - 0.8, y = 8.45, label = "Nominal Schaefer\n79 | 78 border",
           hjust = 1, size = 3.2, color = "#52606D") +
  annotate("text", x = -43 + 0.8, y = 8.45, label = "Participant's true\nvascular border (y = -43 mm)",
           hjust = 0, size = 3.2, fontface = "bold", color = "#172937") +
  geom_point(data = transect_long[method == "Single-voxel"], color = "#A0AEB8", alpha = 0.45, size = 1.2) +
  geom_step(data = truth_curve, aes(x = y, y = truth), color = "#111111", linewidth = 1.05, inherit.aes = FALSE) +
  geom_line(
    data = transect_means[method != "Single-voxel"],
    aes(x = y, y = mean_peak, color = method),
    linewidth = 1.15
  ) +
  scale_color_manual(values = c(
    "Group parcel"                    = col_orange,
    "Graph (hard parcel wall)"        = "#7B5EA7",
    "Graph (soft parcel + bilateral)" = col_teal
  )) +
  coord_cartesian(ylim = c(2.3, 8.9)) +
  labs(
    title    = "Left temporal cortical transect across a shifted subject functional/vascular boundary (Parcels 79 & 78)",
    subtitle = "Grey points: single-voxel estimates | Black step: generating HRF peak | Colored curves: mean cross-fitted HRF peak at each MNI y",
    x = "Anterior-posterior MNI coordinate y (mm)",
    y = "Estimated HRF peak time (s)",
    color = NULL
  ) +
  theme_book

fig2 <- ((p2_a | p2_b | p2_c | p2_d) / p2_transect) +
  plot_layout(heights = c(1.18, 0.92)) +
  plot_annotation(
    title    = "Step 1: Bilateral graph pooling denoises HRF latency without blurring shifted parcel borders",
    subtitle = "Out-of-sample cross-run candidate evidence with soft cross-parcel relaxation and anti-correlation stopping",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15.5, color = "#172937"),
      plot.subtitle = element_text(size = 12, color = "#405362")
    )
  )
save_fig("02-step1-hrf-pooling-modes", fig2, w = 13.8, h = 8.8)

# ==============================================================================
# FIGURE 3: Bilateral Noise Smoothing at Pial-Vein / CSF Cortical Borders
# ==============================================================================
f_rho_true  <- write_vox_map(d$noise_dt$rho_true,   "rho_true")
f_rho_vox   <- write_vox_map(d$noise_dt$rho_voxel,  "rho_vox")
f_rho_par   <- write_vox_map(d$noise_dt$rho_parcel, "rho_par")
f_rho_graph <- write_vox_map(d$noise_dt$rho_graph,  "rho_graph")

rho_scale <- function(nm = "AR(1) rho") {
  scale_fill_viridis_c(name = nm, option = "mako", direction = 1, oob = scales::squish)
}
rho_lim <- c(0.10, 0.88)
rho_brk <- c(0.2, 0.4, 0.6, 0.8)

p3_a <- brain_slice_panel(f_rho_true,  "Generating AR(1) rho",       rho_scale(), rho_lim, breaks = rho_brk, show_legend = FALSE)
p3_b <- brain_slice_panel(f_rho_vox,   "Single-voxel REML",          rho_scale(), rho_lim, breaks = rho_brk, show_legend = FALSE)
p3_c <- brain_slice_panel(f_rho_par,   "Parcel-pooled REML",         rho_scale(), rho_lim, breaks = rho_brk, show_legend = TRUE)
p3_d <- brain_slice_panel(f_rho_graph, "Bilateral graph REML",       rho_scale(), rho_lim, breaks = rho_brk, show_legend = FALSE)

rho_long <- melt(
  d$noise_dt,
  id.vars       = c("voxel", "is_pial_rim", "rho_true"),
  measure.vars  = c("rho_voxel", "rho_parcel", "rho_graph"),
  variable.name = "estimator",
  value.name    = "rho_hat"
)
rho_long[, tissue := fifelse(is_pial_rim, "Pial-vein / CSF rim (true rho = 0.76)", "Interior gray matter (true rho = 0.28)")]
rho_long[, estimator := factor(
  estimator,
  levels = c("rho_voxel", "rho_parcel", "rho_graph"),
  labels = c("Single-voxel REML", "Parcel-pooled REML", "Bilateral graph REML")
)]
rho_ref <- data.table(
  tissue = c("Pial-vein / CSF rim (true rho = 0.76)", "Interior gray matter (true rho = 0.28)"),
  truth  = c(0.76, 0.28)
)

p3_dist <- ggplot(rho_long, aes(x = estimator, y = rho_hat, fill = estimator)) +
  geom_violin(alpha = 0.78, color = "#243B53", linewidth = 0.35, scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", color = "#172937", linewidth = 0.4) +
  geom_hline(data = rho_ref, aes(yintercept = truth), linetype = "dashed", color = "#111111", linewidth = 0.8) +
  facet_wrap(~tissue, nrow = 1) +
  scale_fill_manual(values = c(
    "Single-voxel REML"    = "#9FB3C8",
    "Parcel-pooled REML"   = col_orange,
    "Bilateral graph REML" = col_teal
  )) +
  labs(
    title    = "Tissue-compartment separation of temporal autocorrelation",
    subtitle = "Parcel pooling contaminates deep gray matter with pial-rim autocorrelation; bilateral graph weights preserve both compartments",
    x = NULL, y = "Estimated AR(1) parameter", fill = NULL
  ) +
  theme_book +
  theme(legend.position = "none")

fig3 <- ((p3_a | p3_b | p3_c | p3_d) / p3_dist) +
  plot_layout(heights = c(1.15, 0.85)) +
  plot_annotation(
    title    = "Bilateral REML noise smoothing separates pial-vein borders from interior gray matter",
    subtitle = "Joint kernel in spatial distance, parcel membership, and pilot (rho, log sigma) noise similarity",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15.5, color = "#172937"),
      plot.subtitle = element_text(size = 12, color = "#405362")
    )
  )
save_fig("03-bilateral-pial-rim-noise", fig3, w = 13.8, h = 8.4)

# ==============================================================================
# FIGURE 4: Step 2 SNR-Stratified Regularization (lambda_mode = 'auto' vs 'global')
# ==============================================================================
f_lam_auto <- write_vox_map(log10(d$vox_lag0_maps$lambda_auto), "log10_lambda_auto")
lam_vals   <- log10(d$vox_lag0_maps$lambda_auto)
lam_lim    <- range(lam_vals) + c(-0.05, 0.05)

p4_map <- brain_slice_panel(
  f_lam_auto,
  "Adaptive log10(lambda) map (z = 8 mm)",
  scale_fill_viridis_c(name = "log10(lambda)", option = "cividis"),
  limits      = lam_lim,
  show_legend = TRUE
)

strata_dt <- copy(d$step2_strata)
strata_dt[, tier_label := factor(
  sprintf("Tier %d\n(%s, n=%d)", tier, label, n_voxels),
  levels = sprintf("Tier %d\n(%s, n=%d)", tier, label, n_voxels)
)]

p4_pava <- ggplot(strata_dt, aes(x = tier, y = lambda)) +
  geom_hline(
    yintercept = strata_dt$lambda_global[1L],
    linetype = "dashed", color = col_orange, linewidth = 0.95
  ) +
  annotate(
    "text", x = 1.05, y = strata_dt$lambda_global[1L] * 1.48,
    label = sprintf("Global lambda = %.2f (under-smooths quiet cortex)", strata_dt$lambda_global[1L]),
    hjust = 0, size = 3.35, color = col_orange, fontface = "bold"
  ) +
  geom_line(color = col_blue, linewidth = 1.1) +
  geom_point(size = 3.6, color = col_blue) +
  scale_x_continuous(
    breaks = strata_dt$tier,
    labels = as.character(strata_dt$tier_label),
    expand = expansion(mult = c(0.10, 0.10))
  ) +
  scale_y_log10() +
  labs(
    title    = "Monotone PAVA lambda across OOS task-SNR tiers",
    subtitle = "Quiet background receives heavier smoothing; active foci stay sharp",
    x = "Cross-fitted task-R² tier", y = "Penalty lambda (log scale)"
  ) +
  theme_book

tr_sub <- d$step2_traces[time <= 95]
trace_long <- rbind(
  data.table(time = tr_sub$time, value = tr_sub$null_global, curve = "Global lambda",     panel = "Null background voxel (suppressing high-freq noise)"),
  data.table(time = tr_sub$time, value = tr_sub$null_auto,   curve = "Stratified 'auto'", panel = "Null background voxel (suppressing high-freq noise)"),
  data.table(time = tr_sub$time, value = tr_sub$hot_global,  curve = "Global lambda",     panel = "Focal R Insula hotspot (preserving phasic spikes)"),
  data.table(time = tr_sub$time, value = tr_sub$hot_auto,    curve = "Stratified 'auto'", panel = "Focal R Insula hotspot (preserving phasic spikes)")
)
# Scale true neural impulses onto noise_sd units for visual timing reference in hotspot panel
hot_sig_div <- mean(d$noise_dt[match(vinfo[zone == "Focal hotspot (R Insula core)", voxel], voxel), sig_graph])
spikes_dt <- tr_sub[hot_true > 0, .(
  time  = time,
  value = hot_true / hot_sig_div * 0.22,
  panel = "Focal R Insula hotspot (preserving phasic spikes)"
)]

p4_traces <- ggplot(trace_long, aes(x = time, y = value, color = curve)) +
  geom_hline(yintercept = 0, color = "#CBD2D9", linewidth = 0.4) +
  geom_segment(
    data = spikes_dt,
    aes(x = time, xend = time, y = 0, yend = value),
    inherit.aes = FALSE, color = "#111111", alpha = 0.45, linewidth = 0.8
  ) +
  geom_line(data = trace_long[curve == "Global lambda"], linewidth = 1.15, alpha = 0.82) +
  geom_line(data = trace_long[curve == "Stratified 'auto'"], linewidth = 0.78, alpha = 0.96) +
  facet_wrap(~panel, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c("Global lambda" = col_orange, "Stratified 'auto'" = col_teal)) +
  labs(
    title    = "Deconvolved continuous activity (Run 1, first 95 s)",
    subtitle = "Hotspot curves coincide (Tier 4 keeps lambda = 0.01); Tier 1 suppresses null noise 3x",
    x = "Time (s)", y = "Deconvolved activity (noise-SD units)", color = NULL
  ) +
  theme_book

fig4 <- (p4_map | (p4_pava / p4_traces + plot_layout(heights = c(0.85, 1.35)))) +
  plot_layout(widths = c(1.02, 1.38)) +
  plot_annotation(
    title    = "Step 2: SNR-stratified GCV tuning suppresses null-cortex noise without blunting focal activity",
    subtitle = "deconvolve_crossrun_hrf(..., lambda_mode = 'auto') with cross-run FDR gating and weighted PAVA monotonicity",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15.5, color = "#172937"),
      plot.subtitle = element_text(size = 12, color = "#405362")
    )
  )
save_fig("04-step2-snr-stratified-lambda", fig4, w = 13.6, h = 7.8)

# ==============================================================================
# FIGURE 5: Step 3 Voxelwise Bilateral Graph-EB (tau2_shrinkage)
# ==============================================================================
tau2_cap <- 0.65
f_tau_true <- write_vox_map(pmin(tau2_cap, d$vox_lag0_maps$tau2_true),    "tau2_true")
f_tau_none <- write_vox_map(pmin(tau2_cap, d$vox_lag0_maps$tau2_none),    "tau2_none")
f_tau_eb   <- write_vox_map(pmin(tau2_cap, d$vox_lag0_maps$tau2_eb),      "tau2_eb")
f_tau_geb  <- write_vox_map(pmin(tau2_cap, d$vox_lag0_maps$tau2_grapheb), "tau2_grapheb")

tau_scale <- function(nm = "tau2_bio") {
  scale_fill_viridis_c(name = nm, option = "inferno", direction = 1, oob = scales::squish)
}
tau_lim <- c(0, tau2_cap)
tau_brk <- c(0, 0.2, 0.4, 0.6)

p5_a <- brain_slice_panel(f_tau_true, "Generating tau2_bio",             tau_scale(), tau_lim, breaks = tau_brk, show_legend = FALSE)
p5_b <- brain_slice_panel(f_tau_none, "Raw voxel MoM ('none')",          tau_scale(), tau_lim, breaks = tau_brk, show_legend = FALSE)
p5_c <- brain_slice_panel(f_tau_eb,   "Global EB ('eb')",                tau_scale(), tau_lim, breaks = tau_brk, show_legend = TRUE)
p5_d <- brain_slice_panel(f_tau_geb,  "Bilateral Graph-EB ('graph_eb')", tau_scale(), tau_lim, breaks = tau_brk, show_legend = FALSE)

t3 <- copy(d$step3_tau2_summary)[method != "Euclidean Graph-EB"]
t3[, short_method := factor(
  method,
  levels = c("Raw voxel MoM ('none')", "Global Empirical Bayes ('eb')", "Bilateral Geodesic Graph-EB ('graph_eb')"),
  labels = c("Raw MoM ('none')", "Global EB ('eb')", "Bilateral Graph-EB")
)]

t3_long <- rbind(
  data.table(method = t3$short_method, metric = "Active-cortex zero-collapse fraction\n(lower = fewer false tau2 = 0 holes)", val = t3$active_zero_collapse, target = 0),
  data.table(method = t3$short_method, metric = "Null-cortex mean tau2_bio\n(lower = less positive truncation bias)",        val = t3$null_mean_tau2,       target = 0),
  data.table(method = t3$short_method, metric = "Focal R Insula hotspot mean tau2_bio\n(dashed = true hotspot tau2 = 0.61)", val = t3$hotspot_mean_tau2,    target = 0.78^2)
)
t3_long[, metric := factor(metric, levels = unique(metric))]

p5_metrics <- ggplot(t3_long, aes(x = method, y = val, fill = method)) +
  geom_col(width = 0.58, alpha = 0.90, color = "#172937", linewidth = 0.3) +
  geom_hline(aes(yintercept = target), linetype = "dashed", color = "#111111", linewidth = 0.75) +
  facet_wrap(~metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(
    "Raw MoM ('none')"   = "#9FB3C8",
    "Global EB ('eb')"   = col_orange,
    "Bilateral Graph-EB" = col_teal
  )) +
  labs(
    title    = "Pre-truncation bilateral shrinkage eliminates both active zero-collapse and null inflation",
    subtitle = "Efron-Morris limited translation prevents focal hotspots (R Insula core) from over-shrinking into quiet surround",
    x = NULL, y = NULL, fill = NULL
  ) +
  theme_book +
  theme(legend.position = "none", axis.text.x = element_text(angle = 12, hjust = 1))

fig5 <- ((p5_a | p5_b | p5_c | p5_d) / p5_metrics) +
  plot_layout(heights = c(1.18, 0.88)) +
  plot_annotation(
    title    = "Step 3 (Voxelwise): Bilateral Graph-EB stabilizes latent trial-variance estimation",
    subtitle = "Shrinks scale-free variance ratios prior to zero truncation without smoothing neural timecourses or trial slopes",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15.5, color = "#172937"),
      plot.subtitle = element_text(size = 12, color = "#405362")
    )
  )
save_fig("05-step3-voxelwise-graph-eb-tau2", fig5, w = 13.8, h = 8.5)

# ==============================================================================
# FIGURE 6: Voxelwise (write_nifti_ts) vs Parcel-Wise (backproject_medusa & Step 5)
# ==============================================================================
beta_cap <- 0.95
clamp <- function(x, lim) pmax(-lim, pmin(lim, x))

f_beta_true <- write_vox_map(clamp(d$vox_lag0_maps$beta_true,    beta_cap), "beta_true")
f_beta_vox  <- write_vox_map(clamp(d$vox_lag0_maps$beta_grapheb, beta_cap), "beta_vox")

# Backproject single-subject Step 4 parcel slopes via backproject_medusa()
par_sub_df <- copy(d$c_parcel_lag0)
par_sub_df[, `:=`(parcel_id = as.integer(parcel), beta_clamped = clamp(beta_gls, beta_cap))]
bp_par_sub <- backproject_medusa(
  coef_df    = as.data.frame(par_sub_df),
  brain_mask = ref_obj,
  plot_cols  = "beta_clamped",
  parcel_col = "parcel_id",
  time_col   = "lag",
  effect_col = "term",
  output_dir = nii_dir,
  ncpus      = 1L,
  fill       = NA_real_
)
f_beta_par <- bp_par_sub[1L]

# Backproject Step 5 12-subject cohort mixed_by() parcel slopes & BY-significant outlines
grp_df <- copy(d$group_rpe)
grp_df[, `:=`(
  beta_clamped = clamp(estimate, beta_cap),
  sig_parcel   = fifelse(padj_BY_term <= 0.05, as.numeric(parcel_int), NA_real_)
)]
bp_grp <- backproject_medusa(
  coef_df    = as.data.frame(grp_df),
  brain_mask = ref_obj,
  plot_cols  = c("beta_clamped", "sig_parcel"),
  parcel_col = "parcel_int",
  time_col   = "evt_time",
  effect_col = "term",
  output_dir = file.path(nii_dir, "grp"),
  ncpus      = 1L,
  fill       = NA_real_
)
f_beta_grp <- bp_grp[grep("beta_clamped", bp_grp)]
f_sig_grp  <- bp_grp[grep("sig_parcel",   bp_grp)]

# For signed bisided scales, replace exact 0 in generating map with tiny 0 so only non-zero or all masked voxels show?
# Actually a continuous diverging scale_fill_gradient2 centered at 0 shows null cortex as clean neutral white/ivory!
div_scale <- function(nm = "RPE slope") {
  scale_fill_gradient2(
    name     = nm,
    low      = col_blue,
    mid      = "#F7F7F5",
    high     = col_orange,
    midpoint = 0,
    oob      = scales::squish
  )
}
beta_lim <- c(-beta_cap, beta_cap)
beta_brk <- c(-0.8, -0.4, 0, 0.4, 0.8)

p6_a <- brain_slice_panel(f_beta_true, "1. Generating RPE slope",          div_scale(), beta_lim, breaks = beta_brk, show_legend = FALSE)
p6_b <- brain_slice_panel(f_beta_vox,  "2. Voxelwise Graph-EB (Subj 1)",   div_scale(), beta_lim, breaks = beta_brk, show_legend = FALSE)
p6_c <- brain_slice_panel(f_beta_par,  "3. Parcel panel Step 4 (Subj 1)",  div_scale(), beta_lim, breaks = beta_brk, show_legend = TRUE)
p6_d <- brain_slice_panel(
  f_beta_grp,
  "4. Cohort mixed_by (N=12, BY<=.05)",
  div_scale(), beta_lim, breaks = beta_brk, show_legend = FALSE,
  sig_outline_file = f_sig_grp, sig_color = "#111111"
)

vp_dt <- copy(d$vox_vs_parcel)
vp_dt[, regime_short := factor(
  regime,
  levels = c(
    "Aligned parcel (L Calcarine & L mPFC)",
    "Shifted border parcel (L Temp 78)",
    "Sub-parcel hotspot (R Insula core)",
    "Null cortex (55 background parcels)"
  ),
  labels = c(
    "Aligned homogeneous\nparcels (10 & 89)",
    "Shifted-border\nparcel (78: -/+ RPE)",
    "Sub-parcel focal\nhotspot (155 core)",
    "Null background\ncortex (55 parcels)"
  )
)]

p6_rmse <- ggplot(vp_dt, aes(x = regime_short, y = rmse_beta, fill = pathway)) +
  geom_col(position = position_dodge(width = 0.68), width = 0.60, color = "#172937", linewidth = 0.3) +
  scale_fill_manual(values = c(
    "Voxelwise (Graph + Stratified + Graph-EB)" = col_teal,
    "Parcel-wise (Step 4 multi-stratum panel)"  = col_navy
  )) +
  labs(
    title    = "Slope RMSE across anatomical regimes (lower is better)",
    subtitle = "Parcel pooling wins on aligned parcels; voxelwise wins on shifted & focal",
    x = NULL, y = "RMSE of feedback RPE slope", fill = NULL
  ) +
  theme_book

p6_det <- ggplot(vp_dt, aes(x = regime_short, y = detect_fdr05, fill = pathway)) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "#6C7A86", linewidth = 0.7) +
  geom_col(position = position_dodge(width = 0.68), width = 0.60, color = "#172937", linewidth = 0.3) +
  scale_y_continuous(labels = scales::label_percent(), limits = c(0, 1.02)) +
  scale_fill_manual(values = c(
    "Voxelwise (Graph + Stratified + Graph-EB)" = col_teal,
    "Parcel-wise (Step 4 multi-stratum panel)"  = col_navy
  )) +
  labs(
    title    = "BH-FDR <= 0.05 detection rate (power in active; FPR in null)",
    subtitle = "Opposite-signed sub-domains cancel to zero inside Parcel 78 under parcel averaging",
    x = NULL, y = "Proportion detected (q <= 0.05)", fill = NULL
  ) +
  theme_book

fig6 <- ((p6_a | p6_b | p6_c | p6_d) / (p6_rmse | p6_det + theme(legend.position = "none"))) +
  plot_layout(heights = c(1.18, 0.92)) +
  plot_annotation(
    title    = "Choosing between voxelwise Graph-EB and parcel-wise Step 4 panels",
    subtitle = "MNI152 z = 8 mm feedback RPE maps (lag = 0 s); black bold outlines on Panel 4 mark cohort BY-adjusted q <= 0.05",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15.5, color = "#172937"),
      plot.subtitle = element_text(size = 12, color = "#405362")
    )
  )
save_fig("06-voxelwise-vs-parcelwise-workflows", fig6, w = 14.2, h = 8.9)

cat("Generated all 6 figures in figures/.\n")
