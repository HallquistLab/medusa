# Volumetric MNI152 + Schaefer-200 spatial & parcel pooling simulation for the MEDuSA book.
# Exercises the full NIfTI -> build_spatial_graph -> Step 1 -> Step 2 -> Step 3 / Step 4 / Step 5
# pipeline on real cortical ribbon anatomy (k = 43:45, z = 6, 8, 10 mm).

if (!nzchar(Sys.getenv("KMP_USE_SHM"))) {
  Sys.setenv(KMP_USE_SHM = "0", OMP_NUM_THREADS = "1")
}
.args <- commandArgs(trailingOnly = FALSE)
.farg <- grep("^--file=", .args, value = TRUE)
if (length(.farg) == 1L) setwd(dirname(normalizePath(sub("^--file=", "", .farg))))

library(medusa)
library(RNifti)
library(data.table)

set.seed(20261115)
t_start <- proc.time()[3]

# ------------------------------------------------------------------------------
# 1. Load MNI152 2009c 2mm template & Schaefer-200 7-Network parcellation
# ------------------------------------------------------------------------------
ext_dir <- system.file("extdata", package = "ggbrain")
stopifnot(nzchar(ext_dir))

t1_path    <- file.path(ext_dir, "mni_template_2009c_2mm.nii.gz")
atlas_path <- file.path(ext_dir, "Schaefer_200_7networks_2009c_2mm.nii.gz")
label_path <- file.path(ext_dir, "Schaefer_200_7networks_labels.csv")

t1_nii    <- readNifti(t1_path)
atlas_nii <- readNifti(atlas_path)
atlas_lut <- fread(label_path)

# Build a 3-slice cortical gray-matter ribbon slab (k = 43:45 -> z = 6, 8, 10 mm).
# Restricting to cortical GM intensities (60 <= T1 <= 86) inside Schaefer parcels
# preserves real CSF sulcal/fissure gaps (e.g., interhemispheric & circular sulci)
# as out-of-mask barriers, mirroring a subject ribbon mask in MNI space.
slab_k <- 43:45
mask_arr <- array(0L, dim = dim(t1_nii))
mask_arr[, , slab_k] <- as.integer(
  (atlas_nii[, , slab_k] > 0L) &
  (t1_nii[, , slab_k] >= 60) &
  (t1_nii[, , slab_k] <= 86)
)

# Write mask & atlas to temporary NIfTIs for end-to-end read_nifti_ts() ingestion
tmp_dir <- tempfile("medusa_spatial_sim_")
dir.create(tmp_dir, recursive = TRUE)
on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

mask_file  <- file.path(tmp_dir, "ribbon_mask.nii.gz")
atlas_file <- file.path(tmp_dir, "schaefer200_slab.nii.gz")
writeNifti(updateNifti(mask_arr, template = t1_nii), mask_file)
writeNifti(atlas_nii, atlas_file)

# Pre-extract voxel coordinates to define anatomical ground-truth regimes
lin_idx <- which(mask_arr > 0L)
ijk_mat <- which(mask_arr > 0L, arr.ind = TRUE)
xyz_mat <- voxelToWorld(ijk_mat, t1_nii)
V <- nrow(ijk_mat)

vinfo <- data.table(
  vnum   = seq_len(V),
  i      = ijk_mat[, 1],
  j      = ijk_mat[, 2],
  k      = ijk_mat[, 3],
  x      = xyz_mat[, 1],
  y      = xyz_mat[, 2],
  z      = xyz_mat[, 3],
  parcel = as.integer(atlas_nii[lin_idx])
)

# Identify exterior/pial-vein cortical rim voxels adjacent to dark CSF (T1 < 56)
csf_arr <- t1_nii < 56
is_rim <- logical(V)
for (u in seq_len(V)) {
  ii <- vinfo$i[u]; jj <- vinfo$j[u]; kk <- vinfo$k[u]
  is_rim[u] <- csf_arr[ii - 1L, jj, kk] || csf_arr[ii + 1L, jj, kk] ||
               csf_arr[ii, jj - 1L, kk] || csf_arr[ii, jj + 1L, kk]
}
vinfo[, is_pial_rim := is_rim]

# ------------------------------------------------------------------------------
# 2. Define ground-truth anatomical, vascular, noise, & trial-modulation zones
# ------------------------------------------------------------------------------
tr <- 1.0
K_lib <- hrf_library(tr = tr, preset = "medusa")
K_pars <- attr(K_lib, "parameters")
pick_hrf <- function(delay, disp, ratio = 6) {
  d2 <- (K_pars$peak_delay - delay)^2 +
        4 * (K_pars$dispersion - disp)^2 +
        0.05 * ( ifelse(is.infinite(K_pars$undershoot_ratio), 20, K_pars$undershoot_ratio) - ratio )^2
  which.min(d2)
}

idx_fast  <- pick_hrf(4.0, 0.85, 6)
idx_med   <- pick_hrf(4.5, 0.95, 6)
idx_canon <- pick_hrf(6.0, 1.00, 6)
idx_slow  <- pick_hrf(7.8, 1.20, 6)

kernel_peak_sec <- function(K_mat, tr_val = 1.0) {
  (apply(K_mat, 2L, which.max) - 1L) * tr_val
}

# Default whole-slab background: canonical HRF, silent task modulation
vinfo[, `:=`(
  zone        = "Null background",
  hrf_idx     = idx_canon,
  amp_cue     = 0.0,
  amp_fb      = 0.0,
  beta_rpe    = 0.0,
  tau_bio     = 0.0,
  rho_true    = ifelse(is_pial_rim, 0.76, 0.28),
  sigma_true  = ifelse(is_pial_rim, 0.38, 0.21)
)]

# Zone 1A: Parcel-aligned active cortex (Left Calcarine = 10, Left Sup. Medial PFC = 89)
vinfo[parcel == 10L, `:=`(
  zone     = "Aligned active (L Calcarine)",
  hrf_idx  = idx_fast,
  amp_cue  = 1.15,
  amp_fb   = 0.95,
  beta_rpe = 0.55,
  tau_bio  = 0.46
)]
vinfo[parcel == 89L, `:=`(
  zone     = "Aligned active (L mPFC)",
  hrf_idx  = idx_slow,
  amp_cue  = 1.05,
  amp_fb   = 0.90,
  beta_rpe = 0.50,
  tau_bio  = 0.44
)]

# Zone 1B: Homologous opposing midline banks across CSF gap (Right Calcarine = 109, Right mPFC = 194)
# Cue-responsive with canonical 6.0 s HRF, but zero feedback RPE modulation and zero biological feedback variance.
vinfo[parcel %in% c(109L, 194L), `:=`(
  zone     = "Opposing sulcal null bank",
  hrf_idx  = idx_canon,
  amp_cue  = 1.00,
  amp_fb   = 0.0,
  beta_rpe = 0.0,
  tau_bio  = 0.0
)]

# Zone 2: Shifted subject functional/vascular border inside Left Posterior Temporal (78 & 79)
# Group Schaefer border between 79 & 78 sits around y = -51 mm, but this participant's true
# neurovascular & RPE-sign border sits 8 mm anteriorly at y = -43 mm (splitting parcel 78 in half).
vinfo[parcel %in% c(78L, 79L) & y <= -43, `:=`(
  zone     = "Shifted border: posterior slow (-RPE)",
  hrf_idx  = idx_slow,
  amp_cue  = 1.10,
  amp_fb   = 0.85,
  beta_rpe = -0.62,
  tau_bio  = 0.48
)]
vinfo[parcel %in% c(78L, 79L) & y > -43, `:=`(
  zone     = "Shifted border: anterior fast (+RPE)",
  hrf_idx  = idx_fast,
  amp_cue  = 1.15,
  amp_fb   = 0.90,
  beta_rpe = 0.62,
  tau_bio  = 0.48
)]

# Zone 3: Sub-parcel focal hotspot inside Right Insula (parcel 155; ~11% of parcel voxels)
dist_ins <- sqrt((vinfo$x - 42)^2 + (vinfo$y - 6)^2 + (vinfo$z - 8)^2)
vinfo[parcel == 155L, zone := "Focal parcel: silent surround"]
vinfo[parcel == 155L & dist_ins <= 4.6, `:=`(
  zone     = "Focal hotspot (R Insula core)",
  hrf_idx  = idx_med,
  amp_cue  = 1.45,
  amp_fb   = 1.15,
  beta_rpe = 0.92,
  tau_bio  = 0.78
)]

vinfo[, true_peak_sec := kernel_peak_sec(K_lib[, hrf_idx, drop = FALSE], tr)]
vinfo[, true_tau2_bio := tau_bio^2]

cat(sprintf("Constructed MNI152 slab: V = %d voxels across %d Schaefer parcels.\n",
            V, uniqueN(vinfo$parcel)))
print(vinfo[, .N, by = zone])

# ------------------------------------------------------------------------------
# 3. Simulate 3 runs of 4D NIfTI fMRI with spatially correlated AR(1) innovations
# ------------------------------------------------------------------------------
n_runs   <- 3L
n_time   <- 170L
n_trials <- 22L

# Build geodesic graph on voxel coordinates to induce realistic local spatial noise correlation
g_noise_gen <- build_spatial_graph(
  vinfo[, .(i, j, k, x, y, z)],
  spatial_fwhm   = 4.5,
  spatial_groups = vinfo$parcel,
  metric         = "geodesic"
)

smooth_innovations <- function(U_raw, g, mix = 0.55) {
  # Combine white innovations with graph-averaged neighbor innovations and re-standardize columns
  U_nbr <- matrix(0, nrow(U_raw), ncol(U_raw))
  rp <- g$row_ptr; ci <- g$col_idx
  wt <- exp(-0.5 * (g$dist / g$sigma)^2)
  for (u in seq_len(g$n_voxels)) {
    s <- rp[u]; e <- rp[u + 1L] - 1L
    if (e >= s) {
      idx_n <- ci[s:e]
      w_n   <- wt[s:e] / sum(wt[s:e])
      U_nbr[, u] <- U_raw[, idx_n, drop = FALSE] %*% w_n
    }
  }
  U <- (1 - mix) * U_raw + mix * U_nbr
  sweep(U, 2L, sqrt(colMeans(U^2)), "/")
}

events     <- vector("list", n_runs)
covariates <- vector("list", n_runs)
names(events) <- names(covariates) <- paste0("run", seq_len(n_runs))
run_files <- character(n_runs)
true_neural_runs <- vector("list", n_runs)

active_zones <- setdiff(unique(vinfo$zone), c("Null background", "Opposing sulcal null bank", "Focal parcel: silent surround"))

for (r in seq_len(n_runs)) {
  cue_onsets <- round(cumsum(c(6, runif(n_trials - 1L, 6.4, 7.8))), 1)
  isi_cf     <- round(runif(n_trials, 2.2, 3.4), 1)
  fb_onsets  <- cue_onsets + isi_cf
  rpe_vec    <- as.numeric(scale(rnorm(n_trials)))

  ev_r <- rbind(
    data.frame(onset = cue_onsets, duration = 0, condition = "cue"),
    data.frame(onset = fb_onsets,  duration = 0, condition = "feedback")
  )
  ev_r <- ev_r[order(ev_r$onset), ]
  rownames(ev_r) <- NULL
  events[[r]]     <- ev_r
  covariates[[r]] <- data.frame(rpe = rpe_vec)

  # Continuous neural grid (n_time x V)
  S_neural <- matrix(0, nrow = n_time, ncol = V)
  cue_idx_t <- pmin(n_time, pmax(1L, round(cue_onsets / tr) + 1L))
  fb_idx_t  <- pmin(n_time, pmax(1L, round(fb_onsets  / tr) + 1L))

  # Cue impulses (mild trial-to-trial variability)
  active_cue <- which(vinfo$amp_cue > 0)
  for (ti in seq_len(n_trials)) {
    S_neural[cue_idx_t[ti], active_cue] <- S_neural[cue_idx_t[ti], active_cue] +
      vinfo$amp_cue[active_cue] * (1 + rnorm(length(active_cue), 0, 0.15))
  }

  # Feedback impulses: mean + beta_rpe * rpe + biological trial shock (partially shared within zone)
  for (zn in active_zones) {
    v_zn <- which(vinfo$zone == zn)
    z_shared <- rnorm(n_trials)
    for (ti in seq_len(n_trials)) {
      bio_shock <- vinfo$tau_bio[v_zn] * (0.65 * z_shared[ti] + sqrt(1 - 0.65^2) * rnorm(length(v_zn)))
      S_neural[fb_idx_t[ti], v_zn] <- S_neural[fb_idx_t[ti], v_zn] +
        vinfo$amp_fb[v_zn] + vinfo$beta_rpe[v_zn] * rpe_vec[ti] + bio_shock
    }
  }
  true_neural_runs[[r]] <- S_neural

  # Convolve each HRF stratum with its true HRF kernel
  Y_clean <- matrix(0, nrow = n_time, ncol = V)
  for (h_id in unique(vinfo$hrf_idx)) {
    v_h <- which(vinfo$hrf_idx == h_id)
    k_vec <- K_lib[, h_id]
    L_k <- length(k_vec)
    for (lag_i in seq_len(L_k)) {
      if (lag_i == 1L) {
        Y_clean[, v_h] <- Y_clean[, v_h] + k_vec[lag_i] * S_neural[, v_h, drop = FALSE]
      } else {
        rows_dst <- lag_i:n_time
        rows_src <- 1L:(n_time - lag_i + 1L)
        Y_clean[rows_dst, v_h] <- Y_clean[rows_dst, v_h] + k_vec[lag_i] * S_neural[rows_src, v_h, drop = FALSE]
      }
    }
  }

  # Spatially correlated AR(1) observation noise with pial-rim vs GM heterogeneity
  U_corr <- smooth_innovations(matrix(rnorm(n_time * V), n_time, V), g_noise_gen, mix = 0.55)
  E_noise <- matrix(0, nrow = n_time, ncol = V)
  rho_v   <- vinfo$rho_true
  sig_v   <- vinfo$sigma_true
  inno_sc <- sig_v * sqrt(pmax(1e-6, 1 - rho_v^2))
  E_noise[1L, ] <- sig_v * U_corr[1L, ]
  for (tt in 2L:n_time) {
    E_noise[tt, ] <- rho_v * E_noise[tt - 1L, ] + inno_sc * U_corr[tt, ]
  }

  Y_obs <- 100 + Y_clean + E_noise

  # Write 4D NIfTI run
  arr4d <- array(0, dim = c(dim(t1_nii), n_time))
  vol3  <- array(0, dim = dim(t1_nii))
  for (tt in seq_len(n_time)) {
    vol3[lin_idx] <- Y_obs[tt, ]
    arr4d[, , , tt] <- vol3
  }
  hdr4 <- updateNifti(arr4d, template = t1_nii)
  hdr4$pixdim[5] <- tr
  f_r <- file.path(tmp_dir, sprintf("sub-01_run-%02d_bold.nii.gz", r))
  writeNifti(hdr4, f_r)
  run_files[r] <- f_r
}

# ------------------------------------------------------------------------------
# 4. Ingest via read_nifti_ts() & build spatial graphs (Geodesic vs Euclidean)
# ------------------------------------------------------------------------------
ts <- read_nifti_ts(
  niftis    = run_files,
  mask      = mask_file,
  atlas     = atlas_file,
  tr        = tr,
  run_names = names(events)
)

# Align vinfo strictly to ts$voxel_map order
stopifnot(nrow(ts$voxel_map) == V)
match_idx <- match(ts$voxel_map$linear_index, lin_idx)
stopifnot(!anyNA(match_idx), all(match_idx == seq_len(V)))
vinfo[, voxel := ts$voxel_map$voxel]

g_geo <- build_spatial_graph(ts, spatial_fwhm = 6, parcel_boundary = "soft", metric = "geodesic")
g_euc <- build_spatial_graph(ts, spatial_fwhm = 6, parcel_boundary = "soft", metric = "euclidean")

# Extract cross-sulcal / cross-fissure Euclidean shortcut edges on slice z = 8 mm (k = 44)
shortcut_list <- list()
total_euc_shortcuts <- 0L
for (u in seq_len(V)) {
  se <- g_euc$row_ptr[u]; ee <- g_euc$row_ptr[u + 1L] - 1L
  sg <- g_geo$row_ptr[u]; eg <- g_geo$row_ptr[u + 1L] - 1L
  ne <- if (ee >= se) g_euc$col_idx[se:ee] else integer(0)
  ng <- if (eg >= sg) g_geo$col_idx[sg:eg] else integer(0)
  only_euc <- setdiff(ne, ng)
  total_euc_shortcuts <- total_euc_shortcuts + length(only_euc)
  if (length(only_euc) > 0L && vinfo$k[u] == 44L) {
    in_k44 <- only_euc[vinfo$k[only_euc] == 44L & only_euc > u]
    if (length(in_k44) > 0L) {
      shortcut_list[[length(shortcut_list) + 1L]] <- data.table(
        u = u, v = in_k44,
        x1 = vinfo$x[u], y1 = vinfo$y[u],
        x2 = vinfo$x[in_k44], y2 = vinfo$y[in_k44],
        parcel_u = vinfo$parcel[u], parcel_v = vinfo$parcel[in_k44]
      )
    }
  }
}
shortcuts_z8 <- rbindlist(shortcut_list)
cat(sprintf("Euclidean 6mm graph added %d cross-sulcal/out-of-ribbon directed edges (%d in-plane pairs on z=8mm).\n",
            total_euc_shortcuts, nrow(shortcuts_z8)))

# Also extract local geodesic neighbor edges for two exemplar seed voxels on z = 8 mm
get_ego_edges <- function(g, seed_u) {
  s <- g$row_ptr[seed_u]; e <- g$row_ptr[seed_u + 1L] - 1L
  nbrs <- g$col_idx[s:e]
  nbrs <- nbrs[vinfo$k[nbrs] == 44L]
  data.table(
    seed = seed_u, nbr = nbrs,
    x1 = vinfo$x[seed_u], y1 = vinfo$y[seed_u],
    x2 = vinfo$x[nbrs],   y2 = vinfo$y[nbrs]
  )
}
# Seed on medial edge of Left Calcarine (parcel 10) right next to interhemispheric fissure (x >= -6, z = 8)
seed_calc <- vinfo[parcel == 10L & k == 44L][which.max(x), vnum]
ego_geo_calc <- get_ego_edges(g_geo, seed_calc)
ego_euc_calc <- get_ego_edges(g_euc, seed_calc)

# ------------------------------------------------------------------------------
# 5. Step 1: Compare HRF & Noise Pooling Modes
# ------------------------------------------------------------------------------
cat("Fitting Step 1 HRF & noise models...\n")

# 5a. Pure single-voxel pooling (atlas-free)
hrf_voxel <- estimate_crossrun_hrf(
  bold            = ts$bold,
  events          = events,
  tr              = tr,
  kernel_grid     = K_lib,
  pooling         = "voxel",
  amplitude_model = "trial",
  selection_mode  = "soft_average"
)

# 5b. Default atlas-backed voxel mode (parcel noise pooling + parcel fallback for low-SNR voxels)
hrf_voxel_atlas <- estimate_crossrun_hrf(
  bold            = ts,
  events          = events,
  kernel_grid     = K_lib,
  pooling         = "voxel",
  amplitude_model = "trial",
  selection_mode  = "soft_average"
)

# 5c. Rigid group-parcel pooling
hrf_parcel <- estimate_crossrun_hrf(
  bold            = ts,
  events          = events,
  kernel_grid     = K_lib,
  pooling         = "parcel",
  amplitude_model = "trial",
  selection_mode  = "soft_average"
)

# 5d. Sparse bilateral geodesic graph pooling with soft parcel boundary (recommended)
hrf_graph_soft <- estimate_crossrun_hrf(
  bold            = ts,
  events          = events,
  kernel_grid     = K_lib,
  pooling         = "graph",
  spatial_fwhm    = 6,
  parcel_boundary = "soft",
  amplitude_model = "trial",
  selection_mode  = "soft_average"
)

# 5e. Sparse bilateral geodesic graph pooling with hard parcel wall
hrf_graph_hard <- estimate_crossrun_hrf(
  bold            = ts,
  events          = events,
  kernel_grid     = K_lib,
  pooling         = "graph",
  spatial_fwhm    = 6,
  parcel_boundary = "hard",
  amplitude_model = "trial",
  selection_mode  = "soft_average"
)

# 5f. Euclidean 3D graph pooling ignoring sulcal/parcel barriers
g_euc_ign <- build_spatial_graph(ts, spatial_fwhm = 6, parcel_boundary = "ignore", metric = "euclidean")
hrf_graph_euc <- estimate_crossrun_hrf(
  bold            = ts,
  events          = events,
  kernel_grid     = K_lib,
  pooling         = "graph",
  spatial_graph   = g_euc_ign,
  parcel_boundary = "ignore",
  amplitude_model = "trial",
  selection_mode  = "soft_average"
)

# Compute continuous interpolated peak time (parabolic sub-TR refinement) for smooth evaluation
est_peak_sec <- function(fit) {
  # Average selected kernels across the 3 cross-fitted runs per voxel
  K_avg <- (fit$selected_kernels[[1L]] + fit$selected_kernels[[2L]] + fit$selected_kernels[[3L]]) / 3
  n_k <- nrow(K_avg)
  pk <- numeric(ncol(K_avg))
  for (v in seq_len(ncol(K_avg))) {
    col_v <- K_avg[, v]
    m <- which.max(col_v)
    if (m > 1L && m < n_k) {
      denom <- col_v[m - 1L] - 2 * col_v[m] + col_v[m + 1L]
      delta <- if (abs(denom) > 1e-10) 0.5 * (col_v[m - 1L] - col_v[m + 1L]) / denom else 0
      pk[v] <- (m - 1L + max(-0.5, min(0.5, delta))) * tr
    } else {
      pk[v] <- (m - 1L) * tr
    }
  }
  pk
}

true_pk_cont <- est_peak_sec(list(selected_kernels = replicate(3L, K_lib[, vinfo$hrf_idx, drop = FALSE], simplify = FALSE)))
vinfo[, true_peak_cont := true_pk_cont]

hrf_peaks <- data.table(
  voxel       = vinfo$voxel,
  zone        = vinfo$zone,
  parcel      = vinfo$parcel,
  y           = vinfo$y,
  truth       = true_pk_cont,
  voxel_only  = est_peak_sec(hrf_voxel),
  voxel_atlas = est_peak_sec(hrf_voxel_atlas),
  parcel_pool = est_peak_sec(hrf_parcel),
  graph_soft  = est_peak_sec(hrf_graph_soft),
  graph_hard  = est_peak_sec(hrf_graph_hard),
  graph_euc   = est_peak_sec(hrf_graph_euc)
)

# Noise recovery (average across runs)
avg_noise <- function(fit, col) {
  (fit$noise[[1L]][[col]] + fit$noise[[2L]][[col]] + fit$noise[[3L]][[col]]) / 3
}
noise_dt <- data.table(
  voxel       = vinfo$voxel,
  is_pial_rim = vinfo$is_pial_rim,
  rho_true    = vinfo$rho_true,
  sigma_true  = vinfo$sigma_true,
  rho_voxel   = avg_noise(hrf_voxel, "rho"),
  rho_parcel  = avg_noise(hrf_parcel, "rho"),
  rho_graph   = avg_noise(hrf_graph_soft, "rho"),
  sig_voxel   = avg_noise(hrf_voxel, "noise_sd"),
  sig_parcel  = avg_noise(hrf_parcel, "noise_sd"),
  sig_graph   = avg_noise(hrf_graph_soft, "noise_sd")
)

# Summarize Step 1 HRF & noise accuracy
active_mask <- !(vinfo$zone %in% c("Null background", "Opposing sulcal null bank", "Focal parcel: silent surround"))
shifted_p78 <- vinfo$parcel == 78L
opposing_mid <- vinfo$parcel %in% c(109L, 194L) &
                vinfo$vnum %in% g_euc$col_idx # medial rim within 6mm Euclidean of active L bank
# Identify exact opposing-bank voxels that have Euclidean neighbors in parcel 10 or 89
opp_bleed_vox <- rep(FALSE, V)
for (u in which(vinfo$parcel %in% c(109L, 194L))) {
  se <- g_euc$row_ptr[u]; ee <- g_euc$row_ptr[u + 1L] - 1L
  if (ee >= se && any(vinfo$parcel[g_euc$col_idx[se:ee]] %in% c(10L, 89L))) {
    opp_bleed_vox[u] <- TRUE
  }
}
vinfo[, is_opposing_rim := opp_bleed_vox]

step1_summary <- rbindlist(lapply(
  list(
    list(id = "voxel_only",  label = "Voxel (no pooling)",            fit = hrf_voxel,       rho = noise_dt$rho_voxel,  sig = noise_dt$sig_voxel),
    list(id = "voxel_atlas", label = "Voxel + atlas fallback",        fit = hrf_voxel_atlas, rho = avg_noise(hrf_voxel_atlas, "rho"), sig = avg_noise(hrf_voxel_atlas, "noise_sd")),
    list(id = "parcel_pool", label = "Parcel pooling",                fit = hrf_parcel,      rho = noise_dt$rho_parcel, sig = noise_dt$sig_parcel),
    list(id = "graph_hard",  label = "Geodesic graph (hard parcel)",  fit = hrf_graph_hard,  rho = avg_noise(hrf_graph_hard, "rho"),  sig = avg_noise(hrf_graph_hard, "noise_sd")),
    list(id = "graph_soft",  label = "Geodesic graph (soft parcel)",  fit = hrf_graph_soft,  rho = noise_dt$rho_graph,  sig = noise_dt$sig_graph),
    list(id = "graph_euc",   label = "Euclidean 3D graph (no mask barrier)", fit = hrf_graph_euc, rho = avg_noise(hrf_graph_euc, "rho"), sig = avg_noise(hrf_graph_euc, "noise_sd"))
  ),
  function(m) {
    pk_est <- hrf_peaks[[m$id]]
    n_st <- if (!is.null(m$fit$hrf_strata[[1L]]) && length(m$fit$hrf_strata[[1L]]) > 0L) {
      length(unique(m$fit$hrf_strata[[1L]]))
    } else {
      nrow(unique(t(round(m$fit$selected_kernels[[1L]], 6))))
    }
    data.table(
      method              = m$label,
      mae_active_sec      = mean(abs(pk_est[active_mask] - hrf_peaks$truth[active_mask])),
      mae_shifted_p78_sec = mean(abs(pk_est[shifted_p78] - hrf_peaks$truth[shifted_p78])),
      mae_opposing_rim    = mean(abs(pk_est[opp_bleed_vox] - hrf_peaks$truth[opp_bleed_vox])),
      rmse_rho_rim        = sqrt(mean((m$rho[vinfo$is_pial_rim] - vinfo$rho_true[vinfo$is_pial_rim])^2)),
      rmse_rho_gm         = sqrt(mean((m$rho[!vinfo$is_pial_rim] - vinfo$rho_true[!vinfo$is_pial_rim])^2)),
      rmse_sigma_all      = sqrt(mean((m$sig - vinfo$sigma_true)^2)),
      n_hrf_strata_run1   = as.integer(n_st)
    )
  }
))
cat("\nStep 1 Summary:\n")
print(step1_summary)

# ------------------------------------------------------------------------------
# 6. Step 2: SNR-Stratified vs Global Deconvolution Regularization
# ------------------------------------------------------------------------------
cat("\nFitting Step 2 continuous deconvolution (lambda_mode = 'auto' vs 'global')...\n")

dec_auto <- deconvolve_crossrun_hrf(
  bold        = ts,
  hrf_fit     = hrf_graph_soft,
  scaling     = "noise_sd",
  lambda_mode = "auto"
)

dec_global <- deconvolve_crossrun_hrf(
  bold        = ts,
  hrf_fit     = hrf_graph_soft,
  scaling     = "noise_sd",
  lambda_mode = "global"
)

lambda_vox_auto <- dec_auto$scaling[[1L]]$lambda
step2_strata <- as.data.table(dec_auto$settings$lambda_strata)
step2_strata[, lambda_global := dec_global$settings$lambda[1L]]

# Measure continuous background SD in null voxels & correlation with true neural signal in active/hotspot voxels
null_vox <- which(vinfo$zone == "Null background")
hot_vox  <- which(vinfo$zone == "Focal hotspot (R Insula core)")
act_vox  <- which(active_mask)

bg_sd_auto   <- mean(apply(dec_auto$activity[[1L]][, null_vox, drop = FALSE], 2L, sd))
bg_sd_global <- mean(apply(dec_global$activity[[1L]][, null_vox, drop = FALSE], 2L, sd))

cor_hot_auto   <- mean(vapply(hot_vox, function(v) cor(dec_auto$activity[[1L]][, v],   true_neural_runs[[1L]][, v]), numeric(1)))
cor_hot_global <- mean(vapply(hot_vox, function(v) cor(dec_global$activity[[1L]][, v], true_neural_runs[[1L]][, v]), numeric(1)))

step2_comparison <- data.table(
  mode               = c("Stratified (lambda_mode = 'auto')", "Global (lambda_mode = 'global')"),
  n_lambda_tiers     = c(nrow(step2_strata), 1L),
  lambda_range       = c(sprintf("%.3f - %.3f", min(lambda_vox_auto), max(lambda_vox_auto)),
                         sprintf("%.3f", dec_global$settings$lambda[1L])),
  null_bg_sd         = c(bg_sd_auto, bg_sd_global),
  focal_hotspot_corr = c(cor_hot_auto, cor_hot_global)
)
cat("\nStep 2 Lambda Strata:\n")
print(step2_strata)
print(step2_comparison)

# Save representative run1 time series for 1 focal hotspot voxel & 1 quiet Tier-1 null voxel
v_hot_rep  <- hot_vox[which.max(colSums(true_neural_runs[[1L]][, hot_vox, drop = FALSE]^2))]
v_null_rep <- null_vox[which(lambda_vox_auto[null_vox] > 3.0)[1L]]
step2_traces <- data.table(
  time          = (seq_len(n_time) - 1L) * tr,
  hot_true      = true_neural_runs[[1L]][, v_hot_rep],
  hot_auto      = dec_auto$activity[[1L]][, v_hot_rep],
  hot_global    = dec_global$activity[[1L]][, v_hot_rep],
  null_auto     = dec_auto$activity[[1L]][, v_null_rep],
  null_global   = dec_global$activity[[1L]][, v_null_rep]
)

# ------------------------------------------------------------------------------
# 7. Step 3 (Voxelwise) vs Step 4 (Parcel-Wise) Trial-Epoch Regression
# ------------------------------------------------------------------------------
cat("\nFitting Step 3 voxelwise regressions (graph_eb vs eb vs none) & Step 4 parcel panel...\n")
lags_eval <- -1:3

ep_vox_grapheb <- decon_epoch_regression(
  deconvolution    = dec_auto,
  covariates       = covariates,
  events           = events,
  target_condition = "feedback",
  lags             = lags_eval,
  restore_scale    = TRUE,
  tau2_mode        = "residual",
  tau2_shrinkage   = "graph_eb"
)

ep_vox_eb <- decon_epoch_regression(
  deconvolution    = dec_auto,
  covariates       = covariates,
  events           = events,
  target_condition = "feedback",
  lags             = 0,
  restore_scale    = TRUE,
  tau2_mode        = "residual",
  tau2_shrinkage   = "eb"
)

ep_vox_none <- decon_epoch_regression(
  deconvolution    = dec_auto,
  covariates       = covariates,
  events           = events,
  target_condition = "feedback",
  lags             = 0,
  restore_scale    = TRUE,
  tau2_mode        = "residual",
  tau2_shrinkage   = "none"
)

# Also run voxelwise with Euclidean graph to show cross-fissure tau2 / scale bleed
ep_vox_euc <- decon_epoch_regression(
  deconvolution    = dec_auto,
  covariates       = covariates,
  events           = events,
  target_condition = "feedback",
  lags             = 0,
  restore_scale    = TRUE,
  tau2_mode        = "residual",
  tau2_shrinkage   = "graph_eb",
  spatial_graph    = g_euc_ign
)

# Pathway B: Step 4 Parcel-Wise Aggregation (exact multi-stratum operator + noise_innovations)
ep_parcel <- decon_epoch_regression(
  deconvolution    = dec_auto,
  covariates       = covariates,
  events           = events,
  target_condition = "feedback",
  lags             = lags_eval,
  parcel_groups    = ts$parcel_ids,
  restore_scale    = TRUE,
  tau2_mode        = "residual",
  whiten_panel     = TRUE,
  return_panel     = TRUE
)

# Extract lag = 0 s rpe slopes & tau2_bio across all voxelwise & parcel-wise methods
extract_lag0 <- function(ep) {
  d <- as.data.table(ep$coefficients)[lag == 0 & term == "rpe"]
  d[match(vinfo$voxel, voxel)]
}
c_grapheb <- extract_lag0(ep_vox_grapheb)
c_eb      <- extract_lag0(ep_vox_eb)
c_none    <- extract_lag0(ep_vox_none)
c_euc     <- extract_lag0(ep_vox_euc)

c_parcel_lag0 <- as.data.table(ep_parcel$coefficients)[lag == 0 & term == "rpe"]
c_parcel_lag0[, parcel_int := as.integer(parcel)]
p_match <- match(vinfo$parcel, c_parcel_lag0$parcel_int)

vox_lag0_maps <- data.table(
  voxel          = vinfo$voxel,
  parcel         = vinfo$parcel,
  zone           = vinfo$zone,
  is_pial_rim    = vinfo$is_pial_rim,
  is_opposing    = vinfo$is_opposing_rim,
  beta_true      = vinfo$beta_rpe,
  tau2_true      = vinfo$true_tau2_bio,
  lambda_auto    = lambda_vox_auto,
  # Voxelwise graph_eb
  beta_grapheb   = c_grapheb$beta_gls,
  t_grapheb      = c_grapheb$t_gls,
  se_grapheb     = c_grapheb$se_gls,
  tau2_grapheb   = c_grapheb$tau2_bio,
  # Voxelwise global eb
  beta_eb        = c_eb$beta_gls,
  t_eb           = c_eb$t_gls,
  tau2_eb        = c_eb$tau2_bio,
  # Voxelwise none (raw MoM)
  beta_none      = c_none$beta_gls,
  t_none         = c_none$t_gls,
  se_none        = c_none$se_gls,
  tau2_none      = c_none$tau2_bio,
  # Voxelwise Euclidean graph_eb
  tau2_euc       = c_euc$tau2_bio,
  t_euc          = c_euc$t_gls,
  # Parcel-wise projected onto constituent voxels
  beta_parcel    = c_parcel_lag0$beta_gls[p_match],
  t_parcel       = c_parcel_lag0$t_gls[p_match],
  se_parcel      = c_parcel_lag0$se_gls[p_match],
  tau2_parcel    = c_parcel_lag0$tau2_bio[p_match]
)

# Summarize Step 3 tau2 shrinkage metrics
step3_tau2_summary <- rbindlist(lapply(
  list(
    list(name = "Raw voxel MoM ('none')",            tau2 = vox_lag0_maps$tau2_none,    t = vox_lag0_maps$t_none,    beta = vox_lag0_maps$beta_none),
    list(name = "Global Empirical Bayes ('eb')",     tau2 = vox_lag0_maps$tau2_eb,      t = vox_lag0_maps$t_eb,      beta = vox_lag0_maps$beta_eb),
    list(name = "Euclidean Graph-EB",                tau2 = vox_lag0_maps$tau2_euc,     t = vox_lag0_maps$t_euc,     beta = c_euc$beta_gls),
    list(name = "Bilateral Geodesic Graph-EB ('graph_eb')", tau2 = vox_lag0_maps$tau2_grapheb, t = vox_lag0_maps$t_grapheb, beta = vox_lag0_maps$beta_grapheb)
  ),
  function(m) {
    data.table(
      method                 = m$name,
      active_zero_collapse   = mean(m$tau2[active_mask] < 1e-6),
      null_mean_tau2         = mean(m$tau2[null_vox]),
      opposing_rim_tau2      = mean(m$tau2[opp_bleed_vox]),
      hotspot_mean_tau2      = mean(m$tau2[hot_vox]),
      rmse_tau2_all          = sqrt(mean((m$tau2 - vox_lag0_maps$tau2_true)^2)),
      null_false_pos_p05     = mean(2 * pt(-abs(m$t[null_vox]), df = c_grapheb$df[null_vox]) <= 0.05),
      power_active_fdr05     = {
        p_all <- 2 * pt(-abs(m$t), df = c_grapheb$df)
        q_all <- p.adjust(p_all, method = "BH")
        mean(q_all[active_mask] <= 0.05 & sign(m$beta[active_mask]) == sign(vox_lag0_maps$beta_true[active_mask]))
      }
    )
  }
))
cat("\nStep 3 Voxelwise Tau2 Summary:\n")
print(step3_tau2_summary)

# Summarize Voxelwise (Graph-EB) vs Parcel-Wise (Step 4) across the 4 anatomical regimes
regime_defs <- list(
  "Aligned parcel (L Calcarine & L mPFC)" = which(vinfo$parcel %in% c(10L, 89L)),
  "Shifted border parcel (L Temp 78)"     = which(vinfo$parcel == 78L),
  "Sub-parcel hotspot (R Insula core)"    = hot_vox,
  "Null cortex (55 background parcels)"   = null_vox
)

vox_vs_parcel_summary <- rbindlist(lapply(names(regime_defs), function(rn) {
  idx_r <- regime_defs[[rn]]
  is_null_reg <- grepl("^Null", rn)
  p_vox <- 2 * pt(-abs(vox_lag0_maps$t_grapheb), df = c_grapheb$df)
  q_vox <- p.adjust(p_vox, method = "BH")
  p_par_unique <- c_parcel_lag0$p_value
  q_par_unique <- p.adjust(p_par_unique, method = "BH")
  q_par_vox    <- q_par_unique[p_match]

  det_vox <- if (is_null_reg) {
    mean(q_vox[idx_r] <= 0.05)
  } else {
    mean(q_vox[idx_r] <= 0.05 & sign(vox_lag0_maps$beta_grapheb[idx_r]) == sign(vox_lag0_maps$beta_true[idx_r]))
  }
  det_par <- if (is_null_reg) {
    mean(q_par_vox[idx_r] <= 0.05)
  } else {
    mean(q_par_vox[idx_r] <= 0.05 & sign(vox_lag0_maps$beta_parcel[idx_r]) == sign(vox_lag0_maps$beta_true[idx_r]))
  }

  rbind(
    data.table(
      regime       = rn,
      pathway      = "Voxelwise (Graph + Stratified + Graph-EB)",
      rmse_beta    = sqrt(mean((vox_lag0_maps$beta_grapheb[idx_r] - vox_lag0_maps$beta_true[idx_r])^2)),
      mean_abs_t   = mean(abs(vox_lag0_maps$t_grapheb[idx_r])),
      detect_fdr05 = det_vox
    ),
    data.table(
      regime       = rn,
      pathway      = "Parcel-wise (Step 4 multi-stratum panel)",
      rmse_beta    = sqrt(mean((vox_lag0_maps$beta_parcel[idx_r] - vox_lag0_maps$beta_true[idx_r])^2)),
      mean_abs_t   = mean(abs(vox_lag0_maps$t_parcel[idx_r])),
      detect_fdr05 = det_par
    )
  )
}))
cat("\nVoxelwise vs Parcel-wise Trade-Off Summary:\n")
print(vox_vs_parcel_summary)

# Verify write_nifti_ts() (voxelwise) and backproject_medusa() (parcel-wise) end-to-end
vox_out_nii <- file.path(tmp_dir, "vox_t_grapheb_lag0.nii.gz")
vox_t_vec   <- stats::setNames(vox_lag0_maps$t_grapheb, vinfo$voxel)
write_nifti_ts(vox_t_vec, vox_out_nii, reference = dec_auto, fill = NA_real_)
stopifnot(file.exists(vox_out_nii))

par_bp_df <- subset(ep_parcel$coefficients, lag == 0 & term == "rpe")
bp_files <- backproject_medusa(
  coef_df    = par_bp_df,
  brain_mask = dec_auto,
  plot_cols  = c("beta_gls", "t_gls"),
  output_dir = tmp_dir,
  ncpus      = 1L
)
stopifnot(all(file.exists(bp_files)))

# ------------------------------------------------------------------------------
# 8. Step 5 Cohort Demonstration: 12 Participants -> mixed_by() -> backproject_medusa()
# ------------------------------------------------------------------------------
cat("\nSimulating 12-participant parcel panel cohort for Step 5 mixed_by()...\n")
n_subj <- 12L
parcel_truth_dt <- vinfo[, .(
  beta_true_parcel = mean(beta_rpe),
  n_vox_parcel     = .N
), by = parcel]

subj_sds <- seq(0.16, 0.42, length.out = n_subj) # 2.6x motion/SNR spread across cohort
cohort_list <- vector("list", n_subj)
for (s in seq_len(n_subj)) {
  sid <- sprintf("sub-%02d", s)
  b_rand <- rnorm(nrow(parcel_truth_dt), sd = 0.10)
  dt_s <- CJ(
    run   = c("run1", "run2", "run3"),
    trial = seq_len(n_trials),
    parcel = as.character(parcel_truth_dt$parcel)
  )
  # Generate trial RPE per run
  rpe_s <- dt_s[, .(rpe = as.numeric(scale(rnorm(n_trials))), trial = seq_len(n_trials)), by = run]
  dt_s <- merge(dt_s, rpe_s, by = c("run", "trial"))
  p_idx_s <- match(as.integer(dt_s$parcel), parcel_truth_dt$parcel)
  n_v_p   <- parcel_truth_dt$n_vox_parcel[p_idx_s]
  # Effective parcel measurement SE scales with subject noise SD and parcel size (accounting for spatial corr)
  se_meas <- subj_sds[s] * 4.5 / sqrt(1 + 0.22 * n_v_p)
  tau_lat <- ifelse(abs(parcel_truth_dt$beta_true_parcel[p_idx_s]) > 0.05, 0.28, 0.05)
  b_subj  <- parcel_truth_dt$beta_true_parcel[p_idx_s] + b_rand[p_idx_s]
  dt_s[, `:=`(
    id                      = sid,
    event_index             = 2L * trial,
    evt_time                = 0,
    se_demixed              = se_meas,
    precision_weight_cohort = 1 / (tau_lat^2 + se_meas^2),
    decon_demixed           = b_subj * rpe + rnorm(.N, sd = sqrt(tau_lat^2 + se_meas^2))
  )]
  cohort_list[[s]] <- dt_s
}
group_panel <- rbindlist(cohort_list)

group_fit <- mixed_by(
  data               = group_panel,
  outcomes           = "decon_demixed",
  rhs_model_formulae = ~ rpe + (1 + rpe || id),
  split_on           = c("parcel", "evt_time"),
  weights            = "precision_weight_cohort",
  calculate          = "parameter_estimates_reml",
  padjust_by         = "term",
  padjust_method     = "BY",
  ncores             = 1L
)

group_rpe <- as.data.table(group_fit$coef_df_reml)[term == "rpe"]
group_rpe[, parcel_int := as.integer(parcel)]
group_rpe <- merge(group_rpe, parcel_truth_dt, by.x = "parcel_int", by.y = "parcel")

# ------------------------------------------------------------------------------
# 9. Save compact reproducible artifacts for figures & Quarto chapter
# ------------------------------------------------------------------------------
fwrite(step1_summary,         "step1_hrf_noise_metrics.csv")
fwrite(step2_strata,          "step2_lambda_strata.csv")
fwrite(step2_comparison,      "step2_lambda_comparison.csv")
fwrite(step3_tau2_summary,    "step3_tau2_voxel_metrics.csv")
fwrite(vox_vs_parcel_summary, "voxel_vs_parcel_comparison.csv")

parcel_lag_profiles <- as.data.table(ep_parcel$coefficients)[
  term == "rpe" & as.integer(parcel) %in% c(10L, 78L, 79L, 89L, 109L, 155L),
  .(parcel = as.integer(parcel), lag, beta_gls, se_gls, t_gls, p_value, tau2_bio)
]

saveRDS(
  list(
    slab_k              = slab_k,
    tr                  = tr,
    n_runs              = n_runs,
    n_time              = n_time,
    n_trials            = n_trials,
    voxel_map           = ts$voxel_map,
    geometry            = ts$geometry,
    vinfo               = vinfo,
    shortcuts_z8        = shortcuts_z8,
    total_euc_shortcuts = total_euc_shortcuts,
    seed_calc           = seed_calc,
    ego_geo_calc        = ego_geo_calc,
    ego_euc_calc        = ego_euc_calc,
    hrf_peaks           = hrf_peaks,
    noise_dt            = noise_dt,
    step1_summary       = step1_summary,
    step2_strata        = step2_strata,
    step2_comparison    = step2_comparison,
    step2_traces        = step2_traces,
    vox_lag0_maps       = vox_lag0_maps,
    c_parcel_lag0       = c_parcel_lag0,
    parcel_lag_profiles = parcel_lag_profiles,
    step3_tau2_summary  = step3_tau2_summary,
    vox_vs_parcel       = vox_vs_parcel_summary,
    group_rpe           = group_rpe
  ),
  file = "spatial_pooling_examples.rds"
)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat(sprintf("\nCompleted spatial & parcel pooling simulation in %.1f s.\n",
            proc.time()[3] - t_start))
