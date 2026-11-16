#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importFrom RNifti readNifti
#' @importFrom data.table := data.table is.data.table setorderv
#' @importFrom dplyr %>%
#' @importFrom glue glue
#' @importFrom rlang .data
#' @importFrom stats approx as.formula na.omit
#' @importFrom tidyselect everything
#' @importFrom utils head tail
#' @useDynLib medusa, .registration = TRUE
## usethis namespace: end
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ":=", ".", "..split_on", "..time_column",
    ".align_censored", ".align_epoch_id", ".align_source_id",
    ".column", ".epoch", ".epochs", ".estimate", ".filename", ".gid",
    ".GRP", ".I", ".left", ".member", ".N", ".onset",
    ".result_id", ".right", ".s2_bar", ".SD", ".source", ".support", ".support_status",
    ".parcel_tmp", ".u_ord", ".unit", ".usage", ".w_norm", ".w_raw",
    "atlas", "boundary", "coef_df_ml", "coef_df_reml", "decon_demixed", "dispersion",
    "dt", "dt_split", "effect", "emm", "emp", "emt", "epoch_status", "error_ac1",
    "estimate", "etp", "evt_onset", "evt_time",
    "fit_df", "fitted", "fname", "i", "id", "ii", "isi", "iti_post",
    "lag_noise_vec", "model", "model_name", "name", "outcome", "p.value", "peak_delay",
    "psi_seed_proj", "psi_t0_proj", "rel_seed_trial", "rel_t0_trial",
    "residuals", "rhs", "run",
    "se_demixed", "se_raw", "se_seed", "se_targ", "se_targ_t0",
    "seed_c", "seed_c_raw", "seed_eb_c", "seed_val",
    "shared_noise_loading", "shared_noise_vec",
    "si", "split", "std.error",
    "targ_c", "targ_t0_c", "targ_t0_eb_c", "targ_val", "targ_t0_val",
    "term", "time", "trial", "trial_cor_vec", "undershoot_ratio", "unit",
    "value", "value_source_max_epochs", "vox_ts", "voxel",
    "x", "xm", "y", "ym", "z"
  ))
}
