source_dir <- Sys.getenv("MEDUSA_SOURCE")
if (nzchar(source_dir)) pkgload::load_all(source_dir, compile = TRUE, quiet = TRUE)
knitr::purl("workflow-source.Rmd", output = "workflow-script.R", quiet = TRUE)
source("workflow-script.R")
focal <- subset(group_result$coef_df_reml, term == "rpe_within")
curve <- subset(trajectory$coefficients, parcel == "Parcel_A" & term == "rpe")
stopifnot(nrow(focal) == 1L, nrow(curve) == length(lags_group), nrow(planned) == 720L)
fit_ppi <- function(shared) decon_connectivity_repair(
  group_panel, seed_region = "Parcel_A", target_region = "Parcel_B",
  covariates = "rpe", test_ppi = TRUE, use_whitened = FALSE,
  orthogonalize_seed_x = FALSE, shared_error = shared, method = "mmoms",
  region_col = "parcel"
)
ppi_auto <- fit_ppi("auto")
ppi_zero <- fit_ppi("zero")
out <- list(focal = focal, curve = curve, subjects = length(participant_ids),
            n_trials = nrow(planned), singular = lme4::isSingular(group_fit),
            ppi = rbind(transform(as.data.frame(ppi_auto$group), shared = "auto"),
                        transform(as.data.frame(ppi_zero$group), shared = "zero")),
            diagnostic = subset(ppi_auto$subjects, term == "seed_c:x_c"),
            software = capture.output(sessionInfo()))
saveRDS(out, "example.rds")
print(focal)
print(curve[, .(lag, estimate, p_maxT, n_subjects)])
print(out$ppi[, c("shared", "term", "estimate", "p.value")])
