d <- readRDS('step1_examples.rds')
checks <- list()
for (fit_name in c('soft','hard','voxel','perturbed','signed_condition','signed_trial')) {
 fit <- d[[fit_name]]
 for (r in names(fit$kernel_weights)) {
  w <- fit$kernel_weights[[r]]
  stopifnot(all(is.finite(w)),all(w>=0),max(abs(colSums(w)-1))<1e-12)
  stopifnot(max(abs(d$K%*%w-fit$selected_kernels[[r]]))<1e-12)
 }
 checks[[fit_name]] <- 'Nonnegative unit-sum weights; selected kernels match library-weight products.'
}
stopifnot(d$weight_difference==0,d$kernel_difference==0)
checks$target_exclusion <- 'Exact equality after perturbing only target run 1.'
checks$noise_change <- max(abs(as.matrix(d$soft$noise$run1[,c('rho','psi1','psi2')])-as.matrix(d$perturbed$noise$run1[,c('rho','psi1','psi2')])))
checks$null_status <- table(d$soft$selections$run1$status[17:24])
checks$noise_boundary <- table(d$soft$noise$run1$method)
writeLines(capture.output(checks),'verification.txt')
