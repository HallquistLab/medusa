#' Canonical SPM difference-of-gammas hemodynamic response function
#'
#' @description
#' Evaluates the canonical SPM difference-of-gammas hemodynamic response
#' function at repetition time `TR`, subsampled from 1/16th TR microtime
#' resolution and normalized to unit sum.
#'
#' @param TR Numeric scalar scan repetition time in seconds.
#' @param P Optional numeric parameter vector of length up to 7:
#'   delay of response (6 s), delay of undershoot (16 s), dispersion of
#'   response (1), dispersion of undershoot (1), response-to-undershoot ratio
#'   (6), onset (0 s), and kernel length (32 s).
#'
#' @return A list with elements `hrf` (a single-column numeric matrix containing
#'   the unit-sum HRF kernel on the `TR` grid) and `p` (the 7-element parameter
#'   vector used).
#' @export
spm_hrf <- function(TR, P) {
  fMRI_T <- 16
  p <- c(6, 16, 1, 1, 6, 0, 32)
  if (!missing(P)) {
    p[seq_along(P)] <- P
  }
  dt <- TR / fMRI_T
  u <- 0:(p[7] / dt) - p[6] / dt
  hrf <- stats::dgamma(u, shape = p[1] / p[3], rate = dt / p[3]) -
    stats::dgamma(u, shape = p[2] / p[4], rate = dt / p[4]) / p[5]
  hrf <- hrf[0:(p[7] / TR) * fMRI_T + 1]
  hrf <- hrf / sum(hrf)
  list(hrf = matrix(hrf, ncol = 1L), p = p)
}
