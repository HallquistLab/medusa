library(medusa)
library(data.table)
set.seed(62112026)
nrep <- 180L
n_mc <- 799L
lags <- -1:4
L <- length(lags)
C <- diag(L)
R <- toeplitz(.65^(0:(L-1)))
truth <- c(0, .06, .22, .26, .19, .08)
w <- c(0, -.5, -.5, 0, .5, .5)
subj <- 20L
m <- 20L
r_factor <- chol(R)
Rmat <- function(n, sd) matrix(rnorm(n*L, sd=sd), n) %*% r_factor
ar_path <- function(z, rho) {
  for (j in 2:nrow(z)) z[j,] <- rho*z[j-1L,] + sqrt(1-rho^2)*z[j,]
  z
}
rows <- list()
for (scenario in c("Complete / independent trials", "Gaps / serially dependent trials")) {
  stressed <- startsWith(scenario, "Gaps")
  for (i in seq_len(nrep)) {
    base <- CJ(id = seq_len(subj), run = 1:2, trial = seq_len(m))
    base[, x := if(stressed) as.numeric(ar_path(matrix(rnorm(m),ncol=1),.65)) else rnorm(.N),
         by=.(id,run)]
    base[, onset := cumsum(runif(.N, 3, 9)), by=.(id,run)]
    N <- nrow(base)
    rand <- function() {
      e <- Rmat(N, .65)
      if(stressed) {
        for (start in seq(1,N,by=m)) {
          k <- start:(start+m-1L)
          e[k,] <- ar_path(e[k,,drop=FALSE],.62)
        }
      }
      e
    }
    shared <- rand()
    dt <- rbindlist(lapply(c("Signal","Null A","Null B"), function(p) {
      b <- Rmat(subj,.13) + if(p == "Signal") rep(truth,each=subj) else 0
      e <- (.55 * shared + sqrt(1-.55^2)*rand()) * ifelse(base$run == 1,1,1.6)
      y <- base$x*b[base$id,] + e
      rbindlist(lapply(seq_len(L), function(j) {
        out <- copy(base)
        out[, `:=`(parcel=p,evt_time=lags[j],decon_demixed=y[,j],precision_weight=ifelse(run==1,1,.625^2))]
        if(stressed) {
          out <- out[!trial %in% c(6L,7L,8L,16L)]
          if(j == L) out <- out[id > 3]
        }
        out
      }))
    }))
    f <- medusa_multilag_infer(dt,covariates="x",lag_resolution=C,lag_error_cor=R,
      trajectory_mode="none",use_whitened=FALSE,contrasts=list(window=w),
      n_perm=n_mc,n_boot_latency=0L)
    per <- f$coefficients[, {
      target <- if(parcel[1L] == "Signal") truth else rep(0,L)
      list(maxT=min(p_maxT),tfce=min(p_tfce),omnibus=omnibus_p[1L],
           band=all(conf.low.sim<=target & conf.high.sim>=target),subjects=n_subjects[1L])
    },by=parcel]
    stopifnot(nrow(per)==3L, all(per$subjects == 20L))
    per[, holm_omnibus := p.adjust(omnibus,"holm")]
    con <- f$contrasts[parcel=="Signal"]
    delta <- sum(w*truth)
    rows[[length(rows)+1L]] <- data.table(scenario, replicate=i, subjects=per$subjects[1L],
      null_maxT = per[parcel=="Null A",maxT<=.05], null_tfce = per[parcel=="Null A",tfce<=.05],
      null_omnibus = per[parcel=="Null A",omnibus<=.05],
      null_band = per[parcel=="Null A",band], active_band = per[parcel=="Signal",band],
      raw_epoch_family = any(per[parcel!="Signal",omnibus<=.05]),
      holm_epoch_family = any(per[parcel!="Signal",holm_omnibus<=.05]),
      window_coverage = con$conf.low <= delta & con$conf.high >= delta,
      window_error = con$estimate-delta)
  }
}
draws <- rbindlist(rows)
summary <- melt(draws,id.vars=c("scenario","replicate","subjects","window_error"),
  variable.name="metric",value.name="event")[, {
    ci <- binom.test(sum(event),.N)$conf.int
    list(n=.N, rate=mean(event), low=ci[1L], high=ci[2L])
  },by=.(scenario,metric)]
fwrite(summary,"calibration_summary.csv")
fwrite(draws,"calibration_draws.csv")
print(summary)
