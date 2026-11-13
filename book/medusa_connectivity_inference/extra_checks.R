# Run from this report directory: Rscript extra_checks.R
source('measurement_error_models.R')
n<-8L;d<-list(s=seq_len(n),y=seq_len(n),z=(-3:4),Os=diag(1e8,n),Osy=matrix(0,n,n))
stopifnot(!fit_corrected(list(d,d))$ok)
rank_failure<-tryCatch({likelihood_only(1:10,cbind(1:10,1:10),matrix(1,10,1),diag(10));FALSE},error=function(e)TRUE)
stopifnot(rank_failure)
cat('Weak corrected-moment rejection: PASS\nRank-deficient likelihood rejection: PASS\n')
