# Research prototype, not a production MEDuSA estimator.
# Re-run: Rscript connectivity_prototype.R /path/to/output-directory
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args)) args[1] else "."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(20261110)

# Exact covariance transformation for centered seed, task, and seed*task.
# These inputs are RAW measurement covariances, before any trial df calibration.
error_recipe <- function(Os, Osy, z) {
  n <- length(z); P <- diag(n) - matrix(1/n, n, n); Dz <- diag(z)
  A <- list(P, matrix(0, n, n), P %*% Dz)
  Crows <- array(0, c(n, 3, 3))
  for (k in c(1L,3L)) for (l in c(1L,3L))
    Crows[,k,l] <- rowSums((A[[k]] %*% Os) * A[[l]])
  crows <- matrix(0, n, 3)
  for (k in c(1L,3L)) crows[,k] <- rowSums((A[[k]] %*% Osy) * P)
  list(Crows=Crows, crows=crows,
       C=apply(Crows,c(2,3),sum), c=colSums(crows))
}

# One independent participant is the inference cluster.
# The common-slope estimand is explicit; this is not a random-slope population mean.
fit_corrected <- function(blocks, mode=c("full","diagonal","zero_shared","naive"), row_se=FALSE) {
  mode <- match.arg(mode)
  q <- 3L; A <- matrix(0,q,q); b <- numeric(q)
  prepared <- vector("list",length(blocks))
  for (g in seq_along(blocks)) {
    d <- blocks[[g]]; z <- d$z
    X <- cbind(neural=d$s, task=z, ppi=d$s*z)
    X <- sweep(X,2,colMeans(X),"-"); y <- d$y-mean(d$y)
    Os <- d$Os; Osy <- d$Osy
    if (mode == "diagonal") {Os <- diag(diag(Os)); Osy <- matrix(0,nrow(Os),ncol(Os))}
    if (mode == "zero_shared") Osy[] <- 0
    if (mode == "naive") {Os[] <- 0; Osy[] <- 0}
    recipe <- error_recipe(Os,Osy,z)
    Ag <- crossprod(X)-recipe$C; bg <- as.numeric(crossprod(X,y))-recipe$c
    A <- A+Ag; b <- b+bg
    prepared[[g]] <- list(X=X,y=y,recipe=recipe,A=Ag,b=bg)
  }
  ev <- eigen((A+t(A))/2,symmetric=TRUE,only.values=TRUE)$values
  if (!all(is.finite(ev)) || min(ev)<1e-7 || max(ev)/min(ev)>1e8)
    return(list(ok=FALSE,reason="weak/indefinite corrected moment"))
  theta <- as.numeric(solve(A,b)); bread <- solve(A)
  scores <- do.call(rbind,lapply(prepared,function(d) as.numeric(d$b-d$A%*%theta)))
  if (isTRUE(row_se)) {
    row_scores <- do.call(rbind,lapply(prepared,function(d) {
      ans <- d$X*as.numeric(d$y-d$X%*%theta)-d$recipe$crows
      for (k in 1:3) for (l in 1:3) ans[,k] <- ans[,k]+d$recipe$Crows[,k,l]*theta[l]
      ans
    }))
    meat <- crossprod(row_scores)
  } else {
    G <- length(blocks);meat <- crossprod(scores)*G/(G-1)
  }
  V <- bread%*%meat%*%t(bread)
  if(any(!is.finite(V)) || any(diag(V)<=0))return(list(ok=FALSE,reason="invalid covariance"))
  list(ok=TRUE,estimate=theta,se=sqrt(diag(V)),A=A,scores=scores)
}

n <- 24L
P <- diag(n)-matrix(1/n,n,n)
shift <- diag(n)[,c(n,seq_len(n-1))]
S <- .5*diag(n)+.25*(shift+t(shift))
Ks <- .8*(.3*diag(n)+.7*S)
Kt <- 1.05*(.55*diag(n)+.45*S)
invKs <- solve(Ks);invKt <- solve(Kt)
Sigma <- outer(1:n,1:n,function(i,j).55^abs(i-j))
L <- t(chol(Sigma))
Le <- t(chol(outer(1:n,1:n,function(i,j).65^abs(i-j))))
true <- c(.55,.25,.35)
shock <- function(n,heavy) if (heavy) rt(n,df=5)*sqrt(3/5) else rnorm(n)
make_blocks <- function(G,heavy=FALSE,shared=TRUE,resolution=FALSE,zero=FALSE) {
  rho <- if(shared).45 else 0
  Os <- .65^2*Sigma; Ot <- .55^2*Sigma; Osy <- rho*.65*.55*Sigma
  beta <- if(zero) 0 else true[1]; delta <- if(zero) 0 else true[3]
  lapply(seq_len(G),function(g) {
    z <- as.numeric(scale(rnorm(n)));z <- z-mean(z)
    latent_s <- .35*z+as.numeric(L%*%shock(n,heavy))
    latent_y <- beta*latent_s+.25*z+delta*latent_s*z+as.numeric(Le%*%shock(n,heavy))*.8
    e1 <- shock(n,heavy);e2 <- shock(n,heavy)
    u <- as.numeric(.65*L%*%e1)
    v <- as.numeric(.55*L%*%(rho*e1+sqrt(1-rho^2)*e2))
    s <- if(resolution) as.numeric(Ks%*%latent_s)+u else latent_s+u
    y <- if(resolution) as.numeric(Kt%*%latent_y)+v else latent_y+v
    list(s=s,y=y,z=z,Os=Os,Ot=Ot,Osy=Osy,latent_s=latent_s,latent_y=latent_y)
  })
}
undo_resolution <- function(blocks) lapply(blocks,function(d){
  d$s <- as.numeric(invKs%*%d$s);d$y <- as.numeric(invKt%*%d$y)
  d$Os <- invKs%*%d$Os%*%t(invKs);d$Ot <- invKt%*%d$Ot%*%t(invKt)
  d$Osy <- invKs%*%d$Osy%*%t(invKt);d
})

# Independent Monte Carlo verifies transformed interaction moments, including centering.
z <- as.numeric(scale(seq_len(n)))
Os <- .65^2*Sigma;Osy <- .45*.65*.55*Sigma
r <- error_recipe(Os,Osy,z);B <- 25000L
u <- .65*L%*%matrix(rnorm(n*B),n,B)
v <- .55*(.45*(u/.65)+sqrt(1-.45^2)*(L%*%matrix(rnorm(n*B),n,B)))
uc <- P%*%u;pc <- P%*%(z*u);vc <- P%*%v
empC <- matrix(c(mean(colSums(uc^2)),0,mean(colSums(uc*pc)),0,0,0,
 mean(colSums(uc*pc)),0,mean(colSums(pc^2))),3,3)
empc <- c(mean(colSums(uc*vc)),0,mean(colSums(pc*vc)))
stopifnot(max(abs(empC-r$C))<.12, max(abs(empc-r$c))<.12)
# A moment equation reduces to ordinary within-person OLS when sampling error is zero.
b0 <- make_blocks(10)
b0 <- lapply(b0,function(d){d$s<-d$latent_s;d$y<-d$latent_y;d$Os[]<-0;d$Osy[]<-0;d})
f0 <- fit_corrected(b0)
xx <- do.call(rbind,lapply(b0,function(d){X<-cbind(d$s,d$z,d$s*d$z);sweep(X,2,colMeans(X),"-")}))
yy <- unlist(lapply(b0,function(d)d$y-mean(d$y)))
stopifnot(max(abs(f0$estimate-as.numeric(solve(crossprod(xx),crossprod(xx,yy)))))<1e-10)
# Task multiplication and smoothing do not commute.
commutator <- norm(Ks%*%diag(z)-diag(z)%*%Ks,"F")
stopifnot(commutator>0)
writeLines(c(sprintf("max interaction moment MC discrepancy = %.5f",max(abs(empC-r$C))),
 sprintf("max shared-error moment MC discrepancy = %.5f",max(abs(empc-r$c))),
 sprintf("PPI commutator Frobenius norm = %.5f",commutator),
 "zero-error OLS identity: PASS"),file.path(out_dir,"prototype_checks.txt"))

configs <- data.frame(
 scenario=c("Independent Gaussian","Serial + shared Gaussian","Serial + shared t5", "Different resolution", "Few participants", "Null serial + shared"),
 G=c(60,60,60,60,12,60),heavy=c(FALSE,FALSE,TRUE,FALSE,FALSE,FALSE),
 shared=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),resolution=c(FALSE,FALSE,FALSE,TRUE,FALSE,FALSE),
 zero=c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE))
# "Independent" means measurement correlation is set to identity as well as shared=FALSE.
R <- as.integer(Sys.getenv("MEDUSA_REPS","400"))
rows <- vector("list",nrow(configs)*R*7L*2L);idx<-0L
for (cc in seq_len(nrow(configs))) {
 cfg<-configs[cc,];cat("Scenario:",cfg$scenario,"G=",cfg$G,"reps=",R,"\n")
 Sigma_saved<-Sigma;L_saved<-L;Le_saved<-Le
 if(cc==1L){Sigma<-diag(n);L<-diag(n);Le<-diag(n)}
 for(rep in seq_len(R)) {
  blocks<-make_blocks(cfg$G,cfg$heavy,cfg$shared,cfg$resolution,cfg$zero)
  oracle<-lapply(blocks,function(d){d$s<-d$latent_s;d$y<-d$latent_y;d$Os[]<-0;d$Osy[]<-0;d})
  fits<-list(
   Oracle=fit_corrected(oracle),
   Naive=fit_corrected(blocks,"naive"),
   Diagonal=fit_corrected(blocks,"diagonal"),
   `Zero shared covariance`=fit_corrected(blocks,"zero_shared"),
   `Full moments, row SE`=fit_corrected(blocks,"full",TRUE),
   `Full moments, cluster SE`=fit_corrected(blocks,"full"))
  if(cfg$resolution) fits[["Resolution inverse + full moments"]]<-fit_corrected(undo_resolution(blocks),"full")
  for(method in names(fits))for(k in c(1L,3L)) {
   f<-fits[[method]];truth<-if(cfg$zero)0 else true[k];idx<-idx+1L
   est<-if(f$ok)f$estimate[k] else NA_real_;se<-if(f$ok)f$se[k] else NA_real_
   rows[[idx]]<-data.frame(scenario=cfg$scenario,G=cfg$G,rep=rep,method=method,
    term=if(k==1L)"Connectivity" else "PPI",truth=truth,ok=f$ok,estimate=est,se=se,
    covered=if(f$ok)abs(est-truth)<=qnorm(.975)*se else NA,
    rejected=if(f$ok)abs(est)>qnorm(.975)*se else NA)
  }
 }
 Sigma<-Sigma_saved;L<-L_saved;Le<-Le_saved
}
results<-do.call(rbind,rows[seq_len(idx)])
write.csv(results,file.path(out_dir,"simulation_draws.csv"),row.names=FALSE)
keys<-interaction(results$scenario,results$method,results$term,drop=TRUE)
summary<-do.call(rbind,lapply(split(results,keys),function(d) {
 ok<-d$ok;good<-d[ok,];N<-sum(ok)
 data.frame(scenario=d$scenario[1],G=d$G[1],method=d$method[1],term=d$term[1],truth=d$truth[1],
 attempts=nrow(d),usable=N,failures=sum(!ok),mean=if(N)mean(good$estimate)else NA,
 bias=if(N)mean(good$estimate-good$truth)else NA,empirical_sd=if(N)sd(good$estimate)else NA,
 mean_se=if(N)mean(good$se)else NA,coverage=if(N)mean(good$covered)else NA,
 coverage_mcse=if(N)sqrt(mean(good$covered)*(1-mean(good$covered))/N)else NA,
 rejection=if(N)mean(good$rejected)else NA)
}))
rownames(summary)<-NULL
write.csv(summary,file.path(out_dir,"simulation_summary.csv"),row.names=FALSE)
saveRDS(list(configs=configs,truth=true,n=n,reps=R,session=sessionInfo()),file.path(out_dir,"simulation_settings.rds"))
cat("DONE",nrow(results),"coefficient draws;",sum(!results$ok),"unusable coefficient draws\n")
print(summary[summary$method%in%c("Full moments, cluster SE","Resolution inverse + full moments"),],row.names=FALSE)
