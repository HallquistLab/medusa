# Research helpers for the MEDuSA connectivity exploration.
# Not a validated production API; see report for assumptions and scope.
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


likelihood_only <- function(y,G,Z,Sigma,tol=1e-9) {
  n<-nrow(G);stopifnot(length(y)==n,nrow(Z)==n,all(dim(Sigma)==n))
  W<-solve(t(chol(Sigma)))
  Gw<-W%*%G;Q<-diag(n)
  if(ncol(Z)) {qz<-qr(W%*%Z);Q<-Q-tcrossprod(qr.Q(qz)[,seq_len(qz$rank),drop=FALSE])}
  J<-crossprod(Gw,Q%*%Gw)
  eig<-eigen((J+t(J))/2,symmetric=TRUE,only.values=TRUE)$values
  if(min(eig)<tol*max(eig))stop("Trial basis not identified without regularization")
  A<-solve(J,crossprod(Gw,Q%*%W))
  list(estimate=as.numeric(A%*%y),covariance=A%*%Sigma%*%t(A),operator=A,
       information=J,condition=max(eig)/min(eig))
}
