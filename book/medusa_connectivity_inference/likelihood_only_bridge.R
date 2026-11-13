# Fixed-HRF finite-basis GLS bridge: a research example, not a new package API.
# y = G a + Z gamma + e, Cov(e) = Sigma, with fixed known design/covariance.
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
args<-commandArgs(trailingOnly=TRUE);out<-if(length(args))args[1]else"."
dir.create(out,showWarnings=FALSE,recursive=TRUE)
set.seed(260110)
T<-180L;tr<-1;tt<-0:(T-1);m<-16L;onsets<-cumsum(c(6,runif(m-1,7.2,9.0)))
hrf<-function(t,shape){x<-pmax(t,0);h<-dgamma(x,shape=shape)-dgamma(x,shape=16)/6;h[t<0]<-0;h}
Gs<-vapply(onsets,function(o)hrf(tt-o,5),numeric(T))
Gt<-vapply(onsets,function(o)hrf(tt-o,7),numeric(T))
Z<-cbind(1,seq(-1,1,length.out=T));R<-outer(1:T,1:T,function(i,j).4^abs(i-j))
Ss<-.2^2*R;St<-.25^2*R;Sst<-.35*.2*.25*R
as<-rnorm(m);at<-.55*as+rnorm(m,sd=.8)
fs<-likelihood_only(as.numeric(Gs%*%as),Gs,Z,Ss)
ft<-likelihood_only(as.numeric(Gt%*%at),Gt,Z,St)
stopifnot(max(abs(fs$operator%*%Gs-diag(m)))<1e-10,
 max(abs(ft$operator%*%Gt-diag(m)))<1e-10,
 max(abs(fs$operator%*%Z))<1e-10)
# Posterior-precision subtraction agrees where this reduced basis is identified.
J<-fs$information;prior<-2*diag(m);C<-solve(J+prior)
recovered_J<-solve(C)-prior
stopifnot(max(abs(recovered_J-J))<1e-9)
# The same maneuver cannot identify a full neural time grid with HRF delay.
H<-vapply(1:T,function(j)hrf(tt-tt[j],5),numeric(T))
full_rank<-qr(H,tol=1e-9)$rank
stopifnot(full_rank<T)
# Verify covariance and seed-target measurement covariance by independent noise draws.
B<-10000L;L<-t(chol(R));e1<-L%*%matrix(rnorm(T*B),T,B);e2<-L%*%matrix(rnorm(T*B),T,B)
us<-fs$operator%*%(.2*e1);ut<-ft$operator%*%(.25*(.35*e1+sqrt(1-.35^2)*e2))
Ost<-fs$operator%*%Sst%*%t(ft$operator)
sd_check<-sqrt(rowMeans(us^2)/diag(fs$covariance))
stopifnot(max(abs(sd_check-1))<.04)
cov_check<-us%*%t(ut)/B
scaled_cov_error<-max(abs(cov_check-Ost)/sqrt(outer(diag(fs$covariance),diag(ft$covariance))))
stopifnot(scaled_cov_error<.04)
checks<-data.frame(check=c("Seed A*G identity max error","Target A*G identity max error",
 "Source information condition number","Target information condition number",
 "Full-grid rank","Full-grid columns","MC maximum relative source SD error",
 "MC maximum scaled cross-covariance error"),value=c(max(abs(fs$operator%*%Gs-diag(m))),
 max(abs(ft$operator%*%Gt-diag(m))),fs$condition,ft$condition,full_rank,T,max(abs(sd_check-1)),scaled_cov_error))
write.csv(checks,file.path(out,"likelihood_bridge_checks.csv"),row.names=FALSE)
print(checks,row.names=FALSE)
