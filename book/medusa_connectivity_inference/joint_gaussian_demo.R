# One fitted joint Gaussian latent-seed/PPI measurement model.
# Fixed known noise covariances and resolution, with estimated structural and latent parameters.
args<-commandArgs(TRUE);out<-if(length(args))args[1]else"."
set.seed(90812);n<-24L;G<-60L
R<-outer(1:n,1:n,function(i,j).55^abs(i-j));Re<-outer(1:n,1:n,function(i,j).65^abs(i-j))
S<-diag(n)[,c(n,1:(n-1))];smooth<-.5*diag(n)+.25*(S+t(S))
Ks<-.8*(.3*diag(n)+.7*smooth);Kt<-1.05*(.55*diag(n)+.45*smooth)
K<-matrix(0,2*n,2*n);K[1:n,1:n]<-Ks;K[n+1:n,n+1:n]<-Kt
Om<-rbind(cbind(.65^2*R,.45*.65*.55*R),cbind(.45*.65*.55*R,.55^2*R))
L<-t(chol(R));Le<-t(chol(Re));Lo<-t(chol(Om))
blocks<-lapply(1:G,function(g){z<-as.numeric(scale(rnorm(n)));s<-.35*z+as.numeric(L%*%rnorm(n));
 y<-.55*s+.25*z+.35*z*s+.8*as.numeric(Le%*%rnorm(n));
 list(z=z,observed=as.numeric(K%*%c(s,y)+Lo%*%rnorm(2*n)))})
nll<-function(p) {
 beta<-p[1];gamma<-p[2];delta<-p[3];a<-p[4];Vs<-exp(p[5])*R;Ve<-exp(p[6])*Re
 ans<-0
 for(b in blocks){z<-b$z;d<-beta+delta*z;mu<-c(a*z,(gamma+beta*a)*z+delta*a*z^2)
  # Cov(s,t) = Vs D; Cov(t,t) = D Vs D + Ve.
  Vst<-sweep(Vs,2,d,'*');Vtt<-sweep(sweep(Vs,1,d,'*'),2,d,'*')+Ve
  V<-K%*%rbind(cbind(Vs,Vst),cbind(t(Vst),Vtt))%*%t(K)+Om
  U<-chol(V);r<-b$observed-as.numeric(K%*%mu);w<-backsolve(U,r,transpose=TRUE)
  ans<-ans+sum(log(diag(U)))+.5*sum(w^2)
 }
 ans
}
fit<-optim(c(.4,.2,.2,.2,0,log(.64)),nll,method='BFGS',hessian=TRUE,control=list(maxit=500,reltol=1e-9))
stopifnot(fit$convergence==0,all(eigen(fit$hessian,symmetric=TRUE)$values>0))
se<-sqrt(diag(solve(fit$hessian)))
tab<-data.frame(term=c('Connectivity','Task','PPI','Seed task mean','log seed variance','log structural variance'),
 truth=c(.55,.25,.35,.35,0,log(.64)),estimate=fit$par,se=se)
write.csv(tab,file.path(out,'joint_gaussian_example.csv'),row.names=FALSE)
writeLines(c('optim convergence = 0',sprintf('negative log likelihood = %.6f',fit$value),
 'Conditional Gaussian Hessian SEs; this single dataset does not establish coverage.'),file.path(out,'joint_gaussian_checks.txt'))
print(tab,row.names=FALSE)
