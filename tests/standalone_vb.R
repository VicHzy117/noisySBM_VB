# Invoked by validate.R in a fresh R session outside the replication directory.
args<-commandArgs(TRUE);stopifnot(length(args)==1L)
source(file.path(args[1],"load.R"))
set.seed(17);n<-24L;z<-rep(1:2,each=n/2);ij<-which(upper.tri(matrix(0,n,n)),arr.ind=TRUE)
a<-rbinom(nrow(ij),1,ifelse(z[ij[,1]]==z[ij[,2]],.8,.1));v<-rnorm(length(a),2*a,1)
X<-matrix(0,n,n);X[ij]<-v;X<-X+t(X)
fixed<-vb_fit(X,K=2,seed=11)
selected<-vb_fit(X,K_candidates=1:3,seed=11)
for(f in list(fixed,selected)) {
  g<-vb_test(f,.05)
  stopifnot(is.finite(f$elbo),nrow(f$tau)==n,ncol(f$tau)==f$K,
    max(abs(rowSums(f$tau)-1))<1e-12,length(f$labels)==n,
    all(dim(f$parameters$Q)==c(f$K,f$K)),all(f$parameters$sigma2>0),
    all(dim(g$A)==c(n,n)),all(g$A %in% 0:1),identical(g$A,t(g$A)))
}
cat("PASS: isolated VB directory; fresh R; fixed K, ELBO selection, memberships, parameters, graph tests.\n")
