# Source from the replication root, or through simulation/common.R.
source("data/generate/settings.R")
replicate_seed <- function(setting,n,replicate) {
  tab<-readRDS("data/generate/seeds.rds")
  x<-tab[tab$setting==setting & tab$n==n & tab$replicate==replicate,]
  if(nrow(x)!=1L)stop("Unknown setting/n/replicate")
  x
}
generate_replicate <- function(setting,n,replicate) {
  spec<-simulation_settings()[[setting]]
  if(is.null(spec))stop("Unknown setting")
  s<-replicate_seed(setting,n,replicate)
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(s$data_seed)
  use<-upper.tri(matrix(FALSE,n,n));ij<-which(use,arr.ind=TRUE)
  A<-matrix(0L,n,n);X<-matrix(0,n,n);Z<-NULL;center<-NULL;root_redraws<-0L
  if(spec$family %in% c("two_community","near_duplicate")) {
    K<-spec$K_true;sizes<-rep(n%/%K,K);if(n%%K)sizes[seq_len(n%%K)]<-sizes[seq_len(n%%K)]+1L
    Z<-sample(rep(seq_len(K),sizes),n,replace=FALSE)
    block<-cbind(Z[ij[,1]],Z[ij[,2]])
    av<-rbinom(nrow(ij),1,spec$Q[block]);m<-spec$mu[block]
    if(s$observation_stream=="null_then_signal") {
      xv<-rnorm(length(av));active<-which(av==1L)
      if(length(active))xv[active]<-rnorm(length(active),m[active],sqrt(spec$sigma2[block][active]))
    } else xv<-rnorm(length(av),ifelse(av==1L,m,0),ifelse(av==1L,sqrt(spec$sigma2[block]),1))
    A[use]<-av;X[use]<-xv;A<-A+t(A);X<-X+t(X)
  } else {
    if(spec$family=="bipartite") {
      Z<-sample(rep(1:2,each=n/2))
      U<-runif(nrow(ij));A[use]<-as.integer(Z[ij[,1]]!=Z[ij[,2]] & U<spec$Q[1,2]);A<-A+t(A)
    } else if(spec$family=="star") {
      center<-sample.int(n,1L);A[center,-center]<-A[-center,center]<-1L
    } else {
      r<-as.integer(.2*n);m<-as.integer(floor(.75*r+.5))
      repeat {
        a<-matrix(0L,r,r);a[upper.tri(a)]<-rbinom(choose(r,2),1,.5);a<-a+t(a)
        if(sum(rowSums(a)>0)>=m)break
        root_redraws<-root_redraws+1L
      }
      A[seq_len(r),seq_len(r)]<-a;degree<-rowSums(A)
      for(t in seq.int(r+1L,n)) {
        targets<-sample.int(t-1L,m,replace=FALSE,prob=degree[seq_len(t-1L)])
        A[t,targets]<-A[targets,t]<-1L;degree[targets]<-degree[targets]+1L;degree[t]<-m
      }
      permutation<-sample.int(n);A<-A[permutation,permutation]
    }
    set.seed(s$noise_seed);X[use]<-rnorm(choose(n,2));X<-X+t(X)+spec$mu_scalar*A
  }
  list(setting=setting,n=as.integer(n),replicate=as.integer(replicate),X=X,A=A,Z=Z,
    center=center,data_seed=s$data_seed,noise_seed=s$noise_seed,root_redraws=root_redraws)
}
read_replicate <- function(setting,n,replicate) {
  if(!setting %in% names(simulation_settings()) || !n %in% simulation_sizes || !replicate %in% 1:100)
    stop("Unknown setting/n/replicate")
  x<-readRDS(file.path("data/datasets",setting,sprintf("n%03d.rds",n)))
  x[[as.integer(replicate)]]
}
