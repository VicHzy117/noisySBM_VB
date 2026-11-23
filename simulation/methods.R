source("method/VB/load.R")
source("simulation/sc.R")
.greedy <- new.env(parent=globalenv())
.greedy$loaded<-FALSE
load_method <- function(method) {
  if(method=="VB")load_vb()
  if(method=="VEM") {
    if(!requireNamespace("noisySBM",quietly=TRUE))stop("Install noisySBM")
    if(as.character(packageVersion("noisySBM"))!="0.1.4")
      warning("Paper used noisySBM 0.1.4; a different version may change fits")
  }
  if(method=="Greedy" && !.greedy$loaded) {
    cache<-file.path(tempdir(),paste0("greedy_",Sys.getpid()))
    .greedy$sourceCpp<-function(file,...)Rcpp::sourceCpp(file,env=.greedy,cacheDir=cache,
      rebuild=TRUE,showOutput=FALSE,verbose=FALSE,...)
    old<-getwd();on.exit(setwd(old),add=TRUE);setwd("simulation/Greedy")
    for(f in sort(list.files("R",pattern="[.]R$",full.names=TRUE)))sys.source(f,.greedy)
    setwd(old);sys.source("simulation/greedy_helpers.R",.greedy)
    .greedy$loaded<-TRUE
  }
  invisible(TRUE)
}
fit_method <- function(method,X,K_max,seed,count_pmf_rebuild=TRUE) {
  if(!is.na(seed))set.seed(seed)
  if(method=="VB")return(vb_fit(X,K_candidates=seq_len(K_max),seed=seed,
    count_pmf_rebuild=count_pmf_rebuild,epsilon=.1,kmeans_nstart=25L,n_perturb=2L,
    perturb_fraction=.15,screen_iter=15L,max_iter=300L,tolerance=1e-6,min_iter=3L))
  if(method=="VEM") {
    all<-noisySBM::fitNSBM(X,model="Gauss01",sbmSize=list(Qmin=1L,Qmax=K_max,explor=1.3),
      initParam=list(nbOfTau=2L,nbOfPointsPerTau=2L,maxNbOfPasses=2L,minNbOfPasses=1L),nbCores=1L)
    K<-as.integer(noisySBM::getBestQ(all)[[1L]])
    return(list(K=K,labels=as.integer(all[[K]]$clustering),fit=all[[K]]))
  }
  if(method=="Greedy") {
    f<-.greedy$fit_noisysbmGGM_once(X,NIG=TRUE,Qup=K_max,nbCores=1L,nbOfZ=12L,Nbrepet=5L,
      sigma0=1,verbatim=FALSE)
    return(list(K=as.integer(f$Q),labels=as.integer(f$Z),fit=f))
  }
  if(method=="SC")return(sc_prepare(X[edge_indices(nrow(X))],gamma=.1))
  stop("Unknown method")
}
test_method <- function(method,object,X,alpha) {
  ij<-edge_indices(nrow(X))
  if(method=="VB")return(object$fit$plugin_qvalue_mfdr<alpha)
  if(method=="VEM")return(noisySBM::graphInference(X,object$fit$clustering,object$fit$theta,
    alpha=alpha,modelFamily="Gauss")$A[ij]!=0)
  if(method=="Greedy")return(.greedy$infer_noisysbmGGM_graph(object$fit,alpha)[ij]!=0)
  sc_apply(object,alpha)
}
