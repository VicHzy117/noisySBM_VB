# Source this file; no working-directory assumptions or precompiled objects.
.vb_root <- local({
  paths<-vapply(sys.frames(),function(f)if(is.null(f$ofile))NA_character_ else as.character(f$ofile),character(1))
  dirname(normalizePath(tail(paths[!is.na(paths)],1L),mustWork=TRUE))
})
.vb_backend <- new.env(parent=globalenv())
.vb_backend$root <- .vb_root
.vb_backend$loaded <- FALSE
sys.source(file.path(.vb_root,"R/spectral_initialization.R"),.vb_backend)

load_vb <- function(rebuild=FALSE) {
  if(.vb_backend$loaded && !rebuild) return(invisible(TRUE))
  for(p in c("Rcpp","RcppArmadillo"))
    if(!requireNamespace(p,quietly=TRUE)) stop("Install ",p," first")
  cache <- file.path(tempdir(),paste0("structured_vb_",Sys.getpid()))
  Rcpp::sourceCpp(file.path(.vb_backend$root,"src/nsbm_vi.cpp"),env=.vb_backend,
    cacheDir=cache,rebuild=rebuild,showOutput=FALSE,verbose=FALSE)
  Rcpp::sourceCpp(file.path(.vb_backend$root,"src/spectral_initialization.cpp"),env=.vb_backend,
    cacheDir=cache,rebuild=rebuild,showOutput=FALSE,verbose=FALSE)
  .vb_backend$loaded <- TRUE
  invisible(TRUE)
}
vb_hyperparameters <- function() list(alpha=2,beta=2,alpha_prime=2,beta_prime=2,
  rou=0,kappa=.01,delta=3,xi=2,gamma=1)
vb_block_matrix <- function(x,K) {
  stopifnot(length(x)==K*(K+1)/2)
  a<-matrix(NA_real_,K,K);b<-1L
  for(r in seq_len(K)) for(s in r:K) {a[r,s]<-a[s,r]<-x[b];b<-b+1L}
  a
}
vb_fit <- function(X,K=NULL,K_candidates=NULL,seed=1L,
                   hyperparameters=vb_hyperparameters(),count_pmf_rebuild=TRUE,...) {
  if(!is.null(K) && !is.null(K_candidates)) stop("Specify K or K_candidates, not both")
  load_vb()
  if(!is.null(K)) {
    answer<-.vb_backend$fit_nsbm_vi_spectral(X,K=K,seed=seed,
      hyperparameters=hyperparameters,count_pmf_rebuild=count_pmf_rebuild,...)
    f<-answer$fit;selected_K<-as.integer(K);path<-data.frame(K=K,final_elbo=tail(f$elbo_trace,1),selected=TRUE)
  } else {
    if(is.null(K_candidates)) K_candidates<-seq_len(min(6L,nrow(X)))
    answer<-.vb_backend$fit_nsbm_vi_select_k(X,K_candidates=K_candidates,seed=seed,
      hyperparameters=hyperparameters,count_pmf_rebuild=count_pmf_rebuild,...)
    f<-answer$selected$fit;selected_K<-answer$selected_K;path<-answer$path
  }
  p<-f$plugin_parameters
  structure(list(K=selected_K,labels=max.col(f$tau,ties.method="first"),tau=f$tau,
    parameters=list(Q=vb_block_matrix(p$Q,selected_K),mu=vb_block_matrix(p$mu,selected_K),
      sigma2=vb_block_matrix(p$sigma2,selected_K),pi=p$pi),
    elbo=tail(f$elbo_trace,1),path=path,fit=f,seed=seed,
    count_pmf_rebuild=count_pmf_rebuild),class="structured_vb")
}
vb_test <- function(object,alpha=.05) {
  stopifnot(inherits(object,"structured_vb"),length(alpha)==1L,is.finite(alpha),alpha>=0,alpha<=1)
  n<-nrow(object$tau);q<-object$fit$plugin_qvalue_mfdr
  ij<-cbind(rep.int(seq_len(n-1L),n-seq_len(n-1L)),
    unlist(lapply(seq_len(n-1L),function(i)seq.int(i+1L,n)),use.names=FALSE))
  a<-matrix(0L,n,n);a[ij]<-as.integer(q<alpha);a<-a+t(a)
  list(A=a,alpha=alpha,qvalues=q,lvalues=object$fit$plugin_lvalue,
    pairs=ij,discoveries=sum(q<alpha))
}
# Original entry points remain available, with their original argument names.
fit_nsbm_vi_spectral <- function(...) {load_vb();.vb_backend$fit_nsbm_vi_spectral(...)}
fit_nsbm_vi_select_k <- function(...) {load_vb();.vb_backend$fit_nsbm_vi_select_k(...)}
