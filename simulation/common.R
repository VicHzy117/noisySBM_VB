# All command-line programs resolve their own root before sourcing this file.
source("data/generate/generate.R")
bind_rows <- function(xs) {
  xs<-Filter(function(x)!is.null(x)&&is.data.frame(x)&&nrow(x)>0L,xs)
  if(!length(xs))return(data.frame())
  cols<-unique(unlist(lapply(xs,names)))
  do.call(rbind,lapply(xs,function(x){for(k in setdiff(cols,names(x)))x[[k]]<-NA;x[cols]}))
}
edge_indices <- function(n)cbind(rep.int(seq_len(n-1L),n-seq_len(n-1L)),
  unlist(lapply(seq_len(n-1L),function(i)seq.int(i+1L,n)),use.names=FALSE))
atomic_save <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temp<-tempfile(".writing_",dirname(path));on.exit(unlink(temp),add=TRUE)
  saveRDS(x,temp,compress="gzip",version=3);if(!file.rename(temp,path))stop("Cannot save ",path)
}
elapsed <- function(t)unname((proc.time()-t)["elapsed"])
ari <- function(x,y) {
  t<-table(x,y);c2<-function(v)v*(v-1)/2
  a<-sum(c2(rowSums(t)));b<-sum(c2(colSums(t)));obs<-sum(c2(t));e<-a*b/c2(length(x));den<-(a+b)/2-e
  if(abs(den)<.Machine$double.eps)as.numeric(obs==(a+b)/2)else(obs-e)/den
}
evaluate_calls <- function(calls,A) {
  truth<-A[edge_indices(nrow(A))]!=0
  stopifnot(length(calls)==length(truth),!anyNA(calls))
  tp<-sum(calls & truth);fp<-sum(calls & !truth)
  data.frame(TP=tp,FP=fp,TN=sum(!calls & !truth),FN=sum(!calls & truth),
    discoveries=tp+fp,FDP=fp/max(1,tp+fp),TDR=if(sum(truth))tp/sum(truth)else NA_real_)
}
task_grid <- function() {
  x<-expand.grid(method=simulation_methods,setting=paste0("S",1:14),n=simulation_sizes,
    replicate=1:100,stringsAsFactors=FALSE,KEEP.OUT.ATTRS=FALSE)
  x$task_id<-seq_len(nrow(x));x
}
