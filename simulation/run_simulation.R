# Run from any directory. Example: Rscript simulation/run_simulation.R --setting=S1 --n=50 --replicate=1 --method=VB
file_arg<-grep("^--file=",commandArgs(FALSE),value=TRUE)[1]
script<-normalizePath(sub("^--file=","",file_arg))
setwd(dirname(dirname(script)))
source("simulation/common.R")
args<-commandArgs(TRUE)
option<-function(name,default=NULL) {
  p<-paste0("--",name,"=");v<-args[startsWith(args,p)]
  if(length(v)>1L)stop("Repeated option: ",name)
  if(length(v))substring(v,nchar(p)+1L)else default
}
if(!is.null(option("task-id"))) {
  grid<-task_grid();id<-as.integer(option("task-id"));stopifnot(id %in% grid$task_id)
  task<-grid[id,]
} else task<-data.frame(setting=option("setting","S1"),n=as.integer(option("n","50")),
  replicate=as.integer(option("replicate","1")),method=option("method","VB"))
stopifnot(task$setting %in% names(simulation_settings()),task$n %in% simulation_sizes,
  task$replicate %in% 1:100,task$method %in% simulation_methods)
d<-read_replicate(task$setting,task$n,task$replicate)
s<-replicate_seed(task$setting,task$n,task$replicate)
spec<-simulation_settings()[[task$setting]]
if(!requireNamespace("digest",quietly=TRUE))stop("Install digest")
input_hash<-digest::digest(d[c("X","A","Z")],algo="sha256")
method_seed<-s[[paste0(task$method,"_seed")]]
dest<-option("output","results/data/reruns")
file<-file.path(dest,task$method,task$setting,sprintf("n%03d_rep%03d.rds",task$n,task$replicate))
if(file.exists(file)) {
  x<-readRDS(file)
  stopifnot(isTRUE(x$status$success),x$metadata$input_sha256==input_hash,
    identical(x$metadata$alphas,simulation_alphas),identical(x$metadata$method_seed,method_seed))
  cat("Already complete:",file,"\n");quit(status=0L)
}
source("simulation/methods.R")
backend_start<-proc.time();load_method(task$method);backend_sec<-elapsed(backend_start)
start<-proc.time()
fit<-fit_method(task$method,d$X,spec$K_max,method_seed,s$count_pmf_rebuild)
fit_sec<-elapsed(start)
start<-proc.time();calls<-lapply(simulation_alphas,function(a)test_method(task$method,fit,d$X,a))
inference_sec<-elapsed(start)
key<-task[c("setting","n","replicate","method")]
metrics<-do.call(rbind,lapply(seq_along(simulation_alphas),function(j)
  cbind(key,alpha=simulation_alphas[j],evaluate_calls(calls[[j]],d$A))))
K<-if(task$method=="SC")NA_integer_ else fit$K
labels<-if(task$method=="SC")NULL else fit$labels
if(length(labels) && min(labels)==0L)labels<-labels+1L
fit_row<-cbind(key,data.frame(K_hat=K,K_true=spec$K_true,
  ARI=if(length(labels)&&!is.null(d$Z))ari(labels,d$Z)else NA_real_,
  fit_sec=fit_sec,data_seed=s$data_seed,method_seed=method_seed,timing_batch="new run"))
atomic_save(list(status=list(success=TRUE),metadata=list(input_sha256=input_hash,
  method_seed=method_seed,data_seed=s$data_seed,alphas=simulation_alphas,K_max=spec$K_max,
  count_pmf_rebuild=s$count_pmf_rebuild,backend_sec=backend_sec,inference_sec=inference_sec,
  R=R.version.string,session=capture.output(sessionInfo())),
  diagnostics=fit_row,metrics=metrics,fit=fit,decisions=calls),file)
print(fit_row);print(metrics);cat("Saved:",file,"\n")
