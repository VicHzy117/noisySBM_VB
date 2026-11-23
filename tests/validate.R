script<-normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]))
root<-dirname(dirname(script));setwd(root)
source("simulation/common.R")
Sys.setenv(OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
scratch<-tempfile("paper_replication_check_");dir.create(scratch)
run_R<-function(file,args=character(),cwd=root) {
  old<-getwd();on.exit(setwd(old),add=TRUE);setwd(cwd)
  output<-system2(file.path(R.home("bin"),"Rscript"),c("--vanilla",shQuote(file),shQuote(args)),stdout=TRUE,stderr=TRUE)
  cat(paste(tail(output,10),collapse="\n"),"\n")
  if(!is.null(attr(output,"status")) && attr(output,"status")!=0L)stop("Child R check failed: ",file)
}
file.copy("method/VB",scratch,recursive=TRUE)
run_R(file.path(root,"tests/standalone_vb.R"),file.path(scratch,"VB"),cwd=scratch)
checks<-list(data.frame(check="standalone_VB",passed=TRUE,max_abs_error=0))
maximum<-0
for(s in paste0("S",1:14))for(n in simulation_sizes) {
  archived<-read_replicate(s,n,1L);g<-generate_replicate(s,n,1L)
  err<-max(abs(g$X-archived$X));maximum<-max(maximum,err)
  stopifnot(identical(g$A,archived$A),identical(g$Z,archived$Z),err<=1e-12)
}
checks[[length(checks)+1L]]<-data.frame(check="confirmed_seeds_56_cells",passed=TRUE,max_abs_error=maximum)
published<-readRDS("results/data/metrics.rds")
for(method in simulation_methods) {
  args<-c("--setting=S1","--n=50","--replicate=1",paste0("--method=",method),paste0("--output=",file.path(scratch,"tasks")))
  run_R(file.path(root,"simulation/run_simulation.R"),args)
  p<-file.path(scratch,"tasks",method,"S1/n050_rep001.rds");x<-readRDS(p)
  reference<-published[published$setting=="S1" & published$n==50 & published$replicate==1 & published$method==method,]
  reference<-reference[match(x$metrics$alpha,reference$alpha),]
  difference<-max(abs(x$metrics$TP-reference$TP),abs(x$metrics$FP-reference$FP))
  checks[[length(checks)+1L]]<-data.frame(check=paste0("method_metrics_",method),passed=TRUE,max_abs_error=difference)
  # A repeat must reuse the same task and not append a second record.
  stamp<-file.info(p)$mtime;run_R(file.path(root,"simulation/run_simulation.R"),args)
  stopifnot(identical(stamp,file.info(p)$mtime),nrow(x$metrics)==5L,
    all(x$metrics$FP/pmax(1,x$metrics$discoveries)==x$metrics$FDP))
}
isolated<-file.path(scratch,"results_only");dir.create(isolated)
dir.create(file.path(isolated,"results"));file.copy("results/data",file.path(isolated,"results"),recursive=TRUE)
file.copy("results/code",file.path(isolated,"results"),recursive=TRUE)
dir.create(file.path(isolated,"data/generate"),recursive=TRUE)
file.copy("data/generate/settings.R",file.path(isolated,"data/generate"))
stopifnot(!dir.exists(file.path(isolated,"data/datasets")),!dir.exists(file.path(isolated,"method")))
run_R(file.path(isolated,"results/code/reproduce.R"),cwd=isolated)
v<-read.csv(file.path(isolated,"results/tables/validation.csv"))
stopifnot(nrow(v)==3L,all(v$max_abs_error<1e-12))
checks[[length(checks)+1L]]<-data.frame(check="results_only_figures_and_tables",passed=TRUE,max_abs_error=max(v$max_abs_error))
checks<-do.call(rbind,checks)
write.csv(checks,"tests/validation_results.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"tests/local_session.txt")
cat("\nAll required validation checks passed. Scratch copy:",scratch,"\n");print(checks,row.names=FALSE)
