file_arg<-grep("^--file=",commandArgs(FALSE),value=TRUE)[1]
setwd(dirname(dirname(normalizePath(sub("^--file=","",file_arg)))))
source("simulation/common.R")
args<-commandArgs(TRUE)
input<-if(length(args))args[1]else"results/data/reruns"
output<-if(length(args)>1)args[2]else"results/data/new_summary"
files<-list.files(input,pattern="[.]rds$",recursive=TRUE,full.names=TRUE)
if(!length(files))stop("No task results in ",input)
objects<-lapply(files,readRDS)
stopifnot(all(vapply(objects,function(x)isTRUE(x$status$success),logical(1))))
fits<-bind_rows(lapply(objects,`[[`,"diagnostics"));metrics<-bind_rows(lapply(objects,`[[`,"metrics"))
stopifnot(!anyDuplicated(fits[c("setting","n","replicate","method")]),
  !anyDuplicated(metrics[c("setting","n","replicate","method","alpha")]))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
saveRDS(fits,file.path(output,"fits.rds"));saveRDS(metrics,file.path(output,"metrics.rds"))
cat("Collected",nrow(fits),"fits and",nrow(metrics),"alpha-specific rows into",output,"\n")
