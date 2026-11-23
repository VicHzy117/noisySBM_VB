# Verify the delivered data/code/results after copying or downloading.
script<-normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]))
setwd(dirname(dirname(script)))
if(!requireNamespace("digest",quietly=TRUE))stop("Install digest")
manifest<-"FILES_SHA256.csv"
hash<-function(p)digest::digest(file=p,algo="sha256",serialize=FALSE)
if("--write" %in% commandArgs(TRUE)) {
  files<-sort(list.files(".",recursive=TRUE,all.files=TRUE,full.names=FALSE))
  files<-files[!dir.exists(files) & files!=manifest &
    !grepl("^([.]git/|[.]cache/|results/data/(reruns|new_summary)/)",files) & files!="tests/local_session.txt"]
  stopifnot(!any(grepl("[.](py|pyc|so|dylib|dll|o|rtf|pem)$",files)),
    !any(nzchar(Sys.readlink(files))))
  sizes<-file.info(files)$size
  stopifnot(all(sizes<100*1024^2))
  write.csv(data.frame(path=files,bytes=sizes,sha256=vapply(files,hash,character(1))),manifest,row.names=FALSE)
  cat("Wrote",length(files),"SHA-256 records; all files below 100 MiB.\n")
} else {
  m<-read.csv(manifest,stringsAsFactors=FALSE)
  stopifnot(all(file.exists(m$path)),all(file.info(m$path)$size==m$bytes),
    identical(unname(vapply(m$path,hash,character(1))),m$sha256))
  cat("PASS:",nrow(m),"file hashes.\n")
}
