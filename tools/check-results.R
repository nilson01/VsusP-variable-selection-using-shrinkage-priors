verify_check_results <- function(directory,documentation=FALSE) {
  logs<-list.files(directory,pattern="^00check\\.log$",recursive=TRUE,full.names=TRUE)
  if(length(logs)!=1L)stop("Expected one package-check log.")
  log<-readLines(logs,warn=FALSE)
  if(!any(grepl("^\\* DONE",log)) || any(grepl("\\.\\.\\. (ERROR|WARNING)",log)))
    stop("Package checking is incomplete or reported errors/warnings.")
  rout<-list.files(directory,pattern="^testthat\\.Rout$",recursive=TRUE,full.names=TRUE)
  if(length(rout)!=1L)stop("The installed test suite did not produce its expected report.")
  text<-readLines(rout,warn=FALSE)
  expression<-"\\[ FAIL ([0-9]+) \\| WARN ([0-9]+) \\| SKIP ([0-9]+) \\| PASS ([0-9]+) \\]"
  counts<-regmatches(text,regexec(expression,text));counts<-counts[lengths(counts)==5L]
  if(length(counts)!=1L)stop("Could not verify the installed test counts.")
  counts<-as.integer(counts[[1]][-1]);names(counts)<-c("fail","warn","skip","pass")
  if(any(counts[1:3]!=0L) || counts[4]<1L)
    stop("The suite failed, warned, skipped or ran no assertions.")
  if(documentation) {
    required<-c("PDF version of manual","HTML version of manual","re-building of vignette outputs")
    for(item in required)
      if(!any(grepl(paste0("checking ",item,".*[.][.][.] OK$"),log)))
        stop("Full documentation check did not pass: ",item)
  }
  list(assertions=counts,notes=grep("NOTE",log,value=TRUE))
}

if(sys.nframe()==0L) {
  args<-commandArgs(trailingOnly=TRUE)
  if(length(args)<1L || length(args)>2L || (length(args)==2L && args[2]!="full"))
    stop("Usage: Rscript check-results.R CHECK_DIRECTORY [full]")
  result<-verify_check_results(args[1],documentation=length(args)==2L)
  cat("PASS:",result$assertions["pass"],"installed assertions; no failures, warnings or skips.\n")
}
