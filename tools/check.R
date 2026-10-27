# Run against a built archive, not a source checkout.
main <- function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=3L || !args[2] %in% c("quick","docs","full"))
    stop("Usage: Rscript check.R PACKAGE.tar.gz quick|docs|full RESULTS_DIRECTORY")
  archive <- normalizePath(args[1],mustWork=TRUE)
  if(!grepl("\\.tar\\.gz$",archive)) stop("Supply an R source-package tar.gz archive.")
  mode <- args[2]; output <- args[3]
  if(dir.exists(output) && length(list.files(output,all.files=TRUE,no..=TRUE)))
    stop("Use an empty results directory; previous evidence will not be overwritten.")
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  output <- normalizePath(output,mustWork=TRUE)
  meta <- list(status="incomplete",mode=mode,platform=R.version$platform,
    OS=Sys.info()[["sysname"]],R=as.character(getRversion()),
    archive=basename(archive),archive.MD5=unname(tools::md5sum(archive)),
    remote.checks.enabled=mode=="full",documentation.checks.enabled=mode!="quick")
  on.exit(saveRDS(meta,file.path(output,"result.rds")),add=TRUE)
  tryCatch({
  required <- c("bayesreg","testthat","posterior","knitr","rmarkdown","tinytex","covr","MASS")
  if(mode!="quick")required<-c(required,"V8")
  missing<-required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop("Install check dependencies first: ",paste(missing,collapse=", "))
  meta$dependencies<-vapply(required,function(x)as.character(utils::packageVersion(x)),character(1))
  if(mode!="quick" && !nzchar(Sys.which("pdflatex")))stop("Full checking requires pdflatex on PATH.")
  profile<-file.path(output,"check-profile.R")
  writeLines("options(timeout=10)",profile)
  variables<-c("_R_CHECK_LIMIT_CORES_","_R_CHECK_CRAN_INCOMING_",
               "_R_CHECK_CRAN_INCOMING_REMOTE_","_R_CHECK_SYSTEM_CLOCK_","R_PROFILE_USER",
               "_R_CHECK_FORCE_SUGGESTS_")
  old<-Sys.getenv(variables,unset=NA_character_)
  on.exit({for(i in seq_along(variables)) {
    if(is.na(old[i]))Sys.unsetenv(variables[i]) else
      do.call(Sys.setenv,stats::setNames(list(old[i]),variables[i]))
  }},add=TRUE)
  Sys.setenv(`_R_CHECK_LIMIT_CORES_`="true",`_R_CHECK_FORCE_SUGGESTS_`="true",R_PROFILE_USER=profile)
  for(name in variables[2:4]) do.call(Sys.setenv,stats::setNames(list(if(mode=="full")"true" else "false"),name))
  flags<-"--as-cran"
  if(mode=="quick") flags<-c(flags,"--no-manual","--ignore-vignettes","--no-build-vignettes")
  meta$flags<-flags
  writeLines(capture.output(sessionInfo()),file.path(output,"session.txt"))
  script<-sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)[1])
  script<-gsub("~+~"," ",script,fixed=TRUE) # R encodes spaces in --file arguments.
  helper<-file.path(dirname(normalizePath(script,mustWork=TRUE)),"check-results.R")
  old.wd<-getwd();on.exit(setwd(old.wd),add=TRUE);setwd(output)
  executable<-file.path(R.home("bin"),if(.Platform$OS.type=="windows")"R.exe" else "R")
  status<-system2(executable,c("CMD","check",flags,shQuote(archive)),
    stdout="command.log",stderr="command.log",timeout=240)
  meta$command.status<-status
  if(status!=0L)stop("Package check command failed; inspect command.log.")
  source(helper,local=TRUE)
  verified<-verify_check_results(output,documentation=mode!="quick")
  meta$notes<-verified$notes;meta$assertions<-verified$assertions
  counts<-verified$assertions
  meta$status<-"passed"
  cat("PASS:",mode,"check on",meta$platform,"with",counts[4],"installed assertions.\n")
  cat(length(meta$notes),"note lines retained in result.rds; inspect them before release.\n")
  },error=function(e) {meta$error <<- conditionMessage(e);stop(e)})
}
main()
