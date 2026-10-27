# Reproduce the bounded selection pilot; use a new, empty output directory.
library(VsusP)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L) stop("Usage: Rscript run-selection-pilot.R NEW_OUTPUT_DIRECTORY")
root <- args[1L]
if (dir.exists(root) && length(list.files(root,all.files=TRUE,no..=TRUE))) {
  stop("Output directory must be empty; existing results will not be reused.")
}
dir.create(root,recursive=TRUE,showWarnings=FALSE)
mode <- "pilot"
cells <- expand.grid(effect=c(0,.25,1),rho=c(0,.8))
cells$cell <- seq_len(nrow(cells))
write.csv(cells,file.path(root,"conditions.csv"),row.names=FALSE)
metrics <- function(selected,truth,p) {
  tp <- length(intersect(selected,truth)); fp <- length(setdiff(selected,truth))
  c(empty=as.integer(!length(selected)),fp=fp,fn=length(setdiff(truth,selected)),
    fdp=fp/max(1,length(selected)),exact=as.integer(setequal(selected,truth)),
    sensitivity=if(length(truth))tp/length(truth) else NA_real_)
}
run_case <- function(cell,replicate) {
  seed <- 700000+100*cell+replicate
  n <- 80L;p <- 10L
  set.seed(seed)
  rho <- cells$rho[cell]; effect <- cells$effect[cell]
  Sigma <- outer(seq_len(p),seq_len(p),function(a,b)rho^abs(a-b))
  X <- scale(matrix(rnorm(n*p),n,p) %*% chol(Sigma))
  X <- matrix(as.numeric(X),n,p,dimnames=list(NULL,paste0("x",seq_len(p))))
  beta <- c(effect,-effect,rep(0,p-2L))
  Y <- as.vector(X %*% beta+rnorm(n))
  truth <- which(beta!=0)
  warnings <- character(); started <- proc.time()[[3]]
  tryCatch({
    posterior <- withCallingHandlers(bayesreg::bayesreg(Y~.,data.frame(Y=Y,X),
                 model="gaussian",prior="horseshoe",n.samples=1000,burnin=500,
                 thin=5,n.cores=1),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})
    B <- t(posterior$beta);colnames(B) <- colnames(X)
    stopifnot(nrow(B)==1000L,ncol(B)==p,all(is.finite(B)))
    point <- colMeans(B)
    original <- SAVS(X,point)
    factors <- rep(c(.1,10),p/2)
    changed.X <- sweep(X,2,factors,"*")
    changed.beta <- point/factors
    changed <- SAVS(changed.X,changed.beta)
    recovered <- SAVS(sweep(changed.X,2,factors,"/"),changed.beta*factors)
    stopifnot(identical(recovered$selected,original$selected))
    stopifnot(max(abs(as.vector(X%*%point)-as.vector(changed.X%*%changed.beta)))<1e-10)
    block1 <- SAVS(X,colMeans(B[1:500,,drop=FALSE]))$selected
    block2 <- SAVS(X,colMeans(B[501:1000,,drop=FALSE]))$selected
    instability <- as.integer(!setequal(block1,block2))
    rows <- list()
    add <- function(method,selected) {
      row <- data.frame(cell=cell,replicate=replicate,seed=seed,method=method,
                        effect=effect,rho=rho,t(metrics(selected,truth,p)),
                        block.disagreement=if(method=="SAVS")instability else NA_integer_)
      rows[[length(rows)+1L]] <<- row
    }
    add("SAVS",original$selected);add("SAVS_changed_units",changed$selected)
    add("SAVS_recovered_units",recovered$selected)
    set.seed(seed+1L)
    legacy <- S2MSensitivity(B,c(.1,.5,1),"kmeans","legacy")
    exact <- S2MSensitivity(B,c(.1,.5,1),"exact","refine")
    for (j in seq_len(3)) {
      add(paste0("S2M_legacy_b",legacy$b.i[j]),legacy$selected[[j]])
      add(paste0("S2M_exact_refine_b",exact$b.i[j]),exact$selected[[j]])
    }
    list(status="ok",cell=cell,replicate=replicate,seed=seed,warnings=warnings,
         elapsed=proc.time()[[3]]-started,rows=do.call(rbind,rows),
         X=X,Y=Y,beta=beta,draws=B,legacy=legacy,exact=exact,
         SAVS=original,SAVS.changed=changed,scale.factors=factors)
  },error=function(e)list(status="failed",cell=cell,replicate=replicate,seed=seed,
                         warnings=warnings,error=conditionMessage(e),
                         elapsed=proc.time()[[3]]-started))
}
queue <- if(mode=="timing")data.frame(cell=c(1L,6L),replicate=c(1L,1L)) else
  expand.grid(cell=seq_len(6),replicate=seq_len(10))
for(i in seq_len(nrow(queue))) {
  c <- queue$cell[i];r <- queue$replicate[i]
  file <- file.path(root,sprintf("cell%02d-rep%02d.rds",c,r))
  if(!file.exists(file))saveRDS(run_case(c,r),file)
  result <- readRDS(file)
  cat(sprintf("cell=%d rep=%d status=%s elapsed=%.3f\n",c,r,result$status,result$elapsed))
  flush.console()
}
files <- list.files(root,pattern="^cell[0-9]+-rep[0-9]+[.]rds$",full.names=TRUE)
all <- lapply(files,readRDS)
status <- do.call(rbind,lapply(all,function(a)data.frame(cell=a$cell,replicate=a$replicate,
                    seed=a$seed,status=a$status,elapsed=a$elapsed,
                    warnings=paste(a$warnings,collapse=" | "),
                    error=if(is.null(a$error))"" else a$error)))
write.csv(status,file.path(root,paste0(mode,"-status.csv")),row.names=FALSE)
success <- Filter(function(a)a$status=="ok",all)
if(length(success))write.csv(do.call(rbind,lapply(success,`[[`,"rows")),
                             file.path(root,paste0(mode,"-results.csv")),row.names=FALSE)
cat("Loaded:",find.package("VsusP"),"version",as.character(packageVersion("VsusP")),"\n")
print(sessionInfo())

if (any(status$status != "ok")) stop("Failed cases are retained in pilot-status.csv; no failed case was replaced.")
script.file <- sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[1L])
source(file.path(dirname(script.file),"summarize-selection-pilot.R"))
summarize_selection_pilot(root)
