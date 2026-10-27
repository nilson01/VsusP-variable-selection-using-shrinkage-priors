# Aggregate complete pilot cases; independent datasets are the sampling unit.
summarize_selection_pilot <- function(root) {
  status <- utils::read.csv(file.path(root,"pilot-status.csv"))
  stopifnot(nrow(status)==60L,!anyDuplicated(status[c("cell","replicate")]))
  if (any(status$status != "ok")) stop("Incomplete pilot: inspect all retained failures before summarizing.")
  raw <- utils::read.csv(file.path(root,"pilot-results.csv"),check.names=FALSE)
  stopifnot(nrow(raw)==540L)
  z <- stats::qnorm(.975)
  rows <- list()
  for(cell in seq_len(6)) for(method in sort(unique(raw$method))) {
    x <- raw[raw$cell==cell & raw$method==method,]
    stopifnot(nrow(x)==10L,!anyDuplicated(x$replicate))
    row <- list(cell=cell,effect=x$effect[1],rho=x$rho[1],method=method,
                n_planned=10,n_success=nrow(x),failures=0)
    for(metric in c("empty","fp","fn","fdp","exact","sensitivity","block.disagreement")) {
      a <- x[[metric]];a <- a[!is.na(a)]
      if(!length(a)) next
      n <- length(a);p <- mean(a)
      row[[paste0(metric,".mean")]] <- p
      row[[paste0(metric,".mcse")]] <- stats::sd(a)/sqrt(n)
      if(metric %in% c("empty","exact","block.disagreement")) {
        den <- 1+z^2/n;mid <- (p+z^2/(2*n))/den
        half <- z*sqrt(p*(1-p)/n+z^2/(4*n^2))/den
        row[[paste0(metric,".lower")]] <- max(0,mid-half)
        row[[paste0(metric,".upper")]] <- min(1,mid+half)
        row[[paste0(metric,".failure_lower")]] <- sum(a)/10
        row[[paste0(metric,".failure_upper")]] <- (sum(a)+10-n)/10
      }
    }
    rows[[length(rows)+1L]] <- row
  }
  fields <- unique(unlist(lapply(rows,names)))
  rows <- lapply(rows,function(r){r[setdiff(fields,names(r))]<-NA;as.data.frame(r[fields])})
  result <- do.call(rbind,rows)
  utils::write.csv(result,file.path(root,"summary.csv"),row.names=FALSE)
  invisible(result)
}
