#' Calculate diagnostics with explicit chain identities
#'
#' Delegate rank R-hat, bulk/tail effective sample size and mean Monte Carlo
#' standard error to the optional posterior package.
#' @param draws A native fit, or finite numeric array with iteration, chain,
#'   variable axes. A pooled matrix is not accepted.
#' @param feature.names Unique variable names for an array. A fit supplies its
#'   own names and includes intercept and residual-variance diagnostics.
#' @return Data frame with rhat, ess.bulk, ess.tail, mcse.mean and status. A
#'   single chain has unavailable R-hat. Constant variables have unavailable
#'   diagnostics. Finite diagnostic values do not certify convergence.
#' @details At least four retained iterations per chain are required. Install
#'   posterior separately when using this function. No chain is reconstructed
#'   from a pooled draw matrix. R and Python diagnostics can differ slightly
#'   between backend versions; both use the cited rank-based definitions.
#' @references Vehtari et al. (2021). Rank-normalization, folding, and
#'   localization: An improved R-hat for assessing convergence of MCMC.
#'   \doi{10.1214/20-BA1221}.
#' @export
#' @examples
#' if (requireNamespace("posterior",quietly=TRUE)) {
#'   a <- array(sin(seq_len(800)),c(100,2,4))
#'   ChainDiagnostics(a)
#' }
ChainDiagnostics <- function(draws,feature.names=NULL) {
  a <- vsuspChainArray(draws)
  if(dim(a)[1]<4L) stop("Diagnostics require at least four retained draws per chain.")
  if(!requireNamespace("posterior",quietly=TRUE))
    stop("Install the optional posterior package to use ChainDiagnostics.",call.=FALSE)
  if(inherits(draws,"vsusp_fit")) {
    if(!is.null(feature.names)) stop("A fit already supplies variable names.")
    names <- c(paste0("beta[",draws$feature.names,"]"),"intercept","sigma2")
    expanded <- array(NA_real_,c(dim(a)[1:2],length(names)))
    expanded[,,seq_len(dim(a)[3])] <- a
    expanded[,,length(names)-1L] <- draws$intercept
    expanded[,,length(names)] <- draws$sigma2
    a <- expanded
  } else names <- vsuspNames(feature.names,dim(a)[3])
  out <- data.frame(rhat=rep(NA_real_,dim(a)[3]),ess.bulk=NA_real_,
                    ess.tail=NA_real_,mcse.mean=NA_real_,status="inspect diagnostics",
                    row.names=names)
  for(j in seq_len(dim(a)[3])) {
    x <- matrix(a[,,j],nrow=dim(a)[1],ncol=dim(a)[2])
    if(all(x==x[1,1])) {out$status[j]<-"constant; unavailable";next}
    if(ncol(x)>=2L) out$rhat[j]<-posterior::rhat(x)
    else out$status[j]<-"one chain; R-hat unavailable"
    out$ess.bulk[j]<-posterior::ess_bulk(x)
    out$ess.tail[j]<-posterior::ess_tail(x)
    out$mcse.mean[j]<-posterior::mcse_mean(x)
    values<-unlist(out[j,c("ess.bulk","ess.tail","mcse.mean")],use.names=FALSE)
    if(ncol(x)>=2L)values<-c(values,out$rhat[j])
    if(any(!is.finite(values)))out$status[j]<-paste(out$status[j],"some diagnostics unavailable",sep="; ")
  }
  attr(out,"backend")<-paste("posterior",utils::packageVersion("posterior"))
  out
}

vsuspChainArray <- function(value) {
  if(inherits(value,"vsusp_fit")) {checkVsuspFit(value);return(value$beta)}
  if(!is.array(value) || !is.numeric(value) || is.complex(value) ||
     length(dim(value))!=3L || any(dim(value)<1L) || any(!is.finite(value)))
    stop("draws must retain finite iteration, chain and variable axes.")
  value
}
