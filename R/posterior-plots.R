#' Plot posterior coefficient summaries
#' @param object A native Gaussian fit.
#' @param interval Equal-tail posterior interval probability.
#' @param selected Optional selected predictor indices, using R's one-based
#'   convention. These change display color, not the fitted summaries.
#' @param ... Additional arguments to the base graphics plot.
#' @return Invisibly, the displayed summary data frame. No model is fitted.
#' @section Example figure:
#' Posterior means and 95 percent equal-tail intervals for a six-observation ridge example.
#' \if{html}{\figure{coefficients.png}{Posterior means and 95 percent equal-tail intervals for a six-observation ridge example.}}
#' \if{latex}{\figure{coefficients.pdf}{options: width=5.1in}}
#' @export
#' @examples
#' X <- matrix(c(-1,0,1,2),ncol=1)
#' fit <- FitRidge(X,c(-1,0,1,2),draws=30,chains=2,seed=1)
#' PlotCoefficients(fit)
PlotCoefficients <- function(object,interval=.95,selected=integer(0),...) {
  checkVsuspFit(object)
  values <- summary(object,interval=interval)
  if(!is.numeric(selected) || is.complex(selected) || any(!is.finite(selected)) ||
     any(selected%%1!=0) || any(selected<1 | selected>nrow(values)))
    stop("selected must contain valid one-based predictor indices.")
  colors <- rep("#0072B2",nrow(values));colors[selected]<-"#D55E00"
  vsuspPlotArgs(graphics::plot,list(x=values$mean,y=seq_len(nrow(values)),
    xlim=range(values$lower,values$upper,values$mean),ylim=c(.5,nrow(values)+.5),
    yaxt="n",xlab=expression(beta[j]),ylab="",col=colors,pch=19),list(...))
  graphics::segments(values$lower,seq_len(nrow(values)),values$upper,seq_len(nrow(values)),col=colors)
  graphics::axis(2,at=seq_len(nrow(values)),labels=rownames(values),las=1)
  graphics::abline(v=0,lty=2,col="grey50")
  invisible(values)
}

#' Plot exact spike inclusion frequencies
#' @param object A result from [SpikeSlabSelection()] or a native spike-and-slab fit.
#' @param ... Additional arguments to barplot.
#' @return Invisibly, the displayed frequencies. The dashed line marks the
#'   median-probability-model cutoff; it is not an error-control threshold.
#' @section Example figure:
#' Inclusion frequencies from four fixed coefficient rows. The dashed line marks 0.5.
#' \if{html}{\figure{inclusion.png}{Inclusion frequencies from four fixed coefficient rows. The dashed line marks 0.5.}}
#' \if{latex}{\figure{inclusion.pdf}{options: width=5.1in}}
#' @export
PlotInclusion <- function(object,...) {
  if(inherits(object,"vsusp_fit")) {
    if(object$prior!="spike-slab")stop("Inclusion frequencies require exact spike-and-slab draws.")
    object<-SpikeSlabSelection(PosteriorDraws(object))
  }
  if(!inherits(object,"vsusp_spike"))stop("object must contain exact spike-and-slab summaries.")
  values<-object$inclusion.probability
  vsuspPlotArgs(graphics::barplot,list(height=values,horiz=TRUE,xlim=c(0,1),
    names.arg=vsuspNames(names(values),length(values)),
    col="#0072B2",xlab="Posterior inclusion frequency"),list(...))
  graphics::abline(v=.5,lty=2,col="grey40")
  invisible(values)
}

#' Plot chains without pooling their iterations
#' @param draws Native fit or iteration by chain by variable array.
#' @param variable One-based predictor index, or a fitted predictor name.
#' @param ... Additional arguments to matplot.
#' @return Invisibly, the iteration by chain matrix displayed.
#' @section Example figure:
#' First 300 retained draws from each of four horseshoe chains in the five-observation example. This excerpt is not a convergence assessment.
#' \if{html}{\figure{trace.png}{First 300 retained draws from each of four horseshoe chains in the five-observation example. This excerpt is not a convergence assessment.}}
#' \if{latex}{\figure{trace.pdf}{options: width=5.1in}}
#' @export
PlotTrace <- function(draws,variable=1L,...) {
  values<-vsuspPlotVariable(draws,variable)
  vsuspPlotArgs(graphics::matplot,list(x=seq_len(nrow(values$draws)),y=values$draws,
    type="l",xlab="Retained iteration",ylab=values$label),list(...))
  invisible(values$draws)
}

#' Plot within-chain sample autocorrelations
#' @param draws Native fit or iteration by chain by variable array.
#' @param variable One-based predictor index, or fitted name.
#' @param max.lag Nonnegative integer maximum lag, capped at draws minus one.
#' @param ... Additional arguments to matplot.
#' @return Invisibly, a lag by chain autocorrelation matrix. Constant chains
#'   have NA autocorrelations, rather than evidence of mixing.
#' @section Example figure:
#' Within-chain autocorrelations from all 3000 retained horseshoe draws per chain in the same example.
#' \if{html}{\figure{autocorrelation.png}{Within-chain autocorrelations from all 3000 retained horseshoe draws per chain in the same example.}}
#' \if{latex}{\figure{autocorrelation.pdf}{options: width=5.1in}}
#' @export
PlotAutocorrelation <- function(draws,variable=1L,max.lag=40L,...) {
  values<-vsuspPlotVariable(draws,variable); a<-values$draws
  if(nrow(a)<2L)stop("At least two retained draws are required.")
  checkIntegerScalar(max.lag,"max.lag");lag<-min(max.lag,nrow(a)-1L)
  out<-matrix(NA_real_,lag+1L,ncol(a))
  for(chain in seq_len(ncol(a))) {
    x<-a[,chain];scale<-max(abs(x));if(scale>0)x<-x/scale
    x<-x-mean(x);norm<-sum(x*x)
    if(norm==0)next
    out[,chain]<-vapply(0:lag,function(k)sum(x[seq_len(length(x)-k)]*x[k+seq_len(length(x)-k)])/norm,numeric(1))
  }
  vsuspPlotArgs(graphics::matplot,list(x=0:lag,y=out,type="l",ylim=c(-1,1),
    xlab="Lag",ylab="Autocorrelation"),list(...))
  invisible(out)
}

#' Plot a tuning path
#' @param x An S2M tuning summary.
#' @param ... Additional graphics arguments.
#' @return Invisibly, the displayed thresholds and counts, sorted by threshold.
#' @importFrom graphics plot
#' @section Example figure:
#' Signal counts at two thresholds for fixed seven-coefficient draws. The path does not select a threshold.
#' \if{html}{\figure{tuning.png}{Signal counts at two thresholds for fixed seven-coefficient draws. The path does not select a threshold.}}
#' \if{latex}{\figure{tuning.pdf}{options: width=5.1in}}
#' @export
plot.vsusp_tuning <- function(x,...) {
  i<-order(x$b.i,method="radix")
  values<-data.frame(threshold=x$b.i[i],count=x$H.b.i[i])
  vsuspPlotArgs(graphics::plot,list(x=values$threshold,y=values$count,type="b",
    xlab=expression(b),ylab=expression(H[b])),list(...))
  invisible(values)
}

#' Plot a native Gaussian fit
#' @param x Native Gaussian fit.
#' @param type One of coefficients, trace, acf, or inclusion.
#' @param ... Arguments passed to the selected plotting function.
#' @return Invisibly, the displayed numerical values.
#' @export
plot.vsusp_fit <- function(x,type=c("coefficients","trace","acf","inclusion"),...) {
  type<-match.arg(type)
  switch(type,coefficients=PlotCoefficients(x,...),trace=PlotTrace(x,...),
    acf=PlotAutocorrelation(x,...),inclusion=PlotInclusion(x,...))
}

vsuspPlotVariable <- function(value,variable) {
  a<-vsuspChainArray(value)
  labels<-if(inherits(value,"vsusp_fit"))value$feature.names else vsuspNames(dimnames(a)[[3]],dim(a)[3])
  if(is.character(variable) && length(variable)==1L && !is.na(variable)) variable<-match(variable,labels)
  checkIntegerScalar(variable,"variable",1)
  if(variable>dim(a)[3])stop("variable exceeds the predictor count.")
  list(draws=matrix(a[,,variable],nrow=dim(a)[1],ncol=dim(a)[2]),label=labels[variable])
}

vsuspPlotArgs <- function(fun,args,extras) {
  do.call(fun,utils::modifyList(args,extras))
}
