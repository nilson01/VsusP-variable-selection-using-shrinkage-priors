#' Pool retained coefficient draws for selection
#' @param object A native Gaussian fit returned by [FitGaussian()].
#' @return Numeric matrix with draws in rows and predictors in columns. Chains
#'   are combined only in this returned matrix; the fit retains its chain axes.
#' @export
#' @examples
#' X <- matrix(c(-1,0,1,2), ncol=1)
#' fit <- FitRidge(X,c(-1,0,1,2),draws=20,chains=2,seed=3)
#' dim(PosteriorDraws(fit))
PosteriorDraws <- function(object) {
  checkVsuspFit(object)
  matrix(object$beta,ncol=dim(object$beta)[3],
         dimnames=list(NULL,object$feature.names))
}

#' Summarize native Gaussian posterior draws
#' @param object A native Gaussian fit.
#' @param interval Posterior interval probability, strictly between zero and one.
#' @param ... Further arguments, currently unused.
#' @return Data frame with posterior mean, SD, lower and upper equal-tail
#'   limits. SD is NA when only one draw is available.
#' @export
summary.vsusp_fit <- function(object,interval=.95,...) {
  vsuspInterval(interval); b <- PosteriorDraws(object)
  bounds <- apply(b,2L,stats::quantile,probs=c((1-interval)/2,(1+interval)/2),names=FALSE)
  result <- data.frame(mean=colMeans(b),sd=apply(b,2L,stats::sd),
                      lower=bounds[1,],upper=bounds[2,],row.names=object$feature.names)
  attr(result,"interval") <- interval
  result
}

#' @export
print.vsusp_fit <- function(x,...) {
  checkVsuspFit(x)
  cat("Gaussian",x$prior,"fit;",dim(x$beta)[1],"retained draws per chain;",
      dim(x$beta)[2],"chains\n")
  print(summary(x),...); invisible(x)
}

#' Posterior coefficient means
#' @param object A native Gaussian fit.
#' @param ... Further arguments, currently unused.
#' @return Named numeric vector of posterior coefficient means, excluding
#'   the intercept. Intercept draws remain in object$intercept.
#' @export
coef.vsusp_fit <- function(object,...) colMeans(PosteriorDraws(object))

#' Predict from native Gaussian posterior draws
#' @param object A native Gaussian fit.
#' @param newdata Numeric predictor matrix in the fitted column order. If
#'   column names are supplied, they must match the fitted names and order.
#' @param interval Posterior interval probability.
#' @param include.noise Logical: include independent observation noise in the
#'   interval draws. The reported mean remains the posterior conditional mean.
#' @param seed Optional nonnegative prediction seed; preserves the caller RNG.
#' @param ... Further arguments, currently unused.
#' @return List with mean, lower, upper, interval and kind. Intervals are
#'   pointwise posterior intervals, not simultaneous bands. kind is "mean"
#'   or "observation". Prediction does not refit the model.
#' @importFrom stats predict coef
#' @export
predict.vsusp_fit <- function(object,newdata,interval=.95,include.noise=FALSE,
                              seed=NULL,...) {
  vsuspInterval(interval); b <- PosteriorDraws(object)
  vsuspPredictMatrix(newdata,object$feature.names)
  if (!is.logical(include.noise) || length(include.noise)!=1L || is.na(include.noise))
    stop("include.noise must be TRUE or FALSE.")
  values <- sweep(b%*%t(newdata),1L,as.vector(object$intercept),"+")
  vsuspFinite(values,"Prediction")
  means <- colMeans(values)
  if (include.noise) values <- vsuspWithSeed(seed,function()
    values+matrix(stats::rnorm(length(values),sd=sqrt(as.vector(object$sigma2))),nrow(values)))
  vsuspFinite(values,"Prediction")
  bounds <- apply(values,2L,stats::quantile,
                  probs=c((1-interval)/2,(1+interval)/2),names=FALSE)
  list(mean=means,lower=unname(bounds[1,]),upper=unname(bounds[2,]),
       interval=interval,kind=if(include.noise) "observation" else "mean")
}

#' Summarize sparse coefficient results
#' @param object A result from [SAVS()] or [DSS()].
#' @param ... Further arguments, currently unused.
#' @return Data frame with coefficient values and a selected indicator.
#' @export
summary.vsusp_sparse <- function(object,...) {
  data.frame(coefficient=object$coefficients,selected=object$coefficients!=0,
             row.names=vsuspNames(names(object$coefficients),length(object$coefficients)))
}

#' @export
print.vsusp_sparse <- function(x,...) {
  cat(if(!is.null(x$lambda)) "DSS" else "SAVS",":",length(x$selected),"selected\n")
  print(summary(x),...); invisible(x)
}

#' Predict using sparse coefficient summaries
#' @param object A result from [SAVS()] or [DSS()].
#' @param newdata Numeric matrix in the coefficient order.
#' @param intercept Finite intercept added to the sparse linear prediction.
#' @param ... Further arguments, currently unused.
#' @return Numeric vector. This method returns sparse linear predictions;
#'   it supplies no uncertainty interval and does not refit a model.
#' @export
predict.vsusp_sparse <- function(object,newdata,intercept=0,...) {
  checkNumericScalar(intercept,"intercept")
  labels <- vsuspNames(names(object$coefficients),length(object$coefficients))
  vsuspPredictMatrix(newdata,labels)
  values <- as.vector(newdata%*%object$coefficients)+intercept
  vsuspFinite(values,"Sparse prediction");values
}

#' Summarize an S2M tuning path
#' @param object An S2M sensitivity or coefficient-draw summary.
#' @param ... Further arguments, currently unused.
#' @return Data frame of thresholds and modal signal counts, with modal
#'   frequencies when those were retained. This does not choose a threshold.
#' @export
summary.vsusp_tuning <- function(object,...) {
  out <- data.frame(threshold=object$b.i,count=object$H.b.i)
  if (!is.null(object$mode.proportion)) out$modal.frequency <- object$mode.proportion
  out
}

#' @export
print.vsusp_tuning <- function(x,...) {
  print(summary(x),row.names=FALSE,...);invisible(x)
}

#' Summarize basic two-means selection
#' @param object A result from [TwoMeansBeta()].
#' @param ... Further arguments, currently unused.
#' @return Named selected indices, preserving the fitted ranking.
#' @export
summary.vsusp_selection <- function(object,...) {
  labels <- vsuspNames(names(object$abs.post.median),object$p)
  stats::setNames(object$selected,labels[object$selected])
}

#' @export
print.vsusp_selection <- function(x,...) {
  cat("Two-means:",length(x$selected),"selected\n");print(summary(x),...);invisible(x)
}

#' Summarize exact spike-and-slab inclusion draws
#' @param object A result from [SpikeSlabSelection()].
#' @param ... Further arguments, currently unused.
#' @return Data frame of inclusion frequencies and median-probability-model
#'   membership. Clustering frequencies are not accepted as inclusion draws.
#' @export
summary.vsusp_spike <- function(object,...) {
  probability <- object$inclusion.probability
  data.frame(inclusion.frequency=probability,MPM=probability>=.5,
             row.names=vsuspNames(names(probability),object$p))
}

#' @export
print.vsusp_spike <- function(x,...) {
  print(summary(x),...)
  cat("MPM:",length(x$MPM),"selected; most frequent visited model:",length(x$HPPM),"selected\n")
  invisible(x)
}

vsuspInterval <- function(interval) {
  checkNumericScalar(interval,"interval")
  if(interval<=0 || interval>=1) stop("interval must be strictly between zero and one.")
  invisible(interval)
}

vsuspPredictMatrix <- function(X,labels) {
  checkNumericMatrix(X,"newdata")
  if(ncol(X)!=length(labels)) stop("newdata must preserve the fitted predictor count and order.")
  if(!is.null(colnames(X)) && !identical(colnames(X),labels))
    stop("newdata column names must preserve the fitted predictor order.")
  invisible(X)
}
