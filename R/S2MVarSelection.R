
################################################################
#------------------ S2MVarSelection algorithm -----------------#
################################################################

#'  Variable selection using shrinkage priors :: S2MVarSelection
#'
#' S2MVarSelection function will take S2M: a list obtained from the 2Means.variables function and H: the estimated number of signals obtained from the optimal.H.b.i function. This will give out the important subset of variables for the Gaussian Linear model.
#'
#' @export
#' @param Beta matrix consisting of N posterior samples of p variables that is known either to user or from Sequential2Means function
#' @param H Estimated number of signals obtained from the optimal.b.i function of numeric data type
#'
#' @return a vector containing indices of important subset of variables of dimension H X 1.
#'
#' @references
#'
#' Makalic, E. & Schmidt, D. F.
#' High-Dimensional Bayesian Regularised Regression with the BayesReg Package
#' arXiv:1611.06649, 2016
#'
#' Li, H., & Pati, D.
#' Variable selection using shrinkage priors
#' Computational Statistics & Data Analysis, 107, 107-119.
#'
#' @examples
#'
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- exp(rnorm(p))
#' Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
#' df <- data.frame(X, Y)
#' # Fit a model using gaussian horseshoe+ for 200 samples
#' # # recommended n.samples is 5000 and burning is 2000
#' rv.hs <- bayesreg::bayesreg(Y ~ ., df, "gaussian", "horseshoe+", 200, 100)
#'
#' Beta <- rv.hs$beta
#' H <- 12
#' impVariablesGLM <- S2MVarSelection(Beta, H)
#' impVariablesGLM
#'
S2MVarSelection <- function(Beta, H = 10) {
  # Check for NULL or NaN values in H
  if (is.null(H)) {
    stop(paste("H must not be NULL. "))
  }

  # Check for numeric data type of H
  if (!is.numeric(H)) {
    stop(paste("H must be passed as numeric data type. "))
  }

  # Check for NULL or NaN values in Beta
  if (is.null(Beta) || any(is.na(Beta))) {
    stop(paste("Beta must not be NULL or have NaN values. \n Please check Sequential2Means funnction if Beta is not available."))
  }

  # Check for matrix data type of Beta
  if (!is.matrix(Beta)) {
    stop(paste("Beta must be passed as matrix data type. "))
  }

  # Check for H must be less than total covariates.
  if (H > ncol(Beta)) {
    stop(paste("H must be less than total covariates. "))
  }

  # number of covariates
  p <- ncol(Beta)

  # the medians of the absolute values of the posterior samples of each Beta vector
  abs.post.median <- seq(0, length = p)
  for (i in 1:p) {
    abs.post.median[i] <- stats::median(abs(Beta[, i]))
  }

  # the indices of selected variables
  impVariablesGLM <- order(abs.post.median)[p:(p - H + 1)]

  return(impVariablesGLM)
}



################################################################
#-------------------- OptimalHbi algorithm --------------------#
################################################################


#'  Variable selection using shrinkage priors :: OptimalHbi
#'
#' OptimalHbi function will take b.i and H.b.i as input which comes from the result of TwoMeans function. It will return plot from which you can infer about H: the optimal value of the tuning parameter.
#'
#' @export
#' @param bi a vector holding the values of the tuning parameter specified by the user
#' @param Hbi The estimated number of signals corresponding to each b.i of numeric data type
#'
#' @return the optimal value (numeric) of tuning parameter and the associated H value
#'
#' @references
#'
#' Makalic, E. & Schmidt, D. F.
#' High-Dimensional Bayesian Regularised Regression with the BayesReg Package
#' arXiv:1611.06649, 2016
#'
#' Li, H., & Pati, D.
#' Variable selection using shrinkage priors
#' Computational Statistics & Data Analysis, 107, 107-119.
#'
#' @examples
#'
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n * p), n, p)
#' beta <- exp(rnorm(p))
#' Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
#' df <- data.frame(X, Y)
#' rv.hs <- bayesreg::bayesreg(Y ~ ., df, "gaussian", "horseshoe+", 200, 100)
#'
#' Beta <- t(rv.hs$beta)
#' lower <- 0
#' upper <- 1
#' l <- 5
#' S2Mbeta <- Sequential2MeansBeta(Beta, lower, upper, l)
#'
#' bi <- S2Mbeta$b.i
#' Hbi <- S2Mbeta$H.b.i
#' OptimalHbi(bi, Hbi)
#'
OptimalHbi <- function(bi, Hbi) {
  # Check for NULL or NaN values in b.i
  if (is.null(bi) || any(is.na(bi))) {
    stop(paste("b.i must not be NULL or have NaN values."))
  }

  # Check for vector data type of b.i
  if (!is.vector(bi)) {
    stop(paste("b.i must be passed as vector data type."))
  }

  # Check for the length of b.i and H.b.i
  if (length(bi) != length(Hbi)) {
    stop(paste("b.i and H.b.i must have same length."))
  }

  # plotting tuning parameters Vs number of important variables counts
  plot(bi, Hbi, type = 'b', xlab = "Tuning parameter b.i", ylab = "Estimated number of signals (H.b.i)",
       main = "Plot of b.i vs H.b.i")
}

################################################################
#---------------- S2MVarSelectionV1 algorithm -----------------#
################################################################

#'  Variable selection using shrinkage priors :: S2MVarSelectionV1
#'
#' S2MVarSelectionV1 function will take S2M: a list obtained from the 2Means.variables function and H: the estimated number of signals obtained from the optimal.b.i function. This will give out the important subset of variables for the Gaussian Linear model.
#'
#' @export
#' @param S2M List obtained from the 2Means.variables function
#' @param H Estimated number (numeric) of signals obtained from the optimal.b.i function
#'
#' @return a vector of indices of important subset of variables for the Gaussian Linear modelof shape H X 1
#'
S2MVarSelectionV1 <- function(S2M, H = 10) {
  # number of covariates
  p <- S2M$p

  # the medians of the absolute values of the posterior samples of each Beta vector
  abs.post.median <- S2M$abs.post.median

  # Check for NULL or NaN values in H
  if (is.null(H)) {
    stop(paste("H must not be NULL. "))
  }

  # Check for numeric data type of H
  if (!is.numeric(H)) {
    stop(paste("H must be passed as numeric data type. "))
  }

  # Check for NULL or NaN values in abs.post.median
  if (is.null(abs.post.median) || any(is.na(abs.post.median))) {
    stop(paste("abs.post.median must not be NULL or have NaN values. "))
  }

  # Check for non negative values in abs.post.median
  if (any(abs.post.median < 0)) {
    stop(paste("abs.post.median must not be less than zero. "))
  }

  # Check for matrix data type of abs.post.median
  if (!is.vector(abs.post.median)) {
    stop(paste("abs.post.median must be passed as vector data type. "))
  }

  # Check for H must be less than total covariates.
  if (H > p) {
    stop(paste("H must be less than total covariates. "))
  }


  # the indices of selected variables
  impVariablesGLM <- order(abs.post.median)[p:(p - H + 1)]

  return(impVariablesGLM)
}

