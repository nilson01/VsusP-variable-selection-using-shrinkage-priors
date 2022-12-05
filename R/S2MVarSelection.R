
                ###########################
                ###########################
           ########   S2MVarSelection    #########
                ###########################
                ###########################

#'  Variable selection using shrinkage priors:: Sequential2Means
#'
#' S2MVarSelection function will take S2M: a list obtained from the 2Means.variables function and H: the estimated number of signals obtained from the optimal.b.i function. This will give out the important subset of variables for the Gaussian Linear model.
#'
#' @param Beta matrix consisting of N posterior samples of p variables that is known either to user or from Sequential2Means function
#' @param H Estimated number of signals obtained from the optimal.b.i function
#'
#' @return Indices of important subset of variables
#' @export
#'
#' @examples
#'
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(100), n, p)
#' beta <- exp(rnorm(p))
#' Y <- X %*% beta + rnorm(n, 0, 1)
#' df <- data.frame(X,Y)
#' rv.hs <- bayesreg::bayesreg(Y~. ,df, "gaussian", "horseshoe+", 10000, 2000)
#'
#' Beta = rv.hs$beta
#' H = 10
#' impVariablesGLM = S2MVarSelection(Beta, H)
#' impVariablesGLM
#'
#'
S2MVarSelection <- function(Beta, H) {
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




                ##########################
                ##########################
                ########   V1    #########
                ##########################
                ##########################


#'  Variable selection using shrinkage priors:: Sequential2Means
#'
#' S2MVarSelectionV1 function will take S2M: a list obtained from the 2Means.variables function and H: the estimated number of signals obtained from the optimal.b.i function. This will give out the important subset of variables for the Gaussian Linear model.
#'
#' @param S2M List obtained from the 2Means.variables function
#' @param H Estimated number of signals obtained from the optimal.b.i function
#'
#' @return Indices of important subset of variables for the Gaussian Linear model
#' @export
#'
S2MVarSelectionV1 <- function(S2M, H) {
  # number of covariates
  p <- ncol(S2M$p)

  # the medians of the absolute values of the posterior samples of each Beta vector
  abs.post.median <- S2M$abs.post.median

  # the indices of selected variables
  impVariablesGLM <- order(abs.post.median)[p:(p - H + 1)]

  return(impVariablesGLM)
}
