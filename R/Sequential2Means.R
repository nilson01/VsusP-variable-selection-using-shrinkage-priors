#' Variable selection using shrinkage priors:: Sequential2Means
#'
#' @param Beta.i ith row sample of N by p matrix consisting of N posterior samples of p variables
#' @param b.i_r tuning parameter value from Sequential 2-means (S2M) variable selection algorithm.
#'
#' @return number of noise coefficients
#'
numNoiseCoeff <- function(Beta.i, b.i_r) {
  # perform k means for sequential 2 means to prune the noise coefficients
  fit <- stats::kmeans(abs(Beta.i), 2)
  # perform sequential 2 means on subsequent clusters until threshold is reached

  # check the threshold between two centers
  while (max(fit$centers) - min(fit$centers) > b.i_r) {
    Beta.i <- Beta.i[which(fit$cluster == which.min(fit$centers))]
    if (length(Beta.i) <= 2) {
      break
    }
    fit <- stats::kmeans(abs(Beta.i), 2)
  }

  argminCenter <- which.min(fit$centers)
  return(length(which(fit$cluster == argminCenter)))
}


                                  ###########################
                                  ###########################
                           #########   Sequential2Means    #########
                                  ###########################
                                  ###########################


#' Variable selection using shrinkage priors:: Sequential2Means
#'
#' Sequential2Means function will take as input X: design matrix, Y : response vector, t: vector of tuning parameter values from Sequential 2-means (S2M) variable selection algorithm. The function will return a list S2M which will hold p: the total number of variables, b.i: the values of the tuning parameter, H.b.i : the estimated number of signals corresponding to each b.i, abs.post.median: medians of the absolute values of the posterior samples.
#'
#' @param X Design matrix
#' @param Y Response vector
#' @param b.i Vector of tuning parameter values from Sequential 2-means (S2M) variable selection algorithm.
#' @param prior Shrinkage prior distribution over the Beta. Options: ridge regression("ridge"), lasso regression("lasso"), horseshoe regression ("horseshoe") and horseshoe+ regression ("horseshoe+")
#' @param n.samples Number of posterior samples to generate.
#' @param burnin Number of burn-in samples.
#'
#' @return A list S2M which will hold Beta: N by p matrix consisting of N posterior samples of p variables, b.i: the values of the tuning parameter, H.b.i : the estimated number of signals corresponding to each b.i, abs.post.median: medians of the absolute values of the posterior samples.
#' @export
#'
#' @examples
#'
# Example One Gaussian Model and Horseshoe prior
#'
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' beta <- exp(rnorm(p))
#' Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
#' b.i = seq(0, 1, 0.05)
#' S2M = Sequential2Means(X, Y, b.i, "horseshoe+", 200, 100)
#' Beta = S2M$Beta
#' H.b.i = S2M$H.b.i
#'
#'
# Example Two Gaussian Model and ridge prior
#'
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' beta <- exp(rnorm(p))
#' Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
#' b.i = seq(0, 1, 0.05)
#' S2M = Sequential2Means(X, Y, b.i, "ridge", 200, 100)
#' Beta = S2M$Beta
#' H.b.i = S2M$H.b.i
#'

Sequential2Means <- function(X, Y, b.i, prior="horseshoe+", n.samples = 5000, burnin = 2000) {

  # Check for NULL or NaN values in X
  if(is.null(X) || any(is.na(X))){
    stop(paste("X must not be NULL or have NaN values."))
  }

  # Check for matrix data type of X
  if(! is.matrix(X)){
    stop(paste("X must be passed as matrix data type."))
  }

  # Check for NULL or NaN values in Y
  if(is.null(Y) || any(is.na(Y))){
    stop(paste("Y must not be NULL or have NaN values."))
  }

  # Check for vector data type of Y
  if(! is.vector(Y)){
    stop(paste("Y must be passed as vector data type."))
  }

  # Check for compatibility of dimensions between X and Y
  if (nrow(X) != length(Y)) {
    stop(paste("Dimensions between X and Y are not compatible. "))
  }


  # Check for NULL or NaN values in b.i
  if(is.null(b.i) || any(is.na(b.i))){
    stop(paste("b.i must not be NULL or have NaN values."))
  }

  # Check for vector data type of b.i
  if(! is.vector(b.i)){
    stop(paste("b.i must be passed as vector data type."))
  }

  # Check for prior from available options
  if(! prior %in%  c('ridge', 'lasso', 'horseshoe', 'horseshoe+')){
    stop(paste("Available prior options: ridge regression('ridge'), lasso regression('lasso'), horseshoe regression ('horseshoe') and horseshoe+ regression ('horseshoe+')"))
  }

  # Check for the lower bound of n.samples
  if(n.samples < 100){
    stop(paste("n.samples is recommended to be at least 100"))
  }

  # Check if n.samples is a natural number
  if(n.samples%%1!=0){
    stop(paste("n.samples must be a natural number greater than or equal to 100"))
  }

  # Check for the lower bound of burnin
  if(burnin < 100){
    stop(paste("burnin is recommended to be at least 100"))
  }

  # Check if the number of burnin is a natural number
  if(burnin%%1!=0){
    stop(paste("burnin must be a natural number greater than or equal to 100"))
  }




  # Initializing variables
  ##########################

  # the posterior sample size
  N <- n.samples

  # number of covariates
  p <- ncol(X)

  # number of tuning parameters
  l <- length(b.i)

  # initializing the number of signals corresponding to each b.i
  H.b.i <- rep(0, length = length(b.i))

  # N by p matrix consisting of N posterior samples of p variables

  # Using MCMC for beta sampling
  fit <- bayesreg::bayesreg(Y ~ . ,data.frame(X,Y) , model="gaussian", prior, n.samples, burnin)
  Beta = t(fit$beta)

  # Sequential 2 means algorithm
  for (r in 1:l) {
    impVars.count <- seq(0, length = N)
    for (i in 1:N) {
      impVars.count[i] <- p - numNoiseCoeff(Beta[i, ], b.i[r])
    }
    H.b.i[r] <- as.numeric(names(sort(-table(impVars.count)))[1])
  }

  # list to hold desired output
  S2M <- list(Beta = Beta, H.b.i = H.b.i)

  return(S2M)
}


                          ##########################
                          ##########################
                  ########   Sequential2MeansBeta    #########
                          ##########################
                          ##########################

#' Variable selection using shrinkage prior:: Sequential2Means
#'
#' Sequential2MeansBeta function will take as input Beta : N by p matrix consisting of N posterior samples of p variables, lower : the lower bound of the chosen values of the tuning parameter, upper : the upper bound of the chosen values of the tuning parameter, and l :the number of chosen values of the tuning parameter. The function will return a list S2M which will hold p: the total number of variables, b.i: the values of the tuning parameter, H.b.i : the estimated number of signals corresponding to each b.i, abs.post.median: medians of the absolute values of the posterior samples.
#'
#' @param Beta N by p matrix consisting of N posterior samples of p variables
#' @param lower the lower bound of the chosen values of the tuning parameter
#' @param upper the upper bound of the chosen values of the tuning parameter
#' @param l the number of chosen values of the tuning parameter
#'
#' @return A list S2M which will hold p: the total number of variables, b.i: the values of the tuning parameter, H.b.i : the estimated number of signals corresponding to each b.i, abs.post.median: medians of the absolute values of the posterior samples.
#' @export
#'
#' @examples
#'
#'
# Example One Gaussian Model and Horseshoe prior
#'
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' beta <- exp(rnorm(p))
#' Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
#' df <- data.frame(X,Y)
#' rv.hs <- bayesreg::bayesreg(Y~. ,df, "gaussian", "horseshoe+", 200, 100)
#'
#' Beta = t(rv.hs$beta)
#' lower = 0
#' upper = 1
#' l = 20
#' S2Mbeta = Sequential2MeansBeta( Beta, lower, upper, l )
#' H.b.i = S2Mbeta$H.b.i
#'
# Example Two normal model and lasso prior
#'
#'#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' beta <- exp(rnorm(p))
#' Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
#' df <- data.frame(X,Y)
#' rv.hs <- bayesreg::bayesreg(Y~. ,df, "normal", "lasso", 200, 100)
#'
#' Beta = t(rv.hs$beta)
#' lower = 0
#' upper = 1
#' l = 15
#' S2Mbeta = Sequential2MeansBeta( Beta, lower, upper, l )
#' H.b.i = S2Mbeta$H.b.i
#'

Sequential2MeansBeta <- function(Beta, lower, upper, l) {

  # Check for NULL or NaN values in Beta
  if(is.null(Beta) || any(is.na(Beta))){
    stop(paste("Beta must not be NULL or have NaN values. \n Please check Sequential2Means funnction if Beta is not available."))
  }

  # Check for matrix data type of Beta
  if(! is.matrix(Beta)){
    stop(paste("Beta must be passed as matrix data type. "))
  }

  # Lower bound and Upper bound magnitude comparison check
  if(lower > upper){
    stop(paste("Lower bound is greater than Upper bound. "))
  }

  # Check if l is NULL
  if(is.null(l)){
    stop(paste("l must not be NULL. "))
  }

  # Check if l is natural number greater than zero
  if(l<= 0 | l%%1 !=0){
    stop(paste("l must be an integer greater than or equal to zero. "))
  }


  # Initializing variables
  ##########################

  # number data points
  N <- nrow(Beta)

  # number of covariates
  p <- ncol(Beta)

  # tuning parameters
  b.i <- seq(from=lower, to=upper, length.out = l)

  # initializing the number of signals as zero for corresponding b.i
  H.b.i <- rep(0, length = length(b.i))

  # Sequential 2 means algorithm
  for (r in 1:l) {
    impVars.count <- seq(0, length = N)
    for (i in 1:N) {
      impVars.count[i] <- p - numNoiseCoeff(Beta[i, ], b.i[r])
    }
    H.b.i[r] <- as.numeric(names(sort(-table(impVars.count)))[1])
  }

  # list to hold desired output
  S2Mbeta <- list(p = p, b.i = b.i, H.b.i = H.b.i)
  return(S2Mbeta)
}
