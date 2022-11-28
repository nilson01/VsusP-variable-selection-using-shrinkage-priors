
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
#'
#' @return A list S2M which will hold Beta: N by p matrix consisting of N posterior samples of p variables, b.i: the values of the tuning parameter, H.b.i : the estimated number of signals corresponding to each b.i, abs.post.median: medians of the absolute values of the posterior samples.
#' @export
#'
#' @examples
Sequential2Means <- function(X, Y, b.i) {
  # Initializing variables
  ##########################

  # number data points
  N <- nrow(X)

  # number of covariates
  p <- ncol(X)

  # number of tuning parameters
  l <- length(b.i)

  # initializing the number of signals corresponding to each b.i
  H.b.i <- seq(0, length = length(b.i))

  # N by p matrix consisting of N posterior samples of p variables
  Beta <- matrix(0, nrow = nrow(X), ncol = ncol(Y))

  # .......TODO.......

  # Sequential 2 means algorithm
  for (r in 1:l) {
    impVars.count <- seq(0, length = N)
    for (i in 1:N) {
      impVars.count[i] <- p - numNoiseCoeff(Beta[i, ], b.i[r])
    }
    H.b.i[r] <- as.numeric(names(sort(-table(impVars.count)))[1])
  }

  # list to hold desired output
  S2M <- list(Beta, H.b.i)

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
Sequential2MeansBeta <- function(Beta, lower, upper, l) {
  # Initializing variables
  ##########################

  # number data points
  N <- nrow(Beta)

  # number of covariates
  p <- ncol(Beta)

  # tuning parameters
  b.i <- seq(lower, upper, length = l)

  # initializing the number of signals as zero for corresponding b.i
  H.b.i <- seq(0, length = l)

  # Sequential 2 means algorithm
  for (r in 1:l) {
    impVars.count <- seq(0, length = N)
    for (i in 1:N) {
      impVars.count[i] <- p - numNoiseCoeff(Beta[i, ], b.i[r])
    }
    H.b.i[r] <- as.numeric(names(sort(-table(impVars.count)))[1])
  }

  # list to hold desired output
  S2M <- list(p, b.i, H.b.i)
  return(S2M)
}
