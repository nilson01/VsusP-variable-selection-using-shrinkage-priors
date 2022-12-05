#' Variable selection using shrinkage priors:: OptimalHbi
#'
#' OptimalHbi function will take b.i and H.b.i as input which comes from the result of TwoMeans function. It will return H: the optimal value of the tuning parameter.
#' @param bi The values of the tuning parameter
#' @param Hbi The estimated number of signals corresponding to each b.i
#'
#' @return the optimal value of tuning parameter and the associated H value
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
#' rv.hs <- bayesreg::bayesreg(Y~. ,df , "gaussian", "horseshoe+", 3000, 2000)
#'
#' Beta = rv.hs$beta
#' lower = 0
#' upper = 1
#' l = 0.05
#' S2Mbeta = Sequential2MeansBeta( Beta, lower, upper, l)
#'
#' bi = S2Mbeta$b.i
#' Hbi = S2Mbeta$H.b.i
#' OptimalHbi(bi, Hbi)
#'
OptimalHbi <- function(bi, Hbi) {
  plot(bi,Hbi)
}


