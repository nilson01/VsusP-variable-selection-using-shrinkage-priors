#' Variable selection using shrinkage priors:: Sequential2Means
#'
#' Sequential2Means function will take as input X: design matrix, Y : response vector, t: vector of tuning parameter values from Sequential 2-means (S2M) variable selection algorithm. The function will return a list S2M which will hold p: the total number of variables, b.i: the values of the tuning parameter, H.b.i : the estimated number of signals corresponding to each b.i, abs.post.median: medians of the absolute values of the posterior samples.
#' @param X Design matrix
#' @param Y Response vector
#' @param t Vector of tuning parameter values from Sequential 2-means (S2M) variable selection algorithm.
#'
#' @return A list S2M which will hold p: the total number of variables, b.i: the values of the tuning parameter, H.b.i : the estimated number of signals corresponding to each b.i, abs.post.median: medians of the absolute values of the posterior samples.
#' @export
#'
#' @examples
#' Sequential2Means(X, Y, t)
Sequential2Means <- function(X, Y, t){
  p = 0
  b.i = c()
  H.b.i = c()
  S2M = list(p, b.i, H.b.i)
  return(S2M)
}
