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

  # Initializing variables
  N=dim(Beta)[1]
  p=dim(Beta)[2]
  b.i=seq(lower,upper,length=l)
  H.b.i=NULL

  # Sequential 2 means algorithm
  for(r in 1:l){
    KK=NULL
    for (i in 1:N){
      # perform k means for sequential 2 means to prune the noise coefficients
      fit=kmeans(abs(Beta[i,]),2)
      cen1=min(fit$centers)
      cen2=max(fit$centers)

      # perform sequential 2 meanns on subsequent clusters until threshold is reached
    }
    H.b.i[r]=as.numeric(names(sort(-table(KK)))[1])
  }

  S2M = list(p, b.i, H.b.i)
  return(S2M)
}
