#'  Variable selection using shrinkage priors:: Sequential2Means
#'
#'Sequential2means function will take S2M: a list obtained from the 2Means.variables function and H: the estimated number of signals obtained from the optimal.b.i function. This will give out the important subset of variables for the Gaussian Linear model.
#' @param S2M List obtained from the 2Means.variables function
#' @param H Estimated number of signals obtained from the optimal.b.i function
#'
#' @return Important subset of variables for the Gaussian Linear model
#' @export
#'
#' @examples
#' Sequential2Means(S2M, H)
Sequential2Means <- function(S2M, H){
  impVariablesGLM = c()
  return(impVariablesGLM)
}
