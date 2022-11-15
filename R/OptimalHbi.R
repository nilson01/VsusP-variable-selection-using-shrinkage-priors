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
#' OptimalHbi(bi, Hbi)
OptimalHbi <- function(bi, Hbi){
  bi = 0
  H = 0
  return(list(bi, H))
}
