
################################################################
#------------------ S2MVarSelection algorithm -----------------#
################################################################

#' Select coefficients by posterior median absolute magnitude
#'
#' Return the indices of the H largest column medians of absolute draws.
#' @param Beta Finite numeric matrix with posterior draws in rows and
#'   coefficients in columns. Exclude any intercept column.
#' @param H Integer number of variables to select, from zero to ncol(Beta).
#' @return Integer vector of selected column indices in decreasing order
#'   of median absolute magnitude. Equal medians favor larger column indices.
#'   H = 0 returns integer(0).
#' @details H is supplied by the caller. This function does not estimate H
#'   or test statistical significance.
#' @seealso [S2MVarSelectionV1()], [Sequential2MeansBeta()]
#' @export
#' @references Li, H. and Pati, D. (2017). Variable selection using shrinkage
#'   priors. Computational Statistics & Data Analysis, 107, 107--119.
#'   \doi{10.1016/j.csda.2016.10.008}.
#' @examples
#' Beta <- rbind(c(.01, .02, .03, .3, .31, 10, 11),
#'               c(-.01, -.02, -.03, -.3, -.31, -10, -11))
#' set.seed(1)
#' result <- Sequential2MeansBeta(Beta, .1, .5, 2)
#' S2MVarSelection(Beta, result$H.b.i[2])
#' S2MVarSelection(Beta, 0)
S2MVarSelection <- function(Beta, H = 5) {
  checkNumericMatrix(Beta, "Beta")
  checkSelectionCount(H, ncol(Beta))
  rankSelection(apply(abs(Beta), 2L, stats::median), H)
}



################################################################
#-------------------- OptimalHbi algorithm --------------------#
################################################################


#' Plot the sequential two-means tuning path
#'
#' Display the estimated signal count at each threshold to support a choice
#' of tuning parameter. The function does not choose a threshold or H.
#' @param bi Nonempty vector of finite, nonnegative thresholds.
#' @param Hbi Nonnegative integer counts, with the same length as bi.
#' @return Invisible NULL. A plot is drawn on the current graphics device.
#' @seealso [Sequential2MeansBeta()]
#' @export
#' @references Li, H. and Pati, D. (2017). Variable selection using shrinkage
#'   priors. Computational Statistics & Data Analysis, 107, 107--119.
#'   \doi{10.1016/j.csda.2016.10.008}.
#' @examples
#' Beta <- rbind(c(.01, .02, .03, .3, .31, 10, 11),
#'               c(-.01, -.02, -.03, -.3, -.31, -10, -11))
#' set.seed(1)
#' result <- Sequential2MeansBeta(Beta, .1, .5, 2)
#' OptimalHbi(result$b.i, result$H.b.i)
OptimalHbi <- function(bi, Hbi) {
  checkNumericVector(bi, "bi", nonnegative = TRUE)
  checkNumericVector(Hbi, "Hbi", nonnegative = TRUE)
  if (any(Hbi %% 1 != 0)) stop("Hbi must contain integer counts.")
  if (length(bi) != length(Hbi)) stop("bi and Hbi must have the same length.")
  graphics::plot(bi, Hbi, type = "b", xlab = "Tuning parameter b",
       ylab = "Estimated number of signals", main = "S2M tuning path")
  invisible(NULL)
}

################################################################
#---------------- S2MVarSelectionV1 algorithm -----------------#
################################################################

#' Select coefficients from an S2M summary
#'
#' Select the H largest posterior median absolute coefficients from a
#' summary returned by either S2M fitting function.
#' @param S2M List containing p (a positive integer) and abs.post.median
#'   (a finite, nonnegative numeric vector of length p). Legacy lists with
#'   these two fields are also accepted.
#' @param H Integer number of variables to select, from zero to S2M$p.
#' @return Integer vector of selected column indices, ordered by decreasing
#'   median absolute magnitude. Equal medians favor larger column indices.
#'   H = 0 returns integer(0).
#' @seealso [S2MVarSelection()], [Sequential2MeansBeta()]
#' @export
#' @references Li, H. and Pati, D. (2017). Variable selection using shrinkage
#'   priors. Computational Statistics & Data Analysis, 107, 107--119.
#'   \doi{10.1016/j.csda.2016.10.008}.
#' @examples
#' Beta <- rbind(c(.01, .02, .03, .3, .31, 10, 11),
#'               c(-.01, -.02, -.03, -.3, -.31, -10, -11))
#' set.seed(1)
#' result <- Sequential2MeansBeta(Beta, .1, .5, 2)
#' S2MVarSelectionV1(result, result$H.b.i[2])
S2MVarSelectionV1 <- function(S2M, H = 5) {
  if (!is.list(S2M)) stop("S2M must be a list with p and abs.post.median.")
  checkIntegerScalar(S2M$p, "S2M$p", lower = 1)
  checkNumericVector(S2M$abs.post.median, "S2M$abs.post.median", nonnegative = TRUE)
  if (length(S2M$abs.post.median) != S2M$p) {
    stop("S2M$abs.post.median must have length S2M$p.")
  }
  checkSelectionCount(H, S2M$p)
  rankSelection(S2M$abs.post.median, H)
}
