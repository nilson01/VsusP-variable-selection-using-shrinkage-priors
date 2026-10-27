#' Summarize sensitivity to sequential two-means thresholds
#'
#' Compare signal counts and selected sets over an explicit threshold grid.
#' This function does not choose a threshold.
#' @inheritParams Sequential2MeansBeta
#' @param b.i Nonempty numeric vector of finite, nonnegative thresholds, in the
#'   order to compare. Duplicate thresholds are allowed.
#' @return A vsusp_tuning list with b.i, H.b.i and abs.post.median as in
#'   Sequential2MeansBeta, plus frequencies (a list of per-draw count tables),
#'   selected (a list of selected-index vectors), mode.proportion (the fraction
#'   of draws giving the modal count), and jaccard.previous (selected-set overlap
#'   with the preceding grid entry; NA for the first entry). The shared p field
#'   permits passing this summary to S2MVarSelectionV1.
#' @details The Jaccard similarity is the intersection size divided by the union
#'   size. It describes overlap, not statistical accuracy. Count frequencies
#'   are clustering frequencies, not posterior inclusion probabilities.
#'   A stable tuning path need not recover the true variable set.
#' @seealso [Sequential2MeansBeta()], [OptimalHbi()]
#' @export
#' @examples
#' Beta <- rbind(c(.01,.02,.03,.3,.31,10,11),
#'               c(-.01,-.02,-.03,-.3,-.31,-10,-11))
#' result <- S2MSensitivity(Beta,c(.1,.5),clustering="exact")
#' result$H.b.i
#' result$jaccard.previous
S2MSensitivity <- function(Beta, b.i, clustering = c("kmeans", "exact"),
                           small.set = c("legacy", "refine")) {
  clustering <- matchChoice(clustering, c("kmeans", "exact"), "clustering")
  small.set <- matchChoice(small.set, c("legacy", "refine"), "small.set")
  checkNumericMatrix(Beta, "Beta", min.cols = 2L)
  checkNumericVector(b.i, "b.i", nonnegative = TRUE)
  result <- summarizeS2M(Beta, b.i, clustering, small.set, details = TRUE)
  result$selected <- lapply(result$H.b.i, function(H) rankSelection(result$abs.post.median, H))
  result$mode.proportion <- vapply(result$frequencies,
                                   function(f) max(f$proportion), numeric(1))
  overlap <- rep(NA_real_, length(b.i))
  if (length(b.i) > 1L) {
    for (i in seq.int(2L, length(b.i))) {
      a <- result$selected[[i - 1L]]
      z <- result$selected[[i]]
      u <- length(union(a, z))
      overlap[i] <- if (u == 0L) 1 else length(intersect(a, z)) / u
    }
  }
  result$jaccard.previous <- overlap
  result
}
