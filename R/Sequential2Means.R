#' Count coefficients in the saved S2M noise set
#' @param Beta.i Numeric coefficient vector for one posterior draw.
#' @param b.i_r Nonnegative center-gap threshold.
#' @return Size of the saved lower-mean candidate set.
#' @inheritParams Sequential2MeansBeta
#' @keywords internal
numNoiseCoeff <- function(Beta.i, b.i_r, clustering = "kmeans", small.set = "legacy") {
  # A lower-mean cluster requires at least two distinct magnitudes.
  if (length(unique(abs(Beta.i))) < 2L) {
    stop("Each Beta draw must contain at least two distinct absolute coefficients.")
  }
  # With two distinct observations the two singleton clusters are determined.
  if (length(Beta.i) == 2L) return(1L)
  fit <- splitS2M(Beta.i, clustering)
  # Initialize |D| to the lower-mean cluster (Li and Pati, 2017, section 2.2.1).
  numNoise <- sum(fit$cluster == which.min(fit$centers))

  # perform sequential 2 means on subsequent clusters until threshold is reached
  # between two centers
  while (max(fit$centers) - min(fit$centers) > b.i_r) {
    # finding the noise coefficients from smaller cluster derived from kmeans fit
    Beta.i <- Beta.i[which(fit$cluster == which.min(fit$centers))]
    numNoise <- length(Beta.i)
    # breaking loop for minimum size of beta vector
    if (length(Beta.i) <= 1L ||
        (small.set == "legacy" && length(Beta.i) <= 2L) ||
        length(unique(abs(Beta.i))) < 2L) {
      break
    }
    # clustering the noise coefficients and signals
    fit <- splitS2M(Beta.i, clustering)
  }
  # Return the saved |D|, not the smaller cluster from the stopping split.
  return(numNoise)
}
#--------------------  Sequential2Means  --------------------#
#' Fit Gaussian shrinkage regression and apply sequential two-means
#'
#' Fit a Gaussian linear model with bayesreg, then apply sequential two-means
#' (S2M) to its coefficient draws. The intercept is excluded from selection.
#' @param X Finite numeric design matrix with observations in rows and at
#'   least two predictors in columns.
#' @param Y Finite numeric response vector of length nrow(X).
#' @param b.i Nonempty vector of finite, nonnegative center-gap thresholds.
#' @param prior Shrinkage prior: "ridge" ("rr"), "lasso", "horseshoe"
#'   ("hs"), or "horseshoe+" ("hs+").
#' @param n.samples Integer number of requested posterior draws, at least 100.
#'   bayesreg may round this count up when using parallel workers.
#' @param burnin Integer number of burn-in iterations, at least 100.
#' @param clustering Two-means backend. "kmeans" retains the existing randomized
#'   implementation; "exact" scans all distinct adjacent cuts in one dimension.
#' @param small.set Recursive stopping convention. "legacy" retains a candidate
#'   with at most two elements. "refine" splits a distinct pair and reapplies
#'   the center-gap rule, stopping at a singleton or constant candidate.
#' @return A list with Beta (returned draws in rows, predictors in columns),
#'   p (number of predictors), b.i (thresholds), H.b.i (modal signal counts),
#'   and abs.post.median (column medians of absolute coefficient draws).
#' @details All returned draws are used. Thresholds are on the coefficient
#'   scale, so predictor scaling affects their interpretation. See
#'   [Sequential2MeansBeta()] for stopping conventions and limitations.
#'   A short chain demonstrates the call; assess convergence and effective
#'   sample size before using posterior summaries for inference.
#' @seealso [Sequential2MeansBeta()], [S2MVarSelectionV1()]
#' @export
#' @references Li, H. and Pati, D. (2017). Variable selection using shrinkage
#'   priors. Computational Statistics & Data Analysis, 107, 107--119.
#'   \doi{10.1016/j.csda.2016.10.008}.
#'
#'   Wang, H. and Song, M. (2011). Ckmeans.1d.dp: Optimal k-means clustering
#'   in one dimension by dynamic programming. The R Journal, 3(2), 29--33.
#'   \doi{10.32614/RJ-2011-015}.
#' @examples
#' set.seed(12)
#' X <- matrix(rnorm(120), 40, 3)
#' Y <- 3 * X[, 1] + rnorm(40)
#' result <- Sequential2Means(X, Y, c(.1, .5), prior = "ridge",
#'                            n.samples = 100, burnin = 100)
#' dim(result$Beta)
#' S2MVarSelectionV1(result, result$H.b.i[1])
Sequential2Means <- function(X, Y, b.i, prior = "horseshoe+", n.samples = 5000, burnin = 2000,
                             clustering = c("kmeans", "exact"),
                             small.set = c("legacy", "refine")) {
  clustering <- matchChoice(clustering, c("kmeans", "exact"), "clustering")
  small.set <- matchChoice(small.set, c("legacy", "refine"), "small.set")
  checkNumericMatrix(X, "X", min.cols = 2L)
  checkNumericVector(Y, "Y")
  if (length(Y) != nrow(X)) stop("X and Y must have the same number of observations.")
  checkNumericVector(b.i, "b.i", nonnegative = TRUE)
  prior <- normalizePrior(prior)
  checkIntegerScalar(n.samples, "n.samples", lower = 100)
  checkIntegerScalar(burnin, "burnin", lower = 100)

  # Internal names prevent a predictor named Y from becoming the response.
  data <- as.data.frame(X)
  names(data) <- paste0("X", seq_len(ncol(X)))
  data[[".VsusP.response"]] <- Y
  fit <- fitShrinkageGaussian(.VsusP.response ~ ., data, model = "gaussian",
                           prior = prior, n.samples = n.samples, burnin = burnin)
  Beta <- t(fit$beta)
  if (!is.null(colnames(X))) colnames(Beta) <- colnames(X)
  checkNumericMatrix(Beta, "Beta", min.cols = 2L)
  result <- summarizeS2M(Beta, b.i, clustering, small.set)
  # Use every returned draw: parallel samplers may round the requested count up.
  c(list(Beta = Beta), result)
}



################################################################
#--------------- Sequential2MeansBeta algorithm ---------------#
################################################################

#' Apply sequential two-means to posterior coefficient draws
#'
#' Estimate a signal count at each threshold from supplied coefficient draws,
#' following section 2.2.1 of Li and Pati (2017).
#' @param Beta Finite numeric matrix with posterior draws in rows and
#'   at least two coefficients in columns. Exclude any intercept column.
#'   Each draw must have at least two distinct absolute coefficients.
#' @param lower,upper Finite, nonnegative threshold bounds, with lower <= upper.
#' @param zero.draws Handling of entirely zero draws: `"error"` (default)
#'   rejects an undefined split; `"empty"` assigns zero signals. The latter
#'   is an explicit boundary extension for exact spike-and-slab draws.
#' @param l Positive integer number of equally spaced thresholds.
#' @param clustering Two-means backend. "kmeans" retains the existing randomized
#'   implementation; "exact" scans all distinct adjacent cuts in one dimension.
#' @param small.set Recursive stopping convention. "legacy" retains a candidate
#'   with at most two elements. "refine" splits a distinct pair and reapplies
#'   the center-gap rule, stopping at a singleton or constant candidate.
#' @return A list with p (number of coefficients), b.i (thresholds),
#'   H.b.i (modal signal counts), and abs.post.median (column medians of
#'   absolute coefficient draws).
#' @details For each draw, S2M first splits absolute coefficients into two
#'   clusters. It saves the lower-mean cluster as the noise candidate set.
#'   While the center gap exceeds the threshold, it saves the current
#'   lower-mean cluster and splits that set again. At stopping it uses the
#'   saved set, rather than the lower cluster from the stopping split.
#'
#'   By default, the implementation retains a candidate set of at most two elements, or
#'   a candidate with identical magnitudes, without another split. For an
#'   initial pair of distinct magnitudes, the lower singleton is the noise
#'   set. An initial draw with identical magnitudes has no lower-mean cluster
#'   and is rejected. These are computational boundary conventions. The refine
#'   option can change counts when a retained pair has gap greater than b.
#'
#'   H.b.i is the mode of the per-draw signal counts; a frequency tie selects
#'   the smaller count. Accepted draws always leave at least one signal, so
#'   S2M cannot automatically select an empty model. This clustering summary
#'   is not a posterior inclusion probability or an error-control guarantee.
#'   Results can depend on kmeans initialization; use set.seed() for
#'   reproducibility. The exact backend consumes no random draws. It finds the
#'   minimum within-cluster sum of squares in floating-point arithmetic; equal
#'   losses choose the earliest sorted cut. Numeric rescaling that collapses
#'   distinct values is rejected. An exact clustering optimum does not establish
#'   optimal variable selection. Scaling coefficients changes threshold meaning.
#'   With `zero.draws = "empty"`, all-zero draws contribute zero to the
#'   count frequencies and their row indices are returned as `zero.draws`.
#'   Other constant-magnitude draws remain invalid. This convention does not
#'   make S2M an empty-model selector for continuous-prior draws.
#' @seealso [Sequential2Means()], [OptimalHbi()], [S2MVarSelection()]
#' @export
#' @references Li, H. and Pati, D. (2017). Variable selection using shrinkage
#'   priors. Computational Statistics & Data Analysis, 107, 107--119.
#'   \doi{10.1016/j.csda.2016.10.008}.
#'
#'   Wang, H. and Song, M. (2011). Ckmeans.1d.dp: Optimal k-means clustering
#'   in one dimension by dynamic programming. The R Journal, 3(2), 29--33.
#'   \doi{10.32614/RJ-2011-015}.
#' @examples
#' Beta <- rbind(c(.01, .02, .03, .3, .31, 10, 11),
#'               c(-.01, -.02, -.03, -.3, -.31, -10, -11))
#' set.seed(1)
#' result <- Sequential2MeansBeta(Beta, .1, .5, 2)
#' result$H.b.i
#' S2MVarSelectionV1(result, result$H.b.i[2])
Sequential2MeansBeta <- function(Beta, lower, upper, l,
                                 clustering = c("kmeans", "exact"),
                                 small.set = c("legacy", "refine"),
                                 zero.draws = c("error", "empty")) {
  clustering <- matchChoice(clustering, c("kmeans", "exact"), "clustering")
  small.set <- matchChoice(small.set, c("legacy", "refine"), "small.set")
  zero.draws <- matchChoice(zero.draws, c("error", "empty"), "zero.draws")
  checkNumericMatrix(Beta, "Beta", min.cols = 2L)
  checkNumericScalar(lower, "lower", nonnegative = TRUE)
  checkNumericScalar(upper, "upper", nonnegative = TRUE)
  if (lower > upper) stop("lower must not exceed upper.")
  checkIntegerScalar(l, "l", lower = 1)
  summarizeS2M(Beta, seq(from = lower, to = upper, length.out = l),
               clustering, small.set, zero.draws = zero.draws)
}
