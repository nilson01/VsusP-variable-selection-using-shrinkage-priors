#' Select variables with the signal adaptive variable selector
#'
#' Apply the closed-form signal adaptive variable selector (SAVS) of Ray and
#' Bhattacharya to a coefficient point estimate, such as a posterior mean.
#' The selected set can be empty.
#' @param X Finite numeric design matrix with observations in rows and
#'   coefficients in columns. Exclude the intercept and zero-norm columns.
#' @param beta Finite numeric coefficient vector in the same column order as X.
#'   Supplied names must be unique and nonempty. If both beta and X have
#'   names, their order must agree.
#' @return A vsusp_sparse list with coefficients (sparse estimates), selected (integer column
#'   indices), p (number of coefficients), and log.score (log decision scores).
#'   selected is integer(0) when every coefficient is set to zero.
#' @details A coefficient is retained when its absolute value cubed, multiplied
#'   by the squared norm of its predictor column, exceeds one. Equality gives
#'   zero. An exactly zero input coefficient gives zero without division.
#'   Direct products are used when their factors are representable. Otherwise
#'   logarithms avoid overflow; scores indistinguishable from the threshold
#'   within floating-point rounding error are set to zero.
#'
#'   No centering or scaling is performed. Use the same design matrix and units
#'   that define beta. SAVS depends on units: changing predictor or response
#'   scales can change the selected set. Standardize predictors before fitting
#'   if selection on a common predictor scale is wanted.
#'
#'   This function applies equation (2) and Algorithm 1 of the cited preprint.
#'   It does not solve the full correlated-design penalized optimization
#'   problem. It neither computes inclusion probabilities nor guarantees
#'   false-discovery control. Its behavior depends on the supplied estimate.
#' @references Ray, P. and Bhattacharya, A. (2018). Signal Adaptive Variable
#'   Selector for the Horseshoe Prior. arXiv:1810.09004.
#'   \url{https://arxiv.org/abs/1810.09004}.
#' @section Selection pilot:
#' Empty-selection and exact-support proportions in six Gaussian conditions,
#' with ten datasets per condition and 95 percent Wilson intervals. This small
#' pilot does not establish false-discovery control. See the worked guide for
#' the design and interpretation.
#' \if{html}{\figure{selection-pilot.png}{Empty-selection and exact-support proportions with Wilson intervals.}}
#' \if{latex}{\figure{selection-pilot.pdf}{options: width=5.8in}}
#' @export
#' @examples
#' X <- diag(2)
#' SAVS(X, c(.1, -.1))$selected
#' SAVS(X, c(2, 0))$coefficients
SAVS <- function(X, beta) {
  checkNumericMatrix(X, "X")
  checkNumericVector(beta, "beta")
  if (length(beta) != ncol(X)) stop("beta must have one coefficient per column of X.")
  for (labels in list(colnames(X), names(beta))) {
    if (!is.null(labels) && (anyNA(labels) || any(labels == "") || anyDuplicated(labels))) {
      stop("X column names and beta names must be unique and nonempty when supplied.")
    }
  }
  if (!is.null(names(beta)) && !is.null(colnames(X)) &&
      !identical(names(beta), colnames(X))) {
    stop("beta names must match X column names in the same order.")
  }
  magnitudes <- apply(abs(X), 2L, max)
  if (any(magnitudes == 0)) stop("X must not contain zero-norm columns.")
  # Scale each column before squaring; neither column units nor beta are changed.
  normalized <- sweep(X, 2L, magnitudes, "/")
  log.norm2 <- 2 * log(magnitudes) + log(colSums(normalized^2))
  log.beta <- log(abs(beta))
  log.score <- 3 * log.beta + log.norm2
  # The log sum can round a true equality above zero. Use a direct score when
  # its factors are representable; otherwise retain a conservative roundoff band.
  norm2 <- colSums(X^2)
  cube <- abs(beta)^3
  direct <- which(is.finite(norm2) & norm2 > 0 & is.finite(cube) & cube > 0)
  tolerance <- 8 * .Machine$double.eps * pmax(1, abs(log.norm2) + 3 * abs(log.beta))
  keep <- log.score > tolerance
  if (length(direct)) {
    score <- norm2[direct] * cube[direct]
    keep[direct] <- score > 1
    finite <- which(is.finite(score) & score > 0)
    log.score[direct[finite]] <- log(score[finite])
  }
  coefficients <- beta
  coefficients[] <- 0
  active <- which(keep)
  coefficients[active] <- beta[active] * (-expm1(-log.score[active]))
  labels <- if (!is.null(colnames(X))) colnames(X) else names(beta)
  names(coefficients) <- labels
  names(log.score) <- labels
  structure(list(coefficients = coefficients, selected = unname(which(coefficients != 0)),
       p = ncol(X), log.score = log.score), class = "vsusp_sparse")
}
