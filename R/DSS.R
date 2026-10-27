#' Sparse prediction with the DSS adaptive-lasso approximation
#'
#' Approximate the posterior-mean linear predictor with a sparse coefficient
#' vector using the decoupled shrinkage and selection (DSS) adaptive-L1
#' criterion of Hahn and Carvalho (2015), equation (21).
#' @inheritParams SAVS
#' @param lambda Nonnegative scalar penalty. It must be chosen by the caller;
#'   larger values penalize nonzero coefficients more strongly.
#' @param tol Relative tolerance for the convex optimality (KKT) conditions,
#'   strictly between zero and one. Default 1e-8.
#' @param max.iter Maximum number of coordinate-descent sweeps. Default 10000.
#' @details Here beta represents the posterior mean, and X contains the
#'   prediction points at which the sparse predictor should approximate X beta.
#'   These need not be the observations used to fit the posterior.
#'   The objective is
#'   \deqn{n^{-1}\|X\bar{\beta}-X\gamma\|_2^2 +
#'     \lambda\sum_{j:\bar{\beta}_j\ne0}|\gamma_j|/|\bar{\beta}_j|.}
#'   This is the paper's
#'   adaptive-L1 approximation, not exact minimization of its L0 counting penalty.
#'
#'   For positive lambda, a zero beta entry constrains the corresponding sparse
#'   coefficient to zero, the infinite-weight limit. For lambda = 0, beta itself
#'   is returned as a least-squares minimizer. No intercept, centering or scaling
#'   is added. Exclude the intercept and keep its contribution separately.
#'   The adaptive criterion is unchanged by rescaling predictor columns and
#'   inversely rescaling beta, apart from floating-point error and nonunique
#'   solutions. Changing response units requires changing lambda accordingly.
#'
#'   Cyclic coordinate descent solves the convex surrogate in the transformed
#'   design X times abs(beta). The KKT error is divided by the larger of lambda
#'   and the largest absolute initial gradient. A result is converged when this
#'   ratio is at most tol. The iteration limit produces a warning and
#'   converged = FALSE; inspect this flag before using the selection.
#'   Overflow and detected underflow in transformed columns, norms, sparse
#'   coefficients or squared errors are rejected. With dependent
#'   predictor columns, coefficients and selected sets need not be unique;
#'   coordinates are visited in the supplied column order.
#'
#'   Squared error here measures approximation to the supplied posterior mean
#'   at X, not held-out predictive error. The function does not select lambda,
#'   compute inclusion probabilities or establish false-discovery control.
#' @return A vsusp_sparse list with coefficients, selected (integer column indices, possibly
#'   integer(0)), p, lambda, fitted (X coefficients), mse (mean squared
#'   approximation error), penalty (weighted L1 penalty including lambda),
#'   objective (mse + penalty), iterations, kkt.error (relative KKT error) and
#'   converged. Numerical diagnostics concern the convex surrogate only.
#' @references Hahn, P. R. and Carvalho, C. M. (2015). Decoupling Shrinkage and
#'   Selection in Bayesian Linear Models: A Posterior Summary Perspective.
#'   Journal of the American Statistical Association, 110(509), 435-448.
#'   \doi{10.1080/01621459.2014.993077}.
#' @export
#' @examples
#' X <- diag(2)
#' DSS(X, c(2, .1), lambda = .2)$coefficients
#' DSS(X, c(2, .1), lambda = 10)$selected
DSS <- function(X, beta, lambda, tol = 1e-8, max.iter = 10000L) {
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
  if (any(apply(abs(X), 2L, max) == 0)) stop("X must not contain zero-norm columns.")
  checkNumericScalar(lambda, "lambda", nonnegative = TRUE)
  checkNumericScalar(tol, "tol")
  if (tol <= 0 || tol >= 1) stop("tol must be strictly between zero and one.")
  checkIntegerScalar(max.iter, "max.iter", lower = 1)
  if (max.iter > .Machine$integer.max) stop("max.iter exceeds the supported integer range.")
  finite <- function(x) {
    if (any(!is.finite(x))) {
      stop("DSS arithmetic is outside the representable range; rescale the inputs.",
           call. = FALSE)
    }
    x
  }
  squared <- function(x) {
    result <- finite(x^2)
    if (any(x != 0 & result == 0)) {
      stop("DSS squared errors or norms underflow; rescale the inputs.", call. = FALSE)
    }
    result
  }
  n <- nrow(X)
  target <- finite(as.vector(X %*% beta))
  squared(target) # Detect an unrepresentable null-predictor squared error.
  active <- which(beta != 0)
  Z <- finite(sweep(X[, active, drop = FALSE], 2L, abs(beta[active]), "*"))
  if (length(active) && any(X[, active, drop = FALSE] != 0 & Z == 0)) {
    stop("DSS transformed columns underflow; rescale the inputs.")
  }
  norm2 <- finite(colMeans(squared(Z)))
  if (any(norm2 == 0)) {
    stop("DSS transformed column norms underflow; rescale the inputs.")
  }
  gradient0 <- finite(2 * as.vector(crossprod(Z, target)) / n)
  scale <- max(lambda, abs(gradient0), 0)
  alpha <- numeric(length(active))
  iterations <- 0L
  if (lambda == 0) {
    alpha <- sign(beta[active])
  } else if (length(active) && lambda < max(abs(gradient0))) {
    threshold <- lambda / 2
    if (threshold == 0) stop("lambda/2 underflows; rescale the inputs.")
    residual <- target
    for (iteration in seq_len(max.iter)) {
      for (j in seq_along(active)) {
        partial <- finite(sum(Z[, j] * residual) / n + norm2[j] * alpha[j])
        next.alpha <- finite(sign(partial) * max(abs(partial) - threshold, 0) / norm2[j])
        residual <- finite(residual - Z[, j] * (next.alpha - alpha[j]))
        alpha[j] <- next.alpha
      }
      # Recompute rather than accumulate residual drift across sweeps.
      residual <- finite(target - as.vector(Z %*% alpha))
      gradient <- finite(2 * as.vector(crossprod(Z, residual)) / n)
      violation <- ifelse(alpha != 0, abs(gradient - lambda * sign(alpha)),
                          pmax(abs(gradient) - lambda, 0))
      iterations <- iteration
      if (max(violation) / scale <= tol) break
    }
  }
  coefficients <- beta
  coefficients[] <- 0
  coefficients[active] <- finite(abs(beta[active]) * alpha)
  if (any(alpha != 0 & coefficients[active] == 0)) {
    stop("DSS sparse coefficients underflow; rescale the inputs.")
  }
  labels <- if (!is.null(colnames(X))) colnames(X) else names(beta)
  names(coefficients) <- labels
  fitted <- finite(as.vector(X %*% coefficients))
  residual <- finite(target - fitted)
  effective.alpha <- coefficients[active] / abs(beta[active])
  gradient <- finite(2 * as.vector(crossprod(Z, residual)) / n)
  violation <- ifelse(effective.alpha != 0,
                      abs(gradient - lambda * sign(effective.alpha)),
                      pmax(abs(gradient) - lambda, 0))
  kkt.error <- if (scale == 0) 0 else max(c(0, violation)) / scale
  residual.squared <- squared(residual)
  mse <- finite(mean(residual.squared))
  if (mse == 0 && any(residual.squared != 0)) {
    stop("DSS mean squared error underflows; rescale the inputs.")
  }
  penalty <- finite(lambda * sum(abs(effective.alpha)))
  objective <- finite(mse + penalty)
  converged <- kkt.error <= tol
  if (!converged) warning("DSS did not meet the KKT tolerance; increase max.iter or inspect predictor dependence.",
                          call. = FALSE)
  structure(list(coefficients = coefficients, selected = unname(which(coefficients != 0)),
       p = ncol(X), lambda = lambda, fitted = fitted, mse = mse, penalty = penalty,
       objective = objective, iterations = iterations, kkt.error = kkt.error,
       converged = converged), class = "vsusp_sparse")
}
