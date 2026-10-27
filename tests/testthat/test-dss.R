# Independent small-problem oracle: enumerate faces/signs and solve each
# constrained quadratic's stationary equations, without coordinate descent.
dss_face_oracle <- function(X, beta, lambda) {
  target <- as.vector(X %*% beta)
  p <- ncol(X); n <- nrow(X)
  patterns <- as.matrix(expand.grid(rep(list(c(-1, 0, 1)), p)))
  best <- Inf; answer <- NULL
  for (i in seq_len(nrow(patterns))) {
    signs <- patterns[i, ]; a <- which(signs != 0)
    if (any(beta[a] == 0)) next
    gamma <- numeric(p)
    if (length(a)) {
      XA <- X[, a, drop = FALSE]
      rhs <- crossprod(XA, target) - n * lambda / 2 * signs[a] / abs(beta[a])
      candidate <- tryCatch(as.vector(solve(crossprod(XA), rhs)),
                            error = function(e) NULL)
      if (is.null(candidate) || any(candidate * signs[a] <= 0)) next
      gamma[a] <- candidate
    }
    penalty <- sum(abs(gamma[beta != 0]) / abs(beta[beta != 0]))
    objective <- mean((target - as.vector(X %*% gamma))^2) + lambda * penalty
    if (objective < best) { best <- objective; answer <- gamma }
  }
  list(coefficients = answer, objective = best)
}

test_that("DSS matches the orthogonal solution with the unhalved MSE convention", {
  X <- diag(c(1, 2, 4)); beta <- c(2, -1, .25); lambda <- .2
  expected <- sign(beta) * pmax(abs(beta) - nrow(X) * lambda /
                                (2 * colSums(X^2) * abs(beta)), 0)
  result <- DSS(X, beta, lambda)
  expect_equal(result$coefficients, expected, tolerance = 1e-12)
  expect_true(result$converged)
  expect_lte(result$kkt.error, 1e-8)
  expect_equal(result$objective, result$mse + result$penalty)
  expect_equal(result$fitted, as.vector(X %*% result$coefficients))
  expect_equal(result$penalty, lambda * sum(abs(result$coefficients / beta)))
})

test_that("DSS includes the unpenalized and all-zero boundaries", {
  X <- diag(2)
  expect_identical(DSS(X, c(2, 0), 0)$coefficients, c(2, 0))
  expect_equal(DSS(X, c(2, 0), 0)$mse, 0)
  for (lambda in c(4, 4 + 1e-8, 10)) {
    fit <- DSS(X, c(2, .1), lambda)
    expect_identical(fit$selected, integer(0))
    expect_true(fit$converged)
  }
  expect_length(DSS(X, c(2, .1), 4 - 1e-8)$selected, 1)
  for (lambda in c(0, .1)) {
    fit <- DSS(X, c(0, 0), lambda)
    expect_identical(fit$selected, integer(0))
    expect_identical(fit$objective, 0)
    expect_true(fit$converged)
  }
  expect_equal(DSS(matrix(2, 1, 1), 3, 2)$coefficients, 3 - 1/12)
})

test_that("correlated DSS solutions agree with independently enumerated faces", {
  for (case in seq_len(8)) {
    X <- cbind(c(1, 0, -1, 2, 1), c(.9, .2, -.8, 1.3, 1),
               sin(seq_len(5) + case / 3))
    beta <- c(1, -.7, (-1)^case * .3)
    lambda <- .01 + case / 30
    oracle <- dss_face_oracle(X, beta, lambda)
    fit <- DSS(X, beta, lambda, tol = 1e-10)
    expect_equal(fit$objective, oracle$objective, tolerance = 1e-9)
    expect_equal(fit$coefficients, oracle$coefficients, tolerance = 1e-7)
    expect_true(fit$converged)
    gradient <- 2 * as.vector(crossprod(X, X %*% (beta - fit$coefficients))) / nrow(X)
    a <- fit$selected; inactive <- setdiff(seq_len(ncol(X)), a)
    if (length(a)) expect_equal(gradient[a], lambda * sign(fit$coefficients[a]) /
                                abs(beta[a]), tolerance = 1e-7)
    if (length(inactive)) expect_true(all(abs(gradient[inactive] * beta[inactive]) <= lambda + 1e-8))
  }
})

test_that("DSS reports limits and supports dependent columns and p greater than n", {
  X <- cbind(c(1, 2, 3), c(1, 2, 3))
  fit <- DSS(X, c(1, 1), .5)
  expect_true(fit$converged)
  expect_equal(sum(fit$coefficients), 2 - .5 / (2 * mean(X[, 1]^2)), tolerance = 1e-10)
  X <- rbind(c(1, 2, -1, .5), c(.1, 1, 2, -.5))
  fit <- DSS(X, c(1, -.2, .5, .3), .2)
  expect_true(fit$converged)
  expect_lte(fit$kkt.error, 1e-8)
  X <- cbind(c(1, 2, 3, 4), c(1.1, 1.8, 3.1, 4.2))
  expect_warning(short <- DSS(X, c(1, 1), .01, max.iter = 1), "KKT tolerance")
  expect_false(short$converged)
  expect_gt(short$kkt.error, 1e-8)
})

test_that("DSS respects zero posterior means, names and predictor units", {
  X <- cbind(a = c(1, 0, 2, -1), b = c(.9, .1, 1, -.8), c = c(1, 2, .5, 1))
  beta <- c(a = 1, b = -.5, c = 0)
  fit <- DSS(X, beta, .1, tol = 1e-10)
  expect_identical(fit$coefficients[3], c(c = 0))
  expect_named(fit$coefficients, colnames(X))
  factors <- c(10, -.2, 3)
  changed <- DSS(sweep(X, 2, factors, "*"), beta / factors, .1, tol = 1e-10)
  expect_equal(changed$coefficients * factors, fit$coefficients, tolerance = 1e-9)
  expect_equal(changed$fitted, fit$fitted, tolerance = 1e-9)
  expect_equal(changed$objective, fit$objective, tolerance = 1e-9)
  response <- DSS(X, beta * 3, .1 * 9, tol = 1e-10)
  expect_equal(response$coefficients, fit$coefficients * 3, tolerance = 1e-8)
  expect_equal(response$objective, fit$objective * 9, tolerance = 1e-8)
  set.seed(11); saved <- .Random.seed
  DSS(X, beta, .1)
  expect_identical(.Random.seed, saved)
})

test_that("DSS rejects malformed inputs and nonrepresentable arithmetic", {
  expect_error(DSS(data.frame(a = 1:3), 1, .1), "numeric matrix")
  expect_error(DSS(diag(2), 1, .1), "one coefficient")
  expect_error(DSS(diag(2), c(1, NA), .1), "finite numeric")
  expect_error(DSS(matrix(c(0, 0, 1, 2), 2), c(1, 1), .1), "zero-norm")
  for (lambda in list(-1, NA, Inf, c(1, 2), "a")) {
    expect_error(DSS(diag(2), c(1, 2), lambda))
  }
  for (tol in list(0, -1, 1, NA, Inf, c(.1, .2))) {
    expect_error(DSS(diag(2), c(1, 2), .1, tol = tol))
  }
  for (m in list(0, -1, 1.5, NA, Inf, .Machine$integer.max + 1)) {
    expect_error(DSS(diag(2), c(1, 2), .1, max.iter = m))
  }
  X <- diag(2); colnames(X) <- c("a", "b")
  expect_error(DSS(X, c(b = 1, a = 2), .1), "same order")
  colnames(X) <- c("a", "a")
  expect_error(DSS(X, c(1, 2), .1), "unique and nonempty")
  expect_error(DSS(diag(2), c(a = 1, a = 2), .1), "unique and nonempty")
  expect_error(DSS(matrix(1e200, 1, 1), 1, .1), "representable range")
  expect_error(DSS(matrix(1e-200, 1, 1), 1, .1), "underflow")
  expect_error(DSS(matrix(1e200, 1, 1), 1e200, .1), "representable range")
  expect_error(DSS(matrix(c(1, -1), 1), c(1e-160, 1e-160 * (1 + 1e-15)), 1), "underflow")
})
