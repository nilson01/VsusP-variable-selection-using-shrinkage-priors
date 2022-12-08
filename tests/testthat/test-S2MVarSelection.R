test_that("S2MVarSelection works", {

  # Data Generation

  n <- 100
  p <- 20
  X <- matrix(rnorm(n*p), n, p)
  beta <- exp(rnorm(p))
  Y <- X %*% beta + rnorm(n, 0, 1)
  df <- data.frame(X,Y)
  rv.hs <- bayesreg::bayesreg(Y~. ,df, "gaussian", "horseshoe+", 150, 100)
  Beta = rv.hs$beta
  H = 12

  #--------------------- S2MVarSelection -----------------------#

  # Check for NULL or NaN values of Beta
  testthat::expect_error(S2MVarSelection(NULL, H))

  # Check for matrix data type of Beta
  testthat::expect_error(S2MVarSelection(Beta[,1], H))

  # Check for NULL or NaN values for H
  testthat::expect_error(S2MVarSelection(Beta, NULL))

  # Check for numeric data type of H
  testthat::expect_error(S2MVarSelection(Beta, "H"))

  # Expecting output to be of predefined length for impVariablesGLM
  testthat::expect_length((S2MVarSelection(Beta, H)), H)

  #--------------------- S2MVarSelectionV1 ---------------------#

  # the medians of the absolute values of the posterior samples of each Beta vector
  abs.post.median <- seq(0, length = p)
  for (i in 1:p) {
    abs.post.median[i] <- stats::median(abs(Beta[, i]))
  }

  S2M = list(p = p, abs.post.median=abs.post.median)

  # Check for NULL or NaN values of abs.post.median
  testthat::expect_error(S2MVarSelectionV1(NULL, H))

  # Check for vector data type of abs.post.median
  testthat::expect_error(S2MVarSelection(list(abs.post.median="ABC"), H))

  # Check for NULL or NaN values for H
  testthat::expect_error(S2MVarSelectionV1(S2M, NULL))

  # Check for numeric data type of H
  testthat::expect_error(S2MVarSelectionV1(S2M, "H"))

  # Check for non negative values in abs.post.median
  testthat::expect_error(S2MVarSelectionV1(S2M = list(abs.post.median=c(-1,-2,3)), "H"))

  # Check for vector data type of abs.post.median
  # testthat::expect_length(S2MVarSelectionV1(S2M, H), H)

  #--------------------- OptimalHbi ---------------------#

  # Check for NULL or NaN values in b.i
  testthat::expect_error(OptimalHbi(NULL, Hbi = c(17,12,12,12,12,6,6,6,4,5)))

  # Check for vector data type of b.i
  testthat::expect_error(OptimalHbi(matrix(rnorm(4),2,2), Hbi = c(17,12,12,12,12,6,6,6,4,5)))


})


