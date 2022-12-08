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


  # Check for NULL or NaN values in b.i
  testthat::expect_error(OptimalHbi(NULL, Hbi = c(17,12,12,12,12,6,6,6,4,5)))

  # Check for vector data type of b.i
  testthat::expect_error(OptimalHbi(matrix(rnorm(4),2,2), Hbi = c(17,12,12,12,12,6,6,6,4,5)))


})


