
test_that("Sequential2Means works", {

  # Generating data
  n <- 100
  p <- 20
  X <- matrix(rnorm(n*p), n, p)
  beta <- exp(rnorm(p))
  Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
  b.i = seq(0,1, 0.05)

  #----------------------------------------------#
  ######### Sequential2Means  #########
  #----------------------------------------------#

  ########### INPUT TESTS

  # Check for NULL or NaN values in X
  testthat::expect_error(Sequential2Means(NULL, Y, b.i, "ridge", 200, 100))

  # Check for matrix data type of X
  testthat::expect_error(Sequential2Means(Y, Y, b.i, "ridge", 200, 100))

  # Check for NULL or NaN values in Y
  testthat::expect_error(Sequential2Means(X, NULL, b.i, "ridge", 200, 100))

  # Check for vector data type of Y
  testthat::expect_error(Sequential2Means(X, X, b.i, "ridge", 200, 100))

  # Check for compatibility of dimensions between X and Y
  testthat::expect_error(Sequential2Means(X, Y[-1], b.i, "ridge", 200, 100))

  # Check for NULL or NaN values in b.i
  testthat::expect_error(Sequential2Means(X, Y, NULL, "ridge", 200, 100))

  # Check for vector data type of b.i
  testthat::expect_error(Sequential2Means(X, Y, X, "ridge", 200, 100))

  # Check for priors from available options
  testthat::expect_error(Sequential2Means(X, Y, b.i, "abc", 200, 100))

  # Check for the lower bound of n.samples
  testthat::expect_error(Sequential2Means(X, Y, b.i, "ridge", 10, 200))

  # Check if n.samples is a natural number
  testthat::expect_error(Sequential2Means(X, Y, b.i, "ridge", 200.5213, 100))

  # Check for the lower bound of burnin
  testthat::expect_error(Sequential2Means(X, Y, b.i, "ridge", 200, 10))

  # Check if the number of burnin is a natural number
  testthat::expect_error(Sequential2Means(X, Y, b.i, "ridge", 200, 100.5213))



  ########### OUTPUT TESTS

  # Expecting Beta output to be of specific dimension
  testthat::expect_length(Sequential2Means(X, Y, b.i, "ridge", 200, 100)$Beta[1,], ncol(X))


  testthat::expect_length(Sequential2Means(X, Y, b.i, "ridge", 200, 100)$Beta[,1], 200)

  # Expecting output to be of predefined length for H.b.i
  testthat::expect_length((Sequential2Means(X, Y, b.i, "ridge", 200, 100)$H.b.i), length(b.i))

  # Expecting output to be a list containing Beta, and H.b.i

  testthat::expect_output(str(Sequential2Means(X, Y, b.i, "ridge", 200, 100)), "$ Beta", fixed = TRUE)

  testthat::expect_output(str(Sequential2Means(X, Y, b.i, "ridge", 200, 100)), "$ H.b.i", fixed = TRUE)





  #----------------------------------------------#
  ######### Sequential2MeansBeta   #########
  #----------------------------------------------#

  # Generating data
  n <- 50
  p <- 4
  X <- matrix(rnorm(n*p), n, p)
  beta <- exp(rnorm(p))
  Y <- X %*% beta + rnorm(n, 0, 1)
  df <- data.frame(X,Y)

  # sampling beta from shrinkage prior using gibbs sampling
  rv.hs <- bayesreg::bayesreg(Y~. ,df , "gaussian", "horseshoe+", 2000, 700)
  Beta = t(rv.hs$beta)

  lower = 0
  upper = 1
  l = 7




  ########### INPUT TESTS

  # Check for NULL or NaN values in Beta
  testthat::expect_error(Sequential2MeansBeta( NULL, lower, upper, l))

  # Check for matrix data type of Beta
  testthat::expect_error(Sequential2MeansBeta( Beta[,1], lower, upper, l))

  # Lower bound and Upper bound magnitude comparison check
  testthat::expect_error(Sequential2MeansBeta( Beta, lower = 1, upper= 0, l))

  # Check if l is NULL
  testthat::expect_error(Sequential2MeansBeta( Beta, lower, upper, NULL))

  # Check if l is natural number greater than zero
  testthat::expect_error(Sequential2MeansBeta( Beta, lower, upper, -1))


  ########### OUTPUT TESTS

  # Expecting output to be of predefined length for b.i
  testthat::expect_length((Sequential2MeansBeta( Beta, lower, upper, l)$b.i), l)

  # Expecting output to be of predefined length for H.b.i
  testthat::expect_length((Sequential2MeansBeta( Beta, lower, upper, l)$H.b.i), l)

  # Expecting output to be of fixed length for output p
  testthat::expect_length((Sequential2MeansBeta( Beta, lower, upper, l)$p), 1)

  # Expecting output to be a list containing p, b.i, and H.b.i

  testthat::expect_output(str(Sequential2MeansBeta( Beta, lower, upper, l)), "$ b.i", fixed = TRUE)

  testthat::expect_output(str(Sequential2MeansBeta( Beta, lower, upper, l)), "$ H.b.i", fixed = TRUE)

  testthat::expect_output(str(Sequential2MeansBeta( Beta, lower, upper, l)), "$ p", fixed = TRUE)








})
