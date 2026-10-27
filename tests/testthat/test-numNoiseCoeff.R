# Li and Pati (2017), final publication, section 2.2.1:
# initialize D to the lower-mean cluster; save D before reclustering;
# when the gap no longer exceeds b, return |D| rather than splitting again.
# DOI: 10.1016/j.csda.2016.10.008

test_that("S2M retains the saved noise set at recursive stopping", {
  x <- c(.01, .02, .03, .3, .31, 10, 11)
  # Initial split: five small values versus {10, 11}.
  # Next gap: mean(c(.3, .31)) - mean(c(.01, .02, .03)) = .285.
  set.seed(1)
  expect_equal(numNoiseCoeff(x, .5), 5)
  # Splitting the three smallest values gives gap .015, hence |D| = 3.
  set.seed(1)
  expect_equal(numNoiseCoeff(x, .1), 3)
  set.seed(1)
  expect_equal(numNoiseCoeff(-x, .5), 5)
})

test_that("S2M initializes D and stops when the gap equals b", {
  set.seed(1)
  expect_equal(numNoiseCoeff(c(.01, .02, .03, .3, .31, 10, 11), 20), 5)
  set.seed(1)
  expect_equal(numNoiseCoeff(c(0, 2, 10, 12), 10), 2)
  # Binary fractions give an exactly representable recursive center gap.
  set.seed(1)
  expect_equal(numNoiseCoeff(c(0, .125, .25, 2, 2.125, 16, 18), 1.9375), 5)
})

test_that("S2M preserves the guard for one- and two-member candidates", {
  set.seed(1)
  # Two distinct observations define two singleton clusters.
  expect_equal(numNoiseCoeff(c(.01, 10), .1), 1)
  set.seed(1)
  expect_equal(numNoiseCoeff(c(.01, 10, 11), .1), 1)
  set.seed(1)
  expect_equal(numNoiseCoeff(c(.01, .02, 10), .1), 2)
  # A two-member candidate reached after a recursive split.
  set.seed(1)
  expect_equal(numNoiseCoeff(c(.01, .02, .3, .31, 10, 11), .1), 2)
})

test_that("posterior draws produce the corrected modal signal counts", {
  x <- c(.01, .02, .03, .3, .31, 10, 11)
  Beta <- rbind(x, -x, x, -x)
  set.seed(1)
  result <- Sequential2MeansBeta(Beta, lower = .1, upper = .5, l = 2)
  expect_equal(result$p, 7)
  expect_equal(result$b.i, c(.1, .5))
  expect_equal(result$H.b.i, c(4, 2))
  expect_equal(S2MVarSelection(Beta, result$H.b.i[1]), c(7, 6, 5, 4))
  expect_equal(S2MVarSelection(Beta, result$H.b.i[2]), c(7, 6))
  set.seed(1)
  expect_equal(Sequential2MeansBeta(Beta, 20, 20, 1)$H.b.i, 2)
  # Three draws have two signals; two draws have four. Their mode is two.
  other <- c(.01, .02, .03, 3, 3.1, 10, 11)
  mixed_draws <- rbind(x, -x, x, other, -other)
  set.seed(1)
  expect_equal(Sequential2MeansBeta(mixed_draws, .5, .5, 1)$H.b.i, 2)
})
