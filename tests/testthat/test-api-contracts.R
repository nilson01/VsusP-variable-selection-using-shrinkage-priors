test_that("selection handles empty and complete sets and preserves tie order", {
  B <- rbind(c(1, 1, 3), c(-1, -1, -3))
  expect_identical(S2MVarSelection(B, 0), integer(0))
  expect_identical(S2MVarSelection(B, 3), c(3L, 2L, 1L))
  obj <- list(p = 3, abs.post.median = c(1, 1, 3))
  expect_identical(S2MVarSelectionV1(obj, 0), integer(0))
  expect_identical(S2MVarSelectionV1(obj, 3), c(3L, 2L, 1L))
  for (h in list(-1, .5, NA_real_, Inf, c(1, 2), numeric(0), "1", 4)) {
    expect_error(S2MVarSelection(B, h), "H")
    expect_error(S2MVarSelectionV1(obj, h), "H")
  }
})

test_that("both posterior-draw workflows provide usable selection summaries", {
  B <- rbind(c(.01, .02, .03, .3, .31, 10, 11),
             c(-.01, -.02, -.03, -.3, -.31, -10, -11))
  set.seed(1)
  obj <- Sequential2MeansBeta(B, .5, .5, 1)
  expect_equal(obj$abs.post.median, abs(B[1, ]))
  expect_identical(S2MVarSelectionV1(obj, obj$H.b.i), c(7L, 6L))
  expect_equal(Sequential2MeansBeta(matrix(c(.1, 10), 1), .5, .5, 1)$H.b.i, 1)
})

test_that("degenerate splits have explicit boundary behavior", {
  expect_equal(numNoiseCoeff(c(.1, 10), 0), 1)
  expect_equal(numNoiseCoeff(c(-.1, -10), 20), 1)
  expect_error(Sequential2MeansBeta(matrix(0, 2, 3), 0, 1, 2), "distinct")
  expect_error(Sequential2MeansBeta(matrix(1, 2, 1), 0, 1, 2), "two")
  # An indivisible lower cluster is retained rather than passed to kmeans.
  set.seed(1)
  expect_equal(numNoiseCoeff(c(1, 1, 1, 10, 11), .1), 3)
})

test_that("posterior draws and tuning inputs are checked before clustering", {
  B <- rbind(c(.1, .2, 10), c(-.1, -.2, -10))
  for (bad in list(matrix(Inf, 2, 3), matrix("a", 2, 3),
                   matrix(numeric(0), 0, 3), matrix(1+1i, 2, 3))) {
    expect_error(Sequential2MeansBeta(bad, 0, 1, 2), "Beta")
    expect_error(S2MVarSelection(bad, 1), "Beta")
  }
  for (bad in list(NA_real_, Inf, numeric(0), c(0, 1), "0", -1)) {
    expect_error(Sequential2MeansBeta(B, bad, 1, 2), "lower")
    expect_error(Sequential2MeansBeta(B, 0, bad, 2), "upper")
  }
  for (bad in list(NA_real_, Inf, numeric(0), c(1, 2), "2", 0, 1.5)) {
    expect_error(Sequential2MeansBeta(B, 0, 1, bad), "l")
  }
  expect_error(Sequential2MeansBeta(B, 2, 1, 2), "lower")
})

test_that("legacy summaries and plotting paths reject malformed inputs", {
  for (obj in list(1, list(p=3), list(p=3, abs.post.median=c(1,2)),
                   list(p=3, abs.post.median=c(1,2,Inf)),
                   list(p=3, abs.post.median=c(1,-2,3)),
                   list(p=c(2,3), abs.post.median=c(1,2,3)))) {
    expect_error(S2MVarSelectionV1(obj, 1), "S2M")
  }
  for (bad in list(numeric(0), c(1, Inf), c(-1, 1), c("a", "b"))) {
    expect_error(OptimalHbi(bad, c(1,2)), "bi")
    expect_error(OptimalHbi(c(0,1), bad), "Hbi")
  }
  expect_error(OptimalHbi(c(0,1), c(1,1.5)), "Hbi")
  plotfile <- tempfile(fileext=".pdf")
  grDevices::pdf(plotfile)
  on.exit({grDevices::dev.off(); unlink(plotfile)}, add=TRUE)
  expect_null(OptimalHbi(c(0,1), c(2,1)))
})

test_that("fitting rejects malformed data before starting a sampler", {
  X <- matrix(seq_len(60), 20, 3); Y <- seq_len(20)
  for (bad in list(matrix(Inf,20,3), matrix("a",20,3), matrix(1,20,1))) {
    expect_error(Sequential2Means(bad,Y,0), "X")
  }
  for (bad in list(rep(Inf,20), rep("a",20))) {
    expect_error(Sequential2Means(X,bad,0), "Y")
  }
  for (bad in list(numeric(0), Inf, -1, "0")) {
    expect_error(Sequential2Means(X,Y,bad), "b.i")
  }
  for (bad in list(NA_real_, Inf, c(100,200), numeric(0))) {
    expect_error(Sequential2Means(X,Y,0,n.samples=bad), "n.samples")
    expect_error(Sequential2Means(X,Y,0,burnin=bad), "burnin")
  }
})

test_that("documented prior aliases and fitted summary workflows agree", {
  # Single-worker sampling makes fixed-seed comparisons reproducible.
  sampler <- bayesreg::bayesreg
  local_mocked_bindings(fitShrinkageGaussian = function(...) sampler(..., n.cores = 1),
                        .package = "VsusP")
  set.seed(20)
  X <- matrix(rnorm(120), 40, 3)
  Y <- 3 * X[, 1] + rnorm(40)
  for (pair in list(c("rr","ridge"), c("hs","horseshoe"), c("hs+","horseshoe+"))) {
    set.seed(91)
    alias <- Sequential2Means(X,Y,c(.1,.5),pair[1],101,100)
    set.seed(91)
    canonical <- Sequential2Means(X,Y,c(.1,.5),pair[2],101,100)
    expect_equal(alias,canonical)
    expect_equal(alias$p,3)
    expect_equal(alias$b.i,c(.1,.5))
    expect_equal(alias$abs.post.median,apply(abs(alias$Beta),2,stats::median))
    expect_equal(S2MVarSelectionV1(alias,alias$H.b.i[1]),
                 S2MVarSelection(alias$Beta,alias$H.b.i[1]))
    expect_gte(nrow(alias$Beta),101)
    set.seed(5)
    supplied <- Sequential2MeansBeta(alias$Beta,.1,.5,2)
    expect_equal(alias$H.b.i,supplied$H.b.i)
  }
})

test_that("predictor names do not change the fitted response", {
  sampler <- bayesreg::bayesreg
  local_mocked_bindings(fitShrinkageGaussian = function(...) sampler(..., n.cores = 1),
                        .package = "VsusP")
  set.seed(50)
  X <- matrix(rnorm(120),40,3); Y <- X[,1] * 3 + rnorm(40)
  set.seed(12)
  a <- Sequential2Means(X,Y,.5,"ridge",100,100)
  colnames(X) <- c("Y", ".VsusP.response", "a b")
  set.seed(12)
  b <- Sequential2Means(X,Y,.5,"ridge",100,100)
  expect_equal(unname(a$Beta),unname(b$Beta))
  expect_equal(b$H.b.i,a$H.b.i)
  expect_identical(colnames(b$Beta),colnames(X))
})

test_that("fitting uses rounded-up sampler output when finding the mode", {
  two <- c(.01,.02,.03,.3,.31,10,11)
  four <- c(.01,.02,.03,3,3.1,10,11)
  B <- rbind(matrix(rep(four,51),51,7,byrow=TRUE),
             matrix(rep(two,51),51,7,byrow=TRUE))
  local_mocked_bindings(fitShrinkageGaussian = function(...) list(beta=t(B)),
                        .package="VsusP")
  set.seed(1)
  fit <- Sequential2Means(matrix(seq_len(140),20,7),seq_len(20),.5,
                          "ridge",101,100)
  expect_equal(nrow(fit$Beta),102)
  # All 102 draws tie between counts 2 and 4; the smaller mode is 2.
  # The first 101 alone would instead give mode 4.
  expect_equal(fit$H.b.i,2)
})
