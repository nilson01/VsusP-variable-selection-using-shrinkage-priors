test_that("SAVS agrees with the closed form and returns empty selections", {
  X <- diag(c(2, 1, 4)); b <- c(1, -2, .1)
  # Norms squared are 4,1,16: the first two survive with shrinkage 1/4.
  f <- SAVS(X,b)
  expect_equal(f$coefficients,c(.75,-1.75,0))
  expect_identical(f$selected,c(1L,2L))
  expect_identical(SAVS(diag(2),c(.1,-.1))$selected,integer(0))
  expect_equal(SAVS(diag(2),c(0,0))$coefficients,c(0,0))
  expect_equal(SAVS(diag(2),c(1,-1))$coefficients,c(0,0))
  expect_identical(SAVS(matrix(2,1,1),0)$selected,integer(0))
  # A zero-norm column has no identified coefficient scale.
  expect_error(SAVS(matrix(0,2,1),1),"zero-norm")
})

test_that("SAVS handles names, invalid inputs and extreme finite scales", {
  X <- diag(2);colnames(X) <- c("a","b")
  expect_identical(names(SAVS(X,c(a=2,b=0))$coefficients),c("a","b"))
  expect_error(SAVS(X,c(b=0,a=2)),"names")
  expect_error(SAVS(X,c(1,2,3)),"coefficient")
  expect_error(SAVS(X,c(NA,1)),"beta")
  expect_error(SAVS(X,c(Inf,1)),"beta")
  expect_error(SAVS(matrix("a",2,2),c(1,2)),"X")
  expect_error(SAVS(matrix(Inf,2,2),c(1,2)),"X")
  huge <- SAVS(diag(c(1e200,1e-200)),c(1e-100,1e150))
  expect_equal(huge$coefficients,c(1e-100,1e150))
  expect_identical(huge$selected,c(1L,2L))
  tiny <- SAVS(matrix(1e200,1,1),1e-200)
  expect_equal(tiny$coefficients,0)
})

test_that("SAVS predictor-unit dependence is explicit", {
  X <- diag(2); b <- c(2,0)
  original <- SAVS(X,b)
  scaled.X <- sweep(X,2,c(10,.1),"*")
  scaled.b <- b / c(10,.1)
  expect_equal(as.vector(X %*% b),as.vector(scaled.X %*% scaled.b))
  expect_identical(original$selected,1L)
  expect_identical(SAVS(scaled.X,scaled.b)$selected,integer(0))
  expect_equal(SAVS(sweep(scaled.X,2,c(10,.1),"/"),
                         scaled.b*c(10,.1)),original)
})

test_that("exact two-cluster cuts minimize exhaustive small-input loss", {
  cases <- list(c(0,2),c(0,0,1,2),c(1,2,3),c(0,.01,.03,.3,10,11),
                c(0,0,1,1,2,2,9),c(.1,.2,.3,.4,.5,.6))
  for (x in cases) {
    fit <- exactTwoMeans(x)
    n <- length(x)
    assignments <- expand.grid(rep(list(c(1L,2L)),n-1L))
    loss <- apply(assignments,1,function(a) {
      labels <- c(1L,as.integer(a))
      if (length(unique(labels)) < 2L) return(Inf)
      sum(vapply(split(x,labels),function(v) sum((v-mean(v))^2),numeric(1)))
    })
    got <- sum(vapply(split(x,fit$cluster),function(v)sum((v-mean(v))^2),numeric(1)))
    expect_equal(got,min(loss),tolerance=1e-12)
    expect_true(fit$centers[1] < fit$centers[2])
  }
  # Equal optimal losses: earliest sorted cut, without RNG consumption.
  set.seed(9); saved <- .Random.seed
  expect_equal(exactTwoMeans(c(0,1,2))$cluster,c(1,2,2))
  expect_identical(.Random.seed,saved)
  expect_error(exactTwoMeans(c(1,1,1)),"distinct")
})

test_that("recursive-pair refinement respects the saved candidate and gap", {
  B <- matrix(c(0,.125,10),1)
  # Initial lower set is {0,.125}. Gap .125: retain at equality, split below.
  expect_equal(Sequential2MeansBeta(B,.125,.125,1,"exact","refine")$H.b.i,1)
  expect_equal(Sequential2MeansBeta(B,.1,.1,1,"exact","refine")$H.b.i,2)
  expect_equal(Sequential2MeansBeta(B,.1,.1,1,"exact","legacy")$H.b.i,1)
  expect_equal(Sequential2MeansBeta(matrix(c(.1,10),1),0,0,1,"exact","refine")$H.b.i,1)
  expect_error(Sequential2MeansBeta(B,0,1,2,"bad"),"clustering")
  expect_error(Sequential2MeansBeta(B,0,1,2,"exact","bad"),"small.set")
})

test_that("sensitivity reports frequencies and selected-set overlap", {
  x <- c(.01,.02,.03,.3,.31,10,11)
  z <- c(.01,.02,.03,3,3.1,10,11)
  B <- rbind(x,x,z)
  f <- S2MSensitivity(B,c(.1,.5),"exact")
  expect_equal(f$H.b.i,c(4,2))
  expect_equal(f$mode.proportion,c(1,2/3))
  expect_equal(f$jaccard.previous,c(NA,.5))
  expect_equal(f$frequencies[[2]]$draws,c(2L,1L))
  expect_identical(f$selected[[2]],c(7L,6L))
  expect_equal(S2MVarSelectionV1(f,2),c(7,6))
  one <- S2MSensitivity(B,.5,"exact")
  expect_true(is.na(one$jaccard.previous))
  expect_equal(S2MSensitivity(B,c(.5,.5),"exact")$jaccard.previous,c(NA,1))
  expect_error(S2MSensitivity(B,numeric(0)),"b.i")
})

test_that("the default S2M options preserve existing values and RNG paths", {
  B <- rbind(c(.01,.02,.03,.3,.31,10,11),c(-.01,-.02,-.03,-.3,-.31,-10,-11))
  set.seed(1);a <- Sequential2MeansBeta(B,.1,.5,2); rng <- .Random.seed
  set.seed(1);b <- Sequential2MeansBeta(B,.1,.5,2,"kmeans","legacy")
  expect_equal(a,b)
  expect_identical(.Random.seed,rng)
  expect_equal(a$H.b.i,c(4,2))
})

test_that("SAVS resolves nontrivial exact threshold equality", {
  expect_identical(SAVS(matrix(512,1,1),1/64)$selected,integer(0))
  expect_equal(SAVS(matrix(512,1,1),1/64)$coefficients,0)
  expect_identical(SAVS(matrix(512,1,1),(1/64)*(1+1e-8))$selected,1L)
  expect_identical(SAVS(matrix(512,1,1),(1/64)*(1-1e-8))$selected,integer(0))
})

test_that("SAVS agrees with direct published arithmetic away from its boundary", {
  set.seed(30)
  for (i in seq_len(12)) {
    X <- matrix(rnorm(28),7,4); b <- rnorm(4)
    mu <- 1/abs(b)^2
    expected <- sign(b)/colSums(X^2) * pmax(abs(b)*colSums(X^2)-mu,0)
    expect_equal(SAVS(X,b)$coefficients,expected,tolerance=1e-13)
  }
  X <- diag(2);colnames(X) <- c("a","a")
  expect_error(SAVS(X,c(1,2)),"unique")
  expect_error(exactTwoMeans(c(0,1e-300,1e300)),"collapsed")
  expect_error(exactTwoMeans(c(0,1e-200,1)),"separation")
  expect_identical(SAVS(matrix(2^900,1,1),2^-600)$selected,integer(0))
})

test_that("the fitting wrapper exposes the optional S2M conventions", {
  B <- matrix(c(0,.125,10),1)
  local_mocked_bindings(fitShrinkageGaussian=function(...)list(beta=t(B)),
                        .package="VsusP")
  X <- matrix(seq_len(60),20,3);Y <- seq_len(20)
  f <- Sequential2Means(X,Y,.1,"ridge",100,100,"exact","refine")
  expect_equal(f$H.b.i,2)
  expect_equal(f$H.b.i,Sequential2MeansBeta(B,.1,.1,1,"exact","refine")$H.b.i)
})
