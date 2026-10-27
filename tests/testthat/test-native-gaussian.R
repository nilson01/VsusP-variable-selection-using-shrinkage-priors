test_that("Gaussian sampling block has the analytic conditional mean and covariance", {
  data <- gaussian_example()
  for (X in list(data$X, cbind(data$X, data$X[,1], -data$X[,2], 2*data$X))) {
    X <- scale(X,center=TRUE,scale=FALSE); y <- data$y-mean(data$y)
    d <- seq(.3,1.7,length.out=ncol(X)); p <- ncol(X); n <- nrow(X)
    V <- solve(crossprod(X)+diag(1/d,p))
    mu <- as.vector(V%*%crossprod(X,y)); s2 <- 1.7
    center <- vsuspGaussianBlock(X,y,d,z=numeric(p),delta=numeric(n),sigma2=s2)
    expect_equal(center$beta,mu,tolerance=1e-12)
    # Propagate every independent standard-normal basis vector.
    basis <- diag(p+n)
    L <- vapply(seq_len(p+n),function(k)
      vsuspGaussianBlock(X,y,d,z=basis[seq_len(p),k],
        delta=basis[p+seq_len(n),k],sigma2=s2)$beta-center$beta,numeric(p))
    expect_equal(tcrossprod(L),unname(s2*V),tolerance=1e-12)
    expected.rate <- 1.5+as.numeric(crossprod(y,y)-crossprod(mu,crossprod(X,y)))/2
    expect_equal(center$rate,expected.rate,tolerance=1e-12)
  }
})

test_that("ridge joint draws agree with conjugate posterior moments", {
  d <- gaussian_example(); X <- scale(d$X,center=TRUE,scale=FALSE)
  y <- d$y-mean(d$y); V <- solve(crossprod(X)+diag(2))
  mu <- as.vector(V%*%crossprod(X,y)); shape <- 1.5+(nrow(X)-1)/2
  rate <- 1.5+as.numeric(crossprod(y,y)-crossprod(mu,crossprod(X,y)))/2
  fit <- FitRidge(d$X,d$y,draws=3000,chains=2,seed=921)
  b <- PosteriorDraws(fit); m <- nrow(b); variance <- diag(V)*rate/(shape-1)
  expect_true(all(abs(colMeans(b)-mu)<6*sqrt(variance/m)))
  expect_equal(colMeans(b^2)-colMeans(b)^2,variance,tolerance=.08)
  expect_lt(abs(mean(fit$sigma2)-rate/(shape-1)),
            6*sqrt(rate^2/((shape-1)^2*(shape-2))/m))
  expect_lt(abs(mean(fit$intercept)-mean(d$y)),.04)
})

test_that("native fits retain axes, exact spikes and reproducible separate chains", {
  d <- gaussian_example()
  set.seed(88); before <- .Random.seed
  for (prior in c("ridge","horseshoe","spike-slab")) {
    f <- FitGaussian(d$X,d$y,prior,draws=25,burnin=if(prior=="ridge")0 else 10,
                     chains=2,seed=3)
    g <- FitGaussian(d$X,d$y,prior,draws=25,burnin=if(prior=="ridge")0 else 10,
                     chains=2,seed=3)
    expect_identical(f,g); expect_identical(.Random.seed,before)
    expect_s3_class(f,"vsusp_fit"); expect_equal(dim(f$beta),c(25,2,2))
    expect_equal(dim(PosteriorDraws(f)),c(50,2))
    expect_true(all(is.finite(f$beta))); expect_true(all(f$sigma2>0))
    expect_false(identical(f$beta[,1,],f$beta[,2,]))
    expect_equal(length(unique(f$metadata$chain.seeds)),2)
    expect_equal(coef(f),colMeans(PosteriorDraws(f)))
    expect_output(print(f),"Gaussian")
    if(prior=="horseshoe") {
      expect_true(all(f$extra$lambda2>0));expect_true(all(f$extra$tau2>0))
    }
    if(prior=="spike-slab") {
      expect_identical(f$beta==0,!f$extra$included)
      expect_true(all(f$extra$q>0 & f$extra$q<1));expect_true(all(f$extra$omega2>0))
    }
  }
  expect_identical(FitHorseshoe(d$X,d$y,2,2,1,3),
                   FitGaussian(d$X,d$y,"horseshoe",2,2,1,3))
  expect_identical(FitSpikeSlab(d$X,d$y,2,2,1,3),
                   FitGaussian(d$X,d$y,"spike-slab",2,2,1,3))
  rm(".Random.seed",envir=.GlobalEnv)
  FitRidge(d$X,d$y,draws=2,chains=1,seed=4)
  expect_false(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
  assign(".Random.seed",before,envir=.GlobalEnv)
})

test_that("invalid native inputs fail before sampling", {
  d <- gaussian_example()
  expect_error(FitRidge(cbind(d$X,constant=1),d$y),"Constant")
  expect_error(FitRidge(d$X,d$y[-1]),"match")
  expect_error(FitRidge(d$X,d$y,feature.names=c("a","a")),"unique")
  expect_error(FitRidge(d$X,d$y,draws=2^31),"integer")
  expect_error(FitHorseshoe(d$X,d$y,draws=.Machine$integer.max,burnin=1L),"integer")
  expect_error(FitRidge(d$X,d$y,draws=100000L,chains=100000L),"integer")
  expect_error(FitRidge(d$X,d$y,draws=.5),"integer")
  expect_error(FitGaussian(d$X,d$y,burnin=1),"burnin")
  expect_error(FitGaussian(d$X,d$y,"horseshoe",prior.scale=2),"prior.scale")
  expect_error(FitRidge(d$X,d$y,seed=-1),"seed")
  expect_error(FitRidge(d$X*Inf,d$y),"finite")
  expect_error(FitRidge(d$X+1i,d$y),"finite")
  expect_error(FitRidge(d$X,d$y,chains=Inf),"finite")
  expect_error(FitRidge(d$X,d$y,prior.scale=0),"prior.scale")
})

test_that("posterior summaries and predictions use the retained draws exactly", {
  d <- gaussian_example(); f <- FitRidge(d$X,d$y,4,2,seed=2)
  b <- PosteriorDraws(f); s <- summary(f,interval=.8)
  expect_equal(s$mean,unname(colMeans(b)))
  expect_equal(s$sd,unname(apply(b,2,sd)))
  expect_equal(s$lower,unname(apply(b,2,quantile,.1)))
  v <- sweep(b%*%t(d$X[1:2,]),1,as.vector(f$intercept),"+")
  pred <- predict(f,d$X[1:2,],interval=.8)
  expect_equal(pred$mean,colMeans(v))
  expect_equal(pred$lower,unname(apply(v,2,quantile,.1)))
  expect_equal(pred$upper,unname(apply(v,2,quantile,.9)))
  expect_identical(pred$kind,"mean")
  noisy <- predict(f,d$X[1:2,],include.noise=TRUE,seed=71)
  expect_equal(noisy$mean,pred$mean)
  expect_identical(noisy,predict(f,d$X[1:2,],include.noise=TRUE,seed=71))
  expect_identical(noisy$kind,"observation")
  expect_false(identical(noisy$lower,pred$lower))
  expect_error(predict(f,d$X[,2:1]),"order")
  expect_error(predict(f,matrix(1,2,3)),"count")
  expect_error(predict(f,d$X,include.noise=NA),"TRUE")
  expect_error(summary(f,interval=1),"between")
  f$beta[1,1,1] <- Inf;expect_error(PosteriorDraws(f),"finite")
  one <- FitRidge(d$X,d$y,1,1,seed=5)
  expect_true(all(is.na(summary(one)$sd)))
  expect_equal(predict(one,d$X[1,,drop=FALSE])$lower,
               predict(one,d$X[1,,drop=FALSE])$upper)
})
