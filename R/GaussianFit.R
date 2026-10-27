#' Fit native Gaussian shrinkage regression
#'
#' Fit ridge, horseshoe or exact-zero spike-and-slab regression while retaining
#' separate iteration and chain axes. Existing bayesreg wrappers are unchanged.
#' @param X Finite numeric predictor matrix, without an intercept column.
#'   At least two rows are required. Constant columns are not supported.
#' @param y Finite numeric response vector with one entry per row of X.
#' @param prior One of "ridge", "horseshoe" or "spike-slab".
#' @param draws Positive integer retained iterations per chain.
#' @param burnin Nonnegative integer discarded iterations per chain, additional
#'   to draws. Defaults to zero for ridge and 1000 for the sampled priors.
#' @param chains Positive integer number of separately seeded chains.
#' @param seed Nonnegative integer seed, or NULL to use the current R RNG.
#'   A supplied seed leaves the caller's RNG state unchanged.
#' @param prior.scale Positive ridge prior standard-deviation multiplier.
#'   Other priors require its default value of one.
#' @param feature.names Unique nonempty predictor names. Defaults to column
#'   names, or x1, x2, and so on.
#' @return A vsusp_fit list with beta (iteration by chain by predictor),
#'   intercept and sigma2 (iteration by chain), prior, feature.names, extra
#'   scale/inclusion draws and metadata. Draws count retained iterations.
#' @details The likelihood is Gaussian with a flat intercept prior. Predictors
#'   and response are centered, not standardized. The intercept is integrated
#'   during coefficient/variance updates and then sampled on the original scale.
#'   The residual-variance prior is IG(1.5,1.5), with density proportional to
#'   x^(-a-1) exp(-b/x). Ridge coefficients have conditional variance
#'   sigma2 * prior.scale^2. Horseshoe local and global scales independently
#'   follow half-Cauchy(0,1) priors. For spike-and-slab, inclusion probability
#'   follows Beta(1,15), slab variance multipliers follow IG(1.5,1.5), and
#'   excluded coefficients are exactly zero. R and Python use
#'   different random-number streams. Retained draws do not prove convergence.
#' @references Li, H. and Pati, D. (2017). Variable selection using shrinkage
#'   priors. Computational Statistics & Data Analysis, 107, 107--119.
#'   \doi{10.1016/j.csda.2016.10.008}.
#'
#'   Makalic, E. and Schmidt, D. F. (2016). A simple sampler for the horseshoe
#'   estimator. \doi{10.1109/LSP.2015.2503725}.
#'
#'   Bhattacharya, A., Chakraborty, A. and Mallick, B. K. (2016).
#'   Fast sampling with Gaussian scale-mixture priors in high-dimensional
#'   regression. \doi{10.1093/biomet/asw042}.
#' @seealso [PosteriorDraws()], [ChainDiagnostics()]
#' @export
#' @examples
#' X <- cbind(x1 = c(-2,-1,0,1,2,0), x2 = c(0,1,-1,0,1,-1))
#' fit <- FitRidge(X, c(-4,-2,0,2,4,0), draws=100, chains=2, seed=8)
#' summary(fit)
#' predict(fit, X[1:2,,drop=FALSE])
FitGaussian <- function(X, y, prior = c("ridge", "horseshoe", "spike-slab"),
                        draws = 1000L, burnin = NULL, chains = 4L, seed = NULL,
                        prior.scale = 1, feature.names = colnames(X)) {
  prior <- matchChoice(prior, c("ridge", "horseshoe", "spike-slab"), "prior")
  checkNumericMatrix(X, "X"); checkNumericVector(y, "y")
  if (nrow(X) < 2L || length(y) != nrow(X)) stop("X needs at least two rows and y must match its rows.")
  if (is.null(burnin)) burnin <- if (prior == "ridge") 0L else 1000L
  for (name in c("draws", "burnin", "chains")) {
    value <- get(name); checkIntegerScalar(value, name, if (name == "burnin") 0 else 1)
    if (value > .Machine$integer.max) stop(name, " exceeds the integer limit.")
  }
  draws <- as.double(draws); burnin <- as.double(burnin); chains <- as.double(chains)
  if (draws + burnin > .Machine$integer.max || draws * chains * ncol(X) > .Machine$integer.max)
    stop("Requested draw dimensions exceed the integer limit.")
  checkNumericScalar(prior.scale, "prior.scale")
  if (prior.scale <= 0 || (prior != "ridge" && prior.scale != 1))
    stop("prior.scale must be positive and is configurable only for ridge.")
  if (prior == "ridge" && burnin != 0) stop("Independent ridge draws require burnin=0.")
  feature.names <- vsuspNames(feature.names, ncol(X))
  n <- nrow(X); p <- ncol(X); xm <- colMeans(X); ym <- mean(y)
  Xc <- sweep(X, 2L, xm, "-"); yc <- y - ym; sums <- colSums(Xc^2)
  vsuspFinite(c(Xc, yc, sums, sum(yc^2)), "Centered inputs")
  if (any(sums == 0)) stop("Constant predictor columns are not supported; omit the intercept column.")
  result <- vsuspWithSeed(seed, function() {
    seeds <- sample.int(.Machine$integer.max, chains, replace = FALSE)
    beta <- array(NA_real_, c(draws, chains, p), list(NULL, NULL, feature.names))
    intercept <- sigma.draws <- matrix(NA_real_, draws, chains)
    extra <- list()
    if (prior == "horseshoe") extra <- list(lambda2=beta, tau2=sigma.draws)
    if (prior == "spike-slab") extra <- list(omega2=beta, q=sigma.draws,
                                            included=array(FALSE, dim(beta), dimnames(beta)))
    ridge.factor <- if (prior == "ridge") vsuspGaussianFactor(Xc, yc, rep(prior.scale^2, p)) else NULL
    for (chain in seq_len(chains)) {
      set.seed(seeds[chain]); lambda2 <- nu <- omega2 <- rep(1, p)
      tau2 <- xi <- sigma2 <- 1; q <- 1/16; b <- numeric(p); z <- rep(FALSE, p)
      for (iteration in seq_len(draws + burnin)) {
        if (prior != "spike-slab") {
          d <- if (prior == "ridge") rep(prior.scale^2, p) else tau2 * lambda2
          block <- vsuspGaussianBlock(Xc, yc, d, factor=ridge.factor)
          b <- block$beta; sigma2 <- block$sigma2
          if (prior == "horseshoe") {
            lambda2 <- vsuspIG(rep(1,p), 1/nu + (b/sqrt(sigma2*tau2))^2/2, "lambda2")
            nu <- vsuspIG(rep(1,p), 1+1/lambda2, "nu")
            tau2 <- vsuspIG((p+1)/2, 1/xi + sum((b/sqrt(lambda2))^2)/(2*sigma2), "tau2")
            xi <- vsuspIG(1, 1+1/tau2, "xi")
          }
        } else {
          residual <- as.vector(yc - Xc %*% b)
          for (j in seq_len(p)) {
            r <- residual + Xc[,j]*b[j]; t <- sum(Xc[,j]*r)
            precision <- sums[j] + 1/omega2[j]; variance <- sigma2/precision
            log.bf <- -.5*log1p(omega2[j]*sums[j]) + .5*((t/sqrt(precision))/sqrt(sigma2))^2
            vsuspFinite(c(t, log.bf), "Coordinate conditional")
            vsuspFinite(c(precision,variance), "Coordinate variance", TRUE)
            z[j] <- stats::runif(1) < stats::plogis(stats::qlogis(q)+log.bf)
            b[j] <- if (z[j]) stats::rnorm(1,t/precision,sqrt(variance)) else 0
            if (!is.finite(b[j]) || (z[j] && b[j] == 0)) stop("Included coefficient is nonfinite or rounded to zero.")
            residual <- r - Xc[,j]*b[j]
          }
          omega2 <- vsuspIG(1.5+.5*z, 1.5+.5*(b/sqrt(sigma2))^2, "omega2")
          active <- sum(z); q <- stats::rbeta(1,1+active,15+p-active)
          if (!is.finite(q) || q <= 0 || q >= 1) stop("Inclusion probability reached a numerical boundary.")
          residual <- as.vector(yc-Xc%*%b)
          rate <- 1.5+.5*(sum(residual^2)+sum((b[z]/sqrt(omega2[z]))^2))
          sigma2 <- vsuspIG(1.5+(n-1+active)/2,rate,"sigma2")
        }
        if (iteration > burnin) {
          k <- iteration-burnin; alpha <- stats::rnorm(1,ym-sum(xm*b),sqrt(sigma2/n))
          vsuspFinite(alpha,"Intercept")
          beta[k,chain,] <- b; sigma.draws[k,chain] <- sigma2; intercept[k,chain] <- alpha
          if (prior == "horseshoe") { extra$lambda2[k,chain,] <- lambda2; extra$tau2[k,chain] <- tau2 }
          if (prior == "spike-slab") {
            extra$omega2[k,chain,] <- omega2; extra$q[k,chain] <- q; extra$included[k,chain,] <- z
          }
        }
      }
    }
    structure(list(beta=beta, intercept=intercept, sigma2=sigma.draws,
      prior=prior, feature.names=feature.names, extra=extra,
      metadata=list(n=n,p=p,draws=draws,burnin=burnin,chains=chains,chain.seeds=seeds,
        iteration=seq.int(burnin+1,draws+burnin),x.mean=xm,y.mean=ym,
        centered=TRUE,standardized=FALSE,prior.scale=if(prior=="ridge") prior.scale else NULL,
        sigma2.prior="IG(1.5,1.5)",intercept.prior="flat",
        RNG.kind=RNGkind(),convergence="not established by successful sampling")), class="vsusp_fit")
  })
  checkVsuspFit(result)
  result
}

#' @rdname FitGaussian
#' @export
FitRidge <- function(X,y,draws=1000L,chains=4L,prior.scale=1,seed=NULL,
                     feature.names=colnames(X)) {
  FitGaussian(X,y,"ridge",draws,0L,chains,seed,prior.scale,feature.names)
}

#' @rdname FitGaussian
#' @export
FitHorseshoe <- function(X,y,draws=1000L,burnin=1000L,chains=4L,seed=NULL,
                         feature.names=colnames(X)) {
  FitGaussian(X,y,"horseshoe",draws,burnin,chains,seed,1,feature.names)
}

#' @rdname FitGaussian
#' @export
FitSpikeSlab <- function(X,y,draws=1000L,burnin=1000L,chains=4L,seed=NULL,
                         feature.names=colnames(X)) {
  FitGaussian(X,y,"spike-slab",draws,burnin,chains,seed,1,feature.names)
}

vsuspFinite <- function(x,label,positive=FALSE) {
  if (any(!is.finite(x)) || (positive && any(x<=0)))
    stop(label," exceeded ",if(positive) "positive " else "","finite numeric range.",call.=FALSE)
  invisible(x)
}

vsuspNames <- function(value,p) {
  if (is.null(value)) value <- paste0("x",seq_len(p))
  if (!is.character(value) || length(value)!=p || anyNA(value) || any(!nzchar(value)) || anyDuplicated(value))
    stop("feature.names must have one unique nonempty name per predictor.",call.=FALSE)
  value
}

vsuspWithSeed <- function(seed,fun) {
  if (is.null(seed)) return(fun())
  checkIntegerScalar(seed,"seed"); if (seed>.Machine$integer.max) stop("seed exceeds the integer limit.")
  had <- exists(".Random.seed",.GlobalEnv,inherits=FALSE)
  if(had) old <- get(".Random.seed",.GlobalEnv)
  on.exit(if(had) assign(".Random.seed",old,.GlobalEnv) else
    if(exists(".Random.seed",.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(seed); fun()
}

vsuspIG <- function(shape,rate,label) {
  vsuspFinite(rate,paste(label,"rate"),TRUE)
  gamma <- stats::rgamma(length(rate),shape=shape)
  vsuspFinite(gamma,paste(label,"gamma draw"),TRUE)
  value <- rate/gamma; vsuspFinite(value,label,TRUE); value
}

vsuspGaussianFactor <- function(X,y,d) {
  vsuspFinite(d,"Prior variance",TRUE)
  C <- tcrossprod(sweep(X,2L,sqrt(d),"*"))+diag(nrow(X))
  vsuspFinite(C,"Gaussian covariance")
  U <- chol(C); condition <- rcond(C)
  if (!is.finite(condition) || condition<=1e-12) stop("Gaussian covariance is numerically ill-conditioned.")
  solveC <- function(rhs) backsolve(U,forwardsolve(t(U),rhs))
  quadratic <- sum(y*solveC(y)); vsuspFinite(quadratic,"Marginal quadratic")
  if (quadratic<0) stop("Marginal quadratic became negative.")
  list(U=U,rate=1.5+quadratic/2)
}

vsuspGaussianBlock <- function(X,y,d,z=NULL,delta=NULL,sigma2=NULL,factor=NULL) {
  if (is.null(factor)) factor <- vsuspGaussianFactor(X,y,d)
  if (is.null(sigma2)) sigma2 <- vsuspIG(1.5+(nrow(X)-1)/2,factor$rate,"sigma2")
  if (is.null(z)) z <- stats::rnorm(ncol(X))
  if (is.null(delta)) delta <- stats::rnorm(nrow(X))
  u <- sqrt(d)*z
  rhs <- y/sqrt(sigma2)-as.vector(X%*%u)-delta
  w <- backsolve(factor$U,forwardsolve(t(factor$U),rhs))
  beta <- sqrt(sigma2)*(u+d*as.vector(crossprod(X,w)))
  vsuspFinite(beta,"Gaussian coefficients")
  list(beta=beta,sigma2=sigma2,rate=factor$rate)
}

checkVsuspFit <- function(object) {
  if (!inherits(object,"vsusp_fit")) stop("object must be a vsusp_fit.")
  a <- object$beta; dims <- dim(a)
  if (!is.numeric(a) || is.complex(a) || length(dims)!=3L || any(dims<1) || any(!is.finite(a)))
    stop("beta must retain finite iteration, chain and predictor axes.")
  for (name in c("intercept","sigma2")) {
    b <- object[[name]]
    if (!is.matrix(b) || !is.numeric(b) || is.complex(b) || !identical(dim(b),dims[1:2]) || any(!is.finite(b)))
      stop(name," must match the iteration and chain axes.")
  }
  if (any(object$sigma2<=0)) stop("sigma2 draws must be positive.")
  vsuspNames(object$feature.names,dims[3])
  invisible(object)
}
