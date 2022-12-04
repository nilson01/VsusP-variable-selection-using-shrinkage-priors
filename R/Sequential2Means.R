# sourceCpp("Sequential2MeansC.cpp")

#' Variable selection using shrinkage priors:: Sequential2Means
#'
#' @param Beta.i ith row sample of N by p matrix consisting of N posterior samples of p variables
#' @param b.i_r tuning parameter value from Sequential 2-means (S2M) variable selection algorithm.
#'
#' @return number of noise coefficients
#'
numNoiseCoeff <- function(Beta.i, b.i_r) {
  # perform k means for sequential 2 means to prune the noise coefficients
  fit <- stats::kmeans(abs(Beta.i), 2)
  # perform sequential 2 means on subsequent clusters until threshold is reached

  # check the threshold between two centers
  while (max(fit$centers) - min(fit$centers) > b.i_r) {
    Beta.i <- Beta.i[which(fit$cluster == which.min(fit$centers))]
    if (length(Beta.i) <= 2) {
      break
    }
    fit <- stats::kmeans(abs(Beta.i), 2)
  }

  argminCenter <- which.min(fit$centers)
  return(length(which(fit$cluster == argminCenter)))
}


                                  ###########################
                                  ###########################
                           #########   Sequential2Means    #########
                                  ###########################
                                  ###########################


#' Variable selection using shrinkage priors:: Sequential2Means
#'
#' Sequential2Means function will take as input X: design matrix, Y : response vector, t: vector of tuning parameter values from Sequential 2-means (S2M) variable selection algorithm. The function will return a list S2M which will hold p: the total number of variables, b.i: the values of the tuning parameter, H.b.i : the estimated number of signals corresponding to each b.i, abs.post.median: medians of the absolute values of the posterior samples.
#'
#' @param X Design matrix
#' @param Y Response vector
#' @param b.i Vector of tuning parameter values from Sequential 2-means (S2M) variable selection algorithm.
#'
#' @return A list S2M which will hold Beta: N by p matrix consisting of N posterior samples of p variables, b.i: the values of the tuning parameter, H.b.i : the estimated number of signals corresponding to each b.i, abs.post.median: medians of the absolute values of the posterior samples.
#' @export
#'
#' @examples
Sequential2Means <- function(X, Y, b.i) {
  # Initializing variables
  ##########################

  # number data points
  N <- nrow(X)

  # number of covariates
  R <- ncol(X)

  # number of tuning parameters
  l <- length(b.i)

  iteration = 5000
  p = R

  # initializing the number of signals corresponding to each b.i
  H.b.i <- seq(0, length = length(b.i))

  # N by p matrix consisting of N posterior samples of p variables
  Beta <- matrix(0, nrow = iteration+1, ncol = R)
  Beta[1,] = rep(0.001,R)

  # .......TODO:: USE horse shoe prior and generate beta samples using MCMC.......

  phisq=rep(1,R)
  tausq=1
  yita=rep(1,R)
  n=N

  sigmasq=rep(0,iteration+1)

  #sigmasq~ig(1.5,1.5)
  ta=1.5
  tb=1.5

  sigmasq[1]=1

  for (ii in 1:iteration){
    #generate tausq
    bb=Beta[ii,]
    temp=1/tausq
    u=stats::runif(1,0,1/(1+temp))
    temp1=sum(0.5*(bb^2/(phisq*sigmasq[ii])))
    uu=stats::runif(1,0,stats::pgamma((1-u)/u,(n+1)/2,scale=temp1))
    temp=stats::qgamma(uu,(n+1)/2,scale=temp1)
    tausq=min(1,1/temp)

    # generate phisq
    yita=1/phisq
    u=stats::runif(R,0,1/(1+yita))
    temp=bb^2/(2*tausq*sigmasq[ii])
    uu=stats::runif(R,0,stats::pexp((1-u)/u,temp))
    yita=stats::qexp(uu,temp)
    phisq=1/yita

    # generate beta(j)
    D=sigmasq[ii]*tausq*(phisq)
    Phi=X/sqrt(sigmasq[ii])

    tempu=stats::rnorm(R,0,sd=sqrt((D)))
    tempv=Phi%*%tempu+stats::rnorm(N,0,1)

    tempw=solve((Phi)%*%(D*t(Phi))+diag(N))%*%(Y/sqrt(sigmasq[ii])-tempv)

    Beta[ii+1,]=(tempu+(D*t(Phi))%*%tempw)

    # updating sigma^2
    tempan=ta+(N/2)+(R/2)
    tempbn=tb+0.5*(t(Y-X%*%Beta[ii+1,])%*%(Y-X%*%Beta[ii+1,]))+0.5*sum(Beta[ii+1,]^2/phisq)/tausq
    sigmasq[ii+1]=1/stats::rgamma(1,tempan,rate=tempbn)

  }

  Beta = Beta[2000:5000,]

  # Sequential 2 means algorithm
  for (r in 1:l) {
    impVars.count <- seq(0, length = N)
    for (i in 1:N) {
      impVars.count[i] <- p - numNoiseCoeff(Beta[i, ], b.i[r])
    }
    H.b.i[r] <- as.numeric(names(sort(-table(impVars.count)))[1])
  }

  # list to hold desired output
  S2M <- list(Beta, H.b.i)

  return(S2M)
}



                          ##########################
                          ##########################
                  ########   Sequential2MeansBeta    #########
                          ##########################
                          ##########################

#' Variable selection using shrinkage prior:: Sequential2Means
#'
#' Sequential2MeansBeta function will take as input Beta : N by p matrix consisting of N posterior samples of p variables, lower : the lower bound of the chosen values of the tuning parameter, upper : the upper bound of the chosen values of the tuning parameter, and l :the number of chosen values of the tuning parameter. The function will return a list S2M which will hold p: the total number of variables, b.i: the values of the tuning parameter, H.b.i : the estimated number of signals corresponding to each b.i, abs.post.median: medians of the absolute values of the posterior samples.
#'
#' @param Beta N by p matrix consisting of N posterior samples of p variables
#' @param lower the lower bound of the chosen values of the tuning parameter
#' @param upper the upper bound of the chosen values of the tuning parameter
#' @param l the number of chosen values of the tuning parameter
#'
#' @return A list S2M which will hold p: the total number of variables, b.i: the values of the tuning parameter, H.b.i : the estimated number of signals corresponding to each b.i, abs.post.median: medians of the absolute values of the posterior samples.
#' @export
#'
Sequential2MeansBeta <- function(Beta, lower, upper, l) {
  # Initializing variables
  ##########################

  # number data points
  N <- nrow(Beta)

  # number of covariates
  p <- ncol(Beta)

  # tuning parameters
  b.i <- seq(lower, upper, length = l)

  # initializing the number of signals as zero for corresponding b.i
  H.b.i <- seq(0, length = l)

  # Sequential 2 means algorithm
  for (r in 1:l) {
    impVars.count <- seq(0, length = N)
    for (i in 1:N) {
      impVars.count[i] <- p - numNoiseCoeff(Beta[i, ], b.i[r])
    }
    H.b.i[r] <- as.numeric(names(sort(-table(impVars.count)))[1])
  }

  # list to hold desired output
  S2M <- list(p, b.i, H.b.i)
  return(S2M)
}
