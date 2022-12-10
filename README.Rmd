# VsusP

## Introduction

In the context of Gaussian linear models, the package proposes a solution of the variable selection problem using Shrinkage priors. First obtain a posterior distribution of the number of significant variables by clustering the significant variables and the noise coefficients and obtaining samples from this distribution using MCMC. Then, the significant variables are estimated from the posterior median of absolute beta. This approach requires only one tuning parameters and is applicable to continuous shrinkage priors as well.

<!-- badges: start -->
[![R-CMD-check](https://github.com/nilson01/VsusP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nilson01/VsusP/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/nilson01/VsusP/branch/main/graph/badge.svg)](https://app.codecov.io/gh/nilson01/VsusP?branch=main)
<!-- badges: end -->


## Getting Started

These instructions will install the package to your computer.

### Prerequisites

In order to install the package directly from Github, you need to have the **devtools** package:

```
install.packages("devtools")
```

## Installation

To install the package from Github, first load **devtools**:

```
library(devtools)
```

Next, install **VsusP** as follows:

``` r
devtools::install_github("nilson01/VsusP")
```


The package can now be loaded into R and used:

```
library(VsusP)
```



## Example

To provide an example of the **VsusP** package, we will first simulate some data:

```
# install.packages("MASS") # if not yet installed
library(MASS) 
set.seed(20221208)
sim.XY <- function(n, p, beta){
  X <- matrix(rnorm(n*p), n, p)
  Y <- as.vector(X %*% beta + rnorm(n, 0, 1))
  return(list(X=X, Y=Y, beta=beta)  )
}

n = 100
p = 20
beta <- exp(rnorm(p))
data <- sim.XY(n , p, beta)
```

The function `Sequential2Means` is used get the samples of `Beta` and `H.b.i`. `Beta` can be calculated using different shrinkage priors and gibbs sampling technique, whereas H.b.i is the estimated number of signals corresponding to each tuning parameter `b.i` that is computed using S2M algorithm. If the `Beta` matrix is available prior, then H.b.i is recommended to be computed using `Sequential2MeansBeta`. 

```
b.i = seq(0, 1, 0.05)
S2M = Sequential2Means(X = data$X, Y = data$Y, b.i = b.i, prior = "horseshoe+", n.samples = 5000, burnin = 2000)
Beta = S2M$Beta
H.b.i = S2M$H.b.i
```

Then, using the result from `Sequential2Means` or `Sequential2MeansBeta`, appropriate tuning parameter can be estimated using b.i Vs H.b.i plot from `OptimalHbi` function. 

```
OptimalHbi(bi = b.i, Hbi = H.b.i)

```

Finally, the indices of important subset of variables can be obtained using `S2MVarSelection` function. The input parameter `H` is the optimal H.b.i value from the selected tuning parameter in `OptimalHbi` function. Also, `Beta` is the user defined matrix or the matrix derived from `Sequential2Means` function. 

```
H = 17 
impVariablesGLM = S2MVarSelection(Beta, H)
impVariablesGLM

```

There are several other variants of these functions as per input cases. For example: Sequential2MeansBeta and S2MVarSelectionV1. 
Check the documentation for the different functions on options and details, using e.g., `?Sequential2Means`.


## References

R Core Team (2022). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL <https://www.R-project.org/>.

Li, H. and Pati, D., 2017. Variable selection using shrinkage priors. Computational Statistics & Data Analysis, 107, pp.107-119.
