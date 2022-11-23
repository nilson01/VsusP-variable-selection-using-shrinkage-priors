# VsusP

## Introduction

In the context of Gaussian linear models, the package proposes a solution of the variable selection problem using Shrinkage priors. First obtain a posterior distribution of the number of significant variables by clustering the significant variables and the noise coefficients and obtaining samples from this distribution using MCMC. Then, the significant variables are estimated from the posterior median of absolute beta. This approach requires only one tuning parameters and is applicable to continuous shrinkage priors as well.

## Installation

``` r
devtools::install_github("nilson01/VsusP")
```


