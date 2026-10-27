# VsusP

VsusP applies sequential two-means (S2M) to posterior coefficient draws from
shrinkage regression. The method of Li and Pati (2017) estimates a signal count
at each tuning threshold. Variables are then ranked by their posterior median
absolute coefficients.

<!-- badges: start -->
[![R-CMD-check](https://github.com/nilson01/VsusP-variable-selection-using-shrinkage-priors/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nilson01/VsusP-variable-selection-using-shrinkage-priors/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/nilson01/VsusP-variable-selection-using-shrinkage-priors/branch/main/graph/badge.svg)](https://app.codecov.io/gh/nilson01/VsusP-variable-selection-using-shrinkage-priors?branch=main)
<!-- badges: end -->


## Installation

Install the released version from CRAN:

``` r
install.packages("VsusP")
library(VsusP)
```

To install the development version from GitHub instead:

``` r
install.packages("devtools") # Only needed if devtools is not installed.
devtools::install_github("nilson01/VsusP-variable-selection-using-shrinkage-priors")
library(VsusP)
```

## Example

Use a matrix with posterior draws in rows and coefficients in columns. Exclude
the intercept. These fixed draws illustrate the selection steps:

``` r
Beta <- rbind(c(.01, .02, .03, .3, .31, 10, 11),
              c(-.01, -.02, -.03, -.3, -.31, -10, -11))
colnames(Beta) <- paste0("x", seq_len(ncol(Beta)))
set.seed(1)
result <- Sequential2MeansBeta(Beta, lower = .1, upper = .5, l = 2)
result$H.b.i
#> [1] 4 2
OptimalHbi(result$b.i, result$H.b.i)

H <- result$H.b.i[2]  # Use the second threshold for this illustration.
selected <- S2MVarSelection(Beta, H)
colnames(Beta)[selected]
#> [1] "x7" "x6"
S2MVarSelectionV1(result, H)  # Same indices, using stored summaries.
```

For Gaussian outcomes, `Sequential2Means(X, Y, b.i)` fits a model with
`bayesreg` and applies S2M. It supports ridge, lasso, horseshoe and horseshoe+
priors. For a separate `bayesreg` fit, use `Beta <- t(fit$beta)` before passing
its coefficient draws to VsusP. See `vignette("VsusPVignette")` and
`?Sequential2Means` for the fitting example and input requirements.

## Gaussian fitting and prediction

`FitRidge()`, `FitHorseshoe()` and `FitSpikeSlab()` fit Gaussian regression
and retain separate iteration and chain axes. For example:

``` r
X <- cbind(x1=c(-2,-1,0,1,2,0), x2=c(0,1,-1,0,1,-1))
y <- c(-3,-1,0,1,3,.5)
fit <- FitRidge(X,y,draws=1000,chains=2,seed=8)
summary(fit)
predict(fit,X[1:2,,drop=FALSE])
PlotCoefficients(fit)
Beta <- PosteriorDraws(fit)
path <- S2MSensitivity(Beta,c(.1,.5),clustering="exact")
summary(path)
```

`summary()` reports posterior coefficient means, standard deviations and
pointwise equal-tail intervals. `predict(...,include.noise=TRUE,seed=9)` adds
observation noise to the interval draws. `PosteriorDraws()` combines chains for
selection; the fitted object keeps them separate for diagnostics.

The native fits center inputs without scaling them. Their priors and sampling
controls are in `?FitGaussian`; they differ from the `bayesreg` defaults.
Ridge draws are independent. For horseshoe and spike-and-slab fits, inspect
chains and use adequate burn-in and retained draws. Short examples demonstrate
calls; they do not establish convergence.

| Task | Functions |
|---|---|
| Fit Gaussian regression | `FitRidge`, `FitHorseshoe`, `FitSpikeSlab` |
| Summarize and predict | `summary`, `coef`, `predict`, `PosteriorDraws` |
| Check chains | `ChainDiagnostics` (requires `posterior`) |
| Plot coefficients or chains | `PlotCoefficients`, `PlotTrace`, `PlotAutocorrelation` |
| Summarize exact inclusions | `SpikeSlabSelection`, `PlotInclusion` |
| Select or compare thresholds | `TwoMeansBeta`, `S2MSensitivity`, `plot` |
| Form a sparse predictor | `SAVS`, `DSS`, `predict` |

## Spike-and-slab draws

For spike-and-slab posterior draws, excluded coefficients must be exactly zero.
Use draws in rows and predictors in columns, without the intercept:

``` r
Beta.ss <- rbind(c(1, 0, 0), c(2, 0, 0), c(0, -1, 0), c(0, 0, 0))
colnames(Beta.ss) <- c("x1", "x2", "x3")
selection <- SpikeSlabSelection(Beta.ss)
selection$inclusion.probability
#>   x1   x2   x3
#> 0.50 0.25 0.00
selection$MPM
#> [1] 1
TwoMeansBeta(Beta.ss)$selected
#> [1] 1
path.ss <- Sequential2MeansBeta(Beta.ss, .1, .5, 2,
                                clustering = "exact", zero.draws = "empty")
```

These fixed values illustrate the calls. Replace `Beta.ss` with draws from a
spike-and-slab fit. `MPM` includes coefficients with inclusion frequency at
least 0.5. `HPPM` returns the most frequent complete model visited by the draws;
it does not certify the highest-probability model over all possible models.
Continuous-prior draws cannot supply these inclusion probabilities.

`TwoMeansBeta()` implements the paper's basic 2-M summary. For both clustering
summaries, an all-zero draw is assigned zero signals only by the documented
boundary convention. Other constant-magnitude draws remain invalid.

## Additional selection options

SAVS is a separate selector that can return no variables:

``` r
SAVS(diag(2), c(.1, -.1))$selected
#> integer(0)
```

For fitted Gaussian coefficients, use `SAVS(X, colMeans(fit$Beta))`, with the
same predictor matrix and column order used for fitting. Exclude the intercept.
SAVS depends on predictor and response units; standardize predictors before
fitting when a common predictor scale is needed.

DSS provides a sparse approximation to a posterior-mean linear predictor:

``` r
DSS(diag(2), c(2, .1), lambda = .2)$coefficients
#> [1] 1.9 0.0
```

Use `DSS(X, colMeans(fit$Beta), lambda)` with matching predictor columns.
This implements Hahn and Carvalho's adaptive-L1 approximation, not their exact
variable-count penalty. Choose `lambda` by inspecting the trade-off between
sparsity and prediction; no optimal value is selected automatically. Check
`converged` before using a result. Its `mse` describes approximation to the
supplied mean at `X`, not performance on new outcomes.

S2M also offers optional exact clustering and recursive-pair refinement:

``` r
path <- S2MSensitivity(Beta, c(.1, .5),
                      clustering = "exact", small.set = "refine")
path$H.b.i
path$mode.proportion
path$jaccard.previous
```

The default remains randomized `kmeans` with the legacy small-set rule.
Refinement can change counts. Sensitivity summaries describe stability; they
do not choose an optimal threshold or establish selection accuracy.

## Interpretation

`OptimalHbi()` plots the tuning path; it does not choose an optimal threshold.
The threshold is on the coefficient scale. Assess sensitivity to the threshold,
predictor scaling and posterior sampling before interpreting a selection.

By default, accepted S2M draws always retain at least one signal. They cannot automatically
select an empty model; SAVS and penalized DSS provide separate routes that can do so. `H = 0` is accepted by the two selection functions when
an empty selection is supplied explicitly. A draw with identical absolute
coefficients has no lower-mean cluster and is rejected. With the legacy S2M rule, candidate sets of at
most two elements, or with identical magnitudes, are retained without further
splitting. Clustering counts are not posterior inclusion probabilities or a
guarantee of false-positive control.

## Documentation

- [Reference manual](VsusP.pdf): function inputs, outputs and example figures.
- [Worked guide](vignettes/VsusPVignette.pdf): fitting, selection, prediction and numerical checks.
- In R, open the installed guide with `vignette("VsusPVignette")`.

## References

Li, H. and Pati, D. (2017). Variable selection using shrinkage priors.
*Computational Statistics & Data Analysis*, **107**, 107–119.
[doi:10.1016/j.csda.2016.10.008](https://doi.org/10.1016/j.csda.2016.10.008).

Ray, P. and Bhattacharya, A. (2018). *Signal Adaptive Variable Selector for the
Horseshoe Prior*. [arXiv:1810.09004](https://arxiv.org/abs/1810.09004).

Hahn, P. R. and Carvalho, C. M. (2015). Decoupling Shrinkage and Selection in
Bayesian Linear Models: A Posterior Summary Perspective.
*Journal of the American Statistical Association*, **110**(509), 435–448.
[doi:10.1080/01621459.2014.993077](https://doi.org/10.1080/01621459.2014.993077).

A reproducible, exploratory 60-dataset selection pilot is included under
`inst/validation`; its small sample does not establish false-discovery control.
