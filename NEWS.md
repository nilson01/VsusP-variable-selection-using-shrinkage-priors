# VsusP 1.0.0.9000

## Gaussian fitting and result methods

* Add native `FitRidge()`, `FitHorseshoe()` and `FitSpikeSlab()` samplers with
  explicit priors, intercept draws, retained chain axes and reproducible seeds.
  Existing `bayesreg` fitting calls and prior defaults are unchanged.
* Add `PosteriorDraws()`, coefficient summaries, posterior prediction and sparse
  prediction methods. Existing selection-result fields retain their values;
  added S3 classes provide summaries and plots.
* Add `ChainDiagnostics()` using the suggested `posterior` package for rank
  R-hat, bulk/tail effective sample size and mean Monte Carlo standard error.
* Add coefficient, inclusion, trace and autocorrelation plots. Plot functions
  return the displayed values invisibly.

## Bug fixes

* Correct the S2M noise count at recursive stopping to retain the saved
  candidate set, as specified in Li and Pati (2017), section 2.2.1
  (doi:10.1016/j.csda.2016.10.008). Previously, the smaller cluster from the
  stopping split was counted instead. This could overstate the number of
  selected signals. General-case kmeans calls and the strict center-gap
  comparison are retained.

* Handle an initial pair of distinct coefficients; retain an indivisible
  recursive candidate. Reject constant-magnitude initial draws with a clear
  error rather than assigning an unsupported signal count.
* Return an empty integer vector for H = 0 and validate scalar integer counts,
  finite numeric inputs, matrix dimensions and tuning ranges.
* Return p, b.i and abs.post.median consistently, allowing both public S2M
  results to be passed to S2MVarSelectionV1. Use all returned sampler draws.
* Accept the documented prior aliases rr, hs and hs+.
* Prevent predictor names from changing the fitted response; preserve supplied
  predictor names on the returned coefficient draws.
* Correct the Li--Pati citation, posterior-draw orientation and tuning-plot
  description. Document tie conventions, stopping boundaries and empty-model
  limits. Replace the vignette's unsupported automatic-optimality claims with
  the actual plotting and selection workflow.

## Selection options

* Add `TwoMeansBeta()` for the basic 2-M rule and `SpikeSlabSelection()` for
  exact-zero spike-and-slab inclusion summaries, the median probability model
  and the most frequent visited complete model. Document ties and boundaries.
* Allow `Sequential2MeansBeta(..., zero.draws = "empty")` to assign zero
  signals to all-zero draws. This explicit extension leaves the default
  rejection and all nonzero-draw behavior unchanged.

* Add DSS with the adaptive-L1 surrogate, an explicit penalty and convex
  optimality diagnostics. It supports empty selections, matching prediction
  designs and correlated predictors without another runtime dependency.
  Document zero-mean coefficients, nonunique solutions and numerical limits.

* Add SAVS post-processing of supplied coefficient estimates. This separate
  selector can return no variables. Document predictor-unit dependence,
  exact-zero/equality handling and extreme-value rounding conventions.
* Add optional exact one-dimensional two-cluster fitting and recursive-pair
  refinement. Keep randomized kmeans and legacy small-set stopping as defaults.
* Add S2MSensitivity with modal-count frequencies and selected-set overlap;
  it does not choose an optimal threshold.
* Add reproducible targeted checks for null/weak signals, correlated designs
  and paired changes of predictor units. The small pilot is exploratory and
  does not establish false-discovery control or general selection accuracy.
