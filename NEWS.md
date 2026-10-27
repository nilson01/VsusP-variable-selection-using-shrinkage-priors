# VsusP 1.1.0

## Gaussian fitting and summaries

* Add `FitRidge()`, `FitHorseshoe()` and `FitSpikeSlab()` with explicit priors,
  intercept draws and separate chain axes. Existing `bayesreg` defaults remain.
* Add coefficient summaries, posterior and sparse prediction, `PosteriorDraws()`
  and `ChainDiagnostics()`. Add coefficient, inclusion, trace and autocorrelation
  plots that return the displayed values invisibly.

## Selection methods

* Add `TwoMeansBeta()` for basic two-means selection and `SpikeSlabSelection()`
  for exact-zero inclusion summaries and visited-model frequencies.
* Add SAVS and the adaptive-L1 DSS approximation as separate methods that can
  return an empty selection. Document scaling, tuning and numerical limits.
* Add exact one-dimensional clustering, recursive-pair refinement and
  `S2MSensitivity()`. Randomized kmeans and legacy small-set stopping remain
  the defaults; sensitivity summaries do not choose an optimal threshold.
* Allow all-zero S2M draws with `zero.draws = "empty"`; default rejection remains.

## Fixes

* Retain the saved S2M noise set at recursive stopping, following Li and Pati
  (2017), Section 2.2.1, <doi:10.1016/j.csda.2016.10.008>. The previous count
  used the smaller cluster from the stopping split and could select too many
  signals. The strict center-gap comparison is unchanged.
* Handle distinct initial pairs and indivisible recursive candidates. Reject
  unsupported constant-magnitude draws; accept explicit `H = 0` selections.
* Validate dimensions, finite values, scalar counts and tuning ranges.
* Return consistent selection summaries and use every returned sampler draw.
* Accept prior aliases `rr`, `hs` and `hs+`; preserve predictor names without
  allowing a predictor name to replace the response.
* Correct citations, draw orientation and tuning-plot descriptions. Document
  tie conventions and empty-model limits.

## Documentation and checks

* Add eight reproducible figures and worked fitting, selection and prediction
  examples. Retain the exploratory pilot's outcomes and limitations.
* Remove tracked build products and duplicate documentation. Run one package
  check per CI job and retain strict test-count and coverage verification.
