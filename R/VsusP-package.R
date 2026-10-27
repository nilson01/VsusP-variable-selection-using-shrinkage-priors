#' Variable selection from Bayesian shrinkage priors
#'
#' Select variables from posterior coefficient draws, or fit Gaussian shrinkage
#' regression and summarize its posterior. Selection and posterior fitting are
#' separate operations.
#' @section Gaussian regression:
#' [FitRidge()], [FitHorseshoe()] and [FitSpikeSlab()] retain separate iteration,
#' chain and predictor axes. [PosteriorDraws()] pools chains for selection.
#' Use summary and predict for posterior estimates and pointwise intervals.
#' [ChainDiagnostics()], [PlotTrace()] and [PlotAutocorrelation()] keep chain
#' identities. Successful sampling alone does not establish convergence.
#' @section Selection:
#' [Sequential2MeansBeta()] and [TwoMeansBeta()] summarize supplied draws;
#' [S2MSensitivity()] compares thresholds. [SpikeSlabSelection()] counts exact
#' nonzero coefficients. [SAVS()] and [DSS()] construct sparse coefficient
#' summaries and can return no selected variables. Clustering counts are not
#' inclusion probabilities or false-discovery guarantees.
#' @section Fixed-posterior checks:
#' The following plots use five observations and one predictor. Four chains
#' retain 3000 iterations each after 1000 burn-in iterations. HS denotes
#' horseshoe and SS denotes spike-and-slab. Sampled posterior summaries are
#' compared with independent numerical integration. The audit tolerance is
#' the larger of 0.01 and six estimated Monte Carlo standard errors, computed
#' using 100-draw batches. The shaded range is an audit criterion, not a
#' confidence interval. These checks do not establish performance across
#' repeated datasets.
#' \if{html}{\figure{posterior-agreement.png}{Sampled-minus-integrated posterior summaries divided by the audit tolerance.}}
#' \if{latex}{\figure{posterior-agreement.pdf}{options: width=5.5in}}
#'
#' The spike-and-slab model has two possible inclusion states. Points compare
#' their integrated and sampled probabilities; error bars are plus or minus
#' 1.96 estimated Monte Carlo standard errors. Traces display the first 300
#' retained coefficient draws per chain. All retained draws enter the summaries.
#' \if{html}{\figure{joint-models.png}{Spike-and-slab model probabilities and separate coefficient traces.}}
#' \if{latex}{\figure{joint-models.pdf}{options: width=5.8in}}
#' @seealso [FitGaussian()], [PlotCoefficients()], [PlotInclusion()],
#'   \code{vignette("VsusPVignette")}
#' @references Li, H. and Pati, D. (2017). Variable selection using shrinkage
#'   priors. Computational Statistics & Data Analysis, 107, 107--119.
#'   \doi{10.1016/j.csda.2016.10.008}.
#' @name VsusP
#' @aliases VsusP-package
"_PACKAGE"
