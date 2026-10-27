#' Summarize exact spike-and-slab posterior draws
#'
#' Estimate inclusion probabilities and return the median probability model
#' (MPM) and the most frequent complete inclusion pattern visited by the draws
#' (the empirical HPPM).
#' @param Beta Finite numeric matrix with draws in rows and coefficients in
#'   columns. Exclude any intercept. Excluded coefficients must be exactly
#'   zero; all other values represent inclusion. Column names, if present,
#'   must be nonempty, nonmissing and unique.
#' @return A vsusp_spike list containing `p`, `inclusion.probability`, selected column indices
#'   `MPM` and `HPPM`, logical matrix `visited.models` (one distinct model per
#'   row), `model.frequencies` (model row index, draw count, proportion, and
#'   first draw), and `HPPM.ties` (all model row indices tied for largest count).
#' @details The MPM includes variables with inclusion frequency greater than
#'   or equal to 0.5, as in Li and Pati (2017). The empirical HPPM counts full
#'   joint inclusion patterns, rather than combining marginal probabilities.
#'   It estimates the highest-probability visited model from the supplied
#'   sample; it does not establish the global posterior mode. Model rows are
#'   ordered by first occurrence. If model frequencies tie, `HPPM` uses the
#'   first observed tied model; `HPPM.ties` retains all tied row indices.
#'
#'   Exact zeros are a data contract, not a numerical tolerance. Draws from a
#'   continuous shrinkage prior are inappropriate for this summary: small
#'   nonzero values count as included. An all-zero draw is the empty model,
#'   and an entirely zero matrix returns empty MPM and HPPM selections.
#'   This function fits no model and derives no probabilities from clustering.
#' @export
#' @references Li, H. and Pati, D. (2017). Variable selection using shrinkage
#'   priors. Computational Statistics & Data Analysis, 107, 107--119.
#'   \doi{10.1016/j.csda.2016.10.008}.
#' @examples
#' Beta <- rbind(c(1, 0, 0), c(2, 0, 0), c(0, -1, 0), c(0, -2, 0))
#' SpikeSlabSelection(Beta)
SpikeSlabSelection <- function(Beta) {
  checkPaperDraws(Beta)
  included <- Beta != 0
  keys <- apply(included, 1L, function(x) paste0(as.integer(x), collapse = ""))
  first <- which(!duplicated(keys))
  ids <- match(keys, keys[first])
  counts <- tabulate(ids, nbins = length(first))
  visited <- included[first, , drop = FALSE]
  rownames(visited) <- NULL
  frequencies <- data.frame(model = seq_along(first), draws = counts,
                            proportion = counts / nrow(Beta), first.draw = first)
  ties <- which(counts == max(counts))
  probability <- colMeans(included)
  structure(list(p = ncol(Beta), inclusion.probability = probability,
       MPM = unname(which(probability >= .5)),
       HPPM = unname(which(visited[ties[1L], ])),
       visited.models = visited, model.frequencies = frequencies,
       HPPM.ties = ties), class = "vsusp_spike")
}
