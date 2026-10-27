# Internal validation shared by the paper's posterior-draw summaries.
checkPaperDraws <- function(Beta, min.cols = 1L) {
  checkNumericMatrix(Beta, "Beta", min.cols = min.cols)
  nm <- colnames(Beta)
  if (!is.null(nm) && (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm))) {
    stop("Beta column names must be nonempty, nonmissing and unique.", call. = FALSE)
  }
}

#' Select variables by the paper's two-means rule
#'
#' Cluster each draw's absolute coefficients into two groups. The smaller
#' group size estimates the number of signals for that draw. Select the modal
#' number of variables using posterior medians of absolute coefficients.
#' @param Beta Finite numeric matrix with draws in rows and at least two
#'   coefficients in columns. Exclude any intercept. Column names, if present,
#'   must be nonempty, nonmissing and unique.
#' @param clustering Either `"exact"` (default) for the global one-dimensional
#'   two-means optimum or `"kmeans"` for [stats::kmeans()] with two centers.
#' @return A vsusp_selection list containing `p`, the selected count `H`, per-draw `signal.count`,
#'   `abs.post.median`, selected column indices `selected`, a count table
#'   `frequencies`, the backend `clustering`, and indices `zero.draws` of draws
#'   handled by the zero-model boundary convention.
#' @details This is the 2-M rule in Section 2.1 of Li and Pati (2017).
#'   It uses the smaller cluster size, regardless of which cluster has the
#'   higher mean. It is distinct from the recursive S2M noise-cluster rule.
#'   An equal-size split has that common size. A tied modal count favors the
#'   smaller count. Equal posterior medians favor larger column indices,
#'   matching [S2MVarSelection()]. The exact backend favors the earliest cut
#'   in sorted magnitudes when objective values tie. The kmeans backend uses
#'   R's initialization and RNG; set a seed for reproducibility.
#'
#'   Two-means is undefined for constant magnitudes. As an explicit extension,
#'   an all-zero draw has signal count zero and is not clustered. A constant
#'   nonzero magnitude draw is rejected before consuming clustering RNG.
#'   Neither convention is a published clustering result. `H = 0` returns an
#'   empty selection. This function summarizes supplied draws; it fits no model.
#'   Sparse spike-and-slab draws can give zero medians for every coefficient.
#'   Median ties may then select a coefficient never included in any draw.
#'   Use [SpikeSlabSelection()] for summaries based on inclusion frequencies.
#' @export
#' @references Li, H. and Pati, D. (2017). Variable selection using shrinkage
#'   priors. Computational Statistics & Data Analysis, 107, 107--119.
#'   \doi{10.1016/j.csda.2016.10.008}.
#' @examples
#' Beta <- rbind(c(0, .1, 10, 11, 12), c(0, -.1, -10, -11, -12))
#' TwoMeansBeta(Beta)
#' TwoMeansBeta(matrix(0, 2, 3))
TwoMeansBeta <- function(Beta, clustering = c("exact", "kmeans")) {
  checkPaperDraws(Beta, min.cols = 2L)
  clustering <- matchChoice(clustering, c("exact", "kmeans"), "clustering")
  magnitudes <- abs(Beta)
  constant <- apply(magnitudes, 1L, function(x) all(x == x[1L]))
  zero <- rowSums(magnitudes != 0) == 0L
  if (any(constant & !zero)) {
    stop("Nonzero constant-magnitude draws have no defined two-means split.", call. = FALSE)
  }
  counts <- vapply(seq_len(nrow(Beta)), function(i) {
    if (zero[i]) return(0L)
    split <- splitS2M(Beta[i, ], clustering)
    as.integer(min(tabulate(split$cluster, nbins = 2L)))
  }, integer(1))
  tab <- table(counts)
  frequencies <- data.frame(H = as.integer(names(tab)), draws = as.integer(tab),
                            proportion = as.numeric(tab) / nrow(Beta))
  H <- frequencies$H[which.max(frequencies$draws)]
  medians <- apply(magnitudes, 2L, stats::median)
  structure(list(p = ncol(Beta), H = H, signal.count = counts,
       abs.post.median = medians, selected = rankSelection(medians, H),
       frequencies = frequencies, clustering = clustering,
       zero.draws = which(zero)), class = "vsusp_selection")
}
