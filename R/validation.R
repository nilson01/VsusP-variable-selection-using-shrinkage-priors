# Shared input checks for public selection workflows.
checkNumericVector <- function(x, name, nonnegative = FALSE) {
  if (!is.numeric(x) || is.complex(x) || !is.null(dim(x)) ||
      length(x) == 0L || any(!is.finite(x))) {
    stop(name, " must be a nonempty finite numeric vector.", call. = FALSE)
  }
  if (nonnegative && any(x < 0)) {
    stop(name, " must not contain negative values.", call. = FALSE)
  }
}

checkNumericScalar <- function(x, name, nonnegative = FALSE) {
  checkNumericVector(x, name, nonnegative)
  if (length(x) != 1L) stop(name, " must be a single number.", call. = FALSE)
}

checkIntegerScalar <- function(x, name, lower = 0) {
  checkNumericScalar(x, name)
  if (x < lower || x %% 1 != 0) {
    stop(name, " must be an integer greater than or equal to ", lower, ".", call. = FALSE)
  }
}

checkNumericMatrix <- function(x, name, min.cols = 1L) {
  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
      nrow(x) < 1L || ncol(x) < min.cols || any(!is.finite(x))) {
    stop(name, " must be a finite numeric matrix with at least one row and ",
         if (min.cols == 2L) "two columns." else "one column.", call. = FALSE)
  }
}

checkSelectionCount <- function(H, p) {
  checkIntegerScalar(H, "H")
  if (H > p) stop("H must not exceed the number of variables.", call. = FALSE)
}

rankSelection <- function(medians, H) {
  # Preserve the existing rule: larger column indices win equal-median ties.
  rev(order(medians))[seq_len(H)]
}

normalizePrior <- function(prior) {
  aliases <- c(rr = "ridge", ridge = "ridge", lasso = "lasso",
               hs = "horseshoe", horseshoe = "horseshoe",
               "hs+" = "horseshoe+", "horseshoe+" = "horseshoe+")
  if (!is.character(prior) || length(prior) != 1L || is.na(prior) ||
      !prior %in% names(aliases)) {
    stop("prior must be ridge (rr), lasso, horseshoe (hs), or horseshoe+ (hs+).",
         call. = FALSE)
  }
  unname(aliases[prior])
}

summarizeS2M <- function(Beta, b.i, clustering = "kmeans", small.set = "legacy",
                         details = FALSE, zero.draws = "error") {
  zero <- rowSums(Beta != 0) == 0L
  # Reject undefined initial splits before consuming any clustering RNG draws.
  constant <- apply(abs(Beta), 1L, function(x) length(unique(x)) < 2L)
  if (any(constant & !(zero & zero.draws == "empty"))) {
    stop("Each Beta draw must contain at least two distinct absolute coefficients.")
  }
  counts <- lapply(b.i, function(b) {
    vapply(seq_len(nrow(Beta)), function(i) {
      if (zero[i] && zero.draws == "empty") return(0)
      ncol(Beta) - numNoiseCoeff(Beta[i, ], b, clustering, small.set)
    }, numeric(1))
  })
  frequencies <- lapply(counts, function(x) {
    tab <- table(x)
    data.frame(H = as.numeric(names(tab)), draws = as.integer(tab),
               proportion = as.numeric(tab) / length(x))
  })
  H.b.i <- vapply(frequencies, function(f) f$H[which.max(f$draws)], numeric(1))
  result <- list(p = ncol(Beta), b.i = b.i, H.b.i = H.b.i,
                 abs.post.median = apply(abs(Beta), 2L, stats::median))
  if (details) result$frequencies <- frequencies
  if (zero.draws == "empty") result$zero.draws <- which(zero)
  structure(result, class = "vsusp_tuning")
}

matchChoice <- function(value, choices, name) {
  if (identical(value, choices)) return(choices[1L])
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !value %in% choices) {
    stop(name, " must be one of: ", paste(choices, collapse = ", "), ".", call. = FALSE)
  }
  value
}

# Keep dependency calls behind a local seam so tests mock only this package.
fitShrinkageGaussian <- function(formula, data, model, prior, n.samples, burnin) {
  bayesreg::bayesreg(formula, data, model = model, prior = prior,
                    n.samples = n.samples, burnin = burnin)
}
