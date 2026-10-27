# Two-cluster optimum over distinct adjacent cuts of sorted magnitudes.
# Stable prefix/suffix sums of squared deviations avoid subtracting large sums.
exactTwoMeans <- function(x) {
  o <- order(x)
  z <- x[o]
  if (length(z) < 2L || z[1L] == z[length(z)]) {
    stop("Exact clustering requires two distinct magnitudes.")
  }
  scaled <- z / max(z)
  if (any(diff(z) > 0 & diff(scaled) == 0)) {
    stop("Distinct magnitudes collapsed during rescaling; exact clustering is unavailable at this numeric range.")
  }
  gaps <- diff(scaled)
  if (any(gaps > 0 & gaps^2 == 0)) {
    stop("Magnitude separation is too small for exact clustering at this numeric range.")
  }
  n <- length(z)
  prefix <- suffix <- numeric(n)
  mu <- scaled[1L]
  for (i in 2:n) {
    delta <- scaled[i] - mu
    prefix[i] <- prefix[i - 1L] + delta^2 * (i - 1L) / i
    mu <- mu + delta / i
  }
  mu <- scaled[n]
  for (i in seq.int(n - 1L, 1L)) {
    count <- n - i + 1L
    delta <- scaled[i] - mu
    suffix[i] <- suffix[i + 1L] + delta^2 * (count - 1L) / count
    mu <- mu + delta / count
  }
  candidates <- which(diff(z) > 0)
  loss <- prefix[candidates] + suffix[candidates + 1L]
  cut <- candidates[which.min(loss)]
  cluster <- integer(n)
  cluster[o] <- c(rep.int(1L, cut), rep.int(2L, n - cut))
  list(cluster = cluster,
       centers = c(mean(z[seq_len(cut)]), mean(z[seq.int(cut + 1L, n)])))
}

splitS2M <- function(x, clustering) {
  if (clustering == "exact") return(exactTwoMeans(abs(x)))
  if (length(x) == 2L) {
    o <- order(abs(x))
    cluster <- integer(2L)
    cluster[o] <- 1:2
    return(list(cluster = cluster, centers = abs(x[o])))
  }
  stats::kmeans(abs(x), 2L)
}
