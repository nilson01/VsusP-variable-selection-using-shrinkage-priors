test_that("2-M uses smaller cluster size even when the lower group is smaller", {
  B <- rbind(c(0, .1, 10, 11, 12), c(0, -.1, -10, -11, -12))
  fit <- TwoMeansBeta(B)
  # The optimal split has two small and three large entries. Section 2.1
  # counts the smaller group (2); counting the upper group would give 3.
  expect_identical(fit$signal.count, c(2L, 2L))
  expect_identical(fit$H, 2L)
  expect_identical(fit$selected, c(5L, 4L))
  expect_equal(fit$abs.post.median, c(0, .1, 10, 11, 12))
  expect_equal(fit$frequencies$proportion, 1)
  expect_identical(fit$zero.draws, integer(0))
})

test_that("2-M resolves equal sizes, count modes and median ties explicitly", {
  B <- rbind(c(0, 0, 0, 0, 10), c(0, 0, 10, 10, 10))
  fit <- TwoMeansBeta(B)
  expect_identical(fit$signal.count, c(1L, 2L))
  expect_identical(fit$H, 1L)
  expect_equal(fit$frequencies$proportion, c(.5, .5))
  expect_identical(fit$selected, 5L)
  equal <- TwoMeansBeta(matrix(c(0, 0, 10, 10), 1))
  expect_identical(equal$H, 2L)
  expect_identical(equal$selected, c(4L, 3L))
  # Exact clustering's tied loss uses its earliest sorted cut, no RNG draw.
  set.seed(91)
  saved <- .Random.seed
  expect_identical(TwoMeansBeta(matrix(c(0, 1, 2), 1))$H, 1L)
  expect_identical(.Random.seed, saved)
})

test_that("2-M has a zero-model boundary and rejects nonzero constant magnitudes", {
  fit <- TwoMeansBeta(matrix(0, 3, 4))
  expect_identical(fit$H, 0L)
  expect_identical(fit$selected, integer(0))
  expect_identical(fit$signal.count, c(0L, 0L, 0L))
  expect_identical(fit$zero.draws, 1:3)
  mixed <- TwoMeansBeta(rbind(c(0, 0, 0), c(0, 0, 10)))
  expect_identical(mixed$H, 0L)
  set.seed(92)
  saved <- .Random.seed
  expect_error(TwoMeansBeta(rbind(c(0, 1, 10), c(-2, 2, -2)), "kmeans"),
               "constant-magnitude")
  expect_identical(.Random.seed, saved)
  expect_error(TwoMeansBeta(matrix(2, 1, 3)), "constant-magnitude")
  expect_error(TwoMeansBeta(matrix(0, 1, 1)), "two columns")
})

test_that("2-M exposes the kmeans backend and preserves coefficient names", {
  B <- rbind(c(a = 0, b = .1, c = 10, d = 11, e = 12),
             c(a = 0, b = -.1, c = -10, d = -11, e = -12))
  set.seed(93)
  fit <- TwoMeansBeta(B, "kmeans")
  expect_identical(fit$clustering, "kmeans")
  expect_identical(names(fit$abs.post.median), colnames(B))
  expect_identical(fit$signal.count, c(2L, 2L))
  expect_identical(TwoMeansBeta(matrix(c(-2, 5), 1), "kmeans")$H, 1L)
  expect_error(TwoMeansBeta(B, "unknown"), "clustering")
})

test_that("spike-slab MPM includes equality and HPPM counts joint patterns", {
  B <- rbind(c(a = 1, b = 0, c = 0), c(2, 0, 0),
             c(0, -1, 0), c(0, -2, 0))
  fit <- SpikeSlabSelection(B)
  # MPM {a,b} is never visited. Both visited models occur twice; the
  # empirical HPPM selects the first one and reports the full tie.
  expect_equal(fit$inclusion.probability, c(a = .5, b = .5, c = 0))
  expect_identical(fit$MPM, c(1L, 2L))
  expect_identical(fit$HPPM, 1L)
  expect_identical(fit$HPPM.ties, c(1L, 2L))
  expect_identical(fit$model.frequencies$draws, c(2L, 2L))
  expect_identical(fit$model.frequencies$first.draw, c(1L, 3L))
  expect_equal(fit$model.frequencies$proportion, c(.5, .5))
  expect_equal(unname(fit$visited.models), rbind(c(TRUE, FALSE, FALSE),
                                              c(FALSE, TRUE, FALSE)))
})

test_that("spike-slab summary distinguishes the most frequent joint model", {
  B <- rbind(c(0, 0, 1), c(1, 1, 0), c(1, 1, 0),
             c(1, 0, 1), c(0, 1, 1))
  fit <- SpikeSlabSelection(B)
  expect_equal(fit$inclusion.probability, rep(.6, 3))
  expect_identical(fit$MPM, 1:3)
  expect_identical(fit$HPPM, c(1L, 2L))
  expect_identical(fit$HPPM.ties, 2L)
  expect_identical(fit$model.frequencies$draws, c(1L, 2L, 1L, 1L))
})

test_that("spike-slab exact zero contract includes the empty model", {
  empty <- SpikeSlabSelection(matrix(0, 3, 2))
  expect_identical(empty$MPM, integer(0))
  expect_identical(empty$HPPM, integer(0))
  expect_equal(empty$inclusion.probability, c(0, 0))
  expect_equal(empty$visited.models, matrix(FALSE, 1, 2))
  expect_identical(empty$model.frequencies$draws, 3L)
  tiny <- SpikeSlabSelection(matrix(c(0, 1e-300, -1e-300), 1))
  expect_identical(tiny$MPM, c(2L, 3L))
  expect_identical(tiny$HPPM, c(2L, 3L))
  one <- SpikeSlabSelection(matrix(c(0, -1, 0, -2), 4, 1))
  expect_equal(one$inclusion.probability, .5)
  expect_identical(one$MPM, 1L)
  expect_identical(one$HPPM, integer(0))
  expect_equal(dim(one$visited.models), c(2L, 1L))
})

test_that("paper selection workflows reject malformed shapes, values and names", {
  for (fun in list(TwoMeansBeta, SpikeSlabSelection)) {
    for (bad in list(c(0, 1), data.frame(a = 0, b = 1),
                     matrix(numeric(0), 0, 2), matrix(numeric(0), 2, 0),
                     matrix(c(0, NA), 1), matrix(c(0, Inf), 1),
                     matrix(c(0, NaN), 1), matrix(c("0", "1"), 1),
                     matrix(c(0i, 1i), 1), matrix(c(FALSE, TRUE), 1))) {
      expect_error(fun(bad), "Beta")
    }
    for (nm in list(c("a", "a"), c("", "b"), c(NA_character_, "b"))) {
      B <- matrix(c(0, 1), 1, dimnames = list(NULL, nm))
      expect_error(fun(B), "column names")
    }
  }
})

test_that("S2M handles exact-zero spike draws only by explicit convention", {
  B <- rbind(c(0, 0, 0), c(0, 0, 0), c(0, .1, 10))
  expect_error(Sequential2MeansBeta(B, 1, 2, 2), "distinct")
  out <- Sequential2MeansBeta(B, 1, 2, 2, clustering = "exact",
                              zero.draws = "empty")
  expect_equal(out$H.b.i, c(0, 0))
  expect_identical(out$zero.draws, 1:2)
  expect_identical(S2MVarSelectionV1(out, out$H.b.i[1]), integer(0))
  allzero <- Sequential2MeansBeta(matrix(0, 2, 4), 1, 1, 1,
                                  zero.draws = "empty")
  expect_equal(allzero$H.b.i, 0)
  nonzero <- matrix(c(0, .1, 10), 1)
  standard <- Sequential2MeansBeta(nonzero, 1, 2, 2, clustering = "exact")
  extended <- Sequential2MeansBeta(nonzero, 1, 2, 2, clustering = "exact",
                                   zero.draws = "empty")
  extended$zero.draws <- NULL
  expect_identical(standard, extended)
  expect_error(Sequential2MeansBeta(matrix(1, 2, 3), 1, 1, 1,
                                    zero.draws = "empty"), "distinct")
  expect_error(Sequential2MeansBeta(B, 1, 1, 1, zero.draws = "other"),
               "zero.draws")
})

test_that("median ranking is distinct from spike inclusion ranking", {
  B <- rbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0))
  expect_identical(TwoMeansBeta(B)$H, 1L)
  expect_identical(TwoMeansBeta(B)$selected, 4L)
  expect_identical(SpikeSlabSelection(B)$MPM, integer(0))
  expect_equal(SpikeSlabSelection(B)$inclusion.probability[4], 0)
})
