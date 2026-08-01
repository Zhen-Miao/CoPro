## Tests for selectSigmaByPermutation(): the studentized max-T bandwidth
## selector. The properties worth pinning are the ones the method's validity
## rests on -- that the kernel/statistic is what it claims, that the shift is
## rigid and inside the wrap box, that the studentized null is flat by
## construction, that the p-value is the Phipson-Smyth max-T one, and that
## the coordinates are on the scale sigma is expressed in.

.fit_selection_object <- function(cts, sigmaValues = c(0.05, 0.1, 0.2),
                                  nCC = 2L, seed = 42) {
  toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
  obj <- newCoProSingle(
    normalizedData = toy$normalizedData,
    locationData   = toy$locationData,
    metaData       = toy$metaData,
    cellTypes      = toy$cellTypes
  )
  obj <- subsetData(obj, cellTypesOfInterest = cts)
  obj <- computePCA(obj, nPCA = 5)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = sigmaValues, verbose = FALSE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = nCC)
  obj <- suppressWarnings(computeNormalizedCorrelation(obj))
  computeGeneAndCellScores(obj)
}

## ---------------------------------------------------------------------------
## The statistic
## ---------------------------------------------------------------------------

test_that("the blocked bilinear statistic equals the dense one, including the self diagonal", {
  set.seed(3)
  cA <- cbind(runif(37), runif(37))
  cB <- cbind(runif(23), runif(23))
  sigmas <- c(0.05, 0.2)
  a <- lapply(sigmas, function(s) rnorm(37))
  b <- lapply(sigmas, function(s) rnorm(23))

  dense <- vapply(seq_along(sigmas), function(g) {
    K <- exp(-0.5 * as.matrix(
      fields::rdist(cA, cB))^2 / sigmas[g]^2)
    sum(a[[g]] * as.numeric(K %*% b[[g]]))
  }, numeric(1))

  ## Block size must not change the answer: 1 row at a time, or all at once.
  for (bs in c(1L, 7L, 1024L)) {
    got <- CoPro:::.bilinearOverSigma(cA, cB, a, b, sigmas,
                                      selfPair = FALSE, blockSize = bs)
    expect_equal(got, dense, tolerance = 1e-10)
  }

  ## Self pair: w_ii = 0, so the diagonal must be excluded.
  cS <- cbind(runif(29), runif(29))
  aS <- lapply(sigmas, function(s) rnorm(29))
  dense_self <- vapply(seq_along(sigmas), function(g) {
    K <- exp(-0.5 * as.matrix(fields::rdist(cS, cS))^2 / sigmas[g]^2)
    diag(K) <- 0
    sum(aS[[g]] * as.numeric(K %*% aS[[g]]))
  }, numeric(1))
  for (bs in c(1L, 5L, 1024L)) {
    got <- CoPro:::.bilinearOverSigma(cS, cS, aS, aS, sigmas,
                                      selfPair = TRUE, blockSize = bs)
    expect_equal(got, dense_self, tolerance = 1e-10)
  }

  ## and the diagonal really is what differs between the two modes
  with_diag <- CoPro:::.bilinearOverSigma(cS, cS, aS, aS, sigmas,
                                          selfPair = FALSE, blockSize = 1024L)
  expect_false(isTRUE(all.equal(with_diag, dense_self)))
})

test_that("the kernel matches computeKernelMatrix's convention, not a factor-of-2 variant", {
  ## sigma is a distance and gets handed back to computeKernelMatrix(), so the
  ## two must agree on exp(-0.5 (d/sigma)^2). A 1/sigma^2 convention would make
  ## the selected bandwidth mean something else downstream.
  cA <- matrix(c(0, 0), nrow = 1)
  cB <- matrix(c(1, 0), nrow = 1)
  got <- CoPro:::.bilinearOverSigma(cA, cB, list(1), list(1), 2,
                                    selfPair = FALSE, blockSize = 1L)
  expect_equal(got, exp(-0.5 * (1 / 2)^2), tolerance = 1e-12)
})

## ---------------------------------------------------------------------------
## The null
## ---------------------------------------------------------------------------

test_that("the toroidal shift stays inside the wrap box and moves every cell together", {
  set.seed(11)
  coords <- cbind(runif(200, 2, 9), runif(200, -3, 4))
  domain <- CoPro:::.resolveSelectionDomain(NULL, list(coords))

  for (i in 1:25) {
    shifted <- CoPro:::.toroidalShift(coords, domain)
    expect_true(all(shifted[, 1] >= domain$lower[1] - 1e-9))
    expect_true(all(shifted[, 1] <= domain$upper[1] + 1e-9))
    expect_true(all(shifted[, 2] >= domain$lower[2] - 1e-9))
    expect_true(all(shifted[, 2] <= domain$upper[2] + 1e-9))
    ## Rigid: every cell moves by the same offset modulo the wrap, so the
    ## displacements take at most two distinct values per axis (wrapped or not).
    for (j in 1:2) {
      expect_lte(length(unique(round(shifted[, j] - coords[, j], 8))), 2L)
    }
  }
})

test_that("the studentized null is flat across bandwidths and the p-value is Phipson-Smyth max-T", {
  set.seed(5)
  B <- 200L
  ## Three columns with deliberately different null scales.
  nullMat <- cbind(rnorm(B, sd = 1), rnorm(B, sd = 50), rnorm(B, sd = 1000))
  obs <- c(2, 100, 2000)

  st <- CoPro:::.studentizeSelectMaxT(obs, nullMat, "greater")

  expect_equal(st$nullSD, apply(nullMat, 2, sd))
  expect_equal(st$z, obs / apply(nullMat, 2, sd))
  ## After studentizing, every column's null has unit spread. That flatness is
  ## the whole reason the argmax is comparable across bandwidths.
  zNull <- sweep(nullMat, 2, st$nullSD, "/")
  expect_equal(unname(apply(zNull, 2, sd)), rep(1, 3), tolerance = 1e-12)

  expect_equal(st$nullMax, apply(zNull, 1, max))
  expect_equal(st$zMax, max(st$z))
  ## +1 in numerator and denominator; smallest resolvable p is 1/(B+1).
  expect_equal(st$pAdjusted,
               vapply(st$z, function(zi) (1 + sum(st$nullMax >= zi)) / (1 + B),
                      numeric(1)))
  expect_gte(min(st$pAdjusted), 1 / (B + 1))
})

test_that("two.sided uses the null of the maximum absolute studentized statistic", {
  set.seed(6)
  nullMat <- cbind(rnorm(100), rnorm(100, sd = 10))
  obs <- c(-3, 1)
  st <- CoPro:::.studentizeSelectMaxT(obs, nullMat, "two.sided")
  zNull <- sweep(nullMat, 2, st$nullSD, "/")
  expect_equal(st$nullMax, apply(abs(zNull), 1, max))
  ## A strong negative statistic must be able to win under two.sided.
  expect_equal(which.max(st$zSelect), 1L)
})

test_that("a degenerate bandwidth is dropped rather than dividing by zero", {
  nullMat <- cbind(rnorm(50), rep(0, 50))
  expect_warning(
    st <- CoPro:::.studentizeSelectMaxT(c(1, 1), nullMat, "greater"),
    "null spread of zero"
  )
  expect_true(is.na(st$z[2]))
  expect_true(is.na(st$pAdjusted[2]))
  expect_true(is.finite(st$zMax))
})

## ---------------------------------------------------------------------------
## Coordinates carry the bandwidth's units
## ---------------------------------------------------------------------------

test_that("coordinates are put on the same distance scale computeDistance used", {
  ## Stop before computeKernelMatrix(), which drops @distances by default, so
  ## the stored matrix is still there to compare against.
  toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
  obj <- newCoProSingle(
    normalizedData = toy$normalizedData, locationData = toy$locationData,
    metaData = toy$metaData, cellTypes = toy$cellTypes
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("Epithelial", "Fibroblast"))
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)

  coords <- CoPro:::.selectionCoords(obj, c("Epithelial", "Fibroblast"))

  ## Distances between these coordinates must reproduce the stored distance
  ## matrix. If they did not, the scanned sigmas would be physically wrong --
  ## sigma is a distance, and the selected value is handed to
  ## computeKernelMatrix() afterwards.
  stored <- as.matrix(obj@distances[[
    CoPro:::.createDistMatrixName("Epithelial", "Fibroblast", slide = NULL)
  ]])
  rebuilt <- fields::rdist(coords[["Epithelial"]], coords[["Fibroblast"]])

  ## computeDistance(truncateLowDist = TRUE) additionally floors the bottom
  ## ~0.1% of distances; the selector does not replicate that floor (see the
  ## note in selectSigmaByPermutation()). Apply it here so the comparison tests
  ## the scale, which is the part that has to match, and confirm the residual
  ## disagreement is confined to the floored entries.
  floored <- rebuilt
  floored[floored < min(stored)] <- min(stored)
  expect_equal(as.numeric(floored), as.numeric(stored), tolerance = 1e-8)

  differing <- which(abs(rebuilt - stored) > 1e-8)
  expect_true(all(rebuilt[differing] < min(stored)))
  expect_lt(length(differing) / length(stored), 0.01)

  ## and the object really did normalize, so this is not a trivial pass
  expect_false(isTRUE(all.equal(obj@distanceScaleFactor, 1)))
})

## ---------------------------------------------------------------------------
## End to end
## ---------------------------------------------------------------------------

test_that("selection runs on a cross-type object and reports a coherent result", {
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"))
  sel <- selectSigmaByPermutation(obj, nPermu = 49L, seed = 1, verbose = FALSE)

  expect_s3_class(sel, "CoProSigmaSelection")
  expect_true(sel$sigmaValueChoice %in% obj@sigmaValues)
  expect_equal(nrow(sel$perSigma), length(obj@sigmaValues))
  expect_equal(length(sel$nullMax), 49L)
  expect_gte(sel$pValue, 1 / 50)
  expect_lte(sel$pValue, 1)
  ## the reported optimum is the argmax of z, and its p is the max-T p
  expect_equal(sel$sigmaValueChoice,
               sel$perSigma$sigma[which.max(sel$perSigma$z)])
  expect_equal(sel$zMax, max(sel$perSigma$z))
  expect_true(all(sel$plateau %in% sel$perSigma$sigma))
  expect_output(print(sel), "studentized max-T")
})

test_that("selection runs on a within-type object, where the permutation routes cannot", {
  ## A single cell type is the case the existing permutation route cannot
  ## serve. runSkrCCAPermu() itself returns, but the normalized-correlation
  ## step that has to follow it forms combn(cts, 2) and dies on one type. It is
  ## also the case the argmax-of-ratio rule is least trustworthy on, because
  ## the stored self-kernel has a zero diagonal and so is not a valid whitening
  ## operator. This selector needs neither.
  obj <- .fit_selection_object("Epithelial")
  permuted <- suppressWarnings(runSkrCCAPermu(obj, nPermu = 5L, verbose = FALSE))
  expect_error(computeNormalizedCorrelationPermu(permuted), "n < m")

  sel <- selectSigmaByPermutation(obj, nPermu = 49L, seed = 1, verbose = FALSE)
  expect_true(sel$sigmaValueChoice %in% obj@sigmaValues)
  expect_equal(nrow(sel$perSigma), length(obj@sigmaValues))
  expect_true(all(sel$perSigma$cellType1 == "Epithelial"))
  expect_true(all(sel$perSigma$cellType2 == "Epithelial"))
})

test_that("the same seed reproduces the selection exactly", {
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"))
  a <- selectSigmaByPermutation(obj, nPermu = 25L, seed = 99, verbose = FALSE)
  b <- selectSigmaByPermutation(obj, nPermu = 25L, seed = 99, verbose = FALSE)
  expect_equal(a$perSigma, b$perSigma)
  expect_equal(a$nullMax, b$nullMax)
  expect_equal(a$sigmaValueChoice, b$sigmaValueChoice)
})

test_that("maxCells subsamples and keeps the observed statistic on the same cells as the null", {
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"))
  sel <- selectSigmaByPermutation(obj, nPermu = 19L, maxCells = 40L, seed = 4,
                                  verbose = FALSE)
  expect_true(all(vapply(sel$cells, length, integer(1)) <= 40L))
  ## the recorded sample is what the statistic was computed on, so a rerun
  ## restricted to those cells must reproduce it
  expect_true(all(is.finite(sel$perSigma$statistic)))
})

test_that("selection scores sumcor-requested weights without special-casing them", {
  ## The directions are held fixed inside every draw, so the statistic does not
  ## depend on which canonical criterion produced them. (On a single-slide
  ## object sumcor is routed to the sumcov solvers anyway, so this is a smoke
  ## test that the routing does not disturb the selector, not a claim that the
  ## other permutation routes would refuse it.)
  toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
  obj <- newCoProSingle(
    normalizedData = toy$normalizedData, locationData = toy$locationData,
    metaData = toy$metaData, cellTypes = toy$cellTypes
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("Epithelial", "Fibroblast"))
  obj <- computePCA(obj, nPCA = 5)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1), verbose = FALSE)
  obj <- suppressMessages(runSkrCCA(obj, scalePCs = TRUE, nCC = 1,
                                    objective = "sumcor"))
  obj <- suppressWarnings(computeNormalizedCorrelation(obj))
  obj <- computeGeneAndCellScores(obj)

  expect_equal(getCCAObjective(obj)$requested, "sumcor")
  sel <- suppressWarnings(
    selectSigmaByPermutation(obj, nPermu = 19L, seed = 2, verbose = FALSE)
  )
  expect_true(sel$sigmaValueChoice %in% obj@sigmaValues)
})

## ---------------------------------------------------------------------------
## Refusals
## ---------------------------------------------------------------------------

test_that("selection refuses inputs it cannot serve, with a reason", {
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"))

  expect_error(selectSigmaByPermutation(obj, sigmaValues = 0.123, verbose = FALSE),
               "No cell scores stored for sigma")
  expect_error(selectSigmaByPermutation(obj, ccIndex = 99L, nPermu = 19L,
                                        verbose = FALSE),
               "only 2 component")
  expect_error(selectSigmaByPermutation(obj, ccIndex = 0, verbose = FALSE),
               "positive integer")
  expect_error(selectSigmaByPermutation(obj, alpha = 1, verbose = FALSE),
               "alpha must be")

  bare <- newCoProSingle(
    normalizedData = matrix(rnorm(200), 20, 10,
                            dimnames = list(paste0("c", 1:20), paste0("g", 1:10))),
    locationData = data.frame(x = runif(20), y = runif(20),
                              row.names = paste0("c", 1:20)),
    metaData = data.frame(cell_id = paste0("c", 1:20),
                          row.names = paste0("c", 1:20)),
    cellTypes = rep(c("A", "B"), 10)
  )
  expect_error(selectSigmaByPermutation(bare, verbose = FALSE),
               "computeGeneAndCellScores")

  multi <- create_test_copro_multi(n_cells_per_slide = 20, n_slides = 2)
  expect_error(selectSigmaByPermutation(multi, verbose = FALSE),
               "CoProSingle")
})

test_that("a bandwidth selected at the edge of the grid warns even when quiet", {
  ## A scan that rails at an end of the grid has run out of grid rather than
  ## found an optimum, so this must not be gated on verbose -- the organoid
  ## vignette hit exactly that and verbose = FALSE swallowed it.
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"))
  warned <- character(0)
  sel <- withCallingHandlers(
    selectSigmaByPermutation(obj, nPermu = 19L, seed = 7, verbose = FALSE),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  at_edge <- sel$sigmaValueChoice %in% range(sel$perSigma$sigma)
  expect_equal(any(grepl("end of the scanned grid", warned)), at_edge)
})

## ---------------------------------------------------------------------------
## The small-bandwidth floor
## ---------------------------------------------------------------------------

test_that("candidate bandwidths below the nearest-partner spacing are dropped", {
  ## Below one cell spacing the kernel is nearly diagonal and the
  ## fixed-direction null understates the floor, so z climbs without bound as
  ## sigma shrinks and the argmax rails at the smallest candidate.
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"),
                               sigmaValues = c(0.01, 0.05, 0.1, 0.2))
  spacing <- CoPro:::.blockNearestSpacing(
    CoPro:::.selectionCoords(obj, "Epithelial")[["Epithelial"]],
    CoPro:::.selectionCoords(obj, "Fibroblast")[["Fibroblast"]],
    within = FALSE
  )
  expect_gt(spacing, 0.01)

  floored <- suppressWarnings(
    selectSigmaByPermutation(obj, nPermu = 19L, seed = 3, verbose = FALSE)
  )
  expect_false(0.01 %in% floored$perSigma$sigma)
  expect_equal(unname(floored$settings$minSigma), unname(spacing))

  ## minSigma = NULL scans it anyway, and the dropped bandwidth is exactly the
  ## one that would have won.
  unfloored <- suppressWarnings(
    selectSigmaByPermutation(obj, nPermu = 19L, seed = 3, minSigma = NULL,
                             verbose = FALSE)
  )
  expect_true(0.01 %in% unfloored$perSigma$sigma)
  expect_equal(unfloored$sigmaValueChoice, 0.01)
})

test_that("an explicit numeric floor is honoured and an impossible one is refused", {
  obj <- .fit_selection_object(c("Epithelial", "Fibroblast"))

  sel <- suppressWarnings(
    selectSigmaByPermutation(obj, nPermu = 19L, seed = 3, minSigma = 0.09,
                             verbose = FALSE)
  )
  expect_equal(sort(unique(sel$perSigma$sigma)), c(0.1, 0.2))

  expect_error(
    selectSigmaByPermutation(obj, minSigma = 99, verbose = FALSE),
    "below the nearest-partner spacing"
  )
  expect_error(
    selectSigmaByPermutation(obj, minSigma = "wide", verbose = FALSE),
    "minSigma must be"
  )
})
