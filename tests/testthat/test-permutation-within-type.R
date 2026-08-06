## Within-type (single cell type) permutation testing.
##
## Two defects used to make this path unusable, and the second was the
## dangerous one:
##
##   1. computeNormalizedCorrelationPermu() and runSkrCCAPermu_FairSigma()
##      formed combn(cts, 2), which errors with "n < m" on one cell type.
##   2. .getCellPermu()'s "second_only" rule is `ct_index > 1`, so the only
##      cell type was never permuted. Every draw was the identity, the null
##      equalled the observed statistic, and the p-value was 1 by construction.
##
## Fixing (1) alone would have converted a loud crash into a silently
## meaningless p-value, so the tests below pin the permutation being real, not
## merely the call returning.

.within_type_object <- function(cts = "Epithelial", sigmaValues = c(0.05, 0.1),
                                nCC = 1L) {
  toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
  obj <- newCoProSingle(
    normalizedData = toy$normalizedData, locationData = toy$locationData,
    metaData = toy$metaData, cellTypes = toy$cellTypes
  )
  obj <- subsetData(obj, cellTypesOfInterest = cts)
  obj <- computePCA(obj, nPCA = 5)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = sigmaValues, verbose = FALSE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = nCC)
  suppressWarnings(computeNormalizedCorrelation(obj))
}

test_that("the pair list pairs a lone cell type with itself instead of erroring", {
  expect_equal(CoPro:::.permutationPairTypes("A"),
               matrix(c("A", "A"), nrow = 2, ncol = 1))
  expect_error(combn("A", 2), "n < m")   # what the old code did

  expect_equal(CoPro:::.permutationPairTypes(c("A", "B")), combn(c("A", "B"), 2))
  expect_equal(CoPro:::.permutationPairTypes(c("A", "B", "C")),
               combn(c("A", "B", "C"), 2))
})

test_that("a lone cell type is permuted under every permu_which, not held fixed", {
  obj <- .within_type_object()
  for (which_arg in c("second_only", "first_only", "both")) {
    permu <- CoPro:::.getCellPermu(
      object = obj, permu_method = "bin", nPermu = 5L, cts = "Epithelial",
      permu_which = which_arg, num_bins_x = 5, num_bins_y = 5,
      match_quantile = FALSE
    )
    expect_false(
      CoPro:::.isIdentityPermutation(permu[["Epithelial"]]),
      info = paste("permu_which =", which_arg)
    )
  }
})

test_that("permu_which still selects which type moves when there are two", {
  ## The single-type shortcut must not leak into the multi-type semantics.
  obj <- .within_type_object(cts = c("Epithelial", "Fibroblast"))
  cts <- c("Epithelial", "Fibroblast")
  grab <- function(which_arg) {
    permu <- CoPro:::.getCellPermu(
      object = obj, permu_method = "bin", nPermu = 5L, cts = cts,
      permu_which = which_arg, num_bins_x = 5, num_bins_y = 5,
      match_quantile = FALSE
    )
    vapply(cts, function(ct) CoPro:::.isIdentityPermutation(permu[[ct]]),
           logical(1))
  }
  expect_equal(unname(grab("second_only")), c(TRUE, FALSE))
  expect_equal(unname(grab("first_only")), c(FALSE, TRUE))
  expect_equal(unname(grab("both")), c(FALSE, FALSE))
})

test_that("the within-type null is non-degenerate end to end", {
  obj <- .within_type_object()
  permuted <- suppressWarnings(
    runSkrCCAPermu(obj, nPermu = 30L, verbose = FALSE)
  )
  expect_false(CoPro:::.isIdentityPermutation(permuted@cellPermu[["Epithelial"]]))

  ## Under the old identity permutation every draw reproduced the observed
  ## weights exactly, so |cor| was 1 for all of them.
  observed <- obj@skrCCAOut[[paste0("sigma_", obj@sigmaValueChoice)]][["Epithelial"]][, 1]
  cors <- vapply(permuted@skrCCAPermuOut, function(w) {
    v <- w[["Epithelial"]][, 1]
    abs(sum(v * observed) / sqrt(sum(v^2) * sum(observed^2)))
  }, numeric(1))
  expect_lt(mean(cors > 0.999), 0.5)

  permuted <- suppressWarnings(computeNormalizedCorrelationPermu(permuted))
  null_vals <- vapply(permuted@normalizedCorrelationPermu,
                      function(x) x$normalizedCorrelation[1], numeric(1))
  expect_length(null_vals, 30L)
  expect_true(all(is.finite(null_vals)))
  expect_gt(stats::sd(null_vals), 0)
  expect_gt(length(unique(null_vals)), 25L)

  ## Both sides of the self pair are named for the one cell type.
  row1 <- permuted@normalizedCorrelationPermu[[1]]
  expect_equal(row1$cellType1, "Epithelial")
  expect_equal(row1$cellType2, "Epithelial")

  ## A p-value that is not 1 by construction.
  observed_nc <- max(getNormCorr(obj)$normalizedCorrelation)
  p <- (1 + sum(null_vals >= observed_nc)) / (1 + length(null_vals))
  expect_lt(p, 1)
  expect_gte(p, 1 / 31)
})

test_that("the fair-sigma test runs on one cell type", {
  obj <- .within_type_object()
  out <- suppressWarnings(
    runSkrCCAPermu_FairSigma(obj, nPermu = 20L, verbose = FALSE)
  )
  null_vals <- vapply(out@normalizedCorrelationPermu,
                      function(x) x$normalizedCorrelation[1], numeric(1))
  expect_length(null_vals, 20L)
  expect_gt(stats::sd(null_vals), 0)
  expect_equal(out@normalizedCorrelationPermu[[1]]$cellType1, "Epithelial")
  expect_equal(out@normalizedCorrelationPermu[[1]]$cellType2, "Epithelial")
})

test_that("the conditional test reports the within-type null instead of warning about it", {
  obj <- .within_type_object(nCC = 2L)
  ## It used to warn that second_only gives the identity permutation. It no
  ## longer does, so the warning would now be false.
  warned <- character(0)
  out <- withCallingHandlers(
    runSkrCCAPermu_Conditional(obj, nPermu = 20L, verbose = FALSE),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(any(grepl("identity permutation", warned)))

  per_axis <- out@conditionalPermu$per_axis
  expect_true(all(per_axis$p_raw <= 1))
  expect_true(all(per_axis$p_raw >= 0))
  ## A degenerate null would put every p at exactly 1.
  expect_lt(min(per_axis$p_raw), 1)
})
