# Tests for PCA computation and full workflow

test_that("computePCA works for CoProSingle", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  
  # Check pcaGlobal is populated
  expect_true(length(obj@pcaGlobal) == 2)
  expect_true("CellTypeA" %in% names(obj@pcaGlobal))
  expect_true("CellTypeB" %in% names(obj@pcaGlobal))
  
  # Check PCA results structure
  pca_a <- obj@pcaGlobal[["CellTypeA"]]
  expect_true("x" %in% names(pca_a))  # PC scores
  expect_true("rotation" %in% names(pca_a))  # Loadings
  expect_true("sdev" %in% names(pca_a))  # Standard deviations
  
  # Check dimensions
  expect_equal(ncol(pca_a$x), 10)
})

test_that("sparse PCA uses implicit centering and scaling", {
  set.seed(20260725)
  dense <- matrix(rpois(2400, lambda = 0.15), nrow = 60, ncol = 40)
  rownames(dense) <- paste0("cell", seq_len(nrow(dense)))
  colnames(dense) <- paste0("gene", seq_len(ncol(dense)))
  sparse <- methods::as(Matrix::Matrix(dense, sparse = TRUE), "dgCMatrix")

  set.seed(11)
  fit_sparse <- CoPro:::.run_pca_irlba(sparse, 5, center = TRUE, scale. = TRUE)
  set.seed(11)
  fit_dense <- CoPro:::.run_pca_irlba(dense, 5, center = TRUE, scale. = TRUE)

  expect_s4_class(sparse, "dgCMatrix")
  expect_equal(fit_sparse$sdev, fit_dense$sdev, tolerance = 1e-7)
  expect_equal(abs(fit_sparse$rotation), abs(fit_dense$rotation),
               tolerance = 1e-6)
})

test_that("computePCA chooses one feasible rank across imbalanced cell types", {
  dat <- generate_test_data_single(n_cells = 60, n_genes = 30, seed = 101)
  dat$cellTypes <- c(rep("A", 12), rep("B", 24), rep("C", 24))
  obj <- newCoProSingle(
    normalizedData = dat$normalizedData,
    locationData = dat$locationData,
    metaData = dat$metaData,
    cellTypes = dat$cellTypes
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("A", "B", "C"))

  suppressWarnings(
    obj <- computePCA(obj, nPCA = 15, center = TRUE, scale. = TRUE)
  )
  ranks <- vapply(obj@pcaGlobal, function(x) ncol(x$x), integer(1))
  expect_equal(unname(ranks), rep(obj@nPCA, 3))
  expect_equal(obj@nPCA, 11)
})

test_that("computePCA works for single cell type", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA"))
  
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  
  expect_true(length(obj@pcaGlobal) == 1)
  expect_true("CellTypeA" %in% names(obj@pcaGlobal))
})

test_that("computePCA validates nPCA parameter", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA"))
  
  expect_error(
    computePCA(obj, nPCA = 0),
    "positive integer"
  )
  
  expect_error(
    computePCA(obj, nPCA = -5),
    "positive integer"
  )
})

test_that("computePCA works for CoProMulti", {
  obj <- create_test_copro_multi(n_cells_per_slide = 60, n_slides = 2, 
                                  n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  
  # Check pcaGlobal
  expect_true(length(obj@pcaGlobal) == 2)
  
  # Check pcaResults (slide-specific projections)
  expect_true(length(obj@pcaResults) == 2)  # 2 slides
  expect_true("Slide1" %in% names(obj@pcaResults))
  expect_true("Slide2" %in% names(obj@pcaResults))
})

test_that("runSkrCCA works for CoProSingle with two cell types", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE, normalizeDistance = TRUE)
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  
  # Run skrCCA
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100)
  
  # Check results
  expect_true(length(obj@skrCCAOut) > 0)
  expect_equal(obj@nCC, 2)
  
  # Check weight vectors
  sigma_name <- "sigma_0.1"
  expect_true(sigma_name %in% names(obj@skrCCAOut))
  
  W_A <- obj@skrCCAOut[[sigma_name]][["CellTypeA"]]
  W_B <- obj@skrCCAOut[[sigma_name]][["CellTypeB"]]
  
  expect_true(is.matrix(W_A))
  expect_true(is.matrix(W_B))
  expect_equal(ncol(W_A), 2)  # nCC components
  expect_equal(ncol(W_B), 2)
})

test_that("runSkrCCA works for single cell type", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA"))
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE, normalizeDistance = TRUE)
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  
  # Run skrCCA
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100)

  expect_true(length(obj@skrCCAOut) > 0)
  expect_equal(ncol(obj@skrCCAOut[[1]][["CellTypeA"]]), 2)

  # The same exact all-axis path is used inside permutation fits.
  obj <- suppressMessages(computeNormalizedCorrelation(obj))
  obj <- suppressWarnings(suppressMessages(
    runSkrCCAPermu(obj, nPermu = 2, permu_method = "global",
                   maxIter = 100, verbose = FALSE)
  ))
  expect_equal(length(obj@skrCCAPermuOut), 2)
  expect_true(all(vapply(
    obj@skrCCAPermuOut,
    function(x) ncol(x[["CellTypeA"]]) == 2L,
    logical(1)
  )))
})

test_that("runSkrCCA requires kernel matrices", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA"))
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  
  expect_error(
    runSkrCCA(obj, scalePCs = TRUE, nCC = 2),
    "Kernel matrix is empty"
  )
})

test_that("computeNormalizedCorrelation works", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1), verbose = FALSE, normalizeDistance = TRUE)
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100)
  
  # Compute normalized correlation
  obj <- computeNormalizedCorrelation(obj)
  
  expect_true(length(obj@normalizedCorrelation) > 0)
  expect_true(length(obj@sigmaValueChoice) > 0)
  
  # Check structure of correlation results
  corr_df <- obj@normalizedCorrelation[[1]]
  expect_true("sigmaValues" %in% colnames(corr_df))
  expect_true("cellType1" %in% colnames(corr_df))
  expect_true("cellType2" %in% colnames(corr_df))
  expect_true("normalizedCorrelation" %in% colnames(corr_df))
})

test_that("computeGeneAndCellScores works", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE, normalizeDistance = TRUE)
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100)
  
  # Compute gene and cell scores
  obj <- computeGeneAndCellScores(obj)
  
  expect_true(length(obj@cellScores) > 0)
  expect_true(length(obj@geneScores) > 0)
  
  # Check cell scores structure
  cs_name <- "cellScores|sigma0.1|CellTypeA"
  expect_true(cs_name %in% names(obj@cellScores))
  
  # Check gene scores structure
  gs_name <- "geneScores|sigma0.1|CellTypeA"
  expect_true(gs_name %in% names(obj@geneScores))
  
  # Gene scores should have same number of rows as genes
  expect_equal(nrow(obj@geneScores[[gs_name]]), length(obj@geneList))
})

test_that("getCellScores accessor works", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", 
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE, normalizeDistance = TRUE)
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100)
  obj <- computeGeneAndCellScores(obj)
  
  # Test accessor
  scores <- getCellScores(obj, sigma = 0.1, cellType = "CellTypeA", verbose = FALSE)
  
  expect_true(is.matrix(scores))
  expect_equal(ncol(scores), 2)  # nCC components
})

test_that("Full workflow runs without errors", {
  # This test runs the complete CoPro workflow
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  
  # Run complete pipeline
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  # dropDistances = FALSE so we can assert the distances slot was populated
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1),
                             dropDistances = FALSE, verbose = FALSE, normalizeDistance = TRUE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100)
  obj <- computeNormalizedCorrelation(obj)
  obj <- computeGeneAndCellScores(obj)

  # All key slots should be populated
  expect_true(length(obj@pcaGlobal) > 0)
  expect_true(length(obj@distances) > 0)
  expect_true(length(obj@kernelMatrices) > 0)
  expect_true(length(obj@skrCCAOut) > 0)
  expect_true(length(obj@normalizedCorrelation) > 0)
  expect_true(length(obj@cellScores) > 0)
  expect_true(length(obj@geneScores) > 0)
})

test_that("Multi-slide full workflow runs without errors", {
  obj <- create_test_copro_multi(n_cells_per_slide = 60, n_slides = 2, 
                                  n_cell_types = 2, seed = 42)
  
  # Run complete pipeline
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  # dropDistances = FALSE so we can assert the distances slot was populated
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1),
                             dropDistances = FALSE, verbose = FALSE, normalizeDistance = TRUE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100)
  obj <- computeNormalizedCorrelation(obj)
  obj <- computeGeneAndCellScores(obj)

  # All key slots should be populated
  expect_true(length(obj@pcaGlobal) > 0)
  expect_true(length(obj@pcaResults) > 0)
  expect_true(length(obj@distances) > 0)
  expect_true(length(obj@kernelMatrices) > 0)
  expect_true(length(obj@skrCCAOut) > 0)
  expect_true(length(obj@normalizedCorrelation) > 0)
  expect_true(length(obj@cellScores) > 0)
  expect_true(length(obj@geneScores) > 0)
})

# ===========================================================================
# Tests for computeRegressionGeneScores
# ===========================================================================

test_that("computeRegressionGeneScores works for CoProSingle", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE, normalizeDistance = TRUE)
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100)
  obj <- computeGeneAndCellScores(obj)

  # Save original PCA-based gene scores
  gs_orig <- obj@geneScores

  obj <- computeRegressionGeneScores(obj, verbose = FALSE)

  # geneScoresRegression slot should be populated

  expect_true(length(obj@geneScoresRegression) > 0)

  # Same flat-list keys as geneScores
  expect_equal(sort(names(obj@geneScoresRegression)), sort(names(obj@geneScores)))

  # Check dimensions match
  gs_name <- "geneScores|sigma0.1|CellTypeA"
  expect_true(gs_name %in% names(obj@geneScoresRegression))
  expect_equal(nrow(obj@geneScoresRegression[[gs_name]]), length(obj@geneList))
  expect_equal(ncol(obj@geneScoresRegression[[gs_name]]), 2)  # nCC = 2

  # Row names should be gene names
  expect_equal(rownames(obj@geneScoresRegression[[gs_name]]),
               rownames(obj@geneScores[[gs_name]]))

  # Original PCA-based gene scores should be unchanged
  expect_identical(obj@geneScores, gs_orig)
})

test_that("computeRegressionGeneScores works for CoProMulti", {
  obj <- create_test_copro_multi(n_cells_per_slide = 60, n_slides = 2,
                                  n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE, normalizeDistance = TRUE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100)
  obj <- computeGeneAndCellScores(obj)

  obj <- computeRegressionGeneScores(obj, verbose = FALSE)

  expect_true(length(obj@geneScoresRegression) > 0)
  expect_equal(sort(names(obj@geneScoresRegression)), sort(names(obj@geneScores)))

  gs_name <- "geneScores|sigma0.1|CellTypeA"
  expect_equal(nrow(obj@geneScoresRegression[[gs_name]]), length(obj@geneList))
})

test_that("computeRegressionGeneScores errors without cell scores", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)

  expect_error(
    computeRegressionGeneScores(obj),
    "Cell scores not available"
  )
})

test_that("computeRegressionGeneScores respects sigma parameter", {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1), verbose = FALSE, normalizeDistance = TRUE)
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100)
  obj <- computeGeneAndCellScores(obj)

  # Only compute for sigma = 0.1
  obj <- computeRegressionGeneScores(obj, sigma = 0.1, verbose = FALSE)

  # Should have entries for sigma 0.1 only (2 cell types)
  reg_names <- names(obj@geneScoresRegression)
  expect_true(all(grepl("sigma0.1", reg_names)))
  expect_false(any(grepl("sigma0.05", reg_names)))
})

test_that("uncentered sparse crossproduct equals centered regression scores", {
  set.seed(20260726)
  X <- methods::as(
    Matrix::rsparsematrix(200, 80, density = 0.03), "dgCMatrix"
  )
  score <- rnorm(nrow(X))
  score <- score - mean(score)

  expected <- as.vector(crossprod(
    scale(as.matrix(X), center = TRUE, scale = FALSE), score
  ))
  actual <- as.vector(crossprod(X, score))
  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("computeNormalizedCorrelation works for CoProMulti aggregate mode", {
  obj <- create_test_copro_multi(n_cells_per_slide = 60, n_slides = 2,
                                 n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE, normalizeDistance = TRUE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100)

  slides <- getSlideList(obj)
  cts <- c("CellTypeA", "CellTypeB")
  X_scaled <- CoPro:::.preparePCMatrices(
    pc_data = obj@pcaResults, pca_global = obj@pcaGlobal,
    scalePCs = obj@scalePCs, slides = slides, cts = cts
  )
  weights <- obj@skrCCAOut[[CoPro:::.sigmaName(0.1)]]
  components <- lapply(slides, function(slide) {
    K <- getKernelMatrix(
      obj, 0.1, "CellTypeA", "CellTypeB", slide = slide, verbose = FALSE
    )
    A_w <- X_scaled[[slide]][["CellTypeA"]] %*%
      weights[["CellTypeA"]][, 1, drop = FALSE]
    B_w <- X_scaled[[slide]][["CellTypeB"]] %*%
      weights[["CellTypeB"]][, 1, drop = FALSE]
    list(
      numerator = CoPro:::.kernelXKY(A_w, K, B_w)[1, 1],
      null_sd = sqrt(sum(A_w^2)) * sqrt(sum(B_w^2)) *
        CoPro:::.whitenedFrobNorm(K)
    )
  })
  numerators <- vapply(components, `[[`, numeric(1), "numerator")
  slide_denominators <- vapply(components, `[[`, numeric(1), "null_sd")
  expected_sum <- sum(numerators) / sum(slide_denominators)
  expected_rss <- sum(numerators) / sqrt(sum(slide_denominators^2))

  sum_obj <- computeNormalizedCorrelation(
    obj, calculationMode = "aggregate", normalizer = "unwhitened"
  )
  rss_obj <- computeNormalizedCorrelation(
    obj, calculationMode = "aggregate", normalizer = "unwhitened",
    aggregateDenominator = "rss"
  )

  expect_true(length(sum_obj@normalizedCorrelation) > 0)
  agg_df <- sum_obj@normalizedCorrelation[[1]]
  expect_true(nrow(agg_df) > 0)
  expect_true("sigmaValue" %in% colnames(agg_df))
  expect_true("aggregateCorrelation" %in% colnames(agg_df))
  expect_true(is.numeric(agg_df$sigmaValue))
  observed_sum <- agg_df$aggregateCorrelation[
    agg_df$cellType1 == "CellTypeA" &
      agg_df$cellType2 == "CellTypeB" & agg_df$CC_index == 1
  ]
  rss_df <- rss_obj@normalizedCorrelation[[1]]
  observed_rss <- rss_df$aggregateCorrelation[
    rss_df$cellType1 == "CellTypeA" &
      rss_df$cellType2 == "CellTypeB" & rss_df$CC_index == 1
  ]
  expect_equal(observed_sum, expected_sum, tolerance = 1e-10)
  expect_equal(observed_rss, expected_rss, tolerance = 1e-10)
  expect_identical(
    attr(sum_obj@normalizedCorrelation, "aggregateDenominator"), "sum"
  )
  expect_identical(
    attr(rss_obj@normalizedCorrelation, "aggregateDenominator"), "rss"
  )

  transferred_scores <- setNames(lapply(cts, function(ct) {
    score <- do.call(rbind, lapply(slides, function(slide) {
      X_scaled[[slide]][[ct]] %*% weights[[ct]]
    }))
    colnames(score) <- paste0("CC_", seq_len(ncol(score)))
    score
  }), cts)
  transfer_sum <- getTransferNormCorr(
    obj, transferred_scores, sigma_choice = 0.1,
    calculationMode = "aggregate", normalizer = "unwhitened",
    verbose = FALSE
  )
  transfer_rss <- getTransferNormCorr(
    obj, transferred_scores, sigma_choice = 0.1,
    calculationMode = "aggregate", normalizer = "unwhitened",
    aggregateDenominator = "rss", verbose = FALSE
  )
  expect_equal(
    transfer_sum[[1]]$aggregateCorrelation,
    sum_obj@normalizedCorrelation[[1]]$aggregateCorrelation,
    tolerance = 1e-10
  )
  expect_equal(
    transfer_rss[[1]]$aggregateCorrelation,
    rss_obj@normalizedCorrelation[[1]]$aggregateCorrelation,
    tolerance = 1e-10
  )
  expect_identical(attr(transfer_sum, "aggregateDenominator"), "sum")
  expect_identical(attr(transfer_rss, "aggregateDenominator"), "rss")
})

test_that("column statistics and centering helpers match what they replaced", {
  set.seed(88)
  m <- matrix(rnorm(600 * 25), 600, 25)
  m[sample(length(m), 2000)] <- 0

  # .columnNonzeroFraction() reads a CsparseMatrix's column pointer instead of
  # building a full logical copy of the matrix.
  expect_identical(CoPro:::.columnNonzeroFraction(m),
                   colSums(m != 0) / nrow(m))
  sp <- as(as(as(Matrix::Matrix(m, sparse = TRUE), "generalMatrix"),
              "CsparseMatrix"), "dMatrix")
  expect_identical(as.numeric(CoPro:::.columnNonzeroFraction(sp)),
                   as.numeric(colSums(sp != 0) / nrow(sp)))
  # An explicitly stored zero is not a nonzero, so drop0() has to come first.
  explicit <- sp
  explicit@x[seq_len(15)] <- 0
  expect_identical(as.numeric(CoPro:::.columnNonzeroFraction(explicit)),
                   as.numeric(colSums(explicit != 0) / nrow(explicit)))

  # Centering by sweep() must equal the double-transpose form it replaced.
  expect_identical(CoPro:::.apply_centering_scaling(m, TRUE, FALSE),
                   t(t(m) - colMeans(m)))

  # Scaling without centering must guard degenerate genes the way every other
  # route into PCA does. Bare scale(center = FALSE, scale = TRUE) divides by a
  # zero root-mean-square and hands prcomp_irlba() a column of NaNs.
  deg <- m[, seq_len(4), drop = FALSE]
  colnames(deg) <- paste0("g", seq_len(4))
  deg[, 2] <- 0                        # never detected: divisor would be 0
  deg[, 3] <- 0; deg[seq_len(2), 3] <- 5  # detected in 2/600 cells
  scaled_only <- CoPro:::.apply_centering_scaling(deg, FALSE, TRUE)
  expect_false(anyNA(scaled_only))
  expect_identical(unname(attr(scaled_only, "scaled:scale")[c(2, 3)]), c(1, 1))
  centered_scaled <- CoPro:::center_scale_matrix_opt(deg)
  expect_identical(
    unname(attr(centered_scaled, "scaled:scale")[c(2, 3)]), c(1, 1)
  )
  # Undegenerate columns keep the exact divisor base::scale() would have used.
  expect_equal(
    unname(attr(scaled_only, "scaled:scale")[c(1, 4)]),
    unname(attr(scale(deg, center = FALSE, scale = TRUE),
                "scaled:scale")[c(1, 4)])
  )
  # ...and the guard agrees with the sparse route on the same matrix.
  deg_sp <- as(as(as(Matrix::Matrix(deg, sparse = TRUE), "generalMatrix"),
                  "CsparseMatrix"), "dMatrix")
  expect_equal(
    unname(attr(scaled_only, "scaled:scale")),
    unname(CoPro:::.sparse_pca_parameters(deg_sp, FALSE, TRUE)$scale)
  )
  within <- CoPro:::.withinSlidePCAParameters(
    deg, rep(c("s1", "s2"), each = 300), c("s1", "s2"),
    center = TRUE, scale. = TRUE
  )
  expect_identical(unname(within$scales[, c(2, 3)]),
                   matrix(1, nrow = 2, ncol = 2))

  # .columnSds() falls back to apply() for anything that is not a dense
  # numeric matrix, because matrixStats::colSds() does not accept one.
  expect_identical(CoPro:::.columnSds(sp), apply(sp, 2, stats::sd))
  # Both paths name their result from colnames(), or leave it unnamed.
  named <- m
  colnames(named) <- paste0("g", seq_len(ncol(named)))
  expect_identical(names(CoPro:::.columnSds(named)), colnames(named))
  expect_null(names(CoPro:::.columnSds(m)))
})

test_that("center_scale_matrix_opt() ships .columnSds(), i.e. matrixStats", {
  # This pins *which* SD implementation is in the shipped code path, which the
  # equivalence test below deliberately does not: that one shows the two agree
  # on everything a reader sees, and would pass under either. The two assert
  # different things and both are wanted.
  skip_on_cran()
  skip_if_not_installed("matrixStats")
  q <- function(e) suppressWarnings(suppressMessages(e))

  obj <- q(create_test_copro_single(n_cells = 320, n_genes = 60,
                                    n_cell_types = 2, seed = 11))
  obj <- q(subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB")))
  m <- obj@normalizedDataSub[obj@cellTypesSub == "CellTypeB", , drop = FALSE]

  by_matrix_stats <- matrixStats::colSds(m)
  by_apply <- apply(m, 2, stats::sd)

  # The assertions below can only tell the two implementations apart on input
  # where they actually disagree -- here, 1 ulp on Gene56. Whether they do is a
  # property of the installed matrixStats and BLAS, so establish it rather than
  # assume it, and skip instead of passing vacuously if it ever stops holding.
  skip_if(identical(unname(by_matrix_stats), unname(by_apply)),
          "the two SD implementations agree bit-for-bit on this fixture")

  expect_identical(unname(CoPro:::.columnSds(m)), unname(by_matrix_stats))

  # No column of this fixture is zero-SD or near-empty, so the scale attribute
  # is the raw SD vector with nothing substituted. Assert that, so a change to
  # the substitution rule shows up here rather than being absorbed.
  nz <- colSums(m != 0) / nrow(m)
  expect_length(which(by_matrix_stats < 1e-3 | nz < 0.01), 0L)

  scaled <- CoPro:::center_scale_matrix_opt(m)
  expect_identical(unname(attr(scaled, "scaled:scale")),
                   unname(by_matrix_stats))
})

test_that("the two column-SD implementations give the same science", {
  # matrixStats::colSds() and apply(x, 2, sd) use different variance algorithms
  # and disagree by 1 ulp on some columns. That is enough to flip the sign of a
  # principal component out of prcomp_irlba(), which sounds worse than it is:
  # the CCA weight's coordinate on that PC flips with it and the two cancel.
  #
  # Run the whole pipeline under each implementation and check what a reader of
  # the results would actually see. Component sign is compared with
  # .align_sign() because the sign of a CCA axis is inherently ambiguous -- a
  # different BLAS or R build can flip it under the existing code too -- while
  # the sign-invariant quantities are held to a tight tolerance.
  skip_on_cran()
  skip_if_not_installed("matrixStats")
  q <- function(e) suppressWarnings(suppressMessages(e))

  # Reproduce center_scale_matrix_opt() with a supplied SD function, per cell
  # type, then hand the result to computePCA(center = FALSE, scale. = FALSE).
  # That exercises the real difference between the two implementations without
  # mocking an internal binding.
  prescale <- function(m, sd_fun) {
    col_means <- colMeans(m)
    col_sds <- sd_fun(m)
    col_nz <- colSums(m != 0) / nrow(m)
    safe <- col_sds
    drop <- which(col_sds < 1e-3 | col_nz < 0.01)
    if (length(drop) > 0) safe[drop] <- 1.0
    scale(m, center = col_means, scale = safe)
  }

  run <- function(sd_fun) {
    o <- q(create_test_copro_single(n_cells = 320, n_genes = 60,
                                    n_cell_types = 2, seed = 11))
    o <- q(subsetData(o, cellTypesOfInterest = c("CellTypeA", "CellTypeB")))
    for (ct in c("CellTypeA", "CellTypeB")) {
      rows <- o@cellTypesSub == ct
      o@normalizedDataSub[rows, ] <- prescale(
        o@normalizedDataSub[rows, , drop = FALSE], sd_fun)
    }
    o <- q(computePCA(o, nPCA = 8, center = FALSE, scale. = FALSE))
    o <- q(computeSparseKernel(o, sigmaValues = c(0.05, 0.1), verbose = FALSE, normalizeDistance = TRUE))
    o <- q(runSkrCCA(o, scalePCs = TRUE, nCC = 2))
    o <- q(computeNormalizedCorrelation(o))
    o <- q(computeGeneAndCellScores(o))
    q(computeRegressionGeneScores(o))
  }
  with_matrix_stats <- run(function(x) matrixStats::colSds(x))
  with_apply <- run(function(x) apply(x, 2, stats::sd))

  # Sign-invariant outputs -- the ones the manuscript reports -- must agree
  # tightly, and the selected bandwidth must be the same.
  expect_equal(
    do.call(rbind, with_matrix_stats@normalizedCorrelation)$normalizedCorrelation,
    do.call(rbind, with_apply@normalizedCorrelation)$normalizedCorrelation,
    tolerance = 1e-10
  )
  expect_identical(with_matrix_stats@sigmaValueChoice,
                   with_apply@sigmaValueChoice)

  # Scores and gene weights must agree up to the per-component sign.
  compare_aligned <- function(a, b, label) {
    keys <- intersect(names(a), names(b))
    expect_gt(length(keys), 0L)
    for (k in keys) {
      va <- as.matrix(a[[k]])
      vb <- as.matrix(b[[k]])
      expect_identical(dim(va), dim(vb), label = paste(label, k, "dim"))
      for (cc in seq_len(ncol(va))) {
        expect_equal(unname(.align_sign(va[, cc], vb[, cc])),
                     unname(va[, cc]), tolerance = 1e-8,
                     label = paste(label, k, "component", cc))
      }
    }
  }
  compare_aligned(with_matrix_stats@cellScores, with_apply@cellScores,
                  "cellScores")
  compare_aligned(with_matrix_stats@geneScores, with_apply@geneScores,
                  "geneScores")
  compare_aligned(with_matrix_stats@geneScoresRegression,
                  with_apply@geneScoresRegression, "geneScoresRegression")
})

test_that("per-slide PC scores are stored as views of the global scores", {
  obj <- create_test_copro_multi(n_cells_per_slide = 120, n_slides = 2,
                                 n_cell_types = 2, seed = 91)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- suppressMessages(computePCA(obj, nPCA = 6))

  # Slot shape is unchanged: one entry per slide, keyed by slide name.
  expect_length(obj@pcaResults, 2L)
  expect_setequal(names(obj@pcaResults), c("Slide1", "Slide2"))

  for (sID in c("Slide1", "Slide2")) {
    for (ct in c("CellTypeA", "CellTypeB")) {
      entry <- obj@pcaResults[[sID]][[ct]]
      expect_true(CoPro:::.isPCSlice(entry))
      expect_type(entry$rows, "integer")

      # Resolving it must give exactly what the old code materialized: the
      # rows of the global scores for this slide, labelled with cell IDs.
      resolved <- CoPro:::.resolvePCSlice(entry, obj@pcaGlobal[[ct]]$x)
      keep <- getSlideID(obj)[obj@cellTypesSub == ct] == sID
      expected <- obj@pcaGlobal[[ct]]$x[keep, , drop = FALSE]
      expect_identical(resolved, expected)
      expect_identical(
        rownames(resolved),
        rownames(obj@metaDataSub)[obj@cellTypesSub == ct][keep]
      )
    }
  }

  # The rows partition the cell type exactly once.
  for (ct in c("CellTypeA", "CellTypeB")) {
    all_rows <- sort(unname(unlist(
      lapply(obj@pcaResults, function(s) s[[ct]]$rows))))
    expect_identical(all_rows, seq_len(nrow(obj@pcaGlobal[[ct]]$x)))
  }
})

test_that("within-slide preprocessing keeps shared score views, and legacy matrices still resolve", {
  make <- function(center_per_slide) {
    obj <- create_test_copro_multi(n_cells_per_slide = 120, n_slides = 2,
                                   n_cell_types = 2, seed = 91)
    obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
    suppressMessages(computePCA(obj, nPCA = 6,
                                center_per_slide = center_per_slide))
  }

  # Centering now happens in gene space before the shared projection. Every
  # slide therefore remains a view of one global PC matrix, while its score
  # block is centered without a second, post-PCA transformation.
  centered <- make(TRUE)
  entry <- centered@pcaResults[["Slide1"]][["CellTypeA"]]
  expect_true(CoPro:::.isPCSlice(entry))
  centered_scores <- CoPro:::.resolvePCSlice(
    entry, centered@pcaGlobal[["CellTypeA"]]$x
  )
  expect_equal(colMeans(centered_scores), rep(0, ncol(centered_scores)),
               tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_identical(centered@pcaGlobal[["CellTypeA"]]$preprocessing,
                   "within_slide")

  # The stored shared loading must reproduce every score block after applying
  # that block's stored gene-space affine map.
  ct <- "CellTypeA"
  ct_rows <- centered@cellTypesSub == ct
  X_ct <- centered@normalizedDataSub[ct_rows, , drop = FALSE]
  slide_ct <- getSlideID(centered)[ct_rows]
  pca <- centered@pcaGlobal[[ct]]
  for (sID in getSlideList(centered)) {
    rows <- which(slide_ct == sID)
    Z <- sweep(X_ct[rows, , drop = FALSE], 2L,
               pca$slideCenter[sID, ], "-")
    Z <- sweep(Z, 2L, pca$slideScale[sID, ], "/")
    expect_equal(
      unname(Z %*% pca$rotation),
      unname(pca$x[rows, , drop = FALSE]),
      tolerance = 1e-8,
      info = sID
    )
  }

  # An object saved before slices became views holds matrices. Both
  # representations must prepare to the same thing, for both scalePCs settings.
  lazy <- make(FALSE)
  legacy <- lazy
  for (sID in names(legacy@pcaResults)) {
    for (ct in names(legacy@pcaResults[[sID]])) {
      legacy@pcaResults[[sID]][[ct]] <- CoPro:::.resolvePCSlice(
        legacy@pcaResults[[sID]][[ct]], legacy@pcaGlobal[[ct]]$x)
    }
  }
  expect_true(is.matrix(legacy@pcaResults[["Slide1"]][["CellTypeA"]]))

  for (scale_pcs in c(TRUE, FALSE)) {
    from_view <- CoPro:::.preparePCMatrices(
      pc_data = lazy@pcaResults, pca_global = lazy@pcaGlobal,
      scalePCs = scale_pcs, slides = getSlideList(lazy),
      cts = c("CellTypeA", "CellTypeB"))
    from_matrix <- CoPro:::.preparePCMatrices(
      pc_data = legacy@pcaResults, pca_global = legacy@pcaGlobal,
      scalePCs = scale_pcs, slides = getSlideList(legacy),
      cts = c("CellTypeA", "CellTypeB"))
    for (sID in names(from_view)) {
      for (ct in names(from_view[[sID]])) {
        expect_equal(from_view[[sID]][[ct]], from_matrix[[sID]][[ct]],
                     ignore_attr = TRUE,
                     info = paste(sID, ct, "scalePCs", scale_pcs))
      }
    }
  }
})

test_that("multi-slide PCA default is invariant to per-slide positive affine batch effects", {
  base <- create_test_copro_multi(
    n_cells_per_slide = 140, n_slides = 3, n_genes = 35,
    n_cell_types = 2, seed = 771
  )
  base <- subsetData(base, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  shifted <- base

  slide_id <- getSlideID(shifted)
  p <- ncol(shifted@normalizedDataSub)
  for (k in seq_along(getSlideList(shifted))) {
    rows <- slide_id == getSlideList(shifted)[k]
    multiplier <- seq(0.55, 1.85, length.out = p) * (0.8 + 0.25 * k)
    offset <- seq(-3, 4, length.out = p) + 7 * k
    block <- sweep(shifted@normalizedDataSub[rows, , drop = FALSE],
                   2L, multiplier, "*")
    shifted@normalizedDataSub[rows, ] <- sweep(block, 2L, offset, "+")
  }

  set.seed(91)
  fit_base <- suppressMessages(computePCA(base, nPCA = 8))
  set.seed(91)
  fit_shifted <- suppressMessages(computePCA(shifted, nPCA = 8))

  for (ct in c("CellTypeA", "CellTypeB")) {
    A <- fit_base@pcaGlobal[[ct]]
    B <- fit_shifted@pcaGlobal[[ct]]
    expect_identical(A$preprocessing, "within_slide")
    expect_identical(B$preprocessing, "within_slide")
    expect_equal(A$sdev, B$sdev, tolerance = 1e-9)
    for (pc in seq_len(8L)) {
      sign <- if (sum(A$rotation[, pc] * B$rotation[, pc]) < 0) -1 else 1
      expect_equal(A$rotation[, pc], sign * B$rotation[, pc],
                   tolerance = 1e-8, info = paste(ct, "loading", pc))
      expect_equal(A$x[, pc], sign * B$x[, pc],
                   tolerance = 1e-8, info = paste(ct, "score", pc))
    }
  }
})

test_that("matrix-free sparse within-slide PCA matches dense preprocessing", {
  set.seed(992)
  dense <- matrix(rpois(3600, lambda = 0.7), nrow = 90, ncol = 40,
                  dimnames = list(paste0("c", seq_len(90)),
                                  paste0("g", seq_len(40))))
  sparse <- methods::as(Matrix::Matrix(dense, sparse = TRUE), "dgCMatrix")
  slides <- rep(c("s1", "s2", "s3"), each = 30)

  set.seed(19)
  a <- CoPro:::.run_within_slide_pca(
    dense, slides, unique(slides), 6, center = TRUE, scale. = TRUE
  )
  set.seed(19)
  b <- CoPro:::.run_within_slide_pca(
    sparse, slides, unique(slides), 6, center = TRUE, scale. = TRUE
  )

  expect_equal(a$sdev, b$sdev, tolerance = 1e-8)
  for (pc in seq_len(6L)) {
    sign <- if (sum(a$rotation[, pc] * b$rotation[, pc]) < 0) -1 else 1
    expect_equal(a$rotation[, pc], sign * b$rotation[, pc], tolerance = 2e-5)
    expect_equal(a$x[, pc], sign * b$x[, pc], tolerance = 2e-5)
  }
})

test_that("matrix-free BPCells within-slide PCA matches dense preprocessing", {
  # BPCells remains optional for local installs. The dedicated BPCells CI job
  # installs it and makes this required coverage of the IterableMatrix branch
  # of .run_within_slide_pca(), which the multi-slide default routes
  # out-of-core input through.
  skip_if_not_installed("BPCells")
  set.seed(4242)
  dense <- matrix(rpois(3600, lambda = 0.7), nrow = 90, ncol = 40,
                  dimnames = list(paste0("c", seq_len(90)),
                                  paste0("g", seq_len(40))))
  bp <- BPCells::write_matrix_memory(
    methods::as(Matrix::Matrix(dense, sparse = TRUE), "dgCMatrix"),
    compress = FALSE
  )
  slides <- rep(c("s1", "s2", "s3"), each = 30)

  set.seed(19)
  a <- CoPro:::.run_within_slide_pca(
    dense, slides, unique(slides), 6, center = TRUE, scale. = TRUE
  )
  set.seed(19)
  b <- suppressMessages(CoPro:::.run_within_slide_pca(
    bp, slides, unique(slides), 6, center = TRUE, scale. = TRUE
  ))

  expect_equal(a$sdev, b$sdev, tolerance = 1e-6)
  for (pc in seq_len(6L)) {
    sign <- if (sum(a$rotation[, pc] * b$rotation[, pc]) < 0) -1 else 1
    expect_equal(a$rotation[, pc], sign * b$rotation[, pc], tolerance = 2e-4)
    expect_equal(a$x[, pc], sign * b$x[, pc], tolerance = 2e-4)
  }
})

test_that(".apply_centering_scaling() keeps BPCells input out of core", {
  # computePCA() routes IterableMatrix input through this helper for every
  # combination of center/scale., but only the center = scale. = TRUE branch
  # was ever taught about BPCells. sweep() errored outright ("non-numeric
  # argument to binary operator") and base::scale() silently densified the
  # matrix, which is what BPCells exists to avoid.
  skip_if_not_installed("BPCells")
  set.seed(9)
  dense <- matrix(rpois(800, lambda = 0.6), nrow = 200, ncol = 4,
                  dimnames = list(paste0("c", seq_len(200)),
                                  paste0("g", seq_len(4))))
  dense[, 2] <- 0                    # never detected: zero divisor
  dense[, 3] <- 0; dense[1, 3] <- 5  # detected in 1 of 200 cells: under 1%
  bp <- BPCells::write_matrix_memory(
    methods::as(Matrix::Matrix(dense, sparse = TRUE), "dgCMatrix"),
    compress = FALSE
  )

  for (scale. in c(FALSE, TRUE)) {
    center <- !scale.
    out_bp <- CoPro:::.apply_centering_scaling(bp, center, scale.)
    out_dense <- CoPro:::.apply_centering_scaling(dense, center, scale.)
    expect_true(CoPro:::.is_bpcells(out_bp))
    # BPCells prints a "converting to a dense matrix" notice here. Densifying
    # is exactly what the helper must not do; the test does it deliberately on
    # a 200 x 4 fixture to compare values, and the notice is printed rather
    # than signaled, so it cannot be suppressed.
    expect_equal(as.matrix(out_bp), as.matrix(out_dense),
                 ignore_attr = TRUE, tolerance = 1e-12)
  }

  # .columnNonzeroFraction() feeds every degeneracy guard in the package, and
  # `x != 0` is not defined for an IterableMatrix at all.
  expect_equal(unname(CoPro:::.columnNonzeroFraction(bp)),
               unname(colSums(dense != 0) / nrow(dense)))
  expect_equal(unname(CoPro:::.uncenteredColumnScales(bp)),
               unname(CoPro:::.uncenteredColumnScales(dense)))
  expect_identical(unname(CoPro:::.uncenteredColumnScales(bp)[c(2, 3)]),
                   c(1, 1))
})

test_that("BPCells nonzero fractions count negatives, not stored zeros", {
  skip_if_not_installed("BPCells")
  set.seed(4)
  dense <- cbind(
    pos = abs(rnorm(200)) + 0.5,
    neg = -abs(rnorm(200)) - 0.5,
    mixed = c(rep(-1, 100), rep(1, 100))
  )
  sparse <- methods::as(
    methods::as(Matrix::Matrix(dense, sparse = TRUE), "generalMatrix"),
    "CsparseMatrix"
  )
  bp <- BPCells::write_matrix_memory(sparse, compress = FALSE)

  expected <- rep(1, ncol(dense))
  expect_identical(
    as.numeric(CoPro:::.columnNonzeroFraction(dense)), expected
  )
  expect_identical(
    as.numeric(CoPro:::.columnNonzeroFraction(sparse)), expected
  )
  expect_identical(
    as.numeric(CoPro:::.columnNonzeroFraction(bp)), expected
  )
  expect_equal(
    as.numeric(CoPro:::.uncenteredColumnScales(bp)),
    as.numeric(CoPro:::.uncenteredColumnScales(dense)),
    tolerance = 1e-12
  )
  # The frozen-reference guard rides on the same helper: nothing here is
  # degenerate, so no gene is flagged (binarize() would have flagged `neg`).
  expect_identical(
    unname(CoPro:::.frozen_score_guard(bp, rep("s1", nrow(dense)), FALSE)),
    c(FALSE, FALSE, FALSE)
  )

  # An explicitly stored zero is still a zero. Matrix::Matrix() drops stored
  # zeros on conversion, so build the fixture with sparseMatrix(), which keeps
  # them, and confirm BPCells carries the stored pattern through before
  # trusting the count: matrix_stats(col_stats = "nonzero") reports stored
  # entries, which is exactly why .columnNonzeroFraction() cannot use it.
  stored <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 4, 1, 2, 1, 3),
    j = c(1, 1, 1, 1, 2, 2, 3, 3),
    x = c(-1, -2, 0, 3, 0, 0, 5, 7),
    dims = c(4, 4)
  )
  expect_identical(diff(stored@p), c(4L, 2L, 2L, 0L))
  stored_bp <- BPCells::write_matrix_memory(stored, compress = FALSE)
  stored_counts <- BPCells::matrix_stats(stored_bp, col_stats = "nonzero")
  expect_identical(
    as.numeric(stored_counts$col_stats["nonzero", ]), c(4, 2, 2, 0)
  )

  stored_expected <- c(3, 0, 2, 0) / 4
  expect_identical(
    as.numeric(CoPro:::.columnNonzeroFraction(as.matrix(stored))),
    stored_expected
  )
  expect_identical(
    as.numeric(CoPro:::.columnNonzeroFraction(stored)), stored_expected
  )
  expect_identical(
    as.numeric(CoPro:::.columnNonzeroFraction(stored_bp)), stored_expected
  )
  # Through the frozen-reference guard only the all-zero and the empty column
  # are degenerate; the columns carrying stored zeros are not.
  expect_identical(
    unname(CoPro:::.frozen_score_guard(stored_bp, rep("s1", 4), FALSE)),
    c(FALSE, TRUE, FALSE, TRUE)
  )
})

test_that("the degeneracy guard treats a non-finite or NA scale as unsafe", {
  # Every route into PCA and the frozen-reference guard share one predicate.
  # A finite scale at or above threshold with a nonzero fraction at or above
  # threshold is the only safe case; an NA on either input must come out TRUE,
  # never NA, because `x[NA] <- 1` silently skips that element.
  unsafe <- CoPro:::.unsafeScaleColumns(
    scale_values = c(2, NaN, Inf, NA, 1e-4, 1e-3, 2, 2),
    nonzero_fraction = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.01, 0.001, NA)
  )
  expect_identical(
    unsafe, c(FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
  )
  expect_false(anyNA(unsafe))

  # The non-finite arm is what a one-cell block hits: the sd of a single row
  # is NA on the dense path, and the divisor must still come out as 1.
  one <- matrix(c(1, 2, 3), nrow = 1,
                dimnames = list("c1", c("g1", "g2", "g3")))
  expect_identical(unname(CoPro:::.columnSds(one)), rep(NA_real_, 3))
  dense_one <- CoPro:::center_scale_matrix_opt(one)
  expect_identical(unname(attr(dense_one, "scaled:scale")), c(1, 1, 1))
  expect_identical(as.numeric(dense_one), c(0, 0, 0))

  # A slide holding one cell guards every gene on every slide.
  set.seed(7)
  m <- matrix(rnorm(15), 5, 3, dimnames = list(NULL, c("g1", "g2", "g3")))
  within <- CoPro:::.withinSlidePCAParameters(
    m, c("a", "a", "a", "a", "b"), c("a", "b"), center = TRUE, scale. = TRUE
  )
  expect_identical(unname(within$guarded), c(TRUE, TRUE, TRUE))
  expect_identical(unname(within$scales), matrix(1, nrow = 2, ncol = 3))
})

test_that("the degeneracy guard pins a one-cell BPCells block to unit scale", {
  skip_if_not_installed("BPCells")
  one <- matrix(c(1, 2, 3), nrow = 1,
                dimnames = list("c1", c("g1", "g2", "g3")))
  sparse_one <- methods::as(
    methods::as(Matrix::Matrix(one, sparse = TRUE), "generalMatrix"),
    "CsparseMatrix"
  )
  bp_one <- BPCells::write_matrix_memory(sparse_one, compress = FALSE)
  # The variance of one row is NaN, the non-finite arm of the guard.
  expect_true(all(is.nan(as.numeric(CoPro:::.bpcellsColumnVariances(bp_one)))))
  scaled <- CoPro:::center_scale_matrix_opt(bp_one)
  expect_s4_class(scaled, "IterableMatrix")
  expect_identical(as.numeric(as.matrix(scaled)), c(0, 0, 0))
})

test_that("BPCells column variances never dispatch through colVars()", {
  skip_if_not_installed("BPCells")
  # BPCells::colVars() only reaches its IterableMatrix method when BPCells is
  # attached or MatrixGenerics is installed; the BPCells CI job has neither
  # and failed on every BPCells route through PCA the first time it ran. The
  # helper reads the same streamed statistic through matrix_stats() instead.
  set.seed(11)
  dense <- matrix(rpois(300 * 4, 2) * rnorm(1200), 300, 4,
                  dimnames = list(NULL, paste0("g", 1:4)))
  dense[, 4] <- 0
  sparse <- methods::as(
    methods::as(Matrix::Matrix(dense, sparse = TRUE), "generalMatrix"),
    "CsparseMatrix"
  )
  bp <- BPCells::write_matrix_memory(sparse, compress = FALSE)
  vars <- CoPro:::.bpcellsColumnVariances(bp)
  expect_identical(names(vars), colnames(dense))
  expect_equal(unname(vars), unname(apply(dense, 2, stats::var)),
               tolerance = 1e-12)
  expect_identical(length(CoPro:::.bpcellsColumnVariances(bp[, integer(0)])),
                   0L)
  expect_equal(
    unname(CoPro:::.bpcellsColumnVariances(bp[, 2, drop = FALSE])),
    stats::var(dense[, 2]), tolerance = 1e-12
  )

  # Pin the fix at the source: none of the functions that standardize BPCells
  # input may mention colVars, whichever library this suite runs in.
  for (fn in c(".bpcellsColumnVariances", "center_scale_matrix_opt",
               ".withinSlidePCAParameters", ".frozen_column_sds")) {
    code <- paste(deparse(get(fn, envir = asNamespace("CoPro"))),
                  collapse = "\n")
    expect_false(grepl("colVars", code, fixed = TRUE), label = fn)
  }
})

test_that("the legacy multi-slide combination still runs the legacy path", {
  # NEWS and ?runSkrCCA both promise center_per_slide = FALSE plus
  # objective = "sumcov" reproduces the pre-1.3.0 workflow. Nothing else
  # exercises the two together, and the shared .runSkrCCAUnified() /
  # .applySlideAdequacy() path was substantially refactored around them.
  obj <- suppressWarnings(create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30, n_cell_types = 2,
    seed = 808
  ))
  cts <- c("CellTypeA", "CellTypeB")
  obj <- subsetData(obj, cellTypesOfInterest = cts)

  legacy <- suppressMessages(computePCA(obj, nPCA = 5,
                                        center_per_slide = FALSE))
  expect_identical(legacy@pcaGlobal[[cts[1]]]$preprocessing, "pooled")
  # Pooled preprocessing keeps the raw-unit affine map on the prcomp object and
  # records no per-slide maps.
  expect_null(legacy@pcaGlobal[[cts[1]]]$slideCenter)
  expect_null(legacy@pcaGlobal[[cts[1]]]$slideScale)

  legacy <- suppressMessages(computeDistance(
    legacy, distType = "Euclidean2D", normalizeDistance = TRUE, verbose = FALSE
  ))
  legacy <- suppressMessages(computeKernelMatrix(
    legacy, sigmaValues = 0.1, verbose = FALSE, normalizeDistance = TRUE
  ))
  legacy <- suppressMessages(suppressWarnings(
    runSkrCCA(legacy, nCC = 1, objective = "sumcov")
  ))

  expect_equal(getCCAObjective(legacy)$objective, "sumcov")
  # sumcov fixes the slide weighting by construction, so the recorded
  # slideWeight is a missing value rather than one of the sumcor choices.
  expect_true(is.na(getCCAObjective(legacy)$slideWeight))
  for (ct in cts) {
    w <- legacy@skrCCAOut[["sigma_0.1"]][[ct]]
    expect_equal(nrow(w), 5L)
    expect_equal(sqrt(sum(w[, 1]^2)), 1, tolerance = 1e-8)
  }
})
