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
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE)
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
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE)
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  
  # Run skrCCA
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100)
  
  expect_true(length(obj@skrCCAOut) > 0)
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
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1), verbose = FALSE)
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
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE)
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
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE)
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
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1), verbose = FALSE)
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
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE)
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

test_that("computeNormalizedCorrelation works for CoProMulti aggregate mode", {
  obj <- create_test_copro_multi(n_cells_per_slide = 60, n_slides = 2,
                                 n_cell_types = 2, seed = 42)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE)
  obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2, maxIter = 100)

  obj <- computeNormalizedCorrelation(obj, calculationMode = "aggregate")

  expect_true(length(obj@normalizedCorrelation) > 0)
  agg_df <- obj@normalizedCorrelation[[1]]
  expect_true(nrow(agg_df) > 0)
  expect_true("sigmaValue" %in% colnames(agg_df))
  expect_true("aggregateCorrelation" %in% colnames(agg_df))
  expect_true(is.numeric(agg_df$sigmaValue))
})

