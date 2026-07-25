# Extracted from test-pca-workflow.R:454

# test -------------------------------------------------------------------------
skip_on_cran()
skip_if_not_installed("matrixStats")
q <- function(e) suppressWarnings(suppressMessages(e))
run <- function(sd_fun) {
    testthat::local_mocked_bindings(.columnSds = sd_fun, .package = "CoPro")
    o <- q(create_test_copro_single(n_cells = 320, n_genes = 60,
                                    n_cell_types = 2, seed = 11))
    o <- q(subsetData(o, cellTypesOfInterest = c("CellTypeA", "CellTypeB")))
    o <- q(computePCA(o, nPCA = 8))
    o <- q(computeSparseKernel(o, sigmaValues = c(0.05, 0.1), verbose = FALSE))
    o <- q(runSkrCCA(o, scalePCs = TRUE, nCC = 2))
    o <- q(computeNormalizedCorrelation(o))
    o <- q(computeGeneAndCellScores(o))
    q(computeRegressionGeneScores(o))
  }
