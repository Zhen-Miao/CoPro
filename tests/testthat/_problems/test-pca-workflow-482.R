# Extracted from test-pca-workflow.R:482

# test -------------------------------------------------------------------------
skip_on_cran()
skip_if_not_installed("matrixStats")
q <- function(e) suppressWarnings(suppressMessages(e))
obj <- q(create_test_copro_single(n_cells = 320, n_genes = 60,
                                    n_cell_types = 2, seed = 11))
obj <- q(subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB")))
m <- obj@normalizedDataSub[obj@cellTypesSub == "CellTypeB", , drop = FALSE]
by_apply <- apply(m, 2, stats::sd)
by_matrix_stats <- matrixStats::colSds(m)
expect_false(identical(by_apply, by_matrix_stats))
expect_lt(max(abs(by_apply - by_matrix_stats) / by_apply), 1e-15)
run <- function(perturb) {
    o <- obj
    if (perturb) {
      worst <- which.max(abs(by_apply - by_matrix_stats))
      o@normalizedDataSub[, worst] <-
        o@normalizedDataSub[, worst] * (1 + .Machine$double.eps)
    }
    o <- q(computePCA(o, nPCA = 8))
    o <- q(computeSparseKernel(o, sigmaValues = c(0.05, 0.1), verbose = FALSE))
    o <- q(runSkrCCA(o, scalePCs = TRUE, nCC = 2))
    o <- q(computeNormalizedCorrelation(o))
    o <- q(computeGeneAndCellScores(o))
    q(computeRegressionGeneScores(o))
  }
base <- run(FALSE)
nudged <- run(TRUE)
no_flip <- function(a, b, label) {
    va <- as.numeric(unlist(a))
    vb <- as.numeric(unlist(b))
    expect_lt(max(abs(va - vb)), 1e-8, label = label)
    expect_lt(max(abs(va - vb)), max(abs(va + vb)),
              label = paste(label, "sign"))
  }
