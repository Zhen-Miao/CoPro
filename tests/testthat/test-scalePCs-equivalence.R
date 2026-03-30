# Tests that scalePCs = TRUE and scalePCs = FALSE produce equivalent results.
#
# When scalePCs = TRUE the PCs are whitened and ||w|| = 1 is the constraint.
# When scalePCs = FALSE the PCs are raw and w'Dw = 1 is enforced instead.
# Both formulations solve the same CCA problem, so cell scores Xw and
# gene scores (rotation %*% gene-space weight) must agree (up to sign).

# ---------------------------------------------------------------------------
# Helper: run a full single-slide pipeline and return the object
# ---------------------------------------------------------------------------
.run_single_pipeline <- function(scalePCs, nCC = 2, seed = 42) {
  obj <- create_test_copro_single(n_cells = 100, n_cell_types = 2, seed = seed)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE)
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE,
                    scalePCs = scalePCs)

  suppressWarnings(
    obj <- runSkrCCA(obj, scalePCs = scalePCs, nCC = nCC, maxIter = 200)
  )
  obj <- computeGeneAndCellScores(obj)
  obj
}

# ---------------------------------------------------------------------------
# Helper: run a full multi-slide pipeline and return the object
# ---------------------------------------------------------------------------
.run_multi_pipeline <- function(scalePCs, nCC = 2, seed = 42) {
  obj <- create_test_copro_multi(n_cells_per_slide = 60, n_slides = 2,
                                  n_cell_types = 2, seed = seed)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.1), verbose = FALSE)
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE,
                    scalePCs = scalePCs)

  suppressWarnings(
    obj <- runSkrCCA(obj, scalePCs = scalePCs, nCC = nCC, maxIter = 200)
  )
  obj <- computeGeneAndCellScores(obj)
  obj
}

# ---------------------------------------------------------------------------
# Helper: compare scores across two objects (sign-agnostic)
# Tolerance is 1e-3 because the two optimization paths (whitened unit-norm
# vs raw weighted-norm) converge via different numerical conditioning, so
# small differences (typically < 1e-4) are expected.
# ---------------------------------------------------------------------------
.compare_scores <- function(scores_a, scores_b, label, tol = 1e-3) {
  nms_a <- names(scores_a)
  nms_b <- names(scores_b)
  expect_equal(sort(nms_a), sort(nms_b))

  for (nm in nms_a) {
    mat_a <- scores_a[[nm]]
    mat_b <- scores_b[[nm]]
    expect_equal(dim(mat_a), dim(mat_b), info = paste("dim mismatch for", nm))

    for (cc in seq_len(ncol(mat_a))) {
      a <- mat_a[, cc]
      b <- mat_b[, cc]
      if (sum((a - b)^2) > sum((a + b)^2)) {
        b <- -b
      }
      expect_equal(a, b, tolerance = tol,
                   info = paste(label, "mismatch for", nm, "CC", cc))
    }
  }
}

# ===========================================================================
# Single-slide tests
# ===========================================================================

test_that("scalePCs TRUE vs FALSE produce same scores (single slide, nCC=1)", {
  obj_true  <- .run_single_pipeline(scalePCs = TRUE,  nCC = 1)
  obj_false <- .run_single_pipeline(scalePCs = FALSE, nCC = 1)
  .compare_scores(obj_true@cellScores, obj_false@cellScores, "cell score")
  .compare_scores(obj_true@geneScores, obj_false@geneScores, "gene score")
})

test_that("scalePCs TRUE vs FALSE produce same scores (single slide, nCC=2)", {
  obj_true  <- .run_single_pipeline(scalePCs = TRUE,  nCC = 2)
  obj_false <- .run_single_pipeline(scalePCs = FALSE, nCC = 2)
  .compare_scores(obj_true@cellScores, obj_false@cellScores, "cell score")
  .compare_scores(obj_true@geneScores, obj_false@geneScores, "gene score")
})

test_that("scalePCs TRUE vs FALSE produce same scores (single slide, nCC=3)", {
  obj_true  <- .run_single_pipeline(scalePCs = TRUE,  nCC = 3)
  obj_false <- .run_single_pipeline(scalePCs = FALSE, nCC = 3)
  .compare_scores(obj_true@cellScores, obj_false@cellScores, "cell score")
  .compare_scores(obj_true@geneScores, obj_false@geneScores, "gene score")
})

# ===========================================================================
# Multi-slide tests
# ===========================================================================

test_that("scalePCs TRUE vs FALSE produce same scores (multi slide, nCC=1)", {
  obj_true  <- .run_multi_pipeline(scalePCs = TRUE,  nCC = 1)
  obj_false <- .run_multi_pipeline(scalePCs = FALSE, nCC = 1)
  .compare_scores(obj_true@cellScores, obj_false@cellScores, "cell score")
  .compare_scores(obj_true@geneScores, obj_false@geneScores, "gene score")
})

test_that("scalePCs TRUE vs FALSE produce same scores (multi slide, nCC=2)", {
  obj_true  <- .run_multi_pipeline(scalePCs = TRUE,  nCC = 2)
  obj_false <- .run_multi_pipeline(scalePCs = FALSE, nCC = 2)
  .compare_scores(obj_true@cellScores, obj_false@cellScores, "cell score")
  .compare_scores(obj_true@geneScores, obj_false@geneScores, "gene score")
})

test_that("scalePCs TRUE vs FALSE produce same scores (multi slide, nCC=3)", {
  obj_true  <- .run_multi_pipeline(scalePCs = TRUE,  nCC = 3)
  obj_false <- .run_multi_pipeline(scalePCs = FALSE, nCC = 3)
  .compare_scores(obj_true@cellScores, obj_false@cellScores, "cell score")
  .compare_scores(obj_true@geneScores, obj_false@geneScores, "gene score")
})
