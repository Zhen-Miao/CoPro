# Tests that scalePCs = TRUE and scalePCs = FALSE produce equivalent cell scores.
#
# When scalePCs = TRUE the PCs are whitened and ||w|| = 1 is the constraint.
# When scalePCs = FALSE the PCs are raw and w'Dw = 1 is enforced instead.
# Both formulations solve the same CCA problem, so cell scores Xw must agree
# (up to sign) between the two paths.
#
# NOTE: Gene scores are NOT expected to match because the gene score formula
# uses `w * sdev` (scalePCs=TRUE) vs `w` (scalePCs=FALSE). These represent
# different but valid gene importance measures. Only cell scores (Xw) are
# mathematically guaranteed to be equivalent.

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
# Helper: compare cell scores across two objects (sign-agnostic)
# Tolerance is 1e-3 because the two optimization paths (whitened unit-norm
# vs raw weighted-norm) converge via different numerical conditioning, so
# small differences (typically < 1e-4) are expected.
# ---------------------------------------------------------------------------
.compare_cell_scores <- function(obj_a, obj_b, tol = 1e-3) {
  nms_a <- names(obj_a@cellScores)
  nms_b <- names(obj_b@cellScores)
  expect_equal(sort(nms_a), sort(nms_b))

  for (nm in nms_a) {
    cs_a <- obj_a@cellScores[[nm]]
    cs_b <- obj_b@cellScores[[nm]]
    expect_equal(dim(cs_a), dim(cs_b), info = paste("dim mismatch for", nm))

    # Compare column by column (each CC), allowing sign flip
    for (cc in seq_len(ncol(cs_a))) {
      a <- cs_a[, cc]
      b <- cs_b[, cc]
      # Align sign: pick the sign that minimises residual
      if (sum((a - b)^2) > sum((a + b)^2)) {
        b <- -b
      }
      expect_equal(a, b, tolerance = tol,
                   info = paste("cell score mismatch for", nm, "CC", cc))
    }
  }
}

# ===========================================================================
# Single-slide tests
# ===========================================================================

test_that("scalePCs TRUE vs FALSE produce same cell scores (single slide, nCC=1)", {
  obj_true  <- .run_single_pipeline(scalePCs = TRUE,  nCC = 1)
  obj_false <- .run_single_pipeline(scalePCs = FALSE, nCC = 1)
  .compare_cell_scores(obj_true, obj_false)
})

test_that("scalePCs TRUE vs FALSE produce same cell scores (single slide, nCC=2)", {
  obj_true  <- .run_single_pipeline(scalePCs = TRUE,  nCC = 2)
  obj_false <- .run_single_pipeline(scalePCs = FALSE, nCC = 2)
  .compare_cell_scores(obj_true, obj_false)
})

test_that("scalePCs TRUE vs FALSE produce same cell scores (single slide, nCC=3)", {
  obj_true  <- .run_single_pipeline(scalePCs = TRUE,  nCC = 3)
  obj_false <- .run_single_pipeline(scalePCs = FALSE, nCC = 3)
  .compare_cell_scores(obj_true, obj_false)
})

# ===========================================================================
# Multi-slide tests
# ===========================================================================

test_that("scalePCs TRUE vs FALSE produce same cell scores (multi slide, nCC=1)", {
  obj_true  <- .run_multi_pipeline(scalePCs = TRUE,  nCC = 1)
  obj_false <- .run_multi_pipeline(scalePCs = FALSE, nCC = 1)
  .compare_cell_scores(obj_true, obj_false)
})

test_that("scalePCs TRUE vs FALSE produce same cell scores (multi slide, nCC=2)", {
  obj_true  <- .run_multi_pipeline(scalePCs = TRUE,  nCC = 2)
  obj_false <- .run_multi_pipeline(scalePCs = FALSE, nCC = 2)
  .compare_cell_scores(obj_true, obj_false)
})

test_that("scalePCs TRUE vs FALSE produce same cell scores (multi slide, nCC=3)", {
  obj_true  <- .run_multi_pipeline(scalePCs = TRUE,  nCC = 3)
  obj_false <- .run_multi_pipeline(scalePCs = FALSE, nCC = 3)
  .compare_cell_scores(obj_true, obj_false)
})
