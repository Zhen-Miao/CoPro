# Equivalence tests for the sparse fixed-radius kernel path (method = "sparse")
# against the dense path (method = "dense"). The sparse path must reproduce the
# dense result numerically while returning sparse dgCMatrix kernels.

# Flip the sign of `x` (per column vector) to best match a reference, so weight
# vectors that are equal up to sign compare as equal.
.align_sign <- function(ref, x) {
  if (sum((x - ref)^2) <= sum((x + ref)^2)) x else -x
}

# Compare a sparse kernel against the dense kernel for one (sigma, pair).
.expect_kernel_equal <- function(dense_obj, sparse_obj, sigma, ct1, ct2,
                                 slide = NULL, tol = 1e-8) {
  Kd <- getKernelMatrix(dense_obj, sigma = sigma, cellType1 = ct1,
                        cellType2 = ct2, slide = slide, verbose = FALSE)
  Ks <- getKernelMatrix(sparse_obj, sigma = sigma, cellType1 = ct1,
                        cellType2 = ct2, slide = slide, verbose = FALSE)
  expect_s4_class(Ks, "dgCMatrix")
  expect_equal(dim(Ks), dim(Kd))
  expect_equal(as.matrix(Ks), as.matrix(Kd), tolerance = tol, ignore_attr = TRUE)
}

test_that("sparse kernels match dense kernels entrywise (single-slide, pairwise)", {
  obj <- create_test_copro_single(n_cells = 300, n_genes = 40,
                                  n_cell_types = 2, seed = 1)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))

  # normalizeDistance = TRUE rescales distances to ~0.01 units (small sigmas);
  # FALSE keeps the raw 0-10 coordinate scale (so use larger sigmas).
  sigma_sets <- list(`TRUE` = c(0.05, 0.1, 0.2), `FALSE` = c(0.5, 1, 2))

  for (nd in c(TRUE, FALSE)) {
    sigmas <- sigma_sets[[as.character(nd)]]
    dense <- computeDistance(obj, distType = "Euclidean2D",
                             normalizeDistance = nd, verbose = FALSE)
    dense <- computeKernelMatrix(dense, sigmaValues = sigmas, method = "dense",
                                 normalizeDistance = nd, dropDistances = FALSE,
                                 verbose = FALSE)
    sparse <- computeKernelMatrix(obj, sigmaValues = sigmas, method = "sparse",
                                  distType = "Euclidean2D",
                                  normalizeDistance = nd, verbose = FALSE)

    expect_equal(sort(dense@sigmaValues), sort(sparse@sigmaValues))
    for (s in intersect(dense@sigmaValues, sparse@sigmaValues)) {
      .expect_kernel_equal(dense, sparse, s, "CellTypeA", "CellTypeB")
    }
  }
})

test_that("sparse kernels match dense kernels for all normalization options", {
  obj <- create_test_copro_single(n_cells = 250, n_genes = 30,
                                  n_cell_types = 2, seed = 2)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  sigmas <- c(0.1)

  norm_opts <- list(
    list(normalizeKernel = TRUE,  rowNormalizeKernel = FALSE, colNormalizeKernel = FALSE),
    list(normalizeKernel = FALSE, rowNormalizeKernel = TRUE,  colNormalizeKernel = FALSE),
    list(normalizeKernel = FALSE, rowNormalizeKernel = FALSE, colNormalizeKernel = TRUE)
  )

  for (opt in norm_opts) {
    dense <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
    dense <- do.call(computeKernelMatrix, c(
      list(dense, sigmaValues = sigmas, method = "dense",
           dropDistances = FALSE, verbose = FALSE), opt))
    sparse <- do.call(computeKernelMatrix, c(
      list(obj, sigmaValues = sigmas, method = "sparse",
           distType = "Euclidean2D", verbose = FALSE), opt))
    .expect_kernel_equal(dense, sparse, 0.1, "CellTypeA", "CellTypeB")
  }
})

test_that("sparse kernels match dense kernels for a single (within-type) cell type", {
  obj <- create_test_copro_single(n_cells = 200, n_genes = 30,
                                  n_cell_types = 2, seed = 3)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA"))
  sigmas <- c(0.05, 0.1)

  dense <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
  dense <- computeKernelMatrix(dense, sigmaValues = sigmas, method = "dense",
                               dropDistances = FALSE, verbose = FALSE)
  sparse <- computeKernelMatrix(obj, sigmaValues = sigmas, method = "sparse",
                                distType = "Euclidean2D", verbose = FALSE)

  for (s in intersect(dense@sigmaValues, sparse@sigmaValues)) {
    Ks <- getKernelMatrix(sparse, sigma = s, cellType1 = "CellTypeA",
                          cellType2 = "CellTypeA", verbose = FALSE)
    Kd <- getKernelMatrix(dense, sigma = s, cellType1 = "CellTypeA",
                          cellType2 = "CellTypeA", verbose = FALSE)
    # within-type diagonal must be exactly zero (mirrors diag = Inf -> K = 0)
    expect_true(all(Matrix::diag(Ks) == 0))
    expect_equal(as.matrix(Ks), as.matrix(Kd), tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("sparse kernels match dense kernels for multi-slide objects", {
  obj <- create_test_copro_multi(n_cells_per_slide = 150, n_slides = 2,
                                 n_genes = 30, n_cell_types = 2, seed = 5)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  sigmas <- c(0.1)

  dense <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
  dense <- computeKernelMatrix(dense, sigmaValues = sigmas, method = "dense",
                               dropDistances = FALSE, verbose = FALSE)
  sparse <- computeKernelMatrix(obj, sigmaValues = sigmas, method = "sparse",
                                distType = "Euclidean2D", verbose = FALSE)

  for (sID in getSlideList(obj)) {
    .expect_kernel_equal(dense, sparse, 0.1, "CellTypeA", "CellTypeB", slide = sID)
  }
})

test_that("sparse path actually produces a sparse kernel (large extent, small sigma)", {
  obj <- create_test_copro_single(n_cells = 400, n_genes = 20,
                                  n_cell_types = 2, seed = 6)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  # Original-unit coordinates span 0-10; sigma = 0.3 keeps support local.
  sparse <- computeKernelMatrix(obj, sigmaValues = c(0.3), method = "sparse",
                                distType = "Euclidean2D",
                                normalizeDistance = FALSE, verbose = FALSE)
  K <- getKernelMatrix(sparse, sigma = 0.3, cellType1 = "CellTypeA",
                       cellType2 = "CellTypeB", verbose = FALSE)
  expect_s4_class(K, "dgCMatrix")
  expect_lt(length(K@x), prod(dim(K)))  # genuinely sparse: fewer nonzeros than n_i*n_j
})

test_that("sparse path yields equivalent skrCCA weights and normalized correlation", {
  obj <- create_test_copro_single(n_cells = 300, n_genes = 40,
                                  n_cell_types = 2, seed = 7)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computePCA(obj, nPCA = 10, center = TRUE, scale. = TRUE)
  sigmas <- c(0.05, 0.1)

  dense <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
  dense <- computeKernelMatrix(dense, sigmaValues = sigmas, method = "dense",
                               dropDistances = FALSE, verbose = FALSE)
  dense <- runSkrCCA(dense, scalePCs = TRUE, nCC = 2, maxIter = 300, tol = 1e-7)
  dense <- computeNormalizedCorrelation(dense)

  sparse <- computeKernelMatrix(obj, sigmaValues = sigmas, method = "sparse",
                                distType = "Euclidean2D", verbose = FALSE)
  sparse <- runSkrCCA(sparse, scalePCs = TRUE, nCC = 2, maxIter = 300, tol = 1e-7)
  sparse <- computeNormalizedCorrelation(sparse)

  surviving <- intersect(dense@sigmaValues, sparse@sigmaValues)
  expect_gt(length(surviving), 0)
  for (s in surviving) {
    sn <- paste0("sigma_", s)
    for (ct in c("CellTypeA", "CellTypeB")) {
      Wd <- dense@skrCCAOut[[sn]][[ct]]
      Ws <- sparse@skrCCAOut[[sn]][[ct]]
      expect_equal(dim(Wd), dim(Ws))
      for (cc in seq_len(ncol(Wd))) {
        # values must match up to sign; weight rownames are cosmetic and the
        # Matrix package drops them in sparse products, so ignore attributes
        expect_equal(unname(.align_sign(Wd[, cc], Ws[, cc])),
                     unname(Wd[, cc]), tolerance = 1e-3)
      }
    }
  }

  # Normalized-correlation magnitudes match (sign is weight-orientation dependent)
  nd <- dense@normalizedCorrelation[[1]]
  ns <- sparse@normalizedCorrelation[[1]]
  expect_equal(sort(abs(nd$normalizedCorrelation)),
               sort(abs(ns$normalizedCorrelation)), tolerance = 1e-3)
})

test_that(".exactLowQuantile reproduces R type-7 quantile from the smallest values", {
  set.seed(11)
  x <- runif(2000)
  k <- 200L
  sorted_small <- sort(x)[seq_len(k)]
  for (p in c(1e-3, 5e-3, 1e-2, 0.05)) {
    expect_equal(
      CoPro:::.exactLowQuantile(sorted_small, length(x), p),
      as.numeric(stats::quantile(x, p)),
      tolerance = 1e-12
    )
  }
})

test_that("sparse-safe kernel normalizations match the dense formulas", {
  set.seed(17)
  Kfull <- matrix(0, nrow = 40, ncol = 35)
  idx <- sample(length(Kfull), 250)
  Kfull[idx] <- runif(250)
  Ksp <- as(Matrix::Matrix(Kfull, sparse = TRUE), "CsparseMatrix")

  # Sinkhorn-Knopp scaling
  expect_equal(as.matrix(CoPro:::sinkhorn_knopp(Ksp)),
               CoPro:::sinkhorn_knopp(Kfull), tolerance = 1e-8, ignore_attr = TRUE)

  # Bidirectional cross-correlation for every normalization mode
  A <- matrix(rnorm(40), ncol = 1)
  B <- matrix(rnorm(35), ncol = 1)
  for (nk in c("none", "row_or_col", "sinkhorn_knopp")) {
    cd <- CoPro:::.computeSpatialCrossCorrelation(A, B, Kfull, normalize_K = nk,
                                                 filter_kernel = FALSE)
    cs <- CoPro:::.computeSpatialCrossCorrelation(A, B, Ksp, normalize_K = nk,
                                                 filter_kernel = FALSE)
    expect_equal(cs, cd, tolerance = 1e-8)
  }
})

test_that("dropDistances default clears distances; FALSE retains them", {
  obj <- create_test_copro_single(n_cells = 150, n_genes = 20,
                                  n_cell_types = 2, seed = 8)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)

  dropped <- computeKernelMatrix(obj, sigmaValues = c(0.1), method = "dense",
                                 verbose = FALSE)  # dropDistances = TRUE default
  expect_equal(length(dropped@distances), 0)

  kept <- computeKernelMatrix(obj, sigmaValues = c(0.1), method = "dense",
                              dropDistances = FALSE, verbose = FALSE)
  expect_gt(length(kept@distances), 0)

  # sparse path leaves distances empty too (it never builds them)
  sp <- computeKernelMatrix(obj, sigmaValues = c(0.1), method = "sparse",
                            distType = "Euclidean2D", verbose = FALSE)
  expect_equal(length(sp@distances), 0)
})

test_that("method = 'auto' selects sparse above threshold and dense below", {
  obj <- create_test_copro_single(n_cells = 150, n_genes = 20,
                                  n_cell_types = 2, seed = 9)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))

  # low threshold -> sparse
  a <- computeKernelMatrix(obj, sigmaValues = c(0.1), method = "auto",
                           autoThreshold = 5, distType = "Euclidean2D",
                           verbose = FALSE)
  Ka <- getKernelMatrix(a, sigma = 0.1, cellType1 = "CellTypeA",
                        cellType2 = "CellTypeB", verbose = FALSE)
  expect_s4_class(Ka, "dgCMatrix")

  # high threshold -> dense (requires distances first)
  obj_d <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
  b <- computeKernelMatrix(obj_d, sigmaValues = c(0.1), method = "auto",
                           autoThreshold = 1e9, dropDistances = FALSE,
                           verbose = FALSE)
  Kb <- getKernelMatrix(b, sigma = 0.1, cellType1 = "CellTypeA",
                        cellType2 = "CellTypeB", verbose = FALSE)
  expect_true(is.matrix(Kb))
})

test_that("sparse path does not require computeDistance and rejects Morphology-Aware", {
  obj <- create_test_copro_single(n_cells = 120, n_genes = 20,
                                  n_cell_types = 2, seed = 10)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))

  # No computeDistance() call beforehand -> still works
  sp <- computeKernelMatrix(obj, sigmaValues = c(0.1), method = "sparse",
                            distType = "Euclidean2D", verbose = FALSE)
  expect_gt(length(sp@kernelMatrices), 0)

  expect_error(
    computeKernelMatrix(obj, sigmaValues = c(0.1), method = "sparse",
                        distType = "Morphology-Aware", verbose = FALSE),
    "Euclidean"
  )
})

test_that("sparse whitened-Frobenius normalizer matches dense calculation", {
  set.seed(701)
  K <- Matrix::rsparsematrix(37, 41, density = 0.18)
  X <- Matrix::rsparsematrix(37, 37, density = 0.12)
  Y <- Matrix::rsparsematrix(41, 41, density = 0.12)
  Rx <- Matrix::crossprod(X) + Matrix::Diagonal(37)
  Ry <- Matrix::crossprod(Y) + Matrix::Diagonal(41)

  sparse_norm <- .whitenedFrobNorm(K, Rx, Ry)
  dense_norm <- .whitenedFrobNorm(as.matrix(K), as.matrix(Rx), as.matrix(Ry))
  expect_equal(sparse_norm, dense_norm, tolerance = 1e-10)

  sparse_unwhitened <- .whitenedFrobNorm(K)
  dense_unwhitened <- .whitenedFrobNorm(as.matrix(K))
  expect_equal(sparse_unwhitened, dense_unwhitened, tolerance = 1e-10)
})
