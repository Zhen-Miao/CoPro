# Equivalence tests for the sparse fixed-radius kernel path (method = "sparse")
# against the dense path (method = "dense"). The sparse path must reproduce the
# dense result numerically while returning sparse dgCMatrix cross-kernels and
# triangularly stored dsCMatrix within-type kernels.

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
    expect_s4_class(Ks, "dsCMatrix")
    expect_true(all(Matrix::diag(Ks) == 0))
    expect_equal(2L * length(Ks@x), Matrix::nnzero(Ks))
    expect_equal(as.matrix(Ks), as.matrix(Kd), tolerance = 1e-8, ignore_attr = TRUE)
  }
})

test_that("within-type sparse storage follows symmetry-preserving normalization", {
  obj <- create_test_copro_single(n_cells = 180, n_genes = 25,
                                  n_cell_types = 2, seed = 31)
  obj <- subsetData(obj, cellTypesOfInterest = "CellTypeA")
  with_dist <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)

  for (row_normalize in c(FALSE, TRUE)) {
    dense <- computeKernelMatrix(
      with_dist, sigmaValues = 0.1, method = "dense",
      normalizeKernel = !row_normalize,
      rowNormalizeKernel = row_normalize,
      dropDistances = FALSE, verbose = FALSE
    )
    sparse <- computeKernelMatrix(
      obj, sigmaValues = 0.1, method = "sparse",
      normalizeKernel = !row_normalize,
      rowNormalizeKernel = row_normalize,
      distType = "Euclidean2D", verbose = FALSE
    )
    Kd <- getKernelMatrix(dense, 0.1, "CellTypeA", "CellTypeA",
                          verbose = FALSE)
    Ks <- getKernelMatrix(sparse, 0.1, "CellTypeA", "CellTypeA",
                          verbose = FALSE)

    expect_s4_class(Ks, if (row_normalize) "dgCMatrix" else "dsCMatrix")
    expect_equal(as.matrix(Ks), as.matrix(Kd), tolerance = 1e-8,
                 ignore_attr = TRUE)
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

test_that(".type7QuantileRepeated matches quantile() on the repeated vector", {
  # A symmetric kernel stores one triangle but its quantile is defined on the
  # represented full matrix, so the helper must reproduce exactly what
  # rep(x, each = repetitions) would have given -- without building it.
  # expect_identical, not expect_equal: the helper mirrors quantile()'s type-7
  # arithmetic statement for statement so the clip threshold cannot drift, and
  # a tolerance would hide a rearrangement that reintroduces last-bit error.
  set.seed(29)
  for (n in c(1L, 2L, 3L, 4L, 5L, 17L, 100L, 1001L)) {
    x <- runif(n)
    for (repetitions in c(1L, 2L, 3L, 4L)) {
      for (p in c(0, 1e-3, 0.1, 0.25, 1 / 3, 0.5, 0.66, 0.85, 1 - 1e-4, 1)) {
        expect_identical(
          CoPro:::.type7QuantileRepeated(x, p, repetitions),
          as.numeric(stats::quantile(rep(x, each = repetitions), p,
                                     names = FALSE)),
          info = paste("n", n, "reps", repetitions, "p", p)
        )
      }
    }
  }

  # Heavy ties: the rank arithmetic must not assume distinct values, and the
  # equal-order-statistic branch must be taken rather than interpolated.
  tied <- c(rep(0.5, 40), runif(20))
  for (repetitions in c(1L, 2L)) {
    for (p in c(0.1, 0.5, 0.85)) {
      expect_identical(
        CoPro:::.type7QuantileRepeated(tied, p, repetitions),
        as.numeric(stats::quantile(rep(tied, each = repetitions), p,
                                   names = FALSE))
      )
    }
  }

  # repetitions = 1 must reduce to the ordinary quantile, and NAs are dropped
  # the way the old na.rm = TRUE call dropped them.
  plain <- runif(500)
  expect_identical(CoPro:::.type7QuantileRepeated(plain, 0.85, 1L),
                   as.numeric(stats::quantile(plain, 0.85, names = FALSE)))
  expect_identical(CoPro:::.type7QuantileRepeated(c(plain, NA), 0.85, 1L),
                   as.numeric(stats::quantile(plain, 0.85, names = FALSE)))
  expect_error(CoPro:::.type7QuantileRepeated(numeric(0), 0.5, 1L),
               "no values")
  expect_error(CoPro:::.type7QuantileRepeated(plain, 0.5, 0L),
               "positive integer")
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

  # The public default method is "auto". At the exact threshold it must take
  # the sparse path without requiring computeDistance() first.
  threshold <- .maxCellTypeCount(obj)
  a <- computeKernelMatrix(obj, sigmaValues = c(0.1),
                           autoThreshold = threshold, distType = "Euclidean2D",
                           verbose = FALSE)
  Ka <- getKernelMatrix(a, sigma = 0.1, cellType1 = "CellTypeA",
                        cellType2 = "CellTypeB", verbose = FALSE)
  expect_s4_class(Ka, "dgCMatrix")

  # high threshold -> dense (requires distances first)
  obj_d <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
  b <- computeKernelMatrix(obj_d, sigmaValues = c(0.1),
                           autoThreshold = 1e9, dropDistances = FALSE,
                           verbose = FALSE)
  Kb <- getKernelMatrix(b, sigma = 0.1, cellType1 = "CellTypeA",
                        cellType2 = "CellTypeB", verbose = FALSE)
  expect_true(is.matrix(Kb))
})

test_that("auto accounts for aggregate blocks and per-slide dimensions", {
  obj <- create_test_copro_multi(n_cells_per_slide = 80, n_slides = 3,
                                 n_genes = 10, n_cell_types = 2, seed = 91)
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))

  # Counts across all slides exceed the largest per-slide block; auto should
  # reason about the dimensions actually materialized by each slide.
  per_slide_max <- max(vapply(getSlideList(obj), function(sID) {
    max(vapply(obj@cellTypesOfInterest, function(ct) {
      .countSlideCellType(obj, slide = sID, cellType = ct)
    }, numeric(1)))
  }, numeric(1)))
  expect_equal(.maxCellTypeCount(obj), as.integer(per_slide_max))
  expect_lt(.maxCellTypeCount(obj),
            max(as.integer(table(obj@cellTypesSub))))

  # No individual block reaches 50 cells, but the aggregate cross-slide dense
  # workload exceeds 50^2 entries, so the default auto method selects sparse.
  expect_lt(.maxCellTypeCount(obj), 50L)
  expect_gte(.denseKernelEntryCount(obj), 50^2)
  out <- computeKernelMatrix(obj, sigmaValues = 0.1,
                             autoThreshold = 50L,
                             distType = "Euclidean2D", verbose = FALSE)
  K <- getKernelMatrix(out, sigma = 0.1, cellType1 = "CellTypeA",
                       cellType2 = "CellTypeB", slide = getSlideList(out)[1],
                       verbose = FALSE)
  expect_s4_class(K, "dgCMatrix")
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

test_that("the blocked whitened-Frobenius inner product matches the direct one", {
  # <Rx K Ry, K> is evaluated as sum((Rx K) * (K Ry)) over column blocks so the
  # heavily filled-in triple product is never materialized. The identity needs
  # Rx and Ry symmetric, which .whitenedFrobNorm() guarantees, and the block
  # boundaries must not change the answer.
  set.seed(31)
  mk <- function(nr, nc, k) Matrix::sparseMatrix(
    i = rep(seq_len(nr), each = k),
    j = pmax(1L, pmin(nc, rep(seq_len(nr), each = k) +
                        sample(-8:8, nr * k, TRUE))),
    x = runif(nr * k), dims = c(nr, nc))

  for (dims in list(c(120, 120), c(90, 140), c(140, 90))) {
    nr <- dims[1]; nc <- dims[2]
    K <- mk(nr, nc, 5)
    Rx <- local({ R <- mk(nr, nr, 5); (R + Matrix::t(R)) / 2 })
    Ry <- local({ R <- mk(nc, nc, 5); (R + Matrix::t(R)) / 2 })
    direct <- as.numeric(sum(((Rx %*% K) %*% Ry) * K))

    # Several block sizes, including one column at a time and one single block.
    for (budget in c(1, 50, 500, 1e9)) {
      expect_equal(
        CoPro:::.sparseWhitenedInner(K, Rx, Ry, block_nnz = budget),
        direct, tolerance = 1e-10,
        info = paste(nr, "x", nc, "block_nnz", budget)
      )
    }
  }
})

test_that(".sparseSumSquares counts represented entries for both storage forms", {
  set.seed(37)
  general <- Matrix::sparseMatrix(
    i = sample(60, 200, TRUE), j = sample(45, 200, TRUE),
    x = runif(200), dims = c(60, 45))
  expect_equal(CoPro:::.sparseSumSquares(general), sum(general * general))
  expect_identical(CoPro:::.sparseSumSquares(general), sum(general@x^2))

  # A dsCMatrix stores one triangle; every stored off-diagonal entry stands for
  # two represented entries, so reading @x alone would undercount.
  symmetric <- Matrix::forceSymmetric(Matrix::crossprod(general))
  expect_s4_class(symmetric, "symmetricMatrix")
  expect_equal(CoPro:::.sparseSumSquares(symmetric),
               sum(symmetric * symmetric))
  expect_gt(CoPro:::.sparseSumSquares(symmetric), sum(symmetric@x^2))
})
