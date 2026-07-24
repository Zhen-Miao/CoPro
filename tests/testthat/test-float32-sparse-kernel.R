.align_float32_component <- function(reference, candidate) {
  if (sum(reference * candidate) < 0) -candidate else candidate
}

test_that("direct float32 kernels match float64 sparse kernels", {
  object <- create_test_copro_single(
    n_cells = 500, n_genes = 40, n_cell_types = 2, seed = 81
  )
  object <- subsetData(
    object, cellTypesOfInterest = c("CellTypeA", "CellTypeB")
  )
  sigmas <- c(0.05, 0.1, 0.2)

  double <- computeKernelMatrix(
    object, sigmaValues = sigmas, method = "sparse",
    distType = "Euclidean2D", verbose = FALSE
  )
  float32 <- computeSparseKernelFloat32(
    object, sigmaValues = sigmas, distType = "Euclidean2D",
    verbose = FALSE
  )

  expect_length(float32@distances, 0L)
  expect_equal(float32@sigmaValues, double@sigmaValues)
  expect_lt(
    as.numeric(object.size(float32@kernelMatrices)),
    as.numeric(object.size(double@kernelMatrices))
  )

  for (sigma in sigmas) {
    encoded <- getKernelMatrix(
      float32, sigma, "CellTypeA", "CellTypeB", verbose = FALSE,
      materialize = FALSE
    )
    reference <- getKernelMatrix(
      double, sigma, "CellTypeA", "CellTypeB", verbose = FALSE
    )
    expect_s3_class(encoded, "CoProFloat32SparseMatrix")
    expect_equal(dim(encoded), dim(reference))
    expect_equal(
      as.matrix(asDoubleSparseMatrix(encoded)),
      as.matrix(reference),
      tolerance = 2e-6,
      ignore_attr = TRUE
    )

    reverse <- getKernelMatrix(
      float32, sigma, "CellTypeB", "CellTypeA", verbose = FALSE,
      materialize = FALSE
    )
    expect_equal(dim(reverse), rev(dim(encoded)))
    expect_true(isTRUE(reverse$transposed))
    expect_equal(
      as.matrix(asDoubleSparseMatrix(reverse)),
      t(as.matrix(reference)),
      tolerance = 2e-6,
      ignore_attr = TRUE
    )

    compatible <- getKernelMatrix(
      float32, sigma, "CellTypeA", "CellTypeB", verbose = FALSE
    )
    expect_s4_class(compatible, "dgCMatrix")
    expect_equal(
      as.matrix(compatible), as.matrix(reference),
      tolerance = 2e-6, ignore_attr = TRUE
    )
  }

  fully_materialized <- materializeFloat32Kernels(
    float32, verbose = FALSE
  )
  expect_true(all(vapply(
    fully_materialized@kernelMatrices,
    inherits, logical(1), what = "sparseMatrix"
  )))
  expect_false(any(vapply(
    fully_materialized@kernelMatrices,
    .isFloat32SparseKernel, logical(1)
  )))
})

test_that("parallel float32 operators agree with decoded sparse operations", {
  set.seed(82)
  coordinate_1 <- matrix(runif(600), ncol = 2)
  coordinate_2 <- matrix(runif(420), ncol = 2)
  built <- float32_csr_gaussian_kernels_cpp(
    coordinate_1, coordinate_2, c(0.05, 0.1),
    percentile = 0.005, scaling_factor = 1,
    lower_limit = 1e-7, upper_quantile = 0.85
  )
  encoded <- .newFloat32SparseKernel(built$kernels[[2L]], NULL)
  decoded <- asDoubleSparseMatrix(encoded)
  x_left <- matrix(rnorm(nrow(coordinate_1) * 8), ncol = 8)
  x_right <- matrix(rnorm(nrow(coordinate_2) * 7), ncol = 7)
  expected_y <- as.matrix(crossprod(x_left, decoded %*% x_right))

  for (threads in c(1L, 2L, 4L)) {
    observed_y <- .kernelXKY(
      x_left, encoded, x_right, n_threads = threads
    )
    expect_equal(
      observed_y, expected_y,
      tolerance = 2e-5, ignore_attr = TRUE
    )
  }

  expect_equal(
    .float32KernelMatMult(encoded, x_right, n_threads = 2L),
    as.matrix(decoded %*% x_right),
    tolerance = 2e-5,
    ignore_attr = TRUE
  )
  expect_equal(
    .float32KernelMatMult(t(encoded), x_left, n_threads = 2L),
    as.matrix(t(decoded) %*% x_left),
    tolerance = 2e-5,
    ignore_attr = TRUE
  )
})

test_that("float32 construction supports global, row, and column normalization", {
  object <- create_test_copro_single(
    n_cells = 420, n_genes = 25, n_cell_types = 2, seed = 821
  )
  object <- subsetData(
    object, cellTypesOfInterest = c("CellTypeA", "CellTypeB")
  )
  modes <- list(
    global = list(
      normalizeKernel = TRUE,
      rowNormalizeKernel = FALSE,
      colNormalizeKernel = FALSE
    ),
    row = list(
      normalizeKernel = FALSE,
      rowNormalizeKernel = TRUE,
      colNormalizeKernel = FALSE
    ),
    column = list(
      normalizeKernel = FALSE,
      rowNormalizeKernel = FALSE,
      colNormalizeKernel = TRUE
    )
  )

  for (mode_name in names(modes)) {
    arguments <- modes[[mode_name]]
    reference <- do.call(
      computeKernelMatrix,
      c(
        list(
          object = object, sigmaValues = 0.1, method = "sparse",
          distType = "Euclidean2D", verbose = FALSE
        ),
        arguments
      )
    )
    candidate <- do.call(
      computeSparseKernelFloat32,
      c(
        list(
          object = object, sigmaValues = 0.1,
          distType = "Euclidean2D", verbose = FALSE
        ),
        arguments
      )
    )
    reference_kernel <- getKernelMatrix(
      reference, 0.1, "CellTypeA", "CellTypeB", verbose = FALSE
    )
    encoded <- getKernelMatrix(
      candidate, 0.1, "CellTypeA", "CellTypeB",
      verbose = FALSE, materialize = FALSE
    )
    decoded <- asDoubleSparseMatrix(encoded)
    expect_equal(
      as.matrix(decoded), as.matrix(reference_kernel),
      tolerance = 2e-6, ignore_attr = TRUE,
      info = mode_name
    )
    expect_equal(
      .whitenedFrobNorm(encoded),
      .whitenedFrobNorm(decoded),
      tolerance = 1e-10,
      info = paste(mode_name, "Frobenius norm")
    )
    if (mode_name == "row") {
      positive <- Matrix::rowSums(decoded) > 1e-4
      expect_equal(
        as.numeric(Matrix::rowSums(decoded)[positive]),
        rep.int(1, sum(positive)),
        tolerance = 2e-6
      )
    } else if (mode_name == "column") {
      positive <- Matrix::colSums(decoded) > 1e-4
      expect_equal(
        as.numeric(Matrix::colSums(decoded)[positive]),
        rep.int(1, sum(positive)),
        tolerance = 2e-6
      )
    }
  }

  expect_error(
    computeSparseKernelFloat32(
      object, sigmaValues = 0.1,
      rowNormalizeKernel = TRUE, colNormalizeKernel = TRUE,
      verbose = FALSE
    ),
    "Cannot do both"
  )
})

test_that("normalized float32 self-kernels expand only when asymmetric", {
  object <- create_test_copro_single(
    n_cells = 360, n_genes = 25, n_cell_types = 2, seed = 822
  )
  object <- subsetData(object, cellTypesOfInterest = "CellTypeA")

  for (mode_name in c("row", "column")) {
    arguments <- list(
      rowNormalizeKernel = mode_name == "row",
      colNormalizeKernel = mode_name == "column"
    )
    reference <- do.call(
      computeKernelMatrix,
      c(
        list(
          object = object, sigmaValues = 0.1, method = "sparse",
          distType = "Euclidean2D", verbose = FALSE
        ),
        arguments
      )
    )
    candidate <- do.call(
      computeSparseKernelFloat32,
      c(
        list(
          object = object, sigmaValues = 0.1,
          distType = "Euclidean2D", verbose = FALSE
        ),
        arguments
      )
    )
    encoded <- getKernelMatrix(
      candidate, 0.1, "CellTypeA", "CellTypeA",
      verbose = FALSE, materialize = FALSE
    )
    reference_kernel <- getKernelMatrix(
      reference, 0.1, "CellTypeA", "CellTypeA", verbose = FALSE
    )
    expect_false(isTRUE(encoded$symmetric))
    expect_equal(
      as.matrix(asDoubleSparseMatrix(encoded)),
      as.matrix(reference_kernel),
      tolerance = 2e-6, ignore_attr = TRUE,
      info = mode_name
    )
  }

  globally_normalized <- computeSparseKernelFloat32(
    object, sigmaValues = 0.1, normalizeKernel = TRUE,
    distType = "Euclidean2D", verbose = FALSE
  )
  encoded_global <- getKernelMatrix(
    globally_normalized, 0.1, "CellTypeA", "CellTypeA",
    verbose = FALSE, materialize = FALSE
  )
  expect_true(isTRUE(encoded_global$symmetric))
})

test_that("float32 Frobenius norms match decoded sparse kernels", {
  set.seed(823)
  coordinate_1 <- matrix(runif(320), ncol = 2)
  coordinate_2 <- matrix(runif(280), ncol = 2)
  cross_built <- float32_csr_gaussian_kernels_cpp(
    coordinate_1, coordinate_2, 0.1,
    percentile = 0.005, scaling_factor = 1,
    lower_limit = 1e-7, upper_quantile = 0.85
  )
  self_1_built <- float32_csr_gaussian_kernels_cpp(
    coordinate_1, coordinate_1, 0.1,
    percentile = 0.005, scaling_factor = 1,
    lower_limit = 1e-7, upper_quantile = 0.85,
    symmetric = TRUE
  )
  self_2_built <- float32_csr_gaussian_kernels_cpp(
    coordinate_2, coordinate_2, 0.1,
    percentile = 0.005, scaling_factor = 1,
    lower_limit = 1e-7, upper_quantile = 0.85,
    symmetric = TRUE
  )
  cross <- .newFloat32SparseKernel(cross_built$kernels[[1L]], NULL)
  self_1 <- .newFloat32SparseKernel(self_1_built$kernels[[1L]], NULL)
  self_2 <- .newFloat32SparseKernel(self_2_built$kernels[[1L]], NULL)
  cross_double <- asDoubleSparseMatrix(cross)
  self_1_double <- asDoubleSparseMatrix(self_1)
  self_2_double <- asDoubleSparseMatrix(self_2)

  expect_equal(
    .whitenedFrobNorm(cross),
    .whitenedFrobNorm(cross_double),
    tolerance = 1e-10
  )
  expect_equal(
    .whitenedFrobNorm(self_1),
    .whitenedFrobNorm(self_1_double),
    tolerance = 1e-10
  )
  # Native float32 whitening is intentionally deferred; this verifies the
  # exact compatibility fallback through temporary double sparse matrices.
  expect_equal(
    .whitenedFrobNorm(cross, self_1, self_2),
    .whitenedFrobNorm(cross_double, self_1_double, self_2_double),
    tolerance = 1e-10
  )
})

test_that("float32 kernels run through skrCCA with stable cell scores", {
  object <- create_test_copro_single(
    n_cells = 600, n_genes = 50, n_cell_types = 2, seed = 83
  )
  object <- subsetData(
    object, cellTypesOfInterest = c("CellTypeA", "CellTypeB")
  )
  object <- computePCA(
    object, nPCA = 12, center = TRUE, scale. = TRUE
  )

  double <- computeKernelMatrix(
    object, sigmaValues = c(0.05, 0.1), method = "sparse",
    distType = "Euclidean2D", verbose = FALSE
  )
  float32 <- computeSparseKernelFloat32(
    object, sigmaValues = c(0.05, 0.1),
    distType = "Euclidean2D", verbose = FALSE
  )

  set.seed(8401)
  double <- runSkrCCA(
    double, scalePCs = TRUE, nCC = 3, maxIter = 300, tol = 1e-7
  )
  set.seed(8401)
  float32 <- runSkrCCA(
    float32, scalePCs = TRUE, nCC = 3, maxIter = 300, tol = 1e-7
  )
  double <- computeNormalizedCorrelation(double)
  float32 <- computeNormalizedCorrelation(float32)
  expect_equal(
    getNormCorr(float32)$normalizedCorrelation,
    getNormCorr(double)$normalizedCorrelation,
    tolerance = 2e-4
  )
  double <- computeGeneAndCellScores(double)
  float32 <- computeGeneAndCellScores(float32)

  for (sigma in intersect(double@sigmaValues, float32@sigmaValues)) {
    for (cell_type in c("CellTypeA", "CellTypeB")) {
      reference <- getCellScores(
        double, sigma, cell_type, verbose = FALSE
      )
      candidate <- getCellScores(
        float32, sigma, cell_type, verbose = FALSE
      )
      for (component in seq_len(ncol(reference))) {
        aligned <- .align_float32_component(
          reference[, component], candidate[, component]
        )
        expect_gt(stats::cor(reference[, component], aligned), 0.99999)
        expect_lt(
          sqrt(mean((reference[, component] - aligned)^2)) /
            stats::sd(reference[, component]),
          5e-3
        )
      }
    }
  }
})

test_that("float32 one-type symmetric kernels match float64 and run end to end", {
  object <- create_test_copro_single(
    n_cells = 400, n_genes = 35, n_cell_types = 2, seed = 84
  )
  object <- subsetData(object, cellTypesOfInterest = "CellTypeA")
  object <- computePCA(
    object, nPCA = 8, center = TRUE, scale. = TRUE
  )
  double <- computeKernelMatrix(
    object, sigmaValues = 0.1, method = "sparse",
    distType = "Euclidean2D", verbose = FALSE
  )
  float32 <- computeSparseKernelFloat32(
    object, sigmaValues = 0.1, distType = "Euclidean2D",
    verbose = FALSE
  )

  encoded <- getKernelMatrix(
    float32, 0.1, "CellTypeA", "CellTypeA",
    verbose = FALSE, materialize = FALSE
  )
  reference <- getKernelMatrix(
    double, 0.1, "CellTypeA", "CellTypeA", verbose = FALSE
  )
  expect_true(isTRUE(encoded$symmetric))
  expect_equal(.float32KernelNnz(encoded), length(reference@x))
  expect_equal(length(encoded$x), 4L * length(reference@x))
  expect_equal(
    as.matrix(asDoubleSparseMatrix(encoded)), as.matrix(reference),
    tolerance = 2e-6, ignore_attr = TRUE
  )
  set.seed(8402)
  rhs <- matrix(rnorm(nrow(reference) * 4L), ncol = 4L)
  expect_equal(
    .float32KernelMatMult(encoded, rhs, n_threads = 2L),
    as.matrix(reference %*% rhs),
    tolerance = 2e-5, ignore_attr = TRUE
  )

  set.seed(8501)
  double <- runSkrCCA(
    double, scalePCs = TRUE, nCC = 2, maxIter = 200, tol = 1e-7
  )
  set.seed(8501)
  float32 <- runSkrCCA(
    float32, scalePCs = TRUE, nCC = 2, maxIter = 200, tol = 1e-7
  )
  double <- computeGeneAndCellScores(double)
  float32 <- computeGeneAndCellScores(float32)
  reference_scores <- getCellScores(
    double, 0.1, "CellTypeA", verbose = FALSE
  )
  candidate_scores <- getCellScores(
    float32, 0.1, "CellTypeA", verbose = FALSE
  )
  for (component in seq_len(ncol(reference_scores))) {
    aligned <- .align_float32_component(
      reference_scores[, component], candidate_scores[, component]
    )
    expect_gt(stats::cor(reference_scores[, component], aligned), 0.99999)
  }

  float32 <- computeNormalizedCorrelation(float32)
  expect_no_error({
    one_type_plot_data <- getCorrOneType(
      float32, "CellTypeA", ccIndex = 1,
      sigmaValueChoice = 0.1
    )
  })
  expect_s3_class(one_type_plot_data, "data.frame")
})

test_that("float32 supports more than two cell types", {
  object <- create_test_copro_single(
    n_cells = 600, n_genes = 40, n_cell_types = 3, seed = 85
  )
  object <- subsetData(
    object,
    cellTypesOfInterest = c("CellTypeA", "CellTypeB", "CellTypeC")
  )
  object <- computePCA(
    object, nPCA = 8, center = TRUE, scale. = TRUE
  )
  float32 <- computeSparseKernelFloat32(
    object, sigmaValues = 0.1, distType = "Euclidean2D",
    verbose = FALSE
  )

  expect_length(float32@kernelMatrices, 3L)
  expect_true(all(vapply(
    float32@kernelMatrices, .isFloat32SparseKernel, logical(1)
  )))
  expect_no_error({
    fitted <- runSkrCCA(
      float32, scalePCs = TRUE, nCC = 2,
      maxIter = 300, tol = 1e-6
    )
    fitted <- computeNormalizedCorrelation(fitted)
    fitted <- computeGeneAndCellScores(fitted)
  })
  expect_length(fitted@cellScores, 3L)
  expect_no_error({
    plot_data <- getCorrTwoTypes(
      fitted, "CellTypeA", "CellTypeB", ccIndex = 1,
      sigmaValueChoice = 0.1
    )
  })
  expect_s3_class(plot_data, "data.frame")

  transferred_scores <- stats::setNames(lapply(
    fitted@cellTypesOfInterest,
    function(cell_type) {
      getCellScores(
        fitted, 0.1, cell_type, verbose = FALSE
      )
    }
  ), fitted@cellTypesOfInterest)
  expect_no_error({
    transfer_result <- getTransferNormCorr(
      fitted, transferred_scores, sigma_choice = 0.1,
      verbose = FALSE
    )
  })
  expect_s3_class(transfer_result[[1L]], "data.frame")
})

test_that("float32 supports multi-slide pairwise and one-type kernels", {
  object <- create_test_copro_multi(
    n_cells_per_slide = 180, n_slides = 2,
    n_genes = 35, n_cell_types = 3, seed = 86
  )
  object <- subsetData(
    object,
    cellTypesOfInterest = c("CellTypeA", "CellTypeB", "CellTypeC")
  )
  double <- computeKernelMatrix(
    object, sigmaValues = 0.1, method = "sparse",
    distType = "Euclidean2D", verbose = FALSE
  )
  float32 <- computeSparseKernelFloat32(
    object, sigmaValues = 0.1, distType = "Euclidean2D",
    verbose = FALSE
  )
  expect_equal(
    sort(names(float32@kernelMatrices)),
    sort(names(double@kernelMatrices))
  )
  for (kernel_name in names(float32@kernelMatrices)) {
    parsed <- .parseKernelMatrixName(kernel_name)
    reference <- getKernelMatrix(
      double, parsed$sigma, parsed$cellType1, parsed$cellType2,
      slide = parsed$slide, verbose = FALSE
    )
    candidate <- getKernelMatrix(
      float32, parsed$sigma, parsed$cellType1, parsed$cellType2,
      slide = parsed$slide, verbose = FALSE
    )
    expect_equal(
      as.matrix(candidate), as.matrix(reference),
      tolerance = 2e-6, ignore_attr = TRUE
    )
  }

  float32 <- computePCA(
    float32, nPCA = 8, center = TRUE, scale. = TRUE
  )
  expect_no_error({
    fitted <- runSkrCCA(
      float32, scalePCs = TRUE, nCC = 2,
      maxIter = 300, tol = 1e-6
    )
    fitted <- computeNormalizedCorrelation(
      fitted, calculationMode = "aggregate"
    )
    fitted <- computeGeneAndCellScores(fitted)
  })
  expect_length(fitted@cellScores, 3L)

  one_type <- create_test_copro_multi(
    n_cells_per_slide = 180, n_slides = 2,
    n_genes = 35, n_cell_types = 3, seed = 87
  )
  one_type <- subsetData(
    one_type, cellTypesOfInterest = "CellTypeA"
  )
  one_type <- computeSparseKernelFloat32(
    one_type, sigmaValues = 0.1, distType = "Euclidean2D",
    verbose = FALSE
  )
  expect_length(one_type@kernelMatrices, 2L)
  expect_true(all(vapply(
    one_type@kernelMatrices,
    function(kernel) isTRUE(kernel$symmetric),
    logical(1)
  )))
  one_type <- computePCA(
    one_type, nPCA = 8, center = TRUE, scale. = TRUE
  )
  expect_no_error({
    one_type_fit <- runSkrCCA(
      one_type, scalePCs = TRUE, nCC = 2,
      maxIter = 200, tol = 1e-6
    )
    one_type_fit <- computeNormalizedCorrelation(
      one_type_fit, calculationMode = "aggregate"
    )
    one_type_fit <- computeGeneAndCellScores(one_type_fit)
  })
  expect_length(one_type_fit@cellScores, 1L)
})
