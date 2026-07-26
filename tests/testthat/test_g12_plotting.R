test_that("plotG12Functions works with valid input", {
  # Skip if required packages not available
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")
  skip_if_not_installed("spatstat.random")
  
  # Create synthetic test data
  set.seed(123)
  n_cells <- 60
  cell_ids <- paste0('Cell_', 1:n_cells)
  
  coords <- data.frame(
    x = runif(n_cells, 0, 100),
    y = runif(n_cells, 0, 100)
  )
  rownames(coords) <- cell_ids
  
  cell_types <- sample(c('TypeA', 'TypeB', 'TypeC'), n_cells, replace = TRUE)
  
  meta_data <- data.frame(cellID = cell_ids)
  rownames(meta_data) <- cell_ids
  
  # Create expression data
  n_genes <- 20
  expr_data <- matrix(rnorm(n_cells * n_genes), nrow = n_cells, ncol = n_genes)
  rownames(expr_data) <- cell_ids
  colnames(expr_data) <- paste0('Gene_', 1:n_genes)
  
  # Create test object
  test_obj <- newCoProSingle(
    normalizedData = expr_data,
    locationData = coords,
    metaData = meta_data,
    cellTypes = cell_types
  )
  
  # Subset data to populate locationDataSub and cellTypesSub
  test_obj <- subsetData(test_obj, cellTypesOfInterest = c('TypeA', 'TypeB', 'TypeC'))
  
  # Test basic functionality
  expect_no_error({
    result <- plotG12Functions(
      test_obj,
      r_um_range = c(5, 25),
      nsim = 10,  # Use fewer simulations for speed
      verbose = FALSE
    )
  })
  
  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("plot", "data", "summary"))

  # result$plot is a stable list(combined, individual); combined is present
  # by default since plot_type defaults to "combined".
  expect_type(result$plot, "list")
  expect_named(result$plot, c("combined", "individual"))
  expect_s3_class(result$plot$combined, "ggplot")
  expect_null(result$plot$individual)
  
  # Check data structure
  expect_s3_class(result$data, "data.frame")
  expect_true("r_um" %in% colnames(result$data))
  expect_true("g12_obs" %in% colnames(result$data))
  expect_true("pair_label" %in% colnames(result$data))
  
  # Check summary structure
  expect_s3_class(result$summary, "data.frame")
  expect_true("mean_g12" %in% colnames(result$summary))
  expect_true("max_g12" %in% colnames(result$summary))
})

test_that("plotG12Functions handles edge cases", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")
  skip_if_not_installed("spatstat.random")
  
  # Create minimal test data
  set.seed(456)
  n_cells <- 10
  cell_ids <- paste0('Cell_', 1:n_cells)
  
  coords <- data.frame(
    x = runif(n_cells, 0, 10),
    y = runif(n_cells, 0, 10)
  )
  rownames(coords) <- cell_ids
  
  cell_types <- rep(c('TypeA', 'TypeB'), each = 5)
  
  meta_data <- data.frame(cellID = cell_ids)
  rownames(meta_data) <- cell_ids
  
  expr_data <- matrix(rnorm(n_cells * 10), nrow = n_cells, ncol = 10)
  rownames(expr_data) <- cell_ids
  colnames(expr_data) <- paste0('Gene_', 1:10)
  
  test_obj <- newCoProSingle(
    normalizedData = expr_data,
    locationData = coords,
    metaData = meta_data,
    cellTypes = cell_types
  )
  
  # Subset data to populate locationDataSub and cellTypesSub
  test_obj <- subsetData(test_obj, cellTypesOfInterest = c('TypeA', 'TypeB'))
  
  # Test with insufficient points (should give error when all pairs are skipped)
  expect_error({
    result <- plotG12Functions(
      test_obj,
      min_points_per_type = 20,  # More than we have (5 per type)
      verbose = FALSE
    )
  }, "No valid g_12\\(r\\) results computed")
})

test_that("plotG12Functions parameter validation works", {
  # Test invalid object type (checked first by the function)
  expect_error(
    plotG12Functions(NULL),
    "object must be a CoProSingle or CoProMulti object"
  )
  
  expect_error(
    plotG12Functions("not_an_object"),
    "object must be a CoProSingle or CoProMulti object"
  )
  
  expect_error(
    plotG12Functions(list(a = 1)),
    "object must be a CoProSingle or CoProMulti object"
  )
})

test_that("plotG12Functions returns stable list shape for all plot_type values", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")
  skip_if_not_installed("spatstat.random")

  # Mirror the fixture from the first happy-path test so spatstat reliably
  # produces non-empty g_12 results.
  set.seed(123)
  n_cells <- 60
  cell_ids <- paste0("Cell_", seq_len(n_cells))
  coords <- data.frame(x = runif(n_cells, 0, 100), y = runif(n_cells, 0, 100))
  rownames(coords) <- cell_ids
  cell_types <- sample(c("TypeA", "TypeB", "TypeC"), n_cells, replace = TRUE)
  meta_data <- data.frame(cellID = cell_ids); rownames(meta_data) <- cell_ids
  expr_data <- matrix(rnorm(n_cells * 20), nrow = n_cells, ncol = 20)
  rownames(expr_data) <- cell_ids
  colnames(expr_data) <- paste0("Gene_", seq_len(20))

  test_obj <- newCoProSingle(
    normalizedData = expr_data, locationData = coords,
    metaData = meta_data, cellTypes = cell_types
  )
  test_obj <- subsetData(
    test_obj, cellTypesOfInterest = c("TypeA", "TypeB", "TypeC")
  )

  run_pt <- function(pt) plotG12Functions(
    test_obj, r_um_range = c(5, 25), nsim = 10, verbose = FALSE, plot_type = pt
  )

  # individual: combined is NULL, individual is a facetted ggplot
  r_ind <- run_pt("individual")
  expect_named(r_ind$plot, c("combined", "individual"))
  expect_null(r_ind$plot$combined)
  expect_s3_class(r_ind$plot$individual, "ggplot")

  # both: combined and individual are both ggplots
  r_both <- run_pt("both")
  expect_named(r_both$plot, c("combined", "individual"))
  expect_s3_class(r_both$plot$combined, "ggplot")
  expect_s3_class(r_both$plot$individual, "ggplot")
})
test_that("getCorrTwoTypes keeps float32 kernels encoded and agrees with the decoded product", {
  # The helper used to take getKernelMatrix()'s materialize = TRUE default,
  # decoding an encoded float32 kernel to a double dgCMatrix (4 -> 12 bytes per
  # nonzero) just to call %*% on it. It now applies the compiled operator to
  # the still-encoded kernel.
  q <- function(e) suppressWarnings(suppressMessages(e))
  obj <- q(create_test_copro_single(n_cells = 300, n_genes = 40,
                                    n_cell_types = 2, seed = 64))
  obj <- q(subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB")))
  obj <- q(computePCA(obj, nPCA = 6))
  obj <- q(computeSparseKernelFloat32(obj, sigmaValues = c(0.05, 0.1),
                                      verbose = FALSE))
  obj <- q(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))
  obj <- q(computeNormalizedCorrelation(obj))
  obj <- q(computeGeneAndCellScores(obj))

  sigma <- obj@sigmaValueChoice
  K <- getKernelMatrix(obj, sigma, "CellTypeA", "CellTypeB",
                       verbose = FALSE, materialize = FALSE)
  expect_true(CoPro:::.isFloat32SparseKernel(K))

  df <- q(getCorrTwoTypes(obj, "CellTypeA", "CellTypeB", ccIndex = 1))
  expect_s3_class(df, "data.frame")
  expect_identical(nrow(df), ncol(K))
  expect_true(all(is.finite(df$AK)))

  # Against the decoded double product the old code performed. float32
  # accumulation makes this ~1e-6 relative, not exact -- and the analysis path
  # (.kernelXKY) already accumulates in float32, so the plot now matches the
  # numbers it is plotting rather than a higher-precision recomputation.
  scores_a <- as.matrix(getCellScores(obj, sigma = sigma,
                                      cellType = "CellTypeA", ccIndex = 1,
                                      verbose = FALSE))
  reference <- as.numeric(t(scores_a) %*% asDoubleSparseMatrix(K))
  expect_lt(max(abs(df$AK - reference)) / max(abs(reference)), 1e-5)

  # On a double kernel the helper must be exact.
  expect_identical(
    CoPro:::.leftMultiplyKernel(scores_a, asDoubleSparseMatrix(K)),
    reference
  )
})

test_that("multi-slide getCorrTwoTypes keeps float32 kernels encoded too", {
  # .getCorrTwoTypesCoreMulti() has its own getKernelMatrix() call site, which
  # was changed to materialize = FALSE alongside the single-slide one. The two
  # are separate code paths and only the single-slide one was covered.
  q <- function(e) suppressWarnings(suppressMessages(e))
  obj <- q(create_test_copro_multi(n_cells_per_slide = 200, n_slides = 2,
                                   n_genes = 40, n_cell_types = 2, seed = 65))
  obj <- q(subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB")))
  obj <- q(computePCA(obj, nPCA = 6))
  obj <- q(computeSparseKernelFloat32(obj, sigmaValues = c(0.05, 0.1),
                                      verbose = FALSE))
  obj <- q(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))
  obj <- q(computeNormalizedCorrelation(obj))
  obj <- q(computeGeneAndCellScores(obj))

  sigma <- obj@sigmaValueChoice
  slides <- getSlideList(obj)
  expect_gt(length(slides), 1L)

  df <- q(getCorrTwoTypes(obj, "CellTypeA", "CellTypeB", ccIndex = 1))
  expect_s3_class(df, "data.frame")
  expect_true("slideID" %in% names(df))
  expect_setequal(unique(df$slideID), slides)
  expect_true(all(is.finite(df$AK)))

  # Per slide, the encoded operator must agree with the decoded double product
  # the pre-change code performed.
  for (s in slides) {
    K <- getKernelMatrix(obj, sigma, "CellTypeA", "CellTypeB", slide = s,
                         verbose = FALSE, materialize = FALSE)
    expect_true(CoPro:::.isFloat32SparseKernel(K),
                info = paste("slide", s, "kernel stayed encoded"))
    scores_a <- as.matrix(getCellScores(obj, sigma = sigma,
                                        cellType = "CellTypeA", slide = s,
                                        ccIndex = 1, verbose = FALSE))
    reference <- as.numeric(t(scores_a) %*% asDoubleSparseMatrix(K))
    got <- df$AK[df$slideID == s]
    expect_identical(length(got), length(reference), label = paste("slide", s))
    expect_lt(max(abs(got - reference)) / max(abs(reference)), 1e-5)
  }
})
