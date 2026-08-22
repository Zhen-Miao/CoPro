make_score_reference_fixture <- function() {
  genes <- paste0("g", 1:3)
  cell_ids <- paste0("train", 1:24)
  cell_type <- rep(c("A", "B"), each = 6, times = 2)
  slide_id <- rep(c("s1", "s2"), each = 12)
  expression <- cbind(
    g1 = c(1:6, 3:8, 9:14, 5:10),
    g2 = c(4:9, 2:7, 12:17, 7:12),
    g3 = c(rep(2, 6), 1:6, 5:10, 6:1)
  )
  rownames(expression) <- cell_ids
  location <- data.frame(
    x = seq_len(nrow(expression)), y = seq_len(nrow(expression)),
    row.names = cell_ids
  )
  metadata <- data.frame(row.names = cell_ids)
  object <- newCoProMulti(
    expression, location, metadata, cell_type, slide_id
  )
  object <- subsetData(object, c("A", "B"))
  weights <- list(
    A = matrix(
      c(1, -1, 0.5, 0.5, 1, -0.5), nrow = 3,
      dimnames = list(genes, c("CC_1", "CC_2"))
    ),
    B = matrix(
      c(-0.5, 0.25, 1, 1, -0.5, 0.25), nrow = 3,
      dimnames = list(genes, c("CC_1", "CC_2"))
    )
  )
  object@geneScores <- list(
    "geneScores|sigma0.1|A" = weights$A,
    "geneScores|sigma0.1|B" = weights$B
  )
  object@sigmaValueChoice <- 0.1
  list(object = object, weights = weights)
}

make_score_target <- function(expression, cell_type) {
  cell_ids <- rownames(expression)
  object <- newCoProSingle(
    normalizedData = expression,
    locationData = data.frame(
      x = seq_len(nrow(expression)), y = seq_len(nrow(expression)),
      row.names = cell_ids
    ),
    metaData = data.frame(row.names = cell_ids),
    cellTypes = cell_type
  )
  subsetData(object, c("A", "B"))
}

test_that("frozen reference prediction matches the declared M1 formula", {
  fixture <- make_score_reference_fixture()
  reference <- fit_score_reference(fixture$object)
  target_x <- matrix(
    seq_len(36), nrow = 12,
    dimnames = list(paste0("target", 1:12), paste0("g", 1:3))
  )
  target <- make_score_target(target_x, rep(c("A", "B"), each = 6))

  scores <- predict(reference, target, chunk_size = 1L)
  for (cell_type in c("A", "B")) {
    rows <- target@cellTypesSub == cell_type
    fitted <- reference$references[[cell_type]]
    expected <- sweep(target_x[rows, , drop = FALSE], 2L, fitted$center, "-")
    expected <- sweep(expected, 2L, fitted$scale, "/") %*% fitted$weights
    expect_equal(scores[[cell_type]], expected, tolerance = 1e-14)
  }

  aggregate <- predict(reference, target, aggregate = TRUE)
  expect_identical(rownames(aggregate), rownames(target_x))
  expect_equal(
    as.numeric(aggregate),
    as.numeric(rbind(scores$A, scores$B)[rownames(target_x), , drop = FALSE]),
    tolerance = 1e-14
  )
})

test_that("frozen scores are target-addition and chunk invariant", {
  fixture <- make_score_reference_fixture()
  reference <- fit_score_reference(fixture$object)
  base_x <- matrix(
    seq_len(36), nrow = 12,
    dimnames = list(paste0("kept", 1:12), paste0("g", 1:3))
  )
  extra_x <- matrix(
    c(100, 1, 50, 2, 80, 3), nrow = 2,
    dimnames = list(c("extraA", "extraB"), paste0("g", 1:3))
  )
  alone <- make_score_target(base_x, rep(c("A", "B"), each = 6))
  together <- make_score_target(
    rbind(base_x, extra_x), c(rep(c("A", "B"), each = 6), "A", "B")
  )

  score_alone <- predict(reference, alone, aggregate = TRUE, chunk_size = 1L)
  score_together <- predict(
    reference, together, aggregate = TRUE, chunk_size = 100L
  )
  expect_equal(
    as.numeric(score_alone),
    as.numeric(score_together[rownames(base_x), , drop = FALSE]),
    tolerance = 0
  )
})

test_that("dense and sparse references and targets agree", {
  fixture <- make_score_reference_fixture()
  dense_reference <- fit_score_reference(fixture$object)
  sparse_object <- fixture$object
  sparse_object@normalizedDataSub <- methods::as(
    sparse_object@normalizedDataSub, "CsparseMatrix"
  )
  sparse_reference <- fit_score_reference(sparse_object)

  target_x <- matrix(
    seq_len(36), nrow = 12,
    dimnames = list(paste0("sparse_target", 1:12), paste0("g", 1:3))
  )
  target <- make_score_target(target_x, rep(c("A", "B"), each = 6))
  sparse_target <- target
  sparse_target@normalizedDataSub <- methods::as(
    sparse_target@normalizedDataSub, "CsparseMatrix"
  )

  dense_scores <- predict(dense_reference, target, aggregate = TRUE)
  sparse_scores <- predict(sparse_reference, sparse_target, aggregate = TRUE)
  expect_equal(
    as.numeric(sparse_scores), as.numeric(dense_scores), tolerance = 1e-14
  )
})

test_that("frozen reference integrates with a fitted CoPro object", {
  fitted <- create_test_copro_single(
    n_cells = 60, n_genes = 15, n_cell_types = 2, seed = 808
  )
  fitted <- subsetData(fitted, c("CellTypeA", "CellTypeB"))
  fitted <- computePCA(fitted, nPCA = 5)
  fitted <- computeKernelMatrix(
    fitted, sigmaValues = 0.5, method = "sparse",
    normalizeDistance = FALSE, verbose = FALSE
  )
  fitted <- runSkrCCA(
    fitted, scalePCs = TRUE, nCC = 1, maxIter = 100
  )
  fitted <- computeGeneAndCellScores(fitted)

  reference <- fit_score_reference(fitted, sigma = 0.5)
  transferred <- predict(reference, fitted)
  for (cell_type in fitted@cellTypesOfInterest) {
    native <- getCellScores(
      fitted, sigma = 0.5, cellType = cell_type, verbose = FALSE
    )
    expect_equal(transferred[[cell_type]], native, tolerance = 1e-8)
  }
})

test_that("equal-slide references weight slide moments equally", {
  fixture <- make_score_reference_fixture()
  reference <- fit_score_reference(
    fixture$object, reference_weight = "equal_slide"
  )
  for (cell_type in c("A", "B")) {
    rows <- fixture$object@cellTypesSub == cell_type
    x <- fixture$object@normalizedDataSub[rows, , drop = FALSE]
    slides <- getSlideID(fixture$object)[rows]
    first <- vapply(unique(slides), function(slide) {
      colMeans(x[slides == slide, , drop = FALSE])
    }, numeric(ncol(x)))
    second <- vapply(unique(slides), function(slide) {
      colMeans(x[slides == slide, , drop = FALSE] ^ 2)
    }, numeric(ncol(x)))
    expected_center <- rowMeans(first)
    expected_scale <- sqrt(rowMeans(second) - expected_center ^ 2)
    expect_equal(reference$references[[cell_type]]$center, expected_center)
    expect_equal(reference$references[[cell_type]]$scale, expected_scale)
  }

  one_gene <- fixture$object
  one_gene@geneScores <- lapply(one_gene@geneScores, function(weights) {
    weights["g1", , drop = FALSE]
  })
  one_gene_reference <- fit_score_reference(
    one_gene, reference_weight = "equal_slide"
  )
  expect_length(one_gene_reference$references$A$center, 1L)
  expect_true(is.finite(one_gene_reference$references$A$scale))
})

test_that("frozen reference validates fitted and target contracts", {
  fixture <- make_score_reference_fixture()
  unfitted <- fixture$object
  unfitted@geneScores <- list()
  expect_error(
    fit_score_reference(unfitted),
    "computeGeneAndCellScores"
  )
  no_sigma <- fixture$object
  no_sigma@sigmaValueChoice <- numeric()
  expect_error(
    fit_score_reference(no_sigma),
    "sigma must be supplied"
  )
  expect_error(
    fit_score_reference(fixture$object, sigma = -1),
    "positive finite"
  )

  reference <- fit_score_reference(fixture$object)
  target_x <- matrix(
    1:24, nrow = 12,
    dimnames = list(paste0("target", 1:12), c("g1", "g2"))
  )
  target <- make_score_target(target_x, rep(c("A", "B"), each = 6))
  expect_error(predict(reference, target), "missing frozen-reference genes")
  expect_error(predict(reference, target, chunk_size = 1.5), "positive integer")

  target@cellTypesOfInterest <- "A"
  expect_error(predict(reference, target), "same cell types")
})

test_that("frozen reference has concise print and provenance", {
  fixture <- make_score_reference_fixture()
  reference <- fit_score_reference(fixture$object)
  expect_s3_class(reference, "CoProScoreReference")
  expect_identical(reference$transform, "frozen_log_center_scale")
  expect_identical(reference$training_slides, c("s1", "s2"))
  printed <- capture.output(returned <- print(reference))
  expect_match(printed[[1L]], "CoPro frozen score reference", fixed = TRUE)
  expect_identical(returned, reference)
})
