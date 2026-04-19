test_that("asCoProSingle() coerces a SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  td <- generate_test_data_single(n_cells = 40, n_genes = 15, seed = 3)

  # CoPro stores data as cells-by-genes, but SCE uses genes-by-cells;
  # transpose for the SCE build.
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = t(td$normalizedData)),
    colData = S4Vectors::DataFrame(
      cell_type = td$cellTypes,
      total_counts = td$metaData$total_counts
    )
  )
  colnames(sce) <- rownames(td$normalizedData)
  SingleCellExperiment::reducedDim(sce, "spatial") <- as.matrix(td$locationData)

  obj <- asCoProSingle(sce, spatialDim = "spatial", cellTypeCol = "cell_type")

  expect_s4_class(obj, "CoProSingle")
  expect_equal(nrow(obj@normalizedData), nrow(td$normalizedData))
  expect_equal(ncol(obj@normalizedData), ncol(td$normalizedData))
  expect_equal(as.character(obj@cellTypes), td$cellTypes)
  expect_equal(colnames(obj@locationData), c("x", "y"))
})

test_that("asCoProMulti() coerces a SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("S4Vectors")

  td <- generate_test_data_multi(
    n_cells_per_slide = 25, n_slides = 2, n_genes = 10, seed = 4
  )

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = t(td$normalizedData)),
    colData = S4Vectors::DataFrame(
      cell_type = td$cellTypes,
      slide_id  = td$slideID
    )
  )
  colnames(sce) <- rownames(td$normalizedData)
  SingleCellExperiment::reducedDim(sce, "spatial") <- as.matrix(td$locationData)

  obj <- asCoProMulti(
    sce,
    spatialDim = "spatial",
    cellTypeCol = "cell_type",
    slideCol = "slide_id"
  )

  expect_s4_class(obj, "CoProMulti")
  expect_equal(length(obj@slideList), 2)
})

test_that("asCoProSingle() errors with a helpful message for unsupported input", {
  expect_error(asCoProSingle(list(a = 1), cellTypeCol = "ct"),
               "SingleCellExperiment or Seurat")
})
