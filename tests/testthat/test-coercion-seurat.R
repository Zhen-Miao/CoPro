test_that("asCoProSingle() coerces a Seurat object", {
  skip_if_not_installed("SeuratObject")

  td <- generate_test_data_single(n_cells = 30, n_genes = 12, seed = 5)

  # Seurat uses genes-by-cells; transpose.
  counts <- t(td$normalizedData)
  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  # Place expression data into the "data" layer without re-normalizing.
  seu <- SeuratObject::SetAssayData(seu, layer = "data", new.data = counts)

  # Add metadata and an embedding we can treat as spatial coordinates.
  seu <- SeuratObject::AddMetaData(seu, metadata = td$cellTypes,
                                   col.name = "cell_type")
  coord_mat <- as.matrix(td$locationData)
  colnames(coord_mat) <- c("spatial_1", "spatial_2")
  seu[["spatial"]] <- SeuratObject::CreateDimReducObject(
    embeddings = coord_mat, key = "spatial_", assay = SeuratObject::DefaultAssay(seu)
  )

  obj <- asCoProSingle(seu, spatialDim = "spatial", cellTypeCol = "cell_type")

  expect_s4_class(obj, "CoProSingle")
  expect_equal(ncol(obj@normalizedData), ncol(td$normalizedData))
  expect_equal(nrow(obj@normalizedData), nrow(td$normalizedData))
  expect_equal(colnames(obj@locationData), c("x", "y"))
  expect_equal(as.matrix(obj@normalizedData), td$normalizedData, tolerance = 1e-10)
  expect_equal(as.character(obj@cellTypes), td$cellTypes)
  expect_false("cell_type" %in% colnames(obj@metaData))
})

test_that("asCoProSingle reorders Seurat embeddings to match colnames(x)", {
  skip_if_not_installed("SeuratObject")

  td <- generate_test_data_single(n_cells = 25, n_genes = 10, seed = 7)
  counts <- t(td$normalizedData)
  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  seu <- SeuratObject::SetAssayData(seu, layer = "data", new.data = counts)
  seu <- SeuratObject::AddMetaData(
    seu, metadata = td$cellTypes, col.name = "cell_type"
  )

  # Deliberately shuffle embedding row order; asCoProSingle must realign to
  # colnames(seu) so locationData rows pair correctly with cells.
  coord_mat <- as.matrix(td$locationData)
  rownames(coord_mat) <- colnames(seu)
  colnames(coord_mat) <- c("spatial_1", "spatial_2")
  set.seed(11)
  shuf_idx <- sample(nrow(coord_mat))
  coord_shuf <- coord_mat[shuf_idx, , drop = FALSE]
  seu[["spatial"]] <- SeuratObject::CreateDimReducObject(
    embeddings = coord_shuf, key = "spatial_",
    assay = SeuratObject::DefaultAssay(seu)
  )

  obj <- asCoProSingle(seu, spatialDim = "spatial", cellTypeCol = "cell_type")

  # Row i of locationData must correspond to cell colnames(seu)[i],
  # not to the shuffled order we passed in.
  expect_equal(rownames(obj@locationData), colnames(seu))
  expect_equal(obj@locationData$x, unname(coord_mat[colnames(seu), "spatial_1"]))
  expect_equal(obj@locationData$y, unname(coord_mat[colnames(seu), "spatial_2"]))
})

test_that("asCoProSingle() errors for bad spatialDim in Seurat", {
  skip_if_not_installed("SeuratObject")

  td <- generate_test_data_single(n_cells = 20, n_genes = 10, seed = 10)
  counts <- t(td$normalizedData)
  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  seu <- SeuratObject::SetAssayData(seu, layer = "data", new.data = counts)
  seu <- SeuratObject::AddMetaData(seu, metadata = td$cellTypes, col.name = "cell_type")

  expect_error(
    asCoProSingle(seu, spatialDim = "nonexistent", cellTypeCol = "cell_type"),
    "nonexistent.*neither.*image.*reduction"
  )
})

test_that("asCoProSingle() errors for bad cellTypeCol in Seurat", {
  skip_if_not_installed("SeuratObject")

  td <- generate_test_data_single(n_cells = 20, n_genes = 10, seed = 11)
  counts <- t(td$normalizedData)
  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  seu <- SeuratObject::SetAssayData(seu, layer = "data", new.data = counts)
  seu <- SeuratObject::AddMetaData(seu, metadata = td$cellTypes, col.name = "cell_type")
  coord_mat <- as.matrix(td$locationData)
  colnames(coord_mat) <- c("spatial_1", "spatial_2")
  seu[["spatial"]] <- SeuratObject::CreateDimReducObject(
    embeddings = coord_mat, key = "spatial_", assay = SeuratObject::DefaultAssay(seu)
  )

  expect_error(
    asCoProSingle(seu, spatialDim = "spatial", cellTypeCol = "wrong_col"),
    "wrong_col.*not found"
  )
})
