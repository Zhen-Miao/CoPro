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
})
