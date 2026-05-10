# Coerce a single-cell object into a CoPro object

One-call coercion from `SingleCellExperiment` or `Seurat` objects into
`CoProSingle` (`asCoProSingle()`) or `CoProMulti` (`asCoProMulti()`)
objects. These helpers extract the relevant expression matrix, spatial
coordinates, cell metadata, and cell-type labels, then delegate to
[`newCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/newCoProSingle.md)
/
[`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md)
so all validation stays single-sourced.

## Usage

``` r
asCoProSingle(
  x,
  spatialDim = "spatial",
  cellTypeCol,
  assay = NULL,
  layer = "data",
  ...
)

asCoProMulti(
  x,
  spatialDim = "spatial",
  cellTypeCol,
  slideCol = NULL,
  assay = NULL,
  layer = "data",
  ...
)

# S4 method for class 'ANY'
asCoProSingle(
  x,
  spatialDim = "spatial",
  cellTypeCol,
  assay = NULL,
  layer = "data",
  ...
)

# S4 method for class 'ANY'
asCoProMulti(
  x,
  spatialDim = "spatial",
  cellTypeCol,
  slideCol = NULL,
  assay = NULL,
  layer = "data",
  ...
)
```

## Arguments

- x:

  A `SingleCellExperiment` or `Seurat` object.

- spatialDim:

  Name of the reducedDim / image / reduction that holds 2-D (or 3-D)
  spatial coordinates. Default `"spatial"`.

- cellTypeCol:

  Name of the column in the object's cell-level metadata that holds cell
  type labels. Required.

- assay:

  Optional assay name. For `SingleCellExperiment` defaults to
  `"logcounts"` when available, otherwise the first assay. For `Seurat`
  defaults to `SeuratObject::DefaultAssay(x)`.

- layer:

  Seurat-only. Layer (formerly "slot") to pull expression data from.
  Default `"data"`.

- ...:

  Additional arguments passed through to
  [`newCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/newCoProSingle.md)
  /
  [`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md)
  (for example for future extensions).

- slideCol:

  (`asCoProMulti()` only) Name of the column in the object's cell-level
  metadata that holds the slide/sample identifier.

## Value

A `CoProSingle` or `CoProMulti` object.

## Details

For `SingleCellExperiment` objects:

- `normalizedData` comes from `SummarizedExperiment::assay(x, assay)`
  (default `"logcounts"` if present, else the first assay), transposed
  into cells-by-genes orientation.

- `locationData` comes from
  `SingleCellExperiment::reducedDim(x, spatialDim)` (expected to have 2
  or 3 columns renamed to `x`, `y`, optionally `z`).

- `metaData` comes from `SummarizedExperiment::colData(x)` as a
  data.frame.

- `cellTypes` comes from the column of `colData(x)` named `cellTypeCol`.

For `Seurat` objects:

- `normalizedData` comes from
  `SeuratObject::GetAssayData(x, layer = layer, assay = assay)`,
  transposed into cells-by-genes orientation.

- `locationData` comes from either `x[[spatialDim]]@coordinates` (when
  `spatialDim` refers to a spatial image/slice) or
  `SeuratObject::Embeddings(x, reduction = spatialDim)` (when it refers
  to a dimensionality reduction).

- `metaData` comes from the object's `[[]]` slot.

- `cellTypes` comes from `x[[cellTypeCol]][, 1]`.

Both conversions require the corresponding package to be installed
(`SingleCellExperiment` / `SummarizedExperiment` or `SeuratObject`) and
will error with an install hint otherwise.

## See also

[`newCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/newCoProSingle.md),
[`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md),
[`CreateCoPro()`](https://zhen-miao.github.io/CoPro/reference/CreateCoPro.md)

Other object-creation:
[`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md),
[`newCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/newCoProSingle.md),
[`subsetData()`](https://zhen-miao.github.io/CoPro/reference/subsetData.md)

## Examples

``` r
toy <- readRDS(
  system.file("extdata", "toy_copro_data.rds", package = "CoPro")
)

## --- From a SingleCellExperiment -----------------------------------------
if (requireNamespace("SingleCellExperiment", quietly = TRUE) &&
    requireNamespace("SummarizedExperiment", quietly = TRUE) &&
    requireNamespace("S4Vectors", quietly = TRUE)) {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = t(toy$normalizedData))
  )
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(
    toy$metaData, cellType = toy$cellTypes
  )
  SingleCellExperiment::reducedDim(sce, "spatial") <-
    as.matrix(toy$locationData)
  obj_sce <- asCoProSingle(
    sce, spatialDim = "spatial", cellTypeCol = "cellType"
  )
}

## --- From a Seurat object -----------------------------------------------
if (requireNamespace("SeuratObject", quietly = TRUE)) {
  expr_gbc <- t(toy$normalizedData)  # genes x cells
  srt <- SeuratObject::CreateSeuratObject(
    counts = expr_gbc, min.cells = 0, min.features = 0
  )
  srt <- SeuratObject::SetAssayData(
    srt, layer = "data", new.data = expr_gbc
  )
  srt[["cellType"]] <- toy$cellTypes
  coords <- as.matrix(toy$locationData)
  rownames(coords) <- colnames(srt)
  colnames(coords) <- c("spatial_1", "spatial_2")
  srt[["spatial"]] <- SeuratObject::CreateDimReducObject(
    embeddings = coords, key = "spatial_", assay = "RNA"
  )
  obj_srt <- asCoProSingle(
    srt, spatialDim = "spatial", cellTypeCol = "cellType"
  )
}
```
