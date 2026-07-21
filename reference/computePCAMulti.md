# Compute PCA on Integrated Multi-Slide Data

Performs PCA on the integrated data stored within the `CoProMulti`
object. Assumes integration has created a common space across slides.

## Usage

``` r
computePCA(
  object,
  nPCA = 40,
  center = TRUE,
  scale. = TRUE,
  scalePCs = TRUE,
  dataUse = "raw",
  center_per_slide = FALSE
)
```

## Arguments

- object:

  A `CoProMulti` object with the `integratedData` slot populated.

- nPCA:

  Number of principal components to compute for each cell type.

- center:

  Whether to center the matrix before PCA

- scale.:

  Whether to scale the matrix before PCA

- scalePCs:

  Whether to scale (whiten) PCs by their standard deviation before
  downstream CCA optimization. Default `TRUE` (recommended). When PCs
  are whitened, the unit-norm constraint `||w|| = 1` is equivalent to
  the standard CCA constraint `||w'X'Xw|| = 1`. Setting this to `FALSE`
  distorts the constraint space and may produce unreliable results.

- dataUse:

  What data to use, choices between "raw" and "integrated". Default is
  "raw". For single slide, this argument is ignored.

- center_per_slide:

  After the global PCA, do we do center per slide again? By default this
  is set to FALSE

## Value

A `CoProMulti` object with the `pcaResults` slot populated. `pcaResults`
structure: `list(slideID = list(cellType = pc_matrix))`.

## See also

[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)

Other spatial-pipeline:
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md),
[`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)

## Examples

``` r
toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
obj <- newCoProSingle(
  normalizedData = toy$normalizedData,
  locationData   = toy$locationData,
  metaData       = toy$metaData,
  cellTypes      = toy$cellTypes
)
obj <- subsetData(obj, cellTypesOfInterest = unique(toy$cellTypes))
obj <- computePCA(obj, nPCA = 10)
#> Input is dense (matrixarray), performing irlba pca...
#> Input is dense (matrixarray), performing irlba pca...
```
