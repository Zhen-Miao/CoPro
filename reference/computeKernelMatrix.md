# Compute Kernel Matrix for CoPro

This method calculates the kernel matrices for pairs of cell types based
on their distances and a range of sigma values. The formula of
calculating kernel matrix is: \$\$K(x, y) =
\exp\left(-\frac{\\x-y\\^2}{2 \sigma^2}\right)\$\$ The matrices are
adjusted by clipping the upper quantile of the values to reduce the
effect of outliers. The results are stored within the object.

## Usage

``` r
computeKernelMatrix(
  object,
  sigmaValues,
  lowerLimit = 1e-07,
  upperQuantile = 0.85,
  normalizeKernel = FALSE,
  minAveCellNeighor = 2,
  rowNormalizeKernel = FALSE,
  colNormalizeKernel = FALSE,
  method = c("auto", "dense", "sparse"),
  dropDistances = TRUE,
  autoThreshold = 5000L,
  distType = NULL,
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  normalizeTarget = 0.01,
  truncateLowDist = TRUE,
  verbose = TRUE
)

# S4 method for class 'CoProSingle'
computeKernelMatrix(
  object,
  sigmaValues,
  lowerLimit = 1e-07,
  upperQuantile = 0.85,
  normalizeKernel = FALSE,
  minAveCellNeighor = 2,
  rowNormalizeKernel = FALSE,
  colNormalizeKernel = FALSE,
  method = c("auto", "dense", "sparse"),
  dropDistances = TRUE,
  autoThreshold = 5000L,
  distType = NULL,
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  normalizeTarget = 0.01,
  truncateLowDist = TRUE,
  verbose = TRUE
)

# S4 method for class 'CoProMulti'
computeKernelMatrix(
  object,
  sigmaValues,
  lowerLimit = 1e-07,
  upperQuantile = 0.85,
  normalizeKernel = FALSE,
  minAveCellNeighor = 2,
  rowNormalizeKernel = FALSE,
  colNormalizeKernel = FALSE,
  method = c("auto", "dense", "sparse"),
  dropDistances = TRUE,
  autoThreshold = 5000L,
  distType = NULL,
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  normalizeTarget = 0.01,
  truncateLowDist = TRUE,
  verbose = TRUE
)
```

## Arguments

- object:

  A `CoPro` object.

- sigmaValues:

  A vector of sigma values used for kernel calculation.

- lowerLimit:

  The lower limit for the kernel function, default is 1e-7.

- upperQuantile:

  The quantile used for clipping the kernel values, default is 0.85.

- normalizeKernel:

  Whether to normalize the kernel matrix? Default = FALSE. Note that
  normalization will not affect any downstream analyses, it is for
  numerical stability and easier interpretation only.

- minAveCellNeighor:

  What is the minimum average number of cell in the neighbor? This step
  is to help set up the expected sparsity of the kernel matrix. If a
  kernel sigma value is too small, this result in too few neighbors for
  most cells, resulting in an overly-sparse matrix that makes the
  parameter estimation hard. Thus, the sigma values that results in an
  overly-sparse matrix will be removed for later analysis.

- rowNormalizeKernel:

  Whether the kernel matrix will be row-wise normalized? Note that row
  or column wise normalization will result in an asymmetric result in
  skrCCA inference.

- colNormalizeKernel:

  Whether the kernel matrix will be column-wise normalized? Note that
  row or column wise normalization will result in an asymmetric result
  in skrCCA inference.

- method:

  One of `"auto"`, `"dense"`, or `"sparse"`. `"dense"` is the classic
  path that reads the distance matrices produced by
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md).
  `"sparse"` is a fused, memory-efficient path
  ([`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md))
  that builds sparse `dgCMatrix` kernels directly from coordinates via a
  fixed-radius neighbor search, never forming a dense `n x n` matrix,
  and does not require
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  to have been run. Results are numerically equivalent. `"auto"`
  (default) picks `"sparse"` when any per-slide cell-type block reaches
  `autoThreshold` cells or when the aggregate dense block workload
  reaches `autoThreshold^2` entries; otherwise it picks `"dense"`.

- dropDistances:

  Logical. If `TRUE` (default), the (potentially large) `@distances`
  slot is cleared after kernels are computed, since the downstream
  pipeline only needs the kernels. Set `FALSE` to keep distances for
  inspection via
  [`getDistMat()`](https://zhen-miao.github.io/CoPro/reference/getDistMat.md)
  or to recompute kernels with new sigma values without rebuilding
  distances.

- autoThreshold:

  Integer cell count at which `method = "auto"` selects the sparse path
  for any kernel-block dimension. The sparse path is also selected when
  aggregate dense block entries reach `autoThreshold^2`. Default 5000
  (about 200 MB of doubles before temporary matrices and copies).

- distType, xDistScale, yDistScale, zDistScale, normalizeDistance,
  normalizeTarget, truncateLowDist:

  Distance options used only by the sparse path (see
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  and
  [`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md)).
  `distType` defaults to `"Euclidean3D"` when the coordinates contain a
  `z` column, otherwise `"Euclidean2D"`.

- verbose:

  Whether to output the progress and related information

## Value

The `CoPro` object with computed kernel matrices added. The kernel
matrices are organized into a three-layer nested list object. The first
layer is indexed by the sigma value, and the second and the third layers
are cell types

## Note

To-do: Shall we include row or column normalization of the kernel?

## See also

[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCAMulti.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)

Other spatial-pipeline:
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCAMulti.md),
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
obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1), verbose = FALSE)
```
