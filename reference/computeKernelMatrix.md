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
  method = c("auto", "dense", "sparse", "float32"),
  dropDistances = TRUE,
  autoThreshold = 5000L,
  denseThreshold = 0.3,
  distType = NULL,
  xDistScale = NULL,
  yDistScale = NULL,
  zDistScale = NULL,
  normalizeDistance = NULL,
  normalizeMethod = NULL,
  normalizeTarget = NULL,
  truncateLowDist = NULL,
  verbose = TRUE
)

# S4 method for class 'CoPro'
computeKernelMatrix(
  object,
  sigmaValues,
  lowerLimit = 1e-07,
  upperQuantile = 0.85,
  normalizeKernel = FALSE,
  minAveCellNeighor = 2,
  rowNormalizeKernel = FALSE,
  colNormalizeKernel = FALSE,
  method = c("auto", "dense", "sparse", "float32"),
  dropDistances = TRUE,
  autoThreshold = 5000L,
  denseThreshold = 0.3,
  distType = NULL,
  xDistScale = NULL,
  yDistScale = NULL,
  zDistScale = NULL,
  normalizeDistance = NULL,
  normalizeMethod = NULL,
  normalizeTarget = NULL,
  truncateLowDist = NULL,
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

  Whether to normalize the kernel matrix (default `FALSE`). With neither
  row nor column normalization selected, divide by the median row sum
  among row sums above the numerical cutoff. This is not spectral- or
  Frobenius-norm normalization and does not guarantee that the fitting
  association ratio is bounded by one. Pair- or slide-specific scaling
  can change relative weights in joint fitting and hence downstream
  results. For the separate scale-normalized evaluation statistic, see
  `computeNormalizedCorrelation(normalizer = "unwhitened")`; that
  argument does not change the fitting objective.

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

  One of `"auto"`, `"dense"`, `"sparse"`, or `"float32"`. `"dense"` is
  the classic path that reads the distance matrices produced by
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md).
  `"sparse"` and `"float32"` are fused, memory-efficient paths
  ([`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md)
  and
  [`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md))
  that build kernels directly from coordinates via a fixed-radius
  neighbor search, never forming a dense `n x n` matrix, and do not
  require
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  to have been run. `"float32"` stores 8 bytes per entry against
  `"sparse"`'s 12 and streams one block at a time instead of caching
  every block's neighbors, so it is the cheaper of the two; `"sparse"`
  keeps float64 values for exactness checks.

  `"auto"` (default) picks `"dense"` for small data, and otherwise
  `"float32"`. Small means every per-slide cell-type block is under
  `autoThreshold` cells and the aggregate dense workload is under
  `autoThreshold^2` entries. Above that, a neighbor probe predicts how
  dense the kernel will actually be and warns when a fixed-radius kernel
  would retain more than `denseThreshold` of each block, since sparse
  storage saves nothing there.

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

- denseThreshold:

  Predicted kernel density above which `method = "auto"` warns that a
  sparse representation is a poor fit and suggests a smaller sigma.
  Default 0.3. A `dgCMatrix` costs 12 bytes per stored entry against 8
  for a dense double, so sparse storage is strictly worse past about
  0.67.

- distType, xDistScale, yDistScale, zDistScale, normalizeDistance,
  normalizeMethod, normalizeTarget, truncateLowDist:

  Distance options used only by the sparse paths (see
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  and
  [`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md)).
  Each defaults to `NULL`, meaning "inherit the geometry recorded by
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)",
  so sparse kernels are built on the same coordinates as the distances
  rather than on this function's own defaults. When nothing has been
  recorded, the fallbacks are `xDistScale` / `yDistScale` / `zDistScale`
  `= 1`, `normalizeDistance = FALSE`, `normalizeMethod = "global"`,
  `normalizeTarget = 0.01`, `truncateLowDist = TRUE`, and a `distType`
  of `"Euclidean3D"` when the coordinates contain a `z` column,
  otherwise `"Euclidean2D"`. Passing a value that contradicts the
  recorded geometry is an error; inspect the record with
  [`getDistanceGeometry()`](https://zhen-miao.github.io/CoPro/reference/getDistanceGeometry.md).

  `normalizeDistance` defaulted to `TRUE` before CoPro 1.2.0. Use
  [`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md)
  to choose sigma in the data's own units instead.

- verbose:

  Whether to output the progress and related information

## Value

The `CoPro` object with computed kernel matrices added. The kernel
matrices are organized into a three-layer nested list object. The first
layer is indexed by the sigma value, and the second and the third layers
are cell types

## See also

[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)

Other spatial-pipeline:
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md),
[`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md),
[`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md),
[`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md),
[`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`runSlideLevelInference()`](https://zhen-miao.github.io/CoPro/reference/runSlideLevelInference.md),
[`selectSigmaByPermutation()`](https://zhen-miao.github.io/CoPro/reference/selectSigmaByPermutation.md)

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
