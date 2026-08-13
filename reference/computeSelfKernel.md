# Compute Self-Kernel Matrices for Multiple Cell Types

This function computes within-cell-type kernel matrices for each cell
type when multiple cell types are present. The default `method = "auto"`
uses existing self-distance matrices for small workloads and builds
sparse kernels directly from coordinates for large workloads or when
self-distances have not been materialized.

## Usage

``` r
computeSelfKernel(
  object,
  sigmaValues,
  lowerLimit = 1e-07,
  upperQuantile = 0.85,
  normalizeKernel = FALSE,
  minAveCellNeighor = 2,
  rowNormalizeKernel = FALSE,
  colNormalizeKernel = FALSE,
  verbose = TRUE,
  overwrite = FALSE,
  method = c("auto", "dense", "sparse", "float32"),
  autoThreshold = 5000L,
  distType = NULL,
  xDistScale = NULL,
  yDistScale = NULL,
  zDistScale = NULL,
  normalizeDistance = NULL,
  normalizeMethod = NULL,
  normalizeTarget = NULL,
  truncateLowDist = NULL
)

# S4 method for class 'CoPro'
computeSelfKernel(
  object,
  sigmaValues,
  lowerLimit = 1e-07,
  upperQuantile = 0.85,
  normalizeKernel = FALSE,
  minAveCellNeighor = 2,
  rowNormalizeKernel = FALSE,
  colNormalizeKernel = FALSE,
  verbose = TRUE,
  overwrite = FALSE,
  method = c("auto", "dense", "sparse", "float32"),
  autoThreshold = 5000L,
  distType = NULL,
  xDistScale = NULL,
  yDistScale = NULL,
  zDistScale = NULL,
  normalizeDistance = NULL,
  normalizeMethod = NULL,
  normalizeTarget = NULL,
  truncateLowDist = NULL
)
```

## Arguments

- object:

  A `CoPro` object with multiple cell types

- sigmaValues:

  A vector of sigma values used for kernel calculation

- lowerLimit:

  The lower limit for the kernel function, default is 1e-7

- upperQuantile:

  The quantile used for clipping the kernel values, default is 0.85

- normalizeKernel:

  Whether to normalize the kernel matrix? Default = FALSE

- minAveCellNeighor:

  Minimum average number of neighbors. Default = 2

- rowNormalizeKernel:

  Whether to row-normalize kernel matrices. Default = FALSE

- colNormalizeKernel:

  Whether to column-normalize kernel matrices. Default = FALSE

- verbose:

  Whether to output progress information

- overwrite:

  Whether to overwrite existing kernel matrices. If FALSE, will add
  self-kernel matrices to existing cross-type kernels. Default = FALSE

- method:

  One of `"auto"`, `"dense"`, `"sparse"`, or `"float32"`. Sparse paths
  construct thresholded self-kernels directly from coordinates and do
  not require
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md).
  `"auto"` uses float32 sparse storage for large workloads; request
  `"sparse"` for float64 equivalence checks.

- autoThreshold:

  Cell-count threshold used by `method = "auto"`. Float32 sparse
  construction is selected when a self-kernel dimension reaches this
  value, aggregate dense self-kernel entries reach its square, or
  required self-distance matrices are absent. Default 5000.

- distType, xDistScale, yDistScale, zDistScale, normalizeDistance,
  normalizeTarget, truncateLowDist:

  Distance options for the sparse path, matching
  [`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md):
  `NULL` inherits the geometry recorded by
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  /
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md),
  and a value that contradicts that record is an error.
  `normalizeDistance` additionally accepts `"inherit"`, which reuses the
  scaling factor
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  recorded rather than deriving one from the self-distances; see
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md).
  Building self-kernels never overwrites a factor already recorded on
  the object.

- normalizeMethod:

  How the reference distance is estimated when
  `normalizeDistance = TRUE`: `"global"` (median nearest-neighbor
  distance over all cells, ignoring type labels), `"spacing"` (median
  nearest-partner distance per block, combined by median), or
  `"percentile"` (pre-1.2.0 behavior). See
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md).

## Value

`CoPro` object with self-kernel matrices added to the kernelMatrices
slot

## Examples

``` r
if (FALSE) { # \dontrun{
# Assume you have a CoPro object with multiple cell types
# First compute cross-type distances and kernels
object <- computeDistance(object)
object <- computeKernelMatrix(object, sigmaValues = c(0.01, 0.05, 0.1))

# Then add self-distances and self-kernels
object <- computeSelfDistance(object)
object <- computeSelfKernel(object, sigmaValues = c(0.01, 0.05, 0.1))

# Now you have both cross-type and self-type kernel matrices
} # }
```
