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
  method = c("auto", "dense", "sparse"),
  autoThreshold = 5000L,
  distType = NULL,
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  normalizeTarget = 0.01,
  truncateLowDist = TRUE
)

# S4 method for class 'CoProSingle'
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
  method = c("auto", "dense", "sparse"),
  autoThreshold = 5000L,
  distType = NULL,
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  normalizeTarget = 0.01,
  truncateLowDist = TRUE
)

# S4 method for class 'CoProMulti'
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
  method = c("auto", "dense", "sparse"),
  autoThreshold = 5000L,
  distType = NULL,
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  normalizeTarget = 0.01,
  truncateLowDist = TRUE
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

  One of `"auto"`, `"dense"`, or `"sparse"`. The sparse path constructs
  exact thresholded self-kernels directly from coordinates and does not
  require
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md).

- autoThreshold:

  Cell-count threshold used by `method = "auto"`. Sparse construction is
  selected when a self-kernel dimension reaches this value, aggregate
  dense self-kernel entries reach its square, or required self-distance
  matrices are absent. Default 5000.

- distType, xDistScale, yDistScale, zDistScale, normalizeDistance,
  normalizeTarget, truncateLowDist:

  Distance options for the sparse path, matching
  [`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md).

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
