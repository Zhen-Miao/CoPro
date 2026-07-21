# Compute sparse Gaussian kernels directly from coordinates

A fused, memory-efficient alternative to
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md) +
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md)
for large datasets. It builds, for every cell-type pair (and
within-type), a sparse `dgCMatrix` Gaussian kernel using a fixed-radius
neighbor search, never forming a dense `n x n` matrix. Results are
numerically equivalent to the dense path (every pair beyond the kernel's
support radius is zero anyway). Distances are not stored.

## Usage

``` r
computeSparseKernel(
  object,
  sigmaValues,
  lowerLimit = 1e-07,
  upperQuantile = 0.85,
  normalizeKernel = FALSE,
  minAveCellNeighor = 2,
  rowNormalizeKernel = FALSE,
  colNormalizeKernel = FALSE,
  distType = c("Euclidean2D", "Euclidean3D"),
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  normalizeTarget = 0.01,
  truncateLowDist = TRUE,
  verbose = TRUE
)

# S4 method for class 'CoProSingle'
computeSparseKernel(
  object,
  sigmaValues,
  lowerLimit = 1e-07,
  upperQuantile = 0.85,
  normalizeKernel = FALSE,
  minAveCellNeighor = 2,
  rowNormalizeKernel = FALSE,
  colNormalizeKernel = FALSE,
  distType = c("Euclidean2D", "Euclidean3D"),
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  normalizeTarget = 0.01,
  truncateLowDist = TRUE,
  verbose = TRUE
)

# S4 method for class 'CoProMulti'
computeSparseKernel(
  object,
  sigmaValues,
  lowerLimit = 1e-07,
  upperQuantile = 0.85,
  normalizeKernel = FALSE,
  minAveCellNeighor = 2,
  rowNormalizeKernel = FALSE,
  colNormalizeKernel = FALSE,
  distType = c("Euclidean2D", "Euclidean3D"),
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

- distType:

  "Euclidean2D" or "Euclidean3D" (Morphology-Aware is not supported by
  the sparse path).

- xDistScale, yDistScale, zDistScale:

  per-axis coordinate scales.

- normalizeDistance, normalizeTarget, truncateLowDist:

  distance-processing options, matching
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md).

- verbose:

  Whether to output the progress and related information

## Value

The `CoPro` object with sparse kernel matrices in `@kernelMatrices`.

## See also

[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)

Other spatial-pipeline:
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCAMulti.md),
[`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
