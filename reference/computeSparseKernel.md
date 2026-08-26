# Compute sparse Gaussian kernels directly from coordinates

A fused, memory-efficient alternative to
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md) +
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md)
for large datasets. It builds, for every cell-type pair (and
within-type), a sparse Gaussian kernel using a fixed-radius neighbor
search, never forming a dense `n x n` matrix. Within-type kernels are
stored with one triangle as symmetric `dsCMatrix` objects; cross-type
kernels, and kernels made asymmetric by row/column normalization, use
`dgCMatrix`. Results are numerically equivalent to the dense path (every
pair beyond the kernel's support radius is zero anyway). Distances are
not stored.

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
  the sparse path). `NULL` inherits the recorded geometry.

- xDistScale, yDistScale, zDistScale:

  Per-axis coordinate scales. `NULL` inherits the recorded geometry.

- normalizeDistance, normalizeTarget, truncateLowDist:

  distance-processing options, matching
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md).
  `NULL` inherits the recorded geometry.

- normalizeMethod:

  How the reference distance is estimated when
  `normalizeDistance = TRUE`. `"global"` uses the median
  nearest-neighbor distance over all cells of interest, ignoring their
  type labels, so the unit is a property of the tissue rather than of
  whichever blocks this call builds. `"spacing"` measures each cell-type
  block and takes the median across blocks. `"percentile"` reproduces
  the pre-1.2.0 behavior: the minimum, across blocks, of a low quantile
  of pairwise distances.

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
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md),
[`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md),
[`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md),
[`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`runSlideLevelInference()`](https://zhen-miao.github.io/CoPro/reference/runSlideLevelInference.md),
[`selectSigmaByPermutation()`](https://zhen-miao.github.io/CoPro/reference/selectSigmaByPermutation.md)
