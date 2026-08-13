# Compute block-streamed float32 sparse Gaussian kernels

A large-data kernel path that never materializes a float64 sparse
kernel. Each cell-type block is enumerated once at the largest requested
radius, distances are retained temporarily as float32, and every
sigma-specific kernel is written directly as a row-compressed float32
byte buffer. The block's neighbor data are released before the next
block.

## Usage

``` r
computeSparseKernelFloat32(
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
  normalizeDistance = FALSE,
  normalizeMethod = "global",
  normalizeTarget = 0.01,
  truncateLowDist = TRUE,
  overwrite = TRUE,
  verbose = TRUE,
  nThreads = NULL
)

# S4 method for class 'CoPro'
computeSparseKernelFloat32(
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
  normalizeDistance = FALSE,
  normalizeMethod = "global",
  normalizeTarget = 0.01,
  truncateLowDist = TRUE,
  overwrite = TRUE,
  verbose = TRUE,
  nThreads = NULL
)
```

## Arguments

- object:

  A `CoProSingle` or `CoProMulti` object.

- sigmaValues:

  Positive Gaussian-kernel scale values.

- lowerLimit:

  Kernel support threshold.

- upperQuantile:

  Quantile at which large kernel values are clipped.

- normalizeKernel:

  Whether to divide all entries by the median positive row sum.

- minAveCellNeighor:

  Minimum average represented neighbors.

- rowNormalizeKernel, colNormalizeKernel:

  Whether to normalize each nonempty row or column to sum to one. They
  cannot both be `TRUE`.

- distType:

  `"Euclidean2D"` or `"Euclidean3D"`.

- xDistScale, yDistScale, zDistScale:

  Per-axis coordinate scales.

- normalizeDistance:

  Whether to rescale distances to a common unit.

- normalizeMethod:

  How the reference distance is estimated when
  `normalizeDistance = TRUE`: `"global"` (median nearest-neighbor
  distance over all cells, ignoring type labels), `"spacing"` (median
  nearest-partner distance per block, combined by median), or
  `"percentile"` (pre-1.2.0 behavior). See
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md).

- normalizeTarget:

  Target low distance percentile after normalization.

- truncateLowDist:

  Whether to floor very small distances.

- overwrite:

  Whether to replace existing kernel matrices.

- verbose:

  Whether to report progress.

- nThreads:

  Worker threads for the compiled kernel builder. Pass a positive
  integer to fix the count, including to raise the ceiling for very
  large jobs. `NULL` (default) falls back to
  `getOption("CoPro.float32Threads")`, the global flag this argument
  replaced; with neither set, the count is resolved from the cores
  actually allocated to this process, capped by `OMP_NUM_THREADS`, by
  CRAN's core limit during `R CMD check`, and by a ceiling of eight
  beyond which these memory-bandwidth-bound operators stop scaling.

## Value

The object with encoded float32 kernels in `@kernelMatrices`.

## Details

The resulting kernels are consumed directly by
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
through a parallel float32 `X1' K X2` operator. Only the small PC-space
result is converted to float64. Unnormalized and globally normalized
one-cell-type kernels store only their strict upper triangle and are
applied as symmetric matrices without expansion. Row- or
column-normalized self-kernels are expanded because those operations
make the represented matrix asymmetric.

Both single- and multi-slide objects are supported, with any positive
number of cell types. Multi-slide construction streams one
slide/cell-type block at a time and schedules the largest block first.
Older consumers can request a temporary standard sparse matrix through
[`getKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/getKernelMatrix.md)
or
[`asDoubleSparseMatrix()`](https://zhen-miao.github.io/CoPro/reference/asDoubleSparseMatrix.md).

Global, row, and column normalization are applied directly to float32
values. The ordinary centered Frobenius objective normalizer is also
evaluated from the encoded float32 kernel without materialization. Exact
whitened Frobenius normalization currently uses the compatibility path
through temporary double-precision sparse matrices when within-type
operators are supplied.

## See also

[`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`asDoubleSparseMatrix()`](https://zhen-miao.github.io/CoPro/reference/asDoubleSparseMatrix.md)

Other spatial-pipeline:
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md),
[`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md),
[`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md),
[`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`runSlideLevelInference()`](https://zhen-miao.github.io/CoPro/reference/runSlideLevelInference.md),
[`selectSigmaByPermutation()`](https://zhen-miao.github.io/CoPro/reference/selectSigmaByPermutation.md)
