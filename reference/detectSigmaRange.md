# Detect a usable sigma range from the data

Chooses Gaussian bandwidths by the only quantity that is comparable
across datasets: how many neighbors a cell is effectively coupled to.
For a cell `a`, the effective neighbor count at bandwidth sigma is the
kernel row sum `m_a(sigma) = sum_b exp(-0.5 (d_ab / sigma)^2)`. This
function finds the sigma values at which the median cell reaches
`minNeighbors` and `maxNeighbors`, for every cell-type pair (and slide),
and reports a shared range plus a recommended grid.

## Usage

``` r
detectSigmaRange(
  object,
  minNeighbors = 5,
  maxNeighbors = 20,
  nSigma = 5L,
  nAnchor = 500L,
  nNeighbor = 128L,
  distType = NULL,
  xDistScale = NULL,
  yDistScale = NULL,
  zDistScale = NULL,
  verbose = TRUE
)

# S4 method for class 'CoPro'
detectSigmaRange(
  object,
  minNeighbors = 5,
  maxNeighbors = 20,
  nSigma = 5L,
  nAnchor = 500L,
  nNeighbor = 128L,
  distType = NULL,
  xDistScale = NULL,
  yDistScale = NULL,
  zDistScale = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  A `CoProSingle` or `CoProMulti` object, after
  [`subsetData()`](https://zhen-miao.github.io/CoPro/reference/subsetData.md).

- minNeighbors, maxNeighbors:

  Target effective neighbor counts bracketing the useful range. Defaults
  of 5 and 20 keep the kernel above the point where cells are
  effectively isolated and below the point where every pair is coupled.

- nSigma:

  Number of bandwidths in the recommended grid, spaced logarithmically
  across the detected range.

- nAnchor:

  Anchor cells sampled per block.

- nNeighbor:

  Nearest partners retained per anchor. Truncation error in the
  effective count is about `exp(-nNeighbor / m)`, so the default of 128
  comfortably covers targets up to a few dozen neighbors.

- distType, xDistScale, yDistScale, zDistScale:

  Coordinate geometry. When omitted, any geometry already recorded on
  the object is used. Anchors are chosen by a deterministic stride, not
  at random, so repeated calls on the same object return the same
  bandwidths and the dense and sparse kernel paths agree exactly on the
  distance scale.

- verbose:

  Whether to report per-block progress.

## Value

An object of class `CoProSigmaRange`: a list with `sigmaValues` (the
recommended grid), `sigmaRange` (lower and upper bound), `feasible`
(whether one range satisfies every block), and `blocks` (a per-block
diagnostic `data.frame` with the bandwidths that block would need, its
median nearest-partner spacing, and its effective neighbor count at the
recommended bounds).

## Details

The returned bandwidths are in the coordinate units of `locationData`,
so they can be passed straight to
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md)
with `normalizeDistance = FALSE`. There is no need to rescale distances
first: selecting sigma from the data is what makes the analysis
unit-independent.

Speed comes from sampling. Only `nAnchor` cells per block are examined,
each against its `nNeighbor` nearest partners found by fixed-radius
search, so cost scales with the sample rather than with the number of
pairs. The estimate is deliberately approximate; the downstream
bandwidth selection in
[`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md)
refines the choice within the range.

## See also

[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md)

Other spatial-pipeline:
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md),
[`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md),
[`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md),
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
rng <- detectSigmaRange(obj, verbose = FALSE)
#> detectSigmaRange: normalizeDistance now defaults to FALSE (changed in CoPro 1.2.0); sigma is interpreted in raw coordinate units.
#>   Use detectSigmaRange() to choose sigma for this dataset, or pass normalizeDistance = TRUE to reproduce earlier results.
rng$sigmaValues
#> [1] 0.0979 0.1200 0.1470 0.1790 0.2190
```
