# computeDistance between pairs of cell types

computeDistance between pairs of cell types

## Usage

``` r
computeDistance(
  object,
  distType = c("Euclidean2D", "Euclidean3D", "Morphology-Aware"),
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  normalizeTarget = 0.01,
  truncateLowDist = TRUE,
  verbose = TRUE,
  knn_k = 10,
  geodesic_threshold = 10,
  geodesic_cutoff = 7
)

# S4 method for class 'CoProSingle'
computeDistance(
  object,
  distType = c("Euclidean2D", "Euclidean3D", "Morphology-Aware"),
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  normalizeTarget = 0.01,
  truncateLowDist = TRUE,
  verbose = TRUE,
  knn_k = 10,
  geodesic_threshold = 10,
  geodesic_cutoff = 7
)

# S4 method for class 'CoProMulti'
computeDistance(
  object,
  distType = c("Euclidean2D", "Euclidean3D", "Morphology-Aware"),
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  normalizeTarget = 0.01,
  truncateLowDist = TRUE,
  verbose = TRUE,
  knn_k = 10,
  geodesic_threshold = 10,
  geodesic_cutoff = 7
)
```

## Arguments

- object:

  A `CoPro` object

- distType:

  Type of distance to compute: "Euclidean2D", "Euclidean3D", or
  "Morphology-Aware"

- xDistScale:

  Scale for x distance

- yDistScale:

  Scale for y distance

- zDistScale:

  Scale for z distance

- normalizeDistance:

  Whether to normalize distance? The normalization will make sure that
  the low-percentile cell-cell distance will become `normalizeTarget`
  (default 0.01), thus ensuring no matter which input scale is used for
  the distance matrix, the output will roughly be in mm^3. This ensures
  that the kernel sizes from 0.001 to 0.1 will make sense. Default =
  TRUE.

- normalizeTarget:

  Numeric scalar. The target value that the low-percentile cell-cell
  distance is rescaled to when `normalizeDistance = TRUE`. Default =
  0.01 (preserves historical behavior). Advanced users can tune this
  alongside `sigmaValues` passed to
  [`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md).

- truncateLowDist:

  Whether to truncate small distances so that the cells that are nearly
  overlapping with each other do not have a super small distance.
  Default = TRUE.

- verbose:

  Whether to print info about the quantile of the distance

- knn_k:

  Number of nearest neighbors for KNN graph construction (used only for
  Morphology-Aware distance). Default = 10.

- geodesic_threshold:

  Geodesic distance threshold for regression fitting. Cell pairs with
  geodesic distance \<= this value are used to fit the regression.
  Default = 10.

- geodesic_cutoff:

  Geodesic distance value at which to evaluate the regression for
  determining the Euclidean distance cutoff. Default = 7.

## Value

`CoPro` object with distance matrix computed

## Details

For "Morphology-Aware" distance: This method addresses the issue where
cells appear spatially close (small Euclidean distance) but are actually
topologically distant due to tissue morphology (e.g., semilunar folds in
colon tissue).

The algorithm:

1.  Compute Euclidean distances (d_E) between all cells

2.  Build a K-nearest neighbor (KNN) graph and compute geodesic
    distances (d_g) as shortest path distances on this graph

3.  Fit a linear regression between d_g and d_E for cell pairs where d_g
    \<= geodesic_threshold (default: 10)

4.  Use the fitted model to predict d_E at d_g = geodesic_cutoff
    (default: 7), this becomes cutoff_d_E

5.  For cell pairs where d_g \> geodesic_threshold AND d_E \<
    cutoff_d_E, set distance to the maximum distance in the matrix

This requires the igraph package to be installed.

## See also

[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCAMulti.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)

Other spatial-pipeline:
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCAMulti.md),
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
```
