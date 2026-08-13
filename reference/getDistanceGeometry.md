# Get the coordinate geometry a CoPro object's distances and kernels use

Returns the record written by whichever step built the object's
coordinate basis –
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md),
[`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md),
or
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md)
on its sparse path. Use it whenever an analysis has to reproduce the
geometry the kernels were fit on, for instance to rebuild a spatial
neighbor graph for a permutation null.

## Usage

``` r
getDistanceGeometry(object)

# S4 method for class 'CoPro'
getDistanceGeometry(object)
```

## Arguments

- object:

  A `CoPro` object.

## Value

A named list with `distType`, `xDistScale`, `yDistScale`, `zDistScale`,
`normalizeDistance`, `normalizeTarget`, `truncateLowDist`, and `source`
(the function that wrote the record). `NULL` when no distance- or
kernel-building step has run, or when `object` predates the
`distanceGeometry` slot.

## Details

The record persists after `computeKernelMatrix(dropDistances = TRUE)`
(the default) has cleared `@distances`, so a stored object still knows
which coordinates its kernels live on.

## See also

[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md)

Other accessors:
[`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md),
[`getCellScoresInSitu()`](https://zhen-miao.github.io/CoPro/reference/getCellScoresInSitu.md),
[`getDistMat()`](https://zhen-miao.github.io/CoPro/reference/getDistMat.md),
[`getKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/getKernelMatrix.md)

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
getDistanceGeometry(obj)
#> $distType
#> [1] "Euclidean2D"
#> 
#> $xDistScale
#> [1] 1
#> 
#> $yDistScale
#> [1] 1
#> 
#> $zDistScale
#> [1] 1
#> 
#> $normalizeDistance
#> [1] FALSE
#> 
#> $normalizeMethod
#> [1] "global"
#> 
#> $normalizeTarget
#> [1] 0.01
#> 
#> $truncateLowDist
#> [1] TRUE
#> 
#> $source
#> [1] "computeDistance"
#> 
```
