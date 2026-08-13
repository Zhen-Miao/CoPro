# What objective produced an object's CCA weights

Reads the provenance record
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
attaches to `@skrCCAOut`. Objects computed before the `objective`
argument existed carry no record and are reported as `"sumcov"`, which
is what they were computed under.

## Usage

``` r
getCCAObjective(object)
```

## Arguments

- object:

  A `CoPro` or `CoProMulti` object.

## Value

A list with `space`, `objective`, `requested`, `slideWeight`, `slides`
and `droppedSlides`.

## See also

[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)

Other scores-and-correlation:
[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md),
[`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md),
[`getNormalizerInfo()`](https://zhen-miao.github.io/CoPro/reference/getNormalizerInfo.md)

## Examples

``` r
# \donttest{
toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
obj <- newCoProSingle(
  normalizedData = toy$normalizedData,
  locationData   = toy$locationData,
  metaData       = toy$metaData,
  cellTypes      = toy$cellTypes
)
obj <- subsetData(obj, cellTypesOfInterest = unique(toy$cellTypes))
obj <- computePCA(obj, nPCA = 10)
#> Input is dense (matrix), performing irlba PCA...
#> Input is dense (matrix), performing irlba PCA...
obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE)
#> Warning: The kernel at sigma = 0.1 is predicted to be 52% dense, so a fixed-radius sparse kernel saves little or no memory here (sparse storage costs more than dense past ~67% density).
#>   Sigma below about 0.0756 would keep it under 30%. Use detectSigmaRange() to pick sigma from the data, or build one sigma at a time to bound peak memory.
obj <- runSkrCCA(obj, nCC = 2)
#> Running skrCCA [1/1] for sigma = 0.1 ...
#> skrCCA finished 1 sigma value(s) in 0.0 s.
#> Optimization succeeded for 1 sigma value(s): sigma_0.1
getCCAObjective(obj)
#> $space
#> [1] "pca"
#> 
#> $objective
#> [1] "sumcov"
#> 
#> $requested
#> [1] "sumcov"
#> 
#> $slideWeight
#> [1] NA
#> 
#> $sweep
#> [1] NA
#> 
#> $slides
#> [1] NA
#> 
#> $droppedSlides
#> character(0)
#> 
# }
```
