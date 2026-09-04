# Report the normalizer actually used

The denominator of the normalized correlation determines how the
criterion behaves across the sigma grid, and it used to be decided
implicitly by whether the object happened to carry within-type kernels.
This records the resolved choice on the object so a stored result can be
interpreted later.

## Usage

``` r
getNormalizerInfo(object)
```

## Arguments

- object:

  A `CoPro` or `CoProMulti` object that has been through
  [`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md).

## Value

A list with the normalizer `mode`, a human-readable `description`, and
the fitted autocorrelation `ranges` when `normalizer = "variogram"` was
used. `NULL` if the object predates this record.

## See also

[`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md)

Other scores-and-correlation:
[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md),
[`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md),
[`getCCADiagnostics()`](https://zhen-miao.github.io/CoPro/reference/getCCADiagnostics.md),
[`getCCAObjective()`](https://zhen-miao.github.io/CoPro/reference/getCCAObjective.md)

## Examples

``` r
# obj <- computeNormalizedCorrelation(obj, normalizer = "variogram")
# getNormalizerInfo(obj)
```
