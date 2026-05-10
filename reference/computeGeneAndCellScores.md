# computeGeneAndCellScores

computeGeneAndCellScores

## Usage

``` r
computeGeneAndCellScores(object)

# S4 method for class 'CoPro'
computeGeneAndCellScores(object)

# S4 method for class 'CoProMulti'
computeGeneAndCellScores(object)
```

## Arguments

- object:

  A `CoPro` or `CoProMulti` object containing CCA results and kernel
  matrices.

## Value

A `CoPro` or `CoProMulti` object with gene and cell scores computed

## See also

[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md),
[`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md),
[`getCellScoresInSitu()`](https://zhen-miao.github.io/CoPro/reference/getCellScoresInSitu.md)

Other scores-and-correlation:
[`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md)
