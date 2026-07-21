# Compute Normalized Correlation (approximation)

This method calculates the normalized correlation between pairs of cell
types based on CCA weights and the respective kernel matrix. It uses the
whitened-Frobenius norm \|\|R_x^(1/2) K_c R_y^(1/2)\|\|\_F of the kernel
for normalization (R_x, R_y = matched-sigma within-type kernels).

## Usage

``` r
computeNormalizedCorrelation(object, tol = 1e-04, calculationMode = "perSlide")

# S4 method for class 'CoPro'
computeNormalizedCorrelation(object, tol = 1e-04)

# S4 method for class 'CoProMulti'
computeNormalizedCorrelation(object, tol = 1e-04, calculationMode = "perSlide")
```

## Arguments

- object:

  A `CoPro` or `CoProMulti` object containing CCA results and kernel
  matrices.

- tol:

  tolerance for approximate SVD calculation

- calculationMode:

  (for CoProMulti only) either "perSlide" or "aggregate", for single
  slide analysis, it is ignored, with default value "perSlide".

## Value

The object with the normalized correlation value between any pair of
cell types added as a new slot, `normalizedCorrelation`.

## See also

[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`computeBidirCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeBidirCorrelation.md),
[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md)

Other scores-and-correlation:
[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md)
