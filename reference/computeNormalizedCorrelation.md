# Compute Normalized Correlation (approximation)

This method calculates the normalized correlation between pairs of cell
types based on CCA weights and the respective kernel matrix. It uses the
spectral norm of the kernel matrix for normalization.

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
