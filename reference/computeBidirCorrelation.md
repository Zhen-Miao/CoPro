# Compute Bidirectional Correlation

Computes the mean of two correlations: cor(t(A_w1) %*% K, B_w2) and
cor(A_w1, K %*% B_w2), where A_w1 and B_w2 are cell scores derived from
PCA matrices multiplied by skrCCA weights for the corresponding cell
types, and K is the kernel matrix between the two cell types.

## Usage

``` r
computeBidirCorrelation(
  object,
  calculationMode = "perSlide",
  normalize_K = c("row_or_col", "sinkhorn_knopp", "none"),
  filter_kernel = TRUE,
  K_row_sum_cutoff = 0.05,
  K_col_sum_cutoff = 0.05
)

# S4 method for class 'CoPro'
computeBidirCorrelation(
  object,
  calculationMode = "perSlide",
  normalize_K = c("row_or_col", "sinkhorn_knopp", "none"),
  filter_kernel = TRUE,
  K_row_sum_cutoff = 0.05,
  K_col_sum_cutoff = 0.05
)

# S4 method for class 'CoProMulti'
computeBidirCorrelation(
  object,
  calculationMode = "perSlide",
  normalize_K = c("row_or_col", "sinkhorn_knopp", "none"),
  filter_kernel = TRUE,
  K_row_sum_cutoff = 0.05,
  K_col_sum_cutoff = 0.05
)
```

## Arguments

- object:

  A `CoPro` or `CoProMulti` object containing skrCCA results, PCA
  results, and kernel matrices.

- calculationMode:

  (CoProMulti only) "perSlide" or "aggregate". Ignored for single-slide
  objects.

- normalize_K:

  whether to normalize the kernel matrix, one of "none",
  "sinkhorn_knopp", "row_or_col"

- filter_kernel:

  whether to filter the kernel matrix, default is TRUE

- K_row_sum_cutoff:

  Numeric; cutoff for row sums when normalizing kernel matrix. Default
  5e-3.

- K_col_sum_cutoff:

  Numeric; cutoff for column sums when normalizing kernel matrix.
  Default 5e-3.

## Value

The input object with `@bidirCorrelation` populated (a list keyed by
sigma names). For single-slide objects, each entry is a data frame with
columns `sigmaValues`, `cellType1`, `cellType2`, `CC_index`,
`bidirCorrelation`. For multi-slide objects, the columns depend on
`calculationMode`.

## Details

The results are stored in the slot `bidirCorrelation`.

For multi-slide (`CoProMulti`) objects, results can be computed per
slide (default, returns a data frame with a `slideID` column) or
aggregated across slides (returns a data frame with an
`aggregateCorrelation` column).

## Examples

``` r
# Assuming `obj` is a prepared CoProSingle with PCA, kernels, and skrCCA:
# obj <- computeBidirCorrelation(obj)

# For CoProMulti per-slide results:
# objm <- computeBidirCorrelation(objm, calculationMode = "perSlide")

# For CoProMulti aggregate results:
# objm <- computeBidirCorrelation(objm, calculationMode = "aggregate")
```
