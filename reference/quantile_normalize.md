# quantile normalize data

quantile normalize data

## Usage

``` r
quantile_normalize(
  A,
  B,
  save_Sparse = FALSE,
  ties_method = "min",
  verbose = TRUE
)
```

## Arguments

- A:

  The reference cell-by-feature matrix Values in this matrix define the
  target distribution for each feature.

- B:

  The cell-by-feature matrix we want to normalize to A

- save_Sparse:

  whether to save as a sparse matrix

- ties_method:

  If there are ties, which method to choose

- verbose:

  whether to print progress updates (default TRUE)

## Value

The normalized B matrix
