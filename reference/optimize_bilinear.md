# SkrCCA optimization function for multiple groups (Single Slide) - First Component Uses flat kernel structure for consistent data access

SkrCCA optimization function for multiple groups (Single Slide) - First
Component Uses flat kernel structure for consistent data access

## Usage

``` r
optimize_bilinear(
  X_list,
  flat_kernels,
  sigma,
  max_iter = 1000,
  tol = 1e-05,
  step_size = 1,
  sdev2_list = NULL
)
```

## Arguments

- X_list:

  Named list of data matrices (cell by PC matrix)

- flat_kernels:

  Flat list of kernel matrices with names like
  "kernel\|sigma0.1\|TypeA\|TypeB"

- sigma:

  Sigma value (numeric)

- max_iter:

  Maximum number of iterations

- tol:

  tolerance of accuracy

- step_size:

  Step size for damped power iteration. Default 1 (standard power
  iteration). Values in (0,1) blend old and new weights for smoother
  convergence, which can help with many cells or many CCs.

- sdev2_list:

  Optional named list of squared standard deviations per cell type, used
  for weighted normalization when `scalePCs = TRUE`. Default `NULL`
  (unweighted).

## Value

Named list `w_list` containing the first weight vector component.
