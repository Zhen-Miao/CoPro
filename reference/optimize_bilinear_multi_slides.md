# Multi-slide SkrCCA optimization - First Component

Handles both standard (multiple cell types) and within (single cell
type) cases Uses flat kernel structure for consistent data access

## Usage

``` r
optimize_bilinear_multi_slides(
  X_list_all,
  flat_kernels,
  sigma,
  slides,
  max_iter = 1000,
  tol = 1e-05,
  n_cores = 1,
  direct_solve = TRUE,
  step_size = 1,
  sdev2_list = NULL
)
```

## Arguments

- X_list_all:

  List of lists of data matrices

- flat_kernels:

  Flat list of kernel matrices

- sigma:

  Sigma value (numeric)

- slides:

  Slide IDs

- max_iter:

  Maximum number of iterations

- tol:

  Convergence tolerance

- n_cores:

  Number of cores for parallel computation

- direct_solve:

  For single cell type, use direct eigenvalue solution

- step_size:

  Step size for damped power iteration (default 1). Values in (0,1)
  blend old and new weights for smoother convergence.

- sdev2_list:

  Optional named list of squared standard deviations per cell type for
  weighted normalization. Default `NULL`.

## Value

Named list of weight vectors (first component)
