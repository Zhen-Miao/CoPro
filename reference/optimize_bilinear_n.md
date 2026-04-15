# Run multi version of skrCCA to detect subsequent components (Single Slide) Uses flat kernel structure for consistent data access

Run multi version of skrCCA to detect subsequent components (Single
Slide) Uses flat kernel structure for consistent data access

## Usage

``` r
optimize_bilinear_n(
  X_list,
  flat_kernels,
  sigma,
  w_list,
  cellTypesOfInterest,
  nCC = 2,
  max_iter = 1000,
  tol = 1e-05,
  step_size = 1,
  sdev2_list = NULL
)
```

## Arguments

- X_list:

  Named list of data matrices (subsetted)

- flat_kernels:

  Flat list of kernel matrices

- sigma:

  Sigma value (numeric)

- w_list:

  A named list of weights (subsetted, matrices with previous components
  as columns)

- cellTypesOfInterest:

  A vector specifying cell type names present in the input lists

- nCC:

  Total number of canonical vectors desired (must be \>= 2)

- max_iter:

  Maximum number of iterations for helper function

- tol:

  Tolerance of accuracy for helper function

- step_size:

  Step size for damped power iteration (default 1)

- sdev2_list:

  Optional named list of squared standard deviations per cell type for
  weighted normalization. Default `NULL`.

## Value

A named list of weights (matrices with components 1 to nCC as columns)
