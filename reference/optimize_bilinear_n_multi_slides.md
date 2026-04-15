# Multi-slide SkrCCA optimization - Multiple Components

Computes components 2 to nCC using deflation Uses flat kernel structure
for consistent data access

## Usage

``` r
optimize_bilinear_n_multi_slides(
  X_list_all,
  flat_kernels,
  sigma,
  slides,
  w_list,
  cellTypesOfInterest,
  nCC = 2,
  max_iter = 1000,
  tol = 1e-05,
  n_cores = 1,
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

- w_list:

  Initial weight list with first component(s)

- cellTypesOfInterest:

  Cell types to process

- nCC:

  Total number of components desired

- max_iter:

  Maximum iterations for refinement

- tol:

  Convergence tolerance

- n_cores:

  Number of cores for parallel computation

- step_size:

  Step size for damped power iteration (default 1)

- sdev2_list:

  Optional named list of squared standard deviations per cell type for
  weighted normalization. Default `NULL`.

## Value

Updated weight list with all components
