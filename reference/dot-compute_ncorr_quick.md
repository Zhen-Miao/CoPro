# Quick Normalized Correlation Computation (Internal)

Helper function to quickly compute normalized correlation for a given
set of PC matrices and weights.

## Usage

``` r
.compute_ncorr_quick(PCmats, w_list, flat_kernels, sigma, cts, tol = 1e-04)
```

## Arguments

- PCmats:

  Named list of PC matrices

- w_list:

  Named list of weight vectors (matrices with 1 column)

- flat_kernels:

  Flat list of kernel matrices from CoPro object

- sigma:

  Sigma value for kernel selection

- cts:

  Cell types of interest

- tol:

  Tolerance for SVD computation

## Value

Numeric value of normalized correlation
