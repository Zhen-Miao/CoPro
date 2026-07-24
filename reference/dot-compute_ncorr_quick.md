# Quick Normalized Correlation Computation (Internal)

Helper function to quickly compute normalized correlation for a given
set of PC matrices and weights.

## Usage

``` r
.compute_ncorr_quick(
  PCmats,
  w_list,
  flat_kernels,
  sigma,
  cts,
  tol = 1e-04,
  kernel_info = NULL,
  Y_resi = NULL,
  grams = NULL
)
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

- kernel_info:

  Optional precomputed list containing `K` and `norm_K12` for this
  sigma.

- Y_resi:

  Optional precomputed PC-space operator from `compute_Y_resi()`. When
  supplied, it is used for the numerator so the kernel-vector product is
  not repeated.

- grams:

  Optional named list of per-cell-type Gram matrices `crossprod(X)` used
  for the score norms in the denominator. Valid whenever the permutation
  is a bijection on cells; see `.permutationGrams()`.

## Value

Numeric value of normalized correlation
