# PCA-space SUMCOR optimization – subsequent components

Computes components `k_start + 1` through `nCC` by weight-space
Gram-Schmidt deflation in the metric implied by `sdev2_list`, so
`scalePCs` remains a pure reparametrization.

## Usage

``` r
optimize_sumcor_pca_n(
  X_list_all,
  flat_kernels,
  sigma,
  slides,
  cell_types,
  w_list,
  nCC = 2,
  slideWeight = "equal",
  sdev2_list = NULL,
  max_iter = 200,
  tol = 1e-06,
  step_size = 1,
  n_cores = 1,
  verbose = FALSE,
  ops = NULL
)
```

## Arguments

- X_list_all:

  `X_list_all[[slide]][[cellType]]` cell-by-PC matrices.

- flat_kernels:

  Flat list of kernel matrices.

- sigma:

  Kernel bandwidth (numeric scalar).

- slides:

  Slide IDs.

- cell_types:

  Cell types to optimize over.

- w_list:

  Weight matrices holding the components already computed.

- nCC:

  Total number of components wanted.

- slideWeight:

  Per-slide weighting: `"equal"` (default) or `"size"`.

- sdev2_list:

  Optional named list of squared standard deviations per cell type,
  supplied when `scalePCs = FALSE`.

- max_iter:

  Maximum projected-gradient iterations.

- tol:

  Convergence tolerance on the projected-gradient norm.

- step_size:

  Damping factor in (0, 1\]. Default 1 (undamped). Values below 1
  shorten every trial step, the same stabilization the SUMCOV damped
  power iteration provides; the Armijo test still runs at the damped
  step, so the iteration stays monotone. Also damps the SUMCOV warm
  start.

- n_cores:

  Cores for the per-slide kernel products.

- verbose:

  Report progress.

- ops:

  Optional precomputed `.computeSlideOperators()` structure.

## Value

Named list of `nPC x nCC` weight matrices, with an `"objectives"`
attribute holding the per-axis objective values.
