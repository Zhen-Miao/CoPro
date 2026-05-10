# Gene-space average per-slide CCA — subsequent components

Computes components 2 through nCC using Gram-Schmidt deflation in weight
space. After computing the gradient update for each cell type, the
update is orthogonalized against all previous CC directions before
normalizing.

## Usage

``` r
optimize_genespace_avg_corr_n(
  C_self_slide,
  C_cross_slide,
  slides,
  cell_types,
  w_list,
  nCC = 2,
  max_iter = 3000,
  tol = 1e-06,
  verbose = TRUE
)
```

## Arguments

- C_self_slide:

  Per-slide self-covariance matrices (same as first component).

- C_cross_slide:

  Per-slide cross-covariance matrices.

- slides:

  Slide IDs.

- cell_types:

  Cell type names.

- w_list:

  Named list of weight matrices from previous components. Each entry is
  a G x k matrix where k = number of components already computed.

- nCC:

  Total number of components desired (must be \> existing components).

- max_iter:

  Maximum iterations per component.

- tol:

  Convergence tolerance.

- verbose:

  Print progress.

## Value

Named list of weight matrices, each G x nCC.
