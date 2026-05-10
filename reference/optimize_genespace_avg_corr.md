# Gene-space average per-slide CCA — first component

Power iteration to find the first canonical component that maximizes the
average per-slide canonical correlation across all slides. Each slide's
contribution is self-normalized by its own score standard deviation.

## Usage

``` r
optimize_genespace_avg_corr(
  C_self_slide,
  C_cross_slide,
  slides,
  cell_types,
  max_iter = 3000,
  tol = 1e-06,
  verbose = TRUE
)
```

## Arguments

- C_self_slide:

  Named list of per-slide self-covariance matrices. Structure:
  `C_self_slide[[slide]][[cell_type]]` = G x G matrix.

- C_cross_slide:

  Named list of per-slide cross-covariance matrices. Structure:
  `C_cross_slide[[slide]][["ctA-ctB"]]` = G x G matrix.

- slides:

  Character vector of slide IDs.

- cell_types:

  Character vector of cell type names.

- max_iter:

  Maximum iterations (default 3000). Must be \>= 1.

- tol:

  Convergence tolerance on max weight change (default 1e-6).

- verbose:

  Print progress every 500 iterations (default TRUE).

## Value

Named list of weight vectors, one per cell type (each a G x 1 matrix).
