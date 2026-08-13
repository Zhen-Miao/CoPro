# PCA-space SUMCOR optimization – first component

Maximizes the per-slide self-normalized objective \$\$f(w) =
\frac{\sum_s \sum\_{i\<j} m\_{ij}^{(s)}\\ w_i' Y\_{ij}^{(s)} w_j /
(\sigma_i^{(s)} \sigma_j^{(s)})}{\sum_s \sum\_{i\<j} m\_{ij}^{(s)}}\$\$
where \\\sigma_i^{(s)} = \sqrt{w_i' X_i^{(s)'} X_i^{(s)} w_i}\\ and
\\m\_{ij}^{(s)}\\ is 1 for `slideWeight = "equal"` or \\\sqrt{n_i^{(s)}
n_j^{(s)}}\\ for `"size"`.

## Usage

``` r
optimize_sumcor_pca(
  X_list_all,
  flat_kernels,
  sigma,
  slides,
  cell_types,
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

Named list of single-column weight matrices, with attributes
`"objective"` and `"slideWeight"`.

## Details

When one slide has coinciding per-pair constants (see the file header),
the existing SUMCOV result is used directly and no SUMCOR iteration is
run.
