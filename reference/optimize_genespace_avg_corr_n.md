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
  step_size = 1,
  verbose = TRUE,
  sweep = c("gauss-seidel", "jacobi"),
  objective = c("sumcor", "sumcov")
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

- step_size:

  Step size for damped power iteration (default 1). See
  [`optimize_genespace_avg_corr`](https://zhen-miao.github.io/CoPro/reference/optimize_genespace_avg_corr.md)
  for details.

- verbose:

  Print progress.

- sweep:

  Which block sweep to use.

  `"gauss-seidel"` (default) updates each cell type using the blocks
  already updated in the current sweep. With `sigma` held fixed this
  makes each block update the exact maximizer over that block – the
  objective is linear in `w_i` once the others are fixed – so
  `w_i' g_i = \|g_i\| >= 0` always and the iteration cannot come to rest
  on a negative-objective solution.

  `"jacobi"` updates every block from the previous iterate. All blocks
  then move off stale values, which is not coordinate ascent: the
  iteration can settle on the *negative* singular pair as a period-2
  sign orbit, which `check_convergence()`'s sign-tolerant test accepts
  as converged. That is what the post-hoc sign flip below exists to
  repair, and it is applied only on this path. The flip is a valid
  repair for two cell types, where flipping one block negates the single
  pairwise term; it is **not** valid for three or more, where flipping
  one block negates only the terms touching it and can lower the
  objective. Retained so results computed before `"gauss-seidel"` became
  the default can be reproduced exactly.

  Note what the Gauss-Seidel guarantee covers: the **sign**, not
  solution quality. Under `objective = "sumcor"` the frozen-`sigma`
  sweep maximizes a surrogate rather than the objective itself. Its
  fixed points need not be stationary points of the ratio objective.
  Neither sweep dominates the other – measured across 8 configurations,
  Jacobi was ahead in one. Gauss-Seidel is the default because it cannot
  produce the pathology the sign repair existed to cover, not because it
  is the better optimizer.

- objective:

  `"sumcor"` (default) divides each slide's cross term by that slide's
  own score scales. `"sumcov"` fixes every scale at 1, giving the plain
  sum of kernel-smoothed cross-covariances – the gene-space counterpart
  of `runSkrCCA(objective = "sumcov")`.

## Value

Named list of weight matrices, each G x nCC.
