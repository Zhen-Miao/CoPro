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
  step_size = 1,
  verbose = TRUE,
  sweep = c("gauss-seidel", "jacobi"),
  objective = c("sumcor", "sumcov")
)
```

## Arguments

- C_self_slide:

  Named list of per-slide self-covariance matrices or matrix-free
  operators.

- C_cross_slide:

  Named list of per-slide cross-covariance matrices or matrix-free
  operators.

- slides:

  Character vector of slide IDs.

- cell_types:

  Character vector of cell type names.

- max_iter:

  Maximum iterations (default 3000). Must be \>= 1.

- tol:

  Convergence tolerance on max weight change (default 1e-6).

- step_size:

  Step size for damped power iteration (default 1). Must be in (0, 1\].
  Lower values blend the new iterate with the previous one:
  `w_new = normalize((1 - step_size) * w_old + step_size * w_update)`,
  which damps oscillation when pure power iteration fails to converge.

- verbose:

  Print progress every 500 iterations (default TRUE).

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
  sweep maximizes a surrogate rather than the objective itself, so which
  local optimum each sweep reaches is data-dependent and neither sweep
  dominates the other – measured across 8 configurations, Jacobi was
  ahead in one. Gauss-Seidel is the default because it cannot produce
  the pathology the sign repair existed to cover, not because it is the
  better optimizer.

- objective:

  `"sumcor"` (default) divides each slide's cross term by that slide's
  own score scales. `"sumcov"` fixes every scale at 1, giving the plain
  sum of kernel-smoothed cross-covariances – the gene-space counterpart
  of
  [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)'s
  default.

## Value

Named list of weight vectors, one per cell type (each a G x 1 matrix).
