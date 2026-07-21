# Fit a single conditional canonical axis (internal kernel)

Internal numeric kernel for the sequential step-down permutation test.
Given (possibly permuted) PC matrices and the FIXED observed lower-axis
weight directions, it returns the leading canonical axis of the residual
after those lower directions have been deflated, together with that
axis' normalized correlation.

## Usage

``` r
.fitConditionalAxis(
  PCmats,
  flat_kernels,
  sigma,
  cts,
  W_lower = NULL,
  k_minus_1 = 0,
  Y_resi = NULL,
  kernel_info = NULL,
  maxIter = 200,
  tol = 1e-05
)
```

## Arguments

- PCmats:

  Named list of (possibly permuted) PC matrices, one per cell type.

- flat_kernels:

  Flat kernel list from the CoPro object.

- sigma:

  Kernel bandwidth (numeric).

- cts:

  Cell types of interest.

- W_lower:

  Named list of observed weight matrices whose first `k_minus_1` columns
  are the fixed deflation directions; ignored when `k_minus_1 = 0`.

- k_minus_1:

  Number of lower axes to deflate (0 for the first axis).

- Y_resi:

  Optional precomputed `compute_Y_resi()` structure for this (PCmats,
  sigma); recomputed when `NULL` and `k_minus_1 >= 1`.

- kernel_info:

  Optional precomputed kernel and normalizer information from
  `.get_ncorr_kernel_info()`.

- maxIter, tol:

  Optimization controls.

## Value

List with `w` (named list of 1-column weight matrices for axis k) and
`ncorr` (its normalized correlation).

## Details

For axis `k = 1` (no lower directions; `k_minus_1 = 0`) this is exactly
the first-component optimization used by
[`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md),
so the `k = 1` conditional test reproduces the fair-sigma CC1 test
bit-for-bit.

For `k >= 2` it deflates the observed CC1..CC(k-1) directions from the
cross-product `Y = t(X_i) K_ij X_j` in feature (PC) space, using the
SAME fixed observed directions on every permutation, then optimizes the
leading residual component. Deflating each permutation by the same
observed directions (rather than by that permutation's own leading axis)
is what makes the higher-axis null exchangeable with the observed
statistic and removes the anti-conservative bias of the naive per-axis
permutation p-value.

Under whitened PCs (`scalePCs = TRUE`, so `X^T X = c I`) this Y-space
deflation is algebraically identical to the Freedman-Lane / ter Braak
residualization that removes the observed lower canonical variates from
the data before recomputing `Y`, because
`(I - u u^T) Y (I - v v^T) = Y - (u^T Y v)\, u v^T`.

The expensive `compute_Y_resi()` can be computed once per (permutation,
sigma) and reused across axes by passing it as `Y_resi`; deflation does
not mutate the supplied structure (copy-on-modify), so the same `Y_resi`
is safe to reuse for every `k`.
