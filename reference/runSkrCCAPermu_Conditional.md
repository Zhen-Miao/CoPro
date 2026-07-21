# Conditional (sequential step-down) permutation test across canonical axes

Performs a permutation test that controls BOTH the sigma-selection
multiplicity (via a fair-sigma max-statistic, as in
[`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md))
AND the canonical-axis multiplicity (via a sequential conditional
step-down test). This is the statistically correct way to assign
p-values to CC1, CC2, ... jointly, and it supersedes computing a
separate per-axis p-value from a full multi-component permutation, which
is anti-conservative for axes `k >= 2` when CC1 is strong.

## Usage

``` r
runSkrCCAPermu_Conditional(
  object,
  nPermu = 100,
  sigma_values = NULL,
  permu_method = "bin",
  permu_which = "second_only",
  num_bins_x = NULL,
  num_bins_y = NULL,
  match_quantile = FALSE,
  alpha = 0.05,
  maxIter = 200,
  tol = 1e-05,
  verbose = TRUE
)
```

## Arguments

- object:

  A CoPro object with
  [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
  and
  [`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md)
  already run.

- nPermu:

  Number of permutations (default 100).

- sigma_values:

  Candidate bandwidths for the fair-sigma maximum. `NULL` uses all
  `object@sigmaValues` that have kernels and weights.

- permu_method:

  Permutation null: "bin" (default), "global", "pc", or "toroidal". See
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md).

- permu_which:

  Which cell types to permute: "second_only" (default), "both", or
  "first_only".

- num_bins_x, num_bins_y:

  Bin grid for `permu_method = "bin"`. Default `NULL` is sigma-aware
  (see
  [`.sigmaAwareBins()`](https://zhen-miao.github.io/CoPro/reference/dot-sigmaAwareBins.md));
  the same grid is shared across the sigma sweep. Pass integers to
  override.

- match_quantile:

  Whether to use quantile matching for bin permutation.

- alpha:

  Family-wise significance level for the step-down rule (default 0.05).

- maxIter, tol:

  Optimization controls passed to the axis optimizer.

- verbose:

  Whether to print progress and a summary (default TRUE).

## Value

The CoPro object with results stored in the `@conditionalPermu` slot, a
list whose `per_axis` element is a data frame of `CC_index`,
`observed_stat`, `observed_sigma`, `p_raw`, `p_stepdown`, `mc_floor`,
and `significant`, plus the full null matrices for diagnostics. Use
[`calculate_pvalue_stepdown()`](https://zhen-miao.github.io/CoPro/reference/calculate_pvalue_stepdown.md)
to read the per-axis table.

## Details

### Why a conditional test

The observed CC2 is obtained by deflating the observed (real) CC1
direction and optimizing the residual. A naive permutation p-value
compares this to the CC2 of fully re-optimized permutations, where each
permutation deflates its OWN leading axis. Those two statistics are
produced by different operators and are not exchangeable under the null,
which biases the CC2 p-value downward (too many false positives). The
conditional test instead deflates every permutation by the SAME fixed
observed CC1..CC(k-1) directions, so the observed and null axis-`k`
statistics share one operator and are exchangeable.

### The statistic and step-down rule

For axis `k`, the statistic is the fair-sigma maximum over the candidate
bandwidths of the residual normalized correlation after deflating the
fixed observed CC1..CC(k-1) directions. The observed value is read from
the stored `normalizedCorrelation` (so it matches
[`getNormCorr()`](https://zhen-miao.github.io/CoPro/reference/getNormalizedCorrelation.md)
exactly); each permutation re-optimizes the residual leading axis on
permuted PCs. The raw p-value uses the Phipson & Smyth (2010) estimator,
which counts the observed configuration as one admissible permutation so
the p-value is never exactly zero:

\$\$p\_{\mathrm{raw}}(k) = \frac{1 + \\\\\mathrm{perm}\_k \ge
\mathrm{obs}\_k\\}{1 + m},\$\$

with Monte-Carlo floor `1 / (m + 1)` (reported as `mc_floor`). Closed
step-down control of the family-wise error rate across ordered axes uses

\$\$p\_{\mathrm{stepdown}}(k) = \max\_{j \le k}
p\_{\mathrm{raw}}(j),\$\$

and testing stops at the first axis with `p_stepdown > alpha`; that axis
and all later ones are declared non-significant. No Bonferroni factor is
needed. This is the closed/fixed-sequence test of Marcus, Peritz &
Gabriel (1976) and the permutation "test of canonical axes" of Legendre,
Oksanen & ter Braak (2011).

### Relationship to data residualization

Deflating `Y` by the fixed observed weight directions is algebraically
identical to the Freedman-Lane / ter Braak residualization that removes
the observed lower canonical *variates* from the data before recomputing
the cross-product, whenever the PCs are whitened (`scalePCs = TRUE`, so
`X^T X = c I`): `(I - u u^T) Y (I - v v^T) = Y - (u^T Y v)\, u v^T`. The
fair-sigma maximum over the bandwidth family is the Westfall-Young
(1993) maxT procedure.

## References

Phipson B, Smyth GK (2010). Permutation P-values should never be zero.
*Stat Appl Genet Mol Biol* 9:Article39. Legendre P, Oksanen J, ter Braak
CJF (2011). Testing the significance of canonical axes in redundancy
analysis. *Methods Ecol Evol* 2:269-277. Westfall PH, Young SS (1993).
*Resampling-based Multiple Testing*.

## See also

[`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md)
(the CC1 / sigma-multiplicity case),
[`calculate_pvalue_stepdown()`](https://zhen-miao.github.io/CoPro/reference/calculate_pvalue_stepdown.md)

## Examples

``` r
if (FALSE) { # \dontrun{
br <- runSkrCCA(br, scalePCs = TRUE, nCC = 3)
br <- computeNormalizedCorrelation(br)
br <- runSkrCCAPermu_Conditional(br, nPermu = 200, permu_method = "bin")
calculate_pvalue_stepdown(br)
} # }
```
