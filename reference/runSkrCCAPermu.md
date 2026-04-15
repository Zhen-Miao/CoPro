# Run Spatial CCA with Permutation Testing

Performs permutation testing to assess the significance of spatial
co-progression detected by CoPro. This generates a null distribution by
permuting cell assignments while optionally preserving spatial
structure.

## Usage

``` r
runSkrCCAPermu(
  object,
  tol = 1e-05,
  nPermu = 20,
  maxIter = 200,
  permu_method = "bin",
  permu_which = "second_only",
  num_bins_x = 10,
  num_bins_y = 10,
  match_quantile = FALSE,
  conservative = FALSE,
  n_cores = 1,
  verbose = TRUE
)
```

## Arguments

- object:

  A `CoPro` object with CCA already computed via
  [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)

- tol:

  Tolerance for CCA optimization convergence (default: 1e-5)

- nPermu:

  Number of permutations to run (default: 20). Increase to 100+ for
  publication-quality p-values.

- maxIter:

  Maximum iterations for CCA optimization (default: 200)

- permu_method:

  Method of permutation:

  - "bin" (default): Bin-wise permutation preserving local spatial
    structure

  - "global": Simple random permutation breaking all spatial structure

  - "pc": PC-space permutation (like DIALOGUE) - shuffles values within
    each PC dimension, breaking cell correlation while preserving PC
    distributions

  - "toroidal": Toroidal shift permutation - perfectly preserves spatial
    autocorrelation by shifting coordinates in wrap-around manner

- permu_which:

  Which cell types to permute:

  - "second_only" (default): Keep first cell type fixed, permute others

  - "both": Permute all cell types independently (more conservative)

  - "first_only": Keep others fixed, permute only the first cell type

- num_bins_x:

  Number of bins in x direction for bin-wise permutation (default: 10).
  Use
  [`diagnose_bin_distribution()`](https://zhen-miao.github.io/CoPro/reference/diagnose_bin_distribution.md)
  to choose appropriate values. **More bins = better preserve local
  structure = lower FPR.**

- num_bins_y:

  Number of bins in y direction for bin-wise permutation (default: 10).
  **More bins = better preserve local structure = lower FPR.**

- match_quantile:

  Logical. If TRUE and `permu_method = "bin"`, matches cells between
  tiles based on their relative (quantile) x/y positions. This better
  preserves within-tile spatial autocorrelation structure. Default:
  FALSE. **Setting TRUE helps reduce FPR by better preserving spatial
  structure.**

- conservative:

  Logical. If TRUE, automatically uses settings that better preserve
  spatial autocorrelation to reduce false positive rate:
  `permu_which = "second_only"`, `num_bins_x = 15`, `num_bins_y = 15`,
  `match_quantile = TRUE`. Default: FALSE.

- n_cores:

  Number of cores for parallel computation (default: 1). Set to higher
  values to speed up permutation testing. Use
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  to find available cores. Parallelization uses the `parallel` package.

- verbose:

  Whether to print progress messages (default: TRUE)

## Value

CoPro object with permutation results stored in `@skrCCAPermuOut`

## Details

### Permutation Methods

The function supports three permutation methods:

**"global"**: Simple random shuffling of cells. This breaks ALL spatial
structure and tests against a null of complete spatial randomness.

**"bin"** (default): Bin-wise shuffling that preserves local spatial
structure. This tests against a null where cells have spatial
autocorrelation within their type, but no coordination across types.

**"pc"**: PC-space permutation (like DIALOGUE). Shuffles values within
each PC dimension across cells, breaking cell-to-cell correlation while
preserving the marginal distribution of each PC. This is the same
approach used by DIALOGUE's internal significance testing. Use this to
compare with DIALOGUE's conservative behavior.

**"toroidal"**: Toroidal (wrap-around) shift permutation. Shifts all
cells' coordinates by a random amount, wrapping at boundaries. This
PERFECTLY preserves spatial autocorrelation within each cell type
because relative positions are unchanged. Best choice for reducing FPR
when spatial autocorrelation is strong.

### Which Cell Types to Permute

The `permu_which` parameter controls which cell types are permuted:

**"second_only"** (default): Keep the first cell type FIXED, permute all
others. This is the standard approach for permutation testing. For two
cell types A and B, only B is permuted while A stays fixed.

**"both"**: Permute ALL cell types independently. Both A and B are
shuffled with different random permutations. This breaks MORE structure
and may lead to higher FPR. Use "second_only" for better FPR control.

**"first_only"**: Keep second+ cell types FIXED, permute only the first.
Useful if you want to test from the opposite direction.

### Controlling False Positive Rate

High FPR under null simulations typically means the permutation is
**breaking too much spatial structure**, making the null distribution
too low, so that even random observed values appear significant.

To **reduce FPR**, you need to **better preserve spatial
autocorrelation**:

1.  **Use MORE bins** (not fewer): `num_bins_x = 15, num_bins_y = 15`
    preserves more local structure than the default 10x10.

2.  **Enable quantile matching**: `match_quantile = TRUE` matches cells
    by their relative positions within tiles, better preserving
    within-tile structure.

3.  **Only permute one cell type**: `permu_which = "second_only"`
    (default) keeps one cell type fixed, preserving more structure than
    "both".

4.  **Avoid global permutation**: `permu_method = "global"` breaks ALL
    spatial structure and will likely have high FPR.

The key insight: if permutation doesn't adequately preserve within-type
spatial autocorrelation, the null distribution will be too low, leading
to inflated significance.

## See also

[`computeNormalizedCorrelationPermu`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelationPermu.md)
to compute normalized correlation from permutation results
[`diagnose_bin_distribution`](https://zhen-miao.github.io/CoPro/reference/diagnose_bin_distribution.md)
to check bin distribution

## Examples

``` r
if (FALSE) { # \dontrun{
# After running standard CoPro analysis
br <- runSkrCCA(br, scalePCs = TRUE)
br <- computeNormalizedCorrelation(br)

# Standard permutation (only permute second cell type)
br <- runSkrCCAPermu(br, nPermu = 100, permu_method = "bin",
                     permu_which = "second_only")

# Conservative permutation (lower FPR - better preserves spatial structure)
br <- runSkrCCAPermu(br, nPermu = 100, conservative = TRUE)

# Manual conservative settings: more bins + quantile matching
br <- runSkrCCAPermu(br, nPermu = 100, permu_method = "bin",
                     permu_which = "second_only",
                     num_bins_x = 15, num_bins_y = 15,
                     match_quantile = TRUE)

# PC-space permutation (like DIALOGUE) - shuffles within PC dimensions
br <- runSkrCCAPermu(br, nPermu = 100, permu_method = "pc",
                     permu_which = "second_only")

br <- computeNormalizedCorrelationPermu(br)

# Calculate p-value
observed <- max(getNormCorr(br)$normalizedCorrelation)
permu_values <- sapply(br@normalizedCorrelationPermu,
                       function(x) x$normalizedCorrelation[1])
p_value <- mean(permu_values >= observed)
} # }
```
