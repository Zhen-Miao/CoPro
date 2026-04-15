# Run Permutation Test with Fair Sigma Selection

Performs permutation testing where BOTH observed and permuted data get
to optimize sigma selection. This is the statistically correct approach
that addresses inflated Type I error caused by sigma selection being
applied only to observed data.

## Usage

``` r
runSkrCCAPermu_FairSigma(
  object,
  nPermu = 100,
  sigma_values = NULL,
  permu_method = "bin",
  permu_which = "second_only",
  num_bins_x = 10,
  num_bins_y = 10,
  match_quantile = FALSE,
  maxIter = 200,
  tol = 1e-05,
  n_cores = 1,
  verbose = TRUE
)
```

## Arguments

- object:

  A CoPro object with CCA already computed via
  [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
  and normalized correlation computed via
  [`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md)

- nPermu:

  Number of permutations to run (default: 100)

- sigma_values:

  Vector of sigma values to test. If NULL, uses all sigma values from
  the original analysis (object@sigmaValues)

- permu_method:

  Method of permutation: "bin", "global", "pc", or "toroidal"

- permu_which:

  Which cell types to permute: "second_only", "both", "first_only"

- num_bins_x:

  Number of bins in x for bin-wise permutation

- num_bins_y:

  Number of bins in y for bin-wise permutation

- match_quantile:

  Whether to use quantile matching for bin permutation

- maxIter:

  Maximum iterations for CCA optimization

- tol:

  Convergence tolerance

- n_cores:

  Number of cores for parallel computation (not yet implemented)

- verbose:

  Whether to print progress messages

## Value

CoPro object with fair permutation results stored in:

- @skrCCAPermuOut: Best weights for each permutation

- @normalizedCorrelationPermu: Best ncorr for each permutation

- @fairSigmaPermu: List with sigma selected for each permutation

## Details

### The Sigma Selection Problem

In standard CoPro analysis, the observed data gets to choose the best
sigma (the one maximizing normalized correlation). However, permutation
data uses this SAME sigma, which may not be optimal for permuted data.
This asymmetry can inflate Type I error.

### The Solution

This function runs CCA at EACH sigma value for EACH permutation, then
selects the best sigma for that permutation. Both observed and permuted
data thus have equal opportunity to optimize sigma selection.

### Computational Cost

This is more computationally expensive (nPermu \* nSigma CCA runs
instead of nPermu runs), but provides statistically correct p-values.

## Examples

``` r
if (FALSE) { # \dontrun{
# After running standard CoPro analysis
br <- runSkrCCA(br, scalePCs = TRUE)
br <- computeNormalizedCorrelation(br)

# Run fair sigma permutation test
br <- runSkrCCAPermu_FairSigma(br, nPermu = 100,
                                permu_method = "toroidal")

# Calculate p-value
result <- calculate_pvalue(br)
} # }
```
