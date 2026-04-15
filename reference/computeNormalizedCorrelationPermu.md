# Compute Normalized Correlation for Permutation Results

Calculates the normalized correlation for each permutation, enabling
p-value calculation by comparing observed values to the null
distribution.

## Usage

``` r
computeNormalizedCorrelationPermu(object, tol = 1e-04)
```

## Arguments

- object:

  A `CoPro` object with permutation results from
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)

- tol:

  Tolerance for approximate SVD calculation (default: 1e-4)

## Value

The `CoPro` object with permutation normalized correlations stored in
`@normalizedCorrelationPermu`

## Details

The normalized correlation for each permutation is calculated using the
same formula as the observed data:

\$\$NC = \frac{s_1^T K\_{12} s_2}{\|\|s_1\|\| \cdot \|\|s_2\|\| \cdot
\|\|K\_{12}\|\|\_2}\$\$

where \\s_1\\ and \\s_2\\ are cell scores, \\K\_{12}\\ is the kernel
matrix, and \\\|\|K\_{12}\|\|\_2\\ is the spectral norm.

## Examples

``` r
if (FALSE) { # \dontrun{
# After running permutation testing
br <- computeNormalizedCorrelationPermu(br)

# Extract permutation values and calculate p-value
permu_values <- sapply(br@normalizedCorrelationPermu,
                       function(x) x$normalizedCorrelation[1])
observed <- max(getNormCorr(br)$normalizedCorrelation)

# One-sided p-value (testing if observed > permutation)
p_value <- mean(permu_values >= observed)
} # }
```
