# Calculate P-value from Permutation Results

Helper function to calculate p-value from permutation testing results.

## Usage

``` r
calculate_pvalue(object, cc_index = 1, alternative = "greater")
```

## Arguments

- object:

  A CoPro object with permutation results

- cc_index:

  Which canonical correlation component to use (default: 1)

- alternative:

  Direction of test: "greater" (default), "less", or "two.sided"

## Value

List with the Phipson & Smyth (2010) permutation p-value (`p_value`,
never exactly zero), the Monte-Carlo floor
`mc_floor = 1 / (n_permu + 1)`, the observed value, and the permutation
distribution.

## References

Phipson B, Smyth GK (2010). Permutation P-values should never be zero.
*Stat Appl Genet Mol Biol* 9:Article39.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- calculate_pvalue(br, cc_index = 1, alternative = "greater")
print(paste("P-value:", result$p_value))
} # }
```
