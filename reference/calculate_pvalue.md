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

List with p-value, observed value, and permutation distribution

## Examples

``` r
if (FALSE) { # \dontrun{
result <- calculate_pvalue(br, cc_index = 1, alternative = "greater")
print(paste("P-value:", result$p_value))
} # }
```
