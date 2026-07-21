# Read the step-down per-axis p-value table

Thin reader for the conditional step-down permutation test produced by
[`runSkrCCAPermu_Conditional()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_Conditional.md).

## Usage

``` r
calculate_pvalue_stepdown(object)
```

## Arguments

- object:

  A CoPro object with `@conditionalPermu` populated.

## Value

A data frame with one row per canonical axis: `CC_index`,
`observed_stat`, `observed_sigma`, `p_raw`, `p_stepdown`, and
`significant`. The number of significant axes, `alpha`, and `nPermu` are
attached as attributes.

## See also

[`runSkrCCAPermu_Conditional()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_Conditional.md)

## Examples

``` r
if (FALSE) { # \dontrun{
br <- runSkrCCAPermu_Conditional(br, nPermu = 200)
calculate_pvalue_stepdown(br)
} # }
```
