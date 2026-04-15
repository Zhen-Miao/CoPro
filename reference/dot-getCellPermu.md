# Generate Cell Permutation Indices

Internal function to generate permutation indices for each cell type.

## Usage

``` r
.getCellPermu(
  object,
  permu_method,
  nPermu,
  cts,
  permu_which = "second_only",
  num_bins_x = 10,
  num_bins_y = 10,
  match_quantile = FALSE
)
```

## Arguments

- object:

  A CoPro object

- permu_method:

  "global", "bin", "pc", or "toroidal"

- nPermu:

  Number of permutations

- cts:

  Cell types to permute

- permu_which:

  Which cell types to permute: "second_only", "both", or "first_only"

- num_bins_x:

  Number of bins in x for bin-wise permutation

- num_bins_y:

  Number of bins in y for bin-wise permutation

- match_quantile:

  Logical. If TRUE and permu_method="bin", matches cells between tiles
  based on their relative (quantile) positions to better preserve
  within-tile spatial structure. Default: FALSE.

## Value

List of permutation matrices (for "global"/"bin"/"toroidal") or list of
permuted PC matrices (for "pc"), one per cell type

## Details

Permutation strategy is controlled by `permu_which`:

- "second_only" (default): Keep first cell type FIXED, permute others

- "both": Permute ALL cell types independently

- "first_only": Keep second+ cell types FIXED, permute only first

The default ("second_only") is the standard approach for permutation
testing: we test whether the relationship between cell types is stronger
than expected by chance, while keeping one cell type as reference.
