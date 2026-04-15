# Generate Toroidal Shift Permutation Indices

Shifts spatial coordinates in a toroidal (wrap-around) manner, perfectly
preserving spatial autocorrelation structure. This is useful for
permutation testing when you want to break cross-type coordination while
preserving within-type spatial patterns.

## Usage

``` r
generate_toroidal_permutations(location_data, n_permu = 100, seed = NULL)
```

## Arguments

- location_data:

  Data frame with x, y, cell_ID columns

- n_permu:

  Number of permutations to generate

- seed:

  Optional random seed for reproducibility

## Value

Matrix of permutation indices (n_cells x n_permu). Each column contains
a permutation of row indices that can be used to reorder cells.

## Details

The toroidal shift works by:

1.  Applying a random shift to all coordinates (wrapping at boundaries)

2.  Matching cells based on their new positions to original positions

This preserves ALL spatial autocorrelation within each cell type because
the relative positions of cells are unchanged - only their absolute
positions are shifted.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create example location data
loc_data <- data.frame(
  x = runif(100, 0, 10),
  y = runif(100, 0, 10),
  cell_ID = paste0("cell_", 1:100)
)

# Generate 100 toroidal permutations
perm_matrix <- generate_toroidal_permutations(loc_data, n_permu = 100)

# Apply first permutation
permuted_cells <- loc_data$cell_ID[perm_matrix[, 1]]
} # }
```
