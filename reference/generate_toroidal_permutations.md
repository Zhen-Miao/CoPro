# Generate Toroidal Shift Permutation Indices

Shifts spatial coordinates in a toroidal (wrap-around) manner to break
cross-type coordination while approximately preserving within-type
spatial autocorrelation. Useful as one of several
autocorrelation-respecting nulls for permutation testing (see also the
sigma-aware `bin` null and Moran Spectral Randomization for irregular
tissue).

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

1.  Applying a random shift to all coordinates (wrapping at the torus
    period)

2.  Re-ranking cells by their shifted positions and matching back to the
    original ordering.

The wrap uses the torus *period*, not the coordinate extent. A lattice
with `n` columns at spacing `s` spans `(n - 1) * s`, so wrapping on the
span would identify the last column with the first: two distinct columns
land on the same shifted coordinate and are interleaved by the rank
match rather than translated. The period is the span plus the one gap
the span omits, estimated by the median gap between adjacent distinct
coordinates – exactly `s` on a regular lattice.

Shifts are drawn uniformly on the whole torus. Restricting them away
from zero would drop the group elements whose configurations sit closest
to the observed one, leaving a reference set that is not a group;
because the statistic varies smoothly with the shift, dropping them
pushes null p-values toward both extremes and inflates the rejection
rate.

Important caveats (the docstring previously over-claimed "perfect"
preservation):

- The position matching is by coordinate *rank*. On a regular lattice
  that is exactly a rigid torus translation: the shifted point set is
  the original point set, so matching sorted positions recovers the
  translation and every pairwise torus distance is preserved. On an
  irregular point cloud it is a monotone rearrangement that preserves
  pairwise distances only approximately, so within-type autocorrelation
  is approximately (not exactly) preserved.

- The torus-translation test (Harms et al. 2001) assumes spatial
  stationarity and periodic wrap-around. Gluing opposite tissue edges
  creates artificial neighbours at the seam and is biologically false
  for non-rectangular / hole-bearing tissue, and the family of distinct
  shifts is effectively small (a coarse null with limited power).

For these reasons toroidal should be benchmarked against, not preferred
over, the sigma-aware patch null and graph-based surrogates; let
calibration choose.

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
