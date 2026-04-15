# Compute Self-Distance Matrices for Multiple Cell Types

This function computes within-cell-type distance matrices for each cell
type when multiple cell types are present. Unlike the standard
computeDistance which only computes cross-type distances for multiple
cell types, this function computes self-distances (within each cell
type).

## Usage

``` r
computeSelfDistance(
  object,
  distType = c("Euclidean2D", "Euclidean3D", "Morphology-Aware"),
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  truncateLowDist = TRUE,
  verbose = TRUE,
  overwrite = FALSE
)

# S4 method for class 'CoProSingle'
computeSelfDistance(
  object,
  distType = c("Euclidean2D", "Euclidean3D", "Morphology-Aware"),
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  truncateLowDist = TRUE,
  verbose = TRUE,
  overwrite = FALSE
)

# S4 method for class 'CoProMulti'
computeSelfDistance(
  object,
  distType = c("Euclidean2D", "Euclidean3D", "Morphology-Aware"),
  xDistScale = 1,
  yDistScale = 1,
  zDistScale = 1,
  normalizeDistance = TRUE,
  truncateLowDist = TRUE,
  verbose = TRUE,
  overwrite = FALSE
)
```

## Arguments

- object:

  A `CoPro` object with multiple cell types

- distType:

  Type of distance to compute: "Euclidean2D", "Euclidean3D", or
  "Morphology-Aware"

- xDistScale:

  Scale for x distance, default = 1

- yDistScale:

  Scale for y distance, default = 1

- zDistScale:

  Scale for z distance, default = 1

- normalizeDistance:

  Whether to normalize distance? The normalization will make sure that
  the 0.01% cell-cell distance will become 0.01, thus ensuring
  consistent scaling across cell types. Default = TRUE

- truncateLowDist:

  Whether to truncate small distances so that cells that are nearly
  overlapping do not have super small distances. Default = TRUE.

- verbose:

  Whether to print info about the computation progress

- overwrite:

  Whether to overwrite existing distance matrices. If FALSE, will add
  self-distance matrices to existing cross-type distances. Default =
  FALSE

## Value

`CoPro` object with self-distance matrices added to the distances slot

## Examples

``` r
if (FALSE) { # \dontrun{
# Assume you have a CoPro object with multiple cell types
# First compute cross-type distances
object <- computeDistance(object)

# Then add self-distances
object <- computeSelfDistance(object)

# Now you have both cross-type and self-type distance matrices
} # }
```
