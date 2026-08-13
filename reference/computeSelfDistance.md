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
  distType = NULL,
  xDistScale = NULL,
  yDistScale = NULL,
  zDistScale = NULL,
  normalizeDistance = NULL,
  normalizeMethod = NULL,
  normalizeTarget = NULL,
  truncateLowDist = NULL,
  verbose = TRUE,
  overwrite = FALSE
)

# S4 method for class 'CoProSingle'
computeSelfDistance(
  object,
  distType = NULL,
  xDistScale = NULL,
  yDistScale = NULL,
  zDistScale = NULL,
  normalizeDistance = NULL,
  normalizeMethod = NULL,
  normalizeTarget = NULL,
  truncateLowDist = NULL,
  verbose = TRUE,
  overwrite = FALSE
)

# S4 method for class 'CoProMulti'
computeSelfDistance(
  object,
  distType = NULL,
  xDistScale = NULL,
  yDistScale = NULL,
  zDistScale = NULL,
  normalizeDistance = NULL,
  normalizeMethod = NULL,
  normalizeTarget = NULL,
  truncateLowDist = NULL,
  verbose = TRUE,
  overwrite = FALSE
)
```

## Arguments

- object:

  A `CoPro` object with multiple cell types

- distType:

  Type of distance to compute: `"Euclidean2D"`, `"Euclidean3D"`, or
  `"Morphology-Aware"`. `NULL` (default) inherits the recorded geometry,
  falling back to the coordinate columns present.

- xDistScale:

  Scale for x distance. `NULL` (default) inherits the recorded geometry,
  falling back to 1.

- yDistScale:

  Scale for y distance. `NULL` (default) inherits the recorded geometry,
  falling back to 1.

- zDistScale:

  Scale for z distance. `NULL` (default) inherits the recorded geometry,
  falling back to 1.

- normalizeDistance:

  How to scale the self-distances. One of:

  `NULL` (default)

  :   Inherit the recorded geometry, falling back to `FALSE` as of CoPro
      1.2.0.

  `FALSE`

  :   No scaling; distances stay in the coordinate units of
      `locationData`.

  `TRUE`

  :   Rescale so that the reference cell-cell distance becomes
      `normalizeTarget`. When the object already carries a scale factor,
      that factor is reused rather than re-derived from the within-type
      blocks, so cross-type and within-type distances stay on one unit;
      the substitution is reported. Pass `overwrite = TRUE` to discard
      the existing blocks and re-derive.

  `"inherit"`

  :   Reuse the scaling factor
      [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
      recorded on the object, and error if it recorded none. Differs
      from `TRUE` only on an object with nothing pinned, where `TRUE`
      derives a unit from the within-type blocks; use it to assert that
      the cross-type scale must exist – notably when the self-kernels
      become whitening operators in
      [`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md).

- normalizeMethod:

  How the reference distance is estimated when
  `normalizeDistance = TRUE` and no scale factor is already pinned on
  the object. `"global"` uses the median nearest-neighbor distance over
  all cells, ignoring type labels, which gives the same unit here as in
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md);
  `"spacing"` measures the within-type blocks instead; `"percentile"`
  reproduces the pre-1.2.0 behavior (the minimum, across blocks, of a
  low quantile of pairwise distances). Not available for
  `distType = "Morphology-Aware"`, which falls back to `"percentile"`.
  `NULL` (default) inherits the recorded geometry, falling back to
  `"global"`. Ignored when `normalizeDistance` is `FALSE` or
  `"inherit"`.

- normalizeTarget:

  Numeric scalar. The value the reference distance is rescaled to.
  `NULL` (default) inherits the recorded geometry, falling back to 0.01.

- truncateLowDist:

  Whether to truncate small distances so that cells that are nearly
  overlapping do not have super small distances. `NULL` (default)
  inherits the recorded geometry, falling back to `TRUE`.

- verbose:

  Whether to print info about the computation progress

- overwrite:

  Whether to overwrite existing distance matrices. If FALSE, will add
  self-distance matrices to existing cross-type distances. Default =
  FALSE

## Value

`CoPro` object with self-distance matrices added to the distances slot

## Details

`computeSelfDistance()` is normally additive: it adds within-type blocks
to an object whose cross-type blocks
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
already built. Its distance arguments therefore default to `NULL`,
meaning "use whatever the existing distances were built on"
([`getDistanceGeometry()`](https://zhen-miao.github.io/CoPro/reference/getDistanceGeometry.md));
supplying a value that contradicts that record is an error rather than a
silent second basis. When nothing is recorded – a self-distance-only
workflow – the package defaults apply.

When `normalizeDistance = TRUE` and the object already carries a scale
factor, that factor is reused rather than re-derived from the
within-type blocks, so cross-type and within-type distances stay on one
unit. Pass `overwrite = TRUE` to discard the existing blocks and
re-derive.

## See also

[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`getDistanceGeometry()`](https://zhen-miao.github.io/CoPro/reference/getDistanceGeometry.md)

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
