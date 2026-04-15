# Match Cells by Within-Tile Quantile Position

Internal helper to match cells from target tile to original tile based
on their relative (quantile) positions within each tile. This preserves
spatial structure better than random sampling.

## Usage

``` r
.match_by_quantile_position(orig_points, candidate_points, n_points)
```

## Arguments

- orig_points:

  Data frame of original bin cells with x, y columns

- candidate_points:

  Data frame of candidate cells to sample from

- n_points:

  Number of cells to sample

## Value

Integer vector of indices into candidate_points
