# Compute sigma-aware bin counts for spatial permutation

Internal helper that chooses the number of bins for bin-wise spatial
permutation from the kernel bandwidth `sigma` instead of a fixed
hard-coded grid. The target patch side is `2 * sigma` on the normalized
distance scale, so each patch is large enough to preserve within-type
spatial autocorrelation (length scale `~sigma`) while shuffling whole
patches still breaks longer-range cross-type coordination. A hard-coded
grid (e.g. the historical 10x10) ignores both the tissue extent and
`sigma`, which is the root cause of mis-calibrated bin-wise permutation.

## Usage

``` r
.sigmaAwareBins(object, sigma, min_bins = 2L, verbose = TRUE)
```

## Arguments

- object:

  A CoPro object (uses `@locationDataSub` and `@distances`).

- sigma:

  Kernel bandwidth (numeric, on the normalized distance scale).

- min_bins:

  Minimum number of bins per axis (default 2).

- verbose:

  Whether to message the chosen grid (default TRUE).

## Value

A list with integer `num_bins_x`, `num_bins_y`, and the recovered
`scale_factor`.
