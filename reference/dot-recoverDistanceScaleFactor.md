# Recover the raw-to-normalized distance scaling factor

Internal helper that recovers the single scalar mapping raw
spatial-coordinate distances to the normalized distance scale on which
`sigmaValues` are defined. `computeDistance(normalizeDistance = TRUE)`
multiplies every stored distance by one constant, so for any entry that
was not affected by low-distance truncation the ratio (normalized
distance) / (raw Euclidean distance) equals that constant. We probe
random off-diagonal entries and take the median ratio, which is robust
to the small fraction of truncated entries and avoids copying the full
distance matrix.

## Usage

``` r
.recoverDistanceScaleFactor(object, n_probe = 200L)
```

## Arguments

- object:

  A CoPro object with `@distances` populated.

- n_probe:

  Number of random cell pairs to probe (default 200).

## Value

Numeric scaling factor (normalized = raw \* factor), or `NA_real_` if it
cannot be recovered.

## Details

This assumes isotropic coordinates (`xDistScale = yDistScale`, the
default). Under anisotropic scaling the recovered value is a
representative, not exact, factor; the `[sigma, 4*sigma]` guardrail in
[`.sigmaAwareBins()`](https://zhen-miao.github.io/CoPro/reference/dot-sigmaAwareBins.md)
will flag a grossly mis-sized grid. When `normalizeDistance = FALSE` the
recovered factor is simply the coordinate scale (1 by default), which is
also correct because `sigma` is then on the raw coordinate scale.
