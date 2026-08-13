# Compute Normalized Correlation (approximation)

This method calculates the normalized correlation between pairs of cell
types based on CCA weights and the respective kernel matrix. It divides
the bilinear statistic \\T = a' K b\\ by the whitened-Frobenius norm
\\\\R_x^{1/2} K_c R_y^{1/2}\\\_F\\, which is the null standard deviation
of \\T\\ when the score vectors carry within-type covariance
proportional to \\R_x\\ and \\R_y\\.

## Usage

``` r
computeNormalizedCorrelation(
  object,
  tol = 1e-04,
  calculationMode = "perSlide",
  normalizer = c("legacy", "unwhitened", "kernel", "variogram"),
  normalizerControl = list()
)

# S4 method for class 'CoPro'
computeNormalizedCorrelation(
  object,
  tol = 1e-04,
  calculationMode = "perSlide",
  normalizer = c("legacy", "unwhitened", "kernel", "variogram"),
  normalizerControl = list()
)

# S4 method for class 'CoProMulti'
computeNormalizedCorrelation(
  object,
  tol = 1e-04,
  calculationMode = "perSlide",
  normalizer = c("legacy", "unwhitened", "kernel", "variogram"),
  normalizerControl = list()
)
```

## Arguments

- object:

  A `CoPro` or `CoProMulti` object containing CCA results and kernel
  matrices.

- tol:

  tolerance for approximate SVD calculation

- calculationMode:

  (for CoProMulti only) either "perSlide" or "aggregate", for single
  slide analysis, it is ignored, with default value "perSlide".

- normalizer:

  Which whitening operators to use in the denominator; one of
  `"legacy"`, `"unwhitened"`, `"kernel"`, `"variogram"`. See details.

- normalizerControl:

  A named list tuning the `"variogram"` normalizer. `distType`,
  `xDistScale`, `yDistScale`, `zDistScale` must match what was passed to
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  (they are not stored on the object). `range` accepts a named numeric
  vector of per-cell-type ranges, in normalized distance units, to skip
  estimation. Remaining entries (`maxCells`, `nBins`, `maxLagQuantile`,
  `minCorrelation`, `minBins`, `lowerLimit`) tune the fit and the
  operator's truncation.

## Value

The object with the normalized correlation value between any pair of
cell types added as a new slot, `normalizedCorrelation`. The resolved
normalizer is attached as an attribute and can be read back with
[`getNormalizerInfo()`](https://zhen-miao.github.io/CoPro/reference/getNormalizerInfo.md).

## Choosing the normalizer

The whitening operators decide how the criterion behaves across the
sigma grid, and therefore which bandwidth `sigmaValueChoice` ends up at.

- `"legacy"` (default):

  The historical behaviour, kept so that stored results reproduce: use
  the matched-sigma within-type kernels when the object happens to
  contain them, and \\R = I\\ otherwise. Because
  [`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md)
  builds only cross-type kernels, this is normally \\\\K_c\\\_F\\; but
  it silently becomes the whitened norm on an object that has been
  through
  [`computeSelfKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSelfKernel.md).
  Which one applied is now reported, and recorded by
  [`getNormalizerInfo()`](https://zhen-miao.github.io/CoPro/reference/getNormalizerInfo.md).

- `"unwhitened"`:

  Force \\R_x = R_y = I\\, i.e. \\\\K_c\\\_F\\. This assumes spatially
  independent scores, so it under-counts noise at large sigma and biases
  the selected bandwidth upward.

- `"kernel"`:

  Force \\R = K(\sigma)\\, the matched-sigma self-kernel, with the unit
  diagonal restored. Errors if the self-kernels are absent rather than
  falling back. This over-counts noise at large sigma and biases the
  selected bandwidth downward, and it carries a second problem:
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md)
  normalizes self-distances by its own scaling factor rather than the
  one
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  recorded, so a self-kernel "at the same sigma" is generally at a
  different physical bandwidth from the cross-kernel. Provided as a
  diagnostic, not for analysis.

- `"variogram"`:

  Estimate one within-type autocorrelation range per cell type from the
  feature-averaged spatial autocorrelation of the PC scores, and use \\R
  = \exp(-d^2/2\ell^2)\\. The range is a property of the score field
  rather than of the bandwidth, so it is fitted once and reused across
  the grid. It is fitted from all PC columns, never from a fitted
  canonical score, which would leak the association under test into the
  denominator.

  **This flattens the null when the scores share one correlation length,
  and should be treated as a diagnostic otherwise.** On a real targeted
  panel the principal components within a cell type generally do *not*
  share one: measured per-component ranges spanned 0.3-4.1 cell-spacings
  on a colon MERFISH dataset, where the flat feature-average returned a
  sub-spacing range (so \\R \approx I\\, reproducing `"unwhitened"`)
  while pooling over the leading five components gave a range 6x longer.
  The selected bandwidth is steeply sensitive to that choice – sweeping
  the assumed range over 0.4-5 cell-spacings moved it from 3.4 to 0.3
  spacings on the same data. Check `getNormalizerInfo()$ranges` against
  what you believe about the tissue, supply `normalizerControl$range` if
  you have a better estimate, and prefer a permutation null when the
  selected bandwidth matters.

## See also

[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`computeBidirCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeBidirCorrelation.md),
[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md),
[`getNormalizerInfo()`](https://zhen-miao.github.io/CoPro/reference/getNormalizerInfo.md)

Other scores-and-correlation:
[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md),
[`getCCAObjective()`](https://zhen-miao.github.io/CoPro/reference/getCCAObjective.md),
[`getNormalizerInfo()`](https://zhen-miao.github.io/CoPro/reference/getNormalizerInfo.md)
