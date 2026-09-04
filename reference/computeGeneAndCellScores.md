# computeGeneAndCellScores

computeGeneAndCellScores

## Usage

``` r
computeGeneAndCellScores(object)

# S4 method for class 'CoPro'
computeGeneAndCellScores(object)

# S4 method for class 'CoProMulti'
computeGeneAndCellScores(object)
```

## Arguments

- object:

  A `CoPro` or `CoProMulti` object containing CCA results and kernel
  matrices.

## Value

A `CoPro` or `CoProMulti` object with gene and cell scores computed

## Details

Gene scores are the PCA back-projection
\\R\\\mathrm{diag}(1/\mathrm{sdev})\\w\\, so they carry the units of the
loading matrix \\R\\: **weights per standardized gene**, not per unit of
raw expression. That was already true under pooled preprocessing. Under
the `CoProMulti` default (`computePCA(center_per_slide = TRUE)`) the
standardization is per (slide, cell type), so when the stored slide
scales differ there is no single equivalent coefficient vector in raw
units; the per-slide maps are kept on `pcaGlobal[[ct]]$slideCenter` and
`$slideScale`. Standard deviations below `1e-5` (and non-finite values)
are treated as one during this inverse scaling, adding a final guard
against unstable weights from a numerically rank-deficient PCA.

For gene weights that do not depend on the PCA coordinate system, prefer
[`computeRegressionGeneScores()`](https://zhen-miao.github.io/CoPro/reference/computeRegressionGeneScores.md),
which regresses cell scores on expression directly. Those are also more
robust to `nPCA` and reproduce better across replicates.

## See also

[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md),
[`computeRegressionGeneScores()`](https://zhen-miao.github.io/CoPro/reference/computeRegressionGeneScores.md),
[`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md),
[`getCellScoresInSitu()`](https://zhen-miao.github.io/CoPro/reference/getCellScoresInSitu.md)

Other scores-and-correlation:
[`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md),
[`getCCADiagnostics()`](https://zhen-miao.github.io/CoPro/reference/getCCADiagnostics.md),
[`getCCAObjective()`](https://zhen-miao.github.io/CoPro/reference/getCCAObjective.md),
[`getNormalizerInfo()`](https://zhen-miao.github.io/CoPro/reference/getNormalizerInfo.md)
