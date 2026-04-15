# Changelog

## CoPro 1.0.0

### Major changes

- First public release accompanying the CoPro manuscript.
- Added
  [`computeRegressionGeneScores()`](https://zhen-miao.github.io/CoPro/reference/computeRegressionGeneScores.md)
  for regression-based gene weights, which avoids collinearity issues
  present in PCA back-projection.
- Added
  [`copro_download_data()`](https://zhen-miao.github.io/CoPro/reference/copro_download_data.md)
  for easy download of example datasets via `piggyback`.
- Corrected PCA back-projection formula in
  [`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md):
  gene weights now use `1/sdev` instead of `sdev`.
- Changed
  [`transfer_scores()`](https://zhen-miao.github.io/CoPro/reference/transfer_scores.md)
  default `gs_weight_threshold` from 0.005 to 0.
- Added support for regression-based score transfer via
  `getTransferCellScores(..., gene_score_type = "regression")`.

### New vignettes

- “Within-cell-type spatial patterns (Organoid)” – single cell type
- “Cross-cell-type co-progression (Brain MERFISH)” – two cell types
- “Cross-cell-type co-progression with orthogonal axes (Colon Day 3)”
- “Multi-slide analysis and score transfer (Colon Day 9)”
- “Supervised detection of spatial gradients (Kidney)”

### Documentation

- Full pkgdown website with grouped function reference.
- Cleaned reproducibility scripts in `scripts/` directory.
- Example datasets available via
  [`copro_download_data()`](https://zhen-miao.github.io/CoPro/reference/copro_download_data.md).

## CoPro 0.6.1

- Internal refactoring of cell scores, distance matrices, and kernel
  matrices for cleaner API.
- Added
  [`getCellScoresInSitu()`](https://zhen-miao.github.io/CoPro/reference/getCellScoresInSitu.md),
  [`getDistMat()`](https://zhen-miao.github.io/CoPro/reference/getDistMat.md),
  [`getKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/getKernelMatrix.md)
  accessor functions.
- Added self-correlation methods:
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md),
  [`computeSelfKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSelfKernel.md),
  [`computeSelfBidirCorr()`](https://zhen-miao.github.io/CoPro/reference/computeSelfBidirCorr.md).
- Added
  [`computeBidirCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeBidirCorrelation.md)
  for bidirectional correlation.
- Added
  [`testGeneGLM()`](https://zhen-miao.github.io/CoPro/reference/testGeneGLM.md)
  for gene-level GLM testing.
- Added
  [`resample_spatial()`](https://zhen-miao.github.io/CoPro/reference/resample_spatial.md)
  and
  [`generate_toroidal_permutations()`](https://zhen-miao.github.io/CoPro/reference/generate_toroidal_permutations.md).
