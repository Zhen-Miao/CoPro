# CoPro 1.0.0

## Major changes

* First public release accompanying the CoPro manuscript.
* Added `computeRegressionGeneScores()` for regression-based gene weights,
  which avoids collinearity issues present in PCA back-projection.
* Added `copro_download_data()` for easy download of example datasets
  via `piggyback`.
* Corrected PCA back-projection formula in `computeGeneAndCellScores()`:
  gene weights now use `1/sdev` instead of `sdev`.
* Changed `transfer_scores()` default `gs_weight_threshold` from 0.005 to 0.
* Added support for regression-based score transfer via
 `getTransferCellScores(..., gene_score_type = "regression")`.

## New vignettes

* "Within-cell-type spatial patterns (Organoid)" -- single cell type
* "Cross-cell-type co-progression (Brain MERFISH)" -- two cell types
* "Cross-cell-type co-progression with orthogonal axes (Colon Day 3)"
* "Multi-slide analysis and score transfer (Colon Day 9)"
* "Supervised detection of spatial gradients (Kidney)"

## Documentation

* Full pkgdown website with grouped function reference.
* Cleaned reproducibility scripts in `scripts/` directory.
* Example datasets available via `copro_download_data()`.

# CoPro 0.6.1

* Internal refactoring of cell scores, distance matrices, and kernel
  matrices for cleaner API.
* Added `getCellScoresInSitu()`, `getDistMat()`, `getKernelMatrix()`
  accessor functions.
* Added self-correlation methods: `computeSelfDistance()`,
  `computeSelfKernel()`, `computeSelfBidirCorr()`.
* Added `computeBidirCorrelation()` for bidirectional correlation.
* Added `testGeneGLM()` for gene-level GLM testing.
* Added `resample_spatial()` and `generate_toroidal_permutations()`.
