# CoPro 1.1.0

## New features

* Added `asCoProSingle()` and `asCoProMulti()` S4 generics for one-call
  coercion from `SingleCellExperiment` and `Seurat` objects into CoPro
  objects. Conversions are gated on their respective packages being
  installed and delegate to the existing `newCoProSingle()` /
  `newCoProMulti()` constructors so validation stays single-sourced.
* Exposed `normalizeTarget` argument on `computeDistance()` for users
  who want to control the target value that the low-percentile cell-cell
  distance is rescaled to. Default preserves existing behavior.

## User experience

* `newCoProSingle()` / `newCoProMulti()` now reject `NA`, `NaN`, and `Inf`
  values in `normalizedData` at construction time with an informative
  error instead of producing cryptic downstream failures.
* `newCoProMulti()` now warns when `metaData` already contains a
  `slideID` column that matches the supplied `slideID` (instead of
  silently overwriting).
* `locationData` column standardization now emits a `message()` so silent
  case-folding of `x`/`y`/`z` headers is visible.
* `subsetData()` error on too-few-matched cells now reports the requested
  cell types and the count that was actually found.
* `computeBidirCorrelation()` guards against empty filtered matrices
  (returns zero correlation with a warning rather than crashing).
* `show()` for CoPro objects now reports approximate object size in MB
  and truncates the metadata field list when there are many columns.
* `runSkrCCA()` optimization loop now reports progress via `message()`
  that can be silenced with `suppressMessages()`.
* `plotG12Functions()` now always returns a stable `list(plot, data,
  summary)` shape where `plot` is always a list with `combined` and
  `individual` elements (one of which may be `NULL`), regardless of
  `plot_type`. Default fallback palette for >8 cell-type pairs now uses
  the colorblind-friendly viridis palette rather than `rainbow()`.
* Replaced several `cat()` progress prints in `computeDistance()` and
  `plotG12Functions()` with `message()` gated on `verbose`, so users can
  suppress output cleanly.

## Documentation

* Added `@examples`, `@family`, and `@seealso` annotations across the
  main pipeline functions so the pkgdown reference auto-cross-links
  related steps.
* New `CONTRIBUTING.md` covering test runs, vignette render workflow,
  and roxygen regeneration.

## Internal

* Renamed `R/80_get_cs_in_situ.r` to `R/80_get_cs_in_situ.R` for
  cross-platform portability on case-sensitive filesystems.
* Bumped `actions/checkout` to v6 in the lint workflow to match the
  other CI jobs.

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
