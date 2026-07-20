# CoPro 1.1.0

## Citation

* The CoPro preprint is now available on bioRxiv:
  Miao Z, Qu Y, Huang S, Laux L, Peters S, Aristel A, Zhang Z,
  Niedernhofer L, McMahon A, Kim J, Zhang NR (2026).
  *Dissecting the coordinated progression of cell states in spatial
  transcriptomics with CoPro.* bioRxiv 2026.04.17.719309.
  doi: [10.64898/2026.04.17.719309](https://doi.org/10.64898/2026.04.17.719309).
  `inst/CITATION` and the README have been updated accordingly.

## New features

* Added a sparse, memory-efficient kernel path for large-scale data.
  `computeKernelMatrix()` gains a `method` argument (`"auto"`, `"dense"`,
  `"sparse"`) defaulting to `"auto"`, which selects the sparse path when any
  cell type exceeds `autoThreshold` (default 20000) cells. The new
  `computeSparseKernel()` generic builds sparse `dgCMatrix` Gaussian kernels
  directly from coordinates via an exact fixed-radius neighbor search, never
  forming a dense `n x n` distance or kernel matrix. Results are numerically
  equivalent to the dense path (every pair beyond the kernel's support radius
  is already zero). The sparse path does not require `computeDistance()` to be
  run first and supports `Euclidean2D` / `Euclidean3D` distances.
* `computeKernelMatrix()` gains `dropDistances` (default `TRUE`), which clears
  the large `@distances` slot after kernels are computed, since the downstream
  pipeline only needs the kernels. Set `dropDistances = FALSE` to retain
  distances for inspection via `getDistMat()` or to recompute kernels with new
  sigma values without rebuilding distances.
* Normalized correlation and bandwidth (`sigma`) selection now normalize by the
  whitened-Frobenius norm `||R_x^{1/2} K_c R_y^{1/2}||_F` of the cross-kernel
  instead of its spectral norm `||K||_2`. Here `K_c` is the double-centered
  cross-kernel and `R_x`, `R_y` are the matched-`sigma` within-type kernels;
  this norm is the distribution-free null standard deviation of the bilinear
  statistic `a' K b` and, unlike the spectral norm, does not rail `sigma`
  selection to the grid floor. Affects `computeNormalizedCorrelation()`, the
  permutation tests (`runSkrCCAPermu*()`), and `getTransfer*()` extrapolation.

## Bug fixes

* Sigma-aware bin sizing no longer silently falls back to a hard-coded 10x10
  grid under the default `dropDistances = TRUE`. The raw-to-normalized distance
  scale factor is now stored in a new `@distanceScaleFactor` slot at
  `computeDistance()` time and recovered from there after `@distances` is
  cleared, so `.sigmaAwareBins()` keeps its bandwidth-aware grid.
* The bidirectional-correlation kernel normalizations (`sinkhorn_knopp`,
  `"row_or_col"`) now operate on sparse kernels without densifying them.

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
* `newCoProMulti()` now validates that when `metaData` already contains
  a `slideID` column, its values match the supplied `slideID` argument
  (errors on mismatch instead of silently overwriting).
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
