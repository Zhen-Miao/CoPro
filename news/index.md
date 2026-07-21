# Changelog

## CoPro 1.1.0

### Citation

- The CoPro preprint is now available on bioRxiv: Miao Z, Qu Y, Huang S,
  Laux L, Peters S, Aristel A, Zhang Z, Niedernhofer L, McMahon A, Kim
  J, Zhang NR (2026). *Dissecting the coordinated progression of cell
  states in spatial transcriptomics with CoPro.* bioRxiv
  2026.04.17.719309. doi:
  [10.64898/2026.04.17.719309](https://doi.org/10.64898/2026.04.17.719309).
  `inst/CITATION` and the README have been updated accordingly.

### New features

- Added a sparse, memory-efficient kernel path for large-scale data.
  [`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md)
  gains a `method` argument (`"auto"`, `"dense"`, `"sparse"`) defaulting
  to `"auto"`, which selects the sparse path when any per-slide
  cell-type block reaches `autoThreshold` (default 5000) cells or the
  aggregate dense workload reaches `autoThreshold^2` entries. The new
  [`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md)
  generic builds sparse `dgCMatrix` Gaussian kernels directly from
  coordinates via an exact fixed-radius neighbor search, never forming a
  dense `n x n` distance or kernel matrix. Results are numerically
  equivalent to the dense path (every pair beyond the kernel’s support
  radius is already zero). The sparse path does not require
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  to be run first and supports `Euclidean2D` / `Euclidean3D` distances.
- [`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md)
  gains `dropDistances` (default `TRUE`), which clears the large
  `@distances` slot after kernels are computed, since the downstream
  pipeline only needs the kernels. Set `dropDistances = FALSE` to retain
  distances for inspection via
  [`getDistMat()`](https://zhen-miao.github.io/CoPro/reference/getDistMat.md)
  or to recompute kernels with new sigma values without rebuilding
  distances.
- Normalized correlation and bandwidth (`sigma`) selection now normalize
  by the whitened-Frobenius norm `||R_x^{1/2} K_c R_y^{1/2}||_F` of the
  cross-kernel instead of its spectral norm `||K||_2`. Here `K_c` is the
  double-centered cross-kernel and `R_x`, `R_y` are the matched-`sigma`
  within-type kernels; this norm is the distribution-free null standard
  deviation of the bilinear statistic `a' K b` and, unlike the spectral
  norm, does not rail `sigma` selection to the grid floor. Affects
  [`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md),
  the permutation tests (`runSkrCCAPermu*()`), and `getTransfer*()`
  extrapolation.

### Performance

- Fair-sigma and conditional permutation tests now cache kernel
  normalizers once per bandwidth and reuse each precomputed PC-space
  operator for fitting and scoring. Sparse whitened-Frobenius
  normalization also stays sparse via an equivalent low-rank centering
  formula instead of materializing dense kernels.
- The exact fixed-radius neighbor search used by sparse kernels now runs
  in a deterministic Rcpp engine, with the original R implementation
  retained as a reference fallback. Bin-wise permutations precompute bin
  memberships and neighbor lookups once per cell type, and normalized
  permutation scoring batches all canonical components instead of
  rebuilding permuted PC matrices inside every pair/component loop.
- [`computeSelfKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSelfKernel.md)
  now supports `method = "auto"`, `"dense"`, or `"sparse"`. Its default
  automatically builds exact sparse multitype self-kernels directly from
  coordinates for large workloads or when dense self-distance matrices
  are unavailable.

### Bug fixes

- Sigma-aware bin sizing no longer silently falls back to a hard-coded
  10x10 grid under the default `dropDistances = TRUE`. The
  raw-to-normalized distance scale factor is now stored in a new
  `@distanceScaleFactor` slot at
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  time and recovered from there after `@distances` is cleared, so
  [`.sigmaAwareBins()`](https://zhen-miao.github.io/CoPro/reference/dot-sigmaAwareBins.md)
  keeps its bandwidth-aware grid.

- The bidirectional-correlation kernel normalizations (`sinkhorn_knopp`,
  `"row_or_col"`) now operate on sparse kernels without densifying them.

- Added
  [`asCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/asCoPro.md)
  and
  [`asCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/asCoPro.md)
  S4 generics for one-call coercion from `SingleCellExperiment` and
  `Seurat` objects into CoPro objects. Conversions are gated on their
  respective packages being installed and delegate to the existing
  [`newCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/newCoProSingle.md)
  /
  [`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md)
  constructors so validation stays single-sourced.

- Exposed `normalizeTarget` argument on
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  for users who want to control the target value that the low-percentile
  cell-cell distance is rescaled to. Default preserves existing
  behavior.

### User experience

- [`newCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/newCoProSingle.md)
  /
  [`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md)
  now reject `NA`, `NaN`, and `Inf` values in `normalizedData` at
  construction time with an informative error instead of producing
  cryptic downstream failures.
- [`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md)
  now validates that when `metaData` already contains a `slideID`
  column, its values match the supplied `slideID` argument (errors on
  mismatch instead of silently overwriting).
- `locationData` column standardization now emits a
  [`message()`](https://rdrr.io/r/base/message.html) so silent
  case-folding of `x`/`y`/`z` headers is visible.
- [`subsetData()`](https://zhen-miao.github.io/CoPro/reference/subsetData.md)
  error on too-few-matched cells now reports the requested cell types
  and the count that was actually found.
- [`computeBidirCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeBidirCorrelation.md)
  guards against empty filtered matrices (returns zero correlation with
  a warning rather than crashing).
- `show()` for CoPro objects now reports approximate object size in MB
  and truncates the metadata field list when there are many columns.
- [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
  optimization loop now reports progress via
  [`message()`](https://rdrr.io/r/base/message.html) that can be
  silenced with
  [`suppressMessages()`](https://rdrr.io/r/base/message.html).
- [`plotG12Functions()`](https://zhen-miao.github.io/CoPro/reference/plotG12Functions.md)
  now always returns a stable `list(plot, data, summary)` shape where
  `plot` is always a list with `combined` and `individual` elements (one
  of which may be `NULL`), regardless of `plot_type`. Default fallback
  palette for \>8 cell-type pairs now uses the colorblind-friendly
  viridis palette rather than
  [`rainbow()`](https://rdrr.io/r/grDevices/palettes.html).
- Replaced several [`cat()`](https://rdrr.io/r/base/cat.html) progress
  prints in
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  and
  [`plotG12Functions()`](https://zhen-miao.github.io/CoPro/reference/plotG12Functions.md)
  with [`message()`](https://rdrr.io/r/base/message.html) gated on
  `verbose`, so users can suppress output cleanly.

### Documentation

- Added `@examples`, `@family`, and `@seealso` annotations across the
  main pipeline functions so the pkgdown reference auto-cross-links
  related steps.
- New `CONTRIBUTING.md` covering test runs, vignette render workflow,
  and roxygen regeneration.

### Internal

- Renamed `R/80_get_cs_in_situ.r` to `R/80_get_cs_in_situ.R` for
  cross-platform portability on case-sensitive filesystems.
- Bumped `actions/checkout` to v6 in the lint workflow to match the
  other CI jobs.

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
