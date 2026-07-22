# CoPro 1.1.2

## Inference

* Conditional CC2+ permutation tests now use full projection of the fixed
  observed lower axes on every permuted operator. The weighted oblique form is
  used with `scalePCs = FALSE`, and the PC-variance metric is now propagated
  through every observed and null fit.
* Permutation results record their tested sigma family and aggregation rule.
  `calculate_pvalue()` now compares fixed-sigma nulls only with the observed
  statistic at that sigma, identifies conditional p-values at a data-selected
  fixed sigma, and retains max-over-sigma inference for fair-sigma tests.
  Permutation defaults increase to 999 (Monte-Carlo floor 0.001).
* Added `runSlideLevelInference()` for `CoProMulti`: weights and sigma are
  learned without each held-out biological replicate, held-out normalized
  correlations are combined with equal replicate weight, and uncertainty is
  summarized by a replicate sign-flip test and replicate bootstrap interval.
  Cell-level permutation functions now reject `CoProMulti` objects rather than
  presenting cell shuffles as replicate-level inference.
* Permutation provenance is bound to the null it describes. Running
  `runSkrCCAPermu_Conditional()` after a base `runSkrCCAPermu()` on the same
  object no longer re-labels the base-path null, so a later `calculate_pvalue()`
  returns the same p-value and sigma-selection warning regardless of what else
  has been run on the object.

## Performance

* One-cell-type skrCCA now solves the symmetric quadratic problem directly
  with an exact symmetric eigendecomposition, selecting the largest algebraic
  eigenvalues and obtaining all requested axes from one factorization.
* Sparse within-cell-type kernels are stored as symmetric `dsCMatrix` objects,
  so only one triangle is retained. Cross-type and asymmetrically normalized
  kernels remain general `dgCMatrix` objects.
* Sparse expression PCA now passes centering and scaling vectors directly to
  `irlba` instead of materializing a dense centered matrix. All cell types use
  one common feasible PCA rank, and multiset optimizers are also dimension-aware.
* Regression gene scores use the identity `X' (s - mean(s))` and no longer
  construct a centered dense expression matrix for every sigma and axis.
* Gene-space CCA now applies self- and cross-covariances as matrix-free
  operators (`Z_i' K_ij (Z_j w)`) instead of storing dense `G x G` matrices.
  Euclidean streaming builds exact sparse fixed-radius kernels without dense
  distance or kernel matrices.
* Kernel normalizers are cached on the CoPro object with matrix signatures for
  safe reuse. Fair-sigma and conditional permutations now honor `n_cores` via
  memory-explicit PSOCK workers. When CoPro is not installed (for example under
  `devtools::load_all()`), `n_cores > 1` now falls back to sequential execution
  with a warning instead of aborting inside the worker.

## Documentation

* New vignette *Handling very large datasets (Xenium, large MERFISH)*, a how-to
  for keeping large runs in memory: sparse/BPCells input, early `subsetData()`,
  small `nPCA` for targeted panels, and (most importantly) skipping
  `computeDistance()` in favor of the sparse kernel path
  (`computeKernelMatrix(method = "auto"/"sparse")` / `computeSparseKernel()`).
* Lowered `nPCA` from 40/30 to 15 in the targeted-panel vignettes (brain
  MERFISH, colon D3, colon D9, organoid) so they follow the documented 10–15
  guidance for imaging panels. Cell scores and normalized correlations are
  unaffected by `nPCA`; the change improves gene-weight reproducibility and
  score transfer, and makes the vignettes consistent with the kidney example.

## Internal

* Removed the fully-commented legacy optimizer file and three unused
  per-iteration update-vector helpers superseded by the precomputed PC-space
  operator (`Y_resi`) path. Also removed a stranded file of unexported,
  never-called soft-deprecation wrappers (`newCoPro`, `newCoProm`,
  `subsetDataOne`, `subsetDataMulti`, `computePCAMulti`). User-facing error
  messages that referenced the old `computePCAMulti` / `subsetDataMulti` names
  now point to `computePCA` / `subsetData`.
* Consolidated the `computePCA` documentation topic, which was previously split
  under a stale `computePCAMulti` help page with a bogus method alias.

# CoPro 1.1.1

## Performance

* Two-cell-type skrCCA now solves all requested canonical axes with one exact
  SVD of the summed PC-space cross-operator. Multi-slide analyses retain the
  existing memory-efficient summation and do not construct a stacked cell
  matrix or block-diagonal spatial kernel.

## Bug fixes

* Higher skrCCA axes for three or more cell types now use full two-sided
  projection deflation, preventing later axes from reusing directions selected
  by earlier axes under the default scaled-PC formulation.

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
  per-slide cell-type block reaches `autoThreshold` (default 5000) cells or the
  aggregate dense workload reaches `autoThreshold^2` entries. The new
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

## Performance

* Fair-sigma and conditional permutation tests now cache kernel normalizers
  once per bandwidth and reuse each precomputed PC-space operator for fitting
  and scoring. Sparse whitened-Frobenius normalization also stays sparse via an
  equivalent low-rank centering formula instead of materializing dense kernels.
* The exact fixed-radius neighbor search used by sparse kernels now runs in a
  deterministic Rcpp engine, with the original R implementation retained as a
  reference fallback. Bin-wise permutations precompute bin memberships and
  neighbor lookups once per cell type, and normalized permutation scoring
  batches all canonical components instead of rebuilding permuted PC matrices
  inside every pair/component loop.
* `computeSelfKernel()` now supports `method = "auto"`, `"dense"`, or
  `"sparse"`. Its default automatically builds exact sparse multitype
  self-kernels directly from coordinates for large workloads or when dense
  self-distance matrices are unavailable.

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
