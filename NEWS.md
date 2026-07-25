# CoPro 1.1.2

## Performance

* The float32 sparse operators pack their dense operand row-major before
  threading. `float32_csr_xky_cpp()` previously read `X` in R's column-major
  layout, so every kernel nonzero touched one cache line per PC; it now reads
  one contiguous run. Measured 3.0-3.5x on the cross-type operators at
  `nPCA = 30`, tapering to ~1.3-1.6x at `nPCA = 10` or on within-type
  (symmetric) kernels. Run-to-run variance also collapses (28.5-54.8 s becomes
  16.6-17.0 s on one within-type case), which matters on a shared node. The
  conversion is the same one the strided loop performed on each access and the
  accumulation order is unchanged, so results are bit-identical.
* Permutation draws for a **held-fixed** cell type are stored as a compact
  marker instead of `replicate(nPermu, 1:n)`. That matrix held the integers
  `1..n` repeated `nPermu` times and carried no information: 799 MB at 200,000
  cells and 999 draws, persisted into the object and into every `saveRDS()`.
  Detecting it no longer allocates a same-sized logical either, and each draw
  now skips the row subset that used to copy the fixed type's whole `n x nPCA`
  matrix back to itself. Every drawn permutation and every p-value is
  unchanged.
* Multi-slide `@pcaResults` stores per-slide *views* of the global PC scores
  (row indices) rather than a second copy of the values, halving what PC scores
  occupy: 72x smaller for the slot itself, and total PC-score memory drops from
  ~128 MB to ~64 MB at 200,000 cells and 40 PCs. `.preparePCMatrices()` also
  whitens the global score matrix once per cell type instead of once per
  (slide, cell type). Slices are still materialized when
  `center_per_slide = TRUE`, where re-centering makes them genuinely different
  data, and objects saved with materialized slices continue to work. The slot's
  shape -- one entry per slide, keyed by slide name -- is unchanged.
* The whitened-Frobenius normalizer no longer materializes `(Rx %*% K) %*% Ry`.
  Those two chained sparse products fill in heavily -- 7x then 11x on a
  40k-cell kernel, extrapolating to roughly 1.5 GB at 200k cells -- to produce
  a single scalar. `<Rx K Ry, K>` is now accumulated as
  `sum((Rx K) * (K Ry))` over column blocks, and the low-rank cross term
  applies the operators to its two-column factor instead of to `K`. Blocking
  turned out to be *faster* as well as smaller, because the intermediates stay
  in cache: on a 150k-cell kernel, 4.09 s / 466 MB against 5.97 s / 1996 MB.
  `sum(K * K)` likewise no longer allocates a second sparse matrix. The
  normalizer is unchanged to ~4e-15 relative (floating-point reassociation);
  weights, cell scores, gene scores and the selected sigma are bit-identical.
* The sparse-kernel upper-quantile clip no longer materializes
  `rep(values, each = 2L)` for symmetric kernels. `.type7QuantileRepeated()`
  reads the two order statistics R's type-7 rule needs directly from the stored
  triangle, by selection rather than a full sort. Measured 2.1x faster and
  600 MB less peak at 30M stored values, scaling to several GB at the
  hundreds-of-millions of nonzeros large panels reach. It mirrors
  `stats::quantile()`'s arithmetic statement for statement, so the clip
  threshold is bit-identical.
* `options(CoPro.compactPermutation = TRUE)` additionally stores the
  *permuted* side as one seed per draw and regenerates it on demand, removing
  the remaining `n * nPermu * 4` bytes for the `"global"` and `"bin"` nulls.
  **Off by default, because it changes which permutations are drawn** -- a
  re-run of a saved analysis will move its p-values within Monte Carlo error.
  Leave it off to reproduce existing results.

## Inference

* skrCCA no longer depends on the RNG state. The block-relaxation path (three
  or more cell types, and every conditional higher axis) initialized its weight
  vectors with `irlba()`. When no starting vector is supplied, `irlba()` draws
  one at random and advances the RNG stream, so the starting direction varied
  between runs, between sessions, and between a sequential run and a PSOCK
  worker. Two of the three initializers were affected on that basis:
  `initialize_next_component()` and `initialize_weights_multi_slide()`, neither
  of which passed a starting vector. (`initialize_weights_svd()` passed a fixed
  one and was already deterministic; it was converted for consistency and
  speed, not to fix a defect.) Because a power iteration converges to whichever
  sign its start points at, this flipped the **sign** of CC2+ weight vectors at
  random and moved converged values at the iteration tolerance (observed: sign
  flips plus 2.6e-5 weight differences across seeds on a three-type run). Gene
  weights are read directionally, so a random sign is not cosmetic. Every
  operator involved is `nPC x nPC` or a Gram matrix of that size, so exact
  LAPACK factorizations are now used instead -- deterministic, and cheaper here
  than a Krylov method. This affects `runSkrCCA()` and
  `optimize_bilinear()`/`optimize_bilinear_n()` with three or more cell types,
  the multi-slide initializer, and conditional higher axes; one- and
  two-cell-type runs already used exact decompositions.
* Two-cell-type conditional higher axes are solved by the same exact SVD that
  produced the observed axis. Deflation leaves an ordinary singular-vector
  problem, so block relaxation was only approaching to `tol` an answer that has
  a closed form, and the null statistic was computed by a different solver than
  the observed statistic it is compared against. On colon D3 the two agree to
  1e-17, so this changes no published number there; it removes the last
  iterative step from the conditional permutation inner loop and makes
  sequential and parallel runs bitwise identical.
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

* Permutation tests no longer rebuild the sparse product `X_i' K_ij X_j` for
  every draw. Under `permu_which = "second_only"` (the default) and
  `"first_only"` one cell type is held fixed across all draws, so the kernel is
  applied to that side once and each draw becomes a small dense product in PC
  space; the per-draw cost falls from `O(nnz(K) * nPCA)` to `O(n * nPCA^2)`.
  The identity is exact, so the null distribution and its p-values are
  unchanged. Pairs with both sides permuted (`permu_which = "both"`, and the
  two-permuted-type pairs of a three-type run) keep the original product. The
  same factorization now serves `runSkrCCAPermu()`,
  `computeNormalizedCorrelationPermu()`, `runSkrCCAPermu_FairSigma()`, and
  `runSkrCCAPermu_Conditional()`. Score norms in the normalized-correlation
  denominator are read from precomputed Gram matrices where the null is a
  genuine bijection on cells; the `"bin"` null draws a spatially matched
  resample and the `"pc"` null shuffles each PC column independently, so both
  keep the direct calculation. Measured against the previous release on a
  5,000-cell pair with ~6M kernel nonzeros and 99 permutations:
  `runSkrCCAPermu()` plus `computeNormalizedCorrelationPermu()` ran 5.0x faster
  (1.23 s to 0.25 s), and `runSkrCCAPermu_Conditional()` over three sigmas and
  two axes 4.5x faster (2.14 s to 0.48 s). The realized factor grows with
  kernel density and draw count, since fixed costs (normalizers, permutation
  generation) do not shrink.
  `options(CoPro.factorizePermutation = FALSE)` restores the unfactorized
  operator for direct comparison. Note that this option isolates the
  factorization only -- it does not restore the previous release's
  normalized-correlation code path, so it is a control for the algebra rather
  than a stand-in for the old timing.
* The factorization is an exact rearrangement of the same triple product in
  real arithmetic, so the null distribution and its p-values are unchanged.
  In floating point the two orderings agree to ~1e-15 relative on double
  kernels. Kernels built by `computeSparseKernelFloat32()` accumulate in single
  precision, where the two paths agree to ~1e-6 relative instead -- far below
  the granularity of a rank-based permutation p-value, but not the ~1e-15 of
  the double path.
* Parallel permutation workers no longer receive the whole kernel list, and no
  longer capture the enclosing `CoPro` object through their closure. Each PSOCK
  worker now gets the precomputed operators, the PC matrices, and only its own
  columns of the permutation index matrix, so peak memory no longer scales with
  `n_cores`. A pair with both sides permuted still needs its kernel, but only
  for the sigma values under test rather than the whole stored list.
  `optimize_bilinear()` and `optimize_bilinear_n()` gained an optional
  `Y_resi` argument for supplying those precomputed operators.
* Added `computeSparseKernelFloat32()` for large single- and multi-slide
  analyses with any number of cell types. It streams one neighbor block at a
  time, retains temporary distances as float32, writes kernels directly into a
  compact float32 CSR representation, stores one-type kernels as a symmetric
  triangle, and uses a row-parallel float32 `X1' K X2` operator without an
  `n_cells x nPCA` intermediate. On the included 200k-cell single-slide
  tiled-colon benchmark, peak RSS fell from 8.37 to 2.83 GiB and the largest
  operator product was 4.57x faster with eight threads while four-component
  cell-score NRMSE remained below 1.2e-6. On an eight-slide 200k-cell
  production-API benchmark, peak RSS fell from 8.86 to 2.34 GiB and kernel
  construction was 2.09x faster. Public kernel accessors temporarily
  materialize standard sparse matrices for legacy plotting and transfer code;
  `materializeFloat32Kernels()` provides a whole-object compatibility escape
  hatch. Global, row, and column kernel normalization now remain in float32,
  including asymmetric normalized self-kernels. The ordinary centered
  Frobenius objective norm is computed directly from encoded value sums;
  exact whitened Frobenius normalization retains a temporary double-sparse
  compatibility fallback. The operator thread count is now determined
  automatically from the cores actually allocated to the process, honoring
  common HPC scheduler variables (`SLURM_CPUS_PER_TASK`, `NSLOTS`,
  `PBS_NUM_PPN`, `LSB_DJOB_NUMPROC`) and `OMP_NUM_THREADS` so a single-core
  allocation no longer oversubscribes a shared node; set
  `options(CoPro.float32Threads=)` to override.
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
* The permutation entry points share one draw-evaluation path in
  `R/D0_permutation_plan.R` instead of each carrying a sequential worker plus a
  near-identical parallel copy. `runSkrCCAPermu()` no longer maintains its own
  PSOCK cluster block, and three copies of the `"pc"` column-shuffling helper
  and the now-unused `.parallelPermutationLapply()` were removed.
* `irlba::irlba()` is no longer used anywhere in the package; the five
  `@importFrom` directives naming it (three of which were already stale) were
  dropped. `irlba::prcomp_irlba()` in `computePCA()` is unaffected, and remains
  the one place where a seed still changes results.

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
