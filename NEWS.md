# CoPro 1.3.0

## Choosing the canonical criterion

* **New `objective` argument on `runSkrCCA()`**, with `slideWeight` to control
  how slides are combined, and `space` to select the feature space. The default
  is `objective = "sumcov"`, which is exactly what earlier versions computed, so
  **every existing result is unchanged** unless you opt in.

  The reason this is a choice rather than a fix is that the two criteria are the
  same problem more often than it looks. For a single slide they are
  *identical*. With whitened PCs `X_i'X_i = (n_i - 1) I`, so the constraint
  `||w_i|| = 1` **is** the unit-variance constraint; SUMCOR is scale invariant in
  each `w_i`, so it attains its maximum on that constraint set, where every
  denominator is 1 and it reduces to SUMCOV. Single-slide skrCCA — including
  with more than two cell types — has therefore always been Kettenring SUMCOR,
  and `objective = "sumcor"` is routed to the exact SUMCOV solvers there rather
  than iterating to the same answer.

  They differ only across slides, and there SUMCOV factors exactly as

  ```
  f_cov(w) = sum_{i<j} sum_s sqrt(n_i^(s) n_j^(s)) * sigma_i^(s) * sigma_j^(s) * rho_ij^(s)
  ```

  so it already sums per-slide correlations — weighted by cell count *and* by
  per-slide score scale, because the norm constraint pins only the *pooled*
  variance. The scale factor is what lets a slide with inflated variance along
  the canonical direction dominate. `objective = "sumcor"` drops it;
  `slideWeight = "size"` (default) keeps the cell-count factor, and
  `slideWeight = "equal"` gives strict Kettenring SUMCOR, matching
  `runGeneSpaceCCA()`.

  Read back what was actually optimized with the new `getCCAObjective()`.

* **`space = "gene"`** forwards `runSkrCCA()` to `runGeneSpaceCCA()` (pass the
  bandwidth as `sigmaChoice`), and `runGeneSpaceCCA()` gained `objective`, so the
  space and the criterion can now be varied independently in either direction.

* The PC-space SUMCOR iteration is cheap: `Y_ij^(s)` and the per-slide Gram
  matrices are `nPCA x nPCA` and are built once per sigma, so each sweep costs
  `O(S * nPCA^2)` with no kernel products, against
  `O(S * nnz(K) * nGenes)` per sweep in gene space. It warm-starts from the
  deterministic SUMCOV solution — no RNG dependence — evaluates the true
  objective every sweep, backtracks on a decrease, and returns the best iterate,
  so it can never return something worse than its starting point.

* **`minCellsPerSlide`** (default 10, shared with `runGeneSpaceCCA()`) drops
  slides too thin to divide by under `sumcor`. Under `sumcov` such slides are
  only *reported*, not dropped: a thin slide simply contributes little to the
  summed operator, and dropping it would change results computed before this
  rule existed.

## Batch effects: the larger lever is centering, not the objective

* `objective = "sumcor"` removes per-slide *scale* sensitivity but **not** a
  per-slide *mean shift*, and the mean shift is the bigger problem. PC scores are
  centered globally, so a shared technical shift leaves
  `u_i^(s) ~ M_i^(s) 1 + eps` and the numerator picks up
  `M_i M_j 1'K1`, positive whenever both cell types shift the same way — so
  slides reinforce the batch axis instead of cancelling it. Dividing by `sigma`
  does not rescue this: for a nonnegative kernel the leading singular pair is
  close to the Perron vector, so a constant score reaches
  `rho ~ sigma_max(K)` — near the maximum available — on every slide.

  Pair multi-slide `sumcor` with `computePCA(..., center_per_slide = TRUE)`.
  Gene space already gets this from its per-(slide, cell type)
  z-standardization; PC space does not, and this applies to `sumcov` too.

## Gene-space optimizer: sweep choice, and the sign repair

* **New `sweep` argument on `runGeneSpaceCCA()`** and the gene-space optimizers,
  defaulting to `"gauss-seidel"`. The previous behaviour is `sweep = "jacobi"`
  and reproduces prior results exactly.

  The gene-space power iteration updated every block from the *previous* iterate
  (Jacobi), which is not coordinate ascent: it can settle on the negative
  singular pair as a period-2 sign orbit that the sign-tolerant convergence test
  accepts as converged. That is what the post-hoc sign flip existed to repair.
  The flip is valid for two cell types, where flipping one block negates the
  single pairwise term — but **not** for three or more, where flipping one block
  negates only the pairs touching it and can lower the objective. (Flipping
  *every* block is a no-op: each pairwise term takes two sign flips.) The
  PC-space optimizer never had this problem because it already read
  already-updated blocks.

  Under `"gauss-seidel"` every block update forces `w_i' g_i = ||g_i|| >= 0`;
  summing over blocks gives `f >= 0` at any fixed point, so a negative objective
  is unreachable and no sign repair is needed. The flip is retained only on the
  Jacobi path, and now warns when it fires with three or more cell types.

  Note this guarantee is about the *sign*, not about solution quality: under
  `sumcor` the frozen-`sigma` sweep maximizes a surrogate rather than the
  objective, so which local optimum each sweep finds is data-dependent and
  neither dominates the other.

## Fixes

* **A one-slide `CoProMulti` with three or more cell types no longer errors.**
  `optimize_bilinear_multi_slides()` delegated such objects to the single-slide
  routine, which looks its kernels up with `slide = NULL`; a `CoProMulti` keys
  them by slide ID, so the lookup could only fail with "Could not find kernel
  matrix for pair". One and two cell types took an earlier exact-solver shortcut
  and were unaffected.

* Permutation tests (`runSkrCCAPermu()` and the `FairSigma` / `Conditional`
  variants) now refuse weights fitted under `objective = "sumcor"` rather than
  compare an observed statistic from one criterion against a null built from
  another.

* `vignettes/large_datasets.Rmd` no longer breaks `R CMD check`. Its illustrative
  chunks reference objects too large to ship, and `eval = FALSE` covers knitting
  but not the tangle-and-source step, so `purl = FALSE` is now set on each chunk.

# CoPro 1.2.0

## Choosing sigma from the data

* **New `detectSigmaRange()`.** `sigma` is a distance, so no recommended value
  survives a change of units — microns, pixels and Visium coordinates all need
  different numbers for the same biology. What does transfer is how many
  neighbors a cell is coupled to: the kernel row sum
  `m_a(sigma) = sum_b exp(-0.5 (d_ab / sigma)^2)` is the effective number of
  partners cell `a` talks to, it is dimensionless, and it increases
  monotonically in sigma. `detectSigmaRange()` inverts it, reporting the sigmas
  at which the median cell reaches `minNeighbors` (default 5) and
  `maxNeighbors` (default 20), for every cell-type pair and slide, plus a
  recommended grid and a per-block diagnostic table.

  Cost is bounded by sampling: anchor cells are compared against their nearest
  partners via fixed-radius search rather than against all pairs, so it runs in
  well under a second on hundreds of thousands of cells. Anchors are chosen by
  a deterministic stride, so results are reproducible without a seed.

  The recommended workflow is now:

  ```r
  rng <- detectSigmaRange(obj)
  obj <- computeKernelMatrix(obj, sigmaValues = rng$sigmaValues)
  ```

  `method = "auto"` no longer selects the dense path when an object carries no
  distance matrices, so this works at any size without a `computeDistance()`
  call.

## Breaking changes

* **`normalizeDistance` now defaults to `FALSE`** in `computeDistance()`,
  `computeSelfDistance()`, `computeKernelMatrix()`, `computeSparseKernel()` and
  `computeSparseKernelFloat32()`. Rescaling distances existed so that one
  recommended sigma could travel between datasets; `detectSigmaRange()` now
  does that job without moving the coordinates results are reported on. The
  first step that falls back to the default announces the change once per
  session. **Pass `normalizeDistance = TRUE` explicitly to reproduce results
  from CoPro 1.1.x.**

* **New `normalizeMethod` argument** controls how the reference distance is
  estimated when `normalizeDistance = TRUE`.

  * `"global"` (new default) uses the median nearest-neighbor distance over all
    cells of interest, ignoring their type labels.
  * `"spacing"` uses the median nearest-partner distance measured per
    cell-type block, combined across blocks by median.
  * `"percentile"` reproduces the pre-1.2.0 behavior: the *minimum*, across
    blocks, of a low quantile of pairwise distances.

  The old rule let the single densest block set the unit for the entire object,
  so adding one tightly packed slide or cell type rescaled every other block.
  Both new methods remove that dependence on tissue extent and on which block
  happens to be densest. The percentile is still used to floor small distances
  when `truncateLowDist = TRUE`; the two jobs are now separate quantities
  rather than one shared number.

  `"global"` goes further, and is the default for a reason that only shows up
  when an object holds both cross-type and within-type blocks. Any per-block
  reference is measured on the blocks a step happens to build: a cross-type
  block measures the gap between two compartments, so it moves with
  colocalization, while a within-type block measures one type's own packing, so
  it moves with abundance — in 2-D, a type with a tenth of the cells sits about
  `sqrt(10)` further from its nearest neighbor. `computeDistance()` and
  `computeSelfDistance()` therefore derived different units for the same
  tissue. On a two-type test object where one type is packed into a tenth of
  the other's area, the cross-derived and self-derived factors differ by 1.3x
  under `"spacing"` and by 17x under `"percentile"`; under `"global"` they are
  identical to machine precision, in either order, with no coordination between
  the steps.

  What `"global"` deliberately does not do is equalize effective neighbor
  counts across cell types — a rare type still couples to fewer partners at a
  given bandwidth. That is a question about sigma rather than about units, and
  `detectSigmaRange()` answers it per block. The scale factor fixes the unit;
  sigma fixes the reach.

  `"global"` and `"spacing"` both fall back to `"percentile"` for
  `distType = "Morphology-Aware"`, which has no Euclidean neighbor graph to
  read a spacing from.

* **`normalizeDistance = "inherit"`** is accepted by `computeSelfDistance()` and
  `computeSelfKernel()` as a third value alongside `TRUE` and `FALSE`. It reuses
  the factor `computeDistance()` recorded for the cross-type distances instead
  of deriving a separate one from the within-type blocks, so `sigma = s` means
  the same physical bandwidth for a self-kernel as for a cross-kernel. Deriving
  a separate factor (`TRUE`) now warns when the two disagree, and building
  self-kernels never overwrites the factor already recorded on the object.

## Memory

* **`method = "auto"` now selects the float32 sparse representation** for large
  data instead of float64. It stores 8 bytes per entry against 12, and streams
  one cell-type block at a time rather than caching every block's neighbor
  list, which is the larger saving. On a 12-slide/3-type synthetic benchmark
  the peak dropped from 2.8 GB to 0.75 GB for the same three sigmas, with
  normalized correlations agreeing to 1e-6. `method = "sparse"` still selects
  float64 for exactness checks.

* **`auto` now predicts kernel density before building, and warns when a sparse
  representation would not help.** A fixed-radius kernel is only sparse while
  its support radius stays well below the tissue scale; the radius grows
  linearly in sigma, so density grows as `sigma^d` and saturates. Past about
  two-thirds density, sparse storage costs *more* than dense. The warning
  reports the predicted density and the sigma that would bring it under
  `denseThreshold` (default 0.3), so an out-of-memory run at a too-large sigma
  is now diagnosed rather than merely fatal.

## Bug fixes

* **The sparse kernel path no longer ignores `computeDistance()`.** A CoPro
  object now records the coordinate geometry its distances were built on in a
  new `@distanceGeometry` slot, readable with `getDistanceGeometry()`.
  `computeKernelMatrix(method = "sparse")` used to rebuild coordinates from its
  own `distType` / `xDistScale` / `yDistScale` / `zDistScale` arguments, which
  defaulted to unscaled. Two consequences, both silent: per-axis scales passed
  to `computeDistance()` were dropped, and an object carrying a `z` column got
  3-D kernels even when `computeDistance(distType = "Euclidean2D")` had asked
  for 2-D. Because `method = "auto"` selects the sparse path above
  `autoThreshold` cells, this switched on at exactly the dataset sizes where
  the kernel is not checked by hand. These arguments now default to `NULL`,
  meaning "inherit the recorded geometry"; passing a value that contradicts the
  record is an error rather than a silent override (#39).
* `computeDistance()` on a `CoProMulti` object now records
  `@distanceScaleFactor`, which it previously left empty. Downstream helpers
  fell back to re-deriving the factor from raw coordinates.
* `.sigmaAwareBins()` sizes the bin-permutation grid per axis using the
  recorded coordinate scales. Under anisotropic scaling one axis of the grid
  was previously off by that factor, and `bin` is the default permutation null,
  so this fed into reported p-values.
* `.recoverDistanceScaleFactor()` rebuilds probe distances on the recorded
  geometry instead of assuming raw, unscaled, 2-D `x,y`.
* **Within-type and cross-type distances no longer end up on different
  units.** With `normalizeDistance = TRUE`, `computeDistance()` derived its
  scale factor from the cross-type blocks and `computeSelfDistance()` derived a
  second one from the within-type blocks. The two references differ whenever
  cell types differ in abundance or in how tightly they colocalize — on a
  two-type object with one type packed 10x more densely, the within-type
  reference exceeded the cross-type one in all 30 seeds tried, by 1.07x to
  4.41x — so an object could hold cross and self distances on two units, with
  `@distanceScaleFactor` describing only the cross ones. The new default
  `normalizeMethod = "global"` removes the divergence at its source by
  measuring the cells rather than the blocks (see *Breaking changes* above).
  For the block-based methods, where no shared reference exists, the first
  normalization now wins: a later step adopts the pinned factor instead of
  deriving its own, and says so. `computeDistance()`, which rebuilds
  `@distances` from scratch, and `overwrite = TRUE`, which discards the blocks
  the pin described, are the two ways to re-derive it.
* `computeSelfDistance()` now writes `@distanceScaleFactor`, which it
  previously left untouched while rescaling the blocks it built. On a
  self-distance-only object every downstream helper that maps analysis
  coordinates back to raw ones was reading a factor of 1 for rescaled
  distances.
* `computeSelfDistance()` gained a `normalizeTarget` argument and its geometry
  arguments now default to `NULL`, inheriting the recorded geometry the way the
  kernel entry points do. It previously hardcoded a target of 0.01, so
  `computeDistance(normalizeTarget = 0.05)` followed by `computeSelfDistance()`
  produced blocks 5x apart with nothing reporting the discrepancy. Contradicting
  the record is now an error. **Behavior change:** on a 3-D object with no prior
  `computeDistance()` call, `computeSelfDistance()` now defaults to
  `distType = "Euclidean3D"` (from the coordinate columns present) rather than
  `"Euclidean2D"`; pass `distType` explicitly to pin it.
* **A normalization that could not happen is no longer recorded as one.** When
  `normalizeDistance = TRUE` but no reference could be measured — every
  cell-type block below the 5-cell threshold, or no finite reference among the
  blocks that were built — the step left `normalizeDistance = TRUE` in the
  geometry record beside an untouched scale factor of 1. Nothing had been
  derived from anything, but every later step read that pair as a legitimate
  pinned unit and adopted the 1. The record now reports what happened rather
  than what was asked for, and the step warns. Affected `computeSelfDistance()`
  on both the single- and multi-slide paths, and `computeDistance()` on its two
  multi-slide paths.
* **`computeSparseKernelFloat32()` no longer erases a pinned scale factor.** It
  wrote `@distanceScaleFactor` unconditionally, so calling it with its default
  `normalizeDistance = FALSE` on an object already normalized by
  `computeDistance()` replaced the real factor with 1, and every helper that
  maps analysis coordinates back to raw ones — including the permutation null's
  neighbor graph — silently used the wrong unit. It now guards the write the
  way `computeSparseKernel()` already did.
* **`computeSelfDistance(overwrite = TRUE)` re-derives the scale instead of
  quietly switching normalization off.** `overwrite` clears the geometry record
  to drop the pin, which also dropped the recorded `normalizeDistance = TRUE`,
  so the argument fell back to the new 1.2.0 default of `FALSE` and the result
  was an unnormalized object with a scale factor of 1. The record is now read
  for defaults before it is cleared, so `overwrite` means "re-derive the unit"
  as documented. Arguments supplied alongside `overwrite = TRUE` still win, so
  the geometry remains freely changeable.
* **`normalizeTarget` is validated wherever it is accepted.** Each entry point
  repeated the check inline and `computeSelfDistance()` never gained one, so
  `computeSelfDistance(normalizeTarget = -0.01)` returned a negative scale
  factor and flipped the sign of every distance in the object. The check now
  lives with the other geometry validation and runs on every path.
* **Objects serialized before CoPro 1.2.0 keep their scale factor.** The pin
  was read only when `@distanceGeometry` was populated, so an object carrying a
  valid factor but no record — every object written by an earlier version —
  re-derived a second unit, which is the bug the pin exists to prevent. Such an
  object now keeps its factor. Relatedly, the guard for the missing slot probed
  `methods::slotNames()`, which reads the *class definition* and so answered
  `TRUE` for exactly the objects it was meant to catch; it now uses
  `methods::.hasSlot()`, which asks the instance.

# CoPro 1.1.3

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
* Per-column standard deviations in `center_scale_matrix_opt()` use
  `matrixStats::colSds()` (~2.5x faster than `apply(x, 2, sd)`) via the new
  `.columnSds()` helper, which falls back to `apply()` for sparse input, and
  the nonzero-fraction check reads a sparse matrix's column pointer instead of
  building a full logical copy. `colSds()` and `sd()` use different variance
  algorithms and can disagree by 1 ulp, which is enough to flip the sign of a
  principal component; the CCA weight's coordinate on that PC flips with it and
  the two cancel. **This is the one change in this release that is not
  bit-identical.** Across the 14-scenario baseline it moves the raw weight
  vectors in `@skrCCAOut` (a per-PC sign convention) in 3 scenarios, while every
  reported quantity -- cell scores, gene scores, regression gene weights,
  normalized correlations, null correlations and p-values -- agrees to
  **5.4e-11 relative with no sign changes**, and the selected sigma is identical
  everywhere. Results remain exactly reproducible run to run.
  `test-pca-workflow.R` asserts both that invariance and that this is the
  implementation actually shipped.
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
  in cache: on a 150k-cell kernel, `.whitenedFrobNorm()` itself goes from
  6.16 s / 2231 MB to 4.35 s / 609 MB (1.42x, 3.7x less peak; minimum of 5
  repetitions in each of two checkouts).
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

## Internal

* The permutation entry points share one draw-evaluation path in
  `R/D0_permutation_plan.R` instead of each carrying a sequential worker plus a
  near-identical parallel copy. `runSkrCCAPermu()` no longer maintains its own
  PSOCK cluster block, and three copies of the `"pc"` column-shuffling helper
  and the now-unused `.parallelPermutationLapply()` were removed.
* `irlba::irlba()` is no longer used anywhere in the package; the five
  `@importFrom` directives naming it (three of which were already stale) were
  dropped. `irlba::prcomp_irlba()` in `computePCA()` is unaffected, and remains
  the one place where a seed still changes results.

# CoPro 1.1.2

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
