# Changelog

## CoPro (development version)

### Multi-sample workflow and numerical diagnostics

- PC-space SUMCOR is the primary documented multi-sample workflow,
  matching the existing `CoProMulti` defaults. Gene-space SUMCOR remains
  available as a frozen-denominator surrogate; its fixed points are not
  guaranteed to be stationary points of the ratio objective.
- [`getCCADiagnostics()`](https://zhen-miao.github.io/CoPro/reference/getCCADiagnostics.md)
  exposes saved PC-space SUMCOR stopping status, projected-gradient
  residuals, objective traces, per-slide Gram conditioning, and
  score-denominator flooring. Supplied axes and equivalent SUMCOV
  reductions are labeled separately. Older objects remain readable
  without diagnostic records. Fitting defaults, objectives, and
  inference are unchanged.
- Documentation correction to two 1.3.0 statements below. The per-slide
  SUMCOR term divides a kernel-smoothed cross term by unsmoothed score
  norms, so it is not bounded by 1 and `slideWeight = "equal"` is not
  strict Kettenring SUMCOR; it gives equal nominal coefficients to the
  slide/pair terms. The term is also not cell-count invariant: uniformly
  replicating cells scales it by `sqrt(r_i r_j)` with an unnormalized
  kernel and by `sqrt(r_i / r_j)` under `normalizeKernel = TRUE`.
  `"size"` adds an explicit `sqrt(n_i n_j)` factor rather than
  “reintroducing” one. No code changed;
  [`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md)
  help now states what `normalizeKernel` does.

### Frozen cross-slide score transfer

- [`fit_score_reference()`](https://zhen-miao.github.io/CoPro/reference/fit_score_reference.md)
  freezes cell-type-specific training means, standard deviations, and
  exact PCA back-projected CoPro weights in a self-contained reference.
  Its [`predict()`](https://rdrr.io/r/stats/predict.html) method scores
  new slides without target quantile normalization or target-derived
  parameters, so scores are invariant to which other target slides are
  supplied. This frozen log-expression workflow is the recommended
  default after an internal benchmark; the existing quantile-normalized
  transfer helpers remain available for cross-platform use. Frozen
  references preserve the PCA low-variance/prevalence scale guard,
  reject gene-space weights, and score sparse targets without densifying
  cells-by-genes chunks.
- The benchmark behind that recommendation was a leave-one-sample-out
  comparison of frozen against target-adaptive scoring maps over four
  spatial datasets (kidney seqFISH, lung Xenium, colon MERFISH, and
  liver Xenium); it did not justify replacing the frozen route. Note
  that the packaged guard is CoPro’s PCA preprocessing rule – low
  variance *or* low prevalence – which is stricter than a bare floor
  against division by zero, and is what makes frozen self-transfer
  reproduce
  [`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md)
  exactly under pooled preprocessing.
- A frozen reference collapses to one center and scale per cell type.
  Under the `CoProMulti` default `computePCA(center_per_slide = TRUE)`,
  which standardizes each (slide, cell type) block separately,
  transferred scores are therefore target-invariant but are not on the
  same affine footing as
  [`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md)
  on the fitted object.
  [`fit_score_reference()`](https://zhen-miao.github.io/CoPro/reference/fit_score_reference.md)
  now says so with a message and records `preprocessing` in the returned
  reference.

### Package-review hardening

- `computePCA(center = FALSE, scale. = TRUE)` no longer divides by an
  unguarded root-mean-square on the dense path. A gene that is never
  detected has divisor zero there, which turned its whole column into
  `NaN` and propagated into the PCA; a near-constant gene had its noise
  amplified. That branch now applies the same degeneracy guard the
  centered dense, sparse, and within-slide paths already used – divisor
  pinned to 1 when a gene’s scale is below `1e-3` or its nonzero
  fraction below 1%. Genes that are not degenerate keep exactly the
  divisor they had. The rule itself now lives in one internal predicate,
  `.unsafeScaleColumns()`, shared by every route into PCA and by the
  frozen score reference, and a missing scale or nonzero fraction counts
  as degenerate instead of slipping through an `NA` in the mask.
- [`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md)
  no longer breaks on out-of-core `BPCells` input outside the
  `center = scale. = TRUE` branch. `.columnNonzeroFraction()` used
  `colSums(x != 0)`, and `IterableMatrix` defines no comparison
  operator, so the *default* multi-slide path
  (`center_per_slide = TRUE`) errored with “comparison (!=) is possible
  only for atomic and list types”; it now streams `binarize(x^2)`, which
  counts negative entries as nonzero and explicitly stored zeros as
  zero, matching dense and sparse matrices.
  `center = TRUE, scale. = FALSE` errored in
  [`sweep()`](https://rdrr.io/r/base/sweep.html) and now broadcasts
  through `add_cols()`, and `center = FALSE, scale. = TRUE` silently
  materialized the matrix densely and now uses `multiply_cols()`.
  BPCells input is also checked for non-finite values during object
  construction instead of letting them propagate into PCA.
- BPCells input no longer needs the Bioconductor package
  `MatrixGenerics` to be installed. `BPCells::colVars()` is an S3
  generic whose `IterableMatrix` method is exported but not registered,
  so called from CoPro it dispatched only when MatrixGenerics was
  present (BPCells then re-registers the method as S4) or BPCells was
  attached; in any other library every BPCells route through
  [`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md)
  and
  [`fit_score_reference()`](https://zhen-miao.github.io/CoPro/reference/fit_score_reference.md)
  stopped with “no applicable method for ‘colVars’”. Column variances
  now come from `BPCells::matrix_stats()`, the same streamed computation
  the method itself performs, so results are unchanged. Caught by the
  BPCells CI job on its first run.
- A dedicated `BPCells-tests` GitHub Actions job installs BPCells, with
  the HDF5 system library it links against, and runs the full test suite
  on pushes to `main`, weekly, and on demand, so the `IterableMatrix`
  branches are exercised in CI rather than only on machines that happen
  to have BPCells. The lint workflow, which had never run, is gone.
- Sparse and float32 kernel construction now honors
  `normalizeDistance = "inherit"` on cross-type paths and rejects
  non-finite coordinates with cell IDs in the error. The R
  nearest-neighbor fallback also uses overflow-safe, vectorized
  character grid keys for very large coordinate ranges.
- Self-kernel calls prune only bandwidths they prove invalid. Requesting
  a narrow self-kernel grid no longer hides otherwise valid
  cross-kernels from a broader existing scan.
  `computeSelfKernel(method = "auto")` now uses float32 sparse storage
  when it selects a fused sparse route, reducing memory at the cost of
  float32 rather than float64 kernel precision.
- Permutation entry points reject `nPermu < 2`, restore the caller’s RNG
  state, account explicitly for failed draws, and use the effective null
  size in Monte-Carlo p-values.
  `calculate_pvalue(alternative = "two.sided")` now errors because the
  null for its optimized, max-aggregated normalized-correlation
  statistic is not symmetric about zero. Signed studentized and
  slide-level statistics retain their two-sided alternatives.
- Optimizers now reject non-finite or invalid `max_iter`, `tol`,
  `step_size`, and CCA metric values consistently. Gene-space
  initialization is deterministic without advancing the caller’s random
  stream, and reducible SUMCOR shortcuts honor the caller’s convergence
  controls.
- Optimizer convergence notices are now conditions emitted by
  [`message()`](https://rdrr.io/r/base/message.html) rather than
  ordinary output from [`print()`](https://rdrr.io/r/base/print.html),
  so [`suppressMessages()`](https://rdrr.io/r/base/message.html)
  silences them. Code that captures those notices should use
  `capture.output(type = "message")`.

### Package options are now function arguments

Four [`options()`](https://rdrr.io/r/base/options.html) flags controlled
behavior that belongs in a call, not in global state. Each is now an
argument. The option is still read to supply that argument’s default, so
scripts that set the option keep their behavior and no result changes.

- `CoPro.factorizePermutation` -\> `factorize` on
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md),
  [`computeNormalizedCorrelationPermu()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelationPermu.md),
  [`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md),
  and
  [`runSkrCCAPermu_Conditional()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_Conditional.md).
- `CoPro.compactPermutation` -\> `compactPermutation` on the three
  `runSkrCCAPermu*()` entry points.
- `CoPro.float32Threads` -\> `nThreads` on
  [`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md).
  `NULL` means “not specified” and is resolved to the option inside
  `.float32KernelThreads()` rather than in a default argument, so that
  `computeKernelMatrix(method = "float32")` — which threads `NULL`
  through explicitly — still honors it.
- `CoPro.useRcppFRNN` -\> `useCompiled` on the internal `.frnnGrid()`.

The gain is that the equivalence tests no longer flip global state
mid-loop to compare two paths; they pass the argument.
`test-sparse-frnn.R` used to set the option, call, set it back, and call
again inside a nested loop.

### One method where there was a pair

Six generics carried a `CoProSingle` method and a `CoProMulti` method
whose bodies were identical, or identical apart from an
`is_multi = TRUE/FALSE` literal the dispatcher can derive from the
object’s own class. Each is now a single method on the virtual `CoPro`
base:
[`getNormCorr()`](https://zhen-miao.github.io/CoPro/reference/getNormalizedCorrelation.md),
[`getCellScoresInSitu()`](https://zhen-miao.github.io/CoPro/reference/getCellScoresInSitu.md),
[`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md),
and
[`computeSelfKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSelfKernel.md).
Dispatch and results are unchanged — `CoProSingle` and `CoProMulti` both
extend `CoPro` — and the generics whose two bodies genuinely differ
([`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md),
[`subsetData()`](https://zhen-miao.github.io/CoPro/reference/subsetData.md),
[`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
and the rest) keep their separate methods.

### Documentation and repository

- [`?CoPro`](https://zhen-miao.github.io/CoPro/reference/CoPro-package.md)
  now states the API naming convention outright: `camelCase` is the
  object pipeline (takes a `CoPro` object, returns one), `snake_case` is
  the engine and utility layer (plain matrices and data frames, no
  object), with the two exports that sit outside the rule named
  explicitly. Also summarized in the README.
- The README Quick Start had drifted to the pre-1.2.0 pipeline. It now
  uses
  [`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md)
  instead of a hardcoded sigma grid, drops the
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  step the sparse routes no longer need, and adds
  [`computeRegressionGeneScores()`](https://zhen-miao.github.io/CoPro/reference/computeRegressionGeneScores.md).
- `R/04_optimization_function_refactored.R` is now
  `R/04_optimization.R`. The `_refactored` suffix on a shipped file read
  as an unfinished refactor.
- Untracked `R/.DS_Store`, `vignettes/.DS_Store`, and the stray
  `..Rcheck/` check output that had been committed. Note that
  `vignettes/*.Rmd.orig` are *not* backups — they are the sources
  `vignettes/precompute.R` knits into the shipped `.Rmd` files, and they
  stay.
- Fixed a detached roxygen block in `R/D0_permutation_plan.R`: the docs
  for `.buildYPlan()` sat above `.permutationPairTypes()`, leaving the
  first undocumented and the second with two glued-together blocks.

## CoPro 1.3.0

### Choosing the canonical criterion

- **New recommended defaults for multi-slide CoPro.**
  [`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md)
  now centers and scales genes within each `(slide, cell type)` block
  and then fits one joint PCA per cell type, so all slides share one
  loading matrix while slide-level location and scale cannot choose the
  retained PC subspace.
  [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
  then defaults to `objective = "sumcor"` with `slideWeight = "equal"`
  for `CoProMulti`. Single-slide objects continue to default to
  `"sumcov"`. Use `center_per_slide = FALSE` and `objective = "sumcov"`
  to reproduce the legacy multi-slide workflow.

- **New `objective` argument on
  [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)**,
  with `slideWeight` to control how slides are combined, and `space` to
  select the feature space.

  The reason this is a choice rather than a fix is that the two criteria
  are the same problem more often than it looks. For a single slide they
  usually coincide. Whitened PCs give `X_i'X_i = (n_i - 1) I`, so on
  `||w_i|| = 1` the denominators are `sigma_i = sqrt(n_i - 1)` — `sigma`
  here is the norm `||X_i w_i||`, not a root-mean-square — and the
  objective is SUMCOV reweighted by the *per-pair* constant
  `m_ij / sqrt((n_i - 1)(n_j - 1))`. A per-pair constant leaves the
  maximizer alone only when every pair gets the same one, so the
  reduction to SUMCOV is exact for **one or two cell types at any cell
  counts**, and for **three or more only when the cell counts are
  equal**. In those cases the SUMCOR implementation uses the direct
  SUMCOV decomposition as an exact computational shortcut.

  With three or more cell types at unequal counts the criteria genuinely
  differ even on one slide. An explicit SUMCOR request now optimizes
  SUMCOR in that case. The implementation checks the one-slide Gram
  matrices directly, so filtering a joint multi-slide PCA down to one
  slide cannot trigger an invalid shortcut.

  Across slides they always differ, and there SUMCOV factors exactly as

      f_cov(w) = sum_{i<j} sum_s sigma_i^(s) * sigma_j^(s) * rho_ij^(s)

  with no `sqrt(n_i n_j)` factor: because `sigma` is a norm,
  `sigma_i sigma_j rho_ij = w_i' Y_ij w_j` is the SUMCOV term already,
  and `rho` is cell-count invariant on its own. So SUMCOV already sums
  per-slide correlations weighted by per-slide score scale, and that
  scale factor is what lets a slide with inflated variance along the
  canonical direction dominate. `objective = "sumcor"` drops it.
  `slideWeight = "equal"` (the multi-slide default) gives strict
  Kettenring SUMCOR and matches
  [`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md);
  `slideWeight = "size"` optionally *reintroduces* the cell-count factor
  `sqrt(n_i n_j)` on its own.

  Read back what was actually optimized with the new
  [`getCCAObjective()`](https://zhen-miao.github.io/CoPro/reference/getCCAObjective.md).

- **`space = "gene"`** forwards
  [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
  to
  [`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md)
  (pass the bandwidth as `sigmaChoice`), and
  [`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md)
  gained `objective`, so the space and the criterion can now be varied
  independently in either direction.

- The PC-space SUMCOR optimizer is cheap: `Y_ij^(s)` and the per-slide
  Gram matrices are `nPCA x nPCA` and are built once per sigma, so each
  iteration costs `O(S * nPCA^2)` with no kernel products. It now
  differentiates the denominator exactly, including the
  `-rho * G_i w_i / sigma_i^2` term that the former frozen-scale
  heuristic omitted. Projected-gradient ascent on the product of spheres
  uses a monotone Armijo line search and stops on the constrained
  gradient norm. It remains deterministic through its SUMCOV warm start
  and applies the same projection against earlier axes for later
  components.

- **`step_size` is honored under both objectives.** It previously did
  nothing under `sumcor` (it was never forwarded to the optimizer), so a
  value chosen for stability was silently lost whenever the objective
  changed — which the new multi-slide `sumcor` default would have made
  the common case. Under `sumcov` it remains the damped power iteration.
  Under `sumcor` a value below 1 replaces the adaptive step with that
  fixed step, and damps the SUMCOV warm start too.

  The two are the same operation. Both the current iterate and the
  retracted candidate lie on the geodesic through `w` in the search
  direction, so blending them and renormalizing *is* a retraction:

      normalize((1 - a) * w + a * R_w(t g)) = R_w(tau * g)
      tau = a * t / ((1 - a) * sqrt(1 + t^2 * |g|^2) + a)

  Damping is therefore a shorter step along the same arc. Applying it to
  the trial step rather than after the line search keeps the Armijo test
  on the point actually returned, so damped SUMCOR runs stay monotone;
  expect more iterations, which is the trade being requested.
  `step_size = 1` leaves the adaptive iteration bit-for-bit unchanged.

- **`minCellsPerSlide`** (default 10, shared with
  [`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md))
  drops slides too thin to divide by under `sumcor`. Under `sumcov` such
  slides are only *reported*, not dropped: a thin slide simply
  contributes little to the summed operator, and dropping it would
  change results computed before this rule existed.

### Batch effects: the larger lever is centering, not the objective

- `objective = "sumcor"` removes per-slide *scale* sensitivity but
  **not** a per-slide *mean shift*, and the mean shift is the bigger
  problem. Under the legacy pooled PCA, PC scores are centered globally,
  so a shared technical shift leaves `u_i^(s) ~ M_i^(s) 1 + eps` and the
  numerator picks up `M_i M_j 1'K1`, positive whenever both cell types
  shift the same way — so slides reinforce the batch axis instead of
  cancelling it. Dividing by `sigma` does not rescue this: for a
  nonnegative kernel the leading singular pair is close to the Perron
  vector, so a constant score reaches `rho ~ sigma_max(K)` — near the
  maximum available — on every slide.

  The multi-slide default `computePCA(..., center_per_slide = TRUE)`
  applies that per-(slide, cell type) standardization **before** PCA and
  then stacks the blocks for one shared loading matrix. Consequently
  each slide’s PC block is centered without a second PC-space
  normalization, and the covariance used to select the truncated PC
  space is batch-shift invariant. Gene space already uses the same
  preprocessing principle.

### Gene-space optimizer: sweep choice, and the sign repair

- **New `sweep` argument on
  [`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md)**
  and the gene-space optimizers, defaulting to `"gauss-seidel"`. The
  previous behaviour is `sweep = "jacobi"` and reproduces prior results
  exactly.

  The gene-space power iteration updated every block from the *previous*
  iterate (Jacobi), which is not coordinate ascent: it can settle on the
  negative singular pair as a period-2 sign orbit that the sign-tolerant
  convergence test accepts as converged. That is what the post-hoc sign
  flip existed to repair. The flip is valid for two cell types, where
  flipping one block negates the single pairwise term — but **not** for
  three or more, where flipping one block negates only the pairs
  touching it and can lower the objective. (Flipping *every* block is a
  no-op: each pairwise term takes two sign flips.) The former PC-space
  heuristic avoided that particular orbit by reading already-updated
  blocks; the new full-gradient PC optimizer does not use a block power
  sweep at all.

  Under `"gauss-seidel"` every block update forces
  `w_i' g_i = ||g_i|| >= 0`; summing over blocks gives `f >= 0` at any
  fixed point, so a negative objective is unreachable and no sign repair
  is needed. The flip is retained only on the Jacobi path, and now warns
  when it fires with three or more cell types.

  Note this guarantee is about the *sign*, not about solution quality:
  under `sumcor` the frozen-`sigma` sweep maximizes a surrogate rather
  than the objective, so which local optimum each sweep finds is
  data-dependent and neither dominates the other.

### Choosing sigma from the data

- **New
  [`selectSigmaByPermutation()`](https://zhen-miao.github.io/CoPro/reference/selectSigmaByPermutation.md).**
  Picks the bandwidth by comparing the co-progression statistic at each
  candidate to *its own* permutation null, instead of taking the largest
  normalized correlation.

  The motivation is that no denominator makes the normalized
  correlation’s null level constant in sigma, so `obj@sigmaValueChoice`
  compares numbers with different floors. The un-whitened `||K_c||_F`
  that ships today ignores within-type spatial autocorrelation; whitened
  variants need a within-type correlation operator that the data do not
  pin down, because the principal components of one cell type do not
  share a single correlation length. On a planted-signal simulation the
  null spread of the raw statistic moved 44x across a 32x sigma grid,
  and even the normalized ratio’s null moved 2.7x and peaked mid-grid —
  so its argmax landed on the wrong bandwidth in 20 of 20 replicates,
  always biased high, while the studentized scan found the
  most-detectable scale in 12 of 20 and the rest one grid step away.

  Rather than model that floor, this measures it. At each sigma the
  observed `T(sigma) = a' K(sigma) b` is compared to the mean and
  standard deviation of its own toroidal-shift null, and the bandwidth
  maximizing `z = (T - mean_null) / sd_null` is selected. Because both
  the location and the scale are read off the null itself, `z` has the
  same null level at every bandwidth by construction and carries no
  tuning constant.

  Centering matters most within a single cell type, where the `w_ii = 0`
  convention takes `sum_i a_i^2 k(x_i, x_i + delta)` off every draw and
  leaves the null a negative mean that grows with sigma — on the
  package’s toy object, from `-0.06 sd_null` at sigma 0.05 to
  `-0.36 sd_null` at sigma 0.2. Across two cell types the null mean is
  zero in expectation and subtracting the estimate costs only
  Monte-Carlo noise.

  One pass of `nPermu` draws evaluates every bandwidth, so the draws are
  coupled across the grid and the same draws give the null of
  `max_sigma z`. A draw is one wrap offset *per cell type*, held across
  every pair that type appears in, so each row of the null is the scan
  statistic of one realizable configuration — with three or more cell
  types a type sits in several pairs at once, and an offset redrawn per
  pair would put the same cells in two places within a single draw. The
  reported p-value is therefore already adjusted for having scanned the
  grid (single-step Westfall–Young max-T) and is not circular — the
  selection is replicated inside the null. Looping
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
  over sigma does *not* give this: it redraws permutations on every call
  and its default bin grid is itself sigma-dependent, so the
  per-bandwidth nulls are neither coupled nor comparable.

  Two guards, both on by default. `minSigma = "spacing"` drops
  candidates below the median nearest-partner distance, where the kernel
  is nearly diagonal and the fixed-direction null understates its floor
  — without it the argmax rails at whatever the smallest candidate
  happens to be. And a selection landing at either end of the grid warns
  regardless of `verbose`, because that is a scan running out of grid
  rather than finding an optimum.

  It covers the **within-type (single cell type) case**, which is where
  the argmax-of-ratio rule is least trustworthy: the stored self-kernel
  has a zero diagonal and so is not a valid whitening operator. (The
  re-optimizing route reaches that case too, but only as of the
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
  fix below — it was unusable before.)

  Directions are held fixed inside each draw rather than re-optimized.
  That is what keeps the *selection* level flat in sigma and what makes
  one `O(B)` pass enough, but it inherits the mild anti-conservativeness
  of a fixed-direction null at small sigma, which is what `minSigma`
  guards. For a re-optimizing test at a chosen bandwidth use
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md),
  or
  [`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md)
  for a re-optimizing max-over-sigma test.

  It requires a **Euclidean** geometry and refuses a
  `"Morphology-Aware"` object rather than rescoring it. The selector
  rebuilds the kernel from coordinates at every candidate bandwidth and
  under every shifted configuration, and a morphology-aware distance is
  not a function of the coordinates alone — its geodesic filter is
  fitted on a k-NN graph of the tissue as observed, which a shift
  invalidates.
  [`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md)
  reuses the stored kernels and so still serves those objects.

  The organoid vignette now selects its bandwidth this way.

### Fixes

- **The SUMCOR permutation null now refits each draw with that draw’s
  own Gram matrices.** SUMCOR divides by `sigma_i = sqrt(w_i' G_i w_i)`,
  and reordering rows leaves `G_i = X_i'X_i` alone — which is why a
  bijective null may build the Grams once and reuse them. The default
  `"bin"` null is not a bijection: it *resamples* cells, so duplicates
  and omissions move `G_i`. `"pc"` shuffles each PC column independently
  and moves it too. Reusing the observed Grams there fitted every draw
  with a numerator from the permuted data and a denominator from the
  observed data, so the null weights did not maximize the criterion the
  observed statistic was measured under, and the resulting p-values were
  not valid. The package already drew this distinction for the score
  norms (`.permutationGrams()`); the SUMCOR optimizer now uses the same
  test and recomputes exactly the Grams a draw invalidates. `"global"`
  and `"toroidal"` are bijections and are unchanged, as is every
  `objective = "sumcov"` run. Reachable since 1.3.0 with three or more
  cell types at unequal counts on one slide under an explicit
  `objective = "sumcor"`.

- **`transferred_weight_1` survives `objective = "sumcor"`.** With one
  slide, where SUMCOR reduces to SUMCOV,
  [`optimize_sumcor_pca_n()`](https://zhen-miao.github.io/CoPro/reference/optimize_sumcor_pca_n.md)
  took an exact all-axis shortcut that ignored the supplied `w_list` and
  recomputed every component. A transferred first axis was therefore
  replaced by the solver’s own, silently, and the later axes were
  conditioned on the wrong direction. The shortcut is now taken only
  when the supplied axes are the ones those solvers would have produced;
  otherwise the run falls through to the sequential path that conditions
  on what was given, matching what the SUMCOV route has always done.

- **`step_size` reaches every axis of the SUMCOV warm start, not just
  the first.** In the one-slide reducible shortcut that warm start *is*
  the returned result, so CC2 and beyond ran undamped however small a
  `step_size` was asked for — the axes a user reaching for damping is
  usually trying to stabilize.

- **Cell-level permutation refuses a gene-space fit explicitly.** The
  resolver reads `@pcaGlobal` and re-optimizes with the PC-space
  solvers, which is the wrong feature space for gene-space weights
  whatever criterion they were fitted under. Not reachable through the
  public API —
  [`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md)
  requires `CoProMulti` and every permutation entry point refuses
  `CoProMulti` first — but the two guards are independent and this one
  no longer leans on the other.

- **[`?runSkrCCA`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
  no longer says permutation rejects every SUMCOR fit.** It now
  describes what 1.3.0 actually does: the null is built with whichever
  criterion the weights were fitted under.

- **Within-type permutation testing worked on paper and not in fact.**
  With a single cell type,
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
  returned normally but
  [`computeNormalizedCorrelationPermu()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelationPermu.md)
  — which has to run next — formed `combn(cts, 2)` and died with a bare
  `n < m`.
  [`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md)
  had the same call.

  The crash was hiding the real defect. `permu_which = "second_only"`
  (the default) permutes every cell type *except the first*, implemented
  as `ct_index > 1`. With one cell type that is never true, so the type
  was never permuted: every draw was the identity, the null distribution
  equalled the observed statistic, and the p-value was 1 by construction
  with nothing in the output to indicate it. Fixing only the `combn`
  call would have turned a loud crash into a silently meaningless
  p-value.

  Both are fixed. `permu_which` has nothing to select between when there
  is one cell type, so it is ignored and that type is permuted, giving
  the within-type null: scores are relabelled against their own
  locations, breaking the association between a cell’s score and where
  it sits while leaving the spatial configuration and the score
  distribution intact. Everything downstream was already built for this
  — `.buildYPlan()` has a `self_pair` mode and the permutation worker
  dispatches to `solve_one_type_eigen()` — so only the entry points
  needed the fix.
  [`runSkrCCAPermu_Conditional()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_Conditional.md)
  previously warned that `"second_only"` gives the identity permutation;
  that warning described the bug and is replaced by a note about the
  behaviour that now happens.

  Multi-type semantics are unchanged: with two or more cell types
  `"second_only"`, `"first_only"` and `"both"` still hold exactly the
  types they always did, and there is a test pinning that.

- **[`getCCAObjective()`](https://zhen-miao.github.io/CoPro/reference/getCCAObjective.md)
  now reports the criterion a gene-space fit actually used.**
  [`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md)
  did not record its own provenance, so the reader’s no-record fallback
  answered `"sumcov"` — the opposite of gene space’s `"sumcor"` default
  — for exactly the call this release tells you to inspect. Because
  weights are merged into `@skrCCAOut` rather than replacing it, an
  earlier
  [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
  record could also survive a later gene-space run and describe the
  wrong fit. Both are fixed, and the record now carries `sweep`.

- **`runSkrCCA(space = "gene", ...)` no longer silently discards
  arguments.** `objective`, `tol`, `maxIter`, `minCellsPerSlide`,
  `scalePCs`, `transferred_weight_1`, `n_cores`, `step_size` and
  `slideWeight` are named formals, so none of them reached `...` and
  none were forwarded — a passed `transferred_weight_1` looked exactly
  like a transfer that ran. The four with a gene-space analogue are now
  forwarded **when supplied** (not when left at a default, since the two
  entry points differ: single-slide PCA defaults to objective `"sumcov"`
  while multi-slide PCA and gene space default to `"sumcor"`; `tol` is
  `1e-5` vs `1e-6`; `maxIter` is `200` vs `3000`; `minCellsPerSlide`
  `10` vs `20`), and the rest are an error.

- **A one-slide `CoProMulti` with three or more cell types no longer
  errors.**
  [`optimize_bilinear_multi_slides()`](https://zhen-miao.github.io/CoPro/reference/optimize_bilinear_multi_slides.md)
  delegated such objects to the single-slide routine, which looks its
  kernels up with `slide = NULL`; a `CoProMulti` keys them by slide ID,
  so the lookup could only fail with “Could not find kernel matrix for
  pair”. One and two cell types took an earlier exact-solver shortcut
  and were unaffected.

- **Behavior change for an explicit single-slide
  `objective = "sumcor"`.** The default for single-slide objects is
  still `"sumcov"`, so a default call is unchanged. But an *explicit*
  one-slide `"sumcor"` request used to be routed to the SUMCOV solvers
  unconditionally (with a warning when that was only an approximation);
  it now runs the full-gradient optimizer whenever the reduction does
  not hold, i.e. three or more cell types at unequal counts. Those calls
  return different weights than in 1.2.x. One and two cell types, and
  equal counts, are unaffected because the criteria coincide there.

- **Permutation tests match the criterion the weights were fitted
  with**, rather than refusing every `objective = "sumcor"` fit.
  Cell-level permutation is single-slide only (`CoProMulti` is directed
  to
  [`runSlideLevelInference()`](https://zhen-miao.github.io/CoPro/reference/runSlideLevelInference.md)
  first, as before), so the SUMCOR reduction test decides:

  - one or two cell types at any counts, or three or more at equal
    counts – the fitted weights *are* the SUMCOV maximizer, so the
    existing SUMCOV null is already the matching null and the test
    proceeds;
  - three or more at unequal counts – the criteria genuinely differ, and
    [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
    re-optimizes every draw under SUMCOR.

  A within-slide label permutation permutes the rows of `X_i`, which
  leaves `G_i = X_i'X_i` unchanged, so the per-slide scales SUMCOR
  divides by are permutation-invariant: they are built once and reused
  across draws, and the existing `Y` operator-reuse factorization is
  untouched.

  [`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md)
  and
  [`runSkrCCAPermu_Conditional()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_Conditional.md)
  are restricted to at most two cell types, where the reduction always
  holds, so they are unaffected; they carry a guard in case that
  restriction is relaxed.

- `vignettes/large_datasets.Rmd` no longer breaks `R CMD check`. Its
  illustrative chunks reference objects too large to ship, and
  `eval = FALSE` covers knitting but not the tangle-and-source step, so
  `purl = FALSE` is now set on each chunk.

### Damped power iteration in gene space

- **[`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md)
  gains `step_size`.** The PC-space optimizer
  ([`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
  / `optimize_bilinear*()`) has offered damped power iteration for some
  time; the gene-space optimizer did not, so a gene-space fit that
  oscillated had no recourse other than raising `max_iter` and watching
  it fail again.
  [`optimize_genespace_avg_corr()`](https://zhen-miao.github.io/CoPro/reference/optimize_genespace_avg_corr.md),
  [`optimize_genespace_avg_corr_n()`](https://zhen-miao.github.io/CoPro/reference/optimize_genespace_avg_corr_n.md)
  and
  [`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md)
  now accept `step_size` in (0, 1\], blending each new iterate with the
  previous one:

      w_new = normalize((1 - step_size) * w_old + step_size * w_update)

  The default is `step_size = 1`, which is exactly the previous pure
  power iteration — existing results do not move. In the deflation path
  the blended iterate is re-projected against the previous components,
  since blending can reintroduce directions that were just deflated out.

  Damping composes with `sweep`: the blend is always taken against the
  previous outer iterate, and under `"gauss-seidel"` the per-slide
  scales are refreshed from the damped weight, so later blocks in the
  sweep divide by the scales that were actually adopted.

  Damping trades iterations for stability: reach for it when a fit
  oscillates rather than converges, not as a routine setting.

## CoPro 1.2.0

### Choosing sigma from the data

- **New
  [`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md).**
  `sigma` is a distance, so no recommended value survives a change of
  units — microns, pixels and Visium coordinates all need different
  numbers for the same biology. What does transfer is how many neighbors
  a cell is coupled to: the kernel row sum
  `m_a(sigma) = sum_b exp(-0.5 (d_ab / sigma)^2)` is the effective
  number of partners cell `a` talks to, it is dimensionless, and it
  increases monotonically in sigma.
  [`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md)
  inverts it, reporting the sigmas at which the median cell reaches
  `minNeighbors` (default 5) and `maxNeighbors` (default 20), for every
  cell-type pair and slide, plus a recommended grid and a per-block
  diagnostic table.

  Cost is bounded by sampling: anchor cells are compared against their
  nearest partners via fixed-radius search rather than against all
  pairs, so it runs in well under a second on hundreds of thousands of
  cells. Anchors are chosen by a deterministic stride, so results are
  reproducible without a seed.

  The recommended workflow is now:

  ``` r

  rng <- detectSigmaRange(obj)
  obj <- computeKernelMatrix(obj, sigmaValues = rng$sigmaValues)
  ```

  `method = "auto"` no longer selects the dense path when an object
  carries no distance matrices, so this works at any size without a
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  call.

### Breaking changes

- **`normalizeDistance` now defaults to `FALSE`** in
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md),
  [`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
  [`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md)
  and
  [`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md).
  Rescaling distances existed so that one recommended sigma could travel
  between datasets;
  [`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md)
  now does that job without moving the coordinates results are reported
  on. The first step that falls back to the default announces the change
  once per session. **Pass `normalizeDistance = TRUE` explicitly to
  reproduce results from CoPro 1.1.x.**

- **New `normalizeMethod` argument** controls how the reference distance
  is estimated when `normalizeDistance = TRUE`.

  - `"global"` (new default) uses the median nearest-neighbor distance
    over all cells of interest, ignoring their type labels.
  - `"spacing"` uses the median nearest-partner distance measured per
    cell-type block, combined across blocks by median.
  - `"percentile"` reproduces the pre-1.2.0 behavior: the *minimum*,
    across blocks, of a low quantile of pairwise distances.

  The old rule let the single densest block set the unit for the entire
  object, so adding one tightly packed slide or cell type rescaled every
  other block. Both new methods remove that dependence on tissue extent
  and on which block happens to be densest. The percentile is still used
  to floor small distances when `truncateLowDist = TRUE`; the two jobs
  are now separate quantities rather than one shared number.

  `"global"` goes further, and is the default for a reason that only
  shows up when an object holds both cross-type and within-type blocks.
  Any per-block reference is measured on the blocks a step happens to
  build: a cross-type block measures the gap between two compartments,
  so it moves with colocalization, while a within-type block measures
  one type’s own packing, so it moves with abundance — in 2-D, a type
  with a tenth of the cells sits about `sqrt(10)` further from its
  nearest neighbor.
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  and
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md)
  therefore derived different units for the same tissue. On a two-type
  test object where one type is packed into a tenth of the other’s area,
  the cross-derived and self-derived factors differ by 1.3x under
  `"spacing"` and by 17x under `"percentile"`; under `"global"` they are
  identical to machine precision, in either order, with no coordination
  between the steps.

  What `"global"` deliberately does not do is equalize effective
  neighbor counts across cell types — a rare type still couples to fewer
  partners at a given bandwidth. That is a question about sigma rather
  than about units, and
  [`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md)
  answers it per block. The scale factor fixes the unit; sigma fixes the
  reach.

  `"global"` and `"spacing"` both fall back to `"percentile"` for
  `distType = "Morphology-Aware"`, which has no Euclidean neighbor graph
  to read a spacing from.

- **`normalizeDistance = "inherit"`** is accepted by
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md)
  and
  [`computeSelfKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSelfKernel.md)
  as a third value alongside `TRUE` and `FALSE`. It reuses the factor
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  recorded for the cross-type distances instead of deriving a separate
  one from the within-type blocks, so `sigma = s` means the same
  physical bandwidth for a self-kernel as for a cross-kernel. Deriving a
  separate factor (`TRUE`) now warns when the two disagree, and building
  self-kernels never overwrites the factor already recorded on the
  object.

### Memory

- **`method = "auto"` now selects the float32 sparse representation**
  for large data instead of float64. It stores 8 bytes per entry against
  12, and streams one cell-type block at a time rather than caching
  every block’s neighbor list, which is the larger saving. On a
  12-slide/3-type synthetic benchmark the peak dropped from 2.8 GB to
  0.75 GB for the same three sigmas, with normalized correlations
  agreeing to 1e-6. `method = "sparse"` still selects float64 for
  exactness checks.

- **`auto` now predicts kernel density before building, and warns when a
  sparse representation would not help.** A fixed-radius kernel is only
  sparse while its support radius stays well below the tissue scale; the
  radius grows linearly in sigma, so density grows as `sigma^d` and
  saturates. Past about two-thirds density, sparse storage costs *more*
  than dense. The warning reports the predicted density and the sigma
  that would bring it under `denseThreshold` (default 0.3), so an
  out-of-memory run at a too-large sigma is now diagnosed rather than
  merely fatal.

### Bug fixes

- **The sparse kernel path no longer ignores
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md).**
  A CoPro object now records the coordinate geometry its distances were
  built on in a new `@distanceGeometry` slot, readable with
  [`getDistanceGeometry()`](https://zhen-miao.github.io/CoPro/reference/getDistanceGeometry.md).
  `computeKernelMatrix(method = "sparse")` used to rebuild coordinates
  from its own `distType` / `xDistScale` / `yDistScale` / `zDistScale`
  arguments, which defaulted to unscaled. Two consequences, both silent:
  per-axis scales passed to
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  were dropped, and an object carrying a `z` column got 3-D kernels even
  when `computeDistance(distType = "Euclidean2D")` had asked for 2-D.
  Because `method = "auto"` selects the sparse path above
  `autoThreshold` cells, this switched on at exactly the dataset sizes
  where the kernel is not checked by hand. These arguments now default
  to `NULL`, meaning “inherit the recorded geometry”; passing a value
  that contradicts the record is an error rather than a silent override
  ([\#39](https://github.com/Zhen-Miao/CoPro/issues/39)).
- [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  on a `CoProMulti` object now records `@distanceScaleFactor`, which it
  previously left empty. Downstream helpers fell back to re-deriving the
  factor from raw coordinates.
- [`.sigmaAwareBins()`](https://zhen-miao.github.io/CoPro/reference/dot-sigmaAwareBins.md)
  sizes the bin-permutation grid per axis using the recorded coordinate
  scales. Under anisotropic scaling one axis of the grid was previously
  off by that factor, and `bin` is the default permutation null, so this
  fed into reported p-values.
- [`.recoverDistanceScaleFactor()`](https://zhen-miao.github.io/CoPro/reference/dot-recoverDistanceScaleFactor.md)
  rebuilds probe distances on the recorded geometry instead of assuming
  raw, unscaled, 2-D `x,y`.
- **Within-type and cross-type distances no longer end up on different
  units.** With `normalizeDistance = TRUE`,
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  derived its scale factor from the cross-type blocks and
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md)
  derived a second one from the within-type blocks. The two references
  differ whenever cell types differ in abundance or in how tightly they
  colocalize — on a two-type object with one type packed 10x more
  densely, the within-type reference exceeded the cross-type one in all
  30 seeds tried, by 1.07x to 4.41x — so an object could hold cross and
  self distances on two units, with `@distanceScaleFactor` describing
  only the cross ones. The new default `normalizeMethod = "global"`
  removes the divergence at its source by measuring the cells rather
  than the blocks (see *Breaking changes* above). For the block-based
  methods, where no shared reference exists, the first normalization now
  wins: a later step adopts the pinned factor instead of deriving its
  own, and says so.
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
  which rebuilds `@distances` from scratch, and `overwrite = TRUE`,
  which discards the blocks the pin described, are the two ways to
  re-derive it.
- [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md)
  now writes `@distanceScaleFactor`, which it previously left untouched
  while rescaling the blocks it built. On a self-distance-only object
  every downstream helper that maps analysis coordinates back to raw
  ones was reading a factor of 1 for rescaled distances.
- [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md)
  gained a `normalizeTarget` argument and its geometry arguments now
  default to `NULL`, inheriting the recorded geometry the way the kernel
  entry points do. It previously hardcoded a target of 0.01, so
  `computeDistance(normalizeTarget = 0.05)` followed by
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md)
  produced blocks 5x apart with nothing reporting the discrepancy.
  Contradicting the record is now an error. **Behavior change:** on a
  3-D object with no prior
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  call,
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md)
  now defaults to `distType = "Euclidean3D"` (from the coordinate
  columns present) rather than `"Euclidean2D"`; pass `distType`
  explicitly to pin it.
- **A normalization that could not happen is no longer recorded as
  one.** When `normalizeDistance = TRUE` but no reference could be
  measured — every cell-type block below the 5-cell threshold, or no
  finite reference among the blocks that were built — the step left
  `normalizeDistance = TRUE` in the geometry record beside an untouched
  scale factor of 1. Nothing had been derived from anything, but every
  later step read that pair as a legitimate pinned unit and adopted
  the 1. The record now reports what happened rather than what was asked
  for, and the step warns. Affected
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md)
  on both the single- and multi-slide paths, and
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  on its two multi-slide paths.
- **[`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md)
  no longer erases a pinned scale factor.** It wrote
  `@distanceScaleFactor` unconditionally, so calling it with its default
  `normalizeDistance = FALSE` on an object already normalized by
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  replaced the real factor with 1, and every helper that maps analysis
  coordinates back to raw ones — including the permutation null’s
  neighbor graph — silently used the wrong unit. It now guards the write
  the way
  [`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md)
  already did.
- **`computeSelfDistance(overwrite = TRUE)` re-derives the scale instead
  of quietly switching normalization off.** `overwrite` clears the
  geometry record to drop the pin, which also dropped the recorded
  `normalizeDistance = TRUE`, so the argument fell back to the new 1.2.0
  default of `FALSE` and the result was an unnormalized object with a
  scale factor of 1. The record is now read for defaults before it is
  cleared, so `overwrite` means “re-derive the unit” as documented.
  Arguments supplied alongside `overwrite = TRUE` still win, so the
  geometry remains freely changeable.
- **`normalizeTarget` is validated wherever it is accepted.** Each entry
  point repeated the check inline and
  [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md)
  never gained one, so `computeSelfDistance(normalizeTarget = -0.01)`
  returned a negative scale factor and flipped the sign of every
  distance in the object. The check now lives with the other geometry
  validation and runs on every path.
- **Objects serialized before CoPro 1.2.0 keep their scale factor.** The
  pin was read only when `@distanceGeometry` was populated, so an object
  carrying a valid factor but no record — every object written by an
  earlier version — re-derived a second unit, which is the bug the pin
  exists to prevent. Such an object now keeps its factor. Relatedly, the
  guard for the missing slot probed
  [`methods::slotNames()`](https://rdrr.io/r/methods/slot.html), which
  reads the *class definition* and so answered `TRUE` for exactly the
  objects it was meant to catch; it now uses
  [`methods::.hasSlot()`](https://rdrr.io/r/methods/slot.html), which
  asks the instance.

## CoPro 1.1.3

### Performance

- The float32 sparse operators pack their dense operand row-major before
  threading. `float32_csr_xky_cpp()` previously read `X` in R’s
  column-major layout, so every kernel nonzero touched one cache line
  per PC; it now reads one contiguous run. Measured 3.0-3.5x on the
  cross-type operators at `nPCA = 30`, tapering to ~1.3-1.6x at
  `nPCA = 10` or on within-type (symmetric) kernels. Run-to-run variance
  also collapses (28.5-54.8 s becomes 16.6-17.0 s on one within-type
  case), which matters on a shared node. The conversion is the same one
  the strided loop performed on each access and the accumulation order
  is unchanged, so results are bit-identical.
- Permutation draws for a **held-fixed** cell type are stored as a
  compact marker instead of `replicate(nPermu, 1:n)`. That matrix held
  the integers `1..n` repeated `nPermu` times and carried no
  information: 799 MB at 200,000 cells and 999 draws, persisted into the
  object and into every
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html). Detecting it no
  longer allocates a same-sized logical either, and each draw now skips
  the row subset that used to copy the fixed type’s whole `n x nPCA`
  matrix back to itself. Every drawn permutation and every p-value is
  unchanged.
- Per-column standard deviations in `center_scale_matrix_opt()` use
  [`matrixStats::colSds()`](https://rdrr.io/pkg/matrixStats/man/rowSds.html)
  (~2.5x faster than `apply(x, 2, sd)`) via the new `.columnSds()`
  helper, which falls back to
  [`apply()`](https://rdrr.io/r/base/apply.html) for sparse input, and
  the nonzero-fraction check reads a sparse matrix’s column pointer
  instead of building a full logical copy. `colSds()` and
  [`sd()`](https://rdrr.io/r/stats/sd.html) use different variance
  algorithms and can disagree by 1 ulp, which is enough to flip the sign
  of a principal component; the CCA weight’s coordinate on that PC flips
  with it and the two cancel. **This is the one change in this release
  that is not bit-identical.** Across the 14-scenario baseline it moves
  the raw weight vectors in `@skrCCAOut` (a per-PC sign convention) in 3
  scenarios, while every reported quantity – cell scores, gene scores,
  regression gene weights, normalized correlations, null correlations
  and p-values – agrees to **5.4e-11 relative with no sign changes**,
  and the selected sigma is identical everywhere. Results remain exactly
  reproducible run to run. `test-pca-workflow.R` asserts both that
  invariance and that this is the implementation actually shipped.
- Multi-slide `@pcaResults` stores per-slide *views* of the global PC
  scores (row indices) rather than a second copy of the values, halving
  what PC scores occupy: 72x smaller for the slot itself, and total
  PC-score memory drops from ~128 MB to ~64 MB at 200,000 cells and 40
  PCs. `.preparePCMatrices()` also whitens the global score matrix once
  per cell type instead of once per (slide, cell type). The recommended
  within-slide centering now happens before the shared PCA, so those
  score blocks also remain views; no post-PCA copy is needed. Objects
  saved with the former materialized, re-centered slices remain
  readable. The slot’s shape – one entry per slide, keyed by slide name
  – is unchanged.
- The whitened-Frobenius normalizer no longer materializes
  `(Rx %*% K) %*% Ry`. Those two chained sparse products fill in heavily
  – 7x then 11x on a 40k-cell kernel, extrapolating to roughly 1.5 GB at
  200k cells – to produce a single scalar. `<Rx K Ry, K>` is now
  accumulated as `sum((Rx K) * (K Ry))` over column blocks, and the
  low-rank cross term applies the operators to its two-column factor
  instead of to `K`. Blocking turned out to be *faster* as well as
  smaller, because the intermediates stay in cache: on a 150k-cell
  kernel, `.whitenedFrobNorm()` itself goes from 6.16 s / 2231 MB to
  4.35 s / 609 MB (1.42x, 3.7x less peak; minimum of 5 repetitions in
  each of two checkouts). `sum(K * K)` likewise no longer allocates a
  second sparse matrix. The normalizer is unchanged to ~4e-15 relative
  (floating-point reassociation); weights, cell scores, gene scores and
  the selected sigma are bit-identical.
- The sparse-kernel upper-quantile clip no longer materializes
  `rep(values, each = 2L)` for symmetric kernels.
  `.type7QuantileRepeated()` reads the two order statistics R’s type-7
  rule needs directly from the stored triangle, by selection rather than
  a full sort. Measured 2.1x faster and 600 MB less peak at 30M stored
  values, scaling to several GB at the hundreds-of-millions of nonzeros
  large panels reach. It mirrors
  [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html)’s
  arithmetic statement for statement, so the clip threshold is
  bit-identical.
- `options(CoPro.compactPermutation = TRUE)` additionally stores the
  *permuted* side as one seed per draw and regenerates it on demand,
  removing the remaining `n * nPermu * 4` bytes for the `"global"` and
  `"bin"` nulls. **Off by default, because it changes which permutations
  are drawn** – a re-run of a saved analysis will move its p-values
  within Monte Carlo error. Leave it off to reproduce existing results.

### Inference

- skrCCA no longer depends on the RNG state. The block-relaxation path
  (three or more cell types, and every conditional higher axis)
  initialized its weight vectors with `irlba()`. When no starting vector
  is supplied, `irlba()` draws one at random and advances the RNG
  stream, so the starting direction varied between runs, between
  sessions, and between a sequential run and a PSOCK worker. Two of the
  three initializers were affected on that basis:
  `initialize_next_component()` and `initialize_weights_multi_slide()`,
  neither of which passed a starting vector. (`initialize_weights_svd()`
  passed a fixed one and was already deterministic; it was converted for
  consistency and speed, not to fix a defect.) Because a power iteration
  converges to whichever sign its start points at, this flipped the
  **sign** of CC2+ weight vectors at random and moved converged values
  at the iteration tolerance (observed: sign flips plus 2.6e-5 weight
  differences across seeds on a three-type run). Gene weights are read
  directionally, so a random sign is not cosmetic. Every operator
  involved is `nPC x nPC` or a Gram matrix of that size, so exact LAPACK
  factorizations are now used instead – deterministic, and cheaper here
  than a Krylov method. This affects
  [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
  and
  [`optimize_bilinear()`](https://zhen-miao.github.io/CoPro/reference/optimize_bilinear.md)/[`optimize_bilinear_n()`](https://zhen-miao.github.io/CoPro/reference/optimize_bilinear_n.md)
  with three or more cell types, the multi-slide initializer, and
  conditional higher axes; one- and two-cell-type runs already used
  exact decompositions.
- Two-cell-type conditional higher axes are solved by the same exact SVD
  that produced the observed axis. Deflation leaves an ordinary
  singular-vector problem, so block relaxation was only approaching to
  `tol` an answer that has a closed form, and the null statistic was
  computed by a different solver than the observed statistic it is
  compared against. On colon D3 the two agree to 1e-17, so this changes
  no published number there; it removes the last iterative step from the
  conditional permutation inner loop and makes sequential and parallel
  runs bitwise identical.
- Conditional CC2+ permutation tests now use full projection of the
  fixed observed lower axes on every permuted operator. The weighted
  oblique form is used with `scalePCs = FALSE`, and the PC-variance
  metric is now propagated through every observed and null fit.
- Permutation results record their tested sigma family and aggregation
  rule.
  [`calculate_pvalue()`](https://zhen-miao.github.io/CoPro/reference/calculate_pvalue.md)
  now compares fixed-sigma nulls only with the observed statistic at
  that sigma, identifies conditional p-values at a data-selected fixed
  sigma, and retains max-over-sigma inference for fair-sigma tests.
  Permutation defaults increase to 999 (Monte-Carlo floor 0.001).
- Added
  [`runSlideLevelInference()`](https://zhen-miao.github.io/CoPro/reference/runSlideLevelInference.md)
  for `CoProMulti`: weights and sigma are learned without each held-out
  biological replicate, held-out normalized correlations are combined
  with equal replicate weight, and uncertainty is summarized by a
  replicate sign-flip test and replicate bootstrap interval. Cell-level
  permutation functions now reject `CoProMulti` objects rather than
  presenting cell shuffles as replicate-level inference.
- Permutation provenance is bound to the null it describes. Running
  [`runSkrCCAPermu_Conditional()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_Conditional.md)
  after a base
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
  on the same object no longer re-labels the base-path null, so a later
  [`calculate_pvalue()`](https://zhen-miao.github.io/CoPro/reference/calculate_pvalue.md)
  returns the same p-value and sigma-selection warning regardless of
  what else has been run on the object.

### Performance

- Permutation tests no longer rebuild the sparse product `X_i' K_ij X_j`
  for every draw. Under `permu_which = "second_only"` (the default) and
  `"first_only"` one cell type is held fixed across all draws, so the
  kernel is applied to that side once and each draw becomes a small
  dense product in PC space; the per-draw cost falls from
  `O(nnz(K) * nPCA)` to `O(n * nPCA^2)`. The identity is exact, so the
  null distribution and its p-values are unchanged. Pairs with both
  sides permuted (`permu_which = "both"`, and the two-permuted-type
  pairs of a three-type run) keep the original product. The same
  factorization now serves
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md),
  [`computeNormalizedCorrelationPermu()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelationPermu.md),
  [`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md),
  and
  [`runSkrCCAPermu_Conditional()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_Conditional.md).
  Score norms in the normalized-correlation denominator are read from
  precomputed Gram matrices where the null is a genuine bijection on
  cells; the `"bin"` null draws a spatially matched resample and the
  `"pc"` null shuffles each PC column independently, so both keep the
  direct calculation. Measured against the previous release on a
  5,000-cell pair with ~6M kernel nonzeros and 99 permutations:
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
  plus
  [`computeNormalizedCorrelationPermu()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelationPermu.md)
  ran 5.0x faster (1.23 s to 0.25 s), and
  [`runSkrCCAPermu_Conditional()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_Conditional.md)
  over three sigmas and two axes 4.5x faster (2.14 s to 0.48 s). The
  realized factor grows with kernel density and draw count, since fixed
  costs (normalizers, permutation generation) do not shrink.
  `options(CoPro.factorizePermutation = FALSE)` restores the
  unfactorized operator for direct comparison. Note that this option
  isolates the factorization only – it does not restore the previous
  release’s normalized-correlation code path, so it is a control for the
  algebra rather than a stand-in for the old timing.
- The factorization is an exact rearrangement of the same triple product
  in real arithmetic, so the null distribution and its p-values are
  unchanged. In floating point the two orderings agree to ~1e-15
  relative on double kernels. Kernels built by
  [`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md)
  accumulate in single precision, where the two paths agree to ~1e-6
  relative instead – far below the granularity of a rank-based
  permutation p-value, but not the ~1e-15 of the double path.
- Parallel permutation workers no longer receive the whole kernel list,
  and no longer capture the enclosing `CoPro` object through their
  closure. Each PSOCK worker now gets the precomputed operators, the PC
  matrices, and only its own columns of the permutation index matrix, so
  peak memory no longer scales with `n_cores`. A pair with both sides
  permuted still needs its kernel, but only for the sigma values under
  test rather than the whole stored list.
  [`optimize_bilinear()`](https://zhen-miao.github.io/CoPro/reference/optimize_bilinear.md)
  and
  [`optimize_bilinear_n()`](https://zhen-miao.github.io/CoPro/reference/optimize_bilinear_n.md)
  gained an optional `Y_resi` argument for supplying those precomputed
  operators.
- Added
  [`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md)
  for large single- and multi-slide analyses with any number of cell
  types. It streams one neighbor block at a time, retains temporary
  distances as float32, writes kernels directly into a compact float32
  CSR representation, stores one-type kernels as a symmetric triangle,
  and uses a row-parallel float32 `X1' K X2` operator without an
  `n_cells x nPCA` intermediate. On the included 200k-cell single-slide
  tiled-colon benchmark, peak RSS fell from 8.37 to 2.83 GiB and the
  largest operator product was 4.57x faster with eight threads while
  four-component cell-score NRMSE remained below 1.2e-6. On an
  eight-slide 200k-cell production-API benchmark, peak RSS fell from
  8.86 to 2.34 GiB and kernel construction was 2.09x faster. Public
  kernel accessors temporarily materialize standard sparse matrices for
  legacy plotting and transfer code;
  [`materializeFloat32Kernels()`](https://zhen-miao.github.io/CoPro/reference/materializeFloat32Kernels.md)
  provides a whole-object compatibility escape hatch. Global, row, and
  column kernel normalization now remain in float32, including
  asymmetric normalized self-kernels. The ordinary centered Frobenius
  objective norm is computed directly from encoded value sums; exact
  whitened Frobenius normalization retains a temporary double-sparse
  compatibility fallback. The operator thread count is now determined
  automatically from the cores actually allocated to the process,
  honoring common HPC scheduler variables (`SLURM_CPUS_PER_TASK`,
  `NSLOTS`, `PBS_NUM_PPN`, `LSB_DJOB_NUMPROC`) and `OMP_NUM_THREADS` so
  a single-core allocation no longer oversubscribes a shared node; set
  `options(CoPro.float32Threads=)` to override.
- One-cell-type skrCCA now solves the symmetric quadratic problem
  directly with an exact symmetric eigendecomposition, selecting the
  largest algebraic eigenvalues and obtaining all requested axes from
  one factorization.
- Sparse within-cell-type kernels are stored as symmetric `dsCMatrix`
  objects, so only one triangle is retained. Cross-type and
  asymmetrically normalized kernels remain general `dgCMatrix` objects.
- Sparse expression PCA now passes centering and scaling vectors
  directly to `irlba` instead of materializing a dense centered matrix.
  All cell types use one common feasible PCA rank, and multiset
  optimizers are also dimension-aware.
- Regression gene scores use the identity `X' (s - mean(s))` and no
  longer construct a centered dense expression matrix for every sigma
  and axis.
- Gene-space CCA now applies self- and cross-covariances as matrix-free
  operators (`Z_i' K_ij (Z_j w)`) instead of storing dense `G x G`
  matrices. Euclidean streaming builds exact sparse fixed-radius kernels
  without dense distance or kernel matrices.
- Kernel normalizers are cached on the CoPro object with matrix
  signatures for safe reuse. Fair-sigma and conditional permutations now
  honor `n_cores` via memory-explicit PSOCK workers. When CoPro is not
  installed (for example under `devtools::load_all()`), `n_cores > 1`
  now falls back to sequential execution with a warning instead of
  aborting inside the worker.

### Internal

- The permutation entry points share one draw-evaluation path in
  `R/D0_permutation_plan.R` instead of each carrying a sequential worker
  plus a near-identical parallel copy.
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
  no longer maintains its own PSOCK cluster block, and three copies of
  the `"pc"` column-shuffling helper and the now-unused
  `.parallelPermutationLapply()` were removed.
- [`irlba::irlba()`](https://rdrr.io/pkg/irlba/man/irlba.html) is no
  longer used anywhere in the package; the five `@importFrom` directives
  naming it (three of which were already stale) were dropped.
  [`irlba::prcomp_irlba()`](https://rdrr.io/pkg/irlba/man/prcomp_irlba.html)
  in
  [`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md)
  is unaffected, and remains the one place where a seed still changes
  results.

## CoPro 1.1.2

### Documentation

- New vignette *Handling very large datasets (Xenium, large MERFISH)*, a
  how-to for keeping large runs in memory: sparse/BPCells input, early
  [`subsetData()`](https://zhen-miao.github.io/CoPro/reference/subsetData.md),
  small `nPCA` for targeted panels, and (most importantly) skipping
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  in favor of the sparse kernel path
  (`computeKernelMatrix(method = "auto"/"sparse")` /
  [`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md)).
- Lowered `nPCA` from 40/30 to 15 in the targeted-panel vignettes (brain
  MERFISH, colon D3, colon D9, organoid) so they follow the documented
  10–15 guidance for imaging panels. Cell scores and normalized
  correlations are unaffected by `nPCA`; the change improves gene-weight
  reproducibility and score transfer, and makes the vignettes consistent
  with the kidney example.

### Internal

- Removed the fully-commented legacy optimizer file and three unused
  per-iteration update-vector helpers superseded by the precomputed
  PC-space operator (`Y_resi`) path. Also removed a stranded file of
  unexported, never-called soft-deprecation wrappers (`newCoPro`,
  `newCoProm`, `subsetDataOne`, `subsetDataMulti`, `computePCAMulti`).
  User-facing error messages that referenced the old `computePCAMulti` /
  `subsetDataMulti` names now point to `computePCA` / `subsetData`.
- Consolidated the `computePCA` documentation topic, which was
  previously split under a stale `computePCAMulti` help page with a
  bogus method alias.

## CoPro 1.1.1

### Performance

- Two-cell-type skrCCA now solves all requested canonical axes with one
  exact SVD of the summed PC-space cross-operator. Multi-slide analyses
  retain the existing memory-efficient summation and do not construct a
  stacked cell matrix or block-diagonal spatial kernel.

### Bug fixes

- Higher skrCCA axes for three or more cell types now use full two-sided
  projection deflation, preventing later axes from reusing directions
  selected by earlier axes under the default scaled-PC formulation.

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
