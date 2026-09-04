# runSkrCCA

runSkrCCA

## Usage

``` r
runSkrCCA(
  object,
  scalePCs = TRUE,
  nCC = 2,
  tol = 1e-05,
  transferred_weight_1 = NULL,
  maxIter = 200,
  sigmaChoice = NULL,
  n_cores = 1,
  step_size = 1,
  space = c("pca", "gene"),
  objective = NULL,
  slideWeight = NULL,
  minCellsPerSlide = 10,
  ...
)

# S4 method for class 'CoPro'
runSkrCCA(
  object,
  scalePCs = TRUE,
  nCC = 2,
  tol = 1e-05,
  transferred_weight_1 = NULL,
  maxIter = 200,
  sigmaChoice = NULL,
  n_cores = 1,
  step_size = 1,
  space = c("pca", "gene"),
  objective = NULL,
  slideWeight = NULL,
  minCellsPerSlide = 10,
  ...
)

# S4 method for class 'CoProMulti'
runSkrCCA(
  object,
  scalePCs = TRUE,
  nCC = 2,
  tol = 1e-05,
  transferred_weight_1 = NULL,
  maxIter = 200,
  sigmaChoice = NULL,
  n_cores = 1,
  step_size = 1,
  space = c("pca", "gene"),
  objective = NULL,
  slideWeight = NULL,
  minCellsPerSlide = 10,
  ...
)
```

## Arguments

- object:

  A CoPro object

- scalePCs:

  Whether to scale each PCs to a uniform variance before running the
  program

- nCC:

  Number of canonical vectors to compute, default = 2

- tol:

  Tolerance for termination, default = 1e-5. Under SUMCOR this is the
  projected-gradient norm; under SUMCOV it is the weight-update change.

- transferred_weight_1:

  If we use cross-slide weight transfer function, the transferred weight
  on each PC. Otherwise, the value should be set to NULL.

- maxIter:

  Maximum optimization iterations.

- sigmaChoice:

  Specific sigma value to use (CoProMulti only, ignored for CoPro)

- n_cores:

  Number of cores for parallel processing (CoProMulti only, ignored for
  CoPro)

- step_size:

  Damping factor in (0, 1\]. Default 1 (undamped). **Honored under both
  objectives**, so a value chosen for stability keeps its effect when
  the objective changes.

  Under `"sumcov"` this is the damped power iteration: values below 1
  blend old and new weights,
  `w \leftarrow \mathrm{normalize}((1-\alpha)w_{old} + \alpha w_{new})`,
  which can help with many cells or many CCs.

  Under `"sumcor"` it shortens every trial step of the
  projected-gradient line search. These are the same operation: because
  both the current iterate and the retracted candidate lie on the
  geodesic through `w` in the search direction, blending them and
  renormalizing *is* a retraction, so damping by \\\alpha\\ equals
  taking a step \\\tau = \alpha t / ((1-\alpha)\sqrt{1 + t^2\\g\\^2} +
  \alpha)\\ along the same arc. Applying it to the trial step rather
  than after the fact keeps the Armijo test on the point actually
  returned, so a damped SUMCOR run stays monotone. It also damps the
  SUMCOV warm start.

- space:

  Which feature space to optimize in. `"pca"` (default) runs the
  PC-space optimizer described here. `"gene"` forwards to
  [`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md),
  which needs a single `sigma` – supply it through `sigmaChoice` – and
  accepts its own arguments (`clip`, `min_prevalence`, `streaming`, ...)
  through `...`. PC-space SUMCOR uses the full ratio gradient. The
  current gene-space SUMCOR solver uses a frozen-denominator surrogate
  and does not inherit the PC-space stationarity guarantee.

  Under `space = "gene"`, `objective`, `tol`, `maxIter` and
  `minCellsPerSlide` are forwarded **only when you supply them**,
  because the two entry points have different defaults (`objective` is
  `"sumcov"` for single-slide PCA space and `"sumcor"` for multi-slide
  PCA and gene space; `tol` is `1e-5` vs `1e-6`; `maxIter` is `200` vs
  `3000`; and `minCellsPerSlide` `10` vs `20`). Forwarding a
  `runSkrCCA()` default would silently make this call differ from the
  [`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md)
  call it claims to be, so an unset argument keeps gene space's own
  default. Arguments with no gene-space meaning (`scalePCs`,
  `transferred_weight_1`, `n_cores`, `step_size`, `slideWeight`) are an
  error rather than silently dropped.

- objective:

  Which canonical criterion to maximize.

  `"sumcov"` maximizes the sum of kernel-smoothed cross-covariances
  \\\sum\_{i\<j} w_i' (\sum_s X_i^{(s)'} K\_{ij}^{(s)} X_j^{(s)}) w_j\\
  under \\\\w_i\\ = 1\\.

  `"sumcor"` maximizes the per-slide self-normalized sum, dividing each
  slide's cross term by that slide's own score scales.

  **With one slide these are often the same problem.** When PCA was fit
  to that slide, whitened PCs give \\X_i'X_i = (n_i-1) I\\, so on
  \\\\w_i\\ = 1\\ the denominators are \\\sigma_i = \sqrt{n_i - 1}\\ –
  note \\\sigma\\ is a norm, not a root-mean-square – and the objective
  is SUMCOV reweighted by the *per-pair* constant \\m\_{ij} /
  \sqrt{(n_i-1)(n_j-1)}\\. A per-pair constant leaves the maximizer
  alone only when every pair gets the same one, so the reduction is
  exact for **one or two cell types at any cell counts**, and for
  **three or more only when the counts are equal**. In exactly those
  cases CoPro uses the existing SUMCOV solver. Outside them it runs the
  full-gradient SUMCOR optimizer, so an explicit `objective = "sumcor"`
  always optimizes the requested criterion. The matrix condition is
  checked directly; this matters when several slides informed the PCA
  but filtering leaves only one slide for CCA.

  Across slides they generally differ. SUMCOV factors exactly as
  \\\sum_s \sigma_i^{(s)} \sigma_j^{(s)} \rho\_{ij}^{(s)}\\ – with no
  \\\sqrt{n_i n_j}\\ factor, since \\\sigma_i \sigma_j \rho\_{ij}\\ is
  the SUMCOV term already – so it already sums per-slide correlations,
  weighted by per-slide score scale. That scale factor is what lets a
  slide with inflated variance along the canonical direction dominate;
  `"sumcor"` removes it, and `slideWeight = "size"` adds the explicit
  cell-count factor \\\sqrt{n_i n_j}\\ on its own. The default is
  `"sumcov"` for single-slide objects and `"sumcor"` for multi-slide
  objects. The latter is paired with the within-slide PCA preprocessing
  that
  [`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md)
  now uses by default.

  The cell-level permutation routines build their null with whichever
  criterion the weights were fitted under, so SUMCOR weights no longer
  have to be re-fit before testing. They remain single-slide only:
  `CoProMulti` is directed to
  [`runSlideLevelInference()`](https://zhen-miao.github.io/CoPro/reference/runSlideLevelInference.md)
  whatever the objective. On one slide, where the SUMCOR reduction holds
  (one or two cell types at any counts, three or more at equal counts)
  the fitted weights *are* the SUMCOV maximizer and the existing SUMCOV
  null already matches; otherwise
  [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
  re-optimizes each draw under SUMCOR.
  [`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md)
  and
  [`runSkrCCAPermu_Conditional()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_Conditional.md)
  are restricted to at most two cell types, where only the reducible
  case can arise, and error rather than mix criteria if that is ever
  relaxed.

- slideWeight:

  Per-slide weighting, only valid with `objective = "sumcor"` (an error
  otherwise, since under `"sumcov"` the weighting is fixed by the
  objective). `"equal"` (default) gives equal nominal coefficients to
  slide/pair terms, matching
  [`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md).
  Kernel scaling can still affect their relative influence. Equal slide
  weights are not equal specimen weights when specimens have different
  numbers of sections. `"size"` weights slide `s` by \\\sqrt{n_i^{(s)}
  n_j^{(s)}}\\, so larger slides count for more without per-slide
  variance re-entering.

- minCellsPerSlide:

  Minimum cells per (slide, cell type). Slides below this are
  **dropped** under `objective = "sumcor"`, which divides by a per-slide
  scale that a near-empty slide drives to its floor. Under `"sumcov"`
  they are only reported, not dropped: a thin slide simply contributes
  little to the summed operator, and dropping it would change results
  computed before this rule existed.

- ...:

  Passed to
  [`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md)
  when `space = "gene"`; ignored otherwise.

## Value

CoPro object with skrCCA results computed. The objective actually used
is recorded on `@skrCCAOut` and readable with
[`getCCAObjective()`](https://zhen-miao.github.io/CoPro/reference/getCCAObjective.md).
PC-space SUMCOR numerical diagnostics are readable with
[`getCCADiagnostics()`](https://zhen-miao.github.io/CoPro/reference/getCCADiagnostics.md).

## Batch effects

`objective = "sumcor"` removes per-slide score-scale sensitivity but not
a per-slide *mean shift* on its own. The multi-slide default
`computePCA(center_per_slide = TRUE)` centers and scales each (slide,
cell type) expression block before the shared PCA, so a technical shift
is excluded from both the PC scores and the covariance that chooses the
retained loading space. With the legacy
`computePCA(center_per_slide = FALSE)`, a shared technical shift leaves
\\u_i^{(s)} \approx M_i^{(s)} \mathbf{1} + \epsilon\\ and the numerator
picks up \\M_i M_j \mathbf{1}'K\mathbf{1}\\, positive whenever both cell
types shift the same way. Dividing by `sigma` does not rescue this: for
a nonnegative kernel the leading singular pair is close to the Perron
vector, so a constant score reaches \\\rho \approx \sigma\_{max}(K)\\ on
every slide. Thus batch-robust multi-slide analysis needs both defaults:
within-slide preprocessing before PCA and the per-slide self-normalized
objective.

## Numerical diagnostics

The recommended multi-slide route is
`space = "pca", objective = "sumcor"` with the default within-slide PCA
preprocessing. It uses the full ratio gradient and an Armijo line
search, or an algebraically equivalent SUMCOV reduction. Inspect
[`getCCADiagnostics()`](https://zhen-miao.github.io/CoPro/reference/getCCADiagnostics.md)
after fitting: it records the stopping status, projected-gradient
residual, accepted objective trace, per-slide Gram conditioning, and
denominator-floor use for each axis. Convergence is a first-order
numerical condition, not a global optimum or a biological-recovery
guarantee. The smooth stationarity argument assumes positive marginal
metrics on the retained space; an active hard floor does not satisfy
that assumption. Diagnostics do not change the objective, regularize the
Gram matrices, or choose a new feature space.

## See also

[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md),
[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md),
[`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md),
[`getCCAObjective()`](https://zhen-miao.github.io/CoPro/reference/getCCAObjective.md),
[`getCCADiagnostics()`](https://zhen-miao.github.io/CoPro/reference/getCCADiagnostics.md)

Other spatial-pipeline:
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md),
[`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md),
[`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md),
[`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md),
[`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md),
[`runSlideLevelInference()`](https://zhen-miao.github.io/CoPro/reference/runSlideLevelInference.md),
[`selectSigmaByPermutation()`](https://zhen-miao.github.io/CoPro/reference/selectSigmaByPermutation.md)

## Examples

``` r
# \donttest{
toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
obj <- newCoProSingle(
  normalizedData = toy$normalizedData,
  locationData   = toy$locationData,
  metaData       = toy$metaData,
  cellTypes      = toy$cellTypes
)
obj <- subsetData(obj, cellTypesOfInterest = unique(toy$cellTypes))
obj <- computePCA(obj, nPCA = 10)
#> Input is dense (matrix), performing irlba PCA...
#> Input is dense (matrix), performing irlba PCA...
obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1), verbose = FALSE)
obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2)
#> Running skrCCA [1/2] for sigma = 0.05 ...
#> Running skrCCA [2/2] for sigma = 0.1 ...
#> skrCCA finished 2 sigma value(s) in 0.0 s.
#> Optimization succeeded for 2 sigma value(s): sigma_0.05, sigma_0.1
# }
```
