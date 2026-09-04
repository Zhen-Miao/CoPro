# Inspect saved PC-space SUMCOR solver diagnostics

Reads the numerical diagnostics saved by
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
with `space = "pca", objective = "sumcor"`. No optimization is rerun.

## Usage

``` r
getCCADiagnostics(object, sigma = NULL)
```

## Arguments

- object:

  A fitted `CoProSingle` or `CoProMulti` object.

- sigma:

  Optional single numeric bandwidth. If omitted, return a named list of
  diagnostics for the stored PCA-space bandwidths.

## Value

With `sigma`, a list containing:

- `components`:

  One row per axis: solver, stopping status, convergence, final
  objective, maximum block projected-gradient norm, tolerance,
  iterations, and whether a score denominator reached its floor.

- `conditioning`:

  Per-slide, per-cell-type Gram eigenvalues, numerical rank, condition
  number, and numerical positive-definiteness. Rank counts eigenvalues
  above `n_features * .Machine$double.eps * max_eigenvalue`; that
  threshold is also returned.

- `score_norms`:

  Unfloored score norms at each returned axis, with `floor_active`
  indicating norms at or below `sigma_floor`.

- `objective_traces`:

  Initial and accepted objective values for each full-gradient axis.
  Empty for supplied axes and SUMCOV reductions.

- `coordinate_system`, `sigma_floor`:

  The diagnostic coordinate system (`"whitened_pc"`) and numerical
  score-norm floor. This floor is unrelated to the spatial bandwidth
  named `sigma`.

Without `sigma`, these lists are named by the stored `sigma_<value>`
keys. A PCA-space fit without saved diagnostics (older objects or
SUMCOV) yields `NULL` for that bandwidth; an object without PCA-space
fits yields an empty list. Gene-space fitting creates no diagnostic
record, but any earlier PCA-space results retained in the same object
remain accessible. A requested bandwidth without a PCA-space fit is an
error.

## Details

Residuals and conditioning are measured in the whitened coordinates used
internally, including when the fit used `scalePCs = FALSE`. For later
axes, the gradient is projected off the previously accepted axes as well
as the current direction. A small residual is a first-order condition,
not a global-optimum or biological-recovery certificate.

`solver = "full_gradient"` records the actual stopping reason:
`"gradient_tolerance"`, `"max_iter"`, or `"line_search_stalled"`.
`floor_encountered` includes the initial point and all evaluated
line-search trials, including rejected trials. The warm-start solver is
not monitored. `solver = "sumcov_reduction"` means that the ratio
problem reduced algebraically to SUMCOV. The reduction test is a
numerical one (Gram matrices proportional to the identity to `1e-8`
relative tolerance), so the returned axis is measured against the SUMCOR
residual explicitly and can report `converged = FALSE` with
`status = "residual_above_tolerance"` without the fit itself warning;
read this column rather than assuming the shortcut met `tol`. Iteration
counts are zero for direct one/two-type solutions and unavailable (`NA`)
for the iterative multiblock reduction. Only the returned point is
checked for denominator flooring in that route. `solver = "supplied"`
has `converged = NA`: the axis was not optimized here. Its residual and
score norms are measured on the unit-norm direction of the supplied
weight in whitened coordinates, so they are comparable with the fitted
axes regardless of how the weight was scaled when it was supplied.

The smooth stationarity argument requires positive marginal metrics on
the retained feature space and no active denominator floor. Numerical
full rank alone does not establish good conditioning or the statistical
assumptions of a recovery theorem. Inspect the smallest eigenvalues,
condition numbers, score norms, and residuals together. These
diagnostics neither add a ridge nor filter additional slides or
features.

Recording the diagnostics never changes whether a fit succeeds: if the
record cannot be built, the bandwidth keeps its fitted weights, the
record is `NULL`, and a separate warning says so.

## See also

[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`getCCAObjective()`](https://zhen-miao.github.io/CoPro/reference/getCCAObjective.md)

Other scores-and-correlation:
[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md),
[`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md),
[`getCCAObjective()`](https://zhen-miao.github.io/CoPro/reference/getCCAObjective.md),
[`getNormalizerInfo()`](https://zhen-miao.github.io/CoPro/reference/getNormalizerInfo.md)

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
obj <- computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE)
#> Warning: The kernel at sigma = 0.1 is predicted to be 52% dense, so a fixed-radius sparse kernel saves little or no memory here (sparse storage costs more than dense past ~67% density).
#>   Sigma below about 0.0756 would keep it under 30%. Use detectSigmaRange() to pick sigma from the data, or build one sigma at a time to bound peak memory.
obj <- runSkrCCA(obj, objective = "sumcor", nCC = 2)
#> Objective: sumcor (per-slide self-normalized), slideWeight = "equal", 1 slides.
#> Running skrCCA [1/1] for sigma = 0.1 ...
#> skrCCA finished 1 sigma value(s) in 0.0 s.
#> Optimization succeeded for 1 sigma value(s): sigma_0.1
diagnostics <- getCCADiagnostics(obj, sigma = 0.1)
diagnostics$components
#>   component           solver             status converged objective
#> 1         1 sumcov_reduction gradient_tolerance      TRUE 1.8382729
#> 2         2 sumcov_reduction gradient_tolerance      TRUE 0.3218882
#>   gradient_norm tolerance iterations floor_encountered
#> 1  1.200960e-15     1e-05          0             FALSE
#> 2  7.699152e-16     1e-05          0             FALSE
diagnostics$conditioning
#>           slide  cell_type n_cells n_features rank min_eigenvalue
#> 1 .single_slide Epithelial     100         10   10             99
#> 2 .single_slide Fibroblast     100         10   10             99
#>   max_eigenvalue eigenvalue_tolerance condition_number positive_definite
#> 1             99         2.198242e-13                1              TRUE
#> 2             99         2.198242e-13                1              TRUE
# }
```
