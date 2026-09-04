# Recommended multi-sample workflow: PC-space SUMCOR

CoPro learns shared co-progression directions from several spatial
samples. The recommended fitting route uses **PC-space SUMCOR**:
standardize genes within each slide and cell type, fit one shared PCA
per cell type, then optimize the sum of per-slide standardized spatial
associations with the full ratio gradient. These are already the
`CoProMulti` defaults. The code below spells them out and shows how to
inspect numerical diagnostics.

This walkthrough uses the three-slide Colon Day 0 dataset. Code is shown
without execution during package builds; run it locally after
downloading the example data. The existing [gene-space colon
example](https://zhen-miao.github.io/CoPro/articles/colon_d0_multi_slide.md)
illustrates the alternative surrogate solver and its spatial plots.

## Prepare the multi-slide object

``` r

library(CoPro)
dat <- readRDS(copro_download_data("colon_d0_multi"))
cell_types <- c("Epithelial", "Fibroblast", "Immune")
obj <- newCoProMulti(
  normalizedData = dat$normalizedData,
  locationData = dat$locationData,
  metaData = dat$metaData,
  cellTypes = dat$cellTypes,
  slideID = dat$slideID
)
obj <- subsetData(obj, cellTypesOfInterest = cell_types)
obj <- computePCA(obj, nPCA = 15, center_per_slide = TRUE)
```

Within-slide centering and scaling happen before the shared PCA. The
feature coordinates are therefore shared standardized-gene coordinates.
Retaining 15 PCs is an illustrative choice for this targeted panel; the
retained space is part of the model, and changing it can change the
fitted directions.

## Construct kernels and fit

``` r

sigma_range <- detectSigmaRange(obj)
sigma_values <- sigma_range$sigmaValues
# One geometry-based bandwidth for this example, in locationData units.
sigma_choice <- sigma_values[ceiling(length(sigma_values) / 2)]
obj <- computeKernelMatrix(obj, sigmaValues = sigma_choice)
obj <- runSkrCCA(
  obj, space = "pca", objective = "sumcor", slideWeight = "equal",
  nCC = 2, maxIter = 1000, tol = 1e-5
)
getCCAObjective(obj)
```

Kernels connect cells within the same slide. `slideWeight = "equal"`
assigns equal nominal coefficients to slide/pair terms. It does not
remove differences in kernel geometry or guarantee equal total weight
per specimen when a specimen has more sections. Record the
slide-to-specimen mapping and the intended weighting when describing a
multi-sample analysis.

The fitting term is `a' K b / (||a|| ||b||)`, with the kernel left at
the
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md)
default (clipped at `upperQuantile`, not normalized), as in the package
quick start. It is invariant to score rescaling, but need not lie
between -1 and 1, and it is not cell-count invariant: replicating cells
of the two types by factors `r_i, r_j` scales it by `sqrt(r_i r_j)`.
`normalizeKernel = TRUE` divides each cross-type block by the median of
its row sums; that changes the scaling to `sqrt(r_i / r_j)` and
reweights the slide/pair terms, so treat it as part of the fitted
objective rather than a cosmetic choice and record it if you use it. The
bandwidth above is chosen from geometry, not claimed to maximize
association.

## Inspect convergence and conditioning

``` r

diagnostics <- getCCADiagnostics(obj, sigma = sigma_choice)
diagnostics$components
diagnostics$conditioning
diagnostics$score_norms
diagnostics$objective_traces$CC1
```

`components` reports the solver, stopping status, final objective,
projected gradient norm, tolerance, iterations, and denominator-floor
use for each axis. For a full-gradient fit, `converged = TRUE` means the
projected gradient is within the requested tolerance. `max_iter` or
`line_search_stalled` records why a fit stopped without meeting it;
inspect the residual before interpreting the axis as numerically
stationary. A larger iteration budget may help an iteration-limited fit,
but does not cure poor conditioning.

`conditioning` reports each slide/type Gram matrix in whitened PC
coordinates, including its smallest eigenvalue, condition number, and
numerical rank. These coordinates are used for diagnostics even with
`scalePCs = FALSE`. Full numerical rank does not ensure good
conditioning. `score_norms` reports the unfloored norm at the returned
direction. A floor encounter can also occur during a rejected
line-search trial, as recorded in `components`.

The smooth stationarity argument assumes positive marginal metrics on
the retained space. A hard denominator floor does not supply that
assumption. If a Gram matrix is singular or poorly conditioned, review
the retained PC space and slide adequacy and document any change. CoPro
does not silently add a ridge through this diagnostic API.

The guarantee is first-order stationarity under the stated conditions.
It does not establish a global maximum or recovery of a population
biological program. Later axes are assessed conditionally on the earlier
accepted axes. An equivalent SUMCOV shortcut is labeled
`sumcov_reduction` and checked for its returned SUMCOR residual. A
supplied axis is labeled `supplied`, with `converged = NA`. Older saved
fits have no diagnostic record until refitted.

## Compute scores and regression gene weights

``` r

obj <- computeNormalizedCorrelation(obj, normalizer = "unwhitened")
obj <- computeGeneAndCellScores(obj)
obj <- computeRegressionGeneScores(obj)
epithelial_scores <- getCellScores(
  obj, sigma = sigma_choice, cellType = "Epithelial"
)
```

Use regression gene weights (`geneScoresRegression`) for manuscript gene
rankings and figures. They describe association of genes with the fitted
cell score; changing the fitting space is separate from choosing these
reported gene weights.

The explicit downstream `unwhitened` normalizer divides by the centered
kernel’s Frobenius norm as well as the score norms. With the default
`aggregateDenominator = "sum"`, scores centered within each slide (the
`center_per_slide = TRUE` PCA above), and a nonzero denominator, the
reported statistic is bounded in absolute value by one on every slide
and in aggregate. `aggregateDenominator = "rss"` is a different
statistic that can exceed one across concordant slides, and a
pooled-centered legacy PCA does not satisfy the within-slide centering
the identity relies on. The normalizer does not change the directions
fitted above. Other downstream normalizers answer different questions;
record which one was used.

## The gene-space alternative

[`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md)
works directly in genes and uses its own filtering, clipping, and
per-slide standardization. Under SUMCOR its current iteration uses
frozen denominators. A converged surrogate fixed point need not be a
stationary point of the ratio objective; the Gauss-Seidel sweep does not
remove this limitation. PC-space fitting is the implemented route for
the full-gradient stationarity claim. This distinction does not
demonstrate that PC-space fitting has better biological recovery, since
changing spaces also changes the model.

This walkthrough covers fitting and numerical assessment. Its
diagnostics do not establish statistical significance or independent
biological replication.
