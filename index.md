# CoPro ![CoPro logo](reference/figures/copro-refined-logo-final.jpg)

**CoPro** (Co-Progression) is an R package for detecting co-progression
between cell types in spatial transcriptomics data. It works in both
supervised and unsupervised settings, enabling:

- **Cross-cell-type co-progression**: Identify coordinated gene
  expression patterns between different cell types based on spatial
  proximity
- **Within-cell-type spatial patterns**: Detect tissue
  structure-associated cellular programs within a single cell type
- **Multi-slide analysis**: Analyze patterns consistently across
  multiple tissue slides

For detailed tutorials, function reference, and worked examples with
real datasets, visit the **[CoPro documentation
website](https://zhen-miao.github.io/CoPro/)**.

## Installation

You can install CoPro from [GitHub](https://github.com/Zhen-Miao/CoPro)
with:

``` r

# install.packages("devtools")
devtools::install_github("Zhen-Miao/CoPro")
```

## Quick Start

``` r

library(CoPro)

# Create a CoPro object from your data
obj <- newCoProSingle(
  normalizedData = your_expression_matrix,  # cells x genes
  locationData = your_location_data,        # data.frame with x, y columns
  metaData = your_metadata,                 # data.frame with cell annotations
  cellTypes = your_cell_type_labels         # character vector
)

# Run the analysis pipeline
obj <- subsetData(obj, cellTypesOfInterest = c("TypeA", "TypeB"))
obj <- computePCA(obj, nPCA = 30)

# Pick sigma from the data. sigma is a distance, so it lives in whatever units
# locationData used; detectSigmaRange() reports the sigmas that give the median
# cell 5-20 effective neighbors.
sigmaRange <- detectSigmaRange(obj)
sigmaRange                       # per-block diagnostics + recommended grid

# No computeDistance() needed: the sparse routes build kernels straight from
# coordinates, and method = "auto" picks float32 for anything large.
obj <- computeKernelMatrix(obj, sigmaValues = sigmaRange$sigmaValues)
obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2)
obj <- computeNormalizedCorrelation(obj)
obj <- computeGeneAndCellScores(obj)
obj <- computeRegressionGeneScores(obj)

# Get results
cell_scores <- getCellScores(obj, sigma = obj@sigmaValueChoice, cellType = "TypeA")
```

Gene weights come in two flavors: `@geneScores` (PCA back-projection,
from
[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md))
and `@geneScoresRegression` (from
[`computeRegressionGeneScores()`](https://zhen-miao.github.io/CoPro/reference/computeRegressionGeneScores.md)).
Prefer the regression weights — they avoid PCA collinearity, are
insensitive to the `nPCA` choice, and reproduce better across
replicates.

For `CoProMulti`, the recommended PC-space workflow centers and scales
genes within each `(slide, cell type)` block before one shared PCA, then
fits
`runSkrCCA(space = "pca", objective = "sumcor", slideWeight = "equal")`.
These are the current defaults. The full ratio gradient supports a
first-order stationarity claim under positive marginal conditioning; it
does not establish global optimality or biological recovery. The
gene-space SUMCOR alternative uses a frozen-denominator surrogate and
does not inherit that guarantee.

Inspect saved numerical diagnostics after fitting:

``` r

diagnostics <- getCCADiagnostics(obj)  # one record per fitted bandwidth
diagnostics[[1]]$components           # status, residual, objective, floor use
diagnostics[[1]]$conditioning         # per-slide/type Gram eigenvalues and rank
```

These diagnostics are available for PC-space SUMCOR fits; older objects,
SUMCOV fits, and gene-space fits have no record. Equal slide weights are
equal nominal slide/pair coefficients, not necessarily equal specimen
weights or removal of implicit kernel-size effects. See the [recommended
multi-sample
walkthrough](https://zhen-miao.github.io/CoPro/articles/multi_sample_workflow.html)
for a complete example. Pass `center_per_slide = FALSE` and
`objective = "sumcov"` to reproduce the legacy pooled workflow.

For comparable scores on new slides, freeze the fitted training
reference and apply it unchanged to a subsetted target object:

``` r

score_reference <- fit_score_reference(obj)
target_obj <- subsetData(target_obj, c("TypeA", "TypeB"))
target_scores <- predict(score_reference, target_obj)
```

This target-invariant frozen log-expression transfer is the recommended
default after an internal benchmark. Use
[`getTransferCellScores()`](https://zhen-miao.github.io/CoPro/reference/getTransferCellScores.md)
when cross-platform quantile normalization is specifically intended.
Under the multi-slide default `center_per_slide = TRUE`, a frozen
reference collapses the per-slide block moments to one center and scale
per cell type, so transferred scores are comparable to each other but
not to
[`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md)
on the fitted object; see
[`?fit_score_reference`](https://zhen-miao.github.io/CoPro/reference/fit_score_reference.md).

See the **[Getting Started
vignettes](https://zhen-miao.github.io/CoPro/articles/)** for complete
walkthroughs with real spatial transcriptomics datasets, including
within-cell-type analysis, cross-cell-type co-progression, and
multi-slide experiments.

## API naming

The two naming styles in the API mark two layers. `camelCase` is the
object pipeline — functions that take a `CoProSingle` or `CoProMulti`
and return one, so they chain
([`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md)).
`snake_case` is the engine and utility layer — functions that work on
plain matrices and data frames with no CoPro object involved
([`optimize_bilinear()`](https://zhen-miao.github.io/CoPro/reference/optimize_bilinear.md),
[`optimize_sumcor_pca()`](https://zhen-miao.github.io/CoPro/reference/optimize_sumcor_pca.md),
[`resample_spatial()`](https://zhen-miao.github.io/CoPro/reference/resample_spatial.md),
[`quantile_normalize()`](https://zhen-miao.github.io/CoPro/reference/quantile_normalize.md)).
Reach for the latter when you want the numerical core without the object
wrapper. See
[`?CoPro`](https://zhen-miao.github.io/CoPro/reference/CoPro-package.md)
for the full rule and documented exceptions, including fitted score
references.

## Citation

If you use CoPro in your research, please cite:

> Miao Z, Qu Y, Huang S, Laux L, Peters S, Aristel A, Zhang Z,
> Niedernhofer L, McMahon A, Kim J, Zhang NR (2026). *Dissecting the
> coordinated progression of cell states in spatial transcriptomics with
> CoPro.* bioRxiv 2026.04.17.719309. doi:
> [10.64898/2026.04.17.719309](https://doi.org/10.64898/2026.04.17.719309)

## Getting Help

- Report bugs and request features on [GitHub
  Issues](https://github.com/Zhen-Miao/CoPro/issues)
- Browse the [documentation website](https://zhen-miao.github.io/CoPro/)
  for function reference and tutorials
