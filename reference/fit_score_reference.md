# Fit a frozen CoPro score-transfer reference

Builds a self-contained reference for transferring fitted CoPro axes to
new slides. For each cell type, the function freezes the
training-expression mean and standard deviation together with the exact
PCA back-projected CoPro weights.
[`predict()`](https://rdrr.io/r/stats/predict.html) then applies
\$\$((X\_{target} - \mu\_{train}) / \sigma\_{train}) W\$\$ without
estimating anything from the target cohort.

## Usage

``` r
fit_score_reference(
  object,
  sigma = NULL,
  reference_weight = c("cell_pooled", "equal_slide")
)

# S3 method for class 'CoProScoreReference'
predict(object, newdata, aggregate = FALSE, chunk_size = 20000L, ...)

# S3 method for class 'CoProScoreReference'
print(x, ...)
```

## Arguments

- object:

  For `fit_score_reference()`, a fitted `CoProSingle` or `CoProMulti`
  object. Run
  [`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md)
  first; the exact PCA back-projected weights in `object@geneScores` are
  used for scoring. For the
  [`predict()`](https://rdrr.io/r/stats/predict.html) method, a
  `CoProScoreReference` returned by `fit_score_reference()`.

- sigma:

  Numeric spatial bandwidth identifying the fitted CoPro weights. By
  default, uses the single value in `object@sigmaValueChoice`.

- reference_weight:

  How training moments are combined. `"cell_pooled"` (default) computes
  ordinary cell-pooled means and sample standard deviations.
  `"equal_slide"` averages each training slide's first and second
  moments with equal weight. Both routes retain CoPro's PCA safety
  guard: scales are pinned to 1 for genes with standard deviation below
  `1e-3` or nonzero fraction below 1 percent.

  The two routes use different variance conventions. `"cell_pooled"`
  divides by `n - 1`; `"equal_slide"` is the population second moment of
  the equal-weight slide mixture. On a single training slide they
  therefore differ by a factor of `sqrt(n / (n - 1))` rather than
  coinciding. Both conventions are the ones the selection benchmark
  used, so they are kept as they are.

- newdata:

  A subsetted `CoProSingle` or `CoProMulti` target object using the same
  normalized-expression representation and genes as the reference. No
  PCA, kernel, or CCA fit is required on the target.

- aggregate:

  If `FALSE` (default), return a named list of score matrices by cell
  type. If `TRUE`, return one matrix in target-cell order.

- chunk_size:

  Positive integer number of target cells scored at once.

- ...:

  Unused.

- x:

  A `CoProScoreReference` object.

## Value

`fit_score_reference()` returns a `CoProScoreReference` object
containing the frozen transform, weights, and training provenance.

`predict.CoProScoreReference()` returns target cell scores as a named
list or an aggregated matrix.

## Details

This frozen log-expression workflow is the recommended default for
cross-slide score transfer after an internal benchmark. It preserves a
fixed out-of-sample map: adding or removing other target slides does not
change an existing cell's score.

The reference and target must use the same normalized-expression
representation and modeled gene panel. This function deliberately does
not quantile-normalize the target or recompute target-specific moments.
For cross-platform quantile normalization, use
[`getTransferCellScores()`](https://zhen-miao.github.io/CoPro/reference/getTransferCellScores.md)
instead.

Regression gene scores are useful for interpretation, but they are not
the canonical scoring map. This function therefore always uses the exact
back-projected weights from
[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md).

**Frozen scores and the fitted object's own scores.** Under the
`CoProMulti` default,
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md)
standardizes each (slide, cell type) block with its own center and scale
(`center_per_slide = TRUE`, recorded as
`preprocessing = "within_slide"`). A frozen reference must collapse to
one center and scale per cell type, because an unseen target slide has
no stored block moments. Transferred scores stay mutually comparable and
target-invariant, but they are not on the same affine footing as
[`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md)
on the fitted object, so the two should not be pooled onto one axis
without care. Applying the frozen map back to the training cells
reproduces
[`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md)
exactly only under pooled preprocessing – `CoProSingle`, or
`computePCA(center_per_slide = FALSE)`.
`reference_weight = "equal_slide"` does not change this, since it still
yields a single shared center and scale. `fit_score_reference()` emits a
message when it detects within-slide preprocessing, and records what it
saw in the returned object's `preprocessing` field.

The scale guard is CoPro's PCA preprocessing rule – low variance *or*
low prevalence – rather than a bare floor against division by zero, so a
frozen reference guards exactly the genes the PCA behind its weights
guarded. That correspondence is what makes self-transfer exact under
pooled preprocessing.

For a sparse target, prediction evaluates the affine map as a sparse
matrix multiplication plus a component-level offset. It does not
materialize a dense cells-by-genes chunk.

## See also

[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md),
[`getTransferCellScores()`](https://zhen-miao.github.io/CoPro/reference/getTransferCellScores.md)

## Examples

``` r
genes <- c("g1", "g2")
training_x <- cbind(g1 = 1:20, g2 = c(2:11, 1:10))
rownames(training_x) <- paste0("train", 1:20)
training_type <- rep(c("A", "B"), each = 10)
training <- newCoProSingle(
  normalizedData = training_x,
  locationData = data.frame(
    x = 1:20, y = 1:20, row.names = rownames(training_x)
  ),
  metaData = data.frame(row.names = rownames(training_x)),
  cellTypes = training_type
)
training <- subsetData(training, c("A", "B"))

# In a real analysis these weights come from computeGeneAndCellScores().
training@geneScores <- list(
  "geneScores|sigma0.1|A" = matrix(
    c(1, -1), ncol = 1, dimnames = list(genes, "CC_1")
  ),
  "geneScores|sigma0.1|B" = matrix(
    c(0.5, 0.5), ncol = 1, dimnames = list(genes, "CC_1")
  )
)
training@sigmaValueChoice <- 0.1
reference <- fit_score_reference(training)

target_x <- training_x[1:10, ] + 1
rownames(target_x) <- paste0("target", 1:10)
target_type <- rep(c("A", "B"), each = 5)
target <- newCoProSingle(
  normalizedData = target_x,
  locationData = data.frame(
    x = 1:10, y = 1:10, row.names = rownames(target_x)
  ),
  metaData = data.frame(row.names = rownames(target_x)),
  cellTypes = target_type
)
target <- subsetData(target, c("A", "B"))
transferred <- predict(reference, target)
```
