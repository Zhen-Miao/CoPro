# Compute PCA on Single- or Multi-Slide Data

Performs PCA on the normalized data stored within a `CoProSingle` or
`CoProMulti` object. Multi-slide data must share the same gene columns;
the default preprocessing removes slide-level location and scale
internally.

## Usage

``` r
computePCA(
  object,
  nPCA = 40,
  center = TRUE,
  scale. = TRUE,
  scalePCs = TRUE,
  dataUse = "raw",
  center_per_slide = TRUE
)

# S4 method for class 'CoProSingle'
computePCA(object, nPCA = 40, center = TRUE, scale. = TRUE, scalePCs = TRUE)

# S4 method for class 'CoProMulti'
computePCA(
  object,
  nPCA = 40,
  center = TRUE,
  scale. = TRUE,
  scalePCs = TRUE,
  dataUse = "raw",
  center_per_slide = TRUE
)
```

## Arguments

- object:

  A `CoProSingle` or `CoProMulti` object.

- nPCA:

  Number of principal components to compute for each cell type.

- center:

  Whether to center the matrix before PCA

- scale.:

  Whether to scale the matrix before PCA

- scalePCs:

  Whether to scale (whiten) PCs by their standard deviation before
  downstream CCA optimization. Default `TRUE` (recommended).

- dataUse:

  What data to use, choices between "raw" and "integrated". Default is
  "raw". For single slide, this argument is ignored.

- center_per_slide:

  Apply `center` and `scale.` within each (slide, cell type) block
  before fitting one shared PCA. Default `TRUE`.

## Value

The input object with `pcaGlobal` populated. For `CoProMulti`,
`pcaResults` has structure
`list(slideID = list(cellType = score-row view))`.

## Details

For cell type \\i\\ and slide \\s\\, the recommended multi-slide path
forms \$\$Z_i^{(s)} = (X_i^{(s)} - 1\mu_i^{(s)'})D_i^{(s)-1},\$\$ stacks
the \\Z_i^{(s)}\\ by rows, and computes one truncated SVD \\Z_i =
U_i\Delta_i V_i'\\. Every slide therefore uses the same loading matrix
\\V_i\\, with score block \\Z_i^{(s)}V_i\\. Because each standardized
gene column has zero mean within its block, those PC scores are also
centered within slide automatically; no post-PCA recentering is needed.

The shared loading is in **within-slide standardized gene coordinates**.
If the stored slide scales differ, there is intentionally no single
equivalent coefficient vector in raw expression units.

Consequently the retained subspace is invariant to any per-slide,
gene-wise affine map with positive multipliers – the batch-effect family
this is meant to absorb – provided `center` and `scale.` are both `TRUE`
and no gene trips the degeneracy guard. That guard pins \\d = 1\\ for
genes whose standard deviation is below `1e-3` or whose nonzero fraction
is below 1%, so that dividing by a near-zero scale cannot amplify noise.
It is evaluated on the **raw** block, so exact affine invariance does
not extend to guarded genes: an additive shift makes every entry nonzero
and can suppress a guard that the unshifted data would trip. A gene
guarded on any one slide is guarded on all of them, which keeps slides
mutually comparable; the alternative – a per-block decision – would
leave a gene standardized on one slide and raw on another, reintroducing
a per-slide scale difference in precisely the low-detection genes whose
detection rate is often itself the batch effect.

## See also

[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)

Other spatial-pipeline:
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md),
[`computeSparseKernelFloat32()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernelFloat32.md),
[`detectSigmaRange()`](https://zhen-miao.github.io/CoPro/reference/detectSigmaRange.md),
[`runGeneSpaceCCA()`](https://zhen-miao.github.io/CoPro/reference/runGeneSpaceCCA.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`runSlideLevelInference()`](https://zhen-miao.github.io/CoPro/reference/runSlideLevelInference.md),
[`selectSigmaByPermutation()`](https://zhen-miao.github.io/CoPro/reference/selectSigmaByPermutation.md)

## Examples

``` r
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
```
