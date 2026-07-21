# Run gene-space canonical correlation analysis

Batch-robust CCA that operates directly in gene space instead of PCA
space. For multi-slide data, uses average per-slide canonical
correlation where each slide's contribution is normalized by its own
score variance, preventing batch-level mean shifts from inflating the
objective.

## Usage

``` r
runGeneSpaceCCA(
  object,
  sigma,
  nCC = 2,
  clip = "quantile",
  min_prevalence = 0.008,
  min_cells = 20,
  max_iter = 3000,
  tol = 1e-06,
  streaming = FALSE,
  distanceArgs = list(),
  kernelArgs = list(),
  verbose = TRUE
)

# S4 method for class 'CoPro'
runGeneSpaceCCA(
  object,
  sigma,
  nCC = 2,
  clip = "quantile",
  min_prevalence = 0.008,
  min_cells = 20,
  max_iter = 3000,
  tol = 1e-06,
  streaming = FALSE,
  distanceArgs = list(),
  kernelArgs = list(),
  verbose = TRUE
)

# S4 method for class 'CoProMulti'
runGeneSpaceCCA(
  object,
  sigma,
  nCC = 2,
  clip = "quantile",
  min_prevalence = 0.008,
  min_cells = 20,
  max_iter = 3000,
  tol = 1e-06,
  streaming = FALSE,
  distanceArgs = list(),
  kernelArgs = list(),
  verbose = TRUE
)
```

## Arguments

- object:

  A CoPro object with kernel matrices computed.

- sigma:

  Sigma value to use (single numeric value). Must match a sigma for
  which kernels were computed.

- nCC:

  Number of canonical components (default 2).

- clip:

  Clipping method for gene expression: `"quantile"` for 98th percentile
  (default), or a numeric value for a fixed threshold (e.g., 1 for count
  data).

- min_prevalence:

  Minimum fraction of cells expressing a gene (default 0.008 = 0.8
  percent).

- min_cells:

  Minimum number of cells expressing a gene (default 20).

- max_iter:

  Maximum iterations per component (default 3000).

- tol:

  Convergence tolerance (default 1e-6).

- streaming:

  Logical. If `TRUE`, fuse distance + kernel + covariance reduction into
  a per-slide loop and free the n x n matrices between slides. Bypasses
  [`computeDistance`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  and
  [`computeKernelMatrix`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md):
  `object@distances` and `object@kernelMatrices` are not populated.
  Default `FALSE`.

- distanceArgs:

  Named list of distance parameters passed through to the streaming path
  (e.g., `distType`, `normalizeDistance`, `normalizeTarget`,
  `truncateLowDist`, scaling factors, morphology-aware parameters).
  Ignored when `streaming = FALSE`.

- kernelArgs:

  Named list of kernel parameters passed through to the streaming path
  (e.g., `lowerLimit`, `upperQuantile`, `normalizeKernel`). Ignored when
  `streaming = FALSE`.

- verbose:

  Print progress messages (default TRUE).

## Value

The CoPro object with gene weights in `geneScores`, cell scores in
`cellScores`, and weight vectors in `skrCCAOut`.

## Details

Requires kernel matrices to already be computed via
[`computeKernelMatrix`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md).
Does NOT require
[`computePCA`](https://zhen-miao.github.io/CoPro/reference/computePCA.md).

The objective maximized is: \$\$f\_{avg}(w) = \frac{1}{S}
\sum\_{s=1}^{S} \sum\_{A\<B} \frac{w_A^\top C\_{AB}^{(s)}
w_B}{\sigma_A^{(s)} \sigma_B^{(s)}}\$\$ where \\\sigma_A^{(s)} =
\sqrt{w_A^\top C\_{AA}^{(s)} w_A}\\ is the per-slide score standard
deviation. Subsequent components use Gram-Schmidt deflation in weight
space. The power iteration uses a frozen-sigma surrogate (sigma values
held fixed at the previous iterate when computing each weight update),
making the algorithm an ALS-style alternating maximization rather than
exact coordinate ascent.

Memory scales as \\O(G^2 \times S \times (C + C(C-1)/2))\\ for
precomputed covariance matrices: \\S\\ slides, each storing \\C\\
self-covariances and \\C(C-1)/2\\ cross-covariances of size \\G \times
G\\. For example G=5000, S=10, C=3 gives \\10 \times 6\\ matrices of
\\200\\ MB each, approximately 12 GB.

Per-slide gene handling: each (slide, cell type) expression matrix is
independently centered and scaled. A gene with zero variance on a
particular (slide, cell type) pair is set to zero on that slide and
contributes nothing from it; the gene is still retained globally if it
passes the prevalence filter, with its weight driven by other slides
where it has variance. Slides where any requested cell type has fewer
than 10 cells are dropped with a warning.

Storage: gene-space CCA results live in `@skrCCAOut` under the
`gscca_sigma_<value>` key (distinct from `runSkrCCA`'s `sigma_<value>`
keys) so the two CCA flavors do not collide.
[`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md)
only operates on `runSkrCCA` (PCA-space) outputs; gene-space CCA already
populates `@geneScores` and `@cellScores` directly, so no further call
is needed.

Streaming mode (`streaming = TRUE`) reduces peak memory to roughly one
slide's pairwise n x n footprint at a time. Distance normalization
defaults to `normalizationScope = "global"` (a single factor across all
slides, matching the slot-based pipeline's semantics under
`normalizeDistance = TRUE`); under a fixed RNG seed the streaming result
is bit-identical to the slot-based path. Pass
`distanceArgs = list(normalizationScope = "per_slide")` to opt into
per-slide factors instead – useful when slides have different coordinate
scales, but be aware the per-slide percentile jitter can perturb
degenerate canonical components on heterogeneous datasets. Use
`streaming = FALSE` when downstream code reads `object@distances` or
`object@kernelMatrices`.

## See also

[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md)

Other spatial-pipeline:
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md),
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md),
[`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md),
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
