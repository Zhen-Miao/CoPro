# Handling very large datasets (Xenium, large MERFISH)

## Overview

Imaging-based platforms such as **10x Xenium** and large **MERFISH**
panels routinely produce hundreds of thousands to millions of cells per
slide. Two things blow up with the number of cells `n`:

1.  **The expression matrix** (`n` cells x genes).
2.  **The spatial distance / kernel matrix**, which is `n x n` per
    cell-type block if built densely — the dominant cost at scale.

CoPro handles each with a memory-efficient path. This guide collects the
settings that keep a large run within memory while producing results
that are numerically identical to the classic dense pipeline. The
scientific interpretation of the output (cell scores, gene scores,
normalized correlation) is unchanged; only how they are computed
differs.

The single most important change for large data is: **do not call
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
— use the sparse kernel path instead.**

## 1. Store expression sparsely or on disk

[`newCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/newCoProSingle.md)
and
[`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md)
accept three matrix types for `normalizedData`:

- a base `matrix` (dense),
- a **sparse** `Matrix::dgCMatrix` (recommended for imaging panels,
  which are mostly zeros), or
- a **BPCells** `IterableMatrix`, which keeps the counts on disk and
  streams them (recommended when the matrix does not fit in RAM).

``` r

library(CoPro)

# Sparse in-memory matrix (genes are columns, cells are rows):
obj <- newCoProSingle(
  normalizedData = sparse_expr,     # a cells x genes dgCMatrix
  locationData   = loc_df,          # data.frame with x, y (and optional z)
  metaData       = meta_df,
  cellTypes      = cell_type_labels
)

# On-disk with BPCells (install with remotes::install_github("bnprks/BPCells/r")):
# `bp_mat` is a BPCells IterableMatrix; computePCA() below will call
# BPCells::svds() so the PCA never materializes a dense matrix.
obj <- newCoProSingle(
  normalizedData = bp_mat,
  locationData   = loc_df,
  metaData       = meta_df,
  cellTypes      = cell_type_labels
)
```

## 2. Subset to the cell types of interest early

A Xenium slide may contain dozens of cell types, but a co-progression
analysis usually targets two (cross-type) or one (within-type).
[`subsetData()`](https://zhen-miao.github.io/CoPro/reference/subsetData.md)
drops all other cells, so every downstream step — PCA, kernels, skrCCA —
works on a much smaller problem. Do it before
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md).

``` r

obj <- subsetData(obj, cellTypesOfInterest = c("Tumor", "Fibroblast"))
```

## 3. PCA: use a small `nPCA` for targeted panels

For targeted imaging panels (Xenium is typically 300–500 genes; MERFISH
300–1,000), use **`nPCA = 10–15`**. Too many PCs on a small panel add
unstable noise dimensions that hurt gene-weight reproducibility and
score transfer, without improving cell scores. (For full-transcriptome
scRNA-seq, use 20–30.) The default `nPCA = 40` is intended for large
panels; lower it for Xenium.

When `normalizedData` is a BPCells `IterableMatrix`,
[`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md)
automatically runs `BPCells::svds()`, so the PCA step is itself
out-of-core.

``` r

obj <- computePCA(obj, nPCA = 15, center = TRUE, scale. = TRUE)
```

## 4. Build kernels the sparse way — skip `computeDistance()`

The classic pipeline is
[`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
→
[`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
which forms a dense `n x n` distance matrix. For 100,000 cells that is
~80 GB per block, which is infeasible.

The Gaussian kernel `exp(-0.5 (d / sigma)^2)` is essentially zero once
`d > sigma * sqrt(-2 log(lowerLimit))`. The **sparse path** enumerates
only the cell pairs inside that radius, using a grid-based fixed-radius
neighbor search, and never forms the dense matrix. The result is
numerically identical to the dense path (every entry beyond the radius
is zero anyway).

You get the sparse path in one of two equivalent ways:

``` r

# Option A (recommended): let computeKernelMatrix choose.
# method = "auto" (the default) picks the sparse path when any cell-type block
# reaches autoThreshold (5000) cells, or the aggregate dense workload reaches
# autoThreshold^2 entries. dropDistances = TRUE (default) frees the @distances
# slot afterward. NOTE: do NOT call computeDistance() first on large data — the
# sparse path builds kernels directly from coordinates and would ignore it.
obj <- computeKernelMatrix(
  obj,
  sigmaValues   = c(0.05, 0.1, 0.2),  # required; no default is derived here
  method        = "auto",             # or force "sparse"
  dropDistances = TRUE
)

# Option B (explicit): call the fused sparse constructor directly.
obj <- computeSparseKernel(obj, sigmaValues = c(0.05, 0.1, 0.2))
```

Notes:

- **Force it** with `method = "sparse"` if you want the sparse path
  regardless of size, or raise/lower `autoThreshold` to move the
  automatic cutoff. Note that `method = "auto"` only skips
  [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  when it actually selects the sparse path; if a block falls below the
  threshold it picks `"dense"` and then needs distances. On a machine
  where the dense matrix would not fit, pass `method = "sparse"`
  explicitly so the choice never depends on the threshold.
- The sparse path supports `distType = "Euclidean2D"` and
  `"Euclidean3D"` (inferred from a `z` coordinate column).
  `"Morphology-Aware"` distances require the dense path.
- `sigmaValues` **must** be supplied — because no distance matrix is
  built, there is nothing to derive a default bandwidth from.
- All other kernel options (`lowerLimit`, `upperQuantile`,
  `normalizeKernel`, `minAveCellNeighor`, distance normalization) behave
  exactly as in the dense path.

## 5. skrCCA is already exact and fast for two cell types

For the ordinary two-cell-type (or single within-type) problem,
[`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
now solves **all** requested canonical axes with one exact SVD of the
small PC-space cross-operator (`nPCA x nPCA`). There is no power
iteration to tune, so `maxIter`, `tol`, and `step_size` do not affect
the two-type result — the cost is dominated by the kernel, not the
optimization.

``` r

obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2)
obj <- computeNormalizedCorrelation(obj)

# Prefer regression gene weights for downstream use / transfer:
obj <- computeRegressionGeneScores(obj)
```

(`maxIter` / `step_size` still matter only when analyzing **three or
more** cell types jointly, or when supplying a transferred first axis,
where the sequential deflation path is used.)

## 6. Multi-slide at scale

[`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md)
shares one set of weight vectors across slides. Kernels are built **per
slide** and their small PC-space operators are summed — CoPro never
materializes a block-diagonal `N x N` kernel across the whole cohort.
Use `method = "auto"`/`"sparse"` exactly as above, and set `n_cores` to
parallelize the per-slide kernel and skrCCA work.

``` r

multi <- newCoProMulti(
  normalizedData = expr_list,   # named list of per-slide cells x genes matrices
  locationData   = loc_list,
  metaData       = meta_list,
  cellTypes      = ct_list
)
multi <- subsetData(multi, cellTypesOfInterest = c("Tumor", "Fibroblast"))
multi <- computePCA(multi, nPCA = 15, center = TRUE, scale. = TRUE)
multi <- computeKernelMatrix(multi, sigmaValues = c(0.05, 0.1, 0.2),
                             method = "auto", dropDistances = TRUE)
multi <- runSkrCCA(multi, scalePCs = TRUE, nCC = 2, n_cores = 4)
multi <- computeNormalizedCorrelation(multi)
multi <- computeRegressionGeneScores(multi)
```

## 7. A complete single-slide Xenium recipe

``` r

library(CoPro)

obj <- newCoProSingle(
  normalizedData = xenium_expr,   # dgCMatrix or BPCells IterableMatrix
  locationData   = xenium_loc,    # x, y columns
  metaData       = xenium_meta,
  cellTypes      = xenium_types
)
obj <- subsetData(obj, cellTypesOfInterest = c("Tumor", "Fibroblast"))
obj <- computePCA(obj, nPCA = 15, center = TRUE, scale. = TRUE)

# No computeDistance(): go straight to sparse kernels.
obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1, 0.2),
                           method = "auto", dropDistances = TRUE)

obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2)
obj <- computeNormalizedCorrelation(obj)
obj <- computeRegressionGeneScores(obj)

cell_scores <- getCellScores(obj, sigma = obj@sigmaValueChoice,
                             cellType = "Tumor")
```

## Quick reference: memory-saving knobs

| Setting | Where | Effect |
|----|----|----|
| sparse `dgCMatrix` / BPCells input | [`newCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/newCoProSingle.md) / [`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md) | Avoids a dense expression matrix; BPCells keeps it on disk |
| [`subsetData()`](https://zhen-miao.github.io/CoPro/reference/subsetData.md) before PCA | pipeline order | Shrinks every later step to the cell types of interest |
| `nPCA = 10–15` | [`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCA.md) | Fewer, stabler dimensions for targeted panels |
| `method = "auto"` / `"sparse"` | [`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md) | Never forms the dense `n x n` kernel; exact and lossless |
| **skip [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)** | pipeline order | The dense distance matrix is the main blow-up; the sparse path does not need it |
| `dropDistances = TRUE` | [`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md) | Clears `@distances` after kernels are built (default) |
| `autoThreshold` | [`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md) | Tunes the cell-count cutoff for the sparse path (default 5000) |
| `n_cores` | [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md) (multi-slide) | Parallelizes per-slide kernel and skrCCA work |

## Verifying equivalence (optional)

On a dataset small enough to run both, you can confirm the sparse and
dense paths agree. The stored kernels differ only in representation
(`dgCMatrix` vs dense) and by floating-point rounding:

``` r

# The dense path needs computeDistance() first; the sparse path does not.
obj_dense  <- computeKernelMatrix(computeDistance(obj), sigmaValues = 0.1,
                                  method = "dense", dropDistances = FALSE)
obj_sparse <- computeKernelMatrix(obj, sigmaValues = 0.1, method = "sparse")

k_dense  <- getKernelMatrix(obj_dense,  sigma = 0.1,
                            cellType1 = "Tumor", cellType2 = "Fibroblast")
k_sparse <- as.matrix(getKernelMatrix(obj_sparse, sigma = 0.1,
                            cellType1 = "Tumor", cellType2 = "Fibroblast"))
max(abs(k_dense - k_sparse))   # ~1e-15
```

## See also

- [`vignette("brain_merfish_two_type")`](https://zhen-miao.github.io/CoPro/articles/brain_merfish_two_type.md)
  — cross-cell-type co-progression on a MERFISH panel.
- [`vignette("organoid_one_type")`](https://zhen-miao.github.io/CoPro/articles/organoid_one_type.md)
  — within-cell-type spatial patterns.
- [`?computeKernelMatrix`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md),
  [`?computeSparseKernel`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md)
  — full kernel options.
