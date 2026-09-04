---
name: CoPro
description: "Guide users through CoPro R package workflows (>= 1.3.0): object creation, data-driven sigma selection, PCA- and gene-space CCA, multi-slide analysis, visualization, gene weights, frozen score transfer, permutation and replicate-level inference, large-data settings, and troubleshooting."
---

# /CoPro

Help users run the CoPro R package for detecting co-progression between cell types in spatial transcriptomics data. CoPro finds coordinated spatial gene expression patterns---either between different cell types or within a single cell type.

This guide tracks **CoPro >= 1.3.0**. Several defaults changed in 1.2.0 and 1.3.0 (summarized in "What changed recently" at the end). When unsure about an argument, check `args(fn)` / `?fn` in the installed package rather than older tutorials.

## API naming

- `camelCase` = the **object pipeline**: takes a `CoProSingle`/`CoProMulti`, returns one (or reads a result out of one). These chain: `newCoProSingle() |> subsetData() |> computePCA() |> ...`.
- `snake_case` = the **engine/utility layer** on plain matrices and data frames, no object involved: `optimize_bilinear()`, `optimize_sumcor_pca()`, `transfer_scores()`, `quantile_normalize()`, `resample_spatial()`, `copro_download_data()`.
- Documented exceptions: `fit_score_reference()` returns a fitted `CoProScoreReference`; score it with `predict()`.

## Choosing the right workflow

Ask the user about their data to determine which workflow fits:

| Scenario | Object | CCA call | Vignette |
|----------|--------|----------|----------|
| Single cell type, single slide | `newCoProSingle` | `runSkrCCA` | `organoid_one_type` |
| Two+ cell types, single slide | `newCoProSingle` | `runSkrCCA` | `brain_merfish_two_type`, `colon_d3_cross_type` |
| Two+ cell types, multiple slides --- PC space (default) | `newCoProMulti` | `computePCA` (per-slide centering) + `runSkrCCA` (equal-slide SUMCOR) | `large_datasets` section 6 |
| Two+ cell types, multiple slides --- gene space | `newCoProMulti` | `runGeneSpaceCCA` (or `runSkrCCA(space = "gene")`) | `colon_d0_multi_slide` |
| Train on one slide, score other slides | `newCoProSingle` (ref) | `runSkrCCA` + `fit_score_reference()`/`predict()` (or `getTransferCellScores()`) | `colon_d9_multi_slide` |
| Supervised with known ordering | `newCoProSingle` | `runSkrCCA(transferred_weight_1 = ...)` | `kidney_guided_gradient` |
| Spatial -> scRNA-seq full-transcriptome | matrices | `transfer_scores()` | `kidney_guided_gradient` |
| Hundreds of thousands of cells (Xenium, large MERFISH) | either | sparse/float32 kernels, sparse or BPCells input | `large_datasets` |

## Pre-run input checks

Before running the pipeline, validate the user's data. Flag problems early---running CoPro on malformed input produces silent garbage.

### 1. Required inputs

Verify all four inputs are present and correctly shaped:

```r
stopifnot(is.matrix(normalizedData) || inherits(normalizedData, "dgCMatrix") ||
          inherits(normalizedData, "IterableMatrix"))          # BPCells, on-disk
stopifnot(nrow(normalizedData) == nrow(locationData))
stopifnot(nrow(normalizedData) == length(cellTypes))
stopifnot(all(c("x", "y") %in% colnames(locationData)))
stopifnot(!anyNA(locationData[, c("x", "y")]), all(is.finite(as.matrix(locationData[, c("x", "y")]))))
```

- `normalizedData`: cells x genes (rows = cells, columns = genes). Dense `matrix`, sparse `Matrix::dgCMatrix` (recommended for imaging panels), or a BPCells `IterableMatrix` (recommended when it does not fit in RAM; `computePCA()` then uses `BPCells::svds()` out-of-core). Must have rownames (cell IDs) and colnames (gene names).
- `locationData`: data.frame with `x` and `y` columns (optionally `z`; no other column names). Rownames must match `normalizedData`. Non-finite coordinates are rejected with the offending cell IDs.
- `metaData`: data.frame with cell-level annotations. Rownames must match.
- `cellTypes`: character vector, same length as `nrow(normalizedData)`.
- Multi-slide: add `slideID`, a vector of the same length naming each cell's slide. `newCoProMulti()` takes the *combined* matrices for all slides plus this vector (not a list of per-slide matrices).

Already in a Seurat or SingleCellExperiment object? `asCoProSingle(x, spatialDim = "spatial", cellTypeCol = "celltype")` and `asCoProMulti(..., slideCol = "sample")` do the extraction in one call (Seurat: `layer = "data"`; SCE: `logcounts` assay by default). `CreateCoPro()` additionally splits slides larger than `maxCell` (default 50,000) into sub-slices automatically.

### 2. Cell counts

CoPro needs enough cells per cell type for stable CCA. Check counts after subsetting:

```r
table(cellTypes[cellTypes %in% cellTypesOfInterest])
```

- **Minimum**: ~200 cells per cell type for basic results
- **Recommended**: 1,000+ cells per cell type for stable gene weights and reproducible transfer
- **Warning sign**: if one cell type has <10% of the cells of another, the CCA may be dominated by the larger group

For multi-slide data, also check per-slide counts with `table(slideID, cellTypes)`. Under the multi-slide default (`objective = "sumcor"`) slides with fewer than `minCellsPerSlide` (default 10) cells of a type are *dropped* (and listed in `getCCAObjective(obj)$droppedSlides`); under `"sumcov"` they are only reported.

### 3. Normalization check

CoPro expects **normalized, log-transformed** expression. Raw counts will produce meaningless results. Check for signs of unnormalized data:

```r
range(normalizedData)                               # should NOT be e.g. 0-50000
all(normalizedData == floor(normalizedData))        # should be FALSE for log-transformed data
summary(Matrix::colMeans(normalizedData))           # large spread suggests no normalization
```

Recommended preprocessing (before CoPro):
- Library-size or DESeq2 size-factor normalization
- Log1p transformation (`log(x + 1)`)
- Optional: cap extreme values at the 95th percentile per gene (as in the organoid vignette)

Do NOT pass Seurat's `ScaleData` output (zero-centered, unit-variance)---CoPro centers and scales genes itself in `computePCA(center = TRUE, scale. = TRUE)`, with a guard that pins the scale to 1 for genes with sd < 1e-3 or < 1% nonzero cells (so near-constant genes cannot be amplified).

### 4. Cell type granularity

**Use broad cell types, not fine-grained subtypes.** CoPro discovers spatial gradients *within* each cell type group. If you pass subtypes (e.g., "PTS1", "PTS2", "PTS3" instead of "Tubular"), each subtype has too few cells, and the gradient CoPro would find is exactly the subtype identity---not a new biological axis.

Good: `c("Epithelial", "Fibroblast", "Immune")`
Bad: `c("Goblet", "Enterocyte_mature", "Enterocyte_progenitor", "Paneth", "TA", ...)`

The kidney vignette is instructive: tubular subtypes (PTS1--DCT-CNT) are used as *validation* of the CoPro axis, but the input cell type is just `"Tubular"`. CoPro discovers the corticomedullary gradient; the subtypes confirm it is biologically correct.

Exception: if two subtypes occupy completely different spatial niches with no co-localization (e.g., cortical vs medullary macrophages), grouping them may dilute the signal. Use domain knowledge.

### 5. Spatial coordinate sanity

```r
plot(locationData$x, locationData$y, pch = ".", asp = 1)
```

- Coordinates should look like tissue, not a random cloud
- Check for obvious batch shifts (e.g., two slides plotted on top of each other)
- For multi-slide: confirm `slideID` correctly partitions cells---`plot(x, y, col = as.factor(slideID))` should show non-overlapping groups. Kernels are always built *within* a slide; cells on different slides are never neighbors.
- **Units matter for sigma, not for CoPro.** `sigma` is a distance in whatever units `locationData` uses (microns, pixels, Visium spots). Never carry a sigma value between datasets; get the grid from `detectSigmaRange()` (next section). Since 1.2.0 distances are **not** rescaled by default (`normalizeDistance = FALSE`); the one rescaling worth recommending is `normalizeDistance = TRUE, normalizeMethod = "global"` (see *normalizeDistance* under Key parameters).

## Pipeline function order

### Standard single-slide pipeline (PCA-based)

```r
library(CoPro)

obj <- newCoProSingle(
  normalizedData = expr_matrix,    # cells x genes: matrix, dgCMatrix, or BPCells IterableMatrix
  locationData   = loc_df,         # data.frame with x, y (optional z) columns
  metaData       = meta_df,        # data.frame with cell annotations
  cellTypes      = cell_type_vec   # character vector, length = nrow(expr_matrix)
)
obj <- subsetData(obj, cellTypesOfInterest = c("TypeA", "TypeB"))
obj <- computePCA(obj, nPCA = 15)                 # 10-15 targeted panels, 20-30 full transcriptome

# sigma grid in the data's own coordinate units: the bandwidths at which the
# median cell is coupled to 5-20 effective neighbors
rng <- detectSigmaRange(obj)
rng                                               # prints per-block diagnostics + recommended grid

# No computeDistance() needed: method = "auto" (default) builds kernels straight
# from coordinates (sparse float32 when large, dense when small) and drops the
# distance slot afterward
obj <- computeKernelMatrix(obj, sigmaValues = rng$sigmaValues)

obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2)   # single slide: objective = "sumcov"
obj <- computeNormalizedCorrelation(obj)          # sets obj@sigmaValueChoice (argmax over the grid)
obj <- computeGeneAndCellScores(obj)              # PCA back-projection gene weights + cell scores
obj <- computeRegressionGeneScores(obj)           # regression gene weights (use for reporting)
```

`runSkrCCA()` defaults: `maxIter = 200`, `tol = 1e-5`, `step_size = 1`. For one or two cell types under `"sumcov"` the solution is an exact SVD/eigendecomposition, so those settings do not matter; they matter for three or more cell types, for a transferred first axis, and under `"sumcor"`.

**Global-normalized variant** (the only rescaling to recommend): same sparse route, but distances are first put on a unit-free scale, so the sigma grid can be written down without looking at the data and reused across datasets.

```r
obj <- computeKernelMatrix(obj, sigmaValues = c(0.01, 0.02, 0.05, 0.1),        # 1, 2, 5, 10 neighbor spacings
                           normalizeDistance = TRUE, normalizeMethod = "global")
```

Under `"global"` the median nearest-neighbor distance over all cells of interest becomes `normalizeTarget = 0.01`, so `sigma = 0.01 * k` reads as "k neighbor spacings". `detectSigmaRange()` still reports raw units: after the first normalized call it warns and prints the factor, and `rng$sigmaValues * obj@distanceScaleFactor` converts. Either choice---raw units (default) or `"global"`---is fine; pick one per object (the geometry is recorded, and a contradicting call errors).

**Legacy dense route** (what the colon and kidney vignettes still do): build a dense distance matrix, then kernels from it.

```r
obj <- computeDistance(obj, distType = "Euclidean2D", normalizeDistance = TRUE)   # normalizeMethod = "global"
obj <- computeKernelMatrix(obj, sigmaValues = c(0.01, 0.02, 0.05, 0.1))          # grid in global units
```

Only for small data (dense `n x n` per block) or to reproduce an existing analysis; the sparse route matches it to ~1e-6 in normalized correlation. `normalizeMethod = "percentile"` reproduces the pre-1.2.0 unit exactly and exists for reproduction only.

### Multi-slide pipeline, PC space (the 1.3.0 default)

```r
multi <- newCoProMulti(
  normalizedData = expr_matrix,    # all slides stacked (cells x genes)
  locationData   = loc_df,
  metaData       = meta_df,
  cellTypes      = cell_type_vec,
  slideID        = slide_id_vec    # which slide each cell is on
)
multi <- subsetData(multi, cellTypesOfInterest = cell_types)
multi <- computePCA(multi, nPCA = 15)                  # center_per_slide = TRUE: z-score genes within each
                                                       # (slide, cell type) block, then ONE shared PCA per type
rng   <- detectSigmaRange(multi)                       # reports every (pair, slide) block
multi <- computeKernelMatrix(multi, sigmaValues = rng$sigmaValues)
multi <- runSkrCCA(multi, scalePCs = TRUE, nCC = 2, n_cores = 4)   # default: objective = "sumcor", slideWeight = "equal"
multi <- computeNormalizedCorrelation(multi)           # calculationMode = "perSlide"
multi <- computeGeneAndCellScores(multi)
multi <- computeRegressionGeneScores(multi)
getCCAObjective(multi)                                 # space / objective / slideWeight / dropped slides actually used
```

One loading matrix and one set of CCA weights are shared across slides; kernels live within each slide. The two batch levers are `center_per_slide = TRUE` (removes per-slide *mean shifts* before PCA---the larger effect) and `objective = "sumcor"` (removes per-slide *scale* dominance). Pass `computePCA(center_per_slide = FALSE)` plus `runSkrCCA(objective = "sumcov")` to reproduce the legacy pooled multi-slide workflow. `slideWeight = "size"` reintroduces cell-count weighting under `"sumcor"`.

### Multi-slide pipeline, gene space

```r
sigma_choice <- rng$sigmaValues[3]                     # gene-space CCA takes ONE bandwidth
multi <- computeKernelMatrix(multi, sigmaValues = sigma_choice)
multi <- runGeneSpaceCCA(multi, sigma = sigma_choice, nCC = 1)     # no computePCA() needed
multi <- computeRegressionGeneScores(multi, sigma = sigma_choice)
# same thing through the main entry point:
# multi <- runSkrCCA(multi, space = "gene", sigmaChoice = sigma_choice, nCC = 1)
```

Gene-space CCA operates directly on expression: it z-standardizes per (slide, cell type), filters genes by prevalence (`min_prevalence = 0.008`, `min_cells = 20`), clips extremes (`clip = "quantile"`), and maximizes equal-slide SUMCOR (`objective = "sumcor"` default; `"sumcov"` available). It writes cell scores itself, so `getCellScoresInSitu()` works right after it. Defaults: `max_iter = 3000`, `tol = 1e-6`, `sweep = "gauss-seidel"` (`"jacobi"` reproduces pre-1.3.0 results), `step_size = 1`. It requires a `CoProMulti` object.

### Reference + transfer pipeline

Two routes and two weight types. **Default and first attempt: PCA back-projection weights.** `fit_score_reference()` freezes exactly those, and `getTransferCellScores()` uses them under `gene_score_type = "PCA"` (its default). Regression weights (`gene_score_type = "regression"`, what the D9 vignette shows) are the second attempt, tried only after the back-projection transfer has been checked. Between routes: the **frozen reference** is the default (internal leave-one-slide-out benchmark); `getTransferCellScores()` is for cross-platform transfer where quantile normalization against the reference is specifically wanted, or when the gene panels differ.

```r
# 1. Run the full standard pipeline on the reference slide, through computeGeneAndCellScores().

# 2a. Frozen reference: freezes per-type training mean/sd and the exact PCA
#     back-projected weights at obj@sigmaValueChoice. The target needs the same
#     normalization and gene panel, and only subsetData() -- no PCA/kernel/CCA.
ref <- fit_score_reference(ref_obj)                    # sigma = NULL -> ref_obj@sigmaValueChoice
tar_obj <- newCoProSingle(...)                         # or newCoProMulti(...) for several targets
tar_obj <- subsetData(tar_obj, cellTypesOfInterest = cell_types)
tar_scores <- predict(ref, tar_obj)                    # named list by cell type (cells x nCC)
tar_mat    <- predict(ref, tar_obj, aggregate = TRUE)  # or one matrix in target-cell order
# Optional: evaluate the transferred normalized correlation on the target (kernels only, no PCA)
tar_obj   <- computeKernelMatrix(tar_obj, sigmaValues = ref_obj@sigmaValueChoice)
tar_ncorr <- getTransferNormCorr(tar_obj, transfer_cell_scores = tar_scores,
                                 sigma_choice = ref_obj@sigmaValueChoice)

# 2b. Quantile-normalized transfer (needs PCA + kernels on the target so the
#     transferred normalized correlation can be evaluated there)
tar_obj <- computePCA(tar_obj, nPCA = 15)
tar_obj <- computeKernelMatrix(tar_obj, sigmaValues = sigma_opt)      # or sigma_choice_tar if target units differ
tar_scores <- getTransferCellScores(ref_obj = ref_obj, tar_obj = tar_obj,
                                    sigma_choice = sigma_opt,
                                    gene_score_type = "PCA")           # first attempt: back-projection (default)
tar_ncorr <- getTransferNormCorr(tar_obj = tar_obj, transfer_cell_scores = tar_scores,
                                 sigma_choice = sigma_opt)
# Second attempt, only if the transfer checks on the PCA transfer disappoint
# (needs computeRegressionGeneScores() on ref_obj):
tar_scores_reg <- getTransferCellScores(ref_obj = ref_obj, tar_obj = tar_obj,
                                        sigma_choice = sigma_opt,
                                        gene_score_type = "regression")
tar_ncorr_reg  <- getTransferNormCorr(tar_obj = tar_obj, transfer_cell_scores = tar_scores_reg,
                                      sigma_choice = sigma_opt)
```

A frozen reference collapses to one center and scale per cell type. On a `CoProMulti` fitted with the default `center_per_slide = TRUE`, transferred scores are target-invariant and mutually comparable but **not** on the same affine footing as `getCellScores()` on the fitted object (the function says so with a message and records `preprocessing` in the reference). Self-transfer reproduces `getCellScores()` exactly only under pooled preprocessing (`CoProSingle`, or `center_per_slide = FALSE`). Gene-space weights are rejected by `fit_score_reference()`.

### Supervised mode

When biological ordering is known for one cell type (e.g., nephron segment order):

```r
# Derive supervised weight via OLS on PCA scores (unit-normalized)
x_with_intercept <- cbind(1, tubular_pca_x)            # obj@pcaGlobal$Tubular$x
w1_raw <- lm.fit(x = x_with_intercept, y = as.numeric(ordered_labels))$coefficients[-1]
w1 <- w1_raw / sqrt(sum(w1_raw^2))

# Pass as transferred_weight_1 to runSkrCCA; later axes are conditioned on it
obj <- runSkrCCA(obj, scalePCs = TRUE, maxIter = 500, nCC = 4,
                 transferred_weight_1 = list(
                   TypeA = matrix(w1_typeA, ncol = 1),
                   TypeB = matrix(w1_typeB, ncol = 1)
                 ))
```

Works under both objectives (a 1.3.x fix: the one-slide SUMCOR shortcut used to silently discard the transferred axis). Not forwarded by `space = "gene"`---that is an error, not a silent drop.

## Key parameter guidance

### nPCA

- **scRNA-seq** (full transcriptome, 15k+ genes): `nPCA = 20-30`
- **seqFISH/MERFISH/Xenium** (targeted panel, 300-2000 genes): `nPCA = 10-15`
- The package default is 40; lower it for targeted panels. Too many PCs on a small panel add unstable noise dimensions that hurt gene-weight reproducibility and transfer. Cell scores and normalized correlations are unaffected by `nPCA`.

### sigma: three tools, used in order

1. **`detectSigmaRange(obj, minNeighbors = 5, maxNeighbors = 20, nSigma = 5)`** --- the grid. Inverts the effective-neighbor count (kernel row sum) so the grid spans 5-20 neighbors for the median cell, in the data's own units. Returns `sigmaValues` (log-spaced grid), `sigmaRange`, `feasible` (one range satisfies every block?), and `blocks` (per-block diagnostics). Deterministic (strided anchors), sub-second on 100k+ cells. If `feasible` is `FALSE`, cell types differ a lot in density---inspect `blocks` and widen the neighbor targets or reconsider the pairing.
2. **`computeNormalizedCorrelation(obj)`** --- the default pick. Sets `obj@sigmaValueChoice` to the argmax of normalized correlation over the grid. Fine for cross-type analyses.
3. **`selectSigmaByPermutation(obj, nPermu = 199, maxCells = 2000)`** --- the studentized pick. Compares the statistic at each sigma against *its own* toroidal-shift null and selects the max `z`; the p-value is already max-T adjusted for scanning the grid. Use it for **within-type (single cell type)** analyses, where the argmax rule is biased high. Returns `sigmaValueChoice`, `pValue`, `zMax`, `perSigma` (plot `z` vs `sigma`), `plateau`. Requires Euclidean geometry; warns if the choice lands on a grid edge (extend the grid). The organoid vignette uses this.

Sigma also drives memory: retained kernel pairs grow as `sigma^d`. `method = "auto"` predicts density and warns (with a suggested sigma) when a sparse kernel would save nothing---the fix is a smaller sigma, not more RAM.

Grids in the shipped vignettes are hand-picked and unit-specific: the colon and kidney vignettes (`c(0.005, 0.01, 0.02, 0.05, 0.1)`) use `normalizeDistance = TRUE`, i.e. global units since 1.2.0; the brain (`c(0.1, 0.14, 0.2, 0.5)`) and organoid (`c(0.01, ..., 0.2)`) vignettes use raw coordinate units. Neither kind transfers to another dataset as-is; use `detectSigmaRange()`, or the global unit.

### normalizeDistance

Two settings are recommended; everything else exists to reproduce old results.

- **`FALSE` (default since 1.2.0; was `TRUE`)** --- no rescaling. Sigma is in raw coordinate units and comes from `detectSigmaRange()`. Use this unless there is a reason not to.
- **`TRUE` with `normalizeMethod = "global"`** (the code default for `normalizeMethod`) --- one tissue-wide reference, the median nearest-neighbor distance over all cells of interest regardless of type, is mapped to `normalizeTarget = 0.01`. Cross-type and within-type steps derive the identical unit in either order, so it is safe on objects that mix `computeKernelMatrix()` and `computeSelfKernel()`. Choose it for a grid that reads in neighbor spacings (`c(0.01, 0.02, 0.05, 0.1)` = 1-10 spacings) and is portable between datasets, or to match the colon/kidney vignettes. Accepted by both the sparse route (`computeKernelMatrix(..., normalizeDistance = TRUE, normalizeMethod = "global")`) and the dense route (`computeDistance()`).

Do not suggest `"spacing"` (per-block reference: the unit then depends on cell-type abundance and colocalization, and cross-type vs within-type steps can disagree) or `"percentile"` (pre-1.2.0 rule: the densest block set the unit for the whole object). Both remain only so that older analyses can be reproduced. The geometry used is recorded on the object (`getDistanceGeometry(obj)`, `@distanceScaleFactor`); later kernel and sigma-detection calls inherit it, and contradicting it is an error.

### nCC

Number of canonical components. CC1 captures the strongest signal. Use `nCC = 2-4` to discover multiple orthogonal axes; test CC2+ with `runSkrCCAPermu_Conditional()` (see Inference).

### objective, slideWeight, space (1.3.0)

- `objective`: `"sumcov"` (single-slide default) maximizes the summed kernel covariance; `"sumcor"` (multi-slide and gene-space default) divides each slide's cross term by that slide's own score scales. On one slide they coincide for 1-2 cell types (any counts) or 3+ at equal counts, and CoPro short-circuits to the exact solvers there; otherwise an explicit `"sumcor"` runs the full-gradient optimizer.
- `slideWeight`: `"equal"` (default; equal nominal coefficients on the slide/pair terms -- not Kettenring SUMCOR, since the kernel-smoothed numerator over an unsmoothed denominator is not bounded by 1 and kernel scale still enters) or `"size"` (add the explicit `sqrt(n_i n_j)` factor); only with `"sumcor"`. Read `getCCADiagnostics()` after a PC-space SUMCOR fit for stopping status, projected-gradient residual, and Gram conditioning.
- `space`: `"pca"` (default) or `"gene"` (forwards to `runGeneSpaceCCA()`; pass the bandwidth as `sigmaChoice`).
- `step_size` in (0, 1] damps the iteration under both objectives (shorter step along the same geodesic). Reach for it only when a fit oscillates.
- Read back what was optimized with `getCCAObjective(obj)`.

### center_per_slide (computePCA, CoProMulti)

`TRUE` by default. Standardizes each (slide, cell type) block before the shared PCA, so slide-level location and scale cannot choose the retained subspace. This is the larger batch lever; `"sumcor"` alone does not remove a mean shift. Genes guarded on any slide (sd < 1e-3 or < 1% nonzero) are guarded on all.

## Visualization functions

Load the `copro-figure-style` skill before producing manuscript panels; the snippets below are the minimal correct calls.

### Core accessors

| Function | Purpose |
|----------|---------|
| `getCellScoresInSitu(obj, sigmaValueChoice, ccIndex = 1)` | Data frame with `x`, `y`, `cellScores`, `cellTypesSub` (and `slideID` for multi); rownames = cell IDs |
| `getCellScores(obj, sigma, cellType, slide = NULL, ccIndex = NULL)` | Cells x nCC score matrix for one type (and slide) |
| `getCorrTwoTypes(obj, cellTypeA, cellTypeB, ccIndex = 1, sigmaValueChoice = NULL)` | Cross-type correlation data (`AK` = kernel-smoothed A scores, `B`; plus `slideID` for multi) |
| `getCorrOneType(obj, cellTypeA, ccIndex = 1, sigmaValueChoice = NULL)` | Within-type spatial correlation (`AK`, `B`) |
| `getNormCorr(obj)` | One row per `sigmaValues` x `cellType1`/`cellType2` (`ct12`) x `CC_index` (and `slideID` for multi) with `normalizedCorrelation` |
| `getCCAObjective(obj)` | Which space/objective/slide weighting produced the weights |
| `getKernelMatrix(obj, sigma, cellType1, cellType2, slide = NULL)` | Inspect a kernel block (materialized as sparse) |

### Standard plot patterns

**Normalized correlation across sigma** (sigma selection diagnostic):
```r
ncorr <- getNormCorr(obj)
ggplot(ncorr, aes(x = sigmaValues, y = normalizedCorrelation)) +
  geom_point() + geom_line() + facet_wrap(~ ct12 + CC_index) + theme_minimal()
```

**In-situ cell scores** (viridis, quantile-clamped; option "D" for CC1, "C" for CC2):
```r
cs <- getCellScoresInSitu(obj, sigmaValueChoice = obj@sigmaValueChoice, ccIndex = 1)
lims <- quantile(cs$cellScores, c(0.025, 0.975), na.rm = TRUE)
ggplot(cs[order(cs$cellScores), ], aes(x = x, y = y, color = cellScores)) +
  geom_point(size = 0.3, shape = 16, stroke = 0) +
  scale_color_viridis_c(option = "D", limits = lims, oob = scales::squish) +
  coord_fixed() + theme_classic()
```

**Cross-type correlation scatter**:
```r
df <- getCorrTwoTypes(obj, cellTypeA = "TypeA", cellTypeB = "TypeB", ccIndex = 1,
                      sigmaValueChoice = obj@sigmaValueChoice)
ggplot(df, aes(x = AK, y = B)) +
  geom_point(size = 0.5, alpha = 0.4) +
  xlab("TypeA CC1 (spatially smoothed)") + ylab("TypeB CC1")
```

**Gene weight bar plot** (light = negative, dark = positive; Epithelial pink/darkred, Fibroblast lightblue/darkblue, Immune lightgreen/darkgreen):
```r
key <- paste0("geneScores|sigma", obj@sigmaValueChoice, "|TypeA")
gs <- obj@geneScoresRegression[[key]][, 1]
top <- head(sort(abs(gs), decreasing = TRUE), 20)
top_df <- data.frame(
  gene = factor(names(top), levels = rev(names(top))),
  weight = gs[names(top)],
  direction = ifelse(gs[names(top)] > 0, "positive", "negative")
)
ggplot(top_df, aes(x = gene, y = weight, fill = direction)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("positive" = "darkred", "negative" = "lightpink")) +
  theme_classic() + theme(legend.position = "none", axis.title.y = element_blank())
```

**Multi-slide panels** (use patchwork, NOT facet_wrap, to preserve per-slide coordinates and a shared color limit):
```r
library(patchwork)
cs <- getCellScoresInSitu(multi, sigmaValueChoice = sigma_choice, ccIndex = 1)   # has a slideID column
lims <- quantile(cs$cellScores, c(0.025, 0.975), na.rm = TRUE)
plots <- lapply(getSlideList(multi), function(sl) {                # one ggplot per slide, shared lims
  d <- cs[cs$slideID == sl, ]
  ggplot(d[order(d$cellScores), ], aes(x, y, color = cellScores)) + geom_point(size = 0.2, shape = 16, stroke = 0) +
    scale_color_viridis_c(limits = lims, oob = scales::squish) + coord_fixed() + ggtitle(sl) + theme_classic()
})
wrap_plots(plots, ncol = 3) + plot_layout(guides = "collect")
```

## Gene weights: PCA vs regression

CoPro provides two gene weight methods:

- `@geneScores` via `computeGeneAndCellScores()` --- PCA back-projection. This **is** the scoring functional (cell score = standardized expression x these weights), so it is what `fit_score_reference()` freezes for transfer.
- `@geneScoresRegression` via `computeRegressionGeneScores()` --- per-gene regression of expression on the cell score (marginal association).

**Report regression weights** (figures, tables, gene lists): they avoid PCA collinearity splitting weight between correlated genes, are insensitive to `nPCA`, and reproduce better across replicates. **Transfer with back-projection weights first** (the default of both transfer routes); regression weights are the second thing to try (`getTransferCellScores(gene_score_type = "regression")`, or `gs_ct` from `@geneScoresRegression` for scRNA-seq), judged on the same transfer checks. Regression weights are not the scoring functional, so even self-transfer with them is only approximately faithful, whereas back-projection self-transfer is exact.

Both slots use the same key format: `"geneScores|sigma{X}|{CellType}"`, each a genes x nCC matrix. `computeRegressionGeneScores(sigma = NULL)` uses `obj@sigmaValueChoice`.

## Score transfer to scRNA-seq

Transfer spatial gene weights to scRNA-seq for full-transcriptome analysis:

```r
key <- paste0("geneScores|sigma", obj@sigmaValueChoice, "|TypeA")
gene_weights <- obj@geneScores[[key]][, 1, drop = FALSE]             # first attempt: back-projection
# gene_weights <- obj@geneScoresRegression[[key]][, 1, drop = FALSE] # second attempt (the kidney vignette uses this)
transferred <- transfer_scores(
  mat_A = spatial_expr,      # reference spatial data (cells x shared genes)
  mat_B = scrna_expr,        # target scRNA-seq (cells x shared genes)
  gs_ct = gene_weights,      # gene weight matrix (shared genes x 1)
  use_quantile_normalization = TRUE,
  gs_weight_threshold = 0,
  verbose = FALSE
)
```

Subset all three to the shared genes first (`intersect()` of colnames and `rownames(gs_ct)`). After transfer, run full-transcriptome regression:
```r
for (gene in colnames(scrna_full)) {
  ct <- cor.test(transferred, scrna_full[, gene])
  # collect estimate, p-value; adjust with p.adjust(method = "BH")
}
```

## Inference: permutation and replicate-level tests

### Single slide: cell-level permutation

```r
obj <- runSkrCCAPermu(obj, nPermu = 999, permu_method = "bin")   # null at the fitted sigma
obj <- computeNormalizedCorrelationPermu(obj)
calculate_pvalue(obj, cc_index = 1)                               # alternative = "greater" only
```

- `nPermu` defaults to 999 (Monte-Carlo floor 0.001); `nPermu < 2` is rejected. Use `nPermu = 5-20` only for smoke tests.
- `calculate_pvalue()` warns that a fixed-sigma p-value at the data-selected `obj@sigmaValueChoice` is not adjusted for sigma selection; either pre-declare `sigma =` in `runSkrCCAPermu()` or use `runSkrCCAPermu_FairSigma()`.
- The default `"bin"` null sizes its grid from the recorded distance scale; if it warns "could not recover distance scale", pass `num_bins_x`/`num_bins_y` explicitly (e.g. 10 x 10).
- The null is built under whichever criterion the weights were fitted with (`sumcov` where the SUMCOR reduction holds, a re-optimized `sumcor` null otherwise). No explicit refit needed.
- `alternative = "two.sided"` is an error for `calculate_pvalue()`: the max-aggregated statistic's null is not symmetric.
- **Within-type** (one cell type) permutation works as of 1.3.x: the single type is permuted against its own locations.
- Multiple axes: `runSkrCCAPermu_Conditional(obj, nPermu = 999)` then `calculate_pvalue_stepdown(obj)` gives step-down adjusted p-values per `CC_index`, conditioning each axis on the observed lower axes.
- Max-over-sigma (when sigma was picked by argmax): `runSkrCCAPermu_FairSigma(obj, nPermu = 999)` re-optimizes at every sigma per draw. Both `_FairSigma` and `_Conditional` support at most two cell types.
- Memory/speed knobs are now arguments, not options: `factorize = TRUE` (reuses one kernel product across draws; exact), `compactPermutation = TRUE` (stores one seed per draw; changes which permutations are drawn, so leave off to reproduce saved results), `n_cores`.
- Permutation entry points restore the caller's RNG state.

### Multiple slides: replicate-level inference

Cell-level permutation **refuses** `CoProMulti` objects (cells are not biological replicates). Use

```r
inf <- runSlideLevelInference(multi, cc_index = 1, replicate_id = donor_by_slide,
                              n_resamples = 9999, seed = 1)      # n_resamples must be >= 99
inf    # equal-replicate effect, bootstrap CI, sign-flip p-value, per-replicate held-out results
```

Run it after `computePCA()` and `computeKernelMatrix()` (no `runSkrCCA()` needed: weights and sigma are learned leaving out each replicate and evaluated on the held-out one). `replicate_id`, named by `getSlideList(multi)`, groups slides from the same donor. Currently supports exactly two cell types; with very few replicates the sign-flip p-value floor is coarse (`2^-k`).

## Large datasets (Xenium, large MERFISH)

See `vignette("large_datasets")`. The short version:

- Store expression as `dgCMatrix` or a BPCells `IterableMatrix`; `subsetData()` early; `nPCA = 10-15`.
- **Do not call `computeDistance()`**. `computeKernelMatrix(method = "auto")` builds fixed-radius sparse kernels from coordinates (float32 for large blocks, exact to ~1e-6 in normalized correlation; `method = "sparse"` keeps float64 for exactness checks; `"dense"` is the classic route and needs distances).
- `dropDistances = TRUE` (default) clears `@distances`; `autoThreshold` (5000 cells) switches to sparse; `denseThreshold` (0.3) is the predicted density above which `"auto"` warns.
- When memory is tight, loop over one sigma at a time: `computeKernelMatrix()` holds every requested sigma simultaneously.
- Threads for the float32 operators: `nThreads` on `computeSparseKernelFloat32()`; `computeKernelMatrix(method = "float32")` inherits `options(CoPro.float32Threads)`; HPC scheduler variables (`SLURM_CPUS_PER_TASK` etc.) are honored automatically.
- `n_cores` in `runSkrCCA()` parallelizes per-slide work for `CoProMulti`.
- `materializeFloat32Kernels(obj)` converts float32 kernels to ordinary sparse matrices for legacy code that needs them.

## Post-run biological sanity checks

After running the pipeline, verify that results reflect biology rather than technical artifacts. Do not skip these---CoPro will always produce numbers, even on garbage input.

### 1. Normalized correlation magnitude

```r
ncorr <- getNormCorr(obj)
print(ncorr[ncorr$CC_index == 1, ])
obj@sigmaValueChoice
```

- **Healthy signal**: normalized correlation 0.05--0.3 for CC1 (dataset-dependent)
- **Suspiciously high** (>0.5): likely a batch effect, library-size artifact, or cell type identity signal masquerading as spatial co-progression. Inspect the in-situ plot---does the pattern align with known biology, or does it just separate two clusters?
- **Near zero across all sigma**: the cell types may not spatially co-vary at the scales tested. Check `detectSigmaRange(obj)$blocks` (are the types actually co-localized at 5-20 neighbors?) before widening the grid.
- **Choice on the grid edge**: the scan ran out of grid. Extend it, or use `selectSigmaByPermutation()`.

### 2. Spatial pattern check

Plot cell scores in situ and ask: does this look like a biologically meaningful spatial gradient, or a technical artifact?

**Biological signals** look like:
- Smooth spatial gradients (crypt-to-villus, cortex-to-medulla)
- Region-specific enrichment matching known tissue architecture
- Consistent patterns across cell types in the same spatial region

**Artifact red flags**:
- Scores correlate with cell density rather than tissue structure (library-size artifact)
- Sharp boundaries that match slide edges or imaging FOV boundaries (batch effect)
- Identical scores across all cell types in a region (the CCA found cell-type identity, not co-progression)
- Random salt-and-pepper pattern with no spatial coherence (noise, or sigma too small)

### 3. Gene weight sanity

Check whether the top genes make biological sense:

```r
key <- paste0("geneScores|sigma", obj@sigmaValueChoice, "|TypeA")
gs <- obj@geneScoresRegression[[key]][, 1]
head(sort(gs, decreasing = TRUE), 10)     # top positive
head(sort(gs), 10)                         # top negative
```

- **Good sign**: top genes include known markers for the expected spatial gradient (e.g., pericentral/periportal markers in liver, proximal/distal markers in kidney)
- **Bad sign**: top genes are all ribosomal (`Rpl*`, `Rps*`), mitochondrial (`mt-*`), or housekeeping genes---these usually indicate a technical quality axis, not biology
- **Bad sign**: top genes are all sex-linked (`Xist`, `Ddx3y`)---CoPro found a sex batch effect
- Cross-reference a few top genes with the literature or gene ontology to confirm the axis has a coherent biological interpretation

### 4. Cross-type correlation consistency

For cross-cell-type analyses, verify the correlation between cell types:

```r
df <- getCorrTwoTypes(obj, cellTypeA = "TypeA", cellTypeB = "TypeB", ccIndex = 1,
                      sigmaValueChoice = obj@sigmaValueChoice)
cor(df$AK, df$B)
if ("slideID" %in% names(df)) tapply(seq_len(nrow(df)), df$slideID, function(i) cor(df$AK[i], df$B[i]))
```

- The scatter should show a clear positive trend (CoPro maximizes this)
- If the correlation is driven by a handful of outlier cells, the result may not be robust
- For multi-slide data, check that the correlation holds per slide, not just in aggregate (aggregate correlation can be inflated by batch-level mean shifts; with `center_per_slide = TRUE` this is mostly removed, but still look)

### 5. Statistical validation

Single slide: `runSkrCCAPermu()` + `computeNormalizedCorrelationPermu()` + `calculate_pvalue()`. Multi-slide: `runSlideLevelInference()`. See the Inference section. If the observed normalized correlation is not clearly above the null, the signal may be spurious.

### 6. Multi-slide consistency

For datasets with biological replicates, the strongest validation is cross-slide consistency:
- Do the same genes rank highly across slides (fit each slide alone with `newCoProSingle()` and compare regression weights)?
- Do frozen-reference scores (`predict()`) reproduce the spatial pattern on held-out slides?
- Are per-slide normalized correlations comparable (`getNormCorr(multi)` has one row per slide under `calculationMode = "perSlide"`)?

Inconsistent results across replicates suggest overfitting to slide-specific technical variation.

## Score transfer verification

After transferring scores to other slides or scRNA-seq, verify the transfer worked. A bad transfer is worse than no transfer---it gives false confidence in a meaningless axis.

### Transfer to other spatial slides

**Checks**:

1. **Transferred normalized correlation**: should be comparable to the reference slide's value. A large drop (e.g., reference 0.15 vs transferred 0.03) means the pattern does not generalize. `predict()` and `getTransferCellScores()` return the same shape (named list by cell type of cells x nCC matrices), so this works for both routes; the target only needs kernels at `sigma_choice`.

```r
tar_obj   <- computeKernelMatrix(tar_obj, sigmaValues = sigma_opt)        # if not built yet
tar_ncorr <- getTransferNormCorr(tar_obj = tar_obj, transfer_cell_scores = tar_scores,
                                 sigma_choice = sigma_opt)                # list keyed "sigma_<value>"
```

2. **Visual comparison**: plot transferred scores in situ on the target slide. The spatial pattern should resemble the reference---same gradient direction, similar regions highlighted.

3. **Gene overlap**: the reference and target must share the modeled gene panel. `fit_score_reference()`/`predict()` require the same genes; `getTransferCellScores()` intersects and reports "retaining N genes". If N is small relative to the panel, the transfer is underpowered.

4. **Self-transfer sanity**: `predict(ref, ref_obj)` must reproduce `getCellScores(ref_obj, ...)` exactly on a `CoProSingle`. If it does not, the target preprocessing differs from the reference's.

5. **Weight type**: the back-projection transfer comes first. If checks 1-2 disappoint, repeat with `gene_score_type = "regression"` (or `@geneScoresRegression` as `gs_ct`) and compare the two on transferred normalized correlation and in-situ pattern; keep whichever passes and state which was used.

### Transfer to scRNA-seq

```r
transferred <- transfer_scores(mat_A = spatial_expr, mat_B = scrna_expr, gs_ct = gene_weights,
                               use_quantile_normalization = TRUE, gs_weight_threshold = 0)
```

**Checks**:

1. **Known marker validation**: if the spatial axis has a known biological interpretation (e.g., corticomedullary), check that known markers separate correctly in the transferred scores. In the kidney vignette, transferred nephron scores increase monotonically from PTS1 to DCT-CNT.

```r
boxplot(transferred_score ~ subtype, data = validation_df)
```

2. **Score distribution**: transferred scores should have a unimodal, roughly continuous distribution. A bimodal distribution with a gap suggests the transfer is capturing cell-type identity rather than a gradient.

3. **Full-transcriptome regression**: after transferring, regress all scRNA-seq genes against the transferred score. Significant genes (FDR < 0.05) should include biologically relevant genes beyond the panel genes that defined the axis.

```r
sum(reg_results$fdr < 0.05)    # expect thousands if the axis is real; only panel genes => circular
```

4. **UMAP/t-SNE sanity**: color the scRNA-seq embedding by transferred score. The gradient should align with known biology (e.g., differentiation trajectory), not with technical variables (batch, library size, percent mitochondrial).

### Common transfer pitfalls

- **Too few shared genes**: if <100 genes overlap between reference and target, the transfer is underpowered. Consider a broader spatial panel or a different reference.
- **Batch effects in scRNA-seq**: quantile normalization helps, but large batch effects can still distort transferred scores. Check whether scores correlate with batch labels.
- **Mixing routes**: frozen-reference scores and quantile-normalized scores are on different scales; do not pool them on one axis. Likewise, frozen scores from a `center_per_slide = TRUE` multi-slide fit are not on the fitted object's own score scale.
- **Comparing weight types**: back-projection and regression scores sit on different scales (regression weights have larger magnitude). Compare the two transfers by transferred normalized correlation and in-situ pattern, never by score values, and do not pool them on one axis.

## Available example datasets

Download via `copro_download_data()` (cached under `tools::R_user_dir("CoPro", "cache")`; 5--30 MB each):

| ID | Description |
|----|-------------|
| `"organoid"` | 72 hr intestinal organoid, single cell type (Epithelial), 140 genes --- within-type |
| `"brain_merfish"` | Brain MERFISH, D1/D2 GABAergic striatal neurons --- two-type |
| `"colon_d3"` | Single colon DSS day-3 slide, Epithelial/Fibroblast/Immune + MU labels, 891 genes --- cross-type, `nCC = 4` |
| `"colon_d0_multi"` | 3 healthy colon slides, 3 cell types, 928 genes --- multi-slide gene-space CCA |
| `"colon_d3_multi"` | 3 colon day-3 slides --- multi-slide `newCoProMulti` (no vignette currently) |
| `"colon_d9"` | 3 colon day-9 slides (inflammation) + MU labels --- reference + transfer |
| `"kidney"` | Kidney seqFISH + scRNA-seq, Tubular/Vascular, 1298 genes --- supervised + scRNA-seq transfer |

A tiny bundled object for smoke tests: `readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))` (200 cells, 80 genes, two types).

## Common issues

- **"No kernel for sigma ..."** / kernel lookup errors: every sigma passed to `runSkrCCA(sigmaChoice=)`, `runGeneSpaceCCA(sigma=)`, `getCellScores(sigma=)`, etc. must be one of the values `computeKernelMatrix()` was given. Use `rng$sigmaValues` consistently, and note `obj@sigmaValueChoice` is a single numeric.
- **`method = "auto"` warns "predicted density ..."**: sigma is too large for a sparse kernel at this cell density. Use the smaller sigma it suggests; do not add RAM.
- **Geometry contradiction error** (`distType`/`xDistScale` disagree with the record): kernel steps inherit what `computeDistance()` recorded. Pass the same values or leave them `NULL`.
- **Results differ from an older analysis**: (i) `normalizeDistance` now defaults to `FALSE`; (ii) `CoProMulti` now defaults to `center_per_slide = TRUE` + `objective = "sumcor"`; (iii) gene-space `sweep` defaults to `"gauss-seidel"`. Pass the legacy values explicitly, and only to reproduce (`normalizeDistance = TRUE, normalizeMethod = "percentile"`, `center_per_slide = FALSE, objective = "sumcov"`, `sweep = "jacobi"`). For a new analysis stay with the default (no rescaling) or `normalizeMethod = "global"`.
- **Permutation on `CoProMulti` errors**: intended. Use `runSlideLevelInference()`.
- **`calculate_pvalue(alternative = "two.sided")` errors**: intended; the statistic's null is one-sided.
- **Convergence notices flood the console**: they are `message()`s now; wrap in `suppressMessages()` or capture with `capture.output(type = "message")`.
- **Fit oscillates / does not converge**: set `step_size = 0.5` (both `runSkrCCA()` and `runGeneSpaceCCA()`), expect more iterations; raise `maxIter`/`max_iter`.
- **`fit_score_reference()` errors "Gene scores are unavailable"**: run `computeGeneAndCellScores()` first; it never uses regression weights and rejects gene-space fits.
- **`predict()` errors "Target is missing frozen-reference genes for ..."**: the frozen route needs every modeled gene present in the target, under the same normalization. Use `getTransferCellScores()` (which intersects panels and quantile-normalizes) when panels differ or cross-platform transfer is intended.
- **BPCells input fails**: needs the BPCells package (`remotes::install_github("bnprks/BPCells/r")`); all `computePCA()` branches stream it as of the current release.
- **Multi-slide `facet_wrap` distortion**: use `patchwork::wrap_plots` with per-slide panels instead.

## What changed recently (for reconciling older scripts)

- **1.2.0**: `detectSigmaRange()`; `normalizeDistance` default `FALSE` with `normalizeMethod = "global"/"spacing"/"percentile"` (recommend only the default or `"global"`); `method = "auto"` uses float32 sparse kernels and predicts density; `@distanceGeometry` record; `computeSelfKernel(normalizeDistance = "inherit")`.
- **1.3.0**: `runSkrCCA(objective=, slideWeight=, space=, minCellsPerSlide=)`; `CoProMulti` defaults `center_per_slide = TRUE` + `"sumcor"`/`"equal"`; `getCCAObjective()`; `runGeneSpaceCCA(sweep=, step_size=, objective=)`; `selectSigmaByPermutation()`; within-type permutation fixed; SUMCOR null refits per draw; `transferred_weight_1` survives `"sumcor"`; permutation null matches the fitted criterion.
- **Development (post-1.3.0)**: `fit_score_reference()` / `predict()` frozen transfer; `options(CoPro.*)` became arguments (`factorize`, `compactPermutation`, `nThreads`); one S4 method per accessor on the `CoPro` base class; convergence notices are messages; `computePCA()` degeneracy guard on every branch and full BPCells support; `computeNormalizedCorrelation(normalizer=)` exposes the denominator (`"legacy"` default).
