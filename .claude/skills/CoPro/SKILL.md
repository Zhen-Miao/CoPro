---
name: CoPro
description: "Guide users through CoPro R package workflows: object creation, pipeline execution, visualization, gene weights, score transfer, and troubleshooting."
---

# /CoPro

Help users run the CoPro R package for detecting co-progression between cell types in spatial transcriptomics data. CoPro finds coordinated spatial gene expression patterns---either between different cell types or within a single cell type.

## Choosing the right workflow

Ask the user about their data to determine which workflow fits:

| Scenario | Object | CCA method | Vignette |
|----------|--------|------------|----------|
| Single cell type, single slide | `newCoProSingle` | `runSkrCCA` | `organoid_one_type` |
| Two+ cell types, single slide | `newCoProSingle` | `runSkrCCA` | `brain_merfish_two_type`, `colon_d3_cross_type` |
| Two+ cell types, multiple slides (joint) | `newCoProMulti` | `runGeneSpaceCCA` | `colon_d0_multi_slide` |
| Train on one slide, transfer to others | `newCoProSingle` (ref) | `runSkrCCA` + `getTransferCellScores` | `colon_d9_multi_slide` |
| Supervised with known ordering | `newCoProSingle` | `runSkrCCA` with `transferred_weight_1` | `kidney_guided_gradient` |

## Pipeline function order

### Standard single-slide pipeline (PCA-based)

```r
library(CoPro)

obj <- newCoProSingle(
  normalizedData = expr_matrix,    # cells x genes, normalized + log-transformed
  locationData   = loc_df,         # data.frame with x, y columns
  metaData       = meta_df,        # data.frame with cell annotations
  cellTypes      = cell_type_vec   # character vector, length = nrow(expr_matrix)
)
obj <- subsetData(obj, cellTypesOfInterest = c("TypeA", "TypeB"))
obj <- computePCA(obj, nPCA = 30, center = TRUE, scale. = TRUE)
obj <- computeDistance(obj, distType = "Euclidean2D")
obj <- computeKernelMatrix(obj, sigmaValues = c(0.01, 0.05, 0.1))
obj <- runSkrCCA(obj, scalePCs = TRUE, maxIter = 500, nCC = 2)
obj <- computeNormalizedCorrelation(obj)
obj <- computeGeneAndCellScores(obj)
obj <- computeRegressionGeneScores(obj)
```

### Multi-slide joint pipeline (gene-space CCA)

```r
multi_obj <- newCoProMulti(
  normalizedData = expr_matrix,
  locationData   = loc_df,
  metaData       = meta_df,
  cellTypes      = cell_type_vec,
  slideID        = slide_id_vec    # which cells belong to which slide
)
multi_obj <- subsetData(multi_obj, cellTypesOfInterest = cell_types)
multi_obj <- computeDistance(multi_obj, distType = "Euclidean2D")
multi_obj <- computeKernelMatrix(multi_obj, sigmaValues = 0.01)
multi_obj <- runGeneSpaceCCA(multi_obj, sigma = 0.01, nCC = 1)
multi_obj <- computeRegressionGeneScores(multi_obj, sigma = 0.01)
```

Gene-space CCA does NOT require `computePCA`. It operates directly on gene expression, handles per-slide standardization internally, and maximizes the average per-slide canonical correlation.

### Reference + transfer pipeline

```r
# 1. Run full pipeline on reference slide (see standard pipeline above)
# 2. Create target objects with PCA, distance, kernel (no CCA needed)
tar_obj <- newCoProSingle(...)
tar_obj <- subsetData(tar_obj, cellTypesOfInterest = cell_types)
tar_obj <- computePCA(tar_obj, nPCA = 40, center = TRUE, scale. = TRUE)
tar_obj <- computeDistance(tar_obj, distType = "Euclidean2D")
tar_obj <- computeKernelMatrix(tar_obj, sigmaValues = sigma_choice)

# 3. Transfer gene weights
tar_scores <- getTransferCellScores(
  ref_obj = ref_obj, tar_obj = tar_obj,
  sigma_choice = 0.01, gene_score_type = "regression"
)

# 4. Assess transfer quality
tar_ncorr <- getTransferNormCorr(
  tar_obj = tar_obj,
  transfer_cell_scores = tar_scores,
  sigma_choice = 0.01
)
```

### Supervised mode

When biological ordering is known for one cell type (e.g., nephron segment order):

```r
# Derive supervised weight via OLS on PCA scores
x_with_intercept <- cbind(1, tubular_pca_x)
w1_raw <- lm.fit(x = x_with_intercept, y = as.numeric(ordered_labels))$coefficients
w1_raw <- w1_raw[-1]
w1 <- w1_raw / sqrt(sum(w1_raw^2))

# Pass as transferred_weight_1 to runSkrCCA
obj <- runSkrCCA(obj, scalePCs = TRUE, maxIter = 500, nCC = 4,
                 transferred_weight_1 = list(
                   TypeA = matrix(w1_typeA, ncol = 1),
                   TypeB = matrix(w1_typeB, ncol = 1)
                 ))
```

## Key parameter guidance

### nPCA

- **scRNA-seq** (full transcriptome, 15k+ genes): `nPCA = 20-30`
- **seqFISH/MERFISH** (targeted panel, 500-2000 genes): `nPCA = 10-15`
- Too many PCs on small panels adds unstable noise that hurts gene weight reproducibility

### sigmaValues

Controls the spatial neighborhood bandwidth. Smaller = local patterns, larger = broad patterns. Test multiple values; CoPro selects the optimum via normalized correlation.

Typical ranges:
- Organoids/small tissues: `c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2)`
- Standard spatial transcriptomics: `c(0.005, 0.01, 0.02, 0.05, 0.1)`
- Brain MERFISH (larger spatial scale): `c(0.1, 0.14, 0.2, 0.5)`

### normalizeDistance

`TRUE` by default. Normalizes distances so the 0.01 percentile = 0.01, making kernel behavior consistent across datasets. Set `FALSE` for datasets where absolute distance matters (e.g., organoids).

### nCC

Number of canonical components. CC1 captures the strongest signal. Use `nCC = 2-4` to discover multiple orthogonal axes.

## Visualization functions

### Core accessors

| Function | Purpose |
|----------|---------|
| `getCellScoresInSitu(obj, sigmaValueChoice, ccIndex)` | Cell scores with x,y for spatial plots |
| `getCorrTwoTypes(obj, sigmaValueChoice, cellTypeA, cellTypeB, ccIndex)` | Cross-type correlation (columns: AK, B, slideID) |
| `getCorrOneType(obj, sigmaValueChoice, cellTypeA, ccIndex)` | Within-type spatial correlation (columns: AK, B) |
| `getNormCorr(obj)` | Normalized correlation across sigma values |
| `getCellScores(obj, sigma, cellType)` | Raw cell score vector |

### Standard plot patterns

**In-situ cell scores** (use viridis, quantile-clamped color scale):
```r
cs <- getCellScoresInSitu(obj, sigmaValueChoice = sigma_opt, ccIndex = 1)
lims <- quantile(cs$cellScores, c(0.025, 0.975), na.rm = TRUE)
ggplot(cs, aes(x = x, y = y, color = cellScores)) +
  geom_point(size = 0.5, stroke = 0) +
  scale_color_viridis_c(limits = lims, oob = scales::squish) +
  coord_fixed() + theme_classic()
```

**Cross-type correlation scatter**:
```r
df <- getCorrTwoTypes(obj, sigmaValueChoice = sigma_opt,
                      cellTypeA = "TypeA", cellTypeB = "TypeB", ccIndex = 1)
ggplot(df, aes(x = AK, y = B)) +
  geom_point(size = 0.5, alpha = 0.4) +
  xlab("TypeA CC1 (spatially smoothed)") + ylab("TypeB CC1")
```

**Gene weight bar plot** (manuscript style):
```r
key <- paste0("geneScores|sigma", sigma_opt, "|CellType")
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

**Multi-slide panels** (use patchwork, NOT facet_wrap, to preserve per-slide coordinates):
```r
library(patchwork)
plots <- lapply(slides, function(sl) { ... })
wrap_plots(plots, ncol = 3) + plot_layout(guides = "collect")
```

## Gene weights: PCA vs regression

CoPro provides two gene weight methods:

- `@geneScores` via `computeGeneAndCellScores()` -- PCA back-projection
- `@geneScoresRegression` via `computeRegressionGeneScores()` -- regression-based

**Prefer regression gene weights.** They avoid collinearity (correlated genes splitting weight), are independent of nPCA choice, and show better cross-replicate reproducibility.

Both slots use the same key format: `"geneScores|sigma{X}|{CellType}"`, each containing a genes x nCC matrix.

## Score transfer to scRNA-seq

Transfer spatial gene weights to scRNA-seq for full-transcriptome analysis:

```r
transferred <- transfer_scores(
  mat_A = seqfish_expr,      # reference spatial data (cells x genes)
  mat_B = scrna_expr,        # target scRNA-seq (cells x genes)
  gs_ct = gene_weights,      # gene weight matrix (genes x 1)
  use_quantile_normalization = TRUE,
  gs_weight_threshold = 0,
  verbose = FALSE
)
```

After transfer, run full-transcriptome regression:
```r
for (gene in colnames(scrna_full)) {
  ct <- cor.test(transferred_scores, scrna_full[, gene])
  # collect beta, r, p-value
}
```

## Permutation testing

```r
obj <- runSkrCCAPermu(obj, nPermu = 100, permu_method = "bin",
                       num_bins_x = 10, num_bins_y = 10)
obj <- computeNormalizedCorrelationPermu(obj)
```

Use `nPermu >= 100` for reliable p-values. `nPermu = 5` is only for quick checks.

## Available example datasets

Download via `copro_download_data()`:

| ID | Description |
|----|-------------|
| `"colon_d0_multi"` | 3 healthy colon slides, 3 cell types, 928 genes |
| `"colon_d3"` | Single colon D3 slide, 3 cell types + MU labels, 891 genes |
| `"colon_d9"` | 3 colon D9 slides (inflammation), 3 cell types + MU labels |
| `"brain_merfish"` | Brain MERFISH, D1/D2 GABAergic neurons in striatum |
| `"kidney"` | Kidney seqFISH + scRNA-seq, tubular/vascular, 1298 genes |
| `"organoid"` | 72hr intestinal organoid, single cell type, 140 genes |

## Common issues

- **Kernel too small**: if `computeKernelMatrix` warns about zeros, increase sigma values or lower `lowerLimit`
- **Gene weights all near zero**: check that `computeGeneAndCellScores()` was called AFTER a recent CoPro version (the `1/sdev` fix from 2026-03-29)
- **Transfer scores look wrong**: use `gs_weight_threshold = 0` (default changed from 0.005)
- **Multi-slide facet_wrap distortion**: use `patchwork::wrap_plots` with per-slide panels instead
