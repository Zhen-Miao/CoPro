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

## Pre-run input checks

Before running the pipeline, validate the user's data. Flag problems early---running CoPro on malformed input produces silent garbage.

### 1. Required inputs

Verify all four inputs are present and correctly shaped:

```r
stopifnot(is.matrix(normalizedData) || is(normalizedData, "dgCMatrix"))
stopifnot(nrow(normalizedData) == nrow(locationData))
stopifnot(nrow(normalizedData) == length(cellTypes))
stopifnot(all(c("x", "y") %in% colnames(locationData)))
```

- `normalizedData`: cells x genes matrix (rows = cells, columns = genes). Must have rownames (cell IDs) and colnames (gene names).
- `locationData`: data.frame with `x` and `y` columns. Rownames must match `normalizedData`.
- `metaData`: data.frame with cell-level annotations. Rownames must match.
- `cellTypes`: character vector, same length as `nrow(normalizedData)`.

### 2. Cell counts

CoPro needs enough cells per cell type for stable CCA. Check counts after subsetting:

```r
table(cellTypes[cellTypes %in% cellTypesOfInterest])
```

- **Minimum**: ~200 cells per cell type for basic results
- **Recommended**: 1,000+ cells per cell type for stable gene weights and reproducible transfer
- **Warning sign**: if one cell type has <10% of the cells of another, the CCA may be dominated by the larger group

For multi-slide data, also check per-slide counts---a slide with very few cells of one type will contribute noise.

### 3. Normalization check

CoPro expects **normalized, log-transformed** expression. Raw counts will produce meaningless results. Check for signs of unnormalized data:

```r
# Raw counts: values are integers, large range
range(normalizedData)       # should NOT be e.g. 0-50000
all(normalizedData == floor(normalizedData))  # should be FALSE for log-transformed data

# Quick check: column means should be roughly comparable across genes
summary(colMeans(normalizedData))  # large spread suggests no normalization
```

Recommended preprocessing (before CoPro):
- Library-size or DESeq2 size-factor normalization
- Log1p transformation (`log(x + 1)`)
- Optional: cap extreme values at the 95th percentile per gene (as in the organoid vignette)

Do NOT pass Seurat's `ScaleData` output (zero-centered, unit-variance)---CoPro handles scaling internally via `computePCA(center = TRUE, scale. = TRUE)`.

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
- For multi-slide: confirm `slideID` correctly partitions cells---`plot(x, y, col = as.factor(slideID))` should show non-overlapping groups
- If coordinates are in pixels, `normalizeDistance = TRUE` (default) will handle scale differences

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

## Post-run biological sanity checks

After running the pipeline, verify that results reflect biology rather than technical artifacts. Do not skip these---CoPro will always produce numbers, even on garbage input.

### 1. Normalized correlation magnitude

```r
ncorr <- getNormCorr(obj)
print(ncorr[ncorr$CC_index == 1, ])
```

- **Healthy signal**: normalized correlation 0.05--0.3 for CC1 (dataset-dependent)
- **Suspiciously high** (>0.5): likely a batch effect, library-size artifact, or cell type identity signal masquerading as spatial co-progression. Inspect the in-situ plot---does the pattern align with known biology, or does it just separate two clusters?
- **Near zero across all sigma**: the cell types may not spatially co-vary at the scales tested. Try different sigma values or check whether the cell types are actually co-localized.

### 2. Spatial pattern check

Plot cell scores in situ and ask: does this look like a biologically meaningful spatial gradient, or a technical artifact?

```r
cs <- getCellScoresInSitu(obj, sigmaValueChoice = sigma_opt, ccIndex = 1)
# Plot in situ (see Visualization section)
```

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
key <- paste0("geneScores|sigma", sigma_opt, "|CellType")
gs <- obj@geneScoresRegression[[key]][, 1]
top_pos <- head(sort(gs, decreasing = TRUE), 10)
top_neg <- head(sort(gs), 10)
print(top_pos)
print(top_neg)
```

- **Good sign**: top genes include known markers for the expected spatial gradient (e.g., pericentral/periportal markers in liver, proximal/distal markers in kidney)
- **Bad sign**: top genes are all ribosomal (`Rpl*`, `Rps*`), mitochondrial (`mt-*`), or housekeeping genes---these usually indicate a technical quality axis, not biology
- **Bad sign**: top genes are all sex-linked (`Xist`, `Ddx3y`)---CoPro found a sex batch effect
- Cross-reference a few top genes with the literature or gene ontology to confirm the axis has a coherent biological interpretation

### 4. Cross-type correlation consistency

For cross-cell-type analyses, verify the correlation between cell types:

```r
df <- getCorrTwoTypes(obj, sigmaValueChoice = sigma_opt,
                      cellTypeA = "TypeA", cellTypeB = "TypeB", ccIndex = 1)
cor(df$AK, df$B)
```

- The scatter should show a clear positive trend (CoPro maximizes this)
- If the correlation is driven by a handful of outlier cells, the result may not be robust
- For multi-slide data, check that the correlation holds per slide, not just in aggregate (aggregate correlation can be inflated by batch-level mean shifts)

### 5. Permutation test (statistical validation)

Run permutation testing to confirm the signal exceeds spatial null:

```r
obj <- runSkrCCAPermu(obj, nPermu = 100, permu_method = "bin",
                       num_bins_x = 10, num_bins_y = 10)
obj <- computeNormalizedCorrelationPermu(obj)
```

If the observed normalized correlation is not clearly above the permutation distribution, the signal may be spurious.

### 6. Multi-slide consistency

For datasets with biological replicates, the strongest validation is cross-slide consistency:
- Do the same genes rank highly across slides?
- Do transferred scores reproduce the spatial pattern?
- Are normalized correlations comparable across slides?

Inconsistent results across replicates suggest overfitting to slide-specific technical variation.

## Score transfer verification

After transferring scores to other slides or scRNA-seq, verify the transfer worked. A bad transfer is worse than no transfer---it gives false confidence in a meaningless axis.

### Transfer to other spatial slides

```r
tar_scores <- getTransferCellScores(
  ref_obj = ref_obj, tar_obj = tar_obj,
  sigma_choice = sigma_opt, gene_score_type = "regression"
)
```

**Checks**:

1. **Transferred normalized correlation**: should be comparable to the reference slide's value. A large drop (e.g., reference 0.15 vs transferred 0.03) means the pattern does not generalize.

```r
tar_ncorr <- getTransferNormCorr(tar_obj = tar_obj,
                                  transfer_cell_scores = tar_scores,
                                  sigma_choice = sigma_opt)
```

2. **Visual comparison**: plot transferred scores in situ on the target slide. The spatial pattern should resemble the reference---same gradient direction, similar regions highlighted.

3. **Gene overlap**: check how many genes were retained during transfer. If the gene panels differ substantially between slides, few genes survive intersection and the transfer is unreliable.

```r
# Check the transfer output messages for "retaining N genes"
# If N is very small relative to the panel size, the transfer is underpowered
```

### Transfer to scRNA-seq

```r
transferred <- transfer_scores(
  mat_A = spatial_expr, mat_B = scrna_expr,
  gs_ct = gene_weights,
  use_quantile_normalization = TRUE,
  gs_weight_threshold = 0
)
```

**Checks**:

1. **Known marker validation**: if the spatial axis has a known biological interpretation (e.g., corticomedullary), check that known markers separate correctly in the transferred scores. For example, in the kidney vignette, transferred nephron scores should monotonically increase from PTS1 to DCT-CNT.

```r
# Boxplot of transferred scores by known subtype
boxplot(transferred_score ~ subtype, data = validation_df)
```

2. **Score distribution**: transferred scores should have a unimodal, roughly continuous distribution. A bimodal distribution with a gap suggests the transfer is capturing cell-type identity rather than a gradient.

3. **Full-transcriptome regression**: after transferring, regress all scRNA-seq genes against the transferred score. Check that significant genes (FDR < 0.05) include biologically relevant genes, not just the panel genes that defined the axis.

```r
# Expect thousands of significant genes if the axis is real
# If only the panel genes are significant, the transfer may be circular
sum(reg_results$fdr < 0.05)
```

4. **UMAP sanity**: color the scRNA-seq UMAP by transferred score. The gradient should align with known biology (e.g., differentiation trajectory), not with technical variables (batch, library size, percent mitochondrial).

### Common transfer pitfalls

- **Too few shared genes**: if <100 genes overlap between reference and target, the transfer is underpowered. Consider using a broader spatial panel or a different reference.
- **Batch effects in scRNA-seq**: quantile normalization helps, but large batch effects can still distort transferred scores. Check whether scores correlate with batch labels.
- **Scale mismatch**: transferred scores from regression gene weights have larger magnitude than PCA-based scores. This does not affect relative cell ordering, but absolute values are not directly comparable between methods.

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
