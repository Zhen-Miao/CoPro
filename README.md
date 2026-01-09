
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CoPro <img src="man/figures/copro-refined-logo-final.jpg" align="right" width="150" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/Zhen-Miao/CoPro/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Zhen-Miao/CoPro/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

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

## Installation

You can install the current version of CoPro from
[GitHub](https://github.com/Zhen-Miao/CoPro) with:

``` r
# install.packages("devtools")
devtools::install_github("Zhen-Miao/CoPro")
```

## Quick Start

Here’s a minimal example to get you started with CoPro:

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
obj <- computeDistance(obj, distType = "Euclidean2D")
obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1, 0.2))
obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2)
obj <- computeNormalizedCorrelation(obj)
obj <- computeGeneAndCellScores(obj)

# Get results
cell_scores <- getCellScores(obj, sigma = obj@sigmaValueChoice, cellType = "TypeA")
```

## Detailed Example

Below is a complete workflow using simulated data to demonstrate CoPro’s
functionality.

### Generate Example Data

``` r
library(CoPro)
set.seed(42)

# Simulate 200 cells with 100 genes
n_cells <- 200
n_genes <- 100

# Generate cell types
cell_types <- sample(c("Epithelial", "Immune"), n_cells, replace = TRUE)
cell_ids <- paste0("cell_", seq_len(n_cells))
gene_names <- paste0("Gene", seq_len(n_genes))

# Generate spatial coordinates
location_data <- data.frame(
  x = runif(n_cells, 0, 10),
  y = runif(n_cells, 0, 10),
  row.names = cell_ids
)

# Generate expression data with some spatial pattern
expr_data <- matrix(
  rnorm(n_cells * n_genes, mean = 2, sd = 1),
  nrow = n_cells, ncol = n_genes,
  dimnames = list(cell_ids, gene_names)
)
expr_data[expr_data < 0] <- 0

# Add spatial pattern to first 20 genes
dist_mat <- as.matrix(dist(location_data[, c("x", "y")]))
kernel_mat <- exp(-0.5 * (dist_mat / 1.5)^2)
kernel_mat <- kernel_mat / rowSums(kernel_mat)
pattern <- as.vector(kernel_mat %*% rnorm(n_cells))
for (g in 1:20) {
  expr_data[, g] <- expr_data[, g] + pattern * 0.5
}

# Create metadata
meta_data <- data.frame(
  cell_id = cell_ids,
  cell_type = cell_types,
  row.names = cell_ids
)
```

### Create CoPro Object and Run Analysis

``` r
# Create CoPro object
copro_obj <- newCoProSingle(
  normalizedData = expr_data,
  locationData = location_data,
  metaData = meta_data,
  cellTypes = cell_types
)

# Subset to cell types of interest
copro_obj <- subsetData(copro_obj, 
                        cellTypesOfInterest = c("Epithelial", "Immune"))

# Compute PCA
copro_obj <- computePCA(copro_obj, nPCA = 20, center = TRUE, scale. = TRUE)

# Compute distance matrix
copro_obj <- computeDistance(copro_obj, 
                             distType = "Euclidean2D",
                             normalizeDistance = TRUE,
                             verbose = FALSE)
#> normalizeDistance is set to TRUE, so distance will be normalized, so that 0.01 percentile distance will be scaled to 0.01
#> The scaling factor for normalizing distance is 0.06342367

# Compute kernel matrices for multiple sigma values
copro_obj <- computeKernelMatrix(copro_obj, 
                                 sigmaValues = c(0.05, 0.1, 0.2),
                                 verbose = FALSE)

# Run skrCCA
copro_obj <- runSkrCCA(copro_obj, 
                       scalePCs = TRUE, 
                       nCC = 2, 
                       maxIter = 200)
#> [1] "Convergence reached at 27 iterations (Max diff = 8.164e-06 )"
#> [1] "Convergence reached at 0 iterations (Max diff = 3.607e-08 )"
#> [1] "Convergence reached at 9 iterations (Max diff = 3.613e-06 )"
#> [1] "Convergence reached at 0 iterations (Max diff = 4.922e-11 )"
#> [1] "Convergence reached at 11 iterations (Max diff = 7.557e-06 )"
#> [1] "Convergence reached at 1 iterations (Max diff = 1.665e-16 )"

# Compute normalized correlation to find optimal sigma
copro_obj <- computeNormalizedCorrelation(copro_obj)
#> Calculating spectral norms,  depending on the data size, this may take a while. 
#> Finished calculating spectral norms

# Compute gene and cell scores
copro_obj <- computeGeneAndCellScores(copro_obj)
```

### Examine Results

``` r
# Get the optimal sigma value
optimal_sigma <- copro_obj@sigmaValueChoice
cat("Optimal sigma value:", optimal_sigma, "\n")
#> Optimal sigma value: 0.05

# Get normalized correlation results
ncorr <- getNormCorr(copro_obj)
print(head(ncorr))
#>              sigmaValues  cellType1 cellType2 CC_index normalizedCorrelation
#> sigma_0.05.1        0.05 Epithelial    Immune        1            0.21353667
#> sigma_0.05.2        0.05 Epithelial    Immune        2            0.17627646
#> sigma_0.1.1          0.1 Epithelial    Immune        1            0.17049168
#> sigma_0.1.2          0.1 Epithelial    Immune        2            0.10052308
#> sigma_0.2.1          0.2 Epithelial    Immune        1            0.10826527
#> sigma_0.2.2          0.2 Epithelial    Immune        2            0.05779051
#>                           ct12
#> sigma_0.05.1 Epithelial-Immune
#> sigma_0.05.2 Epithelial-Immune
#> sigma_0.1.1  Epithelial-Immune
#> sigma_0.1.2  Epithelial-Immune
#> sigma_0.2.1  Epithelial-Immune
#> sigma_0.2.2  Epithelial-Immune
```

### Visualize Cell Scores

``` r
library(ggplot2)

# Get cell scores for visualization
cell_scores_df <- getCellScoresInSitu(copro_obj, sigmaValueChoice = optimal_sigma)

# Plot cell scores in spatial coordinates
ggplot(cell_scores_df, aes(x = x, y = y, color = cellScores)) +
  geom_point(size = 1.5) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(title = "Spatial Cell Scores",
       color = "Score") +
  coord_fixed()
```

## Within-Cell-Type Analysis

CoPro can also detect spatial patterns within a single cell type:

``` r
# Create object for single cell type analysis
copro_single <- newCoProSingle(
  normalizedData = expr_data,
  locationData = location_data,
  metaData = meta_data,
  cellTypes = cell_types
)

# Subset to single cell type
copro_single <- subsetData(copro_single, cellTypesOfInterest = c("Epithelial"))

# Run the same pipeline
copro_single <- computePCA(copro_single, nPCA = 20)
copro_single <- computeDistance(copro_single, distType = "Euclidean2D", verbose = FALSE)
copro_single <- computeKernelMatrix(copro_single, sigmaValues = c(0.1, 0.2), verbose = FALSE)
copro_single <- runSkrCCA(copro_single, scalePCs = TRUE, nCC = 2)
copro_single <- computeNormalizedCorrelation(copro_single)
copro_single <- computeGeneAndCellScores(copro_single)
```

## Multi-Slide Analysis

For experiments with multiple tissue slides:

``` r
# Create CoProMulti object
copro_multi <- newCoProMulti(
  normalizedData = combined_expr_matrix,
  locationData = combined_location_data,
  metaData = combined_metadata,
  cellTypes = combined_cell_types,
  slideID = slide_identifiers  # Vector indicating which slide each cell belongs to
)

# The workflow is the same - CoPro automatically handles multi-slide data
copro_multi <- subsetData(copro_multi, cellTypesOfInterest = c("TypeA", "TypeB"))
copro_multi <- computePCA(copro_multi, nPCA = 30)
copro_multi <- computeDistance(copro_multi, distType = "Euclidean2D", verbose = FALSE)
copro_multi <- computeKernelMatrix(copro_multi, sigmaValues = c(0.1), verbose = FALSE)
copro_multi <- runSkrCCA(copro_multi, scalePCs = TRUE, nCC = 2)
copro_multi <- computeNormalizedCorrelation(copro_multi)
copro_multi <- computeGeneAndCellScores(copro_multi)
```

## Key Parameters

- **`nPCA`**: Number of principal components to retain. Higher values
  capture more variance but may include noise. Typical range: 20-50.

- **`sigmaValues`**: Spatial neighborhood bandwidth parameter(s) for the
  Gaussian kernel. Smaller values capture local patterns; larger values
  capture broader patterns. CoPro can test multiple values and select
  the optimal one.

- **`nCC`**: Number of canonical components to compute. The first
  component captures the strongest spatial co-progression signal.

- **`normalizeDistance`**: When `TRUE` (default), normalizes distances
  so that the 0.01 percentile becomes 0.01, ensuring consistent kernel
  behavior across datasets with different spatial scales.

## Output

CoPro provides several types of output:

1.  **Cell Scores**: Per-cell scores indicating the strength of
    co-progression pattern
2.  **Gene Scores**: Per-gene scores indicating contribution to the
    co-progression pattern  
3.  **Normalized Correlation**: A measure of co-progression strength
    (0-1 scale)
4.  **Optimal Sigma**: The spatial scale at which co-progression is
    strongest

## Citation

If you use CoPro in your research, please cite:

> \[Citation information will be added upon publication\]

## Getting Help

- Report bugs and request features on [GitHub
  Issues](https://github.com/Zhen-Miao/CoPro/issues)
- See the [package documentation](https://zhen-miao.github.io/CoPro/)
  for detailed function references
