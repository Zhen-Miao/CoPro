# Package index

## Object Creation

Create CoPro objects from spatial transcriptomics data

- [`newCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/newCoProSingle.md)
  : Function to create a new object
- [`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md)
  : Create a new CoProMulti object for Multi-Slide Analysis
- [`CreateCoPro()`](https://zhen-miao.github.io/CoPro/reference/CreateCoPro.md)
  : Create a CoPro object, automatically choosing Single vs Multi and
  splitting large slices.
- [`CoPro-class`](https://zhen-miao.github.io/CoPro/reference/CoPro-class.md)
  : CoPro object of spatial transcriptomics data
- [`subsetData()`](https://zhen-miao.github.io/CoPro/reference/subsetData.md)
  : subsetData
- [`show(`*`<CoPro>`*`)`](https://zhen-miao.github.io/CoPro/reference/show.md)
  : Show method for CoPro objects

## Core Pipeline

The main CoPro analysis workflow

- [`computePCA()`](https://zhen-miao.github.io/CoPro/reference/computePCAMulti.md)
  : Compute PCA on Integrated Multi-Slide Data
- [`computePCA(`*`<CoProSingle>`*`)`](https://zhen-miao.github.io/CoPro/reference/computePCA.md)
  [`computePCA(`*`<CoProMulti>`*`)`](https://zhen-miao.github.io/CoPro/reference/computePCA.md)
  : Compute PCA on Single-Slide Data
- [`computeDistance()`](https://zhen-miao.github.io/CoPro/reference/computeDistance.md)
  : computeDistance between pairs of cell types
- [`computeKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/computeKernelMatrix.md)
  : Compute Kernel Matrix for CoPro
- [`runSkrCCA()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCA.md)
  : runSkrCCA
- [`computeNormalizedCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelation.md)
  : Compute Normalized Correlation (approximation)
- [`computeGeneAndCellScores()`](https://zhen-miao.github.io/CoPro/reference/computeGeneAndCellScores.md)
  : computeGeneAndCellScores
- [`computeRegressionGeneScores()`](https://zhen-miao.github.io/CoPro/reference/computeRegressionGeneScores.md)
  : Compute regression-based gene scores

## Results Access

Retrieve analysis results

- [`getCellScores()`](https://zhen-miao.github.io/CoPro/reference/getCellScores.md)
  : Get cell scores from CoPro object
- [`getCellScoresInSitu()`](https://zhen-miao.github.io/CoPro/reference/getCellScoresInSitu.md)
  : Get cell score and location information as a data.frame
- [`getNormCorr()`](https://zhen-miao.github.io/CoPro/reference/getNormalizedCorrelation.md)
  : Get normalized correlation vs Sigma squared values
- [`getCorrOneType()`](https://zhen-miao.github.io/CoPro/reference/getCorrOneType.md)
  : Retrieve the Correlation within one cell types
- [`getCorrTwoTypes()`](https://zhen-miao.github.io/CoPro/reference/getCorrTwoTypes.md)
  : Retrieve the Correlation between two cell types
- [`getDistMat()`](https://zhen-miao.github.io/CoPro/reference/getDistMat.md)
  : Get Distance Matrix
- [`getSelfDistMat()`](https://zhen-miao.github.io/CoPro/reference/getSelfDistMat.md)
  : Get Self-Distance Matrix
- [`getKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/getKernelMatrix.md)
  : Get Kernel Matrix
- [`getSelfKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/getSelfKernelMatrix.md)
  : Get Self-Kernel Matrix
- [`getColocScores()`](https://zhen-miao.github.io/CoPro/reference/getColocScores.md)
  : Compute Colocalization Scores for All Cell Type Pairs
- [`getSlideID()`](https://zhen-miao.github.io/CoPro/reference/getSlideID.md)
  : Get slide IDs from CoPro object
- [`getSlideList()`](https://zhen-miao.github.io/CoPro/reference/getSlideList.md)
  : Get slide list from CoPro object
- [`isMultiSlide()`](https://zhen-miao.github.io/CoPro/reference/isMultiSlide.md)
  : Check if object is multi-slide

## Score Transfer

Transfer learned patterns to new slides or datasets

- [`getTransferCellScores()`](https://zhen-miao.github.io/CoPro/reference/getTransferCellScores.md)
  : Get cell score by transferring gene weights from another slide By
  default, quantile normalization is used to ensure distribution match
- [`getTransferNormCorr()`](https://zhen-miao.github.io/CoPro/reference/getTransferNormCorr.md)
  : Compute Normalized Correlation from Transferred Cell Scores
- [`getTransferBidirCorr()`](https://zhen-miao.github.io/CoPro/reference/getTransferBidirCorr.md)
  : Compute Bidirectional Correlation from Transferred Cell Scores
- [`getTransferSelfBidirCorr()`](https://zhen-miao.github.io/CoPro/reference/getTransferSelfBidirCorr.md)
  : Compute Self-Bidirectional Correlation from Transferred Cell Scores
- [`quantile_normalize()`](https://zhen-miao.github.io/CoPro/reference/quantile_normalize.md)
  : quantile normalize data
- [`transfer_scores()`](https://zhen-miao.github.io/CoPro/reference/transfer_scores.md)
  : Transfer cell scores between matrices using gene weights

## Within-Cell-Type Analysis

Self-correlation and single cell type spatial patterns

- [`computeSelfDistance()`](https://zhen-miao.github.io/CoPro/reference/computeSelfDistance.md)
  : Compute Self-Distance Matrices for Multiple Cell Types
- [`computeSelfKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSelfKernel.md)
  : Compute Self-Kernel Matrices for Multiple Cell Types
- [`computeSelfBidirCorr()`](https://zhen-miao.github.io/CoPro/reference/computeSelfBidirCorr.md)
  : Compute Self-Bidirectional Correlation using skrCCA Results
- [`computeBidirCorrelation()`](https://zhen-miao.github.io/CoPro/reference/computeBidirCorrelation.md)
  : Compute Bidirectional Correlation
- [`ensureBidirCorrelationSlot()`](https://zhen-miao.github.io/CoPro/reference/ensureBidirCorrelationSlot.md)
  : Ensure object has bidirCorrelation slot

## Gene-Level Analysis

Gene scoring, testing, and regression

- [`testGeneGLM()`](https://zhen-miao.github.io/CoPro/reference/testGeneGLM.md)
  : testGeneGLM
- [`smoothCellScoresMatrix()`](https://zhen-miao.github.io/CoPro/reference/smoothCellScoresMatrix.md)
  : Smooth the cell scores based on expression neighborhood

## Permutation Testing

Statistical significance via spatial permutations

- [`runSkrCCAPermu()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu.md)
  : Run Spatial CCA with Permutation Testing
- [`runSkrCCAPermu_FairSigma()`](https://zhen-miao.github.io/CoPro/reference/runSkrCCAPermu_FairSigma.md)
  : Run Permutation Test with Fair Sigma Selection
- [`computeNormalizedCorrelationPermu()`](https://zhen-miao.github.io/CoPro/reference/computeNormalizedCorrelationPermu.md)
  : Compute Normalized Correlation for Permutation Results
- [`resample_spatial()`](https://zhen-miao.github.io/CoPro/reference/resample_spatial.md)
  : Spatial Resampling for Permutation Testing
- [`generate_toroidal_permutations()`](https://zhen-miao.github.io/CoPro/reference/generate_toroidal_permutations.md)
  : Generate Toroidal Shift Permutation Indices
- [`compute_ground_truth_ncorr()`](https://zhen-miao.github.io/CoPro/reference/compute_ground_truth_ncorr.md)
  : Compute Normalized Correlation from Ground Truth Scores
- [`calculate_pvalue()`](https://zhen-miao.github.io/CoPro/reference/calculate_pvalue.md)
  : Calculate P-value from Permutation Results

## Visualization

Plotting functions

- [`plotG12Functions()`](https://zhen-miao.github.io/CoPro/reference/plotG12Functions.md)
  : Plot g_12(r) pair correlation functions for colocalization analysis
- [`diagnose_bin_distribution()`](https://zhen-miao.github.io/CoPro/reference/diagnose_bin_distribution.md)
  : Diagnose Bin Distribution

## Utilities

Helper and optimization functions

- [`copro_download_data()`](https://zhen-miao.github.io/CoPro/reference/copro_download_data.md)
  : Download example datasets for CoPro vignettes
- [`assignDistanceManually()`](https://zhen-miao.github.io/CoPro/reference/assignDistanceManually.md)
  : Assign distance matrix manually
- [`optimize_bilinear()`](https://zhen-miao.github.io/CoPro/reference/optimize_bilinear.md)
  : SkrCCA optimization function for multiple groups (Single Slide) -
  First Component Uses flat kernel structure for consistent data access
- [`optimize_bilinear_multi_slides()`](https://zhen-miao.github.io/CoPro/reference/optimize_bilinear_multi_slides.md)
  : Multi-slide SkrCCA optimization - First Component
- [`optimize_bilinear_n()`](https://zhen-miao.github.io/CoPro/reference/optimize_bilinear_n.md)
  : Run multi version of skrCCA to detect subsequent components (Single
  Slide) Uses flat kernel structure for consistent data access
- [`optimize_bilinear_n_multi_slides()`](https://zhen-miao.github.io/CoPro/reference/optimize_bilinear_n_multi_slides.md)
  : Multi-slide SkrCCA optimization - Multiple Components
