pkgname <- "CoPro"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "CoPro-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('CoPro')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("calculate_pvalue")
### * calculate_pvalue

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calculate_pvalue
### Title: Calculate P-value from Permutation Results
### Aliases: calculate_pvalue

### ** Examples

## Not run: 
##D result <- calculate_pvalue(br, cc_index = 1, alternative = "greater")
##D print(paste("P-value:", result$p_value))
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calculate_pvalue", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("computeBidirCorrelation")
### * computeBidirCorrelation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: computeBidirCorrelation
### Title: Compute Bidirectional Correlation
### Aliases: computeBidirCorrelation computeBidirCorrelation,CoPro-method
###   computeBidirCorrelation,CoProMulti-method

### ** Examples

# Assuming `obj` is a prepared CoProSingle with PCA, kernels, and skrCCA:
# obj <- computeBidirCorrelation(obj)

# For CoProMulti per-slide results:
# objm <- computeBidirCorrelation(objm, calculationMode = "perSlide")

# For CoProMulti aggregate results:
# objm <- computeBidirCorrelation(objm, calculationMode = "aggregate")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("computeBidirCorrelation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("computeNormalizedCorrelationPermu")
### * computeNormalizedCorrelationPermu

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: computeNormalizedCorrelationPermu
### Title: Compute Normalized Correlation for Permutation Results
### Aliases: computeNormalizedCorrelationPermu

### ** Examples

## Not run: 
##D # After running permutation testing
##D br <- computeNormalizedCorrelationPermu(br)
##D 
##D # Extract permutation values and calculate p-value
##D permu_values <- sapply(br@normalizedCorrelationPermu,
##D                        function(x) x$normalizedCorrelation[1])
##D observed <- max(getNormCorr(br)$normalizedCorrelation)
##D 
##D # One-sided p-value (testing if observed > permutation)
##D p_value <- mean(permu_values >= observed)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("computeNormalizedCorrelationPermu", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("computeSelfBidirCorr")
### * computeSelfBidirCorr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: computeSelfBidirCorr
### Title: Compute Self-Bidirectional Correlation using skrCCA Results
### Aliases: computeSelfBidirCorr

### ** Examples

## Not run: 
##D # Assuming you have a CoPro object with skrCCA results and self-kernels
##D object <- runSkrCCA(object)
##D object <- computeSelfDistance(object)
##D object <- computeSelfKernel(object, sigmaValues = c(0.01, 0.05, 0.1))
##D 
##D # Compute self-bidirectional correlation using native skrCCA results
##D self_bidir <- computeSelfBidirCorr(object, sigma_choice = 0.05)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("computeSelfBidirCorr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("computeSelfDistance")
### * computeSelfDistance

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: computeSelfDistance
### Title: Compute Self-Distance Matrices for Multiple Cell Types
### Aliases: computeSelfDistance computeSelfDistance,CoProSingle-method
###   computeSelfDistance,CoProMulti-method

### ** Examples

## Not run: 
##D # Assume you have a CoPro object with multiple cell types
##D # First compute cross-type distances
##D object <- computeDistance(object)
##D 
##D # Then add self-distances
##D object <- computeSelfDistance(object)
##D 
##D # Now you have both cross-type and self-type distance matrices
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("computeSelfDistance", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("computeSelfKernel")
### * computeSelfKernel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: computeSelfKernel
### Title: Compute Self-Kernel Matrices for Multiple Cell Types
### Aliases: computeSelfKernel computeSelfKernel,CoProSingle-method
###   computeSelfKernel,CoProMulti-method

### ** Examples

## Not run: 
##D # Assume you have a CoPro object with multiple cell types
##D # First compute cross-type distances and kernels
##D object <- computeDistance(object)
##D object <- computeKernelMatrix(object, sigmaValues = c(0.01, 0.05, 0.1))
##D 
##D # Then add self-distances and self-kernels
##D object <- computeSelfDistance(object)
##D object <- computeSelfKernel(object, sigmaValues = c(0.01, 0.05, 0.1))
##D 
##D # Now you have both cross-type and self-type kernel matrices
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("computeSelfKernel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute_ground_truth_ncorr")
### * compute_ground_truth_ncorr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_ground_truth_ncorr
### Title: Compute Normalized Correlation from Ground Truth Scores
### Aliases: compute_ground_truth_ncorr

### ** Examples

## Not run: 
##D # After running CoPro analysis on simulated data
##D # Get ground truth scores from metadata
##D meta <- br@metaDataSub
##D gt_scores_A <- meta$smoothed_score[meta$cell_type == "A"]
##D gt_scores_B <- meta$smoothed_score[meta$cell_type == "B"]
##D 
##D # Compute ground truth normalized correlation
##D gt_ncorr <- compute_ground_truth_ncorr(
##D   object = br,
##D   scores_ct1 = gt_scores_A,
##D   scores_ct2 = gt_scores_B,
##D   cellType1 = "A",
##D   cellType2 = "B"
##D )
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_ground_truth_ncorr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ensureBidirCorrelationSlot")
### * ensureBidirCorrelationSlot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ensureBidirCorrelationSlot
### Title: Ensure object has bidirCorrelation slot
### Aliases: ensureBidirCorrelationSlot

### ** Examples

# Upgrade legacy object to include bidirCorrelation slot
# obj <- ensureBidirCorrelationSlot(obj)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ensureBidirCorrelationSlot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generate_toroidal_permutations")
### * generate_toroidal_permutations

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generate_toroidal_permutations
### Title: Generate Toroidal Shift Permutation Indices
### Aliases: generate_toroidal_permutations

### ** Examples

## Not run: 
##D # Create example location data
##D loc_data <- data.frame(
##D   x = runif(100, 0, 10),
##D   y = runif(100, 0, 10),
##D   cell_ID = paste0("cell_", 1:100)
##D )
##D 
##D # Generate 100 toroidal permutations
##D perm_matrix <- generate_toroidal_permutations(loc_data, n_permu = 100)
##D 
##D # Apply first permutation
##D permuted_cells <- loc_data$cell_ID[perm_matrix[, 1]]
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generate_toroidal_permutations", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getCellScores")
### * getCellScores

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getCellScores
### Title: Get cell scores from CoPro object
### Aliases: getCellScores getCellScores,CoProSingle-method
###   getCellScores,CoProMulti-method

### ** Examples

## Not run: 
##D # Get all cell scores for a specific sigma and cell type
##D scores <- getCellScores(object, sigma = 0.1, cellType = "TypeA")
##D 
##D # Get scores for a specific canonical component
##D cc1_scores <- getCellScores(object, sigma = 0.1, cellType = "TypeA", ccIndex = 1)
##D 
##D # Get scores for specific cells
##D specific_scores <- getCellScores(object, sigma = 0.1, cellType = "TypeA", 
##D                                  cells = c("cell_1", "cell_2"))
##D 
##D # For multi-slide object, specify slide
##D slide_scores <- getCellScores(object, sigma = 0.1, cellType = "TypeA", 
##D                               slide = "slide1")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getCellScores", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getColocScores")
### * getColocScores

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getColocScores
### Title: Compute Colocalization Scores for All Cell Type Pairs
### Aliases: getColocScores

### ** Examples

## Not run: 
##D # Basic usage with default parameters
##D coloc_results <- getColocScores(object)
##D 
##D # Custom parameters for high-resolution data
##D coloc_results <- getColocScores(object, 
##D                                r_um_range = c(5, 30),
##D                                pixel_size_um = 0.325,  # Convert from pixels to microns
##D                                cell_diam_um = 8,
##D                                nsim = 99)
##D                                
##D # For multi-slide objects
##D coloc_results <- getColocScores(multi_object)
##D # Results will include slideID column
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getColocScores", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getTransferBidirCorr")
### * getTransferBidirCorr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getTransferBidirCorr
### Title: Compute Bidirectional Correlation from Transferred Cell Scores
### Aliases: getTransferBidirCorr

### ** Examples

# Assuming `tar_obj` is prepared and `trans_scores` was computed with
# getTransferCellScores(..., agg_cell_type = FALSE)
# res <- getTransferBidirCorr(tar_obj, trans_scores, sigma_choice = 2.0)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getTransferBidirCorr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getTransferNormCorr")
### * getTransferNormCorr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getTransferNormCorr
### Title: Compute Normalized Correlation from Transferred Cell Scores
### Aliases: getTransferNormCorr

### ** Examples

# Assuming `tar_obj` is prepared and `trans_scores` was computed with
# getTransferCellScores(..., agg_cell_type = FALSE)
# res <- getTransferNormCorr(tar_obj, trans_scores, sigma_choice = 2.0)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getTransferNormCorr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getTransferSelfBidirCorr")
### * getTransferSelfBidirCorr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getTransferSelfBidirCorr
### Title: Compute Self-Bidirectional Correlation from Transferred Cell
###   Scores
### Aliases: getTransferSelfBidirCorr

### ** Examples

## Not run: 
##D # Assuming you have a CoPro object with multiple cell types
##D # First compute standard workflow
##D object <- computeDistance(object)
##D object <- computeKernelMatrix(object, sigmaValues = c(0.01, 0.05, 0.1))
##D 
##D # Add self-distances and self-kernels
##D object <- computeSelfDistance(object)
##D object <- computeSelfKernel(object, sigmaValues = c(0.01, 0.05, 0.1))
##D 
##D # Compute transferred cell scores
##D trans_scores <- getTransferCellScores(ref_obj, tar_obj, sigma_choice = 0.05, 
##D                                      agg_cell_type = FALSE)
##D 
##D # Compute self-bidirectional correlation from transferred scores
##D self_bidir <- getTransferSelfBidirCorr(tar_obj, trans_scores, sigma_choice = 0.05)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getTransferSelfBidirCorr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotG12Functions")
### * plotG12Functions

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotG12Functions
### Title: Plot g_12(r) pair correlation functions for colocalization
###   analysis
### Aliases: plotG12Functions

### ** Examples

## Not run: 
##D # Basic usage
##D g12_plots <- plotG12Functions(object)
##D 
##D # Custom parameters with individual plots
##D g12_plots <- plotG12Functions(object, 
##D                              r_um_range = c(5, 40),
##D                              plot_type = "individual",
##D                              include_confidence = TRUE)
##D 
##D # Access the plot and data
##D print(g12_plots$plot)
##D head(g12_plots$data)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotG12Functions", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("resample_spatial")
### * resample_spatial

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: resample_spatial
### Title: Spatial Resampling for Permutation Testing
### Aliases: resample_spatial

### ** Examples

## Not run: 
##D # Create example data
##D loc_data <- data.frame(
##D   x = runif(100, 0, 10),
##D   y = runif(100, 0, 10),
##D   cell_ID = paste0("cell_", 1:100)
##D )
##D 
##D # Check bin distribution first
##D diagnose_bin_distribution(loc_data, num_bins_x = 5, num_bins_y = 5)
##D 
##D # Perform spatial resampling (random within tile)
##D resampled <- resample_spatial(loc_data, num_bins_x = 5, num_bins_y = 5)
##D 
##D # Perform spatial resampling (quantile-matched within tile)
##D resampled_matched <- resample_spatial(loc_data, num_bins_x = 5, num_bins_y = 5,
##D                                       match_quantile = TRUE)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("resample_spatial", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("runSkrCCAPermu")
### * runSkrCCAPermu

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: runSkrCCAPermu
### Title: Run Spatial CCA with Permutation Testing
### Aliases: runSkrCCAPermu

### ** Examples

## Not run: 
##D # After running standard CoPro analysis
##D br <- runSkrCCA(br, scalePCs = TRUE)
##D br <- computeNormalizedCorrelation(br)
##D 
##D # Standard permutation (only permute second cell type)
##D br <- runSkrCCAPermu(br, nPermu = 100, permu_method = "bin",
##D                      permu_which = "second_only")
##D 
##D # Conservative permutation (lower FPR - better preserves spatial structure)
##D br <- runSkrCCAPermu(br, nPermu = 100, conservative = TRUE)
##D 
##D # Manual conservative settings: more bins + quantile matching
##D br <- runSkrCCAPermu(br, nPermu = 100, permu_method = "bin",
##D                      permu_which = "second_only",
##D                      num_bins_x = 15, num_bins_y = 15,
##D                      match_quantile = TRUE)
##D 
##D # PC-space permutation (like DIALOGUE) - shuffles within PC dimensions
##D br <- runSkrCCAPermu(br, nPermu = 100, permu_method = "pc",
##D                      permu_which = "second_only")
##D 
##D br <- computeNormalizedCorrelationPermu(br)
##D 
##D # Calculate p-value
##D observed <- max(getNormCorr(br)$normalizedCorrelation)
##D permu_values <- sapply(br@normalizedCorrelationPermu,
##D                        function(x) x$normalizedCorrelation[1])
##D p_value <- mean(permu_values >= observed)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("runSkrCCAPermu", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("runSkrCCAPermu_FairSigma")
### * runSkrCCAPermu_FairSigma

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: runSkrCCAPermu_FairSigma
### Title: Run Permutation Test with Fair Sigma Selection
### Aliases: runSkrCCAPermu_FairSigma

### ** Examples

## Not run: 
##D # After running standard CoPro analysis
##D br <- runSkrCCA(br, scalePCs = TRUE)
##D br <- computeNormalizedCorrelation(br)
##D 
##D # Run fair sigma permutation test
##D br <- runSkrCCAPermu_FairSigma(br, nPermu = 100,
##D                                 permu_method = "toroidal")
##D 
##D # Calculate p-value
##D result <- calculate_pvalue(br)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("runSkrCCAPermu_FairSigma", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
