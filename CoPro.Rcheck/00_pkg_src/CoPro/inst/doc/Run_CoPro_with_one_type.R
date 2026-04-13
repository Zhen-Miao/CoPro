## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)

## ----message=FALSE, warning=FALSE, results="hide"-----------------------------
# library("CoPro")
# library_paths <- c(
#   "fields", "ggplot2", "tictoc", "irlba",
#   "Matrix", "ggsci", "ape", "DESeq2","SummarizedExperiment"
# )
# # Load required libraries
# lapply(library_paths, library, character.only = TRUE)

## ----set_parameters-----------------------------------------------------------
# n_pca = 30
# cell_types <- c('Epithelial')
# n_cell_types <- length(cell_types)
# num_bins_x <- 10
# num_bins_y <- 10
# 

## ----load_data, eval=FALSE----------------------------------------------------
# loc <- "/Users/zhenmiao/Library/CloudStorage/Dropbox/DIALOGUE_plus project/Raj_lab_data/72hr/roi1/"
# 
# ra <- read.csv(paste0(loc, 'output/cell_by_gene/cell_by_gene.csv'))
# ra <- ra[,2:ncol(ra)]
# meta_cell <- read.csv(paste0(loc,'output/attributes/cell_attributes.csv'))
# 
# meta_cell$label <- paste0("cell_", meta_cell$label)
# rownames(ra) <- meta_cell$label
# rownames(meta_cell) <- meta_cell$label
# 

## ----filter-------------------------------------------------------------------
# ## filter expression outliers
# filter_out <- function(x) {
#   upper <- quantile(x, 0.95)
#   x[x > upper] <- upper
#   return(x)
# }
# ra <- apply(ra, MARGIN = 2, FUN = filter_out)

## ----normalize----------------------------------------------------------------
# ## normalize data by DEseq2
# counts_matrix <- t(ra) +1 ## count +1, as DEseq2 will need to take log
# mode(counts_matrix) <- "integer"
# 
# condition <- factor(rep("condition", ncol(counts_matrix)))
# colData <- data.frame(row.names = colnames(counts_matrix),
#                       condition = condition)
# 
# # Create DESeqDataSet object
# dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
#                               colData = colData,
#                               design = ~ 1)
# dds <- estimateSizeFactors(dds)
# normalized_counts2 <- t(ra) / sizeFactors(dds)
# ra_norm <- t(log1p(normalized_counts2))
# 

## ----show_in_situ-------------------------------------------------------------
# ggplot(meta_cell) +
#   geom_point(aes(x = center_x, y = center_y), color = "steelblue", size = 0.8) +
#   ggtitle("72h culture, ROI-1") +
#   theme_classic()

## ----create_object------------------------------------------------------------
# location_data <- meta_cell[, c("center_x", "center_y")]
# location_data$x <- location_data$center_x / 5000
# location_data$y <- location_data$center_y / 5000
# 
# location_data$x_bin <- cut(location_data$x, breaks = num_bins_x, labels = FALSE)
# location_data$y_bin <- cut(location_data$y, breaks = num_bins_y, labels = FALSE)
# meta_cell$"Celltype" <- "Epithelial"
# 
# cell_type_data <- meta_cell[, "Celltype"]
# 
# etr <- newCoProSingle(
#   normalizedData = ra_norm,
#   locationData = location_data,
#   metaData = meta_cell,
#   cellTypes = cell_type_data
# ) ## etr stands for enteroid
# cell_types_single <- "Epithelial"
# etr <- subsetData(etr, cellTypesOfInterest = cell_types_single)

## ----distCompute--------------------------------------------------------------
# etr <- computePCA(etr, nPCA = n_pca, center = TRUE, scale. = TRUE)
# ## By default, we normalize the distance
# etr <- computeDistance(etr, distType = "Euclidean2D", normalizeDistance = FALSE)
# 
# sigma_choice <- c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2)
# 
# ## compute kernel. Note that it is a different function for one cell type model.
# etr <- computeKernelMatrix(etr, sigmaValues = sigma_choice,
#                               upperQuantile = 0.85,
#                               normalizeKernel = FALSE, lowerLimit = 5e-7)

## ----skrCCA-------------------------------------------------------------------
# etr <- runSkrCCA(etr, scalePCs = TRUE, maxIter = 500, nCC = 4)
# etr <- computeNormalizedCorrelation(etr, tol = 1e-3)
# 
# 

## ----geneCellScores-----------------------------------------------------------
# etr <- computeGeneAndCellScores(etr)

## ----sigma_choice-------------------------------------------------------------
# ncorr <- getNormCorr(etr)
# 
# ggplot(data = ncorr, aes(
#   x = sigmaValues,
#   y = normalizedCorrelation, group = 1
# )) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(facets = ~ ct12 + CC_index,
#              labeller = purrr::partial(label_both, sep = " = ")) +
#   xlab("sigmaValues") +
#   ylab("Norm. Corr.") +
#   ggtitle(label = "Norm. Corr. across sigma values", subtitle = "Epi-Epi") +
#   theme_minimal()
# 

## -----------------------------------------------------------------------------
# cell_A <- cell_types_single
# df_crtt <- getCorrOneType(etr,
#                           sigmaValueChoice = 0.1,
#                           cellTypeA = cell_A, ccIndex = 1)
# 
# ggplot(df_crtt) +
#   geom_point(aes(x = AK, y = B), size = 0.5) +
#   ggtitle(paste("Correlation plot between", cell_A, "and", cell_A),
#           subtitle = paste("sigma_value = ", 0.1)) +
#   xlab(paste(cell_A, "%*% Kernal_AB")) +
#   ylab(cell_A) +
#   theme_minimal()
# 

## -----------------------------------------------------------------------------
# ## get the subset of second round cells. These are the cells with positive GFP
# ## expression values
# df_crtt$GFP <- ra[, "GFP"]
# df_crtt_sub <- df_crtt[df_crtt$GFP > 0, ]
# 
# ggplot(df_crtt_sub) +
#   geom_point(aes(x = AK, y = B, color = GFP), size = 0.5) +
#   ggtitle(paste("Correlation plot between", cell_A, "and", cell_A),
#           subtitle = paste("sigma_value = ", 0.1)) +
#   xlab(paste(cell_A, "%*% Kernal_AB")) +
#   ylab(cell_A) +
#   theme_minimal()

## -----------------------------------------------------------------------------
# 
# etr_cs <- getCellScoresInSitu(etr,
#                               sigmaValueChoice = 0.1
# )
# 

## -----------------------------------------------------------------------------
# ggplot(data = etr_cs) +
#   geom_point(aes(
#     x = x,
#     y = y,
#     color = cellScores_b
#   ), size = 0.8) +
#   scale_color_discrete(
#     type =
#       c("darkred", "#ffd8d8", "#1e17a4", "#5b79ff")
#   ) +
#   ggtitle("Raj Lab Epi one cell type score") +
#   coord_fixed() +
#   theme_minimal()

## -----------------------------------------------------------------------------
# ggplot(data = etr_cs) +
#   geom_point(aes(
#     x = x,
#     y = y,
#     color = cellScores
#   ), size = 0.8) +
#   coord_fixed() +
#   theme_minimal()

