## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)

## ----message=FALSE, warning=FALSE, results="hide"-----------------------------
# library("CoPro")
# library_paths <- c(
#   "fields", "ggplot2", "tictoc", "irlba",
#   "Matrix", "ggsci", "ape"
# )
# # Load required libraries
# lapply(library_paths, library, character.only = TRUE)

## -----------------------------------------------------------------------------
# n_pca <- 40 ## number of top PCs to retain
# cell_types <- c("061 STR D1 Gaba", "062 STR D2 Gaba")
# n_cell_types <- length(cell_types)
# 
# horizontal_dist_scale <- 1
# vertical_dist_scale <- 1

## ----eval=FALSE---------------------------------------------------------------
# ra <- readRDS("/Users/zhenmiao/Dropbox/Zhuang-ABCA-1.054_1/Zhuang_ABCA_1.054_subset_data.rds")
# meta <- readRDS("/Users/zhenmiao/Dropbox/Zhuang-ABCA-1.054_1/Zhuang_ABCA_1.054_subset_metadata.rds")

## -----------------------------------------------------------------------------
# subset_cells <- meta$subclass %in% cell_types
# 
# ra <- ra[subset_cells, ]
# meta <- meta[subset_cells, ]
# rownames(meta) <- meta$cell_label
# ## make sure coordinates are numeric
# meta$x <- as.numeric(meta$x)
# meta$y <- as.numeric(meta$y)

## -----------------------------------------------------------------------------
# ggplot(data = meta) +
#   geom_point(aes(x = x, y = y, color = subclass)) +
#   theme_minimal() +
#   coord_fixed(ratio = 1)

## ----create_object------------------------------------------------------------
# location_data <- meta[, c("x", "y")]
# cell_type_data <- meta[, "subclass"]
# br <- newCoProSingle(
#   normalizedData = ra,
#   locationData = location_data,
#   metaData = meta,
#   cellTypes = cell_type_data
# )
# br <- subsetData(br, cellTypesOfInterest = cell_types)

## ----distCompute--------------------------------------------------------------
# br <- computePCA(br, nPCA = n_pca, center = TRUE, scale. = TRUE)
# # sigma_square_choice <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)
# sigma_value_choice <- c(0.1, 0.14, 0.2, 0.5)
# ## By default, we normalize the distance
# br <- computeDistance(br, distType = "Euclidean2D", normalizeDistance = FALSE)
# br <- computeKernelMatrix(br, sigmaValues = sigma_value_choice)

## ----skrCCA-------------------------------------------------------------------
# br <- runSkrCCA(br, scalePCs = TRUE, maxIter = 500)
# br <- computeNormalizedCorrelation(br)

## ----geneCellScores-----------------------------------------------------------
# br <- computeGeneAndCellScores(br)

## ----sigma_choice-------------------------------------------------------------
# ncorr <- getNormCorr(br)
# 
# ggplot(data = ncorr, aes(
#   x = sigmaValues,
#   y = normalizedCorrelation, group = 1
# )) +
#   geom_point() +
#   geom_line() +
#   facet_wrap(vars(ct12,CC_index)) +
#   xlab("Sigma squared") +
#   ylab("Norm. Corr.") +
#   ggtitle(label = "Norm. Corr. across sigma squared values") +
#   theme_minimal()

## -----------------------------------------------------------------------------
# cell_A <- "061 STR D1 Gaba"
# cell_B <- "062 STR D2 Gaba"
# 
# df_crtt <- getCorrTwoTypes(br,
#   sigmaValueChoice = 0.14,
#   cellTypeA = cell_A,
#   cellTypeB = cell_B
# )
# 
# ggplot(df_crtt) +
#   geom_point(aes(x = AK, y = B)) +
#   ggtitle(paste("Correlation plot between", cell_A, "and", cell_B)) +
#   xlab(paste(cell_A, "%*% Kernal_AB")) +
#   ylab(cell_B) +
#   theme_minimal()
# 

## -----------------------------------------------------------------------------
# br_cs <- getCellScoresInSitu(br,
#   sigmaValueChoice = 0.14
# )

## -----------------------------------------------------------------------------
# ggplot(data = br_cs) +
#   geom_point(aes(
#     x = x,
#     y = y,
#     color = cellScores_b
#   ), size = 0.8) +
#   scale_color_discrete(
#     type =
#       c("darkgreen", "#d4f8d4", "darkred", "#ffd8d8")
#   ) +
#   coord_fixed() +
#   theme_minimal()

## -----------------------------------------------------------------------------
# ggplot(data = br_cs) +
#   geom_point(aes(
#     x = x,
#     y = y,
#     color = cellScores
#   ), size = 0.8) +
#   coord_fixed() +
#   theme_minimal()

## -----------------------------------------------------------------------------
# br <- runSkrCCAPermu(br, nPermu = 5L, permu_method = "bin",
#                      num_bins_x = 10, num_bins_y = 10)
# 
# br <- computeNormalizedCorrelationPermu(br, tol = 1e-3)
# 
# nc_permu <- br@normalizedCorrelationPermu
# # names(nc_permu) <- NULL
# nc_permu <- do.call(rbind, nc_permu)
# print(nc_permu)

