\name{getTransferSelfBidirCorr}
\alias{getTransferSelfBidirCorr}
\title{Compute Self-Bidirectional Correlation from Transferred Cell Scores}
\usage{
getTransferSelfBidirCorr(
  tar_obj,
  transfer_cell_scores,
  sigma_choice,
  calculationMode = NULL,
  normalize_K = c("row_or_col", "sinkhorn_knopp", "none"),
  filter_kernel = TRUE,
  K_row_sum_cutoff = 0.005,
  K_col_sum_cutoff = 0.005,
  sigma_choice_tar = NULL,
  verbose = TRUE
)
}
\arguments{
\item{tar_obj}{A \code{CoProSingle} or \code{CoProMulti} object containing self-kernel matrices
and metadata needed for alignment. Must have self-kernel matrices computed using
\code{computeSelfKernel()}.}

\item{transfer_cell_scores}{A named list of matrices, with one entry per cell type
(names must be the cell type names). Each matrix should be cells-by-CCs, where
rows are cell IDs and columns are \code{CC_1}, \code{CC_2}, ..., as returned by
\code{getTransferCellScores(agg_cell_type = FALSE)}.}

\item{sigma_choice}{Numeric scalar specifying the sigma value of the self-kernel to use.}

\item{calculationMode}{For \code{CoProMulti} objects only, either "perSlide" or
"aggregate". Ignored for \code{CoProSingle}. Default "perSlide" if \code{tar_obj}
is multi-slide.}

\item{normalize_K}{Character; method for normalizing the kernel matrix, one of 
"row_or_col", "sinkhorn_knopp", or "none". Default "row_or_col".}

\item{filter_kernel}{Logical; whether to filter the kernel matrix. Default TRUE.}

\item{K_row_sum_cutoff}{Numeric; cutoff for row sums when normalizing kernel matrix.
Default 5e-3.}

\item{K_col_sum_cutoff}{Numeric; cutoff for column sums when normalizing kernel matrix.
Default 5e-3.}

\item{sigma_choice_tar}{Numeric; sigma value for target object self-kernel matrices. 
If NULL (default), uses sigma_choice. Not recommended for general use.}

\item{verbose}{Logical; whether to print progress messages.}
}
\value{
A list with one element named \code{paste0("sigma_", sigma_choice)}, whose
value is a data.frame of results. For single-slide objects, the data.frame has
columns \code{sigmaValue}, \code{cellType}, \code{CC_index}, \code{selfBidirCorrelation}.
For multi-slide objects in \code{perSlide} mode, the data.frame additionally includes
\code{slideID}. For \code{aggregate} mode, the correlation column is named
\code{aggregateSelfCorrelation}.
}
\description{
Given transferred cell scores for each cell type, compute the within-cell-type bidirectional
correlation using self-kernel matrices. This function computes spatial autocorrelation
within each cell type using transferred scores, complementing the cross-type analysis 
provided by \code{getTransferBidirCorr()}.
}
\details{
The self-bidirectional correlation is computed as the mean of two correlations:
cor(t(K) \%*\% A_w, A_w) and cor(A_w, K \%*\% A_w), where A_w is the cell score
vector for a cell type and K is the self-kernel matrix for that cell type.
}
\examples{
\dontrun{
# Assuming you have a CoPro object with multiple cell types
# First compute standard workflow
object <- computeDistance(object)
object <- computeKernelMatrix(object, sigmaValues = c(0.01, 0.05, 0.1))

# Add self-distances and self-kernels
object <- computeSelfDistance(object)
object <- computeSelfKernel(object, sigmaValues = c(0.01, 0.05, 0.1))

# Compute transferred cell scores
trans_scores <- getTransferCellScores(ref_obj, tar_obj, sigma_choice = 0.05, 
                                     agg_cell_type = FALSE)

# Compute self-bidirectional correlation from transferred scores
self_bidir <- getTransferSelfBidirCorr(tar_obj, trans_scores, sigma_choice = 0.05)
}
}
