\name{computeSelfBidirCorr}
\alias{computeSelfBidirCorr}
\title{Compute Self-Bidirectional Correlation using skrCCA Results}
\usage{
computeSelfBidirCorr(
  object,
  sigma_choice,
  calculationMode = "perSlide",
  normalize_K = c("row_or_col", "sinkhorn_knopp", "none"),
  filter_kernel = TRUE,
  K_row_sum_cutoff = 0.005,
  K_col_sum_cutoff = 0.005,
  verbose = TRUE
)
}
\arguments{
\item{object}{A \code{CoProSingle} or \code{CoProMulti} object with skrCCA results and
self-kernel matrices computed using \code{computeSelfKernel()}.}

\item{sigma_choice}{Numeric scalar specifying the sigma value to use.}

\item{calculationMode}{For \code{CoProMulti} objects only, either "perSlide" or
"aggregate". Default "perSlide".}

\item{normalize_K}{Character; method for normalizing the kernel matrix, one of 
"row_or_col", "sinkhorn_knopp", or "none". Default "row_or_col".}

\item{filter_kernel}{Logical; whether to filter the kernel matrix. Default TRUE.}

\item{K_row_sum_cutoff}{Numeric; cutoff for row sums when normalizing kernel matrix.
Default 5e-3.}

\item{K_col_sum_cutoff}{Numeric; cutoff for column sums when normalizing kernel matrix.
Default 5e-3.}

\item{verbose}{Logical; whether to print progress messages.}
}
\value{
A list with one element named \code{paste0("sigma_", sigma_choice)}, whose
value is a data.frame of results with the same structure as \code{getSelfBidirCorr()}.
}
\description{
This function computes self-bidirectional correlation directly from a CoPro object
that has skrCCA results, using the object's own cell scores rather than transferred
scores. This is useful for computing spatial autocorrelation patterns within each
cell type using the object's native skrCCA results.
}
\examples{
\dontrun{
# Assuming you have a CoPro object with skrCCA results and self-kernels
object <- runSkrCCA(object)
object <- computeSelfDistance(object)
object <- computeSelfKernel(object, sigmaValues = c(0.01, 0.05, 0.1))

# Compute self-bidirectional correlation using native skrCCA results
self_bidir <- computeSelfBidirCorr(object, sigma_choice = 0.05)
}
}
