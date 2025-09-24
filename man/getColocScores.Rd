\name{getColocScores}
\alias{getColocScores}
\title{Compute Colocalization Scores for All Cell Type Pairs}
\usage{
getColocScores(
  object,
  r_um_range = c(10, 60),
  pixel_size_um = 1,
  cell_diam_um = 10,
  nsim = 199,
  r_step_um = 2,
  min_points_per_type = 20,
  edge_correction = "translation",
  verbose = TRUE,
  include_self = FALSE
)
}
\arguments{
\item{object}{A \code{CoProSingle} or \code{CoProMulti} object with location data and cell types.}

\item{r_um_range}{Numeric vector of length 2 specifying the distance band in microns
over which to summarize colocalization (default: c(10, 60)).}

\item{pixel_size_um}{Numeric scalar specifying microns per coordinate unit. Use this
to convert coordinate units to microns if needed (default: 1, assumes coordinates
are already in microns).}

\item{cell_diam_um}{Numeric scalar specifying typical cell diameter in microns, used
for intensity smoothing bandwidth (default: 10).}

\item{nsim}{Integer specifying number of random-labelling simulations for null
distribution (default: 199).}

\item{r_step_um}{Numeric scalar specifying step size in microns for the r grid
(default: 2).}

\item{min_points_per_type}{Integer specifying minimum number of points required per
cell type to compute colocalization (default: 20).}

\item{edge_correction}{Character specifying edge correction method for pair correlation
function (default: "translation").}

\item{verbose}{Logical specifying whether to print progress messages (default: TRUE).}

\item{include_self}{Logical specifying whether to include self-type pairs (same cell type).
Default FALSE since colocalization is typically measured between different cell types.}
}
\value{
For \code{CoProSingle} objects: A data.frame with columns:
\itemize{
\item \code{cellType1}: First cell type name
\item \code{cellType2}: Second cell type name  
\item \code{colocScore}: Colocalization score (higher = more colocalized)
\item \code{nCells1}: Number of cells of type 1
\item \code{nCells2}: Number of cells of type 2
\item \code{r_um_min}: Minimum distance in analysis band
\item \code{r_um_max}: Maximum distance in analysis band
}

For \code{CoProMulti} objects: Same columns plus:
\itemize{
\item \code{slideID}: Slide identifier
}
}
\description{
This function computes cross-type colocalization scores for all pairs of cell types
in a CoPro object using the inhomogeneous cross pair-correlation function approach.
The colocalization score is a normalized measure based on the standardized difference
between observed and expected (under random labelling) spatial correlation.
}
\details{
The colocalization score is computed as:
CES_z = mean_r((g12_obs(r) - mean_null(r)) / sd_null(r))

Where:
\itemize{
\item g12_obs(r) is the observed inhomogeneous cross pair-correlation function
\item mean_null(r) and sd_null(r) are the mean and standard deviation from random labelling simulations
\item The average is taken over the distance band [r_um_range[1], r_um_range[2]]
}

Positive scores indicate colocalization, negative scores indicate segregation.
}
\examples{
\dontrun{
# Basic usage with default parameters
coloc_results <- getColocScores(object)

# Custom parameters for high-resolution data
coloc_results <- getColocScores(object, 
                               r_um_range = c(5, 30),
                               pixel_size_um = 0.325,  # Convert from pixels to microns
                               cell_diam_um = 8,
                               nsim = 99)
                               
# For multi-slide objects
coloc_results <- getColocScores(multi_object)
# Results will include slideID column
}
}
