\name{plotG12Functions}
\alias{plotG12Functions}
\title{Plot g_12(r) pair correlation functions for colocalization analysis}
\usage{
plotG12Functions(
  object,
  r_um_range = c(10, 60),
  pixel_size_um = 1,
  cell_diam_um = 10,
  nsim = 199,
  r_step_um = 2,
  min_points_per_type = 20,
  edge_correction = "translation",
  plot_type = c("combined", "individual", "both"),
  include_confidence = TRUE,
  confidence_level = 0.95,
  colors = NULL,
  line_width = 1.2,
  alpha_bands = 0.3,
  verbose = TRUE
)
}
\arguments{
\item{object}{A CoProSingle or CoProMulti object with location and cell type data.}

\item{r_um_range}{Numeric vector of length 2 specifying the distance range in 
microns for analysis (default: c(10, 60)).}

\item{pixel_size_um}{Numeric value specifying the size of each pixel/coordinate 
unit in microns (default: 1).}

\item{cell_diam_um}{Numeric value specifying typical cell diameter in microns 
for intensity smoothing (default: 10).}

\item{nsim}{Integer specifying number of random labelling simulations for 
null distribution (default: 199).}

\item{r_step_um}{Numeric value specifying step size for radius grid in microns 
(default: 2).}

\item{min_points_per_type}{Integer specifying minimum number of points required 
per cell type (default: 20).}

\item{edge_correction}{Character specifying edge correction method 
(default: "translation").}

\item{plot_type}{Character specifying plot type: "combined" for all pairs in 
one plot, "individual" for separate plots per pair, or "both" (default: "combined").}

\item{include_confidence}{Logical; whether to include confidence bands from 
simulations (default: TRUE).}

\item{confidence_level}{Numeric; confidence level for bands (default: 0.95).}

\item{colors}{Character vector of colors for different cell type pairs. If NULL,
uses default color palette.}

\item{line_width}{Numeric; width of the g_12(r) lines (default: 1.2).}

\item{alpha_bands}{Numeric; transparency for confidence bands (default: 0.3).}

\item{verbose}{Logical; whether to print progress messages (default: TRUE).}
}
\value{
A list containing:
\itemize{
\item \code{plot}: ggplot object(s) with the g_12(r) plots
\item \code{data}: data.frame with g_12(r) values, confidence intervals, and metadata
\item \code{summary}: summary statistics for each cell type pair
}
}
\description{
This function creates plots of the inhomogeneous cross pair-correlation function
g_12(r) values across different radii for all cell type pairs. It can generate
both individual plots for each pair and a combined plot with all pairs overlaid.
}
\details{
The function computes the inhomogeneous cross pair-correlation function
g_12(r) for each cell type pair and plots the observed values along with 
confidence bands derived from random labelling simulations.

Values of g_12(r):
\itemize{
\item g_12(r) = 1: Random spatial relationship at distance r
\item g_12(r) > 1: Attraction/colocalization at distance r  
\item g_12(r) < 1: Repulsion/segregation at distance r
}
}
\examples{
\dontrun{
# Basic usage
g12_plots <- plotG12Functions(object)

# Custom parameters with individual plots
g12_plots <- plotG12Functions(object, 
                             r_um_range = c(5, 40),
                             plot_type = "individual",
                             include_confidence = TRUE)

# Access the plot and data
print(g12_plots$plot)
head(g12_plots$data)
}
}
