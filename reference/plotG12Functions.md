# Plot g_12(r) pair correlation functions for colocalization analysis

This function creates plots of the inhomogeneous cross pair-correlation
function g_12(r) values across different radii for all cell type pairs.
It can generate both individual plots for each pair and a combined plot
with all pairs overlaid.

## Usage

``` r
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
```

## Arguments

- object:

  A CoProSingle or CoProMulti object with location and cell type data.

- r_um_range:

  Numeric vector of length 2 specifying the distance range in microns
  for analysis (default: c(10, 60)).

- pixel_size_um:

  Numeric value specifying microns per coordinate unit. This converts
  your input coordinates to microns for analysis. For example:

  - If 1 coordinate unit = 50 microns, set pixel_size_um = 50

  - If 1 coordinate unit = 0.325 microns, set pixel_size_um = 0.325

  - If coordinates are already in microns, set pixel_size_um = 1
    (default)

- cell_diam_um:

  Numeric value specifying typical cell diameter in microns for
  intensity smoothing (default: 10).

- nsim:

  Integer specifying number of random labelling simulations for null
  distribution (default: 199).

- r_step_um:

  Numeric value specifying step size for radius grid in microns
  (default: 2).

- min_points_per_type:

  Integer specifying minimum number of points required per cell type
  (default: 20).

- edge_correction:

  Character specifying edge correction method (default: "translation").

- plot_type:

  Character specifying plot type: "combined" for all pairs in one plot,
  "individual" for separate plots per pair, or "both" (default:
  "combined").

- include_confidence:

  Logical; whether to include confidence bands from simulations
  (default: TRUE).

- confidence_level:

  Numeric; confidence level for bands (default: 0.95).

- colors:

  Character vector of colors for different cell type pairs. If NULL,
  uses default color palette.

- line_width:

  Numeric; width of the g_12(r) lines (default: 1.2).

- alpha_bands:

  Numeric; transparency for confidence bands (default: 0.3).

- verbose:

  Logical; whether to print progress messages (default: TRUE).

## Value

A list containing:

- `plot`: ggplot object(s) with the g_12(r) plots

- `data`: data.frame with g_12(r) values, confidence intervals, and
  metadata

- `summary`: summary statistics for each cell type pair

## Details

The function computes the inhomogeneous cross pair-correlation function
g_12(r) for each cell type pair and plots the observed values along with
confidence bands derived from random labelling simulations.

Values of g_12(r):

- g_12(r) = 1: Random spatial relationship at distance r

- g_12(r) \> 1: Attraction/colocalization at distance r

- g_12(r) \< 1: Repulsion/segregation at distance r

## Examples

``` r
if (FALSE) { # \dontrun{
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
} # }
```
