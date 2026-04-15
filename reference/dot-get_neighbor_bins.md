# Get Neighboring Bins

Internal helper function to find neighboring bins (3x3 grid around
target).

## Usage

``` r
.get_neighbor_bins(bin_coords, tar_bin_id, num_bins_x, num_bins_y)
```

## Arguments

- bin_coords:

  Data frame with bin_id, x_bin, y_bin columns

- tar_bin_id:

  Target bin ID (string like "3_5")

- num_bins_x:

  Number of bins in x direction

- num_bins_y:

  Number of bins in y direction

## Value

Character vector of neighboring bin IDs
