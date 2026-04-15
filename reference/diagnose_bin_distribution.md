# Diagnose Bin Distribution

Helper function to check how cells are distributed across bins. Useful
for choosing appropriate bin numbers.

## Usage

``` r
diagnose_bin_distribution(location_data, num_bins_x = 10, num_bins_y = 10)
```

## Arguments

- location_data:

  Data frame with x, y columns

- num_bins_x:

  Number of bins in x direction

- num_bins_y:

  Number of bins in y direction

## Value

List with bin statistics
