# Internal: partition locations into approximately equal-sized spatial blocks

Uses k-means on coordinates x, y, (and z) and recursively subdivides any
cluster that still exceeds maxCell.

## Usage

``` r
.partitionByLocation(locationData, n, maxCell)
```

## Arguments

- locationData:

  data.frame with columns x,y,(optional z)

- n:

  target number of blocks (\>= 1)

- maxCell:

  maximum cells per block (\>= 1)

## Value

character vector of block labels length nrow(locationData)
