# Get normalized correlation vs Sigma squared values

Get a data.frame with normalized correlation vs Sigma squared values.
This helps to evaluate which sigma squared value to choose for
downstream analyses.

## Usage

``` r
getNormCorr(object)

# S4 method for class 'CoProSingle'
getNormCorr(object)

# S4 method for class 'CoProMulti'
getNormCorr(object)
```

## Arguments

- object:

  A `CoPro` object

## Value

A `data.frame` with correlation information
