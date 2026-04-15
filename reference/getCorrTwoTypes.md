# Retrieve the Correlation between two cell types

Retrieve the Correlation between two cell types

## Usage

``` r
getCorrTwoTypes(
  object,
  cellTypeA,
  cellTypeB,
  ccIndex = 1,
  sigmaValueChoice = NULL
)
```

## Arguments

- object:

  A `CoPro` object

- cellTypeA:

  Cell type label for one cell type

- cellTypeB:

  Cell type label for another cell type

- ccIndex:

  Canonical vector index, default = 1

- sigmaValueChoice:

  A particular sigma squared value for the correlation

## Value

A data.frame with two columns, AK and B, where AK represents the cell
score of cell type A times the kernel matrix, and B represents the cell
score of cell type B. If the object is a `CoProMulti` object, the
data.frame will have an additional column, slideID, indicating the slide
ID.
