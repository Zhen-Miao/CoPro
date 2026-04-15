# Create a new CoProMulti object for Multi-Slide Analysis

Initializes a `CoProMulti` object with combined data from multiple
slides.

## Usage

``` r
newCoProMulti(normalizedData, locationData, metaData, cellTypes, slideID)

# S4 method for class 'ANY,matrixOrDataFrame,data.frame,factorOrCharacter,factorOrCharacter'
newCoProMulti(normalizedData, locationData, metaData, cellTypes, slideID)
```

## Arguments

- normalizedData:

  Combined normalized expression matrix (cells x genes) for all slides.
  Rownames should be unique cell identifiers.

- locationData:

  Combined location data frame (cells x coordinates) for all slides.
  Rownames must match `normalizedData`. Columns 'x', 'y', (and
  optionally 'z') required.

- metaData:

  Combined metadata data frame (cells x annotations) for all slides.
  Rownames must match `normalizedData`.

- cellTypes:

  Combined cell type labels vector for all cells. Length must match
  `nrow(normalizedData)`.

- slideID:

  Combined slide/sample identifier vector for all cells. Length must
  match `nrow(normalizedData)`.

## Value

A `CoProMulti` object.
