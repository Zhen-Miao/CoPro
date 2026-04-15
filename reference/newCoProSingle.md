# Function to create a new object

Function to create a new object

## Usage

``` r
newCoProSingle(normalizedData, locationData, metaData, cellTypes)

# S4 method for class 'ANY,matrixOrDataFrame,data.frame,factorOrCharacter'
newCoProSingle(normalizedData, locationData, metaData, cellTypes)
```

## Arguments

- normalizedData:

  A `matrix` object to store normalized data.

- locationData:

  A `data.frame` object to store the location. It should either contain
  two columns named by "x" and "y", or three columns named by "x", "y",
  and "z". No other names allowed

- metaData:

  A `data.frame` object to store metadata for each cell.

- cellTypes:

  A `vector` object with elements being character. It should match the
  number of cells in the data matrix and each represents a cell type
  label of a cell.

## Value

A `CoPro` object
