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

## See also

[`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md),
[`CreateCoPro()`](https://zhen-miao.github.io/CoPro/reference/CreateCoPro.md),
[`asCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/asCoPro.md)

Other object-creation:
[`asCoPro`](https://zhen-miao.github.io/CoPro/reference/asCoPro.md),
[`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md),
[`subsetData()`](https://zhen-miao.github.io/CoPro/reference/subsetData.md)

## Examples

``` r
toy <- readRDS(system.file("extdata", "toy_copro_data.rds", package = "CoPro"))
obj <- newCoProSingle(
  normalizedData = toy$normalizedData,
  locationData   = toy$locationData,
  metaData       = toy$metaData,
  cellTypes      = toy$cellTypes
)
obj
#> 'CoProSingle' object for spatial coordinated progression detection
#> ------------------------
#> Number of cells: 200
#> Number of genes: 80
#> Number of slides: 1 (single-slide analysis)
#> 
#> Processing steps completed:
#> 
#> Available metadata fields:
#> - cell_id
#> - cell_type
#> 
#> Approx. object size: 0.20 MB
```
