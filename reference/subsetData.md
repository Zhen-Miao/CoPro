# subsetData

Take a subset of the original matrix based on cell types of interest.
The original data are stored without being thrown away

## Usage

``` r
subsetData(object, cellTypesOfInterest, saveOriginal = FALSE)
```

## Arguments

- object:

  A `CoPro` object

- cellTypesOfInterest:

  Input cell types of interest as a vector of characters for subsetting
  the data

- saveOriginal:

  Logical, whether to save the original data in the subsetted object.
  Default is FALSE.

## Value

A `CoPro` object with subset slots

## See also

[`newCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/newCoProSingle.md),
[`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md),
[`CreateCoPro()`](https://zhen-miao.github.io/CoPro/reference/CreateCoPro.md)

Other object-creation:
[`asCoPro`](https://zhen-miao.github.io/CoPro/reference/asCoPro.md),
[`newCoProMulti()`](https://zhen-miao.github.io/CoPro/reference/newCoProMulti.md),
[`newCoProSingle()`](https://zhen-miao.github.io/CoPro/reference/newCoProSingle.md)
