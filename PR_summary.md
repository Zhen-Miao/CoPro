# PR Summary: Auto-K Blocking and BPCells PCA Route

## 1. Auto-K Blocking (`CreateCoPro`)

### Usage

Use `CreateCoPro()` instead of `newCoProSingle()` and `newCoProMulti()`. 
It automatically decides which object type to create and handles large slides
by spatially partitioning them into sub-blocks via k-means.


### Function signature

```r
CreateCoPro(
  normalizedData,
  locationData,
  metaData,
  cellTypes,
  slideID  = NULL,      # if NULL, treated as a single slice
  maxCell  = 50000L,    # maximum cells allowed per (sub)slice
  plot     = FALSE,     # if TRUE, plot the auto-slicing diagnostic plot
  seed     = 1L         # RNG seed for reproducible k-means splitting
)
```

### Behavior

| Scenario | Result |
|---|---|
| Single slide, `ncell <= maxCell` | Returns a `CoProSingle` object |
| Single slide, `ncell > maxCell` | Returns a `CoProMulti` with artificially created spatial blocks as `slideID` |
| Multiple slides, all `<= maxCell` | Returns a `CoProMulti` with original `slideID` |
| Multiple slides, some `> maxCell` | Returns a `CoProMulti` with refined `slideID` (original ID + block suffix, e.g. `"slice1_blk2"`) |


---

## 2. BPCells PCA Route (`computePCA`)

### Usage

Support for [BPCells](https://github.com/bnprks/BPCells) on-disk / out-of-memory matrices
as input to `computePCA()`. When `normalizedData` is a BPCells object, the PCA-related
calculations are replaced by BPCell function calls automatically. 

When `normalizedData` is a regular R matrix object, it uses `irlba` as before.

```r
d = readRDS("some/seurat/obj.rds")
e = as(t(GetAssayData(d, layer="counts")), "IterableMatrix")
m = d@meta.data
obj <- CreateCoPro(
  normalizedData = e,  # cells x genes
  locationData = m[,c("x", "y")],        # data.frame with x, y columns
  metaData = m,                 # data.frame with cell annotations
  cellTypes = m$Tier1,
  slideID = m$Slice_ID # set if there are multiple slides
)
```


### New package dependency

`BPCells` must be added to `Imports` or `Suggests` in `DESCRIPTION` depending on whether
BPCells support should be optional or required. The functions used are:
`colVars`, `binarize`, `add_cols`, `multiply_cols`, `colSums`, `svds`.
