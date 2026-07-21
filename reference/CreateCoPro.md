# Create a CoPro object, automatically choosing Single vs Multi and splitting large slices.

Behavior:

1.  one slice and ncell \< p -\> newCoProSingle()

2.  one slice and ncell \> p -\> newCoProMulti() with artificial blocks

3.  multi slices, all \< p -\> newCoProMulti() with original slideID

4.  multi slices, some \> p -\> newCoProMulti() with refined slideID
    (origin + block)

## Usage

``` r
CreateCoPro(
  normalizedData,
  locationData,
  metaData,
  cellTypes,
  slideID = NULL,
  maxCell = 50000L,
  plot = FALSE,
  seed = 1L
)
```

## Arguments

- normalizedData:

  matrix (cells x genes)

- locationData:

  data.frame (cells x coords; needs x,y, optionally z)

- metaData:

  data.frame (cells x annotations)

- cellTypes:

  vector length nrow(normalizedData)

- slideID:

  optional vector length nrow(normalizedData); if NULL, treated as one
  slice

- maxCell:

  integer max cells per (sub)slice (default 50000)

- plot:

  logical; if TRUE, plot the auto-slicing result (default FALSE)

- seed:

  integer seed for reproducible splitting (default 1)

## Note on spatial binning (approximate)

Splitting an oversized slide into spatial blocks (cases 2 and 4) treats
each block as an independent pseudo-slide. Because kernels are computed
only within a block, this is a **block-diagonal approximation** of the
full spatial kernel: cell pairs that straddle a block boundary are
dropped. The approximation is good when a block's linear size greatly
exceeds the kernel support radius
\\\sigma\sqrt{-2\log(\texttt{lowerLimit})}\\, but the error is not
bounded in general and grows with `sigma`. Two further caveats: the
normalized- and bidirectional-correlation metrics average per-slide
quantities, so the block count affects their normalization; and rare
cell types may be silently dropped from blocks where their count falls
to 5 or fewer. For an **exact** large-scale kernel that needs no
binning, use `computeKernelMatrix(..., method = "sparse")` (see
[`computeSparseKernel()`](https://zhen-miao.github.io/CoPro/reference/computeSparseKernel.md)),
which builds a single sparse kernel per slide via a fixed-radius
neighbor search.
