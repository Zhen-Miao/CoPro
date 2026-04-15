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
