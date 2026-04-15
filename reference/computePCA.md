# Compute PCA on Single-Slide Data

This method performs PCA on the normalized data stored within the
`CoProSingle` object. It assumes that the data has already been
integrated across slides.

This method performs PCA on the normalized data stored within the
`CoProMulti` object. It assumes that the data has already been
integrated across slides.

## Usage

``` r
# S4 method for class 'CoProSingle'
computePCA(object, nPCA = 40, center = TRUE, scale. = TRUE, scalePCs = TRUE)

# S4 method for class 'CoProMulti'
computePCA(
  object,
  nPCA = 40,
  center = TRUE,
  scale. = TRUE,
  scalePCs = TRUE,
  dataUse = "raw",
  center_per_slide = FALSE
)
```

## Arguments

- object:

  A `CoProMulti` object with the `normalizedData` slot populated.

- nPCA:

  Number of principal components to compute for each cell type.

- center:

  Whether to center the matrix before PCA

- scale.:

  Whether to scale the matrix before PCA

- scalePCs:

  Whether to scale (whiten) PCs by their standard deviation before
  downstream CCA optimization. Default `TRUE` (recommended).

- dataUse:

  What data to use, choices between "raw" and "integrated". Default is
  "raw". For single slide, this argument is ignored.

- center_per_slide:

  After the global PCA, do we do center per slide again? By default this
  is set to FALSE
