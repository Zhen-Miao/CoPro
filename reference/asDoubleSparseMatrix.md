# Decode a float32 kernel into a standard sparse Matrix

This diagnostic conversion materializes float64 values and should not be
used in a large-data computation.

## Usage

``` r
asDoubleSparseMatrix(x)
```

## Arguments

- x:

  A `CoProFloat32SparseMatrix`.

## Value

A `Matrix::dgCMatrix`, or a `Matrix::dsCMatrix` for a symmetric one-type
kernel.
