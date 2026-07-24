# Materialize all encoded kernels in a CoPro object

Converts every `CoProFloat32SparseMatrix` in `object@kernelMatrices` to
a standard double-precision sparse `Matrix`. This is a compatibility
escape hatch for third-party or older code that accesses the slot
directly rather than using
[`getKernelMatrix()`](https://zhen-miao.github.io/CoPro/reference/getKernelMatrix.md).
It can temporarily require both encoded and decoded copies, so it should
be avoided in memory-critical analysis steps.

## Usage

``` r
materializeFloat32Kernels(object, verbose = TRUE)
```

## Arguments

- object:

  A `CoPro` object.

- verbose:

  Whether to report how many kernels were converted.

## Value

A copy of `object` whose kernel list contains standard sparse matrices.
