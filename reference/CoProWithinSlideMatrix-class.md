# Matrix-free within-slide standardized expression matrix

Internal representation of \$\$Z_s = (X_s -
1\mu_s')\operatorname{diag}(d_s)^{-1}\$\$ stacked in the original row
order. `irlba()` only needs right and left matrix products, so sparse
inputs can stay sparse even though explicit centering would make every
zero nonzero.

## Usage

``` r
# S4 method for class 'CoProWithinSlideMatrix'
dim(x)

# S4 method for class 'CoProWithinSlideMatrix,numeric'
x %*% y

# S4 method for class 'CoProWithinSlideMatrix,matrix'
x %*% y

# S4 method for class 'numeric,CoProWithinSlideMatrix'
x %*% y

# S4 method for class 'matrix,CoProWithinSlideMatrix'
x %*% y
```

## Arguments

- x:

  A `CoProWithinSlideMatrix` object or conformable numeric operand.

- y:

  A `CoProWithinSlideMatrix` object or conformable numeric operand.

## Value

Matrix dimensions for [`dim()`](https://rdrr.io/r/base/dim.html); a
dense low-dimensional product for `%*%`.

## Slots

- `blocks`:

  Per-slide expression blocks in their original storage class.

- `rows`:

  Original row indices for each block.

- `centers`:

  Matrix of per-block gene means.

- `scales`:

  Matrix of per-block gene scales.

- `dims`:

  Integer vector containing the stacked matrix dimensions.
