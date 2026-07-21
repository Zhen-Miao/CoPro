# Whitened-Frobenius null SD of the bilinear statistic a' K b

Internal. Returns \\\\R_x^{1/2} K_c R_y^{1/2}\\\_F\\, the
distribution-free null standard deviation of \\T = a' K b\\, where
\\K_c\\ is the double-centered cross-kernel and \\R_x, R_y\\ are the
within-type correlation operators. This replaces the spectral norm
\\\\K\\\_2\\, which is scale-blind and rails bandwidth selection to the
noise floor. With `Rx`/`Ry` omitted it degrades to the un-whitened
\\\\K_c\\\_F\\ (i.e. \\R_x = R_y = I\\).

## Usage

``` r
.whitenedFrobNorm(K, Rx = NULL, Ry = NULL)
```

## Arguments

- K:

  Cross-type kernel matrix (double-centered internally).

- Rx, Ry:

  Within-type kernels used as correlation operators; symmetrized
  defensively (a no-op when the kernel is un-normalized). `NULL` to skip
  whitening.

## Value

Scalar whitened-Frobenius norm.
