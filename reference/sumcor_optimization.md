# PCA-space SUMCOR optimization for CoPro

The SUMCOV objective that `optimize_bilinear*()` maximizes,
\\\sum\_{i\<j} w_i' (\sum_s X_i^{(s)'} K\_{ij}^{(s)} X_j^{(s)}) w_j\\
subject to \\\\w_i\\ = 1\\, factors exactly as a *slide-weighted*
SUMCOR:

## Details

\$\$f\_{cov}(w) = \sum\_{i\<j} \sum_s \sigma_i^{(s)} \sigma_j^{(s)}\\
\rho\_{ij}^{(s)}\$\$

with \\\rho\_{ij}^{(s)} = u_i' K u_j / (\\u_i\\ \\u_j\\)\\ and
\\u_i^{(s)} = X_i^{(s)} w_i\\. There is no \\\sqrt{n_i n_j}\\ factor
here: \\\sigma\\ is the norm \\\\X_i^{(s)} w_i\\\\, not a
root-mean-square, so \\\sigma_i \sigma_j \rho\_{ij} = w_i' Y\_{ij} w_j\\
is the SUMCOV term already (see `.sumcorSigma()` and
`.sumcorSlideWeight()`). The norm constraint pins the *pooled* variance,
so per-slide variances stay free and a slide with inflated variance
along the canonical direction gets a proportionally larger vote. That is
the batch-domination mode SUMCOR removes.

Two consequences shape this file. First, the slide weight and the scale
factor are separable: `slideWeight = "size"` reintroduces \\\sqrt{n_i
n_j}\\ so larger slides count for more, while dropping \\\sigma_i
\sigma_j\\ removes the batch-scale sensitivity. Second, for a single
slide the two objectives often – but not always – pose the same problem.
With whitened PCs \\X_i'X_i = (n_i-1) I\\, so on \\\\w_i\\ = 1\\ the
denominators are \\\sigma_i = \sqrt{n_i - 1}\\ and

\$\$f\_{equal}(w) = \sum\_{i\<j} \frac{w_i' Y\_{ij} w_j}{
\sqrt{(n_i-1)(n_j-1)}}\$\$

differs from \\f\_{cov}\\ by a *per-pair* constant. A constant that
varies by pair leaves the argmax alone only when there is a single pair,
or when every \\n_i\\ is equal. So the reduction to SUMCOV is exact for
one or two cell types at any cell counts, and for three or more only
when the counts are equal; `.sumcorReducesToSumcov()` is that test.
Outside it the criteria have genuinely different maximizers, and the
full-gradient optimizer runs rather than short-circuiting. Under
`slideWeight = "size"` the mismatch is \\1 + O(1/n)\\ and usually
immaterial; it is `"equal"`, the multi-slide default with equal nominal
coefficients, where it can be material.
