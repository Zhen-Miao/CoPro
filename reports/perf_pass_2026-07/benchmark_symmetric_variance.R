# Focused re-measurement of the symmetric / nPC=10 / 1-thread cell, which was
# the first timing taken in the sweep and so the most likely to be polluted.
suppressMessages(devtools::load_all(".", quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
label <- if (length(args)) args[[1L]] else "run"

set.seed(7)
n <- 120000
A <- cbind(runif(n) * 1000, runif(n) * 1000)
built <- float32_csr_gaussian_kernels_cpp(
  A, A, sigmas = 8, percentile = 0.5, scaling_factor = 1,
  lower_limit = 1e-7, upper_quantile = 0.85,
  truncate_low_distance = TRUE, symmetric = TRUE, normalization = 0L)
K <- CoPro:::.newFloat32SparseKernel(built$kernels[[1]])
rm(built, A); invisible(gc(FALSE))

for (p in c(10L, 30L)) {
  set.seed(3)
  X <- matrix(rnorm(n * p), n, p)
  invisible(CoPro:::.kernelXKY(X, K, X, n_threads = 1L))  # warm
  t <- replicate(6, system.time(
    CoPro:::.kernelXKY(X, K, X, n_threads = 1L))[["elapsed"]])
  cat(sprintf("%-8s symmetric p=%2d threads=1 : min %.3f  median %.3f  all [%s]\n",
              label, p, min(t), median(t),
              paste(sprintf("%.2f", sort(t)), collapse = " ")))
}
