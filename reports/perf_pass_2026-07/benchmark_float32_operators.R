# A/B benchmark for the float32 sparse operators.
# Run from inside a package checkout: Rscript bench_operators.R <label>
suppressMessages(devtools::load_all(".", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
label <- if (length(args)) args[[1L]] else "run"

make_kernel <- function(n_a, n_b, sigma, symmetric, seed = 7) {
  set.seed(seed)
  side <- 1000
  A <- cbind(runif(n_a) * side, runif(n_a) * side)
  B <- if (symmetric) A else cbind(runif(n_b) * side, runif(n_b) * side)
  built <- float32_csr_gaussian_kernels_cpp(
    A, B, sigmas = sigma, percentile = 0.5, scaling_factor = 1,
    lower_limit = 1e-7, upper_quantile = 0.85,
    truncate_low_distance = TRUE, symmetric = symmetric, normalization = 0L
  )
  CoPro:::.newFloat32SparseKernel(built$kernels[[1]])
}

timeit <- function(expr, reps) {
  force(reps)
  eval(expr); # warm
  min(replicate(3, system.time(for (i in seq_len(reps)) eval(expr))[["elapsed"]])) / reps
}

rows <- list()
for (case in list(
  list(tag = "symmetric_1type", n_a = 120000, n_b = 120000, sym = TRUE,  sigma = 8),
  list(tag = "crosstype_2type", n_a = 100000, n_b = 100000, sym = FALSE, sigma = 6)
)) {
  K <- make_kernel(case$n_a, case$n_b, case$sigma, case$sym)
  nnz <- length(K$j)
  for (p in c(10L, 30L)) {
    set.seed(3)
    XL <- matrix(rnorm(case$n_a * p), case$n_a, p)
    XR <- if (case$sym) XL else matrix(rnorm(case$n_b * p), case$n_b, p)
    for (nt in c(1L, 4L)) {
      t_xky <- timeit(quote(CoPro:::.kernelXKY(XL, K, XR, n_threads = nt)), 3)
      rows[[length(rows) + 1L]] <- data.frame(
        label = label, case = case$tag, n = case$n_a, nnz_per_row = round(nnz / case$n_a),
        nPC = p, threads = nt, op = "kernelXKY", seconds = t_xky)
    }
    t_mm <- timeit(quote(CoPro:::.float32KernelMatMult(K, XR, n_threads = 4L)), 3)
    rows[[length(rows) + 1L]] <- data.frame(
      label = label, case = case$tag, n = case$n_a, nnz_per_row = round(nnz / case$n_a),
      nPC = p, threads = 4L, op = "matMult", seconds = t_mm)
  }
  rm(K); invisible(gc(FALSE))
}

out <- do.call(rbind, rows)
print(out, row.names = FALSE)
write.csv(out, file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
                         paste0("operators_", label, ".csv")), row.names = FALSE)
