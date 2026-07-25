suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(Matrix))

set.seed(6)
n <- 150000
k <- 30
mk <- function(nr, nc, kk) sparseMatrix(
  i = rep(seq_len(nr), each = kk),
  j = pmax(1L, pmin(nc, rep(seq_len(nr), each = kk) +
                      sample(-50:50, nr * kk, TRUE))),
  x = runif(nr * kk), dims = c(nr, nc))
K <- mk(n, n, k)
Rxs <- local({ R <- mk(n, n, k); (R + t(R)) / 2 })
Rys <- local({ R <- mk(n, n, k); (R + t(R)) / 2 })

peak_of <- function(fun) {
  invisible(gc(reset = TRUE, full = TRUE))
  before <- gc()[2, "max used"]
  elapsed <- system.time(value <- fun())[["elapsed"]]
  mb <- (gc()[2, "max used"] - before) * 8 / 1e6
  list(value = value, seconds = elapsed, mb = mb)
}

cat(sprintf("n = %d, nnz(K) = %.2f M\n", n, length(K@x) / 1e6))

ref <- peak_of(function() {
  M <- (Rxs %*% K) %*% Rys
  s <- sum(M * K)
  attr(s, "fill") <- length(M@x) / length(K@x)
  s
})
cat(sprintf("  unblocked (Rx K) Ry : %5.2f s   peak %5.0f MB   fill %.0fx\n",
            ref$seconds, ref$mb, attr(ref$value, "fill")))
invisible(gc(FALSE))

for (bn in c(2e6, 1e6, 2e5, 5e4, 2e4)) {
  local({
    budget <- bn
    r <- peak_of(function()
      CoPro:::.sparseWhitenedInner(K, Rxs, Rys, block_nnz = budget))
    per_col <- max(1, length(K@x) / n)
    blk <- max(1L, min(n, as.integer(budget / per_col)))
    cat(sprintf(
      "  block_nnz = %7.0f   : %5.2f s   peak %5.0f MB   %4d blocks   rel diff %.2g\n",
      budget, r$seconds, r$mb, ceiling(n / blk),
      abs(r$value - as.numeric(ref$value)) / abs(as.numeric(ref$value))))
    invisible(gc(FALSE))
  })
}
