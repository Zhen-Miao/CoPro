# Whitened-Frobenius normalizer: true A/B on the shipped function.
#
# An earlier version of this script timed a benchmark-local rewrite of the old
# code, `sum(((Rx %*% K) %*% Ry) * K)`, against the new inner product. That is
# not the pre-pass code path: .whitenedFrobNorm() also computes the low-rank
# cross terms of the double-centering expansion, which the rewrite omitted, so
# the reported speedup credited this change with work the old path never did.
#
# Instead, run the *shipped* .whitenedFrobNorm() on identical input in two
# checkouts and compare, exactly as benchmark_float32_operators.R does:
#
#   git worktree add <dir> <pre-pass commit>
#   cd <dir> && Rscript reports/perf_pass_2026-07/benchmark_whitened_frobenius.R before <out>
#   cd <repo> && Rscript reports/perf_pass_2026-07/benchmark_whitened_frobenius.R after  <out>
#
# The block-size sweep at the end is a parameter study of the new helper and
# only runs on the "after" arm.
#
# Every timed arm is repeated; single-shot timings on this machine vary by up
# to 1.9x (see the variance section of RESULTS.md).

suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
arm <- if (length(args) >= 1L) args[[1L]] else "after"
outdir <- if (length(args) >= 2L) args[[2L]] else "reports/perf_pass_2026-07"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

REPS <- 5L

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

# Peak is measured as R Vcells above the resting set. That is the right
# instrument here: every allocation in this path is an R-level sparse matrix.
timed <- function(fun, reps = REPS) {
  fun()  # warm up; the first call in a session is not representative
  seconds <- numeric(reps)
  peaks <- numeric(reps)
  value <- NULL
  for (i in seq_len(reps)) {
    invisible(gc(reset = TRUE, full = TRUE))
    before <- gc()[2, "max used"]
    seconds[i] <- system.time(value <- fun())[["elapsed"]]
    peaks[i] <- (gc()[2, "max used"] - before) * 8 / 1e6
  }
  list(value = value, min = min(seconds), median = stats::median(seconds),
       max = max(seconds), mb = max(peaks), all = seconds)
}

cat(sprintf("arm = %s, n = %d, nnz(K) = %.2f M, reps = %d\n",
            arm, n, length(K@x) / 1e6, REPS))

whole <- timed(function() CoPro:::.whitenedFrobNorm(K, Rxs, Rys))
cat(sprintf("  .whitenedFrobNorm()  : min %5.2f s  median %5.2f s  peak %5.0f MB\n",
            whole$min, whole$median, whole$mb))
cat(sprintf("    value  %.17g\n", as.numeric(whole$value)))
cat(sprintf("    all    [%s]\n",
            paste(sprintf("%.2f", whole$all), collapse = " ")))

write.csv(
  data.frame(arm = arm, n = n, nnz = length(K@x), reps = REPS,
             seconds_min = whole$min, seconds_median = whole$median,
             seconds_max = whole$max, peak_mb = whole$mb,
             value = as.numeric(whole$value)),
  file.path(outdir, paste0("frobenius_", arm, ".csv")), row.names = FALSE)

if (arm != "after" || !exists(".sparseWhitenedInner", envir = asNamespace("CoPro"))) {
  quit(save = "no")
}

# Parameter study: how the block budget trades time against peak. The reference
# is the same computation with no blocking at all.
cat("\n  block-size sweep on .sparseWhitenedInner()\n")
ref <- timed(function() {
  M <- (Rxs %*% K) %*% Rys
  sum(M * K)
})
cat(sprintf("    unblocked           : min %5.2f s  median %5.2f s  peak %5.0f MB\n",
            ref$min, ref$median, ref$mb))

rows <- list(data.frame(block_nnz = NA_real_, blocks = 1L,
                        seconds_min = ref$min, seconds_median = ref$median,
                        peak_mb = ref$mb, rel_diff = 0))
for (bn in c(2e6, 1e6, 2e5, 5e4, 2e4)) {
  budget <- bn
  r <- timed(function()
    CoPro:::.sparseWhitenedInner(K, Rxs, Rys, block_nnz = budget))
  per_col <- max(1, length(K@x) / n)
  blk <- max(1L, min(n, as.integer(budget / per_col)))
  rel <- abs(r$value - as.numeric(ref$value)) / abs(as.numeric(ref$value))
  cat(sprintf(
    "    block_nnz = %7.0f : min %5.2f s  median %5.2f s  peak %5.0f MB  %4d blocks  rel diff %.2g\n",
    budget, r$min, r$median, r$mb, ceiling(n / blk), rel))
  rows[[length(rows) + 1L]] <- data.frame(
    block_nnz = budget, blocks = ceiling(n / blk),
    seconds_min = r$min, seconds_median = r$median, peak_mb = r$mb,
    rel_diff = rel)
}
write.csv(do.call(rbind, rows),
          file.path(outdir, "frobenius_block_sweep.csv"), row.names = FALSE)
