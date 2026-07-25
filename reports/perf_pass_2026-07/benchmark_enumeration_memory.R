# Peak memory of the parallel neighbour-enumeration path.
#
# Two instruments, because neither alone is honest:
#
#   * `temporary_bytes` from the builder reports the C++ high-water mark of the
#     edge array plus the private per-thread buffers that feed it. R's gc()
#     cannot see any of this -- it counts Vcells only -- so every "peak MB"
#     figure elsewhere in this folder is blind to it.
#   * peak RSS, read from ps in a sampling loop, catches everything the process
#     touches including the allocator's own behaviour.
#
# Usage:
#   Rscript reports/perf_pass_2026-07/benchmark_enumeration_memory.R [outdir]

suppressMessages(devtools::load_all(".", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[[1L]] else "reports/perf_pass_2026-07"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

N <- 120000L
set.seed(4242)
coords <- cbind(runif(N) * 1000, runif(N) * 1000)

# Sample this process's RSS from a detached shell while the build runs.
with_peak_rss <- function(expr) {
  pid <- Sys.getpid()
  sample_file <- tempfile(fileext = ".txt")
  script <- sprintf(
    "while kill -0 %d 2>/dev/null; do ps -o rss= -p %d >> %s; sleep 0.05; done",
    pid, pid, shQuote(sample_file))
  handle <- system2("bash", c("-c", shQuote(script)), wait = FALSE,
                    stdout = NULL, stderr = NULL)
  on.exit(unlink(sample_file), add = TRUE)
  before <- as.numeric(system2("ps", c("-o", "rss=", "-p", pid), stdout = TRUE))
  value <- force(expr)
  Sys.sleep(0.2)
  samples <- suppressWarnings(as.numeric(readLines(sample_file, warn = FALSE)))
  samples <- samples[is.finite(samples)]
  list(value = value,
       rss_before_mb = before / 1024,
       rss_peak_mb = if (length(samples)) max(samples) / 1024 else NA_real_)
}

rows <- list()
for (threads in c(1L, 4L, 8L)) {
  gc(reset = TRUE, full = TRUE)
  timing <- system.time({
    measured <- with_peak_rss(
      float32_csr_gaussian_kernels_cpp(
        coords, coords, sigmas = c(20), percentile = 0.3, scaling_factor = 1,
        lower_limit = 1e-7, upper_quantile = 0.99,
        truncate_low_distance = TRUE, symmetric = TRUE,
        normalization = 0L, n_threads = threads))
  })
  built <- measured$value
  rows[[length(rows) + 1L]] <- data.frame(
    threads       = threads,
    seconds       = as.numeric(timing[["elapsed"]]),
    edges         = built$candidate_pairs,
    cpp_peak_mb   = built$temporary_bytes / 1024^2,
    rss_peak_mb   = measured$rss_peak_mb,
    rss_delta_mb  = measured$rss_peak_mb - measured$rss_before_mb,
    gc_vcells_mb  = sum(gc()[, "max used"] * c(0, 8)) / 1024^2
  )
  cat(sprintf("threads %d  %.2f s  C++ peak %.0f MB  RSS peak %.0f MB\n",
              threads, timing[["elapsed"]],
              built$temporary_bytes / 1024^2, measured$rss_peak_mb))
  rm(built, measured); gc(full = TRUE)
}

out <- do.call(rbind, rows)
write.csv(out, file.path(outdir, "enumeration_memory.csv"), row.names = FALSE)
print(out)
