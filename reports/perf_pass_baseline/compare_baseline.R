# Diff a freshly captured run against the pre-change baseline.
#
# Usage:
#   Rscript reports/perf_pass_baseline/compare_baseline.R [tolerance]
#
# tolerance = 0 (default) demands bit-identical results, which is the bar for
# plan items 1, 3, 6, 7 and 8. Items 2, 4 and 5 change floating-point
# summation order, so run them with a tolerance (1e-10).

args <- commandArgs(trailingOnly = TRUE)
tolerance <- if (length(args) >= 1L) as.numeric(args[[1L]]) else 0

base_dir <- "reports/perf_pass_baseline"
baseline <- readRDS(file.path(base_dir, "baseline.rds"))

tmp <- file.path(base_dir, "current")
system2("Rscript", c(file.path(base_dir, "capture_baseline.R"), tmp),
        stdout = FALSE, stderr = FALSE)
current <- readRDS(file.path(tmp, "baseline.rds"))

stopifnot(identical(names(baseline), names(current)))

failures <- character()
for (scenario in names(baseline)) {
  cmp <- if (tolerance == 0) {
    if (identical(baseline[[scenario]], current[[scenario]])) TRUE else
      "not bit-identical"
  } else {
    all.equal(baseline[[scenario]], current[[scenario]],
              tolerance = tolerance)
  }
  if (isTRUE(cmp)) {
    cat(sprintf("  OK   %s\n", scenario))
  } else {
    failures <- c(failures, scenario)
    cat(sprintf("  FAIL %s\n", scenario))
    cat(paste0("         ", utils::head(as.character(cmp), 5), "\n"))
  }
}

cat("\n")
if (length(failures) == 0L) {
  cat(sprintf("All %d scenarios match (tolerance = %g).\n",
              length(baseline), tolerance))
} else {
  cat(sprintf("%d of %d scenarios differ: %s\n", length(failures),
              length(baseline), paste(failures, collapse = ", ")))
  quit(status = 1L)
}
