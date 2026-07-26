# Diff a freshly captured run against the pre-change baseline.
#
# Usage:
#   Rscript reports/perf_pass_baseline/compare_baseline.R [tolerance]
#
# tolerance = 0 (default) demands bit-identical results, which is the bar for
# plan items 1, 3, 6, 7 and 8. Items 2, 4 and 5 change floating-point
# summation order, so run them with a tolerance (1e-10).
#
# Bit-identity is the wrong bar for one change in this pass: adopting
# matrixStats::colSds() perturbs the column SDs by 1 ulp, which is enough to
# flip the sign of a principal component and, with it, the CCA weight's
# coordinate on that PC. Those two cancel and every reported quantity is
# invariant -- so the bar that matters is not "no bits moved" but "no reported
# quantity changed sign". The per-slot report below states both, because a
# tolerance alone would hide a sign flip: -x and x differ by a relative 2, which
# no sensible tolerance passes, but a *cancelling* pair of flips leaves the
# tolerance clean while the intermediate looks catastrophically different.

args <- commandArgs(trailingOnly = TRUE)
tolerance <- if (length(args) >= 1L) as.numeric(args[[1L]]) else 0

# Quantities a reader of the results actually sees. @skrCCAOut is deliberately
# excluded: it holds the raw CCA weight vectors, whose per-PC sign is a free
# convention that nothing downstream depends on.
REPORTED <- c("normalizedCorrelation", "sigmaValueChoice", "cellScores",
              "geneScores", "geneScoresRegression", "nullCorrelation",
              "pvalue")

flatten_numeric <- function(x) {
  if (is.null(x)) return(numeric(0))
  if (is.data.frame(x)) {
    return(unlist(lapply(x[vapply(x, is.numeric, logical(1))], as.numeric),
                  use.names = FALSE))
  }
  if (is.list(x)) {
    return(unlist(lapply(x, flatten_numeric), use.names = FALSE))
  }
  if (is.numeric(x)) return(as.numeric(x))
  numeric(0)
}

# Count sign changes among entries large enough for a sign to be meaningful.
sign_flips <- function(a, b) {
  if (length(a) != length(b) || !length(a)) return(0L)
  keep <- abs(a) > 1e-12 & abs(b) > 1e-12
  sum(sign(a[keep]) != sign(b[keep]))
}

max_rel <- function(a, b) {
  if (length(a) != length(b) || !length(a)) return(NA_real_)
  scale <- pmax(abs(a), abs(b))
  keep <- scale > 0
  if (!any(keep)) return(0)
  max(abs(a[keep] - b[keep]) / scale[keep])
}

base_dir <- "reports/perf_pass_baseline"
baseline <- readRDS(file.path(base_dir, "baseline.rds"))

# Capture into a directory that is emptied first. Without this, a capture that
# dies half way through leaves the previous run's baseline.rds in place and the
# comparison below happily reports all-OK against a stale file -- the one
# failure mode that would make this whole harness worse than useless.
tmp <- file.path(base_dir, "current")
unlink(tmp, recursive = TRUE)
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(tmp, "capture.log")
status <- system2("Rscript", c(file.path(base_dir, "capture_baseline.R"), tmp),
                  stdout = log_file, stderr = log_file)
if (!identical(status, 0L)) {
  cat("capture_baseline.R exited with status ", status, ":\n", sep = "")
  writeLines(paste0("  ", readLines(log_file, warn = FALSE)))
  stop("baseline capture failed; not comparing against a stale capture")
}
current_file <- file.path(tmp, "baseline.rds")
if (!file.exists(current_file)) {
  stop("baseline capture produced no ", current_file)
}
current <- readRDS(current_file)

if (!identical(names(baseline), names(current))) {
  stop("scenario sets differ.\n  only in baseline: ",
       paste(setdiff(names(baseline), names(current)), collapse = ", "),
       "\n  only in current:  ",
       paste(setdiff(names(current), names(baseline)), collapse = ", "),
       "\nRe-capture the baseline at the pre-pass commit after editing ",
       "capture_baseline.R.")
}

failures <- character()
flipped <- character()
worst_reported <- 0

for (scenario in names(baseline)) {
  b <- baseline[[scenario]]
  cur_s <- current[[scenario]]

  cmp <- if (tolerance == 0) {
    if (identical(b, cur_s)) TRUE else "not bit-identical"
  } else {
    all.equal(b, cur_s, tolerance = tolerance)
  }
  if (isTRUE(cmp)) {
    cat(sprintf("  OK   %s\n", scenario))
  } else {
    failures <- c(failures, scenario)
    cat(sprintf("  FAIL %s\n", scenario))
    cat(paste0("         ", utils::head(as.character(cmp), 5), "\n"))
  }

  # Independently of the tolerance verdict, report what happened to the
  # quantities a reader sees. A slot present in one capture and absent in the
  # other is a harness bug, not a numeric one, so say so explicitly.
  for (slot in intersect(REPORTED, union(names(b), names(cur_s)))) {
    fa <- flatten_numeric(b[[slot]])
    fb <- flatten_numeric(cur_s[[slot]])
    if (!length(fa) && !length(fb)) next
    if (length(fa) != length(fb)) {
      failures <- c(failures, scenario)
      cat(sprintf("       ! %-22s length %d vs %d\n", slot, length(fa),
                  length(fb)))
      next
    }
    nf <- sign_flips(fa, fb)
    mr <- max_rel(fa, fb)
    worst_reported <- max(worst_reported, mr, na.rm = TRUE)
    if (nf > 0L) {
      flipped <- c(flipped, paste0(scenario, "/", slot))
      cat(sprintf("       ! %-22s %d sign flip(s), max rel %.3e\n",
                  slot, nf, mr))
    }
  }
}

cat("\n")
failures <- unique(failures)
if (length(failures) == 0L) {
  cat(sprintf("All %d scenarios match (tolerance = %g).\n",
              length(baseline), tolerance))
} else {
  cat(sprintf("%d of %d scenarios differ: %s\n", length(failures),
              length(baseline), paste(failures, collapse = ", ")))
}
cat(sprintf("Reported quantities: max relative difference %.3e, %d sign flip(s)%s\n",
            worst_reported, length(flipped),
            if (length(flipped)) paste0(": ", paste(flipped, collapse = ", "))
            else ""))

# A sign flip in a reported quantity is a hard failure at any tolerance.
if (length(flipped) > 0L || length(failures) > 0L) quit(status = 1L)
