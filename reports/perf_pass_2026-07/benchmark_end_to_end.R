# End-to-end pipeline timing and peak memory, for A/B between two checkouts.
#
#   Rscript reports/perf_pass_2026-07/benchmark_end_to_end.R <label> <outdir> [reps]
#
# Run once in a checkout of the pre-pass commit and once on the branch head,
# then compare the two CSVs.
#
# The pipeline is stateful -- every stage consumes the previous stage's object
# -- so a stage cannot be repeated on its own. The whole pipeline is run `reps`
# times instead and the per-stage minimum reported, which is the standard
# defence against a co-tenant process inflating a single sample. Peak memory is
# near-deterministic and is reported as the maximum over repetitions.

suppressMessages(devtools::load_all(".", quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
label <- if (length(args) >= 1L) args[[1L]] else "run"
outdir <- if (length(args) >= 2L) args[[2L]] else "."
REPS <- if (length(args) >= 3L) as.integer(args[[3L]]) else 3L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

N_CELLS <- 60000
N_GENES <- 120
SIGMAS <- c(0.05, 0.1)
NPCA <- 30

q <- function(e) suppressWarnings(suppressMessages(e))

# Peak Vcells allocated during `fun`, above the resting set, in MB.
timed <- function(stage, fun) {
  invisible(gc(reset = TRUE, full = TRUE))
  before <- gc()[2, "max used"]
  elapsed <- system.time(value <- fun())[["elapsed"]]
  mb <- (gc()[2, "max used"] - before) * 8 / 1e6
  rows[[length(rows) + 1L]] <<- data.frame(
    label = label, stage = stage, seconds = elapsed, peak_mb = mb)
  message(sprintf("  %-28s %7.2f s  %8.0f MB", stage, elapsed, mb))
  value
}

set.seed(4)
side <- 1000
coords <- data.frame(x = runif(N_CELLS) * side, y = runif(N_CELLS) * side)
cell_types <- rep(c("CellTypeA", "CellTypeB"), length.out = N_CELLS)
ids <- paste0("cell", seq_len(N_CELLS))
rownames(coords) <- ids
expression <- matrix(pmax(0, rnorm(N_CELLS * N_GENES, 2, 1)),
                     N_CELLS, N_GENES,
                     dimnames = list(ids, paste0("g", seq_len(N_GENES))))
meta <- data.frame(row.names = ids, cellType = cell_types)

message("== ", label, ": ", N_CELLS, " cells, ", N_GENES, " genes, nPCA ",
        NPCA, ", ", length(SIGMAS), " sigmas, ", REPS, " repetitions")

run_pipeline <- function() {
  rows <<- list()
  obj <- q(newCoProSingle(normalizedData = expression, locationData = coords,
                          metaData = meta, cellTypes = cell_types))
  obj <- q(subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB")))

  obj <- timed("computePCA", function() q(computePCA(obj, nPCA = NPCA)))
  obj <- timed("computeSparseKernelFloat32",
               function() q(computeSparseKernelFloat32(
                 obj, sigmaValues = SIGMAS, verbose = FALSE)))
  obj <- timed("runSkrCCA",
               function() q(runSkrCCA(obj, scalePCs = TRUE, nCC = 2)))
  obj <- timed("computeNormalizedCorrelation",
               function() q(computeNormalizedCorrelation(obj)))
  obj <- timed("computeGeneAndCellScores",
               function() q(computeGeneAndCellScores(obj)))
  obj <- timed("runSkrCCAPermu (99 draws)", function() {
    set.seed(11)
    q(runSkrCCAPermu(obj, nPermu = 99, permu_method = "global",
                     verbose = FALSE))
  })
  permu_bytes <- as.numeric(utils::object.size(obj@cellPermu))
  obj <- timed("computeNormalizedCorrelationPermu",
               function() q(computeNormalizedCorrelationPermu(obj)))

  rows[[length(rows) + 1L]] <<- data.frame(
    label = label, stage = "@cellPermu stored size",
    seconds = NA_real_, peak_mb = permu_bytes / 1e6)
  do.call(rbind, rows)
}

all_runs <- list()
for (rep_index in seq_len(REPS)) {
  message("-- repetition ", rep_index, " of ", REPS)
  one <- run_pipeline()
  one$rep <- rep_index
  all_runs[[rep_index]] <- one
  message("   subtotal: ", round(sum(one$seconds, na.rm = TRUE), 1), " s")
}
runs <- do.call(rbind, all_runs)

stages <- unique(runs$stage)
out <- do.call(rbind, lapply(stages, function(s) {
  part <- runs[runs$stage == s, ]
  data.frame(label = label, stage = s,
             seconds = suppressWarnings(min(part$seconds, na.rm = TRUE)),
             seconds_median = stats::median(part$seconds),
             seconds_max = suppressWarnings(max(part$seconds, na.rm = TRUE)),
             peak_mb = max(part$peak_mb), reps = REPS)
}))
out$seconds[!is.finite(out$seconds)] <- NA_real_
out$seconds_max[!is.finite(out$seconds_max)] <- NA_real_

write.csv(out, file.path(outdir, paste0("end_to_end_", label, ".csv")),
          row.names = FALSE)
write.csv(runs, file.path(outdir, paste0("end_to_end_", label, "_runs.csv")),
          row.names = FALSE)
for (i in seq_len(nrow(out))) {
  message(sprintf("  %-28s min %7.2f s  median %7.2f s  %8.0f MB",
                  out$stage[i], out$seconds[i], out$seconds_median[i],
                  out$peak_mb[i]))
}
message("total (min per stage): ", round(sum(out$seconds, na.rm = TRUE), 1), " s")
