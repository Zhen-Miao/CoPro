# Numeric baseline for the performance/memory pass.
#
# Runs the documented Quick Start pipeline on a fixture at single- and
# multi-slide, across the double-sparse and float32 kernel backends, plus
# permutation tests, and saves every downstream numeric result. Re-running this
# after each commit and diffing against the saved baseline is how the pass
# proves it did not change the science.
#
# Usage:  Rscript reports/perf_pass_baseline/capture_baseline.R <outdir>

suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[[1L]] else
  "reports/perf_pass_baseline/current"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

source("tests/testthat/helper-test-data.R")

SIGMAS <- c(0.05, 0.1)
CTS <- c("CellTypeA", "CellTypeB")

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# ---------------------------------------------------------------- single slide

run_single <- function(kernel_backend) {
  obj <- quiet(create_test_copro_single(n_cells = 320, n_genes = 60,
                                        n_cell_types = 2, seed = 11))
  obj <- quiet(subsetData(obj, cellTypesOfInterest = CTS))
  obj <- quiet(computePCA(obj, nPCA = 8))
  obj <- quiet(switch(
    kernel_backend,
    dense = {
      obj <- computeDistance(obj, distType = "Euclidean2D", verbose = FALSE)
      computeKernelMatrix(obj, sigmaValues = SIGMAS, verbose = FALSE)
    },
    sparse  = computeSparseKernel(obj, sigmaValues = SIGMAS, verbose = FALSE),
    float32 = computeSparseKernelFloat32(obj, sigmaValues = SIGMAS,
                                         verbose = FALSE)
  ))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 2))
  obj <- quiet(computeNormalizedCorrelation(obj))
  obj <- quiet(computeGeneAndCellScores(obj))
  obj <- quiet(computeRegressionGeneScores(obj))
  obj
}

# ------------------------------------------------------- single, one cell type

run_within_type <- function() {
  obj <- quiet(create_test_copro_single(n_cells = 300, n_genes = 60,
                                        n_cell_types = 2, seed = 13))
  obj <- quiet(subsetData(obj, cellTypesOfInterest = "CellTypeA"))
  obj <- quiet(computePCA(obj, nPCA = 8))
  obj <- quiet(computeSparseKernelFloat32(obj, sigmaValues = SIGMAS,
                                          verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 2))
  obj <- quiet(computeNormalizedCorrelation(obj))
  obj
}

# ----------------------------------------------------------------- multi slide

run_multi <- function(center_per_slide, scalePCs) {
  obj <- quiet(create_test_copro_multi(n_cells_per_slide = 220, n_slides = 2,
                                       n_genes = 60, n_cell_types = 2,
                                       seed = 23))
  obj <- quiet(subsetData(obj, cellTypesOfInterest = CTS))
  obj <- quiet(computePCA(obj, nPCA = 8, scalePCs = scalePCs,
                          center_per_slide = center_per_slide))
  obj <- quiet(computeSparseKernel(obj, sigmaValues = SIGMAS, verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = scalePCs, nCC = 2))
  obj <- quiet(computeNormalizedCorrelation(obj))
  obj <- quiet(computeGeneAndCellScores(obj))
  obj
}

# ------------------------------------------------------------------ permutation

run_permutation <- function(permu_method, permu_which = "second_only") {
  obj <- quiet(create_test_copro_single(n_cells = 260, n_genes = 50,
                                        n_cell_types = 2, seed = 37))
  obj <- quiet(subsetData(obj, cellTypesOfInterest = CTS))
  obj <- quiet(computePCA(obj, nPCA = 6))
  obj <- quiet(computeSparseKernel(obj, sigmaValues = SIGMAS, verbose = FALSE))
  obj <- quiet(runSkrCCA(obj, scalePCs = TRUE, nCC = 1))
  obj <- quiet(computeNormalizedCorrelation(obj))
  set.seed(101)
  obj <- quiet(runSkrCCAPermu(obj, nPermu = 20, permu_method = permu_method,
                              permu_which = permu_which, verbose = FALSE))
  obj <- quiet(computeNormalizedCorrelationPermu(obj))
  obj
}

# ---------------------------------------------------------------------- capture

has <- function(obj, name) name %in% methods::slotNames(obj)

capture <- function(obj) {
  list(
    skrCCAOut             = obj@skrCCAOut,
    normalizedCorrelation = obj@normalizedCorrelation,
    sigmaValueChoice      = obj@sigmaValueChoice,
    cellScores            = if (has(obj, "cellScores")) obj@cellScores,
    geneScores            = if (has(obj, "geneScores")) obj@geneScores,
    geneScoresRegression  = if (has(obj, "geneScoresRegression"))
      obj@geneScoresRegression
  )
}

results <- list()

for (backend in c("dense", "sparse", "float32")) {
  message("== single / ", backend)
  results[[paste0("single_", backend)]] <- capture(run_single(backend))
}

message("== single / within-type float32")
results[["within_type_float32"]] <- capture(run_within_type())

for (cps in c(FALSE, TRUE)) {
  for (spc in c(TRUE, FALSE)) {
    key <- sprintf("multi_cps%s_scale%s", cps, spc)
    message("== ", key)
    results[[key]] <- capture(run_multi(cps, spc))
  }
}

permu_cases <- list(
  list(method = "global",   which = "second_only"),
  list(method = "global",   which = "both"),
  list(method = "global",   which = "first_only"),
  list(method = "bin",      which = "second_only"),
  list(method = "toroidal", which = "second_only"),
  list(method = "pc",       which = "second_only")
)
for (case in permu_cases) {
  key <- paste0("permu_", case$method, "_", case$which)
  message("== ", key)
  obj <- run_permutation(case$method, case$which)
  results[[key]] <- c(
    capture(obj),
    list(
      nullCorrelation = obj@normalizedCorrelationPermu,
      pvalue          = quiet(calculate_pvalue(obj, cc_index = 1))
    )
  )
}

saveRDS(results, file.path(outdir, "baseline.rds"))
cat("Wrote", file.path(outdir, "baseline.rds"), "with",
    length(results), "scenarios\n")
