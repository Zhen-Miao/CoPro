# Extracted from test-permutation-factorization.R:531

# prequel ----------------------------------------------------------------------
.fact_pipeline <- function(n_cell_types = 2, nCC = 2, seed = 42,
                           scalePCs = TRUE) {
  obj <- create_test_copro_single(n_cells = 120,
                                  n_cell_types = n_cell_types, seed = seed)
  cts <- paste0("CellType", LETTERS[seq_len(n_cell_types)])
  obj <- subsetData(obj, cellTypesOfInterest = cts)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1, 0.2),
                             dropDistances = FALSE, verbose = FALSE)
  suppressWarnings(
    obj <- computePCA(obj, nPCA = 8, center = TRUE, scale. = TRUE,
                      scalePCs = scalePCs)
  )
  suppressWarnings(
    obj <- runSkrCCA(obj, scalePCs = scalePCs, nCC = nCC, maxIter = 200)
  )
  computeNormalizedCorrelation(obj)
}
.unfactorized <- function(expr) {
  old <- options(CoPro.factorizePermutation = FALSE)
  on.exit(options(old), add = TRUE)
  force(expr)
}
.permu_null <- function(obj, seed, ...) {
  set.seed(seed)
  out <- suppressMessages(suppressWarnings(
    runSkrCCAPermu(obj, nPermu = 12, verbose = FALSE, ...)
  ))
  out <- suppressMessages(computeNormalizedCorrelationPermu(out))
  list(
    weights = out@skrCCAPermuOut,
    ncorr = lapply(out@normalizedCorrelationPermu,
                   function(x) x$normalizedCorrelation),
    p = suppressWarnings(calculate_pvalue(out, cc_index = 1)$p_value)
  )
}
.flakyReassemblyCase <- function() {
  obj <- .fact_pipeline()
  cts <- obj@cellTypesOfInterest
  set.seed(2)
  cell_permu <- .getCellPermu(obj, "global", 8, cts)
  first_index <- cell_permu[[cts[2]]][1, ]
  flaky <- function(spec) {
    idx <- spec[[2]][1]
    if (idx %% 2L == 0L) NULL else idx
  }
  list(cts = cts, cell_permu = cell_permu, flaky = flaky,
       expected_null = (first_index %% 2L) == 0L)
}

# test -------------------------------------------------------------------------
skip_on_cran()
obj <- .fact_pipeline()
cts <- obj@cellTypesOfInterest
PCmats <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)
set.seed(5)
cell_permu <- .getCellPermu(obj, "global", 6, cts)
n_fixed <- nrow(PCmats[[cts[1]]])
legacy <- cell_permu
legacy[[cts[1]]] <- replicate(6, seq_len(n_fixed))
expect_true(.isIdentityPermutation(legacy[[cts[1]]]))
for (tt in seq_len(6)) {
    compact_draw <- .applyPermutationSpec(
      PCmats, .permutationDrawSpec(cell_permu, cts, tt))
    legacy_draw <- .applyPermutationSpec(
      PCmats, .permutationDrawSpec(legacy, cts, tt))
    expect_identical(compact_draw, legacy_draw)
    # And the fixed side is returned untouched, not copied through seq_len(n).
    expect_identical(compact_draw[[cts[1]]], PCmats[[cts[1]]])
  }
