# Tests for the factorized permutation operator (R/D0_permutation_plan.R).
#
# When a cell type is held fixed across draws, Y_ij = X_i' K_ij X_j is
# evaluated from a precomputed one-sided product instead of the sparse kernel.
# The identity is exact, so the whole point of these tests is that every
# permutation entry point returns the SAME null -- and therefore the same
# p-value -- as the unfactorized path, which
# options(CoPro.factorizePermutation = FALSE) restores.

.fact_pipeline <- function(n_cell_types = 2, nCC = 2, seed = 42,
                           scalePCs = TRUE) {
  obj <- create_test_copro_single(n_cells = 120,
                                  n_cell_types = n_cell_types, seed = seed)
  cts <- paste0("CellType", LETTERS[seq_len(n_cell_types)])
  obj <- subsetData(obj, cellTypesOfInterest = cts)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1, 0.2),
                             dropDistances = FALSE, verbose = FALSE, normalizeDistance = TRUE)
  suppressWarnings(
    obj <- computePCA(obj, nPCA = 8, center = TRUE, scale. = TRUE,
                      scalePCs = scalePCs)
  )
  suppressWarnings(
    obj <- runSkrCCA(obj, scalePCs = scalePCs, nCC = nCC, maxIter = 200)
  )
  computeNormalizedCorrelation(obj)
}

# Run `expr` with the factorization switched off.
.unfactorized <- function(expr) {
  old <- options(CoPro.factorizePermutation = FALSE)
  on.exit(options(old), add = TRUE)
  force(expr)
}

# ---------------------------------------------------------------------------
# Fixed-type detection
# ---------------------------------------------------------------------------

test_that("identity permutations are detected and non-identities are not", {
  expect_true(.isIdentityPermutation(replicate(5, 1:10)))
  expect_false(.isIdentityPermutation(replicate(5, sample(10))))
  # identical on draw 1 but shuffled later: must not be called fixed
  mixed <- cbind(1:10, sample(10), 1:10)
  expect_false(.isIdentityPermutation(mixed))
  # "pc" permutations are stored as a seed list, never as an index matrix
  expect_false(.isIdentityPermutation(list(type = "pc_permute", seeds = 1:5)))
  expect_false(.isIdentityPermutation(NULL))
})

test_that("permu_which decides which cell types the plan can factorize", {
  skip_on_cran()
  obj <- .fact_pipeline()
  cts <- obj@cellTypesOfInterest

  second <- .getCellPermu(obj, "global", 4, cts, permu_which = "second_only")
  expect_equal(unname(.fixedPermutationTypes(second, cts)), c(TRUE, FALSE))

  first <- .getCellPermu(obj, "global", 4, cts, permu_which = "first_only")
  expect_equal(unname(.fixedPermutationTypes(first, cts)), c(FALSE, TRUE))

  both <- .getCellPermu(obj, "global", 4, cts, permu_which = "both")
  expect_equal(unname(.fixedPermutationTypes(both, cts)), c(FALSE, FALSE))
})

# ---------------------------------------------------------------------------
# The factorized operator equals the sparse triple product
# ---------------------------------------------------------------------------

test_that("planned Y equals compute_Y_resi for every permu_which", {
  skip_on_cran()
  obj <- .fact_pipeline()
  cts <- obj@cellTypesOfInterest
  sigma <- obj@sigmaValueChoice
  PCmats <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)

  for (which in c("second_only", "first_only", "both")) {
    set.seed(3)
    cell_permu <- .getCellPermu(obj, "global", 3, cts, permu_which = which)
    plan <- .buildYPlan(PCmats, obj@kernelMatrices, sigma, cts,
                        .fixedPermutationTypes(cell_permu, cts))

    for (tt in 1:3) {
      local_pc <- .applyPermutationSpec(
        PCmats, .permutationDrawSpec(cell_permu, cts, tt)
      )
      expected <- compute_Y_resi(local_pc, obj@kernelMatrices, sigma, cts)
      planned <- .yResiFromPlan(plan, local_pc)
      for (i in cts) for (j in cts) {
        expect_equal(planned[[i]][[j]], expected[[i]][[j]],
                     tolerance = 1e-10,
                     info = paste(which, "draw", tt, i, j))
      }
    }
  }
})

test_that("planned Y equals compute_Y_resi with three cell types", {
  skip_on_cran()
  # Only the pairs touching the fixed type factorize; the rest fall back.
  obj <- .fact_pipeline(n_cell_types = 3, nCC = 1)
  cts <- obj@cellTypesOfInterest
  sigma <- obj@sigmaValueChoice
  PCmats <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)

  set.seed(9)
  cell_permu <- .getCellPermu(obj, "global", 2, cts,
                              permu_which = "second_only")
  plan <- .buildYPlan(PCmats, obj@kernelMatrices, sigma, cts,
                      .fixedPermutationTypes(cell_permu, cts))
  modes <- vapply(utils::combn(cts, 2L, simplify = FALSE), function(p) {
    plan$ops[[p[1]]][[p[2]]]$mode
  }, character(1))
  expect_true("left" %in% modes)   # pairs involving the fixed first type
  expect_true("full" %in% modes)   # the pair of two permuted types

  local_pc <- .applyPermutationSpec(
    PCmats, .permutationDrawSpec(cell_permu, cts, 1L)
  )
  expected <- compute_Y_resi(local_pc, obj@kernelMatrices, sigma, cts)
  planned <- .yResiFromPlan(plan, local_pc)
  for (i in cts) for (j in cts) {
    expect_equal(planned[[i]][[j]], expected[[i]][[j]], tolerance = 1e-10,
                 info = paste(i, j))
  }
})

test_that("planned Y matches the sparse product on float32 kernels", {
  skip_on_cran()
  # float32_csr_xky_cpp accumulates in single precision, so the factorized and
  # unfactorized orderings agree only to ~1e-6 relative here, not to the ~1e-15
  # of the double path. Large analyses use computeSparseKernelFloat32(), so
  # that path needs its own tolerance rather than inheriting the double one.
  obj <- create_test_copro_single(n_cells = 400, n_cell_types = 2, seed = 42)
  cts <- c("CellTypeA", "CellTypeB")
  obj <- subsetData(obj, cellTypesOfInterest = cts)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  sigma <- 0.1
  f32 <- computeSparseKernelFloat32(obj, sigmaValues = sigma, verbose = FALSE, normalizeDistance = TRUE)
  suppressWarnings(
    f32 <- computePCA(f32, nPCA = 8, center = TRUE, scale. = TRUE,
                      scalePCs = TRUE)
  )
  expect_true(.isFloat32SparseKernel(
    get_kernel_matrix_flat(f32@kernelMatrices, sigma, cts[1], cts[2])
  ))

  PCmats <- .getAllPCMats(f32@pcaGlobal, f32@scalePCs)
  set.seed(3)
  cell_permu <- .getCellPermu(f32, "global", 2, cts)
  plan <- .buildYPlan(PCmats, f32@kernelMatrices, sigma, cts,
                      .fixedPermutationTypes(cell_permu, cts))
  expect_identical(plan$ops[[cts[1]]][[cts[2]]]$mode, "left")

  local_pc <- .applyPermutationSpec(
    PCmats, .permutationDrawSpec(cell_permu, cts, 1L)
  )
  expected <- compute_Y_resi(local_pc, f32@kernelMatrices, sigma, cts)
  planned <- .yResiFromPlan(plan, local_pc)
  for (i in cts) for (j in cts) {
    if (is.null(expected[[i]][[j]])) {
      expect_null(planned[[i]][[j]])
      next
    }
    rel <- max(abs(planned[[i]][[j]] - expected[[i]][[j]])) /
      max(abs(expected[[i]][[j]]))
    expect_lt(rel, 1e-4)
  }
})

test_that("a fully fixed plan is constant and matches the unpermuted operator", {
  skip_on_cran()
  obj <- .fact_pipeline()
  cts <- obj@cellTypesOfInterest
  sigma <- obj@sigmaValueChoice
  PCmats <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)
  fixed <- stats::setNames(rep(TRUE, length(cts)), cts)

  plan <- .buildYPlan(PCmats, obj@kernelMatrices, sigma, cts, fixed)
  expect_identical(plan$ops[[cts[1]]][[cts[2]]]$mode, "const")
  expect_equal(.yResiFromPlan(plan, PCmats),
               compute_Y_resi(PCmats, obj@kernelMatrices, sigma, cts),
               tolerance = 1e-10)
})

test_that("Gram score norms equal the direct product, and 'pc' opts out", {
  skip_on_cran()
  obj <- .fact_pipeline()
  cts <- obj@cellTypesOfInterest
  sigma <- obj@sigmaValueChoice
  PCmats <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)
  Wobs <- obj@skrCCAOut[[paste0("sigma_", sigma)]]
  w1 <- lapply(Wobs, function(w) w[, 1, drop = FALSE])
  Y0 <- compute_Y_resi(PCmats, obj@kernelMatrices, sigma, cts)
  kinfo <- .get_ncorr_kernel_info(obj@kernelMatrices, sigma, cts)

  index_permu <- .getCellPermu(obj, "global", 2, cts, permu_which = "both")
  grams <- .permutationGrams(PCmats, index_permu, cts)
  expect_false(any(vapply(grams, is.null, logical(1))))
  expect_equal(
    .compute_ncorr_quick(PCmats, w1, obj@kernelMatrices, sigma, cts,
                         kernel_info = kinfo, Y_resi = Y0, grams = grams),
    .compute_ncorr_quick(PCmats, w1, obj@kernelMatrices, sigma, cts,
                         kernel_info = kinfo, Y_resi = Y0),
    tolerance = 1e-12
  )

  # A "pc" permutation shuffles each PC column separately, which does NOT
  # preserve X' X, so those types must keep the direct calculation.
  pc_permu <- .getCellPermu(obj, "pc", 2, cts, permu_which = "both")
  expect_true(all(vapply(.permutationGrams(PCmats, pc_permu, cts),
                         is.null, logical(1))))
})

test_that("the 'bin' null resamples, so it gets no Gram shortcut", {
  skip_on_cran()
  # The spatially matched patch null draws a resample: a cell can appear twice
  # in a draw and another not at all, so X[draw, ]' X[draw, ] != X' X. Only the
  # held-fixed type (an exact identity) keeps its Gram.
  obj <- .fact_pipeline()
  cts <- obj@cellTypesOfInterest
  PCmats <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)

  set.seed(4)
  bin_both <- .getCellPermu(obj, "bin", 4, cts, permu_which = "both")
  expect_false(any(vapply(cts, function(ct)
    .isRowPermutationMatrix(bin_both[[ct]]), logical(1))))
  expect_true(all(vapply(.permutationGrams(PCmats, bin_both, cts),
                         is.null, logical(1))))

  bin_second <- .getCellPermu(obj, "bin", 4, cts,
                              permu_which = "second_only")
  gram_ok <- !vapply(.permutationGrams(PCmats, bin_second, cts),
                     is.null, logical(1))
  expect_equal(unname(gram_ok), c(TRUE, FALSE))

  # global and toroidal nulls are genuine permutations
  expect_true(.isRowPermutationMatrix(
    .getCellPermu(obj, "global", 3, cts, permu_which = "both")[[cts[2]]]))
  expect_true(.isRowPermutationMatrix(
    .getCellPermu(obj, "toroidal", 3, cts, permu_which = "both")[[cts[2]]]))
})

# ---------------------------------------------------------------------------
# End-to-end: the null distribution is unchanged
# ---------------------------------------------------------------------------

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

test_that("runSkrCCAPermu null is identical with and without factorization", {
  skip_on_cran()
  obj <- .fact_pipeline()
  for (method in c("bin", "global", "pc")) {
    fast <- .permu_null(obj, 101, permu_method = method)
    slow <- .unfactorized(.permu_null(obj, 101, permu_method = method))
    expect_equal(fast$ncorr, slow$ncorr, tolerance = 1e-10, info = method)
    expect_equal(fast$p, slow$p, info = method)
    # Two cell types are closed-form on both paths, so the only difference is
    # floating-point reassociation (measured ~1e-15). Keep the tolerance close
    # to that: this is the quantity the RNG-driven sign-flip bug showed up in,
    # so it should not pass at the testthat default.
    expect_equal(lapply(fast$weights, function(w) lapply(w, abs)),
                 lapply(slow$weights, function(w) lapply(w, abs)),
                 tolerance = 1e-12, info = method)
  }
})

test_that("permu_which = 'both' and 'first_only' are also unchanged", {
  skip_on_cran()
  obj <- .fact_pipeline()
  for (which in c("both", "first_only")) {
    fast <- .permu_null(obj, 55, permu_method = "global", permu_which = which)
    slow <- .unfactorized(
      .permu_null(obj, 55, permu_method = "global", permu_which = which)
    )
    expect_equal(fast$ncorr, slow$ncorr, tolerance = 1e-10, info = which)
    expect_equal(fast$p, slow$p, info = which)
  }
})

test_that("three-type and unscaled-PC runs are unchanged", {
  skip_on_cran()
  # Three or more types still route through iterative block relaxation
  # (optimize_bilinear, tol = 1e-5) even though Y is now identical between the
  # two paths, so the converged weights can differ by more than reassociation
  # noise. Hence a looser tolerance here than in the closed-form two-type case.
  obj3 <- .fact_pipeline(n_cell_types = 3, nCC = 1)
  fast <- .permu_null(obj3, 77, permu_method = "global")
  slow <- .unfactorized(.permu_null(obj3, 77, permu_method = "global"))
  expect_equal(fast$ncorr, slow$ncorr, tolerance = 1e-8)

  # Two types with scalePCs = FALSE is closed-form on both paths.
  obj_uns <- .fact_pipeline(nCC = 2, scalePCs = FALSE)
  fast_u <- .permu_null(obj_uns, 78, permu_method = "global")
  slow_u <- .unfactorized(.permu_null(obj_uns, 78, permu_method = "global"))
  expect_equal(fast_u$ncorr, slow_u$ncorr, tolerance = 1e-12)
  expect_equal(fast_u$p, slow_u$p)
})

test_that("fair-sigma null is identical with and without factorization", {
  skip_on_cran()
  obj <- .fact_pipeline(nCC = 1)
  run <- function() {
    set.seed(303)
    out <- suppressMessages(suppressWarnings(
      runSkrCCAPermu_FairSigma(obj, nPermu = 12, permu_method = "bin",
                               verbose = FALSE)
    ))
    list(
      ncorr = vapply(out@normalizedCorrelationPermu,
                     function(d) d$normalizedCorrelation[1], numeric(1)),
      sigma = vapply(out@normalizedCorrelationPermu,
                     function(d) d$sigmaValues[1], numeric(1)),
      p = suppressWarnings(calculate_pvalue(out, cc_index = 1)$p_value)
    )
  }
  fast <- run()
  slow <- .unfactorized(run())
  expect_equal(fast$ncorr, slow$ncorr, tolerance = 1e-10)
  expect_equal(fast$sigma, slow$sigma)
  expect_equal(fast$p, slow$p)
})

test_that("conditional step-down null is identical with and without factorization", {
  skip_on_cran()
  obj <- .fact_pipeline(nCC = 2)
  run <- function() {
    set.seed(404)
    out <- suppressMessages(suppressWarnings(
      runSkrCCAPermu_Conditional(obj, nPermu = 12, permu_method = "bin",
                                 verbose = FALSE)
    ))
    list(stats = out@conditionalPermu$perm_stats,
         sigma = out@conditionalPermu$perm_sigma,
         per_axis = calculate_pvalue_stepdown(out))
  }
  fast <- run()
  slow <- .unfactorized(run())
  expect_equal(fast$stats, slow$stats, tolerance = 1e-10)
  expect_equal(fast$sigma, slow$sigma)
  expect_equal(fast$per_axis$p_raw, slow$per_axis$p_raw)
  expect_equal(fast$per_axis$p_stepdown, slow$per_axis$p_stepdown)
})

# ---------------------------------------------------------------------------
# Worker payload: what a PSOCK worker actually receives
# ---------------------------------------------------------------------------

test_that("the permutation worker does not capture the kernels", {
  skip_on_cran()
  # Regression: the old workers closed over the CoPro object (or had the whole
  # kernel list clusterExport()ed), so peak memory scaled with n_cores. The
  # worker must now carry only the plan and the PC matrices.
  obj <- .fact_pipeline()
  cts <- obj@cellTypesOfInterest
  sigma <- obj@sigmaValueChoice
  PCmats <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)
  cell_permu <- .getCellPermu(obj, "global", 50, cts)
  plan <- .buildYPlan(PCmats, obj@kernelMatrices, sigma, cts,
                      .fixedPermutationTypes(cell_permu, cts))
  worker <- .makeSkrCCAPermuWorker(PCmats, plan, cts, nCC = 1,
                                   sdev2_list = NULL, maxIter = 50, tol = 1e-5)

  # The closure's environment is the factory frame and nothing else, so what a
  # worker receives is exactly this list.
  captured <- environment(worker)
  expect_setequal(
    ls(captured),
    c("PCmats", "plan", "cts", "nCC", "sdev2_list", "maxIter", "tol")
  )
  expect_false(any(vapply(ls(captured), function(nm)
    methods::is(get(nm, envir = captured), "CoPro"), logical(1))))

  captured_bytes <- sum(vapply(ls(captured), function(nm)
    as.numeric(utils::object.size(get(nm, envir = captured))), numeric(1)))
  expect_lt(captured_bytes,
            as.numeric(utils::object.size(obj@kernelMatrices)))

  # ...and it must not drag the permutation draws along either. Comparing
  # against object.size(cell_permu) used to stand in for that, but a held-fixed
  # type is now a compact marker rather than an n x nPermu matrix, so
  # cell_permu is legitimately smaller than the PC matrices the worker does
  # need. Assert the invariant the size comparison was proxying for: what the
  # worker captures does not grow with the number of draws.
  worker_10 <- .makeSkrCCAPermuWorker(PCmats, plan, cts, nCC = 1,
                                      sdev2_list = NULL, maxIter = 50,
                                      tol = 1e-5)
  bytes_of <- function(w) sum(vapply(ls(environment(w)), function(nm)
    as.numeric(utils::object.size(get(nm, envir = environment(w)))),
    numeric(1)))
  expect_identical(bytes_of(worker_10), captured_bytes)
  expect_false(any(vapply(ls(captured), function(nm)
    identical(get(nm, envir = captured), cell_permu), logical(1))))
})

# Shared fixture for the two reassembly tests. The worker fails on roughly half
# the draws, chosen from the draw's own content so it stays a pure function of
# its spec and gives the same failures sequentially and in parallel.
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

test_that("failed draws keep their slot when results are reassembled", {
  skip_on_cran()
  # The conditional test records a failed draw as NULL and counts it. If
  # reassembly dropped or shifted those slots, every later draw's statistic
  # would be filed against the wrong permutation.
  case <- .flakyReassemblyCase()
  expect_true(any(case$expected_null) && !all(case$expected_null))

  out <- .runPermutationDraws(case$cell_permu, case$cts, 8, case$flaky,
                              n_cores = 1, verbose = FALSE)
  expect_length(out, 8L)
  expect_equal(vapply(out, is.null, logical(1)), case$expected_null)
})

# R CMD check sets _R_CHECK_LIMIT_CORES_, under which parallel::makeCluster()
# refuses to spawn more than two workers and errors out. The package already
# honours this variable when sizing the C++ thread pool
# (R/12c_float32_sparse_kernel.R); these tests must too. Two workers still
# exercise the split-and-reassemble path, which is what is under test.
.test_worker_count <- function(preferred = 3L) {
  limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if (nzchar(limit) && !identical(limit, "false")) 2L else preferred
}

test_that("parallel reassembly of failed draws matches sequential", {
  skip_on_cran()
  # Gated with skip() rather than a silent `if`: under devtools::load_all()
  # .installedCoProLibrary() is NULL and .runPermutationDraws() falls back to
  # sequential, so an `if` would make this a no-op that still reported as a
  # passing test. The skip keeps the gap visible in the test summary.
  if (is.null(.installedCoProLibrary())) {
    skip("CoPro is not installed; the PSOCK path falls back to sequential.")
  }
  case <- .flakyReassemblyCase()
  out <- .runPermutationDraws(case$cell_permu, case$cts, 8, case$flaky,
                              n_cores = 1, verbose = FALSE)
  par_out <- suppressMessages(
    .runPermutationDraws(case$cell_permu, case$cts, 8, case$flaky,
                         n_cores = .test_worker_count(), verbose = FALSE)
  )
  expect_equal(par_out, out)
})

test_that("parallel and sequential draws agree", {
  skip_on_cran()
  if (is.null(.installedCoProLibrary())) {
    skip("CoPro is not installed; the PSOCK path falls back to sequential.")
  }
  obj <- .fact_pipeline()
  cts <- obj@cellTypesOfInterest
  sigma <- obj@sigmaValueChoice
  PCmats <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)
  set.seed(12)
  cell_permu <- .getCellPermu(obj, "global", 7, cts)
  plan <- .buildYPlan(PCmats, obj@kernelMatrices, sigma, cts,
                      .fixedPermutationTypes(cell_permu, cts))
  worker <- .makeSkrCCAPermuWorker(PCmats, plan, cts, nCC = 1,
                                   sdev2_list = NULL, maxIter = 50, tol = 1e-5)

  seq_out <- .runPermutationDraws(cell_permu, cts, 7, worker, n_cores = 1,
                                  verbose = FALSE)
  par_out <- suppressMessages(
    .runPermutationDraws(cell_permu, cts, 7, worker, n_cores = .test_worker_count(),
                         verbose = FALSE)
  )
  expect_length(par_out, 7L)
  expect_equal(par_out, seq_out, tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# Compact permutation storage
# ---------------------------------------------------------------------------

test_that("held-fixed cell types are stored compactly, not as index matrices", {
  skip_on_cran()
  obj <- .fact_pipeline()
  cts <- obj@cellTypesOfInterest

  for (method in c("global", "bin", "toroidal", "pc")) {
    cell_permu <- suppressWarnings(.getCellPermu(obj, method, 40, cts))
    fixed_entry <- cell_permu[[cts[1]]]

    expect_true(.isIdentityPermuEntry(fixed_entry), info = method)
    expect_true(.isIdentityPermutation(fixed_entry), info = method)
    expect_false(is.matrix(fixed_entry), info = method)
    # A bijection, so the Gram shortcut still applies to the fixed side.
    expect_true(.isRowPermutationMatrix(fixed_entry), info = method)

    # The whole point: the marker does not grow with the number of draws.
    many <- suppressWarnings(.getCellPermu(obj, method, 400, cts))
    expect_identical(many[[cts[1]]], fixed_entry, info = method)
  }
})

test_that("a compact fixed entry yields the same draw as an explicit matrix", {
  skip_on_cran()
  obj <- .fact_pipeline()
  cts <- obj@cellTypesOfInterest
  PCmats <- .getAllPCMats(obj@pcaGlobal, obj@scalePCs)
  set.seed(5)
  cell_permu <- .getCellPermu(obj, "global", 6, cts)

  # An object saved by an older CoPro stores the held-fixed type as an explicit
  # replicate(nPermu, 1:n). Both representations must resolve identically.
  n_fixed <- nrow(PCmats[[cts[1]]])
  legacy <- cell_permu
  legacy[[cts[1]]] <- replicate(6, seq_len(n_fixed))

  expect_true(.isIdentityPermutation(legacy[[cts[1]]]))
  for (tt in seq_len(6)) {
    compact_draw <- .applyPermutationSpec(
      PCmats, .permutationDrawSpec(cell_permu, cts, tt))
    legacy_draw <- .applyPermutationSpec(
      PCmats, .permutationDrawSpec(legacy, cts, tt))

    # Values and labels must agree exactly. They are compared without
    # attributes on purpose: the legacy path reached the fixed matrix through
    # `[`, which drops the "scaled:scale" attribute that .getAllPCMats() adds
    # under scalePCs = TRUE, while the compact path returns the matrix
    # untouched and keeps it. No consumer reads that attribute.
    for (ct in cts) {
      expect_equal(unname(compact_draw[[ct]]), unname(legacy_draw[[ct]]),
                   ignore_attr = TRUE, info = paste(ct, "draw", tt))
      expect_identical(dimnames(compact_draw[[ct]]),
                       dimnames(legacy_draw[[ct]]))
    }
    # And the fixed side is returned untouched, not copied through seq_len(n).
    expect_identical(compact_draw[[cts[1]]], PCmats[[cts[1]]])
  }
})

test_that("worker chunking preserves compact entries and slices seeds", {
  skip_on_cran()
  obj <- .fact_pipeline()
  cts <- obj@cellTypesOfInterest
  cell_permu <- .getCellPermu(obj, "pc", 9, cts)

  seen <- list()
  recorder <- function(spec) {
    seen[[length(seen) + 1L]] <<- spec
    NULL
  }
  invisible(.runPermutationDraws(cell_permu, cts, 9, recorder, n_cores = 1L,
                                 verbose = FALSE))
  expect_length(seen, 9L)
  # Held-fixed type: identity marker on every draw.
  expect_true(all(vapply(seen, function(s) isTRUE(s[[cts[1]]]$identity),
                         logical(1))))
  # Permuted "pc" type: the per-draw seeds arrive in draw order.
  expect_identical(
    vapply(seen, function(s) s[[cts[2]]]$pc_seed, numeric(1)),
    as.numeric(cell_permu[[cts[2]]]$seeds)
  )
})

test_that("compact permutation storage is opt-in and produces valid draws", {
  skip_on_cran()
  obj <- .fact_pipeline()
  cts <- obj@cellTypesOfInterest
  n_permuted <- sum(obj@cellTypesSub == cts[2])

  # Off by default: the permuted side is still an explicit index matrix, so
  # saved analyses reproduce exactly.
  expect_false(.useCompactPermutation())
  default_permu <- .getCellPermu(obj, "global", 5, cts)
  expect_true(is.matrix(default_permu[[cts[2]]]))

  old <- options(CoPro.compactPermutation = TRUE)
  on.exit(options(old), add = TRUE)
  expect_true(.useCompactPermutation())

  compact <- .getCellPermu(obj, "global", 5, cts)
  expect_identical(compact[[cts[2]]]$type, "global_seed")
  expect_length(compact[[cts[2]]]$seeds, 5L)
  expect_true(.isRowPermutationMatrix(compact[[cts[2]]]))

  # Every regenerated draw is a genuine permutation, and re-reading a draw
  # gives the same answer (it is a pure function of the stored seed).
  for (tt in seq_len(5)) {
    spec <- .permutationDrawSpec(compact, cts, tt)
    expect_setequal(spec[[cts[2]]], seq_len(n_permuted))
    expect_identical(spec, .permutationDrawSpec(compact, cts, tt))
  }

  # "bin" is a resample rather than a bijection, so it must not claim the
  # Gram shortcut even under compact storage.
  compact_bin <- suppressWarnings(.getCellPermu(obj, "bin", 4, cts))
  expect_identical(compact_bin[[cts[2]]]$type, "bin_seed")
  expect_false(.isRowPermutationMatrix(compact_bin[[cts[2]]]))
  bin_spec <- .permutationDrawSpec(compact_bin, cts, 2L)
  expect_length(bin_spec[[cts[2]]], n_permuted)
  expect_identical(bin_spec, .permutationDrawSpec(compact_bin, cts, 2L))
})
