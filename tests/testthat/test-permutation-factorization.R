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
    expect_equal(lapply(fast$weights, function(w) lapply(w, abs)),
                 lapply(slow$weights, function(w) lapply(w, abs)),
                 tolerance = 1e-8, info = method)
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
  obj3 <- .fact_pipeline(n_cell_types = 3, nCC = 1)
  fast <- .permu_null(obj3, 77, permu_method = "global")
  slow <- .unfactorized(.permu_null(obj3, 77, permu_method = "global"))
  expect_equal(fast$ncorr, slow$ncorr, tolerance = 1e-8)

  obj_uns <- .fact_pipeline(nCC = 2, scalePCs = FALSE)
  fast_u <- .permu_null(obj_uns, 78, permu_method = "global")
  slow_u <- .unfactorized(.permu_null(obj_uns, 78, permu_method = "global"))
  expect_equal(fast_u$ncorr, slow_u$ncorr, tolerance = 1e-8)
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
  # ...and it must not drag the permutation index matrix along either.
  expect_lt(captured_bytes, as.numeric(utils::object.size(cell_permu)))
})

test_that("failed draws keep their slot when results are reassembled", {
  skip_on_cran()
  # The conditional test records a failed draw as NULL and counts it. If
  # reassembly dropped or shifted those slots, every later draw's statistic
  # would be filed against the wrong permutation.
  obj <- .fact_pipeline()
  cts <- obj@cellTypesOfInterest
  set.seed(2)
  cell_permu <- .getCellPermu(obj, "global", 8, cts)
  first_index <- cell_permu[[cts[2]]][1, ]
  # Fail on roughly half the draws, chosen from the draw's own content so the
  # worker stays a pure function of its spec.
  flaky <- function(spec) {
    idx <- spec[[2]][1]
    if (idx %% 2L == 0L) NULL else idx
  }
  expected_null <- (first_index %% 2L) == 0L
  expect_true(any(expected_null) && !all(expected_null))

  out <- .runPermutationDraws(cell_permu, cts, 8, flaky, n_cores = 1,
                              verbose = FALSE)
  expect_length(out, 8L)
  expect_equal(vapply(out, is.null, logical(1)), expected_null)

  if (!is.null(.installedCoProLibrary())) {
    par_out <- suppressMessages(
      .runPermutationDraws(cell_permu, cts, 8, flaky, n_cores = 3,
                           verbose = FALSE)
    )
    expect_equal(par_out, out)
  }
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
    .runPermutationDraws(cell_permu, cts, 7, worker, n_cores = 3,
                         verbose = FALSE)
  )
  expect_length(par_out, 7L)
  expect_equal(par_out, seq_out, tolerance = 1e-12)
})
