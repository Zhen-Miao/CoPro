.diagnostic_fixture <- function() {
  obj <- create_test_copro_multi(
    n_cells_per_slide = 50, n_slides = 3, n_genes = 20,
    n_cell_types = 2, seed = 901
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computePCA(obj, nPCA = 4)
  obj <- computeDistance(obj, normalizeDistance = TRUE, verbose = FALSE)
  computeKernelMatrix(obj, sigmaValues = c(0.08, 0.12), verbose = FALSE)
}

.diagnostic_ops <- function(Gx = diag(c(1, 4))) {
  A <- matrix(c(1, 1, 0, 0), 2)
  list(
    slides = "s1", cell_types = c("x", "y"),
    G = list(s1 = list(x = Gx, y = diag(2))),
    Y = list(s1 = list(x = list(y = A), y = list(x = t(A)))),
    n = list(s1 = list(x = 10L, y = 10L))
  )
}

test_that("saved diagnostics measure every returned bandwidth and axis", {
  obj <- .diagnostic_fixture()
  fit <- suppressMessages(runSkrCCA(obj, nCC = 2, maxIter = 1000))
  records <- getCCADiagnostics(fit)
  expect_named(records, c("sigma_0.08", "sigma_0.12"))
  cts <- fit@cellTypesOfInterest
  slides <- getSlideList(fit)
  X <- CoPro:::.preparePCMatrices(
    fit@pcaResults, fit@pcaGlobal, scalePCs = TRUE, slides = slides, cts = cts
  )
  for (sigma in obj@sigmaValues) {
    d <- getCCADiagnostics(fit, sigma)
    expect_identical(d$components$component, 1:2)
    expect_true(all(d$components$solver == "full_gradient"))
    expect_true(all(d$components$converged))
    expect_true(all(d$components$status == "gradient_tolerance"))
    expect_true(all(d$conditioning$positive_definite))
    expect_false(any(d$components$floor_encountered))
    expect_false(any(d$score_norms$floor_active))
    weights <- fit@skrCCAOut[[CoPro:::.sigmaName(sigma)]]
    ops <- CoPro:::.computeSlideOperators(
      X, fit@kernelMatrices, sigma, slides, cts
    )
    for (cc in 1:2) {
      w <- lapply(weights, function(m) m[, cc, drop = FALSE])
      prev <- if (cc == 1) NULL else lapply(weights, function(m) m[, 1, drop = FALSE])
      g <- CoPro:::.sumcorTangentGradient(
        CoPro:::.sumcorGradient(w, ops, "equal"), w, prev
      )
      expect_equal(d$components$gradient_norm[cc],
                   max(vapply(g, function(x) sqrt(sum(x^2)), numeric(1))))
      # Reconstruct the actual evaluation from scores and kernels, rather than
      # trusting the stored trace or the diagnostic builder's objective.
      values <- vapply(slides, function(s) {
        a <- X[[s]][[cts[1]]] %*% w[[cts[1]]]
        b <- X[[s]][[cts[2]]] %*% w[[cts[2]]]
        as.numeric(crossprod(w[[cts[1]]],
                            ops$Y[[s]][[cts[1]]][[cts[2]]] %*% w[[cts[2]]])) /
          sqrt(sum(a^2) * sum(b^2))
      }, numeric(1))
      expect_equal(d$components$objective[cc], mean(values), tolerance = 1e-12)
      trace <- d$objective_traces[[paste0("CC", cc)]]
      expect_length(trace, d$components$iterations[cc] + 1L)
      expect_equal(tail(trace, 1), d$components$objective[cc])
      expect_true(all(diff(trace) >= -1e-12))
    }
  }
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(fit, path)
  expect_identical(getCCADiagnostics(readRDS(path)), records)
  # Scores and gene-score computation must preserve the solver records.
  scored <- computeGeneAndCellScores(computeNormalizedCorrelation(fit))
  expect_identical(getCCADiagnostics(scored), records)
})

test_that("diagnostics are comparable across PC parameterizations", {
  obj <- .diagnostic_fixture()
  fit <- function(scalePCs) suppressMessages(runSkrCCA(
    obj, scalePCs = scalePCs, nCC = 2, maxIter = 1000, sigmaChoice = 0.08
  ))
  a <- getCCADiagnostics(fit(TRUE), 0.08)
  b <- getCCADiagnostics(fit(FALSE), 0.08)
  expect_identical(a$coordinate_system, "whitened_pc")
  expect_equal(a$conditioning, b$conditioning, tolerance = 1e-10)
  expect_equal(a$components$objective, b$components$objective, tolerance = 1e-8)
  expect_equal(a$score_norms, b$score_norms, tolerance = 1e-7)
})

test_that("supplied axes cannot inherit another fit's convergence claim", {
  obj <- .diagnostic_fixture()
  first <- suppressMessages(runSkrCCA(obj, nCC = 1, sigmaChoice = 0.08))
  supplied <- first@skrCCAOut[["sigma_0.08"]]
  for (ncc in c(1L, 2L)) {
    fit <- suppressMessages(runSkrCCA(
      obj, transferred_weight_1 = supplied, nCC = ncc,
      sigmaChoice = 0.08, maxIter = 1000
    ))
    d <- getCCADiagnostics(fit, 0.08)
    expect_identical(d$components$solver[1], "supplied")
    expect_identical(d$components$status[1], "supplied")
    expect_true(is.na(d$components$converged[1]))
    expect_length(d$objective_traces$CC1, 0L)
    if (ncc == 2L) expect_identical(d$components$solver[2], "full_gradient")
  }
})

test_that("algebraic SUMCOV reductions carry measured residuals", {
  ops <- .diagnostic_ops(Gx = diag(2))
  first <- optimize_sumcor_pca(NULL, NULL, 0.1, "s1", c("x", "y"), ops = ops)
  both <- optimize_sumcor_pca_n(NULL, NULL, 0.1, "s1", c("x", "y"),
                                first, nCC = 2, ops = ops)
  for (w in list(first, both)) {
    d <- attr(w, "ccaDiagnostics")
    expect_true(all(d$components$solver == "sumcov_reduction"))
    expect_true(all(d$components$converged))
    expect_true(all(d$components$iterations == 0L))
    expect_true(all(d$components$gradient_norm <= 1e-6))
    expect_true(all(lengths(d$objective_traces) == 0L))
  }
  supplied <- first
  attr(supplied, "ccaDiagnostics") <- NULL
  extended <- optimize_sumcor_pca_n(NULL, NULL, 0.1, "s1", c("x", "y"),
                                    supplied, nCC = 2, ops = ops)
  d <- attr(extended, "ccaDiagnostics")
  expect_identical(d$components$solver, c("supplied", "sumcov_reduction"))
  expect_true(is.na(d$components$converged[1]))
  expect_true(d$components$converged[2])
})

test_that("iteration limits and line-search failure are not convergence", {
  ops <- .diagnostic_ops()
  w <- list(x = matrix(c(1, 1) / sqrt(2)), y = matrix(c(1, 0)))
  expect_warning(
    limited <- CoPro:::.sumcorIterate(w, ops, "equal", max_iter = 1, tol = 1e-12),
    "did not converge"
  )
  expect_false(limited$converged)
  expect_identical(limited$stop_reason, "max_iter")
  local_mocked_bindings(
    .sumcorRetract = function(...) list(x = NULL, y = NULL), .package = "CoPro"
  )
  expect_warning(
    stalled <- CoPro:::.sumcorIterate(w, ops, "equal", tol = 1e-12),
    "line search stalled"
  )
  expect_false(stalled$converged)
  expect_identical(stalled$stop_reason, "line_search_stalled")
  expect_identical(stalled$iterations, 0L)
})

test_that("a floored rejected trial is recorded even with healthy final norms", {
  ops <- .diagnostic_ops(Gx = diag(c(0, 1)))
  # The nonzero starting direction has finite score norms. Reject a candidate
  # on the null space by returning a nonfinite objective there.
  w <- list(x = matrix(c(0, 1)), y = matrix(c(1, 0)))
  original_objective <- CoPro:::.sumcorObjective
  local_mocked_bindings(
    .sumcorRetract = function(...) list(x = matrix(c(1, 0)), y = matrix(c(1, 0))),
    .sumcorObjective = function(w_list, ops, slideWeight, sigma_all = NULL) {
      if (w_list$x[2, 1] == 0) return(NA_real_)
      original_objective(w_list, ops, slideWeight, sigma_all)
    }, .package = "CoPro"
  )
  expect_warning(
    fit <- CoPro:::.sumcorIterate(w, ops, "equal", tol = 1e-12),
    "line search stalled"
  )
  d <- CoPro:::.recordSumcorAxis(
    CoPro:::.newSumcorDiagnostics(ops), fit$w_list, ops, "equal", 1L, 1e-12,
    "full_gradient", fit
  )
  expect_true(d$components$floor_encountered)
  expect_false(any(d$score_norms$floor_active))
})

test_that("a floored denominator and singular Gram remain visible", {
  ops <- .diagnostic_ops(Gx = diag(c(0, 1)))
  ops$Y$s1$x$y[] <- 0
  ops$Y$s1$y$x[] <- 0
  w <- list(x = matrix(c(1, 0)), y = matrix(c(1, 0)))
  fit <- CoPro:::.sumcorIterate(w, ops, "equal")
  d <- CoPro:::.recordSumcorAxis(
    CoPro:::.newSumcorDiagnostics(ops), w, ops, "equal", 1L, 1e-6,
    "full_gradient", fit
  )
  expect_true(d$components$floor_encountered)
  expect_identical(d$score_norms$score_norm, c(0, 1))
  expect_identical(d$score_norms$floor_active, c(TRUE, FALSE))
  expect_identical(d$conditioning$positive_definite, c(FALSE, TRUE))
  expect_identical(d$conditioning$rank, c(1L, 2L))
  expect_identical(d$conditioning$condition_number[1], Inf)
})

test_that("missing diagnostics remain missing and refitting clears old ones", {
  obj <- .diagnostic_fixture()
  fit <- suppressMessages(runSkrCCA(obj, nCC = 1, sigmaChoice = 0.08))
  old <- fit
  attr(old@skrCCAOut[["sigma_0.08"]], "ccaDiagnostics") <- NULL
  expect_null(getCCADiagnostics(old, 0.08))
  expect_identical(getCCADiagnostics(old), list(sigma_0.08 = NULL))
  covariance <- suppressMessages(runSkrCCA(fit, nCC = 1, objective = "sumcov"))
  expect_true(all(vapply(getCCADiagnostics(covariance), is.null, logical(1))))
  expect_identical(getCCADiagnostics(obj), setNames(list(), character()))
  expect_error(getCCADiagnostics(fit, 0.4), "No PCA-space fit")
  expect_error(getCCADiagnostics(fit, c(0.08, 0.12)), "single finite")
  expect_error(getCCADiagnostics(list()), "CoPro object")
})
