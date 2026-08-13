# Tests for the PCA-space SUMCOR objective (objective = "sumcor").
#
# The claims being pinned here, in order:
#   1. With one slide, sumcor and sumcov are the same problem only when the
#      per-slide Gram matrices and per-pair constants permit the reduction;
#      otherwise runSkrCCA() optimizes SUMCOR itself.
#   2. SUMCOV is the sigma-weighted special case of the SUMCOR family.
#   3. sumcor is invariant to per-slide score rescaling; sumcov is not.
#   4. The full analytical gradient matches finite differences, and projected
#      gradient ascent is monotone and reaches a constrained stationary point.
#   5. It is deterministic (no RNG).
#   6. scalePCs stays a pure reparametrization.
#   7. Thin slides are dropped under sumcor, reported under sumcov.
#   8. Argument validation.
#   9. getCCAObjective() provenance, including for gene-space fits, and
#      space = "gene" forwarding rather than silently dropping arguments.

# ---------------------------------------------------------------- fixtures ---

.sumcor_fixture <- function(nct = 2L, n_slides = 3L, nPCA = 8L, seed = 11L,
                            n_cells_per_slide = 70L) {
  # A one-slide CoProMulti is a deliberate test configuration (it exercises the
  # is_multi = TRUE, S = 1 route), so the constructor's advisory about it is
  # expected here and muffled rather than left to clutter the report.
  obj <- withCallingHandlers(
    create_test_copro_multi(
      n_cells_per_slide = n_cells_per_slide, n_slides = n_slides,
      n_genes = 40, n_cell_types = max(nct, 2L), seed = seed
    ),
    warning = function(w) {
      if (grepl("only one unique slide ID", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  cts <- paste0("CellType", LETTERS[seq_len(nct)])
  obj <- subsetData(obj, cellTypesOfInterest = cts)
  obj <- computePCA(obj, nPCA = nPCA)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE,
                      normalizeDistance = TRUE)
}

# Pull out the pieces the internal solvers take, so tests can drive them
# directly rather than only through runSkrCCA().
.sumcor_parts <- function(obj, cts, scalePCs = TRUE) {
  slides <- getSlideList(obj)
  sdev2 <- if (scalePCs) NULL else {
    setNames(lapply(cts, function(ct) obj@pcaGlobal[[ct]]$sdev^2), cts)
  }
  X <- CoPro:::.preparePCMatrices(
    pc_data = obj@pcaResults, pca_global = obj@pcaGlobal,
    scalePCs = scalePCs, slides = slides, cts = cts
  )
  ops <- CoPro:::.computeSlideOperators(X, obj@kernelMatrices, 0.1, slides,
                                        cts, n_cores = 1)
  list(X = X, ops = ops, slides = slides, sdev2 = sdev2, cts = cts)
}

.unit_weights <- function(parts, seed = 1L) {
  set.seed(seed)
  setNames(lapply(parts$cts, function(ct) {
    v <- matrix(rnorm(ncol(parts$ops$G[[parts$slides[1]]][[ct]])), ncol = 1)
    CoPro:::normalize_vec_weighted(v, parts$sdev2[[ct]])
  }), parts$cts)
}

# ------------------------------------------------- 1. single-slide identity ---

.single_slide_fixture <- function(nct = 2L, nPCA = 8L, seed = 11L,
                                  counts = NULL) {
  obj <- create_test_copro_single(
    n_cells = if (is.null(counts)) 210 else sum(counts), n_genes = 40,
    n_cell_types = max(nct, 2L), seed = seed
  )
  cts <- paste0("CellType", LETTERS[seq_len(nct)])
  # Deliberately unbalanced cell counts: that is what makes the per-pair
  # constants differ and the one-slide reduction inexact.
  if (!is.null(counts)) obj@cellTypes <- rep(cts, times = counts)
  obj <- subsetData(obj, cellTypesOfInterest = cts)
  obj <- computePCA(obj, nPCA = nPCA)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE,
                      normalizeDistance = TRUE)
}

# A one-slide operator set built directly, so the mathematics can be tested
# without a CoPro object in the way. `G_i = (n_i - 1) I` is exactly what
# whitened PCs give, which is the premise the reduction argument rests on.
.synthetic_one_slide_ops <- function(counts, nPC = 3L, seed = 5L) {
  set.seed(seed)
  cts <- paste0("CellType", LETTERS[seq_along(counts)])
  counts <- setNames(as.numeric(counts), cts)
  s <- "Slide1"

  G <- setNames(lapply(cts, function(ct) diag(counts[[ct]] - 1, nPC)), cts)
  Y <- setNames(lapply(cts, function(a) {
    setNames(lapply(cts, function(b) matrix(0, nPC, nPC)), cts)
  }), cts)
  for (i in seq_along(cts)) {
    for (j in seq_along(cts)) {
      if (i < j) {
        M <- matrix(rnorm(nPC * nPC), nPC, nPC)
        Y[[cts[i]]][[cts[j]]] <- M
        Y[[cts[j]]][[cts[i]]] <- t(M)
      }
    }
  }
  list(Y = setNames(list(Y), s), G = setNames(list(G), s),
       n = setNames(list(counts), s), slides = s, cell_types = cts)
}

test_that("the one-slide reduction to sumcov is exact only when the per-pair constants coincide", {
  # The non-circular counterpart to the routing test below. sigma is a norm, so
  # on ||w_i|| = 1 with whitened PCs the denominators are sqrt(n_i - 1), not 1,
  # and the objective is SUMCOV reweighted by m_ij / sqrt((n_i-1)(n_j-1)). That
  # per-pair constant leaves the argmax alone only when every pair gets the same
  # one: at most one pair, or equal cell counts.
  for (sw in c("equal", "size")) {
    for (case in list(list(k = c(70, 70, 70), exact = TRUE),
                      list(k = c(4000, 60, 40), exact = FALSE),
                      list(k = c(150, 60), exact = TRUE))) {
      ops <- .synthetic_one_slide_ops(case$k)
      lbl <- sprintf("slideWeight=%s counts=%s", sw,
                     paste(case$k, collapse = "/"))
      expect_equal(CoPro:::.sumcorReducesToSumcov(ops, sw), case$exact,
                   info = lbl)

      # Check the prediction it encodes: warm-start from the SUMCOV solution and
      # run the SUMCOR iteration. Exact => there is nothing left to gain.
      cts <- ops$cell_types
      warm <- CoPro:::.sumcorWarmStart(ops, cts, NULL, nCC = 1L)
      w0 <- setNames(lapply(cts, function(ct) warm[[ct]][, 1L, drop = FALSE]),
                     cts)
      f0 <- CoPro:::.sumcorObjective(w0, ops, sw)
      fit <- suppressWarnings(CoPro:::.sumcorIterate(
        w_init = w0, ops = ops, slideWeight = sw,
        max_iter = 300, tol = 1e-8, verbose = FALSE
      ))
      f1 <- CoPro:::.sumcorObjective(fit$w_list, ops, sw)

      if (case$exact) {
        expect_equal(f1, f0, tolerance = 1e-8,
                     info = paste("exact case should not improve:", lbl))
      } else {
        # Strict: the whole point of the inexact case is that the SUMCOV
        # solution is *not* already optimal for SUMCOR, so a run that silently
        # fell back to the shortcut would give f1 == f0 and must fail here.
        expect_gt(f1, f0 + 1e-9)
      }
    }
  }

  # "covariance" is SUMCOV by construction, so it is always exact.
  expect_true(CoPro:::.sumcorReducesToSumcov(
    .synthetic_one_slide_ops(c(4000, 60, 40)), "covariance"
  ))
})

test_that("an inexact one-slide sumcor request optimizes sumcor itself", {
  obj <- .single_slide_fixture(nct = 3L, counts = c(150, 60, 40))
  fit_cor <- suppressMessages(suppressWarnings(
    runSkrCCA(obj, nCC = 1, objective = "sumcor", slideWeight = "equal")
  ))
  fit_cov <- suppressMessages(suppressWarnings(
    runSkrCCA(obj, nCC = 1, objective = "sumcov")
  ))

  cts <- c("CellTypeA", "CellTypeB", "CellTypeC")
  X <- CoPro:::.preparePCMatrices(
    pca_global = obj@pcaGlobal, scalePCs = TRUE, cts = cts
  )
  ops <- CoPro:::.computeSingleSlideOperators(
    X, obj@kernelMatrices, 0.1, cts
  )
  weights <- function(fit) setNames(lapply(cts, function(ct) {
    fit@skrCCAOut[["sigma_0.1"]][[ct]][, 1L, drop = FALSE]
  }), cts)
  f_cor <- CoPro:::.sumcorObjective(weights(fit_cor), ops, "equal")
  f_cov <- CoPro:::.sumcorObjective(weights(fit_cov), ops, "equal")

  # Strict, for the same reason as above: `getCCAObjective()` records "sumcor"
  # on the shortcut path too, so only a genuine improvement distinguishes a run
  # that optimized SUMCOR from one that returned the SUMCOV maximizer.
  expect_gt(f_cor, f_cov + 1e-9)
  expect_equal(getCCAObjective(fit_cor)$objective, "sumcor")
  expect_equal(getCCAObjective(fit_cor)$slideWeight, "equal")
})

test_that("one-slide sumcor and sumcov weights agree in exact reduction cases", {
  # One and two cell types have at most one term, so once the one-slide PC Gram
  # matrices are scalar identities, the two criteria differ only by a constant.
  # Check both object classes and both PC parametrizations.
  for (kind in c("single", "multi")) {
    for (nct in c(1L, 2L)) {
      for (scale_pcs in c(TRUE, FALSE)) {
        obj <- if (kind == "single") {
          .single_slide_fixture(nct = nct)
        } else {
          .sumcor_fixture(nct = nct, n_slides = 1L)
        }
        cts <- paste0("CellType", LETTERS[seq_len(nct)])

        fit_cov <- suppressMessages(suppressWarnings(
          runSkrCCA(obj, scalePCs = scale_pcs, nCC = 2, objective = "sumcov")
        ))
        fit_cor <- suppressMessages(suppressWarnings(
          runSkrCCA(obj, scalePCs = scale_pcs, nCC = 2, objective = "sumcor")
        ))

        for (ct in cts) {
          w_cov <- fit_cov@skrCCAOut[["sigma_0.1"]][[ct]]
          w_cor <- fit_cor@skrCCAOut[["sigma_0.1"]][[ct]]
          for (cc in seq_len(ncol(w_cov))) {
            expect_equal(
              as.numeric(w_cor[, cc]),
              .align_sign(as.numeric(w_cor[, cc]), as.numeric(w_cov[, cc])),
              tolerance = 1e-8,
              info = sprintf("%s nct=%d scalePCs=%s ct=%s cc=%d",
                             kind, nct, scale_pcs, ct, cc)
            )
          }
        }
      }
    }
  }

  # Three types are also exact when their cell counts (and hence pair
  # constants) are equal.
  obj <- .single_slide_fixture(nct = 3L, counts = c(70, 70, 70))
  fit_cov <- suppressMessages(runSkrCCA(obj, nCC = 1, objective = "sumcov"))
  fit_cor <- suppressMessages(runSkrCCA(obj, nCC = 1, objective = "sumcor"))
  for (ct in paste0("CellType", LETTERS[1:3])) {
    a <- fit_cov@skrCCAOut[["sigma_0.1"]][[ct]][, 1L]
    b <- fit_cor@skrCCAOut[["sigma_0.1"]][[ct]][, 1L]
    expect_equal(a, .align_sign(a, b), tolerance = 1e-5)
  }
})

test_that("an exact one-slide shortcut still records the requested sumcor criterion", {
  obj <- .single_slide_fixture(nct = 2L)
  fit <- suppressMessages(runSkrCCA(obj, nCC = 1, objective = "sumcor"))
  # The decomposition is a computational shortcut for this criterion, not a
  # change to the objective the result solves.
  expect_equal(getCCAObjective(fit)$objective, "sumcor")
  expect_equal(getCCAObjective(fit)$requested, "sumcor")
  expect_equal(getCCAObjective(fit)$slideWeight, "equal")
})

# ------------------------------------------ 2. sumcov as a special case ------

test_that("the covariance slide weight reproduces the sumcov objective", {
  # Pins the factorization the whole design rests on: SUMCOV is the member of
  # the SUMCOR family whose slide weight is sigma_i * sigma_j.
  for (nct in c(1L, 2L, 3L)) {
    for (scale_pcs in c(TRUE, FALSE)) {
      obj <- .sumcor_fixture(nct = nct)
      cts <- paste0("CellType", LETTERS[seq_len(nct)])
      parts <- .sumcor_parts(obj, cts, scale_pcs)
      w <- .unit_weights(parts)

      Yagg <- CoPro:::.aggregateSlideOperators(parts$ops)
      direct <- if (nct == 1L) {
        as.numeric(crossprod(w[[cts[1]]], Yagg[[cts[1]]][[cts[1]]] %*% w[[cts[1]]]))
      } else {
        sum(vapply(combn(cts, 2, simplify = FALSE), function(p) {
          as.numeric(crossprod(w[[p[1]]], Yagg[[p[1]]][[p[2]]] %*% w[[p[2]]]))
        }, numeric(1)))
      }

      expect_equal(
        CoPro:::.sumcorObjective(w, parts$ops, "covariance"),
        direct, tolerance = 1e-10,
        info = sprintf("nct=%d scalePCs=%s", nct, scale_pcs)
      )
    }
  }
})

test_that("the analytical SUMCOR gradient matches finite differences", {
  for (nct in c(1L, 3L)) {
    obj <- .sumcor_fixture(nct = nct, nPCA = 5L)
    cts <- paste0("CellType", LETTERS[seq_len(nct)])
    parts <- .sumcor_parts(obj, cts)
    w <- .unit_weights(parts, seed = 818)

    for (weighting in c("equal", "size", "covariance")) {
      gradient <- CoPro:::.sumcorGradient(w, parts$ops, weighting)
      errors <- numeric(0)
      eps <- 1e-6
      for (ct in cts) {
        for (j in seq_len(nrow(w[[ct]]))) {
          plus <- minus <- w
          plus[[ct]][j, 1L] <- plus[[ct]][j, 1L] + eps
          minus[[ct]][j, 1L] <- minus[[ct]][j, 1L] - eps
          numerical <- (
            CoPro:::.sumcorObjective(plus, parts$ops, weighting) -
              CoPro:::.sumcorObjective(minus, parts$ops, weighting)
          ) / (2 * eps)
          errors <- c(errors, abs(numerical - gradient[[ct]][j, 1L]))
        }
      }
      expect_lt(max(errors), 1e-6)
    }
  }
})

# ------------------------------------- 3. per-slide scale (in)variance -------

test_that("sumcor is invariant to per-slide score rescaling but sumcov is not", {
  # Inflating one slide's scores is the batch mode sumcor exists to neutralize.
  obj <- .sumcor_fixture(nct = 2L)
  cts <- c("CellTypeA", "CellTypeB")
  parts <- .sumcor_parts(obj, cts)
  w <- .unit_weights(parts)

  X_inflated <- parts$X
  X_inflated[[parts$slides[1]]] <- lapply(
    X_inflated[[parts$slides[1]]], function(m) m * 5
  )
  ops_inflated <- CoPro:::.computeSlideOperators(
    X_inflated, obj@kernelMatrices, 0.1, parts$slides, cts, n_cores = 1
  )

  for (weighting in c("equal", "size")) {
    expect_equal(
      CoPro:::.sumcorObjective(w, ops_inflated, weighting),
      CoPro:::.sumcorObjective(w, parts$ops, weighting),
      tolerance = 1e-9,
      info = paste("slideWeight =", weighting)
    )
  }

  expect_false(isTRUE(all.equal(
    CoPro:::.sumcorObjective(w, ops_inflated, "covariance"),
    CoPro:::.sumcorObjective(w, parts$ops, "covariance"),
    tolerance = 1e-6
  )))
})

test_that("slideWeight = size tracks cell count while equal does not", {
  obj <- .sumcor_fixture(nct = 2L)
  cts <- c("CellTypeA", "CellTypeB")
  parts <- .sumcor_parts(obj, cts)
  w <- .unit_weights(parts)

  n_first <- parts$ops$n[[parts$slides[1]]]
  expect_equal(
    CoPro:::.sumcorSlideWeight("size", parts$ops, parts$slides[1],
                               cts[1], cts[2], NULL),
    sqrt(as.numeric(n_first[[cts[1]]]) * as.numeric(n_first[[cts[2]]]))
  )
  expect_equal(
    CoPro:::.sumcorSlideWeight("equal", parts$ops, parts$slides[1],
                               cts[1], cts[2], NULL),
    1
  )
})

# ----------------------------------------------- 4. ascent guarantees --------

test_that("sumcor never returns worse than its sumcov warm start", {
  for (nct in c(1L, 2L, 3L)) {
    for (weighting in c("size", "equal")) {
      obj <- .sumcor_fixture(nct = nct)
      cts <- paste0("CellType", LETTERS[seq_len(nct)])
      parts <- .sumcor_parts(obj, cts)

      warm <- CoPro:::.sumcorWarmStart(parts$ops, cts, parts$sdev2, 1L)
      warm1 <- setNames(lapply(cts, function(ct) warm[[ct]][, 1, drop = FALSE]),
                        cts)
      obj_warm <- CoPro:::.sumcorObjective(warm1, parts$ops, weighting)

      fit <- suppressWarnings(CoPro:::.sumcorIterate(
        w_init = warm1, ops = parts$ops, slideWeight = weighting,
        sdev2_list = parts$sdev2, max_iter = 300, tol = 1e-10
      ))

      expect_gte(fit$objective, obj_warm - 1e-12)
    }
  }
})

test_that("projected-gradient SUMCOR is monotone and reaches stationarity", {
  obj <- .sumcor_fixture(nct = 3L)
  cts <- c("CellTypeA", "CellTypeB", "CellTypeC")
  parts <- .sumcor_parts(obj, cts)
  warm <- CoPro:::.sumcorWarmStart(parts$ops, cts, parts$sdev2, 1L)
  w <- setNames(lapply(cts, function(ct) warm[[ct]][, 1, drop = FALSE]), cts)

  fit <- suppressWarnings(CoPro:::.sumcorIterate(
    w_init = w, ops = parts$ops, slideWeight = "equal",
    sdev2_list = parts$sdev2, max_iter = 1000, tol = 1e-6
  ))
  expect_true(all(diff(fit$objective_trace) >= -1e-12),
              info = paste(sprintf("%.10f", fit$objective_trace), collapse = " "))
  expect_lte(fit$gradient_norm, 1e-6)
})

test_that("sumcor weights satisfy the CCA constraint on every block", {
  for (scale_pcs in c(TRUE, FALSE)) {
    obj <- .sumcor_fixture(nct = 3L)
    cts <- c("CellTypeA", "CellTypeB", "CellTypeC")
    fit <- suppressMessages(suppressWarnings(
      runSkrCCA(obj, scalePCs = scale_pcs, nCC = 2, objective = "sumcor")
    ))
    sdev2 <- if (scale_pcs) NULL else {
      setNames(lapply(cts, function(ct) obj@pcaGlobal[[ct]]$sdev^2), cts)
    }
    for (ct in cts) {
      W <- fit@skrCCAOut[["sigma_0.1"]][[ct]]
      for (cc in seq_len(ncol(W))) {
        w <- as.numeric(W[, cc])
        norm2 <- if (scale_pcs) sum(w^2) else sum(w^2 * sdev2[[ct]])
        expect_equal(norm2, 1, tolerance = 1e-8,
                     info = sprintf("scalePCs=%s ct=%s cc=%d",
                                    scale_pcs, ct, cc))
      }
    }
  }
})

# ----------------------------------------------------- 5. determinism -------

test_that("sumcor is deterministic and independent of the RNG stream", {
  obj <- .sumcor_fixture(nct = 3L)

  set.seed(1)
  a <- suppressMessages(suppressWarnings(
    runSkrCCA(obj, nCC = 2, objective = "sumcor")
  ))
  set.seed(99999)
  invisible(runif(37))  # advance the stream
  b <- suppressMessages(suppressWarnings(
    runSkrCCA(obj, nCC = 2, objective = "sumcor")
  ))

  for (ct in names(a@skrCCAOut[["sigma_0.1"]])) {
    expect_identical(a@skrCCAOut[["sigma_0.1"]][[ct]],
                     b@skrCCAOut[["sigma_0.1"]][[ct]],
                     info = ct)
  }
})

# --------------------------------------- 6. scalePCs reparametrization ------

test_that("scalePCs stays a pure reparametrization under sumcor", {
  # Same property test-scalePCs-equivalence.R asserts for sumcov: the canonical
  # scores must not depend on whether the PCs were whitened, only the weights'
  # parametrization should.
  for (nct in c(2L, 3L)) {
    obj <- .sumcor_fixture(nct = nct)
    cts <- paste0("CellType", LETTERS[seq_len(nct)])

    fit_t <- suppressMessages(suppressWarnings(
      runSkrCCA(obj, scalePCs = TRUE, nCC = 2, objective = "sumcor")
    ))
    fit_f <- suppressMessages(suppressWarnings(
      runSkrCCA(obj, scalePCs = FALSE, nCC = 2, objective = "sumcor")
    ))

    for (ct in cts) {
      sdev <- obj@pcaGlobal[[ct]]$sdev
      W_t <- fit_t@skrCCAOut[["sigma_0.1"]][[ct]]
      W_f <- fit_f@skrCCAOut[["sigma_0.1"]][[ct]]
      for (cc in seq_len(ncol(W_t))) {
        # w_unwhitened * sdev is the whitened-space weight.
        mapped <- as.numeric(W_f[, cc]) * sdev
        expect_equal(
          as.numeric(W_t[, cc]),
          .align_sign(as.numeric(W_t[, cc]), mapped),
          tolerance = 1e-6,
          info = sprintf("nct=%d ct=%s cc=%d", nct, ct, cc)
        )
      }
    }
  }
})

test_that("sumcor axes are metric-orthogonal within each cell type", {
  for (scale_pcs in c(TRUE, FALSE)) {
    obj <- .sumcor_fixture(nct = 3L)
    cts <- c("CellTypeA", "CellTypeB", "CellTypeC")
    fit <- suppressMessages(suppressWarnings(
      runSkrCCA(obj, scalePCs = scale_pcs, nCC = 3, objective = "sumcor")
    ))
    for (ct in cts) {
      W <- fit@skrCCAOut[["sigma_0.1"]][[ct]]
      d <- if (scale_pcs) rep(1, nrow(W)) else obj@pcaGlobal[[ct]]$sdev^2
      gram <- t(W) %*% (W * d)
      expect_equal(gram, diag(ncol(W)), tolerance = 1e-8,
                   info = sprintf("scalePCs=%s ct=%s", scale_pcs, ct))
    }
  }
})

# -------------------------------------------------- 7. degenerate slides ----

test_that("thin slides are dropped under sumcor and only reported under sumcov", {
  obj <- .sumcor_fixture(nct = 2L, n_slides = 3L)
  cts <- c("CellTypeA", "CellTypeB")

  # A threshold above every slide's count but one leaves a single usable slide,
  # at which point sumcor legitimately reduces to sumcov.
  counts <- vapply(getSlideList(obj), function(s) {
    idx <- getSlideID(obj) == s
    min(vapply(cts, function(ct) sum(idx & obj@cellTypesSub == ct), integer(1)))
  }, integer(1))
  threshold <- sort(counts)[2]  # drops at least the thinnest slide

  expect_warning(
    suppressMessages(runSkrCCA(obj, nCC = 1, objective = "sumcor",
                               minCellsPerSlide = threshold)),
    "dropped"
  )

  # sumcov must not drop anything: it only reports.
  fit_cov <- suppressWarnings(suppressMessages(
    runSkrCCA(obj, nCC = 1, objective = "sumcov")
  ))
  expect_length(getCCAObjective(fit_cov)$droppedSlides, 0)
  expect_equal(getCCAObjective(fit_cov)$objective, "sumcov")
})

test_that(".dropDegenerateSlides errors when nothing survives", {
  obj <- .sumcor_fixture(nct = 2L, n_slides = 2L)
  cts <- c("CellTypeA", "CellTypeB")
  parts <- .sumcor_parts(obj, cts)
  expect_error(
    suppressWarnings(CoPro:::.dropDegenerateSlides(
      parts$X, cts, minCells = 1e6
    )),
    "no slide has at least"
  )
})

# ------------------------------------------------- 8. argument validation ---

test_that("slideWeight is rejected under sumcov", {
  obj <- .sumcor_fixture(nct = 2L)
  expect_error(
    runSkrCCA(obj, nCC = 1, objective = "sumcov", slideWeight = "size"),
    "applies only to objective"
  )
})

test_that("an unknown objective or slideWeight is rejected", {
  obj <- .sumcor_fixture(nct = 2L)
  expect_error(runSkrCCA(obj, nCC = 1, objective = "sumsquares"))
  expect_error(
    runSkrCCA(obj, nCC = 1, objective = "sumcor", slideWeight = "inverse"),
    "arg"
  )
})

test_that("step_size is accepted under sumcor and 1 stays a no-op", {
  obj <- .sumcor_fixture(nct = 2L, n_slides = 3L)
  cts <- c("CellTypeA", "CellTypeB")
  w <- function(fit) {
    setNames(lapply(cts, function(ct) {
      fit@skrCCAOut[["sigma_0.1"]][[ct]][, 1L]
    }), cts)
  }
  fit <- function(...) suppressMessages(suppressWarnings(
    runSkrCCA(obj, nCC = 1, objective = "sumcor", ...)
  ))

  default <- fit()
  undamped <- fit(step_size = 1)
  damped <- fit(step_size = 0.5)

  # Passing the default explicitly must not perturb the iteration at all.
  expect_identical(w(undamped), w(default))

  # Damping trades iterations for stability; it must not move the maximizer.
  for (ct in cts) {
    cosine <- abs(sum(w(damped)[[ct]] * w(default)[[ct]])) /
      sqrt(sum(w(damped)[[ct]]^2) * sum(w(default)[[ct]]^2))
    expect_equal(cosine, 1, tolerance = 1e-5, info = ct)
  }
  expect_equal(getCCAObjective(damped)$objective, "sumcor")
})

test_that("damping changes the SUMCOR trajectory rather than being ignored", {
  # The guard against a silent plumbing regression: if step_size were dropped
  # anywhere between runSkrCCA() and the line search, or were absorbed by the
  # adaptive step hint, both runs would take the same number of iterations.
  obj <- .sumcor_fixture(nct = 3L, n_slides = 3L)
  cts <- paste0("CellType", LETTERS[seq_len(3L)])
  parts <- .sumcor_parts(obj, cts)

  npc <- ncol(parts$ops$G[[parts$slides[[1L]]]][[cts[[1L]]]])
  w0 <- setNames(lapply(cts, function(ct) {
    matrix(1 / sqrt(npc), nrow = npc, ncol = 1L)
  }), cts)
  run <- function(ss) suppressWarnings(CoPro:::.sumcorIterate(
    w0, parts$ops, "equal", step_size = ss, max_iter = 5000, tol = 1e-7
  ))

  undamped <- run(1)
  damped <- run(0.1)

  expect_gt(damped$iterations, undamped$iterations)
  # Both must actually converge, or the iteration counts compare two caps.
  expect_lte(undamped$gradient_norm, 1e-7)
  expect_lte(damped$gradient_norm, 1e-7)
  # Same stationary point, and monotone getting there.
  expect_equal(damped$objective, undamped$objective, tolerance = 1e-6)
  expect_true(all(diff(damped$objective_trace) >= -1e-12))
})

test_that("out-of-range step_size is still rejected under sumcor", {
  obj <- .sumcor_fixture(nct = 2L)
  for (bad in list(0, -0.1, 1.5, "0.5", c(0.5, 0.5), NA_real_)) {
    expect_error(
      suppressMessages(runSkrCCA(obj, nCC = 1, objective = "sumcor",
                                 step_size = bad)),
      "step_size must be a single numeric value",
      info = paste(utils::capture.output(str(bad)), collapse = " ")
    )
  }
})

test_that("minCellsPerSlide must be a single non-negative number", {
  obj <- .sumcor_fixture(nct = 2L)
  expect_error(
    runSkrCCA(obj, nCC = 1, minCellsPerSlide = -1),
    "non-negative"
  )
  expect_error(
    runSkrCCA(obj, nCC = 1, minCellsPerSlide = c(1, 2)),
    "single non-negative"
  )
})

test_that("space = gene forwards to the gene-space implementation", {
  obj <- .sumcor_fixture(nct = 2L)
  fit <- suppressWarnings(suppressMessages(
    runSkrCCA(obj, space = "gene", nCC = 1, sigmaChoice = 0.1,
              max_iter = 200, tol = 1e-4, verbose = FALSE)
  ))
  # Gene space stores under its own key prefix and fills geneScores directly.
  expect_true("gscca_sigma_0.1" %in% names(fit@skrCCAOut))
  expect_gt(length(fit@geneScores), 0)
})

test_that("space = gene needs an unambiguous bandwidth", {
  obj <- create_test_copro_multi(
    n_cells_per_slide = 60, n_slides = 2, n_genes = 30,
    n_cell_types = 2, seed = 11
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computePCA(obj, nPCA = 8)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  obj <- computeKernelMatrix(obj, sigmaValues = c(0.05, 0.1), verbose = FALSE,
                             normalizeDistance = TRUE)
  expect_error(
    runSkrCCA(obj, space = "gene", nCC = 1),
    "one bandwidth at a time"
  )
})

# --------------------------------------------------------- provenance -------

test_that("CoProMulti defaults to equal-slide sumcor and records provenance", {
  obj <- .sumcor_fixture(nct = 2L)

  fit <- suppressMessages(suppressWarnings(
    runSkrCCA(obj, nCC = 1)
  ))
  record <- getCCAObjective(fit)
  expect_equal(record$space, "pca")
  expect_equal(record$objective, "sumcor")
  expect_equal(record$slideWeight, "equal")
  expect_equal(length(record$slides), length(getSlideList(obj)))

  # An object predating the provenance attribute still reads back as sumcov,
  # which is what legacy objects were computed under.
  stripped <- fit
  attr(stripped@skrCCAOut, "ccaObjective") <- NULL
  expect_equal(getCCAObjective(stripped)$objective, "sumcov")
})

test_that("a permutation null matches the criterion the weights were fitted with", {
  # Cell-level permutation is single-slide only (CoProMulti is rejected before
  # this runs), so the SUMCOR reduction test applies in full.
  fit <- function(obj, objective) {
    suppressMessages(suppressWarnings(
      runSkrCCA(obj, nCC = 1, objective = objective)
    ))
  }
  resolve <- function(f, cts) {
    suppressMessages(CoPro:::.resolvePermutationObjective(f, cts, TRUE))
  }

  # Two cell types: one pair, so the per-pair constant is uniform and SUMCOR's
  # maximizer *is* SUMCOV's. The existing SUMCOV null is already the matching
  # null, and refusing here would be refusing a valid test.
  two <- .single_slide_fixture(nct = 2L, counts = c(120, 80))
  cts2 <- c("CellTypeA", "CellTypeB")
  expect_equal(resolve(fit(two, "sumcor"), cts2)$objective, "sumcov")
  expect_equal(resolve(fit(two, "sumcov"), cts2)$objective, "sumcov")

  # Three cell types at equal counts: still reducible.
  eq <- .single_slide_fixture(nct = 3L, counts = c(80, 80, 80))
  cts3 <- c("CellTypeA", "CellTypeB", "CellTypeC")
  expect_equal(resolve(fit(eq, "sumcor"), cts3)$objective, "sumcov")

  # Three at unequal counts: the criteria genuinely differ, so the null has to
  # be re-optimized under SUMCOR, and it needs the permutation-invariant Gram
  # matrices to do it.
  un <- .single_slide_fixture(nct = 3L, counts = c(150, 60, 40))
  resolved <- resolve(fit(un, "sumcor"), cts3)
  expect_equal(resolved$objective, "sumcor")
  expect_equal(resolved$slideWeight, "equal")
  expect_named(resolved$grams, cts3)
  expect_equal(as.integer(resolved$n_cells), c(150L, 60L, 40L))

  # Tests restricted to at most two cell types can only meet the reducible
  # case; the guard is there in case that restriction is ever relaxed.
  expect_equal(
    suppressMessages(CoPro:::.resolvePermutationObjective(
      fit(two, "sumcor"), cts2, TRUE, supports_sumcor = FALSE
    ))$objective,
    "sumcov"
  )
  expect_error(
    CoPro:::.resolvePermutationObjective(
      fit(un, "sumcor"), cts3, TRUE, supports_sumcor = FALSE
    ),
    "different maximizers"
  )
})

test_that("the SUMCOR permutation null is only used where it is needed", {
  # The permutation-invariance argument: a within-slide label permutation
  # permutes rows of X_i, leaving G_i = X_i'X_i untouched, so the SUMCOR
  # denominators can be built once and reused for every draw.
  obj <- .single_slide_fixture(nct = 3L, counts = c(150, 60, 40))
  fit <- suppressMessages(suppressWarnings(
    runSkrCCA(obj, nCC = 1, objective = "sumcor")
  ))
  cts <- c("CellTypeA", "CellTypeB", "CellTypeC")
  resolved <- suppressMessages(
    CoPro:::.resolvePermutationObjective(fit, cts, TRUE)
  )

  X <- CoPro:::.getAllPCMats(allPCs = fit@pcaGlobal, scalePCs = TRUE)[cts]
  set.seed(4)
  for (ct in cts) {
    shuffled <- X[[ct]][sample(nrow(X[[ct]])), , drop = FALSE]
    expect_equal(crossprod(shuffled), resolved$grams[[ct]], tolerance = 1e-10,
                 info = ct)
  }
})

# ------------------------------- 9. provenance and gene-space arg forwarding ---

test_that("getCCAObjective() reports the criterion a gene-space fit used", {
  # runGeneSpaceCCA() defaults to "sumcor". Before it recorded its own
  # provenance, the reader's no-record fallback answered "sumcov" -- the
  # opposite of the truth -- for exactly the call NEWS tells users to inspect.
  obj <- .sumcor_fixture(nct = 2L, n_slides = 2L)
  fit <- suppressMessages(suppressWarnings(
    runGeneSpaceCCA(obj, sigma = 0.1, nCC = 1, max_iter = 100, tol = 1e-4,
                    verbose = FALSE)
  ))
  rec <- getCCAObjective(fit)
  expect_equal(rec$space, "gene")
  expect_equal(rec$objective, "sumcor")
  expect_equal(rec$sweep, "gauss-seidel")

  fit_cov <- suppressMessages(suppressWarnings(
    runGeneSpaceCCA(obj, sigma = 0.1, nCC = 1, max_iter = 100, tol = 1e-4,
                    verbose = FALSE, objective = "sumcov")
  ))
  expect_equal(getCCAObjective(fit_cov)$objective, "sumcov")
})

test_that("a gene-space fit overwrites an earlier PCA-space provenance record", {
  # Weights are merged into @skrCCAOut rather than replacing it, so without an
  # explicit write the earlier PCA record would survive and describe a fit that
  # is not the most recent one.
  obj <- .sumcor_fixture(nct = 2L, n_slides = 2L)
  pca_fit <- suppressMessages(suppressWarnings(
    runSkrCCA(obj, nCC = 1, objective = "sumcor")
  ))
  expect_equal(getCCAObjective(pca_fit)$space, "pca")

  both <- suppressMessages(suppressWarnings(
    runGeneSpaceCCA(pca_fit, sigma = 0.1, nCC = 1, max_iter = 100, tol = 1e-4,
                    verbose = FALSE)
  ))
  expect_equal(getCCAObjective(both)$space, "gene")
})

test_that("space = 'gene' rejects PC-space-only arguments instead of dropping them", {
  # These are named formals of runSkrCCA(), so they never reach `...`; silently
  # discarding one makes a dropped transfer look like a transfer that ran.
  obj <- .sumcor_fixture(nct = 2L, n_slides = 2L)
  for (arg in list(list(transferred_weight_1 = list(a = 1)),
                   list(scalePCs = FALSE),
                   list(step_size = 0.5),
                   list(n_cores = 2),
                   list(slideWeight = "equal"))) {
    expect_error(
      do.call(runSkrCCA, c(list(obj, space = "gene", nCC = 1,
                                sigmaChoice = 0.1), arg)),
      "does not use",
      info = names(arg)
    )
  }
})

test_that("space = 'gene' forwards tol, maxIter and minCellsPerSlide", {
  # Forwarded only when supplied: the two entry points have different defaults
  # (tol 1e-5 vs 1e-6, maxIter 200 vs 3000), so forwarding a runSkrCCA default
  # would silently change the gene-space result.
  # Gene space initializes from rnorm(), so every comparison here fixes the seed
  # first -- otherwise the axis signs alone would differ between runs.
  obj <- .sumcor_fixture(nct = 2L, n_slides = 2L)
  key <- "gscca_sigma_0.1"
  run <- function(...) {
    set.seed(404)
    suppressMessages(suppressWarnings(runSkrCCA(obj, space = "gene", nCC = 1,
                                                sigmaChoice = 0.1, ...)))
  }
  loose <- run(maxIter = 300, tol = 1e-1)
  tight <- run(maxIter = 300, tol = 1e-12)
  expect_false(identical(loose@skrCCAOut[[key]], tight@skrCCAOut[[key]]))

  # Not supplying an argument must reproduce runGeneSpaceCCA()'s own defaults
  # exactly, rather than quietly substituting runSkrCCA()'s -- tol 1e-5 vs 1e-6,
  # maxIter 200 vs 3000, and (the one that changes the criterion, not just the
  # precision) objective "sumcov" vs "sumcor".
  via_skr <- run()
  set.seed(404)
  direct <- suppressMessages(suppressWarnings(
    runGeneSpaceCCA(obj, sigma = 0.1, nCC = 1, verbose = FALSE)
  ))
  expect_equal(via_skr@skrCCAOut[[key]], direct@skrCCAOut[[key]])
  expect_equal(getCCAObjective(via_skr)$objective, "sumcor")

  # An explicit objective does travel.
  expect_equal(getCCAObjective(run(objective = "sumcov"))$objective, "sumcov")
})

# --------------- 10. transferred axes, damping, and the permutation Grams ---

test_that("the one-slide SUMCOR shortcut keeps a transferred first axis", {
  # runSkrCCA(transferred_weight_1 = ...) supplies CC1 from somewhere else --
  # another slide, another study -- and the later axes are only interpretable
  # conditioned on it. The reducible one-slide shortcut used to run the exact
  # SUMCOV solvers over every axis and hand back its own CC1 instead, which is
  # a different component silently wearing the transferred one's label.
  obj <- .single_slide_fixture(nct = 2L, nPCA = 6L)
  cts <- c("CellTypeA", "CellTypeB")
  transferred <- setNames(lapply(cts, function(ct) {
    v <- matrix(0, nrow = 6L, ncol = 1L)
    v[2L, 1L] <- 1                       # a direction the solver would not pick
    v
  }), cts)

  for (objective in c("sumcov", "sumcor")) {
    fit <- suppressMessages(suppressWarnings(runSkrCCA(
      obj, nCC = 2, objective = objective,
      transferred_weight_1 = transferred
    )))
    w <- fit@skrCCAOut[[1L]]
    for (ct in cts) {
      expect_equal(abs(as.numeric(w[[ct]][, 1L])),
                   abs(as.numeric(transferred[[ct]][, 1L])),
                   tolerance = 1e-8,
                   info = paste(objective, ct))
    }
    # CC2 still has to be a genuine second axis, not a copy of what was given.
    expect_gt(max(vapply(cts, function(ct) {
      min(abs(as.numeric(w[[ct]][, 2L]) - as.numeric(transferred[[ct]][, 1L])))
    }, numeric(1))), 0)
  }
})

test_that("step_size reaches every axis of the SUMCOV warm start", {
  # In the one-slide reducible shortcut the warm start *is* the returned
  # result, so damping that only touched CC1 left CC2+ undamped -- exactly the
  # axes a user reaching for step_size is usually trying to stabilize.
  ops <- .synthetic_one_slide_ops(counts = c(80, 80, 80), nPC = 4L, seed = 21L)
  ops_w <- CoPro:::.whitenSlideOperators(ops, NULL)

  iters <- function(step_size) {
    msg <- capture.output(invisible(suppressWarnings(
      CoPro:::.sumcorWarmStart(ops_w, ops$cell_types, NULL, nCC = 3L,
                               step_size = step_size)
    )), type = "message")
    as.integer(sub(" iterations", "", regmatches(
      msg, regexpr("[0-9]+ iterations", msg)
    )))
  }

  undamped <- iters(1)
  damped <- iters(0.2)
  expect_length(undamped, 3L)
  expect_length(damped, 3L)
  # A shorter step is a slower walk to the same fixed point. Every axis must
  # show it, not just the first.
  expect_true(all(damped > undamped))
})

test_that("a SUMCOR null refits each draw with that draw's own Gram matrices", {
  # SUMCOR divides by sigma_i = sqrt(w_i' G_i w_i). Reordering rows leaves
  # G_i = X_i'X_i alone, so a bijective null may reuse it -- but the default
  # "bin" null resamples cells and "pc" shuffles each PC column independently,
  # and both move G_i. Reusing the unpermuted one there fits the null draws
  # with a numerator from the permuted data and a denominator from the
  # observed data.
  obj <- .single_slide_fixture(nct = 3L, nPCA = 5L, counts = c(150, 60, 40))
  cts <- c("CellTypeA", "CellTypeB", "CellTypeC")
  fit <- suppressMessages(suppressWarnings(
    runSkrCCA(obj, nCC = 1, objective = "sumcor")
  ))
  fit <- suppressMessages(computeNormalizedCorrelation(fit))
  expect_equal(
    suppressMessages(CoPro:::.resolvePermutationObjective(fit, cts, TRUE))$objective,
    "sumcor"
  )

  # Rebuild draw 1 by hand and fit it both ways, then ask which one the
  # package actually produced.
  refit <- function(res, grams_from) {
    PCmats <- CoPro:::.getAllPCMats(allPCs = res@pcaGlobal,
                                    scalePCs = res@scalePCs)
    spec <- CoPro:::.permutationDrawSpec(res@cellPermu, cts, 1L)
    PC_draw <- CoPro:::.applyPermutationSpec(PCmats, spec)[cts]
    plan <- CoPro:::.buildYPlan(
      PCmats = PCmats, flat_kernels = res@kernelMatrices,
      sigma = res@sigmaValueChoice, cts = cts,
      fixed = CoPro:::.fixedPermutationTypes(res@cellPermu, cts)
    )
    src <- if (identical(grams_from, "draw")) PC_draw else PCmats
    suppressMessages(suppressWarnings(CoPro:::.fitSumcorPermutedAxes(
      Y_resi = CoPro:::.yResiFromPlan(plan, PC_draw),
      grams = setNames(lapply(cts, function(ct) crossprod(src[[ct]])), cts),
      n_cells = setNames(vapply(cts, function(ct) nrow(PC_draw[[ct]]),
                                integer(1)), cts),
      cts = cts, nCC = 1L,
      sdev2_list = CoPro:::.permutationSdev2(res, cts),
      slideWeight = "equal"
    )))
  }
  axis_gap <- function(a, b) {
    max(vapply(cts, function(ct) {
      max(abs(abs(as.numeric(a[[ct]][, 1L])) - abs(as.numeric(b[[ct]][, 1L]))))
    }, numeric(1)))
  }

  for (method in c("bin", "pc")) {
    res <- suppressMessages(suppressWarnings(
      runSkrCCAPermu(fit, nPermu = 3, permu_method = method, verbose = FALSE)
    ))
    got <- res@skrCCAPermuOut[[1L]]
    expect_lt(axis_gap(got, refit(res, "draw")), 1e-4)
    # And the two really are different fits here, so the check has teeth.
    expect_gt(axis_gap(refit(res, "draw"), refit(res, "observed")), 1e-3)
  }

  # A bijective null is genuinely invariant, so nothing changes for "global".
  res <- suppressMessages(suppressWarnings(
    runSkrCCAPermu(fit, nPermu = 3, permu_method = "global", verbose = FALSE)
  ))
  expect_lt(axis_gap(refit(res, "draw"), refit(res, "observed")), 1e-4)
  expect_lt(axis_gap(res@skrCCAPermuOut[[1L]], refit(res, "draw")), 1e-4)
})

test_that("cell-level permutation refuses a gene-space fit outright", {
  # Unreachable through the public API today: runGeneSpaceCCA() requires
  # CoProMulti, and every permutation entry point refuses CoProMulti first.
  # The two guards are independent, though, and this one must not lean on the
  # other -- the resolver reads @pcaGlobal and re-optimizes with the PC-space
  # solvers, which is the wrong feature space for gene-space weights whatever
  # criterion they were fitted under.
  obj <- .single_slide_fixture(nct = 2L, nPCA = 6L)
  fit <- suppressMessages(suppressWarnings(runSkrCCA(obj, nCC = 1)))
  attr(fit@skrCCAOut, "ccaObjective") <- list(
    space = "gene", objective = "sumcor", requested = "sumcor",
    slideWeight = "equal", sweep = NA_character_, slides = NA_character_,
    droppedSlides = character(0)
  )
  expect_error(
    CoPro:::.resolvePermutationObjective(fit, c("CellTypeA", "CellTypeB"), TRUE),
    "gene space"
  )
})
