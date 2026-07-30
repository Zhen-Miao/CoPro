# Tests for the PCA-space SUMCOR objective (objective = "sumcor").
#
# The claims being pinned here, in order:
#   1. With one slide, sumcor and sumcov are the same problem.
#   2. SUMCOV is the sigma-weighted special case of the SUMCOR family.
#   3. sumcor is invariant to per-slide score rescaling; sumcov is not.
#   4. The iteration is an ascent method and never returns worse than its
#      sumcov warm start.
#   5. It is deterministic (no RNG).
#   6. scalePCs stays a pure reparametrization.
#   7. Thin slides are dropped under sumcor, reported under sumcov.
#   8. Argument validation.

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

.single_slide_fixture <- function(nct = 2L, nPCA = 8L, seed = 11L) {
  obj <- create_test_copro_single(
    n_cells = 210, n_genes = 40, n_cell_types = max(nct, 2L), seed = seed
  )
  cts <- paste0("CellType", LETTERS[seq_len(nct)])
  obj <- subsetData(obj, cellTypesOfInterest = cts)
  obj <- computePCA(obj, nPCA = nPCA)
  obj <- computeDistance(obj, distType = "Euclidean2D",
                         normalizeDistance = TRUE, verbose = FALSE)
  computeKernelMatrix(obj, sigmaValues = 0.1, verbose = FALSE,
                      normalizeDistance = TRUE)
}

test_that("with one slide sumcor and sumcov give the same weights", {
  # The norm constraint IS the unit-variance constraint for whitened PCs, and
  # SUMCOR is scale invariant, so it attains its maximum on that constraint set
  # where every denominator is 1 and it reduces to SUMCOV. Holds for >2 cell
  # types and under both metrics. Checked on a genuine single-slide CoPro (the
  # is_multi = FALSE branch) and on a one-slide CoProMulti (is_multi = TRUE with
  # S = 1), which reach the short circuit by different routes.
  for (kind in c("single", "multi")) {
    for (nct in c(1L, 2L, 3L)) {
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
})

test_that("a single-slide object reports that sumcor reduced to sumcov", {
  obj <- .single_slide_fixture(nct = 2L)
  fit <- NULL
  suppressWarnings(expect_message(
    fit <- runSkrCCA(obj, nCC = 1, objective = "sumcor"),
    "same optimization problem"
  ))
  # The record must say what was actually run, and also what was asked for.
  expect_equal(getCCAObjective(fit)$objective, "sumcov")
  expect_equal(getCCAObjective(fit)$requested, "sumcor")
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

test_that("the sumcor objective is non-decreasing across sweeps", {
  obj <- .sumcor_fixture(nct = 3L)
  cts <- c("CellTypeA", "CellTypeB", "CellTypeC")
  parts <- .sumcor_parts(obj, cts)
  warm <- CoPro:::.sumcorWarmStart(parts$ops, cts, parts$sdev2, 1L)
  w <- setNames(lapply(cts, function(ct) warm[[ct]][, 1, drop = FALSE]), cts)

  trajectory <- numeric(0)
  for (k in seq_len(12)) {
    fit <- suppressWarnings(CoPro:::.sumcorIterate(
      w_init = w, ops = parts$ops, slideWeight = "size",
      sdev2_list = parts$sdev2, max_iter = k, tol = 1e-14
    ))
    trajectory <- c(trajectory, fit$objective)
  }
  expect_true(all(diff(trajectory) >= -1e-12),
              info = paste(sprintf("%.10f", trajectory), collapse = " "))
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

test_that("getCCAObjective records what was run and defaults to sumcov", {
  obj <- .sumcor_fixture(nct = 2L)

  fit <- suppressMessages(suppressWarnings(
    runSkrCCA(obj, nCC = 1, objective = "sumcor", slideWeight = "equal")
  ))
  record <- getCCAObjective(fit)
  expect_equal(record$space, "pca")
  expect_equal(record$objective, "sumcor")
  expect_equal(record$slideWeight, "equal")
  expect_equal(length(record$slides), length(getSlideList(obj)))

  # An object predating the provenance attribute must read back as sumcov,
  # which is what such objects were computed under.
  stripped <- fit
  attr(stripped@skrCCAOut, "ccaObjective") <- NULL
  expect_equal(getCCAObjective(stripped)$objective, "sumcov")
})

test_that("permutation tests refuse sumcor weights", {
  obj <- .sumcor_fixture(nct = 2L)
  fit <- suppressMessages(suppressWarnings(
    runSkrCCA(obj, nCC = 1, objective = "sumcor")
  ))
  # CoProMulti is rejected earlier for cell-level permutation, so exercise the
  # guard directly -- it is the piece that must not silently mix criteria.
  expect_error(
    CoPro:::.rejectSumcorForPermutation(fit, "runSkrCCAPermu()"),
    "mixes criteria"
  )
  fit_cov <- suppressMessages(suppressWarnings(
    runSkrCCA(obj, nCC = 1, objective = "sumcov")
  ))
  expect_null(CoPro:::.rejectSumcorForPermutation(fit_cov))
})
