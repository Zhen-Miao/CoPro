# =============================================================================
# Permutation operator plan
# -----------------------------------------------------------------------------
# Every permutation draw needs the PC-space operator Y_ij = X_i' K_ij X_j, and
# forming it from the sparse kernel costs O(nnz(K_ij) * nPC) per draw. That
# product is the whole cost of a permutation test: everything downstream of it
# (SVD, eigendecomposition, deflation, normalized correlation) works on
# nPC x nPC matrices.
#
# When a cell type is held fixed across all draws -- which is exactly what
# permu_which = "second_only" (the default) and "first_only" do -- one side of
# the triple product never changes, so the kernel can be applied to it ONCE and
# every draw reduces to a small dense product:
#
#   i fixed:  Y_ij = (K_ij' X_i)' X_j[perm, ]
#   j fixed:  Y_ij = X_i[perm, ]' (K_ij X_j)
#
# Per-draw cost falls from O(nnz * nPC) to O(n * nPC^2), a factor of about
# nnz(K) / (n * nPC). The identity is exact -- reordering the rows of X_j is
# the same reordering applied to the columns of K_ij X_j -- so the null
# distribution and the p-values it produces are unchanged. Pairs with both
# sides permuted (permu_which = "both", or the non-fixed pairs of a
# three-type run) fall back to the original sparse product.
#
# The plan also carries the per-type Gram matrices X' X. A row permutation
# leaves X' X invariant, so the score norms in the normalized-correlation
# denominator can be read out of a precomputed nPC x nPC matrix instead of
# forming X[perm, ] %*% w for every draw. This does NOT hold for
# permu_method = "pc", which shuffles each PC column independently and so
# changes the off-diagonal entries; those types keep the direct calculation.
#
# Pass `factorize = FALSE` to a permutation entry point to disable the
# factorization and route every pair through the original sparse product.
# Results are equivalent; the switch exists so the two paths can be compared
# directly. It used to be the global `options(CoPro.factorizePermutation=)`,
# which still supplies the argument's default so old scripts keep working.
#
# -----------------------------------------------------------------------------
# Compact permutation storage
#
# A draw is an integer index per cell, so an explicit index matrix costs
# n * nPermu * 4 bytes per cell type -- 799 MB at 200k cells and 999 draws,
# persisted into the object and into every saveRDS. For a cell type that is
# *held fixed* (which is every type but one under the default
# permu_which = "second_only") that matrix stores the integers 1..n, 999 times,
# and carries no information at all. `.identityPermutation()` replaces it with
# a marker, which also lets `.applyPermutationSpec()` skip the per-draw row
# subset entirely rather than copying an n x nPC matrix back to itself.
#
# The genuinely permuted side can be compressed the same way -- store one seed
# per draw and re-draw on demand -- but that changes which permutations are
# drawn, so a re-run of a saved analysis would move its p-values within Monte
# Carlo error. That is gated behind `compactPermutation = TRUE` on the
# permutation entry points and off by default. The identity marker changes no
# number at all.
# =============================================================================

#' Default for the `compactPermutation` argument
#'
#' Off by default: enabling it changes the drawn permutations, so p-values from
#' a re-run will differ from previously saved ones within Monte Carlo error.
#' The former global `options(CoPro.compactPermutation=)` is still read here so
#' that scripts which set it keep their behavior; an explicit argument wins.
#' @noRd
.defaultCompactPermutation <- function() {
  isTRUE(getOption("CoPro.compactPermutation", FALSE))
}

#' Default for the `factorize` argument
#'
#' On by default. The former global `options(CoPro.factorizePermutation=)` is
#' still read here so that scripts which set it keep their behavior; an
#' explicit argument wins.
#' @noRd
.defaultFactorizePermutation <- function() {
  isTRUE(getOption("CoPro.factorizePermutation", TRUE))
}

#' Compact marker for a cell type that is held fixed across all draws
#' @param n_cell Number of cells of this type.
#' @noRd
.identityPermutation <- function(n_cell) {
  list(type = "identity", n_cell = as.integer(n_cell))
}

#' Is this `cellPermu` entry the compact identity marker?
#' @noRd
.isIdentityPermuEntry <- function(entry) {
  is.list(entry) && identical(entry$type, "identity")
}

#' Is a stored permutation the identity for every draw?
#'
#' `.getCellPermu()` records a held-out cell type with `.identityPermutation()`.
#' Older objects (and `options(CoPro.compactPermutation = FALSE)` runs of other
#' schemes) store an explicit `replicate(nPermu, 1:n)` matrix, so the matrix
#' branch is retained. Detecting identity here -- rather than re-deriving it
#' from `permu_which` -- keeps the factorization correct for any future
#' permutation scheme.
#'
#' @param perm One element of a `cellPermu` list.
#' @return `TRUE` when `perm` is the identity marker, or an index matrix whose
#'   every column is `seq_len(nrow(perm))`.
#' @noRd
.isIdentityPermutation <- function(perm) {
  if (.isIdentityPermuEntry(perm)) return(TRUE)
  if (!is.matrix(perm) || nrow(perm) == 0L || ncol(perm) == 0L) return(FALSE)
  idx <- seq_len(nrow(perm))
  # Compare one column at a time. `all(perm == idx)` would allocate a logical
  # the size of the whole index matrix -- another 799 MB at 200k x 999 -- and a
  # permuted type almost always fails on its first column anyway.
  for (col in seq_len(ncol(perm))) {
    if (!all(perm[, col] == idx)) return(FALSE)
  }
  TRUE
}

#' Is a stored permutation a genuine bijection on the cells in every draw?
#'
#' The `"bin"` null draws a spatially matched *resample* -- a cell can appear
#' twice in a draw and another not at all -- so `X[draw, ]' X[draw, ]` is not
#' `X' X` there. Only bijections license the Gram shortcut in
#' `.permutationGrams()`; the operator factorization itself is valid either
#' way, because it never assumes anything about the permuted side.
#'
#' @param perm One element of a `cellPermu` list.
#' @return `TRUE` when every draw uses each row of the data exactly once.
#' @noRd
.isRowPermutationMatrix <- function(perm) {
  # Compact entries are classified by construction: identity and a reseeded
  # `sample.int()` are bijections; a "bin" resample and a "pc" column shuffle
  # are not.
  if (is.list(perm)) {
    return(identical(perm$type, "identity") ||
             identical(perm$type, "global_seed"))
  }
  if (!is.matrix(perm) || nrow(perm) == 0L || ncol(perm) == 0L) return(FALSE)
  n <- nrow(perm)
  for (col in seq_len(ncol(perm))) {
    # tabulate() drops out-of-range and NA indices, so a short or repeated
    # column always shows up as a count other than one.
    if (!all(tabulate(perm[, col], nbins = n) == 1L)) return(FALSE)
  }
  TRUE
}

#' Which cell types are held fixed across all permutation draws?
#' @param cell_permu A `cellPermu` list.
#' @param cts Cell types of interest.
#' @return Named logical vector over `cts`.
#' @noRd
.fixedPermutationTypes <- function(cell_permu, cts) {
  stats::setNames(vapply(cts, function(ct) {
    .isIdentityPermutation(cell_permu[[ct]])
  }, logical(1)), cts)
}

#' Apply a kernel to a dense matrix, returning a base matrix
#'
#' Handles encoded float32 kernels, `dgCMatrix`, and `dsCMatrix` uniformly and
#' returns an ordinary dense matrix so the per-draw products stay in BLAS.
#' @noRd
.applyKernelDense <- function(K, X) {
  as.matrix(.float32KernelMatMult(K, X))
}

#' Cell-type pairs a permutation statistic is evaluated over
#'
#' One cell type is the within-type (self) problem and pairs with itself;
#' `combn()` cannot express that, it errors with `n < m`. Every consumer of a
#' pair list must go through this, or the within-type path dies on entry --
#' which is exactly what `computeNormalizedCorrelationPermu()` used to do.
#' @noRd
.permutationPairTypes <- function(cts) {
  if (length(cts) == 1L) {
    matrix(c(cts, cts), nrow = 2L, ncol = 1L)
  } else {
    utils::combn(cts, 2L)
  }
}

#' Build the per-pair plan for evaluating Y across permutation draws
#'
#' @param PCmats Named list of unpermuted cell-by-PC matrices.
#' @param flat_kernels Flat kernel list.
#' @param sigma Kernel bandwidth.
#' @param cts Cell types of interest.
#' @param fixed Named logical vector from `.fixedPermutationTypes()`.
#' @param slide Slide ID, or `NULL` for a single slide.
#' @param factorize Use the fixed-side factorization? `FALSE` routes every pair
#'   through the original sparse product.
#' @return A list with `ops` (nested by cell-type pair), `pairs`, `cts` and
#'   `fixed`, consumable by `.yResiFromPlan()`.
#' @noRd
.buildYPlan <- function(PCmats, flat_kernels, sigma, cts, fixed,
                        slide = NULL,
                        factorize = .defaultFactorizePermutation()) {
  if (!isTRUE(factorize)) {
    fixed <- stats::setNames(rep(FALSE, length(cts)), cts)
  }

  ops <- stats::setNames(vector("list", length(cts)), cts)
  for (ct in cts) {
    ops[[ct]] <- stats::setNames(vector("list", length(cts)), cts)
  }

  pairs <- .permutationPairTypes(cts)

  for (pp in seq_len(ncol(pairs))) {
    ct_i <- pairs[1L, pp]
    ct_j <- pairs[2L, pp]
    self_pair <- identical(ct_i, ct_j)
    K <- get_kernel_matrix_flat(flat_kernels, sigma, ct_i, ct_j, slide)

    op <- if (fixed[[ct_i]] && fixed[[ct_j]]) {
      # Nothing on either side moves: Y is the same for every draw.
      list(mode = "const",
           Y = .kernelXKY(PCmats[[ct_i]], K, PCmats[[ct_j]]))
    } else if (self_pair) {
      # A within-type kernel permutes both of its sides together.
      list(mode = "full", K = K)
    } else if (fixed[[ct_i]]) {
      list(mode = "left", M = .applyKernelDense(t(K), PCmats[[ct_i]]))
    } else if (fixed[[ct_j]]) {
      list(mode = "right", L = .applyKernelDense(K, PCmats[[ct_j]]))
    } else {
      list(mode = "full", K = K)
    }

    op$self_pair <- self_pair
    ops[[ct_i]][[ct_j]] <- op
  }

  list(ops = ops, pairs = pairs, cts = cts, fixed = fixed)
}

#' Evaluate the PC-space operators for one permutation draw
#'
#' Produces exactly the structure `compute_Y_resi()` returns, so every
#' downstream solver is unchanged.
#' @param plan A plan from `.buildYPlan()`.
#' @param PCmats_local Permuted cell-by-PC matrices for this draw.
#' @return Nested `Y_resi` list.
#' @noRd
.yResiFromPlan <- function(plan, PCmats_local) {
  cts <- plan$cts
  Y_resi <- stats::setNames(vector("list", length(cts)), cts)
  for (ct in cts) {
    Y_resi[[ct]] <- stats::setNames(vector("list", length(cts)), cts)
  }

  for (pp in seq_len(ncol(plan$pairs))) {
    ct_i <- plan$pairs[1L, pp]
    ct_j <- plan$pairs[2L, pp]
    op <- plan$ops[[ct_i]][[ct_j]]

    Y_ij <- switch(
      op$mode,
      const = op$Y,
      left  = crossprod(op$M, PCmats_local[[ct_j]]),
      right = crossprod(PCmats_local[[ct_i]], op$L),
      full  = .kernelXKY(PCmats_local[[ct_i]], op$K, PCmats_local[[ct_j]]),
      stop("Unknown permutation operator mode: ", op$mode)
    )

    if (isTRUE(op$self_pair)) {
      # Matches compute_symmetric_Y(): the quadratic form depends only on the
      # symmetric part, and a user-supplied kernel may be slightly asymmetric.
      Y_resi[[ct_i]][[ct_j]] <- (Y_ij + t(Y_ij)) * 0.5
    } else {
      Y_resi[[ct_i]][[ct_j]] <- Y_ij
      Y_resi[[ct_j]][[ct_i]] <- t(Y_ij)
    }
  }

  Y_resi
}

#' Precompute per-type Gram matrices for permutation-invariant score norms
#'
#' `sum((X[perm, ] w)^2) = w' (X' X) w` when `perm` is a bijection, because
#' reordering rows leaves `X' X` unchanged. Two of the nulls break that
#' assumption and get `NULL` (keeping the direct calculation): `"pc"` shuffles
#' each PC column independently, which changes the off-diagonal entries, and
#' `"bin"` draws a spatially matched resample rather than a permutation.
#' @param factorize Cache the Gram matrices? `FALSE` returns all-`NULL` so every
#'   draw recomputes `X[perm, ] w` directly.
#' @noRd
.permutationGrams <- function(PCmats, cell_permu, cts,
                              factorize = .defaultFactorizePermutation()) {
  if (!isTRUE(factorize)) {
    return(stats::setNames(vector("list", length(cts)), cts))
  }
  stats::setNames(lapply(cts, function(ct) {
    if (!.isRowPermutationMatrix(cell_permu[[ct]])) return(NULL)
    crossprod(PCmats[[ct]])
  }), cts)
}

#' Fill in the Gram matrices a single draw could not inherit
#'
#' `.permutationGrams()` leaves a `NULL` wherever the null is not a bijection on
#' cells, meaning `X' X` genuinely changed for that cell type and the cached
#' value does not describe this draw. The SUMCOR optimizer divides by
#' \eqn{\sigma_i = \sqrt{w_i' G_i w_i}}, so it needs a real matrix for every
#' type: recompute exactly the entries that moved, from the draw's own permuted
#' PC matrix, and reuse the cached one everywhere else.
#'
#' Recomputing is the whole cost of a `"bin"` or `"pc"` SUMCOR null, and it is
#' `nPC^2 n` per affected type per draw -- the same order as the `Y` operator
#' the draw already builds.
#'
#' @param grams Per-type Gram matrices from `.permutationGrams()`, with `NULL`
#'   where the draw invalidates the cached value.
#' @param PCmats_local This draw's permuted PC matrices.
#' @param cts Cell types of interest.
#' @return A named list with a Gram matrix for every entry of `cts`.
#' @noRd
.drawGrams <- function(grams, PCmats_local, cts) {
  stats::setNames(lapply(cts, function(ct) {
    if (!is.null(grams[[ct]])) return(grams[[ct]])
    crossprod(PCmats_local[[ct]])
  }), cts)
}

#' Drop the kernel from a scoring-info list before sending it to a worker
#'
#' `.compute_ncorr_quick()` only touches `kernel_info$K` when no precomputed
#' `Y_resi` is available. The permutation paths always supply one, so the
#' kernel -- which can be gigabytes -- must not ride along to every worker.
#' @noRd
.slimKernelInfo <- function(kernel_info) {
  list(K = NULL, norm_K12 = kernel_info$norm_K12)
}

#' Is this `cellPermu` entry a `"pc"` null rather than an index matrix?
#'
#' The two permutation representations are discriminated in three places
#' (draw-spec extraction, worker chunking, and the plan itself). Keeping the
#' test in one function means a future third representation cannot satisfy one
#' copy of the predicate and silently miss another.
#' @noRd
.isPcPermuteEntry <- function(entry) {
  is.list(entry) && identical(entry$type, "pc_permute")
}

#' Shuffle values within each PC column (the DIALOGUE-style "pc" null)
#' @noRd
.permutePCMatrix <- function(pc_mat, seed) {
  rng_state <- .captureRNGState()
  on.exit(.restoreRNGState(rng_state), add = TRUE)
  set.seed(seed)
  permuted <- apply(pc_mat, 2L, function(x) sample(x, length(x)))
  rownames(permuted) <- rownames(pc_mat)
  permuted
}

#' Extract the permutation instructions for a single draw
#'
#' Either an integer index vector per cell type, or a one-element list naming
#' the operation to perform: the identity, or the stored seed for a `"pc"`,
#' `"global"` or `"bin"` draw under compact storage.
#' @noRd
.permutationDrawSpec <- function(cell_permu, cts, tt) {
  stats::setNames(lapply(cts, function(ct) {
    entry <- cell_permu[[ct]]
    if (.isIdentityPermuEntry(entry)) {
      list(identity = TRUE)
    } else if (.isPcPermuteEntry(entry)) {
      list(pc_seed = entry$seeds[[tt]])
    } else if (is.list(entry) && identical(entry$type, "global_seed")) {
      .drawGlobalPermutation(entry$n_cell, entry$seeds[[tt]])
    } else if (is.list(entry) && identical(entry$type, "bin_seed")) {
      .drawBinPermutation(entry, entry$seeds[[tt]])
    } else {
      entry[, tt]
    }
  }), cts)
}

#' Re-draw a compact `"global"` permutation from its stored seed
#' @noRd
.drawGlobalPermutation <- function(n_cell, seed) {
  rng_state <- .captureRNGState()
  on.exit(.restoreRNGState(rng_state), add = TRUE)
  set.seed(seed)
  sample.int(n = n_cell, replace = FALSE)
}

#' Re-draw a compact `"bin"` resample from its stored seed
#' @noRd
.drawBinPermutation <- function(entry, seed) {
  rng_state <- .captureRNGState()
  on.exit(.restoreRNGState(rng_state), add = TRUE)
  set.seed(seed)
  .drawSpatialPermutation(entry$prepared,
                          match_quantile = entry$match_quantile)
}

#' Apply one draw's permutation to the PC matrices
#' @noRd
.applyPermutationSpec <- function(PCmats, spec) {
  for (ct in names(spec)) {
    entry <- spec[[ct]]
    # A draw spec carries either an integer index vector or a one-element list
    # naming an operation. Test for the field explicitly rather than for
    # "is a list".
    if (is.list(entry)) {
      if (isTRUE(entry$identity)) {
        # Subsetting by seq_len(n) would copy the whole n x nPC matrix back to
        # itself, once per draw per held-fixed cell type. Leave it alone.
        #
        # One visible consequence: `[` drops attributes, so the old identity
        # subset also stripped the "scaled:scale" attribute that
        # .getAllPCMats() leaves behind when scalePCs = TRUE, while returning
        # the matrix untouched keeps it. Nothing in the package reads that
        # attribute -- every consumer goes through crossprod(), %*% or
        # .kernelXKY(), which see only the values and dim -- and the end-to-end
        # baseline confirms every downstream number is bit-identical.
        next
      }
      if (!is.null(entry$pc_seed)) {
        PCmats[[ct]] <- .permutePCMatrix(PCmats[[ct]], entry$pc_seed)
        next
      }
    }
    PCmats[[ct]] <- PCmats[[ct]][entry, , drop = FALSE]
  }
  PCmats
}

#' Build the per-draw worker for the fixed-sigma permutation test
#'
#' A factory rather than an inline closure: the returned function's
#' environment is exactly this frame, so a PSOCK worker receives the plan and
#' the PC matrices and nothing else. An inline closure would carry the whole
#' CoPro object, kernels included.
#'
#' @param permu_grams Per-type Gram matrices for the SUMCOR denominators, with a
#'   `NULL` wherever this null changes `X' X` and the draw must recompute its
#'   own. Ignored unless `permu_objective` asks for a SUMCOR null.
#' @noRd
.makeSkrCCAPermuWorker <- function(PCmats, plan, cts, nCC, sdev2_list,
                                   maxIter, tol, permu_objective = NULL,
                                   permu_grams = NULL) {
  force(PCmats); force(plan); force(cts); force(nCC)
  force(sdev2_list); force(maxIter); force(tol); force(permu_objective)
  force(permu_grams)
  sumcor <- identical(permu_objective$objective, "sumcor")

  function(spec) {
    PCmats_local <- .applyPermutationSpec(PCmats, spec)[cts]
    Y_resi <- .yResiFromPlan(plan, PCmats_local)

    # Only reached with three or more cell types: fewer than that makes SUMCOR
    # and SUMCOV the same problem, which .resolvePermutationObjective() detects
    # and routes down the exact decompositions below.
    if (sumcor) {
      # The numerator Y is built from this draw's permuted PCs, so the
      # denominator has to be as well. Only a bijection on cells leaves X' X
      # alone; the default "bin" null resamples and "pc" shuffles each column
      # independently, and both move it.
      return(.fitSumcorPermutedAxes(
        Y_resi = Y_resi, grams = .drawGrams(permu_grams, PCmats_local, cts),
        n_cells = permu_objective$n_cells, cts = cts, nCC = nCC,
        sdev2_list = sdev2_list, slideWeight = permu_objective$slideWeight,
        maxIter = maxIter, tol = tol
      ))
    }

    # One and two cell types have exact decompositions; every requested axis
    # comes from a single factorization of the small operator.
    exact_weights <- .solveExactSumcovAxes(
      Y_resi, cts, nCC = nCC, sdev2_list = sdev2_list
    )
    if (!is.null(exact_weights)) return(exact_weights)

    cca_result <- optimize_bilinear(
      X_list = PCmats_local, flat_kernels = NULL, sigma = NULL,
      max_iter = maxIter, tol = tol, sdev2_list = sdev2_list,
      Y_resi = Y_resi
    )
    names(cca_result) <- cts
    if (nCC == 1L) return(cca_result)

    optimize_bilinear_n(
      X_list = PCmats_local, flat_kernels = NULL, sigma = NULL,
      w_list = cca_result, cellTypesOfInterest = cts, nCC = nCC,
      max_iter = maxIter, tol = tol, sdev2_list = sdev2_list,
      Y_resi = Y_resi
    )
  }
}

#' Build the per-draw worker for the fair-sigma permutation test
#'
#' Fits CC1 at every candidate sigma and keeps the best, so the null statistic
#' faces the same sigma selection as the observed one.
#' @noRd
.makeFairSigmaWorker <- function(PCmats, plans, cts, sigma_values, sdev2_list,
                                 kernel_info, grams, maxIter, tol) {
  force(PCmats); force(plans); force(cts); force(sigma_values)
  force(sdev2_list); force(kernel_info); force(grams)
  force(maxIter); force(tol)

  function(spec) {
    PCmats_local <- .applyPermutationSpec(PCmats, spec)[cts]

    best_ncorr <- -Inf
    best_sigma <- sigma_values[1]
    best_weights <- NULL
    errors <- character()

    for (si in seq_along(sigma_values)) {
      Y0 <- .yResiFromPlan(plans[[si]], PCmats_local)
      fit <- tryCatch(
        .fitConditionalAxis(
          PCmats = PCmats_local, flat_kernels = NULL,
          sigma = sigma_values[si], cts = cts, k_minus_1 = 0, Y_resi = Y0,
          kernel_info = kernel_info[[si]], sdev2_list = sdev2_list,
          grams = grams, maxIter = maxIter, tol = tol
        ),
        error = function(e) {
          errors <<- c(errors, conditionMessage(e))
          NULL
        }
      )
      if (is.null(fit)) next
      if (!is.finite(fit$ncorr)) {
        errors <- c(errors, "non-finite normalized correlation")
        next
      }
      if (fit$ncorr > best_ncorr) {
        best_ncorr <- fit$ncorr
        best_sigma <- sigma_values[si]
        best_weights <- fit$w
      }
    }

    list(
      weights = best_weights,
      sigma = best_sigma,
      ncorr = if (is.null(best_weights)) NA_real_ else best_ncorr,
      errors = errors
    )
  }
}

#' Build the per-draw worker for the conditional step-down permutation test
#'
#' For every axis k the fixed observed CC1..CC(k-1) directions are projected
#' out of the permuted operator, and the statistic is maximized over sigma.
#' @noRd
.makeConditionalWorker <- function(PCmats, plans, cts, nCC, sigma_values,
                                   sigma_names, obs_W, sdev2_list,
                                   kernel_info, grams, maxIter, tol) {
  force(PCmats); force(plans); force(cts); force(nCC)
  force(sigma_values); force(sigma_names); force(obs_W)
  force(sdev2_list); force(kernel_info); force(grams)
  force(maxIter); force(tol)

  function(spec) {
    PCmats_local <- .applyPermutationSpec(PCmats, spec)[cts]

    tryCatch({
      stat_k <- rep(-Inf, nCC)
      sig_k <- rep(sigma_values[1], nCC)

      for (si in seq_along(sigma_values)) {
        ## Y is shared by every axis at this sigma.
        Y0 <- .yResiFromPlan(plans[[si]], PCmats_local)
        for (k in seq_len(nCC)) {
          fit <- .fitConditionalAxis(
            PCmats = PCmats_local, flat_kernels = NULL,
            sigma = sigma_values[si], cts = cts,
            W_lower = obs_W[[sigma_names[si]]], k_minus_1 = k - 1,
            Y_resi = Y0, kernel_info = kernel_info[[si]],
            sdev2_list = sdev2_list, grams = grams,
            maxIter = maxIter, tol = tol
          )
          if (!is.finite(fit$ncorr)) {
            stop("non-finite normalized correlation for CC", k,
                 " at sigma ", sigma_values[si])
          }
          if (is.finite(fit$ncorr) && fit$ncorr > stat_k[k]) {
            stat_k[k] <- fit$ncorr
            sig_k[k] <- sigma_values[si]
          }
        }
      }

      list(stat = stat_k, sigma = sig_k)
    }, error = function(e) list(error = conditionMessage(e)))
  }
}

#' Run one worker's slice of the permutation draws
#'
#' Defined at package top level on purpose. `parallel` serializes a function
#' together with its enclosing environment, so a chunk runner defined inside
#' `.runPermutationDraws()` would ship that frame -- including the full
#' permutation index matrix -- to every worker. A namespace function is
#' serialized by reference instead, and only `payload` crosses the wire.
#' @noRd
.runPermutationChunk <- function(payload) {
  lapply(seq_along(payload$draws), function(i) {
    payload$worker(
      .permutationDrawSpec(payload$permu, payload$cts, i)
    )
  })
}

#' Evaluate a permutation worker over all draws, optionally in parallel
#'
#' Each PSOCK worker receives only its own columns of the permutation index
#' matrix and the (small) closure the worker was built from. Nothing here
#' captures the CoPro object or its kernels.
#'
#' @param cell_permu A `cellPermu` list.
#' @param cts Cell types of interest.
#' @param nPermu Number of draws.
#' @param worker A closure taking one draw spec and returning that draw's
#'   result. Build it with a factory so its environment stays small.
#' @param n_cores Worker count; 1 runs sequentially in this process.
#' @param verbose Whether to report worker setup and progress.
#' @param progress_every Emit a progress line every this many sequential
#'   draws, or `NULL` for no progress output.
#' @return List of per-draw results, in draw order.
#' @noRd
.runPermutationDraws <- function(cell_permu, cts, nPermu, worker,
                                 n_cores = 1L, verbose = TRUE,
                                 progress_every = NULL) {
  if (!is.numeric(n_cores) || length(n_cores) != 1L ||
      n_cores < 1L || n_cores != as.integer(n_cores)) {
    stop("n_cores must be a positive integer.")
  }

  run_sequentially <- function() {
    lapply(seq_len(nPermu), function(tt) {
      result <- worker(.permutationDrawSpec(cell_permu, cts, tt))
      if (verbose && !is.null(progress_every) &&
          (tt %% progress_every == 0L || tt == nPermu)) {
        cat(paste("  Completed", tt, "of", nPermu, "permutations\n"))
      }
      result
    })
  }

  if (n_cores == 1L || nPermu <= 1L) return(run_sequentially())

  copro_library <- .installedCoProLibrary()
  if (is.null(copro_library)) {
    warning("Parallel permutation requires an installed CoPro package ",
            "(devtools::load_all() is not sufficient); falling back to ",
            "sequential execution. Run devtools::install() to enable it.")
    return(run_sequentially())
  }

  n_workers <- min(as.integer(n_cores), nPermu)
  if (verbose) message("Using ", n_workers, " PSOCK permutation workers.")
  cluster <- parallel::makeCluster(n_workers)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterCall(cluster, function(lib) {
    .libPaths(c(lib, .libPaths()))
    suppressPackageStartupMessages(library(CoPro))
    NULL
  }, copro_library)

  chunks <- split(
    seq_len(nPermu),
    cut(seq_len(nPermu), n_workers, labels = FALSE)
  )
  payloads <- lapply(chunks, function(chunk) {
    list(
      draws = chunk,
      cts = cts,
      worker = worker,
      permu = stats::setNames(lapply(cts, function(ct) {
        entry <- cell_permu[[ct]]
        if (.isIdentityPermuEntry(entry)) {
          # Already O(1); nothing to slice, and nothing per-draw to ship.
          entry
        } else if (is.list(entry) && !is.null(entry$seeds)) {
          # Every seed-based representation ("pc", and the compact "global" and
          # "bin" entries) carries one seed per draw plus draw-invariant state.
          entry$seeds <- entry$seeds[chunk]
          entry
        } else {
          entry[, chunk, drop = FALSE]
        }
      }), cts)
    )
  })

  results <- parallel::clusterApply(cluster, payloads, .runPermutationChunk)
  unlist(results, recursive = FALSE, use.names = FALSE)
}
