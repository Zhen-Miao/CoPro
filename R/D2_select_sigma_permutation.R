#' Select the kernel bandwidth by studentized permutation (max-T)
#'
#' Chooses the spatial bandwidth \eqn{\sigma} by comparing the co-progression
#' statistic at each candidate bandwidth to *its own* permutation null, rather
#' than by taking the largest normalized correlation across bandwidths.
#'
#' @details
#' ## Why not just take the largest normalized correlation
#'
#' [computeNormalizedCorrelation()] divides the bilinear statistic
#' \eqn{T(\sigma) = a' K(\sigma) b} by a norm of the kernel, and
#' `object@@sigmaValueChoice` is the bandwidth where that ratio is largest. No
#' available denominator makes the *null level* of that ratio constant in
#' \eqn{\sigma}: the un-whitened \eqn{\|K_c\|_F} ignores within-type spatial
#' autocorrelation and drifts down with \eqn{\sigma}, while whitened variants
#' need a within-type correlation operator that the data do not pin down (the
#' principal components of one cell type do not share a single correlation
#' length). Comparing the ratio across \eqn{\sigma} therefore compares numbers
#' whose null expectations differ, and the argmax is biased toward whichever
#' bandwidth happens to have the highest floor.
#'
#' ## What this function does instead
#'
#' It measures the floor rather than modelling it. For each \eqn{\sigma} on the
#' grid it computes the observed statistic \eqn{T(\sigma)} at the cell scores
#' CoPro already fitted at that \eqn{\sigma}, then estimates the mean
#' \eqn{m(\sigma)} and the standard deviation \eqn{s(\sigma)} of
#' \eqn{T(\sigma)} under a spatial null and selects
#'
#' \deqn{\hat\sigma = \arg\max_\sigma\; \frac{T(\sigma) - m(\sigma)}{s(\sigma)}.}
#'
#' Because both the location and the scale are read off the null itself, the
#' studentized statistic \eqn{z(\sigma)} has the same null level at every
#' bandwidth by construction, and it carries no tuning constant.
#'
#' Centering is not cosmetic. Across two cell types \eqn{m(\sigma)} is zero in
#' expectation --- a wrap-around shift leaves every shifted cell equally likely
#' to land anywhere in the box, so each column of the kernel has the same mean
#' and a centered score vector annihilates it --- and subtracting the estimate
#' costs only Monte-Carlo noise. *Within* one cell type it is not zero: the
#' \eqn{w_{ii} = 0} convention takes \eqn{\sum_i a_i^2\,k(x_i, x_i + \delta)}
#' off every draw, and that term has a negative expectation which grows with
#' \eqn{\sigma}. On the package's toy object (\eqn{B = 4000}) the null mean ran
#' from \eqn{-0.06\,s(\sigma)} at \eqn{\sigma = 0.05} to \eqn{-0.36\,s(\sigma)}
#' at \eqn{\sigma = 0.2} --- a \eqn{\sigma}-dependent floor of exactly the kind
#' this function exists to remove, and one that dividing by \eqn{s(\sigma)}
#' alone leaves in place.
#'
#' ## The null, and why the same draws serve selection and inference
#'
#' The null is a toroidal (rigid wrap-around) shift of the coordinates. A rigid
#' shift preserves each cell type's own spatial autocorrelation --- the
#' structure that a plain label shuffle destroys, and whose destruction is what
#' makes shuffle-based nulls anti-conservative --- while removing the alignment
#' between the score fields. It assumes spatial stationarity and wrap-around,
#' which is false at tissue edges; see [runSkrCCAPermu()] for the bin-wise
#' alternative.
#'
#' One draw is one offset *per cell type*, reused everywhere that type appears.
#' With three or more cell types a type sits in several pairs at once, so
#' drawing its offset afresh per pair would put the same cells in two places
#' inside one draw: that row of statistics would come from no single
#' configuration, and a max-T reference built from such rows is not a null for
#' the scan. The first cell type anchors the frame, which costs no randomization
#' --- the offset *between* any two types is still uniform over the box --- and
#' leaves one fewer cloud cut by the wrap seam. The observed configuration is
#' the zero offset, so it is itself an admissible draw, which is what the
#' \eqn{+1} in the p-value below refers to. With a single cell type that type is
#' its own partner, so the shift is applied to one copy while the other is held:
#' a common offset would move both sides together and leave the statistic where
#' it started.
#'
#' One pass of \eqn{B} draws evaluates the null at *every* bandwidth on the
#' grid, so the draws are coupled across \eqn{\sigma}. That coupling is what
#' makes the accompanying p-value valid: the null distribution of
#' \eqn{\max_\sigma z(\sigma)} is available from the same draws, so comparing
#' the observed \eqn{\max_\sigma z(\sigma)} to it accounts for having scanned
#' the grid (single-step Westfall--Young max-T). The selection is replicated
#' inside the null, so the reported p-value is not circular. Running
#' [runSkrCCAPermu()] once per bandwidth would *not* give this: it redraws
#' permutations on every call and its default bin grid is itself
#' \eqn{\sigma}-dependent, so the per-bandwidth nulls are neither coupled nor
#' comparable.
#'
#' `perSigma$pAdjusted` applies the same max-T reference to each bandwidth
#' individually, giving the scales at which coordination is detectable after
#' adjusting for the scan. `plateau` collects those with `pAdjusted <= alpha`.
#' It answers "where is there detectable signal", not "which bandwidths are
#' indistinguishable from the best" --- \eqn{z(\sigma)} is typically flat-topped
#' and \eqn{\hat\sigma} is weakly identified, so treat it as a representative
#' scale within a band and check that conclusions hold across the band.
#'
#' ## Scope, and an honest limitation
#'
#' The canonical directions are held at their observed values inside each draw;
#' they are not re-optimized. That is what is wanted for *selection* --- it
#' measures the \eqn{\sigma}-shape of the noise floor with the signal held
#' fixed --- and it is what makes one \eqn{O(B)} pass enough. It also means the
#' p-value inherits the mild anti-conservativeness of a fixed-direction null,
#' most visibly at small \eqn{\sigma}. For a re-optimizing test at a chosen
#' bandwidth use [runSkrCCAPermu()], or [runSkrCCAPermu_FairSigma()] for a
#' re-optimizing max-over-sigma test of the raw normalized correlation.
#'
#' That small-\eqn{\sigma} inflation is why `minSigma` defaults to `"spacing"`.
#' Below about one cell spacing a Gaussian kernel has almost no mass off the
#' diagonal, the statistic rests on a handful of pairs, and the fixed-direction
#' null understates the floor --- so \eqn{z} keeps climbing as \eqn{\sigma}
#' shrinks and the argmax rails at whatever the smallest candidate happens to
#' be. Flooring the scan at the median nearest-partner distance removes that
#' regime. If the selected bandwidth still lands at an end of the grid the
#' function warns: the scan found the edge, not an optimum.
#'
#' The kernel used here is the plain Gaussian
#' \eqn{\exp(-\tfrac{1}{2}(d/\sigma)^2)} on the same *Euclidean* distance scale
#' `computeKernelMatrix()` used. Euclidean is a requirement, not a default: the
#' kernel is rebuilt from coordinates at every candidate bandwidth and under
#' every shifted configuration, and a `"Morphology-Aware"` metric is not a
#' function of the coordinates alone --- its geodesic filter rests on a k-NN
#' graph of the unshifted tissue, which a shift invalidates. Such objects are
#' refused rather than silently rescored with a Euclidean kernel. Three shaping
#' steps that the kernel pipeline applies are
#' deliberately *not* replicated: the `upperQuantile` cap, the `lowerLimit`
#' truncation, and `computeDistance(truncateLowDist = TRUE)`'s floor on the
#' bottom ~0.1% of distances. The cap is a data-dependent quantile, so it would
#' be recomputed on every shifted configuration and the statistic would stop
#' meaning the same thing across draws; the other two act only where the
#' Gaussian is flat (at distances far below \eqn{\sigma}, or at kernel values
#' near zero) and so move \eqn{T} negligibly. The selected bandwidth is
#' therefore a bandwidth for the same kernel family on the same distance scale,
#' and can be passed straight back to [computeKernelMatrix()].
#'
#' Because the directions are fixed, the statistic does not depend on which
#' canonical criterion produced them: `sumcov` and `sumcor` weights are scored
#' the same way, and no criterion is re-run inside the null.
#'
#' ## Cell types
#'
#' With one cell type this is the within-type (self) problem: `a` and `b` are
#' the same score vector, the kernel is the within-type kernel with a zero
#' diagonal (a cell is not its own neighbour, the \eqn{w_{ii} = 0} convention),
#' and the null shifts one copy of the coordinates against the other. With two
#' or more, every pair is scanned and the max-T reference is taken over the
#' whole (bandwidth x pair) grid, so `pAdjusted` is adjusted for both.
#'
#' @param object A `CoProSingle` object that has been through [runSkrCCA()] and
#'   [computeGeneAndCellScores()], so per-bandwidth cell scores exist.
#' @param sigmaValues Candidate bandwidths. Defaults to `object@@sigmaValues`.
#'   Every value must have stored cell scores.
#' @param ccIndex Which canonical component to select on. Default 1.
#' @param nPermu Number of toroidal draws \eqn{B}. Default 199, giving a
#'   Monte-Carlo floor of 0.005 under the Phipson--Smyth correction. Must be at
#'   least 2, since a single draw has no spread to studentize by.
#' @param alternative `"greater"` (default, matching the sign the optimizer
#'   fixes) or `"two.sided"`.
#' @param minSigma Floor for the scanned grid. `"spacing"` (default) drops
#'   candidates below the median nearest-partner distance, measured on the full
#'   coordinates. A single positive number sets the floor explicitly; `NULL`
#'   scans the whole grid. See the note on small bandwidths below.
#' @param maxCells Cap on the number of cells sampled per cell type. The cost is
#'   quadratic in this, and the same sample is used for the observed statistic
#'   and for every draw, so the studentization stays internally consistent.
#'   `NULL` uses every cell. Default 2000.
#' @param domain Optional `list(lower = , upper = )` numeric vectors giving the
#'   wrap box per coordinate axis, on the object's distance scale. Defaults to
#'   the range of the sampled coordinates.
#' @param alpha Level used to report `plateau`. Default 0.05.
#' @param blockSize Rows per block when accumulating the bilinear statistic.
#'   Caps peak memory; does not change the result. Default 1024.
#' @param seed Optional integer RNG seed.
#' @param verbose Whether to report progress.
#'
#' @return An object of class `CoProSigmaSelection`: a list with
#'   \itemize{
#'     \item `sigmaValueChoice` --- the selected bandwidth.
#'     \item `pValue` --- the max-T permutation p-value at that bandwidth.
#'     \item `zMax` --- the observed \eqn{\max z}.
#'     \item `plateau` --- bandwidths with `pAdjusted <= alpha`.
#'     \item `perSigma` --- a `data.frame` with `sigma`, `cellType1`,
#'       `cellType2`, `statistic` (\eqn{T}), `nullMean` (\eqn{m}), `nullSD`
#'       (\eqn{s}), `z`, and `pAdjusted`.
#'     \item `nullMax` --- the length-`nPermu` null distribution of the maximum.
#'     \item `cells` --- the sampled cell indices per cell type.
#'     \item `spacing` --- the median nearest-partner distance per pair.
#'     \item `settings` --- `nPermu`, `alternative`, `ccIndex`, `alpha`,
#'       `maxCells`, `null`, and `minSigma` (the floor actually applied).
#'   }
#'
#' @references
#' Meinshausen N, Maathuis MH, Buhlmann P (2011). Asymptotic optimality of the
#' Westfall--Young permutation procedure for multiple testing under dependence.
#' *Annals of Statistics* 39(6):3369--3391.
#'
#' Phipson B, Smyth GK (2010). Permutation p-values should never be zero.
#' *Statistical Applications in Genetics and Molecular Biology* 9(1):39.
#'
#' @seealso [detectSigmaRange()] to build the candidate grid,
#'   [runSkrCCAPermu()] for a re-optimizing test at a fixed bandwidth,
#'   [runSkrCCAPermu_FairSigma()] for a re-optimizing max-over-sigma test.
#' @family spatial-pipeline
#'
#' @examples
#' \dontrun{
#' obj <- computeKernelMatrix(obj, sigmaValues = detectSigmaRange(obj)$sigmaValues)
#' obj <- runSkrCCA(obj, scalePCs = TRUE, nCC = 2)
#' obj <- computeNormalizedCorrelation(obj)
#' obj <- computeGeneAndCellScores(obj)
#'
#' sel <- selectSigmaByPermutation(obj, nPermu = 199, seed = 1)
#' sel
#' sel$sigmaValueChoice
#' }
#'
#' @importFrom stats sd runif
#' @export
selectSigmaByPermutation <- function(object,
                                     sigmaValues = NULL,
                                     ccIndex = 1L,
                                     nPermu = 199L,
                                     alternative = c("greater", "two.sided"),
                                     minSigma = "spacing",
                                     maxCells = 2000L,
                                     domain = NULL,
                                     alpha = 0.05,
                                     blockSize = 1024L,
                                     seed = NULL,
                                     verbose = TRUE) {
  alternative <- match.arg(alternative)

  if (!is(object, "CoPro")) {
    stop("Input must be a CoPro object")
  }
  if (inherits(object, "CoProMulti")) {
    stop("selectSigmaByPermutation() supports CoProSingle objects. A toroidal ",
         "shift needs one wrap domain, and a CoProMulti object carries one per ",
         "slide. Select the bandwidth on a representative slide, or pass a ",
         "predeclared bandwidth to the multi-slide routes.")
  }
  if (length(object@cellScores) == 0L) {
    stop("cellScores slot is empty. Run computeGeneAndCellScores() first.")
  }

  ccIndex <- .checkPositiveScalarInt(ccIndex, "ccIndex")
  nPermu <- .checkPositiveScalarInt(nPermu, "nPermu")
  blockSize <- .checkPositiveScalarInt(blockSize, "blockSize")
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single number in (0, 1).")
  }
  if (!is.null(maxCells)) {
    maxCells <- .checkPositiveScalarInt(maxCells, "maxCells")
  }
  if (nPermu < 2L) {
    ## One draw has no spread, so every column studentizes to NA and the
    ## function would return the first candidate with pValue = NA and
    ## zMax = -Inf -- a result shaped like an answer that is not one.
    stop("nPermu must be at least 2: the null spread is estimated from the ",
         "draws themselves, and a single draw has none. Use nPermu >= 199 for ",
         "a usable p-value.")
  }
  if (nPermu < 19L) {
    warning("nPermu = ", nPermu, " gives a Monte-Carlo floor of ",
            format(1 / (nPermu + 1), digits = 3),
            "; the null SD is also noisy. Consider nPermu >= 199.",
            call. = FALSE)
  }
  if (!is.null(seed)) {
    rng_state <- .captureRNGState()
    on.exit(.restoreRNGState(rng_state), add = TRUE)
    set.seed(seed)
  }

  sigmaValues <- .resolveSelectionSigmas(object, sigmaValues)
  cts <- .resolveSelectionCellTypes(object)

  ## Coordinates on the same scale sigma lives on: the recorded per-axis scales
  ## that computeDistance() applied, times the raw -> normalized factor it
  ## recorded. Without this the grid and the kernel would be in different units.
  coords <- .selectionCoords(object, cts)
  pairs <- .selectionPairs(cts)

  ## Floor the grid at the typical nearest-partner spacing, on the FULL
  ## coordinates -- subsampling would inflate the spacing and weaken the floor.
  spacing <- .selectionSpacing(coords, pairs)
  keep <- .applySigmaFloor(sigmaValues, minSigma, spacing, verbose)
  sigmaValues <- keep$sigmaValues

  cells <- .sampleSelectionCells(coords, maxCells)
  coords <- stats::setNames(
    lapply(cts, function(ct) coords[[ct]][cells[[ct]], , drop = FALSE]), cts
  )

  scores <- .selectionScores(object, cts, sigmaValues, ccIndex, cells)
  domain <- .resolveSelectionDomain(domain, coords)

  nSigma <- length(sigmaValues)
  nPair <- length(pairs)
  if (verbose) {
    cat("Studentized max-T bandwidth selection\n")
    cat(sprintf("  bandwidths : %d (%s)\n", nSigma,
                paste(format(sigmaValues), collapse = ", ")))
    cat(sprintf("  pairs      : %d%s\n", nPair,
                if (length(cts) == 1L) "  [within-type]" else ""))
    cat(sprintf("  cells      : %s\n",
                paste(sprintf("%s=%d", cts, vapply(cells, length, integer(1))),
                      collapse = ", ")))
    cat(sprintf("  draws      : %d (toroidal shift)\n", nPermu))
  }

  ## Observed statistic on the (bandwidth x pair) grid.
  statObs <- numeric(nSigma * nPair)
  for (pp in seq_len(nPair)) {
    pr <- pairs[[pp]]
    statObs[.selectionCols(pp, nSigma)] <- .bilinearOverSigma(
      coordsA = coords[[pr$ct1]], coordsB = coords[[pr$ct2]],
      a = scores[[pr$ct1]], b = scores[[pr$ct2]],
      sigmaValues = sigmaValues, selfPair = pr$self, blockSize = blockSize
    )
  }

  ## One O(B) pass: every draw is evaluated at every bandwidth, so the null is
  ## coupled across the grid and max_sigma z has a null distribution. The offset
  ## is drawn once per cell type per draw and reused across every pair that type
  ## appears in, so a row of statNull comes from one configuration -- which is
  ## what makes its maximum a draw from the null of the scan statistic.
  statNull <- matrix(0, nrow = nPermu, ncol = nSigma * nPair)
  report_every <- max(1L, nPermu %/% 5L)
  for (bb in seq_len(nPermu)) {
    shifted <- .drawSelectionShift(coords, domain)
    for (pp in seq_len(nPair)) {
      pr <- pairs[[pp]]
      statNull[bb, .selectionCols(pp, nSigma)] <- .bilinearOverSigma(
        ## Self pair: the type is on both sides, so one copy is held at the
        ## observed positions and only the other moves.
        coordsA = if (pr$self) coords[[pr$ct1]] else shifted[[pr$ct1]],
        coordsB = shifted[[pr$ct2]],
        a = scores[[pr$ct1]], b = scores[[pr$ct2]],
        sigmaValues = sigmaValues, selfPair = pr$self, blockSize = blockSize
      )
    }
    if (verbose && (bb %% report_every == 0L || bb == nPermu)) {
      cat(sprintf("  draw %d/%d\n", bb, nPermu))
    }
  }

  st <- .studentizeSelectMaxT(statObs, statNull, alternative)

  perSigma <- data.frame(
    sigma = rep(sigmaValues, times = nPair),
    cellType1 = rep(vapply(pairs, `[[`, character(1), "ct1"), each = nSigma),
    cellType2 = rep(vapply(pairs, `[[`, character(1), "ct2"), each = nSigma),
    statistic = statObs,
    nullMean = st$nullMean,
    nullSD = st$nullSD,
    z = st$z,
    pAdjusted = st$pAdjusted,
    stringsAsFactors = FALSE
  )
  best <- which.max(st$zSelect)
  plateau <- sort(unique(perSigma$sigma[perSigma$pAdjusted <= alpha]))

  out <- structure(
    list(
      sigmaValueChoice = perSigma$sigma[best],
      pValue = perSigma$pAdjusted[best],
      zMax = st$zMax,
      plateau = plateau,
      perSigma = perSigma,
      nullMax = st$nullMax,
      cells = cells,
      spacing = stats::setNames(
        spacing, vapply(pairs, function(pr) paste(pr$ct1, pr$ct2, sep = "-"),
                        character(1))
      ),
      settings = list(nPermu = nPermu, alternative = alternative,
                      ccIndex = ccIndex, alpha = alpha, maxCells = maxCells,
                      minSigma = keep$floor, null = "toroidal")
    ),
    class = "CoProSigmaSelection"
  )

  if (verbose) {
    cat("\n")
    print(out)
  }
  ## Not gated on `verbose`: a scan that rails at an end of the grid has not
  ## located an optimum, it has run out of grid. That is a property of the
  ## result, not progress chatter, so it must survive verbose = FALSE.
  if (nSigma > 1L && out$sigmaValueChoice %in% range(sigmaValues)) {
    warning("Selected bandwidth ", format(out$sigmaValueChoice),
            " sits at the ",
            if (out$sigmaValueChoice == min(sigmaValues)) "low" else "high",
            " end of the scanned grid, so the most detectable scale may lie ",
            "outside it. Widen sigmaValues (see detectSigmaRange()) and re-run ",
            "before reporting this bandwidth.", call. = FALSE)
  }
  out
}


#' Print a permutation bandwidth selection
#' @param x A `CoProSigmaSelection` from [selectSigmaByPermutation()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @method print CoProSigmaSelection
#' @keywords internal
#' @export
print.CoProSigmaSelection <- function(x, ...) {
  cat("CoPro bandwidth selection (studentized max-T permutation)\n")
  cat(sprintf("  selected sigma : %.4g\n", x$sigmaValueChoice))
  cat(sprintf("  max-T p-value  : %.4g  (B = %d, floor %.4g, %s)\n",
              x$pValue, x$settings$nPermu, 1 / (x$settings$nPermu + 1),
              x$settings$alternative))
  cat(sprintf("  z at optimum   : %.3f\n", x$zMax))
  if (length(x$plateau) == 0L) {
    cat(sprintf("  plateau        : none detectable at alpha = %g\n",
                x$settings$alpha))
  } else {
    cat(sprintf("  plateau        : %s  (pAdjusted <= %g)\n",
                paste(format(x$plateau), collapse = ", "), x$settings$alpha))
  }
  print(x$perSigma)
  invisible(x)
}


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

#' @noRd
.checkPositiveScalarInt <- function(value, what) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 1 || value != as.integer(value)) {
    stop(what, " must be a single positive integer.")
  }
  as.integer(value)
}

#' Bandwidths to scan, checked against what the object can actually supply.
#' @noRd
.resolveSelectionSigmas <- function(object, sigmaValues) {
  if (is.null(sigmaValues)) sigmaValues <- object@sigmaValues
  if (length(sigmaValues) == 0L) {
    stop("No sigmaValues supplied and none stored on the object.")
  }
  if (!is.numeric(sigmaValues) || any(!is.finite(sigmaValues)) ||
      any(sigmaValues <= 0)) {
    stop("sigmaValues must be positive finite numbers.")
  }
  missing <- setdiff(sigmaValues, object@sigmaValues)
  if (length(missing) > 0L) {
    stop("No cell scores stored for sigma ",
         paste(format(missing), collapse = ", "),
         ". Selection uses the directions CoPro fitted at each bandwidth, so ",
         "every candidate must have been through runSkrCCA() and ",
         "computeGeneAndCellScores().")
  }
  sort(unique(as.numeric(sigmaValues)))
}

#' @noRd
.resolveSelectionCellTypes <- function(object) {
  cts <- object@cellTypesOfInterest
  if (length(cts) == 0L) cts <- unique(as.character(object@cellTypesSub))
  if (length(cts) == 0L) stop("No cell types found on the object.")
  as.character(cts)
}

#' Coordinates per cell type, on the distance scale sigma is expressed in.
#'
#' Mirrors what computeDistance() did: scale each axis by the recorded
#' per-axis factor, then apply the raw -> normalized factor it recorded in
#' `@distanceScaleFactor`. A bandwidth is a distance, so getting this wrong
#' would silently scan the wrong physical scales.
#'
#' Only a Euclidean geometry can be rebuilt this way. A "Morphology-Aware"
#' object records distances that are *not* a function of the coordinates alone
#' -- the geodesic filter is fitted on a k-NN graph of the tissue as observed --
#' so there is nothing to reconstruct, and a shift would invalidate the graph
#' the filter was fitted on even if there were. Rebuilding plain Euclidean
#' coordinates from such an object would quietly score a different metric than
#' the one the weights were fitted under, so it is refused instead.
#' @noRd
.selectionCoords <- function(object, cts) {
  loc <- object@locationDataSub
  geom <- .getDistanceGeometry(object)
  distType <- if (is.null(geom$distType)) {
    "Euclidean2D"
  } else {
    as.character(geom$distType)
  }
  if (!distType %in% c("Euclidean2D", "Euclidean3D")) {
    stop("selectSigmaByPermutation() needs a Euclidean geometry, but this ",
         "object records distType = '", distType, "'. The selector rebuilds ",
         "the kernel from coordinates at every candidate bandwidth and under ",
         "every shifted configuration, and a morphology-aware distance is not ",
         "a function of the coordinates alone -- its geodesic filter is fitted ",
         "on a k-NN graph of the unshifted tissue. Rebuild distances with ",
         "distType = 'Euclidean2D' or 'Euclidean3D' to select a bandwidth this ",
         "way, or choose one with runSkrCCAPermu_FairSigma(), which reuses the ",
         "stored kernels instead of rebuilding them.")
  }
  axes <- .geometryAxisScales(geom)
  use3d <- identical(distType, "Euclidean3D")

  needed <- if (use3d) c("x", "y", "z") else c("x", "y")
  if (!all(needed %in% colnames(loc))) {
    stop("locationDataSub must have columns ",
         paste(needed, collapse = ", "), ".")
  }

  mat <- cbind(as.numeric(loc$x) * axes[["x"]],
               as.numeric(loc$y) * axes[["y"]])
  if (use3d) mat <- cbind(mat, as.numeric(loc$z) * axes[["z"]])

  factorRaw <- object@distanceScaleFactor
  if (length(factorRaw) != 1L || !is.finite(factorRaw) || factorRaw <= 0) {
    factorRaw <- 1
  }
  mat <- mat * factorRaw

  labels <- as.character(object@cellTypesSub)
  stats::setNames(lapply(cts, function(ct) {
    idx <- which(labels == ct)
    if (length(idx) < 2L) {
      stop("Cell type '", ct, "' has fewer than 2 cells in the subset.")
    }
    mat[idx, , drop = FALSE]
  }), cts)
}

#' Median nearest-partner spacing per pair, on the full (unsampled) coordinates.
#'
#' Subsampling would inflate this -- half the cells in the same area sit roughly
#' sqrt(2) further apart -- so it has to be measured before `maxCells` bites.
#' @noRd
.selectionSpacing <- function(coords, pairs) {
  vapply(pairs, function(pr) {
    .blockNearestSpacing(coords[[pr$ct1]], coords[[pr$ct2]], within = pr$self)
  }, numeric(1))
}

#' Drop candidate bandwidths below the nearest-partner spacing.
#'
#' Below roughly one cell spacing a Gaussian kernel has almost no mass off the
#' diagonal, so the statistic is carried by a handful of pairs. The
#' fixed-direction null is known to be inflated there, which pulls the argmax
#' down toward the smallest bandwidth on offer. Restricting the scan to
#' `sigma >= spacing` is the same guard the method's Methods description uses.
#' @noRd
.applySigmaFloor <- function(sigmaValues, minSigma, spacing, verbose) {
  floorValue <- if (is.null(minSigma)) {
    NA_real_
  } else if (identical(minSigma, "spacing")) {
    if (all(is.na(spacing))) NA_real_ else max(spacing, na.rm = TRUE)
  } else if (is.numeric(minSigma) && length(minSigma) == 1L &&
             is.finite(minSigma) && minSigma > 0) {
    minSigma
  } else {
    stop("minSigma must be \"spacing\", a single positive number, or NULL.")
  }

  if (!is.finite(floorValue)) {
    return(list(sigmaValues = sigmaValues, floor = NA_real_))
  }
  keep <- sigmaValues >= floorValue
  if (!any(keep)) {
    stop("Every candidate bandwidth is below the nearest-partner spacing (",
         format(floorValue, digits = 3), "). A kernel that narrow has almost ",
         "no mass off the diagonal. Widen sigmaValues (see detectSigmaRange()), ",
         "or pass minSigma = NULL to scan anyway.")
  }
  if (verbose && any(!keep)) {
    cat(sprintf(
      "  dropped %d bandwidth(s) below the nearest-partner spacing %.4g: %s\n",
      sum(!keep), floorValue,
      paste(format(sigmaValues[!keep]), collapse = ", ")
    ))
  }
  list(sigmaValues = sigmaValues[keep], floor = floorValue)
}

#' @noRd
.sampleSelectionCells <- function(coords, maxCells) {
  stats::setNames(lapply(coords, function(m) {
    n <- nrow(m)
    if (is.null(maxCells) || n <= maxCells) seq_len(n) else sort(sample.int(n, maxCells))
  }), names(coords))
}

#' Per-bandwidth cell scores, centered.
#'
#' Centering is what lets the plain kernel stand in for the double-centered
#' one: a' K b = a' K_c b exactly once both score vectors have mean zero.
#' @noRd
.selectionScores <- function(object, cts, sigmaValues, ccIndex, cells) {
  stats::setNames(lapply(cts, function(ct) {
    lapply(sigmaValues, function(s) {
      key <- .createCellScoresName(s, ct, slide = NULL)
      mat <- object@cellScores[[key]]
      if (is.null(mat)) {
        stop("Cell scores not found for sigma = ", s, ", cellType = '", ct, "'.")
      }
      if (ccIndex > ncol(mat)) {
        stop("ccIndex = ", ccIndex, " but only ", ncol(mat),
             " component(s) were computed for '", ct, "'.")
      }
      v <- as.numeric(mat[cells[[ct]], ccIndex])
      v - mean(v)
    })
  }), cts)
}

#' @noRd
.selectionPairs <- function(cts) {
  if (length(cts) == 1L) {
    return(list(list(ct1 = cts[1], ct2 = cts[1], self = TRUE)))
  }
  combos <- utils::combn(cts, 2)
  lapply(seq_len(ncol(combos)), function(k) {
    list(ct1 = combos[1, k], ct2 = combos[2, k], self = FALSE)
  })
}

#' @noRd
.selectionCols <- function(pairIndex, nSigma) {
  (pairIndex - 1L) * nSigma + seq_len(nSigma)
}

#' @noRd
.resolveSelectionDomain <- function(domain, coords) {
  all_coords <- do.call(rbind, coords)
  if (is.null(domain)) {
    return(list(lower = apply(all_coords, 2, min),
                upper = apply(all_coords, 2, max)))
  }
  if (!is.list(domain) || is.null(domain$lower) || is.null(domain$upper)) {
    stop("domain must be a list with `lower` and `upper` numeric vectors.")
  }
  d <- ncol(all_coords)
  if (length(domain$lower) != d || length(domain$upper) != d) {
    stop("domain$lower and domain$upper must each have ", d, " element(s).")
  }
  if (any(domain$upper <= domain$lower)) {
    stop("domain$upper must exceed domain$lower on every axis.")
  }
  list(lower = as.numeric(domain$lower), upper = as.numeric(domain$upper))
}

#' Rigid wrap-around shift, one independent offset per axis.
#'
#' Rigid means every cell moves by the same vector, so the shifted copy keeps
#' its own spatial autocorrelation intact; only its alignment with the other
#' side is destroyed. That is the property a label shuffle does not have.
#' @noRd
.toroidalShift <- function(coords, domain) {
  extent <- domain$upper - domain$lower
  for (j in seq_len(ncol(coords))) {
    ## A constant axis has zero extent -- planar data recorded with a z column,
    ## say -- and `x %% 0` is NaN, which would turn the whole axis, and then
    ## every statistic computed from it, into NaN. There is nothing to wrap
    ## along an axis of zero width, so leave it where it is.
    if (!is.finite(extent[j]) || extent[j] <= 0) next
    shift <- stats::runif(1, 0, extent[j])
    coords[, j] <- (coords[, j] - domain$lower[j] + shift) %% extent[j] +
      domain$lower[j]
  }
  coords
}

#' One draw of the null: a single wrap offset per cell type.
#'
#' Returns the whole coordinate list, shifted. Drawing per cell type rather than
#' per pair is what makes a draw a *configuration*: with three or more types a
#' type appears in several pairs, and an independent offset in each would put
#' the same cells in two places at once, so the row of statistics it produces
#' would not be jointly realizable and the row maximum would not be a draw from
#' the null of the scan.
#'
#' The first cell type is held at its observed positions. That costs no
#' randomization -- the offset between any two types is still uniform over the
#' box -- and leaves one fewer point cloud cut by the wrap seam. With a single
#' cell type there is no partner to be relative to, so that type is the one that
#' moves; the caller pairs the shifted copy against the held one.
#' @noRd
.drawSelectionShift <- function(coords, domain) {
  if (length(coords) == 1L) {
    coords[[1L]] <- .toroidalShift(coords[[1L]], domain)
    return(coords)
  }
  for (k in seq_along(coords)[-1L]) {
    coords[[k]] <- .toroidalShift(coords[[k]], domain)
  }
  coords
}

#' T(sigma) = a' K(sigma) b for every sigma, accumulated in row blocks.
#'
#' `a` and `b` are lists parallel to `sigmaValues`: the canonical directions
#' CoPro fitted at each bandwidth, so the statistic is evaluated at the scores
#' that bandwidth actually produced.
#'
#' Squared distances for a block are formed once and reused across the whole
#' bandwidth grid, so the grid costs one extra exp() per entry rather than a
#' fresh distance computation. Blocking caps peak memory at
#' `blockSize x nrow(coordsB)` instead of the full product.
#' @noRd
.bilinearOverSigma <- function(coordsA, coordsB, a, b, sigmaValues,
                               selfPair, blockSize) {
  nA <- nrow(coordsA)
  out <- numeric(length(sigmaValues))
  sqB <- rowSums(coordsB^2)

  starts <- seq(1L, nA, by = blockSize)
  for (start in starts) {
    idx <- start:min(start + blockSize - 1L, nA)
    blockA <- coordsA[idx, , drop = FALSE]
    d2 <- outer(rowSums(blockA^2), sqB, "+") - 2 * tcrossprod(blockA, coordsB)
    d2[d2 < 0] <- 0
    for (g in seq_along(sigmaValues)) {
      K <- exp(-0.5 * d2 / sigmaValues[g]^2)
      ## w_ii = 0: a cell is not its own neighbour. Held under the shift too,
      ## so the observed and null statistics sum over the same set of pairs.
      if (selfPair) K[cbind(seq_along(idx), idx)] <- 0
      out[g] <- out[g] + sum(a[[g]][idx] * as.numeric(K %*% b[[g]]))
    }
  }
  out
}

#' Center and scale by the per-column permutation moments, then take the max-T
#' reference.
#'
#' `statNull` is `nPermu x (nSigma * nPair)`; because one draw fills a whole
#' row, the row maxima are draws from the null of the scan statistic.
#'
#' Both moments come from the null, and both are subtracted from the observed
#' statistic and from the null itself, so the comparison stays like-for-like.
#' Centering matters most within a cell type, where the `w_ii = 0` convention
#' gives the null a negative, sigma-dependent mean: scaling alone would leave
#' that tilt in the scan and bias the argmax toward wide bandwidths.
#' @noRd
.studentizeSelectMaxT <- function(statObs, statNull, alternative) {
  nPermu <- nrow(statNull)
  nullMean <- colMeans(statNull)
  nullSD <- apply(statNull, 2, stats::sd)
  degenerate <- !is.finite(nullSD) | nullSD <= 0
  if (all(degenerate)) {
    stop("Every bandwidth-pair combination had a null spread of zero, so ",
         "nothing can be studentized and no bandwidth can be compared to ",
         "another. Check that the scanned bandwidths put mass off the kernel ",
         "diagonal (see detectSigmaRange()), that the coordinates are not all ",
         "identical, and that the cell scores are not constant.")
  }
  if (any(degenerate)) {
    warning(sum(degenerate), " of ", length(nullSD), " bandwidth-pair ",
            "combinations had a null spread of zero and were dropped from the ",
            "scan. This usually means the kernel is empty at that bandwidth.",
            call. = FALSE)
  }
  z <- (statObs - nullMean) / nullSD
  zNull <- sweep(sweep(statNull, 2, nullMean, "-"), 2, nullSD, "/")

  if (alternative == "two.sided") {
    zSelect <- abs(z)
    zNullScan <- abs(zNull)
  } else {
    zSelect <- z
    zNullScan <- zNull
  }
  zSelect[degenerate] <- -Inf
  zNullScan[, degenerate] <- -Inf
  nullMax <- apply(zNullScan, 1, max)
  ## Phipson & Smyth (2010): the observed configuration is itself an admissible
  ## draw, so both counts get a +1 and the smallest resolvable p is 1/(B + 1).
  pAdjusted <- vapply(zSelect, function(zi) {
    (1 + sum(nullMax >= zi)) / (1 + nPermu)
  }, numeric(1))

  list(nullMean = nullMean,
       nullSD = ifelse(degenerate, NA_real_, nullSD),
       z = ifelse(degenerate, NA_real_, z),
       zSelect = zSelect,
       zMax = max(zSelect[!degenerate]),
       nullMax = nullMax,
       pAdjusted = ifelse(degenerate, NA_real_, pAdjusted))
}
