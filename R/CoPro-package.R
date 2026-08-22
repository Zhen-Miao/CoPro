#' @keywords internal
#' @aliases CoPro-package
#'
#' @section Naming convention:
#'
#' The exported API deliberately uses two naming styles, split by layer:
#'
#' \itemize{
#'   \item **`camelCase` -- the object pipeline.** Anything whose first argument
#'     is a `CoProSingle` or `CoProMulti` object: constructors
#'     ([newCoProSingle()], [newCoProMulti()]), the S4 generics that advance an
#'     analysis ([computePCA()], [detectSigmaRange()], [computeKernelMatrix()],
#'     [runSkrCCA()], [computeNormalizedCorrelation()]), and the accessors that
#'     read results back out ([getCellScores()], [getNormCorr()],
#'     [getKernelMatrix()]). These take an object and return an object, so they
#'     chain.
#'   \item **`snake_case` -- the engine and utility layer.** Functions that work
#'     on plain matrices, lists, and data frames with no `CoPro` object
#'     involved: the CCA solvers ([optimize_bilinear()],
#'     [optimize_sumcor_pca()], [optimize_genespace_avg_corr()] and their
#'     `_n` / `_multi_slides` variants), the spatial-null builders
#'     ([resample_spatial()],
#'     [generate_toroidal_permutations()], [diagnose_bin_distribution()]), and
#'     standalone helpers ([quantile_normalize()], [transfer_scores()],
#'     [copro_download_data()]). Call these directly when you want the numerical
#'     core without the object wrapper.
#' }
#'
#' A few exports sit outside that rule:
#'
#' \itemize{
#'   \item [runSkrCCAPermu_FairSigma()] and [runSkrCCAPermu_Conditional()] keep
#'     the `camelCase` stem of [runSkrCCAPermu()] with a `_Variant` suffix, so
#'     the three permutation entry points sort together.
#'   \item [calculate_pvalue()], [calculate_pvalue_stepdown()],
#'     [compute_ground_truth_ncorr()], and [fit_score_reference()] take a
#'     `CoPro` object but are `snake_case`. They derive standalone inference or
#'     transfer results rather than advancing and returning the object pipeline.
#' }
#'
#' @section Package options:
#'
#' Behavior that used to be set through `options()` is now reachable as function
#' arguments; the options are still read to supply those arguments' defaults, so
#' existing scripts keep working.
#'
#' \itemize{
#'   \item `CoPro.factorizePermutation` -- default for the `factorize` argument
#'     of [runSkrCCAPermu()] and friends.
#'   \item `CoPro.compactPermutation` -- default for their `compactPermutation`
#'     argument.
#'   \item `CoPro.float32Threads` -- default for the `nThreads` argument of
#'     [computeSparseKernelFloat32()].
#'   \item `CoPro.useRcppFRNN` -- default for the compiled fixed-radius neighbor
#'     engine used internally by the sparse kernel routes.
#' }
#'
#' @useDynLib CoPro, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
