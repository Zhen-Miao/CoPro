#' Replicate-level inference for multi-slide CoPro analyses
#'
#' Estimates an equal-replicate spatial coordination effect without treating
#' cells as biological replicates. For each replicate, canonical weights and
#' the bandwidth are learned from all other replicates and evaluated only on
#' the held-out replicate. The held-out effects are combined with equal weight.
#'
#' @details
#' `replicate_id` may group several slides from the same donor. Every fold then
#' leaves out all slides belonging to one donor. Sigma selection is performed
#' using training slides only, so the held-out effect is not evaluated at a
#' bandwidth selected from that replicate.
#'
#' The p-value uses sign flips of the replicate-level held-out effects and
#' therefore assumes independent replicates and a symmetric null distribution
#' around zero. Up to 15 replicates all sign patterns are enumerated exactly;
#' larger analyses use `n_resamples` Monte-Carlo sign flips with the
#' Phipson--Smyth correction. The confidence interval is a percentile bootstrap
#' over replicates. Both summaries target the equal-replicate mean.
#'
#' @param object A `CoProMulti` object after [computePCA()] and kernel
#'   construction. Exactly two cell types are currently supported.
#' @param cc_index Canonical axis to evaluate.
#' @param sigma_values Candidate sigma values. Defaults to all fitted kernel
#'   bandwidths.
#' @param replicate_id Optional slide-to-replicate mapping. Supply a vector
#'   named by `getSlideList(object)`; an unnamed vector is matched in slide
#'   order. Defaults to one independent replicate per slide.
#' @param alternative One of `"greater"`, `"less"`, or `"two.sided"`.
#' @param n_resamples Number of Monte-Carlo sign flips and bootstrap resamples.
#' @param conf_level Bootstrap confidence level.
#' @param seed Optional random seed.
#' @param verbose Print fold progress.
#'
#' @return A list of class `CoProSlideInference` containing the equal-replicate
#'   effect, confidence interval, replicate sign-flip p-value, Monte-Carlo
#'   floor, and one held-out result per replicate.
#' @family spatial-pipeline
#' @export
runSlideLevelInference <- function(object, cc_index = 1,
                                   sigma_values = NULL,
                                   replicate_id = NULL,
                                   alternative = "greater",
                                   n_resamples = 9999,
                                   conf_level = 0.95,
                                   seed = NULL,
                                   verbose = TRUE) {
  if (!inherits(object, "CoProMulti")) {
    stop("runSlideLevelInference() requires a CoProMulti object.")
  }
  cts <- object@cellTypesOfInterest
  if (length(cts) != 2L) {
    stop("Replicate-level inference currently requires exactly two cell types.")
  }
  if (length(object@pcaResults) == 0L || length(object@pcaGlobal) == 0L) {
    stop("PCA results are missing. Run computePCA() first.")
  }
  if (length(object@kernelMatrices) == 0L) {
    stop("Kernel matrices are missing.")
  }
  if (!is.numeric(cc_index) || length(cc_index) != 1L ||
      cc_index < 1L || cc_index != as.integer(cc_index)) {
    stop("cc_index must be a positive integer.")
  }
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  if (!is.numeric(n_resamples) || length(n_resamples) != 1L ||
      n_resamples < 99L || n_resamples != as.integer(n_resamples)) {
    stop("n_resamples must be an integer >= 99.")
  }
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be in (0, 1).")
  }

  slides <- getSlideList(object)
  if (is.null(replicate_id)) {
    replicate_id <- stats::setNames(slides, slides)
  } else {
    if (length(replicate_id) != length(slides)) {
      stop("replicate_id must have one value per slide.")
    }
    if (is.null(names(replicate_id))) {
      names(replicate_id) <- slides
    }
    if (!all(slides %in% names(replicate_id))) {
      stop("A named replicate_id must contain every slide name.")
    }
    replicate_id <- replicate_id[slides]
  }
  if (anyNA(replicate_id) || any(!nzchar(as.character(replicate_id)))) {
    stop("replicate_id cannot contain missing or empty values.")
  }
  replicate_id <- as.character(replicate_id)
  names(replicate_id) <- slides
  replicate_slides <- split(slides, replicate_id)
  n_replicates <- length(replicate_slides)
  if (n_replicates < 2L) {
    stop("At least two independent replicates are required.")
  }
  if (n_replicates < 5L) {
    warning("Fewer than five independent replicates: p-value resolution and ",
            "bootstrap coverage will be limited.")
  }

  if (is.null(sigma_values)) sigma_values <- object@sigmaValues
  sigma_values <- unique(as.numeric(sigma_values))
  sigma_values <- sigma_values[sigma_values %in% object@sigmaValues]
  if (length(sigma_values) == 0L) {
    stop("No requested sigma_values are available in the object.")
  }

  X_all <- .preparePCMatrices(
    pc_data = object@pcaResults, pca_global = object@pcaGlobal,
    scalePCs = object@scalePCs, slides = slides, cts = cts
  )
  feature_counts <- vapply(X_all[[slides[1L]]][cts], ncol, integer(1))
  if (cc_index > min(feature_counts)) {
    stop("cc_index exceeds the available PCA rank.")
  }
  sdev2_list <- .permutationSdev2(object, cts)

  # Cache every small PC-space operator and kernel normalizer once. Subsequent
  # folds only add small matrices and multiply held-out score vectors.
  Y_by_sigma <- stats::setNames(vector("list", length(sigma_values)),
                                .formatSigmaValue(sigma_values))
  kernel_info <- stats::setNames(vector("list", length(sigma_values)),
                                 .formatSigmaValue(sigma_values))
  for (sigma in sigma_values) {
    skey <- .formatSigmaValue(sigma)
    Y_by_sigma[[skey]] <- stats::setNames(lapply(slides, function(slide) {
      compute_Y_resi(X_all[[slide]], object@kernelMatrices, sigma, cts,
                     slide = slide)
    }), slides)
    kernel_info[[skey]] <- stats::setNames(lapply(slides, function(slide) {
      .get_ncorr_kernel_info(
        object@kernelMatrices, sigma, cts,
        normalizer_cache = attr(object, "kernelNormalizerCache", exact = TRUE),
        slide = slide
      )
    }), slides)
  }

  score_slide <- function(slide, weights, sigma) {
    info <- kernel_info[[.formatSigmaValue(sigma)]][[slide]]
    s1 <- X_all[[slide]][[cts[1L]]] %*%
      weights[[cts[1L]]][, cc_index, drop = FALSE]
    s2 <- X_all[[slide]][[cts[2L]]] %*%
      weights[[cts[2L]]][, cc_index, drop = FALSE]
    denominator <- sqrt(sum(s1^2)) * sqrt(sum(s2^2)) * info$norm_K12
    if (!is.finite(denominator) || denominator <= 1e-12) return(NA_real_)
    .kernelXKY(s1, info$K, s2)[1, 1] / denominator
  }

  fold_rows <- vector("list", n_replicates)
  replicate_names <- names(replicate_slides)
  for (fold in seq_along(replicate_names)) {
    replicate <- replicate_names[fold]
    heldout <- replicate_slides[[replicate]]
    training <- setdiff(slides, heldout)
    if (length(training) == 0L) stop("A fold has no training slides.")

    candidate <- lapply(sigma_values, function(sigma) {
      skey <- .formatSigmaValue(sigma)
      Y12 <- Reduce(
        "+",
        lapply(training, function(slide) {
          Y_by_sigma[[skey]][[slide]][[cts[1L]]][[cts[2L]]]
        })
      )
      Y_train <- stats::setNames(vector("list", 2L), cts)
      Y_train[[cts[1L]]] <- stats::setNames(vector("list", 2L), cts)
      Y_train[[cts[2L]]] <- stats::setNames(vector("list", 2L), cts)
      Y_train[[cts[1L]]][[cts[2L]]] <- Y12
      Y_train[[cts[2L]]][[cts[1L]]] <- t(Y12)
      weights <- solve_two_type_svd(
        Y_train, cts, nCC = cc_index, sdev2_list = sdev2_list
      )
      training_scores <- vapply(
        training, score_slide, numeric(1), weights = weights, sigma = sigma
      )
      heldout_scores <- vapply(
        heldout, score_slide, numeric(1), weights = weights, sigma = sigma
      )
      list(
        sigma = sigma,
        training_score = mean(training_scores, na.rm = TRUE),
        heldout_score = mean(heldout_scores, na.rm = TRUE),
        slide_scores = heldout_scores
      )
    })
    training_stats <- vapply(candidate, `[[`, numeric(1), "training_score")
    if (all(!is.finite(training_stats))) {
      stop("No finite training statistic in fold for replicate '", replicate, "'.")
    }
    chosen <- candidate[[which.max(training_stats)]]
    fold_rows[[fold]] <- data.frame(
      replicate = replicate,
      heldout_slides = paste(heldout, collapse = ","),
      selected_sigma = chosen$sigma,
      training_statistic = chosen$training_score,
      heldout_effect = chosen$heldout_score,
      stringsAsFactors = FALSE
    )
    if (verbose) {
      message(sprintf(
        "  Held out %s: sigma = %g, effect = %.4f",
        replicate, chosen$sigma, chosen$heldout_score
      ))
    }
  }
  folds <- do.call(rbind, fold_rows)
  effects <- folds$heldout_effect
  if (any(!is.finite(effects))) {
    stop("At least one held-out replicate effect is non-finite.")
  }
  estimate <- mean(effects)

  if (!is.null(seed)) {
    rng_state <- .captureRNGState()
    on.exit(.restoreRNGState(rng_state), add = TRUE)
    set.seed(seed)
  }
  stat_transform <- switch(
    alternative,
    greater = function(x) x,
    less = function(x) -x,
    two.sided = function(x) abs(x)
  )
  observed_test <- stat_transform(estimate)
  if (n_replicates <= 15L) {
    patterns <- 0:(2^n_replicates - 1L)
    null_stats <- vapply(patterns, function(pattern) {
      signs <- ifelse(
        bitwAnd(pattern, bitwShiftL(1L, seq_len(n_replicates) - 1L)) != 0L,
        1, -1
      )
      stat_transform(mean(effects * signs))
    }, numeric(1))
    p_value <- mean(null_stats >= observed_test)
    mc_floor <- 1 / length(null_stats)
    null_method <- "exact_replicate_sign_flip"
  } else {
    null_stats <- replicate(
      n_resamples,
      stat_transform(mean(effects * sample(c(-1, 1), n_replicates,
                                           replace = TRUE)))
    )
    p_value <- (1 + sum(null_stats >= observed_test)) / (n_resamples + 1)
    mc_floor <- 1 / (n_resamples + 1)
    null_method <- "monte_carlo_replicate_sign_flip"
  }

  bootstrap_estimates <- replicate(
    n_resamples, mean(sample(effects, n_replicates, replace = TRUE))
  )
  alpha <- (1 - conf_level) / 2
  conf_int <- as.numeric(stats::quantile(
    bootstrap_estimates, probs = c(alpha, 1 - alpha), names = FALSE,
    type = 8
  ))

  result <- list(
    estimate = estimate,
    conf_int = conf_int,
    conf_level = conf_level,
    p_value = p_value,
    mc_floor = mc_floor,
    alternative = alternative,
    null_method = null_method,
    n_replicates = n_replicates,
    replicate_effects = folds,
    sigma_values = sigma_values,
    cc_index = as.integer(cc_index),
    estimand = "equal-replicate mean held-out normalized correlation",
    selection_adjusted = TRUE
  )
  class(result) <- c("CoProSlideInference", "list")
  result
}
