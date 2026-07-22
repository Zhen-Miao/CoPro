test_that("multi-slide inference uses held-out equal-replicate effects", {
  obj <- create_test_copro_multi(
    n_cells_per_slide = 70, n_slides = 3, n_genes = 30,
    n_cell_types = 2, seed = 20260729
  )
  obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
  obj <- computePCA(obj, nPCA = 8, center = TRUE, scale. = TRUE)
  obj <- computeDistance(
    obj, distType = "Euclidean2D", normalizeDistance = TRUE,
    verbose = FALSE
  )
  obj <- computeKernelMatrix(
    obj, sigmaValues = c(0.05, 0.1), verbose = FALSE
  )

  expect_warning(
    result <- runSlideLevelInference(
      obj, cc_index = 1, n_resamples = 199, seed = 99,
      verbose = FALSE
    ),
    "Fewer than five"
  )

  expect_s3_class(result, "CoProSlideInference")
  expect_equal(result$n_replicates, 3L)
  expect_equal(nrow(result$replicate_effects), 3L)
  expect_equal(result$estimate,
               mean(result$replicate_effects$heldout_effect))
  expect_true(all(result$replicate_effects$selected_sigma %in%
                    obj@sigmaValues))
  expect_true(result$p_value >= result$mc_floor && result$p_value <= 1)
  expect_identical(result$null_method, "exact_replicate_sign_flip")
  expect_true(result$selection_adjusted)
})

test_that("cell-level permutation refuses to masquerade as multi-slide inference", {
  obj <- create_test_copro_multi(
    n_cells_per_slide = 30, n_slides = 2, n_genes = 20,
    n_cell_types = 2, seed = 3
  )
  expect_error(
    runSkrCCAPermu(obj, nPermu = 10, verbose = FALSE),
    "runSlideLevelInference"
  )
})
