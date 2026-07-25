# Extracted from test-pca-workflow.R:491

# test -------------------------------------------------------------------------
obj <- create_test_copro_multi(n_cells_per_slide = 120, n_slides = 2,
                                 n_cell_types = 2, seed = 91)
obj <- subsetData(obj, cellTypesOfInterest = c("CellTypeA", "CellTypeB"))
obj <- suppressMessages(computePCA(obj, nPCA = 6))
expect_length(obj@pcaResults, 2L)
expect_setequal(names(obj@pcaResults), c("Slide1", "Slide2"))
for (sID in c("Slide1", "Slide2")) {
    for (ct in c("CellTypeA", "CellTypeB")) {
      entry <- obj@pcaResults[[sID]][[ct]]
      expect_true(CoPro:::.isPCSlice(entry))
      expect_type(entry$rows, "integer")

      # Resolving it must give exactly what the old code materialized: the
      # rows of the global scores for this slide, labelled with cell IDs.
      resolved <- CoPro:::.resolvePCSlice(entry, obj@pcaGlobal[[ct]]$x)
      keep <- getSlideID(obj)[obj@cellTypesSub == ct] == sID
      expected <- obj@pcaGlobal[[ct]]$x[keep, , drop = FALSE]
      expect_identical(resolved, expected)
      expect_identical(
        rownames(resolved),
        rownames(obj@metaDataSub)[obj@cellTypesSub == ct][keep]
      )
    }
  }
for (ct in c("CellTypeA", "CellTypeB")) {
    all_rows <- sort(unlist(lapply(obj@pcaResults, function(s) s[[ct]]$rows)))
    expect_identical(all_rows, seq_len(nrow(obj@pcaGlobal[[ct]]$x)))
  }
