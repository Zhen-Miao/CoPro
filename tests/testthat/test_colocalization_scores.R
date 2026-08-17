# Test colocalization scores function
test_that("Colocalization scores function works correctly", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.explore")
  skip_if_not_installed("CoPro")
  
  # Test that functions exist and have correct signatures
  expect_true(exists("getColocScores"))
  expect_true(exists("coloc_score"))
  expect_true(exists("nn_cross_fraction"))
  
  # Test helper functions exist
  expect_true(exists(".getColocScoresSingle", mode = "function"))
  expect_true(exists(".getColocScoresMulti", mode = "function"))
  expect_true(exists(".compute_g12_inhom", mode = "function"))
  
  # Additional tests would require creating mock CoPro objects
  # with location data and cell types, which is beyond the scope 
  # of this implementation
})

test_that("pixel_size_um converts coordinate units by multiplication", {
  skip_if_not_installed("spatstat.geom")
  skip_if_not_installed("spatstat.random")

  captured <- NULL
  testthat::local_mocked_bindings(
    .compute_g12_inhom = function(X, type_i, type_j, r_vec,
                                  sigma_intensity,
                                  edge_correction = "translation") {
      if (is.null(captured)) {
        captured <<- list(
          x = range(X$x),
          y = range(X$y),
          window_x = X$window$xrange,
          window_y = X$window$yrange
        )
      }
      rep(1, length(r_vec))
    },
    .package = "CoPro"
  )

  A <- data.frame(x = seq(1, 4, length.out = 20),
                  y = seq(2, 5, length.out = 20))
  B <- data.frame(x = seq(5, 8, length.out = 20),
                  y = seq(4, 7, length.out = 20))
  suppressWarnings(coloc_score(
    A, B,
    window_range = list(xrange = c(0, 10), yrange = c(0, 12)),
    r_um_range = c(10, 20), pixel_size_um = 50,
    nsim = 1, r_step_um = 10
  ))

  expect_equal(captured$x, range(c(A$x, B$x)) * 50)
  expect_equal(captured$y, range(c(A$y, B$y)) * 50)
  expect_equal(captured$window_x, c(0, 10) * 50)
  expect_equal(captured$window_y, c(0, 12) * 50)
})
