# The toroidal null is only worth having if it is what it claims to be: on a
# regular lattice, a rigid translation of the torus. That is a checkable
# property -- every pairwise torus distance must survive every draw -- and it
# is the property the wrap period controls.

lattice_locations <- function(n_x, n_y = n_x, spacing = 1) {
  grid <- expand.grid(ix = seq_len(n_x), iy = seq_len(n_y))
  data.frame(
    x = (grid$ix - 1L) * spacing,
    y = (grid$iy - 1L) * spacing,
    cell_ID = paste0("cell_", seq_len(nrow(grid))),
    stringsAsFactors = FALSE
  )
}

torus_distances <- function(coordinate, period) {
  gap <- abs(outer(coordinate, coordinate, "-"))
  pmin(gap, period - gap)
}

test_that(".torusPeriod adds back the gap the observed span omits", {
  period <- CoPro:::.torusPeriod((seq_len(8L) - 1L) * 0.25)
  expect_equal(period, 8 * 0.25)

  # Irregular coordinates have no true period; the median gap is the estimate.
  set.seed(11)
  values <- sort(runif(200))
  expect_equal(
    CoPro:::.torusPeriod(values),
    diff(range(values)) + median(diff(values))
  )

  expect_true(is.na(CoPro:::.torusPeriod(rep(3, 10))))
})

test_that("every draw is a bijection on cells", {
  locations <- lattice_locations(6L)
  perms <- generate_toroidal_permutations(locations, n_permu = 50, seed = 1)

  expect_equal(dim(perms), c(nrow(locations), 50L))
  for (draw in seq_len(ncol(perms))) {
    expect_setequal(perms[, draw], seq_len(nrow(locations)))
  }
})

test_that("draws preserve every pairwise torus distance on a lattice", {
  n_side <- 6L
  spacing <- 0.5
  locations <- lattice_locations(n_side, n_side, spacing)
  period <- n_side * spacing

  reference_x <- torus_distances(locations$x, period)
  reference_y <- torus_distances(locations$y, period)

  perms <- generate_toroidal_permutations(locations, n_permu = 60, seed = 2)
  worst <- 0
  for (draw in seq_len(ncol(perms))) {
    index <- perms[, draw]
    worst <- max(
      worst,
      max(abs(torus_distances(locations$x[index], period) - reference_x)),
      max(abs(torus_distances(locations$y[index], period) - reference_y))
    )
  }
  # A rigid translation distorts nothing; floating point in the wrap is the
  # only source of error.
  expect_lt(worst, 1e-12)
})

test_that("a non-square lattice is wrapped on its own period per axis", {
  locations <- lattice_locations(5L, 8L, spacing = 0.2)
  perms <- generate_toroidal_permutations(locations, n_permu = 40, seed = 3)

  reference_x <- torus_distances(locations$x, 5 * 0.2)
  reference_y <- torus_distances(locations$y, 8 * 0.2)
  for (draw in seq_len(ncol(perms))) {
    index <- perms[, draw]
    expect_lt(max(abs(torus_distances(locations$x[index], 5 * 0.2) -
                        reference_x)), 1e-12)
    expect_lt(max(abs(torus_distances(locations$y[index], 8 * 0.2) -
                        reference_y)), 1e-12)
  }
})

test_that("the reference set is the whole translation group", {
  # Wrapping on the extent instead of the period cannot reach n_x * n_y
  # distinct maps, and restricting shifts away from zero cannot reach the
  # identity. Both show up here.
  n_side <- 5L
  locations <- lattice_locations(n_side)
  perms <- generate_toroidal_permutations(locations, n_permu = 3000, seed = 4)

  distinct <- unique(split(perms, col(perms)))
  expect_equal(length(distinct), n_side^2)

  identity <- seq_len(nrow(locations))
  expect_true(any(vapply(distinct, function(p) all(p == identity), logical(1))))
})

test_that("the induced maps are closed under composition", {
  # A translation group composed with itself stays inside itself. This is what
  # makes the test exact, and it fails as soon as the wrap is off by a column.
  n_side <- 4L
  locations <- lattice_locations(n_side)
  perms <- generate_toroidal_permutations(locations, n_permu = 400, seed = 5)
  distinct <- unique(split(perms, col(perms)))

  keys <- vapply(distinct, paste, character(1), collapse = ",")
  for (a in distinct[seq_len(min(6L, length(distinct)))]) {
    for (b in distinct[seq_len(min(6L, length(distinct)))]) {
      expect_true(paste(a[b], collapse = ",") %in% keys)
    }
  }
})

test_that("shifts are drawn uniformly over the group", {
  n_side <- 4L
  locations <- lattice_locations(n_side)
  perms <- generate_toroidal_permutations(locations, n_permu = 4000, seed = 6)

  counts <- table(apply(perms, 2L, paste, collapse = ","))
  expect_equal(length(counts), n_side^2)
  # 16 equally likely maps over 4000 draws; a chi-square only fails here if the
  # shift distribution is skewed, which is what clipping the shift range does.
  expect_gt(chisq.test(as.numeric(counts))$p.value, 1e-4)
})

test_that("a seed makes the draws reproducible and does not leak", {
  locations <- lattice_locations(4L)
  set.seed(99)
  before <- runif(1)

  set.seed(99)
  invisible(runif(1))
  first <- generate_toroidal_permutations(locations, n_permu = 5, seed = 42)
  after <- runif(1)

  second <- generate_toroidal_permutations(locations, n_permu = 5, seed = 42)
  expect_identical(first, second)

  # The caller's stream continues as if the function had not drawn at all.
  set.seed(99)
  expect_equal(before, runif(1))
  set.seed(99)
  invisible(runif(1))
  expect_equal(after, runif(1))
})

test_that("a degenerate axis warns and returns identity permutations", {
  flat <- data.frame(x = seq_len(10), y = rep(2, 10),
                     cell_ID = paste0("cell_", seq_len(10)))
  expect_warning(perms <- generate_toroidal_permutations(flat, n_permu = 3),
                 "Spatial extent")
  expect_true(all(perms == seq_len(10)))
})
