test_that("closure and integral match analytic endpoints and quadrature", {
  expect_equal(emax_shape(c(0, 0.5, 1), 0, 1), c(0, 0.5, 1))
  expect_equal(emax_shape(c(0, 0.5, 1), 1, 1), c(0, 1, 1))
  expect_equal(emax_shape_area(0), 0.5)
  expect_equal(emax_shape_area(1), 1)
  for (q in c(0.001, 0.01, 0.25, 0.75, 0.95)) {
    numerical <- integrate(function(u) u / (1 - q + q * u), 0, 1,
      rel.tol = 1e-10)$value
    expect_equal(emax_shape_area(q), numerical, tolerance = 1e-9)
  }
  d <- emax_shape_derivatives(c(0, 0.2, 1), 0.75, 1)
  step <- 1e-6
  numerical <- (emax_shape(c(0, 0.2, 1), 0.75 + step, 1) -
    emax_shape(c(0, 0.2, 1), 0.75 - step, 1)) / (2 * step)
  expect_equal(d$first, numerical, tolerance = 1e-8)
})

test_that("raw and summary interfaces use the same observations and F pivot", {
  x <- rep(c(0, 0.5, 1), c(3, 4, 5))
  y <- c(-0.1, 0, 0.1, 0.3, 0.4, 0.5, 0.6, 0.5, 0.6, 0.7, 0.8, 0.9)
  r <- emax_region(x, y)
  expect_equal(r$n, length(y))
  expect_equal(r$counts, c(3, 4, 5))
  expect_equal(r$calibration_rank, 2)
  expect_equal(r$radius_squared, 2 * r$pooled_variance * qf(0.95, 2, 9))
  m <- emax_region(x, y, calibration = "mean")
  expect_gt(m$radius_squared, r$radius_squared)
  s <- emax_region_from_summary(r$x, r$counts, r$means, r$pooled_variance)
  expect_equal(emax_project(r)$interval, emax_project(s)$interval)
  expect_error(emax_region(x, replace(y, 1, NA_real_)), "finite")
  expect_error(emax_region(x, rep(1, length(y))), "positive")
})

test_that("truth is enclosed at flat, linear, and saturation boundaries", {
  for (q in c(0, 0.75, 1)) {
    x <- c(0, 0.5, 1)
    r <- emax_region_from_summary(x, rep(20, 3), emax_response(x, 0, 1, q, 1), 1)
    fit <- emax_project(r)
    truth <- emax_integrated_response(1, q)
    expect_identical(fit$status, "nonempty")
    expect_true(fit$converged)
    expect_lte(fit$interval[1], truth + 1e-10)
    expect_gte(fit$interval[2], truth - 1e-10)
    expect_lte(max(fit$endpoint_gaps), 0.001)
    expect_false(fit$outward_rounding_certified)
    expect_equal(fit$width_bounds[2], diff(fit$interval))
  }
  flat <- emax_project(emax_region_from_summary(c(0, 0.5, 1), rep(20, 3), c(0, 0, 0), 1))
  expect_equal(flat$interval[1], 0)
})

test_that("empty and budget-limited searches retain scientific status", {
  empty <- emax_project(emax_region_from_summary(c(0, 0.5, 1), rep(20, 3), c(2, 1, 0), 0.001))
  expect_identical(empty$status, "empty_set")
  expect_equal(empty$width, 0)
  expect_true(all(is.na(empty$interval)))
  r <- emax_region_from_summary(c(0, 0.5, 1), rep(20, 3), c(0.2, 0.6, 0.8), 1)
  limited <- emax_project(r, tolerance = 1e-12, max_splits = 0)
  expect_identical(limited$termination, "budget_exhausted")
  expect_false(limited$converged)
  expect_error(emax_project(r, max_splits = 0.5), "integer")
  corrupt <- r
  corrupt$radius_squared <- corrupt$radius_squared * 2
  expect_error(emax_project(corrupt), "inconsistent")
  corrupt$counts <- c(20, 20)
  expect_error(emax_project(corrupt), "ordered")
  expect_error(emax_shape(c(0, 2), 0.5, 1), "belong")
})

test_that("baseline and dose-unit changes preserve the target", {
  x <- c(0, 0.25, 0.5, 1)
  means <- emax_response(x, 0, 2, 0.75, 1)
  r <- emax_region_from_summary(x, rep(20, 4), means, 1)
  shift <- emax_region_from_summary(x * 100, rep(20, 4), means + 5, 1, D = 100)
  expect_equal(emax_project(r)$interval, emax_project(shift)$interval, tolerance = 1e-9)
  expect_equal(emax_shape(x, 0.75, 1), emax_shape(x * 100, 0.75, 100))
})
