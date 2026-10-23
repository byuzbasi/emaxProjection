#' Project a Gaussian Emax region onto the integrated response
#'
#' Uses an adaptive search over the entire closed shape interval, exact
#' fixed-shape quadratic slices, necessary box bounds, and feasible witnesses.
#' The refined enclosure is intersected with the original finite-budget
#' enclosure, using at most 256 reference splits. Tolerance and split-budget
#' defaults match the completed research implementation.
#'
#' Bounds hold by the exact-arithmetic argument and are evaluated in ordinary
#' double precision with roundoff cushions. They are not formal machine
#' certificates: `outward_rounding_certified` is always `FALSE`.
#'
#' An empty model intersection has no confidence interval. A finite search
#' without a feasible witness returns `unresolved_feasibility`, not an empty
#' set. Inspect both `status` and `converged`, and report unresolved results.
#' @param region Object from [emax_region()] or [emax_region_from_summary()].
#' @param tolerance Positive absolute endpoint-gap target, in response units.
#' @param max_splits Nonnegative integer adaptive split budget.
#' @return An `emax_projection` list with the outer `interval`, endpoint
#'   brackets, `width_bounds`, `endpoint_gaps`, `status`, `converged`,
#'   `termination`, `witnesses`, and a diagnostic `trace`. Empty-set width is
#'   zero and interval endpoints are missing, matching the research convention.
#'   An empty set counts as noncoverage in coverage studies.
#' @examples
#' region <- emax_region_from_summary(c(0, 0.5, 1), rep(20, 3),
#'   c(0.2, 0.6, 0.8), pooled_variance = 1)
#' result <- emax_project(region)
#' result
#' @export
emax_project <- function(region, tolerance = 0.001, max_splits = 4096L) {
  .check_region(region)
  .scalar(tolerance, "tolerance", positive = TRUE)
  .scalar(max_splits, "max_splits", 0, .Machine$integer.max)
  if (max_splits != floor(max_splits)) stop("max_splits must be an integer", call. = FALSE)
  result <- cpp_project(region, tolerance, as.integer(max_splits))
  result$calibration <- region$calibration
  result$alpha <- region$alpha
  result$tolerance <- tolerance
  result$max_splits <- as.integer(max_splits)
  class(result) <- "emax_projection"
  result
}

#' Fixed-shape slice of an Emax confidence region
#' @inheritParams emax_project
#' @param q Scalar shape in `[0, 1]`.
#' @return List with `status`, and when feasible the nonnegative `delta`
#'   interval and its integrated-response `interval`.
#' @export
emax_fixed_shape <- function(region, q) {
  .check_region(region); .scalar(q, "q", 0, 1)
  cpp_fixed_slice(region, q)
}

#' Analytic outer envelope of an Emax projection
#'
#' Endpoint and interior-dose contrasts bound the effect and shape without
#' optimization. A nonempty envelope does not establish feasibility of the
#' exact model intersection. Use [emax_project()] for the refined interval.
#' @inheritParams emax_project
#' @return A list with `status`, `interval`, and, when nonempty, the effect
#'   bounds `delta` and shape bounds `q`.
#' @export
emax_outer_interval <- function(region) {
  .check_region(region)
  cpp_outer(region)
}

#' Print a projection and its numerical status
#' @param x An `emax_projection` object.
#' @param digits Number of displayed significant digits.
#' @param ... Reserved for compatibility with generic methods.
#' @return The input object, invisibly.
#' @export
print.emax_projection <- function(x, digits = 6, ...) {
  cat("Integrated Emax response | ", x$calibration, " calibration\n", sep = "")
  cat("Status: ", x$status, "; termination: ", x$termination, "\n", sep = "")
  if (x$status != "empty_set") {
    cat("Outer interval: [", format(x$interval[1], digits = digits), ", ",
      format(x$interval[2], digits = digits), "]\n", sep = "")
    cat("Endpoint gaps: ", paste(format(x$endpoint_gaps, digits = digits), collapse = ", "),
      "; converged: ", x$converged, "\n", sep = "")
  }
  cat("Ordinary double precision; no directed-rounding certificate.\n")
  invisible(x)
}

#' Plot the endpoint enclosure diagnostics
#'
#' Shows outer bounds and feasible inner witnesses through the adaptive search.
#' This is a numerical convergence display, not a dose-response confidence band.
#' @param x An `emax_projection` object.
#' @param ... Graphical parameters passed to `graphics::plot()`.
#' @return The input object, invisibly.
#' @export
plot.emax_projection <- function(x, ...) {
  tr <- x$trace
  if (!nrow(tr)) stop("no diagnostic trace is available", call. = FALSE)
  values <- unlist(tr[c("lower_outer", "lower_inner", "upper_inner", "upper_outer")])
  finite <- values[is.finite(values)]
  if (!length(finite)) stop("no finite trace bounds are available", call. = FALSE)
  graphics::plot(tr$splits, tr$lower_outer, type = "l", ylim = range(finite),
    xlab = "Adaptive splits", ylab = "Integrated response", ...)
  graphics::lines(tr$splits, tr$upper_outer)
  graphics::lines(tr$splits, tr$lower_inner, lty = 2, col = "#167D8D")
  graphics::lines(tr$splits, tr$upper_inner, lty = 2, col = "#167D8D")
  graphics::legend("right", c("Outer enclosure", "Feasible witnesses"),
    lty = c(1, 2), col = c("black", "#167D8D"), bty = "n")
  invisible(x)
}
