.scalar <- function(value, name, lower = -Inf, upper = Inf, positive = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < lower || value > upper || (positive && value <= 0))
    stop(name, if (positive) " must be a positive finite scalar" else
      " must be a finite scalar in its stated domain", call. = FALSE)
  invisible(value)
}
.vector <- function(value, name) {
  if (!is.numeric(value) || !length(value) || any(!is.finite(value)))
    stop(name, " must be a nonempty finite numeric vector", call. = FALSE)
  invisible(value)
}
.shape_inputs <- function(x, q, D) {
  .vector(x, "x"); .scalar(q, "q", 0, 1); .scalar(D, "D", positive = TRUE)
  if (any(x < 0 | x > D)) stop("x must belong to [0, D]", call. = FALSE)
}

#' Closed Emax shape and its normalized integrated response
#'
#' The closed model has mean `a + delta * emax_shape(x, q, D)`, with an
#' unrestricted baseline `a`, nonnegative endpoint effect `delta`, and shape
#' `q` in `[0, 1]`. The normalized integrated excess response is
#' `delta * emax_shape_area(q)`. The endpoint `q = 0` is linear; `q = 1`
#' represents immediate saturation at every positive dose, with zero shape
#' at dose zero. For `0 < q < 1`, ordinary Emax parameters satisfy
#' `K = D * (1 - q) / q` and `b = delta / q`.
#'
#' Integration uses the specified linear dose axis. Arbitrary equally spaced
#' codes do not recover an integral on unknown physical doses.
#' @param x Numeric nonempty dose vector in `[0, D]`.
#' @param q Scalar shape in `[0, 1]`.
#' @param D Positive dose-axis upper endpoint, in the same units as `x`.
#' @param a Finite baseline response.
#' @param delta Nonnegative endpoint excess response, in response units.
#' @return `emax_shape()` returns a numeric vector. `emax_shape_area()` returns
#'   a scalar in `[1/2, 1]`. `emax_response()` returns a mean vector, and
#'   `emax_integrated_response()` a scalar in response units.
#' @examples
#' emax_shape(c(0, 0.5, 1), q = 0.75, D = 1)
#' emax_shape_area(0)
#' emax_integrated_response(delta = 2, q = 1)
#' @export
emax_shape <- function(x, q, D) {
  .shape_inputs(x, q, D)
  as.numeric(cpp_shape(x, q, D))
}

#' @rdname emax_shape
#' @export
emax_shape_area <- function(q) {
  .scalar(q, "q", 0, 1)
  cpp_area(q)
}

#' @rdname emax_shape
#' @export
emax_response <- function(x, a, delta, q, D) {
  .scalar(a, "a"); .scalar(delta, "delta", 0)
  a + delta * emax_shape(x, q, D)
}

#' @rdname emax_shape
#' @export
emax_integrated_response <- function(delta, q) {
  .scalar(delta, "delta", 0)
  delta * emax_shape_area(q)
}

#' Shape derivatives with respect to q
#'
#' Derivatives are for the observed-dose shape at fixed `x` and `D`.
#' At `q = 1` these are the one-sided derivatives at each fixed positive dose;
#' zero and maximum doses have zero derivatives. They are not derivatives
#' of the integrated response, whose saturation sensitivity diverges.
#' @inheritParams emax_shape
#' @return List with numeric vectors `first` and `second`.
#' @export
emax_shape_derivatives <- function(x, q, D) {
  .shape_inputs(x, q, D)
  d <- cpp_derivatives(x, q, D)
  lapply(d, as.numeric)
}

.quadratic_hull <- function(coefficients, limits, slack = 0) {
  .vector(coefficients, "coefficients"); .vector(limits, "limits")
  .scalar(slack, "slack", 0)
  if (length(coefficients) != 3L || length(limits) != 2L ||
      limits[1] < 0 || limits[1] > limits[2]) stop("invalid quadratic hull inputs")
  cpp_quadratic(coefficients, limits, slack)
}
