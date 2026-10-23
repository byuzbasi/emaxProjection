#' Construct an exact Gaussian confidence region from replicated doses
#'
#' Assumes fixed distinct doses including zero and `D`, independent Gaussian
#' observations with a common unknown positive variance, at least three dose
#' levels, and positive pure-error degrees of freedom. No observations are
#' removed, imputed, binned, recoded, or transformed. Within-dose variation
#' supplies the pooled variance estimate.
#'
#' Contrast calibration removes the unrestricted intercept and has rank
#' `m - 1`; mean calibration uses the full `m`-dimensional mean region.
#' Both use an exact F pivot under the stated assumptions. Projection onto
#' the integrated response can be conservative. The standard contrast pivot
#' itself is not a claim of a new statistical model.
#' @param x Nonempty numeric vector of observed doses.
#' @param y Corresponding numeric response vector, without missing values.
#' @param D Positive upper dose endpoint. Defaults to the maximum supplied dose.
#' @param alpha Error probability strictly between zero and one.
#' @param calibration Either `"contrast"` (primary method) or `"mean"`
#'   (full mean-region benchmark).
#' @return An `emax_region` object containing the unchanged sufficient
#'   statistics, calibration rank, and squared F-region radius.
#' @examples
#' x <- rep(c(0, 0.5, 1), each = 4)
#' y <- rep(c(0.2, 0.6, 0.8), each = 4) + rep(c(-0.1, 0.1, -0.2, 0.2), 3)
#' region <- emax_region(x, y)
#' region$calibration_rank
#' @export
emax_region <- function(x, y, D = max(x), alpha = 0.05,
                        calibration = c("contrast", "mean")) {
  .vector(x, "x"); .vector(y, "y"); .scalar(D, "D", positive = TRUE)
  if (length(x) != length(y)) stop("x and y lengths differ", call. = FALSE)
  if (any(x < 0 | x > D)) stop("x must belong to [0, D]", call. = FALSE)
  doses <- sort(unique(x))
  groups <- lapply(doses, function(d) y[x == d])
  counts <- vapply(groups, length, integer(1))
  means <- vapply(groups, mean, numeric(1))
  df <- length(y) - length(doses)
  if (df <= 0) stop("positive pure-error degrees of freedom required", call. = FALSE)
  sse <- sum(vapply(groups, function(v) sum((v - mean(v))^2), numeric(1)))
  emax_region_from_summary(doses, counts, means, sse / df, D, alpha, calibration)
}

#' Construct a region from already computed sufficient statistics
#'
#' Uses the same assumptions and calibrations as [emax_region()]. The variance
#' must be the pooled within-dose variance with `sum(counts) - length(x)`
#' degrees of freedom, not a standard error, model residual variance, or
#' population variance. Counts are positive integers. Different replication
#' counts at different doses are supported.
#' @inheritParams emax_region
#' @param counts Replication counts, one per ordered distinct dose.
#' @param means Dose-specific sample means, in the same order as `x`.
#' @param pooled_variance Positive finite pooled within-dose variance.
#' @return An `emax_region` object.
#' @export
emax_region_from_summary <- function(x, counts, means, pooled_variance,
    D = max(x), alpha = 0.05, calibration = c("contrast", "mean")) {
  .vector(x, "x"); .vector(counts, "counts"); .vector(means, "means")
  .scalar(D, "D", positive = TRUE)
  .scalar(alpha, "alpha", 0, 1)
  .scalar(pooled_variance, "pooled_variance", positive = TRUE)
  if (alpha == 0 || alpha == 1) stop("alpha must be strictly between 0 and 1", call. = FALSE)
  calibration <- match.arg(calibration)
  m <- length(x)
  if (m < 3L || length(counts) != m || length(means) != m ||
      x[1] != 0 || x[m] != D || any(diff(x) <= 0) ||
      any(counts < 1 | counts != floor(counts)))
    stop("need ordered distinct doses from 0 to D and positive integer counts", call. = FALSE)
  n <- sum(counts)
  if (!is.finite(n) || n <= m) stop("positive pure-error degrees of freedom required", call. = FALSE)
  rank <- if (calibration == "contrast") m - 1L else m
  radius <- rank * pooled_variance * stats::qf(1 - alpha, rank, n - m)
  if (!is.finite(radius) || radius <= 0)
    stop("F-region radius is not representable as a positive finite double", call. = FALSE)
  structure(list(x = x, counts = counts, means = means, n = n, m = m,
    df = n - m, pooled_variance = pooled_variance, D = D, alpha = alpha,
    radius_squared = radius, calibration = calibration, calibration_rank = rank,
    construction = if (calibration == "contrast")
      "exact normal pure-error F region for centered contrasts" else
      "exact normal pure-error F ellipsoid for dose means"), class = "emax_region")
}

.check_region <- function(region) {
  if (!inherits(region, "emax_region"))
    stop("region must be constructed by emax_region() or emax_region_from_summary()", call. = FALSE)
  fields <- c("x", "counts", "means", "pooled_variance", "D", "alpha", "calibration",
    "radius_squared", "calibration_rank", "n", "m", "df")
  if (!all(fields %in% names(region))) stop("region is missing required fields", call. = FALSE)
  validated <- emax_region_from_summary(region$x, region$counts, region$means,
    region$pooled_variance, region$D, region$alpha, region$calibration)
  for (field in c("radius_squared", "calibration_rank", "n", "m", "df")) {
    .scalar(region[[field]], field)
    if (!isTRUE(all.equal(region[[field]], validated[[field]], tolerance = 1e-14)))
      stop("region has inconsistent ", field, "; reconstruct it from its inputs", call. = FALSE)
  }
  invisible(region)
}
