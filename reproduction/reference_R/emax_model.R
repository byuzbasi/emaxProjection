# Emax feasibility functions, v1. Base R only; no fitting or simulation on source().
# See docs/research_protocol_v1.md and manuscript/theory_v1.tex.

.emax_scalar <- function(value, name, lower = -Inf, upper = Inf,
                         strict_lower = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < lower || value > upper || (strict_lower && value == lower)) {
    stop(name, " must be one finite number in its stated domain", call. = FALSE)
  }
  invisible(value)
}

.emax_vector <- function(value, name) {
  if (!is.numeric(value) || !length(value) || any(!is.finite(value))) {
    stop(name, " must be a nonempty finite numeric vector", call. = FALSE)
  }
  invisible(value)
}

emax_mean <- function(x, a, b, K) {
  .emax_vector(x, "x")
  if (any(x < 0)) stop("x must be nonnegative", call. = FALSE)
  .emax_scalar(a, "a")
  .emax_scalar(b, "b", 0)
  .emax_scalar(K, "K", 0, strict_lower = TRUE)
  # plogis avoids overflow in K + x. x = 0 gives exactly the baseline.
  a + b * stats::plogis(log(x) - log(K))
}

emax_saturation_deficit <- function(b, K, D) {
  .emax_scalar(b, "b", 0)
  .emax_scalar(K, "K", 0, strict_lower = TRUE)
  .emax_scalar(D, "D", 0, strict_lower = TRUE)
  z <- log(D) - log(K)
  # b - average effect, evaluated without subtraction of nearly equal numbers.
  if (z >= 0) return(b * exp(-z) * (z + log1p(exp(-z))))
  t <- exp(z)
  if (t == 0) return(b)
  b * log1p(t) / t
}

emax_average_effect <- function(b, K, D) {
  .emax_scalar(b, "b", 0)
  .emax_scalar(K, "K", 0, strict_lower = TRUE)
  .emax_scalar(D, "D", 0, strict_lower = TRUE)
  z <- log(D) - log(K)
  if (z < log(1e-4)) {
    t <- exp(z)
    # Integrated geometric series; omitted term <= t^5/6 for 0 <= t < 1e-4.
    return(b * t * (1/2 + t * (-1/3 + t * (1/4 - t/5))))
  }
  b - emax_saturation_deficit(b, K, D)
}

emax_shape <- function(x, q, D) {
  .emax_scalar(q, "q", 0, 1)
  .emax_scalar(D, "D", 0, strict_lower = TRUE)
  .emax_vector(x, "x")
  if (any(x < 0 | x > D)) stop("x must belong to [0, D]", call. = FALSE)
  u <- x / D
  if (q == 1) return(as.numeric(x > 0))
  u / (1 - q + q * u)
}

emax_shape_area <- function(q) {
  .emax_scalar(q, "q", 0, 1)
  if (q == 1) return(1)
  if (q < 0.01) {
    # G(q) = sum_{k >= 0} q^k/((k+1)(k+2)); remainder < 1.2e-18.
    k <- 0:7
    return(sum(q^k / ((k + 1) * (k + 2))))
  }
  (q + (1 - q) * log1p(-q)) / q^2
}

emax_anchored_mean <- function(x, a, delta, q, D) {
  .emax_scalar(a, "a")
  .emax_scalar(delta, "delta", 0)
  a + delta * emax_shape(x, q, D)
}

emax_anchored_area <- function(delta, q) {
  .emax_scalar(delta, "delta", 0)
  delta * emax_shape_area(q)
}

emax_summarize <- function(x, y, D) {
  .emax_vector(x, "x")
  .emax_vector(y, "y")
  .emax_scalar(D, "D", 0, strict_lower = TRUE)
  if (length(x) != length(y)) stop("x and y lengths differ", call. = FALSE)
  if (any(x < 0 | x > D)) stop("x must belong to [0, D]", call. = FALSE)
  doses <- sort(unique(x))
  if (length(doses) < 3 || doses[1] != 0 || tail(doses, 1) != D) {
    stop("need dose 0, dose D, and at least one distinct interior dose", call. = FALSE)
  }
  # No binning, omission, rounding, or recoding of dose levels.
  groups <- lapply(doses, function(d) y[x == d])
  counts <- vapply(groups, length, integer(1))
  means <- vapply(groups, mean, numeric(1))
  df <- length(y) - length(doses)
  if (df <= 0) stop("positive pure-error residual degrees of freedom required", call. = FALSE)
  sse <- sum(vapply(groups, function(v) sum((v - mean(v))^2), numeric(1)))
  if (sse <= 0) stop("positive pooled within-dose variance required", call. = FALSE)
  list(x = doses, means = means, counts = counts, n = length(y),
       m = length(doses), df = df, pooled_variance = sse / df, D = D)
}

emax_mean_region <- function(x, y, D, alpha) {
  .emax_scalar(alpha, "alpha", 0, 1, strict_lower = TRUE)
  if (alpha == 1) stop("alpha must be smaller than 1", call. = FALSE)
  result <- emax_summarize(x, y, D)
  result$alpha <- alpha
  result$radius_squared <- result$m * result$pooled_variance *
    stats::qf(1 - alpha, result$m, result$df)
  result$construction <- "exact normal pure-error F ellipsoid for dose means"
  class(result) <- "emax_mean_region"
  result
}

.emax_check_region <- function(region) {
  if (!inherits(region, "emax_mean_region")) {
    stop("region must be produced by emax_mean_region()", call. = FALSE)
  }
  invisible(region)
}

emax_region_contains <- function(region, a, delta, q) {
  .emax_check_region(region)
  mu <- emax_anchored_mean(region$x, a, delta, q, region$D)
  sum(region$counts * (region$means - mu)^2) <= region$radius_squared
}

emax_fixed_shape_projection <- function(region, q) {
  .emax_check_region(region)
  h <- emax_shape(region$x, q, region$D)
  w <- region$counts
  hc <- h - weighted.mean(h, w)
  yc <- region$means - weighted.mean(region$means, w)
  shh <- sum(w * hc^2)
  delta_hat <- sum(w * hc * yc) / shh
  remaining <- region$radius_squared - sum(w * (yc - delta_hat * hc)^2)
  if (remaining < 0) return(list(status = "empty", interval = c(NA_real_, NA_real_)))
  delta_range <- delta_hat + c(-1, 1) * sqrt(remaining / shh)
  if (delta_range[2] < 0) return(list(status = "empty", interval = c(NA_real_, NA_real_)))
  delta_range[1] <- max(0, delta_range[1])
  list(status = "fixed_shape_projection", delta = delta_range,
       interval = delta_range * emax_shape_area(q), q = q)
}

emax_area_outer_interval <- function(region) {
  .emax_check_region(region)
  m <- region$m
  nr <- region$counts
  yy <- region$means
  radius <- sqrt(region$radius_squared)
  delta_hat <- yy[m] - yy[1]
  delta_error <- radius * sqrt(1/nr[m] + 1/nr[1])
  dlo <- max(0, delta_hat - delta_error)
  dhi <- delta_hat + delta_error
  empty <- function(reason) list(status = "empty_outer_set", interval = c(NA_real_, NA_real_),
                                  reason = reason, alpha = region$alpha)
  if (dhi < 0) return(empty("endpoint effect constrained to be negative"))
  qlo <- 0
  qhi <- 1
  for (j in 2:(m - 1)) {
    u <- region$x[j] / region$D
    aj <- yy[j] - yy[1]
    ej <- radius * sqrt(1/nr[j] + 1/nr[1])
    alo <- aj - ej
    ahi <- aj + ej
    if (ahi < 0) return(empty("an interior effect constrained to be negative"))
    if (dhi == 0) {
      if (alo > 0) return(empty("zero endpoint effect incompatible with interior effect"))
      next
    }
    rlo <- max(u, alo / dhi)
    rhi <- if (dlo > 0) min(1, ahi / dlo) else 1
    if (rlo > rhi) return(empty("interior effect ratio incompatible with Emax shape"))
    qlo <- max(qlo, (rlo - u) / (rlo * (1 - u)))
    qhi <- min(qhi, (rhi - u) / (rhi * (1 - u)))
  }
  if (qlo > qhi) return(empty("interior dose constraints have disjoint shape ranges"))
  # Clip roundoff at the parameter-domain endpoints, never dose observations.
  qlo <- min(1, max(0, qlo))
  qhi <- min(1, max(0, qhi))
  list(status = "analytic_outer_interval", alpha = region$alpha,
       interval = c(dlo * emax_shape_area(qlo), dhi * emax_shape_area(qhi)),
       delta = c(dlo, dhi), q = c(qlo, qhi),
       construction = "analytic outer envelope of F-region projection; no optimizer",
       note = "nonempty envelope does not certify nonempty exact model intersection")
}

emax_pair_diagnostics <- function(x, counts, b, kappa, sigma, D) {
  .emax_vector(x, "x")
  .emax_vector(counts, "counts")
  .emax_scalar(kappa, "kappa", 0, strict_lower = TRUE)
  .emax_scalar(sigma, "sigma", 0, strict_lower = TRUE)
  .emax_scalar(b, "b", 0)
  .emax_scalar(D, "D", 0, strict_lower = TRUE)
  if (length(x) != length(counts) || any(counts < 1 | counts != floor(counts)) ||
      any(x < 0 | x > D)) stop("invalid fixed design", call. = FALSE)
  # Stable difference h_K - h_2K; avoids cancellation near saturation.
  difference <- b * stats::plogis(log(x) - log(kappa)) *
    stats::plogis(log(2 * kappa) - log(x)) / 2
  area_gap <- emax_saturation_deficit(b, 2 * kappa, D) -
    emax_saturation_deficit(b, kappa, D)
  list(n = sum(counts), kappa = kappa,
       kl = sum(counts * difference^2) / (2 * sigma^2),
       area_gap = area_gap, hidden_point = kappa, hidden_gap = b / 6)
}
