# Additive precision refinement for the full-model Emax projection.
# Requires emax_model.R and emax_projection_v1.R. Sourcing performs no fits.
# Bounds are exact-arithmetic arguments evaluated in ordinary double precision;
# they are not directed-rounding certificates.

emax_shape_derivatives_v2 <- function(x, q, D) {
  .emax_scalar(q, "q", 0, 1)
  .emax_scalar(D, "D", 0, strict_lower = TRUE)
  .emax_vector(x, "x")
  if (any(x < 0 | x > D)) stop("x must belong to [0, D]", call. = FALSE)
  u <- x / D
  first <- second <- numeric(length(u))
  interior <- u > 0 & u < 1
  denominator <- 1 - q * (1 - u[interior])
  first[interior] <- u[interior] * (1 - u[interior]) / denominator^2
  second[interior] <- 2 * u[interior] * (1 - u[interior])^2 / denominator^3
  list(first = first, second = second)
}

emax_quadratic_hull_v2 <- function(coefficients, limits, slack = 0) {
  .emax_vector(coefficients, "coefficients")
  .emax_vector(limits, "limits")
  .emax_scalar(slack, "slack", 0)
  if (length(coefficients) != 3L || length(limits) != 2L ||
      limits[1] < 0 || limits[1] > limits[2]) stop("invalid quadratic hull inputs")
  lower <- limits[1]
  upper <- limits[2]
  a <- coefficients[1]
  b <- coefficients[2]
  c <- coefficients[3] - slack
  value <- function(x) (a * x + b) * x + c
  scale <- max(1, abs(coefficients), abs(c), abs(limits),
               abs(a) * upper^2, abs(b) * upper)
  cushion <- 1024 * .Machine$double.eps * scale
  if (lower == upper) return(if (value(lower) <= cushion) limits else NULL)

  if (a == 0) {
    if (abs(b) * max(1, upper) <= cushion)
      return(if (c <= cushion) limits else NULL)
    root <- -c / b
    root_cushion <- 1024 * .Machine$double.eps * max(1, abs(root), abs(limits))
    answer <- if (b > 0) c(lower, min(upper, root + root_cushion)) else
      c(max(lower, root - root_cushion), upper)
    if (answer[1] > answer[2]) return(NULL)
    return(answer)
  }

  discriminant <- b^2 - 4 * a * c
  discriminant_cushion <- 4096 * .Machine$double.eps *
    max(1, b^2, abs(4 * a * c))
  if (discriminant < -discriminant_cushion) {
    return(if (a < 0) limits else NULL)
  }
  if (discriminant < 0) return(limits)
  root_discriminant <- sqrt(discriminant)
  if (root_discriminant == 0) {
    root <- -b / (2 * a)
    if (a > 0) {
      root_cushion <- 2048 * .Machine$double.eps * max(1, abs(root), abs(limits))
      answer <- c(max(lower, root - root_cushion), min(upper, root + root_cushion))
      if (answer[1] > answer[2]) return(NULL)
      return(answer)
    }
    return(limits)
  }
  stable <- -0.5 * (b + if (b >= 0) root_discriminant else -root_discriminant)
  roots <- if (stable == 0) c((-b - root_discriminant) / (2 * a),
                              (-b + root_discriminant) / (2 * a)) else
    c(stable / a, c / stable)
  roots <- sort(roots)
  root_cushion <- 4096 * .Machine$double.eps * max(1, abs(roots), abs(limits))
  if (a > 0) {
    answer <- c(max(lower, roots[1] - root_cushion),
                min(upper, roots[2] + root_cushion))
    if (answer[1] > answer[2]) return(NULL)
    return(answer)
  }

  pieces <- list()
  left <- c(lower, min(upper, roots[1] + root_cushion))
  right <- c(max(lower, roots[2] - root_cushion), upper)
  if (left[1] <= left[2]) pieces[[length(pieces) + 1L]] <- left
  if (right[1] <= right[2]) pieces[[length(pieces) + 1L]] <- right
  if (!length(pieces)) return(NULL)
  c(min(vapply(pieces, `[`, numeric(1), 1L)),
    max(vapply(pieces, `[`, numeric(1), 2L)))
}

emax_projection_box_v2 <- function(setup, lo, hi, parent = NULL) {
  stopifnot(0 <= lo, lo <= hi, hi <= 1)
  s <- setup
  legacy_box <- emax_projection_box_v1(s, lo, hi)
  if (is.null(legacy_box)) return(NULL)
  middle <- lo + (hi - lo) / 2
  half_width <- (hi - lo) / 2
  shape <- function(q) emax_shape(s$region$x, q, s$region$D)
  h <- shape(middle)
  H <- s$center(h)
  derivative <- emax_shape_derivatives_v2(s$region$x, middle, s$region$D)
  endpoint_derivative <- emax_shape_derivatives_v2(s$region$x, hi, s$region$D)
  H_first <- s$center(derivative$first)

  endpoint_epsilon <- max(
    sqrt(sum(s$w * (shape(lo) - h)^2)),
    sqrt(sum(s$w * (shape(hi) - h)^2)))
  second_norm <- sqrt(sum(s$w * endpoint_derivative$second^2))
  taylor_epsilon <- half_width * sqrt(sum(H_first^2)) +
    0.5 * half_width^2 * second_norm
  epsilon <- min(endpoint_epsilon, taylor_epsilon)

  A <- sum(H^2)
  b <- sum(s$y * H)
  y_squared <- sum(s$y^2)
  global_limits <- s$outer$delta
  scale <- max(1, y_squared, s$region$radius_squared,
               A * global_limits[2]^2, 2 * abs(b) * global_limits[2])
  slack <- 4096 * .Machine$double.eps * scale

  self_scaled <- emax_quadratic_hull_v2(
    c(A - epsilon^2, -2 * (b + s$radius * epsilon),
      y_squared - s$region$radius_squared), global_limits, slack)
  if (is.null(self_scaled)) return(NULL)

  first_norm <- sqrt(sum(s$w * endpoint_derivative$first^2))
  H_bound <- sqrt(A) + epsilon
  A_first <- 2 * sum(H * H_first)
  A_second_bound <- 2 * first_norm^2 + 2 * H_bound * second_norm
  b_first <- sum(s$y * H_first)
  b_second_bound <- sum(abs(s$y) * s$root_w * endpoint_derivative$second)
  A_error <- abs(A_first) * half_width +
    0.5 * A_second_bound * half_width^2
  b_error <- abs(b_first) * half_width +
    0.5 * b_second_bound * half_width^2
  A_lower <- max(0, A - A_error)
  b_upper <- b + b_error
  coefficient_scaled <- emax_quadratic_hull_v2(
    c(A_lower, -2 * b_upper, y_squared - s$region$radius_squared),
    global_limits, slack)
  if (is.null(coefficient_scaled)) return(NULL)

  delta <- c(max(global_limits[1], legacy_box$delta_lower,
                 self_scaled[1], coefficient_scaled[1]),
             min(global_limits[2], legacy_box$delta_upper,
                 self_scaled[2], coefficient_scaled[2]))
  if (!is.null(parent)) {
    delta <- c(max(delta[1], parent$delta_lower),
               min(delta[2], parent$delta_upper))
  }
  if (delta[1] > delta[2]) return(NULL)
  raw <- c(delta[1] * emax_shape_area(lo),
           delta[2] * emax_shape_area(hi))
  area_cushion <- 512 * .Machine$double.eps *
    max(1, abs(raw), abs(s$outer$interval))
  interval <- c(max(s$outer$interval[1], legacy_box$lower,
                    raw[1] - area_cushion),
                min(s$outer$interval[2], legacy_box$upper,
                    raw[2] + area_cushion))
  if (!is.null(parent)) {
    interval <- c(max(interval[1], parent$lower),
                  min(interval[2], parent$upper))
  }
  if (interval[1] > interval[2]) return(NULL)
  list(lo = lo, hi = hi, lower = interval[1], upper = interval[2],
       delta_lower = delta[1], delta_upper = delta[2], epsilon = epsilon,
       endpoint_epsilon = endpoint_epsilon, taylor_epsilon = taylor_epsilon,
       A_lower = A_lower, b_upper = b_upper)
}

emax_area_projection_v2 <- function(region, tolerance = 0.001,
                                    max_splits = 4096L) {
  .emax_scalar(tolerance, "tolerance", 0, strict_lower = TRUE)
  .emax_scalar(max_splits, "max_splits", 0)
  if (max_splits != as.integer(max_splits)) stop("max_splits must be an integer")
  s <- emax_projection_setup_v1(region)
  empty <- function(reason, splits = 0L, evaluations = 0L,
                    trace = data.frame()) list(
    status = "empty_set", interval = c(NA_real_, NA_real_), width = 0,
    lower_endpoint = c(NA_real_, NA_real_),
    upper_endpoint = c(NA_real_, NA_real_), width_bounds = c(0, 0),
    endpoint_gaps = c(0, 0), converged = TRUE, termination = reason,
    splits = splits, evaluations = evaluations, witness_count = 0L,
    witnesses = data.frame(), trace = trace,
    bound_in_exact_arithmetic = TRUE, outward_rounding_certified = FALSE,
    engine = "emax_projection_v2")
  if (s$outer$status == "empty_outer_set")
    return(empty("empty_original_outer_envelope"))

  cache <- new.env(parent = emptyenv())
  best <- c(Inf, -Inf)
  witnesses <- list()
  evaluations <- 0L
  sample_slice <- function(q) {
    key <- sprintf("%.17g", q)
    if (exists(key, cache, inherits = FALSE))
      return(get(key, cache, inherits = FALSE))
    evaluations <<- evaluations + 1L
    value <- emax_fixed_shape_projection(region, q)
    if (value$status == "empty") {
      assign(key, FALSE, cache)
      return(FALSE)
    }
    H <- s$center(emax_shape(region$x, q, region$D))
    residual <- vapply(value$delta,
      function(delta) sum((s$y - delta * H)^2), numeric(1))
    cushion <- 512 * .Machine$double.eps *
      max(1, region$radius_squared, sum(s$y^2))
    if (any(residual > region$radius_squared + cushion))
      stop("fixed-shape witness failed its residual check")
    best <<- c(min(best[1], value$interval[1]),
               max(best[2], value$interval[2]))
    witnesses[[length(witnesses) + 1L]] <<- data.frame(
      q = q, delta_lower = value$delta[1], delta_upper = value$delta[2],
      area_lower = value$interval[1], area_upper = value$interval[2],
      residual_lower = residual[1], residual_upper = residual[2])
    assign(key, TRUE, cache)
    TRUE
  }
  enclose <- function(lo, hi, parent = NULL) {
    feasible <- vapply(unique(c(lo, lo + (hi - lo) / 2, hi)),
                       sample_slice, logical(1))
    box <- emax_projection_box_v2(s, lo, hi, parent)
    if (is.null(box) && any(feasible))
      stop("v2 shape box excluded a feasible fixed-q slice")
    box
  }

  first <- enclose(s$outer$q[1], s$outer$q[2])
  boxes <- if (is.null(first)) list() else list(first)
  traces <- list()
  splits <- 0L
  bounds <- function() c(
    min(c(best[1], vapply(boxes, `[[`, numeric(1), "lower"))),
    max(c(best[2], vapply(boxes, `[[`, numeric(1), "upper"))))
  termination <- "budget_exhausted"
  repeat {
    if (all(is.finite(best)) && length(boxes))
      boxes <- Filter(function(box)
        box$lower < best[1] || box$upper > best[2], boxes)
    if (!length(boxes) && !all(is.finite(best)))
      return(empty("all_shape_intervals_excluded", splits, evaluations,
        if (length(traces)) do.call(rbind, traces) else data.frame()))
    outer_bounds <- bounds()
    gaps <- if (all(is.finite(best)))
      pmax(0, c(best[1] - outer_bounds[1], outer_bounds[2] - best[2])) else
      c(Inf, Inf)
    traces[[length(traces) + 1L]] <- data.frame(
      splits = splits, active_boxes = length(boxes),
      lower_outer = outer_bounds[1], lower_inner = best[1],
      upper_inner = best[2], upper_outer = outer_bounds[2],
      lower_gap = gaps[1], upper_gap = gaps[2], evaluations = evaluations)
    if (max(gaps) <= tolerance) {
      termination <- "endpoint_tolerance"
      break
    }
    if (splits >= max_splits) break
    index <- if (gaps[1] >= gaps[2])
      which.min(vapply(boxes, `[[`, numeric(1), "lower")) else
      which.max(vapply(boxes, `[[`, numeric(1), "upper"))
    box <- boxes[[index]]
    middle <- box$lo + (box$hi - box$lo) / 2
    if (middle == box$lo || middle == box$hi) {
      termination <- "floating_point_resolution"
      break
    }
    left <- enclose(box$lo, middle, box)
    right <- enclose(middle, box$hi, box)
    boxes <- c(boxes[-index], if (!is.null(left)) list(left),
               if (!is.null(right)) list(right))
    splits <- splits + 1L
  }

  known_feasible <- all(is.finite(best))
  interval <- c(max(s$outer$interval[1], outer_bounds[1]),
                min(s$outer$interval[2], outer_bounds[2]))
  # Intersect the refined search with the independently valid v1 finite-budget
  # enclosure. This preserves the approved v2-within-v1 validation invariant.
  reference_splits <- min(256L, as.integer(max_splits))
  reference <- emax_area_projection_v1(region, tolerance, reference_splits)
  if (reference$status == "empty_set") {
    if (known_feasible) stop("v1 reference declared empty despite a v2 witness")
  } else {
    interval <- c(max(interval[1], reference$interval[1]),
                  min(interval[2], reference$interval[2]))
  }
  cushion <- 1024 * .Machine$double.eps *
    max(1, abs(interval), abs(best[is.finite(best)]))
  if (known_feasible &&
      (interval[1] > best[1] + cushion || interval[2] < best[2] - cushion))
    stop("v2 outer projection bounds excluded a feasible witness")
  gaps <- if (known_feasible)
    pmax(0, c(best[1] - interval[1], interval[2] - best[2])) else c(Inf, Inf)
  if (nrow(traces[[length(traces)]])) {
    traces[[length(traces)]]$lower_outer <- interval[1]
    traces[[length(traces)]]$upper_outer <- interval[2]
    traces[[length(traces)]]$lower_gap <- gaps[1]
    traces[[length(traces)]]$upper_gap <- gaps[2]
  }
  list(status = if (known_feasible) "nonempty" else "unresolved_feasibility",
    interval = interval, width = diff(interval),
    lower_endpoint = if (known_feasible) c(interval[1], best[1]) else
      c(NA_real_, NA_real_),
    upper_endpoint = if (known_feasible) c(best[2], interval[2]) else
      c(NA_real_, NA_real_),
    width_bounds = c(if (known_feasible) max(0, diff(best)) else 0,
                     diff(interval)), endpoint_gaps = gaps,
    converged = known_feasible && max(gaps) <= tolerance,
    termination = termination, splits = splits, evaluations = evaluations,
    witness_count = length(witnesses),
    witnesses = if (length(witnesses)) do.call(rbind, witnesses) else data.frame(),
    trace = do.call(rbind, traces), bound_in_exact_arithmetic = TRUE,
    outward_rounding_certified = FALSE, engine = "emax_projection_v2",
    v1_reference_splits = reference_splits)
}
