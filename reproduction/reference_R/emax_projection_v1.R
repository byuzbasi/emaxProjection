# Additive full-model projection audit. Requires unchanged emax_model.R and
# profile_lr_curvature_v2.R (plr2_ball_interval only). No draws/fits on source().
# Exact-arithmetic enclosures are NOT outward-rounded machine certificates.

emax_projection_setup_v1 <- function(region) {
  .emax_check_region(region)
  w <- region$counts
  stopifnot(length(w) == length(region$x), length(w) == length(region$means),
    all(is.finite(c(w, region$x, region$means, region$radius_squared))),
    all(w > 0), region$radius_squared > 0, region$x[1] == 0,
    tail(region$x, 1) == region$D, all(diff(region$x) > 0))
  root_w <- sqrt(w)
  center <- function(v) root_w * (v - weighted.mean(v, w))
  list(region = region, w = w, root_w = root_w, center = center,
       y = center(region$means), radius = sqrt(region$radius_squared),
       outer = emax_area_outer_interval(region))
}

emax_projection_box_v1 <- function(setup, lo, hi) {
  stopifnot(0 <= lo, lo <= hi, hi <= 1)
  s <- setup
  mid <- lo + (hi-lo)/2
  shape <- function(q) emax_shape(s$region$x, q, s$region$D)
  h <- shape(mid)
  # Weighted centering is an orthogonal projection in these coordinates.
  # Monotonicity bounds every component before centering, including q=1.
  displacement <- pmax(h-shape(lo), shape(hi)-h)
  epsilon <- sqrt(sum(s$w * displacement^2))
  H <- s$center(h)
  expanded_radius <- s$radius + s$outer$delta[2] * epsilon
  delta <- plr2_ball_interval(H, s$y, expanded_radius, s$outer$delta)
  if (is.null(delta)) return(NULL)
  raw <- c(delta[1]*emax_shape_area(lo), delta[2]*emax_shape_area(hi))
  # Roundoff cushions enlarge bounds; they do not certify directed rounding.
  cushion <- 128*.Machine$double.eps*max(1, abs(raw), abs(s$outer$interval))
  interval <- c(max(s$outer$interval[1], raw[1]-cushion),
                min(s$outer$interval[2], raw[2]+cushion))
  if (interval[1] > interval[2]) return(NULL)
  list(lo = lo, hi = hi, lower = interval[1], upper = interval[2],
       delta_lower = delta[1], delta_upper = delta[2], epsilon = epsilon)
}

emax_area_projection_v1 <- function(region, tolerance = 0.001, max_splits = 256L) {
  .emax_scalar(tolerance, "tolerance", 0, strict_lower = TRUE)
  .emax_scalar(max_splits, "max_splits", 0)
  if (max_splits != as.integer(max_splits)) stop("max_splits must be an integer")
  s <- emax_projection_setup_v1(region)
  empty <- function(reason, splits = 0L, evaluations = 0L, trace = data.frame()) list(
    status = "empty_set", interval = c(NA_real_, NA_real_), width = 0,
    lower_endpoint = c(NA_real_, NA_real_), upper_endpoint = c(NA_real_, NA_real_),
    width_bounds = c(0, 0), endpoint_gaps = c(0, 0), converged = TRUE,
    termination = reason, splits = splits, evaluations = evaluations,
    witness_count = 0L, witnesses = data.frame(), trace = trace,
    bound_in_exact_arithmetic = TRUE, outward_rounding_certified = FALSE)
  if (s$outer$status == "empty_outer_set") return(empty("empty_original_outer_envelope"))
  cache <- new.env(parent = emptyenv())
  best <- c(Inf, -Inf)
  witnesses <- list()
  evaluations <- 0L
  sample_slice <- function(q) {
    key <- sprintf("%.17g", q)
    if (exists(key, cache, inherits = FALSE)) return(get(key, cache, inherits = FALSE))
    evaluations <<- evaluations + 1L
    value <- emax_fixed_shape_projection(region, q)
    if (value$status == "empty") { assign(key, FALSE, cache); return(FALSE) }
    H <- s$center(emax_shape(region$x, q, region$D))
    residual <- vapply(value$delta, function(d) sum((s$y-d*H)^2), numeric(1))
    cushion <- 128*.Machine$double.eps*max(1, region$radius_squared, sum(s$y^2))
    if (any(residual > region$radius_squared+cushion))
      stop("fixed-shape witness failed its residual check")
    best <<- c(min(best[1], value$interval[1]), max(best[2], value$interval[2]))
    witnesses[[length(witnesses)+1L]] <<- data.frame(q = q,
      delta_lower = value$delta[1], delta_upper = value$delta[2],
      area_lower = value$interval[1], area_upper = value$interval[2],
      residual_lower = residual[1], residual_upper = residual[2])
    assign(key, TRUE, cache)
    TRUE
  }
  enclose <- function(lo, hi, parent = NULL) {
    feasible <- vapply(unique(c(lo, lo+(hi-lo)/2, hi)), sample_slice, logical(1))
    b <- emax_projection_box_v1(s, lo, hi)
    if (!is.null(b) && !is.null(parent)) {
      b$lower <- max(b$lower, parent$lower)
      b$upper <- min(b$upper, parent$upper)
      if (b$lower > b$upper) b <- NULL
    }
    if (is.null(b) && any(feasible)) stop("shape box excluded a feasible fixed-q slice")
    b
  }
  first <- enclose(s$outer$q[1], s$outer$q[2])
  boxes <- if (is.null(first)) list() else list(first)
  traces <- list()
  splits <- 0L
  bounds <- function() {
    c(min(c(best[1], vapply(boxes, `[[`, numeric(1), "lower"))),
      max(c(best[2], vapply(boxes, `[[`, numeric(1), "upper"))))
  }
  termination <- "budget_exhausted"
  repeat {
    # A box wholly inside the witnessed hull cannot determine either extreme.
    if (all(is.finite(best)) && length(boxes)) boxes <- Filter(function(b)
      b$lower < best[1] || b$upper > best[2], boxes)
    if (!length(boxes) && !all(is.finite(best)))
      return(empty("all_shape_intervals_excluded", splits, evaluations,
                   if (length(traces)) do.call(rbind, traces) else data.frame()))
    ab <- bounds()
    gaps <- if (all(is.finite(best))) pmax(0, c(best[1]-ab[1], ab[2]-best[2])) else c(Inf, Inf)
    traces[[length(traces)+1L]] <- data.frame(splits = splits, active_boxes = length(boxes),
      lower_outer = ab[1], lower_inner = best[1], upper_inner = best[2], upper_outer = ab[2],
      lower_gap = gaps[1], upper_gap = gaps[2], evaluations = evaluations)
    if (max(gaps) <= tolerance) { termination <- "endpoint_tolerance"; break }
    if (splits >= max_splits) break
    j <- if (gaps[1] >= gaps[2]) which.min(vapply(boxes, `[[`, numeric(1), "lower")) else
      which.max(vapply(boxes, `[[`, numeric(1), "upper"))
    b <- boxes[[j]]
    mid <- b$lo+(b$hi-b$lo)/2
    if (mid == b$lo || mid == b$hi) { termination <- "floating_point_resolution"; break }
    left <- enclose(b$lo, mid, b)
    right <- enclose(mid, b$hi, b)
    boxes <- c(boxes[-j], if (!is.null(left)) list(left), if (!is.null(right)) list(right))
    splits <- splits+1L
  }
  known_feasible <- all(is.finite(best))
  # If no feasible slice was found, retain the possible hull and explicit status.
  # Never equate an unsuccessful finite search with an empty model intersection.
  interval <- c(max(s$outer$interval[1], ab[1]), min(s$outer$interval[2], ab[2]))
  cushion <- 256*.Machine$double.eps*max(1, abs(interval), abs(best[is.finite(best)]))
  if (known_feasible && (interval[1] > best[1]+cushion || interval[2] < best[2]-cushion))
    stop("outer projection bounds excluded a feasible witness")
  list(status = if (known_feasible) "nonempty" else "unresolved_feasibility",
    interval = interval, width = diff(interval),
    lower_endpoint = if (known_feasible) c(interval[1], best[1]) else c(NA_real_, NA_real_),
    upper_endpoint = if (known_feasible) c(best[2], interval[2]) else c(NA_real_, NA_real_),
    width_bounds = c(if (known_feasible) max(0, diff(best)) else 0, diff(interval)),
    endpoint_gaps = gaps, converged = known_feasible && max(gaps) <= tolerance,
    termination = termination, splits = splits, evaluations = evaluations,
    witness_count = length(witnesses), witnesses = if (length(witnesses)) do.call(rbind, witnesses) else data.frame(),
    trace = do.call(rbind, traces), bound_in_exact_arithmetic = TRUE,
    outward_rounding_certified = FALSE)
}
