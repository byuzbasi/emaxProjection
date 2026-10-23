# New wrappers reuse R/emax_model.R without changing its interface or formulas.
simulation_grid_v1 <- function(cfg) {
  rows <- list()
  add <- function(design, n, delta, q, phase, regime = "") {
    hit <- which(vapply(rows, function(z) z$design == design && z$n == n &&
                         z$delta == delta && z$q == q, logical(1)))
    if (length(hit)) {
      rows[[hit]]$rate_regime <<- regime
      return(invisible(NULL))
    }
    rows[[length(rows) + 1L]] <<- data.frame(
      cell_id = sprintf("c%03d", length(rows) + 1L), design = design, n = n,
      m = length(cfg$designs[[design]]), a = cfg$a, D = cfg$D, sigma = cfg$sigma,
      delta = delta, q = q, truth = emax_anchored_area(delta, q),
      phase = phase, rate_regime = regime, stringsAsFactors = FALSE)
  }
  for (design in names(cfg$designs)) for (n in cfg$main_n) {
    add(design, n, 0, cfg$flat_q, "comparison")
    for (delta in cfg$effects) for (q in cfg$shapes)
      add(design, n, delta, q, "comparison")
  }
  for (n in cfg$rate_n) {
    add("six", n, 0, cfg$flat_q, "rates", "flat")
    for (q in cfg$rate_shapes)
      add("six", n, cfg$rate_delta, q, "rates", paste0("fixed_q_", q))
    for (cc in cfg$local_c) {
      k <- cc / sqrt(n)
      q <- cfg$D / (cfg$D + k)
      add("six", n, cfg$local_b * q, q, "rates", paste0("finite_K_c_", cc))
    }
  }
  grid <- do.call(rbind, rows)
  stopifnot(nrow(grid) == 93L, sum(grid$phase == "comparison") == 66L,
            sum(nzchar(grid$rate_regime)) == 30L, !anyDuplicated(grid$cell_id),
            all(grid$n %% grid$m == 0), all(grid$n > grid$m),
            all(is.finite(grid$truth)), all(grid$delta >= 0),
            all(grid$q >= 0 & grid$q <= 1))
  grid
}

simulation_plan_v1 <- function(cfg, mode) {
  stopifnot(mode %in% c("smoke", "production"))
  all_cells <- simulation_grid_v1(cfg)
  selected <- rep(TRUE, nrow(all_cells))
  if (mode == "smoke") {
    z <- all_cells
    selected <- (z$design == "six" & z$n == 60 & z$delta == 0) |
      (z$design == "three" & z$n == 60 & z$delta == 0.5 & z$q == 0) |
      (z$design == "six" & z$n == 60 & z$delta == 2 & z$q %in% c(0.75, 0.95, 1)) |
      (z$n == 240 & z$rate_regime == "finite_K_c_0.5")
    stopifnot(sum(selected) == 6L)
  }
  repetitions <- if (mode == "smoke") cfg$smoke_repetitions else cfg$repetitions
  chunk <- if (mode == "smoke") 1L else cfg$chunk_size
  do.call(RNGkind, as.list(cfg$rng))
  set.seed(cfg$seed)
  cell_seed <- .Random.seed
  tasks <- list()
  for (i in seq_len(nrow(all_cells))) {
    if (selected[i]) {
      seed <- cell_seed
      for (first in seq.int(1L, repetitions, by = chunk)) {
        last <- min(first + chunk - 1L, repetitions)
        tasks[[length(tasks) + 1L]] <- list(
          task_id = sprintf("%s_r%04d_%04d", all_cells$cell_id[i], first, last),
          cell_id = all_cells$cell_id[i], phase = all_cells$phase[i],
          first = first, last = last, seed = seed)
        for (unused in seq.int(first, last)) seed <- parallel::nextRNGSubStream(seed)
      }
    }
    cell_seed <- parallel::nextRNGStream(cell_seed)
  }
  list(mode = mode, cells = all_cells[selected, ], tasks = tasks,
       repetitions = repetitions, methods = cfg$methods)
}

simulation_region_v1 <- function(x, counts, means, variance, D, alpha) {
  stopifnot(length(x) == length(counts), length(x) == length(means),
            all(is.finite(c(x, counts, means, variance))), variance > 0,
            all(counts >= 1 & counts == as.integer(counts)),
            length(x) >= 3, x[1] == 0, tail(x, 1) == D, all(diff(x) > 0),
            sum(counts) > length(x), alpha > 0, alpha < 1)
  m <- length(x)
  result <- list(x = x, counts = counts, means = means, n = sum(counts), m = m,
                 df = sum(counts) - m, pooled_variance = variance, D = D,
                 alpha = alpha, radius_squared = m * variance *
                   stats::qf(1 - alpha, m, sum(counts) - m),
                 construction = "exact normal pure-error F ellipsoid for dose means")
  class(result) <- "emax_mean_region"
  result
}

simulation_draw_v1 <- function(cell, cfg, seed) {
  assign(".Random.seed", seed, envir = .GlobalEnv)
  doses <- cfg$designs[[cell$design]]
  counts <- rep(as.integer(cell$n / length(doses)), length(doses))
  mu <- emax_anchored_mean(doses, cell$a, cell$delta, cell$q, cell$D)
  means <- mu + cell$sigma * stats::rnorm(length(doses)) / sqrt(counts)
  variance <- cell$sigma^2 * stats::rchisq(1, cell$n - length(doses)) /
    (cell$n - length(doses))
  simulation_region_v1(doses, counts, means, variance, cell$D, cfg$alpha)
}
