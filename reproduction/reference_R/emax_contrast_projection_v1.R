# Additive contrast-region calibration. Existing Emax functions are unchanged.
emax_contrast_region_v1 <- function(region) {
  .emax_check_region(region)
  stopifnot(region$m >= 3L, region$df > 0, region$pooled_variance > 0,
            region$alpha > 0, region$alpha < 1)
  result <- region
  result$radius_squared <- (region$m - 1) * region$pooled_variance *
    stats::qf(1 - region$alpha, region$m - 1, region$df)
  result$construction <- "exact normal pure-error F region for centered contrasts"
  result$calibration_rank <- region$m - 1L
  result
}

# Private dependency injection reuses the unchanged checkpointed runner.
# Only this environment dispatches the projection to the new calibration.
contrast_runner_environment_v1 <- function() {
  env <- new.env(parent = globalenv())
  paths <- c("R/emax_model.R", "config/simulation_v1.R",
    "R/simulation_design_v1.R", "R/profile_lr_curvature_v2.R",
    "R/emax_projection_v1.R", "R/emax_projection_v2.R",
    "R/simulation_io_v1.R", "R/simulation_runner_v1.R",
    "config/projection_reanalysis_v1.R", "R/projection_reanalysis_runner_v1.R",
    "config/emax_contrast_reanalysis_v1.R")
  for (path in paths) sys.source(path, envir = env)
  # Resolve the region adapter against the private, unchanged model utilities.
  adapter <- emax_contrast_region_v1
  environment(adapter) <- env
  engine <- env$emax_area_projection_v2
  env$emax_area_projection_v2 <- function(region, tolerance = 0.001,
                                         max_splits = 4096L) {
    result <- engine(adapter(region), tolerance, max_splits)
    result$engine <- "emax_contrast_projection_v1"
    result
  }
  env$emax_original_projection_v2 <- engine
  env$emax_contrast_region_v1 <- adapter
  env$projection_reanalysis_sources_v1 <- env$contrast_reanalysis_sources_v1
  env
}
