# Small deterministic reproduction; no simulation or random draws.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the local installed-package library")
.libPaths(c(normalizePath(args[1], mustWork = TRUE), .libPaths()))
library(emaxProjection)
z <- read.delim(system.file("extdata", "clinical_dose_summary.tsv", package="emaxProjection"))
expected <- read.delim("reported_results/clinical_intervals.tsv")
expected <- expected[expected$method == "contrast_projection_v1", ]
variance <- sum((z$n - 1) * z$sd^2) / (sum(z$n) - nrow(z))
region <- emax_region_from_summary(z$dose, z$n, z$mean, variance, D=4)
fit <- emax_project(region)
stopifnot(sum(z$n) == 369, fit$converged,
          max(abs(fit$interval - c(expected$lower, expected$upper))) < 1e-8)
reference <- new.env(parent=globalenv())
for (p in c("emax_model.R", "profile_lr_curvature_v2.R", "emax_projection_v1.R",
            "emax_projection_v2.R", "simulation_design_v1.R", "emax_contrast_projection_v1.R"))
  sys.source(file.path("reference_R", p), reference)
old <- reference$simulation_region_v1(z$dose, z$n, z$mean, variance, 4, 0.05)
old <- reference$emax_contrast_region_v1(old)
reference_fit <- reference$emax_area_projection_v2(old)
stopifnot(identical(fit$status, reference_fit$status),
          max(abs(fit$interval - reference_fit$interval)) < 1e-8)
print(fit)
cat("PASS: 369-observation aggregate example; R and C++ agree; zero random draws.\n")
