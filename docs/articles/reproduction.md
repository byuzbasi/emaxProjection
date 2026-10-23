# Reproduction and validation scope

## What is reproduced here

The package examples evaluate the mathematical shape, project small
deterministic regions, and reproduce one clinical scalar interval from
aggregate summaries. They do not rerun the manuscript’s completed Monte
Carlo study or recompute its benchmark methods. The completed study used
93 design cells and 5,000 replications per cell, for 465,000 datasets.
Porting the numerical core does not replace those completed research
outputs with new results.

## Check the package

With declared dependencies installed and a C++17 compiler available:

``` sh
R CMD build emaxProjection
R CMD check --no-manual emaxProjection_0.1.0.tar.gz
```

These commands build the package, its small vignettes, and its tests.
The source package does not launch production simulations, bootstrap
studies, benchmark sweeps, or model training.

## Equivalence validation

The local port was compared with the preserved R research engine on 32
deterministic end-to-end cases, covering both calibrations, three- and
six-dose inputs, flat and immediate-saturation boundaries, narrow
intervals, empty intersections, exhausted budgets, and unresolved
feasibility. The clinical interval was also compared with the saved
manuscript result. This is a targeted equivalence gate, not exhaustive
validation of all floating-point inputs or operating systems.

`VALIDATION.md` states the measured agreement and tested platform. The
historical R implementation remains the reference for the completed
study. Package tests additionally verify analytic endpoints, numerical
quadrature, F calibration, input validation, and baseline/dose-unit
invariance.

## Record a new analysis

Keep the original observations, input units, variance definition,
calibration, alpha, tolerance, split budget, package version,
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html), result
status, endpoint gaps, and result object. Save to a new location rather
than overwriting an earlier analysis. If stochastic sampling is part of
a separate study, freeze its seed and record that design separately: the
projection functions themselves do not draw random numbers.

## Scientific limits

The assumptions and finite-budget diagnostics remain essential. Neither
the package nor its test suite proves coverage under heteroscedastic or
dependent errors. Ordinary floating-point cushions do not establish a
directed-rounding machine certificate. No speed advantage is claimed
without a measured study.

## Independent reference implementation

The [repository reproduction
folder](https://github.com/byuzbasi/emaxProjection/tree/main/reproduction)
contains the original R projection functions, completed summary tables
and a small clinical comparison script. It does not launch the
production study.
