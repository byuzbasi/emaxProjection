# Validation scope: 0.1.0

The RcppArmadillo implementation was validated locally on Apple arm64
with R 4.6.0, Apple clang 17, Rcpp 1.1.1.1.1, and RcppArmadillo
15.2.7.1.

## Targeted equivalence gate

- 32 deterministic end-to-end comparisons against the preserved R
  engine.
- Both centered contrast and full mean-region calibration.
- Three- and six-dose stored validation inputs and hand-specified
  boundary, narrow, empty, budget-exhausted, and unresolved-feasibility
  cases.
- The maximum absolute difference between interval endpoints was
  **8.881784197001249e-16**; the required comparison tolerance was 1e-9.
- Status, convergence, and termination agreed in every comparison.
- Analytic envelopes and fixed-shape slices agreed within the required
  tolerance.
- Shape, integrated response, derivatives, and quadratic hull primitives
  were checked against the R reference, including the closed endpoints.
- The saved 369-observation clinical interval was reproduced without
  changing the data or fitting any comparator.
- Projection calls preserved the random-number state; no random draws
  were used.
- Historical source and selected input/result seals were unchanged.

## Package tests

Tests cover independent numerical quadrature, analytic model closures,
the Gaussian F calibration, raw/summary consistency, missing-input
rejection, corrupted-region rejection, target containment at boundaries,
numerical status, and invariance to a baseline shift and a dose-unit
change.

The small tests and documentation examples are not new Monte Carlo
evidence. No production run or speed benchmark was performed. Windows
and Linux compilation have not yet been tested. The enclosing
inequalities are exact-arithmetic arguments evaluated in ordinary double
precision; no formal directed-rounding certificate or independent
external proof review is claimed.
