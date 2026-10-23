# Project a Gaussian Emax region onto the integrated response

Uses an adaptive search over the entire closed shape interval, exact
fixed-shape quadratic slices, necessary box bounds, and feasible
witnesses. The refined enclosure is intersected with the original
finite-budget enclosure, using at most 256 reference splits. Tolerance
and split-budget defaults match the completed research implementation.

## Usage

``` r
emax_project(region, tolerance = 0.001, max_splits = 4096L)
```

## Arguments

- region:

  Object from \[emax_region()\] or \[emax_region_from_summary()\].

- tolerance:

  Positive absolute endpoint-gap target, in response units.

- max_splits:

  Nonnegative integer adaptive split budget.

## Value

An \`emax_projection\` list with the outer \`interval\`, endpoint
brackets, \`width_bounds\`, \`endpoint_gaps\`, \`status\`,
\`converged\`, \`termination\`, \`witnesses\`, and a diagnostic
\`trace\`. Empty-set width is zero and interval endpoints are missing,
matching the research convention. An empty set counts as noncoverage in
coverage studies.

## Details

Bounds hold by the exact-arithmetic argument and are evaluated in
ordinary double precision with roundoff cushions. They are not formal
machine certificates: \`outward_rounding_certified\` is always
\`FALSE\`.

An empty model intersection has no confidence interval. A finite search
without a feasible witness returns \`unresolved_feasibility\`, not an
empty set. Inspect both \`status\` and \`converged\`, and report
unresolved results.

## Examples

``` r
region <- emax_region_from_summary(c(0, 0.5, 1), rep(20, 3),
  c(0.2, 0.6, 0.8), pooled_variance = 1)
result <- emax_project(region)
result
#> Integrated Emax response | contrast calibration
#> Status: nonempty; termination: endpoint_tolerance
#> Outer interval: [0, 1.16684]
#> Endpoint gaps: 0.000000000, 0.000637259; converged: TRUE
#> Ordinary double precision; no directed-rounding certificate.
```
