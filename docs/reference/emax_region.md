# Construct an exact Gaussian confidence region from replicated doses

Assumes fixed distinct doses including zero and \`D\`, independent
Gaussian observations with a common unknown positive variance, at least
three dose levels, and positive pure-error degrees of freedom. No
observations are removed, imputed, binned, recoded, or transformed.
Within-dose variation supplies the pooled variance estimate.

## Usage

``` r
emax_region(
  x,
  y,
  D = max(x),
  alpha = 0.05,
  calibration = c("contrast", "mean")
)
```

## Arguments

- x:

  Nonempty numeric vector of observed doses.

- y:

  Corresponding numeric response vector, without missing values.

- D:

  Positive upper dose endpoint. Defaults to the maximum supplied dose.

- alpha:

  Error probability strictly between zero and one.

- calibration:

  Either \`"contrast"\` (primary method) or \`"mean"\` (full mean-region
  benchmark).

## Value

An \`emax_region\` object containing the unchanged sufficient
statistics, calibration rank, and squared F-region radius.

## Details

Contrast calibration removes the unrestricted intercept and has rank
\`m - 1\`; mean calibration uses the full \`m\`-dimensional mean region.
Both use an exact F pivot under the stated assumptions. Projection onto
the integrated response can be conservative. The standard contrast pivot
itself is not a claim of a new statistical model.

## Examples

``` r
x <- rep(c(0, 0.5, 1), each = 4)
y <- rep(c(0.2, 0.6, 0.8), each = 4) + rep(c(-0.1, 0.1, -0.2, 0.2), 3)
region <- emax_region(x, y)
region$calibration_rank
#> [1] 2
```
