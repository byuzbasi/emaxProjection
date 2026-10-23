# Construct a region from already computed sufficient statistics

Uses the same assumptions and calibrations as \[emax_region()\]. The
variance must be the pooled within-dose variance with \`sum(counts) -
length(x)\` degrees of freedom, not a standard error, model residual
variance, or population variance. Counts are positive integers.
Different replication counts at different doses are supported.

## Usage

``` r
emax_region_from_summary(
  x,
  counts,
  means,
  pooled_variance,
  D = max(x),
  alpha = 0.05,
  calibration = c("contrast", "mean")
)
```

## Arguments

- x:

  Nonempty numeric vector of observed doses.

- counts:

  Replication counts, one per ordered distinct dose.

- means:

  Dose-specific sample means, in the same order as \`x\`.

- pooled_variance:

  Positive finite pooled within-dose variance.

- D:

  Positive upper dose endpoint. Defaults to the maximum supplied dose.

- alpha:

  Error probability strictly between zero and one.

- calibration:

  Either \`"contrast"\` (primary method) or \`"mean"\` (full mean-region
  benchmark).

## Value

An \`emax_region\` object.
