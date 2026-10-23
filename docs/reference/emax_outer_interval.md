# Analytic outer envelope of an Emax projection

Endpoint and interior-dose contrasts bound the effect and shape without
optimization. A nonempty envelope does not establish feasibility of the
exact model intersection. Use \[emax_project()\] for the refined
interval.

## Usage

``` r
emax_outer_interval(region)
```

## Arguments

- region:

  Object from \[emax_region()\] or \[emax_region_from_summary()\].

## Value

A list with \`status\`, \`interval\`, and, when nonempty, the effect
bounds \`delta\` and shape bounds \`q\`.
