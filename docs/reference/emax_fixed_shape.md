# Fixed-shape slice of an Emax confidence region

Fixed-shape slice of an Emax confidence region

## Usage

``` r
emax_fixed_shape(region, q)
```

## Arguments

- region:

  Object from \[emax_region()\] or \[emax_region_from_summary()\].

- q:

  Scalar shape in \`\[0, 1\]\`.

## Value

List with \`status\`, and when feasible the nonnegative \`delta\`
interval and its integrated-response \`interval\`.
