# Shape derivatives with respect to q

Derivatives are for the observed-dose shape at fixed \`x\` and \`D\`. At
\`q = 1\` these are the one-sided derivatives at each fixed positive
dose; zero and maximum doses have zero derivatives. They are not
derivatives of the integrated response, whose saturation sensitivity
diverges.

## Usage

``` r
emax_shape_derivatives(x, q, D)
```

## Arguments

- x:

  Numeric nonempty dose vector in \`\[0, D\]\`.

- q:

  Scalar shape in \`\[0, 1\]\`.

- D:

  Positive dose-axis upper endpoint, in the same units as \`x\`.

## Value

List with numeric vectors \`first\` and \`second\`.
