# Closed Emax shape and its normalized integrated response

The closed model has mean \`a + delta \* emax_shape(x, q, D)\`, with an
unrestricted baseline \`a\`, nonnegative endpoint effect \`delta\`, and
shape \`q\` in \`\[0, 1\]\`. The normalized integrated excess response
is \`delta \* emax_shape_area(q)\`. The endpoint \`q = 0\` is linear;
\`q = 1\` represents immediate saturation at every positive dose, with
zero shape at dose zero. For \`0 \< q \< 1\`, ordinary Emax parameters
satisfy \`K = D \* (1 - q) / q\` and \`b = delta / q\`.

## Usage

``` r
emax_shape(x, q, D)

emax_shape_area(q)

emax_response(x, a, delta, q, D)

emax_integrated_response(delta, q)
```

## Arguments

- x:

  Numeric nonempty dose vector in \`\[0, D\]\`.

- q:

  Scalar shape in \`\[0, 1\]\`.

- D:

  Positive dose-axis upper endpoint, in the same units as \`x\`.

- a:

  Finite baseline response.

- delta:

  Nonnegative endpoint excess response, in response units.

## Value

\`emax_shape()\` returns a numeric vector. \`emax_shape_area()\` returns
a scalar in \`\[1/2, 1\]\`. \`emax_response()\` returns a mean vector,
and \`emax_integrated_response()\` a scalar in response units.

## Details

Integration uses the specified linear dose axis. Arbitrary equally
spaced codes do not recover an integral on unknown physical doses.

## Examples

``` r
emax_shape(c(0, 0.5, 1), q = 0.75, D = 1)
#> [1] 0.0 0.8 1.0
emax_shape_area(0)
#> [1] 0.5
emax_integrated_response(delta = 2, q = 1)
#> [1] 2
```
