# Package index

## Model and target

The closed Emax shape and normalized integrated excess response.

- [`emax_shape()`](https://byuzbasi.github.io/emaxProjection/reference/emax_shape.md)
  [`emax_shape_area()`](https://byuzbasi.github.io/emaxProjection/reference/emax_shape.md)
  [`emax_response()`](https://byuzbasi.github.io/emaxProjection/reference/emax_shape.md)
  [`emax_integrated_response()`](https://byuzbasi.github.io/emaxProjection/reference/emax_shape.md)
  : Closed Emax shape and its normalized integrated response
- [`emax_shape_derivatives()`](https://byuzbasi.github.io/emaxProjection/reference/emax_shape_derivatives.md)
  : Shape derivatives with respect to q

## Confidence regions

Exact Gaussian F calibration for fixed replicated doses.

- [`emax_region()`](https://byuzbasi.github.io/emaxProjection/reference/emax_region.md)
  : Construct an exact Gaussian confidence region from replicated doses
- [`emax_region_from_summary()`](https://byuzbasi.github.io/emaxProjection/reference/emax_region_from_summary.md)
  : Construct a region from already computed sufficient statistics

## Projection and diagnostics

Analytic envelopes, fixed-shape slices, and adaptive endpoint
enclosures.

- [`emax_project()`](https://byuzbasi.github.io/emaxProjection/reference/emax_project.md)
  : Project a Gaussian Emax region onto the integrated response
- [`emax_fixed_shape()`](https://byuzbasi.github.io/emaxProjection/reference/emax_fixed_shape.md)
  : Fixed-shape slice of an Emax confidence region
- [`emax_outer_interval()`](https://byuzbasi.github.io/emaxProjection/reference/emax_outer_interval.md)
  : Analytic outer envelope of an Emax projection
- [`print(`*`<emax_projection>`*`)`](https://byuzbasi.github.io/emaxProjection/reference/print.emax_projection.md)
  : Print a projection and its numerical status
- [`plot(`*`<emax_projection>`*`)`](https://byuzbasi.github.io/emaxProjection/reference/plot.emax_projection.md)
  : Plot the endpoint enclosure diagnostics
