# Methods, boundaries, and numerical meaning

## Model and integrated response

Let the fixed dose axis be $`[0,D]`$, where $`D>0`$ has dose units. The
baseline $`a`$ is unrestricted, the endpoint excess response $`\delta`$
is nonnegative, and the dimensionless shape $`q`$ lies in $`[0,1]`$. For
normalized dose $`u=x/D`$, the mean response is defined by

``` math
\mu(x)=a+\delta h_q(x/D),\qquad h_q(u)=\frac{u}{1-q+qu}\quad(q<1).
```

Here $`h_q`$ is dimensionless and
[`emax_shape()`](https://byuzbasi.github.io/emaxProjection/reference/emax_shape.md)
evaluates it. The closure $`h_1(0)=0`$, $`h_1(u)=1`$ for $`u>0`$
represents immediate saturation. At $`q=0`$ the shape is linear. These
endpoint definitions are part of the model, not finite-dose
approximations.

The average integrated excess response $`T`$ is defined on this linear
axis by

``` math
T=D^{-1}\int_0^D\{\mu(x)-a\}\,dx=\delta G(q),\qquad
G(q)=\int_0^1 h_q(u)\,du.
```

The change of variables $`u=x/D`$ explains why normalization removes the
factor $`D`$. Thus $`T`$ and $`\delta`$ have response units.
[`emax_shape_area()`](https://byuzbasi.github.io/emaxProjection/reference/emax_shape.md)
computes $`G`$ with its endpoint values $`G(0)=1/2`$ and $`G(1)=1`$,
using the stable series and logarithmic formula from the research
implementation.

## Gaussian confidence region

Suppose $`m\geq3`$ distinct fixed doses satisfy $`0=x_1<\cdots<x_m=D`$,
with positive integer counts $`n_j`$ and total $`n=\sum_{j=1}^m n_j>m`$.
Let $`\bar Y`$ be the vector of dose means,
$`W=\operatorname{diag}(n_j)`$, and $`S^2`$ the pooled within-dose
variance with $`\nu=n-m`$ degrees of freedom. With
$`v=W^{1/2}\mathbf 1/\sqrt n`$ and $`P=I-vv^\mathsf{T}`$, define
centered coordinates $`y=PW^{1/2}\bar Y`$ and $`H(q)=PW^{1/2}h_q`$. Here
$`h_q`$ denotes the vector of values $`h_q(x_j/D)`$ at dose levels
indexed by $`j=1,\ldots,m`$. The identity matrix $`I`$ has size $`m`$,
and $`\mathbf 1`$ is the length-$`m`$ vector of ones.

Under independent Gaussian errors of common variance, the exact centered
F region accepts parameters satisfying

``` math
\|y-\delta H(q)\|^2\leq(m-1)S^2 F_{m-1,\nu}(1-\alpha).
```

The norm is Euclidean, $`\alpha\in(0,1)`$ is the chosen error
probability, and $`F_{r,\nu}(1-\alpha)`$ denotes an F quantile with
numerator degrees of freedom $`r`$. The rank $`m-1`$ follows from
removing the unrestricted baseline.
[`emax_region()`](https://byuzbasi.github.io/emaxProjection/reference/emax_region.md)
constructs this region;
[`emax_project()`](https://byuzbasi.github.io/emaxProjection/reference/emax_project.md)
encloses the hull of its accepted values of $`\delta G(q)`$. The region
pivot is standard Gaussian inference. The accompanying manuscript
develops the integrated-target boundary theory and computational
construction; no new Emax model is claimed.

`calibration = "mean"` instead uses the full mean-region radius
$`m S^2F_{m,\nu}(1-\alpha)`$. The larger radius provides the
manuscript’s full mean-region benchmark after minimizing over the
intercept. Both are finite-sample constructions under the same Gaussian
assumptions.

## Boundary behavior

At a flat curve, $`\delta=0`$ and shape need not be identifiable. At
positive immediate saturation, the integrated target has logarithmic
sensitivity near $`q=1`$, even though fixed positive-dose shapes have
finite one-sided derivatives. The package includes both closures. It
does not replace their inference with an interior-only normal
approximation.

For unknown physical doses, equally spaced ordinal codes describe an
illustrative code-axis integral. They do not establish physical-dose
area, clinical benefit, or a treatment recommendation.

## Numerical enclosure and status

The adaptive search uses fixed-shape quadratic slices, analytic outer
bounds, necessary shape-box inequalities, and feasible witnesses. The
refined search also intersects the original finite-budget enclosure.
Defaults are endpoint tolerance `0.001` and a budget of `4096` splits;
no cap on the statistical effect parameter is introduced.

`endpoint_gaps` reports the outer-to-witness differences. The result
includes the numerical trace and termination reason. The
exact-arithmetic reasoning does not make ordinary floating-point results
formal directed-rounding certificates. `outward_rounding_certified` is
always false.

An empty model intersection has zero reported width and no endpoints. An
unresolved search is not an empty set. Coverage studies must retain
empty sets as noncoverage and report unresolved computations separately.

## Citation and scope

The scientific reference is Bahadır Yüzbaşı’s unpublished manuscript,
*Confidence-region projection for integrated Emax responses: boundary
behavior and finite-sample inference* (2026). No accepted-paper status
or DOI is claimed. The package ports its validated R construction to
RcppArmadillo; it does not establish robustness to heteroscedasticity,
dependence, or non-Gaussian errors.
