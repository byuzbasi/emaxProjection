#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

namespace {
constexpr double eps = std::numeric_limits<double>::epsilon();
constexpr double inf = std::numeric_limits<double>::infinity();

// Sequential reductions avoid an additional BLAS parallelism layer. Products
// are doubles, accumulated in long double as in R's sum on supported platforms.
double total(const arma::vec& x) {
  long double s = 0;
  for (double v : x) s += v;
  return static_cast<double>(s);
}
double dot(const arma::vec& a, const arma::vec& b) {
  return total(arma::vec(a % b));
}
double sq(const arma::vec& x) { return dot(x, x); }
double norm(const arma::vec& x) { return std::sqrt(sq(x)); }
double scale(std::initializer_list<double> values) {
  double ans = 1;
  for (double v : values) ans = std::max(ans, std::abs(v));
  return ans;
}
arma::vec shape(const arma::vec& x, double q, double D) {
  arma::vec ans(x.n_elem);
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    const double u = x[j] / D;
    ans[j] = q == 1 ? static_cast<double>(x[j] > 0) : u / (1 - q + q * u);
  }
  return ans;
}
double area(double q) {
  if (q == 1) return 1;
  if (q < 0.01) {
    arma::vec terms(8);
    for (int k = 0; k < 8; ++k)
      terms[k] = std::pow(q, k) / ((k + 1.) * (k + 2.));
    return total(terms);
  }
  return (q + (1 - q) * std::log1p(-q)) / (q * q);
}
struct Derivatives { arma::vec first, second; };
Derivatives derivatives(const arma::vec& x, double q, double D) {
  Derivatives d{arma::vec(x.n_elem, arma::fill::zeros),
                arma::vec(x.n_elem, arma::fill::zeros)};
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    const double u = x[j] / D;
    if (u > 0 && u < 1) {
      const double den = 1 - q * (1 - u);
      d.first[j] = u * (1 - u) / std::pow(den, 2);
      d.second[j] = 2 * u * std::pow(1 - u, 2) / std::pow(den, 3);
    }
  }
  return d;
}
struct Interval { double lo = 0, hi = 0; bool valid = false; };
Interval interval(double lo, double hi) { return {lo, hi, lo <= hi}; }

Interval quadratic(const arma::vec& co, Interval limits, double slack) {
  const double a = co[0], b = co[1], c = co[2] - slack;
  const double lower = limits.lo, upper = limits.hi;
  auto value = [&](double x) { return (a * x + b) * x + c; };
  const double cushion = 1024 * eps * scale({a, b, co[2], c, lower, upper,
    std::abs(a) * upper * upper, std::abs(b) * upper});
  if (lower == upper) return value(lower) <= cushion ? limits : Interval{};
  if (a == 0) {
    if (std::abs(b) * std::max(1., upper) <= cushion)
      return c <= cushion ? limits : Interval{};
    const double root = -c / b;
    const double rc = 1024 * eps * scale({root, lower, upper});
    return b > 0 ? interval(lower, std::min(upper, root + rc)) :
      interval(std::max(lower, root - rc), upper);
  }
  const double disc = b * b - 4 * a * c;
  const double dc = 4096 * eps * scale({b * b, 4 * a * c});
  if (disc < -dc) return a < 0 ? limits : Interval{};
  if (disc < 0) return limits;
  const double rd = std::sqrt(disc);
  if (rd == 0) {
    if (a <= 0) return limits;
    const double root = -b / (2 * a);
    const double rc = 2048 * eps * scale({root, lower, upper});
    return interval(std::max(lower, root - rc), std::min(upper, root + rc));
  }
  const double stable = -0.5 * (b + (b >= 0 ? rd : -rd));
  double r1 = stable == 0 ? (-b - rd) / (2 * a) : stable / a;
  double r2 = stable == 0 ? (-b + rd) / (2 * a) : c / stable;
  if (r1 > r2) std::swap(r1, r2);
  const double rc = 4096 * eps * scale({r1, r2, lower, upper});
  if (a > 0) return interval(std::max(lower, r1 - rc), std::min(upper, r2 + rc));
  Interval left = interval(lower, std::min(upper, r1 + rc));
  Interval right = interval(std::max(lower, r2 - rc), upper);
  if (!left.valid) return right;
  if (!right.valid) return left;
  return interval(std::min(left.lo, right.lo), std::max(left.hi, right.hi));
}
struct Region {
  arma::vec x, w, means, rw, y;
  double D, r2, radius, alpha;
  explicit Region(const Rcpp::List& r) :
    x(Rcpp::as<arma::vec>(r["x"])), w(Rcpp::as<arma::vec>(r["counts"])),
    means(Rcpp::as<arma::vec>(r["means"])), rw(arma::sqrt(w)),
    D(Rcpp::as<double>(r["D"])), r2(Rcpp::as<double>(r["radius_squared"])),
    radius(std::sqrt(r2)), alpha(Rcpp::as<double>(r["alpha"])) { y = center(means); }
  arma::vec center(const arma::vec& z) const {
    return rw % (z - dot(z, w) / total(w));
  }
};
struct Outer {
  Interval delta, q, target;
  bool valid = false;
  std::string reason;
};
Outer outer(const Region& s) {
  const arma::uword m = s.x.n_elem;
  const double dhat = s.means[m - 1] - s.means[0];
  const double error = s.radius * std::sqrt(1 / s.w[m - 1] + 1 / s.w[0]);
  const double dlo = std::max(0., dhat - error), dhi = dhat + error;
  Outer ans;
  if (dhi < 0) { ans.reason = "endpoint effect constrained to be negative"; return ans; }
  double qlo = 0, qhi = 1;
  for (arma::uword j = 1; j + 1 < m; ++j) {
    const double u = s.x[j] / s.D, aj = s.means[j] - s.means[0];
    const double ej = s.radius * std::sqrt(1 / s.w[j] + 1 / s.w[0]);
    const double alo = aj - ej, ahi = aj + ej;
    if (ahi < 0) { ans.reason = "an interior effect constrained to be negative"; return ans; }
    if (dhi == 0) {
      if (alo > 0) { ans.reason = "zero endpoint effect incompatible with interior effect"; return ans; }
      continue;
    }
    const double rlo = std::max(u, alo / dhi);
    const double rhi = dlo > 0 ? std::min(1., ahi / dlo) : 1;
    if (rlo > rhi) { ans.reason = "interior effect ratio incompatible with Emax shape"; return ans; }
    qlo = std::max(qlo, (rlo - u) / (rlo * (1 - u)));
    qhi = std::min(qhi, (rhi - u) / (rhi * (1 - u)));
  }
  if (qlo > qhi) { ans.reason = "interior dose constraints have disjoint shape ranges"; return ans; }
  qlo = std::min(1., std::max(0., qlo));
  qhi = std::min(1., std::max(0., qhi));
  ans.delta = interval(dlo, dhi); ans.q = interval(qlo, qhi);
  ans.target = interval(dlo * area(qlo), dhi * area(qhi)); ans.valid = true;
  return ans;
}
struct Slice { Interval delta, target; };
Slice slice(const Region& s, double q) {
  const arma::vec h = shape(s.x, q, s.D);
  const arma::vec hc = h - dot(h, s.w) / total(s.w);
  const arma::vec yc = s.means - dot(s.means, s.w) / total(s.w);
  const double shh = total(arma::vec(s.w % hc % hc));
  const double dhat = total(arma::vec(s.w % hc % yc)) / shh;
  const arma::vec residual = yc - dhat * hc;
  const double remaining = s.r2 - total(arma::vec(s.w % arma::square(residual)));
  if (remaining < 0) return {};
  const double half = std::sqrt(remaining / shh);
  if (dhat + half < 0) return {};
  Interval d = interval(std::max(0., dhat - half), dhat + half);
  return {d, interval(d.lo * area(q), d.hi * area(q))};
}
Interval ball(const arma::vec& direction, const arma::vec& target,
              double radius, Interval limits) {
  const double den = sq(direction), c = dot(direction, target) / den;
  const double slack = radius * radius - sq(arma::vec(target - c * direction));
  const double cushion = 64 * eps * scale({radius * radius, sq(target)});
  if (slack < -cushion) return {};
  const double half = std::sqrt(std::max(0., slack + cushion) / den);
  return interval(std::max(limits.lo, c - half), std::min(limits.hi, c + half));
}
struct Box {
  double lo = 0, hi = 0, lower = 0, upper = 0, dlo = 0, dhi = 0;
  bool valid = false;
};
Box legacy_box(const Region& s, const Outer& o, double lo, double hi) {
  const double middle = lo + (hi - lo) / 2;
  const arma::vec h = shape(s.x, middle, s.D);
  const arma::vec hlo = shape(s.x, lo, s.D), hhi = shape(s.x, hi, s.D);
  arma::vec displacement(h.n_elem);
  for (arma::uword j = 0; j < h.n_elem; ++j)
    displacement[j] = std::max(h[j] - hlo[j], hhi[j] - h[j]);
  const double epsilon = std::sqrt(total(arma::vec(s.w % arma::square(displacement))));
  const Interval d = ball(s.center(h), s.y, s.radius + o.delta.hi * epsilon, o.delta);
  if (!d.valid) return {};
  const double rawlo = d.lo * area(lo), rawhi = d.hi * area(hi);
  const double cushion = 128 * eps * scale({rawlo, rawhi, o.target.lo, o.target.hi});
  const double lower = std::max(o.target.lo, rawlo - cushion);
  const double upper = std::min(o.target.hi, rawhi + cushion);
  return {lo, hi, lower, upper, d.lo, d.hi, lower <= upper};
}
Box refined_box(const Region& s, const Outer& o, double lo, double hi,
                const Box* parent) {
  Box legacy = legacy_box(s, o, lo, hi);
  if (!legacy.valid) return {};
  const double middle = lo + (hi - lo) / 2, hw = (hi - lo) / 2;
  const arma::vec h = shape(s.x, middle, s.D), H = s.center(h);
  const Derivatives d = derivatives(s.x, middle, s.D), ed = derivatives(s.x, hi, s.D);
  const arma::vec Hfirst = s.center(d.first);
  const arma::vec difflo = shape(s.x, lo, s.D) - h, diffhi = shape(s.x, hi, s.D) - h;
  const double endpoint_epsilon = std::max(
    std::sqrt(total(arma::vec(s.w % arma::square(difflo)))),
    std::sqrt(total(arma::vec(s.w % arma::square(diffhi)))));
  const double second_norm = std::sqrt(total(arma::vec(s.w % arma::square(ed.second))));
  const double te = hw * norm(Hfirst) + 0.5 * hw * hw * second_norm;
  const double epsilon = std::min(endpoint_epsilon, te);
  const double A = sq(H), b = dot(s.y, H), ys = sq(s.y);
  const double slack = 4096 * eps * scale({ys, s.r2,
    A * o.delta.hi * o.delta.hi, 2 * std::abs(b) * o.delta.hi});
  const Interval self = quadratic(arma::vec({A - epsilon * epsilon,
    -2 * (b + s.radius * epsilon), ys - s.r2}), o.delta, slack);
  if (!self.valid) return {};
  const double first_norm = std::sqrt(total(arma::vec(s.w % arma::square(ed.first))));
  const double Hb = std::sqrt(A) + epsilon;
  const double Afirst = 2 * dot(H, Hfirst);
  const double Asecond = 2 * first_norm * first_norm + 2 * Hb * second_norm;
  const double bfirst = dot(s.y, Hfirst);
  const double bsecond = total(arma::vec(arma::abs(s.y) % s.rw % ed.second));
  const double Alower = std::max(0., A - std::abs(Afirst) * hw - 0.5 * Asecond * hw * hw);
  const double bupper = b + std::abs(bfirst) * hw + 0.5 * bsecond * hw * hw;
  const Interval coefficient = quadratic(arma::vec({Alower, -2 * bupper, ys - s.r2}), o.delta, slack);
  if (!coefficient.valid) return {};
  double dlo = std::max({o.delta.lo, legacy.dlo, self.lo, coefficient.lo});
  double dhi = std::min({o.delta.hi, legacy.dhi, self.hi, coefficient.hi});
  if (parent) { dlo = std::max(dlo, parent->dlo); dhi = std::min(dhi, parent->dhi); }
  if (dlo > dhi) return {};
  const double rawlo = dlo * area(lo), rawhi = dhi * area(hi);
  const double cushion = 512 * eps * scale({rawlo, rawhi, o.target.lo, o.target.hi});
  double lower = std::max({o.target.lo, legacy.lower, rawlo - cushion});
  double upper = std::min({o.target.hi, legacy.upper, rawhi + cushion});
  if (parent) { lower = std::max(lower, parent->lower); upper = std::min(upper, parent->upper); }
  return {lo, hi, lower, upper, dlo, dhi, lower <= upper};
}
using Row = std::vector<double>;
Rcpp::DataFrame frame(const std::vector<Row>& rows, const std::vector<std::string>& names) {
  if (rows.empty()) return Rcpp::DataFrame::create();
  Rcpp::List columns(names.size());
  for (std::size_t j = 0; j < names.size(); ++j) {
    Rcpp::NumericVector v(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) v[i] = rows[i][j];
    columns[j] = v;
  }
  columns.attr("names") = Rcpp::wrap(names);
  columns.attr("class") = "data.frame";
  columns.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -static_cast<int>(rows.size()));
  return Rcpp::DataFrame(columns);
}
const std::vector<std::string> trace_names = {"splits", "active_boxes", "lower_outer",
  "lower_inner", "upper_inner", "upper_outer", "lower_gap", "upper_gap", "evaluations"};
const std::vector<std::string> witness_names = {"q", "delta_lower", "delta_upper",
  "area_lower", "area_upper", "residual_lower", "residual_upper"};
struct Result {
  bool empty = false, feasible = false;
  Interval target;
  double bestlo = inf, besthi = -inf, gaplo = inf, gaphi = inf;
  int splits = 0, evaluations = 0;
  std::string termination = "budget_exhausted";
  std::vector<Row> witnesses, trace;
};
Result search(const Region& s, const Outer& o, double tolerance, int budget, bool refined) {
  Result r;
  if (!o.valid) { r.empty = true; r.termination = "empty_original_outer_envelope"; return r; }
  std::map<double, bool> cache;
  auto sample = [&](double q) {
    auto it = cache.find(q);
    if (it != cache.end()) return it->second;
    ++r.evaluations;
    Slice value = slice(s, q);
    if (!value.delta.valid) { cache[q] = false; return false; }
    const arma::vec H = s.center(shape(s.x, q, s.D));
    const double rl = sq(arma::vec(s.y - value.delta.lo * H));
    const double rh = sq(arma::vec(s.y - value.delta.hi * H));
    const double cushion = (refined ? 512 : 128) * eps * scale({s.r2, sq(s.y)});
    if (rl > s.r2 + cushion || rh > s.r2 + cushion)
      Rcpp::stop("fixed-shape witness failed its residual check");
    r.bestlo = std::min(r.bestlo, value.target.lo);
    r.besthi = std::max(r.besthi, value.target.hi);
    r.witnesses.push_back({q, value.delta.lo, value.delta.hi,
      value.target.lo, value.target.hi, rl, rh});
    cache[q] = true; return true;
  };
  auto enclose = [&](double lo, double hi, const Box* parent) {
    const double middle = lo + (hi - lo) / 2;
    bool feasible = sample(lo);
    if (middle != lo) feasible = sample(middle) || feasible;
    if (hi != lo && hi != middle) feasible = sample(hi) || feasible;
    Box b = refined ? refined_box(s, o, lo, hi, parent) : legacy_box(s, o, lo, hi);
    if (!refined && parent && b.valid) {
      b.lower = std::max(b.lower, parent->lower); b.upper = std::min(b.upper, parent->upper);
      b.valid = b.lower <= b.upper;
    }
    if (!b.valid && feasible) Rcpp::stop("shape box excluded a feasible fixed-q slice");
    return b;
  };
  const Box first = enclose(o.q.lo, o.q.hi, nullptr);
  std::vector<Box> boxes;
  if (first.valid) boxes.push_back(first);
  while (true) {
    Rcpp::checkUserInterrupt();
    r.feasible = std::isfinite(r.bestlo) && std::isfinite(r.besthi);
    if (r.feasible) boxes.erase(std::remove_if(boxes.begin(), boxes.end(),
      [&](const Box& b) { return b.lower >= r.bestlo && b.upper <= r.besthi; }), boxes.end());
    if (boxes.empty() && !r.feasible) {
      r.empty = true; r.termination = "all_shape_intervals_excluded"; return r;
    }
    double lower = r.bestlo, upper = r.besthi;
    for (const Box& b : boxes) { lower = std::min(lower, b.lower); upper = std::max(upper, b.upper); }
    r.target = interval(lower, upper);
    r.gaplo = r.feasible ? std::max(0., r.bestlo - lower) : inf;
    r.gaphi = r.feasible ? std::max(0., upper - r.besthi) : inf;
    r.trace.push_back({static_cast<double>(r.splits), static_cast<double>(boxes.size()),
      lower, r.bestlo, r.besthi, upper, r.gaplo, r.gaphi, static_cast<double>(r.evaluations)});
    if (std::max(r.gaplo, r.gaphi) <= tolerance) { r.termination = "endpoint_tolerance"; break; }
    if (r.splits >= budget) break;
    std::size_t index = 0;
    for (std::size_t j = 1; j < boxes.size(); ++j)
      if ((r.gaplo >= r.gaphi && boxes[j].lower < boxes[index].lower) ||
          (r.gaplo < r.gaphi && boxes[j].upper > boxes[index].upper)) index = j;
    Box b = boxes[index];
    const double middle = b.lo + (b.hi - b.lo) / 2;
    if (middle == b.lo || middle == b.hi) { r.termination = "floating_point_resolution"; break; }
    const Box left = enclose(b.lo, middle, &b), right = enclose(middle, b.hi, &b);
    boxes.erase(boxes.begin() + index);
    if (left.valid) boxes.push_back(left);
    if (right.valid) boxes.push_back(right);
    ++r.splits;
  }
  r.target.lo = std::max(o.target.lo, r.target.lo);
  r.target.hi = std::min(o.target.hi, r.target.hi);
  return r;
}
Rcpp::List output(Result r, double tolerance, int reference_budget) {
  const Rcpp::NumericVector nas = Rcpp::NumericVector::create(NA_REAL, NA_REAL);
  const Rcpp::NumericVector zeros = Rcpp::NumericVector::create(0., 0.);
  Rcpp::List ans;
  if (r.empty) {
    ans = Rcpp::List::create(Rcpp::Named("status") = "empty_set",
      Rcpp::Named("interval") = nas, Rcpp::Named("width") = 0.,
      Rcpp::Named("lower_endpoint") = nas, Rcpp::Named("upper_endpoint") = nas,
      Rcpp::Named("width_bounds") = zeros, Rcpp::Named("endpoint_gaps") = zeros,
      Rcpp::Named("converged") = true);
  } else {
    ans = Rcpp::List::create(Rcpp::Named("status") = r.feasible ? "nonempty" : "unresolved_feasibility",
      Rcpp::Named("interval") = Rcpp::NumericVector::create(r.target.lo, r.target.hi),
      Rcpp::Named("width") = r.target.hi - r.target.lo,
      Rcpp::Named("lower_endpoint") = r.feasible ? Rcpp::NumericVector::create(r.target.lo, r.bestlo) : nas,
      Rcpp::Named("upper_endpoint") = r.feasible ? Rcpp::NumericVector::create(r.besthi, r.target.hi) : nas,
      Rcpp::Named("width_bounds") = Rcpp::NumericVector::create(
        r.feasible ? std::max(0., r.besthi - r.bestlo) : 0., r.target.hi - r.target.lo),
      Rcpp::Named("endpoint_gaps") = Rcpp::NumericVector::create(r.gaplo, r.gaphi),
      Rcpp::Named("converged") = r.feasible && std::max(r.gaplo, r.gaphi) <= tolerance);
  }
  ans["termination"] = r.termination; ans["splits"] = r.splits;
  ans["evaluations"] = r.evaluations; ans["witness_count"] = static_cast<int>(r.witnesses.size());
  ans["witnesses"] = r.empty ? Rcpp::DataFrame::create() : frame(r.witnesses, witness_names);
  ans["trace"] = frame(r.trace, trace_names);
  ans["bound_in_exact_arithmetic"] = true;
  ans["outward_rounding_certified"] = false;
  ans["engine"] = "RcppArmadillo projection";
  ans["v1_reference_splits"] = reference_budget;
  return ans;
}
} // namespace

// [[Rcpp::export(rng = false)]]
arma::vec cpp_shape(const arma::vec& x, double q, double D) { return shape(x, q, D); }

// [[Rcpp::export(rng = false)]]
double cpp_area(double q) { return area(q); }

// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_derivatives(const arma::vec& x, double q, double D) {
  Derivatives d = derivatives(x, q, D);
  return Rcpp::List::create(Rcpp::Named("first") = d.first, Rcpp::Named("second") = d.second);
}

// [[Rcpp::export(rng = false)]]
SEXP cpp_quadratic(const arma::vec& co, const arma::vec& limits, double slack) {
  Interval ans = quadratic(co, interval(limits[0], limits[1]), slack);
  if (!ans.valid) return R_NilValue;
  return Rcpp::NumericVector::create(ans.lo, ans.hi);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_fixed_slice(Rcpp::List region, double q) {
  Slice ans = slice(Region(region), q);
  if (!ans.delta.valid) return Rcpp::List::create(Rcpp::Named("status") = "empty",
    Rcpp::Named("interval") = Rcpp::NumericVector::create(NA_REAL, NA_REAL));
  return Rcpp::List::create(Rcpp::Named("status") = "fixed_shape_projection",
    Rcpp::Named("delta") = Rcpp::NumericVector::create(ans.delta.lo, ans.delta.hi),
    Rcpp::Named("interval") = Rcpp::NumericVector::create(ans.target.lo, ans.target.hi),
    Rcpp::Named("q") = q);
}

// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_outer(Rcpp::List region) {
  Region s(region); Outer ans = outer(s);
  if (!ans.valid) return Rcpp::List::create(Rcpp::Named("status") = "empty_outer_set",
    Rcpp::Named("interval") = Rcpp::NumericVector::create(NA_REAL, NA_REAL),
    Rcpp::Named("reason") = ans.reason, Rcpp::Named("alpha") = s.alpha);
  return Rcpp::List::create(Rcpp::Named("status") = "analytic_outer_interval",
    Rcpp::Named("alpha") = s.alpha,
    Rcpp::Named("interval") = Rcpp::NumericVector::create(ans.target.lo, ans.target.hi),
    Rcpp::Named("delta") = Rcpp::NumericVector::create(ans.delta.lo, ans.delta.hi),
    Rcpp::Named("q") = Rcpp::NumericVector::create(ans.q.lo, ans.q.hi),
    Rcpp::Named("construction") = "analytic outer envelope of F-region projection; no optimizer",
    Rcpp::Named("note") = "nonempty envelope does not certify nonempty exact model intersection");
}

// [[Rcpp::export(rng = false)]]
Rcpp::List cpp_project(Rcpp::List region, double tolerance, int max_splits) {
  Region s(region); Outer o = outer(s);
  Result r = search(s, o, tolerance, max_splits, true);
  const int rb = std::min(256, max_splits);
  if (!r.empty) {
    Result reference = search(s, o, tolerance, rb, false);
    if (reference.empty && r.feasible) Rcpp::stop("v1 reference declared empty despite a v2 witness");
    if (!reference.empty) {
      r.target.lo = std::max(r.target.lo, reference.target.lo);
      r.target.hi = std::min(r.target.hi, reference.target.hi);
    }
    const double cushion = 1024 * eps * scale({r.target.lo, r.target.hi,
      r.feasible ? r.bestlo : 0, r.feasible ? r.besthi : 0});
    if (r.feasible && (r.target.lo > r.bestlo + cushion || r.target.hi < r.besthi - cushion))
      Rcpp::stop("outer projection bounds excluded a feasible witness");
    r.gaplo = r.feasible ? std::max(0., r.bestlo - r.target.lo) : inf;
    r.gaphi = r.feasible ? std::max(0., r.target.hi - r.besthi) : inf;
    if (!r.trace.empty()) {
      Row& t = r.trace.back(); t[2] = r.target.lo; t[5] = r.target.hi;
      t[6] = r.gaplo; t[7] = r.gaphi;
    }
  }
  return output(r, tolerance, rb);
}
