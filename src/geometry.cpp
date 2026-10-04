#include "tessellation/private/internal.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tessellation/tessellation.h"

namespace iko::internal {
std::string Number(double x) {
  std::ostringstream out;
  out << std::setprecision(17) << x;
  return out.str();
}
std::string BlockKey(Rect r) {
  // Parameter identity names use one representation for mathematical zero.
  // OCCT can supply -0 at reversed boundaries; it is the same knot boundary.
  for (auto& boundary_coordinate : r) {
    if (boundary_coordinate == 0.) {
      boundary_coordinate = 0.;
    }
  }
  return Number(r[0]) + "," + Number(r[1]) + "," + Number(r[2]) + "," +
         Number(r[3]);
}
std::vector<double> Basis(double x, const std::vector<double>& k, int d,
                          bool local, bool extrapolate) {
  int n = static_cast<int>(k.size()) - d - 1;
  if (d < 0 || n <= 0) {
    throw std::invalid_argument("Knot dimensions");
  }
  if (!extrapolate) {
    if (x < k[d] - 1e-12 || x > k[n] + 1e-12) {
      throw std::invalid_argument("Parameter outside knot domain");
    }
    x = std::clamp(x, k[d], k[n]);
  }
  if (local) {
    int span = std::clamp(
        static_cast<int>(std::upper_bound(k.begin(), k.end(), x) - k.begin()) -
            1,
        d, n - 1);
    std::vector<double> left(d + 1);
    std::vector<double> right(d + 1);
    std::vector<double> b(d + 1);
    std::vector<double> out(n);
    b[0] = 1;
    for (int recurrence_degree = 1; recurrence_degree <= d;
         ++recurrence_degree) {
      left[recurrence_degree] = x - k[span + 1 - recurrence_degree];
      right[recurrence_degree] = k[span + recurrence_degree] - x;
      double saved = 0;
      for (int basis_index = 0; basis_index < recurrence_degree;
           ++basis_index) {
        double den =
            right[basis_index + 1] + left[recurrence_degree - basis_index];
        double t = std::abs(den) >= std::numeric_limits<float>::min()
                       ? b[basis_index] / den
                       : 0;
        b[basis_index] = saved + right[basis_index + 1] * t;
        saved = left[recurrence_degree - basis_index] * t;
      }
      b[recurrence_degree] = saved;
    }
    for (int local_basis_index = 0; local_basis_index <= d;
         ++local_basis_index) {
      out[span - d + local_basis_index] = b[local_basis_index];
    }
    return out;
  }
  std::vector<double> b(k.size() - 1);
  for (std::size_t knot_span_index = 0; knot_span_index < b.size();
       ++knot_span_index) {
    b[knot_span_index] =
        (x >= k[knot_span_index] && x < k[knot_span_index + 1]) ? 1. : 0.;
  }
  for (int basis_degree = 1; basis_degree <= d; ++basis_degree) {
    std::vector<double> next(k.size() - basis_degree - 1);
    for (std::size_t basis_index = 0; basis_index < next.size();
         ++basis_index) {
      double den = k[basis_index + basis_degree] - k[basis_index];
      if (den) {
        next[basis_index] += (x - k[basis_index]) / den * b[basis_index];
      }
      den = k[basis_index + basis_degree + 1] - k[basis_index + 1];
      if (den) {
        next[basis_index] +=
            (k[basis_index + basis_degree + 1] - x) / den * b[basis_index + 1];
      }
    }
    b = std::move(next);
  }
  if (k[n] == k.back() && x == k[n]) {
    std::fill(b.begin(), b.end(), 0);
    b.back() = 1;
  }
  return b;
}
// numpy's contiguous reduction uses eight independent accumulators and
// pairwise blocks; Surface.point uses that reduction for its denominator.
static double PairSum(const double* x, std::size_t n) {
  if (n < 8) {
    double r = -0.;
    for (std::size_t element_index = 0; element_index < n; ++element_index) {
      r += x[element_index];
    }
    return r;
  }
  if (n <= 128) {
    double r[8];
    for (int accumulator_index = 0; accumulator_index < 8;
         ++accumulator_index) {
      r[accumulator_index] = x[accumulator_index];
    }
    std::size_t input_index = 8;
    for (; input_index < n - n % 8; input_index += 8) {
      for (int accumulator_index = 0; accumulator_index < 8;
           ++accumulator_index) {
        r[accumulator_index] += x[input_index + accumulator_index];
      }
    }
    double total =
        ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
    for (; input_index < n; ++input_index) {
      total += x[input_index];
    }
    return total;
  }
  std::size_t half = n / 2;
  half -= half % 8;
  return PairSum(x, half) + PairSum(x + half, n - half);
}
void Surface::Validate() const {
  if (p < 1 || q < 1 || nu <= p || nv <= q ||
      static_cast<int>(ku.size()) != nu + p + 1 ||
      static_cast<int>(kv.size()) != nv + q + 1 ||
      static_cast<int>(control.size()) != nu * nv ||
      weights.size() != control.size()) {
    throw std::invalid_argument("Control or knot dimensions");
  }
  for (const auto& axis_knots : {ku, kv}) {
    if (!std::is_sorted(axis_knots.begin(), axis_knots.end()) ||
        !std::all_of(axis_knots.begin(), axis_knots.end(),
                     [](double knot) { return std::isfinite(knot); })) {
      throw std::invalid_argument("Invalid knots");
    }
  }
  for (double weight : weights) {
    if (!(weight > 0)) {
      throw Unsupported("Nonpositive rational weights");
    }
  }
  for (double weight : weights) {
    if (!std::isfinite(weight)) {
      throw std::invalid_argument("Nonfinite weights");
    }
  }
  for (auto control_point : control) {
    if (!Finite(control_point)) {
      throw std::invalid_argument("Nonfinite control coordinates");
    }
  }
  if (!Finite(domain) || domain[0] >= domain[1] || domain[2] >= domain[3]) {
    throw std::invalid_argument("Invalid parameter domain");
  }
}
double Surface::BboxSpan() const {
  Vec3 lo = control.front();
  Vec3 hi = lo;
  for (auto control_point : control) {
    for (int coordinate_index = 0; coordinate_index < 3; ++coordinate_index) {
      lo[coordinate_index] =
          std::min(lo[coordinate_index], control_point[coordinate_index]);
      hi[coordinate_index] =
          std::max(hi[coordinate_index], control_point[coordinate_index]);
    }
  }
  return std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});
}
Vec3 Surface::Point(Vec2 uv) const {
  auto bu = Basis(uv[0], ku, p);
  auto bv = Basis(uv[1], kv, q);
  Vec3 h{};
  std::vector<double> w(nu * nv);
  for (int control_v_index = 0; control_v_index < nv; ++control_v_index) {
    for (int control_u_index = 0; control_u_index < nu; ++control_u_index) {
      int control_point_index = control_v_index * nu + control_u_index;
      w[control_point_index] = (bv[control_v_index] * bu[control_u_index]) *
                               weights[control_point_index];
      for (int coordinate_index = 0; coordinate_index < 3; ++coordinate_index) {
        h[coordinate_index] += w[control_point_index] *
                               control[control_point_index][coordinate_index];
      }
    }
  }
  return Multiply(h, 1. / PairSum(w.data(), w.size()));
}
Vec3 Surface::Evaluate(Vec2 uv) const {
  auto bu = Basis(uv[0], ku, p);
  auto bv = Basis(uv[1], kv, q);
  Vec3 h{};
  double w = 0;
  for (int control_v_index = 0; control_v_index < nv; ++control_v_index) {
    for (int control_u_index = 0; control_u_index < nu; ++control_u_index) {
      int control_point_index = control_v_index * nu + control_u_index;
      double b = (bu[control_u_index] * bv[control_v_index]) *
                 weights[control_point_index];
      w += b;
      for (int coordinate_index = 0; coordinate_index < 3; ++coordinate_index) {
        h[coordinate_index] +=
            b * control[control_point_index][coordinate_index];
      }
    }
  }
  for (auto& weighted_coordinate : h) {
    weighted_coordinate /= w;
  }
  return h;
}
using HomogeneousPoint = std::array<double, 4>;
static HomogeneousPoint Tensor(const std::vector<double>& bu,
                               const std::vector<double>& bv,
                               const std::vector<HomogeneousPoint>& cp, int nu,
                               int nv) {
  HomogeneousPoint value{};
  for (int control_v_index = 0; control_v_index < nv; ++control_v_index) {
    HomogeneousPoint row{};
    for (int control_u_index = 0; control_u_index < nu; ++control_u_index) {
      row = Add(row, Multiply(cp[control_v_index * nu + control_u_index],
                              bu[control_u_index]));
    }
    value = Add(value, Multiply(row, bv[control_v_index]));
  }
  return value;
}
Vec3 Surface::EvaluateControl(Vec2 uv, bool ext) const {
  auto bu = Basis(uv[0], ku, p, true, ext);
  auto bv = Basis(uv[1], kv, q, true, ext);
  std::vector<HomogeneousPoint> h(nu * nv);
  for (int control_point_index = 0; control_point_index < nu * nv;
       ++control_point_index) {
    double w = explicit_weights ? weights[control_point_index] : 1.;
    h[control_point_index] = {control[control_point_index][0] * w,
                              control[control_point_index][1] * w,
                              control[control_point_index][2] * w, w};
  }
  auto value = Tensor(bu, bv, h, nu, nv);
  Vec3 xyz{value[0], value[1], value[2]};
  return explicit_weights ? Multiply(xyz, 1. / value[3]) : xyz;
}
std::array<Vec3, 2> Surface::Derivatives(
    Vec2 uv, const std::array<std::vector<double>, 2>* bk, bool ext) const {
  const auto& bku = bk ? (*bk)[0] : ku;
  const auto& bkv = bk ? (*bk)[1] : kv;
  auto bu = Basis(uv[0], bku, p, true, ext);
  auto bv = Basis(uv[1], bkv, q, true, ext);
  std::vector<HomogeneousPoint> h(nu * nv);
  for (int control_point_index = 0; control_point_index < nu * nv;
       ++control_point_index) {
    h[control_point_index] = {
        control[control_point_index][0] * weights[control_point_index],
        control[control_point_index][1] * weights[control_point_index],
        control[control_point_index][2] * weights[control_point_index],
        weights[control_point_index]};
  }
  auto value = Tensor(bu, bv, h, nu, nv);
  double iw = 1. / value[3];
  Vec3 xyz{value[0] * iw, value[1] * iw, value[2] * iw};
  std::array<Vec3, 2> out;
  for (int axis = 0; axis < 2; ++axis) {
    int degree = axis ? q : p;
    const auto& k = axis ? kv : ku;
    int nx = nu - (axis == 0);
    int ny = nv - (axis == 1);
    std::vector<HomogeneousPoint> net(nx * ny);
    for (int derivative_v_index = 0; derivative_v_index < ny;
         ++derivative_v_index) {
      for (int derivative_u_index = 0; derivative_u_index < nx;
           ++derivative_u_index) {
        int index = axis ? derivative_v_index : derivative_u_index;
        double den = k[index + degree + 1] - k[index + 1];
        double factor = den ? static_cast<double>(degree) / den : 0;
        net[derivative_v_index * nx + derivative_u_index] =
            Multiply(Subtract(h[(derivative_v_index + (axis == 1)) * nu +
                                derivative_u_index + (axis == 0)],
                              h[derivative_v_index * nu + derivative_u_index]),
                     factor);
      }
    }
    auto du =
        axis == 0
            ? Basis(uv[0], std::vector<double>(bku.begin() + 1, bku.end() - 1),
                    p - 1, true, ext)
            : bu;
    auto dv =
        axis == 1
            ? Basis(uv[1], std::vector<double>(bkv.begin() + 1, bkv.end() - 1),
                    q - 1, true, ext)
            : bv;
    auto d = Tensor(du, dv, net, nx, ny);
    Vec3 a{d[0], d[1], d[2]};
    out[axis] =
        explicit_weights ? Multiply(Subtract(a, Multiply(xyz, d[3])), iw) : a;
  }
  return out;
}
Vec3 Curve::Evaluate(double t, bool observed) const {
  bool polynomial =
      prepared && std::all_of(weights.begin(), weights.end(),
                              [](double weight) { return weight == 1.; });
  observed |= prepared;
  auto b = Basis(t, knots, degree, observed);
  Vec3 h{};
  double w = 0;
  for (std::size_t control_point_index = 0;
       control_point_index < control.size(); ++control_point_index) {
    double c = b[control_point_index] * weights[control_point_index];
    w += c;
    for (int coordinate_index = 0; coordinate_index < 3; ++coordinate_index) {
      h[coordinate_index] +=
          observed ? b[control_point_index] *
                         (control[control_point_index][coordinate_index] *
                          weights[control_point_index])
                   : c * control[control_point_index][coordinate_index];
    }
  }
  if (polynomial) {
    return h;
  }
  if (observed) {
    if (std::abs(w) >= kFloatEpsilon) {
      h = Multiply(h, 1. / w);
    }
  } else {
    for (auto& weighted_coordinate : h) {
      weighted_coordinate /= w;
    }
  }
  return h;
}
std::vector<Rect> KnotBlocks(const Surface& s, bool midpoint, bool c0) {
  std::array<std::vector<double>, 2> axes;
  for (int axis = 0; axis < 2; ++axis) {
    auto& values = axes[axis];
    double lo = s.domain[2 * axis];
    double hi = s.domain[2 * axis + 1];
    int d = axis ? s.q : s.p;
    const auto& k = axis ? s.kv : s.ku;
    values = {lo, hi};
    for (double knot : k) {
      if (knot > lo && knot < hi &&
          (!c0 || std::count(k.begin(), k.end(), knot) >= d)) {
        values.push_back(knot);
      }
    }
    values = Unique(values);
    if (midpoint && d > 1 && values.size() == 2) {
      values.insert(values.begin() + 1, (lo + hi) * .5);
    }
  }
  std::vector<Rect> out;
  for (std::size_t v_span_index = 0; v_span_index + 1 < axes[1].size();
       ++v_span_index) {
    for (std::size_t u_span_index = 0; u_span_index + 1 < axes[0].size();
         ++u_span_index) {
      out.push_back({axes[0][u_span_index], axes[0][u_span_index + 1],
                     axes[1][v_span_index], axes[1][v_span_index + 1]});
    }
  }
  return out;
}
Vec3 AnalyticNormal(const Surface& s, Vec2 uv, bool ext) {
  auto d = s.Derivatives(uv, nullptr, ext);
  auto n = Cross(d[0], d[1]);
  double length = std::sqrt((n[0] * n[0] + n[1] * n[1]) + n[2] * n[2]);
  if (std::isfinite(length) && length > kFloatEpsilon) {
    return Multiply(n, 1. / length);
  }
  double nan = std::numeric_limits<double>::quiet_NaN();
  return {nan, nan, nan};
}
std::vector<Vec3> MeshNormals(const Surface& s, const std::vector<Vec2>& uv,
                              const Rect* rect) {
  auto samples = uv;
  bool ext = false;
  for (int axis = 0; axis < 2; ++axis) {
    double lo = s.domain[2 * axis];
    double hi = s.domain[2 * axis + 1];
    double tol =
        32 * kDoubleEpsilon * std::max({1., std::abs(lo), std::abs(hi)});
    int d = axis ? s.q : s.p;
    const auto& k = axis ? s.kv : s.ku;
    int spans = 0;
    for (int knot_span_index = d;
         knot_span_index < static_cast<int>(k.size()) - d - 1;
         ++knot_span_index) {
      spans += k[knot_span_index + 1] > k[knot_span_index];
    }
    bool extension = spans < 5 && lo == k[d] && hi == k[k.size() - d - 1] &&
                     (k.front() < lo || k.back() > hi);
    double allowed_lo = extension ? k.front() : lo;
    double allowed_hi = extension ? k.back() : hi;
    for (auto& parameter_point : samples) {
      if (parameter_point[axis] < allowed_lo - tol ||
          parameter_point[axis] > allowed_hi + tol) {
        throw std::invalid_argument("Normal query outside surface domain");
      }
      parameter_point[axis] =
          std::clamp(parameter_point[axis], allowed_lo, allowed_hi);
      ext |= parameter_point[axis] < lo || parameter_point[axis] > hi;
      if (rect) {
        for (int side_index = 0; side_index < 2; ++side_index) {
          double edge = (*rect)[2 * axis + side_index];
          if (lo < edge && edge < hi &&
              std::abs(parameter_point[axis] - edge) <= tol) {
            parameter_point[axis] += side_index ? -1e-6 : 1e-6;
          }
        }
      }
    }
  }
  std::vector<Vec3> out;
  out.reserve(samples.size());
  for (auto parameter_point : samples) {
    auto n = AnalyticNormal(s, parameter_point, ext);
    if (!Finite(n)) {
      Vec3 total{};
      for (auto neighbor_offset :
           {std::pair<int, int>{0, -1}, {0, 1}, {1, -1}, {1, 1}}) {
        auto q = parameter_point;
        int a = neighbor_offset.first;
        q[a] += neighbor_offset.second *
                (.01 * (s.domain[2 * a + 1] - s.domain[2 * a]));
        if (q[a] < s.domain[2 * a] || q[a] > s.domain[2 * a + 1]) {
          continue;
        }
        auto neighbor = AnalyticNormal(s, q, ext);
        if (Finite(neighbor)) {
          total = Add(total, neighbor);
        }
      }
      double sq =
          (total[0] * total[0] + total[1] * total[1]) + total[2] * total[2];
      if (std::isfinite(sq) && sq > std::numeric_limits<double>::min()) {
        n = Multiply(total, 1. / std::sqrt(sq));
      }
    }
    out.push_back(n);
  }
  return out;
}
std::optional<std::array<double, 2>> AlignedDomain(
    const std::vector<double>& k, int d, std::array<double, 2> domain) {
  if (d < 1 || k.size() < static_cast<std::size_t>(2 * (d + 1)) ||
      !Finite(domain) || !std::is_sorted(k.begin(), k.end())) {
    return {};
  }
  std::array<double, 2> active{k[d], k[k.size() - d - 1]};
  double shortest = std::min(active[1] - active[0], domain[1] - domain[0]);
  if (!(shortest > 0)) {
    return {};
  }
  for (std::size_t knot_index = 1; knot_index < k.size(); ++knot_index) {
    double step = k[knot_index] - k[knot_index - 1];
    if (!std::isfinite(step)) {
      return {};
    }
    if (step > 0) {
      shortest = std::min(shortest, step);
    }
  }
  double limit = shortest * (8 * kDoubleEpsilon);
  for (int endpoint_index = 0; endpoint_index < 2; ++endpoint_index) {
    if (domain[endpoint_index] == active[endpoint_index]) {
      continue;
    }
    if (std::abs(domain[endpoint_index] - active[endpoint_index]) > limit) {
      return {};
    }
    double x = domain[endpoint_index];
    for (int rounding_step = 0;
         rounding_step < 4 && x != active[endpoint_index]; ++rounding_step) {
      x = std::nextafter(x, active[endpoint_index]);
    }
    if (x != active[endpoint_index]) {
      return {};
    }
  }
  return active;
}
Sampling Sample(const std::function<Vec3(double)>& point, double start,
                double end, double chord, double minimum, double maximum,
                bool relative, int initial_segments,
                const std::function<Vec3(double)>& initial_point) {
  if (initial_segments < 1) {
    throw std::invalid_argument("Initial segment count must be positive");
  }
  static const double kCurveSamplingProbes[] = {
      .5210802, .4720207, .6111612, .351051, .7512508, .1921505, .927};
  Sampling result;
  int count = 0;
  std::function<void(double, double, Vec3, Vec3, int)> recurse =
      [&](double a, double b, Vec3 pa, Vec3 pb, int depth) {
        if (++count > 100000) {
          throw std::runtime_error(
              "Independent curve sampling work budget reached");
        }
        auto accept = [&] {
          result.parameters.push_back(a);
          result.points.push_back(pa);
        };
        if (depth > 10) {
          accept();
          return;
        }
        double l2 = Dot(Subtract(pb, pa), Subtract(pb, pa));
        if (l2 < 1e-16) {
          accept();
          return;
        }
        double scale = relative ? l2 : 1.;
        for (int probe_index = 0; probe_index < 7; ++probe_index) {
          double alpha = kCurveSamplingProbes[probe_index];
          double t = (1 - alpha) * a + alpha * b;
          auto p = point(t);
          auto delta =
              Subtract(Add(Multiply(pa, 1 - alpha), Multiply(pb, alpha)), p);
          double error2 = Dot(delta, delta);
          bool split = false;
          if (probe_index == 0 && !relative) {
            auto da = Subtract(pa, p);
            auto db = Subtract(pb, p);
            double half = .5 * (Dot(da, da) + Dot(db, db));
            if (minimum * minimum > half) {
              accept();
              return;
            }
            split = half > maximum * maximum;
          }
          if (split || error2 > chord * chord * scale) {
            recurse(a, t, pa, p, depth + 1);
            recurse(t, b, p, pb, depth + 1);
            return;
          }
        }
        accept();
      };
  auto original = initial_point ? initial_point : point;
  if (initial_segments == 1) {
    auto delta = Subtract(original(end), original(start));
    if (Dot(delta, delta) < 1e-8) {
      initial_segments = 2;
    }
  }
  auto pa = point(start);
  for (int segment_index = 0; segment_index < initial_segments;
       ++segment_index) {
    double a = start + segment_index * ((end - start) / initial_segments);
    double b =
        segment_index + 1 == initial_segments
            ? end
            : start + (segment_index + 1) * ((end - start) / initial_segments);
    auto pb = point(b);
    recurse(a, b, pa, pb, 0);
    pa = pb;
  }
  result.parameters.push_back(end);
  result.points.push_back(pa);
  return result;
}
}  // namespace iko::internal
