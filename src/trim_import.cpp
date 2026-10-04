#include "tessellation/private/repair.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

#include "tessellation/private/bezier.h"
#include "tessellation/private/exact.h"
#include "tessellation/tessellation.h"
#include "tessellation/private/internal.h"

namespace iko::internal {
static bool BezierKnots(const std::vector<double>& k, int d,
                        std::array<double, 2> domain) {
  return k.size() == static_cast<std::size_t>(2 * (d + 1)) &&
         std::all_of(k.begin(), k.begin() + d + 1,
                     [&](double knot) { return knot == k.front(); }) &&
         std::all_of(k.begin() + d + 1, k.end(),
                     [&](double knot) { return knot == k.back(); }) &&
         AlignedDomain(k, d, domain).has_value();
}
void ValidateProvidedTrim(const Face& face) {
  const auto& s = face.surface;
  if (s.p != 1 || s.q != 1 || s.closed_axes[0] || s.closed_axes[1] ||
      !BezierKnots(s.ku, s.p, {s.domain[0], s.domain[1]}) ||
      !BezierKnots(s.kv, s.q, {s.domain[2], s.domain[3]}) ||
      !std::all_of(s.weights.begin(), s.weights.end(),
                   [](double weight) { return weight == 1.; }) ||
      !face.outer || face.outer->size() < 3 || !face.inner.empty()) {
    throw Unsupported(
        "Supplied UV/XYZ wire requires an open bilinear support, one closed "
        "outer wire and no holes");
  }
  const auto& wire = *face.outer;
  for (std::size_t curve_index = 0; curve_index < wire.size(); ++curve_index) {
    const auto& c = wire[curve_index];
    if ((c.degree != 1 && c.degree != 2) ||
        !BezierKnots(c.knots, c.degree, c.domain) ||
        !std::all_of(c.control.begin(), c.control.end(),
                     [](auto control_point) {
                       return Finite(control_point) && control_point[2] == 0.;
                     }) ||
        !std::all_of(c.weights.begin(), c.weights.end(), [](double weight) {
          return std::isfinite(weight) && weight > 0.;
        })) {
      throw Unsupported(
          "Supplied UV/XYZ wire requires finite planar linear/quadratic "
          "Bezier UV curves with positive weights");
    }
    auto delta =
        Subtract(c.control.back(),
                 wire[(curve_index + 1) % wire.size()].control.front());
    if (std::max({std::abs(delta[0]), std::abs(delta[1]), std::abs(delta[2])}) >
        1e-12) {
      throw Unsupported("Open supplied UV/XYZ wire healing is not reproduced");
    }
  }
}
struct Converted {
  Curve source;
  std::vector<Vec3> h;
  std::vector<double> w;
  bool rational = false;
  int identity = -1;
};
static Converted ConvertTrimCurve(const Curve& c, int id) {
  Converted out;
  out.source = c;
  out.identity = id;
  bool constant =
      std::all_of(c.weights.begin(), c.weights.end(),
                  [&](double weight) { return weight == c.weights[0]; });
  out.w = c.weights;
  if (constant) {
    std::fill(out.w.begin(), out.w.end(), 1.);
  }
  out.rational = !std::all_of(out.w.begin(), out.w.end(),
                              [](double weight) { return weight == 1.; });
  for (std::size_t control_point_index = 0;
       control_point_index < c.control.size(); ++control_point_index) {
    out.h.push_back(
        Multiply(c.control[control_point_index], out.w[control_point_index]));
  }
  return out;
}
void RepairTrimWire(Face& face, std::vector<std::string>& trace) {
  auto& s = face.surface;
  if (!BezierKnots(s.ku, s.p, {s.domain[0], s.domain[1]}) ||
      !BezierKnots(s.kv, s.q, {s.domain[2], s.domain[3]})) {
    throw Unsupported(
        "IGES compatibility currently requires a full-domain Bezier support "
        "surface");
  }
  for (int axis = 0; axis < 2; ++axis) {
    auto domain = AlignedDomain(axis ? s.kv : s.ku, axis ? s.q : s.p,
                                {s.domain[2 * axis], s.domain[2 * axis + 1]})
                      .value();
    if (s.domain[2 * axis] != domain[0] ||
        s.domain[2 * axis + 1] != domain[1]) {
      trace.push_back("parameter_domain: canonicalize_full_domain support");
    }
    s.domain[2 * axis] = domain[0];
    s.domain[2 * axis + 1] = domain[1];
  }
  if (!face.outer || !face.inner.empty() || face.outer->size() < 3) {
    throw Unsupported(
        "IGES compatibility currently requires one closed outer wire and no "
        "holes");
  }
  for (double weight : s.weights) {
    if (!std::isfinite(weight) || weight < .001 || weight > 100.) {
      throw Unsupported(
          "IGES support surface weight conversion is outside the verified "
          "range");
    }
  }
  auto& source = *face.outer;
  std::vector<Converted> curves;
  for (std::size_t source_curve_index = 0; source_curve_index < source.size();
       ++source_curve_index) {
    auto& c = source[source_curve_index];
    if (c.weights.empty() || c.control.size() != c.weights.size() ||
        !BezierKnots(c.knots, c.degree, c.domain)) {
      throw Unsupported(
          "IGES compatibility requires finite planar full-domain Bezier "
          "trim curves");
    }
    for (auto control_point : c.control) {
      if (!Finite(control_point) || control_point[2] != 0) {
        throw Unsupported(
            "IGES compatibility requires finite planar full-domain Bezier "
            "trim curves");
      }
    }
    for (double weight : c.weights) {
      if (!std::isfinite(weight)) {
        throw Unsupported(
            "IGES compatibility requires finite planar full-domain Bezier "
            "trim curves");
      }
    }
    auto mm = std::minmax_element(c.weights.begin(), c.weights.end());
    double extent = *mm.second - *mm.first;
    if (extent > 0 && extent < 1e-10) {
      throw Unsupported(
          "Near-constant OCCT rationality classification is not "
          "reproduced");
    }
    auto domain = AlignedDomain(c.knots, c.degree, c.domain).value();
    if (c.domain != domain) {
      trace.push_back("parameter_domain: canonicalize_full_domain curve " +
                      std::to_string(source_curve_index));
    }
    c.domain = domain;
    auto delta = Subtract(
        c.control.back(),
        source[(source_curve_index + 1) % source.size()].control.front());
    if (std::max({std::abs(delta[0]), std::abs(delta[1]), std::abs(delta[2])}) >
        1e-12) {
      throw Unsupported(
          "Open IGES wire may be healed by OCCT before conversion; not "
          "reproduced");
    }
    curves.push_back(ConvertTrimCurve(c, static_cast<int>(source_curve_index)));
  }
  std::vector<Vec2> points;
  for (const auto& curve : source) {
    for (std::size_t control_point_index = 0;
         control_point_index + 1 < curve.control.size();
         ++control_point_index) {
      points.push_back({curve.control[control_point_index][0],
                        curve.control[control_point_index][1]});
    }
  }
  Vec2 lo = points[0];
  Vec2 hi = lo;
  for (auto parameter_point : points) {
    for (int axis = 0; axis < 2; ++axis) {
      lo[axis] = std::min(lo[axis], parameter_point[axis]);
      hi[axis] = std::max(hi[axis], parameter_point[axis]);
    }
  }
  double span = std::max({1., hi[0] - lo[0], hi[1] - lo[1]});
  double guard = 32 * kDoubleEpsilon * span * span;
  bool positive = false;
  bool negative = false;
  for (std::size_t polygon_vertex_index = 0;
       polygon_vertex_index < points.size(); ++polygon_vertex_index) {
    auto a = Subtract(points[(polygon_vertex_index + 1) % points.size()],
                      points[polygon_vertex_index]);
    auto b = Subtract(points[(polygon_vertex_index + 2) % points.size()],
                      points[(polygon_vertex_index + 1) % points.size()]);
    double cross = a[0] * b[1] - a[1] * b[0];
    positive |= cross > guard;
    negative |= cross < -guard;
  }
  if (positive == negative) {
    throw Unsupported(
        "IGES wire orientation currently requires a nondegenerate convex "
        "control polygon");
  }
  constexpr double rejection = 1e-7 / 100.;
  bool rejected = false;
  for (const auto& curve : source) {
    for (double weight : curve.weights) {
      rejected |= weight < rejection;
    }
  }
  if (rejected) {
    auto result = ConnectRejected(face, rejection, trace);
    if (result.wire.empty()) {
      face.outer.reset();
      trace.push_back("occt_face_transfer: discarded_wire_natural_support");
      return;
    }
    negative = negative != result.clockwise;
    if (result.wire.size() == 1 && negative) {
      throw Unsupported(
          "Clockwise single-survivor healing requires the shared-edge "
          "bending "
          "state; midpoint closure is not equivalent");
    }
    curves.clear();
    for (const auto& wire_item : result.wire) {
      curves.push_back(ConvertTrimCurve(wire_item.curve, wire_item.source));
    }
  }
  if (negative) {
    for (auto& converted_curve : curves) {
      converted_curve.source = ReverseCurve(converted_curve.source);
      std::reverse(converted_curve.h.begin(), converted_curve.h.end());
      std::reverse(converted_curve.w.begin(), converted_curve.w.end());
    }
    std::reverse(curves.begin(), curves.end());
    trace.push_back("occt_wire_orientation: reverse_clockwise_outer");
  }
  std::size_t curve_index = 0;
  while (curves.size() > 1 && curve_index < curves.size()) {
    std::size_t next_curve_index = (curve_index + 1) % curves.size();
    auto& a = curves[curve_index];
    auto& b = curves[next_curve_index];
    auto delta = Subtract(a.h.back(), b.h.front());
    double wd = a.rational && b.rational ? a.w.back() - b.w.front() : 0.;
    double d2 = wd * wd + ((delta[1] * delta[1] + delta[0] * delta[0]) +
                           delta[2] * delta[2]);
    if (d2 > 1e-6) {
      Curve line;
      line.degree = 1;
      line.knots = {0, 0, 1, 1};
      line.domain = {0, 1};
      line.weights = {1, 1};
      auto start = a.h.back();
      auto end = b.h.front();
      for (auto& coordinate : start) {
        coordinate /= a.w.back();
      }
      for (auto& coordinate : end) {
        coordinate /= b.w.front();
      }
      line.control = {start, end};
      curves.insert(curves.begin() + next_curve_index,
                    ConvertTrimCurve(line, -1));
      ++curve_index;
      trace.push_back("import_wire: insert_connector");
    } else {
      auto mean = Multiply(Add(a.h.back(), b.h.front()), .5);
      a.h.back() = mean;
      b.h.front() = mean;
      if (a.rational && b.rational) {
        double w = (a.w.back() + b.w.front()) * .5;
        a.w.back() = w;
        b.w.front() = w;
      }
      trace.push_back("import_wire: average_shared_homogeneous_endpoint");
    }
    ++curve_index;
  }
  for (std::size_t curve_index = 0; curve_index < curves.size();
       ++curve_index) {
    auto& a = curves[curve_index];
    auto& b = curves[(curve_index + 1) % curves.size()];
    auto start = a.h.back();
    auto end = b.h.front();
    for (auto& coordinate : start) {
      coordinate /= a.w.back();
    }
    for (auto& coordinate : end) {
      coordinate /= b.w.front();
    }
    auto delta = Subtract(end, start);
    double gap = std::sqrt((delta[1] * delta[1] + delta[0] * delta[0]) +
                           delta[2] * delta[2]);
    if (gap >= 1e-6) {
      auto xyz1 = s.Evaluate({start[0], start[1]});
      auto xyz2 = s.Evaluate({end[0], end[1]});
      if (gap > .01 && Norm(Subtract(xyz2, xyz1)) < .001) {
        throw Unsupported(
            "Periodic/surface-space constructor connector remains "
            "unresolved");
      }
      auto midpoint = Multiply(Add(start, end), .5);
      a.h.back() = midpoint;
      b.h.front() = midpoint;
      trace.push_back("surface_constructor: midpoint_written");
    }
  }
  std::vector<Curve> result;
  for (auto& converted_curve : curves) {
    Curve c = converted_curve.source;
    c.weights.clear();
    c.control.clear();
    std::vector<double> w;
    for (double weight : converted_curve.w) {
      w.push_back(std::clamp(weight, .001, 100.));
    }
    auto mm = std::minmax_element(w.begin(), w.end());
    double extent = *mm.second - *mm.first;
    if (extent > 0 && extent < 1e-10) {
      throw Unsupported(
          "Near-constant constructor weight classification is not "
          "reproduced");
    }
    bool constant = extent == 0;
    for (std::size_t control_point_index = 0; control_point_index < w.size();
         ++control_point_index) {
      auto point = converted_curve.h[control_point_index];
      for (auto& coordinate : point) {
        coordinate /= w[control_point_index];
      }
      c.control.push_back(point);
      c.weights.push_back(constant ? 1. : w[control_point_index]);
    }
    if (w != converted_curve.w) {
      trace.push_back(
          "surface_constructor: clamp_weights_preserve_homogeneous");
    }
    result.push_back(std::move(c));
  }
  if (result.size() == 1) {
    const auto& c = result[0];
    if (c.degree != 3 || c.control.front() != c.control.back()) {
      throw Unsupported(
          "Single surviving IGES curve requires a nondegenerate closed "
          "cubic; "
          "native degenerate output remains unresolved");
    }
    Rat ax(c.control[1][0]);
    Rat ay(c.control[1][1]);
    Rat bx(c.control[2][0]);
    Rat by(c.control[2][1]);
    Rat ox(c.control[0][0]);
    Rat oy(c.control[0][1]);
    Rat area = (ax - ox) * (by - oy) - (ay - oy) * (bx - ox);
    if (exact::Abs(area) <= Rat(guard)) {
      throw Unsupported(
          "Single surviving cubic is degenerate or outside the verified "
          "orientation margin");
    }
    trace.push_back("surface_constructor: single_cubic_endpoint_closure");
  }
  face.outer = std::move(result);
  trace.push_back("surface_constructor: effective_wire");
}
}  // namespace iko::internal
