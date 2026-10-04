#include "tessellation/private/internal.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "tessellation/tessellation.h"

namespace iko::internal {
namespace {
// A representation check, separate from geometric mesh comparison tolerance.
// OCCT's pcurve confusion permits tiny endpoint roundoff on natural boundaries.
constexpr double kParameterConfusion = 1e-9;

bool NaturalRectangle(const Face& f, bool allow_closed = false) {
  if (!f.outer || (f.support_type == "Geom_Plane" && f.outer->size() != 4) ||
      !f.inner.empty() ||
      (!allow_closed &&
       (f.geometrically_closed_axes[0] || f.geometrically_closed_axes[1]))) {
    return false;
  }
  const auto& d = f.surface.domain;
  std::array<std::vector<std::array<double, 2>>, 4> ranges;
  for (const auto& curve : *f.outer) {
    int side = -1;
    for (int boundary_side_index = 0; boundary_side_index < 4;
         ++boundary_side_index) {
      bool fixed = std::all_of(
          curve.control.begin(), curve.control.end(), [&](auto control_point) {
            return std::abs(control_point[boundary_side_index / 2] -
                            d[boundary_side_index]) <= kParameterConfusion;
          });
      if (fixed) {
        side = boundary_side_index;
        break;
      }
    }
    if (side < 0) {
      return false;
    }
    int axis = 1 - side / 2;
    double lo = curve.control.front()[axis];
    double hi = lo;
    for (auto control_point : curve.control) {
      lo = std::min(lo, control_point[axis]);
      hi = std::max(hi, control_point[axis]);
    }
    if (lo < d[2 * axis] - kParameterConfusion ||
        hi > d[2 * axis + 1] + kParameterConfusion) {
      return false;
    }
    ranges[side].push_back({lo, hi});
  }
  for (int boundary_side_index = 0; boundary_side_index < 4;
       ++boundary_side_index) {
    auto& spans = ranges[boundary_side_index];
    if (spans.empty()) {
      return false;
    }
    std::sort(spans.begin(), spans.end());
    int axis = 1 - boundary_side_index / 2;
    double covered = d[2 * axis];
    for (auto span : spans) {
      if (span[0] > covered + kParameterConfusion) {
        return false;
      }
      covered = std::max(covered, span[1]);
    }
    if (std::abs(covered - d[2 * axis + 1]) > kParameterConfusion) {
      return false;
    }
  }
  return true;
}

void ElevateEditedLine(Curve& c) {
  if (c.degree != 1 || c.control.size() != 2) {
    return;
  }
  // Editing an analytic line before conversion retains its old midpoint.
  // A supplied B-spline is already converted and keeps its degree.
  if (!c.analytic_line) {
    return;
  }
  const auto midpoint = Multiply(Add(c.control[0], c.control[1]), .5);
  c.control.insert(c.control.begin() + 1, midpoint);
  c.weights = {1., 1., 1.};
  c.degree = 2;
  c.knots = {c.domain[0], c.domain[0], c.domain[0],
             c.domain[1], c.domain[1], c.domain[1]};
}

void ConvertEditedCurve(Curve& c) {
  if (!c.endpoint_edit_conversion) {
    ElevateEditedLine(c);
    return;
  }
  auto first = c.control.front();
  auto last = c.control.back();
  auto converted = c.endpoint_edit_conversion();
  converted.control.front() = first;
  converted.control.back() = last;
  converted.weights.front() = c.weights.front();
  converted.weights.back() = c.weights.back();
  c = std::move(converted);
}

}  // namespace

void JoinOcctConversionEndpoints(std::vector<Curve>& curves,
                                 std::vector<int>& edge_ids,
                                 std::vector<std::string>& trace,
                                 const std::string& identity) {
  for (auto& curve : curves) {
    if (curve.weights.empty()) {
      throw Unsupported("OCCT adaptation has an empty curve");
    }
    bool constant = std::all_of(
        curve.weights.begin(), curve.weights.end(),
        [&](double weight) { return weight == curve.weights.front(); });
    if (constant) {
      std::fill(curve.weights.begin(), curve.weights.end(), 1.);
    }
  }
  std::size_t curve_index = 0;
  while (!curves.empty() && curve_index < curves.size()) {
    std::size_t next_curve_index = (curve_index + 1) % curves.size();
    auto& a = curves[curve_index];
    auto& b = curves[next_curve_index];
    auto ah = Multiply(a.control.back(), a.weights.back());
    auto bh = Multiply(b.control.front(), b.weights.front());
    auto delta = Subtract(ah, bh);
    bool ar = !std::all_of(a.weights.begin(), a.weights.end(),
                           [](double weight) { return weight == 1.; });
    bool br = !std::all_of(b.weights.begin(), b.weights.end(),
                           [](double weight) { return weight == 1.; });
    double dw = ar && br ? a.weights.back() - b.weights.front() : 0.;
    double distance = dw * dw + ((delta[1] * delta[1] + delta[0] * delta[0]) +
                                 delta[2] * delta[2]);
    if (distance > 1e-6) {
      // A lone curve is closed by the wireframe/surface constructor,
      // rather than by inserting a separate conversion connector.
      // Retain sub-threshold averaging for closed analytic curves.
      if (curves.size() == 1) {
        break;
      }
      Curve line;
      line.control = {a.control.back(), b.control.front()};
      line.weights = {1., 1.};
      line.knots = {0., 0., 1., 1.};
      line.domain = {0., 1.};
      curves.insert(curves.begin() + next_curve_index, std::move(line));
      edge_ids.insert(edge_ids.begin() + next_curve_index,
                      0);  // conversion connector, not an OCCT topological edge
      trace.push_back(identity + ": compatible conversion connector");
      ++curve_index;
    } else if (distance != 0.) {
      auto mean = Multiply(Add(ah, bh), .5);
      double aw = a.weights.back();
      double bw = b.weights.front();
      if (ar && br) {
        aw = bw = (aw + bw) * .5;
      }
      // Analytic endpoint edits use their conversion threshold; already
      // converted two-pole splines retain the pcurve confusion threshold.
      if (Norm(Subtract(Multiply(mean, 1. / aw), a.control.back())) >
          ((a.analytic_line || a.endpoint_edit_conversion)
               ? 1e-8
               : kParameterConfusion)) {
        ConvertEditedCurve(a);
      }
      if (Norm(Subtract(Multiply(mean, 1. / bw), b.control.front())) >
          ((b.analytic_line || b.endpoint_edit_conversion)
               ? 1e-8
               : kParameterConfusion)) {
        ConvertEditedCurve(b);
      }
      a.weights.back() = aw;
      b.weights.front() = bw;
      // Preserve the native homogeneous divide operation, not reciprocal
      // multiply.
      auto ap = mean;
      auto bp = mean;
      for (int coordinate_index = 0; coordinate_index < 3; ++coordinate_index) {
        ap[coordinate_index] /= aw;
        bp[coordinate_index] /= bw;
      }
      a.control.back() = ap;
      b.control.front() = bp;
      trace.push_back(identity + ": compatible homogeneous endpoint average");
    }
    ++curve_index;
  }
  // The later TrimmedSurface constructor joins Euclidean endpoints, but
  // writes their midpoint to homogeneous slots while retaining the weights.
  for (std::size_t curve_index = 0; curve_index < curves.size();
       ++curve_index) {
    auto& a = curves[curve_index];
    auto& b = curves[(curve_index + 1) % curves.size()];
    auto delta = Subtract(a.control.back(), b.control.front());
    double gap = std::sqrt((delta[1] * delta[1] + delta[0] * delta[0]) +
                           delta[2] * delta[2]);
    if (gap >= 1e-6) {
      auto mean = Multiply(Add(a.control.back(), b.control.front()), .5);
      auto ap = mean;
      auto bp = mean;
      for (int coordinate_index = 0; coordinate_index < 3; ++coordinate_index) {
        ap[coordinate_index] /= a.weights.back();
        bp[coordinate_index] /= b.weights.front();
      }
      a.control.back() = ap;
      b.control.front() = bp;
      trace.push_back(identity +
                      ": compatible surface-constructor endpoint join");
    }
  }
  // LucidShape's NCurve constructor clamps stored weights while preserving
  // HXYZ.
  for (auto& curve : curves) {
    auto old = curve.weights;
    for (auto& weight : curve.weights) {
      weight = std::clamp(weight, .001, 100.);
    }
    if (old != curve.weights) {
      for (std::size_t control_point_index = 0;
           control_point_index < old.size(); ++control_point_index) {
        for (auto& coordinate : curve.control[control_point_index]) {
          coordinate = coordinate * old[control_point_index] /
                       curve.weights[control_point_index];
        }
      }
      trace.push_back(identity + ": compatible constructor weight clamp");
    }
    if (std::all_of(
            curve.weights.begin(), curve.weights.end(),
            [&](double weight) { return weight == curve.weights.front(); })) {
      std::fill(curve.weights.begin(), curve.weights.end(), 1.);
    }
  }
}

namespace {
void PlaneDomain(Face& f) {
  const auto old = f.surface;
  Vec2 lo{INFINITY, INFINITY};
  Vec2 hi{-INFINITY, -INFINITY};
  auto include = [&](const auto& loop) {
    for (const auto& curve : loop) {
      // Analytic conic bounds survive conversion to rational control points.
      // Straight segments retain imported bounds for sub-confusion edits;
      // larger endpoint changes use their updated control envelope.
      bool imported_bounds =
          curve.import_bounds && !curve.linear_conversion && curve.degree != 1;
      if (curve.import_bounds && curve.degree == 1) {
        Rect envelope{INFINITY, -INFINITY, INFINITY, -INFINITY};
        for (auto control_point : curve.control) {
          for (int axis = 0; axis < 2; ++axis) {
            envelope[2 * axis] =
                std::min(envelope[2 * axis], control_point[axis]);
            envelope[2 * axis + 1] =
                std::max(envelope[2 * axis + 1], control_point[axis]);
          }
        }
        double movement = 0.;
        for (int boundary_coordinate_index = 0; boundary_coordinate_index < 4;
             ++boundary_coordinate_index) {
          movement = std::max(
              movement,
              std::abs(envelope[boundary_coordinate_index] -
                       (*curve.import_bounds)[boundary_coordinate_index]));
        }
        imported_bounds = movement <= kParameterConfusion;
      }
      if (imported_bounds) {
        auto b = *curve.import_bounds;
        for (int axis = 0; axis < 2; ++axis) {
          lo[axis] = std::min(lo[axis], b[2 * axis]);
          hi[axis] = std::max(hi[axis], b[2 * axis + 1]);
        }
      } else {
        for (auto control_point : curve.control) {
          for (int axis = 0; axis < 2; ++axis) {
            lo[axis] = std::min(lo[axis], control_point[axis]);
            hi[axis] = std::max(hi[axis], control_point[axis]);
          }
        }
      }
    }
  };
  include(*f.outer);
  for (const auto& loop : f.inner) {
    include(loop);
  }
  Rect domain{lo[0], hi[0], lo[1], hi[1]};
  lo = {domain[0], domain[2]};
  hi = {domain[1], domain[3]};
  if (domain == old.domain) {
    return;
  }
  f.surface.ku = {lo[0], lo[0], hi[0], hi[0]};
  f.surface.kv = {lo[1], lo[1], hi[1], hi[1]};
  f.surface.domain = domain;
  f.surface.control.clear();
  for (double v_parameter : {lo[1], hi[1]}) {
    for (double u_parameter : {lo[0], hi[0]}) {
      f.surface.control.push_back(
          old.EvaluateControl({u_parameter, v_parameter}, true));
    }
  }
}
}  // namespace

void AdaptOcctFace(Face& face, std::vector<std::string>& trace) {
  if (!face.outer) {
    throw Unsupported("OCCT face adaptation requires an explicit outer wire");
  }
  const auto& s = face.surface;
  bool nonclamped = s.ku.front() != s.domain[0] || s.ku.back() != s.domain[1] ||
                    s.kv.front() != s.domain[2] || s.kv.back() != s.domain[3];
  if (nonclamped) {
    bool rectangle = NaturalRectangle(face, true);
    face.surface = ClampActiveDomain(face.surface);
    // This conversion retains OCCT's wire; it does not recreate a source-file
    // wire.
    face.constructed_rectangle = face.support_type == "Geom_BSplineSurface" &&
                                 rectangle && face.outer->size() == 4;
    trace.push_back(face.id +
                    ": compatible active-domain endpoint knot insertion; OCCT "
                    "wire retained");
  }
  bool declared_closed_natural =
      face.source_natural_boundary &&
      (face.source_closed_axes[0] || face.source_closed_axes[1]);
  if ((!nonclamped || declared_closed_natural) &&
      (face.support_type == "Geom_BSplineSurface" ||
       face.support_type == "Geom_Plane") &&
      NaturalRectangle(face, declared_closed_natural)) {
    face.outer.reset();
    face.constructed_rectangle = false;
    trace.push_back(face.id +
                    ": eliminate natural rectangular trim; topology retained");
    return;
  }
  JoinOcctConversionEndpoints(*face.outer, face.boundary_edges[0], trace,
                              face.id);
  for (std::size_t inner_loop_index = 0; inner_loop_index < face.inner.size();
       ++inner_loop_index) {
    JoinOcctConversionEndpoints(face.inner[inner_loop_index],
                                face.boundary_edges[inner_loop_index + 1],
                                trace, face.id);
  }
  if (face.support_type == "Geom_Plane") {
    PlaneDomain(face);
  }
  face.surface.Validate();
}
}  // namespace iko::internal
