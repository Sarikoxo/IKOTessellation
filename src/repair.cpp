#include "tessellation/private/repair.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "tessellation/private/bezier.h"
#include "tessellation/tessellation.h"
#include "tessellation/private/internal.h"

namespace iko::internal {
Curve ReverseCurve(Curve c) {
  double origin = c.knots.front() + c.knots.back();
  std::reverse(c.knots.begin(), c.knots.end());
  for (auto& knot : c.knots) {
    knot = origin - knot;
  }
  std::reverse(c.control.begin(), c.control.end());
  std::reverse(c.weights.begin(), c.weights.end());
  c.domain = {c.knots[c.degree], c.knots[c.knots.size() - c.degree - 1]};
  return c;
}
std::vector<Vec2> SampledWire(const std::vector<Curve>& cs, int count) {
  std::vector<Vec2> points;
  for (const auto& curve : cs) {
    double step = (curve.domain[1] - curve.domain[0]) / (count - 1);
    for (int sample_index = 0; sample_index < count - 1; ++sample_index) {
      auto p = curve.Evaluate(curve.domain[0] + sample_index * step);
      points.push_back({p[0], p[1]});
    }
  }
  return points;
}
double SignedArea(const std::vector<Vec2>& p) {
  double area = 0.;
  for (std::size_t vertex_index = 0; vertex_index < p.size(); ++vertex_index) {
    auto a = p[vertex_index];
    auto b = p[(vertex_index + 1) % p.size()];
    area += a[0] * b[1] - a[1] * b[0];
  }
  return area * .5;
}
std::array<double, 2> ResolutionBounds(const Surface& s) {
  if (s.nu != s.nv) {
    throw Unsupported(
        "Bending resolution for rectangular control nets is not verified");
  }
  std::array<double, 2> out;
  for (int axis = 0; axis < 2; ++axis) {
    bool rational = false;
    for (int control_v_index = 0; control_v_index < s.nv - (axis == 1);
         ++control_v_index) {
      for (int control_u_index = 0; control_u_index < s.nu - (axis == 0);
           ++control_u_index) {
        int control_point_index = control_v_index * s.nu + control_u_index;
        int next_control_point_index = control_point_index + (axis ? s.nu : 1);
        rational |= s.weights[control_point_index] !=
                    s.weights[next_control_point_index];
      }
    }
    double bound = 0.;
    for (int control_v_index = 0; control_v_index < s.nv - (axis == 1);
         ++control_v_index) {
      for (int control_u_index = 0; control_u_index < s.nu - (axis == 0);
           ++control_u_index) {
        int control_point_index = control_v_index * s.nu + control_u_index;
        int next_control_point_index = control_point_index + (axis ? s.nu : 1);
        if (rational) {
          for (auto control_point : s.control) {
            auto a = Multiply(
                Subtract(control_point, s.control[next_control_point_index]),
                s.weights[next_control_point_index]);
            auto b = Multiply(
                Subtract(control_point, s.control[control_point_index]),
                s.weights[control_point_index]);
            auto delta = Subtract(a, b);
            double l1 =
                (std::abs(delta[0]) + std::abs(delta[1])) + std::abs(delta[2]);
            bound = std::max(bound, l1);
          }
        } else {
          auto delta = Subtract(s.control[next_control_point_index],
                                s.control[control_point_index]);
          double l1 =
              (std::abs(delta[0]) + std::abs(delta[1])) + std::abs(delta[2]);
          bound = std::max(bound, l1);
        }
      }
    }
    if (rational) {
      bound /= *std::min_element(s.weights.begin(), s.weights.end());
    }
    int degree = axis ? s.q : s.p;
    double span = s.domain[2 * axis + 1] - s.domain[2 * axis];
    out[axis] = bound * degree * std::sqrt(2.) / span;
  }
  if (std::min(out[0], out[1]) <= 0) {
    throw Unsupported("Degenerate support resolution in IGES bending");
  }
  return out;
}
double CandidateTolerance(const Surface& s, const Curve& original,
                          const Curve& candidate) {
  double maximum = 0.;
  double first = 0.;
  double last = 0.;
  double step = (original.domain[1] - original.domain[0]) / 22.;
  for (int sample_index = 0; sample_index < 23; ++sample_index) {
    double t = sample_index == 22 ? original.domain[1]
                                  : original.domain[0] + sample_index * step;
    auto a = original.Evaluate(t);
    auto b = candidate.Evaluate(t);
    double distance =
        Norm(Subtract(s.Evaluate({a[0], a[1]}), s.Evaluate({b[0], b[1]})));
    maximum = std::max(maximum, distance);
    if (sample_index == 0) {
      first = distance;
    }
    if (sample_index == 22) {
      last = distance;
    }
  }
  double endpoints = std::max(first, last);
  if (maximum > endpoints + 1e-10) {
    throw Unsupported(
        "Interior-dominated 3D curve approximation in IGES bending is not "
        "reproduced");
  }
  return 1.00001 * endpoints;
}
RepairResult RepairBending(const Surface& s, const std::vector<WireItem>& wire,
                           std::vector<std::string>& trace) {
  if (wire.size() < 3) {
    throw Unsupported(
        "Bending collapse of fewer than three surviving curves is not "
        "reproduced");
  }
  auto bounds = ResolutionBounds(s);
  auto curves = wire;
  std::vector<double> tolerances(wire.size(), 1e-7);
  for (std::size_t curve_index = 0; curve_index < curves.size();
       ++curve_index) {
    std::size_t previous = (curve_index + curves.size() - 1) % curves.size();
    const auto& a = curves[previous].curve;
    const auto& b = curves[curve_index].curve;
    auto start = a.control.back();
    auto end = b.control.front();
    auto gap = Subtract(end, start);
    Vec2 gap2{gap[0], gap[1]};
    double length = Norm(gap2);
    if (length <= 1e-12) {
      continue;
    }
    auto xyz1 = s.Evaluate({start[0], start[1]});
    auto xyz2 = s.Evaluate({end[0], end[1]});
    if (Norm(Subtract(xyz2, xyz1)) <= .0001) {
      throw Unsupported(
          "Small nonzero rejected-curve gaps require tolerance/bending "
          "healing");
    }
    double required = .5 * std::min(bounds[0], bounds[1]) * length;
    auto midpoint = Multiply(Add(start, end), .5);
    auto ca = a;
    auto cb = b;
    ca.control.back() = midpoint;
    cb.control.front() = midpoint;
    double ta = CandidateTolerance(s, wire[previous].curve, ca);
    double tb = CandidateTolerance(s, wire[curve_index].curve, cb);
    auto tangent_a =
        Subtract(a.control.back(), a.control[a.control.size() - 2]);
    auto tangent_b = Subtract(b.control[1], b.control.front());
    Vec2 tangents[] = {{tangent_a[0], tangent_a[1]},
                       {tangent_b[0], tangent_b[1]}};
    double norms[] = {Norm(tangents[0]), Norm(tangents[1])};
    if (std::min(norms[0], norms[1]) < 1e-12) {
      throw Unsupported("Zero endpoint tangent in IGES bending");
    }
    double cosines[] = {Dot(gap2, tangents[0]) / (length * norms[0]),
                        Dot(gap2, tangents[1]) / (length * norms[1])};
    double threshold = std::cos(.9 * 3.14159265358979323846);
    bool zigzag = cosines[0] < threshold || cosines[1] < threshold;
    bool existing =
        ta + .001 < tolerances[previous] && tb + .001 < tolerances[curve_index];
    auto original1 = s.Evaluate({wire[previous].curve.control.back()[0],
                                 wire[previous].curve.control.back()[1]});
    auto original2 = s.Evaluate({wire[curve_index].curve.control.front()[0],
                                 wire[curve_index].curve.control.front()[1]});
    double tol1 =
        std::max(tolerances[previous], Norm(Subtract(original1, xyz1)));
    double tol2 =
        std::max(tolerances[curve_index], Norm(Subtract(original2, xyz2)));
    double distance = Norm(Subtract(original2, original1));
    double add_threshold = std::sqrt(1.25) * (tol1 + tol2);
    if (!existing &&
        ((!zigzag && distance + .001 >= add_threshold) || required > 1.001)) {
      return RepairCombined(s, wire, trace);
    }
    if (!existing && !(std::max(ta, tb) + .001 < required && required < .999)) {
      throw Unsupported(
          "IGES bending tolerance decision is outside the verified margin");
    }
    if (std::min(std::abs(cosines[0] - threshold),
                 std::abs(cosines[1] - threshold)) < 1e-5) {
      throw Unsupported(
          "IGES bending tangent is close to the zigzag threshold");
    }
    curves[previous].curve = ca;
    curves[curve_index].curve = cb;
    tolerances[previous] = std::max(tolerances[previous], ta);
    tolerances[curve_index] = std::max(tolerances[curve_index], tb);
    trace.push_back("occt_fix_lacking: bend_endpoints " +
                    std::to_string(previous) + " " +
                    std::to_string(curve_index));
  }
  std::vector<Curve> cs;
  for (const auto& wire_item : curves) {
    cs.push_back(wire_item.curve);
  }
  std::optional<CrossingCertificate> certificate;
  std::array<double, 2> areas;
  std::vector<Vec2> points;
  int round = 0;
  for (int sample_count : {65, 129}) {
    points = SampledWire(cs, sample_count);
    if (OuterIntersects(points) && !certificate) {
      certificate = CertifyCrossingWire(cs);
    }
    if (certificate) {
      MatchCrossings(cs, sample_count, *certificate);
    }
    areas[round++] = SignedArea(points);
  }
  Vec2 lo = points[0];
  Vec2 hi = lo;
  for (auto parameter_point : points) {
    for (int axis = 0; axis < 2; ++axis) {
      lo[axis] = std::min(lo[axis], parameter_point[axis]);
      hi[axis] = std::max(hi[axis], parameter_point[axis]);
    }
  }
  double scale = std::max({1., hi[0] - lo[0], hi[1] - lo[1]});
  scale *= scale;
  if (std::min(std::abs(areas[0]), std::abs(areas[1])) < 1e-8 * scale ||
      areas[0] * areas[1] <= 0) {
    throw Unsupported("Ambiguous orientation of a bent IGES wire");
  }
  if (certificate) {
    auto previous_points = SampledWire(cs, 65);
    if (!OuterIntersects(previous_points)) {
      throw Unsupported(
          "Bent sampled boundaries disagree on crossing topology");
    }
    trace.push_back(
        "intersection_and_orientation: retain_certified_crossing_pcurves");
  }
  return {curves, areas[1] < 0};
}
RepairResult ConnectRejected(const Face& face, double rejection,
                             std::vector<std::string>& trace) {
  std::vector<WireItem> wire;
  const auto& source = *face.outer;
  for (std::size_t source_curve_index = 0; source_curve_index < source.size();
       ++source_curve_index) {
    const auto& c = source[source_curve_index];
    if (std::any_of(c.weights.begin(), c.weights.end(),
                    [&](double weight) { return weight < rejection; })) {
      trace.push_back("occt_curve_conversion: reject_curve " +
                      std::to_string(source_curve_index));
      continue;
    }
    bool constant =
        std::all_of(c.weights.begin(), c.weights.end(),
                    [&](double weight) { return weight == c.weights[0]; });
    if (!(c.weights.front() == 1. && c.weights.back() == 1.) && !constant) {
      throw Unsupported(
          "Rejected-curve healing with weighted survivor endpoints is not "
          "reproduced");
    }
    WireItem item{
        c,
        static_cast<int>(source_curve_index),
        1,
        {face.surface.Evaluate({c.control.front()[0], c.control.front()[1]}),
         face.surface.Evaluate({c.control.back()[0], c.control.back()[1]})}};
    if (!wire.empty()) {
      auto first = wire.front().endpoints[0];
      auto last = wire.back().endpoints[1];
      double th = Norm(Subtract(item.endpoints[0], last));
      double tt = Norm(Subtract(item.endpoints[1], last));
      double ht = Norm(Subtract(item.endpoints[1], first));
      double hh = Norm(Subtract(item.endpoints[0], first));
      bool tail_reverse = th > tt;
      bool head_reverse = ht > hh;
      double d1 = tail_reverse ? tt : th;
      double d2 = head_reverse ? hh : ht;
      int choice = d1 > d2 ? (head_reverse ? 3 : 2) : (tail_reverse ? 1 : 0);
      if (std::min(d1, d2) > 1.) {
        throw Unsupported(
            "Rejected-curve transfer exceeds observed maximum connection "
            "tolerance; fallback face not reproduced");
      }
      bool rev_prefix =
          th > .0001 && tt > .0001 && (choice == 2 || choice == 3);
      bool rev_next = rev_prefix ? choice == 2 : tail_reverse;
      if (rev_prefix) {
        std::reverse(wire.begin(), wire.end());
        for (auto& wire_item : wire) {
          wire_item.curve = ReverseCurve(wire_item.curve);
          wire_item.direction = -wire_item.direction;
          std::swap(wire_item.endpoints[0], wire_item.endpoints[1]);
        }
      }
      if (rev_next) {
        item.curve = ReverseCurve(item.curve);
        item.direction = -1;
        std::swap(item.endpoints[0], item.endpoints[1]);
      }
      trace.push_back("occt_connect_next: connected_prefix " +
                      std::to_string(source_curve_index));
    }
    wire.push_back(std::move(item));
  }
  if (wire.empty()) {
    return {};
  }
  if (wire.size() == 1) {
    return {wire, false};
  }
  if (wire.size() < 3) {
    throw Unsupported(
        "Fewer than three surviving IGES curves require additional face/wire "
        "healing");
  }
  bool mixed =
      std::any_of(wire.begin(), wire.end(), [&](const auto& wire_item) {
        return wire_item.direction != wire[0].direction;
      });
  if (mixed) {
    auto result = RepairBending(face.surface, wire, trace);
    std::vector<Vec2> points;
    for (const auto& curve : source) {
      for (std::size_t control_point_index = 0;
           control_point_index + 1 < curve.control.size();
           ++control_point_index) {
        points.push_back({curve.control[control_point_index][0],
                          curve.control[control_point_index][1]});
      }
    }
    result.clockwise = result.clockwise != (SignedArea(points) < 0);
    return result;
  }
  std::vector<WireItem> result;
  for (std::size_t wire_item_index = 0; wire_item_index < wire.size();
       ++wire_item_index) {
    const auto& c = wire[wire_item_index];
    const auto& previous =
        wire[(wire_item_index + wire.size() - 1) % wire.size()];
    auto a = previous.curve.control.back();
    auto b = c.curve.control.front();
    auto delta = Subtract(b, a);
    double length = std::sqrt(delta[0] * delta[0] + delta[1] * delta[1]);
    if (length > 1e-12) {
      double gap = Norm(Subtract(previous.endpoints[1], c.endpoints[0]));
      if (gap <= .0001) {
        throw Unsupported(
            "Small nonzero rejected-curve gaps require tolerance/bending "
            "healing");
      }
      auto end = delta;
      for (auto& direction_coordinate : end) {
        direction_coordinate = direction_coordinate / length * length;
      }
      end = Add(a, end);
      Curve line;
      line.degree = 1;
      line.knots = {0, 0, length, length};
      line.domain = {0, length};
      line.control = {a, end};
      line.weights = {1, 1};
      result.push_back({line, -1, 1, {previous.endpoints[1], c.endpoints[0]}});
      trace.push_back("occt_fix_lacking: insert_length_parameterized_line");
    }
    result.push_back(c);
  }
  return {result, wire.front().direction < 0};
}
}  // namespace iko::internal
