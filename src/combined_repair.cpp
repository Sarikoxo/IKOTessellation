#include "tessellation/private/repair.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "tessellation/private/bezier.h"
#include "tessellation/tessellation.h"
#include "tessellation/private/internal.h"

namespace iko::internal {
static Vec3 EndpointD1(const Curve& c, bool end) {
  int d = c.degree;
  if ((d != 2 && d != 3) ||
      c.control.size() != static_cast<std::size_t>(d + 1)) {
    throw Unsupported("Endpoint D1 requires a full quadratic/cubic Bezier");
  }
  std::vector<Vec3> poles;
  for (int control_point_index = 0; control_point_index <= d;
       ++control_point_index) {
    poles.push_back(
        {c.control[control_point_index][0] * c.weights[control_point_index],
         c.control[control_point_index][1] * c.weights[control_point_index],
         c.weights[control_point_index]});
  }
  double span = c.domain[1] - c.domain[0];
  double reciprocal = 1. / span;
  for (int derivative_order = 0; derivative_order < d; ++derivative_order) {
    for (int control_point_index = d - 1;
         control_point_index >= derivative_order; --control_point_index) {
      poles[control_point_index + 1] = Multiply(
          Subtract(poles[control_point_index + 1], poles[control_point_index]),
          reciprocal);
    }
  }
  double coefficient = end ? span : 0.;
  for (int evaluation_stage = 0; evaluation_stage < d; ++evaluation_stage) {
    for (int control_point_index = evaluation_stage; control_point_index >= 0;
         --control_point_index) {
      poles[control_point_index] =
          Add(poles[control_point_index],
              Multiply(poles[control_point_index + 1], coefficient));
    }
  }
  return {poles[0][0] * (1. / poles[0][2]), poles[0][1] * (1. / poles[0][2]),
          0};
}
struct Vertex {
  Vec3 position;
  double tolerance;
};
struct Record {
  Curve basis;
  Curve original;
  int source = -1;
  int direction = 1;
  std::array<double, 2> range, face_range;
  std::array<int, 2> vertices;
  bool range_alias = true;
};
static bool EndpointContacts(const Record& a, const Record& b) {
  for (int first_endpoint_index = 0; first_endpoint_index < 2;
       ++first_endpoint_index) {
    for (int second_endpoint_index = 0; second_endpoint_index < 2;
         ++second_endpoint_index) {
      double ta = a.basis.domain[first_endpoint_index];
      double tb = b.basis.domain[second_endpoint_index];
      if (a.range[0] <= ta && ta <= a.range[1] && b.range[0] <= tb &&
          tb <= b.range[1] &&
          a.basis.control[first_endpoint_index ? a.basis.control.size() - 1
                                               : 0] ==
              b.basis.control[second_endpoint_index ? b.basis.control.size() - 1
                                                    : 0]) {
        return true;
      }
    }
  }
  return false;
}
static void UnionVertices(std::vector<Record>& records,
                          std::vector<Vertex>& vertices, int first_record_index,
                          int second_record_index) {
  auto& a = records[first_record_index];
  auto& b = records[second_record_index];
  double d[4];
  for (int first_endpoint_index = 0; first_endpoint_index < 2;
       ++first_endpoint_index) {
    for (int second_endpoint_index = 0; second_endpoint_index < 2;
         ++second_endpoint_index) {
      d[2 * first_endpoint_index + second_endpoint_index] =
          Norm(Subtract(vertices[a.vertices[first_endpoint_index]].position,
                        vertices[b.vertices[second_endpoint_index]].position));
    }
  }
  int nearest_endpoint_pair_index = d[0] < std::min({d[1], d[2], d[3]}) ? 0
                                    : d[1] < std::min(d[2], d[3])       ? 1
                                    : d[2] < d[3]                       ? 2
                                                                        : 3;
  int va = a.vertices[nearest_endpoint_pair_index / 2];
  int vb = b.vertices[nearest_endpoint_pair_index % 2];
  double tol = std::max(vertices[va].tolerance, vertices[vb].tolerance);
  if (va == vb || d[nearest_endpoint_pair_index] >= tol) {
    return;
  }
  vertices[va].tolerance = tol;
  for (int record_index :
       {second_record_index,
        (second_record_index - 1 + static_cast<int>(records.size())) %
            static_cast<int>(records.size()),
        (second_record_index + 1) % static_cast<int>(records.size())}) {
    auto& r = records[record_index];
    if (r.vertices[0] != vb && r.vertices[1] != vb) {
      continue;
    }
    for (auto& vertex_id : r.vertices) {
      if (vertex_id == vb) {
        vertex_id = va;
      }
    }
    r.range_alias = false;
  }
}
RepairResult RepairCombined(const Surface& s, const std::vector<WireItem>& wire,
                            std::vector<std::string>& trace) {
  std::vector<Curve> original;
  std::vector<Curve> curves;
  for (const auto& wire_item : wire) {
    auto c = wire_item.curve;
    if ((c.degree != 2 && c.degree != 3) ||
        c.control.size() != static_cast<std::size_t>(c.degree + 1)) {
      throw Unsupported(
          "Combined added-edge cuts require full quadratic/cubic Bezier "
          "survivors");
    }
    ExactBezier checked(c);
    original.push_back(c);
    curves.push_back(c);
  }
  std::vector<double> tolerances(wire.size(), 1e-7);
  std::map<int, Curve> lines;
  std::vector<bool> common(curves.size());
  auto bounds = ResolutionBounds(s);
  for (int curve_index = 0; curve_index < static_cast<int>(curves.size());
       ++curve_index) {
    int previous_curve_index =
        (curve_index - 1 + curves.size()) % curves.size();
    auto a = curves[previous_curve_index];
    auto b = curves[curve_index];
    auto end_a = EndpointD1(a, true);
    auto start_b = EndpointD1(b, false);
    auto gap = Subtract(start_b, end_a);
    double length = std::sqrt(gap[0] * gap[0] + gap[1] * gap[1]);
    if (length <= 1e-12) {
      continue;
    }
    auto midpoint = Multiply(Add(end_a, start_b), .5);
    auto ca = a;
    auto cb = b;
    ca.control.back() = midpoint;
    cb.control.front() = midpoint;
    double ta = CandidateTolerance(s, original[previous_curve_index], ca);
    double tb = CandidateTolerance(s, original[curve_index], cb);
    auto tangent_a =
        Subtract(a.control.back(), a.control[a.control.size() - 2]);
    auto tangent_b = Subtract(b.control[1], b.control.front());
    Vec2 tangents[] = {{tangent_a[0], tangent_a[1]},
                       {tangent_b[0], tangent_b[1]}};
    Vec2 gap2{gap[0], gap[1]};
    double n0 = Norm(tangents[0]);
    double n1 = Norm(tangents[1]);
    if (std::min(n0, n1) < 1e-12) {
      throw Unsupported("Zero endpoint tangent in combined IGES bending");
    }
    double cosines[] = {Dot(gap2, tangents[0]) / (length * n0),
                        Dot(gap2, tangents[1]) / (length * n1)};
    double threshold = std::cos(.9 * 3.14159265358979323846);
    if (std::min(std::abs(cosines[0] - threshold),
                 std::abs(cosines[1] - threshold)) < 1e-5) {
      throw Unsupported(
          "Combined bending tangent is close to the zigzag threshold");
    }
    auto xyz1 = s.Evaluate({a.control.back()[0], a.control.back()[1]});
    auto xyz2 = s.Evaluate({b.control.front()[0], b.control.front()[1]});
    if (Norm(Subtract(xyz2, xyz1)) <= .0001) {
      throw Unsupported(
          "Small nonzero combined gaps require additional tolerance healing");
    }
    auto ox1 = s.Evaluate({original[previous_curve_index].control.back()[0],
                           original[previous_curve_index].control.back()[1]});
    auto ox2 = s.Evaluate({original[curve_index].control.front()[0],
                           original[curve_index].control.front()[1]});
    double tol1 =
        std::max(tolerances[previous_curve_index], Norm(Subtract(ox1, xyz1)));
    double tol2 = std::max(tolerances[curve_index], Norm(Subtract(ox2, xyz2)));
    double required = .5 * std::min(bounds[0], bounds[1]) * length;
    bool existing = ta + .001 < tolerances[previous_curve_index] &&
                    tb + .001 < tolerances[curve_index];
    bool zigzag = cosines[0] < threshold || cosines[1] < threshold;
    bool add_line =
        !existing && !zigzag &&
        Norm(Subtract(ox2, ox1)) > std::sqrt(1.25) * (tol1 + tol2) + .001;
    bool common_line = !existing && !add_line && required > 1.001;
    common[curve_index] = common_line;
    if (add_line || common_line) {
      auto delta = gap;
      for (auto& direction_coordinate : delta) {
        direction_coordinate = direction_coordinate / length * length;
      }
      auto end = Add(end_a, delta);
      Curve line;
      line.degree = 1;
      line.knots = {0, 0, length, length};
      line.control = {end_a, end};
      line.weights = {1, 1};
      line.domain = {0, length};
      lines[curve_index] = std::move(line);
      trace.push_back(common_line ? "fix_lacking: add_common_vertex_line"
                                  : "fix_lacking: add_length_line");
    } else {
      if (!existing &&
          !(std::max(ta, tb) + .001 < required && required < .999)) {
        throw Unsupported(
            "Mixed bending decision outside the existing guard margin");
      }
      curves[previous_curve_index] = ca;
      curves[curve_index] = cb;
      tolerances[previous_curve_index] =
          std::max(tolerances[previous_curve_index], ta);
      tolerances[curve_index] = std::max(tolerances[curve_index], tb);
      trace.push_back("fix_lacking: bend_endpoints");
    }
  }
  std::vector<Vertex> vertices;
  std::vector<std::array<int, 2>> junctions;
  for (int junction_index = 0; junction_index < static_cast<int>(curves.size());
       ++junction_index) {
    int previous_curve_index =
        (junction_index - 1 + curves.size()) % curves.size();
    Vec3 p0 = wire[previous_curve_index].endpoints[1];
    Vec3 p1 = wire[junction_index].endpoints[0];
    std::array<int, 2> ids;
    bool distinct = lines.count(junction_index) && !common[junction_index];
    if (!distinct) {
      p0 = Multiply(Add(p0, p1), .5);
      int id = static_cast<int>(vertices.size());
      vertices.push_back({p0, std::max(tolerances[previous_curve_index],
                                       tolerances[junction_index])});
      ids = {id, id};
    } else {
      ids = {static_cast<int>(vertices.size()),
             static_cast<int>(vertices.size()) + 1};
      vertices.push_back({p0, tolerances[previous_curve_index]});
      vertices.push_back({p1, tolerances[junction_index]});
    }
    junctions.push_back(ids);
  }
  std::vector<Record> records;
  for (int curve_index = 0; curve_index < static_cast<int>(curves.size());
       ++curve_index) {
    if (lines.count(curve_index)) {
      const auto& c = lines.at(curve_index);
      records.push_back(
          {c, c, -1, 1, c.domain, c.domain, junctions[curve_index]});
    }
    auto c = curves[curve_index];
    records.push_back({c,
                       original[curve_index],
                       wire[curve_index].source,
                       wire[curve_index].direction,
                       c.domain,
                       c.domain,
                       {junctions[curve_index][1],
                        junctions[(curve_index + 1) % curves.size()][0]}});
  }
  if (lines.empty()) {
    throw Unsupported(
        "Probe is restricted to combinations containing added lines");
  }
  std::vector<Curve> cs;
  for (auto record : records) {
    cs.push_back(record.basis);
  }
  double area = SignedArea(SampledWire(cs, 65));
  if (std::abs(area) < 1e-8) {
    throw Unsupported("Ambiguous mixed-wire orientation");
  }
  std::vector<int> order(records.size());
  std::iota(order.begin(), order.end(), 0);
  if (area < 0) {
    std::reverse(order.begin(), order.end());
  }
  double maximum =
      std::min(*std::max_element(tolerances.begin(), tolerances.end()), 1.);
  for (int first_order_index = 0;
       first_order_index < static_cast<int>(order.size());
       ++first_order_index) {
    for (int second_order_index = first_order_index + 2;
         second_order_index < static_cast<int>(order.size());
         ++second_order_index) {
      if (first_order_index == 0 &&
          second_order_index == static_cast<int>(order.size()) - 1) {
        continue;
      }
      int first_record_index = order[first_order_index];
      int second_record_index = order[second_order_index];
      auto& a = records[first_record_index];
      auto& b = records[second_record_index];
      if (a.basis.degree > 1 && b.basis.degree > 1) {
        if (!CertifyDisjoint(a.basis, b.basis, true)) {
          throw Unsupported(
              "Nonlinear/nonlinear contact outside the mixed-line "
              "hypothesis");
        }
        if (EndpointContacts(a, b)) {
          UnionVertices(records, vertices, first_record_index,
                        second_record_index);
        }
        continue;
      }
      if (a.basis.degree == 1 && b.basis.degree == 1 &&
          CertifyDisjoint(a.basis, b.basis, true)) {
        if (EndpointContacts(a, b)) {
          UnionVertices(records, vertices, first_record_index,
                        second_record_index);
        }
        continue;
      }
      bool first = a.basis.degree == 1;
      auto contacts = first ? LineContacts(a.basis, b.basis)
                            : LineContacts(b.basis, a.basis);
      int best = -1;
      std::array<double, 2> selected{};
      for (auto contact_parameters : contacts) {
        if (!first) {
          std::swap(contact_parameters[0], contact_parameters[1]);
        }
        if (contact_parameters[0] < a.range[0] - 1e-10 ||
            contact_parameters[0] > a.range[1] + 1e-10 ||
            contact_parameters[1] < b.range[0] - 1e-10 ||
            contact_parameters[1] > b.range[1] + 1e-10) {
          continue;
        }
        int score = (a.range[0] + 1e-10 < contact_parameters[0] &&
                             contact_parameters[0] < a.range[1] - 1e-10
                         ? 1
                         : 0) +
                    (b.range[0] + 1e-10 < contact_parameters[1] &&
                             contact_parameters[1] < b.range[1] - 1e-10
                         ? 2
                         : 0);
        if (score > best) {
          best = score;
          selected = contact_parameters;
        }
      }
      if (best < 0) {
        continue;
      }
      if (best == 0) {
        UnionVertices(records, vertices, first_record_index,
                      second_record_index);
        continue;
      }
      if (best != 3) {
        throw Unsupported(
            "Endpoint/interior split outside the mixed-line hypothesis");
      }
      for (int contact_curve_index = 0; contact_curve_index < 2;
           ++contact_curve_index) {
        auto& r = contact_curve_index ? b : a;
        double t = selected[contact_curve_index];
        const auto& old = r.original;
        auto point = old.Evaluate(t);
        auto p3 = s.Evaluate({point[0], point[1]});
        auto v0 = s.Evaluate({old.control.front()[0], old.control.front()[1]});
        auto v1 = s.Evaluate({old.control.back()[0], old.control.back()[1]});
        double d0 = Norm(Subtract(v0, p3));
        double d1 = Norm(Subtract(v1, p3));
        if (std::min(d0, d1) + .001 >= maximum || std::abs(d0 - d1) < .001) {
          throw Unsupported(
              "3D approximation/vertex cut decision outside the guard "
              "margin");
        }
        double lo = r.range[0];
        double hi = r.range[1];
        double pend = std::abs(lo - t) > std::abs(hi - t) ? lo : hi;
        if (std::abs(pend - t) < 1e-8) {
          throw Unsupported("Nearly collapsed mixed cut");
        }
        r.range = {std::min(pend, t), std::max(pend, t)};
        if (r.range_alias) {
          r.face_range = r.range;
        }
        trace.push_back("intersection_cut: keep_farther_parameter_endpoint");
      }
    }
  }
  std::vector<WireItem> converted;
  for (const auto& record : records) {
    auto c = record.basis;
    if (c.degree == 1) {
      double lo = record.face_range[0];
      double hi = record.face_range[1];
      double span = c.domain[1] - c.domain[0];
      auto unit = Subtract(c.control[1], c.control[0]);
      for (auto& direction_coordinate : unit) {
        direction_coordinate /= span;
      }
      auto a = Add(c.control[0], Multiply(unit, lo - c.domain[0]));
      auto b = Add(c.control[0], Multiply(unit, hi - c.domain[0]));
      c.control = {a, b};
      c.knots = {lo, lo, hi, hi};
      c.domain = {lo, hi};
    }
    converted.push_back({c, record.source, record.direction, {}});
  }
  trace.push_back("occt_conversion: preserve_basis_and_active_range");
  return {converted, area < 0};
}
}  // namespace iko::internal
