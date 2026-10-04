#include "tessellation/private/internal.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <ios>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "tessellation/tessellation.h"

namespace iko::internal {
bool Refine(const Surface& s, const Options& o, const std::array<Vec2, 3>& uv) {
  // Refinement uses the prepared U-then-V control-net evaluation, including
  // the polynomial/rational flag. Dense normalized basis summation can change
  // edge-threshold decisions even when the geometric point is equivalent.
  Vec3 p[3];
  for (int triangle_vertex_index = 0; triangle_vertex_index < 3;
       ++triangle_vertex_index) {
    p[triangle_vertex_index] = s.EvaluateControl(uv[triangle_vertex_index]);
  }
  auto a = Subtract(p[1], p[0]);
  auto b = Subtract(p[2], p[0]);
  auto c = Subtract(p[2], p[1]);
  double lengths[] = {Dot(a, a), Dot(b, b), Dot(c, c)};
  if (*std::min_element(lengths, lengths + 3) <
      o.surface.min_edge * o.surface.min_edge) {
    return false;
  }
  double max_edge = std::max(o.surface.max_edge, .001 * s.BboxSpan());
  if (*std::max_element(lengths, lengths + 3) > max_edge * max_edge) {
    return true;
  }
  if (std::min(lengths[0], lengths[1]) <= std::numeric_limits<double>::min()) {
    return false;
  }
  Vec2 centre{(uv[0][0] + uv[1][0] + uv[2][0]) / 3.,
              (uv[0][1] + uv[1][1] + uv[2][1]) / 3.};
  auto q = s.EvaluateControl(centre);
  double scalar = Dot(Cross(a, b), Subtract(p[0], q));
  double error = scalar * scalar / (lengths[0] * lengths[1]);
  double span =
      .2 * std::min(s.domain[1] - s.domain[0], s.domain[3] - s.domain[2]);
  double factor = std::min(1., std::max(.04, span * span));
  return error > o.surface.chord * o.surface.chord * factor;
}
void Append(Mesh& dest, Mesh source) {
  if (dest.xyz.size() + source.xyz.size() >
      std::numeric_limits<std::uint32_t>::max()) {
    throw std::runtime_error("Mesh index limit reached");
  }
  auto offset = static_cast<std::uint32_t>(dest.xyz.size());
  for (auto triangle : source.triangles) {
    for (auto& vertex_index : triangle) {
      vertex_index += offset;
    }
    dest.triangles.push_back(triangle);
  }
  dest.xyz.insert(dest.xyz.end(), source.xyz.begin(), source.xyz.end());
  dest.uv.insert(dest.uv.end(), source.uv.begin(), source.uv.end());
  dest.normals.insert(dest.normals.end(), source.normals.begin(),
                      source.normals.end());
  dest.identities.insert(dest.identities.end(), source.identities.begin(),
                         source.identities.end());
  dest.unpaired_uv_tail.insert(dest.unpaired_uv_tail.end(),
                               source.unpaired_uv_tail.begin(),
                               source.unpaired_uv_tail.end());
}
static Mesh Records(const Face& f, const std::vector<Vec2>& uv,
                    const std::vector<Triangle>& tris, const Rect* r,
                    const std::string& identity) {
  Mesh out;
  out.triangles = tris;
  out.normals = MeshNormals(f.surface, uv, r);
  out.uv = uv;
  for (auto parameter_point : uv) {
    out.xyz.push_back(f.surface.Evaluate(parameter_point));
  }
  for (auto& vertex_position : out.xyz) {
    for (auto& coordinate : vertex_position) {
      coordinate = Quantize(coordinate);
    }
  }
  for (auto& parameter_point : out.uv) {
    for (auto& coordinate : parameter_point) {
      coordinate = Quantize(coordinate);
    }
  }
  for (auto& vertex_normal : out.normals) {
    for (auto& coordinate : vertex_normal) {
      coordinate = Quantize(coordinate);
    }
  }
  out.identities.assign(uv.size(), identity);
  return out;
}
static std::vector<Vec2> Edge(const Surface& s, const Options& o, int axis,
                              double fixed, const std::vector<double>& values) {
  std::vector<Vec2> out;
  for (std::size_t parameter_span_index = 0;
       parameter_span_index + 1 < values.size(); ++parameter_span_index) {
    auto uv = [=](double t) { return axis ? Vec2{fixed, t} : Vec2{t, fixed}; };
    auto initial = [&](double t) {
      auto p = uv(t);
      return Vec3{p[0], p[1], 0};
    };
    auto r = Sample([&](double t) { return s.Point(uv(t)); },
                    values[parameter_span_index],
                    values[parameter_span_index + 1], o.curve.chord,
                    o.curve.min_edge, o.curve.max_edge, false, 1, initial);
    for (std::size_t sample_index = 0; sample_index + 1 < r.parameters.size();
         ++sample_index) {
      out.push_back(uv(r.parameters[sample_index]));
    }
  }
  out.push_back(axis ? Vec2{fixed, values.back()} : Vec2{values.back(), fixed});
  return out;
}
static PlanarStraightLineGraph Network(const Surface& s, const Options& o,
                                       const std::vector<Rect>& blocks) {
  std::vector<double> u;
  std::vector<double> v;
  for (auto knot_block : blocks) {
    u.push_back(knot_block[0]);
    u.push_back(knot_block[1]);
    v.push_back(knot_block[2]);
    v.push_back(knot_block[3]);
  }
  u = Unique(u);
  v = Unique(v);
  auto south = Edge(s, o, 0, v.front(), u);
  auto east = Edge(s, o, 1, u.back(), v);
  std::vector<Vec2> north;
  std::vector<Vec2> west;
  for (auto south_vertex_iterator = south.rbegin();
       south_vertex_iterator != south.rend(); ++south_vertex_iterator) {
    north.push_back({(*south_vertex_iterator)[0], v.back()});
  }
  for (auto east_vertex_iterator = east.rbegin();
       east_vertex_iterator != east.rend(); ++east_vertex_iterator) {
    west.push_back({u.front(), (*east_vertex_iterator)[1]});
  }
  PlanarStraightLineGraph p;
  std::map<Vec2, int> lookup;
  auto vertex = [&](Vec2 x) {
    auto found = lookup.find(x);
    if (found != lookup.end()) {
      return found->second;
    }
    int id = static_cast<int>(p.uv.size());
    lookup[x] = id;
    p.uv.push_back(x);
    return id;
  };
  std::vector<int> loop;
  for (const auto& boundary_points : {south, east, north, west}) {
    for (std::size_t boundary_vertex_index = 0;
         boundary_vertex_index + 1 < boundary_points.size();
         ++boundary_vertex_index) {
      loop.push_back(vertex(boundary_points[boundary_vertex_index]));
    }
  }
  for (std::size_t loop_vertex_index = 0; loop_vertex_index < loop.size();
       ++loop_vertex_index) {
    p.segments.push_back(
        {loop[loop_vertex_index], loop[(loop_vertex_index + 1) % loop.size()]});
  }
  for (std::size_t u_knot_index = 1; u_knot_index + 1 < u.size();
       ++u_knot_index) {
    p.segments.push_back({vertex({u[u_knot_index], v.front()}),
                          vertex({u[u_knot_index], v.back()})});
  }
  for (std::size_t v_knot_index = 1; v_knot_index + 1 < v.size();
       ++v_knot_index) {
    p.segments.push_back({vertex({u.back(), v[v_knot_index]}),
                          vertex({u.front(), v[v_knot_index]})});
  }
  p.internal_constraints = static_cast<int>(u.size() + v.size()) - 4;
  return p;
}
Mesh Generate(const Face& face, const Options& o) {
  const auto& s = face.surface;
  s.Validate();
  if (face.outer || !face.inner.empty()) {
    return BuildTrim(face, o);
  }
  if (o.smooth_poles) {
    auto pole = DetectPole(s);
    if (pole) {
      return BuildPole(face, o, *pole);
    }
  }
  ArcMap map(s);
  auto blocks = KnotBlocks(s);
  Mesh out;
  out.face_id = face.id;
  out.identity_mode = "knot_blocks";
  out.explicit_weights = s.explicit_weights;
  out.branch = "ordinary";
  auto run = [&](PlanarStraightLineGraph p) {
    for (auto& parameter_point : p.uv) {
      parameter_point = map.Apply(parameter_point);
    }
    auto result = TriangleRun(p, o, [&](const std::array<Vec2, 3>& uv) {
      std::array<Vec2, 3> original;
      for (int triangle_vertex_index = 0; triangle_vertex_index < 3;
           ++triangle_vertex_index) {
        original[triangle_vertex_index] =
            map.Apply(uv[triangle_vertex_index], true);
      }
      return Refine(s, o, original);
    });
    for (auto& parameter_point : result.uv) {
      parameter_point = map.Apply(parameter_point, true);
    }
    return result;
  };
  auto result = run(Network(s, o, blocks));
  std::vector<int> coverage(result.triangles.size());
  for (auto knot_block : blocks) {
    std::vector<Triangle> tris;
    std::set<std::uint32_t> indices;
    for (std::size_t triangle_index = 0;
         triangle_index < result.triangles.size(); ++triangle_index) {
      auto t = result.triangles[triangle_index];
      auto a = result.uv[t[0]];
      auto b = result.uv[t[1]];
      auto c = result.uv[t[2]];
      Vec2 centre{(a[0] + b[0] + c[0]) / 3., (a[1] + b[1] + c[1]) / 3.};
      if (centre[0] > knot_block[0] && centre[0] < knot_block[1] &&
          centre[1] > knot_block[2] && centre[1] < knot_block[3]) {
        ++coverage[triangle_index];
        tris.push_back(t);
        indices.insert(t.begin(), t.end());
      }
    }
    if (indices.empty()) {
      continue;
    }
    std::map<std::uint32_t, std::uint32_t> remap;
    std::vector<Vec2> uv;
    for (auto vertex_index : indices) {
      remap[vertex_index] = static_cast<std::uint32_t>(uv.size());
      uv.push_back(result.uv[vertex_index]);
    }
    for (auto& triangle : tris) {
      for (auto& vertex_index : triangle) {
        vertex_index = remap.at(vertex_index);
      }
    }
    Append(out, Records(face, uv, tris, &knot_block,
                        face.id + "|" + BlockKey(knot_block)));
  }
  if (std::any_of(coverage.begin(), coverage.end(),
                  [](int coverage_count) { return coverage_count != 1; })) {
    throw std::runtime_error(
        "Global mesh contains unassigned/ambiguous knot-block triangles");
  }
  return out;
}
static std::string HexDouble(double v) {
  if (v == 0) {
    return "0x0.0p+0";
  }
  // All component bounds have already been quantized to binary32. Format
  // the binary64 mantissa with Python float.hex's fixed 13 fractional digits.
  int exponent;
  double m = std::frexp(std::abs(v), &exponent) * 2;
  --exponent;
  auto fraction = static_cast<std::uint64_t>((m - 1.) * 4503599627370496.);
  std::ostringstream out;
  if (std::signbit(v)) {
    out << '-';
  }
  out << "0x1." << std::hex << std::setw(13) << std::setfill('0') << fraction
      << "p" << std::dec << (exponent >= 0 ? "+" : "") << exponent;
  return out.str();
}
void ComponentIdentities(Mesh& m, const std::string& face) {
  std::vector<std::uint32_t> parent(m.uv.size());
  std::iota(parent.begin(), parent.end(), 0);
  auto root = [&](std::uint32_t vertex_index) {
    while (parent[vertex_index] != vertex_index) {
      parent[vertex_index] = parent[parent[vertex_index]];
      vertex_index = parent[vertex_index];
    }
    return vertex_index;
  };
  for (auto triangle : m.triangles) {
    auto a = root(triangle[0]);
    parent[root(triangle[1])] = a;
    parent[root(triangle[2])] = a;
  }
  std::map<std::uint32_t, std::vector<std::uint32_t>> groups;
  for (std::uint32_t vertex_index = 0; vertex_index < m.uv.size();
       ++vertex_index) {
    groups[root(vertex_index)].push_back(vertex_index);
  }
  m.identities.resize(m.uv.size());
  for (const auto& component_group : groups) {
    auto p = m.uv[component_group.second.front()];
    Rect bounds{p[0], p[1], p[0], p[1]};
    for (auto vertex_index : component_group.second) {
      for (int axis = 0; axis < 2; ++axis) {
        bounds[axis] = std::min(bounds[axis], m.uv[vertex_index][axis]);
        bounds[axis + 2] = std::max(bounds[axis + 2], m.uv[vertex_index][axis]);
      }
    }
    std::string label = face + "|index-component-uv-box:";
    for (int bound_index = 0; bound_index < 4; ++bound_index) {
      double v = Quantize(bounds[bound_index]);
      if (std::abs(v) <= 1e-8) {
        v = 0.;
      }
      if (bound_index) {
        label += ',';
      }
      label += HexDouble(v);
    }
    for (auto vertex_index : component_group.second) {
      m.identities[vertex_index] = label;
    }
  }
}
}  // namespace iko::internal
