#include "tessellation/private/internal.h"

#include <map>
#include <set>
#include <utility>

namespace iko::internal {
namespace {
// This pipeline never feeds canonical boundaries through the record-exact
// trim cleanup, endpoint editing, or independent face boundary refinement.
constexpr double kParameterMerge = 1e-10;
struct SampleRecord {
  double parameter;
  Vec3 xyz;
  std::uint64_t vertex;
};
struct EdgeRecord {
  const ImportedEdge* source = nullptr;
  std::vector<double> requests;
  std::vector<SampleRecord> samples;
};
struct EdgeUse {
  const Curve* curve;
  int signed_edge;
  // Preserve the exact pcurve parameter of knot/curve sampling requests.
  std::vector<std::pair<double, double>> requests;
};
struct FaceRecord {
  const Face* source;
  std::vector<std::vector<EdgeUse>> loops;
  std::array<std::vector<double>, 2> axes;
};
struct VertexRecord {
  Vec3 xyz;
  std::uint64_t vertex;
};
using Registry = std::map<int, EdgeRecord>;

double Parameter(const std::array<double, 2>& domain, double t) {
  return (1. - t) * domain[0] + t * domain[1];
}
Vec2 UV(const EdgeUse& use, double t) {
  auto p = use.curve->Evaluate(Parameter(use.curve->domain, t));
  if (!Finite(p)) {
    throw Unsupported("Topology edge use has a nonfinite pcurve");
  }
  return {p[0], p[1]};
}
Vec3 XYZ(const ImportedEdge& edge, double t) {
  auto p = edge.point(Parameter(edge.domain, t));
  if (!Finite(p)) {
    throw Unsupported("Topology edge has a nonfinite XYZ curve");
  }
  return p;
}
double Distance2(Vec3 a, Vec3 b) {
  auto delta = Subtract(a, b);
  return Dot(delta, delta);
}
double Tolerance(const Face& face, const ImportedEdge& edge, const Options& options) {
  // A loose imported tolerance must not authorize moving a boundary farther
  // than the user's requested geometry accuracy.
  return std::min({options.curve.chord, options.surface.chord,
                   std::max(edge.tolerance, 1e-8 * std::max(1., face.surface.BboxSpan()))});
}
// Do not assume that a pcurve and the XYZ curve have the same range or speed.
// Search separated local minima, reject multiple valid branches, and validate
// the residual. Endpoint occurrences are resolved by topology, not projection.
// A finite search is deliberately conservative: failure is Unsupported.
double Project(const std::function<Vec3(double)>& curve, Vec3 target,
               double tolerance) {
  constexpr int count = 128;
  std::array<double, count + 1> values;
  for (int i = 0; i <= count; ++i) {
    values[i] = Distance2(curve(static_cast<double>(i) / count), target);
    if (!std::isfinite(values[i])) {
      throw Unsupported("Nonfinite topology parameter mapping");
    }
  }
  std::vector<std::pair<double, double>> candidates;
  candidates.push_back({0., values.front()});
  candidates.push_back({1., values.back()});
  constexpr double ratio = .6180339887498948482;
  for (int i = 0; i <= count; ++i) {
    if ((i > 0 && values[i] > values[i - 1]) ||
        (i < count && values[i] > values[i + 1])) {
      continue;
    }
    double lo = static_cast<double>(std::max(0, i - 1)) / count;
    double hi = static_cast<double>(std::min(count, i + 1)) / count;
    double a = hi - ratio * (hi - lo);
    double b = lo + ratio * (hi - lo);
    double fa = Distance2(curve(a), target);
    double fb = Distance2(curve(b), target);
    for (int iteration = 0; iteration < 70; ++iteration) {
      if (fa <= fb) {
        hi = b; b = a; fb = fa;
        a = hi - ratio * (hi - lo);
        fa = Distance2(curve(a), target);
      } else {
        lo = a; a = b; fa = fb;
        b = lo + ratio * (hi - lo);
        fb = Distance2(curve(b), target);
      }
    }
    double t = fa <= fb ? a : b;
    candidates.push_back({t, std::min(fa, fb)});
  }
  auto best = std::min_element(candidates.begin(), candidates.end(),
                              [](auto a, auto b) { return a.second < b.second; });
  if (best->second > tolerance * tolerance) {
    throw Unsupported("Topology XYZ/pcurve mapping exceeds edge tolerance");
  }
  for (auto candidate : candidates) {
    // Endpoints close to an interior solution are not a second branch.
    if (candidate.second <= tolerance * tolerance &&
        std::abs(candidate.first - best->first) > .5 / count) {
      throw Unsupported("Ambiguous topology XYZ/pcurve parameter mapping");
    }
  }
  return best->first;
}
std::vector<double> Merge(std::vector<double> values) {
  std::sort(values.begin(), values.end());
  std::vector<double> out;
  for (double t : values) {
    if (out.empty() || t - out.back() > kParameterMerge) {
      out.push_back(t);
    }
  }
  if (!out.empty()) {
    out.front() = 0.;
    out.back() = 1.;
  }
  return out;
}
std::vector<double> Breaks(const Curve& c) {
  std::vector<double> values{0., 1.};
  for (double knot : c.knots) {
    if (knot > c.domain[0] && knot < c.domain[1]) {
      values.push_back((knot - c.domain[0]) / (c.domain[1] - c.domain[0]));
    }
  }
  return Unique(values);
}
std::vector<double> KnotCrossings(const EdgeUse& use, int axis, double value) {
  auto breaks = Breaks(*use.curve);
  // Bracket on each knot span. Tangencies do not cross a knot constraint;
  // near coincidences are rejected later if they cannot form a valid PSLG.
  std::vector<double> roots;
  for (std::size_t span = 1; span < breaks.size(); ++span) {
    constexpr int probes = 64;
    double a = breaks[span - 1];
    double fa = UV(use, a)[axis] - value;
    for (int i = 1; i <= probes; ++i) {
      double b = breaks[span - 1] +
                 (breaks[span] - breaks[span - 1]) * i / probes;
      double fb = UV(use, b)[axis] - value;
      if (fa == 0.) roots.push_back(a);
      if (fb == 0.) roots.push_back(b);
      if (fa * fb < 0.) {
        double lo = a, hi = b, flo = fa;
        for (int iteration = 0; iteration < 55; ++iteration) {
          double mid = .5 * (lo + hi);
          double f = UV(use, mid)[axis] - value;
          if (flo * f <= 0.) hi = mid;
          else { lo = mid; flo = f; }
        }
        roots.push_back(.5 * (lo + hi));
      }
      a = b; fa = fb;
    }
  }
  return Unique(roots);
}
void Collect(FaceRecord& face, Registry& registry, const Options& options) {
  const auto& f = *face.source;
  for (auto& loop : face.loops) {
    for (auto& use : loop) {
      auto& record = registry.at(std::abs(use.signed_edge));
      const auto& edge = *record.source;
      auto physical = [&](double t) { return f.surface.EvaluateControl(UV(use, t)); };
      double tolerance = Tolerance(f, edge, options);
      // Check both directed endpoints against the topological edge curve.
      for (int end = 0; end < 2; ++end) {
        double t = use.signed_edge > 0 ? end : 1 - end;
        if (Distance2(physical(end), XYZ(edge, t)) > tolerance * tolerance) {
          throw Unsupported("Topology edge-use endpoints disagree with XYZ edge");
        }
      }
      auto breaks = Breaks(*use.curve);
      double span = .2 * std::min(f.surface.domain[1] - f.surface.domain[0],
                                f.surface.domain[3] - f.surface.domain[2]);
      double factor = std::min(1., std::max(.04, span * span));
      // Refine uses a domain-scaled surface chord. Request a safety margin on
      // every incident face before freezing any shared boundary.
      double boundary_chord = std::min(options.curve.chord,
          .5 * options.surface.chord * std::sqrt(factor));
      std::vector<double> local = breaks;
      for (std::size_t i = 1; i < breaks.size(); ++i) {
        auto sampled = Sample(physical, breaks[i - 1], breaks[i],
                              boundary_chord,
                              std::min(options.curve.min_edge, options.surface.min_edge),
                              std::min(options.curve.max_edge, .45 * options.surface.max_edge),
                              false, 4);
        local.insert(local.end(), sampled.parameters.begin(), sampled.parameters.end());
        auto uv_samples = Sample([&](double t) {
          auto uv = UV(use, t); return Vec3{uv[0], uv[1], 0.};
        }, breaks[i - 1], breaks[i], options.relative_uv_chord, 0., 1., true, 4);
        local.insert(local.end(), uv_samples.parameters.begin(), uv_samples.parameters.end());
      }
      for (int axis = 0; axis < 2; ++axis) {
        for (double knot : face.axes[axis]) {
          auto crossings = KnotCrossings(use, axis, knot);
          local.insert(local.end(), crossings.begin(), crossings.end());
        }
      }
      local = Unique(local);
      double previous = -1.;
      for (double q : local) {
        double t = q == 0. || q == 1.
                       ? (use.signed_edge > 0 ? q : 1. - q)
                       : Project([&](double parameter) { return XYZ(edge, parameter); },
                                 physical(q), tolerance);
        double directed = use.signed_edge > 0 ? t : 1. - t;
        if (directed + kParameterMerge < previous) {
          throw Unsupported("Nonmonotone topology edge-use parameter mapping");
        }
        previous = directed;
        use.requests.push_back({t, q});
        record.requests.push_back(t);
      }
    }
  }
}

bool InsideLoop(Vec2 p, const std::vector<Vec2>& loop) {
  bool inside = false;
  for (std::size_t i = 0, j = loop.size() - 1; i < loop.size(); j = i++) {
    auto a = loop[i], b = loop[j];
    if ((a[1] > p[1]) != (b[1] > p[1]) &&
        p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0]) {
      inside = !inside;
    }
  }
  return inside;
}
bool Inside(Vec2 p, const std::vector<std::vector<Vec2>>& loops) {
  if (!InsideLoop(p, loops.front())) return false;
  for (std::size_t i = 1; i < loops.size(); ++i) {
    if (InsideLoop(p, loops[i])) return false;
  }
  return true;
}
using Segment = std::array<int, 2>;
Segment Ordered(int a, int b) { return a < b ? Segment{a, b} : Segment{b, a}; }

Mesh Build(const FaceRecord& face, const Registry& registry,
           const Options& options, std::uint64_t& next) {
  const auto& f = *face.source;
  PlanarStraightLineGraph graph;
  std::vector<VertexRecord> canonical;
  std::vector<std::vector<Vec2>> loops;
  std::set<Segment> boundary;
  double uv_tolerance = 1e-8 * std::max({1., f.surface.domain[1] - f.surface.domain[0],
                                      f.surface.domain[3] - f.surface.domain[2]});
  for (const auto& loop : face.loops) {
    int offset = static_cast<int>(graph.uv.size());
    std::vector<Vec2> points;
    Vec2 previous_end{};
    std::uint64_t previous_vertex = 0;
    for (std::size_t use_index = 0; use_index < loop.size(); ++use_index) {
      const auto& use = loop[use_index];
      const auto& edge = registry.at(std::abs(use.signed_edge));
      auto ordered = edge.samples;
      if (use.signed_edge < 0) std::reverse(ordered.begin(), ordered.end());
      double previous_q = -1.;
      for (std::size_t i = 0; i < ordered.size(); ++i) {
        const auto& sample = ordered[i];
        double q;
        if (i == 0 || i + 1 == ordered.size()) {
          q = i == 0 ? 0. : 1.;
        } else {
          auto request = std::find_if(use.requests.begin(), use.requests.end(),
                                     [&](auto r) { return std::abs(r.first - sample.parameter) <= kParameterMerge; });
          q = request != use.requests.end()
                  ? request->second
                  : Project([&](double t) { return f.surface.EvaluateControl(UV(use, t)); },
                            sample.xyz, Tolerance(f, *edge.source, options));
        }
        if (q <= previous_q) {
          throw Unsupported("Topology edge-use samples are not strictly ordered");
        }
        previous_q = q;
        auto uv = UV(use, q);
        for (int axis = 0; axis < 2; ++axis) {
          for (double knot : face.axes[axis]) {
            if (std::abs(uv[axis] - knot) < uv_tolerance) uv[axis] = knot;
          }
        }
        if (Distance2(f.surface.EvaluateControl(uv), sample.xyz) >
            std::pow(Tolerance(f, *edge.source, options), 2)) {
          throw Unsupported("Canonical topology sample lies outside its face tolerance");
        }
        if (i == 0 && use_index &&
            (sample.vertex != previous_vertex || Norm(Subtract(uv, previous_end)) > uv_tolerance)) {
          throw Unsupported("Topology wire has a disconnected corner or UV gap");
        }
        if (i + 1 == ordered.size()) {
          previous_end = uv;
          previous_vertex = sample.vertex;
        } else {
          points.push_back(uv);
          canonical.push_back({sample.xyz, sample.vertex});
        }
      }
    }
    if (points.size() < 3 || Norm(Subtract(points.front(), previous_end)) > uv_tolerance ||
        canonical.at(offset).vertex != previous_vertex) {
      throw Unsupported("Topology wire is not a closed nondegenerate UV loop");
    }
    // No proximity welding: even coincident UVs require the same explicit
    // occurrence connectivity. Touching/self-intersecting loops are rejected.
    std::set<Vec2> unique(points.begin(), points.end());
    if (unique.size() != points.size()) {
      throw Unsupported("Topology wire has repeated UV vertices");
    }
    graph.uv.insert(graph.uv.end(), points.begin(), points.end());
    for (int i = 0; i < static_cast<int>(points.size()); ++i) {
      int a = offset + i, b = offset + (i + 1) % static_cast<int>(points.size());
      graph.segments.push_back({a, b});
      graph.markers.push_back(1);
      boundary.insert(Ordered(a, b));
    }
    if (!loops.empty()) {
      auto polygon = Polygon(points);
      const auto& t = polygon.triangles.front();
      graph.holes.push_back(Multiply(Add(Add(polygon.uv[t[0]], polygon.uv[t[1]]), polygon.uv[t[2]]), 1. / 3.));
    }
    loops.push_back(std::move(points));
  }
  std::size_t boundary_points = graph.uv.size();
  // Add knot intersections once. These are local geometric vertices, not
  // topological edge samples. Constraint intersections with the boundary were
  // already requested globally, before any face is triangulated.
  for (double u : face.axes[0]) {
    for (double v : face.axes[1]) {
      Vec2 uv{u, v};
      if (Inside(uv, loops) && std::find(graph.uv.begin(), graph.uv.end(), uv) == graph.uv.end()) {
        graph.uv.push_back(uv);
      }
    }
  }
  std::set<Segment> constraints = boundary;
  for (int axis = 0; axis < 2; ++axis) {
    for (double knot : face.axes[axis]) {
      std::vector<int> on_line;
      for (int i = 0; i < static_cast<int>(graph.uv.size()); ++i) {
        if (graph.uv[i][axis] == knot) on_line.push_back(i);
      }
      std::sort(on_line.begin(), on_line.end(), [&](int a, int b) {
        return graph.uv[a][1 - axis] < graph.uv[b][1 - axis];
      });
      for (std::size_t i = 1; i < on_line.size(); ++i) {
        int a = on_line[i - 1], b = on_line[i];
        auto midpoint = Multiply(Add(graph.uv[a], graph.uv[b]), .5);
        auto segment = Ordered(a, b);
        if (Inside(midpoint, loops) && constraints.insert(segment).second) {
          graph.segments.push_back(segment);
          graph.markers.push_back(2);
          ++graph.internal_constraints;
        }
      }
    }
  }
  ArcMap arc(f.surface);
  for (auto& uv : graph.uv) uv = arc.Apply(uv);
  for (auto& hole : graph.holes) hole = arc.Apply(hole);
  // Y prohibits local boundary splitting, while internal knot constraints may
  // refine. This is safe only after registry union and is verified below.
  const std::size_t original_points = graph.uv.size();
  const std::size_t insertion_budget = options.steiner_limit
                                          ? options.steiner_limit : 100000;
  std::string switches = "pq20zuQYS" + std::to_string(insertion_budget) +
      (options.steiner_limit ? "" : "L" + std::to_string(20 * (graph.internal_constraints + 1))) + "C";
  auto needs_refinement = [&](const std::array<Vec2, 3>& input) {
    auto uv = input;
    for (auto& p : uv) p = arc.Apply(p, true);
    return Refine(f.surface, options, uv);
  };
  auto result = TriangleRun(graph, options, needs_refinement, switches.c_str());
  // Triangle's circumcenter can encroach on an immutable boundary and be
  // discarded even with an unspent Steiner budget. Interior centroids are
  // safe: retriangulate the same PSLG, retaining every previous vertex. No
  // boundary samples are added, removed, or locally reparameterized here.
  for (;;) {
    if (result.uv.size() < graph.uv.size()) {
      throw std::runtime_error("Topology refinement lost input vertices");
    }
    for (std::size_t i = 0; i < graph.uv.size(); ++i) {
      if (result.uv[i] != graph.uv[i]) {
        throw std::runtime_error("Topology refinement changed an input vertex");
      }
    }
    std::vector<Vec2> insertions;
    std::set<Vec2> existing(result.uv.begin(), result.uv.end());
    bool incomplete = false;
    for (auto triangle : result.triangles) {
      std::array<Vec2, 3> points{result.uv[triangle[0]], result.uv[triangle[1]], result.uv[triangle[2]]};
      if (!needs_refinement(points)) continue;
      incomplete = true;
      auto centre = Multiply(Add(Add(points[0], points[1]), points[2]), 1. / 3.);
      if (existing.insert(centre).second) insertions.push_back(centre);
    }
    if (!incomplete) break;
    if (insertions.empty() || result.uv.size() + insertions.size() > original_points + insertion_budget) {
      throw std::runtime_error("Topology refinement budget exhausted or fixed boundary prevents requested refinement");
    }
    graph.uv = result.uv;
    graph.uv.insert(graph.uv.end(), insertions.begin(), insertions.end());
    result = TriangleRun(graph, options, {}, "pQYzS0C");
  }
  if (result.uv.size() < graph.uv.size() || result.triangles.empty()) {
    throw std::runtime_error("Topology triangulation lost input vertices or all triangles");
  }
  for (std::size_t i = 0; i < graph.uv.size(); ++i) {
    if (result.uv[i] != graph.uv[i]) {
      throw std::runtime_error("Topology triangulation changed an input vertex");
    }
  }
  std::map<Segment, int> incidence;
  for (auto triangle : result.triangles) {
    for (int i = 0; i < 3; ++i) {
      ++incidence[Ordered(triangle[i], triangle[(i + 1) % 3])];
    }
  }
  for (auto segment : boundary) {
    if (incidence[segment] != 1) {
      throw Unsupported("Topology boundary was split, removed, or is self-intersecting at segment " +
                        std::to_string(segment[0]) + ":" + std::to_string(segment[1]) +
                        " (incidence " + std::to_string(incidence[segment]) + ")");
    }
  }
  for (auto item : incidence) {
    if (item.second == 1 && !boundary.count(item.first)) {
      throw Unsupported("Topology triangulation introduced an unregistered boundary");
    }
    if (item.second > 2) throw Unsupported("Nonmanifold face triangulation");
  }
  for (auto& uv : result.uv) {
    uv = arc.Apply(uv, true);
    for (int axis = 0; axis < 2; ++axis) {
      for (double knot : face.axes[axis]) {
        double roundoff = 64 * kDoubleEpsilon * std::max({1., std::abs(uv[axis]), std::abs(knot)});
        if (std::abs(uv[axis] - knot) <= roundoff) uv[axis] = knot;
      }
    }
  }
  // Restore original boundary UV rather than retain inverse-map roundoff.
  std::size_t index = 0;
  for (const auto& loop : loops) for (auto uv : loop) result.uv[index++] = uv;
  for (auto triangle : result.triangles) {
    std::array<Vec2, 3> uv{result.uv[triangle[0]], result.uv[triangle[1]], result.uv[triangle[2]]};
    if (Refine(f.surface, options, uv)) {
      throw std::runtime_error("Topology refinement budget exhausted or fixed boundary prevents requested refinement");
    }
  }
  for (std::size_t i = boundary_points; i < result.uv.size(); ++i) {
    canonical.push_back({f.surface.EvaluateControl(result.uv[i]), next++});
  }
  for (auto triangle : result.triangles) {
    const auto& a = canonical.at(triangle[0]);
    const auto& b = canonical.at(triangle[1]);
    const auto& c = canonical.at(triangle[2]);
    auto area = Cross(Subtract(b.xyz, a.xyz), Subtract(c.xyz, a.xyz));
    if (a.vertex == b.vertex || b.vertex == c.vertex || c.vertex == a.vertex ||
        !Finite(area) || !(Norm(area) > 0.)) {
      throw Unsupported("Topology triangulation contains a collapsed geometric facet");
    }
  }
  Mesh out;
  out.face_id = f.id;
  out.identity_mode = "topology_knot_blocks";
  out.branch = "topology_consistent";
  out.explicit_weights = f.surface.explicit_weights;
  auto blocks = KnotBlocks(f.surface, false, true);
  std::vector<int> coverage(result.triangles.size());
  for (const auto& block : blocks) {
    std::map<std::uint32_t, std::uint32_t> remap;
    std::vector<Vec2> uv;
    std::vector<Triangle> triangles;
    for (std::size_t i = 0; i < result.triangles.size(); ++i) {
      auto triangle = result.triangles[i];
      auto centre = Multiply(Add(Add(result.uv[triangle[0]], result.uv[triangle[1]]), result.uv[triangle[2]]), 1. / 3.);
      if (centre[0] <= block[0] || centre[0] >= block[1] || centre[1] <= block[2] || centre[1] >= block[3]) continue;
      for (auto vertex : triangle) {
        const auto& point = result.uv[vertex];
        if (point[0] < block[0] - uv_tolerance || point[0] > block[1] + uv_tolerance ||
            point[1] < block[2] - uv_tolerance || point[1] > block[3] + uv_tolerance) {
          throw Unsupported("Topology triangle crosses an unresolved knot constraint");
        }
      }
      ++coverage[i];
      for (auto& vertex : triangle) {
        auto found = remap.find(vertex);
        if (found == remap.end()) {
          auto id = static_cast<std::uint32_t>(uv.size());
          uv.push_back(result.uv[vertex]);
          found = remap.emplace(vertex, id).first;
        }
        vertex = found->second;
      }
      triangles.push_back(triangle);
    }
    if (triangles.empty()) continue;
    Mesh part;
    part.uv = uv;
    part.normals = MeshNormals(f.surface, uv, &block);
    part.triangles = std::move(triangles);
    part.xyz.resize(uv.size());
    part.geometric_vertices.resize(uv.size());
    part.identities.assign(uv.size(), f.id + "|" + BlockKey(block));
    for (auto entry : remap) {
      part.xyz[entry.second] = canonical.at(entry.first).xyz;
      part.geometric_vertices[entry.second] = canonical.at(entry.first).vertex;
    }
    // Keep full precision in the opt-in path. Quantizing independent UVs or
    // interiors can collapse valid short edges without updating connectivity.
    for (auto normal : part.normals) {
      if (!Finite(normal)) throw Unsupported("Topology mode cannot resolve a singular surface normal");
    }
    Append(out, std::move(part));
  }
  if (std::any_of(coverage.begin(), coverage.end(), [](int n) { return n != 1; })) {
    throw Unsupported("Topology triangles cross an unresolved surface knot constraint");
  }
  return out;
}
}  // namespace

std::vector<Mesh> GenerateTopology(const Model& model,
                                  const std::vector<const Face*>& selected,
                                  const Options& options) {
  Registry registry;
  std::vector<FaceRecord> faces;
  std::map<int, const ImportedEdge*> imported;
  for (const auto& edge : model.edges) {
    if (edge.id <= 0 || !imported.emplace(edge.id, &edge).second) {
      throw Unsupported("Invalid or duplicate topological edge ID");
    }
  }
  // IDs use topology identity only, never proximity/coordinate equality.
  std::uint64_t next = model.vertex_positions.size() + 1;
  for (const auto* f : selected) {
    f->surface.Validate();
    if (!f->outer || f->boundary_edges.size() != f->inner.size() + 1) {
      throw Unsupported("Topology mode requires explicit unadapted OCCT wires");
    }
    FaceRecord face{f, {}, {}};
    for (const auto& block : KnotBlocks(f->surface, false, true)) {
      for (int axis = 0; axis < 2; ++axis) {
        for (int side = 0; side < 2; ++side) {
          double value = block[2 * axis + side];
          if (value > f->surface.domain[2 * axis] && value < f->surface.domain[2 * axis + 1]) {
            face.axes[axis].push_back(value);
          }
        }
      }
    }
    for (auto& axis : face.axes) axis = Unique(axis);
    for (std::size_t i = 0; i < f->boundary_edges.size(); ++i) {
      const auto& curves = i == 0 ? *f->outer : f->inner[i - 1];
      const auto& ids = f->boundary_edges[i];
      if (curves.size() != ids.size() || curves.empty()) {
        throw Unsupported("Topology edge-use provenance does not match wire curves");
      }
      std::vector<EdgeUse> loop;
      for (std::size_t j = 0; j < curves.size(); ++j) {
        auto found = imported.find(std::abs(ids[j]));
        if (ids[j] == 0 || found == imported.end()) {
          throw Unsupported("Topology mode cannot use synthetic/unidentified boundary edges");
        }
        const auto& edge = *found->second;
        if (edge.degenerate) throw Unsupported("Topology mode does not yet support degenerate pole edges");
        if (!edge.point || !(edge.domain[1] > edge.domain[0])) {
          throw Unsupported("Topology edge has no evaluable bounded XYZ curve");
        }
        if (!(curves[j].domain[1] > curves[j].domain[0])) {
          throw Unsupported("Topology edge use has an invalid pcurve interval");
        }
        for (int id : edge.vertices) {
          if (id <= 0 || static_cast<std::size_t>(id) > model.vertex_positions.size() ||
              !Finite(model.vertex_positions[id - 1])) {
            throw Unsupported("Topology edge has no canonical endpoint vertex");
          }
        }
        if (!registry.count(edge.id)) {
          auto sampled = Sample([&](double t) { return XYZ(edge, t); }, 0., 1.,
                                std::min(options.curve.chord, options.surface.chord),
                                std::min(options.curve.min_edge, options.surface.min_edge),
                                std::min(options.curve.max_edge, options.surface.max_edge), false, 8);
          for (double knot : edge.knots) {
            if (knot > edge.domain[0] && knot < edge.domain[1]) {
              sampled.parameters.push_back((knot - edge.domain[0]) / (edge.domain[1] - edge.domain[0]));
            }
          }
          registry.emplace(edge.id, EdgeRecord{&edge, std::move(sampled.parameters), {}});
        }
        loop.push_back({&curves[j], ids[j], {}});
      }
      face.loops.push_back(std::move(loop));
    }
    try {
      Collect(face, registry, options);
    } catch (const Unsupported& error) {
      throw Unsupported("Topology face " + f->id + ": " + error.what());
    }
    faces.push_back(std::move(face));
  }
  for (auto& item : registry) {
    auto& record = item.second;
    const auto& edge = *record.source;
    auto parameters = Merge(std::move(record.requests));
    // The legacy sampler's maximum is a half-chord heuristic. A locked
    // topological segment longer than the surface maximum cannot be fixed by
    // inserting interior points. Subdivide such requests globally, including
    // canonical vertex endpoint offsets, before building any face.
    auto position = [&](double t) {
      if (t == 0. || t == 1.) return model.vertex_positions.at(edge.vertices[t == 1.] - 1);
      return XYZ(edge, t);
    };
    std::vector<double> bounded;
    std::function<void(double, double, int)> split = [&](double a, double b, int depth) {
      if (bounded.size() >= 100000 || depth > 40) {
        throw std::runtime_error("Topology boundary subdivision work budget exhausted");
      }
      if (Norm(Subtract(position(b), position(a))) > .9 * options.surface.max_edge) {
        double midpoint = .5 * (a + b);
        if (midpoint == a || midpoint == b) {
          throw Unsupported("Topology edge parameter precision cannot resolve boundary subdivision");
        }
        split(a, midpoint, depth + 1);
        split(midpoint, b, depth + 1);
      } else {
        bounded.push_back(a);
      }
    };
    for (std::size_t i = 1; i < parameters.size(); ++i) split(parameters[i - 1], parameters[i], 0);
    bounded.push_back(1.);
    for (double t : bounded) {
      bool endpoint = t == 0. || t == 1.;
      auto id = endpoint ? static_cast<std::uint64_t>(edge.vertices[t == 1.]) : next++;
      record.samples.push_back({t, endpoint ? model.vertex_positions.at(id - 1) : XYZ(edge, t), id});
    }
  }
  std::vector<Mesh> out;
  for (const auto& face : faces) {
    try {
      out.push_back(Build(face, registry, options, next));
    } catch (const Unsupported& error) {
      throw Unsupported("Topology face " + face.source->id + ": " + error.what());
    } catch (const std::exception& error) {
      throw std::runtime_error("Topology face " + face.source->id + ": " + error.what());
    }
  }
  return out;
}
}  // namespace iko::internal
