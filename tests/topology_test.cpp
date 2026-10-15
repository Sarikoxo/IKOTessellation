#include "tessellation/private/internal.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using iko::Mesh;
using iko::Options;
using iko::Vec2;
using iko::Vec3;
using iko::internal::Curve;
using iko::internal::Face;
using iko::internal::ImportedEdge;
using iko::internal::Model;
using iko::internal::Rect;
using iko::internal::Surface;
using Id = std::uint64_t;
using Segment = std::array<Id, 2>;

void Require(bool value, const std::string& message) {
  if (!value) {
    throw std::runtime_error(message);
  }
}

bool Near(double a, double b, double tolerance = 2e-6) {
  return std::abs(a - b) <= tolerance;
}

bool Near(Vec3 a, Vec3 b) {
  return Near(a[0], b[0]) && Near(a[1], b[1]) && Near(a[2], b[2]);
}

Vec3 Lerp(Vec3 a, Vec3 b, double t) {
  return iko::internal::Add(a, iko::internal::Multiply(
                                   iko::internal::Subtract(b, a), t));
}

Curve Pcurve(Vec2 a, Vec2 b, std::array<double, 2> domain,
             bool nonlinear = false) {
  Curve result;
  result.degree = nonlinear ? 2 : 1;
  result.domain = domain;
  result.knots.assign(result.degree + 1, domain[0]);
  result.knots.insert(result.knots.end(), result.degree + 1, domain[1]);
  result.control.push_back({a[0], a[1], 0.});
  if (nonlinear) {
    // Still a straight geometric edge, but t does not interpolate distance.
    result.control.push_back({a[0] + .2 * (b[0] - a[0]),
                              a[1] + .2 * (b[1] - a[1]), 0.});
  }
  result.control.push_back({b[0], b[1], 0.});
  result.weights.assign(result.control.size(), 1.);
  return result;
}

// These fixtures use no OCCT objects or files. Topological sharing is explicit:
// equal coordinates alone never cause the builder to reuse a vertex or edge.
struct Fixture {
  Model model;

  int Vertex(Vec3 point) {
    model.vertex_positions.push_back(point);
    model.vertex_count = static_cast<int>(model.vertex_positions.size());
    return model.vertex_count;
  }

  int Edge(int first, int second) {
    for (const auto& edge : model.edges) {
      if (edge.vertices == std::array<int, 2>{first, second}) {
        return edge.id;
      }
      if (edge.vertices == std::array<int, 2>{second, first}) {
        return -edge.id;
      }
    }
    ImportedEdge edge;
    edge.id = static_cast<int>(model.edges.size()) + 1;
    edge.vertices = {first, second};
    edge.domain = {2., 7.};
    edge.tolerance = 1e-8;
    edge.curve_type = "fixture_line";
    edge.degenerate = first == second;
    const auto a = model.vertex_positions.at(first - 1);
    const auto b = model.vertex_positions.at(second - 1);
    edge.point = [a, b](double t) { return Lerp(a, b, (t - 2.) / 5.); };
    model.edges.push_back(std::move(edge));
    return model.edges.back().id;
  }

  void Rectangle(const std::string& name, std::array<int, 4> vertices,
                 Rect domain = {0., 1., 0., 1.},
                 std::vector<double> v_breaks = {},
                 int nonlinear_side = -1) {
    Face face;
    face.id = name;
    face.source_de = 101 + 2 * static_cast<int>(model.faces.size());
    face.support_type = "fixture_bilinear";
    auto& surface = face.surface;
    surface.domain = domain;
    surface.nu = 2;
    surface.ku = {domain[0], domain[0], domain[1], domain[1]};
    std::vector<double> rows{domain[2]};
    rows.insert(rows.end(), v_breaks.begin(), v_breaks.end());
    rows.push_back(domain[3]);
    surface.nv = static_cast<int>(rows.size());
    surface.kv = {domain[2]};
    surface.kv.insert(surface.kv.end(), rows.begin(), rows.end());
    surface.kv.push_back(domain[3]);
    const auto a = model.vertex_positions.at(vertices[0] - 1);
    const auto b = model.vertex_positions.at(vertices[1] - 1);
    const auto c = model.vertex_positions.at(vertices[2] - 1);
    const auto d = model.vertex_positions.at(vertices[3] - 1);
    for (double v : rows) {
      double fraction = (v - domain[2]) / (domain[3] - domain[2]);
      surface.control.push_back(Lerp(a, d, fraction));
      surface.control.push_back(Lerp(b, c, fraction));
    }
    surface.weights.assign(surface.control.size(), 1.);
    const std::array<Vec2, 4> corners{{
        {domain[0], domain[2]}, {domain[1], domain[2]},
        {domain[1], domain[3]}, {domain[0], domain[3]}}};
    face.outer.emplace();
    face.boundary_edges.emplace_back();
    for (int side = 0; side < 4; ++side) {
      // Distinct pcurve parameter ranges exercise geometric inversion rather
      // than accidental equality with the canonical edge's [2, 7] range.
      face.outer->push_back(Pcurve(corners[side], corners[(side + 1) % 4],
                                    {11. + side, 19. + side},
                                    side == nonlinear_side));
      face.boundary_edges.front().push_back(
          Edge(vertices[side], vertices[(side + 1) % 4]));
    }
    surface.Validate();
    model.faces.push_back(std::move(face));
  }
};

Fixture Adjacent(bool knots = false, bool crease = false,
                 bool nonlinear = false) {
  Fixture fixture;
  int a = fixture.Vertex({0., 0., 0.});
  int b = fixture.Vertex({1., 0., 0.});
  int c = fixture.Vertex({1., 1., 0.});
  int d = fixture.Vertex({0., 1., 0.});
  int e = fixture.Vertex(crease ? Vec3{1., 0., 1.} : Vec3{2., 0., 0.});
  int f = fixture.Vertex(crease ? Vec3{1., 1., 1.} : Vec3{2., 1., 0.});
  fixture.Rectangle("left", {a, b, c, d}, {0., 1., 0., 1.},
                    knots ? std::vector<double>{.3, .65}
                          : std::vector<double>{});
  fixture.Rectangle("right", {b, e, f, c}, {10., 14., -3., 5.}, {},
                    nonlinear ? 3 : -1);
  return fixture;
}

Options TopologyOptions() {
  Options options;
  options.topology_consistent = true;
  options.curve.min_edge = .005;
  options.curve.max_edge = .28;
  options.curve.chord = .001;
  options.surface.min_edge = .005;
  options.surface.max_edge = .6;
  options.surface.chord = .001;
  return iko::internal::Effective(options);
}

std::vector<Mesh> Generate(const Model& model,
                           const Options& options = TopologyOptions()) {
  std::vector<const Face*> selected;
  for (const auto& face : model.faces) {
    selected.push_back(&face);
  }
  auto result = iko::internal::GenerateTopology(model, selected, options);
  Require(result.size() == model.faces.size(), "Selected face count changed");
  std::map<Id, Vec3> coordinates;
  for (const auto& mesh : result) {
    Require(!mesh.triangles.empty(), "Topology face has no triangles");
    Require(mesh.xyz.size() == mesh.uv.size() &&
                mesh.xyz.size() == mesh.normals.size() &&
                mesh.xyz.size() == mesh.geometric_vertices.size(),
            "Topology record arrays are not aligned");
    for (std::size_t i = 0; i < mesh.xyz.size(); ++i) {
      Require(iko::internal::Finite(mesh.xyz[i]) &&
                  iko::internal::Finite(mesh.uv[i]) &&
                  iko::internal::Finite(mesh.normals[i]),
              "Topology records are nonfinite");
      auto inserted = coordinates.emplace(mesh.geometric_vertices[i],
                                           mesh.xyz[i]);
      Require(inserted.second || inserted.first->second == mesh.xyz[i],
              "One geometric ID has differing stored XYZ coordinates");
    }
    for (const auto& triangle : mesh.triangles) {
      for (auto index : triangle) {
        Require(index < mesh.xyz.size(), "Triangle index exceeds records");
      }
      Require(triangle[0] != triangle[1] && triangle[0] != triangle[2] &&
                  triangle[1] != triangle[2],
              "Topology generated repeated-index triangle");
      const auto& ids = mesh.geometric_vertices;
      Require(ids[triangle[0]] != ids[triangle[1]] &&
                  ids[triangle[0]] != ids[triangle[2]] &&
                  ids[triangle[1]] != ids[triangle[2]],
              "Triangle repeats one geometric vertex across distinct records");
      auto a = iko::internal::Subtract(mesh.xyz[triangle[1]],
                                       mesh.xyz[triangle[0]]);
      auto b = iko::internal::Subtract(mesh.xyz[triangle[2]],
                                       mesh.xyz[triangle[0]]);
      Require(iko::internal::Norm(iko::internal::Cross(a, b)) > 0.,
              "Regular topology fixture contains a zero-area XYZ triangle");
    }
  }
  return result;
}

bool OnSharedEdge(Vec3 point) {
  return Near(point[0], 1.) && Near(point[2], 0.) &&
         point[1] >= -2e-6 && point[1] <= 1. + 2e-6;
}

std::map<Id, std::size_t> SharedRecords(const Mesh& mesh) {
  std::map<Id, std::size_t> result;
  for (std::size_t i = 0; i < mesh.xyz.size(); ++i) {
    if (OnSharedEdge(mesh.xyz[i])) {
      result.emplace(mesh.geometric_vertices[i], i);
    }
  }
  return result;
}

std::map<Segment, int> SharedSegments(const Mesh& mesh) {
  std::map<Segment, int> result;
  for (const auto& triangle : mesh.triangles) {
    for (int side = 0; side < 3; ++side) {
      auto a = triangle[side];
      auto b = triangle[(side + 1) % 3];
      if (OnSharedEdge(mesh.xyz[a]) && OnSharedEdge(mesh.xyz[b])) {
        Id first = mesh.geometric_vertices[a];
        Id second = mesh.geometric_vertices[b];
        Require(first != second, "Shared-edge triangle is collapsed");
        result[{std::min(first, second), std::max(first, second)}] +=
            first < second ? 1 : -1;
      }
    }
  }
  return result;
}

void CheckAdjacent(const std::vector<Mesh>& meshes, bool crease = false) {
  Require(meshes.size() == 2, "Adjacent fixture needs two faces");
  auto first = SharedRecords(meshes[0]);
  auto second = SharedRecords(meshes[1]);
  Require(first.size() > 2 && first.size() == second.size(),
          "Shared edge does not have matching refined sample counts");
  for (const auto& record : first) {
    auto found = second.find(record.first);
    Require(found != second.end(), "Shared edge has a face-local geometric ID");
    const auto& a = meshes[0];
    const auto& b = meshes[1];
    Require(a.xyz[record.second] == b.xyz[found->second],
            "Shared edge XYZ is not exactly equal");
    Require(a.uv[record.second] != b.uv[found->second],
            "Independent UV charts were welded");
    if (crease) {
      Require(Near(a.normals[record.second], {0., 0., 1.}) &&
                  Near(b.normals[found->second], {-1., 0., 0.}),
              "Sharp-edge face normals were lost or averaged");
    }
  }
  auto a = SharedSegments(meshes[0]);
  auto b = SharedSegments(meshes[1]);
  Require(!a.empty() && a.size() == b.size(),
          "Shared constrained-segment counts differ");
  Require(a.size() + 1 == first.size(),
          "Shared edge contains a hanging or unused geometric sample");
  for (const auto& segment : a) {
    auto found = b.find(segment.first);
    Require(std::abs(segment.second) == 1 && found != b.end() &&
                found->second == -segment.second,
            "Incident faces do not use opposite, matching edge segments");
  }
}

Id At(const Mesh& mesh, Vec3 point) {
  bool found = false;
  Id result = 0;
  for (std::size_t i = 0; i < mesh.xyz.size(); ++i) {
    if (Near(mesh.xyz[i], point)) {
      Require(!found || result == mesh.geometric_vertices[i],
              "Duplicate corner records disagree on geometric identity");
      found = true;
      result = mesh.geometric_vertices[i];
    }
  }
  Require(found, "Expected constrained point is absent");
  return result;
}

void MismatchedChartsAndDirection() {
  auto fixture = Adjacent();
  Require(fixture.model.faces[0].boundary_edges[0][1] ==
              -fixture.model.faces[1].boundary_edges[0][3],
          "Fixture must exercise opposite edge-use directions");
  CheckAdjacent(Generate(fixture.model));
}

void NonlinearPcurveParameterization() {
  auto fixture = Adjacent(false, false, true);
  CheckAdjacent(Generate(fixture.model));
}

void UnilateralKnotSplits() {
  auto fixture = Adjacent(true);
  auto meshes = Generate(fixture.model);
  CheckAdjacent(meshes);
  for (double y : {.3, .65}) {
    Require(At(meshes[0], {1., y, 0.}) == At(meshes[1], {1., y, 0.}),
            "A knot intersection was not propagated to the incident face");
  }
}

void NearEndpointKnotSplits() {
  auto fixture = Adjacent();
  int shared = fixture.model.faces.front().boundary_edges[0][1];
  fixture.model.edges.at(shared - 1).knots = {2.005, 6.995};
  auto check = [&] {
    auto meshes = Generate(fixture.model);
    CheckAdjacent(meshes);
    for (double y : {.001, .999}) {
      Require(At(meshes[0], {1., y, 0.}) == At(meshes[1], {1., y, 0.}),
              "Near-endpoint knot intersection was lost during parameter mapping");
    }
  };
  // First require the canonical XYZ curve's knots alone to propagate. Then
  // repeat with matching knots present in only one incident surface.
  check();
  auto& surface = fixture.model.faces.front().surface;
  surface.nv = 4;
  surface.kv = {0., 0., .001, .999, 1., 1.};
  surface.control = {{0., 0., 0.}, {1., 0., 0.},
                      {0., .001, 0.}, {1., .001, 0.},
                      {0., .999, 0.}, {1., .999, 0.},
                      {0., 1., 0.}, {1., 1., 0.}};
  surface.weights.assign(surface.control.size(), 1.);
  check();
}

void SharpEdgeNormals() {
  auto fixture = Adjacent(false, true);
  CheckAdjacent(Generate(fixture.model), true);
}

void CommonCorner() {
  Fixture fixture;
  int common = fixture.Vertex({0., 0., 0.});
  fixture.Rectangle("xy", {common, fixture.Vertex({1., 0., 0.}),
                            fixture.Vertex({1., 1., 0.}),
                            fixture.Vertex({0., 1., 0.})});
  fixture.Rectangle("xz", {common, fixture.Vertex({-1., 0., 0.}),
                            fixture.Vertex({-1., 0., 1.}),
                            fixture.Vertex({0., 0., 1.})});
  fixture.Rectangle("yz", {common, fixture.Vertex({0., -1., 0.}),
                            fixture.Vertex({0., -1., -1.}),
                            fixture.Vertex({0., 0., -1.})});
  auto meshes = Generate(fixture.model);
  auto id = At(meshes.front(), {0., 0., 0.});
  for (const auto& mesh : meshes) {
    Require(At(mesh, {0., 0., 0.}) == id,
            "Incident edges did not share their topological corner vertex");
  }
}

void UnrelatedCoincidentFaces() {
  Fixture fixture;
  for (int face = 0; face < 2; ++face) {
    fixture.Rectangle("coincident_" + std::to_string(face),
                      {fixture.Vertex({0., 0., 0.}),
                       fixture.Vertex({1., 0., 0.}),
                       fixture.Vertex({1., 1., 0.}),
                       fixture.Vertex({0., 1., 0.})});
  }
  auto meshes = Generate(fixture.model);
  std::set<Id> first(meshes[0].geometric_vertices.begin(),
                     meshes[0].geometric_vertices.end());
  for (Id id : meshes[1].geometric_vertices) {
    Require(first.count(id) == 0,
            "Unrelated coincident geometry was position-welded");
  }
  Require(At(meshes[0], {0., 0., 0.}) != At(meshes[1], {0., 0., 0.}),
          "Unrelated coincident corners share a geometric ID");
}

Model CylinderModel() {
  Model model;
  model.vertex_positions = {{1., 0., 0.}, {1., 0., 1.}};
  model.vertex_count = 2;
  Face face;
  face.id = "closed_cylinder";
  face.surface.p = 2;
  face.surface.nu = 9;
  face.surface.nv = 2;
  face.surface.ku = {0., 0., 0., .25, .25, .5, .5,
                     .75, .75, 1., 1., 1.};
  face.surface.kv = {0., 0., 1., 1.};
  face.surface.domain = {0., 1., 0., 1.};
  face.surface.explicit_weights = true;
  face.surface.closed_axes = {true, false};
  face.geometrically_closed_axes = {true, false};
  const std::array<Vec2, 9> circle{{
      {1., 0.}, {1., 1.}, {0., 1.}, {-1., 1.}, {-1., 0.},
      {-1., -1.}, {0., -1.}, {1., -1.}, {1., 0.}}};
  for (double z : {0., 1.}) {
    for (std::size_t i = 0; i < circle.size(); ++i) {
      face.surface.control.push_back({circle[i][0], circle[i][1], z});
      face.surface.weights.push_back(i % 2 ? std::sqrt(.5) : 1.);
    }
  }
  face.surface.Validate();
  face.outer = std::vector<Curve>{
      Pcurve({0., 0.}, {1., 0.}, {0., 1.}),
      Pcurve({1., 0.}, {1., 1.}, {0., 1.}),
      Pcurve({1., 1.}, {0., 1.}, {0., 1.}),
      Pcurve({0., 1.}, {0., 0.}, {0., 1.})};
  face.boundary_edges = {{1, 2, -3, -2}};
  const auto surface = face.surface;
  for (int id = 1; id <= 3; ++id) {
    ImportedEdge edge;
    edge.id = id;
    edge.vertices = id == 1 ? std::array<int, 2>{1, 1}
                    : id == 2 ? std::array<int, 2>{1, 2}
                              : std::array<int, 2>{2, 2};
    edge.domain = {0., 1.};
    edge.tolerance = 1e-8;
    edge.seam = id == 2;
    // A closed circle is not a degenerate edge even though its endpoint
    // vertex is repeated. The two vertical uses share one imported edge.
    edge.point = [surface, id](double t) {
      return id == 2 ? Vec3{1., 0., t}
                     : surface.EvaluateControl({t, id == 1 ? 0. : 1.});
    };
    model.edges.push_back(std::move(edge));
  }
  model.faces.push_back(std::move(face));
  return model;
}

void RepeatedSeamPreservesOccurrences() {
  auto model = CylinderModel();
  auto meshes = Generate(model);
  const auto& mesh = meshes.front();
  std::array<std::map<Id, std::size_t>, 2> occurrences;
  std::array<std::map<Segment, int>, 2> segments;
  for (std::size_t i = 0; i < mesh.uv.size(); ++i) {
    for (int side = 0; side < 2; ++side) {
      if (Near(mesh.uv[i][0], side)) {
        occurrences[side].emplace(mesh.geometric_vertices[i], i);
      }
    }
  }
  Require(occurrences[0].size() > 2 &&
              occurrences[0].size() == occurrences[1].size(),
          "Repeated seam occurrences have different sample counts");
  for (const auto& item : occurrences[0]) {
    auto found = occurrences[1].find(item.first);
    Require(found != occurrences[1].end() && item.second != found->second &&
                mesh.xyz[item.second] == mesh.xyz[found->second],
            "Repeated seam lost separate records or exact shared geometry");
  }
  for (const auto& triangle : mesh.triangles) {
    bool at_start = false;
    bool at_end = false;
    for (auto index : triangle) {
      at_start |= Near(mesh.uv[index][0], 0.);
      at_end |= Near(mesh.uv[index][0], 1.);
    }
    Require(!(at_start && at_end),
            "A triangle bridges distinct seam occurrences in parameter space");
    for (int edge = 0; edge < 3; ++edge) {
      auto a = triangle[edge];
      auto b = triangle[(edge + 1) % 3];
      for (int side = 0; side < 2; ++side) {
        if (Near(mesh.uv[a][0], side) && Near(mesh.uv[b][0], side)) {
          Id first = mesh.geometric_vertices[a];
          Id second = mesh.geometric_vertices[b];
          Require(first != second, "Cylinder seam segment is collapsed");
          segments[side][{std::min(first, second), std::max(first, second)}] +=
              first < second ? 1 : -1;
        }
      }
    }
  }
  Require(segments[0].size() + 1 == occurrences[0].size() &&
              segments[0].size() == segments[1].size(),
          "Repeated seam occurrences have unmatched constrained segments");
  for (const auto& segment : segments[0]) {
    auto found = segments[1].find(segment.first);
    Require(std::abs(segment.second) == 1 && found != segments[1].end() &&
                found->second == -segment.second,
            "Repeated seam directions or segment identities differ");
  }
}

void SharedCurvedEdge() {
  auto model = CylinderModel();
  model.vertex_positions.push_back({1., 0., 2.});
  model.vertex_count = 3;
  auto face = model.faces.front();
  face.id = "upper_cylinder";
  for (auto& point : face.surface.control) {
    point[2] += 1.;
  }
  face.boundary_edges = {{3, 4, -5, -4}};
  const auto surface = face.surface;
  for (int id : {4, 5}) {
    ImportedEdge edge;
    edge.id = id;
    edge.vertices = id == 4 ? std::array<int, 2>{2, 3}
                           : std::array<int, 2>{3, 3};
    edge.domain = {0., 1.};
    edge.tolerance = 1e-8;
    edge.seam = id == 4;
    edge.point = [surface, id](double t) {
      return id == 4 ? Vec3{1., 0., 1. + t}
                     : surface.EvaluateControl({t, 1.});
    };
    model.edges.push_back(std::move(edge));
  }
  model.faces.push_back(std::move(face));
  auto meshes = Generate(model);
  std::array<std::map<Id, Vec3>, 2> points;
  std::array<std::map<Segment, int>, 2> segments;
  for (std::size_t face_index = 0; face_index < meshes.size(); ++face_index) {
    const auto& mesh = meshes[face_index];
    for (std::size_t i = 0; i < mesh.xyz.size(); ++i) {
      if (Near(mesh.xyz[i][2], 1.)) {
        points[face_index].emplace(mesh.geometric_vertices[i], mesh.xyz[i]);
      }
    }
    for (const auto& triangle : mesh.triangles) {
      for (int side = 0; side < 3; ++side) {
        auto a = triangle[side];
        auto b = triangle[(side + 1) % 3];
        if (Near(mesh.xyz[a][2], 1.) && Near(mesh.xyz[b][2], 1.)) {
          Id first = mesh.geometric_vertices[a];
          Id second = mesh.geometric_vertices[b];
          segments[face_index][{std::min(first, second),
                                 std::max(first, second)}] +=
              first < second ? 1 : -1;
        }
      }
    }
  }
  Require(points[0].size() > 4 && points[0] == points[1],
          "Shared curved edge does not have identical geometric IDs and XYZ");
  Require(segments[0].size() == points[0].size() &&
              segments[0].size() == segments[1].size(),
          "Shared closed curved edge is not a matching complete sample cycle");
  for (const auto& segment : segments[0]) {
    auto found = segments[1].find(segment.first);
    Require(std::abs(segment.second) == 1 && found != segments[1].end() &&
                found->second == -segment.second,
            "Shared curved edge has incompatible triangle connectivity");
  }
}

void ExplicitSteinerBudgetFailsAtomically() {
  auto model = CylinderModel();
  auto options = TopologyOptions();
  options.steiner_limit = 1;
  std::vector<Mesh> published;
  bool exhausted = false;
  try {
    published = iko::internal::GenerateTopology(
        model, {&model.faces.front()}, options);
  } catch (const iko::internal::Unsupported&) {
    throw std::runtime_error(
        "Steiner exhaustion incorrectly classified valid geometry as unsupported");
  } catch (const std::runtime_error& error) {
    exhausted = true;
    Require(std::string(error.what()).find("budget") != std::string::npos,
            "Steiner exhaustion did not provide a work-budget diagnostic");
  }
  Require(exhausted && published.empty(),
          "Exhausted Steiner budget silently published an incomplete mesh");
  // A failed request must leave the serialized Triangle backend ready for
  // another independent request with the normal automatic work allowance.
  auto recovered = Generate(model);
  Require(recovered.size() == 1 && !recovered.front().triangles.empty(),
          "Topology backend did not recover after work-budget exhaustion");
}

void IntraFaceKnotPreservesNormals() {
  Model model;
  model.vertex_positions = {{0., 0., 0.}, {1., 0., .5},
                            {1., 1., .5}, {0., 1., 0.}};
  model.vertex_count = 4;
  Face face;
  face.id = "folded_surface";
  face.surface.nu = 3;
  face.surface.nv = 2;
  face.surface.ku = {0., 0., .5, 1., 1.};
  face.surface.kv = {0., 0., 1., 1.};
  face.surface.domain = {0., 1., 0., 1.};
  face.surface.control = {{0., 0., 0.}, {.5, 0., 0.}, {1., 0., .5},
                          {0., 1., 0.}, {.5, 1., 0.}, {1., 1., .5}};
  face.surface.weights.assign(6, 1.);
  face.surface.Validate();
  const std::array<Vec2, 4> corners{{
      {0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}}};
  face.outer.emplace();
  face.boundary_edges = {{1, 2, 3, 4}};
  const auto surface = face.surface;
  for (int side = 0; side < 4; ++side) {
    auto a = corners[side];
    auto b = corners[(side + 1) % 4];
    face.outer->push_back(Pcurve(a, b, {0., 1.}));
    ImportedEdge edge;
    edge.id = side + 1;
    edge.vertices = {side + 1, (side + 1) % 4 + 1};
    edge.domain = {0., 1.};
    edge.tolerance = 1e-8;
    edge.point = [surface, a, b](double t) {
      return surface.EvaluateControl(
          {a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])});
    };
    model.edges.push_back(std::move(edge));
  }
  model.faces.push_back(std::move(face));
  auto meshes = Generate(model);
  const auto& mesh = meshes.front();
  std::map<Id, std::vector<Vec3>> normals;
  for (std::size_t i = 0; i < mesh.uv.size(); ++i) {
    if (Near(mesh.uv[i][0], .5)) {
      normals[mesh.geometric_vertices[i]].push_back(mesh.normals[i]);
    }
  }
  Require(normals.size() >= 2, "C0 constraint endpoints are missing");
  const double unit = std::sqrt(.5);
  for (const auto& item : normals) {
    bool left = false;
    bool right = false;
    for (auto normal : item.second) {
      left |= Near(normal, {0., 0., 1.});
      right |= Near(normal, {-unit, 0., unit});
    }
    Require(left && right,
            "C0 knot geometric vertex lacks its separate one-sided normals");
  }
  std::map<Segment, std::array<int, 2>> segments;
  for (const auto& triangle : mesh.triangles) {
    bool left = false;
    bool right = false;
    for (auto index : triangle) {
      left |= mesh.uv[index][0] < .5 - 2e-6;
      right |= mesh.uv[index][0] > .5 + 2e-6;
    }
    Require(!(left && right), "Triangle crosses the C0 knot constraint");
    for (int side = 0; side < 3; ++side) {
      auto a = triangle[side];
      auto b = triangle[(side + 1) % 3];
      if (Near(mesh.uv[a][0], .5) && Near(mesh.uv[b][0], .5)) {
        Id first = mesh.geometric_vertices[a];
        Id second = mesh.geometric_vertices[b];
        auto& counts = segments[{std::min(first, second),
                                 std::max(first, second)}];
        ++counts[0];
        counts[1] += first < second ? 1 : -1;
      }
    }
  }
  Require(segments.size() + 1 == normals.size(),
          "C0 constraint contains unmatched geometric samples");
  for (const auto& segment : segments) {
    Require(segment.second[0] == 2 && segment.second[1] == 0,
            "C0 knot block records have incompatible shared connectivity");
  }
}

template <class Function>
void MustBeUnsupported(Function function, const std::string& name) {
  try {
    function();
  } catch (const iko::internal::Unsupported& error) {
    Require(!std::string(error.what()).empty(), name + " has no diagnostic");
    return;
  }
  throw std::runtime_error(name + " did not explicitly report Unsupported");
}

void AmbiguousSeamIsRejected() {
  auto fixture = Adjacent();
  fixture.model.faces.resize(1);
  auto& face = fixture.model.faces.front();
  int east = face.boundary_edges[0][1];
  // A seam-like repeated use whose other pcurve lies on a different physical
  // edge cannot be reconciled. It must not silently produce a cracked mesh.
  fixture.model.edges.at(east - 1).seam = true;
  face.boundary_edges[0][3] = -east;
  MustBeUnsupported([&] { Generate(fixture.model); }, "Ambiguous seam");
}

void DegeneratePoleIsRejected() {
  Fixture fixture;
  int pole = fixture.Vertex({0., 0., 0.});
  fixture.Rectangle("pole", {pole, pole, fixture.Vertex({1., 1., 0.}),
                              fixture.Vertex({-1., 1., 0.})});
  MustBeUnsupported([&] { Generate(fixture.model); }, "Degenerate pole");
}

void LegacyGenerationIsUnchanged() {
  auto fixture = Adjacent();
  auto face = fixture.model.faces.front();
  // Exercise the ordinary legacy path, independent of topology metadata.
  face.outer.reset();
  face.boundary_edges.clear();
  Options defaults;
  Require(!defaults.topology_consistent, "Topology consistency is not opt-in");
  auto implicit = iko::internal::Generate(face, defaults);
  Options disabled = defaults;
  disabled.topology_consistent = false;
  auto explicit_off = iko::internal::Generate(face, disabled);
  Require(implicit.geometric_vertices.empty() &&
              explicit_off.geometric_vertices.empty(),
          "Legacy generation unexpectedly contains geometric identities");
  Require(implicit.xyz == explicit_off.xyz && implicit.uv == explicit_off.uv &&
              implicit.normals == explicit_off.normals &&
              implicit.triangles == explicit_off.triangles &&
              implicit.identities == explicit_off.identities &&
              implicit.unpaired_uv_tail == explicit_off.unpaired_uv_tail &&
              implicit.face_id == explicit_off.face_id &&
              implicit.identity_mode == explicit_off.identity_mode &&
              implicit.branch == explicit_off.branch &&
              implicit.explicit_weights == explicit_off.explicit_weights,
          "Explicitly disabling topology changed legacy records");
  Require(implicit.xyz.size() == 4 && implicit.triangles.size() == 2 &&
              implicit.identity_mode == "knot_blocks" &&
              implicit.branch == "ordinary",
          "Legacy bilinear-plane baseline changed");
}
}  // namespace

int main() {
  try {
    const std::vector<std::pair<const char*, void (*)()>> tests{
        {"mismatched UV charts and opposite edge uses", MismatchedChartsAndDirection},
        {"nonlinear pcurve parameterization", NonlinearPcurveParameterization},
        {"unilateral knot intersections", UnilateralKnotSplits},
        {"near-endpoint knot intersections", NearEndpointKnotSplits},
        {"sharp-edge normal records", SharpEdgeNormals},
        {"common topological corner", CommonCorner},
        {"unrelated coincident faces", UnrelatedCoincidentFaces},
        {"valid repeated seam occurrences", RepeatedSeamPreservesOccurrences},
        {"shared curved topological edge", SharedCurvedEdge},
        {"explicit Steiner budget and backend recovery", ExplicitSteinerBudgetFailsAtomically},
        {"intra-face C0 normal records", IntraFaceKnotPreservesNormals},
        {"ambiguous seam rejection", AmbiguousSeamIsRejected},
        {"degenerate pole rejection", DegeneratePoleIsRejected},
        {"legacy default compatibility", LegacyGenerationIsUnchanged}};
    std::size_t failures = 0;
    for (const auto& test : tests) {
      std::cout << "Checking " << test.first << '\n' << std::flush;
      try {
        test.second();
      } catch (const std::exception& error) {
        ++failures;
        std::cerr << "FAILED " << test.first << ": " << error.what() << '\n';
      }
    }
    std::cout << tests.size() - failures << '/' << tests.size()
              << " topology cases passed\n";
    return failures ? 1 : 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
