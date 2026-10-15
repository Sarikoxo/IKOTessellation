#include "tessellation/tessellation.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>

namespace {
void Require(bool condition, const std::string& message) {
  if (!condition) throw std::runtime_error(message);
}
void VerifySolid(const std::filesystem::path& path, bool require_seam) {
  iko::Options options;
  options.topology_consistent = true;
  auto result = iko::Tessellate(path, options);
  Require(result.Succeeded() && result.complete && result.whole_model,
          path.filename().string() + ": " + result.message);
  Require(result.requested_options.topology_consistent &&
              result.effective_options.topology_consistent,
          "Topology option was not retained");
  const auto& mesh = result.mesh;
  Require(mesh.xyz.size() == mesh.geometric_vertices.size() &&
              mesh.xyz.size() == mesh.uv.size() &&
              mesh.xyz.size() == mesh.normals.size(),
          "Imported topology record arrays are not aligned");
  std::map<std::uint64_t, iko::Vec3> positions;
  std::map<std::uint64_t, std::set<std::size_t>> faces;
  std::map<std::array<std::uint64_t, 2>, int> incidence;
  bool seam = false;
  for (std::size_t f = 0; f < result.faces.size(); ++f) {
    const auto& range = result.faces[f];
    Require(range.vertex_count && range.triangle_count,
            "An imported face was silently omitted");
    std::map<std::uint64_t, iko::Vec2> local_uv;
    for (std::size_t i = range.vertex_offset;
         i < range.vertex_offset + range.vertex_count; ++i) {
      auto id = mesh.geometric_vertices.at(i);
      auto point = positions.emplace(id, mesh.xyz.at(i));
      Require(point.second || point.first->second == mesh.xyz[i],
              "Shared imported vertex has different XYZ records");
      faces[id].insert(f);
      auto uv = local_uv.emplace(id, mesh.uv[i]);
      seam |= !uv.second && uv.first->second != mesh.uv[i];
    }
  }
  for (auto triangle : mesh.triangles) {
    for (int i = 0; i < 3; ++i) {
      auto a = mesh.geometric_vertices.at(triangle[i]);
      auto b = mesh.geometric_vertices.at(triangle[(i + 1) % 3]);
      Require(a != b, "Imported solid has a collapsed geometric facet");
      ++incidence[{std::min(a, b), std::max(a, b)}];
    }
  }
  Require(std::any_of(faces.begin(), faces.end(),
                     [](const auto& item) { return item.second.size() > 1; }),
          "Imported solid did not share any geometric vertices");
  for (const auto& edge : incidence) {
    Require(edge.second == 2, "Imported closed solid has an unmatched/nonmanifold edge");
  }
  Require(!require_seam || seam, "Imported cylinder lost separate seam UV occurrences");
  std::cout << "Verified " << path.filename().string() << ": "
            << result.faces.size() << " faces, " << mesh.triangles.size()
            << " triangles, closed geometric incidence\n";
}
}  // namespace

int main(int argc, char** argv) {
  try {
    Require(argc == 2, "topology import test requires fixture directory");
    std::filesystem::path root = argv[1];
    VerifySolid(root / "prismatic_brep.igs", false);
    VerifySolid(root / "topology_cylinder.igs", true);
    iko::Options limited;
    limited.topology_consistent = true;
    limited.surface.chord = limited.curve.chord = .001;
    limited.steiner_limit = 1;
    auto failed = iko::Tessellate(root / "topology_cylinder.igs", limited);
    Require(failed.status == iko::Status::kComputationFailed && !failed.complete &&
                failed.mesh.xyz.empty() && failed.mesh.triangles.empty() && failed.faces.empty(),
            "Failed topology refinement published a partial model");
    auto legacy = iko::Tessellate(root / "two_planes.igs");
    Require(legacy.Succeeded() && legacy.complete && legacy.mesh.geometric_vertices.empty(),
            "Legacy import changed or backend failed to recover");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
