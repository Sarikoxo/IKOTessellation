#include "tessellation/tessellation.h"

#include <algorithm>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <utility>
#if defined(_MSC_VER) && defined(_DEBUG)
#include <crtdbg.h>
#endif

namespace {
void Require(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}
}  // namespace

int main(int argc, char** argv) {
  try {
    Require(argc == 3, "Expected input and output directories");
    const auto input = std::filesystem::u8path(argv[1]) / "two_planes.igs";
    const auto out = std::filesystem::u8path(argv[2]);
    std::filesystem::create_directories(out);
    auto listed = iko::ListSurfaces(input);
    Require(listed.status == iko::Status::kSuccess && listed.surfaces.size() == 2,
            "Installed DLL surface listing");
    // DLL-allocated strings, STL debug proxies, nested vectors and mesh arrays
    // must support ordinary caller-side copy/move/reallocation/destruction.
    for (int iteration = 0; iteration < 4; ++iteration) {
      auto result = iko::Tessellate(input);
      Require(result.Succeeded() && result.complete && result.faces.size() == 2 &&
                  result.model.required_faces == 2 && !result.mesh.triangles.empty(),
              "Installed DLL complete-model tessellation");
      auto copy = result;
      auto moved = std::move(copy);
      moved.message.append(4096, 'x');
      moved.mesh.xyz.reserve(moved.mesh.xyz.size() + 4096);
      for (auto& identity : moved.mesh.identities) identity.append(256, 'i');
      for (auto& face : moved.faces) {
        face.surface.instance_id.append(256, 'f');
        face.surface.source_entries.reserve(256);
        face.boundary_edges.emplace_back(256, 0);
      }
      moved.repair_trace.emplace_back(4096, 'r');
      for (auto& root : moved.model.roots) root.face_ids.emplace_back(4096, 's');
      auto surface_copy = listed;
      surface_copy.message.append(4096, 'l');
      surface_copy.surfaces.front().source_entries.reserve(256);
      Require(result.mesh.xyz.size() == moved.mesh.xyz.size(), "Value ownership");
      // Caller-allocated storage is consumed by the DLL and remains caller-owned.
      const auto obj = iko::WriteObj(result.mesh, out / "model.obj");
      const auto stl = iko::WriteBinaryStl(result.mesh, out / "model.stl");
      Require(obj.status == iko::Status::kSuccess &&
                  stl.status == iko::Status::kSuccess,
              "Installed DLL OBJ/STL export");
      iko::Mesh invalid = result.mesh;
      invalid.triangles.push_back({UINT32_MAX, 0, 1});
      auto rejected = iko::WriteObj(invalid, out / "invalid.obj");
      Require(rejected.status == iko::Status::kInputProblem, "Export error result");
      rejected.message.append(4096, 'e');
      auto missing = iko::Tessellate(out / "missing-input.igs");
      Require(missing.status == iko::Status::kInputProblem && !missing.complete,
              "Tessellation error result");
      missing.message.append(4096, 'e');
      Require(std::strcmp(iko::StatusName(result.status), "success") == 0,
              "Borrowed status name");
    }
#if defined(_MSC_VER) && defined(_DEBUG)
    Require(_CrtCheckMemory() != 0, "Shared Debug CRT heap consistency");
#endif
    std::cout << "Installed DLL API and STL ownership passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
