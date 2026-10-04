#include <iostream>

#include <tessellation/tessellation.h>

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: tessellation_example model.igs\n";
    return 2;
  }
  iko::Options options;
  options.surface = {0.05, 5.0, 0.1};
  options.curve = {0.05, 5.0, 0.1};
  auto result = iko::Tessellate(argv[1], options);
  if (!result.Succeeded() || !result.complete) {
    std::cerr << result.message << '\n';
    return static_cast<int>(result.status);
  }
  std::cout << result.faces.size() << "/" << result.model.required_faces
            << " complete faces, " << result.model.solids << " solids, "
            << result.mesh.xyz.size() << " vertices, "
            << result.mesh.triangles.size() << " triangles\n";
  // Every vertex has XYZ, UV, normal and identity; keep its original index.
  return 0;
}
