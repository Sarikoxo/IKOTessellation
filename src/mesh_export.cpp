#include "tessellation/tessellation.h"

#include <cstdint>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <limits>
#include <ostream>
#include <stdexcept>

#include "tessellation/private/internal.h"

namespace iko {
static void ValidMesh(const Mesh& m) {
  if (m.xyz.size() != m.uv.size() || m.xyz.size() != m.normals.size() ||
      m.xyz.size() != m.identities.size()) {
    throw std::invalid_argument("Vertex fields mismatch");
  }
  for (auto vertex_position : m.xyz) {
    if (!internal::Finite(vertex_position)) {
      throw std::invalid_argument("Nonfinite vertex coordinates");
    }
  }
  for (auto triangle : m.triangles) {
    for (auto vertex_index : triangle) {
      if (vertex_index >= m.xyz.size()) {
        throw std::invalid_argument("Invalid triangle indices");
      }
    }
  }
}
ExportResult WriteObj(const Mesh& m, const std::filesystem::path& path) {
  try {
    ValidMesh(m);
    std::ofstream out(path);
    if (!out) {
      throw std::runtime_error("Cannot open OBJ output");
    }
    out << std::setprecision(17)
        << "# iko 0.1; original vertex and directed facet records\n";
    for (auto vertex_position : m.xyz) {
      out << "v " << vertex_position[0] << ' ' << vertex_position[1] << ' '
          << vertex_position[2] << '\n';
    }
    for (auto parameter_point : m.uv) {
      out << "vt " << parameter_point[0] << ' ' << parameter_point[1] << '\n';
    }
    for (auto vertex_normal : m.normals) {
      out << "vn " << vertex_normal[0] << ' ' << vertex_normal[1] << ' '
          << vertex_normal[2] << '\n';
    }
    for (auto triangle : m.triangles) {
      out << 'f';
      for (auto vertex_index : triangle) {
        auto n = static_cast<std::uint64_t>(vertex_index) + 1;
        out << ' ' << n << '/' << n << '/' << n;
      }
      out << '\n';
    }
    out.close();
    if (!out) {
      throw std::runtime_error("OBJ output failed");
    }
    return {Status::kSuccess, "OBJ exported"};
  } catch (const std::invalid_argument& e) {
    return {Status::kInputProblem, e.what()};
  } catch (const std::exception& e) {
    return {Status::kComputationFailed, e.what()};
  }
}
static void Word(std::ostream& out, std::uint32_t v, int bytes = 4) {
  for (int byte_index = 0; byte_index < bytes; ++byte_index) {
    out.put(static_cast<char>((v >> (8 * byte_index)) & 255));
  }
}
static void Scalar(std::ostream& out, double v) {
  float x = static_cast<float>(v);
  std::uint32_t bits;
  std::memcpy(&bits, &x, 4);
  Word(out, bits);
}
ExportResult WriteBinaryStl(const Mesh& m, const std::filesystem::path& path) {
  try {
    ValidMesh(m);
    if (m.triangles.size() > std::numeric_limits<std::uint32_t>::max()) {
      throw std::invalid_argument("STL facet limit exceeded");
    }
    std::ofstream out(path, std::ios::binary);
    if (!out) {
      throw std::runtime_error("Cannot open STL output");
    }
    char header[80] = {};
    const char label[] =
        "iko 0.1 binary STL; directed facets, including degeneracies";
    std::memcpy(header, label, sizeof(label) - 1);
    out.write(header, 80);
    Word(out, static_cast<std::uint32_t>(m.triangles.size()));
    for (auto triangle : m.triangles) {
      auto a = m.xyz[triangle[0]];
      auto b = m.xyz[triangle[1]];
      auto c = m.xyz[triangle[2]];
      auto n =
          internal::Cross(internal::Subtract(b, a), internal::Subtract(c, a));
      double length = internal::Norm(n);
      n = length > 0 ? internal::Multiply(n, 1. / length) : Vec3{};
      for (double normal_coordinate : n) {
        Scalar(out, normal_coordinate);
      }
      for (auto vertex_index : triangle) {
        for (double vertex_coordinate : m.xyz[vertex_index]) {
          Scalar(out, vertex_coordinate);
        }
      }
      Word(out, 0, 2);
    }
    out.close();
    if (!out) {
      throw std::runtime_error("STL output failed");
    }
    return {Status::kSuccess, "Binary STL exported"};
  } catch (const std::invalid_argument& e) {
    return {Status::kInputProblem, e.what()};
  } catch (const std::exception& e) {
    return {Status::kComputationFailed, e.what()};
  }
}
}  // namespace iko
