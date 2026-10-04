#include "tessellation/tessellation.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tessellation/private/internal.h"

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>  // IWYU pragma: keep

#include <winnls.h>
#endif
namespace iko {
namespace {
static const char kHelp[] = R"(Usage: tessellate model.igs [options]
  --min-edge X          Set curve and surface minimum edge (default 0.05)
  --max-edge X          Set curve and surface maximum edge (default 5)
  --chord X             Set curve and surface chord height (default 0.1)
  --surface-min X       Override surface minimum edge
  --surface-max X       Override surface maximum edge
  --surface-chord X     Override surface chord height
  --curve-min X         Override curve minimum edge
  --curve-max X         Override curve maximum edge
  --curve-chord X       Override curve chord height
  --relative-uv-chord X Relative UV chord tolerance (default 0.01)
  --surface DE          Diagnostic selection by IGES directory entry
  --list                List effective model faces
  --no-smooth-poles     Disable the baseline polar-band path
  --steiner-limit N     Insertion budget; 0 uses automatic protection
  --out PATH            Optional .obj or binary .stl export
  --help                Show this help
Exit codes: 0 success, 2 input/arguments, 3 unsupported,
            4 computation failure, 5 export failure.
)";
static int Run(int argc, char** argv) {
  Options o;
  std::filesystem::path file;
  std::filesystem::path out;
  bool list = false;
  try {
    for (int argument_index = 1; argument_index < argc; ++argument_index) {
      std::string key = argv[argument_index];
      if (key == "--help" || key == "-h") {
        std::cout << kHelp;
        return 0;
      }
      if (key == "--list" || key == "--list-surfaces") {
        list = true;
        continue;
      }
      if (key == "--no-smooth-poles") {
        o.smooth_poles = false;
        continue;
      }
      if (key.rfind("--", 0) != 0) {
        if (!file.empty()) {
          throw std::invalid_argument("Only one geometry file is accepted");
        }
        file = std::filesystem::u8path(argv[argument_index]);
        continue;
      }
      if (argument_index + 1 >= argc) {
        throw std::invalid_argument("Missing value for " + key);
      }
      std::string value = argv[++argument_index];
      if (key == "--out") {
        out = std::filesystem::u8path(value);
        continue;
      }
      std::size_t used;
      double x = std::stod(value, &used);
      if (used != value.size() || !std::isfinite(x)) {
        throw std::invalid_argument("Invalid numeric option " + key);
      }
      if (key == "--min-edge") {
        o.surface.min_edge = o.curve.min_edge = x;
      } else if (key == "--max-edge") {
        o.surface.max_edge = o.curve.max_edge = x;
      } else if (key == "--chord") {
        o.surface.chord = o.curve.chord = x;
      } else if (key == "--surface-min") {
        o.surface.min_edge = x;
      } else if (key == "--surface-max") {
        o.surface.max_edge = x;
      } else if (key == "--surface-chord") {
        o.surface.chord = x;
      } else if (key == "--curve-min") {
        o.curve.min_edge = x;
      } else if (key == "--curve-max") {
        o.curve.max_edge = x;
      } else if (key == "--curve-chord") {
        o.curve.chord = x;
      } else if (key == "--relative-uv-chord") {
        o.relative_uv_chord = x;
      } else if (key == "--surface" || key == "--steiner-limit") {
        if (x != std::trunc(x) || x < 0 ||
            x > std::numeric_limits<std::int32_t>::max() ||
            (key == "--surface" && x == 0)) {
          throw std::invalid_argument("Invalid integer option " + key);
        }
        if (key == "--surface") {
          o.surface_de = static_cast<int>(x);
        } else {
          o.steiner_limit = static_cast<int>(x);
        }
      } else {
        throw std::invalid_argument("Unknown option " + key);
      }
    }
    if (file.empty()) {
      throw std::invalid_argument("Geometry file is required; use --help");
    }
    if (!out.empty()) {
      auto ext = out.extension().string();
      std::transform(
          ext.begin(), ext.end(), ext.begin(),
          [](unsigned char extension_character) {
            return static_cast<char>(std::tolower(extension_character));
          });
      if (ext != ".obj" && ext != ".stl") {
        throw std::invalid_argument("Output must end in .obj or .stl");
      }
    }
    internal::Effective(
        o);  // Share argument validation with SDK generation, including --list.
  } catch (const std::exception& e) {
    std::cerr << "input_problem: " << e.what() << '\n';
    return 2;
  }
  if (list) {
    auto r = ListSurfaces(file);
    if (r.status != Status::kSuccess) {
      std::cerr << StatusName(r.status) << ": " << r.message << '\n';
      return static_cast<int>(r.status);
    }
    for (const auto& surface : r.surfaces) {
      std::cout << "DE" << surface.directory_entry << " degree "
                << surface.degree[0] << ',' << surface.degree[1]
                << " trimmed=" << surface.trimmed << " weights="
                << (surface.explicit_weights ? "explicit" : "omitted")
                << " instance=" << surface.instance_id << '\n';
    }
    return 0;
  }
  auto r = Tessellate(file, o);
  if (!r.Succeeded()) {
    std::cerr << StatusName(r.status) << ": " << r.message << '\n';
    return static_cast<int>(r.status);
  }
  std::vector<bool> used(r.mesh.xyz.size());
  std::size_t degenerate = 0;
  for (auto triangle : r.mesh.triangles) {
    for (auto vertex_index : triangle) {
      used[vertex_index] = true;
    }
    auto a = r.mesh.xyz[triangle[0]];
    auto b = r.mesh.xyz[triangle[1]];
    auto c = r.mesh.xyz[triangle[2]];
    double u[3];
    double v[3];
    for (int coordinate_index = 0; coordinate_index < 3; ++coordinate_index) {
      u[coordinate_index] = b[coordinate_index] - a[coordinate_index];
      v[coordinate_index] = c[coordinate_index] - a[coordinate_index];
    }
    double n[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                   u[0] * v[1] - u[1] * v[0]};
    degenerate += (n[0] == 0 && n[1] == 0 && n[2] == 0);
  }
  std::cout << "success: " << r.faces.size() << " faces, " << r.mesh.xyz.size()
            << " vertices, " << r.mesh.triangles.size() << " triangles, "
            << std::count(used.begin(), used.end(), false)
            << " unused vertices, " << degenerate << " degenerate facets, "
            << r.mesh.unpaired_uv_tail.size() << " unpaired UV records\n"
            << r.message << '\n';
  if (!out.empty()) {
    auto ext = out.extension().string();
    std::transform(
        ext.begin(), ext.end(), ext.begin(),
        [](unsigned char extension_character) {
          return static_cast<char>(std::tolower(extension_character));
        });
    auto e =
        ext == ".obj" ? WriteObj(r.mesh, out) : WriteBinaryStl(r.mesh, out);
    if (e.status != Status::kSuccess) {
      std::cerr << "export_failed: " << e.message << '\n';
      return 5;
    }
    std::cout << e.message << ": " << out.u8string() << '\n';
  }
  return 0;
}
}  // namespace
}  // namespace iko

#ifdef _WIN32
int wmain(int argc, wchar_t** argv) {
  std::vector<std::string> arguments;
  for (int argument_index = 0; argument_index < argc; ++argument_index) {
    int count = WideCharToMultiByte(CP_UTF8, 0, argv[argument_index], -1,
                                    nullptr, 0, nullptr, nullptr);
    std::string value(count, '\0');
    WideCharToMultiByte(CP_UTF8, 0, argv[argument_index], -1, value.data(),
                        count, nullptr, nullptr);
    value.pop_back();
    arguments.push_back(std::move(value));
  }
  std::vector<char*> pointers;
  for (auto& argument : arguments) {
    pointers.push_back(argument.data());
  }
  return iko::Run(argc, pointers.data());
}
#else
int main(int argc, char** argv) {
  return iko::Run(argc, argv);
}
#endif
