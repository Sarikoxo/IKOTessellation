#ifndef IKO_INCLUDE_IKO_TESSELLATION_H_
#define IKO_INCLUDE_IKO_TESSELLATION_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

// STL storage is owned on both sides of the DLL boundary. Use matching MSVC
// toolsets, dynamic CRTs, configurations and default STL ABI settings.
#if defined(_MSC_VER)
#if !defined(_DLL)
#error iko tessellation requires /MD (Release) or /MDd (Debug)
#endif
#if defined(IKO_TESSELLATION_MSC_VER) && _MSC_VER != IKO_TESSELLATION_MSC_VER
#error iko tessellation requires the MSVC toolset used to build the SDK
#endif
#if defined(IKO_TESSELLATION_ITERATOR_DEBUG_LEVEL) && \
    _ITERATOR_DEBUG_LEVEL != IKO_TESSELLATION_ITERATOR_DEBUG_LEVEL
#error iko tessellation requires matching Release/Debug and iterator debug level
#endif
#endif

#if defined(_WIN32)
#if defined(IKO_TESSELLATION_BUILD)
#define IKO_TESSELLATION_API __declspec(dllexport)
#else
#define IKO_TESSELLATION_API __declspec(dllimport)
#endif
#else
#define IKO_TESSELLATION_API __attribute__((visibility("default")))
#endif

namespace iko {
using Vec2 = std::array<double, 2>;
using Vec3 = std::array<double, 3>;
using Triangle = std::array<std::uint32_t, 3>;
enum class Status {
  kSuccess = 0,
  kInputProblem = 2,
  kUnsupported = 3,
  kComputationFailed = 4,
};
struct Tolerances {
  double min_edge = 0.05;
  double max_edge = 5.0;
  double chord = 0.1;
};
struct Options {
  Tolerances surface;
  Tolerances curve;
  double relative_uv_chord = 0.01;
  bool smooth_poles = true;
  // Zero uses the automatic work limit.
  std::int32_t steiner_limit = 0;
  // Zero generates the complete model. A positive IGES directory entry selects
  // one surface for diagnostics.
  int surface_de = 0;
};
struct Mesh {
  std::vector<Vec3> xyz;
  std::vector<Vec2> uv;
  std::vector<Vec3> normals;
  // Directed, zero-based triangle indices. Vertex records are never welded.
  std::vector<Triangle> triangles;
  std::vector<std::string> identities;
  std::vector<Vec2> unpaired_uv_tail;
  std::string face_id;
  std::string identity_mode;
  std::string branch;
  bool explicit_weights = false;
};
struct SurfaceInfo {
  int directory_entry = 0;
  std::string instance_id;
  std::string label;
  std::array<int, 2> degree{};
  std::array<double, 4> domain{};
  bool trimmed = false;
  bool explicit_weights = false;
  bool reversed = false;
  int solid_id = 0;
  int shell_id = 0;
  std::string support_type;
  std::vector<int> source_entries;
};
struct FaceRange {
  SurfaceInfo surface;
  std::size_t vertex_offset = 0;
  std::size_t vertex_count = 0;
  std::size_t triangle_offset = 0;
  std::size_t triangle_count = 0;
  std::size_t extra_uv_offset = 0;
  std::size_t extra_uv_count = 0;
  std::string identity_mode;
  std::string branch;
  // Signed OCCT edge IDs. Zero identifies an adaptation connector.
  std::vector<std::vector<int>> boundary_edges;
};
struct SourceRoot {
  int directory_entry = 0;
  int entity_type = 0;
  std::vector<std::string> face_ids;
};
struct ModelInfo {
  std::size_t imported_faces = 0;
  std::size_t required_faces = 0;
  int solids = 0;
  int shells = 0;
  int edges = 0;
  int vertices = 0;
  std::string source_units;
  double source_scale = 1.;
  double unit_to_mm = 1.;
  double source_resolution = 0.;
  std::vector<SourceRoot> roots;
};
struct Result {
  Status status = Status::kComputationFailed;
  std::string message;
  Mesh mesh;
  Options requested_options;
  Options effective_options;
  std::vector<std::string> repair_trace;
  std::vector<FaceRange> faces;
  ModelInfo model;
  bool complete = false;
  bool whole_model = true;
  bool Succeeded() const noexcept {
    return status == Status::kSuccess;
  }
};
struct SurfaceList {
  Status status = Status::kComputationFailed;
  std::string message;
  std::vector<SurfaceInfo> surfaces;
};
struct ExportResult {
  Status status = Status::kComputationFailed;
  std::string message;
};
// Imports IGES geometry and generates a mesh. Failures are returned in status
// and message; callers should check Succeeded() and complete before exporting.
IKO_TESSELLATION_API Result Tessellate(const std::filesystem::path& geometry_file,
                  const Options& options = {});
// Lists the faces produced by importing and adapting the geometry.
IKO_TESSELLATION_API SurfaceList ListSurfaces(const std::filesystem::path& geometry_file);
// Exports existing vertex and directed triangle records without welding them.
IKO_TESSELLATION_API ExportResult WriteObj(const Mesh& mesh, const std::filesystem::path& output);
IKO_TESSELLATION_API ExportResult WriteBinaryStl(const Mesh& mesh,
                            const std::filesystem::path& output);
// Borrowed immutable DLL storage; do not free it or use it after DLL unload.
IKO_TESSELLATION_API const char* StatusName(Status status) noexcept;
}  // namespace iko

#endif  // IKO_INCLUDE_IKO_TESSELLATION_H_
