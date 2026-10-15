#ifndef IKO_SRC_INTERNAL_H_
#define IKO_SRC_INTERNAL_H_

#include <algorithm>
#include <array>
#include <cmath>  // IWYU pragma: keep
#include <cstddef>
#include <filesystem>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "tessellation/tessellation.h"

namespace iko::internal {
using Rect = std::array<double, 4>;
constexpr double kDoubleEpsilon = std::numeric_limits<double>::epsilon();
constexpr double kFloatEpsilon = std::numeric_limits<float>::epsilon();
struct Unsupported : std::runtime_error {
  using std::runtime_error::runtime_error;
};
template <std::size_t dimensions>
std::array<double, dimensions> Add(const std::array<double, dimensions>& a,
                                   const std::array<double, dimensions>& b) {
  auto c = a;
  for (std::size_t coordinate_index = 0; coordinate_index < dimensions;
       ++coordinate_index) {
    c[coordinate_index] += b[coordinate_index];
  }
  return c;
}
template <std::size_t dimensions>
std::array<double, dimensions> Subtract(
    const std::array<double, dimensions>& a,
    const std::array<double, dimensions>& b) {
  auto c = a;
  for (std::size_t coordinate_index = 0; coordinate_index < dimensions;
       ++coordinate_index) {
    c[coordinate_index] -= b[coordinate_index];
  }
  return c;
}
template <std::size_t dimensions>
std::array<double, dimensions> Multiply(const std::array<double, dimensions>& a,
                                        double t) {
  auto c = a;
  for (auto& coordinate : c) {
    coordinate *= t;
  }
  return c;
}
template <std::size_t dimensions>
double Dot(const std::array<double, dimensions>& a,
           const std::array<double, dimensions>& b) {
  double s = 0;
  for (std::size_t coordinate_index = 0; coordinate_index < dimensions;
       ++coordinate_index) {
    s += a[coordinate_index] * b[coordinate_index];
  }
  return s;
}
inline Vec3 Cross(Vec3 a, Vec3 b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}
template <std::size_t dimensions>
double Norm(const std::array<double, dimensions>& a) {
  return std::sqrt(Dot(a, a));
}
template <std::size_t dimensions>
bool Finite(const std::array<double, dimensions>& a) {
  return std::all_of(a.begin(), a.end(), [](double coordinate) {
    return std::isfinite(coordinate);
  });
}
inline double Quantize(double x) {
  return static_cast<double>(static_cast<float>(x));
}
inline std::vector<double> Unique(std::vector<double> x) {
  std::sort(x.begin(), x.end());
  x.erase(std::unique(x.begin(), x.end()), x.end());
  return x;
}
std::string Number(double x);
std::string BlockKey(Rect r);
std::vector<double> Basis(double x, const std::vector<double>& knots,
                          int degree, bool local = false,
                          bool extrapolate = false);
struct Surface {
  int p = 1;
  int q = 1;
  int nu = 0;
  int nv = 0;
  std::vector<double> ku;
  std::vector<double> kv;
  std::vector<double> weights;
  std::vector<Vec3> control;
  Rect domain{};
  std::string label;
  bool explicit_weights = false;
  std::array<bool, 2> closed_axes{};
  void Validate() const;
  double BboxSpan() const;
  Vec3 Point(Vec2 uv) const;
  Vec3 Evaluate(Vec2 uv) const;
  Vec3 EvaluateControl(Vec2 uv, bool extrapolate = false) const;
  std::array<Vec3, 2> Derivatives(
      Vec2 uv, const std::array<std::vector<double>, 2>* basis_knots = nullptr,
      bool extrapolate = false) const;
};
struct Curve {
  int degree = 1;
  bool prepared = false;
  bool analytic_line =
      false;  // OCCT line conversion; needed for endpoint editing semantics
  bool linear_conversion = false;
  std::optional<Rect> import_bounds;
  std::vector<double> knots;
  std::vector<double> weights;
  std::vector<Vec3> control;
  std::array<double, 2> domain{};
  // Input-derived conversion, invoked only if an analytic endpoint is edited.
  std::function<Curve()> endpoint_edit_conversion;
  Vec3 Evaluate(double t, bool observed = false) const;
};
struct Face {
  std::string id;
  Surface surface;
  std::optional<std::vector<Curve>> outer;
  std::vector<std::vector<Curve>> inner;
  int source_de = 0;
  bool constructed_rectangle = false;
  bool reversed = false;
  int solid_id = 0;
  int shell_id = 0;
  std::vector<int> source_entries;
  std::vector<std::vector<int>> boundary_edges;
  std::string support_type;
  std::array<bool, 2> geometrically_closed_axes{};
  // Source declarations differ from geometric/periodic closure after OCCT.
  bool source_natural_boundary = false;
  std::array<bool, 2> source_closed_axes{};
};
struct ImportedEdge {
  int id = 0;
  int source_de = 0;
  std::array<int, 2> vertices{};
  bool degenerate = false;
  bool same_parameter = false;
  bool seam = false;
  double tolerance = 0.;
  std::string curve_type;
  std::array<double, 2> domain{};
  // XYZ curve knot breaks, in the same parameter range.
  std::vector<double> knots;
  // World-space, forward topological edge parameterization.
  std::function<Vec3(double)> point;
};
struct ImportCoverage {
  int source_de = 0;
  int source_type = 0;
  std::vector<int> faces;
};
struct Model {
  std::vector<Face> faces;
  std::vector<std::string> repair_trace;
  std::vector<ImportedEdge> edges;
  // OCCT vertex ID minus one indexes a canonical world-space position.
  std::vector<Vec3> vertex_positions;
  std::vector<ImportCoverage> roots;
  int shell_count = 0;
  int solid_count = 0;
  int vertex_count = 0;
  std::string source_units;
  double source_scale = 1.;
  double unit_to_mm = 1.;
  double source_resolution = 0.;
};
// OCCT import and shape healing, then adaptation; no legacy file repair or
// fallback.
Model ReadOcct(const std::filesystem::path& path);
Model ReadOcctRaw(const std::filesystem::path& path,
                  bool topology = false);  // explicit unadapted edge uses
void AdaptOcctFace(Face& face, std::vector<std::string>& trace);
void JoinOcctConversionEndpoints(std::vector<Curve>& curves,
                                 std::vector<int>& edge_ids,
                                 std::vector<std::string>& trace,
                                 const std::string& identity);
std::string ContentSha256(const std::string& bytes);
void RepairTrimWire(Face& face, std::vector<std::string>& trace);
void ValidateProvidedTrim(const Face& face);
Surface ClampActiveDomain(const Surface& source);
void ValidateNaturalRectangle(const Face& face);
Vec3 RectangleSamplingPoint(const Surface& surface, Vec2 uv);
Options Effective(Options options);
std::vector<Rect> KnotBlocks(const Surface& s, bool insert_midpoint = true,
                             bool c0_only = true);
std::vector<Vec3> MeshNormals(const Surface& s, const std::vector<Vec2>& uv,
                              const Rect* rect = nullptr);
Vec3 AnalyticNormal(const Surface& s, Vec2 uv, bool extrapolate = false);
struct Sampling {
  std::vector<double> parameters;
  std::vector<Vec3> points;
};
Sampling Sample(const std::function<Vec3(double)>& point, double start,
                double end, double chord, double minimum, double maximum,
                bool relative = false, int initial_segments = 1,
                const std::function<Vec3(double)>& initial_point = {});
struct AxisMap {
  std::vector<double> source;
  std::vector<double> target;
  double Apply(double x, bool inverse = false) const;
};
struct ArcMap {
  std::array<AxisMap, 2> maps;
  explicit ArcMap(const Surface& s);
  Vec2 Apply(Vec2 uv, bool inverse = false) const {
    return {maps[0].Apply(uv[0], inverse), maps[1].Apply(uv[1], inverse)};
  }
};
struct PlanarStraightLineGraph {
  std::vector<Vec2> uv;
  std::vector<Vec2> holes;
  std::vector<std::array<int, 2>> segments;
  std::vector<int> markers;
  bool disable_refinement = false;
  int internal_constraints = 0;
};
struct BackendMesh {
  std::vector<Vec2> uv;
  std::vector<Triangle> triangles;
  std::vector<std::array<int, 2>> segments;
  std::vector<int> markers;
  std::array<int, 7> auto_stats{};
};
BackendMesh TriangleRun(
    const PlanarStraightLineGraph& input, const Options& options,
    const std::function<bool(const std::array<Vec2, 3>&)>& callback = {},
    const char* switches = nullptr);
BackendMesh Polygon(const std::vector<Vec2>& loop);
bool PrepareLoop(std::vector<Vec2>& loop);
bool Refine(const Surface& s, const Options& options,
            const std::array<Vec2, 3>& uv);
Mesh Generate(const Face& face, const Options& effective_options);
std::vector<Mesh> GenerateTopology(
    const Model& model, const std::vector<const Face*>& selected,
    const Options& options);
Mesh BuildTrim(const Face& face, const Options& options);
Mesh BuildPole(const Face& face, const Options& options, int pole);
struct PoleAxes {
  std::vector<double> angular;
  std::vector<double> radial;
  std::vector<double> lengths;
  double angle = 0.;
};
PoleAxes ComputePoleAxes(const Surface& surface, const Options& options,
                         int pole);
std::optional<int> DetectPole(const Surface& s);
void Append(Mesh& destination, Mesh source);
void ComponentIdentities(Mesh& mesh, const std::string& face_id);
std::optional<std::array<double, 2>> AlignedDomain(
    const std::vector<double>& knots, int degree, std::array<double, 2> domain);
}  // namespace iko::internal

#endif  // IKO_SRC_INTERNAL_H_
