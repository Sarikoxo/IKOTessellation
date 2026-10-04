#include "tessellation/private/internal.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <functional>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#include "tessellation/tessellation.h"

extern "C" {
#define REAL double
#define VOID void
#include "triangle.h"
}

namespace iko::internal {
extern "C" {
extern int triangle_auto_stats[7];
int TriangleExecute(char*, triangulateio*, triangulateio*);
}
namespace {
std::mutex triangle_mutex;
const std::function<bool(const std::array<iko::Vec2, 3>&)>* active_callback =
    nullptr;
std::exception_ptr callback_error;
void Release(triangulateio& m) {
  free(m.pointlist);
  free(m.pointattributelist);
  free(m.pointmarkerlist);
  free(m.trianglelist);
  free(m.triangleattributelist);
  free(m.trianglearealist);
  free(m.neighborlist);
  free(m.segmentlist);
  free(m.segmentmarkerlist);
  free(m.edgelist);
  free(m.edgemarkerlist);
  free(m.normlist);
}
}  // namespace
extern "C" int TriangleUnsuitable(double* a, double* b, double* c, double) {
  if (!active_callback || !*active_callback) {
    return 0;
  }
  try {
    return (*active_callback)({iko::Vec2{a[0], a[1]}, iko::Vec2{b[0], b[1]},
                               iko::Vec2{c[0], c[1]}})
               ? 1
               : 0;
  } catch (...) {
    if (!callback_error) {
      callback_error = std::current_exception();
    }
    return 0;
  }
}
BackendMesh TriangleRun(
    const PlanarStraightLineGraph& p, const Options& o,
    const std::function<bool(const std::array<Vec2, 3>&)>& cb,
    const char* switches) {
  if (p.uv.size() < 3 ||
      p.uv.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument(
        "Triangle requires at least three finite 2D points");
  }
  for (auto parameter_point : p.uv) {
    if (!Finite(parameter_point)) {
      throw std::invalid_argument("Nonfinite PSLG point");
    }
  }
  for (auto segment : p.segments) {
    for (int vertex_index : segment) {
      if (vertex_index < 0 ||
          static_cast<std::size_t>(vertex_index) >= p.uv.size()) {
        throw std::invalid_argument("Invalid PSLG segment");
      }
    }
  }
  if (!p.markers.empty() && p.markers.size() != p.segments.size()) {
    throw std::invalid_argument("Segment marker count mismatch");
  }
  std::string option =
      switches ? switches
      : p.disable_refinement
          ? "QS0pz"
          : "pq20zuQ" +
                (o.steiner_limit
                     ? "S" + std::to_string(o.steiner_limit)
                     : "L" +
                           std::to_string(20 * (p.internal_constraints + 1))) +
                "C";
  std::vector<char> command(option.begin(), option.end());
  command.push_back(0);
  triangulateio input{};
  triangulateio out{};
  input.pointlist = const_cast<double*>(p.uv.front().data());
  input.numberofpoints = static_cast<int>(p.uv.size());
  input.numberofsegments = static_cast<int>(p.segments.size());
  input.segmentlist = p.segments.empty()
                          ? nullptr
                          : const_cast<int*>(p.segments.front().data());
  input.numberofholes = static_cast<int>(p.holes.size());
  input.holelist =
      p.holes.empty() ? nullptr : const_cast<double*>(p.holes.front().data());
  input.segmentmarkerlist =
      p.markers.empty() ? nullptr : const_cast<int*>(p.markers.data());
  std::lock_guard<std::mutex> guard(triangle_mutex);
  active_callback = &cb;
  callback_error = nullptr;
  try {
    int failed = TriangleExecute(command.data(), &input, &out);
    active_callback = nullptr;
    if (failed) {
      throw std::runtime_error("Triangle computation failed (backend status " +
                               std::to_string(failed - 1) + ")");
    }
    BackendMesh result;
    for (int point_index = 0; point_index < out.numberofpoints; ++point_index) {
      result.uv.push_back(
          {out.pointlist[2 * point_index], out.pointlist[2 * point_index + 1]});
    }
    for (int triangle_index = 0; triangle_index < out.numberoftriangles;
         ++triangle_index) {
      result.triangles.push_back(
          {static_cast<std::uint32_t>(out.trianglelist[3 * triangle_index]),
           static_cast<std::uint32_t>(out.trianglelist[3 * triangle_index + 1]),
           static_cast<std::uint32_t>(
               out.trianglelist[3 * triangle_index + 2])});
    }
    if (out.segmentlist) {
      for (int segment_index = 0; segment_index < out.numberofsegments;
           ++segment_index) {
        result.segments.push_back({out.segmentlist[2 * segment_index],
                                   out.segmentlist[2 * segment_index + 1]});
        if (out.segmentmarkerlist) {
          result.markers.push_back(out.segmentmarkerlist[segment_index]);
        }
      }
    }
    std::copy(triangle_auto_stats, triangle_auto_stats + 7,
              result.auto_stats.begin());
    Release(out);
    out = {};
    if (callback_error) {
      std::rethrow_exception(callback_error);
    }
    return result;
  } catch (...) {
    active_callback = nullptr;
    Release(out);
    throw;
  }
}
BackendMesh Polygon(const std::vector<Vec2>& loop) {
  PlanarStraightLineGraph p;
  p.uv = loop;
  PrepareLoop(p.uv);
  for (int vertex_index = 0; vertex_index < static_cast<int>(loop.size());
       ++vertex_index) {
    p.segments.push_back(
        {vertex_index, (vertex_index + 1) % static_cast<int>(loop.size())});
  }
  auto result = TriangleRun(p, Options{}, {}, "QS0pz");
  if (result.triangles.empty()) {
    throw std::invalid_argument("Hole polygon has no interior triangles");
  }
  return result;
}
}  // namespace iko::internal
