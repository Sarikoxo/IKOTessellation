#include "tessellation/private/internal.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "tessellation/tessellation.h"

namespace iko::internal {
static constexpr double kQuadratureNodes[] = {-.93247,  -.6612094, -.2386192,
                                              .2386192, .6612094,  .93247};
static constexpr double kQuadratureWeights[] = {.1713245, .3607616, .4679139,
                                                .4679139, .3607616, .1713245};
static double Length(const Surface& s, int axis, double fixed, double lo,
                     double hi) {
  double a0 = std::min(lo, hi);
  double b0 = std::max(lo, hi);
  std::vector<double> breaks{a0, b0};
  for (double knot : axis ? s.kv : s.ku) {
    if (a0 < knot && knot < b0) {
      breaks.push_back(knot);
    }
  }
  breaks = Unique(breaks);
  double total = 0;
  for (std::size_t integration_span_index = 0;
       integration_span_index + 1 < breaks.size(); ++integration_span_index) {
    double a = breaks[integration_span_index];
    double b = breaks[integration_span_index + 1];
    if (a + 1e-8 > b) {
      continue;
    }
    double value = 0;
    for (int quadrature_index = 0; quadrature_index < 6; ++quadrature_index) {
      double t =
          (b - a) * .5 * kQuadratureNodes[quadrature_index] + (a + b) * .5;
      auto d = s.Derivatives(axis ? Vec2{fixed, t} : Vec2{t, fixed})[axis];
      value += std::sqrt((d[0] * d[0] + d[1] * d[1]) + d[2] * d[2]) *
               kQuadratureWeights[quadrature_index];
    }
    total += (b - a) * .5 * value;
  }
  return total;
}
std::optional<int> DetectPole(const Surface& s) {
  for (int pole : {0, 2, 1, 3}) {
    int axis = pole % 2;
    double fixed = s.domain[2 * (1 - axis) + (pole >= 2)];
    if (Length(s, axis, fixed, s.domain[2 * axis], s.domain[2 * axis + 1]) <
        .001) {
      return pole;
    }
  }
  return {};
}
static std::vector<double> Midpoints(const std::vector<double>& values) {
  std::vector<double> out;
  for (std::size_t parameter_index = 0; parameter_index < values.size();
       ++parameter_index) {
    out.push_back(values[parameter_index]);
    if (parameter_index + 1 < values.size()) {
      out.push_back((values[parameter_index] + values[parameter_index + 1]) *
                    .5);
    }
  }
  return out;
}
PoleAxes ComputePoleAxes(const Surface& s, const Options& o, int pole) {
  int angular_axis = pole % 2;
  int radial_axis = 1 - angular_axis;
  auto uv = [&](double a, double r) {
    return angular_axis ? Vec2{r, a} : Vec2{a, r};
  };
  auto blocks = KnotBlocks(s, false, false);
  std::array<std::vector<double>, 2> original;
  for (auto knot_block : blocks) {
    for (int parameter_axis = 0; parameter_axis < 2; ++parameter_axis) {
      for (int endpoint_index = 0; endpoint_index < 2; ++endpoint_index) {
        original[parameter_axis].push_back(
            knot_block[2 * parameter_axis + endpoint_index]);
      }
    }
  }
  for (auto& axis_parameters : original) {
    axis_parameters = Unique(axis_parameters);
  }
  auto angular = original[angular_axis];
  auto radial = original[radial_axis];
  if (pole >= 2) {
    std::reverse(radial.begin(), radial.end());
  }
  double origin_radial = original[radial_axis].front();
  double far_radial = original[radial_axis].back();
  auto origin = s.EvaluateControl(uv(angular.front(), origin_radial));
  auto first =
      Subtract(s.EvaluateControl(uv(angular.front(), far_radial)), origin);
  auto last =
      Subtract(s.EvaluateControl(uv(angular.back(), far_radial)), origin);
  double cosine = Dot(first, last) /
                  (std::sqrt(Dot(first, first)) * std::sqrt(Dot(last, last)));
  double angle =
      std::acos(std::clamp(cosine, -1., 1.)) * 180. / 3.14159265358979323846;
  if (angle >= 1e-8) {
    auto query =
        uv(angular.front() + .01, origin_radial + (pole >= 2 ? -.01 : .01));
    bool outside = query[0] < s.domain[0] || query[0] > s.domain[1] ||
                   query[1] < s.domain[2] || query[1] > s.domain[3];
    Vec3 n;
    if (outside) {
      if (query[angular_axis] < angular.front() ||
          query[angular_axis] > angular.back() ||
          query[radial_axis] > original[radial_axis].back()) {
        throw Unsupported(
            "Narrow parameter domain requires separate normal wrapping "
            "validation");
      }
      n = AnalyticNormal(s, query, true);
      if (!Finite(n)) {
        throw Unsupported(
            "Invalid extrapolated polar opening normal needs native fallback");
      }
    } else {
      n = MeshNormals(s, {query})[0];
    }
    double orientation = n[1] * last[0] * first[2];
    orientation -= n[0] * last[1] * first[2];
    orientation += n[0] * first[1] * last[2];
    orientation -= n[1] * first[0] * last[2];
    orientation -= first[1] * last[0] * n[2];
    orientation += first[0] * last[1] * n[2];
    if (orientation < 0) {
      angle = 360. - angle;
    }
  }
  if (angle < 1e-8) {
    angle = 360.;
  }
  while (angular.size() < angle / 20.) {
    angular = Midpoints(angular);
  }
  std::vector<double> lengths;
  for (double angular_parameter : angular) {
    lengths.push_back(Length(s, radial_axis, angular_parameter, radial.back(),
                             radial.front()));
  }
  auto chosen = std::max_element(lengths.begin(), lengths.end());
  double fixed = angular[chosen - lengths.begin()];
  double chord = std::min(std::max(.2 * o.surface.chord, .0001 * s.BboxSpan()),
                          o.curve.chord);
  double minimum = std::min(o.surface.min_edge, o.curve.min_edge);
  double maximum = std::max(o.surface.max_edge, .01 * s.BboxSpan());
  std::vector<double> values;
  for (std::size_t radial_span_index = 0; radial_span_index + 1 < radial.size();
       ++radial_span_index) {
    double a = radial[radial_span_index];
    double b = radial[radial_span_index + 1];
    auto point = [&](double t) {
      return s.Point(uv(fixed, (1. - t) * a + t * b));
    };
    auto initial = [&](double t) {
      auto p = uv(fixed, (1. - t) * a + t * b);
      return Vec3{p[0], p[1], 0};
    };
    auto r = Sample(point, 0., 1., chord, minimum, maximum, true, 1, initial);
    for (std::size_t sample_index = 0; sample_index + 1 < r.parameters.size();
         ++sample_index) {
      double t = r.parameters[sample_index];
      values.push_back((1. - t) * a + t * b);
    }
  }
  values.push_back(radial.back());
  double target =
      .75 * std::min(150, 2 * static_cast<int>(*chosen / maximum + 1.));
  while (values.size() < target) {
    values = Midpoints(values);
  }
  radial = std::move(values);
  return {angular, radial, lengths, angle};
}
Mesh BuildPole(const Face& face, const Options& o, int pole) {
  const auto& s = face.surface;
  int angular_axis = pole % 2;
  int radial_axis = 1 - angular_axis;
  auto uv = [&](double a, double r) {
    return angular_axis ? Vec2{r, a} : Vec2{a, r};
  };
  auto axes = ComputePoleAxes(s, o, pole);
  const auto& angular = axes.angular;
  const auto& radial = axes.radial;
  const auto& lengths = axes.lengths;
  Mesh out;
  out.face_id = face.id;
  out.identity_mode = "index_components";
  out.branch = "observed_polar_bands";
  out.explicit_weights = s.explicit_weights;
  for (std::size_t angular_cell_index = 0;
       angular_cell_index + 1 < angular.size(); ++angular_cell_index) {
    double a = angular[angular_cell_index];
    double b = angular[angular_cell_index + 1];
    double average =
        (lengths[angular_cell_index] + lengths[angular_cell_index + 1]) * .5;
    std::vector<int> counts;
    for (double radial_parameter : radial) {
      double arc = Length(s, angular_axis, radial_parameter, a, b);
      counts.push_back(
          arc < 1e-5 ? 1 : static_cast<int>(radial.size() * arc / average) + 2);
    }
    for (std::size_t radial_band_index = 0;
         radial_band_index + 1 < radial.size(); ++radial_band_index) {
      double r0 = radial[radial_band_index];
      double r1 = radial[radial_band_index + 1];
      PlanarStraightLineGraph p;
      std::vector<Vec2> queries;
      for (int boundary_row_index = 0; boundary_row_index < 2;
           ++boundary_row_index) {
        int n = counts[radial_band_index + boundary_row_index];
        double r = radial[radial_band_index + boundary_row_index];
        for (int angular_sample_index = 0; angular_sample_index < n;
             ++angular_sample_index) {
          double fraction =
              n == 1 ? 0. : static_cast<double>(angular_sample_index) / (n - 1);
          double varying =
              n == 1 ? b * .5 + a * .5 : (1. - fraction) * a + fraction * b;
          auto point = uv(varying, r);
          auto query = point;
          query[radial_axis] += (r1 > r0 ? 1. : -1.) *
                                (boundary_row_index == 0 ? 1. : -1.) * 1e-6;
          if (angular_sample_index == 0) {
            query[angular_axis] += 1e-6;
          }
          if (n > 1 && angular_sample_index == n - 1) {
            query[angular_axis] -= 1e-6;
          }
          p.uv.push_back(point);
          queries.push_back(query);
        }
      }
      auto result = TriangleRun(p, o, {}, "Qz");
      if (result.uv != p.uv) {
        throw std::runtime_error(
            "Unexpected Delaunay vertex insertion or reordering");
      }
      Mesh local;
      local.uv = result.uv;
      local.triangles = result.triangles;
      local.normals = MeshNormals(s, queries);
      for (auto parameter_point : local.uv) {
        local.xyz.push_back(s.Evaluate(parameter_point));
      }
      for (auto& vertex_position : local.xyz) {
        for (auto& coordinate : vertex_position) {
          coordinate = Quantize(coordinate);
        }
      }
      for (auto& parameter_point : local.uv) {
        for (auto& coordinate : parameter_point) {
          coordinate = Quantize(coordinate);
        }
      }
      for (auto& vertex_normal : local.normals) {
        for (auto& coordinate : vertex_normal) {
          coordinate = Quantize(coordinate);
        }
        double n = Norm(vertex_normal);
        for (auto& coordinate : vertex_normal) {
          coordinate /= n;
        }
      }
      Append(out, std::move(local));
    }
  }
  ComponentIdentities(out, face.id);
  return out;
}
}  // namespace iko::internal
