#include "tessellation/private/internal.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "tessellation/tessellation.h"

namespace iko::internal {
using HomogeneousPoint = std::array<double, 4>;
using Net = std::vector<std::vector<HomogeneousPoint>>;
static void InsertOnce(std::vector<double>& knots, Net& control, int degree,
                       double value) {
  int insertion_span_index =
      static_cast<int>(std::upper_bound(knots.begin(), knots.end(), value) -
                       knots.begin()) -
      1;
  int knot_multiplicity =
      static_cast<int>(std::count(knots.begin(), knots.end(), value));
  if (!(knots[degree] <= value && value <= knots[knots.size() - degree - 1])) {
    throw std::invalid_argument("Insertion outside active domain");
  }
  if (knot_multiplicity > degree) {
    throw std::invalid_argument("Endpoint already has full multiplicity");
  }
  Net q(control.size() + 1,
        std::vector<HomogeneousPoint>(control.front().size()));
  for (int control_index = 0; control_index < insertion_span_index - degree + 1;
       ++control_index) {
    q[control_index] = control[control_index];
  }
  for (int control_index = insertion_span_index - knot_multiplicity;
       control_index < static_cast<int>(control.size()); ++control_index) {
    q[control_index + 1] = control[control_index];
  }
  for (int control_index = insertion_span_index - degree + 1;
       control_index < insertion_span_index - knot_multiplicity + 1;
       ++control_index) {
    double alpha = (value - knots[control_index]) /
                   (knots[control_index + degree] - knots[control_index]);
    for (std::size_t transverse_index = 0;
         transverse_index < control[control_index].size(); ++transverse_index) {
      for (int homogeneous_channel = 0; homogeneous_channel < 4;
           ++homogeneous_channel) {
        // Match separate NumPy products followed by addition, including
        // the weight channel. Strict floating point disallows contraction.
        q[control_index][transverse_index][homogeneous_channel] =
            alpha *
                control[control_index][transverse_index][homogeneous_channel] +
            (1. - alpha) * control[control_index - 1][transverse_index]
                                  [homogeneous_channel];
      }
    }
  }
  knots.insert(knots.begin() + insertion_span_index + 1, value);
  control = std::move(q);
}
static void ClampAxis(std::vector<double>& knots, Net& control, int degree) {
  double lo = knots[degree];
  double hi = knots[knots.size() - degree - 1];
  for (double endpoint_parameter : {lo, hi}) {
    while (std::count(knots.begin(), knots.end(), endpoint_parameter) <
           degree + 1) {
      InsertOnce(knots, control, degree, endpoint_parameter);
    }
  }
  auto start = std::find(knots.begin(), knots.end(), lo) - knots.begin();
  auto end =
      std::find(knots.rbegin(), knots.rend(), hi).base() - knots.begin() - 1;
  knots = std::vector<double>(knots.begin() + start, knots.begin() + end + 1);
  // Form the final offset first; the knot index can exceed the control count.
  control = Net(control.begin() + start, control.begin() + (end - degree));
}
Surface ClampActiveDomain(const Surface& source) {
  Surface s = source;
  // Python first multiplies XYZ by weight, then clamps U, then clamps V.
  Net u(s.nu, std::vector<HomogeneousPoint>(s.nv));
  for (int control_u_index = 0; control_u_index < s.nu; ++control_u_index) {
    for (int control_v_index = 0; control_v_index < s.nv; ++control_v_index) {
      int control_point_index = control_v_index * s.nu + control_u_index;
      auto cp = s.control[control_point_index];
      double w = s.weights[control_point_index];
      u[control_u_index][control_v_index] = {cp[0] * w, cp[1] * w, cp[2] * w,
                                             w};
    }
  }
  ClampAxis(s.ku, u, s.p);
  s.nu = static_cast<int>(u.size());
  Net v(s.nv, std::vector<HomogeneousPoint>(s.nu));
  for (int control_v_index = 0; control_v_index < s.nv; ++control_v_index) {
    for (int control_u_index = 0; control_u_index < s.nu; ++control_u_index) {
      v[control_v_index][control_u_index] = u[control_u_index][control_v_index];
    }
  }
  ClampAxis(s.kv, v, s.q);
  s.nv = static_cast<int>(v.size());
  s.weights.resize(s.nu * s.nv);
  s.control.resize(s.nu * s.nv);
  for (int control_v_index = 0; control_v_index < s.nv; ++control_v_index) {
    for (int control_u_index = 0; control_u_index < s.nu; ++control_u_index) {
      auto h = v[control_v_index][control_u_index];
      int control_point_index = control_v_index * s.nu + control_u_index;
      s.weights[control_point_index] = h[3];
      s.control[control_point_index] = {h[0] / h[3], h[1] / h[3], h[2] / h[3]};
    }
  }
  s.label += "_clamped";
  s.closed_axes = {};
  return s;
}
static std::array<Vec3, 4> RectangleCorners(const Surface& s) {
  const auto& d = s.domain;
  return {
      {{d[0], d[3], 0.}, {d[0], d[2], 0.}, {d[1], d[2], 0.}, {d[1], d[3], 0.}}};
}
void ValidateNaturalRectangle(const Face& face) {
  if (!face.inner.empty() || !face.outer || face.outer->size() != 4) {
    throw std::invalid_argument(
        "Natural rectangle sampling requires its constructed natural "
        "rectangle");
  }
  auto corners = RectangleCorners(face.surface);
  for (int edge_index = 0; edge_index < 4; ++edge_index) {
    const auto& c = (*face.outer)[edge_index];
    if (c.degree != 1 || c.weights != std::vector<double>{1., 1.} ||
        c.control != std::vector<Vec3>{corners[edge_index],
                                       corners[(edge_index + 1) % 4]}) {
      throw std::invalid_argument(
          "Natural rectangle sampling requires unmodified "
          "constructed polynomial rectangle edges");
    }
  }
}
Vec3 RectangleSamplingPoint(const Surface& s, Vec2 uv) {
  if (!Finite(uv)) {
    throw std::invalid_argument("Nonfinite natural rectangle sample");
  }
  uv = {std::clamp(uv[0], s.domain[0], s.domain[1]),
        std::clamp(uv[1], s.domain[2], s.domain[3])};
  return s.Evaluate(uv);
}
}  // namespace iko::internal
