#include "tessellation/private/internal.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <utility>
#include <vector>

#include "tessellation/tessellation.h"

namespace iko::internal {
static constexpr double kQuadratureNodes[] = {-.93247,  -.6612094, -.2386192,
                                              .2386192, .6612094,  .93247};
static constexpr double kQuadratureWeights[] = {.1713245, .3607616, .4679139,
                                                .4679139, .3607616, .1713245};
double AxisMap::Apply(double x, bool inverse) const {
  const auto& a = inverse ? target : source;
  const auto& b = inverse ? source : target;
  int span_index = std::clamp(
      static_cast<int>(std::lower_bound(a.begin(), a.end(), x) - a.begin()) - 1,
      0, static_cast<int>(a.size()) - 2);
  double scale =
      (b[span_index + 1] - b[span_index]) / (a[span_index + 1] - a[span_index]);
  double shift =
      (b[span_index] * a[span_index + 1] - b[span_index + 1] * a[span_index]) /
      (a[span_index + 1] - a[span_index]);
  return scale * x + shift;
}
static std::vector<int> ActiveSpans(const std::vector<double>& k, int d) {
  std::vector<int> indices;
  for (int knot_span_index = d;
       knot_span_index < static_cast<int>(k.size()) - d - 1;
       ++knot_span_index) {
    if (k[knot_span_index + 1] > k[knot_span_index]) {
      indices.push_back(knot_span_index);
    }
  }
  return indices;
}
static std::vector<double> StencilBasis(double t, const std::vector<double>& k,
                                        int d, int span,
                                        const std::vector<double>& prefix) {
  auto knot = [&](int knot_index) {
    if (knot_index >= 0) {
      return knot_index < static_cast<int>(k.size()) ? k[knot_index] : 0.;
    }
    return knot_index >= -static_cast<int>(prefix.size())
               ? prefix[prefix.size() + knot_index]
               : 0.;
  };
  std::vector<double> b(d + 1);
  std::vector<double> left(d + 1);
  std::vector<double> right(d + 1);
  b[0] = 1;
  for (int recurrence_degree = 1; recurrence_degree <= d; ++recurrence_degree) {
    left[recurrence_degree] = t - knot(span + 1 - recurrence_degree);
    right[recurrence_degree] = knot(span + recurrence_degree) - t;
    double saved = 0;
    for (int basis_index = 0; basis_index < recurrence_degree; ++basis_index) {
      double den =
          right[basis_index + 1] + left[recurrence_degree - basis_index];
      double temp = std::abs(den) >= std::numeric_limits<float>::min()
                        ? b[basis_index] / den
                        : 0.;
      b[basis_index] = saved + right[basis_index + 1] * temp;
      saved = left[recurrence_degree - basis_index] * temp;
    }
    b[recurrence_degree] = saved;
  }
  return b;
}
// The padded slabs preserve the original row stride, channels and U-then-V
// knot allocation. Slots outside the input allocation are exactly zero.
struct Stencil {
  const Surface& s;
  int axis;
  int channels;
  std::vector<double> source;
  std::vector<double> derivative;
  Stencil(const Surface& surface, int a, bool rational)
      : s(surface), axis(a), channels(rational ? 4 : 3) {
    int size = s.nu * s.nv;
    source.resize(channels * size);
    derivative.resize(channels * size);
    for (int channel_index = 0; channel_index < channels; ++channel_index) {
      for (int control_point_index = 0; control_point_index < size;
           ++control_point_index) {
        source[channel_index * size + control_point_index] =
            channel_index == 3
                ? s.weights[control_point_index]
                : s.control[control_point_index][channel_index] *
                      (rational ? s.weights[control_point_index] : 1.);
      }
    }
    const auto& knots = axis ? s.kv : s.ku;
    int degree = axis ? s.q : s.p;
    for (int channel_index = 0; channel_index < channels; ++channel_index) {
      for (int control_v_index = 0; control_v_index < s.nv - (axis == 1);
           ++control_v_index) {
        for (int control_u_index = 0; control_u_index < s.nu - (axis == 0);
             ++control_u_index) {
          int axis_control_index = axis ? control_v_index : control_u_index;
          double den = knots[axis_control_index + degree + 1] -
                       knots[axis_control_index + 1];
          double factor = den ? static_cast<double>(degree) / den : 0.;
          int control_point_index = control_v_index * s.nu + control_u_index;
          int offset = axis ? s.nu : 1;
          derivative[channel_index * size + control_point_index + offset] =
              (source[channel_index * size + control_point_index + offset] -
               source[channel_index * size + control_point_index]) *
              factor;
        }
      }
    }
  }
  Vec3 Evaluate(Vec2 uv, const std::array<std::vector<double>, 2>& bk) const {
    int spans[2];
    for (int parameter_axis = 0; parameter_axis < 2; ++parameter_axis) {
      const auto& k = bk[parameter_axis];
      auto indices = ActiveSpans(k, parameter_axis ? s.q : s.p);
      int span = static_cast<int>(
                     std::upper_bound(k.begin(), k.end(), uv[parameter_axis]) -
                     k.begin()) -
                 1;
      if (indices.size() < 5) {
        span = std::clamp(span, indices.front(), indices.back());
      } else if (span < 0) {
        span = indices.front();
      } else if (span >= static_cast<int>(k.size()) - 1) {
        span = indices.back();
      }
      spans[parameter_axis] = span;
    }
    auto eval = [&](const std::vector<double>& flat, int p, int q) {
      auto bu = StencilBasis(uv[0], bk[0], p, spans[0], {});
      auto bv = StencilBasis(uv[1], bk[1], q, spans[1], bk[0]);
      std::array<double, 4> value{};
      int size = s.nu * s.nv;
      for (int local_v_index = 0; local_v_index <= q; ++local_v_index) {
        std::array<double, 4> subtotal{};
        for (int local_u_index = 0; local_u_index <= p; ++local_u_index) {
          int index = (spans[1] - q + local_v_index) * s.nu + spans[0] - p +
                      local_u_index;
          for (int channel_index = 0; channel_index < channels;
               ++channel_index) {
            int pos = channel_index * size + index;
            double coef = pos >= 0 && pos < static_cast<int>(flat.size())
                              ? flat[pos]
                              : 0.;
            subtotal[channel_index] += bu[local_u_index] * coef;
          }
        }
        for (int channel_index = 0; channel_index < channels; ++channel_index) {
          value[channel_index] += bv[local_v_index] * subtotal[channel_index];
        }
      }
      return value;
    };
    auto d = eval(derivative, s.p - (axis == 0), s.q - (axis == 1));
    Vec3 result{d[0], d[1], d[2]};
    if (channels == 4) {
      auto point = eval(source, s.p, s.q);
      if (!Finite(point) || std::abs(point[3]) < kFloatEpsilon) {
        throw Unsupported(
            "Bounded rational arc stencil has unresolved homogeneous "
            "denominator");
      }
      double iw = 1. / point[3];
      Vec3 xyz{point[0] * iw, point[1] * iw, point[2] * iw};
      result = Multiply(Subtract(result, Multiply(xyz, d[3])), iw);
    }
    return result;
  }
};
static AxisMap FullMap(const Surface& s, int axis, bool padded, bool rational) {
  const auto& k = axis ? s.kv : s.ku;
  const auto& other = axis ? s.ku : s.kv;
  auto indices = ActiveSpans(k, axis ? s.q : s.p);
  auto original = k;
  for (auto& knot : original) {
    knot -= k.front();
  }
  auto rewritten = original;
  std::array<std::vector<double>, 2> bk{s.ku, s.kv};
  bk[axis] = original;
  Stencil stencil(s, axis, rational);
  double scale = 1.;
  int previous = -1;
  auto all = indices;
  all.push_back(static_cast<int>(k.size()) - 1);
  for (std::size_t span_position = 0; span_position < all.size();
       ++span_position) {
    int index = all[span_position];
    int start = previous + 1;
    if (span_position) {
      for (int knot_index = start; knot_index <= index; ++knot_index) {
        rewritten[knot_index] *= scale;
      }
    }
    if (span_position < indices.size()) {
      double a = original[start];
      double b = original[index + 1];
      double total = 0;
      for (double transverse_fraction : {0., .25, .5, .75, 1.}) {
        double c = (1 - transverse_fraction) * other.front() +
                   transverse_fraction * other.back();
        double value = 0.;
        for (int quadrature_index = 0; quadrature_index < 6;
             ++quadrature_index) {
          double varying =
              (b - a) * .5 * kQuadratureNodes[quadrature_index] + (a + b) * .5;
          Vec2 uv = axis ? Vec2{c, varying} : Vec2{varying, c};
          auto d = padded || rational ? stencil.Evaluate(uv, bk)
                                      : s.Derivatives(uv, &bk, true)[axis];
          value += std::sqrt((d[0] * d[0] + d[1] * d[1]) + d[2] * d[2]) *
                   kQuadratureWeights[quadrature_index];
        }
        total += (b - a) * .5 * value;
      }
      if (!std::isfinite(total)) {
        throw Unsupported("Nonclamped native arc integral is nonfinite");
      }
      double length = std::max(1., total / 5.);
      if (std::abs(b) > 1e-8) {
        scale = (scale * a + length) / b;
      }
    }
    previous = index;
  }
  double centre = rewritten.back() * .5;
  for (auto& mapped_knot : rewritten) {
    mapped_knot -= centre;
  }
  indices.push_back(indices.back() + 1);
  AxisMap map;
  for (int knot_index : indices) {
    map.source.push_back(k[knot_index]);
    map.target.push_back(rewritten[knot_index]);
  }
  for (std::size_t mapped_knot_index = 0; mapped_knot_index < map.target.size();
       ++mapped_knot_index) {
    if (!std::isfinite(map.target[mapped_knot_index]) ||
        (mapped_knot_index &&
         map.target[mapped_knot_index] <= map.target[mapped_knot_index - 1])) {
      throw Unsupported(
          "Nonclamped native arc map is nonfinite or nonmonotone");
    }
  }
  return map;
}
static bool Uniform(const std::vector<double>& k) {
  double magnitude = 1.;
  for (double knot : k) {
    magnitude = std::max(magnitude, std::abs(knot));
  }
  double roundoff = 16 * kDoubleEpsilon * magnitude;
  double first = k[1] - k[0];
  for (std::size_t knot_index = 1; knot_index < k.size(); ++knot_index) {
    if (!(k[knot_index] - k[knot_index - 1] > 0) ||
        std::abs(k[knot_index] - k[knot_index - 1] - first) > roundoff) {
      return false;
    }
  }
  return true;
}
static bool Repeated(const Surface& s, int axis, int count, bool weights) {
  for (int control_v_index = 0; control_v_index < (axis ? count : s.nv);
       ++control_v_index) {
    for (int control_u_index = 0; control_u_index < (axis ? s.nu : count);
         ++control_u_index) {
      int a = control_v_index * s.nu + control_u_index;
      int b = axis ? (s.nv - count + control_v_index) * s.nu + control_u_index
                   : control_v_index * s.nu + s.nu - count + control_u_index;
      if (s.control[a] != s.control[b] ||
          (weights && s.weights[a] != s.weights[b])) {
        return false;
      }
    }
  }
  return true;
}
ArcMap::ArcMap(const Surface& s) {
  bool nc[2];
  bool short_or_clamped[2];
  bool complete = true;
  for (int axis = 0; axis < 2; ++axis) {
    const auto& k = axis ? s.kv : s.ku;
    int d = axis ? s.q : s.p;
    double lo = k[d];
    double hi = k[k.size() - d - 1];
    nc[axis] = !(k.front() == lo && k.back() == hi);
    short_or_clamped[axis] = !nc[axis] || ActiveSpans(k, d).size() < 5;
    complete &= s.domain[2 * axis] == lo && s.domain[2 * axis + 1] == hi;
  }
  bool poly = !s.explicit_weights;
  bool unit = std::all_of(s.weights.begin(), s.weights.end(),
                          [](double weight) { return weight == 1.; });
  auto strict = [](const std::vector<double>& k) {
    for (std::size_t knot_index = 1; knot_index < k.size(); ++knot_index) {
      if (k[knot_index] <= k[knot_index - 1]) {
        return false;
      }
    }
    return true;
  };
  bool quadratic = s.p == 2 && s.q == 2 && poly && (nc[0] != nc[1]) &&
                   !(short_or_clamped[0] && short_or_clamped[1]) &&
                   (!nc[0] || strict(s.ku)) && (!nc[1] || strict(s.kv));
  bool cubic_v = s.p == 2 && s.q == 3 && poly && unit && !nc[0] && nc[1] &&
                 !short_or_clamped[1] && s.ku.size() == 6 && s.nu == 3 &&
                 Uniform(s.kv) && Repeated(s, 1, 3, false);
  bool cubic_u = s.p == 3 && s.q == 2 && poly && unit && nc[0] && !nc[1] &&
                 !short_or_clamped[0] && s.kv.size() == 6 && s.nv == 3 &&
                 Uniform(s.ku) && Repeated(s, 0, 3, false);
  auto minmax = std::minmax_element(s.weights.begin(), s.weights.end());
  bool rational = s.p == 2 && s.q == 2 && !poly &&
                  *minmax.second > *minmax.first && nc[0] && !nc[1] &&
                  !short_or_clamped[0] && s.kv.size() == 6 && s.nv == 3 &&
                  Uniform(s.ku) && Repeated(s, 0, 2, true);
  bool padded = quadratic || cubic_v || cubic_u;
  if ((nc[0] || nc[1]) &&
      ((short_or_clamped[0] && short_or_clamped[1]) || padded || rational) &&
      complete) {
    for (int axis = 0; axis < 2; ++axis) {
      maps[axis] = FullMap(s, axis, padded, rational);
    }
    return;
  }
  for (int axis = 0; axis < 2; ++axis) {
    const auto& knots = axis ? s.kv : s.ku;
    double lo = s.domain[2 * axis];
    double hi = s.domain[2 * axis + 1];
    std::vector<double> breaks{lo, hi};
    for (double knot : knots) {
      if (lo < knot && knot < hi) {
        breaks.push_back(knot);
      }
    }
    breaks = Unique(breaks);
    auto integration = breaks;
    for (auto& integration_parameter : integration) {
      integration_parameter -= knots.front();
    }
    std::array<std::vector<double>, 2> bk{s.ku, s.kv};
    for (auto& basis_knot : bk[axis]) {
      basis_knot -= knots.front();
    }
    std::vector<double> target{0.};
    double previous = 1;
    for (std::size_t integration_span_index = 0;
         integration_span_index + 1 < integration.size();
         ++integration_span_index) {
      double a = integration[integration_span_index];
      double b = integration[integration_span_index + 1];
      double total = 0.;
      for (double transverse_fraction : {0., .25, .5, .75, 1.}) {
        int other = 1 - axis;
        double c = (1 - transverse_fraction) * s.domain[2 * other] +
                   transverse_fraction * s.domain[2 * other + 1];
        double value = 0;
        for (int quadrature_index = 0; quadrature_index < 6;
             ++quadrature_index) {
          double varying =
              (b - a) * .5 * kQuadratureNodes[quadrature_index] + (a + b) * .5;
          auto d = s.Derivatives(axis ? Vec2{c, varying} : Vec2{varying, c},
                                 &bk)[axis];
          value += std::sqrt((d[0] * d[0] + d[1] * d[1]) + d[2] * d[2]) *
                   kQuadratureWeights[quadrature_index];
        }
        total += (b - a) * .5 * value;
      }
      double length = std::max(1., total / 5.);
      previous = (previous * a + length) / b;
      target.push_back(previous * b);
    }
    double centre = target.back() * .5;
    for (auto& mapped_parameter : target) {
      mapped_parameter -= centre;
    }
    maps[axis] = {std::move(breaks), std::move(target)};
  }
}
}  // namespace iko::internal
