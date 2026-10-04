#include "tessellation/private/internal.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "tessellation/tessellation.h"

namespace iko::internal {
static Curve PrepareCurve(Curve curve) {
  static constexpr double kQuadratureNodes[] = {
      -.98799251802049, -.93727339240071, -.8482065834104, -.72441773136017,
      -.57097217260854, -.39415134707756, -.2011940939974, 0.,
      .20119409399743,  .3941513470776,   .57097217260854, .72441773136017,
      .84820658341043,  .93727339240071,  .98799251802049};
  static constexpr double kQuadratureWeights[] = {
      .030753241996117, .07036604748811, .1071592204672,  .13957067792615,
      .16626920581699,  .1861610000156,  .19843148532711, .2025782419256,
      .1984314853271,   .18616100001556, .166269205817,   .13957067792615,
      .10715922046717,  .0703660474881,  .03075324199612};
  if (!AlignedDomain(curve.knots, curve.degree, curve.domain)) {
    throw Unsupported(
        "Restricted trim-curve parameter domains are not aligned");
  }
  auto source_knots = curve.knots;
  double origin = curve.knots.front();
  for (auto& knot : curve.knots) {
    knot -= origin;
  }
  int d = curve.degree;
  int n = static_cast<int>(curve.weights.size());
  double a = curve.knots[d];
  double b = curve.knots[curve.knots.size() - d - 1];
  if (d < 1 || a != 0. || b != curve.knots.back() ||
      !std::is_sorted(curve.knots.begin(), curve.knots.end())) {
    throw Unsupported("Unclamped or invalid trim-curve knots are not aligned");
  }
  Surface s;
  s.p = d;
  s.q = 1;
  s.nu = n;
  s.nv = 2;
  // Curve derivative coefficients are initialized before the knot-origin
  // shift. Retain those intervals and shift only the evaluation basis.
  s.ku = source_knots;
  s.kv = {0, 0, 1, 1};
  s.domain = {source_knots[d], source_knots[source_knots.size() - d - 1], 0, 1};
  s.control = curve.control;
  s.control.insert(s.control.end(), curve.control.begin(), curve.control.end());
  s.weights = curve.weights;
  s.weights.insert(s.weights.end(), curve.weights.begin(), curve.weights.end());
  s.explicit_weights = !std::all_of(curve.weights.begin(), curve.weights.end(),
                                    [](double weight) { return weight == 1.; });
  std::array<std::vector<double>, 2> basis_knots{curve.knots, s.kv};
  auto original = curve.knots;
  auto rewritten = original;
  std::vector<int> spans;
  for (int knot_span_index = 0;
       knot_span_index + 1 < static_cast<int>(original.size());
       ++knot_span_index) {
    if (original[knot_span_index + 1] > original[knot_span_index]) {
      spans.push_back(knot_span_index);
    }
  }
  auto all = spans;
  all.push_back(static_cast<int>(original.size()) - 1);
  double scale = 1.;
  int previous = -1;
  for (std::size_t span_position = 0; span_position < all.size();
       ++span_position) {
    int index = all[span_position];
    int start = previous + 1;
    if (span_position) {
      for (int knot_index = start; knot_index <= index; ++knot_index) {
        rewritten[knot_index] *= scale;
      }
    }
    if (span_position < spans.size()) {
      double lo = original[start];
      double hi = original[index + 1];
      double total = 0;
      for (int quadrature_index = 0; quadrature_index < 15;
           ++quadrature_index) {
        double t = (hi - lo) * .5 * kQuadratureNodes[quadrature_index] +
                   (lo + hi) * .5;
        auto der = s.Derivatives({t, 0}, &basis_knots)[0];
        total +=
            std::sqrt((der[0] * der[0] + der[1] * der[1]) + der[2] * der[2]) *
            kQuadratureWeights[quadrature_index];
      }
      double length = std::max(1., (hi - lo) * .5 * total);
      if (std::abs(hi) > 1e-8) {
        scale = (scale * lo + length) / hi;
      }
    }
    previous = index;
  }
  double centre = rewritten.back() * .5;
  for (auto& mapped_knot : rewritten) {
    mapped_knot -= centre;
  }
  curve.knots = std::move(rewritten);
  curve.domain = {curve.knots[d], curve.knots[curve.knots.size() - d - 1]};
  curve.prepared = true;
  return curve;
}
static Vec3 SamplingPoint(const Surface& s, Vec2 uv) {
  bool outside = false;
  for (int axis = 0; axis < 2; ++axis) {
    const auto& k = axis ? s.kv : s.ku;
    int d = axis ? s.q : s.p;
    double lo = k[d];
    double hi = k[k.size() - d - 1];
    if (uv[axis] < lo || uv[axis] > hi) {
      outside = true;
      if (s.closed_axes[axis]) {
        throw Unsupported(
            "Off-domain sampling on a declared closed axis is not "
            "reproduced");
      }
      std::vector<double> active;
      for (double knot : k) {
        if (knot >= lo && knot <= hi) {
          active.push_back(knot);
        }
      }
      if (Unique(active).size() - 1 >= 5) {
        throw Unsupported(
            "Off-domain native span selection with five active spans is "
            "not reproduced");
      }
    }
  }
  auto result = outside ? s.EvaluateControl(uv, true) : s.Evaluate(uv);
  if (!Finite(result)) {
    throw std::invalid_argument("Nonfinite extrapolated trim sampling point");
  }
  return result;
}
static std::vector<double> Crossings(const Curve& c, int axis, double value,
                                     const Options& o) {
  auto r =
      Sample([&](double t) { return c.Evaluate(t); }, c.domain[0], c.domain[1],
             std::max(1e-5, o.curve.chord), std::max(0., o.curve.min_edge),
             std::max(.01, o.curve.max_edge));
  std::vector<double> roots;
  for (std::size_t sample_interval_index = 0;
       sample_interval_index + 1 < r.parameters.size();
       ++sample_interval_index) {
    double a = r.parameters[sample_interval_index];
    double b = r.parameters[sample_interval_index + 1];
    double fa = r.points[sample_interval_index][axis] - value;
    double fb = r.points[sample_interval_index + 1][axis] - value;
    if (fa * fb > 1e-16) {
      continue;
    }
    if (std::abs(fa) < 1e-8) {
      roots.push_back(a);
      continue;
    }
    if (std::abs(fb) < 1e-8) {
      roots.push_back(b);
      continue;
    }
    bool found = false;
    double t = 0;
    double fc = INFINITY;
    for (int root_iteration = 0; root_iteration < 20; ++root_iteration) {
      if (std::abs(fb - fa) <= 1e-8) {
        break;
      }
      t = a - (b - a) * fa / (fb - fa);
      found = true;
      fc = c.Evaluate(t)[axis] - value;
      if (fc * fb < 0) {
        a = b;
        fa = fb;
      } else {
        fa *= std::abs(fc + fb) > 1e-8 ? fb / (fc + fb) : 1.;
      }
      b = t;
      fb = fc;
      if (std::abs(fc) < 1e-10) {
        break;
      }
    }
    if (found && std::abs(fc) < 1e-7) {
      roots.push_back(t);
    }
  }
  return roots;
}
static std::vector<Vec2> FilterLoop(const std::vector<Vec2>& points,
                                    const std::set<int>& junctions) {
  if (points.empty()) {
    throw std::invalid_argument("Empty trim loop");
  }
  Vec2 lo = points[0];
  Vec2 hi = lo;
  for (auto parameter_point : points) {
    for (int axis = 0; axis < 2; ++axis) {
      lo[axis] = std::min(lo[axis], parameter_point[axis]);
      hi[axis] = std::max(hi[axis], parameter_point[axis]);
    }
  }
  auto scaled = Multiply(Subtract(hi, lo), 1e-4);
  double threshold = scaled[1] * scaled[1] + scaled[0] * scaled[0];
  std::vector<int> kept{0};
  for (int vertex_step = 1; vertex_step <= static_cast<int>(points.size());
       ++vertex_step) {
    int source = vertex_step % static_cast<int>(points.size());
    auto delta = Subtract(points[source], points[kept.back()]);
    double d2 = delta[1] * delta[1] + delta[0] * delta[0];
    if (d2 > threshold || ((junctions.count(vertex_step) ||
                            vertex_step == static_cast<int>(points.size())) &&
                           d2 > 1e-12)) {
      kept.push_back(source);
    }
  }
  kept.pop_back();
  std::vector<Vec2> out;
  for (int kept_vertex_index : kept) {
    out.push_back(points[kept_vertex_index]);
  }
  return out;
}
static std::vector<std::array<int, 2>> Cleanup(
    std::vector<std::vector<Vec2>>& loops,
    const std::array<std::vector<double>, 2>& axes) {
  for (int axis = 0; axis < 2; ++axis) {
    for (auto& loop : loops) {
      double low = loop[0][axis];
      double high = low;
      for (auto parameter_point : loop) {
        low = std::min(low, parameter_point[axis]);
        high = std::max(high, parameter_point[axis]);
      }
      double r = std::sqrt(high - low);
      double small = 1e-8 * r;
      double gate = 1e-5 * r;
      for (double knot_parameter : axes[axis]) {
        if (knot_parameter < low + small || knot_parameter > high - small) {
          continue;
        }
        for (int loop_vertex_index = 0;
             loop_vertex_index < static_cast<int>(loop.size());
             ++loop_vertex_index) {
          double previous =
              loop[(loop_vertex_index - 1 + loop.size()) % loop.size()][axis];
          double current = loop[loop_vertex_index][axis];
          double next = loop[(loop_vertex_index + 1) % loop.size()][axis];
          if (!((previous < current && current < next) ||
                (previous > current && current > next))) {
            continue;
          }
          if (!(std::min(previous, next) + small < knot_parameter &&
                knot_parameter < std::max(previous, next) - small)) {
            continue;
          }
          double delta = knot_parameter - current;
          bool accepted = -gate <= delta && delta <= gate;
          if (!accepted) {
            accepted = -50 * gate < delta && delta < 50 * gate &&
                       std::abs(current - previous) > 500 * gate &&
                       std::abs(current - next) > 500 * gate;
          }
          if (accepted) {
            loop[loop_vertex_index][axis] = knot_parameter;
          }
        }
      }
    }
  }
  std::vector<Vec2> points;
  std::vector<std::array<int, 2>> ranges;
  std::vector<std::array<int, 2>> internal;
  std::vector<double> limits;
  for (const auto& loop : loops) {
    int start = static_cast<int>(points.size());
    points.insert(points.end(), loop.begin(), loop.end());
    points.push_back(loop[0]);
    ranges.push_back({start, static_cast<int>(points.size()) - 1});
    Vec2 lo = loop[0];
    Vec2 hi = lo;
    for (auto parameter_point : loop) {
      for (int axis = 0; axis < 2; ++axis) {
        lo[axis] = std::min(lo[axis], parameter_point[axis]);
        hi[axis] = std::max(hi[axis], parameter_point[axis]);
      }
    }
    limits.push_back(.01 * std::max(hi[0] - lo[0], hi[1] - lo[1]));
  }
  for (int axis = 0; axis < 2; ++axis) {
    for (double knot_parameter : axes[axis]) {
      std::vector<int> ids;
      for (int vertex_index = 0; vertex_index < static_cast<int>(points.size());
           ++vertex_index) {
        if (std::abs(points[vertex_index][axis] - knot_parameter) < 1e-11 &&
            std::none_of(ranges.begin(), ranges.end(), [&](auto loop_range) {
              return loop_range[1] == vertex_index;
            })) {
          ids.push_back(vertex_index);
        }
      }
      bool odd = false;
      for (auto loop_range : ranges) {
        odd |= std::count_if(ids.begin(), ids.end(),
                             [&](int vertex_index) {
                               return loop_range[0] <= vertex_index &&
                                      vertex_index < loop_range[1];
                             }) %
                   2 !=
               0;
      }
      if (odd) {
        continue;
      }
      std::stable_sort(ids.begin(), ids.end(),
                       [&](int first_vertex_index, int second_vertex_index) {
                         return points[first_vertex_index][1 - axis] >
                                points[second_vertex_index][1 - axis];
                       });
      for (std::size_t intersection_pair_index = 0;
           intersection_pair_index + 1 < ids.size();
           intersection_pair_index += 2) {
        internal.push_back(
            {ids[intersection_pair_index], ids[intersection_pair_index + 1]});
      }
    }
  }
  struct Candidate {
    int neighbor;
    int owner;
    bool enabled = true;
  };
  std::vector<Candidate> candidates;
  for (auto segment : internal) {
    for (int owner_vertex_index : segment) {
      if (owner_vertex_index == 0 ||
          owner_vertex_index >= static_cast<int>(points.size()) - 1) {
        continue;
      }
      auto found =
          std::find_if(ranges.begin(), ranges.end(), [&](auto loop_range) {
            return loop_range[0] <= owner_vertex_index &&
                   owner_vertex_index <= loop_range[1];
          });
      double limit = limits[found - ranges.begin()];
      for (int neighbor_vertex_index :
           {owner_vertex_index - 1, owner_vertex_index + 1}) {
        if (std::abs(points[neighbor_vertex_index][0] -
                     points[owner_vertex_index][0]) < limit &&
            std::abs(points[neighbor_vertex_index][1] -
                     points[owner_vertex_index][1]) < limit) {
          candidates.push_back({neighbor_vertex_index, owner_vertex_index});
        }
      }
    }
  }
  std::stable_sort(candidates.begin(), candidates.end(),
                   [](auto first_candidate, auto second_candidate) {
                     return first_candidate.neighbor >
                            second_candidate.neighbor;
                   });
  for (auto& candidate : candidates) {
    if (!candidate.enabled) {
      continue;
    }
    int neighbor = candidate.neighbor;
    int owner = candidate.owner;
    bool protected_vertex =
        std::any_of(ranges.begin(), ranges.end(),
                    [&](auto loop_range) {
                      return loop_range[0] == neighbor ||
                             loop_range[1] == neighbor;
                    }) ||
        std::any_of(internal.begin(), internal.end(), [&](auto segment) {
          return segment[0] == neighbor || segment[1] == neighbor;
        });
    if (neighbor >= static_cast<int>(points.size()) || protected_vertex) {
      int internal_segment_index = 0;
      int internal_segment_scan_limit = static_cast<int>(internal.size()) / 2;
      while (internal_segment_index < internal_segment_scan_limit) {
        if (internal[internal_segment_index][0] == owner ||
            internal[internal_segment_index][1] == owner) {
          internal.erase(internal.begin() + internal_segment_index);
          --internal_segment_scan_limit;
        } else {
          ++internal_segment_index;
        }
      }
      for (auto& other_candidate : candidates) {
        if (other_candidate.owner == owner) {
          other_candidate.enabled = false;
        }
      }
      continue;
    }
    points.erase(points.begin() + neighbor);
    for (auto& loop_range : ranges) {
      for (auto& range_endpoint_index : loop_range) {
        if (range_endpoint_index > neighbor) {
          --range_endpoint_index;
        }
      }
    }
    for (auto& segment : internal) {
      for (auto& vertex_index : segment) {
        if (vertex_index > neighbor) {
          --vertex_index;
        }
      }
    }
    for (auto& candidate : candidates) {
      if (candidate.owner > neighbor) {
        --candidate.owner;
      }
    }
  }
  std::vector<int> mapping(points.size(), -1);
  int unique_vertex_count = 0;
  for (int vertex_index = 0; vertex_index < static_cast<int>(points.size());
       ++vertex_index) {
    if (std::none_of(ranges.begin(), ranges.end(), [&](auto loop_range) {
          return loop_range[1] == vertex_index;
        })) {
      mapping[vertex_index] = unique_vertex_count++;
    }
  }
  loops.clear();
  for (auto loop_range : ranges) {
    loops.emplace_back(points.begin() + loop_range[0],
                       points.begin() + loop_range[1]);
  }
  for (auto& edge : internal) {
    for (auto& vertex_index : edge) {
      vertex_index = mapping[vertex_index];
    }
  }
  return internal;
}
bool PrepareLoop(std::vector<Vec2>& p) {
  Vec2 low{INFINITY, INFINITY};
  Vec2 high{-INFINITY, -INFINITY};
  bool previous_repaired = false;
  double twice_area = 0.;
  for (int loop_vertex_index = 0;
       loop_vertex_index < static_cast<int>(p.size()); ++loop_vertex_index) {
    auto current = p[loop_vertex_index];
    auto previous = p[(loop_vertex_index - 1 + p.size()) % p.size()];
    auto next = p[(loop_vertex_index + 1) % p.size()];
    for (int axis = 0; axis < 2; ++axis) {
      low[axis] = std::min(low[axis], current[axis]);
      high[axis] = std::max(high[axis], current[axis]);
    }
    auto a = Subtract(current, previous);
    auto b = Subtract(next, current);
    for (auto* edge_direction : {&a, &b}) {
      double l2 = (*edge_direction)[1] * (*edge_direction)[1] +
                  (*edge_direction)[0] * (*edge_direction)[0];
      if (l2 > 1e-16) {
        *edge_direction = Multiply(*edge_direction, 1. / std::sqrt(l2));
      }
    }
    bool reversal = (b[1] * a[1] + b[0] * a[0]) + 1. < 1e-5;
    if (reversal) {
      if (previous_repaired) {
        return false;
      }
      current = Multiply(Add(next, previous), .5);
      p[loop_vertex_index] = current;
    }
    previous_repaired = reversal;
    twice_area += next[1] * current[0] - current[1] * next[0];
  }
  double area = std::abs(twice_area * .5);
  double bbox = (high[1] - low[1]) * (high[0] - low[0]);
  return bbox >= 1e-6 && area / bbox >= .001;
}
static bool Equal(double a, double b) {
  return std::abs(a - b) <=
         kDoubleEpsilon * std::max({1., std::abs(a), std::abs(b)});
}
static bool Equal(Vec2 a, Vec2 b) {
  return Equal(a[0], b[0]) && Equal(a[1], b[1]);
}
static bool Less(Vec2 a, Vec2 b) {
  if (!Equal(a[0], b[0])) {
    return a[0] < b[0];
  }
  return !Equal(a[1], b[1]) && a[1] < b[1];
}
static int Side(Vec2 a, Vec2 b, Vec2 c) {
  if (Equal(a, b) || Equal(a, c) || Equal(b, c)) {
    return 0;
  }
  if (Less(c, a)) {
    if (Less(c, b)) {
      auto old = a;
      a = c;
      c = b;
      b = old;
    } else {
      auto old = a;
      a = b;
      b = c;
      c = old;
    }
  } else if (!Less(a, b)) {
    auto old = a;
    a = b;
    b = c;
    c = old;
  }
  double dx = b[0] - a[0];
  double dy = b[1] - a[1];
  double ex = c[0] - a[0];
  double ey = c[1] - a[1];
  double det = ey * dx - ex * dy;
  if (std::abs(det) <=
      kDoubleEpsilon * std::max({1., std::abs(dx), std::abs(dy), std::abs(ex),
                                 std::abs(ey)})) {
    return 0;
  }
  return det > 0 ? 1 : -1;
}
std::array<int, 4> SegmentSides(Vec2 a, Vec2 b, Vec2 c, Vec2 d) {
  return {Side(a, b, c), Side(a, b, d), Side(c, d, a), Side(c, d, b)};
}
bool OuterIntersects(const std::vector<Vec2>& p) {
  for (std::size_t first_segment_index = 0; first_segment_index < p.size();
       ++first_segment_index) {
    auto a = p[first_segment_index];
    auto b = p[(first_segment_index + 1) % p.size()];
    for (std::size_t second_segment_index = first_segment_index + 2;
         second_segment_index < p.size() - (first_segment_index == 0);
         ++second_segment_index) {
      auto c = p[second_segment_index];
      auto d = p[(second_segment_index + 1) % p.size()];
      bool separated = false;
      for (int axis = 0; axis < 2; ++axis) {
        double lower =
            std::max(std::min(a[axis], b[axis]), std::min(c[axis], d[axis]));
        double upper =
            std::min(std::max(a[axis], b[axis]), std::max(c[axis], d[axis]));
        separated |=
            lower - upper >
            kDoubleEpsilon * std::max({1., std::abs(lower), std::abs(upper)});
      }
      if (!separated && Side(a, b, c) * Side(a, b, d) <= 0 &&
          Side(c, d, a) * Side(c, d, b) <= 0) {
        return true;
      }
    }
  }
  return false;
}
static bool Quality(const std::vector<std::vector<Vec2>>& loops) {
  if (loops.empty()) {
    return false;
  }
  Vec2 low{INFINITY, INFINITY};
  Vec2 high{-INFINITY, -INFINITY};
  std::vector<double> areas;
  for (const auto& loop : loops) {
    if (loop.size() < 3) {
      return false;
    }
    double area = 0.;
    for (std::size_t loop_vertex_index = 0; loop_vertex_index < loop.size();
         ++loop_vertex_index) {
      auto a = loop[loop_vertex_index];
      auto b = loop[(loop_vertex_index + 1) % loop.size()];
      if (!Finite(a)) {
        return false;
      }
      for (int axis = 0; axis < 2; ++axis) {
        low[axis] = std::min(low[axis], a[axis]);
        high[axis] = std::max(high[axis], a[axis]);
      }
      area += (a[1] - b[1]) * (b[0] + a[0]);
    }
    areas.push_back(area * .5);
  }
  double total = 0;
  for (std::size_t inner_loop_index = 1; inner_loop_index < areas.size();
       ++inner_loop_index) {
    total += areas[inner_loop_index];
  }
  total += areas[0];
  if (high[0] - low[0] < 1e-6 || high[1] - low[1] < 1e-6 ||
      std::abs(total) < 1e-6) {
    return false;
  }
  auto& p = loops[0];
  for (std::size_t first_segment_index = 0; first_segment_index < p.size();
       ++first_segment_index) {
    auto a = p[first_segment_index];
    auto b = p[(first_segment_index + 1) % p.size()];
    for (std::size_t second_segment_index = first_segment_index + 2;
         second_segment_index < p.size() - (first_segment_index == 0);
         ++second_segment_index) {
      auto c = p[second_segment_index];
      auto d = p[(second_segment_index + 1) % p.size()];
      bool separated = false;
      for (int axis = 0; axis < 2; ++axis) {
        double lower =
            std::max(std::min(a[axis], b[axis]), std::min(c[axis], d[axis]));
        double upper =
            std::min(std::max(a[axis], b[axis]), std::max(c[axis], d[axis]));
        separated |=
            lower - upper >
            kDoubleEpsilon * std::max({1., std::abs(lower), std::abs(upper)});
      }
      if (separated) {
        continue;
      }
      if (Side(a, b, c) * Side(a, b, d) <= 0 &&
          Side(c, d, a) * Side(c, d, b) <= 0) {
        return false;
      }
    }
  }
  return true;
}
static PlanarStraightLineGraph TrimNetwork(const Face& face, const Options& o,
                                           const ArcMap& arc) {
  if (!face.outer) {
    throw Unsupported("Implicit outer boundary with holes is not aligned");
  }
  const auto& s = face.surface;
  if (face.constructed_rectangle) {
    ValidateNaturalRectangle(face);
  }
  auto blocks = KnotBlocks(s, false, true);
  std::array<std::vector<double>, 2> axes;
  for (auto knot_block : blocks) {
    for (int axis = 0; axis < 2; ++axis) {
      for (int endpoint_index = 0; endpoint_index < 2; ++endpoint_index) {
        axes[axis].push_back(knot_block[2 * axis + endpoint_index]);
      }
    }
  }
  for (auto& axis_parameters : axes) {
    axis_parameters = Unique(axis_parameters);
    axis_parameters.erase(axis_parameters.begin());
    axis_parameters.pop_back();
  }
  std::vector<std::vector<Curve>> source{*face.outer};
  source.insert(source.end(), face.inner.begin(), face.inner.end());
  std::vector<std::vector<Vec2>> loops;
  for (std::size_t source_loop_index = 0; source_loop_index < source.size();
       ++source_loop_index) {
    const auto& loop = source[source_loop_index];
    if (loop.empty()) {
      throw std::invalid_argument("Empty trim loop");
    }
    std::vector<Vec2> local;
    std::set<int> junctions;
    for (std::size_t curve_index = 0; curve_index < loop.size();
         ++curve_index) {
      const auto& original = loop[curve_index];
      for (double weight : original.weights) {
        if (!(weight > 0)) {
          throw Unsupported("Nonpositive trim weights");
        }
      }
      auto end = original.Evaluate(original.domain[1], true);
      auto start = loop[(curve_index + 1) % loop.size()].Evaluate(
          loop[(curve_index + 1) % loop.size()].domain[0], true);
      auto delta = Subtract(start, end);
      double gap = std::sqrt((delta[1] * delta[1] + delta[0] * delta[0]) +
                             delta[2] * delta[2]);
      bool retain = gap > 1e-4;
      auto c = PrepareCurve(original);
      std::vector<double> knots{c.domain[0], c.domain[1]};
      for (double knot : Unique(c.knots)) {
        if (knot > c.domain[0] && knot < c.domain[1] &&
            std::count(c.knots.begin(), c.knots.end(), knot) >= c.degree) {
          knots.push_back(knot);
        }
      }
      knots = Unique(knots);
      std::vector<double> parameters;
      for (std::size_t knot_span_index = 0; knot_span_index + 1 < knots.size();
           ++knot_span_index) {
        auto result = Sample(
            [&](double t) {
              auto p = c.Evaluate(t);
              // Retain prepared polynomial UV arithmetic for the network.
              // Only physical samples of the verified B rectangle project.
              return face.constructed_rectangle
                         ? RectangleSamplingPoint(s, {p[0], p[1]})
                         : SamplingPoint(s, {p[0], p[1]});
            },
            knots[knot_span_index], knots[knot_span_index + 1], o.curve.chord,
            o.curve.min_edge, o.curve.max_edge, false, 2 * (c.degree + 1),
            [&](double t) { return c.Evaluate(t); });
        parameters.insert(parameters.end(), result.parameters.begin(),
                          result.parameters.end() - 1);
      }
      parameters.push_back(c.domain[1]);
      std::map<double, std::vector<std::pair<int, double>>> intersections;
      for (int axis = 0; axis < 2; ++axis) {
        for (double knot_parameter : axes[axis]) {
          for (double crossing_parameter :
               Crossings(c, axis, knot_parameter, o)) {
            parameters.push_back(crossing_parameter);
            intersections[crossing_parameter].push_back({axis, knot_parameter});
          }
        }
      }
      parameters = Unique(parameters);
      double factor = parameters.size() < 10
                          ? .01
                          : 1. / (parameters.size() * parameters.size());
      double spacing = (parameters.back() - parameters.front()) * factor;
      std::vector<double> kept{parameters.front()};
      for (std::size_t sample_parameter_index = 1;
           sample_parameter_index + 1 < parameters.size();
           ++sample_parameter_index) {
        double delta = parameters[sample_parameter_index] - kept.back();
        if (!intersections.count(parameters[sample_parameter_index]) &&
            delta < spacing) {
          continue;
        }
        if (delta < 1e-6) {
          continue;
        }
        kept.push_back(parameters[sample_parameter_index]);
      }
      kept.push_back(parameters.back());
      std::vector<Vec2> points;
      for (double sample_parameter : kept) {
        auto xyz = c.Evaluate(sample_parameter);
        Vec2 p{xyz[0], xyz[1]};
        auto hits = intersections.find(sample_parameter);
        if (hits != intersections.end()) {
          for (auto axis_intersection : hits->second) {
            p[axis_intersection.first] = axis_intersection.second;
          }
        }
        for (int axis = 0; axis < 2; ++axis) {
          for (double knot_parameter : axes[axis]) {
            if (std::abs(p[axis] - knot_parameter) < 1e-11) {
              p[axis] = knot_parameter;
            }
          }
        }
        points.push_back(p);
      }
      junctions.insert(static_cast<int>(local.size()));
      local.insert(local.end(), points.begin(),
                   retain ? points.end() : points.end() - 1);
    }
    for (auto& parameter_point : local) {
      if (!Finite(parameter_point)) {
        throw std::invalid_argument("Nonfinite trim loop");
      }
      for (int axis = 0; axis < 2; ++axis) {
        if (!s.closed_axes[axis]) {
          parameter_point[axis] =
              std::clamp(parameter_point[axis], s.domain[2 * axis],
                         s.domain[2 * axis + 1]);
        }
      }
    }
    if (source_loop_index) {
      std::reverse(local.begin(), local.end());
    }
    local = FilterLoop(local, junctions);
    double area = 0.;
    for (std::size_t loop_vertex_index = 0; loop_vertex_index < local.size();
         ++loop_vertex_index) {
      auto a = local[loop_vertex_index];
      auto b = local[(loop_vertex_index + 1) % local.size()];
      area += a[0] * b[1] - a[1] * b[0];
    }
    if (local.size() < 3 || area == 0) {
      throw std::invalid_argument("Degenerate trim loop");
    }
    loops.push_back(std::move(local));
  }
  auto internal = Cleanup(loops, axes);
  for (auto& loop : loops) {
    for (auto& parameter_point : loop) {
      parameter_point = arc.Apply(parameter_point);
    }
  }
  PlanarStraightLineGraph p;
  p.disable_refinement = !Quality(loops);
  for (auto& loop : loops) {
    p.disable_refinement |= !PrepareLoop(loop);
  }
  for (std::size_t loop_index = 0; loop_index < loops.size(); ++loop_index) {
    const auto& loop = loops[loop_index];
    int offset = static_cast<int>(p.uv.size());
    p.uv.insert(p.uv.end(), loop.begin(), loop.end());
    for (int loop_vertex_index = 0;
         loop_vertex_index < static_cast<int>(loop.size());
         ++loop_vertex_index) {
      p.segments.push_back(
          {offset + loop_vertex_index,
           offset + (loop_vertex_index + 1) % static_cast<int>(loop.size())});
      p.markers.push_back(1);
    }
    if (loop_index) {
      auto result = Polygon(loop);
      auto t = result.triangles[result.triangles.size() / 2];
      auto seed = Add(Add(result.uv[t[0]], result.uv[t[1]]), result.uv[t[2]]);
      for (auto& seed_coordinate : seed) {
        seed_coordinate /= 3.;
      }
      p.holes.push_back(seed);
    }
  }
  p.internal_constraints = static_cast<int>(internal.size());
  for (std::size_t internal_segment_index = 0;
       internal_segment_index < internal.size(); ++internal_segment_index) {
    p.segments.push_back(internal[internal_segment_index]);
    p.markers.push_back(static_cast<int>(internal_segment_index) + 2);
  }
  return p;
}
Mesh BuildTrim(const Face& face, const Options& o) {
  const auto& s = face.surface;
  ArcMap arc(s);
  auto p = TrimNetwork(face, o, arc);
  auto result = TriangleRun(p, o, [&](const std::array<Vec2, 3>& points) {
    std::array<Vec2, 3> uv;
    for (int triangle_vertex_index = 0; triangle_vertex_index < 3;
         ++triangle_vertex_index) {
      uv[triangle_vertex_index] =
          arc.Apply(points[triangle_vertex_index], true);
    }
    return Refine(s, o, uv);
  });
  int active_vertex_count = static_cast<int>(result.uv.size());
  std::vector<int> vertex_split_flags(active_vertex_count);
  for (std::size_t segment_index = 0; segment_index < result.segments.size();
       ++segment_index) {
    if (result.markers[segment_index] <= 1) {
      continue;
    }
    auto edge = result.segments[segment_index];
    if (std::abs(result.uv[edge[0]][0] - result.uv[edge[1]][0]) < 1e-5) {
      vertex_split_flags[edge[0]] |= 1;
      vertex_split_flags[edge[1]] |= 1;
    }
    if (std::abs(result.uv[edge[0]][1] - result.uv[edge[1]][1]) < 1e-5) {
      vertex_split_flags[edge[0]] |= 2;
      vertex_split_flags[edge[1]] |= 2;
    }
  }
  int reserved_vertex_count = active_vertex_count;
  for (int vertex_split_flags : vertex_split_flags) {
    reserved_vertex_count +=
        vertex_split_flags == 3 ? 3 : vertex_split_flags != 0;
  }
  result.uv.resize(reserved_vertex_count);
  vertex_split_flags.resize(reserved_vertex_count);
  std::vector<std::vector<int>> adjacent(reserved_vertex_count);
  for (int triangle_index = 0;
       triangle_index < static_cast<int>(result.triangles.size());
       ++triangle_index) {
    for (auto vertex_index : result.triangles[triangle_index]) {
      adjacent[vertex_index].push_back(triangle_index);
    }
  }
  const std::map<int, std::pair<int, int>> transitions{
      {1, {1, 2}}, {3, {10, 11}}, {9, {3, 6}}, {10, {4, 7}}, {11, {5, 8}}};
  for (int axis = 0; axis < 2; ++axis) {
    int vertices_to_split = active_vertex_count;
    for (int vertex_index = 0; vertex_index < vertices_to_split;
         ++vertex_index) {
      int vertex_split_flag = vertex_split_flags[vertex_index];
      if (axis == 0 && vertex_split_flag == 2) {
        vertex_split_flags[vertex_index] = 9;
        continue;
      }
      if (axis == 0 ? (vertex_split_flag != 1 && vertex_split_flag != 3)
                    : (vertex_split_flag != 9 && vertex_split_flag != 10 &&
                       vertex_split_flag != 11)) {
        continue;
      }
      std::vector<int> kept;
      std::vector<int> moved;
      for (int adjacent_triangle_index : adjacent[vertex_index]) {
        auto& row = result.triangles[adjacent_triangle_index];
        double a = result.uv[row[0]][axis];
        double b = result.uv[row[1]][axis];
        double c = result.uv[row[2]][axis];
        double centroid = ((a + b) + c) / 3.;
        auto slot = std::find(row.begin(), row.end(), vertex_index);
        if (result.uv[vertex_index][axis] < centroid && slot != row.end()) {
          *slot = active_vertex_count;
          moved.push_back(adjacent_triangle_index);
        } else {
          kept.push_back(adjacent_triangle_index);
        }
      }
      if (!moved.empty()) {
        adjacent[vertex_index] = std::move(kept);
        adjacent[active_vertex_count] = std::move(moved);
        vertex_split_flags[vertex_index] =
            transitions.at(vertex_split_flag).first;
        vertex_split_flags[active_vertex_count] =
            transitions.at(vertex_split_flag).second;
        result.uv[active_vertex_count] = result.uv[vertex_index];
        ++active_vertex_count;
      } else {
        vertex_split_flags[vertex_index] = 0;
      }
    }
  }
  Mesh out;
  out.face_id = face.id;
  out.identity_mode = "index_components";
  out.branch = "trimmed";
  out.explicit_weights = s.explicit_weights;
  out.triangles = std::move(result.triangles);
  std::vector<Vec2> samples;
  for (int vertex_index = 0; vertex_index < active_vertex_count;
       ++vertex_index) {
    auto uv = arc.Apply(result.uv[vertex_index], true);
    auto sample_uv = uv;
    int normal_side_flag = vertex_split_flags[vertex_index];
    if (normal_side_flag == 1 || normal_side_flag == 4 ||
        normal_side_flag == 7) {
      sample_uv[0] -= 1e-6;
    }
    if (normal_side_flag == 2 || normal_side_flag == 5 ||
        normal_side_flag == 8) {
      sample_uv[0] += 1e-6;
    }
    if (normal_side_flag == 3 || normal_side_flag == 4 ||
        normal_side_flag == 5) {
      sample_uv[1] -= 1e-6;
    }
    if (normal_side_flag == 6 || normal_side_flag == 7 ||
        normal_side_flag == 8) {
      sample_uv[1] += 1e-6;
    }
    samples.push_back(sample_uv);
    out.uv.push_back(uv);
    out.xyz.push_back(s.Evaluate(uv));
  }
  out.normals = MeshNormals(s, samples);
  for (int reserved_vertex_index = active_vertex_count;
       reserved_vertex_index < reserved_vertex_count; ++reserved_vertex_index) {
    out.unpaired_uv_tail.push_back(
        arc.Apply(result.uv[reserved_vertex_index], true));
  }
  for (auto& vertex_position : out.xyz) {
    for (auto& coordinate : vertex_position) {
      coordinate = Quantize(coordinate);
    }
  }
  for (auto& parameter_point : out.uv) {
    for (auto& coordinate : parameter_point) {
      coordinate = Quantize(coordinate);
    }
  }
  for (auto& vertex_normal : out.normals) {
    for (auto& coordinate : vertex_normal) {
      coordinate = Quantize(coordinate);
    }
  }
  for (auto& unpaired_parameter_point : out.unpaired_uv_tail) {
    for (auto& coordinate : unpaired_parameter_point) {
      coordinate = Quantize(coordinate);
    }
  }
  ComponentIdentities(out, face.id);
  return out;
}
}  // namespace iko::internal
