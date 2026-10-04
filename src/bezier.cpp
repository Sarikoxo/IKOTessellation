#include "tessellation/private/bezier.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <set>
#include <utility>
#include <vector>

#include "tessellation/private/exact.h"
#include "tessellation/tessellation.h"
#include "tessellation/private/internal.h"

namespace iko::internal {
static int Comb(int population_size, int selection_count) {
  int c = 1;
  for (int factor_index = 1; factor_index <= selection_count; ++factor_index) {
    c = c * (population_size - factor_index + 1) / factor_index;
  }
  return c;
}
static std::vector<Rat> Product(const std::vector<Rat>& a,
                                const std::vector<Rat>& b) {
  int m = static_cast<int>(a.size()) - 1;
  int n = static_cast<int>(b.size()) - 1;
  std::vector<Rat> out(m + n + 1);
  for (int left_coefficient_index = 0; left_coefficient_index <= m;
       ++left_coefficient_index) {
    for (int right_coefficient_index = 0; right_coefficient_index <= n;
         ++right_coefficient_index) {
      out[left_coefficient_index + right_coefficient_index] +=
          a[left_coefficient_index] * b[right_coefficient_index] *
          Rat(Comb(m, left_coefficient_index) *
                  Comb(n, right_coefficient_index),
              Comb(m + n, left_coefficient_index + right_coefficient_index));
    }
  }
  return out;
}
ExactBezier::ExactBezier(const Curve& c) {
  int d = c.degree;
  if (c.weights.size() != static_cast<std::size_t>(d + 1) ||
      c.knots.size() != static_cast<std::size_t>(2 * (d + 1)) ||
      c.domain[0] >= c.domain[1]) {
    throw Unsupported(
        "Full-domain single Bezier span required; no silent domain clipping");
  }
  for (int knot_index = 0; knot_index <= d; ++knot_index) {
    if (c.knots[knot_index] != c.domain[0] ||
        c.knots[d + 1 + knot_index] != c.domain[1]) {
      throw Unsupported(
          "Full-domain single Bezier span required; no silent domain "
          "clipping");
    }
  }
  for (int control_point_index = 0; control_point_index <= d;
       ++control_point_index) {
    if (c.weights[control_point_index] <= 0) {
      throw Unsupported(
          "Positive homogeneous coefficients required for these bounds");
    }
    Rat w(c.weights[control_point_index]);
    h.push_back({Rat(c.control[control_point_index][0]) * w,
                 Rat(c.control[control_point_index][1]) * w,
                 Rat(c.control[control_point_index][2]) * w, w});
  }
  range = {Rat(c.domain[0]), Rat(c.domain[1])};
}
std::vector<ExactPoint> ExactBezier::Control() const {
  std::vector<ExactPoint> out;
  for (const auto& homogeneous_control_point : h) {
    out.push_back(
        {homogeneous_control_point[0] / homogeneous_control_point[3],
         homogeneous_control_point[1] / homogeneous_control_point[3]});
  }
  return out;
}
std::array<ExactRange, 2> ExactBezier::Box() const {
  auto cp = Control();
  std::array<ExactRange, 2> box{{{cp[0][0], cp[0][0]}, {cp[0][1], cp[0][1]}}};
  for (const auto& control_point : cp) {
    for (int axis = 0; axis < 2; ++axis) {
      box[axis][0] = std::min(box[axis][0], control_point[axis]);
      box[axis][1] = std::max(box[axis][1], control_point[axis]);
    }
  }
  return box;
}
std::pair<ExactBezier, ExactBezier> ExactBezier::Split() const {
  auto level = h;
  std::vector<std::array<Rat, 4>> left;
  std::vector<std::array<Rat, 4>> right;
  while (!level.empty()) {
    left.push_back(level.front());
    right.push_back(level.back());
    std::vector<std::array<Rat, 4>> next;
    for (std::size_t control_point_index = 0;
         control_point_index + 1 < level.size(); ++control_point_index) {
      std::array<Rat, 4> p;
      for (int homogeneous_channel = 0; homogeneous_channel < 4;
           ++homogeneous_channel) {
        p[homogeneous_channel] =
            (level[control_point_index][homogeneous_channel] +
             level[control_point_index + 1][homogeneous_channel]) /
            Rat(2);
      }
      next.push_back(p);
    }
    level = std::move(next);
  }
  std::reverse(right.begin(), right.end());
  Rat mid = (range[0] + range[1]) / Rat(2);
  return {ExactBezier(left, {range[0], mid}),
          ExactBezier(right, {mid, range[1]})};
}
static std::vector<Rat> DerivativeNumerator(const ExactBezier& c, int axis) {
  int d = c.Degree();
  std::vector<Rat> x;
  std::vector<Rat> w;
  std::vector<Rat> dx;
  std::vector<Rat> dw;
  for (auto homogeneous_control_point : c.h) {
    x.push_back(homogeneous_control_point[axis]);
    w.push_back(homogeneous_control_point[3]);
  }
  for (int coefficient_index = 0; coefficient_index < d; ++coefficient_index) {
    dx.push_back(Rat(d) * (x[coefficient_index + 1] - x[coefficient_index]));
    dw.push_back(Rat(d) * (w[coefficient_index + 1] - w[coefficient_index]));
  }
  auto a = Product(dx, w);
  auto b = Product(x, dw);
  for (std::size_t coefficient_index = 0; coefficient_index < a.size();
       ++coefficient_index) {
    a[coefficient_index] -= b[coefficient_index];
  }
  return a;
}
bool ExactBezier::Monotone() const {
  for (int axis = 0; axis < 2; ++axis) {
    auto n = DerivativeNumerator(*this, axis);
    bool nonnegative = true;
    bool nonpositive = true;
    bool positive = false;
    bool negative = false;
    for (auto derivative_coefficient : n) {
      nonnegative &= derivative_coefficient >= Rat(0);
      nonpositive &= derivative_coefficient <= Rat(0);
      positive |= derivative_coefficient > Rat(0);
      negative |= derivative_coefficient < Rat(0);
    }
    if ((nonnegative && positive) || (nonpositive && negative)) {
      return true;
    }
  }
  return false;
}
static std::vector<ExactContact> Contacts(const ExactBezier& a,
                                          const ExactBezier& b) {
  auto ac = a.Control();
  auto bc = b.Control();
  std::vector<ExactContact> out;
  for (int first_endpoint_index : {0, a.Degree()}) {
    for (int second_endpoint_index : {0, b.Degree()}) {
      if (ac[first_endpoint_index] == bc[second_endpoint_index]) {
        out.push_back({a.range[first_endpoint_index == 0 ? 0 : 1],
                       b.range[second_endpoint_index == 0 ? 0 : 1]});
      }
    }
  }
  return out;
}
static bool Member(const std::vector<ExactContact>& allowed,
                   const ExactContact& p) {
  return std::find(allowed.begin(), allowed.end(), p) != allowed.end();
}
static bool Separated(const ExactBezier& a, const ExactBezier& b,
                      const std::vector<ExactContact>& allowed) {
  auto ba = a.Box();
  auto bb = b.Box();
  for (int axis = 0; axis < 2; ++axis) {
    if (ba[axis][1] < bb[axis][0] || bb[axis][1] < ba[axis][0]) {
      return true;
    }
  }
  if (allowed.empty()) {
    return false;
  }
  std::array<ExactRange, 2> overlap;
  for (int axis = 0; axis < 2; ++axis) {
    overlap[axis] = {std::max(ba[axis][0], bb[axis][0]),
                     std::min(ba[axis][1], bb[axis][1])};
  }
  if (overlap[0][0] == overlap[0][1] && overlap[1][0] == overlap[1][1] &&
      a.Monotone() && b.Monotone()) {
    ExactPoint point{overlap[0][0], overlap[1][0]};
    auto ac = a.Control();
    auto bc = b.Control();
    std::vector<ExactContact> found;
    for (int first_endpoint_index : {0, a.Degree()}) {
      for (int second_endpoint_index : {0, b.Degree()}) {
        if (ac[first_endpoint_index] == point &&
            bc[second_endpoint_index] == point) {
          found.push_back({a.range[first_endpoint_index == 0 ? 0 : 1],
                           b.range[second_endpoint_index == 0 ? 0 : 1]});
        }
      }
    }
    if (!found.empty() && std::all_of(found.begin(), found.end(),
                                      [&](const auto& endpoint_contact) {
                                        return Member(allowed,
                                                      endpoint_contact);
                                      })) {
      return true;
    }
  }
  for (int axis = 0; axis < 2; ++axis) {
    if ((ba[axis][1] == bb[axis][0] || bb[axis][1] == ba[axis][0]) &&
        ba[axis][0] < ba[axis][1] && bb[axis][0] < bb[axis][1]) {
      auto found = Contacts(a, b);
      if (std::all_of(found.begin(), found.end(),
                      [&](const auto& endpoint_contact) {
                        return Member(allowed, endpoint_contact);
                      })) {
        return true;
      }
    }
  }
  return false;
}
static Rat Width(const ExactBezier& a) {
  auto box = a.Box();
  return std::max(box[0][1] - box[0][0], box[1][1] - box[1][0]);
}
struct PairNode {
  ExactBezier a;
  ExactBezier b;
  int depth;
};
bool CertifyDisjoint(const Curve& ca, const Curve& cb, bool allow) {
  ExactBezier a(ca);
  ExactBezier b(cb);
  auto allowed = allow ? Contacts(a, b) : std::vector<ExactContact>{};
  std::vector<PairNode> pending{{a, b, 0}};
  int nodes = 0;
  while (!pending.empty() && nodes < 4096) {
    auto row = std::move(pending.back());
    pending.pop_back();
    ++nodes;
    if (Separated(row.a, row.b, allowed)) {
      continue;
    }
    if (row.depth >= 20) {
      return false;
    }
    if (Width(row.a) >= Width(row.b)) {
      auto parts = row.a.Split();
      pending.push_back({std::move(parts.first), row.b, row.depth + 1});
      pending.push_back({std::move(parts.second), row.b, row.depth + 1});
    } else {
      auto parts = row.b.Split();
      pending.push_back({row.a, std::move(parts.first), row.depth + 1});
      pending.push_back({row.a, std::move(parts.second), row.depth + 1});
    }
  }
  return pending.empty();
}
static ExactRange IntervalAdd(const ExactRange& a, const ExactRange& b) {
  return {a[0] + b[0], a[1] + b[1]};
}
static ExactRange IntervalNeg(const ExactRange& a) {
  return {-a[1], -a[0]};
}
static ExactRange IntervalMul(const ExactRange& a, const ExactRange& b) {
  Rat values[] = {a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]};
  return {*std::min_element(values, values + 4),
          *std::max_element(values, values + 4)};
}
static ExactRange Constant(Rat a) {
  return {a, a};
}
static std::pair<ExactPoint, ExactPoint> Middle(const ExactBezier& c) {
  auto row = c.h;
  std::array<Rat, 4> derivative{};
  while (row.size() > 1) {
    if (row.size() == 2) {
      for (int homogeneous_channel = 0; homogeneous_channel < 4;
           ++homogeneous_channel) {
        derivative[homogeneous_channel] =
            Rat(c.Degree()) *
            (row[1][homogeneous_channel] - row[0][homogeneous_channel]);
      }
    }
    std::vector<std::array<Rat, 4>> next;
    for (std::size_t control_point_index = 0;
         control_point_index + 1 < row.size(); ++control_point_index) {
      std::array<Rat, 4> p;
      for (int homogeneous_channel = 0; homogeneous_channel < 4;
           ++homogeneous_channel) {
        p[homogeneous_channel] =
            (row[control_point_index][homogeneous_channel] +
             row[control_point_index + 1][homogeneous_channel]) /
            Rat(2);
      }
      next.push_back(p);
    }
    row = std::move(next);
  }
  auto h = row[0];
  ExactPoint point;
  ExactPoint der;
  for (int axis = 0; axis < 2; ++axis) {
    point[axis] = h[axis] / h[3];
    der[axis] =
        (derivative[axis] * h[3] - h[axis] * derivative[3]) / (h[3] * h[3]);
  }
  return {point, der};
}
static std::array<ExactRange, 2> DerivativeBox(const ExactBezier& c) {
  Rat wmin = c.h[0][3];
  Rat wmax = wmin;
  for (auto homogeneous_control_point : c.h) {
    wmin = std::min(wmin, homogeneous_control_point[3]);
    wmax = std::max(wmax, homogeneous_control_point[3]);
  }
  ExactRange reciprocal{Rat(1) / (wmax * wmax), Rat(1) / (wmin * wmin)};
  std::array<ExactRange, 2> out;
  for (int axis = 0; axis < 2; ++axis) {
    auto n = DerivativeNumerator(c, axis);
    out[axis] = IntervalMul({*std::min_element(n.begin(), n.end()),
                             *std::max_element(n.begin(), n.end())},
                            reciprocal);
  }
  return out;
}
static int RootBox(const ExactBezier& a, const ExactBezier& b,
                   const std::vector<ExactContact>& allowed) {
  auto ma = Middle(a);
  auto mb = Middle(b);
  std::array<std::array<Rat, 2>, 2> jacobian{
      {{ma.second[0], -mb.second[0]}, {ma.second[1], -mb.second[1]}}};
  Rat determinant =
      jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];
  if (determinant == Rat(0)) {
    return 0;
  }
  std::array<std::array<Rat, 2>, 2> r{
      {{jacobian[1][1] / determinant, -jacobian[0][1] / determinant},
       {-jacobian[1][0] / determinant, jacobian[0][0] / determinant}}};
  auto ba = DerivativeBox(a);
  auto bb = DerivativeBox(b);
  std::array<std::array<ExactRange, 2>, 2> ji{
      {{ba[0], IntervalNeg(bb[0])}, {ba[1], IntervalNeg(bb[1])}}};
  std::array<std::array<ExactRange, 2>, 2> e;
  Rat contraction;
  for (int matrix_row = 0; matrix_row < 2; ++matrix_row) {
    Rat sum;
    for (int matrix_column = 0; matrix_column < 2; ++matrix_column) {
      e[matrix_row][matrix_column] = IntervalAdd(
          Constant(Rat(matrix_row == matrix_column)),
          IntervalNeg(IntervalAdd(
              IntervalMul(Constant(r[matrix_row][0]), ji[0][matrix_column]),
              IntervalMul(Constant(r[matrix_row][1]), ji[1][matrix_column]))));
      sum += std::max(exact::Abs(e[matrix_row][matrix_column][0]),
                      exact::Abs(e[matrix_row][matrix_column][1]));
    }
    contraction = std::max(contraction, sum);
  }
  if (contraction >= Rat(1)) {
    return 0;
  }
  for (auto endpoint_contact : Contacts(a, b)) {
    if (Member(allowed, endpoint_contact)) {
      return 1;
    }
  }
  ExactPoint f{ma.first[0] - mb.first[0], ma.first[1] - mb.first[1]};
  for (int parameter_axis = 0; parameter_axis < 2; ++parameter_axis) {
    Rat center =
        Rat(1, 2) - (r[parameter_axis][0] * f[0] + r[parameter_axis][1] * f[1]);
    auto image = IntervalAdd(
        Constant(center),
        IntervalAdd(
            IntervalMul(e[parameter_axis][0], {Rat(-1, 2), Rat(1, 2)}),
            IntervalMul(e[parameter_axis][1], {Rat(-1, 2), Rat(1, 2)})));
    if (!(Rat(0) < image[0] && image[0] <= image[1] && image[1] < Rat(1))) {
      return 0;
    }
  }
  return 2;
}
static std::pair<std::vector<std::array<ExactRange, 2>>, int> PairRoots(
    const ExactBezier& a, const ExactBezier& b,
    const std::vector<ExactContact>& allowed) {
  std::vector<PairNode> pending{{a, b, 0}};
  int nodes = 0;
  std::vector<std::array<ExactRange, 2>> roots;
  while (!pending.empty() && nodes < 4096) {
    auto row = std::move(pending.back());
    pending.pop_back();
    ++nodes;
    if (Separated(row.a, row.b, allowed)) {
      continue;
    }
    int proof = RootBox(row.a, row.b, allowed);
    if (proof) {
      if (proof == 2) {
        roots.push_back({row.a.range, row.b.range});
      }
      continue;
    }
    if (row.depth >= 28) {
      throw Unsupported(
          "Bent Bezier contact is unresolved within the subdivision bound");
    }
    if (Width(row.a) >= Width(row.b)) {
      auto parts = row.a.Split();
      pending.push_back({parts.first, row.b, row.depth + 1});
      pending.push_back({parts.second, row.b, row.depth + 1});
    } else {
      auto parts = row.b.Split();
      pending.push_back({row.a, parts.first, row.depth + 1});
      pending.push_back({row.a, parts.second, row.depth + 1});
    }
  }
  if (!pending.empty()) {
    throw Unsupported("Bent Bezier contact exceeded the work bound");
  }
  return {roots, nodes};
}
static void Simple(const ExactBezier& c) {
  std::vector<std::pair<ExactBezier, int>> pending{{c, 0}};
  std::vector<ExactBezier> pieces;
  int nodes = 0;
  while (!pending.empty()) {
    auto row = std::move(pending.back());
    pending.pop_back();
    if (row.first.Monotone()) {
      pieces.push_back(std::move(row.first));
      if (pieces.size() > 32) {
        throw Unsupported("Bent simplicity exceeded the piece bound");
      }
      continue;
    }
    if (row.second >= 8) {
      throw Unsupported("Bent Bezier simplicity is unresolved");
    }
    auto parts = row.first.Split();
    pending.push_back({std::move(parts.second), row.second + 1});
    pending.push_back({std::move(parts.first), row.second + 1});
  }
  for (std::size_t first_piece_index = 0; first_piece_index < pieces.size();
       ++first_piece_index) {
    for (std::size_t second_piece_index = first_piece_index + 1;
         second_piece_index < pieces.size(); ++second_piece_index) {
      std::vector<ExactContact> allowed;
      if (pieces[first_piece_index].range[1] ==
          pieces[second_piece_index].range[0]) {
        allowed.push_back({pieces[first_piece_index].range[1],
                           pieces[second_piece_index].range[0]});
      }
      auto result = PairRoots(pieces[first_piece_index],
                              pieces[second_piece_index], allowed);
      nodes += result.second;
      if (!result.first.empty()) {
        throw Unsupported("Self-intersecting individual bent Bezier curve");
      }
      if (nodes > 8192) {
        throw Unsupported("Bent simplicity exceeded the work bound");
      }
    }
  }
}
CrossingCertificate CertifyCrossingWire(const std::vector<Curve>& curves) {
  std::vector<ExactBezier> exact;
  for (const auto& curve : curves) {
    if (curve.degree < 1 || curve.degree > 3) {
      throw Unsupported(
          "Bent crossing verification is limited to degrees 1 through 3");
    }
    exact.emplace_back(curve);
  }
  for (std::size_t curve_index = 0; curve_index < exact.size(); ++curve_index) {
    if (exact[curve_index].Control().back() !=
        exact[(curve_index + 1) % exact.size()].Control().front()) {
      throw Unsupported("Bent crossing wire is not exactly joined");
    }
  }
  for (auto& exact_curve : exact) {
    Simple(exact_curve);
  }
  CrossingCertificate certificate;
  for (std::size_t first_curve_index = 0; first_curve_index < exact.size();
       ++first_curve_index) {
    for (std::size_t second_curve_index = first_curve_index + 1;
         second_curve_index < exact.size(); ++second_curve_index) {
      std::vector<ExactContact> allowed;
      if (second_curve_index == first_curve_index + 1) {
        allowed.push_back({exact[first_curve_index].range[1],
                           exact[second_curve_index].range[0]});
      }
      if (first_curve_index == 0 && second_curve_index == exact.size() - 1) {
        allowed.push_back({exact[first_curve_index].range[0],
                           exact[second_curve_index].range[1]});
      }
      auto found =
          Contacts(exact[first_curve_index], exact[second_curve_index]);
      for (auto found_contact : found) {
        if (!Member(allowed, found_contact)) {
          throw Unsupported(
              "Unexpected coincident endpoints in bent crossing wire");
        }
      }
      for (auto allowed_contact : allowed) {
        if (!Member(found, allowed_contact)) {
          throw Unsupported(
              "Unexpected coincident endpoints in bent crossing wire");
        }
      }
      auto roots = PairRoots(exact[first_curve_index],
                             exact[second_curve_index], allowed);
      for (const auto& root_ranges : roots.first) {
        certificate.crossings.push_back({{static_cast<int>(first_curve_index),
                                          static_cast<int>(second_curve_index)},
                                         root_ranges});
      }
    }
  }
  if (certificate.crossings.empty()) {
    throw Unsupported(
        "Sampled bent crossing has no certified whole-curve counterpart");
  }
  return certificate;
}
static Rat ExactCross(ExactPoint a, ExactPoint b) {
  return a[0] * b[1] - a[1] * b[0];
}
void MatchCrossings(const std::vector<Curve>& curves, int count,
                    const CrossingCertificate& certificate) {
  std::vector<Vec2> points;
  std::vector<std::vector<double>> times;
  for (const auto& curve : curves) {
    std::vector<double> ts;
    double step = (curve.domain[1] - curve.domain[0]) / (count - 1);
    for (int sample_index = 0; sample_index < count; ++sample_index) {
      ts.push_back(sample_index == count - 1
                       ? curve.domain[1]
                       : curve.domain[0] + sample_index * step);
    }
    for (int sample_index = 0; sample_index < count - 1; ++sample_index) {
      auto p = curve.Evaluate(ts[sample_index]);
      points.push_back({p[0], p[1]});
    }
    times.push_back(ts);
  }
  std::set<int> seen;
  for (std::size_t first_segment_index = 0; first_segment_index < points.size();
       ++first_segment_index) {
    auto a = points[first_segment_index];
    auto b = points[(first_segment_index + 1) % points.size()];
    for (std::size_t second_segment_index = first_segment_index + 2;
         second_segment_index < points.size() - (first_segment_index == 0);
         ++second_segment_index) {
      auto c = points[second_segment_index];
      auto d = points[(second_segment_index + 1) % points.size()];
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
      auto sides = SegmentSides(a, b, c, d);
      if (sides[0] * sides[1] > 0 || sides[2] * sides[3] > 0) {
        continue;
      }
      if (std::any_of(sides.begin(), sides.end(), [](int orientation_sign) {
            return orientation_sign == 0;
          })) {
        throw Unsupported(
            "Bent sampled boundary has a touching or overlapping contact");
      }
      ExactPoint ab{Rat(b[0]) - Rat(a[0]), Rat(b[1]) - Rat(a[1])};
      ExactPoint cd{Rat(d[0]) - Rat(c[0]), Rat(d[1]) - Rat(c[1])};
      ExactPoint ca{Rat(c[0]) - Rat(a[0]), Rat(c[1]) - Rat(a[1])};
      Rat den = ExactCross(ab, cd);
      Rat alpha = ExactCross(ca, cd) / den;
      Rat beta = ExactCross(ca, ab) / den;
      std::array<int, 2> pair{
          static_cast<int>(first_segment_index) / (count - 1),
          static_cast<int>(second_segment_index) / (count - 1)};
      int k1 = static_cast<int>(first_segment_index) % (count - 1);
      int k2 = static_cast<int>(second_segment_index) % (count - 1);
      std::array<Rat, 2> params{
          Rat(times[pair[0]][k1]) +
              alpha * (Rat(times[pair[0]][k1 + 1]) - Rat(times[pair[0]][k1])),
          Rat(times[pair[1]][k2]) +
              beta * (Rat(times[pair[1]][k2 + 1]) - Rat(times[pair[1]][k2]))};
      std::vector<int> matches;
      for (int crossing_index = 0;
           crossing_index < static_cast<int>(certificate.crossings.size());
           ++crossing_index) {
        const auto& r = certificate.crossings[crossing_index];
        if (r.curves == pair && r.ranges[0][0] < params[0] &&
            params[0] < r.ranges[0][1] && r.ranges[1][0] < params[1] &&
            params[1] < r.ranges[1][1]) {
          matches.push_back(crossing_index);
        }
      }
      if (matches.size() != 1 || seen.count(matches[0])) {
        throw Unsupported(
            "Bent sampled crossing does not uniquely match its certified "
            "parameter box");
      }
      seen.insert(matches[0]);
    }
  }
  if (seen.size() != certificate.crossings.size()) {
    throw Unsupported("Bent sampling missed a certified crossing");
  }
}
}  // namespace iko::internal
