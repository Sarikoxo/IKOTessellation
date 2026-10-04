#include "tessellation/private/bezier.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>
#include <utility>
#include <vector>

#include "tessellation/tessellation.h"
#include "tessellation/private/internal.h"

namespace iko::internal {
using Polynomial = std::vector<Rat>;
static Polynomial Trim(Polynomial polynomial) {
  while (!polynomial.empty() && polynomial.back() == Rat(0)) {
    polynomial.pop_back();
  }
  return polynomial;
}
static Rat Value(const Polynomial& p, Rat x) {
  Rat result;
  for (auto coefficient_iterator = p.rbegin(); coefficient_iterator != p.rend();
       ++coefficient_iterator) {
    result = result * x + *coefficient_iterator;
  }
  return result;
}
static std::pair<Polynomial, Polynomial> Divide(Polynomial dividend,
                                                Polynomial divisor) {
  dividend = Trim(std::move(dividend));
  divisor = Trim(std::move(divisor));
  if (divisor.empty()) {
    throw std::runtime_error("Zero Sturm divisor");
  }
  Polynomial quotient(dividend.size() >= divisor.size()
                          ? dividend.size() - divisor.size() + 1
                          : 0);
  while (dividend.size() >= divisor.size()) {
    int offset = static_cast<int>(dividend.size() - divisor.size());
    Rat coefficient = dividend.back() / divisor.back();
    quotient[offset] = coefficient;
    for (std::size_t divisor_coefficient_index = 0;
         divisor_coefficient_index < divisor.size();
         ++divisor_coefficient_index) {
      dividend[offset + divisor_coefficient_index] -=
          coefficient * divisor[divisor_coefficient_index];
    }
    dividend = Trim(std::move(dividend));
  }
  return {Trim(quotient), dividend};
}
static Polynomial Derivative(const Polynomial& p) {
  Polynomial out;
  for (std::size_t coefficient_index = 1; coefficient_index < p.size();
       ++coefficient_index) {
    out.push_back(Rat(static_cast<int>(coefficient_index)) *
                  p[coefficient_index]);
  }
  return out;
}
static std::vector<Polynomial> Sturm(const Polynomial& p) {
  std::vector<Polynomial> sturm_sequence{p, Derivative(p)};
  while (sturm_sequence.back().size() > 1) {
    auto remainder =
        Divide(sturm_sequence[sturm_sequence.size() - 2], sturm_sequence.back())
            .second;
    if (remainder.empty()) {
      throw Unsupported("Multiple line/cubic root");
    }
    for (auto& remainder_coefficient : remainder) {
      remainder_coefficient = -remainder_coefficient;
    }
    sturm_sequence.push_back(std::move(remainder));
  }
  return sturm_sequence;
}
static std::vector<double> UnitRoots(Polynomial p) {
  p = Trim(std::move(p));
  if (p.size() > 4) {
    throw Unsupported("Mixed polynomial degree exceeds cubic");
  }
  if (p.empty()) {
    throw Unsupported("Overlapping line/cubic");
  }
  if (p.size() == 1) {
    return {};
  }
  auto der = Derivative(p);
  Sturm(p);
  std::vector<double> roots;
  for (int unit_endpoint : {0, 1}) {
    if (Value(p, Rat(unit_endpoint)) == Rat(0)) {
      roots.push_back(static_cast<double>(unit_endpoint));
      p = Divide(p, {-Rat(unit_endpoint), Rat(1)}).first;
    }
  }
  if (p.size() > 1) {
    auto sturm_sequence = Sturm(p);
    auto variations = [&](Rat x) {
      std::vector<int> signs;
      for (const auto& sturm_polynomial : sturm_sequence) {
        auto v = Value(sturm_polynomial, x);
        if (v != Rat(0)) {
          signs.push_back(v > Rat(0) ? 1 : -1);
        }
      }
      int n = 0;
      for (std::size_t sign_index = 1; sign_index < signs.size();
           ++sign_index) {
        n += signs[sign_index] != signs[sign_index - 1];
      }
      return n;
    };
    auto count = [&](Rat a, Rat b) { return variations(a) - variations(b); };
    struct Interval {
      Rat lo;
      Rat hi;
      int count;
      int depth;
    };
    std::vector<Interval> pending{{Rat(0), Rat(1), count(Rat(0), Rat(1)), 0}};
    int nodes = 0;
    while (!pending.empty()) {
      auto row = std::move(pending.back());
      pending.pop_back();
      ++nodes;
      if (!row.count) {
        continue;
      }
      if (row.depth >= 64 || nodes > 512) {
        throw Unsupported("Unresolved cubic root bounds");
      }
      Rat middle = (row.lo + row.hi) / Rat(2);
      if (row.count == 1 &&
          (Value(p, middle) == Rat(0) || (row.hi - row.lo).Value() < 1e-16)) {
        roots.push_back(middle.Value());
        continue;
      }
      int left = count(row.lo, middle);
      pending.push_back({middle, row.hi, row.count - left, row.depth + 1});
      pending.push_back({row.lo, middle, left, row.depth + 1});
    }
  }
  std::sort(roots.begin(), roots.end());
  for (std::size_t root_index = 1; root_index < roots.size(); ++root_index) {
    if (roots[root_index] - roots[root_index - 1] < 1e-8) {
      throw Unsupported("Nearly multiple line/cubic roots");
    }
  }
  double scale = 0;
  for (auto derivative_coefficient : der) {
    scale = std::max(scale, std::abs(derivative_coefficient.Value()));
  }
  if (!std::isfinite(scale) || scale == 0) {
    throw Unsupported("Unresolved cubic derivative scale");
  }
  for (double root_parameter : roots) {
    if (std::abs(Value(der, Rat(root_parameter)).Value()) < 1e-10 * scale) {
      throw Unsupported("Nearly tangential line/cubic root");
    }
  }
  return roots;
}
static std::pair<ExactPoint, ExactPoint> LineData(const Curve& c) {
  ExactBezier verified(c);
  if (c.degree != 1 ||
      !std::all_of(c.weights.begin(), c.weights.end(),
                   [](double weight) { return weight == 1.; })) {
    throw Unsupported("Polynomial full-span line required for mixed contacts");
  }
  ExactPoint origin{Rat(c.control[0][0]), Rat(c.control[0][1])};
  ExactPoint delta{Rat(c.control[1][0]) - origin[0],
                   Rat(c.control[1][1]) - origin[1]};
  if (delta[0] == Rat(0) && delta[1] == Rat(0)) {
    throw Unsupported("Zero line in mixed contacts");
  }
  return {origin, delta};
}
static Rat Determinant(ExactPoint x, ExactPoint y) {
  return x[0] * y[1] - x[1] * y[0];
}
std::vector<std::array<double, 2>> LineContacts(const Curve& line,
                                                const Curve& curve) {
  auto data = LineData(line);
  ExactBezier verified(curve);
  std::vector<std::array<double, 2>> contacts;
  if (curve.degree == 1) {
    auto other = LineData(curve);
    ExactPoint difference{other.first[0] - data.first[0],
                          other.first[1] - data.first[1]};
    auto det = Determinant(data.second, other.second);
    if (det == Rat(0)) {
      if (Determinant(difference, data.second) == Rat(0)) {
        throw Unsupported(
            "Collinear mixed lines require overlap classification");
      }
      return {};
    }
    double u = (Determinant(difference, other.second) / det).Value();
    double v = (Determinant(difference, data.second) / det).Value();
    if (u < -1e-10 || u > 1. + 1e-10 || v < -1e-10 || v > 1. + 1e-10) {
      return {};
    }
    Vec2 la{line.control[1][0] - line.control[0][0],
            line.control[1][1] - line.control[0][1]};
    Vec2 lb{curve.control[1][0] - curve.control[0][0],
            curve.control[1][1] - curve.control[0][1]};
    if (std::abs(det.Value()) <= 1e-12 * (Norm(la) * Norm(lb))) {
      throw Unsupported("Nearly parallel mixed line contact");
    }
    return {{{line.domain[0] +
                  std::clamp(u, 0., 1.) * (line.domain[1] - line.domain[0]),
              curve.domain[0] + std::clamp(v, 0., 1.) *
                                    (curve.domain[1] - curve.domain[0])}}};
  }
  std::vector<Rat> bernstein;
  for (std::size_t control_point_index = 0;
       control_point_index < curve.control.size(); ++control_point_index) {
    ExactPoint delta{
        Rat(curve.control[control_point_index][0]) - data.first[0],
        Rat(curve.control[control_point_index][1]) - data.first[1]};
    bernstein.push_back(Rat(curve.weights[control_point_index]) *
                        Determinant(data.second, delta));
  }
  std::vector<double> roots;
  if (curve.degree == 2) {
    auto c = bernstein[0];
    auto b = Rat(2) * (bernstein[1] - bernstein[0]);
    auto a = bernstein[0] - Rat(2) * bernstein[1] + bernstein[2];
    if (a == Rat(0)) {
      if (b == Rat(0)) {
        if (c == Rat(0)) {
          throw Unsupported("Overlapping line/quadratic mixed cut");
        }
        return {};
      }
      roots = {(-c / b).Value()};
    } else {
      Rat disc = b * b - Rat(4) * a * c;
      if (disc < Rat(0)) {
        return {};
      }
      double scale = std::max({std::abs((b * b).Value()),
                               std::abs((Rat(4) * a * c).Value()), 1e-300});
      if (disc == Rat(0) || disc.Value() <= 1e-12 * scale) {
        throw Unsupported("Tangential or nearly multiple mixed cut");
      }
      double q =
          -.5 * (b.Value() + std::copysign(std::sqrt(disc.Value()), b.Value()));
      roots = {q / a.Value(), c.Value() / q};
    }
    std::sort(roots.begin(), roots.end());
  } else if (curve.degree == 3) {
    auto a = bernstein[0];
    auto b = bernstein[1];
    auto c = bernstein[2];
    auto d = bernstein[3];
    roots = UnitRoots({a, Rat(3) * (b - a), Rat(3) * (a - Rat(2) * b + c),
                       -a + Rat(3) * b - Rat(3) * c + d});
  } else {
    throw Unsupported(
        "Full positive rational cubic required for mixed contacts");
  }
  Vec2 delta{line.control[1][0] - line.control[0][0],
             line.control[1][1] - line.control[0][1]};
  for (double curve_fraction : roots) {
    if (curve_fraction < -1e-10 || curve_fraction > 1. + 1e-10) {
      continue;
    }
    curve_fraction = std::clamp(curve_fraction, 0., 1.);
    double t =
        curve.domain[0] + curve_fraction * (curve.domain[1] - curve.domain[0]);
    auto point = curve.Evaluate(t);
    Vec2 difference{point[0] - line.control[0][0],
                    point[1] - line.control[0][1]};
    double v = Dot(difference, delta) / Dot(delta, delta);
    if (v >= -1e-10 && v <= 1. + 1e-10) {
      v = std::clamp(v, 0., 1.);
      contacts.push_back(
          {line.domain[0] + v * (line.domain[1] - line.domain[0]), t});
    }
  }
  return contacts;
}
}  // namespace iko::internal
