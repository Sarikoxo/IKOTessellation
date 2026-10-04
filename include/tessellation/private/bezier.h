#ifndef IKO_SRC_BEZIER_H_
#define IKO_SRC_BEZIER_H_

#include <array>
#include <utility>
#include <vector>

#include "tessellation/private/exact.h"
#include "tessellation/tessellation.h"
#include "tessellation/private/internal.h"

namespace iko::internal {
using exact::Rat;
using ExactPoint = exact::Point;
using ExactRange = exact::Range;
using ExactContact = exact::Contact;
struct ExactBezier {
  std::vector<std::array<Rat, 4>> h;
  ExactRange range;
  explicit ExactBezier(const Curve& curve);
  ExactBezier(std::vector<std::array<Rat, 4>> net, ExactRange r)
      : h(std::move(net)), range(std::move(r)) {}
  int Degree() const {
    return static_cast<int>(h.size()) - 1;
  }
  std::vector<ExactPoint> Control() const;
  std::array<ExactRange, 2> Box() const;
  std::pair<ExactBezier, ExactBezier> Split() const;
  bool Monotone() const;
};
bool CertifyDisjoint(const Curve& a, const Curve& b,
                     bool allow_endpoints = true);
struct CrossingBox {
  std::array<int, 2> curves;
  std::array<ExactRange, 2> ranges;
};
struct CrossingCertificate {
  std::vector<CrossingBox> crossings;
};
CrossingCertificate CertifyCrossingWire(const std::vector<Curve>& curves);
void MatchCrossings(const std::vector<Curve>& curves, int count,
                    const CrossingCertificate& certificate);
std::vector<std::array<double, 2>> LineContacts(const Curve& line,
                                                const Curve& curve);
std::array<int, 4> SegmentSides(Vec2 a, Vec2 b, Vec2 c, Vec2 d);
bool OuterIntersects(const std::vector<Vec2>& points);
}  // namespace iko::internal

#endif  // IKO_SRC_BEZIER_H_
