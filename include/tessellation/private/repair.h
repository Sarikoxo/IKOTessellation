#ifndef IKO_SRC_REPAIR_H_
#define IKO_SRC_REPAIR_H_

#include <array>
#include <string>
#include <vector>

#include "tessellation/tessellation.h"
#include "tessellation/private/internal.h"

namespace iko::internal {
struct WireItem {
  Curve curve;
  int source = -1;
  int direction = 1;
  std::array<Vec3, 2> endpoints;
};
struct RepairResult {
  std::vector<WireItem> wire;
  bool clockwise = false;
};
Curve ReverseCurve(Curve curve);
RepairResult ConnectRejected(const Face& face, double rejection_threshold,
                             std::vector<std::string>& trace);
RepairResult RepairBending(const Surface& surface,
                           const std::vector<WireItem>& wire,
                           std::vector<std::string>& trace);
RepairResult RepairCombined(const Surface& surface,
                            const std::vector<WireItem>& wire,
                            std::vector<std::string>& trace);
std::array<double, 2> ResolutionBounds(const Surface& surface);
double CandidateTolerance(const Surface& surface, const Curve& original,
                          const Curve& candidate);
std::vector<Vec2> SampledWire(const std::vector<Curve>& curves, int count);
double SignedArea(const std::vector<Vec2>& points);
}  // namespace iko::internal

#endif  // IKO_SRC_REPAIR_H_
