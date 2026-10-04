#include "tessellation/private/internal.h"

namespace iko::internal {
// Shared argument validation for the DLL and CLI; no geometry dependencies.
Options Effective(Options o) {
  if (o.surface_de < 0) {
    throw std::invalid_argument("Invalid surface selection");
  }
  for (auto tolerances : {o.surface, o.curve}) {
    if (!std::isfinite(tolerances.min_edge) ||
        !std::isfinite(tolerances.max_edge) ||
        !std::isfinite(tolerances.chord) ||
        !(0 < tolerances.min_edge &&
          tolerances.min_edge < tolerances.max_edge) ||
        tolerances.chord <= 0) {
      throw std::invalid_argument(
          "Invalid tolerances: require finite 0 < min-edge < max-edge and "
          "chord > 0");
    }
  }
  if (!std::isfinite(o.relative_uv_chord) || o.relative_uv_chord <= 0 ||
      o.steiner_limit < 0) {
    throw std::invalid_argument("Invalid UV tolerance or Steiner limit");
  }
  for (auto tolerances : {&o.surface, &o.curve}) {
    tolerances->max_edge = std::max(.01, tolerances->max_edge);
    tolerances->chord = std::max(1e-5, tolerances->chord);
  }
  return o;
}
}  // namespace iko::internal
