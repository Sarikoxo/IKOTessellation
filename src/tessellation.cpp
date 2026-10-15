#include "tessellation/tessellation.h"

#include <algorithm>
#include <exception>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tessellation/private/internal.h"

namespace iko {
namespace {
SurfaceInfo Info(const internal::Face& f) {
  return {f.source_de,
          f.id,
          f.surface.label,
          {f.surface.p, f.surface.q},
          f.surface.domain,
          (f.outer || !f.inner.empty()),
          f.surface.explicit_weights,
          f.reversed,
          f.solid_id,
          f.shell_id,
          f.support_type,
          f.source_entries};
}
}  // namespace
const char* StatusName(Status s) noexcept {
  switch (s) {
    case Status::kSuccess:
      return "success";
    case Status::kInputProblem:
      return "input_problem";
    case Status::kUnsupported:
      return "unsupported";
    default:
      return "computation_failed";
  }
}
Result Tessellate(const std::filesystem::path& file, const Options& o) {
  Result r;
  r.requested_options = o;
  r.whole_model = !o.surface_de;
  bool computing = false;
  std::string current_face;
  try {
    r.effective_options = internal::Effective(o);
    auto model = o.topology_consistent ? internal::ReadOcctRaw(file, true)
                                       : internal::ReadOcct(file);
    r.repair_trace = model.repair_trace;
    if (model.faces.empty()) {
      throw internal::Unsupported("OCCT model contains no effective faces");
    }
    r.model.imported_faces = model.faces.size();
    r.model.solids = model.solid_count;
    r.model.shells = model.shell_count;
    r.model.edges = static_cast<int>(model.edges.size());
    r.model.vertices = model.vertex_count;
    r.model.source_units = model.source_units;
    r.model.source_scale = model.source_scale;
    r.model.unit_to_mm = model.unit_to_mm;
    r.model.source_resolution = model.source_resolution;
    for (const auto& root : model.roots) {
      SourceRoot value{root.source_de, root.source_type, {}};
      for (int face_id : root.faces) {
        value.face_ids.push_back(model.faces.at(face_id - 1).id);
      }
      r.model.roots.push_back(std::move(value));
    }
    std::vector<const internal::Face*> selected;
    for (const auto& face : model.faces) {
      if (!o.surface_de || face.source_de == o.surface_de ||
          std::find(face.source_entries.begin(), face.source_entries.end(),
                    o.surface_de) != face.source_entries.end()) {
        selected.push_back(&face);
      }
    }
    if (selected.empty()) {
      throw std::invalid_argument(
          "Selected surface directory entry does not exist");
    }
    r.model.required_faces = selected.size();
    computing = true;
    // Publish atomically: any required-face failure leaves no partial mesh.
    Mesh complete;
    std::vector<FaceRange> ranges;
    std::vector<Mesh> topology_parts;
    if (o.topology_consistent) {
      topology_parts = internal::GenerateTopology(model, selected, r.effective_options);
    }
    std::size_t face_index = 0;
    for (const auto* face : selected) {
      current_face = face->id;
      auto part = o.topology_consistent
                      ? std::move(topology_parts.at(face_index++))
                      : internal::Generate(*face, r.effective_options);
      // LucidShape getMeshSurface exports the support's parametric orientation.
      // Keep the OCCT topological use direction in FaceRange, without editing
      // the exported triangle/normal records (verified on reversed BRep faces).
      FaceRange range{Info(*face),
                      complete.xyz.size(),
                      part.xyz.size(),
                      complete.triangles.size(),
                      part.triangles.size(),
                      complete.unpaired_uv_tail.size(),
                      part.unpaired_uv_tail.size(),
                      part.identity_mode,
                      part.branch,
                      face->boundary_edges};
      ranges.push_back(std::move(range));
      if (selected.size() == 1) {
        complete = std::move(part);
      } else {
        internal::Append(complete, std::move(part));
      }
    }
    if (selected.size() > 1) {
      complete.identity_mode = "face_ranges";
      complete.branch = "whole_model";
    }
    r.mesh = std::move(complete);
    r.faces = std::move(ranges);
    r.complete = true;
    r.status = Status::kSuccess;
    r.message = "Complete mesh generated for " +
                std::to_string(selected.size()) + " required faces";
  } catch (const internal::Unsupported& e) {
    r.status = Status::kUnsupported;
    r.message = e.what();
  } catch (const std::invalid_argument& e) {
    r.status = computing ? Status::kComputationFailed : Status::kInputProblem;
    r.message = e.what();
  } catch (const std::out_of_range& e) {
    r.status = computing ? Status::kComputationFailed : Status::kInputProblem;
    r.message = e.what();
  } catch (const std::exception& e) {
    r.status = Status::kComputationFailed;
    r.message = e.what();
  }
  if (!r.complete && !current_face.empty()) {
    r.message = "Incomplete model at " + current_face + ": " + r.message;
  }
  return r;
}
SurfaceList ListSurfaces(const std::filesystem::path& file) {
  SurfaceList result;
  try {
    auto model = internal::ReadOcct(file);
    for (const auto& face : model.faces) {
      result.surfaces.push_back(Info(face));
    }
    if (result.surfaces.empty()) {
      throw internal::Unsupported("OCCT model contains no effective faces");
    }
    result.status = Status::kSuccess;
    result.message = "Surfaces listed";
  } catch (const internal::Unsupported& e) {
    result.status = Status::kUnsupported;
    result.message = e.what();
  } catch (const std::invalid_argument& e) {
    result.status = Status::kInputProblem;
    result.message = e.what();
  } catch (const std::out_of_range& e) {
    result.status = Status::kInputProblem;
    result.message = e.what();
  } catch (const std::exception& e) {
    result.status = Status::kComputationFailed;
    result.message = e.what();
  }
  return result;
}
}  // namespace iko
