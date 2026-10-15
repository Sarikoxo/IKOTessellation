#include "tessellation/private/internal.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iterator>
#include <map>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Adaptor3d_CurveOnSurface.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include <BRepTools.hxx>
#include <BRepTools_NurbsConvertModification.hxx>
#include <BRepTools_ReShape.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BndLib_Add2dCurve.hxx>
#include <Bnd_Box2d.hxx>
#include <DEIGES_Parameters.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2dConvert.hxx>
#include <Geom2dConvert_ApproxCurve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_Circle.hxx>
#include <Geom2d_Ellipse.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAbs_Shape.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomConvert.hxx>
#include <GeomConvert_ApproxCurve.hxx>
#include <GeomLib_CheckCurveOnSurface.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Ellipse.hxx>
#include <Geom_Line.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <IGESControl_Controller.hxx>
#include <IGESControl_Reader.hxx>
#include <IGESData_IGESEntity.hxx>
#include <IGESData_IGESModel.hxx>
#include <IGESGeom_BSplineSurface.hxx>
#include <IGESGeom_TrimmedSurface.hxx>
#include <Interface_Static.hxx>
#include <Message.hxx>
#include <Message_Messenger.hxx>
#include <NCollection_Sequence.hxx>
#include <Precision.hxx>
#include <ProjLib_ComputeApproxOnPolarSurface.hxx>
#include <ShapeFix_Wireframe.hxx>
#include <ShapeProcess.hxx>
#include <Standard_Failure.hxx>
#include <Standard_Handle.hxx>
#include <Standard_Type.hxx>
#include <Standard_TypeDef.hxx>
#include <Standard_Version.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <TopAbs_Orientation.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_DataMapOfShapeInteger.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <TransferBRep.hxx>
#include <XSControl_TransferReader.hxx>
#include <XSControl_WorkSession.hxx>

#include "tessellation/tessellation.h"

static_assert(OCC_VERSION_MAJOR == 8 && OCC_VERSION_MINOR == 0 &&
                  OCC_VERSION_MAINTENANCE == 1,
              "The importer requires the pinned OCCT 8.0.1 release");
static_assert(OCCT_TESSELLATION_PROJECTION_POLICY == 1,
              "Build OCCT with dependencies/build_occt.ps1 for the required "
              "projection policy");

namespace iko::internal {
namespace {
std::mutex import_mutex;  // IGES uses process-wide Interface_Static state.

struct ProjectionPolicyGuard {
  bool original = Tessellation_ProjectionPolicy::Active;
  ProjectionPolicyGuard() {
    Tessellation_ProjectionPolicy::Active = true;
  }
  ~ProjectionPolicyGuard() {
    Tessellation_ProjectionPolicy::Active = original;
  }
  ProjectionPolicyGuard(const ProjectionPolicyGuard&) = delete;
  ProjectionPolicyGuard& operator=(const ProjectionPolicyGuard&) = delete;
};

struct SettingsGuard {
  std::vector<std::pair<std::string, std::string>> original;
  NCollection_Sequence<occ::handle<Message_Printer>> printers;
  SettingsGuard() : printers(Message::DefaultMessenger()->Printers()) {
    IGESControl_Controller::Init();
    // Pin the upstream defaults rather than inheriting a caller's IGES
    // settings.
    const std::pair<const char*, const char*> settings[] = {
        {"read.iges.bspline.continuity", "1"},
        {"read.iges.bspline.approxd1.mode", "Off"},
        {"read.iges.onlyvisible", "Off"},
        {"read.iges.faulty.entities", "Off"},
        {"read.precision.mode", "File"},
        {"read.precision.val", "0.001"},
        {"read.maxprecision.mode", "Preferred"},
        {"read.maxprecision.val", "1"},
        {"read.surfacecurve.mode", "Default"},
        {"read.encoderegularity.angle", "0.01"},
        {"xstep.cascade.unit", "MM"}};
    try {
      for (auto setting : settings) {
        if (!Interface_Static::IsPresent(setting.first)) {
          throw std::runtime_error(
              std::string("Missing OCCT import parameter: ") + setting.first);
        }
        original.emplace_back(setting.first,
                              Interface_Static::CVal(setting.first));
        if (!Interface_Static::SetCVal(setting.first, setting.second)) {
          throw std::runtime_error(
              std::string("Cannot set OCCT import parameter: ") +
              setting.first);
        }
      }
      Message::DefaultMessenger()->ChangePrinters().Clear();
    } catch (...) {
      Restore();
      throw;
    }
  }
  void Restore() noexcept {
    for (const auto& setting : original) {
      Interface_Static::SetCVal(setting.first.c_str(), setting.second.c_str());
    }
    Message::DefaultMessenger()->ChangePrinters() = printers;
  }
  ~SettingsGuard() {
    Restore();
  }
};

template <class S>
std::vector<double> Knots(const S& s, bool u) {
  std::vector<double> out;
  int n = u ? s->NbUKnots() : s->NbVKnots();
  for (int knot_index = 1; knot_index <= n; ++knot_index) {
    int count = u ? s->UMultiplicity(knot_index) : s->VMultiplicity(knot_index);
    double value = u ? s->UKnot(knot_index) : s->VKnot(knot_index);
    out.insert(out.end(), count, value);
  }
  return out;
}

Surface Support(const TopoDS_Face& face, std::string& type) {
  TopLoc_Location location;
  auto geometry = BRep_Tool::Surface(face, location);
  if (geometry.IsNull()) {
    throw Unsupported("OCCT face has no support surface");
  }
  type = geometry->DynamicType()->Name();
  auto spline = occ::handle<Geom_BSplineSurface>::DownCast(geometry);
  if (!spline.IsNull()) {
    spline = occ::handle<Geom_BSplineSurface>::DownCast(spline->Copy());
  } else {
    double u0;
    double u1;
    double v0;
    double v1;
    BRepTools::UVBounds(face, u0, u1, v0, v1);
    Rect natural;
    geometry->Bounds(natural[0], natural[1], natural[2], natural[3]);
    Rect bounded{u0, u1, v0, v1};
    // Close analytic endpoints to their natural parameter values within the
    // kernel's pcurve confusion. Preserve genuinely partial revolutions.
    for (int bound_index = 0; bound_index < 4; ++bound_index) {
      if (!Precision::IsInfinite(natural[bound_index]) &&
          std::abs(bounded[bound_index] - natural[bound_index]) <=
              Precision::PConfusion()) {
        bounded[bound_index] = natural[bound_index];
      }
    }
    u0 = bounded[0];
    u1 = bounded[1];
    v0 = bounded[2];
    v1 = bounded[3];
    if (!Finite(Rect{u0, u1, v0, v1}) || u0 >= u1 || v0 >= v1) {
      throw Unsupported("OCCT analytic support requires finite face bounds");
    }
    // Exact analytic-to-NURBS conversion, not a sampled approximation.
    spline = GeomConvert::SurfaceToBSplineSurface(
        new Geom_RectangularTrimmedSurface(geometry, u0, u1, v0, v1));
  }
  if (spline.IsNull()) {
    throw Unsupported("OCCT support conversion to NURBS failed");
  }
  if (spline->IsUPeriodic()) {
    spline->SetUNotPeriodic();
  }
  if (spline->IsVPeriodic()) {
    spline->SetVNotPeriodic();
  }
  Surface out;
  out.p = spline->UDegree();
  out.q = spline->VDegree();
  out.nu = spline->NbUPoles();
  out.nv = spline->NbVPoles();
  out.ku = Knots(spline, true);
  out.kv = Knots(spline, false);
  spline->Bounds(out.domain[0], out.domain[1], out.domain[2], out.domain[3]);
  out.closed_axes = {spline->IsUClosed(), spline->IsVClosed()};
  out.explicit_weights = spline->IsURational() || spline->IsVRational();
  for (int control_v_index = 1; control_v_index <= out.nv; ++control_v_index) {
    for (int control_u_index = 1; control_u_index <= out.nu;
         ++control_u_index) {
      auto p = spline->Pole(control_u_index, control_v_index)
                   .Transformed(location.Transformation());
      out.control.push_back({p.X(), p.Y(), p.Z()});
      // OCCT's polynomial surface evaluation ignores its dormant weight
      // array. Some 7.3 imported polynomial arrays contain zeroes; they are
      // not weights.
      double w = out.explicit_weights
                     ? spline->Weight(control_u_index, control_v_index)
                     : 1.;
      if (!(w > 0.)) {
        throw Unsupported(
            "OCCT support weight at (" + std::to_string(control_u_index) + "," +
            std::to_string(control_v_index) + ") is " + std::to_string(w));
      }
      out.weights.push_back(w);
    }
  }
  out.Validate();
  return out;
}

Curve Pcurve(const TopoDS_Edge& edge, const TopoDS_Face& face,
             bool endpoint_edit = false,
             const occ::handle<Geom2d_Curve>& replacement = {},
             bool retain_basis = false) {
  double first;
  double last;
  auto geometry = replacement;
  if (geometry.IsNull()) {
    geometry = BRep_Tool::CurveOnSurface(edge, face, first, last);
  } else {
    first = geometry->FirstParameter();
    last = geometry->LastParameter();
  }
  if (geometry.IsNull()) {
    throw Unsupported("OCCT boundary edge has no pcurve");
  }
  if (!std::isfinite(first) || !std::isfinite(last) || first >= last) {
    throw Unsupported("OCCT boundary edge has an invalid pcurve range");
  }
  auto spline = occ::handle<Geom2d_BSplineCurve>::DownCast(geometry);
  if (endpoint_edit) {
    // An edited analytic circle retains its angular parameterization.
    // The ordinary rational conversion uses a different parameter map.
    Geom2dConvert_ApproxCurve conversion(
        new Geom2d_TrimmedCurve(geometry, first, last), 1e-8, GeomAbs_C1, 100,
        8);
    if (!conversion.IsDone()) {
      throw Unsupported("Edited analytic pcurve approximation failed");
    }
    spline = conversion.Curve();
  } else if (!spline.IsNull()) {
    spline = occ::handle<Geom2d_BSplineCurve>::DownCast(spline->Copy());
    // The native B-spline conversion retains its complete basis after
    // wireframe range cuts. Newly projected/joint-converted curves instead
    // use their active range, just as analytic curves do.
    if (!retain_basis) {
      if (std::abs(first - spline->FirstParameter()) <=
          Precision::Confusion()) {
        first = spline->FirstParameter();
      }
      if (std::abs(last - spline->LastParameter()) <= Precision::Confusion()) {
        last = spline->LastParameter();
      }
      spline->Segment(first, last);
    }
  } else {
    spline = Geom2dConvert::CurveToBSplineCurve(
        new Geom2d_TrimmedCurve(geometry, first, last));
  }
  if (spline.IsNull()) {
    throw Unsupported("OCCT pcurve conversion failed");
  }
  if (spline->IsPeriodic()) {
    spline->SetNotPeriodic();
  }
  if (edge.Orientation() == TopAbs_REVERSED) {
    spline->Reverse();
  }
  Curve out;
  Bnd_Box2d bounds;
  BndLib_Add2dCurve::Add(geometry, first, last, 0., bounds);
  double bu0;
  double bv0;
  double bu1;
  double bv1;
  bounds.Get(bu0, bv0, bu1, bv1);
  out.import_bounds = Rect{bu0, bu1, bv0, bv1};
  auto basis = geometry;
  while (auto trimmed = occ::handle<Geom2d_TrimmedCurve>::DownCast(basis)) {
    basis = trimmed->BasisCurve();
  }
  out.analytic_line = !occ::handle<Geom2d_Line>::DownCast(basis).IsNull();
  out.degree = spline->Degree();
  out.linear_conversion = out.degree == 1 && spline->NbPoles() == 2;
  for (int knot_index = 1; knot_index <= spline->NbKnots(); ++knot_index) {
    out.knots.insert(out.knots.end(), spline->Multiplicity(knot_index),
                     spline->Knot(knot_index));
  }
  out.domain = {spline->FirstParameter(), spline->LastParameter()};
  for (int control_point_index = 1; control_point_index <= spline->NbPoles();
       ++control_point_index) {
    auto p = spline->Pole(control_point_index);
    out.control.push_back({p.X(), p.Y(), 0.});
    out.weights.push_back(
        spline->IsRational() ? spline->Weight(control_point_index) : 1.);
  }
  return out;
}

Rect CurveEnvelope(const Curve& c) {
  Rect box{INFINITY, -INFINITY, INFINITY, -INFINITY};
  for (auto control_point : c.control) {
    for (int axis = 0; axis < 2; ++axis) {
      box[2 * axis] = std::min(box[2 * axis], control_point[axis]);
      box[2 * axis + 1] = std::max(box[2 * axis + 1], control_point[axis]);
    }
  }
  return box;
}

Rect ConicEndpointEnvelope(const Curve& c,
                           const occ::handle<Geom2d_Curve>& geometry,
                           double first, double last) {
  auto circle = occ::handle<Geom2d_Circle>::DownCast(geometry);
  auto ellipse = occ::handle<Geom2d_Ellipse>::DownCast(geometry);
  if (circle.IsNull() && ellipse.IsNull()) {
    throw Unsupported("Edited conic envelope requires a circle or ellipse");
  }
  auto axes = circle.IsNull() ? ellipse->Position() : circle->Position();
  double major = circle.IsNull() ? ellipse->MajorRadius() : circle->Radius();
  double minor = circle.IsNull() ? ellipse->MinorRadius() : circle->Radius();
  Rect box{INFINITY, -INFINITY, INFINITY, -INFINITY};
  auto include = [&](Vec3 p) {
    for (int axis = 0; axis < 2; ++axis) {
      box[2 * axis] = std::min(box[2 * axis], p[axis]);
      box[2 * axis + 1] = std::max(box[2 * axis + 1], p[axis]);
    }
  };
  include(c.control.front());
  include(c.control.back());
  constexpr double pi = 3.1415926535897932384626433832795;
  for (int axis = 0; axis < 2; ++axis) {
    double root = std::atan2(minor * axes.YDirection().Coord(axis + 1),
                             major * axes.XDirection().Coord(axis + 1));
    auto lower = std::ceil((first - root) / pi);
    auto upper = std::floor((last - root) / pi);
    for (double half_turn_index = lower; half_turn_index <= upper;
         ++half_turn_index) {
      double parameter = root + half_turn_index * pi;
      if (parameter > first && parameter < last) {
        auto p = geometry->Value(parameter);
        include({p.X(), p.Y(), 0.});
      }
    }
  }
  return box;
}

struct ConversionBounds {
  Rect domain{INFINITY, -INFINITY, INFINITY, -INFINITY};
  std::map<int, Rect> edges;
};

ConversionBounds ComputeConversionBounds(
    const TopoDS_Face& original,
    const TopTools_DataMapOfShapeInteger& edge_ids) {
  ConversionBounds out;
  auto forward = original;
  forward.Orientation(TopAbs_FORWARD);
  for (TopoDS_Iterator wire_iterator(forward); wire_iterator.More();
       wire_iterator.Next()) {
    if (wire_iterator.Value().ShapeType() != TopAbs_WIRE) {
      continue;
    }
    std::vector<TopoDS_Edge> edges;
    std::vector<Curve> curves;
    std::vector<int> ids;
    for (BRepTools_WireExplorer edge_explorer(
             TopoDS::Wire(wire_iterator.Value()), forward);
         edge_explorer.More(); edge_explorer.Next()) {
      edges.push_back(edge_explorer.Current());
      curves.push_back(Pcurve(edge_explorer.Current(), forward));
      ids.push_back(static_cast<int>(edges.size()));
    }
    if (curves.empty()) {
      throw Unsupported("OCCT conversion bounds have an empty wire");
    }
    std::vector<bool> dirty(curves.size(), false);
    for (std::size_t curve_index = 0; curve_index < curves.size();
         ++curve_index) {
      std::size_t next_curve_index = (curve_index + 1) % curves.size();
      const auto& a = curves[curve_index];
      const auto& b = curves[next_curve_index];
      double gap = 0.;
      for (int coordinate_index = 0; coordinate_index < 3; ++coordinate_index) {
        gap = std::max(
            gap,
            std::abs(a.control.back()[coordinate_index] * a.weights.back() -
                     b.control.front()[coordinate_index] * b.weights.front()));
      }
      // Envelope invalidation tests each coordinate before the endpoint mean.
      // It is distinct from the representation-promotion and connector tests.
      if (gap > Precision::PConfusion()) {
        dirty[curve_index] = dirty[next_curve_index] = true;
      }
    }
    std::vector<std::string> trace;
    // This copy supplies bounds in the original parameterization. Actual output
    // pcurves are joined after joint NURBS conversion by the same numerical
    // rule.
    JoinOcctConversionEndpoints(curves, ids, trace, "");
    for (std::size_t curve_index = 0; curve_index < curves.size();
         ++curve_index) {
      const auto& c = curves[curve_index];
      auto hull = CurveEnvelope(c);
      auto box = c.import_bounds.value_or(hull);
      if (ids[curve_index]) {
        std::size_t source = ids[curve_index] - 1;
        double first;
        double last;
        auto geometry =
            BRep_Tool::CurveOnSurface(edges[source], forward, first, last);
        while (auto trimmed_curve =
                   occ::handle<Geom2d_TrimmedCurve>::DownCast(geometry)) {
          geometry = trimmed_curve->BasisCurve();
        }
        if (dirty[source] &&
            !occ::handle<Geom2d_BSplineCurve>::DownCast(geometry).IsNull()) {
          box = hull;
        }
        if (c.linear_conversion && c.degree != 1) {
          box = hull;
        }
        if (c.degree == 1) {
          double movement = 0.;
          for (int bound_index = 0; bound_index < 4; ++bound_index) {
            movement = std::max(movement,
                                std::abs(hull[bound_index] - box[bound_index]));
          }
          if (movement > Precision::PConfusion()) {
            box = hull;
          }
        }
        if (dirty[source] &&
            (!occ::handle<Geom2d_Circle>::DownCast(geometry).IsNull() ||
             !occ::handle<Geom2d_Ellipse>::DownCast(geometry).IsNull())) {
          box = ConicEndpointEnvelope(c, geometry, first, last);
        }
        int id = (edges[source].Orientation() == TopAbs_REVERSED ? -1 : 1) *
                 edge_ids.Find(edges[source]);
        out.edges[id] = box;
      }
      for (int bound_index = 0; bound_index < 4; ++bound_index) {
        out.domain[bound_index] =
            bound_index % 2
                ? std::max(out.domain[bound_index], box[bound_index])
                : std::min(out.domain[bound_index], box[bound_index]);
      }
    }
  }
  if (!Finite(out.domain) || out.domain[0] >= out.domain[1] ||
      out.domain[2] >= out.domain[3]) {
    throw Unsupported("OCCT conversion bounds are invalid");
  }
  return out;
}

class BoundsModification : public BRepTools_NurbsConvertModification {
 public:
  explicit BoundsModification(Rect domain) : bounds_(domain) {}
  // OCCT requires this name for its virtual conversion interface.
  Standard_Boolean NewSurface(const TopoDS_Face& face,
                              occ::handle<Geom_Surface>& surface,
                              TopLoc_Location& location,
                              Standard_Real& tolerance,
                              Standard_Boolean& reverse_wires,
                              Standard_Boolean& reverse_face) override {
    // An isolated bounds face passes the input-derived envelope to upstream
    // NewSurface. OCCT still owns natural/periodic limits, conversion and
    // knot reparameterization. No transferred topology or pcurve is mutated.
    BRep_Builder builder;
    TopoDS_Face domain_face;
    TopLoc_Location placement;
    auto support = BRep_Tool::Surface(face, placement);
    builder.MakeFace(domain_face, support, placement,
                     BRep_Tool::Tolerance(face));
    TopoDS_Wire wire;
    builder.MakeWire(wire);
    const auto& d = bounds_;
    Vec2 corners[] = {{d[0], d[2]}, {d[1], d[2]}, {d[1], d[3]}, {d[0], d[3]}};
    for (int corner_index = 0; corner_index < 4; ++corner_index) {
      TColgp_Array1OfPnt2d poles(1, 2);
      poles(1) = gp_Pnt2d(corners[corner_index][0], corners[corner_index][1]);
      poles(2) = gp_Pnt2d(corners[(corner_index + 1) % 4][0],
                          corners[(corner_index + 1) % 4][1]);
      TColStd_Array1OfReal knots(1, 2);
      knots(1) = 0.;
      knots(2) = 1.;
      TColStd_Array1OfInteger mults(1, 2);
      mults(1) = mults(2) = 2;
      occ::handle<Geom2d_BSplineCurve> curve =
          new Geom2d_BSplineCurve(poles, knots, mults, 1);
      TopoDS_Edge edge;
      builder.MakeEdge(edge);
      builder.UpdateEdge(edge, curve, support, placement, 0.);
      builder.Range(edge, support, placement, 0., 1.);
      builder.Add(wire, edge);
    }
    builder.Add(domain_face, wire);
    return BRepTools_NurbsConvertModification::NewSurface(
        domain_face, surface, location, tolerance, reverse_wires, reverse_face);
  }

 private:
  Rect bounds_;
};

class BoundsConversion : public BRepBuilderAPI_NurbsConvert {
 public:
  BoundsConversion(const TopoDS_Face& face, Rect domain) {
    myModification = new BoundsModification(domain);
    Perform(face);
  }
};

occ::handle<Geom2d_BSplineCurve> VertexPcurve(const TopoDS_Edge& edge,
                                              const TopoDS_Face& original,
                                              const TopoDS_Face& converted) {
  if (BRep_Tool::Degenerated(edge) ||
      BRepTools::IsReallyClosed(edge, original)) {
    return {};
  }
  double first;
  double last;
  auto geometry = BRep_Tool::Curve(edge, first, last);
  if (geometry.IsNull()) {
    throw Unsupported("OCCT boundary has no XYZ curve for endpoint conversion");
  }
  auto basis = geometry;
  while (auto trim = occ::handle<Geom_TrimmedCurve>::DownCast(basis)) {
    basis = trim->BasisCurve();
  }
  // Straight XYZ edges retain their joint conversion and the established UV
  // endpoint rule. Moving their XYZ poles changes matching line
  // representations.
  if (!occ::handle<Geom_Line>::DownCast(basis).IsNull()) {
    return {};
  }
  auto forward = edge;
  forward.Orientation(TopAbs_FORWARD);
  TopoDS_Vertex va;
  TopoDS_Vertex vb;
  TopExp::Vertices(forward, va, vb);
  if (va.IsNull() || vb.IsNull()) {
    throw Unsupported("OCCT boundary has no endpoint vertex");
  }
  auto a = BRep_Tool::Pnt(va);
  auto b = BRep_Tool::Pnt(vb);
  bool edit_a = a.Distance(geometry->Value(first)) > 1e-8;
  bool edit_b = b.Distance(geometry->Value(last)) > 1e-8;
  if (!edit_a && !edit_b) {
    return {};
  }
  ProjectionPolicyGuard projection_policy;
  auto spline = occ::handle<Geom_BSplineCurve>::DownCast(geometry);
  if (!spline.IsNull()) {
    spline = occ::handle<Geom_BSplineCurve>::DownCast(spline->Copy());
    if (std::abs(first - spline->FirstParameter()) <= Precision::Confusion()) {
      first = spline->FirstParameter();
    }
    if (std::abs(last - spline->LastParameter()) <= Precision::Confusion()) {
      last = spline->LastParameter();
    }
    spline->Segment(first, last);
  } else {
    if (!occ::handle<Geom_Circle>::DownCast(basis).IsNull() ||
        !occ::handle<Geom_Ellipse>::DownCast(basis).IsNull()) {
      GeomConvert_ApproxCurve approximation(
          new Geom_TrimmedCurve(geometry, first, last), 1e-8, GeomAbs_C1, 100,
          8);
      if (!approximation.IsDone()) {
        throw Unsupported("Edited analytic XYZ curve approximation failed");
      }
      spline = approximation.Curve();
    } else {
      spline = GeomConvert::CurveToBSplineCurve(
          new Geom_TrimmedCurve(geometry, first, last));
    }
  }
  if (spline.IsNull()) {
    throw Unsupported("Edited XYZ curve conversion failed");
  }
  if (spline->IsPeriodic()) {
    spline->SetNotPeriodic();
  }
  if (edit_a) {
    spline->SetPole(1, a);
  }
  if (edit_b) {
    spline->SetPole(spline->NbPoles(), b);
  }

  double f2;
  double l2;
  auto seed = BRep_Tool::CurveOnSurface(edge, original, f2, l2);
  if (seed.IsNull()) {
    throw Unsupported("Edited XYZ boundary has no original pcurve");
  }
  while (auto trim = occ::handle<Geom2d_TrimmedCurve>::DownCast(seed)) {
    seed = trim->BasisCurve();
  }
  if (!seed->IsPeriodic()) {
    if (seed->FirstParameter() - f2 > Precision::PConfusion()) {
      f2 = seed->FirstParameter();
    }
    if (l2 - seed->LastParameter() > Precision::PConfusion()) {
      l2 = seed->LastParameter();
    }
  }
  occ::handle<Geom2dAdaptor_Curve> curve2 =
      new Geom2dAdaptor_Curve(new Geom2d_TrimmedCurve(seed, f2, l2));
  occ::handle<GeomAdaptor_Curve> curve3 =
      new GeomAdaptor_Curve(spline, first, last);
  // Recompute the XYZ/UV discrepancy after the endpoint change. A face-wide
  // constant discards real curve error; the old edge tolerance includes an
  // endpoint discrepancy that may have disappeared. No ShapeFix runs here.
  double check_first = std::max(
      first, std::max(spline->FirstParameter(), seed->FirstParameter()));
  double check_last =
      std::min(last, std::min(spline->LastParameter(), seed->LastParameter()));
  occ::handle<GeomAdaptor_Curve> check_curve =
      new GeomAdaptor_Curve(spline, check_first, check_last);
  occ::handle<Geom2dAdaptor_Curve> check_uv =
      new Geom2dAdaptor_Curve(seed, check_first, check_last);
  occ::handle<GeomAdaptor_Surface> check_surface =
      new GeomAdaptor_Surface(BRep_Tool::Surface(original));
  GeomLib_CheckCurveOnSurface check(check_curve);
  // The former Perform(seed, true) argument disabled multithreading.
  check.SetParallel(false);
  check.Perform(new Adaptor3d_CurveOnSurface(check_uv, check_surface));
  if (!check.IsDone()) {
    throw Unsupported("Edited XYZ/UV discrepancy evaluation failed");
  }
  double tolerance = std::max(Precision::Confusion(), check.MaxDistance());
  auto surface = BRep_Tool::Surface(converted);
  Rect bounds;
  surface->Bounds(bounds[0], bounds[1], bounds[2], bounds[3]);
  Vec2 margin{(bounds[1] - bounds[0]) * .1, (bounds[3] - bounds[2]) * .1};
  for (int axis = 0; axis < 2; ++axis) {
    bool periodic = axis ? surface->IsVPeriodic() : surface->IsUPeriodic();
    if (periodic) {
      double period = axis ? surface->VPeriod() : surface->UPeriod();
      double width = bounds[2 * axis + 1] - bounds[2 * axis];
      if (period < width + 2 * margin[axis]) {
        margin[axis] = period <= width ? 0. : (period - width) * .5;
      }
    }
  }
  occ::handle<GeomAdaptor_Surface> support = new GeomAdaptor_Surface(
      surface, bounds[0] - margin[0], bounds[1] + margin[0],
      bounds[2] - margin[1], bounds[3] + margin[1]);
  ProjLib_ComputeApproxOnPolarSurface projection(curve2, curve3, support,
                                                 tolerance);
  if (!projection.IsDone() || projection.BSpline().IsNull()) {
    throw Unsupported("Edited XYZ boundary projection failed");
  }
  return projection.BSpline();
}

std::vector<Curve> Wire(const TopoDS_Wire& w, const TopoDS_Face& face,
                        const TopTools_DataMapOfShapeInteger& edges,
                        std::vector<int>& ids, bool retain_basis) {
  std::size_t count = 0;
  for (TopoDS_Iterator edge_iterator(w); edge_iterator.More();
       edge_iterator.Next()) {
    if (edge_iterator.Value().ShapeType() == TopAbs_EDGE) {
      ++count;
    }
  }
  std::vector<Curve> curves;
  for (BRepTools_WireExplorer edge_explorer(w, face); edge_explorer.More();
       edge_explorer.Next()) {
    curves.push_back(
        Pcurve(edge_explorer.Current(), face, false, {}, retain_basis));
    ids.push_back(
        (edge_explorer.Current().Orientation() == TopAbs_REVERSED ? -1 : 1) *
        edges.Find(edge_explorer.Current()));
  }
  if (curves.size() != count || curves.empty()) {
    throw Unsupported("OCCT wire traversal does not cover every boundary edge");
  }
  return curves;
}

int AncestorId(const TopoDS_Face& face,
               const TopTools_IndexedDataMapOfShapeListOfShape& relations,
               const TopTools_IndexedMapOfShape& ancestors) {
  if (!relations.Contains(face) || relations.FindFromKey(face).IsEmpty()) {
    return 0;
  }
  if (relations.FindFromKey(face).Extent() != 1) {
    throw Unsupported(
        "A face belongs to multiple ancestor shapes; occurrence mapping "
        "required");
  }
  return ancestors.FindIndex(relations.FindFromKey(face).First());
}
}  // namespace

Model ReadOcctRaw(const std::filesystem::path& path, bool topology) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::invalid_argument("Cannot read geometry file");
  }
  std::string bytes((std::istreambuf_iterator<char>(in)), {});
  if (bytes.empty()) {
    throw std::invalid_argument("Empty IGES file");
  }
  const auto sha = ContentSha256(bytes);
  std::lock_guard<std::mutex> lock(import_mutex);
  try {
    SettingsGuard settings;
    IGESControl_Reader reader;
    // 8.0.1 replaced the old resource/sequence defaults with typed parameters.
    // Run the upstream IGES ShapeFix operation once, without external
    // resources.
    reader.SetShapeFixParameters(
        DEIGES_Parameters::GetDefaultShapeFixParameters());
    ShapeProcess::OperationsFlags operations;
    operations.set(ShapeProcess::Operation::FixShape);
    reader.SetShapeProcessFlags(operations);
    reader.SetReadVisible(Standard_False);
    if (reader.ReadFile(path.u8string().c_str()) != IFSelect_RetDone) {
      throw std::invalid_argument("OCCT could not read IGES file");
    }
    auto source = reader.IGESModel();
    if (source.IsNull() || !source->NbEntities()) {
      throw std::invalid_argument("OCCT returned no IGES entities");
    }
    int roots = reader.NbRootsForTransfer();
    if (!roots) {
      throw Unsupported("IGES has no transferable roots");
    }
    int transferred = reader.TransferRoots();
    if (transferred != roots) {
      throw Unsupported("OCCT did not transfer every IGES root (" +
                        std::to_string(transferred) + "/" +
                        std::to_string(roots) + ")");
    }
    auto whole = reader.OneShape();
    if (whole.IsNull()) {
      throw Unsupported("OCCT returned an empty model");
    }
    // Native OCC import performs a wireframe stage after root transfer.
    // Direct NURBS supports need its active-range cuts and endpoint edits.
    // Analytic supports already use the joint conversion/endpoint adapter
    // below, which preserves their observed edit/conversion thresholds.
    // Retain the resulting hierarchy and apply face replacements to source
    // binders as well as roots, so coverage describes the effective model.
    occ::handle<BRepTools_ReShape> frame_changes = new BRepTools_ReShape;
    TopTools_IndexedMapOfShape initial_faces;
    TopExp::MapShapes(whole, TopAbs_FACE, initial_faces);
    int changed_frames = 0;
    for (int initial_face_index = 1;
         initial_face_index <= initial_faces.Extent(); ++initial_face_index) {
      const auto& original = initial_faces(initial_face_index);
      auto geometry = BRep_Tool::Surface(TopoDS::Face(original));
      while (
          auto trimmed =
              occ::handle<Geom_RectangularTrimmedSurface>::DownCast(geometry)) {
        geometry = trimmed->BasisSurface();
      }
      if (!geometry->IsKind(STANDARD_TYPE(Geom_BSplineSurface)) &&
          std::string(geometry->DynamicType()->Name()) !=
              "Geom_BezierSurface") {
        continue;
      }
      occ::handle<ShapeFix_Wireframe> fixer = new ShapeFix_Wireframe(original);
      fixer->SetPrecision(1e-8);
      fixer->SetMaxTolerance(1e-8);
      fixer->ModeDropSmallEdges() = true;
      bool small = fixer->FixSmallEdges();
      bool gaps = fixer->FixWireGaps();
      auto repaired = fixer->Shape();
      if (repaired.IsNull() || repaired.ShapeType() != TopAbs_FACE) {
        throw Unsupported("OCCT wireframe adaptation did not retain a face");
      }
      if (small || gaps) {
        ++changed_frames;
      }
      if (!repaired.IsSame(original)) {
        frame_changes->Replace(original, repaired);
      }
    }
    whole = frame_changes->Apply(whole);
    TopTools_IndexedMapOfShape face_map;
    TopTools_IndexedMapOfShape edges;
    TopTools_IndexedMapOfShape vertices;
    TopTools_IndexedMapOfShape shells;
    TopTools_IndexedMapOfShape solids;
    TopExp::MapShapes(whole, TopAbs_FACE, face_map);
    TopExp::MapShapes(whole, TopAbs_EDGE, edges);
    TopExp::MapShapes(whole, TopAbs_VERTEX, vertices);
    TopExp::MapShapes(whole, TopAbs_SHELL, shells);
    TopExp::MapShapes(whole, TopAbs_SOLID, solids);
    if (face_map.IsEmpty()) {
      throw Unsupported("OCCT model contains no faces");
    }
    int occurrences = 0;
    for (TopExp_Explorer face_explorer(whole, TopAbs_FACE);
         face_explorer.More(); face_explorer.Next()) {
      ++occurrences;
    }
    if (occurrences != face_map.Extent()) {
      throw Unsupported(
          "Repeated face occurrences require occurrence-aware assembly");
    }
    TopTools_DataMapOfShapeInteger edge_ids;
    for (int edge_index = 1; edge_index <= edges.Extent(); ++edge_index) {
      edge_ids.Bind(edges(edge_index), edge_index);
    }
    TopTools_IndexedDataMapOfShapeListOfShape face_shells;
    TopTools_IndexedDataMapOfShapeListOfShape face_solids;
    TopExp::MapShapesAndAncestors(whole, TopAbs_FACE, TopAbs_SHELL,
                                  face_shells);
    TopExp::MapShapesAndAncestors(whole, TopAbs_FACE, TopAbs_SOLID,
                                  face_solids);
    Model model;
    model.repair_trace.push_back(
        "OCCT direct-NURBS wireframe adaptation: FixSmallEdges then "
        "FixWireGaps; precision/max "
        "tolerance 1e-8 mm; " +
        std::to_string(changed_frames) +
        " changed faces; analytic endpoint handling retained in joint "
        "conversion");
    const auto& header = source->GlobalSection();
    model.source_units =
        header.UnitName().IsNull() ? "" : header.UnitName()->ToCString();
    model.source_scale = header.Scale();
    model.unit_to_mm = header.UnitValue();
    model.source_resolution = header.Resolution();
    model.shell_count = shells.Extent();
    model.solid_count = solids.Extent();
    model.vertex_count = vertices.Extent();
    for (int i = 1; i <= vertices.Extent(); ++i) {
      auto p = BRep_Tool::Pnt(TopoDS::Vertex(vertices(i)));
      model.vertex_positions.push_back({p.X(), p.Y(), p.Z()});
    }
    auto process = reader.WS()->TransferReader()->TransientProcess();
    std::map<int, std::vector<int>> face_sources;
    for (int entity_index = 1; entity_index <= source->NbEntities();
         ++entity_index) {
      auto ent = source->Entity(entity_index);
      if (ent->TypeNumber() != 144 && ent->TypeNumber() != 510) {
        continue;
      }
      auto shape = TransferBRep::ShapeResult(process, ent);
      if (shape.IsNull()) {
        throw Unsupported("Source face DE" + std::to_string(source->DNum(ent)) +
                          " has no effective OCCT shape");
      }
      shape = frame_changes->Apply(shape);
      for (TopExp_Explorer face_explorer(shape, TopAbs_FACE);
           face_explorer.More(); face_explorer.Next()) {
        int id = face_map.FindIndex(face_explorer.Current());
        if (!id) {
          // Resolve an unlocated transfer binder before the wireframe
          // stage replaces its complete-model TShape. Follow that exact
          // face replacement while retaining the model's placement.
          std::vector<int> partners;
          for (int initial_face_index = 1;
               initial_face_index <= initial_faces.Extent();
               ++initial_face_index) {
            if (face_explorer.Current().IsPartner(
                    initial_faces(initial_face_index))) {
              partners.push_back(initial_face_index);
            }
          }
          if (partners.size() == 1) {
            id = face_map.FindIndex(
                frame_changes->Apply(initial_faces(partners.front())));
            if (id) {
              model.repair_trace.push_back(
                  "Source face DE" + std::to_string(source->DNum(ent)) +
                  ": unique pre-wireframe TShape correspondence; "
                  "complete-model "
                  "placement retained");
            }
          }
        }
        if (!id) {
          // OCCT 7.3 MergeTransferInfo may replace a located binder face
          // with its unlocated ShapeFix result. The complete model still
          // owns the placed face. Match only a unique identical TShape,
          // including its support and full boundary topology; retain the
          // complete model's placement rather than the binder's location.
          std::vector<int> partners;
          for (int face_index = 1; face_index <= face_map.Extent();
               ++face_index) {
            if (face_explorer.Current().IsPartner(face_map(face_index))) {
              partners.push_back(face_index);
            }
          }
          if (partners.size() == 1) {
            id = partners.front();
            model.repair_trace.push_back("Source face DE" +
                                         std::to_string(source->DNum(ent)) +
                                         ": unique TShape correspondence; "
                                         "complete-model placement retained");
          }
        }
        if (!id) {
          TopLoc_Location source_location;
          auto source_surface = BRep_Tool::Surface(
              TopoDS::Face(face_explorer.Current()), source_location);
          while (auto trimmed =
                     occ::handle<Geom_RectangularTrimmedSurface>::DownCast(
                         source_surface)) {
            source_surface = trimmed->BasisSurface();
          }
          std::vector<int> matches;
          for (int face_index = 1; face_index <= face_map.Extent();
               ++face_index) {
            TopLoc_Location actual_location;
            auto actual_surface = BRep_Tool::Surface(
                TopoDS::Face(face_map(face_index)), actual_location);
            while (auto trimmed =
                       occ::handle<Geom_RectangularTrimmedSurface>::DownCast(
                           actual_surface)) {
              actual_surface = trimmed->BasisSurface();
            }
            if (source_surface == actual_surface &&
                source_location == actual_location) {
              matches.push_back(face_index);
            }
          }
          if (matches.size() == 1) {
            id = matches.front();
          } else {
            model.repair_trace.push_back("Unresolved source face DE" +
                                         std::to_string(source->DNum(ent)) +
                                         ": " + std::to_string(matches.size()) +
                                         " shared-support candidates");
          }
        }
        if (id) {
          face_sources[id].push_back(source->DNum(ent));
        }
      }
    }
    for (int root_index = 1; root_index <= roots; ++root_index) {
      auto ent = occ::handle<IGESData_IGESEntity>::DownCast(
          reader.RootForTransfer(root_index));
      auto shape =
          frame_changes->Apply(TransferBRep::ShapeResult(process, ent));
      ImportCoverage coverage{source->DNum(ent), ent->TypeNumber(), {}};
      for (TopExp_Explorer face_explorer(shape, TopAbs_FACE);
           face_explorer.More(); face_explorer.Next()) {
        int id = face_map.FindIndex(face_explorer.Current());
        if (!id) {
          throw Unsupported(
              "Transferred root has a face outside the complete model");
        }
        coverage.faces.push_back(id);
        if (face_sources[id].empty()) {
          face_sources[id].push_back(source->DNum(ent));
        }
      }
      model.roots.push_back(std::move(coverage));
    }
    std::map<int, int> source_counts;
    for (const auto& face_source_entry : face_sources) {
      if (!face_source_entry.second.empty()) {
        ++source_counts[*std::min_element(face_source_entry.second.begin(),
                                          face_source_entry.second.end())];
      }
    }
    std::map<int, int> source_parts;
    for (int edge_index = 1; edge_index <= edges.Extent(); ++edge_index) {
      auto e = TopoDS::Edge(edges(edge_index));
      e.Orientation(TopAbs_FORWARD);
      ImportedEdge info;
      info.id = edge_index;
      info.degenerate = BRep_Tool::Degenerated(e);
      info.same_parameter = BRep_Tool::SameParameter(e);
      info.tolerance = BRep_Tool::Tolerance(e);
      TopoDS_Vertex a;
      TopoDS_Vertex b;
      TopExp::Vertices(e, a, b);
      info.vertices = {vertices.FindIndex(a), vertices.FindIndex(b)};
      double first;
      double last;
      TopLoc_Location loc;
      auto c = BRep_Tool::Curve(e, loc, first, last);
      if (!c.IsNull()) {
        info.curve_type = c->DynamicType()->Name();
        info.domain = {first, last};
        auto basis = c;
        while (auto trimmed = occ::handle<Geom_TrimmedCurve>::DownCast(basis)) {
          basis = trimmed->BasisCurve();
        }
        if (auto spline = occ::handle<Geom_BSplineCurve>::DownCast(basis)) {
          for (int knot = 1; knot <= spline->NbKnots(); ++knot) {
            double t = spline->Knot(knot);
            if (t > first && t < last) info.knots.push_back(t);
          }
        }
        info.point = [c, loc](double t) {
          auto p = c->Value(t).Transformed(loc.Transformation());
          return Vec3{p.X(), p.Y(), p.Z()};
        };
      }
      model.edges.push_back(std::move(info));
    }
    for (int face_index = 1; face_index <= face_map.Extent(); ++face_index) {
      auto original = TopoDS::Face(face_map(face_index));
      auto f = original;
      std::string original_type =
          BRep_Tool::Surface(original)->DynamicType()->Name();
      TopTools_DataMapOfShapeInteger converted_edge_ids;
      std::optional<ConversionBounds> bounds;
      bool convert = original_type != "Geom_Plane" &&
                     original_type != "Geom_BSplineSurface" &&
                     original_type != "Geom_BezierSurface";
      if (original_type == "Geom_Plane") {
        for (TopExp_Explorer edge_explorer(original, TopAbs_EDGE);
             edge_explorer.More(); edge_explorer.Next()) {
          double first;
          double last;
          auto c = BRep_Tool::Curve(TopoDS::Edge(edge_explorer.Current()),
                                    first, last);
          if (c.IsNull()) {
            continue;
          }
          while (auto trimmed = occ::handle<Geom_TrimmedCurve>::DownCast(c)) {
            c = trimmed->BasisCurve();
          }
          std::string type = c->DynamicType()->Name();
          // Existing NURBS XYZ/pcurves already use a compatible parameter
          // representation; analytic conics need the joint conversion.
          convert |= type != "Geom_Line" && type != "Geom_BSplineCurve" &&
                     type != "Geom_BezierCurve";
        }
      }
      if (convert) {
        // Convert analytic curved support and its boundary parameterization
        // together. Retain supplied plane/NURBS UV geometry untouched.
        bounds = ComputeConversionBounds(original, edge_ids);
        ProjectionPolicyGuard projection_policy;
        BoundsConversion conversion(original, bounds->domain);
        if (!conversion.IsDone()) {
          throw Unsupported("OCCT analytic NURBS conversion failed");
        }
        f = TopoDS::Face(conversion.ModifiedShape(original));
        TopTools_IndexedMapOfShape local_edges;
        TopExp::MapShapes(original, TopAbs_EDGE, local_edges);
        for (int local_edge_index = 1; local_edge_index <= local_edges.Extent();
             ++local_edge_index) {
          converted_edge_ids.Bind(
              conversion.ModifiedShape(local_edges(local_edge_index)),
              edges.FindIndex(local_edges(local_edge_index)));
        }
      }
      f.Orientation(TopAbs_FORWARD);
      Face face;
      face.source_entries = face_sources[face_index];
      std::sort(face.source_entries.begin(), face.source_entries.end());
      face.source_entries.erase(
          std::unique(face.source_entries.begin(), face.source_entries.end()),
          face.source_entries.end());
      if (!face.source_entries.empty()) {
        face.source_de = face.source_entries.front();
      }
      face.id = sha + ":DE" + std::to_string(face.source_de);
      if (!face.source_de) {
        face.id += "/face" + std::to_string(face_index);
      } else if (source_counts[face.source_de] > 1) {
        face.id += "/part" + std::to_string(++source_parts[face.source_de]);
      }
      face.reversed = original.Orientation() == TopAbs_REVERSED;
      if (original.Orientation() != TopAbs_FORWARD && !face.reversed) {
        throw Unsupported(
            "Internal/external OCCT face orientation is not adapted");
      }
      face.shell_id = AncestorId(original, face_shells, shells);
      face.solid_id = AncestorId(original, face_solids, solids);
      try {
        face.surface = Support(f, face.support_type);
      } catch (const std::exception& e) {
        throw Unsupported("Source DE" + std::to_string(face.source_de) +
                          " support: " + e.what());
      }
      face.support_type = original_type;
      face.geometrically_closed_axes = face.surface.closed_axes;
      // Geometric closure is not the NSurface wrapping flag. OCC-imported
      // native NSurfaces export both optional flags as zero; preserve seam
      // topology separately and use their ordinary NURBS evaluation.
      face.surface.closed_axes = {false, false};
      if (face.source_de) {
        auto ent = source->Entity((face.source_de + 1) / 2);
        if (ent->HasName()) {
          face.surface.label = ent->NameValue()->ToCString();
        }
        auto trimmed = occ::handle<IGESGeom_TrimmedSurface>::DownCast(ent);
        if (!trimmed.IsNull()) {
          face.source_natural_boundary = trimmed->OuterBoundaryType() == 0 &&
                                         trimmed->NbInnerContours() == 0;
          ent = trimmed->Surface();
        }
        auto bspline = occ::handle<IGESGeom_BSplineSurface>::DownCast(ent);
        if (!bspline.IsNull()) {
          face.source_closed_axes = {bspline->IsClosedU(),
                                     bspline->IsClosedV()};
        }
      }
      auto outer = BRepTools::OuterWire(f);
      if (outer.IsNull()) {
        throw Unsupported("OCCT face has no outer boundary");
      }
      face.boundary_edges.emplace_back();
      const auto& actual_edge_ids = convert ? converted_edge_ids : edge_ids;
      face.outer =
          Wire(outer, f, actual_edge_ids, face.boundary_edges.back(), !convert && !topology);
      for (TopoDS_Iterator wire_iterator(f); wire_iterator.More();
           wire_iterator.Next()) {
        if (wire_iterator.Value().ShapeType() != TopAbs_WIRE ||
            wire_iterator.Value().IsSame(outer)) {
          continue;
        }
        face.boundary_edges.emplace_back();
        face.inner.push_back(Wire(TopoDS::Wire(wire_iterator.Value()), f,
                                  actual_edge_ids, face.boundary_edges.back(),
                                  !convert && !topology));
      }
      if (convert && original_type == "Geom_Plane") {
        // Retain each conic's analytic envelope and input-derived angular
        // conversion for an endpoint edit; keep unedited joint conversion.
        auto retain_plane_pcurves = [&](auto& curves, const auto& ids) {
          for (std::size_t curve_index = 0; curve_index < curves.size();
               ++curve_index) {
            auto original_edge =
                TopoDS::Edge(edges(std::abs(ids[curve_index])));
            original_edge.Orientation(ids[curve_index] < 0 ? TopAbs_REVERSED
                                                           : TopAbs_FORWARD);
            auto c = Pcurve(original_edge, original);
            curves[curve_index].import_bounds =
                bounds->edges.at(ids[curve_index]);
            curves[curve_index].linear_conversion = c.linear_conversion;
            curves[curve_index].analytic_line = c.analytic_line;
            double first;
            double last;
            auto original_pc =
                BRep_Tool::CurveOnSurface(original_edge, original, first, last);
            while (
                auto trimmed_curve =
                    occ::handle<Geom2d_TrimmedCurve>::DownCast(original_pc)) {
              original_pc = trimmed_curve->BasisCurve();
            }
            if (!occ::handle<Geom2d_Circle>::DownCast(original_pc).IsNull()) {
              curves[curve_index].endpoint_edit_conversion = [original_edge,
                                                              original] {
                try {
                  return Pcurve(original_edge, original, true);
                } catch (const Standard_Failure& e) {
                  throw Unsupported(std::string("Edited circle conversion: ") +
                                    e.GetMessageString());
                }
              };
            }
          }
        };
        retain_plane_pcurves(*face.outer, face.boundary_edges[0]);
        for (std::size_t inner_loop_index = 0;
             inner_loop_index < face.inner.size(); ++inner_loop_index) {
          retain_plane_pcurves(face.inner[inner_loop_index],
                               face.boundary_edges[inner_loop_index + 1]);
        }
      }
      if (convert && original_type != "Geom_Plane") {
        auto align_vertex_pcurves = [&](auto& curves, const auto& ids) {
          for (std::size_t curve_index = 0; curve_index < curves.size();
               ++curve_index) {
            auto original_edge =
                TopoDS::Edge(edges(std::abs(ids[curve_index])));
            original_edge.Orientation(ids[curve_index] < 0 ? TopAbs_REVERSED
                                                           : TopAbs_FORWARD);
            auto replacement = VertexPcurve(original_edge, original, f);
            if (!replacement.IsNull()) {
              curves[curve_index] =
                  Pcurve(original_edge, f, false, replacement);
              model.repair_trace.push_back(
                  face.id +
                  ": topological XYZ endpoint conversion and projection");
            }
          }
        };
        align_vertex_pcurves(*face.outer, face.boundary_edges[0]);
        for (std::size_t inner_loop_index = 0;
             inner_loop_index < face.inner.size(); ++inner_loop_index) {
          align_vertex_pcurves(face.inner[inner_loop_index],
                               face.boundary_edges[inner_loop_index + 1]);
        }
      }
      for (const auto& boundary_loop : face.boundary_edges) {
        for (int edge_id : boundary_loop) {
          if (BRep_Tool::IsClosed(TopoDS::Edge(edges(std::abs(edge_id))),
                                  original)) {
            model.edges[std::abs(edge_id) - 1].seam = true;
          }
        }
      }
      model.faces.push_back(std::move(face));
    }
    model.repair_trace.push_back(
        "OCCT 8.0.1 IGES TransferRoots: fixed settings, default "
        "ShapeFix once; no legacy file repair");
    model.repair_trace.push_back(
        "Analytic joint/edited-endpoint projection: scoped 24-point "
        "quadrature, ascending degrees and original normal-root "
        "solve; OCCT release defaults restored outside this stage");
    return model;
  } catch (const Standard_Failure& e) {
    throw Unsupported(std::string("OCCT import/adaptation: ") +
                      e.GetMessageString());
  }
}
Model ReadOcct(const std::filesystem::path& path) {
  auto model = ReadOcctRaw(path);
  for (auto& face : model.faces) {
    AdaptOcctFace(face, model.repair_trace);
  }
  return model;
}
}  // namespace iko::internal
