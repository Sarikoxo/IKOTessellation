// Validation executable only; never linked into the installed CLI/library.
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>
#include "tessellation/private/exact.h"
#include "tessellation/private/internal.h"
#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif
using namespace iko;
using namespace iko::internal;
static void Quote(std::ostream& out, const std::string& x) {
  out << '"';
  for (unsigned char c : x) {
    if (c == '"' || c == '\\')
      out << '\\';
    if (c == '\n')
      out << "\\n";
    else if (c == '\r')
      out << "\\r";
    else if (c == '\t')
      out << "\\t";
    else
      out << c;
  }
  out << '"';
}
static void Value(std::ostream& out, double x) {
  if (std::isnan(x))
    out << "NaN";
  else if (std::isinf(x))
    out << (x > 0 ? "Infinity" : "-Infinity");
  else
    out << std::setprecision(17) << x;
}
template <std::size_t N>
static void Array(std::ostream& out, const std::array<double, N>& v) {
  out << '[';
  for (std::size_t i = 0; i < N; ++i) {
    if (i)
      out << ',';
    Value(out, v[i]);
  }
  out << ']';
}
template <std::size_t N>
static void Arrays(std::ostream& out,
                   const std::vector<std::array<double, N>>& vs) {
  out << '[';
  for (std::size_t i = 0; i < vs.size(); ++i) {
    if (i)
      out << ',';
    Array(out, vs[i]);
  }
  out << ']';
}
namespace iko::internal {
void DumpMesh(const Mesh& m, std::ostream& out) {
  out << "{\"xyz\":";
  Arrays(out, m.xyz);
  out << ",\"uv\":";
  Arrays(out, m.uv);
  out << ",\"normals\":";
  Arrays(out, m.normals);
  out << ",\"triangles\":[";
  for (std::size_t i = 0; i < m.triangles.size(); ++i) {
    if (i)
      out << ',';
    auto t = m.triangles[i];
    out << '[' << t[0] << ',' << t[1] << ',' << t[2] << ']';
  }
  out << "],\"identities\":[";
  for (std::size_t i = 0; i < m.identities.size(); ++i) {
    if (i)
      out << ',';
    Quote(out, m.identities[i]);
  }
  out << "],\"metadata\":{\"identity_mode\":";
  Quote(out, m.identity_mode);
  out << ",\"unpaired_uv_tail\":";
  Arrays(out, m.unpaired_uv_tail);
  out << "}}";
}
}  // namespace iko::internal
static Curve ReadCurve(std::istream& in) {
  Curve c;
  int n;
  in >> c.degree >> n >> c.domain[0] >> c.domain[1];
  if (n <= 0 || n > 1000000)
    throw std::invalid_argument("Invalid test curve dimensions");
  c.knots.resize(n + c.degree + 1);
  for (auto& x : c.knots)
    in >> x;
  c.weights.resize(n);
  for (auto& x : c.weights)
    in >> x;
  c.control.resize(n);
  for (auto& p : c.control)
    for (auto& x : p)
      in >> x;
  return c;
}
static void DumpGeometry(const Face& f) {
  const auto& s = f.surface;
  std::cout << "{\"domain\":";
  Array(std::cout, s.domain);
  std::cout << ",\"control\":";
  Arrays(std::cout, s.control);
  std::cout << ",\"weights\":[";
  for (std::size_t i = 0; i < s.weights.size(); ++i) {
    if (i)
      std::cout << ',';
    Value(std::cout, s.weights[i]);
  }
  std::cout << "],\"knots_u\":[";
  for (std::size_t i = 0; i < s.ku.size(); ++i) {
    if (i)
      std::cout << ',';
    Value(std::cout, s.ku[i]);
  }
  std::cout << "],\"knots_v\":[";
  for (std::size_t i = 0; i < s.kv.size(); ++i) {
    if (i)
      std::cout << ',';
    Value(std::cout, s.kv[i]);
  }
  std::cout << "],\"face_id\":";
  Quote(std::cout, f.id);
  std::cout << ",\"explicit_weights\":"
            << (s.explicit_weights ? "true" : "false")
            << ",\"constructed_rectangle\":"
            << (f.constructed_rectangle ? "true" : "false") << ",\"outer\":[";
  if (f.outer)
    for (std::size_t ci = 0; ci < f.outer->size(); ++ci) {
      if (ci)
        std::cout << ',';
      const auto& c = (*f.outer)[ci];
      std::cout << "{\"degree\":" << c.degree << ",\"domain\":";
      Array(std::cout, c.domain);
      std::cout << ",\"control\":";
      Arrays(std::cout, c.control);
      std::cout << ",\"knots\":[";
      for (std::size_t i = 0; i < c.knots.size(); ++i) {
        if (i)
          std::cout << ',';
        Value(std::cout, c.knots[i]);
      }
      std::cout << "],\"weights\":[";
      for (std::size_t i = 0; i < c.weights.size(); ++i) {
        if (i)
          std::cout << ',';
        Value(std::cout, c.weights[i]);
      }
      std::cout << "]}";
    }
  std::cout << "],\"has_outer\":" << (f.outer ? "true" : "false") << "}\n";
}
static Face ReadFace(std::istream& in, Options& o, std::vector<Vec2>& queries) {
  Face f;
  in >> std::quoted(f.id);
  auto& s = f.surface;
  in >> s.p >> s.q >> s.nu >> s.nv >> s.explicit_weights >> s.closed_axes[0] >>
      s.closed_axes[1];
  for (auto& x : s.domain)
    in >> x;
  if (s.nu <= 0 || s.nv <= 0 || s.nu > 1000000 / s.nv)
    throw std::invalid_argument("Invalid test dimensions");
  s.ku.resize(s.nu + s.p + 1);
  s.kv.resize(s.nv + s.q + 1);
  s.weights.resize(s.nu * s.nv);
  s.control.resize(s.nu * s.nv);
  for (auto& x : s.ku)
    in >> x;
  for (auto& x : s.kv)
    in >> x;
  for (auto& x : s.weights)
    in >> x;
  for (auto& p : s.control)
    for (auto& x : p)
      in >> x;
  int loops;
  in >> loops;
  for (int i = 0; i < loops; ++i) {
    int n;
    in >> n;
    std::vector<Curve> cs;
    for (int j = 0; j < n; ++j)
      cs.push_back(ReadCurve(in));
    if (i == 0)
      f.outer = std::move(cs);
    else
      f.inner.push_back(std::move(cs));
  }
  in >> o.surface.min_edge >> o.surface.max_edge >> o.surface.chord >>
      o.curve.min_edge >> o.curve.max_edge >> o.curve.chord >>
      o.relative_uv_chord >> o.smooth_poles >> o.steiner_limit;
  int count;
  in >> count;
  queries.resize(count);
  for (auto& p : queries)
    for (auto& x : p)
      in >> x;
  if (!in)
    throw std::invalid_argument("Truncated test fixture");
  s.Validate();
  return f;
}
template <class V>
void ImportArray(std::ostream& out, const V& values) {
  out << '[';
  bool first = true;
  for (const auto& x : values) {
    if (!first)
      out << ',';
    first = false;
    out << x;
  }
  out << ']';
}
template <class V>
void ImportArrays(std::ostream& out, const V& values) {
  out << '[';
  bool first = true;
  for (const auto& x : values) {
    if (!first)
      out << ',';
    first = false;
    ImportArray(out, x);
  }
  out << ']';
}
void DumpCurve(std::ostream& out, const Curve& c) {
  out << "{\"degree\":" << c.degree << ",\"analytic_line\":" << c.analytic_line
      << ",\"knots\":";
  ImportArray(out, c.knots);
  out << ",\"domain\":";
  ImportArray(out, c.domain);
  out << ",\"control\":";
  ImportArrays(out, c.control);
  out << ",\"weights\":";
  ImportArray(out, c.weights);
  out << '}';
}
void DumpLoop(std::ostream& out, const std::vector<Curve>& cs) {
  out << '[';
  for (std::size_t i = 0; i < cs.size(); ++i) {
    if (i)
      out << ',';
    DumpCurve(out, cs[i]);
  }
  out << ']';
}

static void DumpImportModel(const Model& m, std::ostream& out) {
  out << std::setprecision(17) << "{\"units\":" << std::quoted(m.source_units)
      << ",\"scale\":" << m.source_scale << ",\"unit_to_mm\":" << m.unit_to_mm
      << ",\"resolution\":" << m.source_resolution
      << ",\"shells\":" << m.shell_count << ",\"solids\":" << m.solid_count
      << ",\"vertices\":" << m.vertex_count << ",\"roots\":[";
  for (std::size_t i = 0; i < m.roots.size(); ++i) {
    if (i)
      out << ',';
    const auto& r = m.roots[i];
    out << "{\"de\":" << r.source_de << ",\"type\":" << r.source_type
        << ",\"faces\":";
    ImportArray(out, r.faces);
    out << '}';
  }
  out << "],\"edges\":[";
  for (std::size_t i = 0; i < m.edges.size(); ++i) {
    if (i)
      out << ',';
    const auto& e = m.edges[i];
    out << "{\"id\":" << e.id << ",\"vertices\":";
    ImportArray(out, e.vertices);
    out << ",\"degenerate\":" << e.degenerate
        << ",\"same_parameter\":" << e.same_parameter << ",\"seam\":" << e.seam
        << ",\"tolerance\":" << e.tolerance
        << ",\"type\":" << std::quoted(e.curve_type) << '}';
  }
  out << "],\"faces\":[";
  for (std::size_t i = 0; i < m.faces.size(); ++i) {
    if (i)
      out << ',';
    const auto& f = m.faces[i];
    const auto& s = f.surface;
    out << "{\"index\":" << i + 1 << ",\"id\":" << std::quoted(f.id)
        << ",\"de\":" << f.source_de << ",\"sources\":";
    ImportArray(out, f.source_entries);
    out << ",\"reversed\":" << f.reversed << ",\"shell\":" << f.shell_id
        << ",\"solid\":" << f.solid_id
        << ",\"support_type\":" << std::quoted(f.support_type)
        << ",\"boundary_edges\":";
    ImportArrays(out, f.boundary_edges);
    out << ",\"label\":" << std::quoted(s.label) << ",\"degree\":[" << s.p
        << ',' << s.q << "],\"shape\":[" << s.nv << ',' << s.nu
        << "],\"domain\":";
    ImportArray(out, s.domain);
    out << ",\"knots_u\":";
    ImportArray(out, s.ku);
    out << ",\"knots_v\":";
    ImportArray(out, s.kv);
    out << ",\"control\":";
    ImportArrays(out, s.control);
    out << ",\"weights\":";
    ImportArray(out, s.weights);
    out << ",\"explicit_weights\":" << s.explicit_weights;
    out << ",\"closed_axes\":";
    ImportArray(out, s.closed_axes);
    out << ",\"geometrically_closed_axes\":";
    ImportArray(out, f.geometrically_closed_axes);
    out << ",\"source_natural_boundary\":" << f.source_natural_boundary;
    out << ",\"source_closed_axes\":";
    ImportArray(out, f.source_closed_axes);
    out << ",\"constructed_rectangle\":" << f.constructed_rectangle;
    out << ",\"has_outer\":" << bool(f.outer) << ",\"outer\":";
    DumpLoop(out, f.outer ? *f.outer : std::vector<Curve>{});
    out << ",\"inner\":[";
    for (std::size_t j = 0; j < f.inner.size(); ++j) {
      if (j)
        out << ',';
      DumpLoop(out, f.inner[j]);
    }
    out << "]}";
  }
  out << "],\"trace\":[";
  for (std::size_t i = 0; i < m.repair_trace.size(); ++i) {
    if (i)
      out << ',';
    out << std::quoted(m.repair_trace[i]);
  }
  out << "]}\n";
}
static int Run(int argc, char** argv) {
  if (argc == 2 && std::string(argv[1]) == "--describe") {
    std::cout
        << "{\"protocol_revision\":1,\"sdk_project\":\"iko\","
           "\"modes\":[\"mesh\",\"evaluate\",\"api_model\",\"model_"
           "geometry\",\"exact\"]}\n";
    return 0;
  }
  if (argc < 3) {
    std::cerr
        << "usage: tessellation_probe fixture (mesh|evaluate|api_model|exact) "
           "[IGS DE]; or IGS model_geometry [raw]\n";
    return 2;
  }
  try {
    std::cout << std::setprecision(17);
    if (std::string(argv[2]) == "model_geometry") {
      if (argc > 4 || (argc == 4 && std::string(argv[3]) != "raw"))
        throw std::invalid_argument(
            "model_geometry accepts only an optional raw flag");
      const auto path = std::filesystem::u8path(argv[1]);
      DumpImportModel(argc == 4 ? ReadOcctRaw(path) : ReadOcct(path),
                      std::cout);
      return 0;
    }
    Options o;
    std::string mode = argv[2];
    Face f;
    std::vector<Vec2> queries;
    if (mode == "exact") {
      std::ifstream input(std::filesystem::u8path(argv[1]));
      if (!input)
        throw std::invalid_argument("Cannot read exact arithmetic input");
      double a, b;
      std::cout << '[';
      bool first = true;
      while (input >> a >> b) {
        if (!first)
          std::cout << ',';
        first = false;
        exact::Rat x(a), y(b);
        std::cout << '[';
        Value(std::cout, (x + y).Value());
        std::cout << ',';
        Value(std::cout, (x - y).Value());
        std::cout << ',';
        Value(std::cout, (x * y).Value());
        std::cout << ',';
        Value(std::cout, (x / y).Value());
        std::cout << ']';
      }
      std::cout << "]\n";
      return 0;
    }
    if (mode == "iges") {
      auto model = ReadOcct(std::filesystem::u8path(argv[1]));
      int de = argc > 3 ? std::stoi(argv[3]) : 0;
      auto it = std::find_if(
          model.faces.begin(), model.faces.end(),
          [&](const Face& f) { return de == 0 || f.source_de == de; });
      if (it == model.faces.end())
        throw std::invalid_argument("Missing test face");
      f = *it;
    } else {
      std::ifstream in(std::filesystem::u8path(argv[1]));
      f = ReadFace(in, o, queries);
    }
    if (mode == "import_mesh" || mode == "import_geometry") {
      if (argc < 5)
        throw std::invalid_argument("Missing IGES test source or DE");
      auto model = ReadOcct(std::filesystem::u8path(argv[3]));
      int de = std::stoi(argv[4]);
      auto it = std::find_if(
          model.faces.begin(), model.faces.end(),
          [&](const Face& f) { return de == 0 || f.source_de == de; });
      if (it == model.faces.end())
        throw std::invalid_argument("Missing test face");
      f = *it;
      if (mode == "import_geometry") {
        DumpGeometry(f);
        return 0;
      }
    }
    if (mode == "api_mesh" || mode == "api_model") {
      if (argc < 5)
        throw std::invalid_argument("Missing IGES test source or DE");
      o.surface_de = std::stoi(argv[4]);
      auto result = iko::Tessellate(std::filesystem::u8path(argv[3]), o);
      if (!result.Succeeded() || !result.complete) {
        std::cerr << result.message << '\n';
        return result.Succeeded() ? 4 : int(result.status);
      }
      if (mode == "api_model") {
        std::cout << "{\"complete\":" << (result.complete ? "true" : "false")
                  << ",\"whole_model\":"
                  << (result.whole_model ? "true" : "false")
                  << ",\"imported_faces\":" << result.model.imported_faces
                  << ",\"required_faces\":" << result.model.required_faces
                  << ",\"solids\":" << result.model.solids
                  << ",\"shells\":" << result.model.shells
                  << ",\"edges\":" << result.model.edges
                  << ",\"vertices\":" << result.model.vertices
                  << ",\"source_scale\":" << result.model.source_scale
                  << ",\"unit_to_mm\":" << result.model.unit_to_mm
                  << ",\"source_resolution\":" << result.model.source_resolution
                  << ",\"source_units\":";
        Quote(std::cout, result.model.source_units);
        std::cout << ",\"roots\":[";
        for (std::size_t i = 0; i < result.model.roots.size(); ++i) {
          if (i)
            std::cout << ',';
          const auto& r = result.model.roots[i];
          std::cout << "{\"de\":" << r.directory_entry
                    << ",\"type\":" << r.entity_type << ",\"faces\":[";
          for (std::size_t j = 0; j < r.face_ids.size(); ++j) {
            if (j)
              std::cout << ',';
            Quote(std::cout, r.face_ids[j]);
          }
          std::cout << "]}";
        }
        std::cout << ']' << ",\"faces\":[";
        for (std::size_t i = 0; i < result.faces.size(); ++i) {
          if (i)
            std::cout << ',';
          const auto& f = result.faces[i];
          std::cout << "{\"de\":" << f.surface.directory_entry << ",\"id\":";
          Quote(std::cout, f.surface.instance_id);
          std::cout << ",\"reversed\":"
                    << (f.surface.reversed ? "true" : "false")
                    << ",\"solid\":" << f.surface.solid_id
                    << ",\"shell\":" << f.surface.shell_id
                    << ",\"source_entries\":[";
          for (std::size_t j = 0; j < f.surface.source_entries.size(); ++j) {
            if (j)
              std::cout << ',';
            std::cout << f.surface.source_entries[j];
          }
          std::cout << "],\"boundary_edges\":[";
          for (std::size_t j = 0; j < f.boundary_edges.size(); ++j) {
            if (j)
              std::cout << ',';
            std::cout << '[';
            for (std::size_t k = 0; k < f.boundary_edges[j].size(); ++k) {
              if (k)
                std::cout << ',';
              std::cout << f.boundary_edges[j][k];
            }
            std::cout << ']';
          }
          std::cout << ']' << ",\"vertex_offset\":" << f.vertex_offset
                    << ",\"vertex_count\":" << f.vertex_count
                    << ",\"triangle_offset\":" << f.triangle_offset
                    << ",\"triangle_count\":" << f.triangle_count
                    << ",\"extra_uv_offset\":" << f.extra_uv_offset
                    << ",\"extra_uv_count\":" << f.extra_uv_count
                    << ",\"identity_mode\":";
          Quote(std::cout, f.identity_mode);
          std::cout << '}';
        }
        std::cout << "],\"mesh\":";
      }
      DumpMesh(result.mesh, std::cout);
      if (mode == "api_model")
        std::cout << '}';
      std::cout << '\n';
      return 0;
    }
    if (mode == "pole_axes") {
      auto pole = DetectPole(f.surface);
      if (!pole)
        throw std::invalid_argument("No pole");
      auto axes = ComputePoleAxes(f.surface, Effective(o), *pole);
      std::cout << "{\"pole\":" << *pole
                << ",\"angle\":" << std::setprecision(17) << axes.angle;
      for (const auto& field :
           {std::pair<const char*, const std::vector<double>*>{"angular",
                                                               &axes.angular},
            {"radial", &axes.radial},
            {"lengths", &axes.lengths}}) {
        std::cout << ",\"" << field.first << "\":[";
        for (std::size_t i = 0; i < field.second->size(); ++i) {
          if (i)
            std::cout << ',';
          Value(std::cout, (*field.second)[i]);
        }
        std::cout << ']';
      }
      std::cout << "}\n";
    } else if (mode == "evaluate") {
      ArcMap arc(f.surface);
      std::vector<Vec3> points, evaluation, control, du, dv;
      for (auto uv : queries) {
        points.push_back(f.surface.Point(uv));
        evaluation.push_back(f.surface.Evaluate(uv));
        control.push_back(f.surface.EvaluateControl(uv));
        auto d = f.surface.Derivatives(uv);
        du.push_back(d[0]);
        dv.push_back(d[1]);
      }
      std::cout << "{\"point\":";
      Arrays(std::cout, points);
      std::cout << ",\"evaluate\":";
      Arrays(std::cout, evaluation);
      std::cout << ",\"control\":";
      Arrays(std::cout, control);
      std::cout << ",\"du\":";
      Arrays(std::cout, du);
      std::cout << ",\"dv\":";
      Arrays(std::cout, dv);
      std::cout << ",\"normals\":";
      Arrays(std::cout, MeshNormals(f.surface, queries));
      std::cout << ",\"maps\":[";
      for (int a = 0; a < 2; ++a) {
        if (a)
          std::cout << ',';
        std::cout << "{\"source\":[";
        for (std::size_t i = 0; i < arc.maps[a].source.size(); ++i) {
          if (i)
            std::cout << ',';
          Value(std::cout, arc.maps[a].source[i]);
        }
        std::cout << "],\"target\":[";
        for (std::size_t i = 0; i < arc.maps[a].target.size(); ++i) {
          if (i)
            std::cout << ',';
          Value(std::cout, arc.maps[a].target[i]);
        }
        std::cout << "]}";
      }
      std::cout << "]}\n";
    } else if (mode == "mesh" || mode == "iges" || mode == "import_mesh") {
      auto mesh = Generate(f, Effective(o));
      DumpMesh(mesh, std::cout);
      std::cout << '\n';
    } else {
      throw std::invalid_argument("Unknown probe mode: " + mode);
    }
    return 0;
  } catch (const Unsupported& e) {
    std::cerr << e.what() << '\n';
    return 3;
  } catch (const std::invalid_argument& e) {
    std::cerr << e.what() << '\n';
    return 2;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 4;
  }
}

#ifdef _WIN32
int wmain(int argc, wchar_t** argv) {
  std::vector<std::string> arguments;
  for (int argument_index = 0; argument_index < argc; ++argument_index) {
    int count = WideCharToMultiByte(CP_UTF8, 0, argv[argument_index], -1,
                                    nullptr, 0, nullptr, nullptr);
    std::string argument(count, '\0');
    WideCharToMultiByte(CP_UTF8, 0, argv[argument_index], -1, argument.data(),
                        count, nullptr, nullptr);
    argument.pop_back();
    arguments.push_back(std::move(argument));
  }
  std::vector<char*> pointers;
  for (auto& argument : arguments)
    pointers.push_back(argument.data());
  return Run(argc, pointers.data());
}
#else
int main(int argc, char** argv) {
  return Run(argc, argv);
}
#endif
