#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <future>
#include <iostream>
#include <set>
#if defined(_MSC_VER) && defined(_DEBUG)
#include <crtdbg.h>
#endif
#include "tessellation/private/internal.h"
extern "C" {
#define REAL double
#define VOID void
#include "triangle.h"
int TriangleExecute(char*, triangulateio*, triangulateio*);
}
using namespace iko;
static void Require(bool value, const char* message) {
  if (!value)
    throw std::runtime_error(message);
}
static std::uint32_t Word(std::istream& in) {
  unsigned char bytes[4];
  in.read(reinterpret_cast<char*>(bytes), 4);
  return std::uint32_t(bytes[0]) | (std::uint32_t(bytes[1]) << 8) |
         (std::uint32_t(bytes[2]) << 16) | (std::uint32_t(bytes[3]) << 24);
}
static float Scalar(std::istream& in) {
  auto bits = Word(in);
  float x;
  std::memcpy(&x, &bits, 4);
  return x;
}
int main(int argc, char** argv) {
#if defined(_MSC_VER) && defined(_DEBUG)
  // A CRT assertion must fail an unattended test, rather than wait in a dialog.
  _set_abort_behavior(0, _WRITE_ABORT_MSG | _CALL_REPORTFAULT);
  _CrtSetReportHook([](int type, char* message, int*) -> int {
    if (type == _CRT_ERROR || type == _CRT_ASSERT) {
      std::fputs(message, stderr);
      std::fflush(stderr);
      std::abort();
    }
    return 0;
  });
#endif
  try {
    Require(argc == 3, "contract requires test data and output directories");
    std::filesystem::path root = argv[1], out = argv[2];
    std::filesystem::create_directories(out);
    std::cout << "Checking two-plane API, selection and concurrency\n"
              << std::flush;
    // Two independent square faces, one in XY and one in YZ.
    auto file = root / "two_planes.igs";
    auto listed = ListSurfaces(file);
    Require(listed.status == Status::kSuccess && listed.surfaces.size() == 2,
            "Surface listing");
    auto whole = Tessellate(file);
    Require(whole.Succeeded() && whole.complete && whole.whole_model &&
                whole.faces.size() == listed.surfaces.size() &&
                whole.model.required_faces == 2 &&
                std::all_of(whole.faces.begin(), whole.faces.end(),
                            [](const auto& face) {
                              return face.vertex_count > 0 &&
                                     face.triangle_count > 0;
                            }),
            "Default invocation generates all required model faces");
    Options o;
    o.surface_de = listed.surfaces.front().directory_entry;
    auto r = Tessellate(file, o);
    Require(r.Succeeded() && r.complete && !r.whole_model &&
                r.faces.size() == 1 && !r.mesh.triangles.empty(),
            r.message.c_str());
    Require(r.mesh.face_id == listed.surfaces.front().instance_id,
            "Source face identity");
    o.surface_de = listed.surfaces.back().directory_entry;
    auto second = Tessellate(file, o);
    Require(second.Succeeded() && second.complete && !second.whole_model &&
                second.faces.size() == 1 && !second.mesh.triangles.empty() &&
                second.mesh.face_id == listed.surfaces.back().instance_id,
            "The second source DE selects only the second plane");
    o.surface_de = 2147483647;
    Require(Tessellate(file, o).status == Status::kInputProblem,
            "Missing selected DE");
    Require(Tessellate(out / "missing.igs").status == Status::kInputProblem,
            "Missing input");
    std::ofstream(out / "invalid.igs") << "invalid";
    Require(Tessellate(out / "invalid.igs").status == Status::kInputProblem,
            "Malformed IGES");
    o = {};
    o.surface.chord = 0;
    Require(Tessellate(file, o).status == Status::kInputProblem,
            "Invalid tolerances");
    o = {};
    o.surface_de = listed.surfaces.front().directory_entry;
    std::vector<std::future<Result>> calls;
    for (int i = 0; i < 4; ++i)
      calls.push_back(
          std::async(std::launch::async, [&] { return Tessellate(file, o); }));
    for (auto& call : calls) {
      auto other = call.get();
      Require(other.Succeeded() && other.mesh.xyz == r.mesh.xyz &&
                  other.mesh.triangles == r.mesh.triangles &&
                  other.mesh.identities == r.mesh.identities,
              "Concurrent calls retain independent results");
    }

    // ShapeFix can return unlocated transfer binders for placed source faces.
    // Each independently declared face must retain its individual source DE,
    // so selecting it cannot accidentally select a group of unrelated faces.
    std::cout << "Checking placed source faces\n" << std::flush;
    auto placed_file = root / "motor_faces.igs";
    auto placed_list = ListSurfaces(placed_file);
    Require(placed_list.status == Status::kSuccess &&
                placed_list.surfaces.size() == 223,
            "Complete placed face listing");
    std::set<int> placed_sources;
    for (const auto& surface : placed_list.surfaces) {
      Require(surface.directory_entry > 1 &&
                  surface.source_entries.size() == 1 &&
                  surface.source_entries.front() == surface.directory_entry,
              "Placed face has its individual source entry");
      placed_sources.insert(surface.directory_entry);
    }
    Require(placed_sources.size() == placed_list.surfaces.size(),
            "Placed source face identities are unique");
    Options placed_options;
    placed_options.surface_de = placed_list.surfaces.front().directory_entry;
    auto placed_result = Tessellate(placed_file, placed_options);
    Require(placed_result.Succeeded() && placed_result.complete &&
                !placed_result.whole_model && placed_result.faces.size() == 1 &&
                placed_result.faces.front().surface.instance_id ==
                    placed_list.surfaces.front().instance_id,
            "Placed source DE selects its own face");

    std::cout << "Checking periodic surfaces and domain boundaries\n"
              << std::flush;
    auto periodic_file = root / "periodic_eight.igs";
    auto periodic_list = ListSurfaces(periodic_file);
    Require(periodic_list.status == Status::kSuccess &&
                periodic_list.surfaces.size() == 1 &&
                periodic_list.surfaces[0].trimmed,
            "Unified listing exposes the constructed natural rectangle");
    Options periodic_options;
    periodic_options.surface = {.05, 3.13, .063};
    periodic_options.curve = {.05, 5.27, .093};
    auto periodic_result = Tessellate(periodic_file, periodic_options);
    Require(periodic_result.Succeeded() &&
                !periodic_result.repair_trace.empty() &&
                periodic_result.mesh.face_id ==
                    periodic_list.surfaces[0].instance_id &&
                periodic_result.requested_options.surface.max_edge == 3.13 &&
                periodic_result.effective_options.curve.chord == .093 &&
                periodic_result.mesh.explicit_weights &&
                std::any_of(
                    periodic_result.repair_trace.begin(),
                    periodic_result.repair_trace.end(),
                    [](const auto& s) {
                      return s.find("active-domain endpoint knot insertion") !=
                             std::string::npos;
                    }),
            "Periodic public SDK call, independent tolerances, adaptation and "
            "source identity");
    periodic_options.surface_de = periodic_list.surfaces[0].directory_entry;
    auto selected = Tessellate(periodic_file, periodic_options);
    Require(selected.Succeeded() &&
                selected.mesh.xyz == periodic_result.mesh.xyz &&
                selected.mesh.uv == periodic_result.mesh.uv &&
                selected.mesh.normals == periodic_result.mesh.normals &&
                selected.mesh.triangles == periodic_result.mesh.triangles &&
                selected.mesh.identities == periodic_result.mesh.identities,
            "Default and explicit single-face selection agree");
    auto seven = Tessellate(root / "periodic_seven.igs", periodic_options);
    Require(seven.Succeeded() &&
                std::any_of(
                    seven.repair_trace.begin(), seven.repair_trace.end(),
                    [](const auto& s) {
                      return s.find("active-domain endpoint knot insertion") !=
                             std::string::npos;
                    }),
            "Seven-span adaptation from OCCT geometry");
    auto b_source = internal::ReadOcct(periodic_file).faces[0];
    internal::ValidateNaturalRectangle(b_source);
    auto bad_rectangle = b_source;
    (*bad_rectangle.outer)[0].control[1][0] += .01;
    bool rejected = false;
    try {
      internal::ValidateNaturalRectangle(bad_rectangle);
    } catch (const std::invalid_argument&) {
      rejected = true;
    }
    Require(rejected, "Rectangle projection accepted a changed boundary");
    const auto& bs = b_source.surface;
    Vec2 outside{std::nextafter(bs.domain[0], -INFINITY),
                 std::nextafter(bs.domain[3], INFINITY)};
    Require(internal::RectangleSamplingPoint(bs, outside) ==
                bs.Evaluate({bs.domain[0], bs.domain[3]}),
            "Natural rectangle physical sampling projection ordering");
    rejected = false;
    try {
      internal::RectangleSamplingPoint(bs, {NAN, 0.});
    } catch (const std::invalid_argument&) {
      rejected = true;
    }
    Require(rejected, "Nonfinite Natural rectangle physical sample admission");

    // A Triangle fatal input must return across the C boundary. The next
    // ordinary call proves that the adapter and callback state recovered.
    std::cout << "Checking backend recovery and mesh export\n" << std::flush;
    double points[] = {0., 0., 1., 0.};
    triangulateio input{}, output{};
    input.pointlist = points;
    input.numberofpoints = 2;
    char switches[] = "Qz";
    Require(TriangleExecute(switches, &input, &output) != 0,
            "Triangle fatal input did not return failure");
    internal::PlanarStraightLineGraph p;
    p.uv = {{0., 0.}, {1., 0.}, {0., 1.}};
    Require(internal::TriangleRun(p, Options{}, {}, "Qz").triangles.size() == 1,
            "Backend recovery after fatal input");
    p.segments = {{0, 1}, {1, 2}, {2, 0}};
    bool propagated = false;
    try {
      internal::TriangleRun(p, Options{}, [](const auto&) -> bool {
        throw std::runtime_error("callback failure");
      });
    } catch (const std::runtime_error&) {
      propagated = true;
    }
    Require(propagated, "Callback exception propagation");

    Mesh m;
    m.xyz = {
        {0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}, {9., 9., 9.}};
    m.uv = {{0., 0.}, {1., 0.}, {0., 1.}, {0., 0.}, {.5, .5}};
    m.normals = std::vector<Vec3>(5, {0., 0., 1.});
    m.identities = {"seamA", "b", "c", "seamB", "unused"};
    m.triangles = {{0, 1, 2}, {0, 1, 2}, {2, 1, 0}, {0, 3, 1}};
    m.unpaired_uv_tail = {{.75, .25}};
    auto original = m;
    Require(WriteObj(m, out / "records.obj").status == Status::kSuccess,
            "OBJ export");
    Require(WriteBinaryStl(m, out / "records.stl").status == Status::kSuccess,
            "STL export");
    Require(m.xyz == original.xyz && m.uv == original.uv &&
                m.normals == original.normals &&
                m.triangles == original.triangles &&
                m.identities == original.identities &&
                m.unpaired_uv_tail == original.unpaired_uv_tail,
            "Export changed standard mesh");
    std::ifstream obj(out / "records.obj");
    std::string line;
    int v = 0, vt = 0, vn = 0, f = 0;
    std::vector<std::string> facets;
    while (std::getline(obj, line)) {
      if (line.rfind("v ", 0) == 0)
        ++v;
      else if (line.rfind("vt ", 0) == 0)
        ++vt;
      else if (line.rfind("vn ", 0) == 0)
        ++vn;
      else if (line.rfind("f ", 0) == 0) {
        ++f;
        facets.push_back(line);
      }
    }
    Require(v == 5 && vt == 5 && vn == 5 && f == 4 && facets[0] == facets[1] &&
                facets[2] == "f 3/3/3 2/2/2 1/1/1" &&
                facets[3] == "f 1/1/1 4/4/4 2/2/2",
            "OBJ lost duplicates, unused records, winding or degeneracy");
    Require(std::filesystem::file_size(out / "records.stl") ==
                84 + 50 * m.triangles.size(),
            "Binary STL size");
    std::ifstream stl(out / "records.stl", std::ios::binary);
    stl.seekg(80);
    Require(Word(stl) == 4, "STL facet count");
    for (std::size_t ti = 0; ti < m.triangles.size(); ++ti) {
      Vec3 n;
      for (auto& x : n)
        x = Scalar(stl);
      Require(n == (ti == 2   ? Vec3{0., 0., -1.}
                    : ti == 3 ? Vec3{}
                              : Vec3{0., 0., 1.}),
              "STL normal/winding");
      for (auto i : m.triangles[ti])
        for (double x : m.xyz[i])
          Require(Scalar(stl) == float(x), "STL facet coordinates/order");
      char attribute[2];
      stl.read(attribute, 2);
      Require(attribute[0] == 0 && attribute[1] == 0, "STL attribute word");
    }
    Require(WriteObj(m, out / "absent/mesh.obj").status ==
                Status::kComputationFailed,
            "Export I/O status");
    m.triangles[0][0] = 100;
    Require(WriteBinaryStl(m, out / "bad.stl").status == Status::kInputProblem,
            "Invalid mesh status");
    std::cout << "SDK calls, selection, concurrency, backend failure/recovery "
                 "and lossless "
                 "mesh export contract passed\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
