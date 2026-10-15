# IKO | Tessellation

Tessellate IGS models and export OBJ or binary STL files.

## Use the SDK

Add this to your `CMakeLists.txt`:

```cmake
find_package(iko CONFIG REQUIRED COMPONENTS tessellation)
target_link_libraries(my_app PRIVATE iko::tessellation)
```

C++ example:

```cpp
#include <tessellation/tessellation.h>
auto result = iko::Tessellate("model.igs");
if (result.Succeeded() && result.complete) {
  iko::WriteObj(result.mesh, "model.obj");
}
```

Check `result.status` and `result.message` for errors. See the [public header](include/tessellation/tessellation.h) for the full API.

To try the example, open a VS x64 developer shell at the SDK repository root, set the installation path, and run:

```powershell
$IKO_INSTALL_ROOT = "<absolute path to the iko installation directory>"
cmake -S examples -B .runs/consumer -G Ninja `
  -DCMAKE_BUILD_TYPE=Release `
  "-DCMAKE_PREFIX_PATH=$IKO_INSTALL_ROOT"
cmake --build .runs/consumer
./.runs/consumer/tessellation_consumer.exe model.igs
```

Replace `model.igs` with your input file. For Debug, replace `Release` with `Debug` and `.runs/consumer` with `.runs/consumer-debug`.

## Build and install

Install the Visual Studio C++ build tools, Windows SDK, Git, and CMake/Ninja. Open PowerShell at the SDK repository root. Before your first build, download the OCCT source:

```powershell
git clone --depth 1 --branch V8_0_1 --single-branch `
  https://github.com/Open-Cascade-SAS/OCCT.git `
  .runs/deps/occt-8.0.1
```

If you already have the source, start here. These commands build, test, and install both Release and Debug:

```powershell
$IKO_INSTALL_ROOT = "<absolute path to the iko installation directory>"
./dependencies/build_occt.ps1
./build.ps1 -PrefixPath .runs/occt-8.0.1-install-dynamic `
  -RunTests -InstallRoot $IKO_INSTALL_ROOT

./dependencies/build_occt.ps1 -BuildType Debug
./build.ps1 -BuildType Debug `
  -PrefixPath .runs/occt-8.0.1-install-dynamic `
  -RunTests -InstallRoot $IKO_INSTALL_ROOT
```

The dependency script patches OCCT 8.0.1 and keeps its build files under the SDK repository's `.runs` directory. SDK builds go to `.runs/build` for Release and `.runs/build-debug` for Debug. If you already have a patched OCCT build for the required configuration, skip the dependency step and pass its installation directory to `-PrefixPath`.

`-InstallRoot` is the final installation directory, corresponding to CMake's `IKO_INSTALL_ROOT`. Omit it to build without installing. To skip tests, replace `-RunTests` with `-NoTests`.

## Installation layout

Release and Debug can share an installation directory. CMake selects the files for the current configuration. The `.lib` files are DLL import libraries.

```text
<IKO_INSTALL_ROOT>/
  bin/release/   # iko_tessellation.dll, tessellate.exe
  bin/debug/     # iko_tessellationd.dll, tessellate.exe
  lib/release/   # iko_tessellation.lib
  lib/debug/     # iko_tessellationd.lib
  lib/cmake/     # CMake packages for iko/ and iko_tessellation/
  include/tessellation/tessellation.h
  share/iko/tessellation/   # Third-party licenses, source, and patches
```

Only public headers are installed. Reinstalling removes obsolete files recorded for the current module and configuration while preserving other modules, configurations, user files, and shared files still in use. `-Clean` only cleans build artifacts for the current configuration.

See the [test guide](tests/README.md) for tests and C++ interface requirements, and the [validation record](tests/validation.json) for verified results. Third-party materials are listed in the [license notices](THIRD_PARTY_NOTICES.md) and [OCCT patch](dependencies/occt_projection-8.0.1.patch).

## Topology-consistent meshing (opt-in)

The default remains record-exact: faces are sampled independently and vertex
records are never welded. To coordinate an imported BRep's shared edges:

```cpp
iko::Options options;
options.topology_consistent = true;
auto result = iko::Tessellate("model.igs", options);
```

The CLI equivalent is `tessellate model.igs --topology-consistent --out model.obj`.
This mode builds one ordered XYZ sample registry per OCCT topological edge,
unions sampling and surface-knot requests from the selected incident faces,
and maps the resulting samples onto each directed pcurve independently. Shared
edge endpoints use the canonical OCCT vertex position. It retains explicit
wires, including natural rectangles, instead of using the legacy trim cleanup.
The triangulator cannot independently insert boundary vertices; boundary
segments are checked before any result is published. If a fixed boundary blocks
Triangle's circumcenter insertion, bounded interior-centroid refinement retries
the same constraints while retaining all previous vertices. Automatic mode caps
insertions at 100,000 per face; an explicit `steiner_limit` covers all backend and
manual interior insertions.

UVs, normals, and triangle indices remain face-local records. In this mode,
`mesh.geometric_vertices` supplies one ID per vertex record: equal IDs represent
the same geometric vertex, including copies on adjacent faces, seam occurrences,
and surface-knot normal discontinuities. IDs are scoped to that result and are
not stable across calls/options. OBJ/STL exports still preserve separate records;
applications needing geometric connectivity can use these IDs without merging
UVs or sharp normals. Topology-mode coordinates retain double precision.

Only existing OCCT edge/vertex identities are shared. Coincident independent
surfaces are not automatically sewn. A diagnostic `surface_de` selection
coordinates only the selected faces. The support's parametric normal and winding
convention is unchanged; `FaceRange.surface.reversed` remains separate metadata.

This is a conservative, opt-in implementation. Degenerate pole edges, missing or
synthetic topology, ambiguous/nonmonotone XYZ-to-pcurve mappings, disconnected
wires, unresolved knot constraints, and singular normals return an explicit
failure with no partial mesh. Valid repeated seam uses retain their distinct UV
occurrences. Sampling is adaptive and numerical, using the existing tolerance
and minimum-edge heuristics plus finite mapping probes; it is not a proof of a
strict global chord-error bound. Fixed boundaries can prevent requested interior
refinement; exhaustion is reported rather than silently returning an incomplete
refinement. Rebuild SDK consumers when updating these public value-type layouts.
