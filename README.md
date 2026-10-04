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
