# Topology-mode validation, 2026-10-15

This records checks for the opt-in topology implementation. It does not replace
`validation.json`, which describes the earlier Windows baseline.

## Production build and tests

Built the SDK on Linux with GCC, C++17, a static/PIC OCCT 8.0.1 installation
including the repository's projection-policy patch, and the existing CMake
configuration. No production hash/platform code was changed.

Commands, from the repository root (CMake/CTest must be on `PATH`):

```sh
cmake -S . -B .runs/build-linux -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PWD/.runs/occt-8.0.1-install-linux"
cmake --build .runs/build-linux --parallel 4
ctest --test-dir .runs/build-linux \
  -R '^(sdk_topology|sdk_exact)$' --output-on-failure
```

- Production SDK, CLI, probe, and test executables compiled successfully.
- `sdk_topology`: all 14 synthetic cases passed, including curved shared edges,
  repeated seams, nonlinear parameter maps, near-end knots, sharp normals,
  limited-budget failure/recovery, and legacy default compatibility.
- `sdk_exact`: passed.
- The complete Linux CTest run passed 2 of 6 tests. `sdk_topology_import`,
  `sdk_contract`, and `sdk_cli` encounter the existing Windows-only SHA-256
  guard. `sdk_package` assumes Windows `.lib` cleanup artifacts and fails its
  obsolete-flat-layout check on Linux. These are not a full supported-platform
  validation. Windows DLL/ABI and complete Windows CTest checks remain unrun.

## Isolated real-import checks

A separate, untracked Linux executable linked the actual OCCT importer and
meshing objects, substituting only a clearly labeled constant source-identity
function for the platform-blocked hash. This does **not** validate SHA-256 or
the installed DLL boundary and is not a production workaround.

The same approach ran the exact committed `topology_import_test.cpp`, linked
against `src/tessellation.cpp`, `src/mesh_export.cpp`, the core archive, and OCCT:

```sh
c++ -O1 -std=c++17 -ffp-contract=off -fno-fast-math -Iinclude \
  tests/topology_import_test.cpp src/tessellation.cpp src/mesh_export.cpp \
  .runs/topology/identity_stub.cpp \
  .runs/build-linux/tests/libiko_tessellation_test_support.a \
  -Wl,--start-group .runs/occt-8.0.1-install-linux/lib/*.a -Wl,--end-group \
  -pthread -ldl -lrt -lm -o .runs/topology/public_import_identity_stub
.runs/topology/public_import_identity_stub tests/data
```

The identity stub above is a validation-only local file, deliberately not
included in the SDK or test suite. The normal `sdk_topology_import` test uses
the real public SDK and production hash on the supported Windows build.

Verified results:

- `prismatic_brep.igs`: 10 faces, 3,254 triangles, 376 cross-face geometric IDs,
  zero open or nonmanifold geometric segments.
- `topology_cylinder.igs`: 3 faces, 1,116 triangles, 238 cross-face geometric IDs,
  zero open or nonmanifold segments; distinct seam UV occurrences retained.
- `two_planes.igs`: 2 faces, 104 triangles; unrelated topology remains separate.
- The committed import test's limited-budget atomic failure and subsequent
  legacy recovery checks also passed in that isolated executable.

Triangle counts are observations of this build, not new exact-output contracts
for the opt-in mode. Default record-exact behavior remains a separate contract.

## Sanitizers and review

An independent build of the production topology sources and Triangle with
AddressSanitizer and UndefinedBehaviorSanitizer passed all 14 synthetic cases
without diagnostics. LeakSanitizer cannot run under this environment's ptrace
restriction; leak detection was disabled with `ASAN_OPTIONS=detect_leaks=0`.
This is not a leak-test pass. Independent code review and `git diff --check`
found no remaining confirmed defects at submission.
