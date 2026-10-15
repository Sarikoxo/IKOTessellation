# Run the tests

Build the dependencies for the required configuration as described in the [main README](../README.md), then run these commands at the SDK repository root:

```powershell
./build.ps1 -PrefixPath .runs/occt-8.0.1-install-dynamic `
  -RunTests
./build.ps1 -BuildType Debug `
  -PrefixPath .runs/occt-8.0.1-install-dynamic `
  -RunTests
```

The default Windows build has 7 CTest tests. Output goes to `tests/output` inside each SDK build directory. See [validation.json](validation.json) for verified configurations and results.

## Coverage

- `sdk_topology_import`: public SDK integration on the prismatic BRep and a cylinder with shared circular edges/repeated seam; closed geometric edge incidence, canonical XYZ, atomic failure, and legacy recovery.
- `sdk_topology`: 14 synthetic cases covering shared straight/curved edges, mismatched/nonlinear UV parameterizations, opposite uses, unilateral/near-endpoint knots, shared corners, sharp and intra-face C0 normals, valid repeated seams, unrelated coincident faces, explicit ambiguous-seam/pole rejection, explicit insertion-budget failure/recovery, and default compatibility.

- `sdk_contract`: models, face selection, concurrent calls, periodic surfaces, source identities, exports, and error recovery.
- `sdk_exact`: double-precision addition, subtraction, multiplication, and division.
- `sdk_cli`: command-line arguments, paths containing Chinese characters, and OBJ/STL exports.
- `sdk_package`: installed consumers, objects returned by the DLL, component lookup, and installation cleanup. Source builds still reject unpatched OCCT.
- `sdk_dll`: only iko interfaces are exported, no OCCT DLL is required, and the `.lib` is an import library.

## C++ interface requirements

The SDK, its dependencies, and the application must use the same architecture, MSVC toolset, and default STL and structure layouts. Match the settings below:

| Configuration | Runtime | `_ITERATOR_DEBUG_LEVEL` |
| --- | --- | --- |
| Release | `/MD` | 0 |
| Debug | `/MDd` | 2 |

`/MT` and `/MTd` are unsupported. Do not mix Release and Debug.

Returned strings and containers can be copied, modified, and destroyed normally. The application owns input paths and meshes. The application and SDK should use consistent memory allocation and deallocation methods. Before unloading the DLL, wait for all calls to finish and destroy returned objects. Do not free the pointer returned by `StatusName` or use it after unloading the DLL.

The `topology_cylinder.igs` fixture was generated with patched OCCT 8.0.1,
using `IGESControl_Writer("MM", 1)` and
`BRepPrimAPI_MakeCylinder(2., 3.).Shape()` (radius 2 mm, height 3 mm).
The public import tests require the supported Windows build because the SDK's
production SHA-256 implementation intentionally rejects other platforms.
The historical `validation.json` records the pre-topology Windows baseline; it
is not a claim that the new mode has passed that entire platform matrix.

See [topology validation](topology_validation.md) for the new mode's actual
checks, commands, and supported-platform verification limits.
