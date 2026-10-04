# Triangle 1.6

The SDK statically compiles `third_party/triangle/triangle.c`.
The source retains Jonathan Richard Shewchuk's complete copyright, license,
disclaimer and notices. The existing Win64, queue and work-limit changes are
marked in that file. `src/triangle_backend.c` adds a separate C failure and
allocation boundary without changing the backend's numerical operations.

The applicable license text from the source is:

```text
Copyright 1993, 1995, 1997, 1998, 2002, 2005
Jonathan Richard Shewchuk
2360 Woolsey #H
Berkeley, California 94705-1927
jrs@cs.berkeley.edu

This program may be freely redistributed under the condition that the
copyright notices (including this entire header and the copyright
notice printed when the `-h' switch is selected) are not removed, and
no compensation is received. Private, research, and institutional
use is free. You may distribute modified versions of this code UNDER
THE CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE TO IT IN THE
SAME FILE REMAIN UNDER COPYRIGHT OF THE ORIGINAL AUTHOR, BOTH SOURCE
AND OBJECT CODE ARE MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR
NOTICE IS GIVEN OF THE MODIFICATIONS. Distribution of this code as
part of a commercial system is permissible ONLY BY DIRECT ARRANGEMENT
WITH THE AUTHOR. (If you are not directly supplying this code to a
customer, and you are instead telling them how they can obtain it for
free, then you are not required to make any arrangement with me.)

Neither I nor Carnegie Mellon warrant this code in any way
whatsoever. This code is provided "as-is". Use at your own risk.
```

Keep this notice and the complete backend source header when distributing.
This migration does not grant a new commercial license for Triangle.

# Open CASCADE Technology 8.0.1

The importer links the official OCCT V8_0_1 sources, commit
`b8f597c677811d1f9f4d8a97f5ae2825c0353a42`, with the explicit compatibility
patch `dependencies/occt_projection-8.0.1.patch`.
Source: https://github.com/Open-Cascade-SAS/OCCT/tree/V8_0_1

OCCT is licensed under GNU LGPL 2.1 with the OCCT exception. Complete texts
are preserved in `third_party/occt/` and installed with this notice. This build
statically links the required OCCT toolkit dependency closure into the iko DLL,
using the shared MSVC runtime (/MD Release, /MDd Debug). Only the iko public API
is exported; no OCCT object crosses that API. Distributors must retain
the notices and meet the applicable source/relinking terms. The pinned upstream
source is built using `dependencies/build_occt.ps1`. No OCCT DLL is required by this verified build.

The patch adds a thread-local projection policy, disabled by default. The SDK
enables it only during analytic NURBS conversion and edited-endpoint projection.
It retains 24-point quadrature, ascending fitting degrees, the existing interval
cutting limit and the original normal-root solve for boundary compatibility.
IGES reading, ShapeFix and all other OCCT calls retain the release defaults.
Distribute this patch with the pinned OCCT source and its license texts.
