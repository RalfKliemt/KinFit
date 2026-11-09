# Installation

## Requirements

- C++17 compiler (the C++ standard is taken over from the ROOT build)
- CMake ≥ 3.10
- ROOT 6 with the components Core, Matrix, Physics, MathCore, RIO and Tree.
  The build is tested with ROOT 6.40.

## Build and install

```bash
git clone https://github.com/KinFit/KinFit.git
cd KinFit
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/your/install/prefix
cmake --build build -j
cmake --install build
```

This produces `libKinFit` (`.so` on Linux, `.dylib` on macOS) together with its
ROOT dictionary (`libKinFit_rdict.pcm`, `libKinFit.rootmap`). The library is
installed to `<prefix>/lib`, the headers to `<prefix>/include`, and the CMake
package files to `<prefix>/lib/cmake/KinFit`.

### Build options

| Option | Default | Effect |
| --- | --- | --- |
| `KINFIT_BUILD_EXAMPLES` | `OFF` | Build `toy_lambda_fit` and `vertex_constB` from `docs/examples/` |
| `CMAKE_INSTALL_PREFIX` | system default | Installation prefix |

### After a ROOT upgrade

The CMake cache stores absolute paths into the ROOT installation. If ROOT is
updated (for example by a package manager), configuration fails with errors
like `.../root/6.40.00/include/root/vdt/vdtMath.h cannot be read`. Delete the
build directory and configure again.

## Using KinFit from another CMake project

```cmake
find_package(ROOT REQUIRED COMPONENTS Core Physics Tree)
find_package(KinFit REQUIRED)        # add <prefix> to CMAKE_PREFIX_PATH

add_executable(my_analysis my_analysis.cxx)
target_link_libraries(my_analysis PRIVATE KinFit::KinFit)
```

The `KinFit::KinFit` target carries the include directory and the ROOT
libraries it needs.

## Using KinFit in a ROOT session or macro

Load the library once, for example in your `rootlogon.C`:

```cpp
gSystem->Load("/your/install/prefix/lib/libKinFit");
```

The dictionary makes all KinFit classes available at the prompt and in
interpreted macros without any `#include`. If a macro should also compile
as a standalone program, guard the KinFit includes so that the interpreter
does not see them twice:

```cpp
#ifndef __CLING__
#include "KinFitter.h"
#endif
```

The examples in [examples/](examples/) follow this pattern.
