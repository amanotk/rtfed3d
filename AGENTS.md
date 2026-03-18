# AGENTS.md

This file is for coding agents working in `/home/amano/rtfed3d`.

## Project Summary

- `rtfed3d` is a C++17 scientific computing code for relativistic two-fluid electrodynamics.
- The build system is CMake.
- MPI and parallel HDF5 are required.
- `external/mdspan/` contains a vendored upstream `kokkos/mdspan` tree.
- `external/` also contains vendored single-header third-party dependencies.
- `examples/` contains runnable problem setups, not unit tests.
- `tests/` contains regression infrastructure and golden HDF5 snapshots.

## Agent Rules Sources

- There is currently no repo-local `AGENTS.md` other than this file.
- There is no `.cursorrules` file.
- There is no `.cursor/rules/` directory.
- There is no `.github/copilot-instructions.md` file.
- If any of those are added later, treat them as higher-priority repo guidance and merge them into your behavior.

## Repository Layout

- `include/`: project headers
- `src/`: shared solver implementation
- `examples/`: simulation drivers and `.cfg` inputs
- `tests/`: regression helpers and reference HDF5 outputs
- `external/`: vendored third-party headers
- `external/mdspan/`: vendored upstream mdspan headers and metadata
- `out/`: generated outputs; ignored by git
- `build/`, `build-*`: generated build trees; ignored by git

## Build Commands

Use the repo root as the working directory unless you have a good reason not to.

### Standard configure/build

```bash
cmake -S . -B build
cmake --build build
```

### Optimized build

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release
```

### Build a single target

```bash
cmake --build build --target riemann1d
cmake --build build --target diffusion
```

### Reconfigure after CMake changes

```bash
cmake -S . -B build
```

## Test Commands

CTest is enabled from the main `CMakeLists.txt`.

### Run the full test suite

```bash
ctest --test-dir build --output-on-failure
```

### Run a single test by exact or regex name

```bash
ctest --test-dir build -R diffusion_small_regression --output-on-failure
ctest --test-dir build -R diffusion_small_mpi4_regression --output-on-failure
ctest --test-dir build -R riemann1d_help --output-on-failure
```

### List available tests

```bash
ctest --test-dir build -N
```

### Re-run only failed tests

```bash
ctest --test-dir build --rerun-failed --output-on-failure
```

## Regression Test Notes

- The current golden regression uses `examples/diffusion_small.cfg`.
- The reference snapshot is `tests/data/diffusion_small_reference.h5`.
- The same golden file is used for both the single-process and 4-rank MPI regression.
- `riemann1d_help` is a smoke test only; it is not the current golden physics regression.
- If you change solver numerics, array layout, MPI decomposition, or HDF5 output semantics, re-run the tests before touching any golden files.
- Do not refresh golden snapshots casually. Only do so when the behavior change is intentional and understood.

## Running Examples Manually

Example executables are emitted into `build/examples/` or `build-release/examples/`.

Typical workflow:

```bash
cmake -S . -B build
cmake --build build --target riemann1d
cd build/examples
./riemann1d
```

Override the config path explicitly when needed:

```bash
./riemann1d -c /absolute/path/to/examples/riemann1d.cfg
```

For MPI runs, use the executable's domain flags:

```bash
mpiexec -n 4 ./diffusion -c /absolute/path/to/examples/diffusion_small.cfg -x 4
```

## Lint / Formatting

- There is no repo-wide `lint` target in CMake.
- There is no committed clang-format or clang-tidy config for the project sources.
- Do not invent a formatting regime and do not mass-reformat files.
- Preserve the existing local style in touched files.
- Python helper scripts in `tests/` follow modern PEP 8-ish formatting, but there is no committed formatter config.

## C++ Style Guidelines

Follow the style already present in the codebase.

### Headers and includes

- Use quoted includes for project headers, e.g. `#include "global.hpp"`.
- Keep project support headers in use through the configured include paths rather than relative paths.
- Do not add unnecessary includes.
- Prefer keeping include order stable with nearby files rather than applying an external convention.

### Formatting and layout

- Use the existing brace style: opening braces typically go on the next line for classes/functions.
- Keep the current spacing style, e.g. `if(this != &other)` and `for(int i=0; i < n ;i++)`.
- Do not reformat legacy loops just for style.
- Match the file's existing comment density.
- Preserve the file footer style if it already has local variable comments.

### Types

- Use project aliases such as `float64`, `int32`, `T_vector`, `T_tensor`, and related existing types.
- For core arrays, prefer the project wrapper types defined through `include/decls.hpp` rather than raw `std::vector` or raw mdspan in calling code.
- Keep new array/storage changes compatible with the current mdspan-backed `Array` wrapper unless intentionally refactoring that layer.
- Use `const` where the surrounding code already does.

### Naming

- Match existing naming rather than modernizing names opportunistically.
- Classes use `CamelCase`, e.g. `Global`, `PeriodicBoundary`, `MC2`.
- Methods and free functions use lower-case or mixed legacy names, e.g. `energy`, `get_dt_factor`, `xflux_fluid_hll`.
- Local loop indices are typically `ix`, `iy`, `iz`, `k`, `n`, `ibc`.
- Do not rename widely used physics variables just for readability.

### Error handling

- Existing C++ code mostly does not use exceptions for control flow.
- Prefer explicit checks and existing failure mechanisms.
- For configure-time dependency problems, use clear CMake `message(FATAL_ERROR ...)` diagnostics.
- For Python regression helpers, raising `AssertionError` or standard exceptions with precise messages is the established pattern.

### Numerical/physics changes

- Be conservative with any solver or boundary-condition change.
- Small indexing changes can alter the physics result or MPI consistency.
- Validate with `ctest` after touching solver, boundary, integrator, array, or HDF5 code.
- If a change affects outputs, inspect whether the difference is expected before updating snapshots.

## Python Style Guidelines

- Keep helper scripts simple and dependency-light.
- Use `pathlib` for paths when editing existing modern scripts such as `tests/run_and_compare_hdf5.py`.
- Prefer explicit error messages over silent failures.
- Do not introduce heavyweight test frameworks unless there is a clear reason.

## Git and Change Hygiene

- Never commit generated `build/`, `build-*`, or `out/` artifacts.
- Avoid editing vendored code under `external/` unless the task is specifically about updating the vendored dependency.
- If vendoring is updated, also update `external/mdspan/VENDORED.md`.
- Keep commits scoped: build-system changes, regression updates, and solver changes should usually be separate commits.

## What to Check Before Finishing

- Did CMake still configure successfully?
- Did the touched target(s) rebuild?
- Did `ctest --test-dir build --output-on-failure` pass if behavior-affecting code changed?
- Did you avoid modifying ignored/generated artifacts?
- Did you preserve the project's legacy style where possible instead of rewriting it?
