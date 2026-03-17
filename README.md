# Relativistic Two-Fluid Electrodynamics 3D

## Requirements

- Parallel HDF5 (compiled with `--enable-parallel`)
- MPI C++ compiler wrapper
- CMake 3.20 or newer
- C++17 compiler

## Repository Layout

- `include/`: project headers
- `src/`: shared solver implementation
- `examples/`: runnable problem drivers and example configuration files
- `base/`: vendored support headers used by the solver and examples
- `external/mdspan/`: vendored upstream `kokkos/mdspan` headers

The vendored `mdspan` source tree records its upstream revision in
`external/mdspan/VENDORED.md`.

## Build

### Clone
```bash
$ git clone git@github.com:amanotk/rtfed3d.git
$ cd rtfed3d
```

### Configure
If HDF5 lives in a non-standard location, pass the usual CMake hints for
MPI/HDF5 during configure.

```bash
$ cmake -S . -B build
```

### Compile
```bash
$ cmake --build build
```

Example executables are written to `build/examples/`.

### Test
```bash
$ ctest --test-dir build --output-on-failure
```

The regression suite currently uses a small `diffusion` example for both
single-process and 4-process MPI snapshot checks.

### Run
Each example expects its `.cfg` file in the current working directory. The build
copies matching config files next to each executable, so the simplest workflow is:

```bash
$ cd build/examples
$ ./riemann1d
```

If you want an optimized build, you can add `-DCMAKE_BUILD_TYPE=Release` during
configure, but it is optional.

You can also override the config file explicitly:

```bash
$ ./riemann1d -c /path/to/riemann1d.cfg
```

## Reference
- Amano, T. 2016. A second-order divergence-constrained multidimensional numerical scheme for relativistic two-fluid electrodynamics. The Astrophysical Journal, 831(1), 100. https://doi.org/10.3847/0004-637X/831/1/100
