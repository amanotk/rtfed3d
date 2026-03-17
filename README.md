# Relativistic Two-Fluid Electrodynamics 3D

## Requirements

- Parallel HDF5 (compiled with `--enable-parallel`)
- Blitz++
  - https://github.com/blitzpp/blitz/releases
- MPI C++ compiler wrapper
- CMake 3.20 or newer

## Repository Layout

- `include/`: project headers
- `src/`: shared solver implementation
- `examples/`: runnable problem drivers and example configuration files
- `base/`: vendored support headers used by the solver and examples

## Build

### Clone
```bash
$ git clone git@github.com:amanotk/rtfed3d.git
$ cd rtfed3d
```

### Configure
If Blitz++ or HDF5 live in non-standard locations, pass `-DBLITZ_INCLUDE_DIR=...`,
`-DBLITZ_LIBRARY=...`, and the usual CMake hints for MPI/HDF5 during configure.

```bash
$ cmake --preset release
```

### Compile
```bash
$ cmake --build --preset release
```

Example executables are written to `build/release/examples/`.

### Run
Each example expects its `.cfg` file in the current working directory. The build
copies matching config files next to each executable, so the simplest workflow is:

```bash
$ cd build/release/examples
$ ./riemann1d
```

You can also override the config file explicitly:

```bash
$ ./riemann1d -c /path/to/riemann1d.cfg
```

## Reference
- Amano, T. 2016. A second-order divergence-constrained multidimensional numerical scheme for relativistic two-fluid electrodynamics. The Astrophysical Journal, 831(1), 100. https://doi.org/10.3847/0004-637X/831/1/100
