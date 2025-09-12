# Relativistic Two-Fluid Electrodynamics 3D

## Requirements

- Parallel HDF5 (compiled with `--enable-parallel`)
- blitz
  - https://github.com/blitzpp/blitz/releases

## Sample
### Clone
```bash
$ git clone git@github.com:amanotk/rtfed3d.git
$ cd rtfed3d
$ git submodule update -i
```

### Compilation
Copy the default compiler configuration file to `Makefile.inc` and edit it as you like
```
$ cp Makefile-MPI.inc Makefile.inc
```
You may need to specify the include and library directories for blitz and hdf5 if they are installed in non-standard directories.

## Reference
- Amano, T. 2016. A second-order divergence-constrained multidimensional numerical scheme for relativistic two-fluid electrodynamics. The Astrophysical Journal, 831(1), 100. https://doi.org/10.3847/0004-637X/831/1/100
