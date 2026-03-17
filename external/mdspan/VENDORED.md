# Vendored mdspan

- Upstream: `https://github.com/kokkos/mdspan`
- Upstream branch: `stable`
- Upstream revision at vendoring time: `80fc772eb812b45097c28fc0a46d8ff006138d69`

This directory contains the full upstream source tree as vendored project code.
Local project code should include headers from `external/mdspan/include/` through
the main project CMake configuration rather than building the upstream mdspan
project itself.
