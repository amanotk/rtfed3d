// -*- C++ -*-
#ifndef _BC_CONDUCTING_HPP_
#define _BC_CONDUCTING_HPP_

///
/// @brief Conducting Wall Boundary Condition
///
/// $Id: bc_conducting.hpp,v 99bd1b55086c 2015/09/02 19:54:43 amano $
///
#include "global.hpp"
#include "mpiutil.hpp"

///
/// Implementation of conducting wall boundary condition
///
class ConductingBoundary
{
public:
  ConductingBoundary()
  {
  }

  void set_field_x(Global& g, T_vector& eb, const int Nb)
  {
    if (mpiutil::isLowerBoundary(2)) {
      // normal
      for (int ibc = 1; ibc < Nb; ibc++) {
        for (int iz = 0; iz < g.Mz; iz++) {
          for (int iy = 0; iy < g.My; iy++) {
            eb(iz, iy, g.Lbx - 1 - ibc, 0) = eb(iz, iy, g.Lbx + ibc - 1, 0);
            eb(iz, iy, g.Lbx - 1 - ibc, 3) = eb(iz, iy, g.Lbx + ibc - 1, 3);
          }
        }
      }
      // transverse
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iz = 0; iz < g.Mz; iz++) {
          for (int iy = 0; iy < g.My; iy++) {
            eb(iz, iy, g.Lbx - 1 - ibc, 1) = eb(iz, iy, g.Lbx + ibc, 1);
            eb(iz, iy, g.Lbx - 1 - ibc, 2) = eb(iz, iy, g.Lbx + ibc, 2);
            eb(iz, iy, g.Lbx - 1 - ibc, 4) = eb(iz, iy, g.Lbx + ibc, 4);
            eb(iz, iy, g.Lbx - 1 - ibc, 5) = eb(iz, iy, g.Lbx + ibc, 5);
          }
        }
      }
    }

    if (mpiutil::isUpperBoundary(2)) {
      // normal
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iz = 0; iz < g.Mz; iz++) {
          for (int iy = 0; iy < g.My; iy++) {
            eb(iz, iy, g.Ubx + 1 + ibc, 0) = eb(iz, iy, g.Ubx - ibc - 1, 0);
            eb(iz, iy, g.Ubx + 1 + ibc, 3) = eb(iz, iy, g.Ubx - ibc - 1, 3);
          }
        }
      }
      // transverse
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iz = 0; iz < g.Mz; iz++) {
          for (int iy = 0; iy < g.My; iy++) {
            eb(iz, iy, g.Ubx + 1 + ibc, 1) = eb(iz, iy, g.Ubx - ibc, 1);
            eb(iz, iy, g.Ubx + 1 + ibc, 2) = eb(iz, iy, g.Ubx - ibc, 2);
            eb(iz, iy, g.Ubx + 1 + ibc, 4) = eb(iz, iy, g.Ubx - ibc, 4);
            eb(iz, iy, g.Ubx + 1 + ibc, 5) = eb(iz, iy, g.Ubx - ibc, 5);
          }
        }
      }
    }
  }

  void set_field_y(Global& g, T_vector& eb, const int Nb)
  {
    if (mpiutil::isLowerBoundary(1)) {
      // normal
      for (int ibc = 1; ibc < Nb; ibc++) {
        for (int iz = 0; iz < g.Mz; iz++) {
          for (int ix = 0; ix < g.Mx; ix++) {
            eb(iz, g.Lby - 1 - ibc, ix, 1) = eb(iz, g.Lby + ibc - 1, ix, 1);
            eb(iz, g.Lby - 1 - ibc, ix, 4) = eb(iz, g.Lby + ibc - 1, ix, 4);
          }
        }
      }
      // transverse
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iz = 0; iz < g.Mz; iz++) {
          for (int ix = 0; ix < g.Mx; ix++) {
            eb(iz, g.Lby - 1 - ibc, ix, 2) = eb(iz, g.Lby + ibc, ix, 2);
            eb(iz, g.Lby - 1 - ibc, ix, 0) = eb(iz, g.Lby + ibc, ix, 0);
            eb(iz, g.Lby - 1 - ibc, ix, 5) = eb(iz, g.Lby + ibc, ix, 5);
            eb(iz, g.Lby - 1 - ibc, ix, 3) = eb(iz, g.Lby + ibc, ix, 3);
          }
        }
      }
    }

    if (mpiutil::isUpperBoundary(1)) {
      // normal
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iz = 0; iz < g.Mz; iz++) {
          for (int ix = 0; ix < g.Mx; ix++) {
            eb(iz, g.Uby + 1 + ibc, ix, 1) = eb(iz, g.Uby - ibc - 1, ix, 1);
            eb(iz, g.Uby + 1 + ibc, ix, 4) = eb(iz, g.Uby - ibc - 1, ix, 4);
          }
        }
      }
      // transverse
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iz = 0; iz < g.Mz; iz++) {
          for (int ix = 0; ix < g.Mx; ix++) {
            eb(iz, g.Uby + 1 + ibc, ix, 2) = eb(iz, g.Uby - ibc, ix, 2);
            eb(iz, g.Uby + 1 + ibc, ix, 0) = eb(iz, g.Uby - ibc, ix, 0);
            eb(iz, g.Uby + 1 + ibc, ix, 5) = eb(iz, g.Uby - ibc, ix, 5);
            eb(iz, g.Uby + 1 + ibc, ix, 3) = eb(iz, g.Uby - ibc, ix, 3);
          }
        }
      }
    }
  }

  void set_field_z(Global& g, T_vector& eb, const int Nb)
  {
    if (mpiutil::isLowerBoundary(0)) {
      // normal
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iy = 0; iy < g.My; iy++) {
          for (int ix = 0; ix < g.Mx; ix++) {
            eb(g.Lbz - 1 - ibc, iy, ix, 2) = eb(g.Lbz + ibc - 1, iy, ix, 2);
            eb(g.Lbz - 1 - ibc, iy, ix, 5) = eb(g.Lbz + ibc - 1, iy, ix, 5);
          }
        }
      }
      // transverse
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iy = 0; iy < g.My; iy++) {
          for (int ix = 0; ix < g.Mx; ix++) {
            eb(g.Lbz - 1 - ibc, iy, ix, 0) = eb(g.Lbz + ibc, iy, ix, 0);
            eb(g.Lbz - 1 - ibc, iy, ix, 1) = eb(g.Lbz + ibc, iy, ix, 1);
            eb(g.Lbz - 1 - ibc, iy, ix, 3) = eb(g.Lbz + ibc, iy, ix, 3);
            eb(g.Lbz - 1 - ibc, iy, ix, 4) = eb(g.Lbz + ibc, iy, ix, 4);
          }
        }
      }
    }

    if (mpiutil::isUpperBoundary(0)) {
      // normal
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iy = 0; iy < g.My; iy++) {
          for (int ix = 0; ix < g.Mx; ix++) {
            eb(g.Ubz + 1 + ibc, iy, ix, 2) = eb(g.Ubz - ibc - 1, iy, ix, 2);
            eb(g.Ubz + 1 + ibc, iy, ix, 5) = eb(g.Ubz - ibc - 1, iy, ix, 5);
          }
        }
      }
      // transverse
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iy = 0; iy < g.My; iy++) {
          for (int ix = 0; ix < g.Mx; ix++) {
            eb(g.Ubz + 1 + ibc, iy, ix, 0) = eb(g.Ubz - ibc, iy, ix, 0);
            eb(g.Ubz + 1 + ibc, iy, ix, 1) = eb(g.Ubz - ibc, iy, ix, 1);
            eb(g.Ubz + 1 + ibc, iy, ix, 3) = eb(g.Ubz - ibc, iy, ix, 3);
            eb(g.Ubz + 1 + ibc, iy, ix, 4) = eb(g.Ubz - ibc, iy, ix, 4);
          }
        }
      }
    }
  }

  void set_fluid_x(Global& g, T_vector& uf, T_vector& eb, const int Nb)
  {
    const static float64 A[10] = {+1.0, -1.0, +1.0, +1.0, +1.0, +1.0, -1.0, +1.0, +1.0, +1.0};
    const static float64 B[6]  = {-1.0, -1.0, -1.0, +1.0, +1.0, +1.0};

    if (mpiutil::isLowerBoundary(2)) {
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iz = 0; iz < g.Mz; iz++) {
          for (int iy = 0; iy < g.My; iy++) {
            for (int k = 0; k < 10; k++) {
              uf(iz, iy, g.Lbx - 1 - ibc, k) = A[k] * uf(iz, iy, g.Lbx + ibc, k);
            }
            for (int k = 0; k < 6; k++) {
              eb(iz, iy, g.Lbx - 1 - ibc, k) = B[k] * eb(iz, iy, g.Lbx + ibc, k);
            }
          }
        }
      }
    }

    if (mpiutil::isUpperBoundary(2)) {
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iz = 0; iz < g.Mz; iz++) {
          for (int iy = 0; iy < g.My; iy++) {
            for (int k = 0; k < 10; k++) {
              uf(iz, iy, g.Ubx + 1 + ibc, k) = A[k] * uf(iz, iy, g.Ubx - ibc, k);
            }
            for (int k = 0; k < 6; k++) {
              eb(iz, iy, g.Ubx + 1 + ibc, k) = B[k] * eb(iz, iy, g.Ubx - ibc, k);
            }
          }
        }
      }
    }
  }

  void set_fluid_y(Global& g, T_vector& uf, T_vector& eb, const int Nb)
  {
    const static float64 A[10] = {+1.0, +1.0, -1.0, +1.0, +1.0, +1.0, +1.0, -1.0, +1.0, +1.0};
    const static float64 B[6]  = {-1.0, -1.0, -1.0, +1.0, +1.0, +1.0};

    if (mpiutil::isLowerBoundary(1)) {
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iz = 0; iz < g.Mz; iz++) {
          for (int ix = 0; ix < g.Mx; ix++) {
            for (int k = 0; k < 10; k++) {
              uf(iz, g.Lby - 1 - ibc, ix, k) = A[k] * uf(iz, g.Lby + ibc, ix, k);
            }
            for (int k = 0; k < 6; k++) {
              eb(iz, g.Lby - 1 - ibc, ix, k) = B[k] * eb(iz, g.Lby + ibc, ix, k);
            }
          }
        }
      }
    }

    if (mpiutil::isUpperBoundary(1)) {
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iz = 0; iz < g.Mz; iz++) {
          for (int ix = 0; ix < g.Mx; ix++) {
            for (int k = 0; k < 10; k++) {
              uf(iz, g.Uby + 1 + ibc, ix, k) = A[k] * uf(iz, g.Uby - ibc, ix, k);
            }
            for (int k = 0; k < 6; k++) {
              eb(iz, g.Uby + 1 + ibc, ix, k) = B[k] * eb(iz, g.Uby - ibc, ix, k);
            }
          }
        }
      }
    }
  }

  void set_fluid_z(Global& g, T_vector& uf, T_vector& eb, const int Nb)
  {
    const static float64 A[10] = {+1.0, +1.0, +1.0, -1.0, +1.0, +1.0, +1.0, +1.0, -1.0, +1.0};
    const static float64 B[6]  = {-1.0, -1.0, -1.0, +1.0, +1.0, +1.0};

    if (mpiutil::isLowerBoundary(0)) {
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iy = 0; iy < g.My; iy++) {
          for (int ix = 0; ix < g.Mx; ix++) {
            for (int k = 0; k < 10; k++) {
              uf(g.Lbz - 1 - ibc, iy, ix, k) = A[k] * uf(g.Lbz + ibc, iy, ix, k);
            }
            for (int k = 0; k < 6; k++) {
              eb(g.Lbz - 1 - ibc, iy, ix, k) = B[k] * eb(g.Lbz + ibc, iy, ix, k);
            }
          }
        }
      }
    }

    if (mpiutil::isUpperBoundary(0)) {
      for (int ibc = 0; ibc < Nb; ibc++) {
        for (int iy = 0; iy < g.My; iy++) {
          for (int ix = 0; ix < g.Mx; ix++) {
            for (int k = 0; k < 10; k++) {
              uf(g.Ubz + 1 + ibc, iy, ix, k) = A[k] * uf(g.Ubz - ibc, iy, ix, k);
            }
            for (int k = 0; k < 6; k++) {
              eb(g.Ubz + 1 + ibc, iy, ix, k) = B[k] * eb(g.Ubz - ibc, iy, ix, k);
            }
          }
        }
      }
    }
  }
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
