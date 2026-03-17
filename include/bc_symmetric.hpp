// -*- C++ -*-
#ifndef _BC_SYMMETRIC_HPP_
#define _BC_SYMMETRIC_HPP_

///
/// @brief Symmetric Boundary Condition
///
/// $Id: bc_symmetric.hpp,v 99bd1b55086c 2015/09/02 19:54:43 amano $
///
#include "global.hpp"
#include "mpiutil.hpp"

///
/// Implementation of symmetric boundary condition
///
class SymmetricBoundary
{
public:
  SymmetricBoundary() {}

  void set_field_x(Global &g, T_vector &eb, const int Nb)
  {
    if( mpiutil::isLowerBoundary(2) ) {
      // normal
      for(int ibc=1; ibc < Nb ;ibc++) {
        for(int iz=0; iz < g.Mz ;iz++) {
          for(int iy=0; iy < g.My ;iy++) {
            eb(iz,iy,g.Lbx-1-ibc,0) = eb(iz,iy,g.Lbx+ibc-1,0);
            eb(iz,iy,g.Lbx-1-ibc,3) = eb(iz,iy,g.Lbx+ibc-1,3);
          }
        }
      }
      // transverse
      set_lower_x(g, eb, 1, Nb);
      set_lower_x(g, eb, 2, Nb);
      set_lower_x(g, eb, 4, Nb);
      set_lower_x(g, eb, 5, Nb);
    }

    if( mpiutil::isUpperBoundary(2) ) {
      // normal
      for(int ibc=0; ibc < Nb ;ibc++) {
        for(int iz=0; iz < g.Mz ;iz++) {
          for(int iy=0; iy < g.My ;iy++) {
            eb(iz,iy,g.Ubx+1+ibc,0) = eb(iz,iy,g.Ubx-ibc-1,0);
            eb(iz,iy,g.Ubx+1+ibc,3) = eb(iz,iy,g.Ubx-ibc-1,3);
          }
        }
      }
      // transverse
      set_upper_x(g, eb, 1, Nb);
      set_upper_x(g, eb, 2, Nb);
      set_upper_x(g, eb, 4, Nb);
      set_upper_x(g, eb, 5, Nb);
    }
  }

  void set_field_y(Global &g, T_vector &eb, const int Nb)
  {
    if( mpiutil::isLowerBoundary(1) ) {
      // normal
      for(int ibc=1; ibc < Nb ;ibc++) {
        for(int iz=0; iz < g.Mz ;iz++) {
          for(int ix=0; ix < g.Mx ;ix++) {
            eb(iz,g.Lby-1-ibc,ix,1) = eb(iz,g.Lby+ibc-1,ix,1);
            eb(iz,g.Lby-1-ibc,ix,4) = eb(iz,g.Lby+ibc-1,ix,4);
          }
        }
      }
      // transverse
      set_lower_y(g, eb, 2, Nb);
      set_lower_y(g, eb, 0, Nb);
      set_lower_y(g, eb, 5, Nb);
      set_lower_y(g, eb, 3, Nb);
    }

    if( mpiutil::isUpperBoundary(1) ) {
      // normal
      for(int ibc=0; ibc < Nb ;ibc++) {
        for(int iz=0; iz < g.Mz ;iz++) {
          for(int ix=0; ix < g.Mx ;ix++) {
            eb(iz,g.Uby+1+ibc,ix,1) = eb(iz,g.Uby-ibc-1,ix,1);
            eb(iz,g.Uby+1+ibc,ix,4) = eb(iz,g.Uby-ibc-1,ix,4);
          }
        }
      }
      // transverse
      set_upper_y(g, eb, 2, Nb);
      set_upper_y(g, eb, 0, Nb);
      set_upper_y(g, eb, 5, Nb);
      set_upper_y(g, eb, 3, Nb);
    }
  }

  void set_field_z(Global &g, T_vector &eb, const int Nb)
  {
    if( mpiutil::isLowerBoundary(0) ) {
      // normal
      for(int ibc=1; ibc < Nb ;ibc++) {
        for(int iy=0; iy < g.My ;iy++) {
          for(int ix=0; ix < g.Mx ;ix++) {
            eb(g.Lbz-1-ibc,iy,ix,2) = eb(g.Lbz+ibc-1,iy,ix,2);
            eb(g.Lbz-1-ibc,iy,ix,5) = eb(g.Lbz+ibc-1,iy,ix,5);
          }
        }
      }
      // transverse
      set_lower_z(g, eb, 0, Nb);
      set_lower_z(g, eb, 1, Nb);
      set_lower_z(g, eb, 3, Nb);
      set_lower_z(g, eb, 4, Nb);
    }

    if( mpiutil::isUpperBoundary(0) ) {
      // normal
      for(int ibc=0; ibc < Nb ;ibc++) {
        for(int iy=0; iy < g.My ;iy++) {
          for(int ix=0; ix < g.Mx ;ix++) {
            eb(g.Ubz+1+ibc,iy,ix,2) = eb(g.Ubz-ibc-1,iy,ix,2);
            eb(g.Ubz+1+ibc,iy,ix,5) = eb(g.Ubz-ibc-1,iy,ix,5);
          }
        }
      }
      // transverse
      set_upper_z(g, eb, 0, Nb);
      set_upper_z(g, eb, 1, Nb);
      set_upper_z(g, eb, 3, Nb);
      set_upper_z(g, eb, 4, Nb);
    }
  }

  void set_fluid_x(Global &g, T_vector &uf, T_vector &eb, const int Nb)
  {
    if( mpiutil::isLowerBoundary(2) ) {
      set_lower_x(g, uf, Nb);
      set_lower_x(g, eb, Nb);
    }

    if( mpiutil::isUpperBoundary(2) ) {
      set_upper_x(g, uf, Nb);
      set_upper_x(g, eb, Nb);
    }
  }

  void set_fluid_y(Global &g, T_vector &uf, T_vector &eb, const int Nb)
  {
    if( mpiutil::isLowerBoundary(1) ) {
      set_lower_y(g, uf, Nb);
      set_lower_y(g, eb, Nb);
    }

    if( mpiutil::isUpperBoundary(1) ) {
      set_upper_y(g, uf, Nb);
      set_upper_y(g, eb, Nb);
    }
  }

  void set_fluid_z(Global &g, T_vector &uf, T_vector &eb, const int Nb)
  {
    if( mpiutil::isLowerBoundary(0) ) {
      set_lower_z(g, uf, Nb);
      set_lower_z(g, eb, Nb);
    }

    if( mpiutil::isUpperBoundary(0) ) {
      set_upper_z(g, uf, Nb);
      set_upper_z(g, eb, Nb);
    }
  }

  void set_lower_x(Global &g, T_vector &x, int Ic, const int Nb)
  {
    for(int ix=0; ix < Nb ;ix++) {
      for(int iz=0; iz < g.Mz ;iz++) {
        for(int iy=0; iy < g.My ;iy++) {
          x(iz,iy,g.Lbx-Nb+ix,Ic) = x(iz,iy,g.Lbx+Nb-1-ix,Ic);
        }
      }
    }
  }

  void set_lower_y(Global &g, T_vector &x, int Ic, const int Nb)
  {
    for(int iy=0; iy < Nb ;iy++) {
      for(int iz=0; iz < g.Mz ;iz++) {
        for(int ix=0; ix < g.Mx ;ix++) {
          x(iz,g.Lby-Nb+iy,ix,Ic) = x(iz,g.Lby+Nb-1-iy,ix,Ic);
        }
      }
    }
  }

  void set_lower_z(Global &g, T_vector &x, int Ic, const int Nb)
  {
    for(int iz=0; iz < Nb ;iz++) {
      for(int iy=0; iy < g.My ;iy++) {
        for(int ix=0; ix < g.Mx ;ix++) {
          x(g.Lbz-Nb+iz,iy,ix,Ic) = x(g.Lbz+Nb-1-iz,iy,ix,Ic);
        }
      }
    }
  }

  void set_upper_x(Global &g, T_vector &x, int Ic, const int Nb)
  {
    for(int ix=0; ix < Nb ;ix++) {
      for(int iz=0; iz < g.Mz ;iz++) {
        for(int iy=0; iy < g.My ;iy++) {
          x(iz,iy,g.Ubx+Nb-ix,Ic) = x(iz,iy,g.Ubx-Nb+1+ix,Ic);
        }
      }
    }
  }

  void set_upper_y(Global &g, T_vector &x, int Ic, const int Nb)
  {
    for(int iy=0; iy < Nb ;iy++) {
      for(int iz=0; iz < g.Mz ;iz++) {
        for(int ix=0; ix < g.Mx ;ix++) {
          x(iz,g.Uby+Nb-iy,ix,Ic) = x(iz,g.Uby-Nb+1+iy,ix,Ic);
        }
      }
    }
  }

  void set_upper_z(Global &g, T_vector &x, int Ic, const int Nb)
  {
    for(int iz=0; iz < Nb ;iz++) {
      for(int iy=0; iy < g.My ;iy++) {
        for(int ix=0; ix < g.Mx ;ix++) {
          x(g.Ubz+Nb-iz,iy,ix,Ic) = x(g.Ubz-Nb+1+iz,iy,ix,Ic);
        }
      }
    }
  }

  void set_lower_x(Global &g, T_vector &x, const int Nb)
  {
    for(int ic=0; ic < x.extent(3) ;ic++) {
      set_lower_x(g, x, ic, Nb);
    }
  }

  void set_lower_y(Global &g, T_vector &x, const int Nb)
  {
    for(int ic=0; ic < x.extent(3) ;ic++) {
      set_lower_y(g, x, ic, Nb);
    }
  }

  void set_lower_z(Global &g, T_vector &x, const int Nb)
  {
    for(int ic=0; ic < x.extent(3) ;ic++) {
      set_lower_z(g, x, ic, Nb);
    }
  }

  void set_upper_x(Global &g, T_vector &x, const int Nb)
  {
    for(int ic=0; ic < x.extent(3) ;ic++) {
      set_upper_x(g, x, ic, Nb);
    }
  }

  void set_upper_y(Global &g, T_vector &x, const int Nb)
  {
    for(int ic=0; ic < x.extent(3) ;ic++) {
      set_upper_y(g, x, ic, Nb);
    }
  }

  void set_upper_z(Global &g, T_vector &x, const int Nb)
  {
    for(int ic=0; ic < x.extent(3) ;ic++) {
      set_upper_z(g, x, ic, Nb);
    }
  }
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
