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
private:
  blitz::Range II;

public:
  SymmetricBoundary()
  {
    II = blitz::Range::all();
  }

  void set_field_x(Global &g, T_vector &eb, const int Nb)
  {
    if( mpiutil::isLowerBoundary(2) ) {
      // normal
      for(int ibc=1; ibc < Nb ;ibc++) {
        eb(II,II,g.Lbx-1-ibc,0) = eb(II,II,g.Lbx+ibc-1,0);
        eb(II,II,g.Lbx-1-ibc,3) = eb(II,II,g.Lbx+ibc-1,3);
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
        eb(II,II,g.Ubx+1+ibc,0) = eb(II,II,g.Ubx-ibc-1,0);
        eb(II,II,g.Ubx+1+ibc,3) = eb(II,II,g.Ubx-ibc-1,3);
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
        eb(II,g.Lby-1-ibc,II,1) = eb(II,g.Lby+ibc-1,II,1);
        eb(II,g.Lby-1-ibc,II,4) = eb(II,g.Lby+ibc-1,II,4);
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
        eb(II,g.Uby+1+ibc,II,1) = eb(II,g.Uby-ibc-1,II,1);
        eb(II,g.Uby+1+ibc,II,4) = eb(II,g.Uby-ibc-1,II,4);
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
        eb(g.Lbz-1-ibc,II,II,2) = eb(g.Lbz+ibc-1,II,II,2);
        eb(g.Lbz-1-ibc,II,II,5) = eb(g.Lbz+ibc-1,II,II,5);
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
        eb(g.Ubz+1+ibc,II,II,2) = eb(g.Ubz-ibc-1,II,II,2);
        eb(g.Ubz+1+ibc,II,II,5) = eb(g.Ubz-ibc-1,II,II,5);
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
      set_lower_x(g, uf, II, Nb);
      set_lower_x(g, eb, II, Nb);
    }

    if( mpiutil::isUpperBoundary(2) ) {
      set_upper_x(g, uf, II, Nb);
      set_upper_x(g, eb, II, Nb);
    }
  }

  void set_fluid_y(Global &g, T_vector &uf, T_vector &eb, const int Nb)
  {
    if( mpiutil::isLowerBoundary(1) ) {
      set_lower_y(g, uf, II, Nb);
      set_lower_y(g, eb, II, Nb);
    }

    if( mpiutil::isUpperBoundary(1) ) {
      set_upper_y(g, uf, II, Nb);
      set_upper_y(g, eb, II, Nb);
    }
  }

  void set_fluid_z(Global &g, T_vector &uf, T_vector &eb, const int Nb)
  {
    if( mpiutil::isLowerBoundary(0) ) {
      set_lower_z(g, uf, II, Nb);
      set_lower_z(g, eb, II, Nb);
    }

    if( mpiutil::isUpperBoundary(0) ) {
      set_upper_z(g, uf, II, Nb);
      set_upper_z(g, eb, II, Nb);
    }
  }

  template <class T_range>
  void set_lower_x(Global &g, T_vector &x, T_range Ic, const int Nb)
  {
    blitz::Range R1(g.Lbx-Nb, g.Lbx-1);
    blitz::Range R2(g.Lbx+Nb-1, g.Lbx);
    x(II,II,R1,Ic) = x(II,II,R2,Ic);
  }

  template <class T_range>
  void set_lower_y(Global &g, T_vector &x, T_range Ic, const int Nb)
  {
    blitz::Range R1(g.Lby-Nb, g.Lby-1);
    blitz::Range R2(g.Lby+Nb-1, g.Lby);
    x(II,R1,II,Ic) = x(II,R2,II,Ic);
  }

  template <class T_range>
  void set_lower_z(Global &g, T_vector &x, T_range Ic, const int Nb)
  {
    blitz::Range R1(g.Lbz-Nb, g.Lbz-1);
    blitz::Range R2(g.Lbz+Nb-1, g.Lbz);
    x(R1,II,II,Ic) = x(R2,II,II,Ic);
  }

  template <class T_range>
  void set_upper_x(Global &g, T_vector &x, T_range Ic, const int Nb)
  {
    blitz::Range R1(g.Ubx+Nb, g.Ubx+1, -1);
    blitz::Range R2(g.Ubx-Nb+1, g.Ubx);
    x(II,II,R1,Ic) = x(II,II,R2,Ic);
  }

  template <class T_range>
  void set_upper_y(Global &g, T_vector &x, T_range Ic, const int Nb)
  {
    blitz::Range R1(g.Uby+Nb, g.Uby+1, -1);
    blitz::Range R2(g.Uby-Nb+1, g.Uby);
    x(II,R1,II,Ic) = x(II,R2,II,Ic);
  }

  template <class T_range>
  void set_upper_z(Global &g, T_vector &x, T_range Ic, const int Nb)
  {
    blitz::Range R1(g.Ubz+Nb, g.Ubz+1, -1);
    blitz::Range R2(g.Ubz-Nb+1, g.Ubz);
    x(R1,II,II,Ic) = x(R2,II,II,Ic);
  }
};

// Local Variables:
// c-file-style   : "gnu"
// c-file-offsets : ((innamespace . 0) (inline-open . 0))
// End:
#endif
