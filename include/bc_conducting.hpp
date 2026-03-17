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
private:
  blitz::Range II;

public:
  ConductingBoundary()
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
      for(int ibc=0; ibc < Nb ;ibc++) {
        eb(II,II,g.Lbx-1-ibc,1) = eb(II,II,g.Lbx+ibc  ,1);
        eb(II,II,g.Lbx-1-ibc,2) = eb(II,II,g.Lbx+ibc  ,2);
        eb(II,II,g.Lbx-1-ibc,4) = eb(II,II,g.Lbx+ibc  ,4);
        eb(II,II,g.Lbx-1-ibc,5) = eb(II,II,g.Lbx+ibc  ,5);
      }
    }

    if( mpiutil::isUpperBoundary(2) ) {
      // normal
      for(int ibc=0; ibc < Nb ;ibc++) {
        eb(II,II,g.Ubx+1+ibc,0) = eb(II,II,g.Ubx-ibc-1,0);
        eb(II,II,g.Ubx+1+ibc,3) = eb(II,II,g.Ubx-ibc-1,3);
      }
      // transverse
      for(int ibc=0; ibc < Nb ;ibc++) {
        eb(II,II,g.Ubx+1+ibc,1) = eb(II,II,g.Ubx-ibc  ,1);
        eb(II,II,g.Ubx+1+ibc,2) = eb(II,II,g.Ubx-ibc  ,2);
        eb(II,II,g.Ubx+1+ibc,4) = eb(II,II,g.Ubx-ibc  ,4);
        eb(II,II,g.Ubx+1+ibc,5) = eb(II,II,g.Ubx-ibc  ,5);
      }
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
      for(int ibc=0; ibc < Nb ;ibc++) {
        eb(II,g.Lby-1-ibc,II,2) = eb(II,g.Lby+ibc  ,II,2);
        eb(II,g.Lby-1-ibc,II,0) = eb(II,g.Lby+ibc  ,II,0);
        eb(II,g.Lby-1-ibc,II,5) = eb(II,g.Lby+ibc  ,II,5);
        eb(II,g.Lby-1-ibc,II,3) = eb(II,g.Lby+ibc  ,II,3);
      }
    }

    if( mpiutil::isUpperBoundary(1) ) {
      // normal
      for(int ibc=0; ibc < Nb ;ibc++) {
        eb(II,g.Uby+1+ibc,II,1) = eb(II,g.Uby-ibc-1,II,1);
        eb(II,g.Uby+1+ibc,II,4) = eb(II,g.Uby-ibc-1,II,4);
      }
      // transverse
      for(int ibc=0; ibc < Nb ;ibc++) {
        eb(II,g.Uby+1+ibc,II,2) = eb(II,g.Uby-ibc  ,II,2);
        eb(II,g.Uby+1+ibc,II,0) = eb(II,g.Uby-ibc  ,II,0);
        eb(II,g.Uby+1+ibc,II,5) = eb(II,g.Uby-ibc  ,II,5);
        eb(II,g.Uby+1+ibc,II,3) = eb(II,g.Uby-ibc  ,II,3);
      }
    }
  }

  void set_field_z(Global &g, T_vector &eb, const int Nb)
  {
    if( mpiutil::isLowerBoundary(0) ) {
      // normal
      for(int ibc=0; ibc < Nb ;ibc++) {
        eb(g.Lbz-1-ibc,II,II,2) = eb(g.Lbz+ibc-1,II,II,2);
        eb(g.Lbz-1-ibc,II,II,5) = eb(g.Lbz+ibc-1,II,II,5);
      }
      // transverse
      for(int ibc=0; ibc < Nb ;ibc++) {
        eb(g.Lbz-1-ibc,II,II,0) = eb(g.Lbz+ibc  ,II,II,0);
        eb(g.Lbz-1-ibc,II,II,1) = eb(g.Lbz+ibc  ,II,II,1);
        eb(g.Lbz-1-ibc,II,II,3) = eb(g.Lbz+ibc  ,II,II,3);
        eb(g.Lbz-1-ibc,II,II,4) = eb(g.Lbz+ibc  ,II,II,4);
      }
    }

    if( mpiutil::isUpperBoundary(0) ) {
      // normal
      for(int ibc=0; ibc < Nb ;ibc++) {
        eb(g.Ubz+1+ibc,II,II,2) = eb(g.Ubz-ibc-1,II,II,2);
        eb(g.Ubz+1+ibc,II,II,5) = eb(g.Ubz-ibc-1,II,II,5);
      }
      // transverse
      for(int ibc=0; ibc < Nb ;ibc++) {
        eb(g.Ubz+1+ibc,II,II,0) = eb(g.Ubz-ibc  ,II,II,0);
        eb(g.Ubz+1+ibc,II,II,1) = eb(g.Ubz-ibc  ,II,II,1);
        eb(g.Ubz+1+ibc,II,II,3) = eb(g.Ubz-ibc  ,II,II,3);
        eb(g.Ubz+1+ibc,II,II,4) = eb(g.Ubz-ibc  ,II,II,4);
      }
    }
  }

  void set_fluid_x(Global &g, T_vector &uf, T_vector &eb, const int Nb)
  {
    const static float64 A[10] = { +1.0, -1.0, +1.0, +1.0, +1.0,
                                   +1.0, -1.0, +1.0, +1.0, +1.0};
    const static float64 B[6]  = { -1.0, -1.0, -1.0, +1.0, +1.0, +1.0};

    if( mpiutil::isLowerBoundary(2) ) {
      for(int ibc=0; ibc < Nb ;ibc++) {
        // fluid
        for(int k=0; k < 10 ;k++) {
          uf(II,II,g.Lbx-1-ibc,k) = A[k] * uf(II,II,g.Lbx+ibc,k);
        }
        // field
        for(int k=0; k < 6 ;k++) {
          eb(II,II,g.Lbx-1-ibc,k) = B[k] * eb(II,II,g.Lbx+ibc,k);
        }
      }
    }

    if( mpiutil::isUpperBoundary(2) ) {
      for(int ibc=0; ibc < Nb ;ibc++) {
        // fluid
        for(int k=0; k < 10 ;k++) {
          uf(II,II,g.Ubx+1+ibc,k) = A[k] * uf(II,II,g.Ubx-ibc,k);
        }
        // field
        for(int k=0; k < 6 ;k++) {
          eb(II,II,g.Ubx+1+ibc,k) = B[k] * eb(II,II,g.Ubx-ibc,k);
        }
      }
    }
  }

  void set_fluid_y(Global &g, T_vector &uf, T_vector &eb, const int Nb)
  {
    const static float64 A[10] = { +1.0, +1.0, -1.0, +1.0, +1.0,
                                   +1.0, +1.0, -1.0, +1.0, +1.0};
    const static float64 B[6]  = { -1.0, -1.0, -1.0, +1.0, +1.0, +1.0};

    if( mpiutil::isLowerBoundary(1) ) {
      for(int ibc=0; ibc < Nb ;ibc++) {
        // fluid
        for(int k=0; k < 10 ;k++) {
          uf(II,g.Lby-1-ibc,II,k) = A[k] * uf(II,g.Lby+ibc,II,k);
        }
        // field
        for(int k=0; k < 6 ;k++) {
          eb(II,g.Lby-1-ibc,II,k) = B[k] * eb(II,g.Lby+ibc,II,k);
        }
      }
    }

    if( mpiutil::isUpperBoundary(1) ) {
      for(int ibc=0; ibc < Nb ;ibc++) {
        // fluid
        for(int k=0; k < 10 ;k++) {
          uf(II,g.Uby+1+ibc,II,k) = A[k] * uf(II,g.Uby-ibc,II,k);
        }
        // field
        for(int k=0; k < 6 ;k++) {
          eb(II,g.Uby+1+ibc,II,k) = B[k] * eb(II,g.Uby-ibc,II,k);
        }
      }
    }
  }

  void set_fluid_z(Global &g, T_vector &uf, T_vector &eb, const int Nb)
  {
    const static float64 A[10] = { +1.0, +1.0, +1.0, -1.0, +1.0,
                                   +1.0, +1.0, +1.0, -1.0, +1.0};
    const static float64 B[6]  = { -1.0, -1.0, -1.0, +1.0, +1.0, +1.0};

    if( mpiutil::isLowerBoundary(0) ) {
      for(int ibc=0; ibc < Nb ;ibc++) {
        // fluid
        for(int k=0; k < 10 ;k++) {
          uf(g.Lbz-1-ibc,II,II,k) = A[k] * uf(g.Lbz+ibc,II,II,k);
        }
        // field
        for(int k=0; k < 6 ;k++) {
          eb(g.Lbz-1-ibc,II,II,k) = B[k] * eb(g.Lbz+ibc,II,II,k);
        }
      }
    }

    if( mpiutil::isUpperBoundary(0) ) {
      for(int ibc=0; ibc < Nb ;ibc++) {
        // fluid
        for(int k=0; k < 10 ;k++) {
          uf(g.Ubz+1+ibc,II,II,k) = A[k] * uf(g.Ubz-ibc,II,II,k);
        }
        // field
        for(int k=0; k < 6 ;k++) {
          eb(g.Ubz+1+ibc,II,II,k) = B[k] * eb(g.Ubz-ibc,II,II,k);
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
